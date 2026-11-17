"""Optical drives: how the absorbed power density p_abs(t) is obtained.

Three levels, in increasing coupling to the electromagnetic problem:

``PowerDrive``      p_abs(t) prescribed.  Used for tests and for controlled B0/B1 comparisons.
``FieldDrive``      the local (in-film) field envelope |E(t)| is prescribed and
                    p_abs = (omega/2) eps0 eps''[f](omega) |E|^2          ([14] Eq. (5), [15] Eq. (6)).
                    This is the WP1 dataset setting ("in-film peak intensity").
``ThinFilmDrive``   the *incident* intensity envelope is prescribed and the local field follows from a
                    transfer-matrix solution of the layer stack with the instantaneous eps[f]
                    (adiabatic self-consistent scheme of [15] Sec. II B).  WP3 replaces this by FDTD.

All drives expose ``power(t, eps)`` returning p_abs in eV/(nm^3 fs), where ``eps`` is the complex
permittivity of the film at the pump frequency evaluated from the current distribution.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

from . import constants as c
from .tmm import layer_average_field_sq


def gaussian_envelope(t, t0: float, fwhm: float):
    """Intensity envelope exp(-4 ln2 ((t-t0)/fwhm)^2), peak 1."""
    return np.exp(-4.0 * math.log(2.0) * ((np.asarray(t, dtype=float) - t0) / fwhm) ** 2)


@dataclass
class Pulse:
    """Gaussian pulse: intensity FWHM ``fwhm`` (fs), centre ``t0`` (fs), vacuum wavelength (nm)."""

    wavelength: float
    fwhm: float
    t0: float = 0.0

    @property
    def omega(self) -> float:
        return 2.0 * math.pi * c.C_LIGHT / self.wavelength

    @property
    def photon_energy(self) -> float:
        return c.HBAR * self.omega

    def envelope(self, t):
        return gaussian_envelope(t, self.t0, self.fwhm)


@dataclass
class PowerDrive:
    pulse: Pulse
    peak_power: float  # eV/(nm^3 fs)
    needs_eps: bool = False

    def power(self, t: float, eps: complex | None = None) -> float:
        return float(self.peak_power * self.pulse.envelope(t))

    @classmethod
    def from_absorbed_energy(cls, pulse: Pulse, u_abs_j_m3: float) -> "PowerDrive":
        """Drive that deposits a total energy density u_abs (J/m^3)."""
        u = u_abs_j_m3 * c.J_PER_M3_IN_EV_PER_NM3
        area = pulse.fwhm * math.sqrt(math.pi / (4.0 * math.log(2.0)))
        return cls(pulse, u / area)


@dataclass
class FieldDrive:
    """Prescribed in-film intensity I_film(t) = I_peak * envelope; |E|^2 = 2 Z0 I / n with n = 1 convention.

    The "in-film intensity" is defined through the local field amplitude, I = |E_loc|^2 / (2 Z0), i.e. the
    intensity a vacuum plane wave of the same field amplitude would carry.  This avoids the ambiguity of
    dividing by a refractive index that vanishes at the ENZ point.
    """

    pulse: Pulse
    peak_intensity_w_cm2: float
    needs_eps: bool = True

    def field_sq(self, t: float) -> float:
        e0 = c.intensity_to_field(self.peak_intensity_w_cm2)
        return float(e0**2 * self.pulse.envelope(t))

    def power(self, t: float, eps: complex) -> float:
        return 0.5 * self.pulse.omega * c.EPS0 * max(eps.imag, 0.0) * self.field_sq(t)

    def local_intensity(self, t: float, eps: complex | None = None) -> float:
        """In-film intensity in W/cm^2 (used by intensity-driven interband terms)."""
        return float(self.peak_intensity_w_cm2 * self.pulse.envelope(t))


@dataclass
class ThinFilmDrive:
    """Self-consistent adiabatic drive: incident plane wave on a layer stack containing the active film.

    ``layers_before`` / ``layers_after`` are sequences of (eps, thickness_nm) on the incidence side and the
    exit side of the active film; ``eps_in`` / ``eps_out`` are the semi-infinite media.  The absorbed power
    is the film average of (omega/2) eps0 eps'' |E|^2 as in [15] (the film is much thinner than the
    wavelength, so [15] uses the spatial average as well).
    """

    pulse: Pulse
    peak_intensity_w_cm2: float  # incident, in the incidence medium
    thickness: float  # active film, nm
    angle_deg: float = 0.0  # angle of incidence in the incidence medium
    polarization: str = "p"
    eps_in: complex = 1.0
    eps_out: complex = 1.0
    layers_before: Sequence = ()
    layers_after: Sequence = ()
    needs_eps: bool = True
    n_sub: int = 16

    def enhancement(self, eps: complex) -> float:
        """<|E|^2>_film / |E_inc|^2 for the current film permittivity."""
        layers = list(self.layers_before) + [(eps, self.thickness)] + list(self.layers_after)
        return layer_average_field_sq(
            layers,
            index=len(self.layers_before),
            wavelength=self.pulse.wavelength,
            angle_deg=self.angle_deg,
            polarization=self.polarization,
            eps_in=self.eps_in,
            eps_out=self.eps_out,
            n_sub=self.n_sub,
        )

    def power(self, t: float, eps: complex) -> float:
        n_in = math.sqrt(abs(complex(self.eps_in).real))
        e_inc_sq = c.intensity_to_field(self.peak_intensity_w_cm2, n_in) ** 2 * float(self.pulse.envelope(t))
        return 0.5 * self.pulse.omega * c.EPS0 * max(eps.imag, 0.0) * self.enhancement(eps) * e_inc_sq

    def local_intensity(self, t: float, eps: complex) -> float:
        """|E_loc|^2 / (2 Z0) averaged over the film, in W/cm^2."""
        n_in = math.sqrt(abs(complex(self.eps_in).real))
        return float(self.enhancement(eps) * self.peak_intensity_w_cm2 * self.pulse.envelope(t) / n_in)

    def absorptance(self, eps: complex) -> float:
        """Fraction of the incident power absorbed in the film (per unit area, at the given angle)."""
        n_in = math.sqrt(abs(complex(self.eps_in).real))
        cos_t = math.cos(math.radians(self.angle_deg))
        # absorbed per area = (omega/2) eps0 eps'' <|E|^2> d ; incident flux along z = n |E|^2 cos / (2 Z0)
        k0 = self.pulse.omega / c.C_LIGHT
        return float(k0 * max(eps.imag, 0.0) * self.enhancement(eps) * self.thickness / (n_in * cos_t))
