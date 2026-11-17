"""Transfer-matrix method for planar multilayers (s and p polarization, complex permittivities).

Convention: fields ~ exp(i(k_x x + k_z z - omega t)); light is incident from z < 0 in a medium of
permittivity ``eps_in`` at angle ``angle_deg``; layers are listed in order of propagation; the exit
medium has ``eps_out``.  Im(k_z) >= 0 selects decaying waves.

Only what WP1 needs is implemented: reflection/transmission coefficients and the average of |E|^2 over one
layer (normalised to the incident |E|^2), which sets the absorbed power density in the adiabatic
self-consistent scheme of ref. [15].  WP3 replaces this module by a time-domain solver.
"""
from __future__ import annotations

import cmath
import math
from typing import Sequence

import numpy as np

from . import constants as c


def _kz(eps: complex, k0: float, kx: complex) -> complex:
    kz = cmath.sqrt(eps * k0 * k0 - kx * kx)
    if kz.imag < 0 or (kz.imag == 0 and kz.real < 0):
        kz = -kz
    return kz


def _solve(layers, wavelength, angle_deg, polarization, eps_in, eps_out):
    """Return (r, t, amplitudes) where amplitudes[j] = (A_j, B_j) are forward/backward tangential-field
    amplitudes at the *left* interface of medium j (j = 0 is the incidence medium, last is the exit medium).

    The tangential field is E_y for s polarization and H_y for p polarization.
    """
    k0 = 2.0 * math.pi / wavelength
    eps_list = [complex(eps_in)] + [complex(e) for e, _ in layers] + [complex(eps_out)]
    d_list = [0.0] + [float(d) for _, d in layers] + [0.0]
    kx = cmath.sqrt(eps_list[0]) * k0 * math.sin(math.radians(angle_deg))
    kz = [_kz(e, k0, kx) for e in eps_list]
    # admittance-like factor: continuity of (F, q dF/dz-ish):  s: q = kz ;  p: q = kz/eps
    q = [kz_j if polarization == "s" else kz_j / e for kz_j, e in zip(kz, eps_list)]
    n = len(eps_list)
    # propagate from the exit medium backwards: in the exit medium only a forward wave (amplitude 1)
    amps = [None] * n
    amps[-1] = (1.0 + 0j, 0.0 + 0j)
    for j in range(n - 2, -1, -1):
        a_next, b_next = amps[j + 1]
        # field and derivative-like quantity at the interface between j and j+1 (left side of j+1)
        f = a_next + b_next
        g = q[j + 1] * (a_next - b_next)
        # amplitudes in medium j at that interface (its right side)
        a_r = 0.5 * (f + g / q[j])
        b_r = 0.5 * (f - g / q[j])
        # move to the left interface of medium j
        ph = cmath.exp(1j * kz[j] * d_list[j])
        amps[j] = (a_r / ph, b_r * ph)
    a0, b0 = amps[0]
    r = b0 / a0
    t = 1.0 / a0
    amps = [(a / a0, b / a0) for a, b in amps]
    return r, t, amps, kz, kx, eps_list, d_list


def reflection_transmission(layers: Sequence, wavelength: float, angle_deg: float = 0.0, polarization: str = "p",
                            eps_in: complex = 1.0, eps_out: complex = 1.0):
    """Power reflectance R and transmittance T (T is meaningful for a lossless exit medium)."""
    r, t, amps, kz, kx, eps_list, _ = _solve(layers, wavelength, angle_deg, polarization, eps_in, eps_out)
    R = abs(r) ** 2
    if polarization == "s":
        T = abs(t) ** 2 * (kz[-1].real / kz[0].real)
    else:
        T = abs(t) ** 2 * ((kz[-1] / eps_list[-1]).real / (kz[0] / eps_list[0]).real)
    return float(R), float(T), r, t


def layer_average_field_sq(layers: Sequence, index: int, wavelength: float, angle_deg: float = 0.0,
                           polarization: str = "p", eps_in: complex = 1.0, eps_out: complex = 1.0,
                           n_sub: int = 16) -> float:
    """<|E|^2> averaged over layer ``index`` (0-based within ``layers``), in units of the incident |E|^2."""
    r, t, amps, kz, kx, eps_list, d_list = _solve(layers, wavelength, angle_deg, polarization, eps_in, eps_out)
    j = index + 1
    a, b = amps[j]
    d = d_list[j]
    z = (np.arange(n_sub) + 0.5) * d / n_sub
    fwd = a * np.exp(1j * kz[j] * z)
    bwd = b * np.exp(-1j * kz[j] * z)
    if polarization == "s":
        # tangential field is E_y itself; incident amplitude 1
        e2 = np.abs(fwd + bwd) ** 2
    else:
        # tangential field is H_y;  E_x = (kz/(omega eps0 eps)) (fwd - bwd),  E_z = -(kx/(omega eps0 eps)) (fwd + bwd)
        # incident |E|^2 = |H|^2 |k_0n|^2 / (omega eps0 eps_in)^2  ->  normalise by k_in^2/eps_in^2
        k0 = 2.0 * math.pi / wavelength
        eps_j = eps_list[j]
        ex = kz[j] / eps_j * (fwd - bwd)
        ez = -kx / eps_j * (fwd + bwd)
        norm = abs(cmath.sqrt(eps_list[0]) * k0 / eps_list[0]) ** 2
        e2 = (np.abs(ex) ** 2 + np.abs(ez) ** 2) / norm
    return float(np.mean(e2))


def absorptance_of_layer(layers: Sequence, index: int, wavelength: float, angle_deg: float = 0.0,
                         polarization: str = "p", eps_in: complex = 1.0, eps_out: complex = 1.0,
                         n_sub: int = 64) -> float:
    """Fraction of incident power absorbed in layer ``index``: k0 eps'' <|E|^2> d / (n_in cos theta)."""
    eps_j, d = complex(layers[index][0]), float(layers[index][1])
    k0 = 2.0 * math.pi / wavelength
    n_in = cmath.sqrt(complex(eps_in)).real
    avg = layer_average_field_sq(layers, index, wavelength, angle_deg, polarization, eps_in, eps_out, n_sub)
    return float(k0 * eps_j.imag * avg * d / (n_in * math.cos(math.radians(angle_deg))))
