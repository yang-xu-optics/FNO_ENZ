"""Calibration of the two linear-optics parameters that refs. [14, 15] do not print.

``eps_inf`` and ``gamma_imp`` are fixed by one complex number: the permittivity of the cold film at one
frequency.  The default target is the value plotted in [15] (Fig. 1 inset): eps = 0.172 + 0.399i at
230 THz and 300 K.  For a measured sample, pass its ellipsometry value instead.
"""
from __future__ import annotations

from scipy.optimize import fsolve

from . import constants as c
from .material import Material
from .model import KineticModel

EPS_UN2023 = 0.172 + 0.399j
F_UN2023_THZ = 230.0


def calibrate_linear_optics(mat: Material, N: int = 128, f_thz: float = F_UN2023_THZ,
                            eps_target: complex = EPS_UN2023, T: float = 300.0) -> Material:
    """Return ``mat`` with eps_inf and gamma_imp such that eps(f_thz, T) = eps_target on an N-point grid."""
    w = c.frequency_thz_to_omega(f_thz)

    def resid(x):
        m = mat.with_(eps_inf=float(x[0]), gamma_imp=float(x[1]))
        model = KineticModel(m, N=N)
        f0, _ = model.split(model.initial_state(T))
        e = complex(model.permittivity(f0, T, w))
        return [e.real - eps_target.real, e.imag - eps_target.imag]

    x = fsolve(resid, [mat.eps_inf, max(mat.gamma_imp, 0.02)])
    return mat.with_(eps_inf=float(x[0]), gamma_imp=float(x[1]))
