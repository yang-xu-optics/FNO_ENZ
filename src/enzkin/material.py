"""Material parameter sets.

The default :func:`ito_sarkar2023` reproduces the indium tin oxide parameters of Sarkar, Un & Sivan
(ref. [14]) and Un, Sarkar & Sivan (ref. [15]) so that results compare directly with those papers.
Every field is a plain float so parameter sets can be sampled (WP1 dataset) and inferred (WP4).
"""
from __future__ import annotations

from dataclasses import dataclass, replace

from . import constants as c
from .band import KaneBand


@dataclass(frozen=True)
class Material:
    # --- band and carriers -------------------------------------------------------------
    n_e: float  # conduction-electron density, 1/nm^3   (1 nm^-3 = 1e27 m^-3)
    m0_rel: float  # band-edge effective mass / m_e
    C: float  # Kane non-parabolicity, 1/eV
    E_max: float  # top of the modelled conduction band, eV
    # --- dielectric -----------------------------------------------------------------------
    eps_b: float  # background permittivity screening the e-e Coulomb interaction
    eps_inf: float  # high-frequency (interband) permittivity entering eps(omega)
    # --- electron-electron ---------------------------------------------------------------
    c_ee: float  # dimensionless prefactor of the screened-Coulomb collision integral (see ee.py)
    # --- electron-phonon (acoustic deformation potential, Debye spectrum) ------------------
    D_def: float  # deformation potential, eV
    rho_mass: float  # mass density, kg/m^3
    v_sound: float  # sound velocity, m/s
    E_debye: float  # Debye energy, eV
    C_ph: float  # lattice heat capacity, J/(m^3 K)
    # --- momentum relaxation not affecting f(E) --------------------------------------------
    gamma_imp: float  # electron-impurity collision rate entering the damping only, 1/fs
    name: str = "material"
    # --- optional energy-dependent impurity damping (Brooks-Herring, collisions/eimp.py) ---------
    n_imp: float = 0.0  # ionised-impurity density, 1/nm^3 (0 switches the term off)
    Z_imp: float = 1.0  # impurity charge

    # ---- derived -------------------------------------------------------------------------
    @property
    def band(self) -> KaneBand:
        return KaneBand(self.m0_rel, self.C)

    @property
    def hbar_v_sound(self) -> float:
        """hbar * v_sound in eV nm."""
        return c.HBAR * self.v_sound * 1e-6  # m/s -> nm/fs is 1e9/1e15

    @property
    def rho_mass_working(self) -> float:
        return self.rho_mass * c.KG_PER_M3_IN_WORKING

    @property
    def C_ph_working(self) -> float:
        """Lattice heat capacity in eV/(nm^3 K)."""
        return self.C_ph * c.J_PER_M3_IN_EV_PER_NM3

    @property
    def kappa_eph(self) -> float:
        """Electron-phonon pair-kernel constant  m0^2 D^2 / (4 pi^3 hbar^3 rho (hbar v_s)^4),  1/(eV^4 fs nm^3).

        With it, the number-conserving form of [14, Eq. (9)] reads
        d(rho f)/dt |_E = kappa (1+2CE) int d eps eps^2 (1+2C(E +- eps)) [...]  (see collisions/eph.py),
        and the T -> 0 coupling is G = kappa (1+2C mu)^2 k_B (2 hbar v_s k_F)^4 / 4 = [14, Eq. (27)].
        """
        import math

        m0 = self.band.m0
        return (m0**2 * self.D_def**2) / (
            4.0 * math.pi**3 * c.HBAR**3 * self.rho_mass_working * self.hbar_v_sound**4
        )

    def with_(self, **kw) -> "Material":
        return replace(self, **kw)


def ito_sarkar2023(**overrides) -> Material:
    """ITO as in refs. [14, 15].

    Values and where they come from:

    ========== ================= ==================================================================
    n_e        1.5 nm^-3         [14] Sec. II  (1.5e27 m^-3, giving E_F = 0.88 eV)
    m0_rel     0.3964            [14] Eq. (1) text, from Liu et al. 2014
    C          0.4191 eV^-1      [14] Eq. (1) text
    E_max      4.0 eV            [14] Sec. II ("total size" of the conduction band).  A 4 eV band saturates
                                 above ~30 GW/cm^2 in the geometry of [15]; E_max = 6 eV (with N scaled to keep
                                 dE) reproduces that paper's peak temperatures to 5 %  (docs/WORKLOG.md 1.5)
    eps_b      4.0               [14] below Eq. (14)
    eps_inf    3.45              not printed in [14, 15]; set so that Re eps = 0 at 225 THz and 300 K,
                                 the ENZ point quoted in [15] Sec. III  (see docs/theory/wp1_kinetic_model.md)
    c_ee       3/16              makes the prefactor equal to [14] Eq. (15) as printed; c_ee = 1 is the
                                 spin-summed direct-term golden-rule value (see collisions/ee.py)
    D_def      17.2 eV           [14] below Eq. (5)
    rho_mass   7120 kg/m^3       [14] below Eq. (5)
    v_sound    6400 m/s          [14] Sec. II B
    E_debye    0.0603 eV         [14] Sec. II B (0.06 eV; 60.3 meV in the source comments; q_D = 14.3 nm^-1)
    C_ph       2.54e6 J/m^3/K    [14] Sec. III B, [15] App. B
    gamma_imp  0.0 fs^-1         [15] includes an e-impurity rate in the damping but prints no value;
                                 treated as a per-sample fit parameter (ellipsometry / Hall mobility)
    ========== ================= ==================================================================
    """
    base = dict(
        n_e=1.5,
        m0_rel=0.3964,
        C=0.4191,
        E_max=4.0,
        eps_b=4.0,
        eps_inf=3.45,
        c_ee=3.0 / 16.0,
        D_def=17.2,
        rho_mass=7120.0,
        v_sound=6400.0,
        E_debye=0.0603,
        C_ph=2.54e6,
        gamma_imp=0.0,
        name="ITO (Sarkar/Un/Sivan 2023)",
    )
    base.update(overrides)
    return Material(**base)


def ito_un2026(**overrides) -> Material:
    """ITO as in Table S1 of Un, Sarkar & Sivan, Nanophotonics 15, e70272 (2026) (ref. [62]).

    Same band, density and lattice as :func:`ito_sarkar2023`, but the authors' later damping budget:

    ========== ================= ==================================================================
    D_def      8.6 eV            Table S1, "obtained by fitting the ITO permittivity at the bulk limit and
                                 at room temperature" (Wang et al. 2019, plan ref. [7]); [14] used 17.2 eV
                                 with the same e-ph formula, so G is 4 times smaller here
    n_imp      0.05 nm^-3        Table S1, Brooks-Herring ionised-impurity scattering (collisions/eimp.py)
    Z_imp      1                 Table S1
    gamma_imp  0.0 fs^-1         replaced by the energy-dependent Brooks-Herring rate
    eps_b      4.0               not in Table S1; (S6), (S9), (S18) screen with eps_b, (S3)-(S4) with eps_inf
    eps_inf    3.45              not printed; calibrate per sample (enzkin.calibrate)
    ========== ================= ==================================================================

    With eps_inf = 3.45 this set gives eps(230 THz, 300 K) = 0.14 + 0.17i, against 0.172 + 0.399i plotted in
    [15]: it describes a lower-loss film than [15] (scripts/check_un2026.py).
    """
    base = dict(D_def=8.6, n_imp=0.05, Z_imp=1.0, gamma_imp=0.0, name="ITO (Un/Sarkar/Sivan 2026, Table S1)")
    base.update(overrides)
    return ito_sarkar2023(**base)
