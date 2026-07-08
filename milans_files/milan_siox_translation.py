"""
python translation of milan's matlab files:
  - SiOxOCs.m
  - SiOxmultilayers.m

this is intentionally close to the matlab, including its sign conventions,
matlab 1-based indexing choices, and matrix ordering. it is not a cleaned-up
or corrected thin-film implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.io import loadmat


@dataclass
class MilanResult:
    wavenumber_cm: np.ndarray
    absorbance_ra: np.ndarray
    r_stack: np.ndarray
    r_bare_gold: np.ndarray
    n_siox: np.ndarray
    n_gold: np.ndarray
    experimental_wavenumber_cm: Optional[np.ndarray] = None
    experimental_shifted: Optional[np.ndarray] = None


def _load_first_nonmeta(path: str | Path, preferred: str | None = None) -> np.ndarray:
    """load a matlab .mat array by preferred name, else by sole non-metadata variable."""
    mat = loadmat(path)
    keys = [k for k in mat.keys() if not k.startswith("__")]
    if preferred is not None and preferred in mat:
        return mat[preferred]
    if len(keys) == 1:
        return mat[keys[0]]
    raise KeyError(f"could not choose variable in {path}; found {keys}")


def milan_gold_index(wavenumber_cm: np.ndarray) -> np.ndarray:
    """
    milan's empirical gold model from SiOxmultilayers.m.

    note: K below is the spectral wavenumber in cm^-1.
    kappa_au is the extinction coefficient of gold, not the wavevector.
    """
    K = np.asarray(wavenumber_cm, dtype=np.complex128)
    kappa_au = 0.95 * 14000.0 / (K ** 0.75) - 10.0
    n_au = 20.0 * kappa_au / (24.0 + 0.1 * K)
    return n_au + 1j * kappa_au


def interface_coefficients_milan(
    Ni: complex,
    Nf: complex,
    ui: complex,
    uf: complex,
    polar: int,
) -> tuple[complex, complex]:
    """
    fresnel-like coefficients exactly as in milan's matlab.

    Ni, Nf are complex refractive indices of incident/final media.
    ui = sqrt(Ni^2 - n_air^2 sin(theta0)^2) = Ni*cos(theta_i)
    uf = sqrt(Nf^2 - n_air^2 sin(theta0)^2) = Nf*cos(theta_f)

    polar=1 -> p polarization, polar=0 -> s polarization.
    """
    Ni2 = Ni ** 2
    Nf2 = Nf ** 2
    if polar == 1:  # p polarization
        r = -(ui * Nf2 - uf * Ni2) / (ui * Nf2 + uf * Ni2)
        t = 2.0 * Ni * Nf * uf / (Nf2 * ui + Ni2 * uf)
    elif polar == 0:  # s polarization
        r = (ui - uf) / (ui + uf)
        t = 2.0 * uf / (ui + uf)
    else:
        raise ValueError("polar must be 1 for p or 0 for s")
    return r, t

#"/mnt/data",
def simulate_siox_multilayers(
    data_dir: str | Path = ".",
    x: float = 1.15,
    td_nm_like: float = 10.0,
    no_column_matlab: int = 2,
    n_air: float = 1.0,
    d_unit: float = 1e-9,
    L: int = 1,
    angle_deg: float = 80.0,
    polar: int = 1,
) -> MilanResult:
    """
    close translation of SiOxmultilayers.m.

    important matlab-to-python notes:
    - no_column_matlab=2 means column index 1 in python.
    - q(1660, No) in matlab is q[1659, no_column_matlab-1] in python.
    - the code propagates through the SiOx layer inside the loop and then
      again immediately before the gold interface. with L=1 this produces
      two propagation factors through SiOx, i.e. effective thickness 2*td.
    - td = d_unit * td_nm_like, exactly as in matlab. no unit correction is
      made here.
    """
    data_dir = Path(data_dir)
    Z = _load_first_nonmeta(data_dir / "SiOxRI.mat", "SiOxRI")
    # the uploaded file contains variable SN, while the matlab script uses exspa
    Qexp = _load_first_nonmeta(data_dir / "exspa.mat", "exspa")
    Coeff = _load_first_nonmeta(data_dir / "coeff.mat", "coeff")

    n12 = n_air ** 2
    td = d_unit * td_nm_like
    angle = np.deg2rad(angle_deg)
    s = np.sin(angle)
    s2 = s * s

    # matlab: for p=1:N-1, K=Z(p+1,1). skip header row.
    K_values = np.real(Z[1:, 0]).astype(float)
    n_points = len(K_values)

    # empirical gold n+i*kappa model.
    N_AU = milan_gold_index(K_values)

    # matlab experimental pre-processing. this is computed but not used later.
    col = no_column_matlab - 1
    exp_wavenumber = None
    exp_shifted = None
    if Qexp is not None and Qexp.ndim == 2 and Qexp.shape[1] > col and Qexp.shape[0] >= 1660:
        temp = Qexp[1659, col]  # matlab Q(1660, No)
        exp_wavenumber = Qexp[:, 0].astype(float)
        exp_shifted = Qexp[:, col].astype(float) - float(temp)

    A = np.empty(n_points, dtype=float)
    r_stack = np.empty(n_points, dtype=np.complex128)
    r_bare_gold = np.empty(n_points, dtype=np.complex128)
    n_siox = np.empty(n_points, dtype=np.complex128)

    for p, K in enumerate(K_values):
        # interpolate/extrapolate SiOx optical constant at composition x.
        # matlab: Pol=[Coeff(p,1) Coeff(p,2) Coeff(p,3)]; Nf=polyval(Pol,x)
        Pol = Coeff[p, :3]
        Nf = np.polyval(Pol, x)
        n_siox[p] = Nf

        # initial interface: air -> SiOx
        Ni = complex(n_air)
        Ni2 = Ni ** 2
        Nf2 = Nf ** 2
        ui = np.sqrt(Ni2 - n12 * s2 + 0j)
        uf = np.sqrt(Nf2 - n12 * s2 + 0j)
        r, t = interface_coefficients_milan(Ni, Nf, ui, uf, polar)
        ti = 1.0 / t
        M = np.array([[ti, -r * ti], [-r * ti, ti]], dtype=np.complex128)

        # matlab loop: for ss=1:L. for L=1 and same Nf on both sides, this is
        # basically propagation through a SiOx slab.
        for _ in range(L):
            Ni = Nf
            Ni2 = Nf2
            Nf = np.polyval(Pol, x)
            Nf2 = Nf ** 2
            ui = uf
            uf = np.sqrt(Nf2 - n12 * s2 + 0j)
            D = np.exp(2.0 * np.pi * 1j * K * td * ui)
            Di = 1.0 / D
            r, t = interface_coefficients_milan(Ni, Nf, ui, uf, polar)
            ti = 1.0 / t
            Q = np.array([[D * ti, -r * ti * Di], [-r * ti * D, Di * ti]], dtype=np.complex128)
            M = Q @ M

        # final block: another SiOx propagation factor, then interface SiOx -> Au.
        Ni = Nf
        Ni2 = Nf2
        Nf = N_AU[p]
        Nf2 = Nf ** 2
        ui = uf
        uf = np.sqrt(Nf2 - n12 * s2 + 0j)
        D = np.exp(2.0 * np.pi * 1j * K * td * ui)
        Di = 1.0 / D
        r, t = interface_coefficients_milan(Ni, Nf, ui, uf, polar)
        ti = 1.0 / t
        Qmat = np.array([[D * ti, -r * ti * Di], [-r * ti * D, Di * ti]], dtype=np.complex128)
        M = Qmat @ M

        # bare air -> gold reference reflection amplitude r14.
        Ni = complex(n_air)
        Ni2 = n12
        ui = np.sqrt(n12 - n12 * s2 + 0j)
        # Nf, Nf2, uf are still gold from above.
        if polar == 1:
            r14 = -(ui * Nf2 - uf * Ni2) / (ui * Nf2 + uf * Ni2)
        else:
            r14 = (ui - uf) / (ui + uf)

        # matlab: A(p)=-2*log10(abs(M(2,1)/(M(2,2)*r14)))
        r_multilayer = M[1, 0] / M[1, 1]
        A[p] = -2.0 * np.log10(abs(r_multilayer / r14))
        r_stack[p] = r_multilayer
        r_bare_gold[p] = r14

    return MilanResult(
        wavenumber_cm=K_values,
        absorbance_ra=A,
        r_stack=r_stack,
        r_bare_gold=r_bare_gold,
        n_siox=n_siox,
        n_gold=N_AU,
        experimental_wavenumber_cm=exp_wavenumber,
        experimental_shifted=exp_shifted,
    )

#"/mnt/data",
def siox_ocs_curves(
    data_dir: str | Path = '.',
    x_values: np.ndarray | None = None,
) -> tuple[np.ndarray, dict[float, np.ndarray]]:
    """
    python version of the commented plotting section in SiOxOCs.m.

    returns wavenumber and real(polyval(coeff, x)) for each x value.
    """
    data_dir = Path(data_dir)
    Z = _load_first_nonmeta(data_dir / "SiOxRI.mat", "SiOxRI")
    Coeff = _load_first_nonmeta(data_dir / "coeff.mat", "coeff")
    K_values = np.real(Z[1:, 0]).astype(float)
    if x_values is None:
        # matlab commented loop: pp=1:11; x=1+(pp-1)*0.1
        x_values = 1.0 + 0.1 * np.arange(11)
    curves: dict[float, np.ndarray] = {}
    for x in x_values:
        Nf = np.array([np.polyval(Coeff[p, :3], x) for p in range(len(K_values))])
        curves[float(x)] = np.real(Nf)
    return K_values, curves


if __name__ == "__main__":
    result = simulate_siox_multilayers()
    print("computed", len(result.wavenumber_cm), "points")
    print("wavenumber range:", result.wavenumber_cm[0], "to", result.wavenumber_cm[-1], "cm^-1")
    print("absorbance_ra min/max:", np.nanmin(result.absorbance_ra), np.nanmax(result.absorbance_ra))
    print("first five rows: wavenumber_cm, absorbance_ra")
    for K, A in zip(result.wavenumber_cm[:5], result.absorbance_ra[:5]):
        print(f"{K:10.3f}  {A: .9e}")
