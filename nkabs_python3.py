"""modernized implementation of the rocha-pilling nkabs transmission inversion.

this module preserves the published fixed-point algorithm while making the
following deliberate changes:

1. python 3 and in-memory numpy/pandas calculations.
2. film thickness is supplied directly in cm.
3. wavenumber is supplied in cm^-1.
4. the ambient medium is explicitly air/vacuum, n_ambient = 1.
5. the visible film index is used only as the kramers-kronig anchor.
6. descending ftir wavenumber data are sorted internally and restored on output.
7. the maclaurin kramers-kronig grid spacing uses n_points - 1 intervals.
8. all constants use numpy values rather than rounded literals.
9. a maximum iteration count and numerical checks prevent silent infinite loops.

this remains a normal-incidence transmission model for one isotropic film on a
transparent, nonabsorbing substrate. it is not yet the gold-reflection model.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any
import warnings

import numpy as np
import pandas as pd

try:
    from numba import njit
except ImportError:  # numba is optional; the numpy fallback remains valid
    njit = None


# explicit mathematical constants used by the rocha-pilling equations
LN_10 = np.log(10.0)
TWO_PI = 2.0 * np.pi
FOUR_PI = 4.0 * np.pi
N_AMBIENT = 1.0  # air/vacuum at normal incidence


def load_absorbance_file(
    path: str | Path = "example1-spectrum.txt",
    *,
    comment: str = "#",
) -> pd.DataFrame:
    """load a two-column whitespace-delimited absorbance spectrum.

    parameters
    ----------
    path
        path relative to the current working directory unless an absolute path
        is supplied. column 1 must be wavenumber in cm^-1 and column 2 must be
        absorbance.
    comment
        lines beginning with this character are ignored.

    returns
    -------
    pandas.dataframe
        columns are ``wavenumber_cm1`` and ``absorbance``.
    """

    path = Path(path)
    data = pd.read_csv(
        path,
        sep=r"\s+",
        comment=comment,
        header=None,
        usecols=[0, 1],
        names=["wavenumber_cm1", "absorbance"],
    )

    if data.empty:
        raise ValueError(f"no spectral data were found in {path!s}")

    return data


def _validate_and_sort_spectrum(
    wavenumber_cm1: np.ndarray,
    absorbance: np.ndarray,
    *,
    spacing_rtol: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """validate the spectrum and sort it into ascending wavenumber order.

    the maclaurin principal-value quadrature used by nkabs requires a uniformly
    spaced grid. ftir exports commonly run from high to low wavenumber, which is
    fine: the function sorts internally and later restores the user's order.
    """

    nu_input = np.asarray(wavenumber_cm1, dtype=float)
    absorbance_input = np.asarray(absorbance, dtype=float)

    if nu_input.ndim != 1 or absorbance_input.ndim != 1:
        raise ValueError("wavenumber and absorbance must be one-dimensional arrays")
    if nu_input.size != absorbance_input.size:
        raise ValueError("wavenumber and absorbance must contain the same number of points")
    if nu_input.size < 4:
        raise ValueError("at least four spectral points are required")
    if not np.all(np.isfinite(nu_input)) or not np.all(np.isfinite(absorbance_input)):
        raise ValueError("wavenumber and absorbance must contain only finite values")
    if np.any(nu_input <= 0.0):
        raise ValueError("all wavenumbers must be positive and expressed in cm^-1")

    # stable sorting handles either ordinary ascending data or the usual ftir
    # convention in which the largest wavenumber appears first.
    sort_index = np.argsort(nu_input, kind="mergesort")
    nu = nu_input[sort_index]
    absorbance_sorted = absorbance_input[sort_index]

    spacing = np.diff(nu)
    if np.any(spacing <= 0.0):
        raise ValueError("duplicate wavenumbers are not allowed")

    # equation (12) of rocha and pilling assumes a constant interval h.
    h = float(np.mean(spacing))
    if not np.allclose(spacing, h, rtol=spacing_rtol, atol=0.0):
        max_fractional_deviation = float(np.max(np.abs(spacing - h)) / h)
        raise ValueError(
            "the maclaurin kramers-kronig method requires uniformly spaced "
            f"wavenumbers; maximum fractional spacing deviation is "
            f"{max_fractional_deviation:.3e}. interpolate to a uniform grid first."
        )

    return nu, absorbance_sorted, sort_index, h



def _kk_maclaurin_python(
    nu: np.ndarray,
    k: np.ndarray,
    n_anchor: float,
    h: float,
) -> np.ndarray:
    """pure-python/numpy fallback for the alternating-point kk sum."""

    n = np.empty_like(k)
    indices = np.arange(nu.size)

    for i, nu_i in enumerate(nu):
        use = (indices % 2) != (i % 2)
        nu_j = nu[use]
        k_j = k[use]
        integrand = 0.5 * k_j * (
            1.0 / (nu_j - nu_i) + 1.0 / (nu_j + nu_i)
        )
        n[i] = n_anchor + (2.0 / np.pi) * (2.0 * h) * np.sum(integrand)

    return n


if njit is not None:

    @njit(cache=True)
    def _kk_maclaurin_numba(
        nu: np.ndarray,
        k: np.ndarray,
        n_anchor: float,
        h: float,
    ) -> np.ndarray:
        """compiled version of the O(N^2) maclaurin principal-value sum."""

        n_points = nu.size
        n = np.empty_like(k)
        prefactor = (2.0 / np.pi) * (2.0 * h)

        for i in range(n_points):
            total = 0.0
            # use source indices of parity opposite to the target index
            first_j = 1 if i % 2 == 0 else 0
            for j in range(first_j, n_points, 2):
                total += 0.5 * k[j] * (
                    1.0 / (nu[j] - nu[i]) + 1.0 / (nu[j] + nu[i])
                )
            n[i] = n_anchor + prefactor * total

        return n


def kramers_kronig_maclaurin(
    wavenumber_cm1: np.ndarray,
    k: np.ndarray,
    n_anchor: float,
    *,
    grid_spacing_cm1: float | None = None,
) -> np.ndarray:
    r"""calculate n from k using the alternating-point maclaurin quadrature.

    this evaluates the finite-window relation used in the nkabs paper,

    .. math::

       n(\nu_i) = n_\mathrm{anchor}
       + \frac{2}{\pi}(2h)\sum_j
       \frac{1}{2}\left[
       \frac{k_j}{\nu_j-\nu_i}+\frac{k_j}{\nu_j+\nu_i}\right],

    where only indices j of parity opposite to i are included. this avoids the
    singular point j = i without evaluating it directly.

    notes
    -----
    this is a truncated kramers-kronig integral over the measured interval, not
    a full zero-to-infinity transform. edge errors and missing-band effects are
    therefore still present.
    """

    nu = np.asarray(wavenumber_cm1, dtype=float)
    k = np.asarray(k, dtype=float)

    if nu.ndim != 1 or k.ndim != 1 or nu.size != k.size:
        raise ValueError("wavenumber and k must be one-dimensional arrays of equal length")

    if grid_spacing_cm1 is None:
        spacing = np.diff(nu)
        h = float(np.mean(spacing))
    else:
        h = float(grid_spacing_cm1)

    # the maclaurin sum is O(N^2). use numba automatically when available;
    # otherwise fall back to the transparent numpy implementation above.
    if njit is not None:
        return _kk_maclaurin_numba(nu, k, float(n_anchor), h)

    return _kk_maclaurin_python(nu, k, float(n_anchor), h)


def transmission_model(
    wavenumber_cm1: np.ndarray,
    thickness_cm: float,
    n_film: np.ndarray,
    k_film: np.ndarray,
    n_substrate: float,
) -> tuple[np.ndarray, np.ndarray]:
    """calculate normalized theoretical transmittance and fresnel correction.

    geometry at normal incidence:

        air/vacuum (material 0) | film (material 1) | substrate (material 2)

    the returned transmittance is normalized to the bare air-substrate interface,
    matching the t01*t12/t02 structure used by rocha and pilling.
    """

    nu = np.asarray(wavenumber_cm1, dtype=float)
    m_film = np.asarray(n_film, dtype=float) + 1j * np.asarray(k_film, dtype=float)
    m0 = complex(N_AMBIENT, 0.0)
    m2 = complex(float(n_substrate), 0.0)

    # normal-incidence fresnel amplitude coefficients
    t01 = 2.0 * m0 / (m0 + m_film)
    t12 = 2.0 * m_film / (m_film + m2)
    t02 = 2.0 * m0 / (m0 + m2)
    r01 = (m0 - m_film) / (m0 + m_film)
    r12 = (m_film - m2) / (m_film + m2)

    # one-pass complex phase delta = 2*pi*nu*d*m_film.
    # exp(2j*delta) is the round-trip phase/attenuation factor.
    delta = TWO_PI * nu * thickness_cm * m_film
    round_trip = np.exp(2j * delta)

    fresnel_amplitude = (t01 * t12 / t02) / (1.0 + r01 * r12 * round_trip)
    fresnel_correction = np.abs(fresnel_amplitude)

    # |exp(i*delta)|^2 = exp(-4*pi*nu*d*k), the single-pass intensity loss.
    lambert_factor = np.exp(-FOUR_PI * nu * thickness_cm * k_film)
    transmittance = lambert_factor * fresnel_correction**2

    return transmittance, fresnel_correction


def nkabs_transmission(
    wavenumber_cm1: np.ndarray | pd.Series,
    absorbance: np.ndarray | pd.Series,
    *,
    thickness_cm: float,
    n_visible_film: float,
    n_substrate: float,
    mape_tolerance_percent: float = 1.0e-3,
    max_iterations: int = 100,
    absorbance_scale: str = "decadic",
    spacing_rtol: float = 5.0e-3,
    output_csv: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """recover film n and k from a normal-incidence transmission absorbance.

    parameters
    ----------
    wavenumber_cm1, absorbance
        one-dimensional spectral arrays. wavenumber must be in cm^-1. either
        ascending or descending order is accepted.
    thickness_cm
        physical film thickness in cm. no unit conversion is performed.
    n_visible_film
        real visible/high-frequency refractive index used as the additive
        kramers-kronig anchor for the film. it is not used as the ambient index.
    n_substrate
        one real, wavelength-independent substrate refractive index. this is the
        next assumption to relax when a dispersive substrate is introduced.
    mape_tolerance_percent
        convergence threshold in percent, using the mape definition from the
        nkabs paper with theoretical transmittance in the denominator.
    max_iterations
        hard stop preventing an infinite fixed-point loop.
    absorbance_scale
        ``"decadic"`` for A = log10(I0/I), matching the original python code,
        or ``"natural"`` for optical depth tau = ln(I0/I), matching the equation
        as printed in the paper.
    spacing_rtol
        relative tolerance used to decide whether the spectral grid is uniform.
    output_csv
        optional final csv path. relative paths use the notebook's current
        working directory. no intermediate files are written.
    verbose
        print one concise line per iteration.

    returns
    -------
    dict
        ``data`` is the final dataframe in the same wavenumber order supplied by
        the user; ``history`` records convergence by iteration; ``summary`` gives
        parameters and final diagnostics.
    """

    thickness_cm = float(thickness_cm)
    n_visible_film = float(n_visible_film)
    n_substrate = float(n_substrate)
    tolerance_fraction = float(mape_tolerance_percent) / 100.0

    if thickness_cm <= 0.0:
        raise ValueError("thickness_cm must be positive")
    if n_visible_film <= 0.0 or n_substrate <= 0.0:
        raise ValueError("refractive indices must be positive")
    if tolerance_fraction <= 0.0:
        raise ValueError("mape_tolerance_percent must be positive")
    if max_iterations < 1:
        raise ValueError("max_iterations must be at least one")

    input_wavenumber = np.asarray(wavenumber_cm1, dtype=float)
    input_absorbance = np.asarray(absorbance, dtype=float)

    nu, absorbance_sorted, sort_index, h = _validate_and_sort_spectrum(
        input_wavenumber,
        input_absorbance,
        spacing_rtol=spacing_rtol,
    )

    scale = absorbance_scale.strip().lower()
    if scale == "decadic":
        # base-10 ftir absorbance A = log10(I0/I)
        optical_depth = LN_10 * absorbance_sorted
    elif scale == "natural":
        # natural-log absorbance, often called optical depth tau = ln(I0/I)
        optical_depth = absorbance_sorted.copy()
    else:
        raise ValueError("absorbance_scale must be 'decadic' or 'natural'")

    experimental_transmittance = np.exp(-optical_depth)

    # first iteration uses a unit fresnel/interference correction, as in nkabs.
    correction = np.ones_like(nu)
    history_rows: list[dict[str, float | int]] = []
    converged = False

    for iteration in range(1, max_iterations + 1):
        if np.any(correction <= 0.0) or not np.all(np.isfinite(correction)):
            raise FloatingPointError("invalid fresnel correction encountered")

        # rearranged lambert-fresnel expression, equation (8) in the paper.
        alpha_cm1 = (optical_depth + 2.0 * np.log(correction)) / thickness_cm
        k_film = alpha_cm1 / (FOUR_PI * nu)

        # no clipping is imposed: negative k values are retained so that baseline
        # or model problems remain visible rather than being silently hidden.
        n_film = kramers_kronig_maclaurin(
            nu,
            k_film,
            n_visible_film,
            grid_spacing_cm1=h,
        )

        theoretical_transmittance, new_correction = transmission_model(
            nu,
            thickness_cm,
            n_film,
            k_film,
            n_substrate,
        )

        if not np.all(np.isfinite(theoretical_transmittance)):
            raise FloatingPointError("nonfinite theoretical transmittance encountered")

        residual = experimental_transmittance - theoretical_transmittance
        safe_theoretical = np.maximum(
            theoretical_transmittance,
            np.finfo(float).tiny,
        )

        # preserve the paper's definitions, while labeling chi-square honestly.
        mape_fraction = float(np.mean(np.abs(residual / safe_theoretical)))
        rocha_chi_square = float(np.sum(residual**2 / safe_theoretical))
        rmse = float(np.sqrt(np.mean(residual**2)))

        history_rows.append(
            {
                "iteration": iteration,
                "mape_percent": 100.0 * mape_fraction,
                "rocha_chi_square": rocha_chi_square,
                "rmse_transmittance": rmse,
                "minimum_k": float(np.min(k_film)),
                "maximum_k": float(np.max(k_film)),
            }
        )

        if verbose:
            print(
                f"iteration {iteration:3d}: "
                f"mape = {100.0 * mape_fraction:.6g}%  "
                f"rmse = {rmse:.6g}"
            )

        correction = new_correction

        if mape_fraction <= tolerance_fraction:
            converged = True
            break

    if not converged:
        warnings.warn(
            "nkabs did not reach the requested mape tolerance before "
            f"max_iterations={max_iterations}",
            RuntimeWarning,
            stacklevel=2,
        )

    # all arrays above are in ascending order. restore the user's original order.
    inverse_sort = np.empty_like(sort_index)
    inverse_sort[sort_index] = np.arange(sort_index.size)

    final_residual = experimental_transmittance - theoretical_transmittance
    data_sorted = pd.DataFrame(
        {
            "wavenumber_cm1": nu,
            "absorbance_input": absorbance_sorted,
            "experimental_transmittance": experimental_transmittance,
            "n_film": n_film,
            "k_film": k_film,
            "theoretical_transmittance": theoretical_transmittance,
            "fresnel_correction_amplitude": correction,
            "transmittance_residual": final_residual,
        }
    )
    data = data_sorted.iloc[inverse_sort].reset_index(drop=True)

    history = pd.DataFrame(history_rows)
    final_metrics = history.iloc[-1]
    summary = {
        "converged": converged,
        "iterations": int(final_metrics["iteration"]),
        "mape_percent": float(final_metrics["mape_percent"]),
        "rocha_chi_square": float(final_metrics["rocha_chi_square"]),
        "rmse_transmittance": float(final_metrics["rmse_transmittance"]),
        "thickness_cm": thickness_cm,
        "n_visible_film": n_visible_film,
        "n_ambient": N_AMBIENT,
        "n_substrate": n_substrate,
        "absorbance_scale": scale,
        "grid_spacing_cm1": h,
        "number_of_points": int(nu.size),
        "input_wavenumber_descending": bool(input_wavenumber[0] > input_wavenumber[-1]),
    }

    if output_csv is not None:
        output_path = Path(output_csv)
        data.to_csv(output_path, index=False)
        summary["output_csv"] = str(output_path.resolve())

    return {"data": data, "history": history, "summary": summary}
