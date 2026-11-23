#!/usr/bin/env python3
# Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
# SPDX-License-Identifier: MIT
"""
hrtf2hrir.py
Convert HRTF complex spectra (BEM/MLFMM output) to windowed HRIRs.

Usage:
  python hrtf2hrir.py config_L.txt [options]

The right-ear result file is derived automatically from the left-ear result_file
in the config: replaces the last _L with _R (or _G with _D for legacy naming).
Use --result-R to override.

Spatial convention: x=front, z=up, trigonometric sense (az=atan2(y,x)).

Result file format (Fortran output):
  N_ext lines per frequency block, N_freq blocks.
  Columns: Re(H)  Im(H)  iter  time

Frequency grid:
  f[k] = k * fe/(2*N_freq),  k = 1..N_freq
  f_min = fe/(2*N_freq),  f_max = fe/2,  N_time = 2*N_freq
"""

import argparse
import csv
import numpy as np
from pathlib import Path
from scipy.io import wavfile
from scipy.signal.windows import tukey


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def parse_config(config_file):
    """
    Parse config.txt ('!' comment line followed by value line).
    Returns (params dict, base directory Path).
    """
    config_file = Path(config_file)
    lines = [l.strip() for l in config_file.read_text().splitlines() if l.strip()]

    values = []
    i = 0
    while i < len(lines):
        if lines[i].startswith('!'):
            if i + 1 < len(lines):
                values.append(lines[i + 1])
            i += 2
        else:
            i += 1

    # first value lines of the case file, the same for Ray and Diff (format of 2026-09;
    # the 2024 format had level_max / level_min after the frequency line)
    keys = ['N_freq', 'fe', 'f_line', 'mesh_file', 'pext_file', 'result_file']
    params = dict(zip(keys, values))
    params['N_freq'] = int(params['N_freq'])
    params['fe']     = float(params['fe'])
    params['f_min']  = float(params['f_line'].replace(',', ' ').split()[0])
    return params, config_file.parent


def derive_result_r(stem):
    """
    Derive right-ear result stem from left-ear stem.
    Tries _L -> _R first, then _G -> _D (legacy French naming).
    Falls back to appending _R if neither pattern is found.
    """
    for src, dst in [('_L', '_R'), ('_G', '_D')]:
        idx = stem.rfind(src)
        if idx >= 0:
            return stem[:idx] + dst + stem[idx + len(src):]
    return stem + '_R'


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_pext(pext_file):
    """Load external points -> (N_ext, 3) Cartesian. Format: N_ext on line 1, then x y z."""
    return np.loadtxt(pext_file, skiprows=1)


def load_hrtf(result_file, N_freq, N_ext):
    """
    Load result file -> complex (N_freq, N_ext).
    Structure: N_ext lines (points) per frequency, N_freq blocks.
    Columns: Re(H)  Im(H)  iter  time
    """
    data = np.loadtxt(result_file)
    expected = N_freq * N_ext
    if data.shape[0] != expected:
        raise ValueError(
            f"{result_file}: expected {expected} lines "
            f"(N_freq={N_freq} x N_ext={N_ext}), got {data.shape[0]}."
        )
    return (data[:, 0] + 1j * data[:, 1]).reshape(N_freq, N_ext)


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def cartesian_to_spherical(pext):
    """
    (N_ext, 3) Cartesian -> az (deg), el (deg), r (m).
    Convention: x=front, z=up, trigonometric sense.
      az = atan2(y, x)   0 = front, 90 = left
      el = asin(z / r)   0 = horizontal, +90 = zenith
    """
    x, y, z = pext[:, 0], pext[:, 1], pext[:, 2]
    r  = np.sqrt(x**2 + y**2 + z**2)
    az = np.degrees(np.arctan2(y, x))
    el = np.degrees(np.arcsin(np.clip(z / r, -1.0, 1.0)))
    return az, el, r


# ---------------------------------------------------------------------------
# Processing
# ---------------------------------------------------------------------------

def diffuse_field_response(H_L, H_R, weights=None):
    """
    DFR(f) = sqrt( sum_i w_i * (|H_L(f,i)|^2 + |H_R(f,i)|^2) / 2 )

    weights: (N_ext,) quadrature weights of the directions (Lebedev, Fibonacci...), sum=1.
             If None, uniform arithmetic mean.
    Returns shape (N_freq,).
    """
    N_ext = H_L.shape[1]
    if weights is None:
        w = np.ones(N_ext) / N_ext
    else:
        w = np.asarray(weights, dtype=float)
        if w.shape != (N_ext,):
            raise ValueError(f"weights: expected {N_ext} values, got {w.shape}.")
        w = w / w.sum()

    power = 0.5 * (np.abs(H_L)**2 + np.abs(H_R)**2)  # (N_freq, N_ext)
    return np.sqrt(power @ w)                           # (N_freq,)


def asymmetric_window(N, n_fade_in, n_fade_out):
    """
    Asymmetric window: cosine fade-in over n_fade_in samples,
    flat centre, cosine fade-out over n_fade_out samples.
    Suited to HRIRs whose peak is at the start of the array.
    """
    w = np.ones(N)
    if n_fade_in > 0:
        w[:n_fade_in] = 0.5 * (1 - np.cos(np.pi * np.arange(n_fade_in) / n_fade_in))
    if n_fade_out > 0:
        w[-n_fade_out:] = 0.5 * (1 + np.cos(np.pi * np.arange(n_fade_out) / n_fade_out))
    return w


def hrtf_to_hrir(H, fe, f_min, N_out, n_fade_in, n_fade_out, tukey_alpha=None, delay=32):
    """
    One-sided HRTF spectrum -> windowed HRIR via IFFT.

    H          : complex (N_freq, N_ext), frequencies f_min..fe/2
    fe         : sampling rate (Hz)
    f_min      : first computed frequency (Hz), expected = fe/(2*N_freq)
    N_out      : output HRIR length in samples
    n_fade_in  : fade-in length in samples (default 4)
    n_fade_out : fade-out length in samples (default 8)
    tukey_alpha: if given, replaces asymmetric window with symmetric Tukey
    delay      : pure delay added (samples). The HRTF is normalised by the
                 free-field pressure at the head centre: the ipsilateral ear
                 receives the sound BEFORE t = 0 (up to ~0.3 ms) and the
                 contralateral one after (up to ~0.7 ms). Without a delay the
                 ipsilateral onset wraps to the end of the IFFT buffer and is cut.

    Returns real HRIR (N_out, N_ext).
    """
    N_freq, N_ext = H.shape
    N_time = 2 * N_freq

    df_expected = fe / N_time
    if abs(f_min - df_expected) > 1.0:
        print(f"  Warning: f_min={f_min:.4f} Hz != expected df={df_expected:.4f} Hz")

    if N_out > N_time:
        raise ValueError(
            f"N_out={N_out} > N_time={N_time} (=2*N_freq). "
            f"Reduce --n-out or increase N_freq."
        )

    # The solver uses the e^{-i omega t} convention; numpy's irfft assumes
    # e^{+i omega t}: conjugate, otherwise the response is reversed in time
    # (2026-11-20: before this fix the contralateral ear was lost).
    H = np.conj(H)
    dc = H[0].real                 # DC: constant extrapolation from f_min (Im=0), before the delay
    f = f_min * np.arange(1, N_freq + 1)
    H = H * np.exp(-2j * np.pi * f * delay / fe)[:, np.newaxis]

    # Full one-sided spectrum: DC (0) + positive frequencies + Nyquist
    H_full = np.empty((N_freq + 1, N_ext), dtype=complex)
    H_full[0]        = dc
    H_full[1:N_freq] = H[:-1]      # positive frequencies excluding Nyquist
    H_full[N_freq]   = H[-1].real  # Nyquist: Im forced to 0

    # IFFT -> real signal, peak at start of array
    hrir_full = np.fft.irfft(H_full, n=N_time, axis=0)  # (N_time, N_ext)

    # Windowing
    if tukey_alpha is not None:
        win = tukey(N_out, alpha=tukey_alpha)
    else:
        win = asymmetric_window(N_out, n_fade_in, n_fade_out)

    return hrir_full[:N_out] * win[:, np.newaxis]  # (N_out, N_ext)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def save_npz(out_file, hrir_L, hrir_R, pext, az, el, r, fe, equalized):
    np.savez(out_file,
             hrir_L=hrir_L, hrir_R=hrir_R,
             pext=pext, az=az, el=el, r=r,
             fe=np.array(fe), equalized=np.array(equalized))
    print(f"NPZ: {out_file}  {hrir_L.shape}")


def save_wav_dir(wav_dir, hrir_L, hrir_R, pext, az, el, r, fe):
    """
    Save float32 WAV files + index CSV.
      wav_dir/index.csv    - index, az, el, dist, x, y, z
      wav_dir/XXXX_L.wav   - left ear
      wav_dir/XXXX_R.wav   - right ear
    """
    wav_dir = Path(wav_dir)
    wav_dir.mkdir(parents=True, exist_ok=True)

    fe_int = int(fe)
    N_ext  = hrir_L.shape[1]

    with open(wav_dir / 'index.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['index', 'az_deg', 'el_deg', 'dist_m', 'x', 'y', 'z'])
        for i in range(N_ext):
            stem = f"{i + 1:04d}"
            wavfile.write(wav_dir / f"{stem}_L.wav", fe_int,
                          hrir_L[:, i].astype(np.float32))
            wavfile.write(wav_dir / f"{stem}_R.wav", fe_int,
                          hrir_R[:, i].astype(np.float32))
            writer.writerow([
                i + 1,
                f"{az[i]:.4f}", f"{el[i]:.4f}", f"{r[i]:.4f}",
                f"{pext[i, 0]:.6f}", f"{pext[i, 1]:.6f}", f"{pext[i, 2]:.6f}",
            ])

    print(f"WAV: {N_ext * 2} float32 files + index.csv -> {wav_dir}/")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Convert HRTF spectra (BEM/MLFMM) to windowed HRIRs.')
    parser.add_argument('config', help='Left-ear config file (config_L.txt)')
    parser.add_argument('--result-R', default=None,
                        help='Right-ear result file stem (overrides auto-convention). '
                             'Auto: replaces last _L with _R in result_file, '
                             'or last _G with _D (legacy French naming).')
    parser.add_argument('--no-equalization', action='store_true',
                        help='Disable diffuse-field equalisation')
    parser.add_argument('--n-out', type=int, default=256,
                        help='HRIR length in samples (default: 256)')
    parser.add_argument('--fade-in', type=int, default=4,
                        help='Asymmetric fade-in in samples (default: 4)')
    parser.add_argument('--fade-out', type=int, default=8,
                        help='Asymmetric fade-out in samples (default: 8)')
    parser.add_argument('--delay', type=int, default=32,
                        help='Pure delay added before the onset, in samples (default: 32 = '
                             '0.67 ms at 48 kHz; the ipsilateral ear leads t = 0 by up to ~0.3 ms)')
    parser.add_argument('--tukey-alpha', type=float, default=None,
                        help='If given, replaces asymmetric window with symmetric Tukey '
                             '[0=rect, 1=Hann] (e.g. 0.1)')
    parser.add_argument('--output', '-o', default=None,
                        help='Output .npz file (default: hrir_<stem>.npz)')
    parser.add_argument('--wav-dir', default=None,
                        help='Output directory for WAV files + index CSV (optional)')
    parser.add_argument('--weights', default=None,
                        help='quadrature weights of the directions (N_ext values, sum=1), '
                             'for the diffuse-field equalisation: weights_lebedev_N.txt or '
                             'weights_fibonacci_N.txt (generate_lebedev.py / generate_fibonacci.py), '
                             'or any file. Default: uniform weights (arithmetic mean).')
    parser.add_argument('--no-normalize', action='store_true',
                        help='Disable global peak normalisation')
    args = parser.parse_args()

    # --- Config (left ear) ---
    params, base = parse_config(args.config)

    N_freq = params['N_freq']
    fe     = params['fe']
    N_time = 2 * N_freq
    # f_min from config is only meaningful for single-frequency runs (N_freq=1)
    f_min  = params['f_min'] if N_freq == 1 else fe / N_time

    print(f"fe={fe} Hz | N_freq={N_freq} | f_min={f_min:.4f} Hz | "
          f"df={fe/N_time:.4f} Hz")
    print(f"N_time={N_time} ({N_time/fe*1e3:.2f} ms)  ->  "
          f"N_out={args.n_out} ({args.n_out/fe*1e3:.2f} ms)")

    # --- Directions (pext from left-ear config) ---
    pext  = load_pext(base / params['pext_file'])
    N_ext = pext.shape[0]
    az, el, r = cartesian_to_spherical(pext)
    print(f"N_ext={N_ext} directions")

    # --- Result file stems ---
    stem_L = params['result_file']
    stem_R = args.result_R if args.result_R else derive_result_r(stem_L)
    print(f"Result files: L={stem_L}.txt  R={stem_R}.txt")

    # --- Load HRTFs ---
    H_L = load_hrtf(base / (stem_L + '.txt'), N_freq, N_ext)
    H_R = load_hrtf(base / (stem_R + '.txt'), N_freq, N_ext)
    print(f"HRTFs loaded (L and R): {H_L.shape}")

    # --- quadrature weights of the directions (optional) ---
    weights = None
    if args.weights:
        weights = np.loadtxt(args.weights)
        print(f"Quadrature weights loaded: {args.weights} ({len(weights)} values, "
              f"sum={weights.sum():.6f})")
    else:
        print("Weights: uniform arithmetic mean (no --weights given)")

    # --- Diffuse-field equalisation ---
    equalized = not args.no_equalization
    if equalized:
        DFR = diffuse_field_response(H_L, H_R, weights)
        DFR = np.maximum(DFR, 1e-12)
        H_L = H_L / DFR[:, np.newaxis]
        H_R = H_R / DFR[:, np.newaxis]
        print("Diffuse-field equalisation applied")
    else:
        print("No equalisation (--no-equalization)")

    # --- HRTF -> HRIR ---
    if args.tukey_alpha is not None:
        print(f"Symmetric Tukey window alpha={args.tukey_alpha}, N_out={args.n_out}")
    else:
        print(f"Asymmetric window: fade-in={args.fade_in} fade-out={args.fade_out} "
              f"samples, N_out={args.n_out}")
    print(f"Delay: {args.delay} samples ({1e3 * args.delay / fe:.3f} ms)")
    hrir_L = hrtf_to_hrir(H_L, fe, f_min, args.n_out,
                           args.fade_in, args.fade_out, args.tukey_alpha, args.delay)
    hrir_R = hrtf_to_hrir(H_R, fe, f_min, args.n_out,
                           args.fade_in, args.fade_out, args.tukey_alpha, args.delay)

    # --- Global normalisation (same factor for L+R, all directions) ---
    if not args.no_normalize:
        peak = max(np.abs(hrir_L).max(), np.abs(hrir_R).max())
        if peak > 0:
            scale = 0.99 / peak
            hrir_L *= scale
            hrir_R *= scale
            print(f"Global normalisation: factor={scale:.6f} (peak={peak:.6f})")

    # --- NPZ ---
    if args.output is None:
        out_file = base / f"hrir_{stem_L}.npz"
    else:
        out_file = Path(args.output)

    save_npz(out_file, hrir_L, hrir_R, pext, az, el, r, fe, equalized)

    # --- WAV (optional) ---
    if args.wav_dir:
        save_wav_dir(args.wav_dir, hrir_L, hrir_R, pext, az, el, r, fe)


if __name__ == '__main__':
    main()
