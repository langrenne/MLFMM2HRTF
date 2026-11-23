#!/usr/bin/env python3
# Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
# SPDX-License-Identifier: MIT
"""
frequency_grid.py
Frequency list for an HRTF run, regular or MIXED (2026-10-02).

  regular : f[k] = k * fe/(2*N),  k = 1..N           (N = 512 for fe = 48 kHz: 46.875 Hz step)
  mixed   : the regular grid up to --f-lin (default 10 kHz), then a logarithmic grid
            at 1/--oct octave (default 1/96) up to fe/2, fe/2 included.
            With the defaults: 335 frequencies instead of 512 (-35 % of computing time).
            The missing frequencies are rebuilt afterwards by interpolate_hrtf.py.

Interpolation error of the mixed grid (2026-11-18, head and torso HT07_20, 337 k elements,
precise setting, 335 -> 512 frequencies against a run on the 512 frequencies): p99 0.004 dB at 10-16 kHz, 0.04 dB at 16-24 kHz, 0.06 dB contralateral;
worst case 4 dB at the bottom of a notch 58 dB below the median level (narrower than the step).

Usage:
  python frequency_grid.py --mixed                      -> prints NF, then the list
  python frequency_grid.py --mixed --config config.txt  -> writes NF and the list into
                                                           the case file (lines 2 and 6)
Options: --fe 48000 --n 512 --f-lin 10000 --oct 96
"""
import argparse
import math
import numpy as np


def regular_grid(fe, n):
    return np.arange(1, n + 1) * fe / (2.0 * n)


def mixed_grid(fe, n, f_lin, n_oct):
    df = fe / (2.0 * n)
    k_lin = int(math.floor(f_lin / df))
    lin = np.arange(1, k_lin + 1) * df
    f0 = lin[-1]
    m = int(math.floor(n_oct * math.log2((fe / 2.0) / f0)))
    log = f0 * 2.0 ** (np.arange(1, m + 1) / n_oct)
    grid = np.concatenate([lin, log])
    if grid[-1] < fe / 2.0 - 1e-6:
        grid = np.append(grid, fe / 2.0)
    return grid


def write_into_config(config, freqs):
    """Replace the number of frequencies (2nd line) and the frequency line (6th line)."""
    with open(config) as f:
        lines = f.read().splitlines()
    lines[1] = str(len(freqs))
    lines[5] = ' '.join(f'{x:.6f}' for x in freqs)
    if len(lines[5]) > 20000:
        raise SystemExit('frequency line longer than 20 000 characters (solver limit)')
    with open(config, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def main():
    p = argparse.ArgumentParser(description='Regular or mixed frequency list for the HRTF solver.')
    p.add_argument('--mixed', action='store_true', help='mixed grid (regular, then 1/oct octave)')
    p.add_argument('--fe', type=float, default=48000.0, help='sampling rate (Hz)')
    p.add_argument('--n', type=int, default=512, help='number of frequencies of the full grid')
    p.add_argument('--f-lin', type=float, default=10000.0, help='end of the regular part (Hz)')
    p.add_argument('--oct', type=int, default=96, help='logarithmic step: 1/oct octave')
    p.add_argument('--config', default=None, help='case file to update (lines 2 and 6)')
    a = p.parse_args()
    f = mixed_grid(a.fe, a.n, a.f_lin, a.oct) if a.mixed else regular_grid(a.fe, a.n)
    if a.config:
        write_into_config(a.config, f)
        print(f'{a.config}: {len(f)} frequencies, {f[0]:.3f} to {f[-1]:.3f} Hz')
    else:
        print(len(f))
        print(' '.join(f'{x:.6f}' for x in f))


if __name__ == '__main__':
    main()
