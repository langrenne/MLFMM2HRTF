#!/usr/bin/env python3
# Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
# SPDX-License-Identifier: MIT
"""
interpolate_hrtf.py
Rebuild the full regular frequency grid (default 512 frequencies) from a run
computed on a MIXED grid (frequency_grid.py --mixed), 2026-10-02.

Method (validated on the m7 and m2 meshes, 2026-10-01/02, "method B"):
  1. remove the geometric delay of each direction, tau = (|P| - |P - S|) / c,
     S = source (Diff monopole, or centroid of the Ray source triangle);
     its sign is the one that leaves the slowest phase rotation with frequency;
  2. natural cubic spline on the real and imaginary parts of the residual;
  3. put the delay back.
Frequencies already computed are copied as they are (no interpolation).

Input : the case file of the mixed run (config_RAY / config_DIFF format) and its
        result file <stem>.txt (stem read from the case file).
Output: <stem>_<N>.txt   same format as the solver (N_ext lines per frequency:
                         Re, Im, GMRES iterations, time ; iterations = 0 and
                         time = 0 on interpolated frequencies) -> readable by
                         hrtf2hrir.py and by MATLAB (see MANUAL.md)
        <stem>_<N>.npz   f (N), H (N x N_ext complex), computed (N bool), Pext
        config_<stem>_<N>.txt  copy of the case file with NF = N and a single
                         frequency (regular grid), stem <stem>_<N>, for hrtf2hrir.py

Usage: python interpolate_hrtf.py config_mixed.txt [--n 512] [--source X Y Z]
"""
import argparse
import os
import re
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def cubic_spline(x, y):
    """Natural cubic spline (numpy only); y of shape (n, m). Returns a function of xi."""
    x = np.asarray(x, float)
    y = np.asarray(y)
    n = len(x)
    h = np.diff(x)
    A = np.zeros((n, n))
    r = np.zeros((n,) + y.shape[1:], dtype=y.dtype)
    A[0, 0] = A[-1, -1] = 1.0
    for i in range(1, n - 1):
        A[i, i - 1], A[i, i], A[i, i + 1] = h[i - 1], 2 * (h[i - 1] + h[i]), h[i]
        r[i] = 3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])
    c = np.linalg.solve(A, r)
    b = (y[1:] - y[:-1]) / h[:, None] - h[:, None] * (2 * c[:-1] + c[1:]) / 3
    d = (c[1:] - c[:-1]) / (3 * h[:, None])

    def f(xi):
        k = np.clip(np.searchsorted(x, xi) - 1, 0, n - 2)
        t = (xi - x[k])[:, None]
        return y[k] + b[k] * t + c[k] * t ** 2 + d[k] * t ** 3
    return f


def read_case(config):
    """Values of a case file (lines not starting with '!'), Ray or Diff.

    Order (config.md): NF, fe, frequencies, mesh, Pext, stem, source, normalisation,
    verbose, log, c, technical file, [skin impedance]. Legacy Diff order (before
    2026-11-20): the skin impedance is value 9, before the verbose (as in the solver,
    the new order is recognised by value 9 being a lone integer).
    Ray: the source is a triangle number; Diff: a point x y z.
    """
    with open(config) as f:
        raw = f.read().splitlines()
    val = [l.strip() for l in raw if l.strip() and not l.lstrip().startswith('!')]
    if len(val) not in (12, 13):
        raise SystemExit(f'{config}: {len(val)} values, expected 12 or 13 (see config.md)')
    code = 'ray' if len(val[6].split()) == 1 else 'diff'
    new_order = re.fullmatch(r'[+-]?\d+', val[8].split()[0] if val[8].split() else '') is not None
    c = float(val[10] if new_order else val[11])
    d = dict(code=code, raw=raw, nf=int(val[0]), fe=float(val[1]), freq=val[2], mesh=val[3],
             pext=val[4], stem=val[5], source=val[6], c=c)
    return d


def frequencies(case):
    v = [float(x) for x in case['freq'].replace('d', 'e').split()]
    if len(v) > 1:
        return np.array(v)
    if case['nf'] == 1:
        return np.array(v)
    return np.arange(1, case['nf'] + 1) * case['fe'] / (2.0 * case['nf'])


def source_point(case, here, forced):
    if forced is not None:
        return np.array(forced, float)
    v = case['source'].replace('d', 'e').split()
    if case['code'] == 'diff':
        return np.array([float(x) for x in v[:3]])
    from maillage_msh import lire_maillage          # Ray: centroid of the source triangle (.msh or .stl)
    P, T = lire_maillage(os.path.join(here, case['mesh']))
    return P[T[int(v[0]) - 1]].mean(axis=0)


def main():
    p = argparse.ArgumentParser(description='Mixed-grid HRTF run -> full regular grid.')
    p.add_argument('config', help='case file of the mixed run')
    p.add_argument('--n', type=int, default=512, help='number of frequencies of the full grid')
    p.add_argument('--source', type=float, nargs=3, default=None, help='source point (m), overrides the case file')
    a = p.parse_args()

    here = os.path.dirname(os.path.abspath(a.config))
    case = read_case(a.config)
    fc = frequencies(case)
    pext = np.loadtxt(os.path.join(here, case['pext']), skiprows=1)
    data = np.loadtxt(os.path.join(here, case['stem'] + '.txt'))
    nx = pext.shape[0]
    if data.shape[0] != len(fc) * nx:
        raise SystemExit(f'{case["stem"]}.txt: {data.shape[0]} lines, expected {len(fc)} x {nx}')
    H = (data[:, 0] + 1j * data[:, 1]).reshape(len(fc), nx)
    it = data[:, 2].reshape(len(fc), nx)[:, 0]

    fg = np.arange(1, a.n + 1) * case['fe'] / (2.0 * a.n)
    if fc[0] > fg[0] + 1e-6 or fc[-1] < fg[-1] - 1e-6:
        raise SystemExit(f'computed frequencies {fc[0]:.3f}-{fc[-1]:.3f} Hz do not cover {fg[0]:.3f}-{fg[-1]:.3f} Hz')

    S = source_point(case, here, a.source)
    tau = (np.linalg.norm(pext, axis=1) - np.linalg.norm(pext - S, axis=1)) / case['c']

    def residual(sgn):
        return H * np.exp(-sgn * 1j * 2 * np.pi * fc[:, None] * tau[None, :])
    rot = {s: np.median(np.abs(np.angle(residual(s)[1:] / residual(s)[:-1]))) for s in (+1, -1)}
    sgn = min(rot, key=rot.get)
    R = residual(sgn)
    Ri = cubic_spline(fc, R.real)(fg) + 1j * cubic_spline(fc, R.imag)(fg)
    Hg = Ri * np.exp(sgn * 1j * 2 * np.pi * fg[:, None] * tau[None, :])
    iters = np.zeros(a.n)
    computed = np.zeros(a.n, bool)
    for i, f in enumerate(fc):                      # computed frequencies: copied exactly
        j = np.argmin(np.abs(fg - f))
        if abs(fg[j] - f) < 1e-6:
            Hg[j] = H[i]
            iters[j] = it[i]
            computed[j] = True

    out = case['stem'] + f'_{a.n}'
    with open(os.path.join(here, out + '.txt'), 'w') as f:
        for j in range(a.n):
            for x in range(nx):
                f.write(f'{Hg[j, x].real:.16E} {Hg[j, x].imag:.16E} {int(iters[j]):d} 0.0\n')
    np.savez(os.path.join(here, out + '.npz'), f=fg, H=Hg, computed=computed, Pext=pext)

    raw = list(case['raw'])                         # case file of the full grid, for hrtf2hrir.py
    idx = [i for i, l in enumerate(raw) if l.strip() and not l.lstrip().startswith('!')]
    raw[idx[0]] = str(a.n)
    raw[idx[2]] = f'{fg[0]:.6f}'
    raw[idx[5]] = out
    with open(os.path.join(here, f'config_{out}.txt'), 'w') as f:
        f.write('\n'.join(raw) + '\n')

    print(f'{case["code"]}: {len(fc)} computed frequencies -> {a.n} ({computed.sum()} copied, '
          f'{a.n - computed.sum()} interpolated); |tau| max {np.abs(tau).max() * 1e6:.0f} us, sign {sgn:+d}')
    print(f'written: {out}.txt, {out}.npz, config_{out}.txt')


if __name__ == '__main__':
    main()
