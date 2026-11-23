# Python tools

Preparation of a computation (mesh checks, sources, directions, frequencies)
and post-processing (interpolation, HRIR, SOFA). The step-by-step guide is
[`../MANUAL.md`](../MANUAL.md); [`../run_subject.sh`](../run_subject.sh) chains
the solver and the post-processing for one subject.

```bash
pip install -r requirements.txt          # numpy, scipy, sofar
```

| Script | Step | Needs |
|--------|------|-------|
| [`diag_msh.py`](#diag_mshpy) | check a mesh | numpy |
| [`msh_flip.py`](#msh_flippy) | reverse the normals | — |
| [`source_oreille.py`](#source_oreillepy) | source triangle (Ray) and point (Diff) | numpy |
| [`generate_lebedev.py`](#generate_lebedevpy--generate_fibonaccipy) | directions + quadrature weights | numpy |
| [`generate_fibonacci.py`](#generate_lebedevpy--generate_fibonaccipy) | directions, uniform weights | numpy |
| [`frequency_grid.py`](#frequency_gridpy) | regular or mixed frequency list | numpy |
| [`interpolate_hrtf.py`](#interpolate_hrtfpy) | mixed grid → 512 frequencies | numpy |
| [`hrtf2hrir.py`](#hrtf2hrirpy) | HRTFs → windowed HRIRs (NPZ, WAV) | numpy, scipy |
| [`npz2sofa.py`](#npz2sofapy) | HRIRs → SOFA file | sofar |
| `maillage_msh.py` | library: mesh reading, distances, ray intersections | numpy |

```
head mesh ──► diag_msh.py (msh_flip.py) ──► source_oreille.py ──► case files L, R
generate_lebedev.py or generate_fibonacci.py ──► pext_<grid>_N.txt (case files), weights_<grid>_N.txt
frequency_grid.py --mixed --config ──► frequency line of the case files
solver L, solver R ──► <stem>.txt ──► interpolate_hrtf.py ──► <stem>_512.txt
──► hrtf2hrir.py ──► hrir_<subject>.npz (+ WAV) ──► npz2sofa.py ──► <subject>.sofa
```

Conventions: x front, y left, z up (right ear at y < 0); azimuth
`atan2(y, x)` (0° front, 90° left), elevation `asin(z / r)`; in SOFA, azimuth
in [0°, 360°[. Time dependence e^{−iωt} in the solver results.

Mesh formats read by all scripts: gmsh 2.2 ASCII or binary (`.msh`), STL binary
or ASCII (`.stl`, millimetres converted to metres when the extent exceeds 2,
as in the solvers).

---

## `diag_msh.py`

```bash
python diag_msh.py head.msh            # or head.stl
```

Checks, before any computation: flat triangles only; closed surface (every
edge shared by two triangles); extent and units; edge lengths (the octree is
limited to 10 levels: avoid edges below ~0.5 mm); minimum angle of the
triangles; orientation of the normals (signed volume: they must point out of
the head, into the air); the Oy axis crossing the mesh on both sides with the
origin inside the head. Ends with `summary: no problem found` (exit code 0) or
the list of problems (exit code 1). `--tube` adds the checks of the tube solver.

## `msh_flip.py`

```bash
python msh_flip.py head.msh head_ok.msh      # gmsh 2.2 / 4 ASCII; 2.2 binary -> 2.2 ASCII
python msh_flip.py head.stl head_ok.stl      # STL binary or ASCII
```

Reverses all surface triangles. The order of the triangles is kept, so a Ray
source triangle number stays valid; the rest of the file is copied unchanged.

## `source_oreille.py`

```bash
python source_oreille.py head.msh --ear R                     # right ear (L: left)
python source_oreille.py head.msh --check -0.0005 -0.0870 0   # check a given Diff point
```

- **Ray**: the source triangle (number from 1, in the order of the file): the
  first triangle crossed by the Oy axis from the head centre towards the ear.
- **Diff**: the monopole position (m): the point of the air closest to the
  canal entrance at ≥ `--dmin` + `--margin` mm from the mesh (1 + 1 mm by
  default), then checked (exact distance, winding number ≈ 0: in the air).

| Option | Default | Meaning |
|--------|---------|---------|
| `--ear` | `R` | `R` right (y < 0), `L` left (y > 0) |
| `--dmin` | 1 | minimum distance of the Diff point to the mesh (mm); beyond ~3.5 mm a notch c/(4d) enters the band |
| `--margin` | 1 | margin above `--dmin` (mm) |
| `--step` | 0.5 | search step (mm) |
| `--rlat` | 10 | half-width of the lateral search (mm) |
| `--depth` | 15 | maximum depth from the canal entrance (mm) |
| `--check X Y Z` | — | only check the given Diff point (exit code 1 if it fails) |

For several meshes of the same head, keep the same Diff point and `--check` it
on each mesh. (Former option names `--oreille D/G`, `--verifier`, `--marge`,
`--pas`, `--prof` are still accepted.)

## `generate_lebedev.py`, `generate_fibonacci.py`

```bash
python generate_lebedev.py --n 2702 --dist 3.0 --out-dir run/
python generate_fibonacci.py --n 2000 --dist 3.0 --out-dir run/
```

Write `pext_<grid>_NNNN.txt` (N on the first line, then N lines `x y z` in
metres: the external points file of the case file) and `weights_<grid>_NNNN.txt`
(quadrature weights, sum 1, for the diffuse-field equalisation of
`hrtf2hrir.py`). Lebedev grids exist for 32 sizes only (6 to 5810 points; the
closest is used); Fibonacci accepts any N, with uniform weights.

| Option | Default | Meaning |
|--------|---------|---------|
| `--n` | 2702 | number of points |
| `--dist` | 3.0 | distance of the points (m) |
| `--out-dir` | `.` | output directory |

Lebedev: port of the MATLAB version of R. Parrish, from the Fortran of
C. van Wüllen and the C of D. Laikov. Reference: V. I. Lebedev, D. N. Laikov,
Doklady Mathematics 59(3), 477–481, 1999.

## `frequency_grid.py`

```bash
python frequency_grid.py --mixed --config config_S01_L.txt   # writes NF and the list into the case file
python frequency_grid.py --mixed                             # prints NF, then the list
```

Mixed grid: the regular grid `k·fe/(2N)` up to `--f-lin`, then 1/`--oct` octave
up to fe/2: 335 frequencies instead of 512 with the defaults (−35 % of
computing time). Interpolation error back to 512 frequencies (99th percentile,
head and torso): 0.004 dB at 10–16 kHz, 0.04 dB at 16–24 kHz, 0.06 dB on the
contralateral side; larger only at the bottom of deep notches.

| Option | Default | Meaning |
|--------|---------|---------|
| `--mixed` | — | mixed grid (without it: the regular grid) |
| `--fe` | 48000 | sampling rate (Hz) |
| `--n` | 512 | number of frequencies of the full grid |
| `--f-lin` | 10000 | end of the regular part (Hz) |
| `--oct` | 96 | logarithmic step, 1/oct octave |
| `--config` | — | case file to update |

## `interpolate_hrtf.py`

```bash
python interpolate_hrtf.py config_S01_L.txt          # --n 512 by default
```

Rebuilds the regular grid from a run on a mixed (or any) frequency list: the
geometric delay of each direction (source → external point) is removed, the
real and imaginary parts are interpolated by natural cubic splines, the delay is
put back; computed frequencies are copied unchanged. Writes `<stem>_512.txt`
(solver format; iterations 0 on interpolated frequencies), `<stem>_512.npz`
(`f`, `H`, `computed`, `Pext`) and `config_<stem>_512.txt` (case file of the
full grid, for `hrtf2hrir.py`). `--source X Y Z` overrides the source point
(default: Diff point, or centroid of the Ray triangle).

## `hrtf2hrir.py`

```bash
python hrtf2hrir.py config_S01_L_512.txt --weights weights_lebedev_2702.txt -o hrir_S01.npz
```

Takes the **left-ear** case file (regular grid); the right-ear results are
found by replacing the last `_L` of the stem by `_R` (or `_G` by `_D`), or
given with `--result-R`. Steps: conjugation (solver e^{−iωt} → e^{+iωt}),
pure delay, diffuse-field equalisation
`DFR(f) = sqrt(Σ w_i (|H_L|² + |H_R|²)/2)`, inverse FFT (DC: real part at
f_min; Nyquist: real part), window, global peak normalisation.

| Option | Default | Meaning |
|--------|---------|---------|
| `--result-R` | from the stem | stem of the right-ear results |
| `--weights` | uniform | quadrature weights of the directions: `weights_lebedev_N.txt` or `weights_fibonacci_N.txt` |
| `--no-equalization` | — | no diffuse-field equalisation |
| `--n-out` | 256 | HRIR length (samples) |
| `--delay` | 32 | pure delay added (samples, 0.67 ms at 48 kHz): the HRTF is referred to the head centre, so the ipsilateral ear receives the sound up to ~0.3 ms before t = 0 |
| `--fade-in` | 4 | cosine fade-in (samples) |
| `--fade-out` | 8 | cosine fade-out (samples) |
| `--tukey-alpha` | — | symmetric Tukey window instead of the two fades |
| `--no-normalize` | — | no global peak normalisation (factor 0.99 / max) |
| `--wav-dir` | — | also one float32 WAV file per direction and ear (`NNNN_L.wav`, `NNNN_R.wav`) and `index.csv` (index, azimuth, elevation, distance, x, y, z) |
| `-o` | `hrir_<stem>.npz` | output file: `hrir_L`, `hrir_R` (samples × directions), `pext`, `az`, `el`, `r`, `fe`, `equalized` |

## `npz2sofa.py`

```bash
python npz2sofa.py hrir_S01.npz S01.sofa --subject S01
```

SOFA file AES69-2022, convention SimpleFreeFieldHRIR, `Data.IR` of shape
(directions, 2 ears: 0 left, 1 right, samples). Metadata are read from
`sofa_meta.txt` in the current directory, else from `python/sofa_meta.txt`
(fill it in once: `title`, `author`, `org`, `license`, `db`, `comment`, as
`key = value`); the options `--title`, `--author`, `--org`, `--license`,
`--comment`, `--db`, `--subject`, `--meta FILE` override them.
