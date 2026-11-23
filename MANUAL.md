# User manual — computing individual HRTFs

This manual walks through a complete computation, from a head mesh to HRTFs on
512 frequencies, ready for MATLAB, Python or the SOFA export. The line-by-line
reference of the input files is [`config.md`](config.md); building is described
in [`README.md`](README.md).

```
head mesh (.msh / .stl) ──► checks ──► source (Ray triangle or Diff point) ──► direction grid
      ──► case file (+ optional technical file) ──► RAM estimate ──► solver (Ray or Diff)
      ──► <stem>.txt  ──► (mixed grid: interpolation to 512 frequencies)
      ──► MATLAB / Python, or hrtf2hrir.py ──► npz2sofa.py
```

---

## 1. Two solvers, one ear at a time

Both solvers use acoustic **reciprocity**: the source is placed at the ear and
the pressure is computed at the external points (the directions). The boundary
integral equation is solved with the CHI formulation (Burton–Miller equations
written at points just inside the head), which avoids the spurious resonances
of the interior problem.

| Program | Source | Use |
|---------|--------|-----|
| `main_HRTF_Ray` | one vibrating **triangle** of the mesh at the entrance of the (blocked) ear canal | classical reciprocal HRTF |
| `main_HRTF_Diff` | a **monopole** in the air, ~2 mm in front of the canal entrance | the same source point can be kept for several meshes of the same head: preferred to compare meshes |

Both solvers accept a **skin impedance**: a locally reacting surface, uniform
over the head and frequency dependent (two reflection coefficients after
Ghorbal et al., or a file of reflection, admittance or impedance values; see
`config.md`), given on the last line of the case file; without that line, the
skin is rigid. With Ray, the skin impedance requires the default `kreuzer`
source (the `gumerov` source model does not admit it).

One run computes one ear. For both ears, run twice with two case files
(convention `<subject>_L` / `<subject>_R`).

## 2. Build and environment

```bash
git clone --recursive <repository-url> MLFMM2HRTF && cd MLFMM2HRTF
make                            # builds main_HRTF_Ray and main_HRTF_Diff (options: README.md)
export OPENBLAS_NUM_THREADS=1   # mandatory (see below)
export OMP_NUM_THREADS=24       # MAXIMUM number of threads: the number of physical cores
```

The solver is parallel over the frequencies. `OMP_NUM_THREADS` is an upper
bound: the solver itself decides how many frequencies run at the same time,
from the available RAM (section 7). `OPENBLAS_NUM_THREADS=1` is required:
otherwise OpenBLAS starts its own threads inside each frequency and the cores
are oversubscribed.

Tested on Linux (Ubuntu, gfortran, OpenMP, LAPACK/BLAS). The RAM estimate reads
`/proc/meminfo`: on Windows, use WSL2.

## 3. Preparing the mesh

Requirements:

- gmsh **2.2** format (`.msh`, ASCII or binary, coordinates in **metres**), or
  **STL** (`.stl`, binary or ASCII; if the mesh is larger than 2 units it is
  taken as millimetres and converted to metres); flat triangles only;
- one closed surface (head, possibly with torso), normals pointing **out of the
  head, into the air**;
- the head **aligned on the Oy axis**: Oy goes through both ear canal entrances,
  right ear on the y < 0 side, left ear on y > 0 (needed for the automatic
  source placement); x points to the front, z up;
- element size adapted to the highest frequency: about six elements per
  wavelength, i.e. edges of about **2.4 mm at 24 kHz**. Use finer elements
  (0.5–1 mm) in the pinna of the computed ear and up to ~2 mm elsewhere
  (graded mesh). Coarser elements far from the ear make the computation faster
  but limit the accuracy of the contralateral HRTF at high frequencies.
- smallest elements about **0.5 mm** (the octree is limited to 10 levels; for a
  head with torso, level 10 has cells of about 0.5 mm).

Checks, before any computation:

```bash
python python/diag_msh.py head.msh               # triangles, closed, units, edges, quality, normals, ear axis
python python/msh_flip.py head.msh head_ok.msh   # if the normals point into the head (.msh or .stl)
```

`diag_msh.py` ends with `summary: no problem found` or the list of problems
(exit code 1). `msh_flip.py` keeps the order of the triangles, so a Ray source
triangle number stays valid.

A few degenerate triangles (minimum angle < 1°, edges of a few hundredths of a
millimetre) do not stop the solver but are best removed (merge close vertices).

## 4. Source

```bash
python python/source_oreille.py head.msh --ear R        # right ear (L: left)
```

prints

- **Ray**: the source triangle number (from 1, in the order of the file), i.e.
  the first triangle crossed by the Oy axis from the head centre towards the ear;
- **Diff**: the monopole position (m): the point of the air closest to the
  canal entrance at ≥ 2 mm from the mesh (`--dmin 1` mm plus a 1 mm margin).

To compare several meshes of the same head with Diff, keep **the same point**
and check it on each mesh:

```bash
python python/source_oreille.py other.msh --check -0.00050 -0.06502 0.00100
#   -> distance to the mesh and "in the air: True"
```

A Diff source at a distance d from the skin creates an interference notch near
c/(4d): keep d below about 3.5 mm so that it stays above 24 kHz.

**Choosing the source by hand.** The script is a convenience: any triangle
(Ray) or point in the air (Diff) can be used, for instance picked in a mesh
viewer. The Ray triangle number counts **from 1**, in the order of the
triangles in the file; check the numbering of the viewer: MeshLab, for
instance, numbers the faces from 0, so **add 1** to the index it shows. Check a
hand-picked Diff point with `--check x y z` as above.

## 5. Direction grid (external points)

```bash
python python/generate_lebedev.py   --n 2702 --dist 1.5 --out-dir grid/
python python/generate_fibonacci.py --n 2600 --dist 1.5 --out-dir grid/
```

Pext file format: first line = number of points, then one point per line,
`x y z` in metres. Each script also writes the quadrature weights of its grid
(`weights_lebedev_N.txt`, `weights_fibonacci_N.txt`), used later for the
diffuse-field equalisation (section 9):

| Grid | Sizes | Weights |
|------|-------|---------|
| Lebedev | 32 sizes only (6 … 5810; the closest one is used) | exact quadrature rule, unequal weights |
| Fibonacci | any number of points | uniform, 1/N (nearly uniform covering of the sphere) |

Any other set of directions can be used; give its weights to `hrtf2hrir.py`
(`--weights`), or none for uniform weights.

## 6. Input files

Each run reads a **case file** (`config_RAY.txt`, `config_DIFF.txt`): mesh,
directions, frequencies, source, output stem, speed of sound, and on its last
line the name of an optional **technical file**. Line-by-line reference:
[`config.md`](config.md).

File names are relative to the directory where the solver is run (e.g.
`grid/pext_lebedev_2702.txt`), or absolute.

### Settings: fast or precise

The last line of the case file is either `-` or the name of a technical file
(`key = value`, see [`config.md`](config.md), part 2):

| Last line | Settings |
|-----------|----------|
| `-` | **fast** settings (default) |
| a file containing `reglage = precis` | **precise** settings |

| Setting | Target solver error (FMM + GMRES) | Time |
|---------|------|------|
| fast | below 0.1 dB — for production | 1 |
| precise | below 0.05 dB — reference to judge other settings | × 1.2–1.4 |

Both targets assume a mesh valid up to 24 kHz (elements of at most ~2.4 mm);
they were checked against the exact solution on spheres (99th percentile over
the directions). On a coarser mesh, the error due to the mesh dominates at high
frequencies, whatever the setting. These deviations are far below what the ear
perceives (about 1 dB); they are criteria to compare numerical settings.

### Frequencies

Case file, lines 2 (number of frequencies) and 6 (frequency line):

- **regular grid**: `NF` = 512 and a single value on the frequency line →
  f[k] = k · fe / (2 · NF), 46.875 Hz to 24 kHz for fe = 48 kHz;
- **list**: `NF` values on the frequency line (any frequencies, e.g. a few test
  frequencies);
- **mixed grid** (about half the computing time): regular step up to 10 kHz,
  then 1/96 octave, 335 frequencies, rebuilt on the 512 regular frequencies
  afterwards (section 9):

```bash
python python/frequency_grid.py --mixed --config config_S01_R.txt   # writes lines 2 and 6
```

## 7. RAM: estimate before running

The memory of a frequency grows with the frequency and with the mesh size; on
a large mesh, RAM, not the number of cores, limits the number of frequencies
computed at the same time. Put `estimer = 1` in a technical file: the solver
reads the mesh, prints the RAM per frequency and the number of threads
possible, then **stops without computing**.

During the computation, the solver starts a frequency only when the memory
reserved by the running ones leaves room for it (budget: 85 % of the RAM
available at start, or `ram_max` GB). High and low frequencies are
interleaved; each thread then chains neighbouring frequencies, the solution of
one being the initial guess of the next.

Orders of magnitude (24 cores, 62 GB, idle machine):

| Mesh | Elements | RAM per frequency, 47 Hz → 24 kHz | Frequencies at the same time |
|------|----------|----------------------|------------------|
| head + torso, 0.7–7 mm | 48 000 | 0.8 → 2.5 GB | 24 |
| head + torso, 0.7–3 mm | 158 000 | 2.4 → 5.4 GB | 12 |
| head + torso, 0.7–2 mm | 337 000 | 3.2 → 7.0 GB | 10 |
| head + torso, 0.7–1.6 mm | 484 000 | 5.2 → 9.7 GB | 6 |

Computing time of one ear, fast setting (precise in parentheses), same machine:

| Mesh | Elements | Mixed grid, 335 frequencies | Regular grid, 512 frequencies |
|------|----------|------|------|
| head + torso, 0.7–7 mm | 48 000 | 1 h 18 (1 h 36) | 2 h 24 (3 h 09) |
| head + torso, 0.7–3 mm | 158 000 | 1 h 53 (2 h 19) | 3 h 24 (4 h 24) |
| head + torso, 0.7–2 mm | 337 000 | 3 h 16 (4 h 20) | 5 h 38 (7 h 29) |
| head + torso, 0.7–1.6 mm | 484 000 | 6 h 12 (8 h 02) | — |
| head only, 0.5–1 mm | 409 000 | 3 h 06 (4 h 24) | — |

## 8. Running

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=24
nohup ./main_HRTF_Ray config_S01_R.txt > S01_R.out 2>&1 &
```

Run **one long computation at a time**: two solvers in parallel compete for the
RAM and the memory bandwidth (for both ears, run the second after the first).
**Results are written only at the end** of the run. The log file (named in
the case file) shows the progress; with verbose 1 (production):

```
 MLFMM2HRTF  main_HRTF_Ray  2026-11-23 10:17

 case file    : config_S01_R.txt
 threads      : 24 (OpenMP)
 frequencies  : 335 (list), 46.875 to 24000.000 Hz
 c            : 340.00 m/s
 mesh         : S01_R.msh: 366112 triangles, 183058 nodes, edges 0.62 to 2.45 mm
 technical    : - (none: default settings)
 setting      : fast
 octree       : adaptive, 8 levels, box 495.4 mm, leaf cell 1.935 mm
 skin         : rigid
 source       : triangle 214178 (kreuzer: vibrating triangle, area 0.2718 mm2)
 directions   : pext_lebedev_2702.txt: 2702 points
 normalisation: yes (free field at the head centre)
 preparation  : 11.5 s
 RAM/thread   : 2.81 GB (46.9 Hz) to 7.38 GB (23919.1 Hz)
 memory       : shared 5.35 GB, available 50.2 GB, budget of the threads 42.7 GB
 parallel     : 24 threads at most: 5 at the largest frequency, 15 at the smallest

 frequencies as they finish (335):
       #   freq (Hz)   iter   time (s)   done
       1      46.875     30      101.5   1/335
     322   21933.875     35      563.0   2/335
 ...

 done         : 335 frequencies in 4 h 00 min
 iterations   : 25 to 48 (iter_max 150)
 not converged: none
 results      : S01_R.txt
```

`#` is the number of the frequency in the case file; frequencies finish out of
order (several run in parallel). A frequency that reaches `iter_max` is marked
and listed in the summary, with a `WARNING` line giving its backward error.
`grep -c "/335" S01_R.log` counts the frequencies done.

### Whole chain for one subject: `run_subject.sh`

```bash
./run_subject.sh config_S01_L.txt config_S01_R.txt          # ./run_subject.sh -h for the options
```

runs, in the directory of the case files: the solver for the left ear, then
the right ear (Ray or Diff, chosen from the source line of the case file; an
ear whose result file already exists is skipped), the interpolation to 512
frequencies if the frequency line is a list (mixed grid), `hrtf2hrir.py` and
`npz2sofa.py` (section 9), and writes `S01.sofa`. The quadrature weights are
found from the name of the external points file (`pext_lebedev_N.txt` →
`weights_lebedev_N.txt`, `pext_fibonacci_N.txt` → `weights_fibonacci_N.txt`). Main options: `--wav` (one WAV file
per direction and ear, e.g. for a convolution engine), `--threads N`,
`--post-only` (results already computed), `--force`, `--wait` (wait until no
other solver runs). Progress is written to `run_S01.log`.

## 9. Results

### Result file `<stem>.txt`

ASCII, one line per (frequency, direction); frequencies in the outer loop, in
the order of the frequency list; directions in the order of the Pext file:

```
Re(H)   Im(H)   GMRES_iterations   total_time_s
```

`NF × N_ext` lines. H is the pressure normalised by the free-field pressure of
the same source at the head centre (head absent), if the normalisation line of
the case file is 1.

### Mixed grid → 512 frequencies

```bash
python python/interpolate_hrtf.py config_S01_R.txt          # --n 512 by default
```

The geometric delay of each direction is removed, the real and imaginary parts
are interpolated by cubic splines, and the delay is put back; computed
frequencies are copied unchanged. Accuracy, interpolated 335 → 512 frequencies
against a run on the 512 frequencies (head and torso, 337,000 elements): 99th
percentile 0.004 dB at 10–16 kHz, 0.04 dB at 16–24 kHz, 0.06 dB on the
contralateral side. The largest deviations sit at the bottom of deep notches
(worst case 4 dB, in a notch 58 dB below the median level, narrower than the
frequency step: the regular grid does not resolve it either). Outputs, next to
the case file:

| File | Content |
|------|---------|
| `<stem>_512.txt` | same format as the solver (iterations and time = 0 on interpolated frequencies) |
| `<stem>_512.npz` | `f` (512), `H` (512 × N_ext, complex), `computed` (512, bool), `Pext` |
| `config_<stem>_512.txt` | case file of the full grid, for `hrtf2hrir.py` |

### Reading in MATLAB

```matlab
d    = load('S01_R.txt');                                  % (NF*Next) x 4
P    = readmatrix('pext_lebedev_2702.txt', 'NumHeaderLines', 1);   % Next x 3 (m)
Next = size(P, 1);
NF   = size(d, 1) / Next;
H    = reshape(complex(d(:,1), d(:,2)), Next, NF).';       % NF x Next
it   = d(1:Next:end, 3);                                   % GMRES iterations per frequency
f    = (1:NF)' * 48000 / (2*NF);                           % regular grid (or the list of the case file)
[az, el] = cart2sph(P(:,1), P(:,2), P(:,3));               % x front, y left, z up
```

### Reading in Python

```python
import numpy as np
P = np.loadtxt('pext_lebedev_2702.txt', skiprows=1)
d = np.loadtxt('S01_R.txt')
H = (d[:, 0] + 1j * d[:, 1]).reshape(-1, len(P))   # NF x Next
# or, after interpolation:
z = np.load('S01_R_512.npz'); f, H = z['f'], z['H']
```

### HRIR and SOFA

```bash
python python/hrtf2hrir.py config_S01_L.txt --weights grid/weights_lebedev_2702.txt -o hrir_S01.npz
#                                         or --weights grid/weights_fibonacci_2600.txt
python python/npz2sofa.py hrir_S01.npz S01.sofa --subject S01
```

`hrtf2hrir.py` takes the **left-ear** case file and finds the right-ear results
itself. With a mixed grid, give it the `config_<stem>_512.txt` file written by
`interpolate_hrtf.py`. Options:

| Option | Default | Effect |
|--------|---------|--------|
| `--result-R` | last `_L` of the stem → `_R` | stem of the right-ear results |
| `--weights` | uniform | quadrature weights of the directions, for the diffuse-field equalisation: `weights_lebedev_N.txt` or `weights_fibonacci_N.txt` (section 5) |
| `--no-equalization` | — | no diffuse-field equalisation |
| `--n-out` | 256 | HRIR length (samples) |
| `--delay` | 32 | pure delay added (samples; 0.67 ms at 48 kHz): the HRTF is referred to the head centre, so the ipsilateral ear receives the sound up to ~0.3 ms before t = 0 |
| `--fade-in` | 4 | cosine fade-in (samples) |
| `--fade-out` | 8 | cosine fade-out (samples) |
| `--tukey-alpha` | — | symmetric Tukey window instead of the two fades |
| `--no-normalize` | — | no global peak normalisation |
| `--wav-dir` | — | also write one WAV file per direction and ear, with an index file |
| `-o` | `hrir_<stem>.npz` | output file |

`npz2sofa.py` writes a SOFA file (AES69-2022, SimpleFreeFieldHRIR). Fill in
`python/sofa_meta.txt` once (title, author, organisation, licence; a
`sofa_meta.txt` in the run directory takes precedence); the options
`--title`, `--author`, `--org`, `--license`, `--comment`, `--db` and
`--subject` override it for one file.

Spatial conventions: x front, y left, z up; azimuth `atan2(y, x)` (0° front,
90° left), elevation `asin(z / r)`; in SOFA, azimuth in [0°, 360°[.

## 10. Checks and pitfalls

- **Normals and source**: inverted normals or a source inside the head give a
  result without any error message → `diag_msh.py` and `source_oreille.py
  --check` before every new mesh.
- **Validity of the mesh**: the solver converges whatever the mesh; the
  accuracy at high frequency is set by the element size. Compare two meshes on
  a few frequencies (frequency list) to know the highest valid frequency.
- **Iterations**: a few tens per frequency is normal; a frequency that reaches
  `iter_max` (150) is not converged — check the mesh first.
- **Test run**: `cd test && ./run_test.sh` (2 cases, Ray and Diff, on a small
  sphere, about 15 s) checks that the build gives the reference results.
