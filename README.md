# MLFMM2HRTF

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24473755.svg)](https://doi.org/10.5281/zenodo.24473755)

Computation of individual head-related transfer functions (HRTFs) from a head
mesh, up to 24 kHz, with a broadband multilevel fast multipole method (MLFMM)
and the boundary element method (BEM). Meshes of several hundred thousand
elements, strongly graded (fine at the ear, coarser elsewhere), are computed on
a single workstation.

- **Reciprocity**: the source is at the ear, the pressure is computed at the
  external points (the directions); one run per ear.
- **CHI formulation** (Burton–Miller equations at points just inside the head):
  no spurious resonances, accurate at all frequencies.
- **Broadband MLFMM** with an adaptive octree: each element enters the FMM at
  the level that matches its size.
- **Two settings**: fast (default, solver error below 0.1 dB) and precise
  (below 0.05 dB), validated against the exact solution on spheres.
- **Memory-aware scheduling**: the solver estimates the RAM of each frequency
  and runs as many frequencies in parallel as the memory allows.
- Python tools: mesh checks, source placement, direction grids, mixed
  frequency grid and interpolation, HRIR windowing and diffuse-field
  equalisation, SOFA export (AES69-2022).

Two solvers differ by the source model:

| Program | Source |
|---------|--------|
| `main_HRTF_Ray` | one vibrating triangle of the mesh at the entrance of the blocked ear canal |
| `main_HRTF_Diff` | a monopole in the air ~2 mm from the canal entrance |

Both accept a skin impedance (locally reacting surface, uniform over the head,
frequency dependent: reflection coefficient, admittance or impedance); for Ray,
with the default `kreuzer` source only.

## Documentation

- [`MANUAL.md`](MANUAL.md): step-by-step guide, from the mesh to the SOFA file
  (mesh, source, directions, settings, RAM, running, results, MATLAB/Python).
- [`config.md`](config.md): reference of the input files. Only one is required,
  the **case file** (`config_RAY.txt` / `config_DIFF.txt`, copied and filled in
  for each ear: mesh, directions, frequencies, source, output). Its last line
  names an optional **technical file** (`key = value`, e.g. `reglage = precis`
  or `estimer = 1`), or is `-`: no technical file, default fast settings.
- [`python/README.md`](python/README.md): the Python tools, one by one.

## Dependencies

- gfortran (Fortran 2008) with OpenMP;
- LAPACK and BLAS, development packages: `liblapack-dev libblas-dev` or
  `libopenblas-dev` (Debian, Ubuntu), `lapack-devel blas-devel` or
  `openblas-devel` (Fedora);
- Linux (the RAM estimate reads `/proc/meminfo`; on Windows, use WSL2).
- Python 3 with `numpy`, `scipy`, `sofar` (`pip install -r python/requirements.txt`).

## Building

```bash
git clone https://github.com/langrenne/MLFMM2HRTF.git
cd MLFMM2HRTF
make                             # main_HRTF_Ray and main_HRTF_Diff
cd test && ./run_test.sh         # 2 regression cases (Ray, Diff) on a small sphere, ~15 s
```

Build options:

```bash
make ARCH=-march=native          # processor older than ~2013 or ARM (default: -march=x86-64-v3, AVX2)
make LAPACK="-lopenblas"         # another LAPACK/BLAS (default: -llapack -lblas)
make LAPACK="-L/path/to/lib -llapack -lblas"
```

After changing an option, rebuild everything: `make clean && make ARCH=...`.

## Quick start

```bash
python python/source_oreille.py S01.msh --ear R               # source triangle (Ray) and point (Diff)
python python/generate_lebedev.py --n 2702 --dist 1.5 --out-dir .       # in the run directory
cp config_RAY.txt config_S01_R.txt   # edit: mesh, pext_lebedev_2702.txt, stem S01_R, source triangle
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=24
./main_HRTF_Ray config_S01_R.txt                              # then config_S01_L.txt
python python/hrtf2hrir.py config_S01_L.txt --weights weights_lebedev_2702.txt -o hrir_S01.npz
python python/npz2sofa.py hrir_S01.npz S01.sofa --subject S01
```

Or the whole chain for one subject (both ears, interpolation, HRIR, SOFA):
`./run_subject.sh config_S01_L.txt config_S01_R.txt` (`-h` for the options).

Computing time of one ear on 24 cores and 62 GB, fast setting: about 1 h 20
for a 48 000-element head and torso, 3 h 15 for 337 000 elements (335
frequencies of the mixed grid; see `MANUAL.md`).

## Repository structure

```
mlfmm/                  FMM/octree library (git submodule)
main_HRTF_Ray.f08       solver, triangle source
main_HRTF_Diff.f08      solver, monopole source
routines_Fait_tout.f08  BEM/MLFMM routines shared by both solvers
routines_gmres.f08      GMRES
config_RAY.txt          example case files
config_DIFF.txt
makefile
run_subject.sh          whole chain for one subject: both ears -> SOFA
python/                 tools: mesh checks, sources, grids, interpolation, HRIR, SOFA
test/                   regression tests (run_test.sh)
```

## Conventions

- Frame: x front, y left, z up. Right ear on the y < 0 side.
- Azimuth `atan2(y, x)`: 0° front, 90° left; elevation `asin(z/r)`. In SOFA,
  azimuth in [0°, 360°[.
- Time dependence e^{−iωt}.
- Normals pointing into the air (out of the head).

## Limitations

- Flat triangular elements only.
- The octree is limited to 10 levels: avoid elements smaller than about 0.5 mm
  on a head with torso.
- Tested on Linux only.

## Licence

Code: MIT License, Copyright (c) 2024-2026 Conservatoire national des arts et
métiers (Cnam), Christophe Langrenne; see [`LICENSE`](LICENSE).

Third-party code keeps its own licence (SPARSKIT: LGPL 2.1 or later;
modern_fftpack / FFTPACK 5.1: MIT and UCAR notice; Bessel routines of Zhang
and Jin; Lebedev quadrature after R. Parrish: BSD; Peterson sort: public
domain): see [`THIRD_PARTY.md`](THIRD_PARTY.md).

HRTF data (SOFA files) computed with this software and distributed by the
authors: CC BY 4.0.

If you use MLFMM2HRTF, please cite the software, C. Langrenne, *MLFMM2HRTF*,
Zenodo, [doi:10.5281/zenodo.24473755](https://doi.org/10.5281/zenodo.24473755)
(all versions; the DOI of each version is on that page, and GitHub's "Cite this
repository" gives the reference), the preprint (reference to come) and, for the
Lebedev grids, V. I. Lebedev and D. N. Laikov, Doklady Mathematics 59(3),
477–481, 1999.
