# Third-party code

MLFMM2HRTF is distributed under the licence given in [`LICENSE`](LICENSE). The
parts below come from other authors and keep their own terms.

| Code | Where | Licence |
|------|-------|---------|
| SPARSKIT (ILUT preconditioner, sparse format conversions), Y. Saad | `mlfmm/sparse_kit.f08` | GNU LGPL 2.1 or later |
| modern_fftpack (J. Lo Kim Lin), from FFTPACK 5.1 (UCAR) | `mlfmm/fftpack_*.f90` | MIT + UCAR BSD-type notice |
| Bessel and spherical Bessel routines, S. Zhang and J. Jin | `mlfmm/module_Acoustiques.f08` | free use with acknowledgement |
| Lebedev quadrature, port of the MATLAB code of R. Parrish | `mlfmm/module_Utils.f08`, `python/generate_lebedev.py` | BSD |
| Unique sort, W. M. Peterson | `mlfmm/module_sort.f08` | public domain |

## SPARSKIT

Copyright (C) 2005 the Regents of the University of Minnesota. SPARSKIT is free
software; you can redistribute it and/or modify it under the terms of the GNU
Lesser General Public License as published by the Free Software Foundation
(version 2.1 of the License, or any later version).
<https://www-users.cse.umn.edu/~saad/software/SPARSKIT/>.
The LGPL applies to `mlfmm/sparse_kit.f08` only (routines adapted to complex
arithmetic); it does not extend to the rest of the code.

## modern_fftpack / FFTPACK 5.1

MIT License, Copyright (c) 2016 Jon Lo Kim Lin; FFTPACK 5.1, Copyright (c) 2004
the University Corporation for Atmospheric Research (UCAR). Full texts:
[`mlfmm/LICENSE_fftpack.md`](mlfmm/LICENSE_fftpack.md).

## Zhang and Jin

The Bessel routines are copyrighted by Shanjie Zhang and Jianming Jin, who give
permission to incorporate them into a user program provided that the copyright
is acknowledged. Reference: S. Zhang, J. Jin, *Computation of Special
Functions*, Wiley, 1996, ISBN 0-471-11963-6.

## Lebedev quadrature

Grids of V. I. Lebedev and D. N. Laikov, "A quadrature formula for the sphere of
the 131st algebraic order of accuracy", Doklady Mathematics 59(3), 477–481,
1999 (citation requested by the authors). Code: port of `getLebedevSphere.m`
by Rob Parrish (MATLAB File Exchange 27097, 2010), itself translated from the
Fortran of C. van Wüllen and the C of D. Laikov.

Licence of `getLebedevSphere.m` (applies to the ported code in
`mlfmm/module_Utils.f08` and `python/generate_lebedev.py`):

```
Copyright (c) 2010, Robert Parrish
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## Peterson sort

`unique_sort`, after W. M. Peterson, released into the public domain.

## Quadrature rules

The triangle rules of Dunavant (D. A. Dunavant, "High degree efficient symmetrical
Gaussian quadrature rules for the triangle", Int. J. Numer. Methods Eng. 21,
1129–1148, 1985) are tables of numbers taken from the paper (checked against
J. Burkardt's `dunavant.f90`); no code is copied.
