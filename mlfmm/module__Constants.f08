! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: MIT
! module__Constants : constantes mathematiques du code.
! Reecrit le 2026-11-23 (memes noms et memes valeurs qu'avant).
module module__Constants
	use module__Precision
	implicit none
	private
	public :: pi, jcomplex
	real(prec),    parameter :: pi       = 3.141592653589793_prec
	complex(prec), parameter :: jcomplex = (0.0_prec, 1.0_prec)      ! i
end module module__Constants
