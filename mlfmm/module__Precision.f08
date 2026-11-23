! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: MIT
! module__Precision : sortes (kinds) des reels et des entiers du code.
! Reecrit le 2026-11-23 (memes noms et memes valeurs qu'avant).
module module__Precision
	implicit none
	private
	public :: prec, iprec, iprec4, iprec8
	integer, parameter :: prec   = kind(1.0d0)            ! reels et complexes : double precision
	integer, parameter :: iprec  = kind(1)                ! entiers courants (32 bits)
	integer, parameter :: iprec4 = selected_int_kind(9)   ! entiers 32 bits
	integer, parameter :: iprec8 = selected_int_kind(18)  ! entiers 64 bits
end module module__Precision
