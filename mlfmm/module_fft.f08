! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: MIT
! =============================================================================
! module_fft : FFT complexes multiples pour les translations HF (M2M, L2L) de
! module_MLFMM.
!
! Adaptateur vers modern_fftpack (fichiers fftpack_*.f90, licences dans
! LICENSE_fftpack.md) : memes noms et memes arguments que FFTPACK 5.1
! (cfftmi, cfftmf, cfftmb), mais le tableau c est passe en COMPLEXE comme dans
! notre code (modern_fftpack le declare real c(2,lenc), meme memoire).
!
! Conventions (identiques a FFTPACK 5.1) :
!   cfftmf : c_k <- (1/n) sum_j c_j exp(-2 i pi j k / n)   (normalisee)
!   cfftmb : c_k <-       sum_j c_j exp(+2 i pi j k / n)
! module_MLFMM appelle volontairement cfftmb la ou l'on attendrait une FFT
! directe et inversement (convention e^{-i omega t}) : ne pas permuter.
! Sans etat (pas de SAVE) : appelable depuis plusieurs threads.
! =============================================================================
module module_fft
	use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer
	use complex_transform_routines, only: lib_cfftmi => cfftmi, &
		lib_cfftmf => cfftmf, lib_cfftmb => cfftmb
	implicit none
	private
	public :: cfftmi, cfftmf, cfftmb

	integer, parameter :: dp = kind(1.0d0)

contains

	subroutine cfftmi(n, wsave, lensav, ier)
		integer, intent(in) :: n, lensav
		real(dp), intent(out) :: wsave(lensav)
		integer, intent(out) :: ier
		call lib_cfftmi(n, wsave, lensav, ier)
	end subroutine cfftmi

	subroutine cfftmf(lot, jump, n, inc, c, lenc, wsave, lensav, work, lenwrk, ier)
		integer, intent(in) :: lot, jump, n, inc, lenc, lensav, lenwrk
		complex(dp), intent(inout), target :: c(lenc)
		real(dp), intent(in) :: wsave(lensav)
		real(dp), intent(inout) :: work(lenwrk)
		integer, intent(out) :: ier
		real(dp), pointer :: cr(:,:)
		call c_f_pointer(c_loc(c), cr, [2, lenc])
		call lib_cfftmf(lot, jump, n, inc, cr, lenc, wsave, lensav, work, lenwrk, ier)
	end subroutine cfftmf

	subroutine cfftmb(lot, jump, n, inc, c, lenc, wsave, lensav, work, lenwrk, ier)
		integer, intent(in) :: lot, jump, n, inc, lenc, lensav, lenwrk
		complex(dp), intent(inout), target :: c(lenc)
		real(dp), intent(in) :: wsave(lensav)
		real(dp), intent(inout) :: work(lenwrk)
		integer, intent(out) :: ier
		real(dp), pointer :: cr(:,:)
		call c_f_pointer(c_loc(c), cr, [2, lenc])
		call lib_cfftmb(lot, jump, n, inc, cr, lenc, wsave, lensav, work, lenwrk, ier)
	end subroutine cfftmb

end module module_fft
