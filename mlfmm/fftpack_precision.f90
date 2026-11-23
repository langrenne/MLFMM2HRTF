! -----------------------------------------------------------------------------
! From modern_fftpack (https://github.com/jlokimlin/modern_fftpack), copied
! unchanged into mlfmm except for this header and the file name (prefix
! fftpack_). Only the complex multiple-vector transforms (cfftmi, cfftmf,
! cfftmb) are used, through module_fft.f08. Licences: mlfmm/LICENSE_fftpack.md
!
! modern_fftpack: MIT License, Copyright (c) 2016 Jon Lo Kim Lin.
! Modernization of FFTPACK 5.1 (P. N. Swarztrauber, R. A. Valent), whose
! notice follows.
!
! Copyright (c) 2004 the University Corporation for Atmospheric
! Research ("UCAR"). All rights reserved. Developed by NCAR's
! Computational and Information Systems Laboratory, UCAR,
! www.cisl.ucar.edu.
! Redistribution and use of the Software in source and binary forms,
! with or without modification, is permitted provided that the conditions
! in mlfmm/LICENSE_fftpack.md are met.
! -----------------------------------------------------------------------------
module fftpack_precision

    ! Explicit typing only
    implicit none

    ! Everything is private unless stated otherwise
    private
    public :: wp, ip
    public :: pimach, epmach

    !-----------------------------------------------
    ! Dictionary: precision constants
    !-----------------------------------------------
    integer,   parameter :: ip = selected_int_kind(r=9)
    integer,   parameter :: sp = selected_real_kind(p=6, r=37)
    integer,   parameter :: dp = selected_real_kind(p=15, r=307)
    integer,   parameter :: qp = selected_real_kind(p=33, r=4931)
    integer,   parameter :: wp = dp
    !-----------------------------------------------

contains


    pure function pimach() result (return_value)
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        real (wp) :: return_value
        !-----------------------------------------------

        return_value = 3.141592653589793238462643383279502884197169399375105820974_wp

    end function pimach


    pure function epmach() result (return_value)
        !
        ! Purpose:
        !
        ! Computes an approximate machine epsilon (accuracy), i.e.,
        !
        ! the smallest number x of the kind wp such that 1 + x > 1.
        !
        !-----------------------------------------------
        ! Dummy arguments
        !-----------------------------------------------
        real (wp) :: return_value
        !-----------------------------------------------

        return_value = epsilon(1.0_wp)

    end function epmach


end module fftpack_precision
