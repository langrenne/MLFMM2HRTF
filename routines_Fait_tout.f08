! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: MIT
!
!----------------------------------------------------------------------------
!
! Routine qui calcule la pression sur toute la surface à partir du second membre
! Structure rigide
! Burton/Miller, donc avec alpha et calcule des 4 opérateurs D,M,H,K
!----------------------------------------------------------------------------
    subroutine fait_tout(i1,NL,N,octree,SPH,SPH0,k,dp1,iter,sol2,iguess,x0)
!
!----------------------------------------------------------------------------
!
! Variables in :
! i1 : rang de fréquence
! NL : taille octree
! octree(NL)
! k : nbre d'onde
! dp1(NE) : gradient sur la structure
! SPH : MLFMM variables frequency dependant
! SPH0 : MLFMM variables non frequency dependant
!
! variables out :
! iter : nombre d'iterations pour le résultat
! sol2(NE) : pression sur la structure
! iguess, x0(NE) (2026-10-01) : iguess = 0, depart GMRES a zero ; iguess = 1,
!   depart x0 (solution de la frequence voisine, voir guess dans parametres)
!----------------------------------------------------------------------------
    use omp_lib
    use module__Precision
    use parametres
    use module_Octree
    use module_BEM
    use module_MLFMM
    use module_gmres
    use module_utils
    implicit none
    ! variables d'entrées
    integer :: i1,NL
    real(prec) :: k
    complex(prec) :: dp1(NE)
    type(tree) :: octree(NL)

    ! variables de sorties
    integer,intent(out) :: iter
    complex(prec),intent(out) :: sol2(NE)
    integer, intent(in) :: iguess
    complex(prec), intent(in) :: x0(NE)

    ! variables locales
!	type(FMM), dimension(:),allocatable :: SPH
	type(FMM) :: SPH(NL)
	type(FMM0) :: SPH0

    integer :: i2,N(NL)
	integer :: ierr,iwk,lfil,lwork
!	integer :: jA0(num_K0),iA0(NE+1),jU(NE),jwork(2*NE)
	integer :: jU(NE)
!	complex(prec) :: w(NE+1)
    complex(prec) :: Anear(num_K0),alpha,var1(NE)
	integer, dimension(:), allocatable :: jLU
	complex(prec), dimension(:), allocatable :: ALU

	integer :: irc(5), icntl(8), info(3)
	integer, allocatable :: jA0(:),iA0(:),jwork(:)
	real(prec) :: cntl(5), rinfo(2)

    complex(prec), allocatable :: Bnear(:),Anear0(:),temp3(:),w(:)
	complex(prec), allocatable :: work(:)
    real :: t01,t02
    real(prec) :: eps0_k            ! epsilon0 de la bande de frequence
	! CHIEF (essai) : systeme borde de taille n_gm = NE + nc (nc = 0 : inchange)
	integer :: nc, n_gm, nc_orth
	complex(prec), allocatable :: Cc(:,:), Mc(:,:), dc(:)

!
!---------------------------------------------------------------------------
!
! initialisation MLFMM
! methode donnée dans le config.txt (0=BF ; 1=HF ; 2=HYB)
! précision selon epsilon0 donnée dans config.txt
!---------------------------------------------------------------------------
!
!    write(*,*) 'in Faittout k=',k
!	write(*,*) 'in fait tout'
!    write(*,*) 'SPH0%tab_T_M2M=',SPH0%tab_T_M2M(:)

    t01=omp_get_wtime()
!	allocate(SPH(NL))
	eps0_k = eps0_de_k(k)          ! precision FMM de la bande de frequence (parametres)
	if (nb_bandes > 0 .and. verbose >= 2) write(log_unit,'(A,I5,A,ES9.2,A,ES9.2)') ' #', i1, ' band: epsilon0 =', eps0_k, &
	     '  epsilon2 =', eps2_de_k(k)
	call init_FMM_T0(octree,NL,k,N,SPH,SPH0,FMM_meth,eps0_k)
	t02=omp_get_wtime()
	if (verbose >= 2) write(log_unit,'(A,I5,A,*(I2))') ' #', i1, ' LF levels (1 = LF):', SPH(:)%BF
	if (verbose >= 2) write(log_unit,'(A,I5,A,*(I4))') ' #', i1, ' FMM orders N:', N
	if (verbose >= 2) write(log_unit,'(A,I4,A,*(I5))') 'CHECK[',i1,'] ordre_N    =', N
!	write(*,*) 'après init_FMM'
!    write(*,*) 'SPH0%tab_T_M2M=',SPH0%tab_T_M2M(:)

!
!---------------------------------------------------------------------------
!
! Premier calcul FMM pour second membre
!
!---------------------------------------------------------------------------
!
    alpha=ji*k/(k**2+C_BM)                         ! C_BM = pi par defaut (parametres)
    if (k < k_max_chief) alpha=(0.0_prec,0.0_prec)   ! CHIEF remplace Burton/Miller (cle interne f_max_chief)
	t01=omp_get_wtime()
!	call BM_1_MLFMM_FFT(dp1,octree,NL,N,SPH,alpha,k,var1)
!	call MLFMM_FFT_dnp(dp1,octree,NL,N,SPH,alpha,k,var1)
!	call MLFMM_FFT_p(dp1,octree,NL,N,SPH,alpha,k,var1)
!	write(*,*) 'FFT, var1=',var1(1:3)
	call MLFMM_HYB_dnp(dp1,octree,NL,N,SPH,SPH0,alpha,k,var1)
	t02=omp_get_wtime()
	if (verbose >= 2) write(log_unit,'(A,I5,A,F9.1,A)') ' #', i1, ' right-hand side (FMM)  ', t02-t01, ' s'
	if (verbose >= 2) write(log_unit,'(A,I4,A,ES15.8)') 'CHECK[',i1,'] norm_var1  =', sqrt(sum(abs(var1)**2))
!
!---------------------------------------------------------------------------
!
! calcul BEM de Anear et Bnear
!
!---------------------------------------------------------------------------
!

	t01=omp_get_wtime()
	allocate(Bnear(num_K0))
	call ABnear_FMM_T0(Anear,Bnear,octree,NL,k,alpha)
	t02=omp_get_wtime()
	if (verbose >= 2) write(log_unit,'(A,I5,A,F9.1,A)') ' #', i1, ' near field            ', t02-t01, ' s'
	if (verbose >= 2) then
		write(log_unit,'(A,I4,A,ES15.8)') 'CHECK[',i1,'] norm_Anear =', sqrt(sum(abs(Anear)**2))
		write(log_unit,'(A,I4,A,ES15.8)') 'CHECK[',i1,'] norm_Bnear =', sqrt(sum(abs(Bnear)**2))
	end if

	! Pour second membre et deallocation Bnear
	allocate(temp3(NE))
    temp3(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
	do i2=1,num_K0
	  temp3(iA(i2))=temp3(iA(i2))+Bnear(i2)*dp1(jA(i2))
	  if (isnan(abs(Anear(i2)))) write(log_unit,*) 'ERROR #', i1, ': near-field term Anear(', i2, ') is NaN'
	  if (isnan(abs(Bnear(i2)))) write(log_unit,*) 'ERROR #', i1, ': near-field term Bnear(', i2, ') is NaN'
	end do
	! Ray (Kreuzer) avec impedance (2026-11-06) : dnp = -i k A_red p + q0 ; second membre B.q0 ci-dessus,
	! operateur A + i k A_red B (champ proche ; champ lointain : fmm_T0_adapt, V2 <- Vc)
	if (beta_Vc /= (0.0_prec,0.0_prec)) Anear = Anear + ji*k*beta_Vc*Bnear
    deallocate(Bnear)
	if (verbose >= 2) write(log_unit,'(A,I4,A,ES15.8)') 'CHECK[',i1,'] norm_rhs   =', sqrt(sum(abs(temp3+var1)**2))

!
!---------------------------------------------------------------------------
!
!		decomposition ILUT
!
!---------------------------------------------------------------------------
!
	t01=omp_get_wtime()

!		passage coo en crs de Anear
	allocate(Anear0(num_K0))
	allocate(iA0(NE+1))
	allocate(jA0(num_K0))
	call zcoocsr(NE,num_K0,Anear,iA,jA,Anear0,jA0,iA0)
	! allocations
	lfil=lfil_ilut      ! config HRTF (ligne facultative, defaut 50)
	! taille des facteurs ILUT (2026-09-28) : au plus lfil valeurs par ligne dans L
	! et dans U, plus la diagonale -> borne EXACTE NE*(2*lfil+1), independante du
	! near-field num_K0. coef_ilut >= 1 ne peut donc jamais manquer de place.
	! (ancien : iwk = num_K0*coef_ilut, surdimensionne et dependant de NL)
	iwk=nint(coef_ilut*real(NE,prec)*real(2*lfil+1,prec))
	allocate(ALU(iwk))
	allocate(jLU(iwk))
	allocate(jwork(2*NE))
	allocate(w(NE+1))
! Ilut avec precision epsilon1 donnée dans le fichier config.txt
! si err=-4 augmenter l'allocation de iwk avec un facteur plus grand
	call ilut(NE,num_K0,Anear0,jA0,iA0,lfil,epsilon1,ALU,jLU,jU,iwk,w,jwork,ierr)
	if (ierr.ne.0) then
		  write(log_unit,'(A,I5,A,I0)') ' ERROR #', i1, ': ILUT preconditioner failed, ierr = ', ierr
	end if
	t02=omp_get_wtime()
	if (verbose >= 2) write(log_unit,'(A,I5,A,F9.1,A)') ' #', i1, ' ILUT preconditioner   ', t02-t01, ' s'

	deallocate(Anear0)
	deallocate(jA0)
	deallocate(iA0)
	deallocate(jwork)
	deallocate(w)
!---------------------------------------------------------------------------
!
! GMRES résolution
!
!---------------------------------------------------------------------------
!
    ! Init GMRES
	call gmres_defauts(icntl,cntl)
	! warning display
	icntl(2)=0
	! left preconditioning donné dans le fichier config.txt
	icntl(4) = FMM_pre  !
	! ICGS orthogonalization
	icntl(5) = 0
	! init guess : 0, ou x0 fourni (frequence voisine)
	icntl(6) = merge(1, 0, iguess > 0)
	icntl(8) = 0
	! Maximum number of iterations
	icntl(7) = iter_max            ! parametres (HRTF : cle iter_max, defaut 150)
	if (verbose >= 4) icntl(3) = log_unit   ! historique de convergence GMRES dans le journal (diagnostic, 1 frequence)
	! tolerance donnée dans config.txt
    cntl(1) = eps2_de_k(k)          ! tolerance de la bande de frequence (parametres)

	! CHIEF (essai) : la ou Burton/Miller est coupe, nc lignes interieures en plus
	! (champ rayonne nul a l'interieur : D p - S q = 0, soit C p = S dp1)
	nc = 0
	if (NC_chief > 0 .and. k < k_max_chief .and. k >= k_min_chief) nc = NC_lignes  ! f_min_chief : pas de CHIEF en BF
	n_gm = NE + nc
	if (nc > 0) then
		allocate(Cc(nc,NE), Mc(nc,NE), dc(nc))
		call lignes_chief_comb(k, nc, Cc, Mc)
		dc = matmul(Mc, dp1)
		deallocate(Mc)
		call normalise_chief(nc, Cc, dc)
		if (chief_orth > 0.0_prec) then
			call orthonormalise_chief(nc, Cc, dc, nc_orth)
			Cc = Cc(1:nc_orth,:)     ! reallocation a la taille gardee
			dc = dc(1:nc_orth)
			nc = nc_orth
			n_gm = NE + nc
		endif
	endif
    !write(*,*) 'epsilon2=',epsilon2
    !write(*,*) 'NE=',NE
    !write(*,*) 'restart=',restart
    !write(*,*) 'NP=',NP

	! allocation GMRES
	if (icntl(5)==0 .OR. icntl(5)==1) then
        if (icntl(8)==0) then
            lwork=restart**2+restart*(n_gm+5)+6*n_gm+2  ! ICNTL(5)=0 ou 1 et ICNTL(8)=0
        else
            lwork=restart**2+restart*(n_gm+5)+5*n_gm+2  ! ICNTL(5)=0 ou 1 et ICNTL(8)=1
        end if
    else
        if (icntl(8)==0) then
            lwork=restart**2+restart*(n_gm+5)+6*n_gm+restart+1  ! ICNTL(5)=2 ou 3 et ICNTL(8)=0
        else
            lwork=restart**2+restart*(n_gm+5)+5*n_gm+restart+1  ! ICNTL(5)=2 ou 3 et ICNTL(8)=1
        end if
    end if

    !write(*,*) 'lwork=',lwork
	allocate(work(lwork))

	!Second membre
	work(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
	work(n_gm+1:n_gm+NE)=temp3+var1
	if (nc > 0) work(n_gm+NE+1:2*n_gm)=dc      ! CHIEF : lignes interieures
	if (iguess > 0) work(1:NE) = x0     ! estimation initiale (x dans work(1:n), voir module_gmres)
	deallocate(temp3)

	t01=omp_get_wtime()

	if (nc > 0) then
		call sub_gmres_chief(icntl,cntl,info,rinfo,&
		Anear,ALU,jLU,jU,iwk,octree,NL,N,SPH,SPH0,alpha,k,(0.0_prec,0.0_prec),.false.,nc,Cc,lwork,work)
	else
		call sub_gmres_HYB(icntl,cntl,info,rinfo,&
		Anear,ALU,jLU,jU,iwk,octree,NL,N,SPH,SPH0,alpha,k,lwork,work)
	endif
!	call sub_gmres_FFT(icntl,cntl,info,rinfo,&
!	Anear,ALU,jLU,jU,iwk,octree,NL,N,SPH,alpha,k,lwork,work)

	if (info(1).ne.0) then
		  if (info(1) == -4) then
		    write(log_unit,'(A,I5,A,I0,A,ES9.2,A,ES9.2,A)') ' WARNING #', i1, ': GMRES not converged in ', info(2), &
		      ' iterations (backward error ', rinfo(1), ', tolerance ', cntl(1), ')'
		  else
		    write(log_unit,'(A,I5,A,I0)') ' ERROR #', i1, ': GMRES failed, info = ', info(1)
		  endif
	endif
	t02=omp_get_wtime()
	if (verbose >= 2) write(log_unit,'(A,I5,A,I4,A,F9.1,A)') ' #', i1, ' GMRES', info(2), ' iterations ', t02-t01, ' s'
	if (verbose >= 3) write(log_unit,'(A,I5,A,ES9.2)') ' #', i1, ' GMRES true backward error', rinfo(1)
	iter=info(2)
	sol2=work(1:NE)
	if (nc > 0) then      ! CHIEF : lambda ~ 0 si les lignes sont compatibles (grand : point sur un noeud ?)
		write(log_unit,'(A,I5,A,I4,A,ES10.3)') ' #', i1, ' CHIEF:', nc, ' points, max|lambda| =', maxval(abs(work(NE+1:n_gm)))
		deallocate(Cc, dc)
	endif
	if (verbose >= 2) then
		write(log_unit,'(A,I4,A,I5)')   'CHECK[',i1,'] gmres_iter =', iter
		write(log_unit,'(A,I4,A,ES15.8)') 'CHECK[',i1,'] norm_sol2  =', sqrt(sum(abs(sol2)**2))
	end if

	end subroutine fait_tout
!
!----------------------------------------------------------------------------
!
! Routine qui calcule la pression sur toute la surface à partir du second membre
! Structure rigide
! Burton/Miller, donc avec alpha et calcule des 4 opérateurs D,M,H,K
!----------------------------------------------------------------------------
    subroutine fait_tout_HRTF_DIFF_T0(i1,NL,octree,SPH0,k,p1,dp1,pref,pref0,A_red,iter,sol,iguess,x0,surf)
!
!----------------------------------------------------------------------------
!
! Variables in :
! i1 : rang de fréquence
! NL : taille octree
! octree(NL)
! k : nbre d'onde
! p1(NE) : pression sur la structure
! dp1(NE) : gradient sur la structure
! pref(Next) : pour normalisation par pression au centre
! A_red(complx) : Admittance complexe réduite (A_red=1/Z_red), 0 pour paroi rigide
!
! variables out :
! iter : nombre d'iterations pour le résultat
! sol(Next) : pression aux points Pext
!----------------------------------------------------------------------------
    use omp_lib
    use module__Precision
    use parametres
    use module_Octree
    use module_BEM
    use module_MLFMM
    use module_gmres
    use module_utils
    implicit none
!    INCLUDE 'mpif.h'
!    INCLUDE 'dmumps_struc.h'
!    INCLUDE 'zmumps_struc.h'
    ! variables d'entrées
    integer :: i1,NL
    real(prec) :: k
    complex(prec) :: dp1(NE),p1(NE),pref(Next),pref0(Next),A_red
    type(tree) :: octree(NL)

    ! variables de sorties
    integer,intent(out) :: iter
    complex(prec),intent(inout) :: sol(Next)
    ! estimation initiale GMRES (iguess, x0) et pression sur la surface (surf), comme
    ! fait_tout_HRTF_T0 (guess, Diff : 2026-10-02). iguess = -1 : x0 et surf inutilises
    integer, intent(in) :: iguess
    complex(prec), intent(in) :: x0(NE)
    complex(prec), intent(out) :: surf(NE)

    ! variables locales
!	type(FMM), dimension(:),allocatable :: SPH
	type(FMM) :: SPH(NL)
	type(FMM0) :: SPH0

    integer :: i2,N(NL),LEV
	integer :: ierr,iwk,lfil,lwork
!	integer :: jA0(num_K0),iA0(NE+1),jU(NE),jwork(2*NE)
	integer :: jU(NE)
!	complex(prec) :: w(NE+1)
    complex(prec) :: Apnear(num_K0),alpha,pout(Next),sol2(NE)
	integer, dimension(:), allocatable :: jLU
	complex(prec), dimension(:), allocatable :: ALU

	integer :: irc(5), icntl(8), info(3)
	integer, allocatable :: jA0(:),iA0(:),jwork(:)
	real(prec) :: cntl(5), rinfo(2)

    complex(prec), allocatable :: Anear(:),Bnear(:),Anear0(:),temp3(:),w(:)
	complex(prec), allocatable :: work(:)
    real :: t01,t02
    real(prec) :: eps0_k            ! epsilon0 de la bande de frequence
	! CHIEF (essai) : systeme borde de taille n_gm = NE + nc (nc = 0 : inchange)
	integer :: nc, n_gm, nc_orth
	complex(prec), allocatable :: Cc(:,:), Mc(:,:), dc(:)

	! test bem/mlfmm
	integer :: i3,no
	complex(prec) :: dnsol2(NE),p_test(5),D1,M1,H1,K1,p_test2(Next),var1(NE),var2(NE),temp(NE),pmlfmm2(NE)
	integer, dimension(1:3) :: nPe
	real(prec), dimension(1:3,1:3) :: PQ
	real(prec), dimension(1:3,1:IS_v) :: PG_v
	real(prec) :: PG2(2,IS_v)
	!real(prec), dimension(1:3,1:IS_v) :: PG_v
	!real(prec) :: PG_c2(2,IS_c),PG_v2(2,IS_v),PGc3(IS_c),PGv3(IS_v),Pcin(3),norcin(3)
	real(prec), dimension(1:2,1:IS_s) :: PG_s
	real(prec), dimension(1:IS_s) :: x, wth
	!real(prec), dimension(1:3,1:IS_c) :: F_c, dF1_c, dF2_c
	real(prec), dimension(1:3,1:IS_v) :: F_v, dF1_v, dF2_v
	real(prec), dimension(1:3,1:IS_f) :: F_f, dF1_f, dF2_f
	real(prec) :: Qc(1:3,IS_c), Qv(1:3,IS_v), Qf(1:3,IS_f)

!
!---------------------------------------------------------------------------
!
! initialisation MLFMM
! methode donnée dans le config.txt (0=BF ; 1=HF ; 2=HYB)
! précision selon epsilon0 donnée dans config.txt
!---------------------------------------------------------------------------
!
!    write(*,*) 'in Faittout k=',k
    t01=omp_get_wtime()
!	allocate(SPH(NL))
	eps0_k = eps0_de_k(k)          ! precision FMM de la bande de frequence (parametres)
	if (nb_bandes > 0 .and. verbose >= 2) write(log_unit,'(A,I5,A,ES9.2,A,ES9.2)') ' #', i1, ' band: epsilon0 =', eps0_k, &
	     '  epsilon2 =', eps2_de_k(k)
	beta_Vc = A_red                ! Diff : Vc = V2 fusionne avec V1 (fusion_Vc, module_MLFMM, 2026-11-06)
	call init_FMM_T0(octree,NL,k,N,SPH,SPH0,FMM_meth,eps0_k)
	t02=omp_get_wtime()
	if (verbose >= 2) write(log_unit,'(A,I5,A,*(I4))') ' #', i1, ' FMM orders N:', N
	if (verbose >= 2) write(log_unit,'(A,I4,A,*(I5))') 'CHECK[',i1,'] ordre_N    =', N
!
!---------------------------------------------------------------------------
!
! Premier calcul FMM pour initialisation
!
!---------------------------------------------------------------------------
!
    alpha=ji*k/(k**2+C_BM)                         ! C_BM = pi par defaut (parametres)
    if (k < k_max_chief) alpha=(0.0_prec,0.0_prec)   ! CHIEF remplace Burton/Miller (cle interne f_max_chief)
!	t01=omp_get_wtime()
!	call BM_1_MLFMM_FFT(dp1,octree,NL,N,SPH,alpha,k,var1)
!	call MLFMM_FFT_dnp(dp1,octree,NL,N,SPH,alpha,k,var1)
!	t02=omp_get_wtime()
!	write(*,*) 'passe BM_1 :',i1, ' in ',t02-t01,' sec'
!
!---------------------------------------------------------------------------
!
! calcul BEM de Anear et Bnear
!
!---------------------------------------------------------------------------
!

	t01=omp_get_wtime()
	allocate(Anear(num_K0))
	allocate(Bnear(num_K0))
	call ABnear_FMM_T0(Anear,Bnear,octree,NL,k,alpha)
	t02=omp_get_wtime()
	if (verbose >= 2) write(log_unit,'(A,I5,A,F9.1,A)') ' #', i1, ' near field            ', t02-t01, ' s'
	if (verbose >= 2) then
		write(log_unit,'(A,I4,A,ES15.8)') 'CHECK[',i1,'] norm_Anear =', sqrt(sum(abs(Anear)**2))
		write(log_unit,'(A,I4,A,ES15.8)') 'CHECK[',i1,'] norm_Bnear =', sqrt(sum(abs(Bnear)**2))
	end if

!
!---------------------------------------------------------------------------
!
! correction de Anear pour Admittance de paroi en mesh(3)
!
!---------------------------------------------------------------------------
!
	! IMPEDANCE : Anear.p - Bnear.dnp avec dnp = -i k A_red p
	! CONVENTION (2026-10-06) : normale vers le FLUIDE, e^{-iwt}, Z = p/v_paroi (vitesse entrant dans la paroi)
	! => v.n = -p/Z et dnp = -i k A_red p (NumCalc/Mesh2HRTF identique) ; avant : +i k A_red p (paroi active)
	Apnear=Anear + ji * k * A_red * Bnear
    deallocate(Anear)
    deallocate(Bnear)
	if (verbose >= 2) write(log_unit,'(A,I4,A,ES15.8)') 'CHECK[',i1,'] norm_Apnear=', sqrt(sum(abs(Apnear)**2))

!
!---------------------------------------------------------------------------
!
!		decomposition ILUT
!
!---------------------------------------------------------------------------
!
	t01=omp_get_wtime()

!		passage coo en crs de Anear
	allocate(Anear0(num_K0))
	allocate(iA0(NE+1))
	allocate(jA0(num_K0))
	call zcoocsr(NE,num_K0,Apnear,iA,jA,Anear0,jA0,iA0)
	! allocations
	lfil=lfil_ilut      ! config HRTF (ligne facultative, defaut 50)
	! taille des facteurs ILUT (2026-09-28) : au plus lfil valeurs par ligne dans L
	! et dans U, plus la diagonale -> borne EXACTE NE*(2*lfil+1), independante du
	! near-field num_K0. coef_ilut >= 1 ne peut donc jamais manquer de place.
	! (ancien : iwk = num_K0*coef_ilut, surdimensionne et dependant de NL)
	iwk=nint(coef_ilut*real(NE,prec)*real(2*lfil+1,prec))
	allocate(ALU(iwk))
	allocate(jLU(iwk))
	allocate(jwork(2*NE))
	allocate(w(NE+1))
! Ilut avec precision epsilon1 donnée dans le fichier config.txt
! si err=-4 augmenter l'allocation de iwk avec un facteur plus grand
	call ilut(NE,num_K0,Anear0,jA0,iA0,lfil,epsilon1,ALU,jLU,jU,iwk,w,jwork,ierr)
	if (ierr.ne.0) then
		  write(log_unit,'(A,I5,A,I0)') ' ERROR #', i1, ': ILUT preconditioner failed, ierr = ', ierr
	end if
	t02=omp_get_wtime()
	if (verbose >= 2) write(log_unit,'(A,I5,A,F9.1,A)') ' #', i1, ' ILUT preconditioner   ', t02-t01, ' s'

	deallocate(Anear0)
	deallocate(jA0)
	deallocate(iA0)
	deallocate(jwork)
	deallocate(w)
!---------------------------------------------------------------------------
!
! GMRES résolution
!
!---------------------------------------------------------------------------
!
    ! Init GMRES
	call gmres_defauts(icntl,cntl)
	! warning display
	icntl(2)=0
	! left preconditioning donné dans le fichier config.txt
	icntl(4) = FMM_pre  !
	! ICGS orthogonalization
	icntl(5) = 0
	! init guess
	icntl(6) = merge(1, 0, iguess > 0)   ! init guess : 0, ou x0 fourni (frequence voisine)
	icntl(8) = 0
	! Maximum number of iterations
	icntl(7) = iter_max            ! parametres (HRTF : cle iter_max, defaut 150)
	if (verbose >= 4) icntl(3) = log_unit   ! historique de convergence GMRES dans le journal (diagnostic, 1 frequence)
	! tolerance donnée dans config.txt
    cntl(1) = eps2_de_k(k)          ! tolerance de la bande de frequence (parametres)

	! CHIEF (essai) : la ou Burton/Miller est coupe, nc lignes interieures en plus
	! (champ total nul a l'interieur : p_inc + D p - S q = 0 avec q = -i k A_red p,
	! soit (D + i k A_red S) p = -p_inc, meme combinaison que Apnear)
	nc = 0
	if (NC_chief > 0 .and. k < k_max_chief .and. k >= k_min_chief) nc = NC_lignes  ! f_min_chief : pas de CHIEF en BF
	n_gm = NE + nc
	if (nc > 0) then
		allocate(Cc(nc,NE), Mc(nc,NE), dc(nc))
		call lignes_chief_comb(k, nc, Cc, Mc)
		Cc = Cc + ji*k*A_red*Mc
		deallocate(Mc)
		call combine_chief_vec(pinc_chief(:,i1), dc)
		dc = -dc
		call normalise_chief(nc, Cc, dc)
		if (chief_orth > 0.0_prec) then
			call orthonormalise_chief(nc, Cc, dc, nc_orth)
			Cc = Cc(1:nc_orth,:)     ! reallocation a la taille gardee
			dc = dc(1:nc_orth)
			nc = nc_orth
			n_gm = NE + nc
		endif
	endif

	! allocation GMRES
	if (icntl(5)==0 .OR. icntl(5)==1) then
        if (icntl(8)==0) then
            lwork=restart**2+restart*(n_gm+5)+6*n_gm+2  ! ICNTL(5)=0 ou 1 et ICNTL(8)=0
        else
            lwork=restart**2+restart*(n_gm+5)+5*n_gm+2  ! ICNTL(5)=0 ou 1 et ICNTL(8)=1
        end if
    else
        if (icntl(8)==0) then
            lwork=restart**2+restart*(n_gm+5)+6*n_gm+restart+1  ! ICNTL(5)=2 ou 3 et ICNTL(8)=0
        else
            lwork=restart**2+restart*(n_gm+5)+5*n_gm+restart+1  ! ICNTL(5)=2 ou 3 et ICNTL(8)=1
        end if
    end if

	allocate(work(lwork))

	!Second membre
	work(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
	work(n_gm+1:n_gm+NE)=-p1-alpha*dp1
	if (nc > 0) work(n_gm+NE+1:2*n_gm)=dc      ! CHIEF : lignes interieures
	if (iguess > 0) work(1:NE) = x0     ! estimation initiale (x dans work(1:n), voir module_gmres)
	if (verbose >= 2) write(log_unit,'(A,I4,A,ES15.8)') 'CHECK[',i1,'] norm_rhs   =', sqrt(sum(abs(work(n_gm+1:n_gm+NE))**2))

	t01=omp_get_wtime()

	if (nc > 0) then
		call sub_gmres_chief(icntl,cntl,info,rinfo,&
		Apnear,ALU,jLU,jU,iwk,octree,NL,N,SPH,SPH0,alpha,k,A_red,.true.,nc,Cc,lwork,work)
	else
		call sub_gmres_Z_HYB(icntl,cntl,info,rinfo,&
		Apnear,ALU,jLU,jU,iwk,octree,NL,N,SPH,SPH0,alpha,k,A_red,lwork,work)
	endif

	if (info(1).ne.0) then
		  if (info(1) == -4) then
		    write(log_unit,'(A,I5,A,I0,A,ES9.2,A,ES9.2,A)') ' WARNING #', i1, ': GMRES not converged in ', info(2), &
		      ' iterations (backward error ', rinfo(1), ', tolerance ', cntl(1), ')'
		  else
		    write(log_unit,'(A,I5,A,I0)') ' ERROR #', i1, ': GMRES failed, info = ', info(1)
		  endif
	endif
	t02=omp_get_wtime()
	if (verbose >= 2) write(log_unit,'(A,I5,A,I4,A,F9.1,A)') ' #', i1, ' GMRES', info(2), ' iterations ', t02-t01, ' s'
	if (verbose >= 3) write(log_unit,'(A,I5,A,ES9.2)') ' #', i1, ' GMRES true backward error', rinfo(1)
	iter=info(2)
	sol2=work(1:NE)
	if (nc > 0) then      ! CHIEF : lambda ~ 0 si les lignes sont compatibles (grand : point sur un noeud ?)
		write(log_unit,'(A,I5,A,I4,A,ES10.3)') ' #', i1, ' CHIEF:', nc, ' points, max|lambda| =', maxval(abs(work(NE+1:n_gm)))
		deallocate(Cc, dc)
	endif
	if (iguess >= 0) surf = sol2      ! pression de surface pour la frequence suivante du bloc
	if (verbose >= 2) then
		write(log_unit,'(A,I4,A,I5)')   'CHECK[',i1,'] gmres_iter =', iter
		write(log_unit,'(A,I4,A,ES15.8)') 'CHECK[',i1,'] norm_sol2  =', sqrt(sum(abs(sol2)**2))
	end if

    ! calcul les HRTF correspondant à Pext
	LEV=2
   	t01=omp_get_wtime()
	if (methode_propagation == 1) then      ! parametres (fichier technique MLFMM)
		call propagation_hierarchique(1,2.0_prec,1,k,NL,N,octree,SPH,SPH0,sol2,-ji*k*A_red*sol2,pout)
	else
		call propagation_FMM_only(LEV,k,NL,N,octree,SPH,SPH0,sol2,-ji*k*A_red*sol2,pout)
	endif
	t02=omp_get_wtime()
	if (verbose >= 2) write(log_unit,'(A,I5,A,F9.1,A)') ' #', i1, ' external points       ', t02-t01, ' s'
	if (verbose >= 2) write(log_unit,'(A,I4,A,ES15.8)') 'CHECK[',i1,'] norm_pout  =', sqrt(sum(abs(pout)**2))
	sol=(pref + pout)/pref0
	if (verbose >= 2) write(log_unit,'(A,I4,A,ES15.8)') 'CHECK[',i1,'] norm_HRTF  =', sqrt(sum(abs(sol)**2))


	end subroutine fait_tout_HRTF_DIFF_T0

!
!----------------------------------------------------------------------------
!
! A partir de fait_tout, la pression sur la surface et le gradient, calcule la propagation au point désiré
! pour le calcul d'HRTF, selon Gumerov par réciprocité.
! Maillage en T3, grandeurs physiques (p et dnp) constantes sur l'élément => T0
!
!----------------------------------------------------------------------------
    subroutine fait_tout_HRTF_T0(i1,NL,octree,SPH0,k,dp1,pref0,pref,iter,sol,iguess,x0,surf)
!
!----------------------------------------------------------------------------
!
! A partir de fait_tout, la pression sur la surface et le gradient, calcule la propagation au point désiré
! pour le calcul d'HRTF, selon Gumerov par réciprocité.
! Maillage en T3, grandeurs physiques (p et dnp) constantes sur l'élément => T0
!
!----------------------------------------------------------------------------
    use omp_lib
    use module__Precision
    use parametres
    use module_Octree
    use module_BEM
    use module_MLFMM
    use module_gmres
    use module_utils
    implicit none
    ! variables d'entrées
    integer :: i1,NL
    real(prec) :: k
    complex(prec) :: dp1(NE),pref0(Next),pref(Next)
    type(tree) :: octree(NL)
	type(FMM0) :: SPH0

    ! variables de sortie
    integer, intent(out) :: iter
    complex(prec),intent(out) :: sol(Next)
    ! estimation initiale GMRES (iguess, x0) et pression sur la surface (surf), pour
    ! enchainer les frequences voisines (guess, 2026-10-01). iguess = -1 : x0 et surf
    ! ne sont pas utilises (appel sans guess, tableaux factices de taille 1)
    integer, intent(in) :: iguess
    complex(prec), intent(in) :: x0(NE)
    complex(prec), intent(out) :: surf(NE)

    real :: t01,t02
    complex(prec)  :: sol2(NE),pout(Next),dnp_tot(NE)
    ! variables locales
	type(FMM) :: SPH(NL)
    integer :: i2,N(NL),LEV

!    write(*,*) 'in Faittout_hrtf k =',k
    ! Calcul la pression sur la surface à partir du gradient
!    write(*,*) 'avant fait_tout'
!    write(*,*) 'SPH0%tab_T_M2M=',SPH0%tab_T_M2M(:)

    call fait_tout(i1,NL,N,octree,SPH,SPH0,k,dp1,iter,sol2,iguess,x0)
    if (iguess >= 0) surf = sol2      ! iguess = -1 : ni estimation ni sortie surface (x0, surf factices)
!	write(*,*) 'sol2=',sol2(1:3)
    ! calcul les HRTF correspondant à Pext
    LEV=3
   	t01=omp_get_wtime()
	dnp_tot = dp1 - ji*k*beta_Vc*sol2       ! impedance (Kreuzer) : dnp = -i k A_red p + q0 ; rigide : dp1
	if (methode_propagation == 1) then      ! parametres (fichier technique MLFMM)
		call propagation_hierarchique(1,2.0_prec,1,k,NL,N,octree,SPH,SPH0,sol2,dnp_tot,pout)
	else
		call propagation_FMM_only(LEV,k,NL,N,octree,SPH,SPH0,sol2,dnp_tot,pout)
	endif
	t02=omp_get_wtime()
	if (verbose >= 2) write(log_unit,'(A,I5,A,F9.1,A)') ' #', i1, ' external points       ', t02-t01, ' s'
	if (verbose >= 2) write(log_unit,'(A,I4,A,ES15.8)') 'CHECK[',i1,'] norm_pout  =', sqrt(sum(abs(pout)**2))
	do i2=1,Next
		sol(i2)=(2.0_prec*pref(i2) + pout(i2))/pref0(i2)
	end do
	if (verbose >= 2) write(log_unit,'(A,I4,A,ES15.8)') 'CHECK[',i1,'] norm_HRTF  =', sqrt(sum(abs(sol)**2))

	end subroutine fait_tout_HRTF_T0

!
!----------------------------------------------------------------------------
!
! CHIEF (essai 2026-10-09) : lignes interieures en dense. Pour chaque point
! Q_chief(kc) : Dc(kc,j) = D1 et Mc(kc,j) = M1 de l'element j, memes noyaux et
! meme quadrature (IS_v) que la propagation directe (propagation_hierarchique) :
! champ rayonne en Q = sum_j D1 p_j - M1 q_j.
!
!----------------------------------------------------------------------------
	subroutine lignes_chief(k, nc, Dc, Mc)
	use module__Precision
	use parametres
	use module_BEM
	implicit none
	real(prec), intent(in) :: k
	integer, intent(in) :: nc
	complex(prec), intent(out) :: Dc(nc,NE), Mc(nc,NE)
	call lignes_chief_pts(k, nc, Q_chief(1:nc,:), Dc, Mc)
	end subroutine lignes_chief

! lignes D1, M1 pour une liste de points Qp(npt,3) quelconque
	subroutine lignes_chief_pts(k, npt, Qp, Dc, Mc)
	use module__Precision
	use parametres
	use module_BEM
	implicit none
	real(prec), intent(in) :: k, Qp(npt,3)
	integer, intent(in) :: npt
	complex(prec), intent(out) :: Dc(npt,NE), Mc(npt,NE)
	real(prec) :: PG_v(3,IS_v), F_v(3,IS_v), dF1_v(3,IS_v), dF2_v(3,IS_v), PQ(3,3), Qv(3,IS_v), norP(3)
	complex(prec) :: D1, K1, M1, H1
	integer :: kc, no

	call pgt6(PG_v,IS_v)
	call FORM_T3(PG_v(1:2,:),IS_v,F_v,dF1_v,dF2_v)
	norP=(/0.0_prec,0.0_prec,1.0_prec/)          ! sans effet sur D1 et M1
	do no=1,NE
		PQ=transpose(Ps(E(no,:),:))
		Qv=matmul(PQ,F_v)
		do kc=1,npt
			call DMKH_T0(k,Qv,norc(no,:),Jac0(no),IS_v,PG_v(3,:),Qp(kc,:),norP,D1,K1,M1,H1)
			Dc(kc,no)=D1
			Mc(kc,no)=M1
		end do
	end do
	end subroutine lignes_chief_pts

! CHIEF-B (essai 2026-10-13) : lignes aux NC_chief points puis combinaison G
! (NC_lignes contraintes, fichier chief_comb). Sans combinaison : lignes_chief.
	subroutine lignes_chief_comb(k, nc, Dc, Mc)
	use module__Precision
	use parametres
	use module_BEM
	implicit none
	real(prec), intent(in) :: k
	integer, intent(in) :: nc
	complex(prec), intent(out) :: Dc(nc,NE), Mc(nc,NE)
	complex(prec), allocatable :: Dp(:,:), Mp(:,:)
	integer :: jz, j0, j1
	integer, parameter :: npq = 32
	if (NZ_comb == 0) then
		call lignes_chief(k, nc, Dc, Mc)
		return
	endif
	! par paquets de npq points : memoire npq x NE au lieu de NC_chief x NE
	allocate(Dp(npq,NE), Mp(npq,NE))
	Dc = (0.0_prec,0.0_prec)
	Mc = (0.0_prec,0.0_prec)
	do j0 = 1, NC_chief, npq
		j1 = min(j0+npq-1, NC_chief)
		call lignes_chief_pts(k, j1-j0+1, Q_chief(j0:j1,:), Dp(1:j1-j0+1,:), Mp(1:j1-j0+1,:))
		do jz = 1, NZ_comb
			if (jG_comb(jz) < j0 .or. jG_comb(jz) > j1) cycle
			Dc(iG_comb(jz),:) = Dc(iG_comb(jz),:) + vG_comb(jz)*Dp(jG_comb(jz)-j0+1,:)
			Mc(iG_comb(jz),:) = Mc(iG_comb(jz),:) + vG_comb(jz)*Mp(jG_comb(jz)-j0+1,:)
		end do
	end do
	deallocate(Dp, Mp)
	end subroutine lignes_chief_comb

! CHIEF-B : meme combinaison G pour une grandeur donnee aux points (Diff : champ incident)
	subroutine combine_chief_vec(vp, vl)
	use module__Precision
	use parametres
	implicit none
	complex(prec), intent(in) :: vp(NC_chief)
	complex(prec), intent(out) :: vl(NC_lignes)
	integer :: jz
	if (NZ_comb == 0) then
		vl = vp
		return
	endif
	vl = (0.0_prec,0.0_prec)
	do jz = 1, NZ_comb
		vl(iG_comb(jz)) = vl(iG_comb(jz)) + vG_comb(jz)*vp(jG_comb(jz))
	end do
	end subroutine combine_chief_vec

! CHIEF (essai 2026-10-13) : contraintes orthonormalisees. C C^H = U S U^H (nc x nc) ;
! nouvelles lignes S^-1/2 U^H C (et second membre), en ne gardant que S > chief_orth x max.
! Les contraintes presque dependantes (16 par bloc CHIEF-B) donnaient lambda ~ 1e3.
	subroutine orthonormalise_chief(nc, Cc, dc, ne_)
	use module__Precision
	use parametres
	use module_BEM, only: NE
	implicit none
	integer, intent(in) :: nc
	integer, intent(out) :: ne_
	complex(prec), intent(inout) :: Cc(nc,NE), dc(nc)      ! lignes gardees en tete (1:ne_)
	complex(prec), allocatable :: G(:,:), Cn(:,:), dn(:), wk(:)
	real(prec), allocatable :: s(:), rwk(:)
	integer :: info_l, i, i0
	allocate(G(nc,nc), s(nc), wk(4*nc), rwk(3*nc))
	G = matmul(Cc, transpose(conjg(Cc)))
	call zheev('V', 'U', nc, G, nc, s, wk, 4*nc, rwk, info_l)      ! valeurs propres croissantes
	if (info_l /= 0) write(log_unit,*) 'CHIEF: zheev info =', info_l
	i0 = 1
	do while (s(i0) < chief_orth*s(nc))
		i0 = i0 + 1
	end do
	ne_ = nc - i0 + 1
	allocate(Cn(ne_,NE), dn(ne_))
	do i = i0, nc
		Cn(i-i0+1,:) = matmul(conjg(G(:,i)), Cc) / sqrt(s(i))
		dn(i-i0+1) = dot_product(G(:,i), dc) / sqrt(s(i))
	end do
	if (verbose >= 2) write(log_unit,'(A,I4,A,I4,A,ES9.2)') ' CHIEF:', ne_, ' constraints kept out of', nc, &
		', condition number ', sqrt(s(nc)/s(i0))
	Cc(1:ne_,:) = Cn
	dc(1:ne_) = dn
	deallocate(G, s, wk, rwk, Cn, dn)
	end subroutine orthonormalise_chief

! CHIEF : chaque ligne interieure ramenee a une norme 2 unite (second membre avec),
! pour qu'elle pese comme une ligne de surface (diagonale -1/2) dans la GMRES
	subroutine normalise_chief(nc, Cc, dc)
	use module__Precision
	use module_BEM, only: NE
	implicit none
	integer, intent(in) :: nc
	complex(prec), intent(inout) :: Cc(nc,NE), dc(nc)
	real(prec) :: s
	integer :: kc
	do kc=1,nc
		s=1.0_prec/sqrt(sum(abs(Cc(kc,:))**2))
		Cc(kc,:)=s*Cc(kc,:)
		dc(kc)=s*dc(kc)
	end do
	end subroutine normalise_chief
