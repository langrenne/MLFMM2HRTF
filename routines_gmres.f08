! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: MIT
! ---------------------------------------------------------------------
!
! VERSION FFT
!
! ---------------------------------------------------------------------

	subroutine sub_gmres_FFT(icntl,cntl,info,rinfo,&
	Anear,ALU,jLU,jU,iwk,octree,NL,N,SPH,alpha,k,lwork,work)

	use module_gmres
	use module_Octree
	use module_MLFMM
	use module_BEM
	use parametres
	implicit none

	integer, intent(in) :: lwork,NL,N(NL),iwk,jU(NE),jLU(iwk)
!	integer, intent(in) :: iA(num_K0),jA(num_K0)

	integer :: revcom, colx, coly, colz, nbscal,i1,i2
	integer :: irc(5), icntl(8), info(3)
	type(etat_gmres) :: etat          ! etat du GMRES, local a ce sous-programme (thread)
	double precision :: cntl(5), rinfo(2),k
	double complex :: work(lwork),test(NE),var1(NE),alpha,Anear(num_K0),ALU(iwk),temp(NE)

	type(tree) :: octree(NL)
	type(FMM) :: SPH(NL)
	character(len=1) :: answer


    10	 	continue
		call gmres_rc(NE,restart,lwork,work,irc,icntl,cntl,info,rinfo,etat)
		revcom = irc(1)
		colx = irc(2)
		coly = irc(3)
		colz = irc(4)
		nbscal = irc(5)

		if (revcom.eq.1) then
		  !* perform the matrix vector product
		  ! work(colz) <-- A * work(colx)

		  test=work(colx:colx+NE-1)
		  ! BEM
		  temp=dcmplx(0.0d0,0.0d0)
          do i2=1,num_K0
            temp(iA(i2))=temp(iA(i2))+Anear(i2)*test(jA(i2))
          end do
          work(colz:colz+NE-1)=temp
		 ! FMM en pression
		  !call BM_2_MLFMM_FFT(test,octree,NL,N,SPH,alpha,k,var1)
		  call MLFMM_FFT_p(test,octree,NL,N,SPH,alpha,k,var1)
		  work(colz:colz+NE-1)=work(colz:colz+NE-1)+var1
		  goto 10
		elseif (revcom.eq.2) then
		  !* perform the left preconditioning
		  !* work(colz) <-- M_1^{-1} * work(colx)
		  !write(*,*) 'passe left preconditionning'
		  test=work(colx:colx+NE-1)
		  call lusol(NE, test, var1, ALU, jLU, jU,iwk)
		  work(colz:colz+NE-1)=var1
		  goto 10
		elseif (revcom.eq.3) then
		  !* perform the right preconditioning
		  !* work(colz) <-- M_2^{-1} * work(colx)
		  ! MEME ILUT que la branche gauche. C'etait une simple copie,
		  ! si bien que FMM_pre=2 revenait a tourner SANS preconditionneur.
		  ! Interet du cote droit : le critere d'arret de la GMRES porte
		  ! alors sur le VRAI residu et non sur le residu preconditionne.
		  !* perform the left preconditioning
		  !* work(colz) <-- M_1^{-1} * work(colx)
		  !write(*,*) 'passe left preconditionning'
		  test=work(colx:colx+NE-1)
		  call lusol(NE, test, var1, ALU, jLU, jU,iwk)
		  work(colz:colz+NE-1)=var1
		  goto 10

		endif


	end subroutine sub_gmres_FFT
! ---------------------------------------------------------------------
!
! VERSION HYB
!
! ---------------------------------------------------------------------

	subroutine sub_gmres_HYB(icntl,cntl,info,rinfo,&
	Anear,ALU,jLU,jU,iwk,octree,NL,N,SPH,SPH0,alpha,k,lwork,work)

	use module_gmres
	use module_Octree
	use module_MLFMM
	use module_BEM
	use parametres
	implicit none

	integer, intent(in) :: lwork,NL,N(NL),iwk,jU(NE),jLU(iwk)
!	integer, intent(in) :: iA(num_K0),jA(num_K0)

	integer :: revcom, colx, coly, colz, nbscal,i1,i2
	integer :: irc(5), icntl(8), info(3)
	type(etat_gmres) :: etat          ! etat du GMRES, local a ce sous-programme (thread)
	double precision :: cntl(5), rinfo(2),k
	double complex :: work(lwork),test(NE),var1(NE),alpha,Anear(num_K0),ALU(iwk),temp(NE)

	type(tree) :: octree(NL)
	type(FMM) :: SPH(NL)
	type(FMM0) :: SPH0
	character(len=1) :: answer


    10	 	continue
		call gmres_rc(NE,restart,lwork,work,irc,icntl,cntl,info,rinfo,etat)
		revcom = irc(1)
		colx = irc(2)
		coly = irc(3)
		colz = irc(4)
		nbscal = irc(5)

		if (revcom.eq.1) then
		  !* perform the matrix vector product
		  ! work(colz) <-- A * work(colx)

		  test=work(colx:colx+NE-1)
		  ! BEM
		  temp=dcmplx(0.0d0,0.0d0)
          do i2=1,num_K0
            temp(iA(i2))=temp(iA(i2))+Anear(i2)*test(jA(i2))
          end do
          work(colz:colz+NE-1)=temp
		 ! FMM en pression
		  !call BM_2_MLFMM_FFT(test,octree,NL,N,SPH,alpha,k,var1)
		  call MLFMM_HYB_p(test,octree,NL,N,SPH,SPH0,alpha,k,var1)
		  work(colz:colz+NE-1)=work(colz:colz+NE-1)+var1
		  goto 10
		elseif (revcom.eq.2) then
		  !* perform the left preconditioning
		  !* work(colz) <-- M_1^{-1} * work(colx)
		  !write(*,*) 'passe left preconditionning'
		  test=work(colx:colx+NE-1)
		  call lusol(NE, test, var1, ALU, jLU, jU,iwk)
		  work(colz:colz+NE-1)=var1
		  goto 10
		elseif (revcom.eq.3) then
		  !* perform the right preconditioning
		  !* work(colz) <-- M_2^{-1} * work(colx)
		  ! MEME ILUT que la branche gauche. C'etait une simple copie,
		  ! si bien que FMM_pre=2 revenait a tourner SANS preconditionneur.
		  ! Interet du cote droit : le critere d'arret de la GMRES porte
		  ! alors sur le VRAI residu et non sur le residu preconditionne.
		  !* perform the left preconditioning
		  !* work(colz) <-- M_1^{-1} * work(colx)
		  !write(*,*) 'passe left preconditionning'
		  test=work(colx:colx+NE-1)
		  call lusol(NE, test, var1, ALU, jLU, jU,iwk)
		  work(colz:colz+NE-1)=var1
		  goto 10

		endif


	end subroutine sub_gmres_HYB
! ---------------------------------------------------------------------
!
! VERSION Z sur toute la surface de la structure
!
! ---------------------------------------------------------------------

	subroutine sub_gmres_Z_FFT(icntl,cntl,info,rinfo,&
	Anear,ALU,jLU,jU,iwk,octree,NL,N,SPH,alpha,k,A_red,lwork,work)

	use module_gmres
	use module_Octree
	use module_MLFMM
	use module_BEM
	use parametres
	implicit none

	integer, intent(in) :: lwork,NL,N(NL),iwk,jU(NE),jLU(iwk)
!	integer, intent(in) :: iA(num_K0),jA(num_K0)

	integer :: revcom, colx, coly, colz, nbscal,i1,i2
	integer :: irc(5), icntl(8), info(3)
	type(etat_gmres) :: etat          ! etat du GMRES, local a ce sous-programme (thread)
	double precision :: cntl(5), rinfo(2),k
	double complex :: work(lwork),test(NE),var1(NE),alpha,Anear(num_K0),ALU(iwk),temp(NE),A_red

	type(tree) :: octree(NL)
	type(FMM) :: SPH(NL)
	character(len=1) :: answer

	!write(*,*) 'A_red=',A_red
	!write(*,*) 'num_K0=',num_K0
!	Adm = ji * k * A_red
!	alpha=dcmplx(0.0d0,0.0d0)
!	Adm = (1.0d0-ji)/2.0d0 * k * (0.99985d0*dsqrt(2*lv/k)+(Gam-1.0d0)*dsqrt(2*lh/k)) ! pour cone

    10	 	continue
!		write(*,*) 'lwork=',lwork
!		write(*,*) 'NE=',NE
		!write(*,*) 'i1 =',i1
!		write(*,*) 'restart =',restart
		call gmres_rc(NE,restart,lwork,work,irc,icntl,cntl,info,rinfo,etat)
		revcom = irc(1)
		colx = irc(2)
		coly = irc(3)
		colz = irc(4)
		nbscal = irc(5)
!		write(*,*) 'revcom=',revcom,'colx=',colx,'coly=',coly,'colz=',colz,'nbscal=',nbscal
!		read(*,*) answer

		if (revcom.eq.1) then
		  !* perform the matrix vector product
		  ! work(colz) <-- A * work(colx)
!		  write(*,*) 'passe matvec'
		  test=work(colx:colx+NE-1)
		  ! BEM
		  temp=dcmplx(0.0d0,0.0d0)
          do i2=1,num_K0
!			write(*,*) 'iA(',i2,')=',iA(i2)
!			write(*,*) 'jA(',i2,')=',jA(i2)
!			write(*,*) 'Anear(',i2,')=',Anear(i2)
!			write(*,*) 'temp(',i2,')=',temp(i2)
            temp(iA(i2))=temp(iA(i2))+Anear(i2)*test(jA(i2))
          end do
          work(colz:colz+NE-1)=temp

		  ! FMM en pression
		  !call BM_2_MLFMM_FFT(test,octree,NL,N,SPH,alpha,k,var1)
!		  call MLFMM_FFT_p(test,octree,NL,N,SPH,alpha,k,var1)
!		  work(colz:colz+NE-1)=work(colz:colz+NE-1)+var1

		  ! FMM en gradient
!		  var1=dcmplx(0.0d0,0.0d0)
!		  test=dcmplx(0.0d0,0.0d0)

!			test=work(colx:colx+NE-1) * ji * k * A_red
!		  	call BM_dnp_MLFMM_FFT(test,octree,NL,N,SPH,alpha,k,var1)
!			call MLFMM_FFT_dnp(test,octree,NL,N,SPH,alpha,k,var1)

			temp=-work(colx:colx+NE-1) * ji * k * A_red   ! dnp = -i k A_red p (voir sub_gmres_Z_HYB)
			call MLFMM_FFT_p_dnp(test,temp,octree,NL,N,SPH,alpha,k,var1)
			work(colz:colz+NE-1)=work(colz:colz+NE-1) + var1

!		  write(*,*) 'passe mul :',i1
			goto 10
		elseif (revcom.eq.2) then
		  !* perform the left preconditioning
		  !* work(colz) <-- M_1^{-1} * work(colx)
		  test=work(colx:colx+NE-1)
		  call lusol(NE, test, var1, ALU, jLU, jU,iwk)
		  work(colz:colz+NE-1)=var1
		  goto 10
		elseif (revcom.eq.3) then
		  !* perform the right preconditioning
		  !* work(colz) <-- M_2^{-1} * work(colx)
		  ! MEME ILUT que la branche gauche. C'etait une simple copie,
		  ! si bien que FMM_pre=2 revenait a tourner SANS preconditionneur.
		  ! Interet du cote droit : le critere d'arret de la GMRES porte
		  ! alors sur le VRAI residu et non sur le residu preconditionne.
		  !* perform the left preconditioning
		  !* work(colz) <-- M_1^{-1} * work(colx)
		  test=work(colx:colx+NE-1)
		  call lusol(NE, test, var1, ALU, jLU, jU,iwk)
		  work(colz:colz+NE-1)=var1
		  goto 10

		endif


	end subroutine sub_gmres_Z_FFT
! ---------------------------------------------------------------------
!
! VERSION Z sur toute la surface de la structure
! Version Hybride BF et HF
!
! ---------------------------------------------------------------------

	subroutine sub_gmres_Z_HYB(icntl,cntl,info,rinfo,&
	Anear,ALU,jLU,jU,iwk,octree,NL,N,SPH,SPH0,alpha,k,A_red,lwork,work)

	use module_gmres
	use module_Octree
	use module_MLFMM
	use module_BEM
	use parametres
	implicit none

	integer, intent(in) :: lwork,NL,N(NL),iwk,jU(NE),jLU(iwk)
!	integer, intent(in) :: iA(num_K0),jA(num_K0)

	integer :: revcom, colx, coly, colz, nbscal,i1,i2
	integer :: irc(5), icntl(8), info(3)
	type(etat_gmres) :: etat          ! etat du GMRES, local a ce sous-programme (thread)
	double precision :: cntl(5), rinfo(2),k
	double complex :: work(lwork),test(NE),var1(NE),alpha,Anear(num_K0),ALU(iwk),temp(NE),A_red

	type(tree) :: octree(NL)
	type(FMM) :: SPH(NL)
	type(FMM0) :: SPH0
	character(len=1) :: answer

	!write(*,*) 'A_red=',A_red
	!write(*,*) 'num_K0=',num_K0
!	Adm = ji * k * A_red
!	alpha=dcmplx(0.0d0,0.0d0)
!	Adm = (1.0d0-ji)/2.0d0 * k * (0.99985d0*dsqrt(2*lv/k)+(Gam-1.0d0)*dsqrt(2*lh/k)) ! pour cone

    10	 	continue
!		write(*,*) 'lwork=',lwork
!		write(*,*) 'NE=',NE
		!write(*,*) 'i1 =',i1
!		write(*,*) 'restart =',restart
		call gmres_rc(NE,restart,lwork,work,irc,icntl,cntl,info,rinfo,etat)
		revcom = irc(1)
		colx = irc(2)
		coly = irc(3)
		colz = irc(4)
		nbscal = irc(5)
!		write(*,*) 'revcom=',revcom,'colx=',colx,'coly=',coly,'colz=',colz,'nbscal=',nbscal
!		read(*,*) answer

		if (revcom.eq.1) then
		  !* perform the matrix vector product
		  ! work(colz) <-- A * work(colx)
!		  write(*,*) 'passe matvec'
		  test=work(colx:colx+NE-1)
		  ! BEM
		  temp=dcmplx(0.0d0,0.0d0)
          do i2=1,num_K0
!			write(*,*) 'iA(',i2,')=',iA(i2)
!			write(*,*) 'jA(',i2,')=',jA(i2)
!			write(*,*) 'Anear(',i2,')=',Anear(i2)
!			write(*,*) 'temp(',i2,')=',temp(i2)
            temp(iA(i2))=temp(iA(i2))+Anear(i2)*test(jA(i2))
          end do
          work(colz:colz+NE-1)=temp

		  ! FMM en pression
		  !call BM_2_MLFMM_FFT(test,octree,NL,N,SPH,alpha,k,var1)
!		  call MLFMM_FFT_p(test,octree,NL,N,SPH,alpha,k,var1)
!		  work(colz:colz+NE-1)=work(colz:colz+NE-1)+var1

		  ! FMM en gradient
!		  var1=dcmplx(0.0d0,0.0d0)
!		  test=dcmplx(0.0d0,0.0d0)

!			test=work(colx:colx+NE-1) * ji * k * A_red
!		  	call BM_dnp_MLFMM_FFT(test,octree,NL,N,SPH,alpha,k,var1)
!			call MLFMM_FFT_dnp(test,octree,NL,N,SPH,alpha,k,var1)

			! IMPEDANCE : dnp = -i k A_red p (meme signe que Apnear de fait_tout_HRTF_DIFF_T0)
			! CONVENTION (2026-10-06) : normale vers le FLUIDE, e^{-iwt}, Z = p/v_paroi (vitesse entrant dans la paroi)
			! => v.n = -p/Z et dnp = -i k A_red p (NumCalc/Mesh2HRTF identique) ; avant : +i k A_red p (paroi active)
			temp=-work(colx:colx+NE-1) * ji * k * A_red
			call MLFMM_p_dnp(test,temp,octree,NL,N,SPH,SPH0,alpha,k,var1)
			work(colz:colz+NE-1)=work(colz:colz+NE-1) + var1

!		  write(*,*) 'passe mul :',i1
			goto 10
		elseif (revcom.eq.2) then
		  !* perform the left preconditioning
		  !* work(colz) <-- M_1^{-1} * work(colx)
		  test=work(colx:colx+NE-1)
		  call lusol(NE, test, var1, ALU, jLU, jU,iwk)
		  work(colz:colz+NE-1)=var1
		  goto 10
		elseif (revcom.eq.3) then
		  !* perform the right preconditioning
		  !* work(colz) <-- M_2^{-1} * work(colx)
		  ! MEME ILUT que la branche gauche. C'etait une simple copie,
		  ! si bien que FMM_pre=2 revenait a tourner SANS preconditionneur.
		  ! Interet du cote droit : le critere d'arret de la GMRES porte
		  ! alors sur le VRAI residu et non sur le residu preconditionne.
		  !* perform the left preconditioning
		  !* work(colz) <-- M_1^{-1} * work(colx)
		  test=work(colx:colx+NE-1)
		  call lusol(NE, test, var1, ALU, jLU, jU,iwk)
		  work(colz:colz+NE-1)=var1
		  goto 10

		endif


	end subroutine sub_gmres_Z_HYB
! ---------------------------------------------------------------------
!
! VERSION CHIEF (essai 2026-10-09, branche chief) : systeme borde carre
!
!   [ A  Y ] [ p      ]   [ b ]      A : operateur de surface (Anear + MLFMM),
!   [ C  0 ] [ lambda ] = [ d ]      C : nc lignes interieures (dense), Y = C^H
!
! Taille n = NE + nc. lambda n'a pas de sens physique : ~0 si les lignes CHIEF
! sont compatibles avec la surface ; elles levent la singularite de A aux
! frequences irregulieres. impedance = .false. : Ray (MLFMM_HYB_p) ;
! .true. : Diff (MLFMM_p_dnp avec dnp = -i k A_red p, comme sub_gmres_Z_HYB).
! Preconditionneur par blocs EXACT pour la bordure (M = ILUT de Anear) :
!   [M Y ; C 0]^-1 via le complement de Schur S = C M^-1 Y (nc x nc), calcule une
!   fois : Z = M^-1 Y (nc resolutions ILUT). Pour r = (r1, r2) :
!   y1 = M^-1 r1 ; lambda = S^-1 (C y1 - r2) ; x = y1 - Z lambda.
!
! ---------------------------------------------------------------------

	subroutine sub_gmres_chief(icntl,cntl,info,rinfo,&
	Anear,ALU,jLU,jU,iwk,octree,NL,N,SPH,SPH0,alpha,k,A_red,impedance,nc,Cc,lwork,work)

	use module_gmres
	use module_Octree
	use module_MLFMM
	use module_BEM
	use parametres
	implicit none

	integer, intent(in) :: lwork,NL,N(NL),iwk,jU(NE),jLU(iwk),nc
	logical, intent(in) :: impedance
	double complex, intent(in) :: Cc(nc,NE)

	integer :: revcom, colx, coly, colz, nbscal, i2, ntot
	integer :: irc(5), icntl(8), info(3)
	type(etat_gmres) :: etat          ! etat du GMRES, local a ce sous-programme (thread)
	double precision :: cntl(5), rinfo(2),k
	double complex :: work(lwork),test(NE),var1(NE),var2(NE),alpha,Anear(num_K0),ALU(iwk),temp(NE),A_red

	type(tree) :: octree(NL)
	type(FMM) :: SPH(NL)
	type(FMM0) :: SPH0
	double complex, allocatable :: Z(:,:), S(:,:), rl(:)
	integer, allocatable :: ipiv(:)
	integer :: kc, info_l

	ntot = NE + nc
	! preconditionneur de la bordure : Z = M^-1 C^H, S = C Z (LU)
	allocate(Z(NE,nc), S(nc,nc), rl(nc), ipiv(nc))
	do kc=1,nc
		test=conjg(Cc(kc,:))
		call lusol(NE, test, var1, ALU, jLU, jU,iwk)
		Z(:,kc)=var1
	end do
	S=matmul(Cc,Z)
	call zgetrf(nc,nc,S,nc,ipiv,info_l)
	if (info_l /= 0) write(log_unit,*) 'CHIEF: singular Schur complement, info =', info_l
    10	 	continue
		call gmres_rc(ntot,restart,lwork,work,irc,icntl,cntl,info,rinfo,etat)
		revcom = irc(1)
		colx = irc(2)
		coly = irc(3)
		colz = irc(4)
		nbscal = irc(5)

		if (revcom.eq.1) then
		  ! work(colz) <-- [A Y ; C 0] * work(colx)
		  test=work(colx:colx+NE-1)
		  temp=dcmplx(0.0d0,0.0d0)
          do i2=1,num_K0
            temp(iA(i2))=temp(iA(i2))+Anear(i2)*test(jA(i2))
          end do
		  if (impedance) then
			var2=-test * ji * k * A_red          ! dnp = -i k A_red p (convention de sub_gmres_Z_HYB)
			call MLFMM_p_dnp(test,var2,octree,NL,N,SPH,SPH0,alpha,k,var1)
		  else
			call MLFMM_HYB_p(test,octree,NL,N,SPH,SPH0,alpha,k,var1)
		  endif
		  temp=temp+var1
		  ! Y lambda, Y = C^H
		  temp=temp+matmul(transpose(conjg(Cc)),work(colx+NE:colx+ntot-1))
		  work(colz:colz+NE-1)=temp
		  ! C p
		  work(colz+NE:colz+ntot-1)=matmul(Cc,test)
		  goto 10
		elseif (revcom.eq.2 .or. revcom.eq.3) then
		  ! preconditionneur (gauche ou droite) : [M Y ; C 0]^-1 par le complement de Schur
		  test=work(colx:colx+NE-1)
		  call lusol(NE, test, var1, ALU, jLU, jU,iwk)
		  rl=matmul(Cc,var1)-work(colx+NE:colx+ntot-1)
		  call zgetrs('N',nc,1,S,nc,ipiv,rl,nc,info_l)
		  work(colz:colz+NE-1)=var1-matmul(Z,rl)
		  work(colz+NE:colz+ntot-1)=rl
		  goto 10
		endif
	deallocate(Z, S, rl, ipiv)

	end subroutine sub_gmres_chief
