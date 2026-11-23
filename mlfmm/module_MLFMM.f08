! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: MIT
	module module_MLFMM
	use module__Precision
	use parametres
	use module_Acoustiques
	use module_Octree
	use module_BEM
	use module_fft     ! cfftmi/cfftmf/cfftmb (modern_fftpack)
	use omp_lib
	implicit none

	integer :: maxi_BF
	! Ordre des operations des translations HF entre niveaux d'ORDRES DIFFERENTS
	! (upward_HF, L2L de down_HF), 2026-09-29 :
	!   .false. : ancien ordre (dephasage sur la grille de l'enfant puis
	!             interpolation / filtrage puis dephasage) : repliement, le
	!             produit multipole x exp(ik s.d) depasse la bande de la grille
	!             enfant (propagation depuis le niveau 2 : 2e-5 au lieu de 7e-7)
	!   .true.  : interpoler PUIS dephaser (montee), dephaser PUIS filtrer
	!             (descente), sur la grille du parent : ordre MLFMA classique
	!   Cle interp_solveur (parametres, 2026-10-30) : 1 PAR DEFAUT = ordre
	!   classique aussi dans le solveur (variable partagee, lue avant la boucle
	!   parallele) ; 0 = ancien ordre. HRTF : meme solution a 1e-3 dB, +2 % de temps.
	!   Reglage (2026-09-29) : .false. dans le solveur (gain de precision non
	!   mesurable en T6, surcout 6 a 70 %) ; .true. pendant la montee des
	!   PROPAGATIONS (propagation_FMM_only, propagation_hierarchique), ou le
	!   repliement dominait (2e-5 -> 2,6e-8). THREADPRIVATE : la propagation
	!   tourne dans la boucle parallele en frequence, a cote des GMRES des
	!   autres threads, et ne doit changer que son propre calcul.
	logical :: interp_puis_dephase = .false.
	!$omp threadprivate(interp_puis_dephase)
	! liberer V1 apres le second membre (mode 1) du chemin T0 adaptatif : GMRES n'utilise que V2 (Ray
	! rigide) ; V1 est recalcule pour la propagation (mode 3). Pose par main_HRTF_Ray (2026-11-06).
	logical :: liberer_V1 = .false.
	! Diff (2026-11-06) : le produit n'utilise V1 qu'avec dnp = -i k A_red p ; V2.p - V1.dnp = Vc.p avec
	! Vc = V2 + i k A_red V1 (BF) ou V2 - A_red V1 (HF, normalisation des injections) : un seul tableau.
	! fusion_Vc : pose par main_HRTF_Diff ; beta_Vc = A_red de la frequence (fait_tout_HRTF_DIFF_T0).
	logical :: fusion_Vc = .false.
	complex(prec) :: beta_Vc = (0.0_prec, 0.0_prec)
	!$omp threadprivate(beta_Vc)
	! propagation_hierarchique : cellules acceptees et elements integres
	! directement au dernier appel hors region parallele (diagnostic)
	integer(8) :: n_prop_cellules = 0, n_prop_directs = 0
	! methode_propagation : voir parametres (lire_parametres_MLFMM)
	! variables frequency dependant : à deallocated après traitement
	Type FMM

	real(prec), dimension(:,:), allocatable :: Ps
	real(prec), dimension(:), allocatable  :: wp, theta, IR
	complex(prec), dimension(:), allocatable  :: in,out,V1,V2,outBF,W1,W2,inBF
	complex(prec), dimension(:,:), allocatable  :: tab_M2L
	real(prec), dimension(:,:,:), allocatable :: mat_Q
	integer, dimension(316,3) :: tab_xyz
	integer, dimension(:,:), allocatable  :: tab1_xyz
	integer, dimension(:), allocatable  :: iv,jv,iw,jw
	integer   :: BF
	integer :: nnzv,nlv,ncv,nnzw,nlw,ncw
	real(prec), dimension(:), allocatable :: trans_M2M,trans_L2L
	complex(prec), dimension(:,:), allocatable :: trans_M2L
	integer, dimension(:,:), allocatable :: trans_M2L_i,trans_M2L_j
	integer, dimension(:), allocatable :: trans_M2M_i,trans_M2M_j,rot_num,trans_L2L_i,trans_L2L_j
	real(prec), dimension(15) :: tab_T_M2L
	integer :: tab_T_M2M
	complex(prec), dimension(:,:), allocatable :: rot_M2M,rot_M2L
	integer, dimension(:,:), allocatable :: rot_M2M_i,rot_M2M_j
	integer, dimension(:,:), allocatable :: rot_M2L_i,rot_M2L_j
	real(prec), dimension(290,3) :: tab_R_M2L
	real(prec), dimension(8,3) :: tab_R_M2M

	end type FMM

	! variables non frequency dependant : rot et vecteur translation
	! on ne deallocate pas
	type FMM0

	complex(prec), dimension(:,:), allocatable :: rot_M2M,rot_M2L
	integer, dimension(:,:), allocatable :: rot_M2M_i,rot_M2M_j
	integer, dimension(:,:), allocatable :: rot_M2L_i,rot_M2L_j
	integer, dimension(:), allocatable :: rot_num
	real(prec), dimension(:,:),allocatable :: tab_T_M2L
	real(prec), dimension(:),allocatable :: tab_T_M2M
	real(prec),allocatable :: b(:,:),a(:,:),b_ex(:,:)

	end type FMM0

	!PRIVATE
	PUBLIC :: init_FMM_T0, ordres_FMM_T0, memoire_thread_T0, liberer_V1, fusion_Vc, beta_Vc

	! Interfaces génériques : résolution par présence/absence de SPH0 (type FMM0)
	interface pass_upward
		module procedure pass_upward_FFT   ! sans SPH0
		module procedure pass_upward_HYB   ! avec SPH0
	end interface pass_upward

	interface pass_downward
		module procedure pass_downward_FFT ! sans SPH0
		module procedure pass_downward_HYB ! avec SPH0
	end interface pass_downward

	! Interfaces génériques : calcul des moments (surface BEM)
	interface moments_FMM
		module procedure moments_FMM_HF    ! sans SPH0 (5 args)
		module procedure moments_FMM_BF    ! avec SPH0 (6 args)
	end interface moments_FMM

	interface moments_FMM_propa
		module procedure moments_FMM_propa_HF ! sans SPH0 (5 args)
		module procedure moments_FMM_propa_BF ! avec SPH0 (6 args)
	end interface moments_FMM_propa

	! Interfaces génériques : calcul des expansions locales (surface BEM)
	interface expansions_BM
		module procedure expansions_BM_HF  ! sans SPH0 (6 args)
		module procedure expansions_BM_BF  ! avec SPH0 (7 args)
	end interface expansions_BM

	interface expansions_BM_propa
		module procedure expansions_BM_propa_HF ! sans SPH0 (5 args)
		module procedure expansions_BM_propa_BF ! avec SPH0 (6 args)
	end interface expansions_BM_propa

	! Interface générique : upward (passage montant)
	interface upward
		module procedure upward_HF  ! sans SPH0 (6 args)
		module procedure upward_BF  ! avec SPH0 (6 args, type FMM0)
	end interface upward

	! Interfaces génériques : produit MLFMM (FFT sans SPH0, HYB avec SPH0)
	interface MLFMM_dnp
		module procedure MLFMM_FFT_dnp ! sans SPH0 (8 args)
		module procedure MLFMM_HYB_dnp ! avec SPH0 (9 args)
	end interface MLFMM_dnp

	interface MLFMM_p
		module procedure MLFMM_FFT_p   ! sans SPH0 (8 args)
		module procedure MLFMM_HYB_p   ! avec SPH0 (9 args)
	end interface MLFMM_p

	interface MLFMM_p_dnp
		module procedure MLFMM_FFT_p_dnp ! sans SPH0 (9 args)
		module procedure MLFMM_HYB_p_dnp ! avec SPH0 (10 args)
	end interface MLFMM_p_dnp

	contains

! --------------------------------------------------------------------------------
!
! VERSION T0 pour FFT, BF ou version hybride HYB
! précalcul de ce qu'il faut
! HFBF = 0 ! tout en BF
! HFBF = 1 ! tout en HF
! HFBF = 2 ! Hybride
!
! pour epsil2 donné
!
! --------------------------------------------------------------------------------
	subroutine init_FMM_T0(octree,NL,k,N,SPH,SPH0,HFBF,epsil2)
! initialisation variable pour MLFMM par FFT
	implicit none
! variables entrée
	integer, intent(inout) :: NL
	integer, intent(inout) :: HFBF
	real(prec), intent(inout) :: k,epsil2
	type(tree), intent(inout) :: octree(NL)
! variables sortie
	integer,intent(out) :: N(NL)
	type(FMM), intent(out) :: SPH(NL)
	type(FMM0),intent(inout) :: SPH0
! variable locales
	integer :: i1, num
!	integer :: n2p,n2m,m0,m11,m0p,m0m
!	integer :: n2,m2,n1,m1,Lmax,indL,ind1,ind2,ind3,ind33
	real(prec) :: r,sigma,kr,ka,kb,ka1,kb1
!	real(prec), allocatable :: Pt(:,:),wp(:),theta0(:),PL(:,:),cos_phi(:,:),cos_phi0(:)
	real(prec) :: plow, phigh
	complex(prec), dimension(:), allocatable  :: M2L
	integer :: n1, n2, res1, i2, nx, ny, nz, npu_gro, numT
	real(prec) :: cote
	real(prec), dimension(1:3) :: vt
!	real(prec), allocatable :: Ps_leb(:,:)
!	complex(prec), dimension(:), allocatable :: M2L
    integer, dimension(:), allocatable :: M2L_i,M2L_j
	integer :: BF_niv(NL)

	! DEBUG : character(len=1) :: answer

	maxi_BF=0
	! ordres et choix BF/HF par niveau : routine separee (2026-10-02), partagee avec
	! l'estimation de la RAM par frequence (memoire_thread_T0)
	call ordres_FMM_T0(octree,NL,k,HFBF,epsil2,N,BF_niv)
	SPH(:)%BF=BF_niv

    ! tout en HF
!    SPH(:)%BF=0
    ! tout en BF
!    SPH(:)%BF=1

	! allocations pour formulation diagonale ou non
	do i1=1,NL
!        if (SPH(i1)%BF==0 .or. i1==NL) then
!		write(*,*) 'allocation niveau : ',i1
            allocate(SPH(i1)%Ps((2*N(i1)+1)*(N(i1)+1),3))
            allocate(SPH(i1)%wp((2*N(i1)+1)*(N(i1)+1)))
            allocate(SPH(i1)%theta(N(i1)+1))
!   		write(*,*) 'fin allocation niveau : ',i1
            call gauss_sphere(NL,SPH,N,i1)
!        endif
	enddo

!	% paramêtres de translation Interpolation version diagonale !
	do i1=1,NL-1
!        write(*,*) 'boucle i1=',i1
	    !write(*,*) i1
            if (SPH(i1+1)%BF==0) then

                n1=(N(i1)+1)
                n2=(N(i1+1)+1)
!                call coef_IR_v2(NL,N,n1,n2,SPH,i1)  ! modif 2022
                call coef_IR_v2(NL,N,n1,n2,SPH,i1)   ! on va jusqu'à l'ordre le plus grand et on tronquera dans L2L
	    !write(*,*) SPH(i1)%IR(:)
!	    read(*,*) answer
            endif
	enddo
!
!---------------------------------------------------------------------
!
! maillage Lebedev dernier niveau pour test intégration HF2BF
!
!---------------------------------------------------------------------
!
!    ray=0.1_prec
!    if (SPH(NL)%BF==1) then
!        call getLebedevSphere(N(NL),Leb)
!        ! changement de rayon
!        allocate(Ps_leb(Leb%n,3))
!        call car2sph(Leb%Pc,Ps_leb,Leb%n)
!        Ps_leb(:,1)=ray
!        call sph2car(Leb%Pc,Ps_leb,Leb%n)
!    endif
!
!---------------------------------------------------------------------
!
! précalcul des translations formalisme diagonal
!
!---------------------------------------------------------------------
!
    do i1=1,NL
		if (SPH(i1)%BF==0) then
 	  npu_gro=size(SPH(i1)%Ps,1)  !% nbre de points sur sph unitaire niveau

 	  allocate(SPH(i1)%tab_M2L(npu_gro,0:342))  ! 7^3 au lieu de 316, mais c'est plus simple
	  SPH(i1)%tab_M2L=cmplx(0.0_prec,0.0_prec, kind=prec)
      if (allocated(M2L)) deallocate(M2L)
      allocate(M2L(npu_gro))

	  cote=octree(i1)%cell(1)%length
	  num=0
	  do i2=0,342
        ! passage en base 7 pour nxnynz
        nz=mod(i2,7)
        res1=(i2-nz)/7
        ny=mod(res1,7)
        nx=(res1-ny)/7
        ! on retranche 3 pour aller de -3 à 3
        nx=nx-3
        ny=ny-3
        nz=nz-3
        ! test interaction
        if (abs(nx) > 1 .or. abs(ny) > 1 .or. abs(nz) > 1) then
            num=num+1
            ! calcul du vecteur de translation
            vt=real((/nx,ny,nz/),8)*cote
            !write(*,*) "vt=",vt
            call T_LR(NL,N(i1),i1,k,SPH,vt,M2L)
            SPH(i1)%tab_M2L(:,i2)=M2L
        end if
       end do
       endif
    end do
!
!---------------------------------------------------------------------
!
! précalcul des translations formalisme BF
! Les vecteurs |t| sont déjà calculés par init_BF_MLFMM_v1
!---------------------------------------------------------------------
!
    ! les "1" translations par niveau pour M2M et L2L
    do i1=NL,2,-1
        if (SPH(i1)%BF==1) then
            call M2M_BF(NL,k,N,SPH,SPH0,i1) ! matrice creuse de translation récurrence
        endif
    enddo

    ! les 15 translations par niveaux pour M2L!

    do i1=1,NL
        if (SPH(i1)%BF==1) then
        if (allocated(SPH(i1)%trans_M2L)) deallocate(SPH(i1)%trans_M2L)
        if (allocated(SPH(i1)%trans_M2L_i)) deallocate(SPH(i1)%trans_M2L_i)
        if (allocated(SPH(i1)%trans_M2L_j)) deallocate(SPH(i1)%trans_M2L_j)
        do i2=1,15
            call M2L_BF(M2L,M2L_i,M2L_j,numT,k,N(i1),SPH0%tab_T_M2L(i1,i2)) ! matrice creuse de translation récurrence
            if (i2==1) then
                allocate(SPH(i1)%trans_M2L(numT,15))
                allocate(SPH(i1)%trans_M2L_i(numT,15))
                allocate(SPH(i1)%trans_M2L_j(numT,15))
            endif
            SPH(i1)%trans_M2L(:,i2)=M2L
            SPH(i1)%trans_M2L_i(:,i2)=M2L_i
            SPH(i1)%trans_M2L_j(:,i2)=M2L_j
        enddo
        endif
    enddo

    end subroutine init_FMM_T0

! ------------------------------------------------------------------------
! Ordres de troncature N et choix BF (1) / HF (0) par niveau pour le nombre
! d'onde k (2026-10-02 : extrait de init_FMM_T0, logique INCHANGEE) : regle
! plow/phigh, seuil de bascule seuil_de_k(k), plancher Nmin_ord, N croissant
! vers le haut, BF monotone (BF aux niveaux fins), boost dN_BF sous k_boost.
! Sert a init_FMM_T0 et a l'estimation de la RAM (memoire_thread_T0).
! ------------------------------------------------------------------------
	subroutine ordres_FMM_T0(octree,NL,k,HFBF,epsil2,N,BF)
	implicit none
	integer, intent(in) :: NL, HFBF
	type(tree), intent(in) :: octree(NL)
	real(prec), intent(in) :: k, epsil2
	integer, intent(out) :: N(NL), BF(NL)
	integer :: i1
	real(prec) :: r,sigma,kr,ka,kb,ka1,kb1
	real(prec) :: plow, phigh
	real(prec) :: ph(NL)        ! phigh par niveau (re-test BF/HF apres le boost)
	real(prec) :: Nbasc         ! essai 10/11 : N du test de bascule (reel si bascule_reelle = 1)

	sigma=2.0_prec
    BF(:)=1
!	write(*,*) 'in init BF(',i1,')=', BF(:)
!	write(*,*) 'in init_FMM_T0'
!    write(*,*) 'SPH0%tab_T_M2M=',SPH0%tab_T_M2M(:)

!	write(*,*) 'init ordre'
	do i1=1,NL
!        write(*,*) 'niveau :',i1
        r=octree(i1)%cell(1)%length  / 2.0_prec * sqrt(3.0_prec) ! un demi diametre : rayon de sphère engolobant la cellule.
!        r=octree(i1)%cell(1)%length
        ! indicateurs
        kr=k*r
!        plow=1.0_prec-log(epsil2*(1-1.0_prec/sigma)**(3.0_prec/2.0_prec))
!        write(*,*) 'plow1=',plow
        plow=1.0_prec-plow_coef*log10(epsil2*(1-1.0_prec/sigma)**(3.0_prec/2.0_prec))   ! plow_coef : essai 09/11 (defaut 1)
!        write(*,*) 'plow2=',plow
        !plow=-1.0_prec/log(sigma) * log(epsil2*kr*(1-1.0_prec/sigma)**(3.0_prec/2.0_prec))-1.0_prec
!        plow=-10.0_prec/15000.0_prec*k*c/2.0_prec/pi+10.0_prec;
!        plow=6.0_prec
!        write(*,*) 'plow1=',plow
!        phigh=kr+ (kr)**(1.0_prec/3.0_prec)*0.7_prec*(3.0_prec*log(1.0_prec/(epsil2*sigma)))**(2.0_prec/3.0_prec)
        phigh=kr+ (kr)**(1.0_prec/3.0_prec)*(3.0_prec*log10(1.0_prec/(epsil2)))**(2.0_prec/3.0_prec)
!        phigh=kr+0.5_prec*(kr)**(1.0_prec/3.0_prec)
!        write(*,*) 'phigh=',phigh
        N(i1)=nint((plow**4.0_prec+phigh**4.0_prec)**(1.0_prec/4.0_prec))
        Nbasc=real(N(i1), prec)
        if (bascule_reelle == 1) Nbasc=(plow**4.0_prec+phigh**4.0_prec)**(1.0_prec/4.0_prec)   ! sans arrondi : N - phigh monotone
        ph(i1)=phigh
!	N(i1)=5
!        write(*,*) 'plow=',plow,' phigh=',phigh,' N0(',i1,')=',N(i1)

        ka=0.5_prec*((4.0_prec*(kr)**2+1)**(0.5_prec)-1.0_prec)
        ka1=0.5_prec*(kr)**2
!        ka=3.0_prec/(2.0_prec**(3.0_prec/2.0_prec)*(sigma-1.0_prec)**(3.0_prec/2.0_prec))*log(1.0_prec/(epsil2*sigma))
        kb=0.5_prec*((4.0_prec*(kr*sigma)**2+1)**(0.5_prec)-1.0_prec)
        kb1=0.5_prec*(kr*sigma)**2
!        write(*,*) 'kr=',kr,' ka=',ka,' ka1=',ka1,' kb=',kb,' kb1=',kb1

        ! (supprime) ancienne clause "cap" N=6+nint(phigh) si N-phigh>5 : rendait N
        ! non-monotone (l'ordre diminuait a frequence croissante). Remplacee par le
        ! plancher Nmin_ord applique apres la boucle.

        ! test bf ou hf
 !       write(*,*) 'niveau ',i1,' 2kr=',2.0_prec*kr
 !       if (2.0_prec*k*r > 0.175_prec) then
 !           BF(i1)=0
 !       endif
 !      if (N(i1) > 12) then
 !           BF(i1)=0
 !       endif
!        write(*,*) ' N1(',i1,')=',N(i1)
!        prot=nint(kr*sqrt(2.0_prec)+epsil2*log(kr*sqrt(2.0_prec)+pi))
!        write(*,*) 'prot=',prot
        ! 1er test : bascule HF a kr~0.49 (seuil 2.5 au lieu de 2 Gumerov) = point de coincidence
        ! BF=HF pour N=8/eps=1e-7 (~700 Hz pour L2). HF (diagonal) plus tot = plus rapide.
        ! Seuil parametrable (seuil_HF, parametres, defaut 2.5) depuis le 2026-09-30,
        ! par bande de frequence (HRTF) depuis le 2026-10-02.
        ! (cle kr_hf_min, essai 2026-11-03 : pas de HF sous un k r donne ; x4 a x6 en temps au-dessus
        ! de 10 kHz sur la tete sans gain de precision, retiree le 2026-11-20)
        if (Nbasc-phigh < seuil_de_k(k)) then
            BF(i1)=0
        endif
        ! 2ème test
!        if (real(N(i1), prec) > plow) then
!            BF(i1)=0
!        endif
        ! 2ème test
!        if (kr > 0.2_prec) then
!            BF(i1)=0
!        endif
        !N(i1)=int(phigh)
        !if (N(i1) < 4) N(i1)=4


!            if (SPH(i1)%clebsch==0) then
!                if (N(i1) > maxi_BF) maxi_BF=N(i1)  ! ordre maxi pour calcul rotation
!            endif
!            else
!                    SPH(i1)%clebsch=1  ! pour formulation diagonale
!            endif
!                SPH(i1)%clebsch=0
!                if (N(i1) > maxi_BF) maxi_BF=N(i1)  ! ordre maxi pour calcul rotation
!            endif
	enddo

	! plancher d'ordre : garantit N >= Nmin_ord a tous les niveaux (accuracy basses freq).
	! Applique apres le test BF/HF pour ne pas fausser la decision BF/HF.
	do i1=1,NL
		if (N(i1) < Nmin_ord) N(i1) = Nmin_ord
		N(i1) = N(i1) + dN_ord      ! increment d'ordre (parametres, 0 par defaut) : test de precision
	end do

	! raz tout en HF
	if (HFBF==1) then
        BF(:)=0
	endif
	! raz tout en BF pour les ordres inférieurs à Nmax pré-calculé dans le "main" programme (default 20)
  	if (HFBF==0) then
		do i1=1,NL
			if (N(i1) > 12) then
				BF(i1)=0
			else
				BF(i1)=1
			endif
		enddo
	endif
	!test
!	SPH(NL)%BF=1
!	write(*,*) 'BF=',BF(:)
	! test
	!SPH(1)%BF=0
    !test, N au moins croissant
	do i1=NL,2,-1
            if (N(i1-1)<N(i1)) then
                N(i1-1)=N(i1)
            endif
            ! Niveau BF ou HF
!            if (N(i1)>10) BF(i1)=0  ! les deux ligne supprimer pour test bf hf
	enddo
!    if (N(1)>10) SPH(1)%BF=0  ! la deuxième ligne

    ! BF/HF monotone : BF aux niveaux fins, HF aux niveaux grossiers, une seule transition.
    ! Une fois BF a un niveau, tous les niveaux plus fins le sont aussi (i1 croissant = plus fin).
    ! Corrige un motif non-monotone que le critere par-niveau (N-phigh<2) peut produire.
    do i1=1,NL-1
        if (BF(i1)==1) BF(i1+1)=1
    end do

    ! BOOST BF (2026-09-30) : + dN_BF a l'ordre des niveaux BF seulement. Octree
    ! adaptatif : chaque element est pose au ras de sa cellule (seuil x arete),
    ! ses points de Gauss debordent de la sphere supposee par la regle d'ordre ;
    ! en basses frequences N est au plancher et la troncature domine (id13 :
    ! 0,7 % sur l'operateur a 50-200 Hz). En HF les ordres sont grands et le
    ! debordement ne compte plus : on ne touche pas aux niveaux HF (sauf pour
    ! garder N croissant vers le haut). dN_BF = 0 : inchange. Seulement si
    ! k < k_boost (coupure en frequence : au-dessus l'erreur est deja faible et
    ! le boost couterait sur les niveaux fins, les plus peuples). id13, s1,6,
    ! erreur de l'operateur 200 Hz / 1 kHz / 2 kHz : 6,3e-3 / 1,9e-3 / 8e-4 ->
    ! +4 : 3,2e-4 / 1,1e-4 / 7e-5, produit FMM x2,3.
    if (dN_BF > 0 .and. k < k_boost) then
        do i1=1,NL
            if (BF(i1)==1) N(i1)=N(i1)+dN_BF
        end do
        do i1=NL,2,-1
            if (N(i1-1)<N(i1)) N(i1-1)=N(i1)
        end do
        ! Re-test BF/HF apres le boost (2026-11-02, diag516_20261102) : la monotonie ci-dessus
        ! remonte l'ordre des niveaux HF au-dessus de la transition (m7 serre, 516 Hz : niveau 2,
        ! cellule 62 mm, k r ~ 0,5, N 9 -> 13 ; idem a chaque niveau qui passe en HF : 281, 1031 Hz).
        ! La translation diagonale y perd des chiffres : le vrai residu reste > 1e-5 (epsilon2) alors que celui d'Arnoldi descend a 1e-13, GMRES
        ! s'arrete a iter_max. On applique donc le critere de bascule a l'ordre FINAL : un niveau
        ! HF devenu N - phigh >= seuil repasse en BF (stable a N eleve). Transition toujours
        ! unique (les niveaux sous un niveau BF le sont deja). Hybride seulement.
        ! Mesure (plancher_20261102, m7 serre, 42 f sous 2 kHz, epsilon2 1e-9) : plancher 1,7e-5 au plus
        ! avant, 3,5e-9 apres ; temps par iteration : -20 a +8 %, total -10 %.
        if (HFBF == 2) then
            do i1=1,NL
                if (BF(i1)==0 .and. real(N(i1),prec)-ph(i1) >= seuil_de_k(k)) BF(i1)=1
            end do
            do i1=1,NL-1
                if (BF(i1)==1) BF(i1+1)=1
            end do
        endif
    endif

    ! (essai 2026-11-03, krhf_20261103 : ordre + 1 aux seuls niveaux HF au lieu de kr_hf_min : bosse
    ! 2,9-3,7 kHz corrigee a moitie seulement, +4 % de temps ; c'est la translation HF a k r ~ 0,5,
    ! pas l'ordre. Cle dN_HF retiree.)
	end subroutine ordres_FMM_T0

! ------------------------------------------------------------------------
!
! calcul les poids et points de Gauss pour intégration sur sphère unitaire
!
! ------------------------------------------------------------------------
!
	subroutine gauss_sphere(NL,SPH,N,i1)
! calcul les points de gauss et poids associé pour une sphère à l'ordre NL
	!use parametres
	!use init_param_FMM
	implicit none
! variables d'entrée
	integer, intent(in) :: NL,i1
	integer,dimension(1:NL), intent(in) :: N(NL)
! variables de sorties
	real(prec), dimension(1:(N(i1)+1)*(2*N(i1)+1),1:3) :: Pt
	real(prec), dimension(1:(N(i1)+1)*(2*N(i1)+1)) :: wp
	real(prec), dimension(1:N(i1)+1) :: theta
! variables locales
	real(prec) :: phi(2*N(i1)+1), x(N(i1)+1), wth(N(i1)+1), wphi(2*N(i1)+1)
	real(prec), dimension(1:2*N(i1)+1,1:N(i1)+1) :: wi, wj
	real(prec), dimension(1:(N(i1)+1)*(2*N(i1)+1),1:3) :: Psph
	integer :: ij
	! DEBUG : character(len=1) :: answer
!	type FMM
!		real(prec), dimension(:,:), pointer :: Ps
!		real(prec), dimension(:), pointer :: wp, theta
!	end type FMM
!	type int_FMM
!		type(FMM), pointer, dimension(:) :: level
!	end type int_FMM

	type(FMM) :: SPH(NL)

!	write(*,*) 'niveau : ',i1,' sur ',NL
!	write(*,*) 'debut allocation'
!		allocate(SPH(i1)%Ps((2*N(i1)+1)*(N(i1)+1),3))
!		allocate(SPH(i1)%wp((2*N(i1)+1)*(N(i1)+1)))
!		allocate(SPH(i1)%theta(N(i1)+1))
!	write(*,*) 'fin allocation'

		phi=(/ ( (ij-1.e0_prec)*2.0_prec*pi/(2.0_prec*N(i1)+1.0_prec),ij=1,2*N(i1)+1) /)          !	[0:2*pi/(2*L+1):2*pi-2*pi/(2*L+1)];

!		% calcul des poids
		call gauleg(N(i1)+1, x, wth)
		wphi=spread(2.0_prec*pi/(2.0_prec*N(i1)+1),1,2*N(i1)+1)
		wi=spread(wth,1,2*N(i1)+1)
		wj=spread(wphi,2,N(i1)+1)
		wp(:)=reshape(wi,(/(N(i1)+1)*(2*N(i1)+1)/)) * reshape(wj,(/(N(i1)+1)*(2*N(i1)+1)/))

!		% calcul des points
		theta=acos(x);
		wi=spread(theta,1,2*N(i1)+1)
		wj=spread(phi,2,N(i1)+1)
!		write(*,*) wi(1,:)
		Psph(:,2)=reshape(wi,(/(N(i1)+1)*(2*N(i1)+1)/))
		Psph(:,3)=reshape(wj,(/(N(i1)+1)*(2*N(i1)+1)/))
		Psph(:,1)=1.e0_prec;
		!write(*,*) Psph(1:5,2),Psph(1:5,3)
		!read(*,*) answer
		call sph2car(Pt,Psph,(N(i1)+1)*(2*N(i1)+1))
!		write(*,*) 'debut Pt'
		SPH(i1)%Ps=Pt
!		write(*,*) 'ecriture wp'
	        SPH(i1)%wp=wp
!		write(*,*) 'ecriture theta'
	        SPH(i1)%theta=theta


	end subroutine gauss_sphere
! -----------------------------------------------------------------------------
!
! VERSION v2 pour FFT
!
! -----------------------------------------------------------------------------
	subroutine coef_IR_v2(NL,N,n1,n2,SPH,i1)
!	use init_param_FMM
!	use parametres
	implicit none
! variables d'entrée
	integer, intent(in) :: n1, n2, NL,i1
	integer, dimension(1:NL) :: N
! variable de sortie
	type(FMM),intent(inout) :: SPH(NL)

! variable locales
!	real(prec) :: PL(n1*n2,N(i1+1)+2),cos_phi(n1,n2),cos_phi0(n1*n2)
	real(prec) :: mat_Qnm(n1,n2)
	integer :: i2
	! DEBUG : character(len=1) :: answer
!    write(*,*) 'in coef_IR_v2',i1,n1,n2
    allocate(SPH(i1)%mat_Q(n1,n2,0:n2-1))   ! modif 2022
!    allocate(SPH(i1)%mat_Q(n1,n2,0:n1-1))   ! modif 2022
	!write(*,*) shape(SPH(i1)%mat_Q)
	!read(*,*) answer
	do i2=0,n2-1      ! modif 2022
!	do i2=0,n1-1       ! modif 2022
            !write(*,*) 'n1 = ',n1,' n2 = ',n2
	    call mat_Q(N(i1+1),i2,n2,SPH(i1+1)%theta,n1,SPH(i1)%theta,mat_Qnm)
	    SPH(i1)%mat_Q(:,:,i2)=mat_Qnm
	    !write(*,*) 'Qnm passe(',i2,') = ',mat_Qnm
	    !read(*,*) answer
	enddo

	end subroutine coef_IR_v2
! -----------------------------------------------------------------
!
! mat_Q pour M2M FFT
!
! -----------------------------------------------------------------
!
	subroutine mat_Q(n,m,n1,theta1,n2,theta2,mat_Qnm)
	! calcul la matrice des Qnm de changement d'ordre pour M2M et L2L
	integer :: n,m,n1,n2
	real(prec) :: theta1(n1),theta2(n2)
	real(prec) :: mat_Qnm(n2,n1)

	real(prec) :: Qnm1(m:n,n1),Qnm2(m:n,n2)

	call coef_Q(n,m,n1,theta1,Qnm1)
	call coef_Q(n,m,n2,theta2,Qnm2)

	mat_Qnm=matmul(transpose(Qnm2),Qnm1)
	!write(*,*) 'shape Q =',shape(mat_Qnm)
	!write(*,*) 'mat_Qnm',mat_Qnm

	end subroutine mat_Q
!
! --------------------------------------------------------------------------
!
! Coefs Q pour mat_Q pour M2M FFT
!
! --------------------------------------------------------------------------
	subroutine coef_Q(n,m,n1,theta1,Qnm)
	! calcul la matrice Qnm pour le calcul de la matrice des Qnm pour M2M et L2L
	! avec m positif
!	use parametres
	implicit none
	integer :: n,m,n1
	real(prec) :: theta1(n1)
	real(prec) :: Qnm(m:n,n1)

	real(prec) :: x(n1)
	integer :: i1

	x=cos(theta1)
	! calcul du premier terme
	!Qnm(m,:)=(-1.0_prec)**m * (1-x**2)**real(m/2.0_prec, prec) * sqrt(gamma(real(2*m+1, prec))) / (2.0_prec**m * gamma(real(m+1, prec)))
	Qnm(m,:)=(-1.0_prec)**m * (1.0_prec-x**2)**real(m/2.0_prec, prec) &
		* gamma(real(2*m+1, prec)) / (2.0_prec**m * gamma(real(m+1, prec)))

	do i1=m+1,n
	    if (i1==m+1) then
		!Qnm(i1,:)= sqrt(real(1+2*m, prec)) * x * Qnm(m,:)
		Qnm(i1,:)= (real(1+2*m, prec)) * x * Qnm(m,:)
	    else
!	        Qnm(i1,:)=1.0_prec / sqrt(real(i1**2-m**2, prec)) * ( (2.0_prec*i1-1.0_prec)*x*Qnm(i1-1,:) - sqrt(real((i1-1)**2-m**2, prec))*Qnm(i1-2,:) )
	        Qnm(i1,:)=1.0_prec / (real(i1-m, prec)) * ( (2.0_prec*i1-1.0_prec)*x*Qnm(i1-1,:) - (real(i1+m-1, prec))*Qnm(i1-2,:) )
	    endif
	enddo
	do i1=m,n
!                Qnm(i1,:)=sqrt(real(2*i1+1, prec)/(4.0_prec*pi)*gamma(real(i1-m+1, prec))/gamma(real(i1+m+1, prec))) * Qnm(i1,:)
				Qnm(i1,:)=sqrt(real(2*i1+1, prec)/(4.0_prec*pi) &
					*gamma(real(i1-m+1, prec))/gamma(real(i1+m+1, prec))) * Qnm(i1,:)
	end do
	!write(*,*) 'shape Q =',shape(Qnm),'m =',m,'n = ',n,'n1 =',n1
	!write(*,*) 'Qnm =',Qnm
	end subroutine coef_Q
! --------------------------------------------------------------------------------
!
! initialisation rotations pour version hybride BF
! précalcul des distances pour translations
! pour un ordre Nmax donné
! Variables non dépendants de la fréquence
! A faire dans le "main", avant la boucle en fréquence
! --------------------------------------------------------------------------------
	subroutine init_BF_MLFMM_v1(octree,NL,SPH,Nmax)
! initialisation variable pour MLFMM par FFT
	!use mod_tree_v2 ! pour type tree
	!use parametres  ! pour epsilon
	!use init_param_FMM
	implicit none
! variables entrée
	integer, intent(inout) :: NL
	integer, intent(in) :: Nmax
	type(tree), intent(inout) :: octree(NL)
! variables sortie
!	integer, intent(out) :: maxi_BF
	type(FMM0), intent(inout) :: SPH
! variable locales
	integer :: i1, num, numR, numT, i2, i3, ind
	integer :: n2p,n2m,m0,m11,m0p,m0m
	integer :: npos1,pos1(15)
	! (supprimé) real(prec) :: r, a,sigma,kr
!	real(prec), allocatable :: Pt(:,:),wp(:),theta0(:),PL(:,:),cos_phi(:,:),cos_phi0(:)
	! (supprimé) real(prec) :: test(45,45),plow,phigh
	! (supprimé) complex(prec), dimension(:), allocatable  :: M2L
	integer :: n2, res1, nx, ny, nz, m2
	real(prec) :: cote,R21,D21(3),rm2m1(3)
	real(prec), dimension(1:3) :: vt
    complex(prec), dimension(:),allocatable :: Rot
    integer,dimension(:),allocatable :: ix,iy

	! DEBUG : character(len=1) :: answer

!
!-----------------------------------------------------
! Initialise le nombre de termes dans la matrice de rotation
! Stockage sparse
!-----------------------------------------------------
!
    allocate(SPH%rot_num(Nmax))
    num=1
    do i1=1,Nmax
        num=num+(2*i1+1)**2
        SPH%rot_num(i1)=num
    end do
    if (verbose >= 2) write(log_unit,*) 'init_BF_MLFMM: num_init =', num

!
!--------------------------------------------------------------
! Initalise les 316 rotations, non fréquence dépendant pour M2L
! On initialise à de 0 à 342 (7^3=343) pour le stockage et la facilité de trouver la rotation en passant en base 7
! et on calcule 316 (343-27) intéractions au lieu de 290 normalement différentes, mais c'est plus simple que de faire des tests après pour les retrouver
!--------------------------------------------------------------
!
    allocate(SPH%rot_M2L(num,0:342))
    allocate(SPH%rot_M2L_i(num,0:342))
    allocate(SPH%rot_M2L_j(num,0:342))

	  cote=octree(1)%cell(1)%length
	  num=0
	  do i2=0,342  ! pour avoir toutes les intéractions
        ! passage en base 7 pour nxnynz
        nz=mod(i2,7)
        res1=(i2-nz)/7
        ny=mod(res1,7)
        nx=(res1-ny)/7
        ! on retranche 3 pour aller de -3 à 3
        nx=nx-3
        ny=ny-3
        nz=nz-3
        ! test interaction
        if (abs(nx) > 1 .or. abs(ny) > 1 .or. abs(nz) > 1) then   ! que les intéractions
            num=num+1
            ! calcul du vecteur de translation
            vt=real((/nx,ny,nz/),8)*cote
            R21=sqrt(dot_product(vt,vt)) ! norme
            D21=vt/R21 ! direction
            call ROT_BF(Nmax,D21,Rot,ix,iy,numR)
            SPH%rot_M2L(:,i2)=Rot
            SPH%rot_M2L_i(:,i2)=ix
            SPH%rot_M2L_j(:,i2)=iy
        end if
       end do
!
!-----------------------------------------------------
! Initalise les 8 rotations, non fréquence dépendant pour M2M et/ou L2L
!-----------------------------------------------------
!
    allocate(SPH%rot_M2M(SPH%rot_num(Nmax),8))
    allocate(SPH%rot_M2M_i(SPH%rot_num(Nmax),8))
    allocate(SPH%rot_M2M_j(SPH%rot_num(Nmax),8))
    do i1=0,1
        do i2=0,1
            do i3=0,1
                ind=i3+2*i2+4*i1+1   ! passage base 2 en base 10 : de 0 à 7 +1
                D21=(/i1,i2,i3/)
                D21=(D21-0.5_prec)*2.0_prec/sqrt(3.0_prec)
                call ROT_BF(Nmax,D21,Rot,ix,iy,numR)
                SPH%rot_M2M(:,ind)=Rot
                SPH%rot_M2M_i(:,ind)=ix
                SPH%rot_M2M_j(:,ind)=iy
            enddo
        enddo
    enddo
!    write(*,*) 'pass 8 rotations'
!
!-----------------------------------------------------
! Initalise la translation par niveau, non fréquence dépendant pour M2M et/ou L2L
! initialise juste la norme du vecteur |t| de translation
!-----------------------------------------------------
!
    allocate(SPH%tab_T_M2M(NL))
    do i1=NL,2,-1
       rm2m1=octree(i1-1)%cell(1)%center - octree(i1)%cell(octree(i1-1)%cell(1)%child(1))%center
       SPH%tab_T_M2M(i1)=sqrt(dot_product(rm2m1,rm2m1)) ! norme
    enddo
!    write(*,*) 'SPH0%tab_T_M2M=',SPH%tab_T_M2M(:)

!    write(*,*) 'pass ',NL,' translations'
!
!-----------------------------------------------------
! Initalise les 15 translations par niveaux, non fréquence dépendant pour M2L
! initialise juste la norme du vecteur |t| de translation
!-----------------------------------------------------
!
    allocate(SPH%tab_T_M2L(NL,15))
    sph%tab_T_M2L=0.0_prec

    do i1=1,NL
	  cote=octree(i1)%cell(1)%length
	  num=0
	  numT=0
      i2=0
        do while (numT < 15)  ! 15 maxi
        ! passage en base 7 pour nxnynz
        nz=mod(i2,7)
        res1=(i2-nz)/7
        ny=mod(res1,7)
        nx=(res1-ny)/7
        ! on retranche 3 pour aller de -3 à 3
        nx=nx-3
        ny=ny-3
        nz=nz-3
        ! test interaction
        if (abs(nx) > 1 .or. abs(ny) > 1 .or. abs(nz) > 1) then   ! que les intéractions
            num=num+1
            ! calcul du vecteur de translation
            vt=real((/nx,ny,nz/),8)*cote
            R21=sqrt(dot_product(vt,vt)) ! norme
            ! test norme
           	call FindInVector(15,abs(SPH%tab_T_M2L(i1,:)-R21) < 1e-10_prec,npos1,pos1)
            if (npos1==0) then
                numT=numT+1
                SPH%tab_T_M2L(i1,numT)=R21
            end if
        end if
        i2=i2+1
        end do
    end do

!
!---------------------------------------------------------------------
!
! précalcul des coefficients formalisme BF
! pour le calcul de la dérivée normale de la fonction de Green
! utilisé dans moments_FMM_p_T3_BF
!---------------------------------------------------------------------
!
!    if (allocated(SPH0%a)) deallocate(SPH0%a)
    allocate(SPH%a(2,(Nmax+1)**2))
!    if (allocated(b)) deallocate(b)
    allocate(SPH%b(4,(Nmax+1)**2))
    allocate(SPH%b_ex(4,(Nmax+1)**2))

    num=0
    do n2=0,Nmax
        n2p=n2+1
        n2m=n2-1
        do m2=-n2,n2
            num=num+1
            m0=-m2
            m0p=m0+1
            m0m=m0-1
            m11=-m0-1

            if (m11>=0) then
                SPH%b(1,num)=sqrt(real(n2p-m11-1, prec)*real(n2p-m11, prec)/real(2*n2p-1, prec)/real(2*n2p+1, prec))
            else
                SPH%b(1,num)=-sqrt(real(n2p-m11-1, prec)*real(n2p-m11, prec)/real(2*n2p-1, prec)/real(2*n2p+1, prec))
            endif
            if (-m2-1>=0) then
                SPH%b_ex(1,num)=sqrt(real(n2p-(-m2-1)-1, prec)*real(n2p-(-m2-1), prec)/real(2*n2p-1, prec)/real(2*n2p+1, prec))
            else
                SPH%b_ex(1,num)=-sqrt(real(n2p-(-m2-1)-1, prec)*real(n2p-(-m2-1), prec)/real(2*n2p-1, prec)/real(2*n2p+1, prec))
            endif

            if (abs(m0p)<=n2m) then
                if (m0>=0) then
                    SPH%b(2,num)=sqrt(real(n2-m0-1, prec)*real(n2-m0, prec)/real(2*n2-1, prec)/real(2*n2+1, prec))
                else
                    SPH%b(2,num)=-sqrt(real(n2-m0-1, prec)*real(n2-m0, prec)/real(2*n2-1, prec)/real(2*n2+1, prec))
                endif
            else
                SPH%b(2,num)=0.0_prec
            endif
            if (abs(m2+1)<=n2m) then
                if (m2>=0) then
                    SPH%b_ex(2,num)=sqrt(real(n2-m2-1, prec)*real(n2-m2, prec)/real(2*n2-1, prec)/real(2*n2+1, prec))
                else
                    SPH%b_ex(2,num)=-sqrt(real(n2-m2-1, prec)*real(n2-m2, prec)/real(2*n2-1, prec)/real(2*n2+1, prec))
                endif
            else
                SPH%b_ex(2,num)=0.0_prec
            endif

            if (m0m>=0) then
                SPH%b(3,num)=sqrt(real(n2p-m0m-1, prec)*real(n2p-m0m, prec)/real(2*n2p-1, prec)/real(2*n2p+1, prec))
            else
                SPH%b(3,num)=-sqrt(real(n2p-m0m-1, prec)*real(n2p-m0m, prec)/real(2*n2p-1, prec)/real(2*n2p+1, prec))
            endif
            if ((m2-1)>=0) then
                SPH%b_ex(3,num)=sqrt(real(n2p-(m2-1)-1, prec)*real(n2p-(m2-1), prec)/real(2*n2p-1, prec)/real(2*n2p+1, prec))
            else
                SPH%b_ex(3,num)=-sqrt(real(n2p-(m2-1)-1, prec)*real(n2p-(m2-1), prec)/real(2*n2p-1, prec)/real(2*n2p+1, prec))
            endif

            if (abs(m0m)<=n2m) then
                if ((-m0)>=0) then
                    SPH%b(4,num)=sqrt(real(n2+m0-1, prec)*real(n2+m0, prec)/real(2*n2-1, prec)/real(2*n2+1, prec))
                else
                    SPH%b(4,num)=-sqrt(real(n2+m0-1, prec)*real(n2+m0, prec)/real(2*n2-1, prec)/real(2*n2+1, prec))
                endif
            else
                SPH%b(4,num)=0.0_prec
            endif
            if (abs(m2-1)<=n2m) then
                if ((-m2)>=0) then
                    SPH%b_ex(4,num)=sqrt(real(n2+m2-1, prec)*real(n2+m2, prec)/real(2*n2-1, prec)/real(2*n2+1, prec))
                else
                    SPH%b_ex(4,num)=-sqrt(real(n2+m2-1, prec)*real(n2+m2, prec)/real(2*n2-1, prec)/real(2*n2+1, prec))
                endif
            else
                SPH%b_ex(4,num)=0.0_prec
            endif

            if (n2m>=abs(m0)) then
                SPH%a(1,num)=sqrt(real(n2m+1+m0, prec)*real(n2m+1-m0, prec)/real(2*n2m+1, prec)/real(2*n2m+3, prec))
            else
                SPH%a(1,num)=cmplx(0.0_prec,0.0_prec, kind=prec)
            endif

            SPH%a(2,num)=sqrt(real(n2+1+m0, prec)*real(n2+1-m0, prec)/real(2*n2+1, prec)/real(2*n2+3, prec))

        enddo
    enddo


	end subroutine init_BF_MLFMM_v1

! ---------------------------------------------------------------------------------------------------
!
! Version FFT (diagonale version), avec Burton/Miller
! pour variables constantes sur élément T0
! pour gradient dnp
!
! version déjà initialisée par init_FMM_T0
! la première fois calcul des moments et expansion pour initialiser l'intégration
! ---------------------------------------------------------------------------------------------------
	subroutine MLFMM_FFT_dnp(varin,octree,NL,N,SPH,alpha,k,varout)
	implicit none
! variables d'entrée
	integer, intent(in) :: NL
	integer, dimension (1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	complex(prec), intent(in) :: alpha
	complex(prec), dimension(NE) :: varin  ! version T0
	type(FMM),intent(inout) :: SPH(NL)
	type(tree),intent(inout) :: octree(NL)

! variables de sortie
	complex(prec), intent(out) :: varout(NE)    ! version T0
! variables locales
	integer :: i1
	! DEBUG : real(prec) :: t1,t2
	! (supprimé) complex(prec),allocatable :: temp1
	! DEBUG : character(len=1) :: answer

!   Calcul moment diagonaux
	if (.not. allocated(SPH(1)%V1)) call moments_FMM_HF(octree,NL,N,SPH,k)
	SPH(NL)%in(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
	do i1=1,SPH(1)%nnzv
		SPH(NL)%in(SPH(1)%iv(i1))=SPH(NL)%in(SPH(1)%iv(i1))+SPH(1)%V1(i1)*varin(SPH(1)%jv(i1))*ji*k/(4.0_prec*pi)**2    ! gradient
	end do

  	! pass upward
   	call pass_upward_FFT(NL,N,octree,SPH,k)
!    write(*,*) 'SPH(1)%in',SPH(1)%in(1:3)
!    call HF2BF_up_FMM(octree,NL,N,SPH,1)
!    write(*,*) 'SPH(1)%inBF',SPH(1)%inBF(1:3)
   	! passe downward_v2 : translations initialisées
   	call pass_downward_FFT(octree,NL,N,SPH,k)
!    call HF2BF_down_FMM(octree,NL,N,SPH,NL)
!	open(unit=10,file='dp_outBF_HFBF.txt')
!	do i1=1,(N(NL)+1)**2
!    write(10,*) SPH(NL)%outBF(i1)
!    enddo
!    close(10)
!	open(unit=10,file='dp_out_HF.txt')
!	do i1=1,(N(NL)+1)*(2*N(NL)+1)
!    write(10,*) SPH(NL)%out(i1)
!    enddo
!    close(10)
!    write(*,*) 'SPH(NL)%outBF',SPH(NL)%outBF(1:3)
!    write(*,*) 'SPH(NL)%out',SPH(NL)%out(1:3)


   	! expansion
	varout(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
	if (.not. allocated(SPH(1)%W1)) call expansions_BM_HF(octree,NL,N,SPH,k,alpha)

	do i1=1,SPH(1)%nnzw
		varout(SPH(1)%iw(i1))=varout(SPH(1)%iw(i1))+SPH(1)%W1(i1)*SPH(NL)%out(SPH(1)%jw(i1))
	end do


   	! prise en compte pour pression
!    varout=varout*ji*k/(4.0_prec*pi)**2
	end subroutine MLFMM_FFT_dnp
! ---------------------------------------------------------------------------------------------------
!
! Version Hybride (diagonale+BF), avec Burton/Miller
! pour variables constantes sur l'élément T0
! pour gradient de pression dnp
!
! version déjà initialisée par init_v5
! la première fois calcul des moments et expansion pour initialiser l'intégration
! ---------------------------------------------------------------------------------------------------
	subroutine MLFMM_HYB_dnp(varin,octree,NL,N,SPH,SPH0,alpha,k,varout)
	implicit none
! variables d'entrée
	integer, intent(in) :: NL
	integer, dimension (1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	complex(prec) :: alpha
	complex(prec), dimension(NE) :: varin   ! vraie version
	type(FMM),intent(inout) :: SPH(NL)
	type(FMM0),intent(inout) :: SPH0
	type(tree) :: octree(NL)

! variables de sortie
    complex(prec), intent(inout) :: varout(NE)   ! version T0
! variables locales
    integer :: i1
	! (supprimé) complex(prec), allocatable :: outBF(:)
	! DEBUG : real(prec) :: t1,t2
	! DEBUG : character(len=1) :: answer

	if (fmm_adaptatif) then                  ! octree adaptatif (module_BEM)
		call fmm_T0_adapt(1,varin,varin,octree,NL,N,SPH,SPH0,alpha,k,varout)
		return
	endif
   	if (SPH(NL)%BF==1) then
        if (.not. allocated(SPH(1)%V1)) call moments_FMM_BF(octree,NL,N,SPH,SPH0,k)
		SPH(NL)%inBF(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
		do i1=1,SPH(1)%nnzv
			SPH(NL)%inBF(SPH(1)%iv(i1))=SPH(NL)%inBF(SPH(1)%iv(i1))+SPH(1)%V1(i1)*varin(SPH(1)%jv(i1))*ji*k/(4.0_prec*pi)    ! gradient
		end do
    else
		if (.not. allocated(SPH(1)%V1)) call moments_FMM_HF(octree,NL,N,SPH,k)
		SPH(NL)%in(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
		do i1=1,SPH(1)%nnzv
			SPH(NL)%in(SPH(1)%iv(i1))=SPH(NL)%in(SPH(1)%iv(i1))+SPH(1)%V1(i1)*varin(SPH(1)%jv(i1))*ji*k/(4.0_prec*pi)**2    ! gradient
		end do
    endif
!    write(*,*) 'iv=',SPH(1)%iv(1:3)
!    write(*,*) 'jv=',SPH(1)%jv(1:3)
!    write(*,*) 'V1=',SPH(1)%V1(1:3)

!    write(*,*) 'SPH(NL)%inBF',SPH(NL)%inBF(1:10)
!    write(*,*) 'SPH(NL)%in',SPH(NL)%in(1:10)
   	! pass upward
   	call pass_upward_HYB(NL,N,octree,SPH,SPH0,k)
!  	call BF2HF_up_FMM(octree,NL,N,SPH,1)
!    write(*,*) 'SPH(1)%inBF',SPH(1)%inBF(1:3)
!    write(*,*) 'SPH(1)%in',SPH(1)%in(1:3)
   	! passe downward_v2 : translations initialisées
   	call pass_downward_HYB(octree,NL,N,SPH,SPH0,k)
!	call HF2BF_down_FMM(octree,NL,N,SPH,NL)
!   	if (SPH(NL)%BF==1) then
!		call BF2HF_down_FMM(octree,NL,N,SPH,NL)
!	endif
!	open(unit=10,file='dp_outBF_BF.txt')
!	do i1=1,(N(NL)+1)**2
!    write(10,*) SPH(NL)%outBF(i1)
!    enddo
!    close(10)
!	open(unit=10,file='dp_out_BFHF.txt')
!	do i1=1,(N(NL)+1)*(2*N(NL)+1)
!    write(10,*) SPH(NL)%out(i1)
!    enddo
!    close(10)
!    write(*,*) 'SPH(NL)%outBF',SPH(NL)%outBF(1:3)
!    write(*,*) 'SPH(NL)%out',SPH(NL)%out(1:3)
   	! expansion
   	varout(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
   	if (SPH(NL)%BF==1) then
		if (.not. allocated(SPH(1)%W1)) call expansions_BM_BF(octree,NL,N,SPH,SPH0,k,alpha)
!		write(*,*) 'nnzw=',SPH(1)%nnzw,' size(W1)=',size(SPH(1)%W1)
		do i1=1,SPH(1)%nnzw
			varout(SPH(1)%iw(i1))=varout(SPH(1)%iw(i1))+SPH(1)%W1(i1)*SPH(NL)%outBF(SPH(1)%jw(i1))
!			if (SPH(1)%W1(i1) .ne. SPH(1)%W1(i1)) then
!				write(*,*) 'i1=',i1
!				read(*,*) answer
!			endif
		end do
	else
		if (.not. allocated(SPH(1)%W1)) call expansions_BM_HF(octree,NL,N,SPH,k,alpha)
		do i1=1,SPH(1)%nnzw
			varout(SPH(1)%iw(i1))=varout(SPH(1)%iw(i1))+SPH(1)%W1(i1)*SPH(NL)%out(SPH(1)%jw(i1))
		end do
	endif

	end subroutine MLFMM_HYB_dnp
! ---------------------------------------------------------------------------------------------------
!
! Version Hybride (diagonale+BF), avec Burton/Miller
! pour variables constantes sur l'élément T0
! pour gradient de pression dnp
!
! Version pour propagation en FMM donc sur octree2 avec les Pext en plus dans l'octree
! alpha imposé à 0 pour la HIE (sans la contribution de alpha H dnp)
!
! version déjà initialisée par init_v5
! la première fois calcul des moments et expansion pour initialiser l'intégration
! ---------------------------------------------------------------------------------------------------
	subroutine MLFMM_HYB_propa_dnp(varin,octree,NL,N,SPH,SPH0,k,varout)
	implicit none
! variables d'entrée
	integer :: NL
	integer, dimension (1:NL) :: N
	real(prec) :: k
!	complex(prec) :: alpha
	complex(prec), dimension(NE+Next) :: varin   ! vraie version
	type(FMM) :: SPH(NL)
	type(FMM0) :: SPH0
	type(tree) :: octree(NL)

! variables de sortie
    complex(prec) :: varout(Next)   ! version T0
! variables locales
    integer :: i1
	! (supprimé) complex(prec), allocatable :: outBF(:)
	! DEBUG : real(prec) :: t1,t2
	! DEBUG : character(len=1) :: answer

! 	write(*,*) 'debut propa_dnp'
!   	write(*,*) 'nb_c dans 1387=',octree(4)%cell(1387)%Nb_c
  	if (SPH(NL)%BF==1) then
!        if (.not. allocated(SPH(1)%V1)) call moments_FMM_propa_BF(octree,NL,N,SPH,SPH0,k)
        call moments_FMM_propa_BF(octree,NL,N,SPH,SPH0,k)
		SPH(NL)%inBF(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
		do i1=1,SPH(1)%nnzv
			SPH(NL)%inBF(SPH(1)%iv(i1))=SPH(NL)%inBF(SPH(1)%iv(i1))+SPH(1)%V1(i1)*varin(SPH(1)%jv(i1))*ji*k/(4.0_prec*pi)    ! gradient
		end do
    else
		if (.not. allocated(SPH(1)%V1)) call moments_FMM_propa_HF(octree,NL,N,SPH,k)
		SPH(NL)%in(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
		do i1=1,SPH(1)%nnzv
			SPH(NL)%in(SPH(1)%iv(i1))=SPH(NL)%in(SPH(1)%iv(i1))+SPH(1)%V1(i1)*varin(SPH(1)%jv(i1))*ji*k/(4.0_prec*pi)**2    ! gradient
		end do
    endif
!    write(*,*) 'pass moments in propa_dnp'
!	write(*,*) 'nb_c dans 1387=',octree(4)%cell(1387)%Nb_c
!    write(*,*) 'shape(SPH(NL)%inBF) in propa_dnp',shape(SPH(NL)%inBF)
!    write(*,*) 'iv=',SPH(1)%iv(1:3)
!    write(*,*) 'jv=',SPH(1)%jv(1:3)
!    write(*,*) 'V1=',SPH(1)%V1(1:3)

!    write(*,*) 'SPH(NL)%inBF',SPH(NL)%inBF(1:10)
!    write(*,*) 'SPH(NL)%in',SPH(NL)%in(1:10)
   	! pass upward
   	call pass_upward_HYB(NL,N,octree,SPH,SPH0,k)
!   	write(*,*) 'pass upward_HYB dans propa_dnp'
!   	write(*,*) 'nb_c dans 1387=',octree(4)%cell(1387)%Nb_c
!  	call BF2HF_up_FMM(octree,NL,N,SPH,1)
!    write(*,*) 'SPH(1)%inBF',SPH(1)%inBF(1:3)
!    write(*,*) 'SPH(1)%in',SPH(1)%in(1:3)
   	! passe downward_v2 : translations initialisées
   	call pass_downward_HYB(octree,NL,N,SPH,SPH0,k)
!   	write(*,*) 'pass downward_HYB dans propa_dnp'
!   	write(*,*) 'nb_c dans 1387=',octree(4)%cell(1387)%Nb_c
!	call HF2BF_down_FMM(octree,NL,N,SPH,NL)
!   	if (SPH(NL)%BF==1) then
!		call BF2HF_down_FMM(octree,NL,N,SPH,NL)
!	endif
!	open(unit=10,file='dp_outBF_BF.txt')
!	do i1=1,(N(NL)+1)**2
!    write(10,*) SPH(NL)%outBF(i1)
!    enddo
!    close(10)
!	open(unit=10,file='dp_out_BFHF.txt')
!	do i1=1,(N(NL)+1)*(2*N(NL)+1)
!    write(10,*) SPH(NL)%out(i1)
!    enddo
!    close(10)
!    write(*,*) 'SPH(NL)%outBF',SPH(NL)%outBF(1:3)
!    write(*,*) 'SPH(NL)%out',SPH(NL)%out(1:3)
   	! expansion
   	varout(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
   	if (SPH(NL)%BF==1) then
		if (.not. allocated(SPH(1)%W1)) call expansions_BM_propa_BF(octree,NL,N,SPH,SPH0,k)
!		write(*,*) 'nnzw=',SPH(1)%nnzw,' size(W1)=',size(SPH(1)%W1)
!		write(*,*) 'après expansion'
!		write(*,*) 'nb_c dans 1387=',octree(4)%cell(1387)%Nb_c
!		write(*,*) 'nnzw=',SPH(1)%nnzw,' size(iw)=',size(SPH(1)%iw),' size(W1)=',size(SPH(1)%W1)
!		write(*,*) 'size(jw)=',size(SPH(1)%jw),' size(outBF)=',size(SPH(NL)%outBF)
		do i1=1,SPH(1)%nnzw
			varout(SPH(1)%iw(i1))=varout(SPH(1)%iw(i1))+SPH(1)%W1(i1)*SPH(NL)%outBF(SPH(1)%jw(i1))
!			if (SPH(1)%W1(i1) .ne. SPH(1)%W1(i1)) then
!				write(*,*) 'i1=',i1
!				read(*,*) answer
!			endif
		end do
!		write(*,*) 'après mul mat'
!	   	write(*,*) 'nb_c dans 1387=',octree(4)%cell(1387)%Nb_c
	else
		if (.not. allocated(SPH(1)%W1)) call expansions_BM_propa_HF(octree,NL,N,SPH,k)

		do i1=1,SPH(1)%nnzw
			varout(SPH(1)%iw(i1))=varout(SPH(1)%iw(i1))+SPH(1)%W1(i1)*SPH(NL)%out(SPH(1)%jw(i1))
		end do

	endif
!	write(*,*) 'pass expansions in propa_dnp'
!   	write(*,*) 'nb_c dans 1387=',octree(4)%cell(1387)%Nb_c

!	read(*,*) answer

	end subroutine MLFMM_HYB_propa_dnp
! ---------------------------------------------------------------------------------------------------
!
! Version FFT (diagonale version), avec Burton/Miller
! pour variables constantes sur élément T0
! pour pression p
!
! version déjà initialisée par init_FMM_T0
! la première fois calcul des moments pour initialiser l'intégration
! ---------------------------------------------------------------------------------------------------
	subroutine MLFMM_FFT_p(varin,octree,NL,N,SPH,alpha,k,varout)
	implicit none
! variables d'entrée
	integer, intent(in) :: NL
	integer, dimension (1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	complex(prec), intent(in) :: alpha
	complex(prec), dimension(NE) :: varin  ! version T0
	type(FMM),intent(inout) :: SPH(NL)
	type(tree) :: octree(NL)

! variables de sortie
	complex(prec), intent(inout) :: varout(NE)    ! version T0
! variables locales
	integer :: i1
	! DEBUG : real(prec) :: t1,t2
!	complex(prec) :: temp1(NE)
	! DEBUG : character(len=1) :: answer

!   Calcul moment diagonaux
	if (.not. allocated(SPH(1)%V1)) call moments_FMM_HF(octree,NL,N,SPH,k)
	SPH(NL)%in(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
	do i1=1,SPH(1)%nnzv
		SPH(NL)%in(SPH(1)%iv(i1))=SPH(NL)%in(SPH(1)%iv(i1))+SPH(1)%V2(i1)*varin(SPH(1)%jv(i1))*(k)**2/(4.0_prec*pi)**2 ! pression
	end do

  	! pass upward
   	call pass_upward_FFT(NL,N,octree,SPH,k)
!    write(*,*) 'SPH(1)%in',SPH(1)%in(1:3)
!    call HF2BF_up_FMM(octree,NL,N,SPH,1)
!    write(*,*) 'SPH(1)%inBF',SPH(1)%inBF(1:3)

   	! passe downward_v2 : translations initialisées
   	call pass_downward_FFT(octree,NL,N,SPH,k)
!   call HF2BF_down_FMM(octree,NL,N,SPH,NL)
!    write(*,*) 'SPH(NL)%outBF',SPH(NL)%outBF(1:3)
!    write(*,*) 'SPH(NL)%out',SPH(NL)%out(1:3)

   	! expansion
	varout(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
	if (.not. allocated(SPH(1)%W1)) call expansions_BM_HF(octree,NL,N,SPH,k,alpha)

	do i1=1,SPH(1)%nnzw
		varout(SPH(1)%iw(i1))=varout(SPH(1)%iw(i1))+SPH(1)%W1(i1)*SPH(NL)%out(SPH(1)%jw(i1))
	end do

   	! prise en compte pour pression
!    varout=varout*(k)**2/(4.0_prec*pi)**2
	end subroutine MLFMM_FFT_p
! ---------------------------------------------------------------------------------------------------
!
! Version Hybride (diagonale+BF), avec Burton/Miller
! pour variables constantes sur l'élément T0
! pour pression p
!
! version déjà initialisée par init_v5
! la première fois calcul des moments et expansion pour initialiser l'intégration
! ---------------------------------------------------------------------------------------------------
	subroutine MLFMM_HYB_p(varin,octree,NL,N,SPH,SPH0,alpha,k,varout)
	implicit none
! variables d'entrée
	integer :: NL
	integer, dimension (1:NL) :: N
	real(prec) :: k
	complex(prec) :: alpha
	complex(prec), dimension(NE) :: varin   ! vraie version
	type(FMM) :: SPH(NL)
	type(FMM0) :: SPH0
	type(tree) :: octree(NL)

! variables de sortie
    complex(prec), intent(inout) :: varout(NE)   ! version T0
! variables locales
    integer :: i1
	! (supprimé) complex(prec), allocatable :: outBF(:)
	! DEBUG : real(prec) :: t1,t2
	! DEBUG : character(len=1) :: answer

	if (fmm_adaptatif) then                  ! octree adaptatif (module_BEM)
		call fmm_T0_adapt(2,varin,varin,octree,NL,N,SPH,SPH0,alpha,k,varout)
		return
	endif
   	if (SPH(NL)%BF==1) then
        if (.not. allocated(SPH(1)%V1)) call moments_FMM_BF(octree,NL,N,SPH,SPH0,k)
		SPH(NL)%inBF(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
		do i1=1,SPH(1)%nnzv
			SPH(NL)%inBF(SPH(1)%iv(i1))=SPH(NL)%inBF(SPH(1)%iv(i1))+SPH(1)%V2(i1)*varin(SPH(1)%jv(i1))*(ji*k)/(4.0_prec*pi) ! pression
		end do
    else
		if (.not. allocated(SPH(1)%V1)) call moments_FMM_HF(octree,NL,N,SPH,k)
		SPH(NL)%in(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
		do i1=1,SPH(1)%nnzv
			SPH(NL)%in(SPH(1)%iv(i1))=SPH(NL)%in(SPH(1)%iv(i1))+SPH(1)%V2(i1)*varin(SPH(1)%jv(i1))*(k)**2/(4.0_prec*pi)**2 ! pression
		end do
    endif
   	! pass upward
   	call pass_upward_HYB(NL,N,octree,SPH,SPH0,k)
!   	call BF2HF_up_FMM(octree,NL,N,SPH,1)
!    write(*,*) 'SPH(1)%in',SPH(1)%in(1:3)
!    write(*,*) 'SPH(1)%inBF',SPH(1)%inBF(1:3)
!    write(*,*) 'SPH(1)%in',SPH(1)%in(1:3)
   	! passe downward_v2 : translations initialisées
   	call pass_downward_HYB(octree,NL,N,SPH,SPH0,k)
!	call BF2HF_down_FMM(octree,NL,N,SPH,NL)
!	open(unit=10,file='p_outBF_BF.txt')
!	do i1=1,(N(NL)+1)**2
!    write(10,*) SPH(NL)%outBF(i1)
!    enddo
!    close(10)
!	open(unit=10,file='p_out_BFHF.txt')
!	do i1=1,(N(NL)+1)*(2*N(NL)+1)
!    write(10,*) SPH(NL)%out(i1)
!    enddo
!    close(10)
!    write(*,*) 'SPH(NL)%outBF',SPH(NL)%outBF(1:3)
!    write(*,*) 'SPH(NL)%out',SPH(NL)%out(1:3)
   	! expansion
	varout(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
	if (SPH(NL)%BF==1) then
		if (.not. allocated(SPH(1)%W1)) call expansions_BM_BF(octree,NL,N,SPH,SPH0,k,alpha)
		do i1=1,SPH(1)%nnzw
			varout(SPH(1)%iw(i1))=varout(SPH(1)%iw(i1))+SPH(1)%W1(i1)*SPH(NL)%outBF(SPH(1)%jw(i1))
		end do
	else
		if (.not. allocated(SPH(1)%W1)) call expansions_BM_HF(octree,NL,N,SPH,k,alpha)
		do i1=1,SPH(1)%nnzw
			varout(SPH(1)%iw(i1))=varout(SPH(1)%iw(i1))+SPH(1)%W1(i1)*SPH(NL)%out(SPH(1)%jw(i1))
		end do
	endif
	end subroutine MLFMM_HYB_p
! ---------------------------------------------------------------------------------------------------
!
! Version Hybride (diagonale+BF), avec Burton/Miller
! pour variables constantes sur l'élément T0
! pour pression p
!
! Version pour propagation en FMM donc sur octree2 avec les Pext en plus dans l'octree
! alpha imposé à 0 pour la HIE (sans la contribution de dnP)
!
! version déjà initialisée par init_v5
! la première fois calcul des moments et expansion pour initialiser l'intégration
! ---------------------------------------------------------------------------------------------------
	subroutine MLFMM_HYB_propa_p(varin,octree,NL,N,SPH,SPH0,k,varout)
	implicit none
! variables d'entrée
	integer :: NL
	integer, dimension (1:NL) :: N
	real(prec) :: k
!	complex(prec), intent(in) :: alpha
	complex(prec), dimension(NE+Next) :: varin   ! vraie version
	type(FMM) :: SPH(NL)
	type(FMM0) :: SPH0
	type(tree) :: octree(NL)

! variables de sortie
    complex(prec) :: varout(Next)   ! version T0
! variables locales
    integer :: i1
	! (supprimé) complex(prec), allocatable :: outBF(:)
	! DEBUG : real(prec) :: t1,t2
	! DEBUG : character(len=1) :: answer

!	write(*,*) 'début propa_p'
!   	write(*,*) 'nb_c dans 1387=',octree(4)%cell(1387)%Nb_c
   	if (SPH(NL)%BF==1) then
        if (.not. allocated(SPH(1)%V1)) call moments_FMM_propa_BF(octree,NL,N,SPH,SPH0,k)
        !call moments_FMM_propa_BF(octree,NL,N,SPH,SPH0,k)
		SPH(NL)%inBF(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
		do i1=1,SPH(1)%nnzv
			SPH(NL)%inBF(SPH(1)%iv(i1))=SPH(NL)%inBF(SPH(1)%iv(i1))+SPH(1)%V2(i1)*varin(SPH(1)%jv(i1))*(ji*k)/(4.0_prec*pi) ! pression
		end do
    else
		if (.not. allocated(SPH(1)%V1)) call moments_FMM_propa_HF(octree,NL,N,SPH,k)
		SPH(NL)%in(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
		do i1=1,SPH(1)%nnzv
			SPH(NL)%in(SPH(1)%iv(i1))=SPH(NL)%in(SPH(1)%iv(i1))+SPH(1)%V2(i1)*varin(SPH(1)%jv(i1))*(k)**2/(4.0_prec*pi)**2 ! pression
		end do
    endif
!    write(*,*) 'pass moments in propa_p'

    !write(*,*) 'shape(SPH(NL)%inBF) in propa_p',shape(SPH(NL)%inBF)
   	! pass upward
   	call pass_upward_HYB(NL,N,octree,SPH,SPH0,k)
!   	write(*,*) 'pass upward_HYB dans propa_p'
!   	call BF2HF_up_FMM(octree,NL,N,SPH,1)
!    write(*,*) 'SPH(1)%in',SPH(1)%in(1:3)
!    write(*,*) 'SPH(1)%inBF',SPH(1)%inBF(1:3)
!    write(*,*) 'SPH(1)%in',SPH(1)%in(1:3)
   	! passe downward_v2 : translations initialisées
   	call pass_downward_HYB(octree,NL,N,SPH,SPH0,k)
!	write(*,*) 'pass downward_HYB dans propa_p'
!	call BF2HF_down_FMM(octree,NL,N,SPH,NL)
!	open(unit=10,file='p_outBF_BF.txt')
!	do i1=1,(N(NL)+1)**2
!    write(10,*) SPH(NL)%outBF(i1)
!    enddo
!    close(10)
!	open(unit=10,file='p_out_BFHF.txt')
!	do i1=1,(N(NL)+1)*(2*N(NL)+1)
!    write(10,*) SPH(NL)%out(i1)
!    enddo
!    close(10)
!    write(*,*) 'SPH(NL)%outBF',SPH(NL)%outBF(1:3)
!    write(*,*) 'SPH(NL)%out',SPH(NL)%out(1:3)
   	! expansion
	varout(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
	if (SPH(NL)%BF==1) then
		if (.not. allocated(SPH(1)%W1)) call expansions_BM_propa_BF(octree,NL,N,SPH,SPH0,k)
		do i1=1,SPH(1)%nnzw
			varout(SPH(1)%iw(i1))=varout(SPH(1)%iw(i1))+SPH(1)%W1(i1)*SPH(NL)%outBF(SPH(1)%jw(i1))
		end do
	else
		if (.not. allocated(SPH(1)%W1)) call expansions_BM_propa_HF(octree,NL,N,SPH,k)
		do i1=1,SPH(1)%nnzw
			varout(SPH(1)%iw(i1))=varout(SPH(1)%iw(i1))+SPH(1)%W1(i1)*SPH(NL)%out(SPH(1)%jw(i1))
		end do
	endif
	end subroutine MLFMM_HYB_propa_p
! ---------------------------------------------------------------------------------------------------
!
! Version FFT (diagonale version), avec Burton/Miller
! pour variables constantes sur élément T0
! pour combinaison pression p et gradient dnp
!
! equivalent BEM de (D1+alphaK1) p - (M1+alphaH1) dnp : donc dnp avec signe '+' car le '-' est dans la formulation
!
! version déjà initialisée par init_FMM_T0
! la première fois calcul des moments pour initialiser l'intégration
! ---------------------------------------------------------------------------------------------------
	subroutine MLFMM_FFT_p_dnp(varin1,varin2,octree,NL,N,SPH,alpha,k,varout)
	implicit none
! variables d'entrée
	integer, intent(in) :: NL
	integer, dimension (1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	complex(prec), intent(in) :: alpha
	complex(prec), dimension(NE) :: varin1,varin2  ! version T0
	type(FMM),intent(inout) :: SPH(NL)
	type(tree) :: octree(NL)

! variables de sortie
	complex(prec), intent(inout) :: varout(NE)    ! version T0
! variables locales
	integer :: i1
	! DEBUG : real(prec) :: t1,t2
!	complex(prec) :: temp1(NE)
	! DEBUG : character(len=1) :: answer

!   Calcul moment diagonaux
	if (.not. allocated(SPH(1)%V1)) call moments_FMM_HF(octree,NL,N,SPH,k)
	SPH(NL)%in(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
	do i1=1,SPH(1)%nnzv
		SPH(NL)%in(SPH(1)%iv(i1))=SPH(NL)%in(SPH(1)%iv(i1))+SPH(1)%V2(i1)*varin1(SPH(1)%jv(i1))*(k)**2/(4.0_prec*pi)**2 - & ! pression
			ji*k/(4.0_prec*pi)**2 * SPH(1)%V1(i1)*varin2(SPH(1)%jv(i1))  ! gradient
	end do

  	! pass upward
   	call pass_upward_FFT(NL,N,octree,SPH,k)

   	! passe downward_v2 : translations initialisées
   	call pass_downward_FFT(octree,NL,N,SPH,k)

   	! expansion
	varout(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
	if (.not. allocated(SPH(1)%W1)) call expansions_BM_HF(octree,NL,N,SPH,k,alpha)

	do i1=1,SPH(1)%nnzw
		varout(SPH(1)%iw(i1))=varout(SPH(1)%iw(i1))+SPH(1)%W1(i1)*SPH(NL)%out(SPH(1)%jw(i1))
	end do

   	! prise en compte pour pression
!    varout=varout*(k)**2/(4.0_prec*pi)**2
	end subroutine MLFMM_FFT_p_dnp
! ---------------------------------------------------------------------------------------------------
!
! Version Hybride (diagonale+BF), avec Burton/Miller
! pour variables constantes sur l'élément T0
! pour combinaison de pression p et son gradient dnp
!
! equivalent BEM de (D1+alphaK1) p - (M1+alphaH1) dnp : donc dnp avec signe '+' car le '-' est dans la formulation
!
! version déjà initialisée par init_v5
! la première fois calcul des moments et expansion pour initialiser l'intégration
! ---------------------------------------------------------------------------------------------------
	subroutine MLFMM_HYB_p_dnp(varin1,varin2,octree,NL,N,SPH,SPH0,alpha,k,varout)
	implicit none
! variables d'entrée
	integer, intent(in) :: NL
	integer, dimension (1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	complex(prec), intent(in) :: alpha
	complex(prec), dimension(NE) :: varin1,varin2   ! vraie version
	type(FMM),intent(inout) :: SPH(NL)
	type(FMM0),intent(inout) :: SPH0
	type(tree) :: octree(NL)

! variables de sortie
    complex(prec), intent(inout) :: varout(NE)   ! version T0
! variables locales
    integer :: i1
	! (supprimé) complex(prec), allocatable :: outBF(:)
	! DEBUG : real(prec) :: t1,t2
	! DEBUG : character(len=1) :: answer

	if (fmm_adaptatif) then                  ! octree adaptatif (module_BEM)
		call fmm_T0_adapt(4,varin1,varin2,octree,NL,N,SPH,SPH0,alpha,k,varout)
		return
	endif
   	if (SPH(NL)%BF==1) then
        if (.not. allocated(SPH(1)%V1)) call moments_FMM_BF(octree,NL,N,SPH,SPH0,k)
		SPH(NL)%inBF(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
		do i1=1,SPH(1)%nnzv
			SPH(NL)%inBF(SPH(1)%iv(i1))=SPH(NL)%inBF(SPH(1)%iv(i1))+SPH(1)%V2(i1)*varin1(SPH(1)%jv(i1))*(ji*k)/(4.0_prec*pi) - & ! pression
			ji*k/(4.0_prec*pi) * SPH(1)%V1(i1)*varin2(SPH(1)%jv(i1))  ! gradient
		end do
    else
		if (.not. allocated(SPH(1)%V1)) call moments_FMM_HF(octree,NL,N,SPH,k)
		SPH(NL)%in(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
		do i1=1,SPH(1)%nnzv
			SPH(NL)%in(SPH(1)%iv(i1))=SPH(NL)%in(SPH(1)%iv(i1))+SPH(1)%V2(i1)*varin1(SPH(1)%jv(i1))*(k)**2/(4.0_prec*pi)**2  - & ! pression
			ji*k/(4.0_prec*pi)**2 * SPH(1)%V1(i1)*varin2(SPH(1)%jv(i1))  ! gradient
		end do
    endif
   	! pass upward
   	call pass_upward_HYB(NL,N,octree,SPH,SPH0,k)
   	! passe downward_v2 : translations initialisées
   	call pass_downward_HYB(octree,NL,N,SPH,SPH0,k)
!	call BF2HF_down_FMM(octree,NL,N,SPH,NL)
   	! expansion
	varout(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
   	if (SPH(NL)%BF==1) then
		if (.not. allocated(SPH(1)%W1)) call expansions_BM_BF(octree,NL,N,SPH,SPH0,k,alpha)
		do i1=1,SPH(1)%nnzw
			varout(SPH(1)%iw(i1))=varout(SPH(1)%iw(i1))+SPH(1)%W1(i1)*SPH(NL)%outBF(SPH(1)%jw(i1))
		end do
	else
		if (.not. allocated(SPH(1)%W1)) call expansions_BM_HF(octree,NL,N,SPH,k,alpha)
		do i1=1,SPH(1)%nnzw
			varout(SPH(1)%iw(i1))=varout(SPH(1)%iw(i1))+SPH(1)%W1(i1)*SPH(NL)%out(SPH(1)%jw(i1))
		end do
	endif
	end subroutine MLFMM_HYB_p_dnp
!
! --------------------------------------------------------------------------------------------------
!
! Pass_upward_v1 : version FFT
!
! --------------------------------------------------------------------------------------------------
    subroutine pass_upward_FFT(NL,N,octree,SPH,k)
    implicit none
	integer, intent(in) :: NL
	integer, dimension (1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	type(FMM),intent(inout) :: SPH(NL)
	type(tree) :: octree(NL)

	integer :: i1

	do i1=NL,2,-1
         call upward_HF(NL,N,SPH,octree,k,i1)
	enddo
    end subroutine pass_upward_FFT
!
! --------------------------------------------------------------------------------------------------
!
! Pass_upward_BF : version BF
!
! --------------------------------------------------------------------------------------------------
    subroutine pass_upward_HYB(NL,N,octree,SPH,SPH0,k)
    implicit none
	integer :: NL
	integer, dimension (1:NL) :: N
	real(prec) :: k
	type(FMM) :: SPH(NL)
	type(FMM0) :: SPH0
	type(tree) :: octree(NL)
	! DEBUG : real(prec) :: t1,t2
	integer :: i1

	do i1=NL,2,-1

        if (SPH(i1)%BF==1) then  ! niveau BF
!			write(*,*) 'in upward_BF'
!			write(*,*) 'nb_c dans 1387=',octree(4)%cell(1387)%Nb_c

			call upward_BF(NL,N,SPH,SPH0,octree,i1)
!            write(*,*) 'pass upward_BF'
        else
            call upward_HF(NL,N,SPH,octree,k,i1)
        endif
	enddo
    end subroutine pass_upward_HYB
! --------------------------------------------------------------------------------------------------
!
! UPWARD version diagonale
!
! --------------------------------------------------------------------------------------------------
    subroutine upward_HF(NL,N,SPH,octree,k,i1)
    implicit none
! variables d'entrée
    integer, intent(in) :: NL
    integer, dimension (1:NL), intent(in) :: N
	type(FMM) :: SPH(NL)
	type(tree) :: octree(NL)
	real(prec) :: k
	integer :: i1
! variables de sortie
! variables locales
    integer :: npu_gro, npu_fin, i2,i3
	real(prec), dimension(:), allocatable :: prodsca
	complex(prec), dimension(:), allocatable :: Mq,M_gro
    real(prec), dimension(1:3) :: rm2m1
    integer, dimension(2) :: taille


       ! write(*,*) 'i1=',i1
    npu_fin=(N(i1)+1)*(2*N(i1)+1) ! nbre de points sphère unitaire niveau i1
    npu_gro=(N(i1-1)+1)*(2*N(i1-1)+1) ! nbre de points sphère unitaire niveau i1-1
!    write(*,*) 'in upward, npu_gro(',i1-1,')=',npu_gro
    if (allocated(SPH(i1-1)%in)) deallocate(SPH(i1-1)%in)
    allocate(SPH(i1-1)%in(1:octree(i1-1)%NL*npu_gro))

    ! INTERPOLER PUIS DEPHASER (2026-09-29) : le produit multipole x exp(ik s.d)
    ! a une bande ~ N(i1) + k|d| que la grille de l'ENFANT ne resout pas ; on
    ! interpole donc chaque enfant (exact, bande N(i1)) vers la grille du PARENT,
    ! puis on dephase et on somme sur cette grille. Voir interp_puis_dephase.
    if ((interp_puis_dephase .or. interp_solveur == 1) .and. N(i1-1) /= N(i1)) then
        allocate(Mq(npu_fin), M_gro(npu_gro), prodsca(npu_gro))
        do i2=1,octree(i1-1)%NL
            SPH(i1-1)%in((i2-1)*npu_gro+1:i2*npu_gro)=cmplx(0.0_prec,0.0_prec, kind=prec)
            do i3=1,octree(i1-1)%cell(i2)%Nb_c
                if (.not. octree(i1)%cell(octree(i1-1)%cell(i2)%child(i3))%a_source) cycle   ! adaptatif : multipole nul
                rm2m1=octree(i1-1)%cell(i2)%center - octree(i1)%cell(octree(i1-1)%cell(i2)%child(i3))%center
                Mq=SPH(i1)%in(((octree(i1-1)%cell(i2)%child(i3))-1)*npu_fin+1:octree(i1-1)%cell(i2)%child(i3)*npu_fin)
                call M2M_FFT_v1(Mq,SPH,NL,N,i1-1,M_gro)          ! grille enfant -> grille parent
                prodsca=matmul(SPH(i1-1)%Ps,rm2m1)
                SPH(i1-1)%in((i2-1)*npu_gro+1:i2*npu_gro)=SPH(i1-1)%in((i2-1)*npu_gro+1:i2*npu_gro) + &
                    M_gro*exp(ji*k*prodsca)
            enddo
        enddo
        return
    endif

    do i2=1,octree(i1-1)%NL     ! boucle sur cellules level -1

        if (allocated(Mq)) deallocate(Mq)
            allocate(Mq(npu_fin))
            Mq=cmplx(0.0_prec,0.0_prec, kind=prec)
            do i3=1,octree(i1-1)%cell(i2)%Nb_c  ! boucle sur les childs
            	if (.not. octree(i1)%cell(octree(i1-1)%cell(i2)%child(i3))%a_source) cycle   ! adaptatif : multipole nul
            ! calcul produit scalaire
                rm2m1=octree(i1-1)%cell(i2)%center - octree(i1)%cell(octree(i1-1)%cell(i2)%child(i3))%center
                taille=shape(SPH(i1)%Ps)
                if (allocated(prodsca)) deallocate(prodsca)
                allocate(prodsca(taille(1)))

                prodsca=matmul(SPH(i1)%Ps,rm2m1)
                Mq=Mq+ &
        SPH(i1)%in(((octree(i1-1)%cell(i2)%child(i3))-1)*npu_fin+1:octree(i1-1)%cell(i2)%child(i3)*npu_fin)*exp(ji*k*prodsca)
            enddo
            if (N(i1-1)==N(i1)) then
                SPH(i1-1)%in((i2-1)*npu_gro+1:i2*npu_gro)=Mq
            else
                if (allocated(M_gro)) deallocate(M_gro)
                    allocate(M_gro(npu_gro))
                    call M2M_FFT_v1(Mq,SPH,NL,N,i1-1,M_gro)
                    SPH(i1-1)%in((i2-1)*npu_gro+1:i2*npu_gro)=M_gro
            endif
    enddo

    end subroutine upward_HF
! ------------------------------------------------------------------------------------------------
!
! UPWARD version BF, version pré-calculée
!
! ------------------------------------------------------------------------------------------------
    subroutine upward_BF(NL,N,SPH,SPH0,octree,i1)
    implicit none
! variables d'entrée
    integer :: NL
    integer, dimension (1:NL) :: N
	type(FMM) :: SPH(NL)
	type(FMM0) :: SPH0
	type(tree) :: octree(NL)
	integer :: i1
! variables de sortie
! variables locales
    integer :: nanm_fin, nanm_gro, i2, i3, i4, Ni1, Ni2, ind, Nx, Ny, Nz
	complex(prec), dimension(:), allocatable :: Mq_fin,Mq_gro,Mq_gro2
	real(prec), dimension(1:3) :: rm2m1,D21
	real(prec) :: R21
	! (supprimé) complex(prec), dimension(:), allocatable :: Rot
    ! (supprimé) integer, dimension(:), allocatable :: ix,iy
	! (supprimé) real(prec), dimension(:), allocatable :: M2M
    ! (supprimé) integer, dimension(:), allocatable :: M2M_i,M2M_j

    Ni1=N(i1-1)
    Ni2=N(i1)
	nanm_fin=(Ni2+1)**2              ! nombre de coefs HS correspondant au niveau i1
	nanm_gro=(Ni1+1)**2              ! nombre de coefs HS correspondant au niveau i1-1
	if (allocated(SPH(i1-1)%inBF)) deallocate(SPH(i1-1)%inBF)
	allocate(SPH(i1-1)%inBF(nanm_gro*octree(i1-1)%NL))
	SPH(i1-1)%inBF(:)=cmplx(0.0_prec,0.0_prec, kind=prec)

!	write(*,*) 'Cellules par niveaux maillage total ',octree(:)%NL
!	write(*,*) 'nb_c dans 1387=',octree(4)%cell(1387)%Nb_c
	do i2=1,octree(i1-1)%NL     ! boucle sur cellules level -1
		if (allocated(Mq_gro)) deallocate(Mq_gro)
		allocate(Mq_gro(nanm_gro))
		if (allocated(Mq_gro2)) deallocate(Mq_gro2)
		allocate(Mq_gro2(nanm_gro))
		if (allocated(Mq_fin)) deallocate(Mq_fin)
		allocate(Mq_fin(nanm_fin))
		do i3=1,octree(i1-1)%cell(i2)%Nb_c  ! boucle sur les childs
			if (.not. octree(i1)%cell(octree(i1-1)%cell(i2)%child(i3))%a_source) cycle   ! adaptatif : multipole nul
		! calcul produit scalaire
			!write(*,*) 'center=',octree(i1-1)%cell(i2)%center
!			if (i2==1387) then
!				write(*,*) 'i1=',i1,' i2=',i2,' i3=',i3
!				write(*,*) 'child_center=',octree(i1)%cell(octree(i1-1)%cell(i2)%child(i3))%center
!			endif
			rm2m1=octree(i1-1)%cell(i2)%center - octree(i1)%cell(octree(i1-1)%cell(i2)%child(i3))%center
			R21=sqrt(dot_product(rm2m1,rm2m1)) ! norme
			D21=rm2m1/R21 ! direction
			! Bidouille Rotation déjà calculée, passage base 2
			Nx=nint(D21(1)*sqrt(3.0_prec)/2.0_prec+0.5_prec)
			Ny=nint(D21(2)*sqrt(3.0_prec)/2.0_prec+0.5_prec)
			Nz=nint(D21(3)*sqrt(3.0_prec)/2.0_prec+0.5_prec)
			ind=Nz+2*Ny+4*Nx+1   ! passage base 2 en base 10 : de 0 à 7 +1
			!write(*,*) 'ind_rot=',ind
			! translation M2M coaxiale
			! rotation matrice creuse niveau fin
			Mq_fin=cmplx(0.0_prec,0.0_prec, kind=prec)
!			write(*,*) 'rot(Ni2)=',SPH0%rot_num(Ni2)
			do i4=1,SPH0%rot_num(Ni2)
					Mq_fin(SPH0%rot_M2M_i(i4,ind))=Mq_fin(SPH0%rot_M2M_i(i4,ind))+ &
					SPH(i1)%inBF(((octree(i1-1)%cell(i2)%child(i3))-1)*nanm_fin+SPH0%rot_M2M_j(i4,ind))*SPH0%rot_M2M(i4,ind)
			enddo
			! translation matrice creuse
			Mq_gro2=cmplx(0.0_prec,0.0_prec, kind=prec)
!			write(*,*) 'trans_M2M=',size(SPH(i1)%trans_M2M)
			do i4=1,size(SPH(i1)%trans_M2M)
					Mq_gro2(SPH(i1)%trans_M2M_i(i4))=Mq_gro2(SPH(i1)%trans_M2M_i(i4))+ &
					Mq_fin(SPH(i1)%trans_M2M_j(i4))*SPH(i1)%trans_M2M(i4)
			enddo
			! rotation matrice creuse transposée conjuguée
			Mq_gro=cmplx(0.0_prec,0.0_prec, kind=prec)
!			write(*,*) 'rot(Ni1)=',SPH0%rot_num(Ni1)
			do i4=1,SPH0%rot_num(Ni1)
					Mq_gro(SPH0%rot_M2M_j(i4,ind))=Mq_gro(SPH0%rot_M2M_j(i4,ind))+ &
					Mq_gro2(SPH0%rot_M2M_i(i4,ind))*conjg(SPH0%rot_M2M(i4,ind))
			enddo
			SPH(i1-1)%inBF((i2-1)*nanm_gro+1:i2*nanm_gro)=SPH(i1-1)%inBF((i2-1)*nanm_gro+1:i2*nanm_gro) + Mq_gro
		enddo
	enddo
!	write(*,*) 'pass in upward_BF before BF2HF'
	if (SPH(i1-1)%BF==0) then
		call BF2HF_up_FMM(octree,NL,N,SPH,i1-1)
!		write(*,*) 'pass_BF2HF_up_FMM'
	endif

	end subroutine upward_BF
!
! ----------------------------------------------------------------------------------------
!
! Passage moment BF -> HF, upward, toutes les cellules
! SPH(Lev)%inBF - > SPH(Lev)%in
!
! ----------------------------------------------------------------------------------------
	subroutine BF2HF_up_FMM(octree,NL,N,SPH,Lev)
	implicit none
! variables d'entrée
	integer :: NL,Lev
	integer, dimension(1:NL) :: N
	type(tree) :: octree(NL)

! variables de sortie
	type(FMM) :: SPH(NL)
! variables locales
	integer :: i1,n2,m2,ind
	integer :: npu,nanm
	real(prec), allocatable, dimension(:,:) :: Psph
	complex(prec),allocatable, dimension(:,:) :: Ynm
	! DEBUG : character(len=1) :: answer

	npu=(N(Lev)+1)*(2*N(Lev)+1)  ! nbre de points sur sphère unitaire
	nanm=(N(Lev)+1)**2
	! nbre de coefs correspondants
	allocate(Psph(npu,3))
	call car2sph(SPH(Lev)%Ps,Psph,npu)
!	write(*,*) 'shape(SPH(Lev)%Ps)=',shape(SPH(Lev)%Ps)
	if (allocated(SPH(Lev)%in)) deallocate(SPH(Lev)%in)
	allocate(SPH(Lev)%in(octree(Lev)%NL*npu))
	allocate(Ynm(nanm,npu))
!	write(*,*) 'nanm= ',nanm,' npu= ',npu,' N(Lev)=',N(Lev),' Lev=',Lev,' shape(Ynm)=',shape(Ynm)

	call Ynm_gum(N(Lev),Psph(:,2),Psph(:,3),npu,Ynm)
	SPH(Lev)%in=cmplx(0.0_prec,0.0_prec, kind=prec)
	do n2=0,N(Lev)
		do m2=-n2,n2
			ind=n2**2+n2+1+m2
			do i1=1,octree(Lev)%NL 	! boucle cellules level
				SPH(Lev)%in((i1-1)*npu+1:i1*npu)=SPH(Lev)%in((i1-1)*npu+1:i1*npu)+&
            ji**(-n2)* SPH(Lev)%inBF((i1-1)*nanm+ind) * Ynm(ind,:)
!            ji**(-n2)* SPH(Lev)%outBF((i1-1)*nanm+ind) * Ynm /(ji*k)
			enddo
		enddo
	enddo

	end subroutine BF2HF_up_FMM
!
! ----------------------------------------------------------------------------------------
!
! Passage moment BF -> HF, downward pour test, toutes les cellules
! SPH(Lev)%outBF - > SPH(Lev)%out
!
! ----------------------------------------------------------------------------------------
	subroutine BF2HF_down_FMM(octree,NL,N,SPH,Lev)
	implicit none
! variables d'entrée
	integer, intent(in) :: NL,Lev
	integer, dimension(1:NL), intent(in) :: N
	type(tree) :: octree(NL)

! variables de sortie
	type(FMM),intent(inout) :: SPH(NL)
! variables locales
	integer :: i1,n2,m2,ind
	integer :: npu,nanm
	real(prec), allocatable, dimension(:,:) :: Psph
	complex(prec),allocatable, dimension(:,:) :: Ynm
	! DEBUG : character(len=1) :: answer

	npu=(N(Lev)+1)*(2*N(Lev)+1)  ! nbre de points sur sphère unitaire
	nanm=(N(Lev)+1)**2
	! nbre de coefs correspondants
	allocate(Psph(npu,3))
	call car2sph(SPH(Lev)%Ps,Psph,npu)
	if (allocated(SPH(Lev)%out)) deallocate(SPH(Lev)%out)
	allocate(SPH(Lev)%out(octree(Lev)%NL*npu))
	allocate(Ynm(nanm,npu))
	call Ynm_gum(N(Lev),Psph(:,2),Psph(:,3),npu,Ynm)
	SPH(Lev)%out=cmplx(0.0_prec,0.0_prec, kind=prec)
	do n2=0,N(Lev)
		do m2=-n2,n2
			ind=n2**2+n2+1+m2
			do i1=1,octree(Lev)%NL 	! boucle cellules level
				SPH(Lev)%out((i1-1)*npu+1:i1*npu)=SPH(Lev)%out((i1-1)*npu+1:i1*npu)+&
            ji**(-n2)* SPH(Lev)%outBF((i1-1)*nanm+ind) * Ynm(ind,:)
!            ji**(-n2)* SPH(Lev)%outBF((i1-1)*nanm+ind) * Ynm /(ji*k)
			enddo
		enddo
	enddo

	end subroutine BF2HF_down_FMM
! --------------------------------------------------------------------------------------------------
!
! pass_downward_FFT : version FFT
!
! --------------------------------------------------------------------------------------------------
    subroutine pass_downward_FFT(octree,NL,N,SPH,k)
    implicit none
    ! variables d'entrée
	integer, intent(in) :: NL
	integer, dimension (1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	type(FMM),intent(inout) :: SPH(NL)
	type(tree) :: octree(NL)

	! variables locales
!    complex(prec), allocatable :: L(:),L2(:)
    integer :: i1


!    SPH(1)%tab_R_M2L=0.0_prec

    ! M2L suivi de L2L jusqu'à niveau NL-1, M2L seule au dernier niveau NL
	do i1=1,NL    ! boucle en niveau
            call down_HF(NL,N,SPH,octree,k,i1)
    enddo  ! fin boucle en niveau

    end subroutine pass_downward_FFT
! --------------------------------------------------------------------------------------------------
!
! pass_downward_BF : version BF
!
! --------------------------------------------------------------------------------------------------
    subroutine pass_downward_HYB(octree,NL,N,SPH,SPH0,k)
    implicit none
    ! variables d'entrée
	integer, intent(in) :: NL
	integer, dimension (1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	type(FMM),intent(inout) :: SPH(NL)
	type(FMM0),intent(inout) :: SPH0
	type(tree) :: octree(NL)

	! variables locales
    ! (supprimé) complex(prec), allocatable :: L(:),L2(:)
    integer :: i1
    ! DEBUG : real(prec) :: t1,t2

    SPH(1)%tab_R_M2L=0.0_prec

    ! M2L suivi de L2L jusqu'à niveau NL-1 et M2L en NL
	do i1=1,NL    ! boucle en niveau
!        t1=omp_get_wtime()
        if (i1 < NL) then
            if (SPH(i1)%BF==1 .and. SPH(i1+1)%BF==1) then
                call down_BF(NL,N,SPH,SPH0,octree,i1)
            elseif (SPH(i1)%BF==0 .and. SPH(i1+1)%BF==0) then
                call down_HF(NL,N,SPH,octree,k,i1)
            else
                call down_HYB(NL,N,SPH,SPH0,octree,i1)
            endif
        else
            if (SPH(i1)%BF==1) then
                call down_BF(NL,N,SPH,SPH0,octree,i1)
            else
                call down_HF(NL,N,SPH,octree,k,i1)
            endif
        endif
!        t2=omp_get_wtime()
!        write(*,*) 'downward niveau ',i1,' en ',t2-t1,' sec'

    enddo  ! fin boucle en niveau

    end subroutine pass_downward_HYB
! ------------------------------------------------------------------------------------------------
!
! DOWNWARD M2L+L2L version HF 3, pré-calculée selon init_v5 (classement base 7)
!
! ------------------------------------------------------------------------------------------------
    subroutine down_HF(NL,N,SPH,octree,k,i1)
! variables d'entrée
    integer, intent(in) :: NL
    integer, dimension (1:NL), intent(in) :: N
	type(FMM) :: SPH(NL)
	type(tree) :: octree(NL)
	real(prec) :: k
	integer :: i1
! variables de sortie dans SPH()%in
!	complex(prec),dimension(:) :: L2
! variables locales
    integer :: i2, i3, nxnynz(3), no
    integer :: npu_gro,npu_fin,nc_gro,nc_fin
    real(prec) :: cote
	complex(prec), dimension(:), allocatable :: M2L
	complex(prec), dimension(:), allocatable :: Lp1,Lp2
	real(prec), dimension(1:3) :: rlm,rl1l2
	real(prec), allocatable :: prodsca(:)
	integer, dimension(2) :: taille

    npu_gro=(N(i1)+1)*(2*N(i1)+1) ! nbre de points sphère unitaire niveau i1
    nc_gro=octree(i1)%NL  !% nombre de cellule par niveau ind
	allocate(Lp1(npu_gro))

    if (i1<NL) then
        npu_fin=(N(i1+1)+1)*(2*N(i1+1)+1) ! nbre de points sphère unitaire niveau i1+1
        nc_fin=octree(i1+1)%NL  !% nombre de cellule par niveau ind
        if (allocated(SPH(i1+1)%out)) deallocate(SPH(i1+1)%out)
        allocate(SPH(i1+1)%out(npu_fin*nc_fin))
        SPH(i1+1)%out=cmplx(0.0_prec,0.0_prec, kind=prec)
        allocate(Lp2(npu_fin))
    endif
    if (i1==1) then
        if (allocated(SPH(1)%out)) deallocate(SPH(1)%out)
        allocate(SPH(1)%out(npu_gro*nc_gro))
        SPH(1)%out=cmplx(0.0_prec,0.0_prec, kind=prec)
    endif

	cote=octree(i1)%cell(1)%length  !% dim des cotés des cellules

	do i2=1,octree(i1)%NL   ! boucle en cellules par niveau
		if (.not. octree(i1)%cell(i2)%a_cible) cycle   ! adaptatif : aucun noeud dans la cellule
	do i3=1,octree(i1)%cell(i2)%Nb_i  ! boucle en cellule en interaction avec cellule en cours
		if (.not. octree(i1)%cell(octree(i1)%cell(i2)%interac(i3))%a_source) cycle   ! adaptatif : multipole nul

          rlm=octree(i1)%cell(i2)%center - octree(i1)%cell(octree(i1)%cell(i2)%interac(i3))%center;
	      ! test redondance
	      nxnynz=nint(rlm/cote)
          no= sum((nxnynz+3)*(/49,7,1/))  ! conversion en base 7, cf init_FMM_v2
          M2L=SPH(i1)%tab_M2L(:,no)


             SPH(i1)%out((i2-1)*npu_gro+1:i2*npu_gro)=SPH(i1)%out((i2-1)*npu_gro+1:i2*npu_gro)+M2L * &
	 SPH(i1)%in((octree(i1)%cell(i2)%interac(i3)-1)*npu_gro+1:octree(i1)%cell(i2)%interac(i3)*npu_gro)

	end do   ! fin boucle cellule en interaction


    ! translation L2L vers niveau plus fin
    if (i1<NL .and. (interp_puis_dephase .or. interp_solveur == 1) .and. N(i1) /= N(i1+1)) then
       ! DEPHASER PUIS FILTRER (2026-09-29) : dephasage sur la grille du PARENT
       ! (bande N(i1) suffisante), puis filtrage vers la grille de l'enfant.
       ! Voir interp_puis_dephase.
       if (allocated(prodsca)) deallocate(prodsca)
       allocate(prodsca(npu_gro))
       do i3=1,octree(i1)%cell(i2)%Nb_c
           if (.not. octree(i1+1)%cell(octree(i1)%cell(i2)%child(i3))%a_cible) cycle   ! adaptatif : aucun noeud
           rl1l2=octree(i1+1)%cell(octree(i1)%cell(i2)%child(i3))%center - octree(i1)%cell(i2)%center
           prodsca=matmul(SPH(i1)%Ps,rl1l2)
           Lp1=SPH(i1)%out((i2-1)*npu_gro+1:i2*npu_gro)*exp(ji*k*prodsca)
           Lp2=cmplx(0.0_prec,0.0_prec, kind=prec)
           call L2L_FFT_v1(Lp1,SPH,NL,N,i1,Lp2)             ! grille parent -> grille enfant
           SPH(i1+1)%out((octree(i1)%cell(i2)%child(i3)-1)*npu_fin+1:octree(i1)%cell(i2)%child(i3)*npu_fin)=Lp2
       enddo
    elseif (i1<NL) then
       Lp1=SPH(i1)%out((i2-1)*npu_gro+1:i2*npu_gro)

		if (N(i1)==N(i1+1)) then
			Lp2=Lp1
		else
            Lp2=cmplx(0.0_prec,0.0_prec, kind=prec)
            call L2L_FFT_v1(Lp1,SPH,NL,N,i1,Lp2)
		endif

	do i3=1,octree(i1)%cell(i2)%Nb_c
		if (.not. octree(i1+1)%cell(octree(i1)%cell(i2)%child(i3))%a_cible) cycle   ! adaptatif : aucun noeud

		rl1l2=octree(i1+1)%cell(octree(i1)%cell(i2)%child(i3))%center - octree(i1)%cell(i2)%center;

        taille=shape(SPH(i1+1)%Ps)
		if (allocated(prodsca)) deallocate(prodsca)
		allocate(prodsca(taille(1)))
		call dgemm('N','N',taille(1),1,3,1.0_prec,SPH(i1+1)%Ps,taille(1),rl1l2,3,0.0_prec,prodsca,taille(1))

		SPH(i1+1)%out((octree(i1)%cell(i2)%child(i3)-1)*npu_fin+1:octree(i1)%cell(i2)%child(i3)*npu_fin)= &
            Lp2 * exp(ji*k*prodsca)
	enddo
	endif

	enddo   ! fin boucle cellule par niveau

    end subroutine down_HF
! ------------------------------------------------------------------------------------------------
!
! DOWNWARD M2L + L2L version BF, version pré-calculée
!
! ------------------------------------------------------------------------------------------------
    subroutine down_BF(NL,N,SPH,SPH0,octree,i1)
    implicit none
! variables d'entrée
    integer, intent(in) :: NL
    integer, dimension (1:NL), intent(in) :: N
	type(FMM) :: SPH(NL)
	type(FMM0) :: SPH0
	type(tree) :: octree(NL)
	integer :: i1,ind
! variables de sortie
! variables locales
	integer :: i2, i3, i4, no, Ni1, Ni2, nanm_gro, nanm_fin
	integer :: npos1, Nx, Ny, Nz, nxnynz(3), nc_gro, nc_fin
	real(prec) :: cote
	complex(prec), dimension((N(i1)+1)**2) :: Lp1,Lp2
	complex(prec), dimension(:), allocatable :: Lq_fin,Lq_gro,Lq_fin2

	real(prec), dimension(1:3) :: rlm,D21,rl2l1
	real(prec) :: R21
	! (supprimé) complex(prec), dimension(:), allocatable :: Rot,M2L
	! (supprimé) integer, dimension(:), allocatable :: ix,iy
	! (supprimé) real(prec), dimension(:), allocatable :: L2L
	! (supprimé) integer, dimension(:), allocatable :: M2L_i,M2L_j,L2L_i,L2L_j
	! DEBUG : character(len=1) :: answer
	! DEBUG : real(prec) :: t1,t2

	Ni1=N(i1)
	nanm_gro=(Ni1+1)**2              ! nombre de coefs HS correspondant au niveau i1
	nc_gro=octree(i1)%NL  !% nombre de cellule par niveau ind
	if (i1<NL) Ni2=N(i1+1)   ! lecture gardee (evite N(NL+1) hors bornes, valeur non utilisee si i1=NL)

	if (i1<NL) then
		nanm_fin=(Ni2+1)**2              ! nombre de coefs HS correspondant au niveau i1+1
		nc_fin=octree(i1+1)%NL  !% nombre de cellule par niveau ind
		if (allocated(SPH(i1+1)%outBF)) deallocate(SPH(i1+1)%outBF)
		allocate(SPH(i1+1)%outBF(nanm_fin*nc_fin))
		SPH(i1+1)%outBF=cmplx(0.0_prec,0.0_prec, kind=prec)
		allocate(Lq_gro(nanm_gro))
		allocate(Lq_fin2(nanm_fin))
		allocate(Lq_fin(nanm_fin))
	endif
	if (i1==1) then
		if (allocated(SPH(1)%outBF)) deallocate(SPH(1)%outBF)
		allocate(SPH(1)%outBF(nanm_gro*nc_gro))
		SPH(1)%outBF=cmplx(0.0_prec,0.0_prec, kind=prec)
	endif

   	cote=octree(i1)%cell(1)%length  !% dim des cotés des cellules

	do i2=1,octree(i1)%NL   ! boucle en cellules par niveau
		if (.not. octree(i1)%cell(i2)%a_cible) cycle   ! adaptatif : aucun noeud dans la cellule
		do i3=1,octree(i1)%cell(i2)%Nb_i  ! boucle en cellule en interaction avec cellule en cours
			if (.not. octree(i1)%cell(octree(i1)%cell(i2)%interac(i3))%a_source) cycle   ! adaptatif : multipole nul
			rlm=octree(i1)%cell(octree(i1)%cell(i2)%interac(i3))%center-octree(i1)%cell(i2)%center
			! test redondance
			R21=sqrt(dot_product(rlm,rlm)) ! norme
			D21=rlm/R21 ! direction
			! test Rotation déjà calculée
			nxnynz=nint(rlm/cote)
			no= sum((nxnynz+3)*(/49,7,1/))  ! conversion en base 7, cf init_FMM_v2
			do i4=1,15
				if (abs(SPH0%tab_T_M2L(i1,i4)-R21) < 1e-10_prec) then
					npos1=i4
					exit
				endif
			end do

			! rotation matrice creuse
			Lp1=cmplx(0.0_prec,0.0_prec, kind=prec)
			do i4=1,SPH0%rot_num(Ni1)
				Lp1(SPH0%rot_M2L_i(i4,no))=Lp1(SPH0%rot_M2L_i(i4,no))+SPH0%rot_M2L(i4,no) * &
					SPH(i1)%inBF((octree(i1)%cell(i2)%interac(i3)-1)*nanm_gro+SPH0%rot_M2L_j(i4,no))
			enddo
			! translation matrice creuse
			Lp2=cmplx(0.0_prec,0.0_prec, kind=prec)
			do i4=1,size(SPH(i1)%trans_M2L(:,npos1))
				Lp2(SPH(i1)%trans_M2L_i(i4,npos1))=Lp2(SPH(i1)%trans_M2L_i(i4,npos1))+ &
					Lp1(SPH(i1)%trans_M2L_j(i4,npos1))*SPH(i1)%trans_M2L(i4,npos1)
			enddo
			! rotation matrice creuse transposée conjuguée
			Lp1=cmplx(0.0_prec,0.0_prec, kind=prec)
			do i4=1,SPH0%rot_num(Ni1)
				Lp1(SPH0%rot_M2L_j(i4,no))=Lp1(SPH0%rot_M2L_j(i4,no))+ &
					Lp2(SPH0%rot_M2L_i(i4,no))*conjg(SPH0%rot_M2L(i4,no))
			enddo

			SPH(i1)%outBF((i2-1)*nanm_gro+1:i2*nanm_gro)=SPH(i1)%outBF((i2-1)*nanm_gro+1:i2*nanm_gro)+Lp1

		end do   ! fin boucle cellule en interaction
	end do

	if (i1 < NL) then
		do i2=1,octree(i1)%NL   ! boucle en cellules par niveau
			if (.not. octree(i1)%cell(i2)%a_cible) cycle   ! adaptatif : aucun noeud dans la cellule

			do i3=1,octree(i1)%cell(i2)%Nb_c  ! boucle sur les childs
				if (.not. octree(i1+1)%cell(octree(i1)%cell(i2)%child(i3))%a_cible) cycle   ! adaptatif : aucun noeud
				! calcul produit scalaire
				rl2l1= octree(i1)%cell(i2)%center - octree(i1+1)%cell(octree(i1)%cell(i2)%child(i3))%center
				R21=sqrt(dot_product(rl2l1,rl2l1)) ! norme
				D21=rl2l1/R21 ! direction
				! test Rotation déjà calculée
				! Bidouille Rotation déjà calculée, passage base 2
				Nx=nint(D21(1)*sqrt(3.0_prec)/2.0_prec+0.5_prec)
				Ny=nint(D21(2)*sqrt(3.0_prec)/2.0_prec+0.5_prec)
				Nz=nint(D21(3)*sqrt(3.0_prec)/2.0_prec+0.5_prec)
				ind=Nz+2*Ny+4*Nx+1   ! passage base 2 en base 10 : de 0 à 7 +1

				! translation M2M coaxiale
				! rotation matrice creuse
				Lq_gro=cmplx(0.0_prec,0.0_prec, kind=prec)
				do i4=1,SPH0%rot_num(Ni1)
					Lq_gro(SPH0%rot_M2M_i(i4,ind))=Lq_gro(SPH0%rot_M2M_i(i4,ind))+ &
					SPH(i1)%outBF((i2-1)*nanm_gro+SPH0%rot_M2M_j(i4,ind))*SPH0%rot_M2M(i4,ind)
				enddo
				! translation matrice creuse
				Lq_fin2=cmplx(0.0_prec,0.0_prec, kind=prec)
				do i4=1,size(SPH(i1)%trans_L2L)
					Lq_fin2(SPH(i1)%trans_L2L_i(i4))=Lq_fin2(SPH(i1)%trans_L2L_i(i4))+ &
					Lq_gro(SPH(i1)%trans_L2L_j(i4))*SPH(i1)%trans_L2L(i4)
				enddo
				! rotation matrice creuse transposée conjuguée
				Lq_fin=cmplx(0.0_prec,0.0_prec, kind=prec)
				do i4=1,SPH0%rot_num(Ni2)
					Lq_fin(SPH0%rot_M2M_j(i4,ind))=Lq_fin(SPH0%rot_M2M_j(i4,ind))+ &
					Lq_fin2(SPH0%rot_M2M_i(i4,ind))*conjg(SPH0%rot_M2M(i4,ind))
				enddo
				SPH(i1+1)%outBF((octree(i1)%cell(i2)%child(i3)-1)*nanm_fin+1:octree(i1)%cell(i2)%child(i3)*nanm_fin)= Lq_fin
			enddo
		enddo
	endif

	end subroutine down_BF
!
! ------------------------------------------------------------------------------------------------
!
! DOWNWARD M2L+L2L version HYB 3, pré-calculée selon init_v5 (classement base 7)
! M2L en HF, puis L2L en HF ou BF suivant SPH(Lev)%BF
!
! ------------------------------------------------------------------------------------------------
    subroutine down_HYB(NL,N,SPH,SPH0,octree,i1)
! variables d'entrée
    integer, intent(in) :: NL
    integer, dimension (1:NL), intent(in) :: N
	type(FMM) :: SPH(NL)
	type(FMM0) :: SPH0
	type(tree) :: octree(NL)
	integer :: i1
! variables de sortie dans SPH()%in
! variables locales
    integer :: i2,i3,i4,nxnynz(3),no,Nx,Ny,Nz,ind
    integer :: npu_gro,npu_fin,nc_gro,nc_fin,nanm_fin,nanm_gro
    real(prec) :: cote
	complex(prec), dimension(:), allocatable :: M2L
	complex(prec), dimension(:), allocatable :: Lp1,Lp2
	complex(prec), dimension(:), allocatable :: Lq_fin,Lq_gro,Lq_fin2

	real(prec), dimension(1:3) :: rlm, rl2l1, D21
	real(prec) :: R21
	! (supprimé) real(prec), allocatable :: prodsca(:)
	! (supprimé) integer, dimension(2) :: taille
	! DEBUG : real(prec) :: t1,t2


	npu_gro=(N(i1)+1)*(2*N(i1)+1) ! nbre de points sphère unitaire niveau i1
	nc_gro=octree(i1)%NL  !% nombre de cellule par niveau ind
	allocate(Lp1(npu_gro))

	if (i1<NL) then
		nc_fin=octree(i1+1)%NL  !% nombre de cellule par niveau ind
		if (SPH(i1+1)%BF==0) then
			npu_fin=(N(i1+1)+1)*(2*N(i1+1)+1) ! nbre de points sphère unitaire niveau i1+1
			if (allocated(SPH(i1+1)%out)) deallocate(SPH(i1+1)%out)
			allocate(SPH(i1+1)%out(npu_fin*nc_fin))
			SPH(i1+1)%out=cmplx(0.0_prec,0.0_prec, kind=prec)
			allocate(Lp2(npu_fin))
		else
			nanm_gro=(N(i1)+1)**2              ! nombre de coefs HS correspondant au niveau i1
			nanm_fin=(N(i1+1)+1)**2              ! nombre de coefs HS correspondant au niveau i1
			if (allocated(SPH(i1+1)%outBF)) deallocate(SPH(i1+1)%outBF)
			allocate(SPH(i1+1)%outBF(nanm_fin*nc_fin))
			SPH(i1+1)%outBF=cmplx(0.0_prec,0.0_prec, kind=prec)
			allocate(Lq_gro(nanm_gro))
			allocate(Lq_fin2(nanm_fin))
			allocate(Lq_fin(nanm_fin))
		endif
	endif
	if (i1==1) then
		if (allocated(SPH(1)%out)) deallocate(SPH(1)%out)
		allocate(SPH(1)%out(npu_gro*nc_gro))
		SPH(1)%out=cmplx(0.0_prec,0.0_prec, kind=prec)
	endif

	cote=octree(i1)%cell(1)%length  !% dim des cotés des cellules
    ! translation
	do i2=1,octree(i1)%NL   ! boucle en cellules par niveau
		if (.not. octree(i1)%cell(i2)%a_cible) cycle   ! adaptatif : aucun noeud dans la cellule

		do i3=1,octree(i1)%cell(i2)%Nb_i  ! boucle en cellule en interaction avec cellule en cours
			if (.not. octree(i1)%cell(octree(i1)%cell(i2)%interac(i3))%a_source) cycle   ! adaptatif : multipole nul

			rlm=octree(i1)%cell(i2)%center - octree(i1)%cell(octree(i1)%cell(i2)%interac(i3))%center
			! test redondance
			nxnynz=nint(rlm/cote)
			no= sum((nxnynz+3)*(/49,7,1/))  ! conversion en base 7, cf init_FMM_v2
			M2L=SPH(i1)%tab_M2L(:,no)

			SPH(i1)%out((i2-1)*npu_gro+1:i2*npu_gro)=SPH(i1)%out((i2-1)*npu_gro+1:i2*npu_gro)+M2L * &
	 SPH(i1)%in((octree(i1)%cell(i2)%interac(i3)-1)*npu_gro+1:octree(i1)%cell(i2)%interac(i3)*npu_gro)

		end do   ! fin boucle cellule en interaction
	end do

    ! translation L2L vers niveau plus fin
    ! conversion HF2BF
    call HF2BF_down_FMM(octree,NL,N,SPH,i1)    ! HF2BF de out en outBF
	! translation L2L
	do i2=1,octree(i1)%NL
		if (.not. octree(i1)%cell(i2)%a_cible) cycle   ! adaptatif : aucun noeud dans la cellule
		do i3=1,octree(i1)%cell(i2)%Nb_c  ! boucle sur les childs
			if (.not. octree(i1+1)%cell(octree(i1)%cell(i2)%child(i3))%a_cible) cycle   ! adaptatif : aucun noeud
			! calcul produit scalaire
			rl2l1= octree(i1)%cell(i2)%center - octree(i1+1)%cell(octree(i1)%cell(i2)%child(i3))%center
			R21=sqrt(dot_product(rl2l1,rl2l1)) ! norme
			D21=rl2l1/R21 ! direction
			! test Rotation déjà calculée
			! Bidouille Rotation déjà calculée, passage base 2
			Nx=nint(D21(1)*sqrt(3.0_prec)/2.0_prec+0.5_prec)
			Ny=nint(D21(2)*sqrt(3.0_prec)/2.0_prec+0.5_prec)
			Nz=nint(D21(3)*sqrt(3.0_prec)/2.0_prec+0.5_prec)
			ind=Nz+2*Ny+4*Nx+1   ! passage base 2 en base 10 : de 0 à 7 +1
			! translation M2M coaxiale
			! rotation matrice creuse
			Lq_gro=cmplx(0.0_prec,0.0_prec, kind=prec)
			do i4=1,SPH0%rot_num(N(i1))
				Lq_gro(SPH0%rot_M2M_i(i4,ind))=Lq_gro(SPH0%rot_M2M_i(i4,ind))+ &
					SPH(i1)%outBF((i2-1)*nanm_gro+SPH0%rot_M2M_j(i4,ind))*SPH0%rot_M2M(i4,ind)
			enddo
			! translation matrice creuse
			Lq_fin2=cmplx(0.0_prec,0.0_prec, kind=prec)
			do i4=1,size(SPH(i1)%trans_L2L)
				Lq_fin2(SPH(i1)%trans_L2L_i(i4))=Lq_fin2(SPH(i1)%trans_L2L_i(i4))+ &
					Lq_gro(SPH(i1)%trans_L2L_j(i4))*SPH(i1)%trans_L2L(i4)
			enddo
			! rotation matrice creuse transposée conjuguée
			Lq_fin=cmplx(0.0_prec,0.0_prec, kind=prec)
			do i4=1,SPH0%rot_num(N(i1+1))
				Lq_fin(SPH0%rot_M2M_j(i4,ind))=Lq_fin(SPH0%rot_M2M_j(i4,ind))+ &
					Lq_fin2(SPH0%rot_M2M_i(i4,ind))*conjg(SPH0%rot_M2M(i4,ind))
			enddo
			SPH(i1+1)%outBF((octree(i1)%cell(i2)%child(i3)-1)*nanm_fin+1:octree(i1)%cell(i2)%child(i3)*nanm_fin)= Lq_fin
		enddo
	enddo   ! fin boucle cellule par niveau

	end subroutine down_HYB
!
! ----------------------------------------------------------------------------------------
!
! Passage moment HF -> BF : downward, toutes les cellules
! SPH(Lev)%out -> SPH(Lev)%outBF
!
! ----------------------------------------------------------------------------------------
        subroutine HF2BF_down_FMM(octree,NL,N,SPH,Lev)
        implicit none
! variables d'entrée
	integer, intent(in) :: NL,Lev
	integer, dimension(1:NL), intent(in) :: N
	type(tree) :: octree(NL)

! variables de sortie
	type(FMM),intent(inout) :: SPH(NL)
! variables locales
	integer :: i1,n2,m2,ind,npu,nanm
	real(prec), allocatable, dimension(:,:) :: Psph
	complex(prec),allocatable,dimension(:,:) :: Ynm
	! DEBUG : character(len=1) :: answer

	npu=(N(Lev)+1)*(2*N(Lev)+1)  ! nbre de points sur sphère unitaire
	nanm=(N(Lev)+1)**2              ! nbre de coefs correspondants
	allocate(Psph(npu,3))
	call car2sph(SPH(Lev)%Ps,Psph,npu)

	if (allocated(SPH(Lev)%outBF)) deallocate(SPH(Lev)%outBF)
	allocate(SPH(Lev)%outBF(octree(Lev)%NL*nanm))
	allocate(Ynm(nanm,npu))
	call Ynm_gum(N(Lev),Psph(:,2),Psph(:,3),npu,Ynm)

	do n2=0,N(Lev)
		do m2=-n2,n2
			ind=n2**2+n2+1+m2
			do i1=1,octree(Lev)%NL 	! boucle cellules level
				SPH(Lev)%outBF((i1-1)*nanm+ind)=ji**(n2)*sum(SPH(Lev)%out((i1-1)*npu+1:i1*npu) * &
					conjg(Ynm(ind,:))*SPH(Lev)%wp)
			enddo
		enddo
	enddo

	end subroutine HF2BF_down_FMM
!
! ----------------------------------------------------------------------------------------
!
! Passage moment HF -> BF : upward pour test, toutes les cellules
! SPH(Lev)%in -> SPH(Lev)%inBF
!
! ----------------------------------------------------------------------------------------
        subroutine HF2BF_up_FMM(octree,NL,N,SPH,Lev)
        implicit none
! variables d'entrée
	integer, intent(in) :: NL,Lev
	integer, dimension(1:NL), intent(in) :: N
	type(tree) :: octree(NL)

! variables de sortie
	type(FMM),intent(inout) :: SPH(NL)
! variables locales
	integer :: i1,n2,m2,ind,npu,nanm
	real(prec), allocatable, dimension(:,:) :: Psph
	complex(prec),allocatable,dimension(:,:) :: Ynm
	! DEBUG : character(len=1) :: answer

	npu=(N(Lev)+1)*(2*N(Lev)+1)  ! nbre de points sur sphère unitaire
	nanm=(N(Lev)+1)**2              ! nbre de coefs correspondants
	allocate(Psph(npu,3))
	call car2sph(SPH(Lev)%Ps,Psph,npu)

	if (allocated(SPH(Lev)%inBF)) deallocate(SPH(Lev)%inBF)
	allocate(SPH(Lev)%inBF(octree(Lev)%NL*nanm))
	allocate(Ynm(nanm,npu))
	call Ynm_gum(N(Lev),Psph(:,2),Psph(:,3),npu,Ynm)

	do n2=0,N(Lev)
		do m2=-n2,n2
			ind=n2**2+n2+1+m2
			do i1=1,octree(Lev)%NL 	! boucle cellules level
				SPH(Lev)%inBF((i1-1)*nanm+ind)=ji**(n2)*sum(SPH(Lev)%in((i1-1)*npu+1:i1*npu) * &
					conjg(Ynm(ind,:))*SPH(Lev)%wp)
			enddo
		enddo
	enddo

	end subroutine HF2BF_up_FMM
!
! -----------------------------------------------------------------------------------------------------------------------------
! calcul des T_LR : fonction de transfert pour tab_M2L
! Version FFT
! -----------------------------------------------------------------------------------------------------------------------------
	subroutine T_LR(NL,L,LEV,k,SPH,rlm,T)
	implicit none
! variables d'entrée	
	integer, intent(in) :: NL,L,LEV
	real(prec), intent(in) :: k,rlm(3)
	type(FMM) :: SPH(NL)
! variables de sortie
	complex(prec) :: T(size(SPH(LEV)%Ps,1))
! variables locales
	integer, dimension(0:L) :: indices
	integer :: i1,NM
    real(prec) :: SY(0:250),DY(0:250),SJ(0:250),DJ(0:250)
 	real(prec), dimension(size(SPH(LEV)%Ps,1)) :: ps_sprlm,cos_sr,n_sp
	real(prec), dimension(size(SPH(LEV)%Ps,1),1:L+1) :: PL_L
	real(prec) :: n_rlm
	complex(prec), dimension(0:L) :: H1,coef
	complex(prec) :: j
	! DEBUG : character(len=1) :: answer
	
	j=cmplx(0.0_prec,1.0_prec, kind=prec)
	
!% produit scalaire sp.rlm
	ps_sprlm=matmul(SPH(LEV)%Ps,rlm)
!% norme de sp
	n_sp=sqrt(sum(SPH(LEV)%Ps**2,dim=2))
!% norme de rlm
	n_rlm=sqrt(sum(rlm**2))
	cos_sr=ps_sprlm/n_sp/n_rlm
! calcul des polynomes de legendre
	call poly_leg(L,cos_sr,size(SPH(LEV)%Ps,1),PL_L)
! calcul fonction de hankel
	CALL SPHY(L,k*n_rlm,NM,SY,DY)
	if (k*n_rlm>3000.0_prec) then
		call sph_j1(L,k*n_rlm,SJ)
	else
		CALL SPHJ(L,k*n_rlm,NM,SJ,DJ)
	endif
	H1(:)=SJ(0:L)+j*SY(0:L)

	indices=[(i1,i1=0,L)]
	coef=(2.0_prec*indices+1.0_prec) * (j**indices) * H1
	T=matmul(PL_L,coef)
	!call gemm('N','N',size(SPH(LEV)%Ps,1),1,L+1,cmplx(1.0_prec,0.0_prec, kind=prec),PL_L,size(SPH(LEV)%Ps,1),coef,L+1,cmplx(0.0_prec,0.0_prec, kind=prec),&
	!T,size(SPH(LEV)%Ps,1))
	end subroutine T_LR
!
! ---------------------------------------------------------------------------------------
!
! Translation M2L par recurrence, coaxiale translation en z (positif : module)
! complexité O(p^3)
! Version BF
!
! ---------------------------------------------------------------------------------------
    subroutine M2L_BF(T_M2L,T_M2L_i,T_M2L_j,num,k,N11,Z)
        implicit none
! variables d'entrée
        integer,intent(in) :: N11
        real(prec),intent(inout) :: k,Z
! variables de sortie
        integer, intent(out) :: num
        complex(prec), intent(out),allocatable :: T_M2L(:)
        integer, intent(out),allocatable :: T_M2L_i(:),T_M2L_j(:)


! variables locales
        complex(prec) :: T0_M2L((N11+1)**2 * (N11+1))
        integer :: T0_M2L_i((N11+1)**2 * (N11+1)),T0_M2L_j((N11+1)**2 * (N11+1))
        complex(prec) :: trans_M2L((2*N11)**2+2*N11+1,(2*N11)**2+2*N11+1)
        complex(prec) :: E2_M2M,E3_M2M,E4_M2M
        integer :: NM
        integer :: l2,n1,m1,ind1,ind2,P1
        integer :: n1p, n1m, l2p, l2m, m1p, m0
        integer :: ind3,ind4,ind5,ind6,ind7,ind8,ind01,ind02
        real(prec) :: jn1(0:2*N11),Djn1(0:2*N11)
        real(prec) :: yn(0:2*N11),Dyn(0:2*N11)
        real(prec) :: a1,a2,a3,a4,b1,b2,b3
        complex(prec) :: hn(0:2*N11)
        ! DEBUG : character(len=1) :: answer

        ! récurrence pour tranlation
        P1=N11+1
        CALL SPHJ(2*P1-2,k*Z,NM,jn1,Djn1)
        CALL SPHY(2*P1-2,k*Z,NM,yn,Dyn)
        hn=jn1 + ji * yn
        n1=0
        m1=0
        trans_M2L=cmplx(0.0_prec,0.0_prec, kind=prec)
        ! initialisation 4.83
        do l2=0,2*P1-2
            ind1=n1**2+n1+1+m1
            ind2=l2**2+l2+1+m1
            trans_M2L(ind1,ind2)=(-1.0_prec)**l2*sqrt(real(2*l2+1, prec))*hn(l2)
            ! par symétrie
            if (ind2>1) then
                trans_M2L(ind2,ind1)=1.0_prec/(-1.0_prec)**(n1+l2)*trans_M2L(ind1,ind2)
            endif
        enddo
!        write(*,*) 'pass 2'

        ! on avance à n=m eq 4.84
        do m1=0,P1-1
            do l2=m1+1,2*P1-3-m1
                n1=m1
                n1p=n1+1
                n1m=n1-1
                l2p=l2+1
                l2m=l2-1
                m1p=m1+1
                m0=-m1-1
                b1=-sqrt(real(n1p-m0-1, prec)*real(n1p-m0, prec)/real(2*n1p-1, prec)/real(2*n1p+1, prec))
                b2=-sqrt(real(l2-m0-1, prec)*real(l2-m0, prec)/real(2*l2-1, prec)/real(2*l2+1, prec))
                b3=sqrt(real(l2p-m1-1, prec)*real(l2p-m1, prec)/real(2*l2p-1, prec)/real(2*l2p+1, prec))
                ind1=n1p**2+n1p+1+m1p
                ind2=l2**2+l2+1+m1p
                ind3=m1**2+m1+1+m1
                ind4=l2m**2+l2m+1+m1
                ind5=m1**2+m1+1+m1
                ind6=l2p**2+l2p+1+m1
                trans_M2L(ind1,ind2)=1.0_prec/b1*(b2*trans_M2L(ind3,ind4)-b3*trans_M2L(ind5,ind6))
                !symétrie
                ind01=n1p**2+n1p+1+m0
                ind02=l2**2+l2+1+m0
                trans_M2L(ind01,ind02)=trans_M2L(ind1,ind2)
                if (n1p/=l2) then
                    trans_M2L(ind2,ind1)=1.0_prec/(-1.0_prec)**(n1p+l2)*trans_M2L(ind1,ind2)
                    trans_M2L(ind02,ind01)=trans_M2L(ind2,ind1)
                endif
            enddo
        enddo
!        write(*,*) 'pass 3'

        ! on avance en n2
        do m1=0,P1-1
            do n1=m1,2*P1-3-m1
                do l2=n1+1,2*P1-3-n1
                    l2p=l2+1
                    l2m=l2-1
                    n1m=n1-1
                    n1p=n1+1
                    ind1=n1p**2+n1p+1+m1
                    ind2=l2**2+l2+1+m1
                        a1=sqrt(real(n1+1+m1, prec)*real(n1+1-m1, prec)/real(2*n1+1, prec)/real(2*n1+3, prec))
                        a2=sqrt(real(l2m+1+m1, prec)*real(l2m+1-m1, prec)/real(2*l2m+1, prec)/real(2*l2m+3, prec))
                        a3=sqrt(real(l2+1+m1, prec)*real(l2+1-m1, prec)/real(2*l2+1, prec)/real(2*l2+3, prec))
                        if (n1-1<m1) then
                            a4=0.0_prec
                        else
                            a4=sqrt(real(n1m+1+m1, prec)*real(n1m+1-m1, prec)/real(2*n1m+1, prec)/real(2*n1m+3, prec))
                        endif
                        ind3=n1**2+n1+1+m1
                        ind4=l2m**2+l2m+1+m1
                        ind5=n1**2+n1+1+m1
                        ind6=l2p**2+l2p+1+m1
                        ind7=n1m**2+n1m+1+m1
                        ind8=l2**2+l2+1+m1
                        E2_M2M=trans_M2L(ind3,ind4)
                        E3_M2M=trans_M2L(ind5,ind6)
                        E4_M2M=trans_M2L(ind7,ind8)
                    trans_M2L(ind1,ind2)=1.0_prec/a1*(a2*E2_M2M-a3*E3_M2M+a4*E4_M2M)
                    ! symétrie
                    if (m1>0) then
                        ind01=n1p**2+n1p+1-m1
                        ind02=l2**2+l2+1-m1
                        trans_M2L(ind01,ind02)=trans_M2L(ind1,ind2);
                    endif
                    if (ind1/=ind2) then
                        trans_M2L(ind2,ind1)=1.0_prec/(-1.0_prec)**(n1p+l2)*trans_M2L(ind1,ind2)
                       if (m1>0) then
                        trans_M2L(ind02,ind01)=trans_M2L(ind2,ind1)
                       endif
                    endif
                enddo
            enddo
        enddo

        num=0
        T0_M2L=cmplx(0.0_prec,0.0_prec, kind=prec)
        do n1=0,N11
            do m1=-n1,n1
                ind1=n1**2+n1+1+m1
                do l2=abs(m1),N11
                            num=num+1
                            ind2=l2**2+l2+1+m1
                            T0_M2L(num)=trans_M2L(ind2,ind1)
                            T0_M2L_i(num)=ind2
                            T0_M2L_j(num)=ind1
                enddo
            enddo
        enddo
        if (allocated(T_M2L)) deallocate(T_M2L)
        allocate(T_M2L(num))
        if (allocated(T_M2L_i)) deallocate(T_M2L_i)
        allocate(T_M2L_i(num))
        if (allocated(T_M2L_j)) deallocate(T_M2L_j)
        allocate(T_M2L_j(num))
        T_M2L=T0_M2L(1:num)
        T_M2L_i=T0_M2L_i(1:num)
        T_M2L_j=T0_M2L_j(1:num)

    end subroutine M2L_BF
! -----------------------------------------------------------------------------------------------------------------------------
! calcul des M2M en MLFMM en FFT version diagonale
! -----------------------------------------------------------------------------------------------------------------------------
	subroutine M2M_FFT_v1(Mq,SPH,NL,N,ind,Mq2)
	!use init_param_FMM
	implicit none
! variables d'entrée
	integer, intent(in) :: NL,ind
	integer, dimension (1:NL), intent(in) :: N
!	complex(prec), intent(inout) :: Mq(:)
	complex(prec), intent(inout) :: Mq((N(ind+1)+1)*(2*N(ind+1)+1))
	type(FMM) :: SPH(NL)

! variable de sortie
!	complex(prec),allocatable, intent(inout) :: Mq2(:)
	complex(prec), intent(inout) :: Mq2((N(ind)+1)*(2*N(ind)+1))
! variables locales
	integer :: npu_fin, npu_gro,lensav,ier,lenwrk,jump,inc,m,i1
	integer :: npu_fin2,npu_gro2
	real(prec), allocatable :: wsave(:)
	real(prec), allocatable :: work(:)
!	real(prec), allocatable :: mat(:,:)
!	complex(prec), allocatable :: wF(:)
	! (supprimé) real(prec) :: mat((N(ind)+1),(N(ind+1)+1))
	complex(prec) :: wF((N(ind+1)+1)*(2*N(ind+1)+1))
	! DEBUG : character(len=1) :: answer
	complex(prec) :: Mqtilde(2*N(ind+1)+1,N(ind)+1)
	!complex(prec),dimension(N(ind)+1,-N(ind):N(ind)) :: Mqtilde2
	complex(prec),dimension(2*N(ind)+1,N(ind)+1) :: Mqtilde2

!	write(*,*) 'ind=',ind
	npu_fin=(N(ind+1)+1) ! nbre de points sphère unitaire niveau i1
	npu_fin2=2*N(ind+1)+1
	npu_gro=(N(ind)+1) ! nbre de points sphère unitaire niveau i1-1
	npu_gro2=2*N(ind)+1
!	write(*,*) shape(SPH(ind)%IR),ind,npu_gro,npu_fin
!	write(*,*) SPH(ind)%IR(1:npu_gro)
!	read(*,*) answer
	!allocate(mat(npu_gro,npu_fin))
	!allocate(wF(npu_fin))
	!mat=reshape(SPH(ind)%IR,(/ npu_gro,npu_fin /));
	! multiplication par les poids de Gauss/Legendre sphère
	wF=SPH(ind+1)%wp * Mq
	!write(*,*) 'shape wp =',shape(SPH(ind+1)%wp), 'shape(Mq) = ',shape(Mq),&
	!'shape(wF) = ',shape(wF)

	!write(*,*) wF
	!read(*,*) answer

	! FFTF init
	lensav=2*(npu_fin2) + INT(LOG(REAL(npu_fin2))/LOG(2.0_prec)) + 4
	!write(*,*) 'lensav = ',lensav
	allocate(wsave(lensav))
	call cfftmi ( npu_fin2, wsave, lensav, ier )
	if (ier==2) then
		write(*,*) 'lensav n''est pas assez grand M2M_FFT forward'
	endif
	! FFT Foward
	lenwrk=2*(npu_fin2)*(npu_fin)
	allocate(work(lenwrk))
	jump=npu_fin2
	inc=1
	!write(*,*) 'lot =',npu_fin
	!write(*,*) 'jump =',jump
	!write(*,*) 'N = ',npu_fin2
        !write(*,*) 'inc = ',inc

	!write(*,*) 'lenc=',npu_fin*npu_fin2
	!write(*,*) ' lenc < ',(npu_fin-1) * jump + inc * (npu_fin2-1) +1

	call cfftmb ( npu_fin, jump,npu_fin2, inc, wF, (N(ind+1)+1)*(2*N(ind+1)+1), wsave, lensav, work, &
  lenwrk, ier )
  	if (ier.ne.0) then
	    write(*,*) 'erreur après fft forward : ',ier
	endif
	!write(*,*) 'passe FFT F'
	!write(*,*) wF(1:9)

	!write(*,*) 'npu_fin = ',npu_fin
        !read(*,*) answer

	jump=npu_fin2
	do i1=1,npu_fin2
!	do i2=-npu_fin:npu_fin
            if (i1 < (npu_fin+1)) then
                m=1-i1
            else
                m=-((i1-npu_fin)-npu_fin)
            endif
         !   write(*,*) 'ind = ',ind,'i1 ='    ,i1,' m =',m
         !   write(*,*) 'wf = ',wF((i1-1)+1:(N(ind+1)+1)*(2*N(ind+1)+1):jump)
            !write(*,*) 'mat_Q in =',SPH(ind)%mat_Q(:,:,abs(m))
	   !m=-npu_fin+i1-1
!	   write(*,*) shape(SPH(ind)%mat_Q(:,:,abs(m)))
!	   write(*,*) shape(wF((i1-1)+1:(N(ind+1)+1)*(2*N(ind+1)+1):jump))
!	   write(*,*) shape(Mqtilde(:,i1))
	   !write(*,*) 'mat_Q =',SPH(ind)%mat_Q
           Mqtilde(i1,:)=matmul(SPH(ind)%mat_Q(:,:,abs(m)),wF((i1-1)+1:(N(ind+1)+1)*(2*N(ind+1)+1):jump))
	enddo
	!write(*,*) shape(Mqtilde)
	!write(*,*) shape(SPH(ind)%mat_Q)
	!write(*,*) 'Mqtilde1 =',Mqtilde(1:3,1)
	!write(*,*) 'Mqtilde2 =',Mqtilde(1:3,2)
	!write(*,*) 'Mqtilde3 =',Mqtilde(1:3,3)
	!write(*,*) 'Mqtilde4 =',Mqtilde(1:3,4)
	!write(*,*) 'Mqtilde5 =',Mqtilde(1:3,5)
	!write(*,*) 'Mqtilde6 =',Mqtilde(1:3,6)
	!write(*,*) 'Mqtilde7 =',Mqtilde(1:3,7)
	!write(*,*) 'Mqtilde8 =',Mqtilde(:,8)
	!write(*,*) 'Mqtilde9 =',Mqtilde(:,9)
	!read(*,*) answer
	! FFT backward
	lensav=2*(npu_gro2) + INT(LOG(REAL(npu_gro2))/LOG(2.0_prec)) + 4
	deallocate(wsave)
	allocate(wsave(lensav))
	call cfftmi ( npu_gro2, wsave, lensav, ier )
	if (ier==2) then
	    write(*,*) 'lensav n''est pas assez grand M2M_FFT backward'
	endif

	lenwrk=2*(npu_gro2)*(npu_gro)
	deallocate(work)
	allocate(work(lenwrk))

	jump=npu_gro2
	inc=1

	! zeros-padding
	Mqtilde2=cmplx(0.0_prec,0.0_prec, kind=prec)
	!write(*,*) 'ind1 = ',N(ind+1)+1
	Mqtilde2(1:N(ind+1)+1,:)=Mqtilde(1:N(ind+1)+1,:)
	!write(*,*) 'ind1 = ',N(ind)+1,'ind2 = ',2*N(ind)+1
	Mqtilde2(2*N(ind)+1-N(ind+1)+1:2*N(ind)+1,:)=Mqtilde(N(ind+1)+2:2*N(ind+1)+1,:)
	!write(*,*) 'Mqtilde2 =',Mqtilde2(:,1)
	!Mqtilde2(:,-N(ind+1):N(ind+1))=Mqtilde
	!Mq2=reshape(transpose(Mqtilde2),(/(N(ind)+1)*(2*N(ind)+1)/))
	Mq2=reshape(Mqtilde2,(/(N(ind)+1)*(2*N(ind)+1)/))
	call cfftmf ( npu_gro, jump, npu_gro2, inc, Mq2, (N(ind)+1)*(2*N(ind)+1), wsave, lensav, work, &
  lenwrk, ier )
  	if (ier.ne.0) then
	    write(*,*) 'erreur après fft backward : ',ier
	endif
	!write(*,*) 'Mq2 =',npu_gro2*Mq2(1:10)
	!read(*,*) answer
	Mq2=Mq2*npu_gro2
!	write(*,*) 'shape de Mq2',shape(Mq2)
!	read(*,*) answer
!	deallocate(mat)
!	deallocate(wF)
	end subroutine M2M_FFT_v1
! -----------------------------------------------------------------------------------------------------------------------------
! calcul des L2L en MLFMM par FFT version diagonale
! -----------------------------------------------------------------------------------------------------------------------------
	subroutine L2L_FFT_v1(Mq,SPH,NL,N,ind,Mq2)
!	use init_param_FMM
!	implicit none
! variables d'entrée
	integer, intent(in) :: NL,ind
	integer, dimension (1:NL), intent(in) :: N
!	complex(prec), intent(inout) :: Mq(:)
	complex(prec), intent(inout) :: Mq((N(ind)+1)*(2*N(ind)+1))
	type(FMM) :: SPH(NL)

! variable de sortie
!	complex(prec),allocatable, intent(inout) :: Mq2(:)
	complex(prec), intent(inout) :: Mq2((N(ind+1)+1)*(2*N(ind+1)+1))
! variables locales
	integer :: npu_fin, npu_gro,lensav,ier,lenwrk,jump,inc,m,i1
	integer :: npu_fin2,npu_gro2
	real(prec), allocatable :: wsave(:)
	real(prec), allocatable :: work(:),work2(:)
!	real(prec), allocatable :: mat(:,:)
!	complex(prec), allocatable :: wF(:)
	! (supprimé) real(prec) :: mat((N(ind+1)+1),(N(ind)+1))
	complex(prec) :: wF((N(ind)+1)*(2*N(ind)+1))
	complex(prec) :: wFt((2*N(ind)+1),(N(ind)+1))
	complex(prec) :: wF2t(2*N(ind+1)+1,N(ind)+1)
	! DEBUG : character(len=1) :: answer
	! (supprimé) complex(prec) :: Mqtilde(2*N(ind)+1,N(ind+1)+1)
	!complex(prec),dimension(N(ind)+1,-N(ind):N(ind)) :: Mqtilde2
	complex(prec),dimension(2*N(ind+1)+1,N(ind+1)+1) :: Mqtilde2

 !   write(*,*) 'in L2L'
!	write(*,*) 'ind=',ind
	npu_fin=(N(ind+1)+1) ! nbre de points sphère unitaire niveau i1
	npu_fin2=2*N(ind+1)+1
	npu_gro=(N(ind)+1) ! nbre de points sphère unitaire niveau i1-1
	npu_gro2=2*N(ind)+1

	!write(*,*) 'npu_fin =',npu_fin,' npu_fin2 =',npu_fin2
	!write(*,*) 'npu_gro =',npu_gro,' npu_gro2 =',npu_gro2
	! multiplication par les poids de Gauss/Legendre sphère
	wF=SPH(ind)%wp * Mq


	! FFTF init
	lensav=2*(npu_gro2) + INT(LOG(REAL(npu_gro2))/LOG(2.0_prec)) + 4
	!write(*,*) 'lensav = ',lensav
	allocate(wsave(lensav))
	call cfftmi ( npu_gro2, wsave, lensav, ier )
	if (ier==2) then
	    write(*,*) 'lensav n''est pas assez grand M2M_FFT forward'
	endif
	! FFT Foward
	lenwrk=2*(npu_gro2)*(npu_gro)
	allocate(work(lenwrk))
	jump=npu_gro2
	inc=1
!	write(*,*) 'lot =',npu_fin
!	write(*,*) 'jump =',jump
!	write(*,*) 'N = ',npu_fin2
	!write(*,*) 'inc = ',inc

	!write(*,*) 'lenc=',npu_fin*npu_fin2
	!write(*,*) ' lenc < ',(npu_fin-1) * jump + inc * (npu_fin2-1) +1

	call cfftmb ( npu_gro, jump,npu_gro2, inc, wF, (npu_gro*npu_gro2), wsave, lensav, work, &
  lenwrk, ier )
  	if (ier.ne.0) then
	    write(*,*) 'erreur après fft forward : ',ier
	endif
 	wFt=reshape(wF,(/npu_gro2,npu_gro/))

	! on tronque
    	wF2t(1:npu_fin,:)=wFt(1:npu_fin,:)
    	wF2t(npu_fin+1:npu_fin2,:)=wFt(npu_gro2-npu_fin+2:npu_gro2,:)
!	write(*,*) wF(1:9)

!	write(*,*) 'passe FFT F'
	!write(*,*) 'npu_fin = ',npu_fin
!        read(*,*) answer

	jump=npu_fin2
	!write(*,*) 'shape(Mqtilde) =',shape(Mqtilde)
    	!write(*,*) 'shape(mat_Q) =',shape(SPH(ind)%mat_Q)
	!write(*,*) 'npu_gro2=',npu_gro2
	do i1=1,npu_fin2
!	do i2=-npu_fin:npu_fin
  !          write(*,*) 'i1=',i1
            if (i1 < (npu_fin+1)) then
                m=1-i1
            else
                m=i1-(npu_fin2+1)
            endif
    !        write(*,*) 'ind = ',ind,'i1 ='    ,i1,' m =',m
         !   write(*,*) 'wf = ',wF((i1-1)+1:(N(ind+1)+1)*(2*N(ind+1)+1):jump)
            !write(*,*) 'mat_Q in =',SPH(ind)%mat_Q(:,:,abs(m))
	   !m=-npu_fin+i1-1
	   !write(*,*) shape(transpose(SPH(ind)%mat_Q(:,:,abs(m))))
	   !write(*,*) shape(wF((i1-1)+1:(npu_gro*npu_gro2):jump))
	   !write(*,*) shape(Mqtilde(i1,:))
	   !write(*,*) 'abs(m)=',abs(m)
	   !mat=transpose(SPH(ind)%mat_Q(:,:,abs(m)))
	   !write(*,*) 'mat =',shape(mat)
 !          Mqtilde(i1,:)=matmul(transpose(SPH(ind)%mat_Q(:,:,abs(m))),wF((i1-1)+1:(npu_fin*npu_fin2):jump))
           Mqtilde2(i1,:)=matmul(transpose(SPH(ind)%mat_Q(:,:,abs(m))),wF2t(i1,:))    ! version tronquée
	enddo
	!write(*,*) shape(Mqtilde)
	!write(*,*) shape(SPH(ind)%mat_Q)
	!write(*,*) 'Mqtilde1 =',Mqtilde(1:3,1)
	!write(*,*) 'Mqtilde2 =',Mqtilde(1:3,2)
	!write(*,*) 'Mqtilde3 =',Mqtilde(1:3,3)
	!write(*,*) 'Mqtilde4 =',Mqtilde(1:3,4)
	!write(*,*) 'Mqtilde5 =',Mqtilde(1:3,5)
	!write(*,*) 'Mqtilde6 =',Mqtilde(1:3,6)
	!write(*,*) 'Mqtilde7 =',Mqtilde(1:3,7)
	!write(*,*) 'Mqtilde8 =',Mqtilde(:,8)
	!write(*,*) 'Mqtilde9 =',Mqtilde(:,9)
	!read(*,*) answer
	!write(*,*) 'passe matmul'
	! FFT backward
	lensav=2*(npu_fin2) + INT(LOG(REAL(npu_fin2))/LOG(2.0_prec)) + 4
	deallocate(wsave)
	allocate(wsave(lensav))
	call cfftmi ( npu_fin2, wsave, lensav, ier )
	if (ier==2) then
	    write(*,*) 'lensav n''est pas assez grand M2M_FFT backward'
	endif
 !       write(*,*) 'passe init'
	lenwrk=2*(npu_fin2)*(npu_fin)
!	write(*,*) 'lenwrk=',lenwrk
	!deallocate(work)
!	write(*,*) 'alloc'
	allocate(work2(lenwrk))

	jump=npu_fin2
	inc=1

!	write(*,*) 'passe alloc'
	! troncature
!	Mqtilde2=cmplx(0.0_prec,0.0_prec, kind=prec)
	!write(*,*) 'ind1 = ',N(ind+1)+1
!	Mqtilde2(1:N(ind+1)+1,:)=Mqtilde(1:N(ind+1)+1,:)
	!write(*,*) 'ind1 = ',N(ind)+1,'ind2 = ',2*N(ind)+1
!	Mqtilde2(N(ind+1)+2:2*N(ind+1)+1,:)=Mqtilde(2*N(ind)+1-N(ind+1)+1:2*N(ind)+1,:)
	!write(*,*) 'Mqtilde2 =',Mqtilde2(:,1)
	!Mqtilde2(:,-N(ind+1):N(ind+1))=Mqtilde
	!Mq2=reshape(transpose(Mqtilde2),(/(N(ind)+1)*(2*N(ind)+1)/))
	Mq2=reshape(Mqtilde2,(/(npu_fin*npu_fin2)/))
!	write(*,*) 'passe troncature'
	call cfftmf ( npu_fin, jump, npu_fin2, inc, Mq2, npu_fin*npu_fin2, wsave, lensav, work2, &
  lenwrk, ier )
  	if (ier.ne.0) then
	    write(*,*) 'erreur après fft backward : ',ier
	endif
	!write(*,*) 'Mq2 =',npu_gro2*Mq2(1:13*7)
	!read(*,*) answer
	Mq2=Mq2*npu_fin2
!	write(*,*) 'Mq2 =',Mq2(1:10)
!	read(*,*) answer

	end subroutine L2L_FFT_v1
! ---------------------------------------------------------------------------------------
!
! Translation M2M et L2L par recurrence, coaxiale translation en z (positif : module)
! complexité O(p^3)
!
! Passe d'un ordre N(i1) à N(i1-1) : M2M : fin à gros : N(i1) <= N(i1-1)
! Passe d'un ordre N(i1-1) à N(i1) : L2L : gros à fin : N(i1-1) >= N(i1)
!
! TODO : Supprimer L2L car identique à M2M, gain mémoire
! ---------------------------------------------------------------------------------------
    subroutine M2M_BF(NL,k,N,SPH,SPH0,i1)
        implicit none
! variables d'entrée
        integer,intent(in) :: i1, NL, N(NL)
        real(prec),intent(in) :: k
        type(FMM), intent(inout) :: SPH(NL)
        type(FMM0), intent(inout) :: SPH0
! variables locales
        integer :: N11
        real(prec) :: Z
        real(prec) :: T0_M2M((N(i1-1)+1)**2 * (N(i1-1)+1))
        integer :: T0_M2M_i((N(i1-1)+1)**2 * (N(i1-1)+1)),T0_M2M_j((N(i1-1)+1)**2 * (N(i1-1)+1))
        real(prec) :: T0_L2L((N(i1-1)+1)**2 * (N(i1-1)+1))
        integer :: T0_L2L_i((N(i1-1)+1)**2 * (N(i1-1)+1)),T0_L2L_j((N(i1-1)+1)**2 * (N(i1-1)+1))
        real(prec) :: trans_M2L((2*N(i1-1))**2+2*N(i1-1)+1,(2*N(i1-1))**2+2*N(i1-1)+1)
        real(prec) :: E2_M2M,E3_M2M,E4_M2M
        integer :: NM,numM,numL
        integer :: l2, n1, m1, ind1, ind2, P1
        integer :: n1p, n1m, l2p, l2m, m1p, m0
        integer :: ind3,ind4,ind5,ind6,ind7,ind8,ind01,ind02
        real(prec) :: jn(0:2*N(i1-1)),Djn(0:2*N(i1-1))
        real(prec) :: a1,a2,a3,a4,b1,b2,b3
        ! DEBUG : character(len=1) :: answer


        ! récurrence pour tranlation
        N11=N(i1-1)   ! le plus grand des deux ordres (ou égaux)
        P1=N11+1
!        write(*,*) 'SPH0%tab_T_M2M=',SPH0%tab_T_M2M(:)
        Z=SPH0%tab_T_M2M(i1)
!        write(*,*) 'Z=',Z,' k=',k

        CALL SPHJ(2*P1-2,k*Z,NM,jn,Djn)

        n1=0
        m1=0
        trans_M2L=cmplx(0.0_prec,0.0_prec, kind=prec)

        ! initialisation 4.83
        do l2=0,2*P1-2
            ind1=n1**2+n1+1+m1
            ind2=l2**2+l2+1+m1
            trans_M2L(ind1,ind2)=(-1.0_prec)**l2*sqrt(real(2*l2+1, prec))*jn(l2)
            ! par symétrie
            if (ind2>1) then
                trans_M2L(ind2,ind1)=1.0_prec/(-1.0_prec)**(n1+l2)*trans_M2L(ind1,ind2)
            endif
        enddo
!        write(*,*) 'pass 2'

        ! on avance à n=m eq 4.84
        do m1=0,P1-1
            do l2=m1+1,2*P1-3-m1
                n1=m1
                n1p=n1+1
                n1m=n1-1
                l2p=l2+1
                l2m=l2-1
                m1p=m1+1
                m0=-m1-1
                b1=-sqrt(real(n1p-m0-1, prec)*real(n1p-m0, prec)/real(2*n1p-1, prec)/real(2*n1p+1, prec))
                b2=-sqrt(real(l2-m0-1, prec)*real(l2-m0, prec)/real(2*l2-1, prec)/real(2*l2+1, prec))
                b3=sqrt(real(l2p-m1-1, prec)*real(l2p-m1, prec)/real(2*l2p-1, prec)/real(2*l2p+1, prec))
                ind1=n1p**2+n1p+1+m1p
                ind2=l2**2+l2+1+m1p
                ind3=m1**2+m1+1+m1
                ind4=l2m**2+l2m+1+m1
                ind5=m1**2+m1+1+m1
                ind6=l2p**2+l2p+1+m1

                trans_M2L(ind1,ind2)=1.0_prec/b1*(b2*trans_M2L(ind3,ind4)-b3*trans_M2L(ind5,ind6))

                !symétrie
                ind01=n1p**2+n1p+1+m0
                ind02=l2**2+l2+1+m0
                trans_M2L(ind01,ind02)=trans_M2L(ind1,ind2)
                if (n1p/=l2) then
                    trans_M2L(ind2,ind1)=1.0_prec/(-1.0_prec)**(n1p+l2)*trans_M2L(ind1,ind2)
                    trans_M2L(ind02,ind01)=trans_M2L(ind2,ind1)
                endif
            enddo
        enddo
!        write(*,*) 'pass 3'

        ! on avance en n2
        do m1=0,P1-1
            do n1=m1,2*P1-3-m1
                do l2=n1+1,2*P1-3-n1
                    l2p=l2+1
                    l2m=l2-1
                    n1m=n1-1
                    n1p=n1+1
                    ind1=n1p**2+n1p+1+m1
                    ind2=l2**2+l2+1+m1
                    a1=sqrt(real(n1+1+m1, prec)*real(n1+1-m1, prec)/real(2*n1+1, prec)/real(2*n1+3, prec))
                    a2=sqrt(real(l2m+1+m1, prec)*real(l2m+1-m1, prec)/real(2*l2m+1, prec)/real(2*l2m+3, prec))
                    a3=sqrt(real(l2+1+m1, prec)*real(l2+1-m1, prec)/real(2*l2+1, prec)/real(2*l2+3, prec))
                    if (n1-1<m1) then
                        a4=0.0_prec
                    else
                        a4=sqrt(real(n1m+1+m1, prec)*real(n1m+1-m1, prec)/real(2*n1m+1, prec)/real(2*n1m+3, prec))
                    endif
                    ind3=n1**2+n1+1+m1
                    ind4=l2m**2+l2m+1+m1
                    ind5=n1**2+n1+1+m1
                    ind6=l2p**2+l2p+1+m1
                    ind7=n1m**2+n1m+1+m1
                    ind8=l2**2+l2+1+m1
                    E2_M2M=trans_M2L(ind3,ind4)
                    E3_M2M=trans_M2L(ind5,ind6)
                    E4_M2M=trans_M2L(ind7,ind8)
                    trans_M2L(ind1,ind2)=1.0_prec/a1*(a2*E2_M2M-a3*E3_M2M+a4*E4_M2M)

                    ! symétrie
                    if (m1>0) then
                        ind01=n1p**2+n1p+1-m1
                        ind02=l2**2+l2+1-m1
                        trans_M2L(ind01,ind02)=trans_M2L(ind1,ind2);
                    endif
                    if (ind1/=ind2) then
                       trans_M2L(ind2,ind1)=1.0_prec/(-1.0_prec)**(n1p+l2)*trans_M2L(ind1,ind2)
                       if (m1>0) then
                        trans_M2L(ind02,ind01)=trans_M2L(ind2,ind1)
                       endif
                    endif
                enddo
            enddo
        enddo

        numM=0
        numL=0
        T0_M2M=cmplx(0.0_prec,0.0_prec, kind=prec)
        T0_L2L=cmplx(0.0_prec,0.0_prec, kind=prec)
        do n1=0,N11
            do m1=-n1,n1
                ind1=n1**2+n1+1+m1  ! indice niveau gro
                do l2=abs(m1),N11
                            ind2=l2**2+l2+1+m1   ! indice niveau fin
                            if (ind2 <= (N(i1)+1)**2) then
                                numL=numL+1
                                numM=numM+1
                                T0_L2L(numL)=trans_M2L(ind2,ind1)
                                T0_L2L_i(numL)=ind2
                                T0_L2L_j(numL)=ind1
                                T0_M2M(numM)=trans_M2L(ind2,ind1)
                                T0_M2M_i(numM)=ind1
                                T0_M2M_j(numM)=ind2
                            end if

                enddo
            enddo
        enddo

        if (allocated(SPH(i1-1)%trans_L2L)) deallocate(SPH(i1-1)%trans_L2L)
        allocate(SPH(i1-1)%trans_L2L(numL))
        if (allocated(SPH(i1-1)%trans_L2L_i)) deallocate(SPH(i1-1)%trans_L2L_i)
        allocate(SPH(i1-1)%trans_L2L_i(numL))
        if (allocated(SPH(i1-1)%trans_L2L_j)) deallocate(SPH(i1-1)%trans_L2L_j)
        allocate(SPH(i1-1)%trans_L2L_j(numL))
        if (allocated(SPH(i1-1)%trans_M2M)) deallocate(SPH(i1-1)%trans_M2M)
        allocate(SPH(i1)%trans_M2M(numM))
        if (allocated(SPH(i1-1)%trans_M2M_i)) deallocate(SPH(i1-1)%trans_M2M_i)
        allocate(SPH(i1)%trans_M2M_i(numM))
        if (allocated(SPH(i1-1)%trans_M2M_j)) deallocate(SPH(i1-1)%trans_M2M_j)
        allocate(SPH(i1)%trans_M2M_j(numM))

        SPH(i1-1)%trans_L2L=T0_L2L(1:numL)
        SPH(i1-1)%trans_L2L_i=T0_L2L_i(1:numL)
        SPH(i1-1)%trans_L2L_j=T0_L2L_j(1:numL)
        SPH(i1)%trans_M2M=T0_M2M(1:numM)
        SPH(i1)%trans_M2M_i=T0_M2M_i(1:numM)
        SPH(i1)%trans_M2M_j=T0_M2M_j(1:numM)

    end subroutine M2M_BF
!
! -------------------------------------------------------------------
!
! Rotation pour coaxiale translation
! La direction doit être de norme 1
! avec nouveau calcul de Ynm_gum_1 total
! -------------------------------------------------------------------
	subroutine ROT_BF(N11,dir,Rot,ix,iy,ind)
        implicit none
! variables d'entrée
        integer, intent(in) :: N11
        real(prec), intent(in) :: dir(3)
! variable de sortie
        complex(prec), intent(out),dimension(:),allocatable :: Rot
        integer, intent(out),dimension(:),allocatable :: ix,iy
        integer, intent(out) :: ind
! variables locales
        integer :: P1
        integer :: n2,m2,ind1,ind2,n1,m1
        integer :: n2m,m2p,m1m,m1p,m0
        integer :: ind3,ind4,ind5,ind6,ind7,ind8
        real(prec) :: H0M2M,HM2M((2*N11+1)**2,(2*N11+1)**2)
        real(prec) :: beta_M2M,alpha_M2M,gamma_M2M
        real(prec) :: b1,b2,b3,a1,Ps(3)
        complex(prec) :: Rot_M2M((N11+1)**2,(N11+1)**2)
        complex(prec) :: Ynm((2*N11+1)**2)
        ! DEBUG : character(len=1) :: answer


        ! angles d'Euler
        call car2sph(dir,Ps,1)
        beta_M2M=Ps(2)
        alpha_M2M=Ps(3)
        gamma_M2M=-pi/2.0_prec
        if (abs(beta_M2M-pi)<1e-12 .or. abs(beta_M2M) < 1e-12) then
            alpha_M2M=0.0_prec
            gamma_M2M=0.0_prec
        endif

	P1=N11+1
        call Ynm_gum(2*N11,beta_M2M,0.0_prec,Ynm)
        do n2=0,2*P1-2
            do m1=-n2,n2
                ind1=n2**2+n2+1+m1
                ind2=n2**2+n2+1
                H0M2M=dreal(sqrt(4.0_prec*pi/real(2*n2+1, prec))*Ynm(ind1))
                HM2M(ind1,ind2)=H0M2M;
                if (ind1 <= (N11+1)**2 .and. ind2<=(N11+1)**2) then
                    Rot_M2M(ind1,ind2)=exp(-ji*real(m1, prec)*gamma_M2M)*H0M2M
                endif
            enddo
        enddo

        do m2=0,N11-1
            m2p=m2+1
            do n2=2,2*P1-2-m2
                if (n2-1>=m2p) then

                    n2m=n2-1;
                    do m1=-n2m,n2m
                        m1m=m1-1
                        m1p=m1+1
                        m0=-m1-1
                        b1=sqrt(real(n2-m2-1, prec)*real(n2-m2, prec)/real(2*n2-1, prec)/real(2*n2+1, prec))
                        b2=sign(1,2*m0+1)*sqrt(real(n2-m0-1, prec)*real(n2-m0, prec)/real(2*n2-1, prec)/real(2*n2+1, prec))
                        b3=sign(1,2*m1m+1)*sqrt(real(n2-m1m-1, prec)*real(n2-m1m, prec)/real(2*n2-1, prec)/real(2*n2+1, prec))
                        a1=sqrt(real(n2m+m1+1, prec)*real(n2m+1-m1, prec)/real(2*n2m+1, prec)/real(2*n2m+3, prec))

                        ind1=n2m**2+n2m+1+m1
                        ind2=n2m**2+n2m+1+m2p
                        ind3=n2**2+n2+1+m1p
                        ind4=n2**2+n2+1+m2
                        ind5=n2**2+n2+1+m1m
                        ind6=n2**2+n2+1+m2
                        ind7=n2**2+n2+1+m1
                        ind8=n2**2+n2+1+m2
                        H0M2M=1/b1 * (0.5_prec * (b2*(1.0_prec-cos(beta_M2M))*HM2M(ind3,ind4) - &
            b3*(1.0_prec+cos(beta_M2M))*HM2M(ind5,ind6)) - a1 * sin(beta_M2M) * HM2M(ind7,ind8) )
                        HM2M(ind1,ind2)=H0M2M
                        if (ind1 <= (N11+1)**2 .and. ind2<=(N11+1)**2) then
                            Rot_M2M(ind1,ind2)=exp(ji*real(m2p, prec)*alpha_M2M)*exp(-ji*real(m1, prec)*gamma_M2M)*H0M2M
                        endif

                    enddo
                endif
            enddo
        enddo
        !% size(Rot_M2M)
        do n2=2,N11+1
            n2m=n2-1
            do m1=-n2m,n2m

                !% Rot_M2M(1,1)
                do m2=-n2m,-1
                    ind1=n2m**2+n2m+1+m1
                    ind2=n2m**2+n2m+1+m2
                    ind3=n2m**2+n2m+1-m1
                    ind4=n2m**2+n2m+1-m2
                    if (ind1 <= (N11+1)**2 .and. ind2<=(N11+1)**2) then

                        Rot_M2M(ind1,ind2)=conjg(Rot_M2M(ind3,ind4))

                    endif
                enddo
            enddo
        enddo
        ! création de la matrice creuse tri-bandes
        ind=0
        do n1=0,N11
            ind=ind+(2*n1+1)**2
        enddo
        allocate(rot(ind))
        allocate(ix(ind))
        allocate(iy(ind))
        ind=0
        do n1=0,N11
            do m1=-n1,n1
                ind1=n1**2+n1+1+m1
                n2=n1 ! classement par ordre
                    do m2=-n2,n2
                        ind2=n2**2+n2+1+m2
                        ind=ind+1
                        rot(ind)=Rot_M2M(ind1,ind2)
                        ix(ind)=ind1
                        iy(ind)=ind2
                    enddo
            enddo
        enddo

	end subroutine ROT_BF
!
! -----------------------------------------------------------------------------------------------------------------------------
! calcul des moments version diagonale
! -----------------------------------------------------------------------------------------------------------------------------
	subroutine moments_FMM_HF(octree,NL,N,SPH,k)
	implicit none
! variables d'entrée
	integer, intent(in) :: NL
	integer, dimension(1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	type(FMM),intent(inout) :: SPH(NL)
	type(tree) :: octree(NL)
! variables de sortie!
!	complex(prec), dimension(:), allocatable, intent(inout) :: W1,W2
!	integer, dimension(:), allocatable, intent(inout) :: iw,jw
!	integer, intent(out) :: nlw,nnzw
! variables locales
	integer :: i1,i2,no,nu,j0,num
	real(prec), dimension(1:3,1:IS_f) :: PG_f,Rmq,Q
	real(prec) :: Nor(3)
	real(prec), dimension(1:3,1:IS_f) :: F_f, dF1_f, dF2_f
	real(prec) :: PG2(2,IS_f)
	real(prec), dimension(1:3) :: M
	real(prec), dimension(1:3,1:3) :: E1
	real(prec), dimension(1:size(SPH(NL)%Ps,1),1:IS_f) :: rqs,rns
	real(prec), dimension(1:size(SPH(NL)%Ps,1)) :: rns0
	real(prec), dimension(1:size(SPH(NL)%Ps,1),1:3) :: Pu
!	complex(prec), dimension(size(SPH(NL)%Ps,1)*coef_2*NE) :: Dq,Mq
!	integer, dimension(size(SPH(NL)%Ps,1)*coef_2*NE) :: ii,jj
!	character(len=1) :: answer

!	write(*,*) 'passe moment'
	call pgt6(PG_f,IS_f)  ! points dans l'élément
	PG2=PG_f(1:2,:)
	call FORM_T3(PG2,IS_f,F_f,dF1_f,dF2_f)

	nu=(N(NL)+1)*(2*N(NL)+1) !% nombre de points sur la sphere
	Pu=SPH(NL)%Ps !% points sur la sphère unitaire

	! allocations
	allocate(SPH(1)%V1(1:NE*nu))
	allocate(SPH(1)%V2(1:NE*nu))
	allocate(SPH(1)%iv(1:NE*nu))
	allocate(SPH(1)%jv(1:NE*nu))

	num=0
	do i1=1,octree(NL)%NL 	! boucle cellules level le plus fin
	    M=octree(NL)%cell(i1)%center     ! coordonnées du centre de la cellule

	    do i2=1,octree(NL)%cell(i1)%Nb_c
	        no=octree(NL)%cell(i1)%child(i2)  !% n° du point central
        	E1=transpose(Ps(E(no,:),:))    !% fournit les coordonnées des sommets du triangle no
        	Q=matmul(E1,F_f)  !% coordonnées réel des points de Gauss
        	!% vecteurs
        	Rmq=spread(M,2,IS_f)-Q !% vecteurs sources et centre M
        	!% produit scalaire avec points sur la sphère untiaire
        	rqs=0.0_prec
        	call dgemm('N','N',size(SPH(NL)%Ps,1),IS_f,3,1.0_prec,SPH(NL)%Ps,size(SPH(NL)%Ps,1),Rmq,3,0.0_prec,rqs,&
        	size(SPH(NL)%Ps,1))
        	Nor=norc(no,:)   ! sans jacobien
        	rns0=matmul(SPH(NL)%Ps,Nor)
         	rns=spread(rns0,2,IS_f)

         	num=num+1
            SPH(1)%V2((num-1)*nu+1:num*nu)=matmul(exp(ji*k*rqs)*rns*Jac0(no),PG_f(3,:))
        	SPH(1)%V1((num-1)*nu+1:num*nu)=matmul(exp(ji*k*rqs)*Jac0(no),PG_f(3,:))

        	SPH(1)%iv((num-1)*nu+1:num*nu)=[((i1-1)*nu+j0,j0=1,nu)]
        	SPH(1)%jv((num-1)*nu+1:num*nu)=no

	    end do
	end do
!	SPH(1)%V1=Mq(1:num*nu)
!	SPH(1)%V2=Dq(1:num*nu)
!	SPH(1)%iv=ii(1:num*nu)
!	SPH(1)%jv=jj(1:num*nu)
	SPH(1)%nnzv=NE*nu
	SPH(1)%nlv=nu*octree(NL)%NL
	SPH(1)%ncv=NE
    allocate(SPH(NL)%in(SPH(1)%nlv))

	end subroutine moments_FMM_HF
!
! -----------------------------------------------------------------------------------------------------------------------------
! calcul des moments version diagonale
! version propa octree augmentée des Pexts
! -----------------------------------------------------------------------------------------------------------------------------
	subroutine moments_FMM_propa_HF(octree,NL,N,SPH,k)
	implicit none
! variables d'entrée
	integer, intent(in) :: NL
	integer, dimension(1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	type(FMM),intent(inout) :: SPH(NL)
	type(tree) :: octree(NL)
! variables de sortie!
!	complex(prec), dimension(:), allocatable, intent(inout) :: W1,W2
!	integer, dimension(:), allocatable, intent(inout) :: iw,jw
!	integer, intent(out) :: nlw,nnzw
! variables locales
	integer :: i1,i2,no,nu,j0,num
	real(prec), dimension(1:3,1:IS_f) :: PG_f,Rmq,Q
	real(prec) :: Nor(3)
	real(prec), dimension(1:3,1:IS_f) :: F_f, dF1_f, dF2_f
	real(prec) :: PG2(2,IS_f)
	real(prec), dimension(1:3) :: M
	real(prec), dimension(1:3,1:3) :: E1
	real(prec), dimension(1:size(SPH(NL)%Ps,1),1:IS_f) :: rqs,rns
	real(prec), dimension(1:size(SPH(NL)%Ps,1)) :: rns0
	real(prec), dimension(1:size(SPH(NL)%Ps,1),1:3) :: Pu
!	complex(prec), dimension(size(SPH(NL)%Ps,1)*coef_2*NE) :: Dq,Mq
!	integer, dimension(size(SPH(NL)%Ps,1)*coef_2*NE) :: ii,jj
!	character(len=1) :: answer

!	write(*,*) 'passe moment'
	call pgt6(PG_f,IS_f)  ! points dans l'élément
	PG2=PG_f(1:2,:)
	call FORM_T3(PG2,IS_f,F_f,dF1_f,dF2_f)

	nu=(N(NL)+1)*(2*N(NL)+1) !% nombre de points sur la sphere
	Pu=SPH(NL)%Ps !% points sur la sphère unitaire

	! allocations
	allocate(SPH(1)%V1(1:(NE+Next)*nu))
	allocate(SPH(1)%V2(1:(NE+Next)*nu))
	allocate(SPH(1)%iv(1:(NE+Next)*nu))
	allocate(SPH(1)%jv(1:(NE+Next)*nu))

	num=0
	do i1=1,octree(NL)%NL 	! boucle cellules level le plus fin
	    M=octree(NL)%cell(i1)%center     ! coordonnées du centre de la cellule

	    do i2=1,octree(NL)%cell(i1)%Nb_c
	        no=octree(NL)%cell(i1)%child(i2)  !% n° du point central
			num=num+1
			if (no <= NE) then
				E1=transpose(Ps(E(no,:),:))    !% fournit les coordonnées des sommets du triangle no
				Q=matmul(E1,F_f)  !% coordonnées réel des points de Gauss
				!% vecteurs
				Rmq=spread(M,2,IS_f)-Q !% vecteurs sources et centre M
				!% produit scalaire avec points sur la sphère untiaire
				rqs=0.0_prec
				call dgemm('N','N',size(SPH(NL)%Ps,1),IS_f,3,1.0_prec,SPH(NL)%Ps,size(SPH(NL)%Ps,1),Rmq,3,0.0_prec,rqs,&
				size(SPH(NL)%Ps,1))
				Nor=norc(no,:)   ! sans jacobien
				rns0=matmul(SPH(NL)%Ps,Nor)
				rns=spread(rns0,2,IS_f)

				SPH(1)%V2((num-1)*nu+1:num*nu)=matmul(exp(ji*k*rqs)*rns*Jac0(no),PG_f(3,:))
				SPH(1)%V1((num-1)*nu+1:num*nu)=matmul(exp(ji*k*rqs)*Jac0(no),PG_f(3,:))

				SPH(1)%iv((num-1)*nu+1:num*nu)=[((i1-1)*nu+j0,j0=1,nu)]
				SPH(1)%jv((num-1)*nu+1:num*nu)=no
			else
				SPH(1)%V2((num-1)*nu+1:num*nu)=cmplx(0.0_prec,0.0_prec, kind=prec)
				SPH(1)%V1((num-1)*nu+1:num*nu)=cmplx(0.0_prec,0.0_prec, kind=prec)

				SPH(1)%iv((num-1)*nu+1:num*nu)=[((i1-1)*nu+j0,j0=1,nu)]
				SPH(1)%jv((num-1)*nu+1:num*nu)=no

			endif
		end do
	end do
!	SPH(1)%V1=Mq(1:num*nu)
!	SPH(1)%V2=Dq(1:num*nu)
!	SPH(1)%iv=ii(1:num*nu)
!	SPH(1)%jv=jj(1:num*nu)
	SPH(1)%nnzv=(NE+Next)*nu
	SPH(1)%nlv=nu*octree(NL)%NL
	SPH(1)%ncv=NE+Next
    allocate(SPH(NL)%in(SPH(1)%nlv))

	end subroutine moments_FMM_propa_HF
! -----------------------------------------------------------------------------------------------------------------------------
! calcul des moments init_param_FMM version BF
! -----------------------------------------------------------------------------------------------------------------------------
	subroutine moments_FMM_BF(octree,NL,N,SPH,SPH0,k)
	implicit none
! variables d'entrée
	integer, intent(in) :: NL
	integer, dimension(1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	type(FMM),intent(inout) :: SPH(NL)
	type(FMM0), intent(inout) :: SPH0
	type(tree) :: octree(NL)

! variables locales
	integer :: i1,i2,i4,no,nanm,NM,num
	integer :: m2,n2,n2p,n2m,m0,m0m,m0p,m11
	integer :: acn_nmm,acn_nmmm,acn_nmmp,acn_npm,acn_npmm,acn_npmp,acn_nm
	!real(prec) :: b1,b2,b3,b4,a01,a02
	real(prec), dimension(2,IS_f) :: PG2
	real(prec), dimension(3,IS_f) :: PG_f
	real(prec), dimension(IS_f) :: theta,phi
	real(prec) :: Rmq(3,IS_f),Rmq2(IS_f,3),Pmq(IS_f,3),M(3)
	real(prec), dimension(3,IS_f) :: F_f,dF1_f,dF2_f,Q
	real(prec), dimension(3,3) :: E1
    real(prec) :: jn(0:N(NL)+1),Djn(0:N(NL)+1)
    real(prec) :: jn2(0:N(NL)+1,IS_f)
	complex(prec) ,dimension((N(NL)+2)**2,IS_f) :: Ynm
	complex(prec),dimension(IS_f) :: dx,dy,dz
	complex(prec),dimension(IS_f) :: jY_1,jY_2,jY_3,jY_4,jY_5,jY_6
!	character(len=1) :: answer

	call pgt6(PG_f,IS_f)  ! points dans l'élément
	PG2=PG_f(1:2,:)
	call FORM_T3(PG2,IS_f,F_f,dF1_f,dF2_f)

	nanm=(N(NL)+1)**2 !% nombre de composante HS
	allocate(SPH(1)%V1(nanm*NE))
	allocate(SPH(1)%V2(nanm*NE))
	allocate(SPH(1)%iv(nanm*NE))
	allocate(SPH(1)%jv(nanm*NE))

	num=0
	do i1=1,octree(NL)%NL 	! boucle cellules level le plus fin
		M=octree(NL)%cell(i1)%center     ! coordonnées du centre de la cellule
		do i2=1,octree(NL)%cell(i1)%Nb_c
			no=octree(NL)%cell(i1)%child(i2)  !% n° du point central
			E1=transpose(Ps(E(no,:),:))    !% fournit les coordonnées des sommets du triangle no
			Q=matmul(E1,F_f)  !% coordonnées réel des points de Gauss
        	num=num+1
			Rmq=Q-spread(M,2,IS_f)
			Rmq2=transpose(Rmq)
			call car2sph(Rmq2,Pmq,IS_f)
			theta=Pmq(:,2)
			phi=Pmq(:,3)
			do i4=1,IS_f
				call SPHJ(N(NL)+1,k*Pmq(i4,1),NM,jn,Djn)
				jn2(:,i4)=jn
			enddo
			call Ynm_gum(N(NL)+1,theta,phi,IS_f,Ynm)

			do i4=1,(N(NL)+1)**2
				n2=int(sqrt(real(i4-1)))
				m2=i4-n2**2-n2-1
				n2p=n2+1
				n2m=n2-1
				m0=-m2
				m0p=m0+1
				m0m=m0-1
				m11=-m0-1

				acn_npmp=n2p**2+n2p+m0p+1
				acn_npmm=n2p**2+n2p+m0m+1
				acn_nmmp=n2m**2+n2m+m0p+1
				acn_nmmm=n2m**2+n2m+m0m+1
				acn_nmm =n2m**2+n2m+m0+1
				acn_npm =n2p**2+n2p+m0+1
				acn_nm= n2**2+n2+m0+1

				jY_1=SPH0%b(1,i4) * jn2(n2p,:) * (Ynm(acn_npmp,:))
				jY_3=SPH0%b(3,i4) * jn2(n2p,:) * (Ynm(acn_npmm,:))
				jY_6=SPH0%a(2,i4) * jn2(n2p,:) * (Ynm(acn_npm,:))

				if (abs(SPH0%b(2,i4)) > 1.0e-5_prec) then
					jY_2=SPH0%b(2,i4) * jn2(n2m,:) * (Ynm(acn_nmmp,:))
				else
					jY_2=cmplx(0.0_prec,0.0_prec, kind=prec)
				endif
				if (abs(SPH0%b(4,i4)) > 1.0e-5_prec) then
					jY_4=SPH0%b(4,i4) * jn2(n2m,:) * (Ynm(acn_nmmm,:))
				else
					jY_4=cmplx(0.0_prec,0.0_prec, kind=prec)
				endif
				if (abs(SPH0%a(1,i4)) > 1.0e-5_prec) then
					jY_5=SPH0%a(1,i4) * jn2(n2m,:) * (Ynm(acn_nmm,:))
				else
					jY_5=cmplx(0.0_prec,0.0_prec, kind=prec)
				endif
				dx=k/2.0_prec * (jY_1-jY_2+jY_3-jY_4)
				dy=-ji*k/2.0_prec * (jY_1-jY_2-jY_3+jY_4)
				dz=k * (jY_5 - jY_6)
				! moment G
				SPH(1)%V1((num-1)*nanm+i4)= sum(jn2(n2,:)*(Ynm(acn_nm,:)) * Jac0(no) * PG_f(3,:))  ! le conjugué est dans le m0=-m2
				! moment nq.grad(G)
				! le jacobien est dans la normale
				SPH(1)%V2((num-1)*nanm+i4)= sum( &
					(norc(no,1)*dx+norc(no,2)*dy+norc(no,3)*dz) * Jac0(no) * PG_f(3,:))
				SPH(1)%iv((num-1)*nanm+i4)=(i1-1)*nanm+i4
				SPH(1)%jv((num-1)*nanm+i4)=no

			enddo

	    end do
	end do

	SPH(1)%nnzv=NE*nanm
	SPH(1)%nlv=nanm*octree(NL)%NL
	SPH(1)%ncv=NE
    allocate(SPH(NL)%inBF(SPH(1)%nlv))

	end subroutine moments_FMM_BF
! -----------------------------------------------------------------------------------------------------------------------------
! calcul des moments init_param_FMM version BF
! version propa pour octree augmenté avec Pext
! -----------------------------------------------------------------------------------------------------------------------------
	subroutine moments_FMM_propa_BF(octree,NL,N,SPH,SPH0,k)
	implicit none
! variables d'entrée
	integer :: NL
	integer, dimension(1:NL) :: N
	real(prec) :: k
	type(FMM) :: SPH(NL)
	type(FMM0) :: SPH0
	type(tree) :: octree(NL)

! variables locales
	integer :: i1,i2,i4,no,nanm,NM,num,j0
	integer :: m2,n2,n2p,n2m,m0,m0m,m0p,m11
	integer :: acn_nmm,acn_nmmm,acn_nmmp,acn_npm,acn_npmm,acn_npmp,acn_nm
	!real(prec) :: b1,b2,b3,b4,a01,a02
	real(prec), dimension(2,IS_f) :: PG2
	real(prec), dimension(3,IS_f) :: PG_f
	real(prec), dimension(IS_f) :: theta,phi
	real(prec) :: Rmq(3,IS_f),Rmq2(IS_f,3),Pmq(IS_f,3),M(3)
	real(prec), dimension(3,IS_f) :: F_f,dF1_f,dF2_f,Q
	real(prec), dimension(3,3) :: E1
    real(prec) :: jn(0:N(NL)+1),Djn(0:N(NL)+1)
    real(prec) :: jn2(0:N(NL)+1,IS_f)
	complex(prec) ,dimension((N(NL)+2)**2,IS_f) :: Ynm
	complex(prec),dimension(IS_f) :: dx,dy,dz
	complex(prec),dimension(IS_f) :: jY_1,jY_2,jY_3,jY_4,jY_5,jY_6
!	character(len=1) :: answer

!	write(*,*) 'pass in moment_propa_BF'
	call pgt6(PG_f,IS_f)  ! points dans l'élément
	PG2=PG_f(1:2,:)
	call FORM_T3(PG2,IS_f,F_f,dF1_f,dF2_f)

	nanm=(N(NL)+1)**2 !% nombre de composante HS
	if (allocated(SPH(1)%V1)) deallocate(SPH(1)%V1)
	if (allocated(SPH(1)%V2)) deallocate(SPH(1)%V2)
	if (allocated(SPH(1)%iv)) deallocate(SPH(1)%iv)
	if (allocated(SPH(1)%jv)) deallocate(SPH(1)%jv)
	allocate(SPH(1)%V1(nanm*(NE+Next)))
	allocate(SPH(1)%V2(nanm*(NE+Next)))
	allocate(SPH(1)%iv(nanm*(NE+Next)))
	allocate(SPH(1)%jv(nanm*(NE+Next)))

	num=0
	do i1=1,octree(NL)%NL 	! boucle cellules level le plus fin
		M=octree(NL)%cell(i1)%center     ! coordonnées du centre de la cellule
		do i2=1,octree(NL)%cell(i1)%Nb_c
			no=octree(NL)%cell(i1)%child(i2)  !% n° du point central
			num=num+1
			if (no<=NE) then
				E1=transpose(Ps(E(no,:),:))    !% fournit les coordonnées des sommets du triangle no
				Q=matmul(E1,F_f)  !% coordonnées réel des points de Gauss
				Rmq=Q-spread(M,2,IS_f)
				Rmq2=transpose(Rmq)
				call car2sph(Rmq2,Pmq,IS_f)
				theta=Pmq(:,2)
				phi=Pmq(:,3)
				do i4=1,IS_f
					call SPHJ(N(NL)+1,k*Pmq(i4,1),NM,jn,Djn)
					jn2(:,i4)=jn
				enddo
				!write(*,*) 'NL=',NL
				!write(*,*) 'N(',NL,')+1=',N(NL)+1,' IS_f=',IS_f
				!write(*,*) 'num=',num,' no=',no
				call Ynm_gum(N(NL)+1,theta,phi,IS_f,Ynm)

				do i4=1,(N(NL)+1)**2
					n2=int(sqrt(real(i4-1)))
					m2=i4-n2**2-n2-1
					n2p=n2+1
					n2m=n2-1
					m0=-m2
					m0p=m0+1
					m0m=m0-1
					m11=-m0-1

					acn_npmp=n2p**2+n2p+m0p+1
					acn_npmm=n2p**2+n2p+m0m+1
					acn_nmmp=n2m**2+n2m+m0p+1
					acn_nmmm=n2m**2+n2m+m0m+1
					acn_nmm =n2m**2+n2m+m0+1
					acn_npm =n2p**2+n2p+m0+1
					acn_nm= n2**2+n2+m0+1

					jY_1=SPH0%b(1,i4) * jn2(n2p,:) * (Ynm(acn_npmp,:))
					jY_3=SPH0%b(3,i4) * jn2(n2p,:) * (Ynm(acn_npmm,:))
					jY_6=SPH0%a(2,i4) * jn2(n2p,:) * (Ynm(acn_npm,:))

					if (abs(SPH0%b(2,i4)) > 1.0e-5_prec) then
						jY_2=SPH0%b(2,i4) * jn2(n2m,:) * (Ynm(acn_nmmp,:))
					else
						jY_2=cmplx(0.0_prec,0.0_prec, kind=prec)
					endif
					if (abs(SPH0%b(4,i4)) > 1.0e-5_prec) then
						jY_4=SPH0%b(4,i4) * jn2(n2m,:) * (Ynm(acn_nmmm,:))
					else
						jY_4=cmplx(0.0_prec,0.0_prec, kind=prec)
					endif
					if (abs(SPH0%a(1,i4)) > 1.0e-5_prec) then
						jY_5=SPH0%a(1,i4) * jn2(n2m,:) * (Ynm(acn_nmm,:))
					else
						jY_5=cmplx(0.0_prec,0.0_prec, kind=prec)
					endif
					dx=k/2.0_prec * (jY_1-jY_2+jY_3-jY_4)
					dy=-ji*k/2.0_prec * (jY_1-jY_2-jY_3+jY_4)
					dz=k * (jY_5 - jY_6)
					! moment G
					SPH(1)%V1((num-1)*nanm+i4)= sum(jn2(n2,:)*(Ynm(acn_nm,:)) * Jac0(no) * PG_f(3,:))  ! le conjugué est dans le m0=-m2
					! moment nq.grad(G)
					! le jacobien est dans la normale
					SPH(1)%V2((num-1)*nanm+i4)= sum( &
						(norc(no,1)*dx+norc(no,2)*dy+norc(no,3)*dz) * Jac0(no) * PG_f(3,:))
					SPH(1)%iv((num-1)*nanm+i4)=(i1-1)*nanm+i4
					SPH(1)%jv((num-1)*nanm+i4)=no

				enddo
			else
					SPH(1)%V1((num-1)*nanm+1:num*nanm)= cmplx(0.0_prec,0.0_prec, kind=prec)  ! Pour contribution Pext ==0
					! moment nq.grad(G)
					SPH(1)%V2((num-1)*nanm+1:num*nanm)= cmplx(0.0_prec,0.0_prec, kind=prec)       ! Pour contribution Pext==0
					SPH(1)%iv((num-1)*nanm+1:num*nanm)=[((i1-1)*nanm+j0,j0=1,nanm)]
					SPH(1)%jv((num-1)*nanm+1:num*nanm)=no

			endif
	    end do
	end do

!	write(*,*) 'num=',num, 'NE+Next=',NE+Next

	SPH(1)%nnzv=(NE+Next)*nanm
	SPH(1)%nlv=nanm*octree(NL)%NL
	SPH(1)%ncv=NE+Next
	if (allocated(SPH(NL)%inBF)) deallocate(SPH(NL)%inBF)
    allocate(SPH(NL)%inBF(SPH(1)%nlv))

	end subroutine moments_FMM_propa_BF
!
! -----------------------------------------------------------------------------------------------------------------------------
! calcul des moments version diagonale
! coefficients séparés
! -----------------------------------------------------------------------------------------------------------------------------
	subroutine expansions_FMM_HF(octree,NL,N,SPH,k)
	implicit none
! variables d'entrée
	integer, intent(in) :: NL
	integer, dimension(1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	type(FMM),intent(inout) :: SPH(NL)
	type(tree) :: octree(NL)
! variables de sortie!
!	complex(prec), dimension(:), allocatable, intent(inout) :: W1,W2
!	integer, dimension(:), allocatable, intent(inout) :: iw,jw
!	integer, intent(out) :: ncw,nnzw
! variables locales
	integer :: i1,i2,no,nu,j0,num
	real(prec), dimension(1:3) :: Rlp
	real(prec), dimension(1:3) :: M
	real(prec), dimension(1:size(SPH(NL)%Ps,1)) :: rps,rns
	real(prec), dimension(1:size(SPH(NL)%Ps,1),1:3) :: Pu
!	complex(prec), dimension(size(SPH(NL)%Ps,1)*coef_2*NE) :: Dq,Mq
!	integer, dimension(size(SPH(NL)%Ps,1)*coef_2*NE) :: ii,jj
!	character(len=1) :: answer

!	write(*,*) 'passe expansion'
	nu=(N(NL)+1)*(2*N(NL)+1) !% nombre de points sur la sphere
	Pu=SPH(NL)%Ps !% points sur la sphère unitaire

	! allocations
	allocate(SPH(1)%W1(1:NE*nu))
	allocate(SPH(1)%W2(1:NE*nu))
	allocate(SPH(1)%iw(1:NE*nu))
	allocate(SPH(1)%jw(1:NE*nu))

	num=0
	do i1=1,octree(NL)%NL 	! boucle cellules level le plus fin

        M=octree(NL)%cell(i1)%center     ! coordonnées du centre de la cellule

	    do i2=1,octree(NL)%cell(i1)%Nb_c
	        no=octree(NL)%cell(i1)%child(i2)  !% n° du point central
            Rlp=Pc(no,:)-M
	        rps=matmul(SPH(NL)%Ps,Rlp)
        	rns=matmul(SPH(NL)%Ps,norc(no,:))
        	num=num+1
        	SPH(1)%W1((num-1)*nu+1:num*nu)=exp(ji*k*rps)*SPH(NL)%wp
        	SPH(1)%W2((num-1)*nu+1:num*nu)=exp(ji*k*rps)*rns*SPH(NL)%wp
        	SPH(1)%jw((num-1)*nu+1:num*nu)=[((i1-1)*nu+j0,j0=1,nu)]
        	SPH(1)%iw((num-1)*nu+1:num*nu)=no
	    end do
	end do

!	SPH(1)%W1=Dq(1:num*nu)
!	SPH(1)%W2=Mq(1:num*nu)
!	SPH(1)%iw=ii(1:num*nu)
!	SPH(1)%jw=jj(1:num*nu)
	SPH(1)%nnzw=NE*nu
	SPH(1)%ncw=nu*octree(NL)%NL
	SPH(1)%nlw=NE

	end subroutine expansions_FMM_HF
! -----------------------------------------------------------------------------------------------------------------------------
! calcul des moments version diagonale
! version Burton/Miller regroupée, gain de mémoire
! -----------------------------------------------------------------------------------------------------------------------------
	subroutine expansions_BM_HF(octree,NL,N,SPH,k,alpha)
	implicit none
! variables d'entrée
	integer, intent(in) :: NL
	integer, dimension(1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	complex(prec), intent(in) :: alpha
	type(FMM),intent(inout) :: SPH(NL)
	type(tree) :: octree(NL)
! variables de sortie!
!	complex(prec), dimension(:), allocatable, intent(inout) :: W1,W2
!	integer, dimension(:), allocatable, intent(inout) :: iw,jw
!	integer, intent(out) :: ncw,nnzw
! variables locales
	integer :: i1,i2,no,nu,j0,num
	real(prec), dimension(1:3) :: Rlp
	real(prec), dimension(1:3) :: M
	real(prec), dimension(1:size(SPH(NL)%Ps,1)) :: rps,rns
	real(prec), dimension(1:size(SPH(NL)%Ps,1),1:3) :: Pu
!	complex(prec), dimension(size(SPH(NL)%Ps,1)*coef_2*NE) :: Dq,Mq
!	integer, dimension(size(SPH(NL)%Ps,1)*coef_2*NE) :: ii,jj
!	character(len=1) :: answer

	nu=(N(NL)+1)*(2*N(NL)+1) !% nombre de points sur la sphere
	Pu=SPH(NL)%Ps !% points sur la sphère unitaire

	! allocations
	allocate(SPH(1)%W1(1:NE*nu))
!	allocate(SPH(1)%W2(1:num*nu))
	allocate(SPH(1)%iw(1:NE*nu))
	allocate(SPH(1)%jw(1:NE*nu))


	num=0
	do i1=1,octree(NL)%NL 	! boucle cellules level le plus fin

        M=octree(NL)%cell(i1)%center     ! coordonnées du centre de la cellule

	    do i2=1,octree(NL)%cell(i1)%Nb_c
	        no=octree(NL)%cell(i1)%child(i2)  !% n° du point central
            Rlp=Pt(no,:)-M   ! cible (CHI : point interieur)
	        rps=matmul(SPH(NL)%Ps,Rlp)
        	rns=matmul(SPH(NL)%Ps,norc(no,:))
        	num=num+1
        	SPH(1)%W1((num-1)*nu+1:num*nu)=exp(ji*k*rps)*SPH(NL)%wp + &
				alpha * ji * k * exp(ji*k*rps)*rns*SPH(NL)%wp
!			SPH(1)%W1((num-1)*nu+1:num*nu)=ji * k * exp(ji*k*rps)*rns*SPH(NL)%wp
        	SPH(1)%jw((num-1)*nu+1:num*nu)=[((i1-1)*nu+j0,j0=1,nu)]
        	SPH(1)%iw((num-1)*nu+1:num*nu)=no
	    end do
	end do

!	SPH(1)%W1=Dq(1:num*nu)+ alpha * ji * k * Mq(1:num*nu)
!	SPH(1)%W2=Mq(1:num*nu)
!	SPH(1)%iw=ii(1:num*nu)
!	SPH(1)%jw=jj(1:num*nu)
	SPH(1)%nnzw=NE*nu
	SPH(1)%ncw=nu*octree(NL)%NL
	SPH(1)%nlw=NE

	end subroutine expansions_BM_HF
! -----------------------------------------------------------------------------------------------------------------------------
! calcul des moments version diagonale
! version Burton/Miller regroupée, gain de mémoire
!
! version pour propa avec octree augmenté avec Pext
!
! -----------------------------------------------------------------------------------------------------------------------------
	subroutine expansions_BM_propa_HF(octree,NL,N,SPH,k)
	implicit none
! variables d'entrée
	integer, intent(in) :: NL
	integer, dimension(1:NL), intent(in) :: N
	real(prec), intent(in) :: k
!	complex(prec), intent(in) :: alpha
	type(FMM),intent(inout) :: SPH(NL)
	type(tree) :: octree(NL)
! variables de sortie!
!	complex(prec), dimension(:), allocatable, intent(inout) :: W1,W2
!	integer, dimension(:), allocatable, intent(inout) :: iw,jw
!	integer, intent(out) :: ncw,nnzw
! variables locales
	integer :: i1,i2,no,nu,j0,num
	real(prec), dimension(1:3) :: Rlp
	real(prec), dimension(1:3) :: M
	real(prec), dimension(1:size(SPH(NL)%Ps,1)) :: rps
	real(prec), dimension(1:size(SPH(NL)%Ps,1),1:3) :: Pu
!	complex(prec), dimension(size(SPH(NL)%Ps,1)*coef_2*NE) :: Dq,Mq
!	integer, dimension(size(SPH(NL)%Ps,1)*coef_2*NE) :: ii,jj
!	character(len=1) :: answer

	nu=(N(NL)+1)*(2*N(NL)+1) !% nombre de points sur la sphere
	Pu=SPH(NL)%Ps !% points sur la sphère unitaire

	! allocations
	allocate(SPH(1)%W1(1:Next*nu))
!	allocate(SPH(1)%W2(1:num*nu))
	allocate(SPH(1)%iw(1:Next*nu))
	allocate(SPH(1)%jw(1:Next*nu))

	!write(*,*) 'in expansion_propa HF, NE=',NE
	num=0
	do i1=1,octree(NL)%NL 	! boucle cellules level le plus fin

        M=octree(NL)%cell(i1)%center     ! coordonnées du centre de la cellule
!		write(*,*) 'nchild=',octree(NL)%cell(i1)%Nb_c,' en i1=',i1
!		write(*,*) 'no=',octree(NL)%cell(i1)%child(1)
	    do i2=1,octree(NL)%cell(i1)%Nb_c
	        no=octree(NL)%cell(i1)%child(i2)  !% n° du point central
	        !write(*,*) 'no=',no,' norc=',norc(no,:)
	        if (no>NE) then
				Rlp=P_prop(no,:)-M
				rps=matmul(SPH(NL)%Ps,Rlp)
				! rns=matmul(SPH(NL)%Ps,norc(no,:)) ! pour l'instant pas de W2, donc pas de norc
				num=num+1
				SPH(1)%W1((num-1)*nu+1:num*nu)=exp(ji*k*rps)*SPH(NL)%wp
!				SPH(1)%W1((num-1)*nu+1:num*nu)=exp(ji*k*rps)*SPH(NL)%wp + &
!				alpha * ji * k * exp(ji*k*rps)*rns*SPH(NL)%wp
!			SPH(1)%W1((num-1)*nu+1:num*nu)=ji * k * exp(ji*k*rps)*rns*SPH(NL)%wp
				SPH(1)%jw((num-1)*nu+1:num*nu)=[((i1-1)*nu+j0,j0=1,nu)]
				SPH(1)%iw((num-1)*nu+1:num*nu)=no-NE
			endif
		end do
	end do


!	SPH(1)%W1=Dq(1:num*nu)+ alpha * ji * k * Mq(1:num*nu)
!	SPH(1)%W2=Mq(1:num*nu)
!	SPH(1)%iw=ii(1:num*nu)
!	SPH(1)%jw=jj(1:num*nu)
	SPH(1)%nnzw=num*nu
	SPH(1)%ncw=nu*octree(NL)%NL
	SPH(1)%nlw=num

	end subroutine expansions_BM_propa_HF
! -----------------------------------------------------------------------------------------------------------------------------
! calcul des moments version BF
! version Burton/Miller regroupée, gain de mémoire
! attention FAUX pour l'instant pour le d2, donc pour BM
! -----------------------------------------------------------------------------------------------------------------------------
	subroutine expansions_BM_BF(octree,NL,N,SPH,SPH0,k,alpha)
	implicit none
! variables d'entrée
	integer, intent(in) :: NL
	integer, dimension(1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	complex(prec), intent(in) :: alpha
	type(tree) :: octree(NL)
	type(FMM),intent(inout) :: SPH(NL)
	type(FMM0), intent(inout) :: SPH0
! variables de sortie!
! variables locales
	integer :: i1,i2,i4,no,n2,m2,N1,nanm,NM,num
	integer :: n2p,n2m,m0,m0m,m0p,m11
	integer :: acn_nmm,acn_nmmm,acn_nmmp,acn_npm,acn_npmm,acn_npmp
	real(prec), dimension(1:3) :: Rlp
	real(prec), dimension(1:3) :: M
!	real(prec) :: theta,phi
	complex(prec) :: Ynm((N(NL)+2)**2)
	complex(prec) :: dx,dy,dz
	complex(prec) :: jY_1,jY_2,jY_3,jY_4,jY_5,jY_6
	real(prec), dimension(3) :: Psph
!	character(len=1) :: answer
    real(prec) :: jn(0:N(NL)+1),Djn(0:N(NL)+1)  !,npm,nR

!
! précalcul
!
    N1=N(NL)
	nanm=(N1+1)**2

	allocate(SPH(1)%W1(1:NE*nanm))
	allocate(SPH(1)%iw(1:NE*nanm))
	allocate(SPH(1)%jw(1:NE*nanm))
	SPH(1)%nnzw=NE*nanm
	SPH(1)%ncw=nanm*octree(NL)%NL
	SPH(1)%nlw=NE

    num=0
	do i1=1,octree(NL)%NL 	! boucle cellules level le plus fin
	    M=octree(NL)%cell(i1)%center     ! coordonnées du centre de la cellule
	    do i2=1,octree(NL)%cell(i1)%Nb_c
	        no=octree(NL)%cell(i1)%child(i2)  !% n° du point central
	        num=num+1
            Rlp=Pt(no,:)-M   ! cible (CHI : point interieur)
            !nR=sqrt(sum(Rlp**2))
            !npm=Rlp(1)*norc(no,1)+Rlp(2)*norc(no,2)+Rlp(3)*norc(no,3)
            call car2sph(Rlp,Psph,1)
            call SPHJ(N1+1,k*Psph(1),NM,jn,Djn)
			call Ynm_gum(N(NL)+1,Psph(2),Psph(3),Ynm)
!			write(*,*) 'kr=',k*Psph(1),' theta=',Psph(2),' phi=',Psph(3)
			do i4=1,(N(NL)+1)**2
				n2=int(sqrt(real(i4-1)))
				m2=i4-n2**2-n2-1
				n2p=n2+1
				n2m=n2-1
				m0=m2
				m0p=m0+1
				m0m=m0-1
				m11=-m0-1

				acn_npmp=n2p**2+n2p+m0p+1
				acn_npmm=n2p**2+n2p+m0m+1
				acn_nmmp=n2m**2+n2m+m0p+1
				acn_nmmm=n2m**2+n2m+m0m+1
				acn_nmm =n2m**2+n2m+m0+1
				acn_npm =n2p**2+n2p+m0+1
!				acn_nm=n2**2+n2+m2+1

!				jY_1=SPH0%b(1,i4) * jn(n2p) * (Ynm(acn_npmp))
!				jY_3=SPH0%b(3,i4) * jn(n2p) * (Ynm(acn_npmm))
!				jY_6=SPH0%a(2,i4) * jn(n2p) * (Ynm(acn_npm))

				jY_1=SPH0%b_ex(1,i4) * jn(n2p) * (Ynm(acn_npmp))
				jY_3=SPH0%b_ex(3,i4) * jn(n2p) * (Ynm(acn_npmm))
				jY_6=SPH0%a(2,i4) * jn(n2p) * (Ynm(acn_npm))

				if (abs(SPH0%b_ex(2,i4)) > 1.0e-5_prec) then
					jY_2=SPH0%b_ex(2,i4) * jn(n2m) * (Ynm(acn_nmmp))
				else
					jY_2=cmplx(0.0_prec,0.0_prec, kind=prec)
				endif
				if (abs(SPH0%b_ex(4,i4)) > 1.0e-5_prec) then
					jY_4=SPH0%b_ex(4,i4) * jn(n2m) * (Ynm(acn_nmmm))
				else
					jY_4=cmplx(0.0_prec,0.0_prec, kind=prec)
				endif
				if (abs(SPH0%a(1,i4)) > 1.0e-5_prec) then
					jY_5=SPH0%a(1,i4) * jn(n2m) * (Ynm(acn_nmm))
				else
					jY_5=cmplx(0.0_prec,0.0_prec, kind=prec)
				endif
				dx=k/2.0_prec * (jY_1-jY_2+jY_3-jY_4)
				dy=-ji*k/2.0_prec * (jY_1-jY_2-jY_3+jY_4)
				dz=k * (jY_5 - jY_6)
				! moment G
				SPH(1)%W1((num-1)*nanm+i4)= 4.0_prec*pi*jn(n2)*Ynm(i4) + alpha * &
				! moment nq.grad(G)
					4.0_prec*pi*(norc(no,1)*dx+norc(no,2)*dy+norc(no,3)*dz)
				!SPH(1)%W1((num-1)*nanm+i4)=	4.0_prec*pi*(norc(no,1)*dx+norc(no,2)*dy+norc(no,3)*dz)
				SPH(1)%jw((num-1)*nanm+i4)=(i1-1)*nanm+i4
				SPH(1)%iw((num-1)*nanm+i4)=no

			enddo
	    end do
	end do

	end subroutine expansions_BM_BF
! -----------------------------------------------------------------------------------------------------------------------------
! calcul des moments version BF
! version Burton/Miller regroupée, gain de mémoire
! attention FAUX pour l'instant pour le d2, donc pour BM
!
! version augmenté pour octree + Pext
! -----------------------------------------------------------------------------------------------------------------------------
	subroutine expansions_BM_propa_BF(octree,NL,N,SPH,SPH0,k)
	implicit none
! variables d'entrée
	integer :: NL
	integer, dimension(1:NL) :: N
	real(prec) :: k
!	complex(prec), intent(in) :: alpha
	type(tree) :: octree(NL)
	type(FMM) :: SPH(NL)
	type(FMM0) :: SPH0

! variables de sortie!
! variables locales
	integer :: i1,i2,i4,no,n2,m2,N1,nanm,NM,num
	integer :: n2p,n2m,m0,m0m,m0p,m11
	integer :: acn_nmm,acn_nmmm,acn_nmmp,acn_npm,acn_npmm,acn_npmp
	real(prec), dimension(1:3) :: Rlp
	real(prec), dimension(1:3) :: M
!	real(prec) :: theta,phi
	complex(prec) :: Ynm((N(NL)+2)**2)
	complex(prec) :: dx,dy,dz
	complex(prec) :: jY_1,jY_2,jY_3,jY_4,jY_5,jY_6
	real(prec), dimension(3) :: Psph
!	character(len=1) :: answer
    real(prec) :: jn(0:N(NL)+1),Djn(0:N(NL)+1)  !,npm,nR

!
! précalcul
!
!	write(*,*) 'pass expansion_propa_BF'
    N1=N(NL)
	nanm=(N1+1)**2

	allocate(SPH(1)%W1(1:Next*nanm))
	allocate(SPH(1)%iw(1:Next*nanm))
	allocate(SPH(1)%jw(1:Next*nanm))
	SPH(1)%nnzw=Next*nanm
	SPH(1)%ncw=nanm*octree(NL)%NL
	SPH(1)%nlw=Next

    num=0
	do i1=1,octree(NL)%NL 	! boucle cellules level le plus fin
	    M=octree(NL)%cell(i1)%center     ! coordonnées du centre de la cellule
	    do i2=1,octree(NL)%cell(i1)%Nb_c
	        no=octree(NL)%cell(i1)%child(i2)  !% n° du point central
	        if (no>NE) then
				num=num+1
				Rlp=P_prop(no,:)-M
				!nR=sqrt(sum(Rlp**2))
				!npm=Rlp(1)*norc(no,1)+Rlp(2)*norc(no,2)+Rlp(3)*norc(no,3)
				call car2sph(Rlp,Psph,1)
				call SPHJ(N1+1,k*Psph(1),NM,jn,Djn)
				call Ynm_gum(N(NL)+1,Psph(2),Psph(3),Ynm)
!				write(*,*) 'kr=',k*Psph(1),' theta=',Psph(2),' phi=',Psph(3)
				do i4=1,(N(NL)+1)**2
					n2=int(sqrt(real(i4-1)))
					m2=i4-n2**2-n2-1
					n2p=n2+1
					n2m=n2-1
					m0=m2
					m0p=m0+1
					m0m=m0-1
					m11=-m0-1

					acn_npmp=n2p**2+n2p+m0p+1
					acn_npmm=n2p**2+n2p+m0m+1
					acn_nmmp=n2m**2+n2m+m0p+1
					acn_nmmm=n2m**2+n2m+m0m+1
					acn_nmm =n2m**2+n2m+m0+1
					acn_npm =n2p**2+n2p+m0+1
!					acn_nm=n2**2+n2+m2+1

!					jY_1=SPH0%b(1,i4) * jn(n2p) * (Ynm(acn_npmp))
!					jY_3=SPH0%b(3,i4) * jn(n2p) * (Ynm(acn_npmm))
!					jY_6=SPH0%a(2,i4) * jn(n2p) * (Ynm(acn_npm))

					jY_1=SPH0%b_ex(1,i4) * jn(n2p) * (Ynm(acn_npmp))
					jY_3=SPH0%b_ex(3,i4) * jn(n2p) * (Ynm(acn_npmm))
					jY_6=SPH0%a(2,i4) * jn(n2p) * (Ynm(acn_npm))

					if (abs(SPH0%b_ex(2,i4)) > 1.0e-5_prec) then
						jY_2=SPH0%b_ex(2,i4) * jn(n2m) * (Ynm(acn_nmmp))
					else
						jY_2=cmplx(0.0_prec,0.0_prec, kind=prec)
					endif
					if (abs(SPH0%b_ex(4,i4)) > 1.0e-5_prec) then
						jY_4=SPH0%b_ex(4,i4) * jn(n2m) * (Ynm(acn_nmmm))
					else
						jY_4=cmplx(0.0_prec,0.0_prec, kind=prec)
					endif
					if (abs(SPH0%a(1,i4)) > 1.0e-5_prec) then
						jY_5=SPH0%a(1,i4) * jn(n2m) * (Ynm(acn_nmm))
					else
						jY_5=cmplx(0.0_prec,0.0_prec, kind=prec)
					endif
					dx=k/2.0_prec * (jY_1-jY_2+jY_3-jY_4)
					dy=-ji*k/2.0_prec * (jY_1-jY_2-jY_3+jY_4)
					dz=k * (jY_5 - jY_6)
					! moment G
!					SPH(1)%W1((num-1)*nanm+i4)= 4.0_prec*pi*jn(n2)*Ynm(i4) + alpha * &
!					! moment nq.grad(G)
!					4.0_prec*pi*(norc(no,1)*dx+norc(no,2)*dy+norc(no,3)*dz)
					SPH(1)%W1((num-1)*nanm+i4)= 4.0_prec*pi*jn(n2)*Ynm(i4)

					!SPH(1)%W1((num-1)*nanm+i4)=	4.0_prec*pi*(norc(no,1)*dx+norc(no,2)*dy+norc(no,3)*dz)
					SPH(1)%jw((num-1)*nanm+i4)=(i1-1)*nanm+i4
					SPH(1)%iw((num-1)*nanm+i4)=no-NE
				enddo
			endif
	    end do
	end do
!	write(*,*) 'num=',num, 'Next=',Next

	end subroutine expansions_BM_propa_BF
!
! -----------------------------------------------------------------------------------------------------------------------------
! calcul propagation par ONLY FMM et pour HRTF car points éloignés de la structure
! calcul p=Dps - M dnps : le signe '-' est intégré au calcul, mettre (+) dnps en entrée
! -----------------------------------------------------------------------------------------------------------------------------
!
	subroutine propagation_FMM_only(LEV,k,NL,N,octree,SPH,SPH0,var_p,var_dnp,var_out)
	implicit none
! variables d'entrée
	integer, intent(in) :: NL,LEV
	integer, dimension(1:NL), intent(in) :: N
	type(FMM),intent(inout) :: SPH(NL)
	type(FMM0),intent(inout) :: SPH0
	type(tree) :: octree(NL)
	real(prec) :: k
	complex(prec) :: var_p(NE),var_dnp(NE)
	! variable de sortie
	complex(prec) :: var_out(Next)
	! variables locales
	integer :: i1, i2
	integer :: LEVa   ! niveau de propagation effectif (octree adaptatif : voir plus bas)
	logical :: interp_sauve
	integer, allocatable :: npu_gro
	real(prec) :: rlm(3)
	complex(prec), dimension(:), allocatable :: M2L
!	character(len=1) :: answer


	! montee avec interpoler PUIS dephaser (voir interp_puis_dephase)
	interp_sauve=interp_puis_dephase
	interp_puis_dephase=.true.
	LEVa=LEV
	if (fmm_adaptatif) then
		! octree adaptatif : un element pose AU-DESSUS de LEV n'est pas dans les
		! multipoles du niveau LEV -> on propage depuis le niveau le plus grossier
		! qui les contient tous (points exterieurs lointains : translation valable)
		LEVa=min(LEV,minval(niv_E(1:NE)))
		if (LEVa < LEV .and. .not. omp_in_parallel()) write(*,'(A,I3,A,I3)') &
			' propagation_FMM_only (adaptatif) : niveau de propagation ', LEV, ' -> ', LEVa
		call montee_T0_adapt(octree,NL,N,SPH,SPH0,k,3,var_p,var_dnp)
	else
  	if (SPH(NL)%BF==1) then
		SPH(NL)%inBF(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
		do i1=1,SPH(1)%nnzv
			SPH(NL)%inBF(SPH(1)%iv(i1))=SPH(NL)%inBF(SPH(1)%iv(i1))-SPH(1)%V1(i1)*var_dnp(SPH(1)%jv(i1))*ji*k/(4.0_prec*pi) + &  ! gradient
				SPH(1)%V2(i1)*var_p(SPH(1)%jv(i1))*(ji*k)/(4.0_prec*pi) ! pression
		end do
   else
		SPH(NL)%in(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
		do i1=1,SPH(1)%nnzv
			SPH(NL)%in(SPH(1)%iv(i1))=SPH(NL)%in(SPH(1)%iv(i1))-SPH(1)%V1(i1)*var_dnp(SPH(1)%jv(i1))*ji*k/(4.0_prec*pi)**2 + &  ! gradient
				SPH(1)%V2(i1)*var_p(SPH(1)%jv(i1))*(k)**2/(4.0_prec*pi)**2 ! pression
		end do
   endif
 !   write(*,*) 'iv=',SPH(1)%iv(1:3)
  !  write(*,*) 'jv=',SPH(1)%jv(1:3)
   ! write(*,*) 'V1=',SPH(1)%V1(1:3)

!    write(*,*) 'SPH(NL)%inBF',SPH(NL)%inBF(1:10)
 !   write(*,*) 'SPH(NL)%in',SPH(NL)%in(1:10)
	! pass upward
   	!call pass_upward_FFT(NL,N,octree,SPH,k)
   	call pass_upward_HYB(NL,N,octree,SPH,SPH0,k)
	endif
	interp_puis_dephase=interp_sauve

!	! moment p et dnp
!	var_out(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
!	SPH(NL)%in(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
!	do i1=1,SPH(1)%nnzv
!		SPH(NL)%in(SPH(1)%iv(i1))=SPH(NL)%in(SPH(1)%iv(i1))+k**2/(4.0_prec*pi)**2 * SPH(1)%V2(i1)*var_p(SPH(1)%jv(i1)) - &
!			ji*k/(4.0_prec*pi)**2 * SPH(1)%V1(i1)*var_dnp(SPH(1)%jv(i1))
!	end do
!
 ! 	! pass upward
  ! 	call pass_upward_FFT(NL,N,octree,SPH,k)
!! conversion en HF si besoin
	if (SPH(LEVa)%BF==1) then ! cas BF, on passe en HF
		call BF2HF_up_FMM(octree,NL,N,SPH,LEVa)
	endif



!  	! translation jusqu'aux points Pext
	npu_gro=size(SPH(LEVa)%Ps,1)  !% nbre de points sur sph unitaire niveau
!	write(*,*) 'npu1=',npu_gro
!	npu_gro=(N(1)+1)*(2*N(1)+1)!	if (allocated(SPH(1)%out)) deallocate(SPH(1)%out)
!	write(*,*) 'npu2=',npu_gro
!	allocate(SPH(1)%out(Next*npu_gro))
	if (allocated(M2L)) deallocate(M2L)
	allocate(M2L(npu_gro))
	var_out(:)=cmplx(0.0_prec,0.0_prec, kind=prec)

	do i2=1,Next
		do i1=1,octree(LEVa)%NL
!			write(*,*) 'i1=',i1,' center=',octree(1)%cell(i1)%center
			rlm=Pext(i2,:)-octree(LEVa)%cell(i1)%center

			call T_LR(NL,N(LEVa),LEVa,k,SPH,rlm,M2L)
			var_out(i2)=var_out(i2) + sum( M2L * &
				SPH(LEVa)%in((i1-1)*npu_gro+1:(i1)*npu_gro)*SPH(LEVa)%wp)
		end do
	enddo

	end subroutine propagation_FMM_only
!
! -----------------------------------------------------------------------------------------------------------------------------
! calcul propagation par MLFMM
! -----------------------------------------------------------------------------------------------------------------------------
!

	subroutine propa_FMM_p(NL,octree,SPH0,k,var_p,var_dnp,pout)
! variables d'entrée
	integer :: NL
	type(FMM0) :: SPH0
	type(tree) :: octree(NL)
	real(prec) :: k
	complex(prec) :: var_p(NE),var_dnp(NE),var1(Next),var2(Next)
	! variable de sortie
	complex(prec) :: pout(Next)
	! variables locales
	integer :: i0,i2,i3,i4
	integer, dimension(1:NL) :: N
!	integer, allocatable :: npu_gro
!	real(prec) :: rlm(3)
!	complex(prec), dimension(:), allocatable :: M2L
	integer :: no,no1
	integer, dimension(1:3) :: nPe
	real(prec), dimension(1:3,1:3) :: PQ
	real(prec) :: Qc(1:3,IS_c), Qv(1:3,IS_v)
	real(prec), dimension(1:3,1:IS_c) :: F_c, dF1_c, dF2_c
    real(prec), dimension(1:3,1:IS_v) :: F_v, dF1_v, dF2_v
	real(prec), dimension(1:3,1:IS_c) :: PG_c
	real(prec), dimension(1:3,1:IS_v) :: PG_v
	complex(prec) :: var(NE+Next)
	complex(prec) :: D1,M1
	type(FMM)	:: SPH(NL)
	real(prec) :: eps0_k
	! DEBUG : character(len=1) :: answer

!    write(*,*) 'Cellules par niveaux dans propa',octree(:)%NL
!	write(*,*) 'avant propa_dnp'
!	write(*,*) 'nchild=',octree(NL)%cell(2481)%Nb_c
!	write(*,*) 'no=',octree(NL)%cell(2481)%child(1)
!	write(*,*) 'NL=',NL
!	write(*,*) 'NE=',NE
	pout=cmplx(0.0_prec,0.0_prec, kind=prec)

	! init FMM0 (epsilon0 de la bande de frequence, parametres)
	eps0_k = eps0_de_k(k)
	call init_FMM_T0(octree,NL,k,N,SPH,SPH0,FMM_meth,eps0_k)
	! contribution FMM
	var=cmplx(0.0_prec,0.0_prec, kind=prec)
	var(1:NE)=var_dnp
	call MLFMM_HYB_propa_dnp(var,octree,NL,N,SPH,SPH0,k,var2)
	var=cmplx(0.0_prec,0.0_prec, kind=prec)
	var(1:NE)=var_p
	call MLFMM_HYB_propa_p(var,octree,NL,N,SPH,SPH0,k,var1)

	pout=var1-var2

	! contribution BEM
	call pgt6(PG_c,IS_c)  ! points dans la même cellule
	call pgt6(PG_v,IS_v)  ! points dans la cellule voisine
	call FORM_T3(PG_c(1:2,:),IS_c,F_c,dF1_c,dF2_c)
	call FORM_T3(PG_v(1:2,:),IS_v,F_v,dF1_v,dF2_v)
!! init
	do i0=1,octree(NL)%NL   ! boucle sur cellule dernier niveau
		do i2=1,octree(NL)%cell(i0)%Nb_c
			no=octree(NL)%cell(i0)%child(i2)
			if (no <= NE) then
				nPe=E(no,:) ! Points de l'élément
				PQ=transpose(Ps(nPe,:)) ! Coordonnées des points de l'éléments
				Qc=matmul(PQ,F_c)
				Qv=matmul(PQ,F_v)
				! points dans la cellule
				do i3=1,octree(NL)%cell(i0)%Nb_c
					no1=octree(NL)%cell(i0)%child(i3)
					if (no1 > NE) then
						call DM_T0(k,Qc,norc(no,:),Jac0(no),IS_c,PG_c(3,:),P_prop(no1,:),D1,M1)
						pout(no1-NE)=pout(no1-NE)+D1*var_p(no)-M1*var_dnp(no)
					endif
				end do
!! points dans les cellules voisines
				do i3=1,octree(NL)%cell(i0)%Nb_v
					do i4=1,octree(NL)%cell(octree(NL)%cell(i0)%voisin(i3))%Nb_c
						no1=octree(NL)%cell(octree(NL)%cell(i0)%voisin(i3))%child(i4)
						if (no1 > NE) then
							call DM_T0(k,Qv,norc(no,:),Jac0(no),IS_v,PG_v(3,:),P_prop(no1,:),D1,M1)
							pout(no1-NE)=pout(no1-NE)+D1*var_p(no)-M1*var_dnp(no)
						endif
					end do
				end do
			endif
		end do
	end do




	end subroutine propa_FMM_p

!============================================================================
! OCTREE ADAPTATIF T0 (HRTF), 2026-09-29
!
! Meme principe que le T6 (module_MLFMM_T6, fmm_T6_adapt) : l'element no est
! pose au niveau il = niv_E(no) (module_BEM) ; son moment est pris au centre
! de SA cellule de ce niveau, a l'ordre N(il), au format du niveau (BF ou
! HF) ; la montee ajoute a chaque niveau les moments de ses elements ; les
! cibles (barycentres) sont evaluees dans leur cellule feuille (feuille,
! il < NL) ou au niveau NL (a_cible), voir marquer_sources_cibles.
!
! Formules RECOPIEES de moments_FMM_BF/HF, expansions_BM_BF/HF et des
! injections de MLFMM_HYB_dnp / MLFMM_HYB_p / propagation_FMM_only : avec
! niv_E(:) = NL, resultat identique au bit pres.
!
! NORMALISATION HF : en T0 le facteur est applique a l'INJECTION (x ik/(4pi)^2
! pour V1, k^2/(4pi)^2 pour V2 en HF ; ik/(4pi) en BF). BF2HF_up_FMM d'un
! moment BF x ik/(4pi) donne le moment HF brut x ik/(4pi)^2 : un element
! injecte a un niveau HF prend donc le MEME facteur que si la feuille etait
! HF, quel que soit le format de la feuille. De meme l'expansion HF
! (exp.wp) vaut 4 pi fois l'expansion BF sans facteur (l'expansion BF porte
! deja 4 pi) : aucun facteur supplementaire aux feuilles HF.
!
! Stockage par niveau : SPH(il)%V1/V2 (moments, blocs sans iv/jv),
! SPH(il)%W1 (expansions Burton/Miller, blocs sans iw/jw), SPH(il)%nnzv,
! SPH(il)%nnzw. Drapeaux : allocated(SPH(NL)%V1), allocated(SPH(NL)%W1).
! NE PAS MELANGER avec le chemin uniforme (SPH(1)%V1/W1 + iv/jv/iw/jw).
!============================================================================

!----------------------------------------------------------------------------
! pre_moments_T0_niv : moments V1 (G) et V2 (dnG) des elements du niveau il
! (formules de moments_FMM_BF / moments_FMM_HF), parcours cellules du niveau
! il -> pts d'element de niveau il -> coefficient.
!----------------------------------------------------------------------------
	subroutine pre_moments_T0_niv(octree,NL,il,N,SPH,SPH0,k)
	implicit none
	integer, intent(in) :: NL, il
	integer, dimension(1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	type(FMM),intent(inout) :: SPH(NL)
	type(FMM0), intent(in) :: SPH0
	type(tree), intent(in) :: octree(NL)
	integer :: ic,i2,i4,no,nanm,nu,NM,num,nel
	integer :: m2,n2,n2p,n2m,m0,m0m,m0p
	integer :: acn_nmm,acn_nmmm,acn_nmmp,acn_npm,acn_npmm,acn_npmp,acn_nm
	real(prec), dimension(2,IS_f) :: PG2
	real(prec), dimension(3,IS_f) :: PG_f
	real(prec), dimension(IS_f) :: theta,phi
	real(prec) :: Rmq(3,IS_f),Rmq2(IS_f,3),Pmq(IS_f,3),M(3),Nor(3)
	real(prec), dimension(3,IS_f) :: F_f,dF1_f,dF2_f,Q
	real(prec), dimension(3,3) :: E1
	real(prec), allocatable :: jn(:),Djn(:),jn2(:,:),rqs(:,:),rns(:,:),rns0(:)
	complex(prec), allocatable :: Ynm(:,:)
	complex(prec),dimension(IS_f) :: dx,dy,dz
	complex(prec),dimension(IS_f) :: jY_1,jY_2,jY_3,jY_4,jY_5,jY_6

	call pgt6(PG_f,IS_f)
	PG2=PG_f(1:2,:)
	call FORM_T3(PG2,IS_f,F_f,dF1_f,dF2_f)
	nel=count(niv_E(1:NE)==il)
	if (allocated(SPH(il)%V1)) deallocate(SPH(il)%V1)
	if (allocated(SPH(il)%V2)) deallocate(SPH(il)%V2)

	if (SPH(il)%BF==1) then
		nanm=(N(il)+1)**2
		if (nanm > size(SPH0%a,2)) stop 'pre_moments_T0_niv : ordre BF superieur a celui de init_BF_MLFMM_v1'
		if (fusion_Vc) then                  ! Diff : Vc dans V2, V1 factice
			allocate(SPH(il)%V1(1), SPH(il)%V2(nanm*nel))
		else
			allocate(SPH(il)%V1(nanm*nel), SPH(il)%V2(nanm*nel))
		endif
		allocate(jn(0:N(il)+1),Djn(0:N(il)+1),jn2(0:N(il)+1,IS_f),Ynm((N(il)+2)**2,IS_f))
		num=0
		do ic=1,octree(il)%NL
			M=octree(il)%cell(ic)%center
			do i2=1,octree(il)%cell(ic)%Nb_p
				no=octree(il)%cell(ic)%pts(i2)
				if (no > NE) cycle                ! cibles separees (CHI) : pas de moment
				if (niv_E(no) /= il) cycle
				E1=transpose(Ps(E(no,:),:))
				Q=matmul(E1,F_f)
				num=num+1
				Rmq=Q-spread(M,2,IS_f)
				Rmq2=transpose(Rmq)
				call car2sph(Rmq2,Pmq,IS_f)
				theta=Pmq(:,2)
				phi=Pmq(:,3)
				do i4=1,IS_f
					call SPHJ(N(il)+1,k*Pmq(i4,1),NM,jn,Djn)
					jn2(:,i4)=jn
				enddo
				call Ynm_gum(N(il)+1,theta,phi,IS_f,Ynm)
				do i4=1,(N(il)+1)**2
					n2=int(sqrt(real(i4-1)))
					m2=i4-n2**2-n2-1
					n2p=n2+1
					n2m=n2-1
					m0=-m2
					m0p=m0+1
					m0m=m0-1
					acn_npmp=n2p**2+n2p+m0p+1
					acn_npmm=n2p**2+n2p+m0m+1
					acn_nmmp=n2m**2+n2m+m0p+1
					acn_nmmm=n2m**2+n2m+m0m+1
					acn_nmm =n2m**2+n2m+m0+1
					acn_npm =n2p**2+n2p+m0+1
					acn_nm= n2**2+n2+m0+1
					jY_1=SPH0%b(1,i4) * jn2(n2p,:) * (Ynm(acn_npmp,:))
					jY_3=SPH0%b(3,i4) * jn2(n2p,:) * (Ynm(acn_npmm,:))
					jY_6=SPH0%a(2,i4) * jn2(n2p,:) * (Ynm(acn_npm,:))
					if (abs(SPH0%b(2,i4)) > 1.0e-5_prec) then
						jY_2=SPH0%b(2,i4) * jn2(n2m,:) * (Ynm(acn_nmmp,:))
					else
						jY_2=cmplx(0.0_prec,0.0_prec, kind=prec)
					endif
					if (abs(SPH0%b(4,i4)) > 1.0e-5_prec) then
						jY_4=SPH0%b(4,i4) * jn2(n2m,:) * (Ynm(acn_nmmm,:))
					else
						jY_4=cmplx(0.0_prec,0.0_prec, kind=prec)
					endif
					if (abs(SPH0%a(1,i4)) > 1.0e-5_prec) then
						jY_5=SPH0%a(1,i4) * jn2(n2m,:) * (Ynm(acn_nmm,:))
					else
						jY_5=cmplx(0.0_prec,0.0_prec, kind=prec)
					endif
					dx=k/2.0_prec * (jY_1-jY_2+jY_3-jY_4)
					dy=-ji*k/2.0_prec * (jY_1-jY_2-jY_3+jY_4)
					dz=k * (jY_5 - jY_6)
					if (fusion_Vc) then
						SPH(il)%V2((num-1)*nanm+i4)= sum( &
							(norc(no,1)*dx+norc(no,2)*dy+norc(no,3)*dz) * Jac0(no) * PG_f(3,:)) + &
							ji*k*beta_Vc * sum(jn2(n2,:)*(Ynm(acn_nm,:)) * Jac0(no) * PG_f(3,:))
					else
					SPH(il)%V1((num-1)*nanm+i4)= sum(jn2(n2,:)*(Ynm(acn_nm,:)) * Jac0(no) * PG_f(3,:))
					SPH(il)%V2((num-1)*nanm+i4)= sum( &
						(norc(no,1)*dx+norc(no,2)*dy+norc(no,3)*dz) * Jac0(no) * PG_f(3,:))
					endif
				enddo
			enddo
		enddo
		SPH(il)%nnzv=num*nanm
	else
		nu=(N(il)+1)*(2*N(il)+1)
		if (fusion_Vc) then
			allocate(SPH(il)%V1(1), SPH(il)%V2(nu*nel))
		else
			allocate(SPH(il)%V1(nu*nel), SPH(il)%V2(nu*nel))
		endif
		allocate(rqs(nu,IS_f),rns(nu,IS_f),rns0(nu))
		num=0
		do ic=1,octree(il)%NL
			M=octree(il)%cell(ic)%center
			do i2=1,octree(il)%cell(ic)%Nb_p
				no=octree(il)%cell(ic)%pts(i2)
				if (no > NE) cycle                ! cibles separees (CHI) : pas de moment
				if (niv_E(no) /= il) cycle
				E1=transpose(Ps(E(no,:),:))
				Q=matmul(E1,F_f)
				Rmq=spread(M,2,IS_f)-Q
				rqs=0.0_prec
				call dgemm('N','N',size(SPH(il)%Ps,1),IS_f,3,1.0_prec,SPH(il)%Ps,size(SPH(il)%Ps,1),Rmq,3,0.0_prec,rqs,&
				size(SPH(il)%Ps,1))
				Nor=norc(no,:)
				rns0=matmul(SPH(il)%Ps,Nor)
				rns=spread(rns0,2,IS_f)
				num=num+1
				if (fusion_Vc) then
					SPH(il)%V2((num-1)*nu+1:num*nu)=matmul(exp(ji*k*rqs)*rns*Jac0(no),PG_f(3,:)) - &
						beta_Vc*matmul(exp(ji*k*rqs)*Jac0(no),PG_f(3,:))
				else
				SPH(il)%V2((num-1)*nu+1:num*nu)=matmul(exp(ji*k*rqs)*rns*Jac0(no),PG_f(3,:))
				SPH(il)%V1((num-1)*nu+1:num*nu)=matmul(exp(ji*k*rqs)*Jac0(no),PG_f(3,:))
				endif
			enddo
		enddo
		SPH(il)%nnzv=num*nu
	endif
	end subroutine pre_moments_T0_niv

!----------------------------------------------------------------------------
! injecter_T0_niv : AJOUTE au developpement du niveau il les moments de ses
! elements. mode 1 : dnp (V1), 2 : p (V2), 3 : propagation (V2.p - V1.dnp),
! 4 : p_dnp du GMRES avec impedance (V2.p - V1.dnp, ecriture de
! MLFMM_HYB_p_dnp). Expressions recopiees de MLFMM_HYB_dnp, MLFMM_HYB_p,
! propagation_FMM_only, MLFMM_HYB_p_dnp (bit pres en uniforme).
!----------------------------------------------------------------------------
	subroutine injecter_T0_niv(octree,NL,il,N,SPH,k,mode,varp,vardnp)
	implicit none
	integer, intent(in) :: NL, il, mode
	integer, dimension(1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	type(tree), intent(in) :: octree(NL)
	type(FMM), intent(inout) :: SPH(NL)
	complex(prec), intent(in) :: varp(NE), vardnp(NE)
	integer :: ic, i2, i4, no, ncoef, idx, ico, mode_e

	if (SPH(il)%nnzv == 0) return
	mode_e = mode                              ! Diff, Vc = V2 fusionne : V2.p - V1.dnp = Vc.p
	if (fusion_Vc .and. (mode == 3 .or. mode == 4)) mode_e = 2
	if (fusion_Vc .and. mode == 1) stop 'injecter_T0_niv : mode 1 (dnp seul) incompatible avec fusion_Vc'
	if (SPH(il)%BF==1) then
		ncoef=(N(il)+1)**2
	else
		ncoef=(N(il)+1)*(2*N(il)+1)
	endif
	idx=0
	do ic=1,octree(il)%NL
		do i2=1,octree(il)%cell(ic)%Nb_p
			no=octree(il)%cell(ic)%pts(i2)
			if (no > NE) cycle                ! cibles separees (CHI) : pas de moment
			if (niv_E(no) /= il) cycle
			do i4=1,ncoef
				idx=idx+1
				ico=(ic-1)*ncoef+i4
				if (SPH(il)%BF==1) then
					select case (mode_e)
					case (1)
						SPH(il)%inBF(ico)=SPH(il)%inBF(ico)+SPH(il)%V1(idx)*vardnp(no)*ji*k/(4.0_prec*pi)
					case (2)
						SPH(il)%inBF(ico)=SPH(il)%inBF(ico)+SPH(il)%V2(idx)*varp(no)*(ji*k)/(4.0_prec*pi)
					case (4)
						SPH(il)%inBF(ico)=SPH(il)%inBF(ico)+SPH(il)%V2(idx)*varp(no)*(ji*k)/(4.0_prec*pi) - &
							ji*k/(4.0_prec*pi) * SPH(il)%V1(idx)*vardnp(no)
					case default
						SPH(il)%inBF(ico)=SPH(il)%inBF(ico)-SPH(il)%V1(idx)*vardnp(no)*ji*k/(4.0_prec*pi) + &
							SPH(il)%V2(idx)*varp(no)*(ji*k)/(4.0_prec*pi)
					end select
				else
					select case (mode_e)
					case (1)
						SPH(il)%in(ico)=SPH(il)%in(ico)+SPH(il)%V1(idx)*vardnp(no)*ji*k/(4.0_prec*pi)**2
					case (2)
						SPH(il)%in(ico)=SPH(il)%in(ico)+SPH(il)%V2(idx)*varp(no)*(k)**2/(4.0_prec*pi)**2
					case (4)
						SPH(il)%in(ico)=SPH(il)%in(ico)+SPH(il)%V2(idx)*varp(no)*(k)**2/(4.0_prec*pi)**2  - &
							ji*k/(4.0_prec*pi)**2 * SPH(il)%V1(idx)*vardnp(no)
					case default
						SPH(il)%in(ico)=SPH(il)%in(ico)-SPH(il)%V1(idx)*vardnp(no)*ji*k/(4.0_prec*pi)**2 + &
							SPH(il)%V2(idx)*varp(no)*(k)**2/(4.0_prec*pi)**2
					end select
				endif
			enddo
		enddo
	enddo
	end subroutine injecter_T0_niv

!----------------------------------------------------------------------------
! montee_T0_adapt : moments du niveau fin (developpement remis a zero), puis
! montee (upward_BF/HF inchanges) avec les moments de chaque niveau.
!----------------------------------------------------------------------------
	subroutine montee_T0_adapt(octree,NL,N,SPH,SPH0,k,mode,varp,vardnp)
	implicit none
	integer, intent(in) :: NL, mode
	integer, dimension(1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	type(tree) :: octree(NL)
	type(FMM), intent(inout) :: SPH(NL)
	type(FMM0), intent(inout) :: SPH0
	complex(prec), intent(in) :: varp(NE), vardnp(NE)
	integer :: il, i1

	! moments : drapeau sur V2 (V1 peut avoir ete libere apres le second membre, liberer_V1) ;
	! les modes 1, 3 et 4 ont besoin de V1 : recalcul (identique au bit pres) s'il manque
	if (.not. allocated(SPH(NL)%V2) .or. (mode /= 2 .and. .not. allocated(SPH(NL)%V1))) then
		do il=1,NL
			call pre_moments_T0_niv(octree,NL,il,N,SPH,SPH0,k)
		enddo
	endif
	if (SPH(NL)%BF==1) then
		if (.not. allocated(SPH(NL)%inBF)) allocate(SPH(NL)%inBF((N(NL)+1)**2*octree(NL)%NL))
		SPH(NL)%inBF(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
	else
		if (.not. allocated(SPH(NL)%in)) allocate(SPH(NL)%in((N(NL)+1)*(2*N(NL)+1)*octree(NL)%NL))
		SPH(NL)%in(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
	endif
	call injecter_T0_niv(octree,NL,NL,N,SPH,k,mode,varp,vardnp)
	do i1=NL,2,-1
		if (SPH(i1)%BF==1) then
			call upward_BF(NL,N,SPH,SPH0,octree,i1)
		else
			call upward_HF(NL,N,SPH,octree,k,i1)
		endif
		call injecter_T0_niv(octree,NL,i1-1,N,SPH,k,mode,varp,vardnp)
	enddo
	end subroutine montee_T0_adapt

!----------------------------------------------------------------------------
! pre_expansions_T0_niv : developpements Burton/Miller W1 aux cibles evaluees
! au niveau il (feuille si il < NL, a_cible si il = NL), formules de
! expansions_BM_BF / expansions_BM_HF. SPH(il)%W1, SPH(il)%nnzw.
!----------------------------------------------------------------------------
	subroutine pre_expansions_T0_niv(octree,NL,il,N,SPH,SPH0,k,alpha)
	implicit none
	integer, intent(in) :: NL, il
	integer, dimension(1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	complex(prec), intent(in) :: alpha
	type(tree), intent(in) :: octree(NL)
	type(FMM), intent(inout) :: SPH(NL)
	type(FMM0), intent(in) :: SPH0
	integer :: ic, i2, i4, no, N1, nanm, nu, NM, num, ncib, ip
	integer :: n2, m2, n2p, n2m, m0, m0m, m0p
	integer :: acn_nmm,acn_nmmm,acn_nmmp,acn_npm,acn_npmm,acn_npmp
	real(prec) :: Rlp(3), M(3), Psph(3), Pcib(3), cG, cH
	real(prec), allocatable :: rps(:), rns(:), jn(:), Djn(:)
	complex(prec), allocatable :: Ynm(:)
	complex(prec) :: dx,dy,dz,jY_1,jY_2,jY_3,jY_4,jY_5,jY_6

	ncib=0
	do ic=1,octree(il)%NL
		if (.not. cible_evaluee(octree,NL,il,ic)) cycle
		do i2=1,octree(il)%cell(ic)%Nb_p
			if (.not. cible_active(octree(il)%cell(ic)%pts(i2))) cycle
			ncib=ncib+1
		enddo
	enddo
	if (allocated(SPH(il)%W1)) deallocate(SPH(il)%W1)
	N1=N(il)
	if (SPH(il)%BF==1) then
		nanm=(N1+1)**2
		allocate(SPH(il)%W1(ncib*nanm))
		allocate(jn(0:N1+1),Djn(0:N1+1),Ynm((N1+2)**2))
		num=0
		do ic=1,octree(il)%NL
			if (.not. cible_evaluee(octree,NL,il,ic)) cycle
			M=octree(il)%cell(ic)%center
			do i2=1,octree(il)%cell(ic)%Nb_p
				ip=octree(il)%cell(ic)%pts(i2)
				if (.not. cible_active(ip)) cycle
				call cible_T0(ip,no,Pcib,cG,cH)    ! ligne, point (CHI : interieur), poids (hybride)
				num=num+1
				Rlp=Pcib-M
				call car2sph(Rlp,Psph,1)
				call SPHJ(N1+1,k*Psph(1),NM,jn,Djn)
				call Ynm_gum(N1+1,Psph(2),Psph(3),Ynm)
				do i4=1,(N1+1)**2
					n2=int(sqrt(real(i4-1)))
					m2=i4-n2**2-n2-1
					n2p=n2+1
					n2m=n2-1
					m0=m2
					m0p=m0+1
					m0m=m0-1
					acn_npmp=n2p**2+n2p+m0p+1
					acn_npmm=n2p**2+n2p+m0m+1
					acn_nmmp=n2m**2+n2m+m0p+1
					acn_nmmm=n2m**2+n2m+m0m+1
					acn_nmm =n2m**2+n2m+m0+1
					acn_npm =n2p**2+n2p+m0+1
					jY_1=SPH0%b_ex(1,i4) * jn(n2p) * (Ynm(acn_npmp))
					jY_3=SPH0%b_ex(3,i4) * jn(n2p) * (Ynm(acn_npmm))
					jY_6=SPH0%a(2,i4) * jn(n2p) * (Ynm(acn_npm))
					if (abs(SPH0%b_ex(2,i4)) > 1.0e-5_prec) then
						jY_2=SPH0%b_ex(2,i4) * jn(n2m) * (Ynm(acn_nmmp))
					else
						jY_2=cmplx(0.0_prec,0.0_prec, kind=prec)
					endif
					if (abs(SPH0%b_ex(4,i4)) > 1.0e-5_prec) then
						jY_4=SPH0%b_ex(4,i4) * jn(n2m) * (Ynm(acn_nmmm))
					else
						jY_4=cmplx(0.0_prec,0.0_prec, kind=prec)
					endif
					if (abs(SPH0%a(1,i4)) > 1.0e-5_prec) then
						jY_5=SPH0%a(1,i4) * jn(n2m) * (Ynm(acn_nmm))
					else
						jY_5=cmplx(0.0_prec,0.0_prec, kind=prec)
					endif
					dx=k/2.0_prec * (jY_1-jY_2+jY_3-jY_4)
					dy=-ji*k/2.0_prec * (jY_1-jY_2-jY_3+jY_4)
					dz=k * (jY_5 - jY_6)
					SPH(il)%W1((num-1)*nanm+i4)= cG*4.0_prec*pi*jn(n2)*Ynm(i4) + cH*alpha * &
						4.0_prec*pi*(norc(no,1)*dx+norc(no,2)*dy+norc(no,3)*dz)
				enddo
			enddo
		enddo
		SPH(il)%nnzw=num*nanm
	else
		nu=(N1+1)*(2*N1+1)
		allocate(SPH(il)%W1(ncib*nu), rps(nu), rns(nu))
		num=0
		do ic=1,octree(il)%NL
			if (.not. cible_evaluee(octree,NL,il,ic)) cycle
			M=octree(il)%cell(ic)%center
			do i2=1,octree(il)%cell(ic)%Nb_p
				ip=octree(il)%cell(ic)%pts(i2)
				if (.not. cible_active(ip)) cycle
				call cible_T0(ip,no,Pcib,cG,cH)    ! ligne, point (CHI : interieur), poids (hybride)
				Rlp=Pcib-M
				rps=matmul(SPH(il)%Ps,Rlp)
				rns=matmul(SPH(il)%Ps,norc(no,:))
				num=num+1
				SPH(il)%W1((num-1)*nu+1:num*nu)=cG*exp(ji*k*rps)*SPH(il)%wp + &
					cH*alpha * ji * k * exp(ji*k*rps)*rns*SPH(il)%wp
			enddo
		enddo
		SPH(il)%nnzw=num*nu
	endif
	end subroutine pre_expansions_T0_niv

	! cellule ic du niveau il dont les cibles sont evaluees a ce niveau
	logical function cible_evaluee(octree,NL,il,ic)
	integer, intent(in) :: NL, il, ic
	type(tree), intent(in) :: octree(NL)
	if (il == NL) then
		cible_evaluee = octree(il)%cell(ic)%a_cible
	else
		cible_evaluee = octree(il)%cell(ic)%feuille
	endif
	end function cible_evaluee

!----------------------------------------------------------------------------
! appliquer_expansions_T0_niv : varout(no) += W1 . out (meme ordre que la
! boucle sur iw/jw de MLFMM_HYB_p / MLFMM_HYB_dnp).
!----------------------------------------------------------------------------
	subroutine appliquer_expansions_T0_niv(octree,NL,il,N,SPH,varout)
	implicit none
	integer, intent(in) :: NL, il
	integer, dimension(1:NL), intent(in) :: N
	type(tree), intent(in) :: octree(NL)
	type(FMM), intent(in) :: SPH(NL)
	complex(prec), intent(inout) :: varout(NE)
	integer :: ic, i2, i4, no, ncoef, idx, ip

	if (SPH(il)%nnzw == 0) return
	if (SPH(il)%BF==1) then
		ncoef=(N(il)+1)**2
	else
		ncoef=(N(il)+1)*(2*N(il)+1)
	endif
	idx=0
	do ic=1,octree(il)%NL
		if (.not. cible_evaluee(octree,NL,il,ic)) cycle
		do i2=1,octree(il)%cell(ic)%Nb_p
			ip=octree(il)%cell(ic)%pts(i2)
			if (.not. cible_active(ip)) cycle
			no=ip
			if (cibles_sep .and. ip > NE) no=ip-NE
			do i4=1,ncoef
				idx=idx+1
				if (SPH(il)%BF==1) then
					varout(no)=varout(no)+SPH(il)%W1(idx)*SPH(il)%outBF((ic-1)*ncoef+i4)
				else
					varout(no)=varout(no)+SPH(il)%W1(idx)*SPH(il)%out((ic-1)*ncoef+i4)
				endif
			enddo
		enddo
	enddo
	end subroutine appliquer_expansions_T0_niv

!----------------------------------------------------------------------------
! fmm_T0_adapt : produit adaptatif complet. mode 1 : dnp (second membre,
! MLFMM_HYB_dnp), 2 : p (GMRES, MLFMM_HYB_p), 4 : p et dnp (GMRES avec
! impedance, MLFMM_HYB_p_dnp).
!----------------------------------------------------------------------------
	subroutine fmm_T0_adapt(mode,varp,vardnp,octree,NL,N,SPH,SPH0,alpha,k,varout)
	implicit none
	integer, intent(in) :: NL, mode
	integer, dimension(1:NL), intent(in) :: N
	real(prec), intent(in) :: k
	complex(prec), intent(in) :: alpha
	complex(prec), intent(in) :: varp(NE), vardnp(NE)
	type(tree) :: octree(NL)
	type(FMM), intent(inout) :: SPH(NL)
	type(FMM0), intent(inout) :: SPH0
	complex(prec), intent(inout) :: varout(NE)
	integer :: il

	call montee_T0_adapt(octree,NL,N,SPH,SPH0,k,mode,varp,vardnp)
	call pass_downward_HYB(octree,NL,N,SPH,SPH0,k)
	if (.not. allocated(SPH(NL)%W1)) then
		do il=1,NL
			call pre_expansions_T0_niv(octree,NL,il,N,SPH,SPH0,k,alpha)
		enddo
	endif
	varout(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
	do il=1,NL
		call appliquer_expansions_T0_niv(octree,NL,il,N,SPH,varout)
	enddo
	if (mode == 1 .and. beta_Vc /= (0.0_prec,0.0_prec) .and. .not. fusion_Vc) then
		! Ray (Kreuzer) avec impedance (2026-11-06) : second membre fait avec V1 et V2 bruts ; GMRES
		! voit dnp = -i k A_red p : V2 <- Vc (meme combinaison que fusion_Vc de Diff), V1 libere.
		! La propagation (mode 3) recalcule V1 et V2 bruts (montee_T0_adapt) avec le dnp total.
		do il=1,NL
			if (.not. allocated(SPH(il)%V1)) cycle
			if (size(SPH(il)%V1) == size(SPH(il)%V2)) then
				if (SPH(il)%BF==1) then
					SPH(il)%V2 = SPH(il)%V2 + ji*k*beta_Vc*SPH(il)%V1
				else
					SPH(il)%V2 = SPH(il)%V2 - beta_Vc*SPH(il)%V1
				endif
			endif
			deallocate(SPH(il)%V1)
		enddo
	elseif (mode == 1 .and. liberer_V1) then  ! second membre fait : V1 inutile jusqu'a la propagation
		do il=1,NL
			if (allocated(SPH(il)%V1)) deallocate(SPH(il)%V1)
		enddo
	endif
	end subroutine fmm_T0_adapt

!============================================================================
! PROPAGATION HIERARCHIQUE vers les points exterieurs (SOFA), 2026-09-29
!
! Remplace le balayage de propagation_FMM_only (tous les couples (Pext,
! cellule du niveau LEV), T_LR recalcule a chaque couple) par une DESCENTE :
! pour chaque point P, on part des cellules du niveau il0 ;
!   - cellule acceptee si |P - centre| >= fact x cote : son multipole est
!     evalue en P (methode 0 : HF, T_LR comme propagation_FMM_only ;
!     methode 1 : BF au point, 4 pi somme h_n(kr) Y_nm(r) a_nm) ;
!   - sinon : les elements POSES au niveau de la cellule (niv_E, ou tous au
!     niveau fin en uniforme) sont integres directement (BEM, DMKH_T0,
!     D.p - M.dnp), puis on descend dans les enfants.
!
! Normalisation BF : T_LR = somme (2n+1) i^n h_n P_n(s.r) = 4 pi somme i^n h_n
! Y_nm(r) conj(Y_nm(s)) ; avec a_nm = i^n somme_s wp conj(Y_nm(s)) in(s)
! (meme definition que HF2BF_up_FMM) : evaluation HF = 4 pi somme h_n Y_nm a_nm.
!
! Harmoniques : Ynm_normalise_1 (Legendre entierement normalises, recurrence
! stable). Ynm_gum deborde des N >= 85 (gamma(171)) : ne pas l'utiliser ici.
!
! Compteurs (hors region parallele) : n_prop_cellules, n_prop_directs.
!============================================================================

!----------------------------------------------------------------------------
! Ynm_normalise_1 : Y_nm(theta,phi), n = 0..Nmax, rangement ACN n^2+n+m+1,
! MEME convention que Ynm_gum (pas de phase de Condon-Shortley,
! Y_{n,-m} = conj(Y_nm)). Recurrence sur les Legendre normalises
!   Q_00 = 1/sqrt(4 pi), Q_mm = sqrt((2m+1)/(2m)) sin(theta) Q_{m-1,m-1},
!   Q_{m+1,m} = sqrt(2m+3) cos(theta) Q_mm,
!   Q_nm = a_nm (cos(theta) Q_{n-1,m} - b_nm Q_{n-2,m}),
!   a_nm = sqrt((4n^2-1)/(n^2-m^2)), b_nm = sqrt(((n-1)^2-m^2)/(4(n-1)^2-1)).
!----------------------------------------------------------------------------
	subroutine Ynm_normalise_1(Nmax,theta,phi,Ynm)
	implicit none
	integer, intent(in) :: Nmax
	real(prec), intent(in) :: theta, phi
	complex(prec), intent(out) :: Ynm((Nmax+1)**2)
	real(prec) :: x, sn, Qmm, Q1, Q2, Q0, anm, bnm
	complex(prec) :: eim
	integer :: n, m

	x = cos(theta) ; sn = sin(theta)
	Qmm = 1.0_prec/sqrt(4.0_prec*pi)
	do m = 0, Nmax
		if (m > 0) Qmm = sqrt(real(2*m+1,prec)/real(2*m,prec))*sn*Qmm
		eim = exp(ji*real(m,prec)*phi)
		Q2 = Qmm                                   ! Q_{m,m}
		call poser(m, m, Q2)
		if (m < Nmax) then
			Q1 = sqrt(real(2*m+3,prec))*x*Qmm      ! Q_{m+1,m}
			call poser(m+1, m, Q1)
			do n = m+2, Nmax
				anm = sqrt(real(4*n*n-1,prec)/real(n*n-m*m,prec))
				bnm = sqrt(real((n-1)*(n-1)-m*m,prec)/real(4*(n-1)*(n-1)-1,prec))
				Q0 = anm*(x*Q1 - bnm*Q2)
				call poser(n, m, Q0)
				Q2 = Q1 ; Q1 = Q0
			end do
		endif
	end do
	contains
	subroutine poser(nn, mm, Q)
	integer, intent(in) :: nn, mm
	real(prec), intent(in) :: Q
	Ynm(nn*nn+nn+mm+1) = Q*eim
	if (mm > 0) Ynm(nn*nn+nn-mm+1) = conjg(Q*eim)
	end subroutine poser
	end subroutine Ynm_normalise_1

!----------------------------------------------------------------------------
! propagation_hierarchique : voir l'en-tete de la section.
!   il0     : niveau de depart (indice d'octree, 1 = niveau NLF)
!   fact    : critere d'acceptation |P - centre| >= fact x cote (2 conseille)
!   methode : 0 HF (T_LR), 1 BF au point
!   var_p, var_dnp : p et dp/dn sur les elements ; var_out(Next) : p en Pext
!----------------------------------------------------------------------------
	subroutine propagation_hierarchique(il0,fact,methode,k,NL,N,octree,SPH,SPH0,var_p,var_dnp,var_out)
	implicit none
	integer, intent(in) :: il0, methode, NL
	integer :: ild      ! niveau de depart effectif
	real(prec), intent(in) :: fact, k
	integer, dimension(1:NL), intent(in) :: N
	type(tree) :: octree(NL)
	type(FMM), intent(inout) :: SPH(NL)
	type(FMM0), intent(inout) :: SPH0
	complex(prec), intent(in) :: var_p(NE), var_dnp(NE)
	complex(prec), intent(out) :: var_out(Next)

	type coefs_BF
		complex(prec), allocatable :: a(:,:)       ! (ncoef, cellules) : a_nm par cellule
	end type coefs_BF
	type(coefs_BF), allocatable :: cBF(:)
	integer :: i1, il, ic, ipx, is, j, no, nu, nanm, nst, top, NM, nn, mm, ind
	integer, allocatable :: pile_l(:), pile_c(:)
	real(prec) :: P(3), rlm(3), d, Psph(3), th, ph, norP(3)
	real(prec), allocatable :: SJ(:), DJ(:), SY(:), DY(:)
	complex(prec), allocatable :: M2L(:), Yp(:), hn(:)
	complex(prec) :: acc, D1, K1, M1, H1
	real(prec), dimension(1:3,1:IS_v) :: PG_v, F_v, dF1_v, dF2_v
	real(prec) :: PQ(3,3), Qv(3,IS_v)
	integer(8) :: ncel, ndir
	logical :: interp_sauve

	! 1. montee (moments + M2M), comme propagation_FMM_only, avec interpoler
	!    PUIS dephaser (voir interp_puis_dephase)
	interp_sauve=interp_puis_dephase
	interp_puis_dephase=.true.
	if (fmm_adaptatif) then
		call montee_T0_adapt(octree,NL,N,SPH,SPH0,k,3,var_p,var_dnp)
	else
		if (SPH(NL)%BF==1) then
			if (.not. allocated(SPH(1)%V1)) call moments_FMM_BF(octree,NL,N,SPH,SPH0,k)
			SPH(NL)%inBF(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
			do i1=1,SPH(1)%nnzv
				SPH(NL)%inBF(SPH(1)%iv(i1))=SPH(NL)%inBF(SPH(1)%iv(i1))-SPH(1)%V1(i1)*var_dnp(SPH(1)%jv(i1))*ji*k/(4.0_prec*pi) + &
					SPH(1)%V2(i1)*var_p(SPH(1)%jv(i1))*(ji*k)/(4.0_prec*pi)
			end do
		else
			if (.not. allocated(SPH(1)%V1)) call moments_FMM_HF(octree,NL,N,SPH,k)
			SPH(NL)%in(:)=cmplx(0.0_prec,0.0_prec, kind=prec)
			do i1=1,SPH(1)%nnzv
				SPH(NL)%in(SPH(1)%iv(i1))=SPH(NL)%in(SPH(1)%iv(i1))-SPH(1)%V1(i1)*var_dnp(SPH(1)%jv(i1))*ji*k/(4.0_prec*pi)**2 + &
					SPH(1)%V2(i1)*var_p(SPH(1)%jv(i1))*(k)**2/(4.0_prec*pi)**2
			end do
		endif
		call pass_upward_HYB(NL,N,octree,SPH,SPH0,k)
	endif
	interp_puis_dephase=interp_sauve

	! niveau de depart : un element pose AU-DESSUS de il0 ne serait dans aucune
	! cellule visitee -> on part du niveau du plus gros element s'il est plus haut
	ild=il0
	if (fmm_adaptatif) ild=min(il0,minval(niv_E(1:NE)))

	! 2. developpements des niveaux ild..NL au format de la methode
	allocate(cBF(NL))
	do il=ild,NL
		if (methode==0) then
			if (SPH(il)%BF==1) call BF2HF_up_FMM(octree,NL,N,SPH,il)
		else
			nanm=(N(il)+1)**2
			if (SPH(il)%BF==1) then
				allocate(cBF(il)%a(nanm,octree(il)%NL))
				cBF(il)%a=reshape(SPH(il)%inBF,[nanm,octree(il)%NL])
			else
				! a_nm = i^n somme_s wp conj(Y_nm(s)) in(s), harmoniques stables
				nu=(N(il)+1)*(2*N(il)+1)
				allocate(cBF(il)%a(nanm,octree(il)%NL), Yp(nanm))
				cBF(il)%a=cmplx(0.0_prec,0.0_prec, kind=prec)
				do is=1,nu
					th=acos(max(-1.0_prec,min(1.0_prec,SPH(il)%Ps(is,3))))
					ph=atan2(SPH(il)%Ps(is,2),SPH(il)%Ps(is,1))
					call Ynm_normalise_1(N(il),th,ph,Yp)
					Yp=conjg(Yp)*SPH(il)%wp(is)
					do ic=1,octree(il)%NL
						cBF(il)%a(:,ic)=cBF(il)%a(:,ic)+Yp*SPH(il)%in((ic-1)*nu+is)
					end do
				end do
				do nn=0,N(il)
					do mm=-nn,nn
						ind=nn*nn+nn+1+mm
						cBF(il)%a(ind,:)=cBF(il)%a(ind,:)*ji**nn
					end do
				end do
				deallocate(Yp)
			endif
		endif
	end do

	! 3. descente, point par point
	call pgt6(PG_v,IS_v)
	call FORM_T3(PG_v(1:2,:),IS_v,F_v,dF1_v,dF2_v)
	norP=(/0.0_prec,0.0_prec,1.0_prec/)          ! sans effet sur D1 et M1
	nst=octree(ild)%NL+8*NL+8
	allocate(pile_l(nst), pile_c(nst))
	allocate(M2L(size(SPH(ild)%Ps,1)))
	ncel=0 ; ndir=0
	var_out=cmplx(0.0_prec,0.0_prec, kind=prec)
	do ipx=1,Next
		P=Pext(ipx,:)
		acc=cmplx(0.0_prec,0.0_prec, kind=prec)
		top=0
		do ic=octree(ild)%NL,1,-1
			top=top+1 ; pile_l(top)=ild ; pile_c(top)=ic
		end do
		do while (top > 0)
			il=pile_l(top) ; ic=pile_c(top) ; top=top-1
			if (.not. octree(il)%cell(ic)%a_source) cycle
			rlm=P-octree(il)%cell(ic)%center
			d=sqrt(sum(rlm**2))
			if (d >= fact*octree(il)%cell(ic)%length) then
				ncel=ncel+1
				if (methode==0) then
					nu=size(SPH(il)%Ps,1)
					if (size(M2L) /= nu) then
						deallocate(M2L) ; allocate(M2L(nu))
					endif
					call T_LR(NL,N(il),il,k,SPH,rlm,M2L)
					acc=acc+sum(M2L*SPH(il)%in((ic-1)*nu+1:ic*nu)*SPH(il)%wp)
				else
					nanm=(N(il)+1)**2
					allocate(Yp(nanm), hn(0:N(il)), SJ(0:N(il)), DJ(0:N(il)), SY(0:N(il)), DY(0:N(il)))
					call car2sph(rlm,Psph,1)
					call Ynm_normalise_1(N(il),Psph(2),Psph(3),Yp)
					call SPHY(N(il),k*Psph(1),NM,SY,DY)
					if (k*Psph(1) > 3000.0_prec) then
						call sph_j1(N(il),k*Psph(1),SJ)
					else
						call SPHJ(N(il),k*Psph(1),NM,SJ,DJ)
					endif
					hn=SJ+ji*SY
					do nn=0,N(il)
						acc=acc+4.0_prec*pi*hn(nn)*sum(Yp(nn*nn+1:nn*nn+2*nn+1)*cBF(il)%a(nn*nn+1:nn*nn+2*nn+1,ic))
					end do
					deallocate(Yp, hn, SJ, DJ, SY, DY)
				endif
			else
				! elements poses a ce niveau : integration directe
				do j=1,octree(il)%cell(ic)%Nb_p
					no=octree(il)%cell(ic)%pts(j)
					if (no > NE) cycle             ! cibles separees (CHI)
					if (fmm_adaptatif) then
						if (niv_E(no) /= il) cycle
					else
						if (il /= NL) cycle
					endif
					PQ=transpose(Ps(E(no,:),:))
					Qv=matmul(PQ,F_v)
					call DMKH_T0(k,Qv,norc(no,:),Jac0(no),IS_v,PG_v(3,:),P,norP,D1,K1,M1,H1)
					acc=acc+D1*var_p(no)-M1*var_dnp(no)
					ndir=ndir+1
				end do
				if (il < NL) then
					do j=octree(il)%cell(ic)%Nb_c,1,-1
						top=top+1 ; pile_l(top)=il+1 ; pile_c(top)=octree(il)%cell(ic)%child(j)
					end do
				endif
			endif
		end do
		var_out(ipx)=acc
	end do
	if (.not. omp_in_parallel()) then
		n_prop_cellules=ncel ; n_prop_directs=ndir
	endif
	end subroutine propagation_hierarchique


! ------------------------------------------------------------------------
! RAM d'UN thread HRTF T0 (octets) pour le nombre d'onde k (2026-10-02), sans
! rien allouer : ordres par ordres_FMM_T0 (precision eps0_de_k(k)), puis les
! tableaux de fait_tout_HRTF_T0 et des produits FMM (comptage de
! test/test_memoire_t0.f08) :
!   FMM (toute la frequence) : moments V1+V2 (adaptatif 32 o, uniforme 40 o
!     x ncoef x elements du niveau), expansions W1 (16 / 24 o x ncoef x
!     cibles ; avec cibles separees, CHI, seuls les points > NE), developpements (2 x 16 o x ncoef x cellules), tables M2L HF
!     (16 o x points x 343) et BF (15 x 24 o x (N+1)^3, majorant)
!   puis le maximum des trois phases :
!     champ proche : Anear + Bnear + Apnear       48 o x num_K0
!     ILUT         : Apnear + Anear0 + jA0 (36 o x num_K0) + ALU/jLU
!     GMRES        : Apnear (16 o x num_K0) + ALU/jLU + work
!   ILUT : 20 o x coef_ilut x NE x (2 lfil + 1) ; GMRES : 16 o x lwork.
! liberer_V1 : V1 (16 o x ncoef x elements, adaptatif) retire des trois phases, garde a la propagation.
! ncoef = (N+1)^2 en BF, (N+1)(2N+1) en HF. Detail facultatif (octets) :
! det = (/ FMM, champ proche, ILUT, GMRES /).
! ------------------------------------------------------------------------
	real(prec) function memoire_thread_T0(octree,NL,k,HFBF,det) result(m)
	implicit none
	integer, intent(in) :: NL, HFBF
	type(tree), intent(in) :: octree(NL)
	real(prec), intent(in) :: k
	real(prec), intent(out), optional :: det(4)
	integer :: N(NL), BF(NL), il, ic, p
	integer(8) :: nel, ncib
	real(prec) :: nc, m_fmm, m_ilut, m_gm, K0, m_v1
	call ordres_FMM_T0(octree,NL,k,HFBF,eps0_de_k(k),N,BF)
	m_fmm = 0.0_prec ; m_v1 = 0.0_prec
	do il = 1, NL
		if (BF(il) == 1) then
			nc = real((N(il)+1)**2, prec)
		else
			nc = real((N(il)+1)*(2*N(il)+1), prec)
		endif
		if (fmm_adaptatif) then
			nel = count(niv_E == il)
			ncib = 0
			if (size(octree(il)%ipt) == NE) then     ! T0 : les NE barycentres sont les cibles
				do ic = 1, octree(il)%NL
					if (il == NL) then
						if (octree(il)%cell(ic)%a_cible) ncib = ncib + octree(il)%cell(ic)%Nb_p
					else
						if (octree(il)%cell(ic)%feuille) ncib = ncib + octree(il)%cell(ic)%Nb_p
					endif
				end do
			else                                      ! cibles separees (CHI, 2 NE points) : Nb_p compte
				do p = NE+1, size(octree(il)%ipt)     ! aussi les sources ; seuls les points > NE (2026-10-29)
					ic = octree(il)%ipt(p)
					if (il == NL) then
						if (octree(il)%cell(ic)%a_cible) ncib = ncib + 1
					else
						if (octree(il)%cell(ic)%feuille) ncib = ncib + 1
					endif
				end do
			endif
			m_fmm = m_fmm + 32.0_prec*nc*real(nel,prec) + 16.0_prec*nc*real(ncib,prec)
			m_v1 = m_v1 + 16.0_prec*nc*real(nel,prec)
		elseif (il == NL) then
			m_fmm = m_fmm + 40.0_prec*nc*real(NE,prec) + 24.0_prec*nc*real(NE,prec)
		endif
		m_fmm = m_fmm + 2.0_prec*16.0_prec*nc*real(octree(il)%NL,prec)
		if (BF(il) == 0) then
			m_fmm = m_fmm + 16.0_prec*nc*343.0_prec
		else
			m_fmm = m_fmm + 15.0_prec*24.0_prec*real(N(il)+1,prec)**3
		endif
	end do
	K0 = real(num_K0,prec)
	m_ilut = 20.0_prec*coef_ilut*real(NE,prec)*real(2*lfil_ilut+1,prec)
	m_gm = 16.0_prec*(real(restart,prec)**2 + real(restart,prec)*real(NE+5,prec) + 6.0_prec*real(NE,prec))
	m = m_fmm + max(48.0_prec*K0, 36.0_prec*K0 + m_ilut, 16.0_prec*K0 + m_ilut + m_gm)
	! liberer_V1 (2026-11-06) : V1 absent pendant champ proche, ILUT et GMRES, present a la propagation
	! (champ proche et GMRES liberes a ce moment) : maximum des deux situations
	if (liberer_V1 .and. fmm_adaptatif) m = max(m - m_v1, m_fmm)
	if (fusion_Vc .and. fmm_adaptatif) m = m - m_v1         ! Diff : V1 jamais stocke
	if (present(det)) det = (/ m_fmm, 48.0_prec*K0, m_ilut, m_gm /)
	end function memoire_thread_T0
	end module module_MLFMM
