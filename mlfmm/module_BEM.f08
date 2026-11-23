! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: MIT
	module module_BEM
! Convention des noms d'elements (2026-11-20) :
!   T0 : triangle plan, inconnue CONSTANTE sur l'element, collocation au barycentre
!        (HRTF : *_T0, ex. DMKH_T0, DMKH_sing_T0, DMKH_proche_T0 ; anciennement *_const_T3) ;
!   T3 : triangle plan, valeurs NODALES lineaires (tuyau : module_BEM_T3, module_MLFMM_T3) ;
!   T6 : triangle quadratique a 6 noeuds, valeurs nodales (tuyau : module_BEM_T6, module_MLFMM_T6).
	use module__Precision
	use omp_lib
	use module_octree
	use module_utils
!	use sparse_kit
!	use unique_sort_mod
	implicit none

! Maillage
	integer :: NE,NP,Next,NS,NE_prop
	integer, dimension(:,:), allocatable :: E,Eext
	integer, dimension(:), allocatable :: Ino,Ino1
!	integer, dimension(:), allocatable :: noP,noC
!	integer, dimension(:), allocatable :: iAp,jAp,noCp,noPp,iA,jA,Part,ind_Part,Identic
	integer, dimension(:), allocatable :: iAp,jAp,iA,jA,Part,ind_Part,Identic
	real(prec), dimension(:,:), allocatable :: Ps, Pext, norext,P_prop
	real(prec), dimension(:,:), allocatable :: Pc, norc, Pc1, norc1
	real(prec), dimension(:,:), allocatable :: Pt   ! points cibles (lignes) : Pc, ou interieurs (CHI)
	! CHI + octree adaptatif : octree sur 2 NE points, sources = barycentres Pc (1..NE), cibles =
	! points interieurs Pt (NE+1..2NE, ligne p-NE), chacun dans SA cellule (2026-10-14)
	logical :: cibles_sep = .false.
	real(prec), dimension(:,:), allocatable :: base_a1,base_a2
	real(prec), dimension(:), allocatable :: Jac0

	! Angle solide et terme libre -- voir angles_solides() pour les conventions.
	!   type_probleme : 0 = interieur (tube ferme/ferme), 1 = exterieur (tube ouvert, HRTF)
	!   somme_dG0     : oint d_nQ G0 dS, calcule par Angle_solide / Angle_solide_global_T6
	!   Omega_moins   : angle solide / 4pi cote du volume ferme par le maillage, dans [0,1]
	!   Omega_plus    : angle solide / 4pi cote de l'infini = 1 - Omega_moins
	!   terme_libre   : angle solide / 4pi cote FLUIDE, pose sur la diagonale de Dnear
	integer :: type_probleme = 0
	real(prec), dimension(:), allocatable :: somme_dG0, Omega_moins, Omega_plus, terme_libre
	integer :: num_K0,num_K1,num_M0,num_Kp,num_Mp
	integer :: num_K0_child  ! T3 : frontiere dans Ino/Ino1 : [1,num_K0_child]=same-cell (IS_c), [num_K0_child+1,num_K0]=voisin (IS_v)
	! T6 SANS DOUBLONS (2026-09-28) : iA/jA/ind_Part/num_K1 decrivent les couples UNIQUES
	! (ligne, colonne, partition) ; pos_near(num_K1_brut) donne, pour chaque contribution
	! elementaire (paire (i1,jn) -> (i1-1)*6+jn), l'indice du couple unique ou l'accumuler.
	integer, dimension(:), allocatable :: pos_near
	integer :: num_K1_brut = 0
	! Anear (double couche) ne depend pas de la partition : il est stocke sur les
	! couples uniques (ligne, colonne) au format CSR (iAc(NP+1), jAc(num_Ka)), lu
	! tel quel par l'ILUT. posA_near(num_K1_brut) : contribution brute -> indice CSR.
	integer, dimension(:), allocatable :: iAc, jAc, posA_near
	integer :: num_Ka = 0
	! octree ADAPTATIF (2026-09-29) : niveau de chaque element (indice il de
	! octree(il), 1..Nlev), donne par niveaux_elements (module_Octree)
	integer, dimension(:), allocatable :: niv_E
	! FMM adaptative active (tuyau : ligne NL de la config a -1) : les produits
	! MLFMM_HYB_(p_)dnp_T6 passent par leur version _adapt
	logical :: fmm_adaptatif = .false.
	! T0 adaptatif (HRTF) : type de chaque paire (iA cible, jA source) du champ
	! proche, construit par init_Anear_T0_adapt : 0 singuliere (cible = source),
	! 1 meme cellule (quadrature IS_c), 2 cellule voisine (IS_v)
	integer, dimension(:), allocatable :: type_paire_T0

	! impédance en une zone (T3 tuyau : IZ = indices near-field de la zone impedante)
	integer, dimension(:), allocatable :: IZ,IP
	integer :: num_Z0, num_P0
	complex(prec), dimension(:), allocatable :: PZ

!	integer, parameter :: IS_c=7, IS_v=16, IS_s=9, coef=200, coef_2=2, IS_f=16, IS_e=16
	integer, parameter :: IS_c=33, IS_v=33, IS_s=12, coef=200, coef_2=2, IS_f=16, IS_e=16
	! IS_s3 : nb de points de la quadrature singuliere T3 (auto-terme simple couche, DM_sing_T3).
	! DOIT etre un carre parfait car pgc4_gauss remplit IS05**2 = int(sqrt(IS_s3))**2 points
	! (16 = 4x4). NE PAS confondre avec IS_s (=12) utilise par gauleg dans le near-field T0/HRTF.
	integer, parameter :: IS_s3=49

	! partitionnement du maillage
	Type Bord
		integer, dimension(:),allocatable :: curve
		integer	:: Nb_c
	end type Bord
	Type BEM
!        type(Bord), allocatable :: Lim(:)
        integer, allocatable :: Lim(:)
		integer, dimension(:,:), allocatable :: E,Bord
		real(prec), allocatable :: Ps(:,:)
		integer, dimension(:), allocatable :: table_E,table_Ps,NB0
		integer	:: NE,NP,NL,NB
	end type BEM

	type(BEM), allocatable :: mesh(:)

	contains

!--------------------------------------------------------------------
! Angles solides et terme libre, deduits de somme_dG0.
!
! CONVENTION : la normale SORT TOUJOURS DU DOMAINE FLUIDE.
!   - interieur (tube ferme/ferme) : le fluide est DANS le volume ferme
!     par le maillage ; normale sortante de ce volume.
!   - exterieur (tube ferme/ouvert, HRTF) : le fluide est AUTOUR ;
!     normale dirigee vers la structure, donc entrante dans le volume.
! C'est l'inverse de formulation_isoparametrique.pdf, qui ne sert donc
! pas de reference ici.
!
! Angles solides divises par 4pi, toujours POSITIFS :
!   Omega_moins = fraction de sphere vue cote volume ferme
!   Omega_plus  = fraction vue cote infini,   Omega_moins + Omega_plus = 1
! (1/2 chacun sur une surface lisse).
!
! Par Gauss sur le volume ferme :
!   normale sortante du volume (interieur) : somme_dG0 = -Omega_moins
!   normale entrante (exterieur)           : somme_dG0 = +Omega_moins
!
! La representation integrale, normale sortant du fluide, est la meme
! dans les deux cas :   c p = PV oint (G dnp - p dnG)   d'ou
!       (D + c) p = S dnp,   c = angle solide cote FLUIDE
!   interieur : c = Omega_moins      exterieur : c = Omega_plus
!
! CONTROLE : Omega_moins hors de [0,1] => normales dans la mauvaise
! convention pour ce type de probleme -> arret (rattrapage sans
! remailler : test/msh_flip.py).
!--------------------------------------------------------------------
! couples_uniques_near(nne) : structure creuse du champ proche SANS
! DOUBLONS, commune au T3 (nne = 3) et au T6 (nne = 6). 2026-09-29.
!
! Entree : les paires (element Ino(i1), point de collocation Ino1(i1)),
! i1 = 1..num_K0, construites par init_Anear_part_T3/T6_CPU_v3. Chaque
! paire donne nne contributions elementaires, une par noeud de l'element :
! la contribution brute (i1, jn) porte le numero (i1-1)*nne + jn.
!
! Un meme couple (ligne no1, colonne j) revient autant de fois que j
! appartient a des elements du champ proche de no1 (~2,5 fois en T6, les
! milieux d'arete appartenant a 2 elements et les sommets a ~6). On ne
! stocke que les couples UNIQUES :
!   - Bnear (simple couche, multiplie dp/dn, qui differe d'une partition a
!     l'autre sur les jantes) : cle (ligne, colonne, partition)
!       -> iA, jA, ind_Part, num_K1 ; pos_near(brut) = indice du couple
!   - Anear (double couche, multiplie p, continu) : cle (ligne, colonne),
!     au format CSR lu tel quel par l'ILUT (colonnes non triees : l'ILUT
!     n'en a pas besoin, mais les exige UNIQUES)
!       -> iAc(NP+1), jAc, num_Ka ; posA_near(brut) = indice CSR
!   - Identic(no1) = indice CSR du couple (no1, no1), ou poser le terme libre
! L'assemblage (DMnear_FMM_part_T3/T6_v3) accumule chaque contribution
! brute a sa place : Dnear(posA_near(r)) += ..., Mnear(pos_near(r)) += ...
!
! Dedoublonnage ligne par ligne en temps lineaire, marqueur (colonne,
! partition) prive a chaque thread ; ordre de premiere apparition.
!--------------------------------------------------------------------
	subroutine couples_uniques_near(nne)
	implicit none
	integer, intent(in) :: nne
	integer :: i1, no, no1, jn, loc, i, r, raw, j, pp, c, ca, u, npart
	integer, allocatable :: rp(:), ridx(:), fill(:), nu(:), up(:), nua(:), mark(:,:), markA(:)

	if (allocated(iA))        deallocate(iA)
	if (allocated(jA))        deallocate(jA)
	if (allocated(ind_Part))  deallocate(ind_Part)
	if (allocated(Identic))   deallocate(Identic)
	if (allocated(pos_near))  deallocate(pos_near)
	if (allocated(posA_near)) deallocate(posA_near)
	if (allocated(iAc))       deallocate(iAc)
	if (allocated(jAc))       deallocate(jAc)
	num_K1_brut = nne*num_K0
	npart = maxval(Part(1:NE))

	! contributions brutes regroupees par ligne (point de collocation no1)
	allocate(rp(NP+1), fill(NP), ridx(num_K1_brut))
	rp = 0
	do i1 = 1, num_K0
		rp(Ino1(i1)+1) = rp(Ino1(i1)+1) + nne
	end do
	rp(1) = 1
	do i = 1, NP
		rp(i+1) = rp(i) + rp(i+1)
	end do
	fill = rp(1:NP)
	do i1 = 1, num_K0
		no1 = Ino1(i1)
		do jn = 1, nne
			ridx(fill(no1)) = (i1-1)*nne + jn
			fill(no1) = fill(no1) + 1
		end do
	end do
	deallocate(fill)

	! passe 1 : nombre de couples uniques par ligne, pour Bnear (colonne,
	! partition) et pour Anear (colonne)
	allocate(nu(NP), up(NP+1), nua(NP), iAc(NP+1))
	!$omp parallel private(i,r,raw,i1,jn,no,j,pp,c,ca,mark,markA)
	allocate(mark(NP,npart), markA(NP)) ; mark = 0 ; markA = 0
	!$omp do schedule(dynamic,64)
	do i = 1, NP
		c = 0 ; ca = 0
		do r = rp(i), rp(i+1)-1
			raw = ridx(r) ; i1 = (raw-1)/nne + 1 ; jn = raw - (i1-1)*nne
			no = Ino(i1) ; j = E(no,jn) ; pp = Part(no)
			if (mark(j,pp) == 0) then
				c = c + 1 ; mark(j,pp) = c
			endif
			if (markA(j) == 0) then
				ca = ca + 1 ; markA(j) = ca
			endif
		end do
		nu(i) = c ; nua(i) = ca
		do r = rp(i), rp(i+1)-1
			raw = ridx(r) ; i1 = (raw-1)/nne + 1 ; jn = raw - (i1-1)*nne
			no = Ino(i1) ; mark(E(no,jn),Part(no)) = 0 ; markA(E(no,jn)) = 0
		end do
	end do
	!$omp end do
	deallocate(mark, markA)
	!$omp end parallel
	up(1) = 1 ; iAc(1) = 1
	do i = 1, NP
		up(i+1)  = up(i)  + nu(i)
		iAc(i+1) = iAc(i) + nua(i)
	end do
	num_K1 = up(NP+1) - 1
	num_Ka = iAc(NP+1) - 1

	! passe 2 : couples uniques et correspondance brute -> unique
	allocate(iA(num_K1), jA(num_K1), ind_Part(num_K1), pos_near(num_K1_brut), Identic(NP))
	allocate(jAc(num_Ka), posA_near(num_K1_brut))
	!$omp parallel private(i,r,raw,i1,jn,no,j,pp,c,ca,u,mark,markA)
	allocate(mark(NP,npart), markA(NP)) ; mark = 0 ; markA = 0
	!$omp do schedule(dynamic,64)
	do i = 1, NP
		c = 0 ; ca = 0
		do r = rp(i), rp(i+1)-1
			raw = ridx(r) ; i1 = (raw-1)/nne + 1 ; jn = raw - (i1-1)*nne
			no = Ino(i1) ; j = E(no,jn) ; pp = Part(no)
			if (mark(j,pp) == 0) then
				c = c + 1 ; mark(j,pp) = c
				u = up(i) + c - 1
				iA(u) = i ; jA(u) = j ; ind_Part(u) = pp
			endif
			pos_near(raw) = up(i) + mark(j,pp) - 1
			if (markA(j) == 0) then
				ca = ca + 1 ; markA(j) = ca
				jAc(iAc(i) + ca - 1) = j
			endif
			posA_near(raw) = iAc(i) + markA(j) - 1
		end do
		do r = rp(i), rp(i+1)-1
			raw = ridx(r) ; i1 = (raw-1)/nne + 1 ; jn = raw - (i1-1)*nne
			no = Ino(i1) ; mark(E(no,jn),Part(no)) = 0 ; markA(E(no,jn)) = 0
		end do
	end do
	!$omp end do
	deallocate(mark, markA)
	!$omp end parallel
	deallocate(rp, ridx, nu, up, nua)

	! terme libre : couple (no1, no1) dans Anear (CSR). Identic(no1) = 0 si aucun
	! element de no1 ne le voit dans son champ proche (octree trop profond : la
	! cellule est plus petite que les elements) ; DMnear_* s'arrete alors.
	Identic = 0
	do i1 = 1, num_K0
		no  = Ino(i1)
		no1 = Ino1(i1)
		loc = 0
		do jn = 1, nne
			if (E(no,jn) == no1) then
				loc = jn ; exit
			endif
		end do
		if (loc /= 0) Identic(no1) = posA_near((i1-1)*nne + loc)
	enddo
	if (any(Identic == 0)) write(*,'(A,I8,A)') ' ATTENTION couples_uniques_near : ', count(Identic == 0), &
	     ' noeud(s) sans auto-terme (cellule plus petite que les elements : NL trop grand)'
	write(*,'(A,I2,A,I11,A,I11,A,I11)') ' champ proche T', nne, ' sans doublon : contributions brutes', &
	     num_K1_brut, ' -> Bnear', num_K1, ' , Anear CSR', num_Ka
	end subroutine couples_uniques_near

!--------------------------------------------------------------------
! init_Anear_adapt(octree,Nlev,nne) : paires (element, noeud) du champ
! proche pour l'octree ADAPTATIF (2026-09-29), T3 (nne=3) et T6 (nne=6).
!
! L'element no est vu a SON niveau il = niv_E(no), dans sa cellule
! ic = octree(il)%ipt(no) :
!   bloc child  (quadrature IS_c) : noeuds de ic
!   bloc voisin (quadrature IS_v) : noeuds des voisines de ic au niveau il
! Meme rangement que init_Anear_part_T3/T6_CPU_v3 : [1,num_K0_child] child,
! [num_K0_child+1,num_K0] voisin. Puis couples_uniques_near(nne) : tout
! l'aval (DMnear, ILUT) est inchange.
!
! Parcours : niveaux il croissants -> cellules -> pts (elements de niveau
! il). Au niveau le plus fin pts = child (completer_octree) : avec
! niv_E(:) = Nlev on retrouve Ino/Ino1 de init_Anear_part_T6_CPU_v3 dans
! le MEME ORDRE, donc Anear/Bnear au bit pres.
! Necessite completer_octree et niv_E.
!--------------------------------------------------------------------
	subroutine init_Anear_adapt(octree,Nlev,nne)
	implicit none
	integer, intent(in) :: Nlev, nne
	type(tree), intent(in) :: octree(Nlev)
	integer :: il, ic, i2, i3, i4, cv, no, no1
	integer :: num_K, num_K_child, num_K_voisin, idx_c, idx_v
	real :: t01, t02

	t01 = omp_get_wtime()
	num_K_child = 0 ; num_K_voisin = 0
	!$OMP PARALLEL DO private(il,ic,i3,i4,cv) shared(octree,niv_E,NE) &
	!$OMP REDUCTION(+:num_K_child,num_K_voisin)
	do no=1,NE
		il = niv_E(no)
		ic = octree(il)%ipt(no)
		do i3=1,octree(il)%cell(ic)%Nb_p
			if (octree(il)%cell(ic)%pts(i3) > NE) num_K_child = num_K_child+1
		enddo
		do i3=1,octree(il)%cell(ic)%Nb_v
			cv = octree(il)%cell(ic)%voisin(i3)
			do i4=1,octree(il)%cell(cv)%Nb_p
				if (octree(il)%cell(cv)%pts(i4) > NE) num_K_voisin = num_K_voisin+1
			end do
		end do
	enddo
	!$OMP END PARALLEL DO
	num_K = num_K_child + num_K_voisin
	write(*,*) 'pass 1 adaptatif : num_K=',num_K,' (child=',num_K_child,' voisin=',num_K_voisin,')'
	num_K0 = num_K
	num_K0_child = num_K_child
	if (allocated(Ino))  deallocate(Ino)
	if (allocated(Ino1)) deallocate(Ino1)
	allocate(Ino(num_K),Ino1(num_K))

	idx_c = 0 ; idx_v = num_K_child
	do il=1,Nlev
		do ic=1,octree(il)%NL
			do i2=1,octree(il)%cell(ic)%Nb_p
				no = octree(il)%cell(ic)%pts(i2)
				if (no > NE) cycle
				if (niv_E(no) /= il) cycle
				do i3=1,octree(il)%cell(ic)%Nb_p
					no1 = octree(il)%cell(ic)%pts(i3)
					if (no1 > NE) then
						idx_c = idx_c+1
						Ino1(idx_c) = no1-NE ; Ino(idx_c) = no
					endif
				enddo
				do i3=1,octree(il)%cell(ic)%Nb_v
					cv = octree(il)%cell(ic)%voisin(i3)
					do i4=1,octree(il)%cell(cv)%Nb_p
						no1 = octree(il)%cell(cv)%pts(i4)
						if (no1 > NE) then
							idx_v = idx_v+1
							Ino1(idx_v) = no1-NE ; Ino(idx_v) = no
						endif
					end do
				end do
			enddo
		enddo
	enddo

	call couples_uniques_near(nne)
	t02 = omp_get_wtime()
	write(*,*) 'init_Anear adaptatif en',t02-t01,' sec ; num_K1=',num_K1, &
	           ' couples uniques (contributions brutes',num_K1_brut,') ; Anear CSR',num_Ka

	end subroutine init_Anear_adapt

!--------------------------------------------------------------------
	subroutine angles_solides(arret)
	implicit none
	logical, intent(in), optional :: arret   ! .false. : signale sans arreter (diagnostic)
	real(prec), parameter :: tol = 1.0e-2_prec
	integer :: nfaux

	if (allocated(Omega_moins)) deallocate(Omega_moins)
	if (allocated(Omega_plus))  deallocate(Omega_plus)
	if (allocated(terme_libre)) deallocate(terme_libre)
	allocate(Omega_moins(NP), Omega_plus(NP), terme_libre(NP))

	select case (type_probleme)
	case (0)
		Omega_moins = -somme_dG0
		Omega_plus  = 1.0_prec - Omega_moins
		terme_libre = Omega_moins
	case (1)
		Omega_moins = somme_dG0
		Omega_plus  = 1.0_prec - Omega_moins
		terme_libre = Omega_plus
	case default
		write(*,*) 'angles_solides : type_probleme inconnu (', type_probleme, ') : 0 ou 1'
		stop 1
	end select

	write(*,'(A,F8.4,A,F8.4)') ' Omega-/4pi min=', minval(Omega_moins), ' max=', maxval(Omega_moins)
	write(*,'(A,F8.4,A,F8.4)') ' Omega+/4pi min=', minval(Omega_plus),  ' max=', maxval(Omega_plus)

	nfaux = count(Omega_moins < -tol .or. Omega_moins > 1.0_prec + tol)
	if (nfaux > 0) then
		write(*,'(A,I8,A,I8,A)') ' ERREUR : ', nfaux, ' noeuds sur ', NP, ' ont Omega-/4pi hors de [0,1]'
		if (type_probleme == 0) then
			write(*,*) ' Probleme INTERIEUR : la normale doit SORTIR du volume ferme (du fluide).'
		else
			write(*,*) ' Probleme EXTERIEUR : la normale doit pointer VERS la structure.'
		endif
		write(*,*) ' Verifier le type de probleme en tete de config, ou retourner le maillage (test/msh_flip.py).'
		if (present(arret)) then
			if (.not. arret) return
		endif
		stop 1
	endif

	end subroutine angles_solides

	subroutine pre_BEM_T0()
!-----------------------------------------------------------------------------------------------
!
! pré-calcul BEM après récupération du maillage GMSH
! Calculs les points centraux, normales correspondantes et Jacobiens de chaque éléments T0
!
! base_a1 et base_a2
! Pc et norc
! Jac0
! A faire avant octree
!
!-----------------------------------------------------------------------------------------------
!
	implicit none
	integer :: i1,nPe(3)
	real(prec) :: a1(3), a2(3), Nor(3),PQ(3,3)


	allocate(Pc(NE,3))
	allocate(Jac0(1:NE))
	allocate(base_a1(NE,3))  !
	allocate(base_a2(NE,3))  !
	allocate(norc(NE,3))

	do i1=1,NE
		nPe=E(i1,:) ! Points de l'élément
		PQ=transpose(Ps(nPe,:)) ! Coordonnées des points de l'éléments
		Pc(i1,:)=matmul(PQ,(/ 1.0_prec/3.0_prec,1.0_prec/3.0_prec,1.0_prec/3.0_prec /))
		a1=matmul(PQ,(/-1,1,0/))
		a2=matmul(PQ,(/-1,0,1/))
		base_a1(i1,:)=a1
		base_a2(i1,:)=a2
		Nor(1)=a1(2)*a2(3)-a1(3)*a2(2)
		Nor(2)=a1(3)*a2(1)-a1(1)*a2(3)
		Nor(3)=a1(1)*a2(2)-a1(2)*a2(1)
		Jac0(i1)=sqrt(Nor(1)**2+Nor(2)**2+Nor(3)**2)
		norc(i1,:)=Nor/Jac0(i1)
	enddo
	end subroutine pre_BEM_T0

!
!-----------------------------------------------------------------------------------------------
! Init_Anear
! Initialisation de base pour calculer la matrice creuse Anear, définition de la taille
! sans partitionnement, en T0
! à faire après octree
!
!-----------------------------------------------------------------------------------------------
!
	subroutine init_Anear(num_K,octree,NL)
!
! num_K (out) : taille de la matrice creuse Anear/Bnear
! octree(NL)
!
	implicit none

! Variables d'entrée
	integer, intent(in) :: NL
	type(tree), intent(in) :: octree(NL)
! Variables de sortie
	integer, intent(out) :: num_K

! Variables locales

	integer :: i1, i2, i3

! init
	num_K=0
	do i1=1,octree(NL)%NL   ! boucle sur cellule dernier niveau
		do i2=1,octree(NL)%cell(i1)%Nb_c
		  num_K=num_K+octree(NL)%cell(i1)%Nb_c
  		  do i3=1,octree(NL)%cell(i1)%Nb_v
			num_K=num_K+octree(NL)%cell(octree(NL)%cell(i1)%voisin(i3))%Nb_c
		  end do
		end do
	end do

	allocate(iA(num_K0))
	allocate(jA(num_K0))

	end subroutine init_Anear
!--------------------------------------------------------------------
! Init_Anear_Z
! Initialisation pour matrice creuse Anear avec impédance
!
!--------------------------------------------------------------------

	subroutine init_Anear_Z(octree,NL)
!	calcul la matrice creuse Anear pour FMM
! et aussi pour le calcul des moments
	use parametres
	implicit none

! Variables d'entrée
	integer, intent(in) :: NL
	type(tree), intent(in) :: octree(NL)
! Variables de sortie
	integer :: num_K

! Variables locales

	integer :: i1, i2, i3, i4

! init
	num_K=0
    do i1=1,octree(NL)%NL   ! boucle sur cellule dernier niveau

        do i2=1,octree(NL)%cell(i1)%Nb_c
		  num_K=num_K+octree(NL)%cell(i1)%Nb_c

  		  do i3=1,octree(NL)%cell(i1)%Nb_v
			num_K=num_K+octree(NL)%cell(octree(NL)%cell(i1)%voisin(i3))%Nb_c
		  end do
		end do
	end do

	num_K0=num_K
	num_K=0
!	num_M=0
	allocate(iA(num_K0))
	allocate(jA(num_K0))
	allocate(PZ(num_K0))
	PZ=0.0_prec

	do i1=1,octree(NL)%NL   ! boucle sur cellule dernier niveau

		do i2=1,octree(NL)%cell(i1)%Nb_c

		! points dans la cellule
		  do i3=1,octree(NL)%cell(i1)%Nb_c

				num_K=num_K+1
		        iA(num_K)=octree(NL)%cell(i1)%child(i3)
                jA(num_K)=octree(NL)%cell(i1)%child(i2)
                if (Part(jA(num_K)) == 3) then
					PZ(num_K)=1.0_prec
				endif
          end do

		  ! points dans les cellules voisines
  		  do i3=1,octree(NL)%cell(i1)%Nb_v

			do i4=1,octree(NL)%cell(octree(NL)%cell(i1)%voisin(i3))%Nb_c
				num_K=num_K+1
				iA(num_K)=octree(NL)%cell(octree(NL)%cell(i1)%voisin(i3))%child(i4)
				jA(num_K)=octree(NL)%cell(i1)%child(i2)

				if (Part(jA(num_K)) == 3) then
					PZ(num_K)=1.0_prec
				endif
			end do
		  end do
		end do
	end do

	end subroutine init_Anear_Z
!
!--------------------------------------------------------------------------------------------
!
! Calcul Anear et Bnear pour FMM, BEM méthode Burton/Miller
!
!--------------------------------------------------------------------------------------------
!
	subroutine ABnear_FMM_T0(Anear,Bnear,octree,NL,k,alpha)
!	calcul la matrice creuse Anear et Bnear pour FMM
	use module_octree
	use parametres
	implicit none

! Variables d'entrée
	integer, intent(in) :: NL
	real(prec), intent(in) :: k
	complex(prec), intent(in) :: alpha
	type(tree), intent(in) :: octree(NL)
! Variables de sortie
	complex(prec), intent(inout) :: Anear(num_K0),Bnear(num_K0)
! Variables locales
!	integer, parameter :: IS_c=28, IS_v=28, IS_s=12, coef=200
	integer :: no,num_K,no1
	integer, dimension(1:3) :: nPe
	real(prec), dimension(1:3,1:3) :: PQ
	real(prec), dimension(1:3,1:IS_c) :: PG_c
	real(prec), dimension(1:3,1:IS_v) :: PG_v
	real(prec), dimension(1:2,1:IS_s) :: PG_s
	real(prec), dimension(1:IS_s) :: x, wth
	real(prec), dimension(1:3,1:IS_c) :: F_c, dF1_c, dF2_c
    real(prec), dimension(1:3,1:IS_v) :: F_v, dF1_v, dF2_v
	real(prec) :: Qc(1:3,IS_c), Qv(1:3,IS_v)
	complex(prec) :: D1, K1, M1, H1
	integer :: i1, i2, i3, i4
	if (fmm_adaptatif) then                ! octree adaptatif : liste de init_Anear_T0_adapt
		call ABnear_FMM_T0_adapt(Anear,Bnear,k,alpha)
		return
	endif
	call verifier_points_cibles()
	call pgt6(PG_c,IS_c)  ! points dans la même cellule
	call pgt6(PG_v,IS_v)  ! points dans la cellule voisine
	call gauleg(IS_s, x, wth)
	PG_s(1,:)=x		! points singuliers
	PG_s(2,:)=wth
	call FORM_T3(PG_c(1:2,:),IS_c,F_c,dF1_c,dF2_c)
	call FORM_T3(PG_v(1:2,:),IS_v,F_v,dF1_v,dF2_v)

! init
	num_K=0
	do i1=1,octree(NL)%NL   ! boucle sur cellule dernier niveau
		do i2=1,octree(NL)%cell(i1)%Nb_c
		  no=octree(NL)%cell(i1)%child(i2)
		  nPe=E(no,:) ! Points de l'élément
		  PQ=transpose(Ps(nPe,:)) ! Coordonnées des points de l'éléments
		  Qc=matmul(PQ,F_c)
		  Qv=matmul(PQ,F_v)
! points dans la cellule
		  do i3=1,octree(NL)%cell(i1)%Nb_c
			no1=octree(NL)%cell(i1)%child(i3)
			if (i2.eq.i3) then    ! point singulier
					num_K=num_K+1
					iA(num_K)=octree(NL)%cell(i1)%child(i3)
					jA(num_K)=octree(NL)%cell(i1)%child(i2)
				if (chi_d > 0.0_prec) then   ! CHI : point interieur, pas de saut 1/2
					call DMKH_proche_T0(k,PQ,norc(no,:),Pt(no1,:),norc(no1,:),D1,K1,M1,H1)
					Anear(num_K)=D1+alpha*K1
					Bnear(num_K)=M1+alpha*H1
				else
				call DMKH_sing_T0(no1,PQ,k,PG_s,IS_s,K1,M1)
					Anear(num_K)=alpha*K1-cmplx(0.5_prec,0.0_prec,kind=prec) ! prise en compte de l'angle solide+ BM
					Bnear(num_K)=M1+alpha*cmplx(0.5_prec,0.0_prec,kind=prec) ! prise en compte de l'angle solide+ BM
				endif
			else			! point non singulier
				call DMKH_T0(k,Qc,norc(no,:),Jac0(no),IS_c,PG_c(3,:),Pt(no1,:),norc(no1,:),D1,K1,M1,H1)
				num_K=num_K+1
				iA(num_K)=octree(NL)%cell(i1)%child(i3)
				jA(num_K)=octree(NL)%cell(i1)%child(i2)
				Anear(num_K)=D1+alpha*K1
				Bnear(num_K)=M1+alpha*H1
			endif
		  end do
! points dans les cellules voisines
  		  do i3=1,octree(NL)%cell(i1)%Nb_v
			do i4=1,octree(NL)%cell(octree(NL)%cell(i1)%voisin(i3))%Nb_c
				no1=octree(NL)%cell(octree(NL)%cell(i1)%voisin(i3))%child(i4)
				call DMKH_T0(k,Qv,norc(no,:),Jac0(no),IS_v,PG_v(3,:),Pt(no1,:),norc(no1,:),D1,K1,M1,H1)
				num_K=num_K+1
				iA(num_K)=octree(NL)%cell(octree(NL)%cell(i1)%voisin(i3))%child(i4)
				jA(num_K)=octree(NL)%cell(i1)%child(i2)
				Anear(num_K)=D1+alpha*K1
				Bnear(num_K)=M1+alpha*H1
			end do
		  end do
		end do
	end do
	end subroutine ABnear_FMM_T0
!
!--------------------------------------------------------------------------------------------
!
! OCTREE ADAPTATIF T0 (HRTF), 2026-09-29. Octree construit sur les seuls
! barycentres (make_octree_all(NE,Pc,...)) : chaque point est a la fois
! source (element) et cible (collocation au barycentre).
!
! init_Anear_T0_adapt : paires (cible iA, source jA) du champ proche, prises
! AU NIVEAU DE L'ELEMENT SOURCE il = niv_E(no), dans sa cellule
! ic = octree(il)%ipt(no) : cibles de ic (type 0 si cible = source, 1
! sinon) puis cibles des voisines de ic au niveau il (type 2). Rangement par
! element, dans l'ordre de ABnear_FMM_T0 (cellules -> elements -> cibles de
! la cellule -> cibles des voisines) : avec niv_E(:) = Nlev on retrouve ses
! iA/jA dans le MEME ORDRE. Construit une fois (hors boucle en frequence).
! Necessite completer_octree et niv_E.
!
!--------------------------------------------------------------------------------------------
!
	subroutine init_Anear_T0_adapt(octree,Nlev)
	implicit none
	integer, intent(in) :: Nlev
	type(tree), intent(in) :: octree(Nlev)
	integer :: il, ic, i2, i3, i4, cv, no, num_K, p, l

	! l = 1 : comptage, l = 2 : remplissage. cibles_sep : seules les cibles p > NE (ligne p - NE)
	do l=1,2
	num_K = 0
	do il=1,Nlev
		do ic=1,octree(il)%NL
			do i2=1,octree(il)%cell(ic)%Nb_p
				no = octree(il)%cell(ic)%pts(i2)
				if (no > NE) cycle
				if (niv_E(no) /= il) cycle
				do i3=1,octree(il)%cell(ic)%Nb_p          ! cibles de la cellule
					p = octree(il)%cell(ic)%pts(i3)
					if (cibles_sep) then
						if (p <= NE) cycle
						p = p - NE
					endif
					num_K = num_K+1
					if (l == 1) cycle
					iA(num_K) = p
					jA(num_K) = no
					if (p == no) then
						type_paire_T0(num_K) = 0
					else
						type_paire_T0(num_K) = 1
					endif
				enddo
				do i3=1,octree(il)%cell(ic)%Nb_v          ! cibles des voisines
					cv = octree(il)%cell(ic)%voisin(i3)
					do i4=1,octree(il)%cell(cv)%Nb_p
						p = octree(il)%cell(cv)%pts(i4)
						if (cibles_sep) then
							if (p <= NE) cycle
							p = p - NE
						endif
						num_K = num_K+1
						if (l == 1) cycle
						iA(num_K) = p
						jA(num_K) = no
						type_paire_T0(num_K) = 2
						if (p == no) type_paire_T0(num_K) = 0   ! cibles separees : propre dans une voisine
					enddo
				enddo
			enddo
		enddo
	enddo
	if (l == 1) then
		num_K0 = num_K
		if (allocated(iA)) deallocate(iA)
		if (allocated(jA)) deallocate(jA)
		if (allocated(type_paire_T0)) deallocate(type_paire_T0)
		allocate(iA(num_K0), jA(num_K0), type_paire_T0(num_K0))
	endif
	enddo
	if (count(type_paire_T0 == 0) /= NE) then
		write(*,*) 'init_Anear_T0_adapt : ', count(type_paire_T0 == 0), ' termes propres pour ', NE, ' elements'
		stop 'init_Anear_T0_adapt : point cible hors du voisinage de son element (CHI trop profond ?)'
	endif
	end subroutine init_Anear_T0_adapt
!
!--------------------------------------------------------------------------------------------
!
! ABnear_FMM_T0_adapt : coefficients Burton/Miller du champ proche T0 sur la
! liste de init_Anear_T0_adapt (iA, jA, type_paire_T0). Memes noyaux, memes
! quadratures et memes formules que ABnear_FMM_T0 (qui, elle, reconstruit
! iA/jA a chaque appel) : resultat identique au bit pres en uniforme.
!
!--------------------------------------------------------------------------------------------
!
	subroutine ABnear_FMM_T0_adapt(Anear,Bnear,k,alpha)
	use parametres
	implicit none
	real(prec), intent(in) :: k
	complex(prec), intent(in) :: alpha
	complex(prec), intent(inout) :: Anear(num_K0),Bnear(num_K0)
	integer :: i, no, no1, no_prec
	integer, dimension(1:3) :: nPe
	real(prec), dimension(1:3,1:3) :: PQ
	real(prec), dimension(1:3,1:IS_c) :: PG_c
	real(prec), dimension(1:3,1:IS_v) :: PG_v
	real(prec), dimension(1:2,1:IS_s) :: PG_s
	real(prec), dimension(1:IS_s) :: x, wth
	real(prec), dimension(1:3,1:IS_c) :: F_c, dF1_c, dF2_c
	real(prec), dimension(1:3,1:IS_v) :: F_v, dF1_v, dF2_v
	real(prec) :: Qc(1:3,IS_c), Qv(1:3,IS_v)
	complex(prec) :: D1, K1, M1, H1
	! quadrature selon la distance (quad_proche = 1, 2026-10-30) : regle choisie par
	! S = Hmax/D (Hmax plus grand cote de l'element, D distance cible - centre) et k Hmax,
	! comme l'indice de severite de Rezayat, Shippy & Rizzo (CMAME 55, 1986) ; seuils
	! calibres sur m3 (quadrature_20261030/calibre.py) pour une erreur normalisee <= 1e-7.
	! Essai 2026-11-02 (bm_quad_20261102) : Burton-Miller avec les couples proches en polaires +
	! sinh 16 x 16 (quasi exact, DMKH_proche_T0) : p99 a 469 Hz -7 a -23 % (pire sur u1,5), champ
	! proche x8 a x17 : l'erreur de BM vient de la formulation, pas de cette regle.
	! Essai 2026-11-02 (quadbm_20261102) : 48 points sur tous les couples = cette regle a 4e-7 dB pres
	! (sphere g1,5, BM, 2,1 / 5,2 / 9,4 kHz) ; 33 partout = quad_proche 0. Au-dessus de 3 kHz, cote oppose,
	! l'ancienne regle a 33 points compensait en partie l'erreur de BM (p99 0,08 au lieu de 0,45 dB).
	integer, parameter :: NR = 6
	integer, parameter :: nr_pts(NR) = [7, 16, 19, 27, 33, 48]
	real(prec) :: PGr(3,48,NR), Fr(3,48,NR), dFr1(3,48), dFr2(3,48), Qr(3,48,NR)
	real(prec) :: Hmax, S, kH
	integer :: r, compte(NR)

	call verifier_points_cibles()
	PGr = 0.0_prec ; Fr = 0.0_prec ; compte = 0
	do r = 1, NR
		call pgt6(PGr(:,1:nr_pts(r),r), nr_pts(r))
		call FORM_T3(PGr(1:2,1:nr_pts(r),r), nr_pts(r), Fr(:,1:nr_pts(r),r), dFr1(:,1:nr_pts(r)), dFr2(:,1:nr_pts(r)))
	enddo
	call pgt6(PG_c,IS_c)
	call pgt6(PG_v,IS_v)
	call gauleg(IS_s, x, wth)
	PG_s(1,:)=x
	PG_s(2,:)=wth
	call FORM_T3(PG_c(1:2,:),IS_c,F_c,dF1_c,dF2_c)
	call FORM_T3(PG_v(1:2,:),IS_v,F_v,dF1_v,dF2_v)

	no_prec=0
	do i=1,num_K0
		no=jA(i)
		if (no /= no_prec) then        ! paires rangees par element source
			nPe=E(no,:)
			PQ=transpose(Ps(nPe,:))
			if (quad_proche == 1) then
				do r = 1, NR
					Qr(:,1:nr_pts(r),r) = matmul(PQ, Fr(:,1:nr_pts(r),r))
				enddo
				Hmax = max(norm2(PQ(:,2)-PQ(:,1)), norm2(PQ(:,3)-PQ(:,2)), norm2(PQ(:,1)-PQ(:,3)))
				kH = k*Hmax
			else
				Qc=matmul(PQ,F_c)
				Qv=matmul(PQ,F_v)
			endif
			no_prec=no
		endif
		no1=iA(i)
		if (quad_proche == 1 .and. type_paire_T0(i) /= 0) then
			if (chi_d > 0.0_prec .and. proche_chi(no1,no,PQ)) then   ! CHI : element quasi singulier
				call DMKH_proche_T0(k,PQ,norc(no,:),Pt(no1,:),norc(no1,:),D1,K1,M1,H1)
			else
				S = Hmax/norm2(Pt(no1,:)-Pc(no,:))
				if (S < 0.15_prec .and. kH < 0.5_prec) then
					r = 1                                        ! 7 points
				elseif (kH >= 2.0_prec) then                     ! hors calibration : 33 (ancien) ou 48 tres pres
					r = 5 ; if (S >= 0.7_prec) r = 6
				elseif (S < 0.4_prec) then
					r = 2                                        ! 16
				elseif (S < 0.5_prec) then
					r = 3                                        ! 19
				elseif (S < 0.7_prec) then
					r = 4                                        ! 27
				else
					r = 6                                        ! 48
				endif
				compte(r) = compte(r) + 1
				call DMKH_T0(k,Qr(:,1:nr_pts(r),r),norc(no,:),Jac0(no),nr_pts(r),PGr(3,1:nr_pts(r),r), &
					Pt(no1,:),norc(no1,:),D1,K1,M1,H1)
			endif
			Anear(i)=D1+alpha*K1
			Bnear(i)=M1+alpha*H1
			cycle
		endif
		select case (type_paire_T0(i))
		case (0)
			if (chi_d > 0.0_prec) then   ! CHI : point interieur, pas de saut 1/2
				call DMKH_proche_T0(k,PQ,norc(no,:),Pt(no1,:),norc(no1,:),D1,K1,M1,H1)
				Anear(i)=D1+alpha*K1
				Bnear(i)=M1+alpha*H1
			else
			call DMKH_sing_T0(no1,PQ,k,PG_s,IS_s,K1,M1)
			Anear(i)=alpha*K1-cmplx(0.5_prec,0.0_prec,kind=prec)
			Bnear(i)=M1+alpha*cmplx(0.5_prec,0.0_prec,kind=prec)
			endif
		case (1)
			if (chi_d > 0.0_prec .and. proche_chi(no1,no,PQ)) then   ! CHI : element quasi singulier
				call DMKH_proche_T0(k,PQ,norc(no,:),Pt(no1,:),norc(no1,:),D1,K1,M1,H1)
			else
			call DMKH_T0(k,Qc,norc(no,:),Jac0(no),IS_c,PG_c(3,:),Pt(no1,:),norc(no1,:),D1,K1,M1,H1)
			endif
			Anear(i)=D1+alpha*K1
			Bnear(i)=M1+alpha*H1
		case default
			if (chi_d > 0.0_prec .and. proche_chi(no1,no,PQ)) then   ! CHI : element quasi singulier
				call DMKH_proche_T0(k,PQ,norc(no,:),Pt(no1,:),norc(no1,:),D1,K1,M1,H1)
			else
			call DMKH_T0(k,Qv,norc(no,:),Jac0(no),IS_v,PG_v(3,:),Pt(no1,:),norc(no1,:),D1,K1,M1,H1)
			endif
			Anear(i)=D1+alpha*K1
			Bnear(i)=M1+alpha*H1
		end select
	enddo
	if (quad_proche == 1 .and. verbose >= 2) then
		!$omp critical (crit_quad_proche)
		write(log_unit,'(A,F9.1,A,6(I0,A))') ' near-field quadrature, k = ', k, ' m-1: pairs with 7/16/19/27/33/48 points: ', &
			compte(1),' / ',compte(2),' / ',compte(3),' / ',compte(4),' / ',compte(5),' / ',compte(6),''
		!$omp end critical (crit_quad_proche)
	endif
	end subroutine ABnear_FMM_T0_adapt
!
!--------------------------------------------------------------------------------------------
!
! CHI (Cunefare 1989, essai 2026-10-14) : points cibles Pt = Pc - chi_d sqrt(Jac0) norc
! (dans le corps, la normale va vers le fluide). chi_d = 0 : Pt = Pc
	subroutine init_points_cibles()
	! Profondeur par element : chi_i = min(chi_d, 0,5 x epaisseur locale / h_i), h_i = sqrt(Jac0),
	! pour que le point reste a mi-epaisseur dans les parties minces (pavillon). Epaisseur =
	! distance le long de -n jusqu'a la paroi opposee, cherchee jusqu'a 2 chi_d h_i (au-dela,
	! chi_i = chi_d) par lancer de rayon (Moller-Trumbore) sur les triangles d'une grille
	! uniforme (chaque triangle range dans les cellules que couvre sa boite ; un seul tableau
	! contigu rempli en deux passes, comptage puis remplissage : Ericson 2005, 7.1.4).
	! Limite par la courbure (2026-10-15) : profondeur <= chi_courbure / kappa_i, kappa_i = plus grande
	! courbure convexe vue des voisins j (|Pc_j - Pc_i| < 1,5 h_i) : kappa_ij = (n_j - n_i).(Pc_j - Pc_i) / |Pc_j - Pc_i|^2.
	! Deux points voisins a la profondeur d sur un rayon R se rapprochent du facteur 1 - d/R.
	use parametres
	implicit none
	integer :: i, j, l, ix, iy, iz, ix0, ix1, iy0, iy1, iz0, iz1, nx, ny, nz, nred, ncourb
	integer, allocatable :: cnt(:), lst(:), pos(:)
	real(prec) :: xmin(3), xmax(3), cs, hmoy, hi, ls, t, tmin, b0(3), b1(3), chimin, dv(3), dd, kap, kmax
	real(prec), allocatable :: chi_el(:)
	logical, allocatable :: courbe(:)
	if (allocated(Pt)) deallocate(Pt)
	allocate(Pt(NE,3))
	if (chi_d <= 0.0_prec) then
		Pt = Pc
		return
	endif
	allocate(chi_el(NE), courbe(NE))
	! grille : cellule ~ longueur des rayons (2 chi_d h moyen), au plus 2e7 cellules
	xmin = minval(Ps(1:NP,:),1) - 1.0e-9_prec
	xmax = maxval(Ps(1:NP,:),1) + 1.0e-9_prec
	hmoy = sum(sqrt(Jac0(1:NE)))/NE
	cs = 2.0_prec*chi_d*hmoy
	do while (product(int((xmax-xmin)/cs, 8)+1_8) > 20000000_8)
		cs = 1.5_prec*cs
	enddo
	nx = int((xmax(1)-xmin(1))/cs)+1; ny = int((xmax(2)-xmin(2))/cs)+1; nz = int((xmax(3)-xmin(3))/cs)+1
	allocate(cnt(nx*ny*nz+1))
	cnt = 0
	do l=1,2                       ! 1 : comptage, 2 : remplissage
		if (l == 2) then
			do i=2,nx*ny*nz+1
				cnt(i) = cnt(i) + cnt(i-1)
			enddo
			allocate(lst(cnt(nx*ny*nz+1)), pos(nx*ny*nz))
			pos = cnt(1:nx*ny*nz)
		endif
		do j=1,NE
			b0 = minval(Ps(E(j,1:3),:),1); b1 = maxval(Ps(E(j,1:3),:),1)
			call cellules(b0, b1, ix0, ix1, iy0, iy1, iz0, iz1)
			do iz=iz0,iz1; do iy=iy0,iy1; do ix=ix0,ix1
				if (l == 1) then
					cnt(ix+(iy-1)*nx+(iz-1)*nx*ny+1) = cnt(ix+(iy-1)*nx+(iz-1)*nx*ny+1) + 1
				else
					pos(ix+(iy-1)*nx+(iz-1)*nx*ny) = pos(ix+(iy-1)*nx+(iz-1)*nx*ny) + 1
					lst(pos(ix+(iy-1)*nx+(iz-1)*nx*ny)) = j
				endif
			enddo; enddo; enddo
		enddo
	enddo
	!$omp parallel do private(i,j,l,ix,iy,iz,ix0,ix1,iy0,iy1,iz0,iz1,hi,ls,t,tmin,b0,b1,dv,dd,kap,kmax) schedule(dynamic,256)
	do i=1,NE
		hi = sqrt(Jac0(i))
		ls = 2.0_prec*chi_d*hi
		b0 = min(Pc(i,:), Pc(i,:)-ls*norc(i,:)); b1 = max(Pc(i,:), Pc(i,:)-ls*norc(i,:))
		call cellules(b0, b1, ix0, ix1, iy0, iy1, iz0, iz1)
		tmin = ls
		do iz=iz0,iz1; do iy=iy0,iy1; do ix=ix0,ix1
			do l=cnt(ix+(iy-1)*nx+(iz-1)*nx*ny)+1, cnt(ix+(iy-1)*nx+(iz-1)*nx*ny+1)
				j = lst(l)
				if (j == i) cycle
				t = rayon_triangle(Pc(i,:), -norc(i,:), Ps(E(j,1),:), Ps(E(j,2),:), Ps(E(j,3),:))
				if (t > 1.0e-9_prec*hi .and. t < tmin) tmin = t
			enddo
		enddo; enddo; enddo
		chi_el(i) = min(chi_d, 0.5_prec*tmin/hi)
		if (k_chi_max > 0.0_prec) then          ! limite par la longueur d'onde (2026-10-29)
			if (k_chi_max*hi > 2.0_prec*pi/3.0_prec) then
				chi_el(i) = min(chi_el(i), 0.5_prec)
			elseif (k_chi_max*hi > 2.0_prec*pi/6.0_prec) then
				chi_el(i) = min(chi_el(i), 0.75_prec)
			endif
		endif
		courbe(i) = .false.
		if (chi_courbure > 0.0_prec) then
			kmax = 0.0_prec
			b0 = Pc(i,:) - 1.5_prec*hi; b1 = Pc(i,:) + 1.5_prec*hi
			call cellules(b0, b1, ix0, ix1, iy0, iy1, iz0, iz1)
			do iz=iz0,iz1; do iy=iy0,iy1; do ix=ix0,ix1
				do l=cnt(ix+(iy-1)*nx+(iz-1)*nx*ny)+1, cnt(ix+(iy-1)*nx+(iz-1)*nx*ny+1)
					j = lst(l)
					if (j == i) cycle
					dv = Pc(j,:) - Pc(i,:)
					dd = dot_product(dv, dv)
					if (dd >= (1.5_prec*hi)**2 .or. dd <= 0.0_prec) cycle
					kap = dot_product(norc(j,:) - norc(i,:), dv)/dd
					if (kap > kmax) kmax = kap
				enddo
			enddo; enddo; enddo
			if (kmax*chi_el(i)*hi > chi_courbure) then
				chi_el(i) = chi_courbure/(kmax*hi)
				courbe(i) = .true.
			endif
		endif
		Pt(i,:) = Pc(i,:) - chi_el(i)*hi*norc(i,:)
	enddo
	!$omp end parallel do
	nred = count(chi_el < chi_d*(1.0_prec-1.0e-12_prec))
	chimin = minval(chi_el)
	ncourb = count(courbe)
	if (verbose >= 2) write(log_unit,'(A,F6.3,A,I8,A,I8,A,F6.3,A,F6.3)') ' CHI: depth', chi_d, ' x h, reduced on', nred, &
		' elements (curvature:', ncourb, '), min', chimin, ', chi_courbure', chi_courbure
	if (k_chi_max > 0.0_prec .and. verbose >= 2) write(log_unit,'(A,F8.1,A,I8,A,I8,A)') ' CHI: wavelength limit, k max =', &
		k_chi_max, ' m-1: chi <= 0.75 on', &
		count(k_chi_max*sqrt(Jac0(1:NE)) > 2.0_prec*pi/6.0_prec .and. k_chi_max*sqrt(Jac0(1:NE)) <= 2.0_prec*pi/3.0_prec), &
		' elements (h > lambda/6), chi <= 0.5 on', count(k_chi_max*sqrt(Jac0(1:NE)) > 2.0_prec*pi/3.0_prec), ' (h > lambda/3)'
	deallocate(cnt, lst, pos, chi_el, courbe)
	contains
	subroutine cellules(a0, a1, jx0, jx1, jy0, jy1, jz0, jz1)
	real(prec), intent(in) :: a0(3), a1(3)
	integer, intent(out) :: jx0, jx1, jy0, jy1, jz0, jz1
	jx0 = max(1, min(nx, int((a0(1)-xmin(1))/cs)+1)); jx1 = max(1, min(nx, int((a1(1)-xmin(1))/cs)+1))
	jy0 = max(1, min(ny, int((a0(2)-xmin(2))/cs)+1)); jy1 = max(1, min(ny, int((a1(2)-xmin(2))/cs)+1))
	jz0 = max(1, min(nz, int((a0(3)-xmin(3))/cs)+1)); jz1 = max(1, min(nz, int((a1(3)-xmin(3))/cs)+1))
	end subroutine cellules
	end subroutine init_points_cibles

	! Moller-Trumbore : parametre t du point O + t d sur le triangle (A,B,C), -1 si pas d'intersection
	pure real(prec) function rayon_triangle(O, d, A, B, C)
	real(prec), intent(in) :: O(3), d(3), A(3), B(3), C(3)
	real(prec) :: e1(3), e2(3), hv(3), sv(3), qv(3), det, u, v
	rayon_triangle = -1.0_prec
	e1 = B-A; e2 = C-A
	hv = [d(2)*e2(3)-d(3)*e2(2), d(3)*e2(1)-d(1)*e2(3), d(1)*e2(2)-d(2)*e2(1)]
	det = dot_product(e1, hv)
	if (abs(det) < 1.0e-30_prec) return
	sv = O-A
	u = dot_product(sv, hv)/det
	if (u < 0.0_prec .or. u > 1.0_prec) return
	qv = [sv(2)*e1(3)-sv(3)*e1(2), sv(3)*e1(1)-sv(1)*e1(3), sv(1)*e1(2)-sv(2)*e1(1)]
	v = dot_product(d, qv)/det
	if (v < 0.0_prec .or. u+v > 1.0_prec) return
	rayon_triangle = dot_product(e2, qv)/det
	end function rayon_triangle

	! point Q du triangle (A,B,C) le plus proche de X (Ericson, Real-Time Collision Detection, 5.1.5)
	pure subroutine point_plus_proche(X, A, B, C, Q)
	real(prec), intent(in) :: X(3), A(3), B(3), C(3)
	real(prec), intent(out) :: Q(3)
	real(prec) :: ab(3), ac(3), ap(3), bp(3), cp(3), d1, d2, d3, d4, d5, d6, va, vb, vc, v, w, den
	ab = B-A; ac = C-A; ap = X-A
	d1 = dot_product(ab,ap); d2 = dot_product(ac,ap)
	if (d1 <= 0 .and. d2 <= 0) then; Q = A; goto 10; endif
	bp = X-B; d3 = dot_product(ab,bp); d4 = dot_product(ac,bp)
	if (d3 >= 0 .and. d4 <= d3) then; Q = B; goto 10; endif
	vc = d1*d4 - d3*d2
	if (vc <= 0 .and. d1 >= 0 .and. d3 <= 0) then; Q = A + d1/(d1-d3)*ab; goto 10; endif
	cp = X-C; d5 = dot_product(ab,cp); d6 = dot_product(ac,cp)
	if (d6 >= 0 .and. d5 <= d6) then; Q = C; goto 10; endif
	vb = d5*d2 - d1*d6
	if (vb <= 0 .and. d2 >= 0 .and. d6 <= 0) then; Q = A + d2/(d2-d6)*ac; goto 10; endif
	va = d3*d6 - d5*d4
	if (va <= 0 .and. (d4-d3) >= 0 .and. (d5-d6) >= 0) then
		Q = B + (d4-d3)/((d4-d3)+(d5-d6))*(C-B); goto 10
	endif
	den = 1.0_prec/(va+vb+vc); v = vb*den; w = vc*den
	Q = A + ab*v + ac*w
10	continue
	end subroutine point_plus_proche

	pure real(prec) function distance_point_triangle(X, A, B, C)
	real(prec), intent(in) :: X(3), A(3), B(3), C(3)
	real(prec) :: Q(3)
	call point_plus_proche(X, A, B, C, Q)
	distance_point_triangle = sqrt(sum((X-Q)**2))
	end function distance_point_triangle

	! CHI : le point cible de l'element no1 est-il a moins de 0,5 h de l'element source no ?
	! (paroi voisine ou opposee : la regle a 33 points y est fausse, voir quad_proche.py)
	logical function proche_chi(no1, no, PQ)
	integer, intent(in) :: no1, no
	real(prec), intent(in) :: PQ(3,3)
	real(prec) :: dpl
	dpl = abs(dot_product(Pt(no1,:)-PQ(:,1), norc(no,:)))
	proche_chi = .false.
	if (dpl < 1.0e-12_prec*sqrt(Jac0(no))) return      ! point dans le plan : Ma inapplicable
	proche_chi = distance_point_triangle(Pt(no1,:), PQ(:,1), PQ(:,2), PQ(:,3)) < 0.5_prec*sqrt(Jac0(no))
	end function proche_chi

	! point ip de l'octree evalue comme cible ? T0 : barycentres ; CHI : points interieurs
	logical function cible_active(ip)
	integer, intent(in) :: ip
	cible_active = (cibles_sep .eqv. ip > NE)
	end function cible_active

	! cible ip de l'octree : ligne no, point P, poids de l'equation simple (cG, noyau G) et de
	! l'equation derivee (cH, multiplie par alpha) : 1 et 1 (BM et CHI ; 0/1 servaient au CHI
	! hybride, essai retire le 2026-11-20)
	subroutine cible_T0(ip, no, P, cG, cH)
	integer, intent(in) :: ip
	integer, intent(out) :: no
	real(prec), intent(out) :: P(3), cG, cH
	cG = 1.0_prec ; cH = 1.0_prec
	no = ip
	if (cibles_sep .and. ip > NE) no = ip-NE
	P = Pt(no,:)
	end subroutine cible_T0

	! filet de securite (programmes de test) : Pt = Pc si personne ne l'a initialise
	subroutine verifier_points_cibles()
	if (allocated(Pt)) return
	!$omp critical (crit_points_cibles)
	if (.not. allocated(Pt)) call init_points_cibles()
	!$omp end critical (crit_points_cibles)
	end subroutine verifier_points_cibles

	! CHI : D, K, M, H de l'element (sommets PQ) au point P proche (element propre, ou paroi
	! voisine / opposee a moins de 0,5 h). Polaires autour d'un point c de l'element, triangle
	! coupe en 3 depuis c (sous-triangles degeneres sautes), Duffy en w, et sur chaque rayon
	! c + rho e : |x - P|^2 = a^2 + (rho + b)^2.
	! c = point de l'element le plus proche de P (Ericson 5.1.5 : le pied si P se projette dans
	! l'element), b = e.(c-P) (0 si c est le pied) ; transformation sinh (Johnston & Elliott)
	! dans la variable decalee s = rho + b : s = a sinh(nu), ds = g dnu. Polaires depuis c et
	! sous-triangles : schema PART de Hayami (Publ. RIMS 41, 2005), qui compare log L2, log L1
	! et (rho+d)^(-1/5). Jusqu'au 2026-10-30, pied DANS l'element : Ma & Kamiya (EABE 26, 2002)
	! eta = log sqrt(d^2 + rho^2) (= log L2 de Hayami). Comparaison a 16 x 16 points, 24 kHz
	! (quadrature_20261030/proche2.py) : Ma est meilleur seulement pour d < 0,05 h (4e-8 contre
	! 4e-7 a 5e-6 en sinh), cas exclu par les plafonds de profondeur CHI (d/h >= 0,2 environ) ;
	! il echoue quand le pied est plus pres d'une arete que d (critere de Hayami, 77 % des
	! couples reels de m3) : max 5,9e-5 contre 3e-7 en sinh ; log L1 et (rho+d)^(-1/5) ne font
	! pas mieux que sinh. Hors de l'element, log converge mal (b -> 0 le long de l'arete :
	! racine a la borne), d'ou sinh (chiefb_20261013/quad_signe.py).
	! Essai 2026-11-02 (profondeur_20261102) : Ma retabli selon le critere de Hayami (pied dans
	! l'element, plus loin de chaque arete que d), chi = 0,05 et 0,25 sur 3 spheres : p99 identique
	! a sinh a 4 chiffres ; retire.
	subroutine DMKH_proche_T0(k,PQ,norQ,P,norP,D1,K1,M1,H1)
	use parametres
	implicit none
	real(prec), intent(in) :: k, PQ(3,3), norQ(3), P(3), norP(3)
	complex(prec), intent(out) :: D1, K1, M1, H1
	! Nombre de points (quad_chi = 1, 2026-11-03) selon d/h, d distance de P a l'element, h = sqrt(2 A) :
	! Gauss 48 sans decoupage si d/h >= 0,9, sinh 8 x 8 / 10 x 10 / 12 x 12 aux seuils 0,75 / 0,5 / 0,25,
	! 16 x 16 en dessous ; erreur <= 6e-7 (couples reels de m3, 24 kHz, quadrature_20261030/proche3.py
	! et proche4.py). quad_chi = 0 : 16 x 16 partout.
	integer, parameter :: NMAX = 16
	real(prec) :: xu(NMAX), wu(NMAX), cpt(3), ea(3), eb(3), ec(3), er(3), eu(3), dh
	real(prec) :: Q(3,NMAX*NMAX), wq(NMAX*NMAX), lr, e0, e1, eta, rr, jac, w, surf, surf0, bb, a2, aa
	real(prec) :: PG48(3,48), F48(3,48), dF1(3,48), dF2(3,48), Q48(3,48)
	complex(prec) :: d, kk, m, h
	integer :: is, iu, iw, n, nn, npt
	ec = cross3(PQ(:,2)-PQ(:,1), PQ(:,3)-PQ(:,1))
	surf0 = dot_product(ec, norQ)
	call point_plus_proche(P, PQ(:,1), PQ(:,2), PQ(:,3), cpt)   ! le pied si P se projette sur l'element
	D1 = 0; K1 = 0; M1 = 0; H1 = 0
	nn = NMAX
	if (quad_chi == 1) then
		dh = norm2(P-cpt)/sqrt(surf0)
		if (dh >= 0.9_prec) then              ! assez loin : regle de Gauss a 48 points, sans decoupage
			call pgt6(PG48, 48)
			call FORM_T3(PG48(1:2,:), 48, F48, dF1, dF2)
			Q48 = matmul(PQ, F48)
			call DMKH_T0(k,Q48,norQ,surf0,48,PG48(3,:),P,norP,D1,K1,M1,H1)
			return
		elseif (dh >= 0.75_prec) then
			nn = 8
		elseif (dh >= 0.5_prec) then
			nn = 10
		elseif (dh >= 0.25_prec) then
			nn = 12
		endif
	endif
	npt = nn*nn
	call gauleg(nn, xu(1:nn), wu(1:nn))
	do is=1,3
		ea = PQ(:,is) - cpt
		eb = PQ(:,mod(is,3)+1) - PQ(:,is)
		surf = dot_product(cross3(ea, eb), norQ)
		if (surf <= 1.0e-12_prec*surf0) cycle         ! c sur ce cote : sous-triangle degenere
		n = 0
		do iw=1,nn
			w = (xu(iw)+1.0_prec)/2.0_prec
			er = ea + w*eb                    ! rayon de c vers l'arete
			lr = sqrt(sum(er**2))
			eu = er/lr
			bb = dot_product(eu, cpt-P)
			a2 = max(sum((cpt-P)**2) - bb**2, 1.0e-30_prec)
			aa = sqrt(a2)
			e0 = asinh(bb/aa)
			e1 = asinh((lr+bb)/aa)
			do iu=1,nn
				n = n+1
				eta = e0 + (e1-e0)*(xu(iu)+1.0_prec)/2.0_prec
				rr = aa*sinh(eta) - bb
				jac = rr*aa*cosh(eta)*(e1-e0)/2.0_prec            ! rr drr = rr a cosh(nu) dnu
				Q(:,n) = cpt + rr*eu
				wq(n) = wu(iu)*jac*surf/lr**2 * wu(iw)/2.0_prec    ! dA = surf/L^2 rr drr dw
			enddo
		enddo
		call DMKH_T0(k,Q(:,1:npt),norQ,1.0_prec,npt,wq(1:npt),P,norP,d,kk,m,h)
		D1 = D1+d; K1 = K1+kk; M1 = M1+m; H1 = H1+h
	enddo
	contains
	pure function cross3(u, v)
	real(prec), intent(in) :: u(3), v(3)
	real(prec) :: cross3(3)
	cross3 = [u(2)*v(3)-u(3)*v(2), u(3)*v(1)-u(1)*v(3), u(1)*v(2)-u(2)*v(1)]
	end function cross3
	end subroutine DMKH_proche_T0
!
!--------------------------------------------------------------------------------------------
!
! Calcul Anear pour FMM, BEM méthode Burton/Miller
! pour diffraction sur structure rigide Bnear =0
!
!--------------------------------------------------------------------------------------------
!
	subroutine Anear_FMM_T0(Anear,octree,NL,k,alpha)
!	calcul la matrice creuse Anear et Bnear pour FMM
	use module_octree
	use parametres
	implicit none

! Variables d'entrée
	integer, intent(in) :: NL
	real(prec), intent(in) :: k
	complex(prec), intent(in) :: alpha
	type(tree), intent(in) :: octree(NL)
! Variables de sortie
	complex(prec), intent(inout) :: Anear(num_K0)
! Variables locales
	integer :: no,num_K,no1
	integer, dimension(1:3) :: nPe
	real(prec), dimension(1:3,1:3) :: PQ
	real(prec), dimension(1:3,1:IS_c) :: PG_c
	real(prec), dimension(1:3,1:IS_v) :: PG_v
	real(prec), dimension(1:2,1:IS_s) :: PG_s
	real(prec), dimension(1:IS_s) :: x, wth
	real(prec), dimension(1:3,1:IS_c) :: F_c, dF1_c, dF2_c
    real(prec), dimension(1:3,1:IS_v) :: F_v, dF1_v, dF2_v
	real(prec) :: Qc(1:3,IS_c), Qv(1:3,IS_v)
	complex(prec) :: D1, K1
	integer :: i1, i2, i3, i4
	call pgt6(PG_c,IS_c)  ! points dans la même cellule
	call pgt6(PG_v,IS_v)  ! points dans la cellule voisine
	call gauleg(IS_s, x, wth)
	PG_s(1,:)=x		! points singuliers
	PG_s(2,:)=wth
	call FORM_T3(PG_c(1:2,:),IS_c,F_c,dF1_c,dF2_c)
	call FORM_T3(PG_v(1:2,:),IS_v,F_v,dF1_v,dF2_v)

! init
	num_K=0
	do i1=1,octree(NL)%NL   ! boucle sur cellule dernier niveau
		do i2=1,octree(NL)%cell(i1)%Nb_c
		  no=octree(NL)%cell(i1)%child(i2)
		  nPe=E(no,:) ! Points de l'élément
		  PQ=transpose(Ps(nPe,:)) ! Coordonnées des points de l'éléments
		  Qc=matmul(PQ,F_c)
		  Qv=matmul(PQ,F_v)
! points dans la cellule
		  do i3=1,octree(NL)%cell(i1)%Nb_c
			no1=octree(NL)%cell(i1)%child(i3)
			if (i2.eq.i3) then    ! point singulier
				call DK_sing_T0(no1,PQ,k,PG_s,IS_s,K1)
					num_K=num_K+1
					iA(num_K)=octree(NL)%cell(i1)%child(i3)
					jA(num_K)=octree(NL)%cell(i1)%child(i2)
					Anear(num_K)=alpha*K1-cmplx(0.5_prec,0.0_prec,kind=prec) ! prise en compte de l'angle solide+ BM
			else			! point non singulier
				call DK_T0(k,Qc,norc(no,:),Jac0(no),IS_c,PG_c(3,:),Pc(no1,:),norc(no1,:),D1,K1)
				num_K=num_K+1
				iA(num_K)=octree(NL)%cell(i1)%child(i3)
				jA(num_K)=octree(NL)%cell(i1)%child(i2)
				Anear(num_K)=D1+alpha*K1
			endif
		  end do
! points dans les cellules voisines
  		  do i3=1,octree(NL)%cell(i1)%Nb_v
			do i4=1,octree(NL)%cell(octree(NL)%cell(i1)%voisin(i3))%Nb_c
				no1=octree(NL)%cell(octree(NL)%cell(i1)%voisin(i3))%child(i4)
				call DK_T0(k,Qv,norc(no,:),Jac0(no),IS_v,PG_v(3,:),Pc(no1,:),norc(no1,:),D1,K1)
				num_K=num_K+1
				iA(num_K)=octree(NL)%cell(octree(NL)%cell(i1)%voisin(i3))%child(i4)
				jA(num_K)=octree(NL)%cell(i1)%child(i2)
				Anear(num_K)=D1+alpha*K1
			end do
		  end do
		end do
	end do
	end subroutine Anear_FMM_T0
!------------------------------------------------------------------
!
! Intégration singulière DMKH_sing_T0
!
!------------------------------------------------------------------
	subroutine DMKH_sing_T0(no,PQ,k,PG,IS,K2,M2)
! calcul le coef K pour BEM dans le cas singulier
	use parametres
	implicit none
! variables d'entrée
	integer, intent(in) :: IS,no
	real(prec), intent(in) :: PQ(1:3,1:3), k
	real(prec), intent(in) :: PG(1:2,1:IS)
! variables de sortie
	complex(prec), intent(inout) :: K2,M2
! variables locales
	integer,dimension (1:3), parameter :: permut=(/ 2, 3, 1 /)
	real(prec) :: prodsca,theta
	real(prec) :: Alpha(1:3,1:size(PG,2)), r(1:3,1:size(PG,2)), Jac1(1:3)
	real(prec) :: R0(1:3,1:3),Rd(1:3,1:3)
	real(prec) :: nr0(1:3),nrd(1:3), h(1:3), alpha1(1:3), alpham(1:3), alpha2(1:3),test
	complex(prec) :: K0(1:3,1:size(PG,2)), K1
	complex(prec) :: M0(1:3,1:size(PG,2)), M1
	integer :: nPG
	character(len=1) :: answer


!	one(1,:)=1.0_prec
!	one1(1,:)=1.0_prec

	R0=spread(Pc(no,:),2,3)-PQ
	Rd=PQ(:,permut)-PQ
	nr0=sqrt(sum(R0**2,1))
	nrd=sqrt(sum(Rd**2,1))

! Triangle 1
	prodsca=dot_product(R0(:,1),Rd(:,1))
	theta=acos(prodsca/nr0(1)/nrd(1))
	h(1)=sin(theta)*nr0(1)
	alpha1(1)=0
	alpham(1)=pi/2-theta
	test=h(1)/nr0(2)
	if (test.gt.1) test=1.0_prec
	alpha2(1)=alpham(1)+acos(test)
	if (alpha2(1) /= alpha2(1)) then
	  write(*,*) 'alpha2(1) isnan = ',test,acos(test)
	  read(*,*) answer
	end if


! Triangle 2
	prodsca=dot_product(R0(:,2),Rd(:,2))
	theta=acos(prodsca/nr0(2)/nrd(2))
	h(2)=sin(theta)*nr0(2)
	alpha1(2)=0
	alpham(2)=pi/2-theta
	test=h(2)/nr0(3)
	if (test.gt.1) test=1.0_prec
	alpha2(2)=alpham(2)+acos(test)

! Triangle 3
	prodsca=dot_product(R0(:,3),Rd(:,3))
	theta=acos(prodsca/nr0(3)/nrd(3))
	h(3)=sin(theta)*nr0(3)
	alpha1(3)=0
	alpham(3)=pi/2-theta
	test=h(3)/nr0(1)
	if (test.gt.1) test=1.0_prec
	alpha2(3)=alpham(3)+acos(test)

	nPG=size(PG,2)
	Alpha=(spread(alpha2-alpha1,dim=2,ncopies=nPG)*spread(PG(1,:),1,3)+ spread(alpha1+alpha2,2,nPG))/2
	r=spread(h,2,nPG) / cos(Alpha-spread(alpham,2,nPG))
	Jac1=(alpha2-alpha1)/2

	K0=exp(ji*k*r)/(4*pi*r)*(spread(Jac1,2,nPG))
	M0=exp(ji*k*r)*(spread(Jac1,2,nPG))
!	write(*,*) 'M0=',M0
!	write(*,*) 'K0=',K0

	K1=sum(matmul(K0,PG(2,:)))
	M1=sum(matmul(M0,PG(2,:)))
!	write(*,*) 'M1=',M1,' K1=',K1

	K2=ji*k/2.0_prec - K1
	M2=ji/(2.0_prec*k)*(1-M1/(2.0_prec*pi))
!	write(*,*) 'M2=',M2,' K2=',K2
!	write(*,*) 'k=',k,' ji=',ji

	end subroutine DMKH_sing_T0
!------------------------------------------------------------------
!
! Intégration singulière DMKH_sing_T0
!
!------------------------------------------------------------------
	subroutine DK_sing_T0(no,PQ,k,PG,IS,K2)
! calcul le coef K pour BEM dans le cas singulier
	use parametres
	implicit none
! variables d'entrée
	integer, intent(in) :: IS,no
	real(prec), intent(in) :: PQ(1:3,1:3), k
	real(prec), intent(in) :: PG(1:2,1:IS)
! variables de sortie
	complex(prec), intent(inout) :: K2
! variables locales
	integer,dimension (1:3), parameter :: permut=(/ 2, 3, 1 /)
	real(prec) :: prodsca,theta
	real(prec) :: Alpha(1:3,1:size(PG,2)), r(1:3,1:size(PG,2)), Jac1(1:3)
	real(prec) :: R0(1:3,1:3),Rd(1:3,1:3)
	real(prec) :: nr0(1:3),nrd(1:3), h(1:3), alpha1(1:3), alpham(1:3), alpha2(1:3),test
	complex(prec) :: K0(1:3,1:size(PG,2)), K1
	integer :: nPG
	character(len=1) :: answer


!	one(1,:)=1.0_prec
!	one1(1,:)=1.0_prec

	R0=spread(Pc(no,:),2,3)-PQ
	Rd=PQ(:,permut)-PQ
	nr0=sqrt(sum(R0**2,1))
	nrd=sqrt(sum(Rd**2,1))

! Triangle 1
	prodsca=dot_product(R0(:,1),Rd(:,1))
	theta=acos(prodsca/nr0(1)/nrd(1))
	h(1)=sin(theta)*nr0(1)
	alpha1(1)=0
	alpham(1)=pi/2-theta
	test=h(1)/nr0(2)
	if (test.gt.1) test=1.0_prec
	alpha2(1)=alpham(1)+acos(test)
	if (alpha2(1) /= alpha2(1)) then
	  write(*,*) 'alpha2(1) isnan = ',test,acos(test)
	  read(*,*) answer
	end if


! Triangle 2
	prodsca=dot_product(R0(:,2),Rd(:,2))
	theta=acos(prodsca/nr0(2)/nrd(2))
	h(2)=sin(theta)*nr0(2)
	alpha1(2)=0
	alpham(2)=pi/2-theta
	test=h(2)/nr0(3)
	if (test.gt.1) test=1.0_prec
	alpha2(2)=alpham(2)+acos(test)

! Triangle 3
	prodsca=dot_product(R0(:,3),Rd(:,3))
	theta=acos(prodsca/nr0(3)/nrd(3))
	h(3)=sin(theta)*nr0(3)
	alpha1(3)=0
	alpham(3)=pi/2-theta
	test=h(3)/nr0(1)
	if (test.gt.1) test=1.0_prec
	alpha2(3)=alpham(3)+acos(test)

	nPG=size(PG,2)
	Alpha=(spread(alpha2-alpha1,dim=2,ncopies=nPG)*spread(PG(1,:),1,3)+ spread(alpha1+alpha2,2,nPG))/2
	r=spread(h,2,nPG) / cos(Alpha-spread(alpham,2,nPG))
	Jac1=(alpha2-alpha1)/2

	K0=exp(ji*k*r)/(4*pi*r)*(spread(Jac1,2,nPG))
	K1=sum(matmul(K0,PG(2,:)))
	K2=ji*k/2.0_prec - K1

	end subroutine DK_sing_T0

!------------------------------------------------------------------
!
! Intégration DMKH_T0
!
!------------------------------------------------------------------
	subroutine DMKH_T0(k,Q,norQ,jac1,IS,wq,P,norP,D1,K1,M1,H1)
! calcul les coefs D, M, H et K pour BEM
	use parametres
	implicit none
! variables d'entrée
	integer, intent(in) :: IS
	real(prec), intent(in) :: Q(1:3,1:IS), k,jac1
	real(prec), intent(in) :: norQ(1:3),norP(3),wq(IS),P(3)
! variables de sortie
	complex(prec), intent(out) :: D1, K1, M1, H1
! variables locales
	real(prec) :: R(1:3,IS)
	real(prec) :: nR(IS), NQR(IS), NPR(IS), Nn
	complex(prec) :: D0(IS), K0(IS), M0(IS), H0(IS)
	complex(prec) :: E(IS)   ! exp(ikr), calcule une seule fois (2026-10-30 : 4 fois avant, ~50 % du temps HRTF ici)
	! DEBUG : character(len=1) :: answer


! calcul des vecteurs rayons jusqu'aux points de Gauss
!
! Calcul de R et norme de R
!
	R=Q-spread(P,2,IS)   ! vecteur R
	nR=sqrt(sum(R**2,1))               ! norme de R

!produits scalaires (rayon.Normale*Jacobien)
	NQR=matmul(transpose(R),norQ)
	NPR=matmul(transpose(R),norP)
	Nn=dot_product(norP,norQ)
	!Jac1=sqrt(sum(norQ**2))

! fonction à intégrer
	E=exp(ji*k*nR)
	M0=E /(4*pi*nR) * Jac1
	D0=-(1-ji*k*nR)*E*NQR /(4*pi*(nR**3)) *Jac1  !% dNQ green
	H0=(1-ji*k*nR)*E*NPR /(4*pi*(nR**3)) * Jac1   !% dNP green
	K0= (1-ji*k*nR)*E * Nn / (4*pi*(nR**3)) * Jac1 - ((3 - 3 *ji*k*nR + (ji*k*nR)**2)*E) * NQR * NPR /&
		(4*pi*(nR**5)) * Jac1     !% derivée seconde de G


! intégration

 	D1=sum(wq*D0)
 	K1=sum(wq*K0)
 	M1=sum(wq*M0)
 	H1=sum(wq*H0)

	end subroutine DMKH_T0
!------------------------------------------------------------------
!
! Intégration DM_T0
!
!------------------------------------------------------------------
	subroutine DM_T0(k,Q,norQ,jac1,IS,wq,P,D1,M1)
! calcul les coefs D, M, H et K pour BEM
	use parametres
	implicit none
! variables d'entrée
	integer, intent(in) :: IS
	real(prec), intent(in) :: Q(1:3,1:IS), k,jac1
	real(prec), intent(in) :: norQ(1:3),wq(IS),P(3)
! variables de sortie
	complex(prec), intent(out) :: D1,M1
! variables locales
	real(prec) :: R(1:3,IS)
	real(prec) :: nR(IS), NQR(IS)
	complex(prec) :: D0(IS), M0(IS), E(IS)
	! DEBUG : character(len=1) :: answer


! calcul des vecteurs rayons jusqu'aux points de Gauss
!
! Calcul de R et norme de R
!
	R=Q-spread(P,2,IS)   ! vecteur R
	nR=sqrt(sum(R**2,1))               ! norme de R

!produits scalaires (rayon.Normale*Jacobien)
	NQR=matmul(transpose(R),norQ)
!	NPR=matmul(transpose(R),norP)
!	Nn=dot_product(norP,norQ)
	!Jac1=sqrt(sum(norQ**2))

! fonction à intégrer
	E=exp(ji*k*nR)   ! une seule fois (2026-11-03)
	M0=E /(4*pi*nR) * Jac1
	D0=-(1-ji*k*nR)*E*NQR /(4*pi*(nR**3)) *Jac1  !% dNQ green
!	H0=(1-ji*k*nR)*exp(ji*k*nR)*NPR /(4*pi*(nR**3)) * Jac1   !% dNP green
!	K0= (1-ji*k*nR)*exp(ji*k*nR) * Nn / (4*pi*(nR**3)) * Jac1 - ((3 - 3 *ji*k*nR + (ji*k*nR)**2)*exp(ji*k*nR)) * NQR * NPR /&
!		(4*pi*(nR**5)) * Jac1     !% derivée seconde de G


! intégration

 	D1=sum(wq*D0)
 !	K1=sum(wq*K0)
 	M1=sum(wq*M0)
 !	H1=sum(wq*H0)

	end subroutine DM_T0
!------------------------------------------------------------------
!
! Intégration DK_T0
!
!------------------------------------------------------------------
	subroutine DK_T0(k,Q,norQ,jac1,IS,wq,P,norP,D1,K1)
! calcul les coefs D, M, H et K pour BEM
	use parametres
	implicit none
! variables d'entrée
	integer, intent(in) :: IS
	real(prec), intent(in) :: Q(1:3,1:IS), k,jac1
	real(prec), intent(in) :: norQ(1:3),norP(3),wq(IS),P(3)
! variables de sortie
	complex(prec), intent(out) :: D1, K1
! variables locales
	real(prec) :: R(1:3,IS)
	real(prec) :: nR(IS), NQR(IS), NPR(IS), Nn
	complex(prec) :: D0(IS), K0(IS)
	! DEBUG : character(len=1) :: answer


! calcul des vecteurs rayons jusqu'aux points de Gauss
!
! Calcul de R et norme de R
!
	R=Q-spread(P,2,IS)   ! vecteur R
	nR=sqrt(sum(R**2,1))               ! norme de R

!produits scalaires (rayon.Normale*Jacobien)
	NQR=matmul(transpose(R),norQ)
	NPR=matmul(transpose(R),norP)
	Nn=dot_product(norP,norQ)
	!Jac1=sqrt(sum(norQ**2))

! fonction à intégrer
	D0=-(1-ji*k*nR)*exp(ji*k*nR)*NQR /(4*pi*(nR**3)) *Jac1  !% dNQ green
	K0= (1-ji*k*nR)*exp(ji*k*nR) * Nn / (4*pi*(nR**3)) * Jac1 - ((3 - 3 *ji*k*nR + (ji*k*nR)**2)*exp(ji*k*nR)) * NQR * NPR /&
		(4*pi*(nR**5)) * Jac1     !% derivée seconde de G


! intégration

 	D1=sum(wq*D0)
 	K1=sum(wq*K0)

	end subroutine DK_T0
!-------------------------------------------------------------------------------------
!
! ABnear_FMM_v3 : calcul D et M : pas de Burton/Miller
!
!-------------------------------------------------------------------------------------
    subroutine ABnear0_FMM_T0(Anear,Bnear,octree,NL,k)
!	calcul la matrice creuse Anear et Bnear pour FMM
	use parametres
	implicit none

! Variables d'entrée
	integer, intent(in) :: NL
!	integer, intent(in) :: E(1:NE,1:3)
!	real(prec), intent(in) :: Pc(1:NE,1:3), Ps(1:NP,1:3), nor(1:NE,1:3)
!	real(prec), dimension(1:NP,1:3), intent(in) :: Ps
	real(prec), intent(in) :: k
!	complex(prec), intent(in) :: alpha
	type(tree), intent(in) :: octree(NL)
! Variables de sortie

!	integer, intent(inout) :: iA(num_K0), jA(num_K0)
	complex(prec), intent(inout) :: Anear(num_K0),Bnear(num_K0)

! Variables locales
!	integer, parameter :: IS_c=28, IS_v=28, IS_s=12, coef=200
	integer :: no,num_K,no1
	integer, dimension(1:3) :: nPe
	real(prec), dimension(1:3,1:3) :: PQ
	real(prec), dimension(1:3,1:IS_c) :: PG_c
	real(prec), dimension(1:3,1:IS_v) :: PG_v
	real(prec) :: PG_c2(2,IS_c),PG_v2(2,IS_v),PGc3(IS_c),PGv3(IS_v),Pcin(3),norcin(3)
	real(prec), dimension(1:2,1:IS_s) :: PG_s
	real(prec), dimension(1:IS_s) :: x, wth
	real(prec), dimension(1:3,1:IS_c) :: F_c, dF1_c, dF2_c
	real(prec), dimension(1:3,1:IS_v) :: F_v, dF1_v, dF2_v
	real(prec) :: Qc(1:3,IS_c), Qv(1:3,IS_v)

!	integer, dimension(1:coef*NE),allocatable :: iK0, jK0, iD0, jD0
!	complex(prec), dimension(1:coef*NE),allocatable :: Dnear0, Knear0
	complex(prec) :: D1, M1
	integer :: i1, i2, i3, i4

! points pour intégration
	call pgt6(PG_c,IS_c)  ! points dans la même cellule
	call pgt6(PG_v,IS_v)  ! points dans la cellule voisine
	call gauleg(IS_s, x, wth)
	PG_s(1,:)=x		! points singuliers
	PG_s(2,:)=wth
	PG_c2=PG_c(1:2,:)
	PG_v2=PG_v(1:2,:)
	PGc3=PG_c(3,:)
	PGv3=PG_v(3,:)
	call FORM_T3(PG_c2,IS_c,F_c,dF1_c,dF2_c)
	call FORM_T3(PG_v2,IS_v,F_v,dF1_v,dF2_v)

! init
	!num_D=0
	num_K=0
	do i1=1,octree(NL)%NL   ! boucle sur cellule dernier niveau

		do i2=1,octree(NL)%cell(i1)%Nb_c

			no=octree(NL)%cell(i1)%child(i2)
			nPe=E(no,:) ! Points de l'élément
			PQ=transpose(Ps(nPe,:)) ! Coordonnées des points de l'éléments
			Qc=matmul(PQ,F_c)
			Qv=matmul(PQ,F_v)

			! points dans la cellule
			do i3=1,octree(NL)%cell(i1)%Nb_c
				no1=octree(NL)%cell(i1)%child(i3)
				if (i2.eq.i3) then    ! point singulier
					call DM_sing_T0(no1,PQ,k,PG_s,IS_s,M1)
					num_K=num_K+1
					Anear(num_K)=-cmplx(0.5_prec,0.0_prec,kind=prec) ! prise en compte de l'angle solide
					Bnear(num_K)=M1
				else			! point non singulier
					Pcin=Pc(no1,:)
					norcin=norc(no1,:)
					call DM_T0(k,Qc,norc(no,:),Jac0(no),IS_c,PGc3,Pcin,D1,M1)
					num_K=num_K+1
					Anear(num_K)=D1
					Bnear(num_K)=M1
				endif
			end do
			!points dans les cellules voisines
			do i3=1,octree(NL)%cell(i1)%Nb_v
				do i4=1,octree(NL)%cell(octree(NL)%cell(i1)%voisin(i3))%Nb_c
					no1=octree(NL)%cell(octree(NL)%cell(i1)%voisin(i3))%child(i4)
					Pcin=Pc(no1,:)
					norcin=norc(no1,:)
					call DM_T0(k,Qv,norc(no,:),Jac0(no),IS_v,PGv3,Pcin,D1,M1)
					num_K=num_K+1
					Anear(num_K)=D1
					Bnear(num_K)=M1
				end do
			end do
		enddo
	end do

	end subroutine ABnear0_FMM_T0
!------------------------------------------------------------------
!
! Intégration singulière DM_sing_T0
!
!------------------------------------------------------------------
	subroutine DM_sing_T0(no,PQ,k,PG,IS,M2)
! calcul le coef K pour BEM dans le cas singulier
	use parametres
	implicit none
! variables d'entrée
	integer, intent(in) :: IS,no
	real(prec), intent(in) :: PQ(1:3,1:3), k
	real(prec), intent(in) :: PG(1:2,1:IS)
! variables de sortie
	complex(prec), intent(inout) :: M2
!	complex(prec), intent(inout) :: K2,M2
! variables locales
	integer,dimension (1:3), parameter :: permut=(/ 2, 3, 1 /)
	real(prec) :: prodsca,theta
	real(prec) :: Alpha(1:3,1:size(PG,2)), r(1:3,1:size(PG,2)), Jac(1:3)
	real(prec) :: R0(1:3,1:3),Rd(1:3,1:3)
	real(prec) :: nr0(1:3),nrd(1:3), h(1:3), alpha1(1:3), alpham(1:3), alpha2(1:3),test
!	complex(prec) :: K0(1:3,1:size(PG,2)), K1
	complex(prec) :: M0(1:3,1:size(PG,2)), M1
	integer :: nPG
	character(len=1) :: answer


	R0=spread(Pc(no,:),2,3)-PQ
	Rd=PQ(:,permut)-PQ
	nr0=sqrt(sum(R0**2,1))
	nrd=sqrt(sum(Rd**2,1))

! Triangle 1
	prodsca=dot_product(R0(:,1),Rd(:,1))
	theta=acos(prodsca/nr0(1)/nrd(1))
	h(1)=sin(theta)*nr0(1)
	alpha1(1)=0
	alpham(1)=pi/2-theta
	test=h(1)/nr0(2)
	if (test.gt.1) test=1.0_prec
	alpha2(1)=alpham(1)+acos(test)
	if (isnan(alpha2(1))) then
	  write(*,*) 'alpha2(1) isnan = ',test,acos(test)
	  read(*,*) answer
	end if


! Triangle 2
	prodsca=dot_product(R0(:,2),Rd(:,2))
	theta=acos(prodsca/nr0(2)/nrd(2))
	h(2)=sin(theta)*nr0(2)
	alpha1(2)=0
	alpham(2)=pi/2-theta
	test=h(2)/nr0(3)
	if (test.gt.1) test=1.0_prec
	alpha2(2)=alpham(2)+acos(test)

! Triangle 3
	prodsca=dot_product(R0(:,3),Rd(:,3))
	theta=acos(prodsca/nr0(3)/nrd(3))
	h(3)=sin(theta)*nr0(3)
	alpha1(3)=0
	alpham(3)=pi/2-theta
	test=h(3)/nr0(1)
	if (test.gt.1) test=1.0_prec
	alpha2(3)=alpham(3)+acos(test)

	nPG=size(PG,2)
	Alpha=(spread(alpha2-alpha1,dim=2,ncopies=nPG)*spread(PG(1,:),1,3)+ spread(alpha1+alpha2,2,nPG))/2
	r=spread(h,2,nPG) / cos(Alpha-spread(alpham,2,nPG))
	Jac=(alpha2-alpha1)/2

	M0=exp(ji*k*r)*(spread(Jac,2,nPG))
	M1=sum(matmul(M0,PG(2,:)))

	M2=ji/(2.0_prec*k)*(1-M1/(2.0_prec*pi))

	end subroutine DM_sing_T0
!------------------------------------------------------------------
!
! Intégration DM_T0
!
!------------------------------------------------------------------
!	subroutine DM_T0(k,Q,norQ,Jac1,IS,wq,P,norP,D1,M1)
!! calcul les coefs D, M, H et K pour BEM
!	use parametres
!	implicit none
!! variables d'entrée
!	integer, intent(in) :: IS
!	real(prec), intent(in) :: Q(1:3,1:IS), k,Jac1
!	real(prec), intent(in) :: norQ(1:3),norP(3),wq(IS),P(3)
!! variables de sortie
!!	complex(prec), intent(out) :: D1, K1, M1, H1
!	complex(prec), intent(out) :: D1, M1
!! variables locales
!	real(prec) :: R(1:3,IS), Jac
!	real(prec) :: nR(IS), NQR(IS), NPR(IS), Nn
!	complex(prec) :: D0(IS), M0(IS)
!!	complex(prec) :: D0(IS), K0(IS), M0(IS), H0(IS)
!	character(len=1) :: answer
!
!
!! calcul des vecteurs rayons jusqu'aux points de Gauss
!!
!! Calcul de R et norme de R
!!
!	R=Q-spread(P,2,IS)   ! vecteur R
!	nR=sqrt(sum(R**2,1))               ! norme de R
!
!!produits scalaires (rayon.Normale*Jacobien)
!	NQR=matmul(transpose(R),norQ)
!
!! fonction à intégrer
!	M0=exp(ji*k*nR) /(4*pi*nR) * Jac1
!	D0=-(1-ji*k*nR)*exp(ji*k*nR)*NQR /(4*pi*(nR**3)) * Jac1   !% dNQ green
!
!! intégration
!
!	D1=sum(wq*D0)
!	M1=sum(wq*M0)
!
!	end subroutine DM_T0

!-----------------------------------------------------------------------
!
! BEM propa : pout=D pin - M dpin
!
!-----------------------------------------------------------------------
	subroutine bem_const_propa_tot(k,pout,pin,dpin)
!	calcul la matrice creuse Anear pour FMM
	use parametres
	implicit none

! Variables d'entrée
	real(prec), intent(in) :: k
	complex(prec) :: pin(NE),dpin(NE)
! Variables de sortie

	complex(prec), intent(inout) :: pout(Next)

! Variables locales
	integer, dimension(1:3) :: nPe
	real(prec), dimension(1:3,1:3) :: PQ
	real(prec), dimension(1:3,1:IS_e) :: PG_e
	real(prec), dimension(1:3,1:IS_e) :: F_e, dF1_e, dF2_e
	real(prec) :: Qe(1:3,IS_e)

	complex(prec) :: D1,M1
	integer :: i1, i2

! points pour intégration
	call pgt6(PG_e,IS_e)  ! points dans la même cellule
	call FORM_T3(PG_e(1:2,:),IS_e,F_e,dF1_e,dF2_e)
	pout=cmplx(0.0_prec,0.0_prec,kind=prec) ! raz pout

	do i1=1,NE   ! boucle sur point surface
		  nPe=E(i1,:) ! Points de l'élément
		  PQ=transpose(Ps(nPe,:)) ! Coordonnées des points de l'éléments
		  !% calcul de la base tangente naturelle en tous les points de Gauss
		  !a1=matmul(PQ,(/-1,1,0/))
		  !a2=matmul(PQ,(/-1,0,1/))
		  !call vect_product(a1,a2,Nor)
		  Qe=matmul(PQ,F_e)
! points dans la cellule
		  do i2=1,Next ! boucle sur point extérieur

				call DM_T0_propa(k,Qe,norc(i1,:),Jac0(i1),IS_e,PG_e(3,:),Pext(i2,:),D1,M1)
				pout(i2)=pout(i2)+D1*pin(i1)-M1*dpin(i1)
		  end do
	end do
	end subroutine bem_const_propa_tot
!------------------------------------------------------------------------------
!
! Pour BEM propa
!
!------------------------------------------------------------------------------
	subroutine DM_T0_propa(k,Q,norQ,Jac1,IS,wq,P,D1,M1)
! calcul les coefs D et M pour BEM
	use parametres
	implicit none
! variables d'entrée
	integer, intent(in) :: IS
	real(prec), intent(in) :: Q(1:3,IS), k,P(3),norQ(3)
	real(prec), intent(in) :: wq(1:IS),Jac1
! variables de sortie
	complex(prec), intent(out) :: D1, M1
! variables locales
!	complex(prec) :: I
	real(prec) :: R(1:3,IS)
	real(prec) :: nR(IS), NQR(IS) !, NPR(1:size(PG,2)), Nn
	complex(prec) :: D0(IS), M0(IS)

! calcul des vecteurs rayons jusqu'aux points de Gauss
!
! Calcul de R et norme de R
!
	R=Q-spread(P,2,IS)   ! vecteur R
    nR=sqrt(sum(R**2,1))               ! norme de R

!produits scalaires (rayon.Normale*Jacobien)
	NQR=matmul(transpose(R),norQ)

! fonction à intégrer
    M0=exp(ji*k*nR) /(4*pi*(nR)) * Jac1
	D0=-(1-ji*k*nR)*exp(ji*k*nR)*NQR /(4*pi*(nR**3)) * Jac1   !% dN green

! intégration

 	D1=sum(wq*D0)
 	M1=sum(wq*M0)

	end subroutine DM_T0_propa

!------------------------------------------------------------------
!
! points et poids de Gauss Legendre pour intégration
!
!------------------------------------------------------------------
	SUBROUTINE  gauleg(ngp, xabsc, weig)

	implicit none
        INTEGER  i, j, m
        real(prec)  p1, p2, p3, pp, z, z1,EPS,M_PI
        INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
        real(prec), INTENT(OUT) :: xabsc(ngp), weig(ngp)

	EPS=3.0e-15_prec
	M_PI=2*acos(0.0)
	   m = (ngp + 1) / 2
!* Roots are symmetric in the interval - so only need to find half of them  */

	   do i = 1, m				! Loop over the desired roots */

     		z = cos( M_PI * (i-0.25_prec) / (ngp+0.5_prec) )
!*   Starting with the above approximation to the ith root,
!*          we enter the main loop of refinement by NEWTON'S method   */
100     	p1 = 1.0_prec
        	p2 = 0.0_prec
!*  Loop up the recurrence relation to get the Legendre
!*  polynomial evaluated at z                 */

        	do j = 1, ngp
           	p3 = p2
           	p2 = p1
           	p1 = ((2.0_prec*j-1.0_prec) * z * p2 - (j-1.0_prec)*p3) / j
        	enddo

!* p1 is now the desired Legendre polynomial. We next compute pp,
!* its derivative, by a standard relation involving also p2, the
!* polynomial of one lower order.      */
        	pp = ngp*(z*p1-p2)/(z*z-1.0_prec)
        	z1 = z
        	z = z1 - p1/pp             ! Newton's Method  */

        	if (abs(z-z1) .gt. EPS) GOTO  100

      	xabsc(i) =  - z                    	! Roots will be bewteen -1.0 & 1.0 */
      	xabsc(ngp+1-i) =  + z                	! and symmetric about the origin  */
      	weig(i) = 2.0_prec/((1.0_prec-z*z)*pp*pp) ! Compute the weight and its       */
      	weig(ngp+1-i) = weig(i)               ! symmetric counterpart         */

	end do     ! i loop

	End subroutine gauleg
!
!------------------------------------------------------------------
!
! Fonction de forme T3 et dérivées
!
!------------------------------------------------------------------
	subroutine FORM_T6(NP,n1,F,dF1,dF2)
!
! Fonctions de forme T6 quadratiques et leurs derivees, aux points de Gauss NP.
!
! CONVENTION DES NOEUDS : celle du Matlab (Forme6.m / d1Forme.m), qui tourne
! dans le sens trigonometrique en alternant sommet et milieu d'arete :
!
!   1 = sommet A   (0  ,0  )      2 = milieu AB  (1/2,0  )
!   3 = sommet B   (1  ,0  )      4 = milieu BC  (1/2,1/2)
!   5 = sommet C   (0  ,1  )      6 = milieu CA  (0  ,1/2)
!
! GMSH (type 9) ordonne differemment : A B C mAB mBC mCA. La permutation
! E_matlab = E_gmsh(/1,4,2,5,3,6/) est faite UNE fois a la lecture, dans
! lecture_gmsh ; a partir de E tout est dans la convention ci-dessus.
!
	implicit none
	integer, intent(in) :: n1
	real(prec), dimension(1:2,1:n1), intent(in) :: NP
	real(prec), dimension(1:6,1:n1), intent(out) :: F, dF1, dF2
	real(prec), dimension(1:n1) :: x1, x2, x3

	x1 = NP(1,:)
	x2 = NP(2,:)
	x3 = 1.0_prec - x1 - x2

	! fonctions de forme
	F(1,:) = -x3*(1.0_prec-2.0_prec*x3)
	F(2,:) =  4.0_prec*x1*x3
	F(3,:) = -x1*(1.0_prec-2.0_prec*x1)
	F(4,:) =  4.0_prec*x1*x2
	F(5,:) = -x2*(1.0_prec-2.0_prec*x2)
	F(6,:) =  4.0_prec*x2*x3

	! derivee / x1
	dF1(1,:) =  1.0_prec-4.0_prec*x3
	dF1(2,:) =  4.0_prec*(x3-x1)
	dF1(3,:) = -1.0_prec+4.0_prec*x1
	dF1(4,:) =  4.0_prec*x2
	dF1(5,:) =  0.0_prec
	dF1(6,:) = -4.0_prec*x2

	! derivee / x2
	dF2(1,:) =  1.0_prec-4.0_prec*x3
	dF2(2,:) = -4.0_prec*x1
	dF2(3,:) =  0.0_prec
	dF2(4,:) =  4.0_prec*x1
	dF2(5,:) = -1.0_prec+4.0_prec*x2
	dF2(6,:) =  4.0_prec*(x3-x2)

	end subroutine FORM_T6

	subroutine FORM_T3(NP,n1,F,dF1,dF2)
!
! Calcul les fonctions de forme T3 aux points de Gauss NP
!
! Variables entrées
	implicit none
	integer, intent(in)  :: n1
	real(prec), dimension(1:2,1:n1), intent(in) :: NP
!
! Variables sorties
	real(prec), dimension(1:3,1:n1), intent(out) :: F, dF1, dF2
!
! Variables locales
	real(prec), dimension(1:n1) :: x1, x2, x3
!	write(*,*) NP(1,:)
!	write(*,*) NP(2,:)
	x1=NP(1,:)
	x2=NP(2,:)
	x3=1-x1-x2

! Fonctions de forme
	F(1,:)=1-x1-x2
	F(2,:)=x1
	F(3,:)=x2

! Dérivée 1
	dF1(1,:)=-1
	dF1(2,:)=1
	dF1(3,:)=0
! Dérivée 2
	dF2(1,:)=-1
	dF2(2,:)=0
	dF2(3,:)=1

	end subroutine FORM_T3

!
! -----------------------------------------------------------------
!
! Points de Gauss pour intégration sur un triangle
!
!------------------------------------------------------------------
	subroutine pgt6(PtG,IS)
!pgt6(PtG,IS)
!donne les points de Gauss en coordonn�es param�triques
!ainsi que les fonctions de poids (3� ligne)
!d'un élément triangulaire à six noeuds en fonction de
!l'indicateur de sévérité
	implicit none
! variable d'entrée
	integer, intent(in) :: IS
! variable de sortie
	real(prec), dimension(1:3,1:IS), intent(out) :: PtG
! variable locale
	real(prec), dimension(1:4,1:IS) :: P

! GARDE-FOU. Les regles 28, 64, 76 et 112 TELLES QUE SAISIES ICI sont fausses. Elles viennent
! du livre de M. Bonnet (Boundary Integral Equation Methods for Solids and Fluids, 1999), dont
! les tables sont a priori justes : c'est la saisie (coefficients, permutation des orbites ;
! 76 = regle 19 recopiee sur 4 sous-triangles) qui est en cause. Reverifie le 2026-11-02.
! testees par exactitude sur les monomes x^a y^b (2026-09-23, voir
! test_moment_t6), elles donnent
!       IS= 28 : exacte au degre 1 seulement, puis 1e-9 -> 6e-3
!       IS= 64 : exacte au degre 8, alors que 64 points devraient donner bien
!                plus ; erreur 1.5e-6 des le degre 9
!       IS= 76 : exacte au degre 0, erreur 2e-3 des le degre 1
!       IS=112 : exacte au degre 1, puis 2.5e-10 -> 1.8e-4
! Les regles saines sont : 3 (degre 2), 4 (3), 7 (5), 16 (8), 19 (9),
! 27 (11), 33 (12), 37 (13), 48 (15).
! Aucune regle cassee n'est utilisee en production (IS_c=IS_v=33, IS_f=IS_e=16),
! mais les choisir "pour affiner" degraderait silencieusement le resultat.
! On refuse plutot que de renvoyer des poids faux.
	if (IS.eq.28 .or. IS.eq.64 .or. IS.eq.76 .or. IS.eq.112) then
	   write(*,*) 'pgt6 : la regle a ',IS,' points est mal saisie (exactitude polynomiale non verifiee).'
	   write(*,*) '       Regles saines : 3, 4, 7, 16, 19, 27, 33, 37, 48.'
	   stop 1
	endif

	if (IS.eq.3) then
		PtG(1,1)=1.0_prec/6.0_prec
		PtG(2,1)=1.0_prec/6.0_prec
		PtG(3,1)=1.0_prec/6.0_prec
		PtG(1,2)=2.0_prec/3.0_prec
		PtG(2,2)=1.0_prec/6.0_prec
		PtG(3,2)=1.0_prec/6.0_prec
		PtG(1,3)=1.0_prec/6.0_prec
		PtG(2,3)=2.0_prec/3.0_prec
		PtG(3,3)=1.0_prec/6.0_prec
	elseif (IS.eq.4) then
		PtG(1,1)=1.0_prec/3.0_prec
		PtG(2,1)=1.0_prec/3.0_prec
		PtG(3,1)=-27.0_prec/96.0_prec
		PtG(1,2)=1.0_prec/5.0_prec
		PtG(2,2)=1.0_prec/5.0_prec
		PtG(3,2)=25.0_prec/96.0_prec
		PtG(1,3)=3.0_prec/5.0_prec
		PtG(2,3)=1.0_prec/5.0_prec
		PtG(3,3)=25.0_prec/96.0_prec
		PtG(1,4)=1.0_prec/5.0_prec
		PtG(2,4)=3.0_prec/5.0_prec
		PtG(3,4)=25.0_prec/96.0_prec

	elseif (IS.eq.7) then
		P(1,1)=0.470142064105115_prec
		P(2,1)=0.470142064105115_prec
		P(3,1)=1.0_prec-P(1,1)-P(2,1)
   		P(4,1)=0.066197076394253_prec
		P(1,2)=P(2,1)
		P(2,2)=P(3,1)
		P(3,2)=P(1,1)
		P(4,2)=P(4,1)
		P(1,3)=P(2,2)
		P(2,3)=P(3,2)
		P(3,3)=P(1,2)
		P(4,3)=P(4,1)
		P(1,4)=0.101286507323456_prec
		P(2,4)=0.101286507323456_prec
		P(3,4)=1.0_prec-P(1,4)-P(2,4)
		P(4,4)=0.062969590272414_prec
		P(1,5)=P(2,4)
		P(2,5)=P(3,4)
		P(3,5)=P(1,4)
		P(4,5)=P(4,4)
		P(1,6)=P(2,5)
		P(2,6)=P(3,5)
		P(3,6)=P(1,5)
		P(4,6)=P(4,4)
		P(1,7)=1.0_prec/3.0_prec
		P(2,7)=1.0_prec/3.0_prec
		P(3,7)=1.0_prec/3.0_prec
		P(4,7)=0.1125_prec
		PtG(1,:)=P(1,:)
		PtG(2,:)=P(2,:)
		PtG(3,:)=P(4,:)

	elseif (IS.eq.16) then
	   P(1,1)=0.081414823414554_prec
	   P(2,1)=0.459292588292722_prec
	   P(3,1)=1.0_prec-P(1,1)-P(2,1)
	   P(4,1)=0.047545817133642_prec
	   P(1,2)=P(2,1)
	   P(2,2)=P(3,1)
	   P(3,2)=P(1,1)
	   P(4,2)=P(4,1)
	   P(1,3)=P(2,2)
	   P(2,3)=P(3,2)
	   P(3,3)=P(1,2)
	   P(4,3)=P(4,1)
	   P(1,4)=0.898905543365937_prec
	   P(2,4)=0.050547228317031_prec
	   P(3,4)=1.0_prec-P(1,4)-P(2,4)
	   P(4,4)=0.016229248811599_prec
	   P(1,5)=P(2,4)
	   P(2,5)=P(3,4)
	   P(3,5)=P(1,4)
	   P(4,5)=P(4,4)
	   P(1,6)=P(2,5)
	   P(2,6)=P(3,5)
	   P(3,6)=P(1,5)
	   P(4,6)=P(4,4)
	   P(1,7)=1.0_prec/3.0_prec
	   P(2,7)=1.0_prec/3.0_prec
	   P(3,7)=1.0_prec/3.0_prec
	   P(4,7)=0.072157803838893_prec
	   P(1,8)=0.658861384496479_prec
	   P(2,8)=0.170569307751760_prec
	   P(3,8)=1.0_prec-P(1,8)-P(2,8)
	   P(4,8)=0.051608685267359_prec
	   P(1,9)=P(2,8)
	   P(2,9)=P(3,8)
	   P(3,9)=P(1,8)
	   P(4,9)=P(4,8)
	   P(1,10)=P(2,9)
	   P(2,10)=P(3,9)
	   P(3,10)=P(1,9)
	   P(4,10)=P(4,8)
	   P(1,11)=0.008394777409957_prec  !%multiplicit� 6 � discuter
	   P(2,11)=0.728492392955404_prec
	   P(3,11)=1.0_prec-P(1,11)-P(2,11)
	   P(4,11)=0.013615157087217_prec
	   P(1,12)=P(2,11)
	   P(2,12)=P(3,11)
	   P(3,12)=P(1,11)
	   P(4,12)=P(4,11)
	   P(1,13)=P(2,12)
	   P(2,13)=P(3,12)
	   P(3,13)=P(1,12)
	   P(4,13)=P(4,11)
	   P(1,14)=P(2,11)
	   P(2,14)=P(1,11)
	   P(3,14)=P(3,11)
	   P(4,14)=P(4,11)
	   P(1,15)=P(1,11)
	   P(2,15)=P(3,11)
	   P(3,15)=P(2,11)
	   P(4,15)=P(4,11)
	   P(1,16)=P(3,11)
	   P(2,16)=P(2,11)
	   P(3,16)=P(1,11)
	   P(4,16)=P(4,11)
	   PtG(1,:)=P(1,:)
	   PtG(2,:)=P(2,:)
	   PtG(3,:)=P(4,:)

	elseif (IS.eq.19) then ! % modif 06-2004
	   P(1,1)=0.020634961602525_prec
	   P(2,1)=0.489682519198737_prec
	   P(3,1)=1.0_prec-P(1,1)-P(2,1)
	   P(4,1)=0.015667350113570_prec
	   P(1,2)=P(2,1)
	   P(2,2)=P(3,1)
	   P(3,2)=P(1,1)
	   P(4,2)=P(4,1)
	   P(1,3)=P(2,2)
	   P(2,3)=P(3,2)
	   P(3,3)=P(1,2)
	   P(4,3)=P(4,1)
	   P(1,4)=0.910540973211094_prec
	   P(2,4)=0.044729513394452_prec
	   P(3,4)=1.0_prec-P(1,4)-P(2,4)
	   P(4,4)=0.012788837829349_prec
	   P(1,5)=P(2,4)
	   P(2,5)=P(3,4)
	   P(3,5)=P(1,4)
	   P(4,5)=P(4,4)
	   P(1,6)=P(2,5)
   	   P(2,6)=P(3,5)
	   P(3,6)=P(1,5)
	   P(4,6)=P(4,4)
	   P(1,7)=1.0_prec/3.0_prec
	   P(2,7)=1.0_prec/3.0_prec
	   P(3,7)=1.0_prec/3.0_prec
	   P(4,7)=0.048567898141398_prec
	   P(1,8)=0.125820817014129_prec
	   P(2,8)=0.437089591492935_prec
	   P(3,8)=1.0_prec-P(1,8)-P(2,8)
	   P(4,8)=0.038913770502388_prec
	   P(1,9)=P(2,8)
	   P(2,9)=P(3,8)
	   P(3,9)=P(1,8)
	   P(4,9)=P(4,8)
	   P(1,10)=P(2,9)
	   P(2,10)=P(3,9)
	   P(3,10)=P(1,9)
	   P(4,10)=P(4,8)
	   ! FIX : P11-13 etaient des doublons exacts de l'orbite a=0.437090 (P8-10).
	   ! Remplaces par l'orbite MANQUANTE a=0.188203535619033 (Dunavant degre-9) -> poids somment a 0.5.
	   P(1,11)=0.623592928761935_prec
	   P(2,11)=0.188203535619032_prec
	   P(3,11)=1.0_prec-P(1,11)-P(2,11)
	   P(4,11)=0.039823869463604_prec
	   P(1,12)=P(2,11)
	   P(2,12)=P(3,11)
	   P(3,12)=P(1,11)
	   P(4,12)=P(4,11)
	   P(1,13)=P(2,12)
	   P(2,13)=P(3,12)
	   P(3,13)=P(1,12)
	   P(4,13)=P(4,11)

	   P(1,14)=0.036838412054736_prec
	   P(2,14)=0.741198598784498_prec
	   P(3,14)=1.0_prec-P(1,14)-P(2,14)
	   P(4,14)=0.021641769688645_prec
	   P(1,15)=P(2,14)
	   P(2,15)=P(3,14)
	   P(3,15)=P(1,14)
	   P(4,15)=P(4,14)
	   P(1,16)=P(2,15)
	   P(2,16)=P(3,15)
	   P(3,16)=P(1,15)
	   P(4,16)=P(4,14)
	   P(1,17)=P(3,14)
	   P(2,17)=P(2,14)
	   P(3,17)=P(1,14)
	   P(4,17)=P(4,14)
	   P(1,18)=P(2,14)
	   P(2,18)=P(1,14)
	   P(3,18)=P(3,14)
	   P(4,18)=P(4,14)
	   P(1,19)=P(1,14)
	   P(2,19)=P(3,14)
	   P(3,19)=P(2,14)
	   P(4,19)=P(4,14)
	   ! (l'ancienne tentative commentee de l'orbite 0.188204 est maintenant integree en P11-13 ci-dessus)
	   PtG(1,:)=P(1,:)
	   PtG(2,:)=P(2,:)
	   PtG(3,:)=P(4,:)

	elseif (IS.eq.28) then
	   P(1,1)=1.0_prec/3.0_prec
	   P(2,1)=1.0_prec/3.0_prec
	   P(3,1)=1.0_prec/3.0_prec
	   P(4,1)=0.08797730116222190_prec/2.0_prec
	   P(1,2)=0.9480217181434233_prec
	   P(2,2)=0.0258914092828833_prec
	   P(3,2)=1.0_prec-P(1,2)-P(2,2)
	   P(4,2)=0.02623293466120857_prec/6.0_prec
	   P(1,3)=0.8114249947041546_prec
	   P(2,3)=0.09428750264792270_prec
	   P(3,3)=1.0_prec-P(1,3)-P(2,3)
	   P(4,3)=0.1142447159818060_prec/6.0_prec
	   P(1,4)=0.01072644996557060_prec
	   P(2,4)=0.4946367750172147_prec
	   P(3,4)=1.0_prec-P(1,4)-P(2,4)
	   P(4,4)=0.05656634416839376_prec/6.0_prec
	   P(1,5)=0.5853132347709715_prec
	   P(2,5)=0.2073433826145142_prec
	   P(3,5)=1.0_prec-P(1,5)-P(2,5)
	   P(4,5)=0.2164790926342230_prec/6.0_prec
	   P(1,6)=0.1221843885990187_prec
	   P(2,6)=0.4389078057004907_prec
	   P(3,6)=1.0_prec-P(1,6)-P(2,6)
	   P(4,6)=0.2079874161166116_prec/6.0_prec
	   P(1,7)=0.0_prec
	   P(2,7)=0.8588702812826364_prec
	   P(3,7)=1.0_prec-P(1,7)-P(2,7)
	   P(4,7)=0.04417430269980344_prec/12.0_prec
	   P(1,8)=0.04484167758913055_prec
	   P(2,8)=0.6779376548825902_prec
	   P(3,8)=1.0_prec-P(1,8)-P(2,8)
	   P(4,8)=0.2463378925757316_prec/12.0_prec
	   P(1,9)=P(2,2)
	   P(2,9)=P(3,2)
	   P(3,9)=P(1,2)
	   P(4,9)=P(4,2)
	   P(1,10)=P(2,9)
	   P(2,10)=P(3,9)
	   P(3,10)=P(1,9)
	   P(4,10)=P(4,9)
	   P(1,11)=P(2,3)
	   P(2,11)=P(3,3)
	   P(3,11)=P(1,3)
	   P(4,11)=P(4,3)
	   P(1,12)=P(2,11)
	   P(2,12)=P(3,11)
	   P(3,12)=P(1,11)
	   P(4,12)=P(4,11)
	   P(1,13)=P(2,4)
	   P(2,13)=P(3,4)
	   P(3,13)=P(1,4)
	   P(4,13)=P(4,4)
	   P(1,14)=P(2,13)
	   P(2,14)=P(3,13)
	   P(3,14)=P(1,13)
	   P(4,14)=P(4,13)
	   P(1,15)=P(2,5)
	   P(2,15)=P(3,5)
	   P(3,15)=P(1,5)
	   P(4,15)=P(4,5)
	   P(1,16)=P(2,15)
	   P(2,16)=P(3,15)
	   P(3,16)=P(1,15)
	   P(4,16)=P(4,15)
	   P(1,17)=P(2,6)
	   P(2,17)=P(3,6)
	   P(3,17)=P(1,6)
	   P(4,17)=P(4,6)
	   P(1,18)=P(2,17)
	   P(2,18)=P(3,17)
	   P(3,18)=P(1,17)
	   P(4,18)=P(4,17)
	   P(1,19)=P(2,7)
	   P(2,19)=P(3,7)
	   P(3,19)=P(1,7)
	   P(4,19)=P(4,7)
	   P(1,20)=P(3,7)
	   P(2,20)=P(1,7)
	   P(3,20)=P(2,7)
	   P(4,20)=P(4,7)
	   P(1,21)=P(3,7)
	   P(2,21)=P(2,7)
	   P(3,21)=P(1,7)
	   P(4,21)=P(4,7)
	   P(1,22)=P(2,7)
	   P(2,22)=P(1,7)
	   P(3,22)=P(3,7)
	   P(4,22)=P(4,7)
	   P(1,23)=P(1,7)
	   P(2,23)=P(3,7)
	   P(3,23)=P(2,7)
	   P(4,23)=P(4,7)
	   P(1,24)=P(2,8)
	   P(2,24)=P(3,8)
	   P(3,24)=P(1,8)
	   P(4,24)=P(4,8)
	   P(1,25)=P(3,8)
	   P(2,25)=P(1,8)
	   P(3,25)=P(2,8)
	   P(4,25)=P(4,8)
	   P(1,26)=P(3,8)
	   P(2,26)=P(2,8)
	   P(3,26)=P(1,8)
	   P(4,26)=P(4,8)
	   P(1,27)=P(2,8)
	   P(2,27)=P(1,8)
	   P(3,27)=P(3,8)
	   P(4,27)=P(4,8)
	   P(1,28)=P(1,8)
	   P(2,28)=P(3,8)
	   P(3,28)=P(2,8)
	   P(4,28)=P(4,8)
	   PtG(1,:)=P(1,:)
	   PtG(2,:)=P(2,:)
	   PtG(3,:)=P(4,:)

	elseif (IS.eq.27) then
!	   %27 points de Gauss, regle de Dunavant degre 11 (source : Burkardt, dunavant.f90,
!	   %rule 11 ; verifiee independamment par exactitude des monomes jusqu'au degre 11)
	   P(1,1)=0.534611048270758_prec
	   P(2,1)=0.534611048270758_prec
	   P(3,1)=1.0_prec-P(1,1)-P(2,1)
	   P(4,1)=0.000463503164481_prec
	   P(1,2)=0.398969302965855_prec
	   P(2,2)=0.398969302965855_prec
	   P(3,2)=1.0_prec-P(1,2)-P(2,2)
	   P(4,2)=0.038574767457406_prec
	   P(1,3)=0.203309900431282_prec
	   P(2,3)=0.203309900431282_prec
	   P(3,3)=1.0_prec-P(1,3)-P(2,3)
	   P(4,3)=0.029661488690387_prec
	   P(1,4)=0.119350912282581_prec
	   P(2,4)=0.119350912282581_prec
	   P(3,4)=1.0_prec-P(1,4)-P(2,4)
	   P(4,4)=0.018092270251709_prec
	   P(1,5)=0.032364948111276_prec
	   P(2,5)=0.032364948111276_prec
	   P(3,5)=1.0_prec-P(1,5)-P(2,5)
	   P(4,5)=0.006829865501339_prec
	   P(1,6)=0.050178138310495_prec
	   P(2,6)=0.356620648261293_prec
	   P(3,6)=1.0_prec-P(1,6)-P(2,6)
	   P(4,6)=0.026168555981102_prec
	   P(1,7)=0.021022016536166_prec
	   P(2,7)=0.171488980304042_prec
	   P(3,7)=1.0_prec-P(1,7)-P(2,7)
	   P(4,7)=0.010353829819571_prec
	   P(1,8)=P(2,1)
	   P(2,8)=P(3,1)
	   P(3,8)=P(1,1)
	   P(4,8)=P(4,1)
	   P(1,9)=P(2,8)
	   P(2,9)=P(3,8)
	   P(3,9)=P(1,8)
	   P(4,9)=P(4,8)
	   P(1,10)=P(2,2)
	   P(2,10)=P(3,2)
	   P(3,10)=P(1,2)
	   P(4,10)=P(4,2)
	   P(1,11)=P(2,10)
	   P(2,11)=P(3,10)
	   P(3,11)=P(1,10)
	   P(4,11)=P(4,10)
	   P(1,12)=P(2,3)
	   P(2,12)=P(3,3)
	   P(3,12)=P(1,3)
	   P(4,12)=P(4,3)
	   P(1,13)=P(2,12)
	   P(2,13)=P(3,12)
	   P(3,13)=P(1,12)
	   P(4,13)=P(4,12)
	   P(1,14)=P(2,4)
	   P(2,14)=P(3,4)
	   P(3,14)=P(1,4)
	   P(4,14)=P(4,4)
	   P(1,15)=P(2,14)
	   P(2,15)=P(3,14)
	   P(3,15)=P(1,14)
	   P(4,15)=P(4,14)
	   P(1,16)=P(2,5)
	   P(2,16)=P(3,5)
	   P(3,16)=P(1,5)
	   P(4,16)=P(4,5)
	   P(1,17)=P(2,16)
	   P(2,17)=P(3,16)
	   P(3,17)=P(1,16)
	   P(4,17)=P(4,16)
	   P(1,18)=P(2,6)
	   P(2,18)=P(3,6)
	   P(3,18)=P(1,6)
	   P(4,18)=P(4,6)
	   P(1,19)=P(3,6)
	   P(2,19)=P(1,6)
	   P(3,19)=P(2,6)
	   P(4,19)=P(4,6)
	   P(1,20)=P(3,6)
	   P(2,20)=P(2,6)
	   P(3,20)=P(1,6)
	   P(4,20)=P(4,6)
	   P(1,21)=P(2,6)
	   P(2,21)=P(1,6)
	   P(3,21)=P(3,6)
	   P(4,21)=P(4,6)
	   P(1,22)=P(1,6)
	   P(2,22)=P(3,6)
	   P(3,22)=P(2,6)
	   P(4,22)=P(4,6)
	   P(1,23)=P(2,7)
	   P(2,23)=P(3,7)
	   P(3,23)=P(1,7)
	   P(4,23)=P(4,7)
	   P(1,24)=P(3,7)
	   P(2,24)=P(1,7)
	   P(3,24)=P(2,7)
	   P(4,24)=P(4,7)
	   P(1,25)=P(3,7)
	   P(2,25)=P(2,7)
	   P(3,25)=P(1,7)
	   P(4,25)=P(4,7)
	   P(1,26)=P(2,7)
	   P(2,26)=P(1,7)
	   P(3,26)=P(3,7)
	   P(4,26)=P(4,7)
	   P(1,27)=P(1,7)
	   P(2,27)=P(3,7)
	   P(3,27)=P(2,7)
	   P(4,27)=P(4,7)
	   PtG(1,:)=P(1,:)
	   PtG(2,:)=P(2,:)
	   PtG(3,:)=P(4,:)

	elseif (IS.eq.33) then
!	   %33 points de Gauss, regle de Dunavant degre 12 (source : Burkardt, dunavant.f90,
!	   %rule 12 ; verifiee independamment par exactitude des monomes jusqu'au degre 12)
	   P(1,1)=0.488217389773805_prec
	   P(2,1)=0.488217389773805_prec
	   P(3,1)=1.0_prec-P(1,1)-P(2,1)
	   P(4,1)=0.012865533220227_prec
	   P(1,2)=0.439724392294460_prec
	   P(2,2)=0.439724392294460_prec
	   P(3,2)=1.0_prec-P(1,2)-P(2,2)
	   P(4,2)=0.021846272269019_prec
	   P(1,3)=0.271210385012116_prec
	   P(2,3)=0.271210385012116_prec
	   P(3,3)=1.0_prec-P(1,3)-P(2,3)
	   P(4,3)=0.031429112108943_prec
	   P(1,4)=0.127576145541586_prec
	   P(2,4)=0.127576145541586_prec
	   P(3,4)=1.0_prec-P(1,4)-P(2,4)
	   P(4,4)=0.017398056465355_prec
	   P(1,5)=0.021317350453210_prec
	   P(2,5)=0.021317350453210_prec
	   P(3,5)=1.0_prec-P(1,5)-P(2,5)
	   P(4,5)=0.003083130525780_prec
	   P(1,6)=0.115343494534698_prec
	   P(2,6)=0.275713269685514_prec
	   P(3,6)=1.0_prec-P(1,6)-P(2,6)
	   P(4,6)=0.020185778883191_prec
	   P(1,7)=0.022838332222257_prec
	   P(2,7)=0.281325580989940_prec
	   P(3,7)=1.0_prec-P(1,7)-P(2,7)
	   P(4,7)=0.011178386601152_prec
	   P(1,8)=0.025734050548330_prec
	   P(2,8)=0.116251915907597_prec
	   P(3,8)=1.0_prec-P(1,8)-P(2,8)
	   P(4,8)=0.008658115554329_prec
	   P(1,9)=P(2,1)
	   P(2,9)=P(3,1)
	   P(3,9)=P(1,1)
	   P(4,9)=P(4,1)
	   P(1,10)=P(2,9)
	   P(2,10)=P(3,9)
	   P(3,10)=P(1,9)
	   P(4,10)=P(4,9)
	   P(1,11)=P(2,2)
	   P(2,11)=P(3,2)
	   P(3,11)=P(1,2)
	   P(4,11)=P(4,2)
	   P(1,12)=P(2,11)
	   P(2,12)=P(3,11)
	   P(3,12)=P(1,11)
	   P(4,12)=P(4,11)
	   P(1,13)=P(2,3)
	   P(2,13)=P(3,3)
	   P(3,13)=P(1,3)
	   P(4,13)=P(4,3)
	   P(1,14)=P(2,13)
	   P(2,14)=P(3,13)
	   P(3,14)=P(1,13)
	   P(4,14)=P(4,13)
	   P(1,15)=P(2,4)
	   P(2,15)=P(3,4)
	   P(3,15)=P(1,4)
	   P(4,15)=P(4,4)
	   P(1,16)=P(2,15)
	   P(2,16)=P(3,15)
	   P(3,16)=P(1,15)
	   P(4,16)=P(4,15)
	   P(1,17)=P(2,5)
	   P(2,17)=P(3,5)
	   P(3,17)=P(1,5)
	   P(4,17)=P(4,5)
	   P(1,18)=P(2,17)
	   P(2,18)=P(3,17)
	   P(3,18)=P(1,17)
	   P(4,18)=P(4,17)
	   P(1,19)=P(2,6)
	   P(2,19)=P(3,6)
	   P(3,19)=P(1,6)
	   P(4,19)=P(4,6)
	   P(1,20)=P(3,6)
	   P(2,20)=P(1,6)
	   P(3,20)=P(2,6)
	   P(4,20)=P(4,6)
	   P(1,21)=P(3,6)
	   P(2,21)=P(2,6)
	   P(3,21)=P(1,6)
	   P(4,21)=P(4,6)
	   P(1,22)=P(2,6)
	   P(2,22)=P(1,6)
	   P(3,22)=P(3,6)
	   P(4,22)=P(4,6)
	   P(1,23)=P(1,6)
	   P(2,23)=P(3,6)
	   P(3,23)=P(2,6)
	   P(4,23)=P(4,6)
	   P(1,24)=P(2,7)
	   P(2,24)=P(3,7)
	   P(3,24)=P(1,7)
	   P(4,24)=P(4,7)
	   P(1,25)=P(3,7)
	   P(2,25)=P(1,7)
	   P(3,25)=P(2,7)
	   P(4,25)=P(4,7)
	   P(1,26)=P(3,7)
	   P(2,26)=P(2,7)
	   P(3,26)=P(1,7)
	   P(4,26)=P(4,7)
	   P(1,27)=P(2,7)
	   P(2,27)=P(1,7)
	   P(3,27)=P(3,7)
	   P(4,27)=P(4,7)
	   P(1,28)=P(1,7)
	   P(2,28)=P(3,7)
	   P(3,28)=P(2,7)
	   P(4,28)=P(4,7)
	   P(1,29)=P(2,8)
	   P(2,29)=P(3,8)
	   P(3,29)=P(1,8)
	   P(4,29)=P(4,8)
	   P(1,30)=P(3,8)
	   P(2,30)=P(1,8)
	   P(3,30)=P(2,8)
	   P(4,30)=P(4,8)
	   P(1,31)=P(3,8)
	   P(2,31)=P(2,8)
	   P(3,31)=P(1,8)
	   P(4,31)=P(4,8)
	   P(1,32)=P(2,8)
	   P(2,32)=P(1,8)
	   P(3,32)=P(3,8)
	   P(4,32)=P(4,8)
	   P(1,33)=P(1,8)
	   P(2,33)=P(3,8)
	   P(3,33)=P(2,8)
	   P(4,33)=P(4,8)
	   PtG(1,:)=P(1,:)
	   PtG(2,:)=P(2,:)
	   PtG(3,:)=P(4,:)

	elseif (IS.eq.37) then
!	   %37 points de Gauss, regle de Dunavant degre 13 (source : Burkardt, dunavant.f90,
!	   %rule 13 ; verifiee independamment par exactitude des monomes jusqu'au degre 13)
	   P(1,1)=1.0_prec/3.0_prec
	   P(2,1)=1.0_prec/3.0_prec
	   P(3,1)=1.0_prec/3.0_prec
	   P(4,1)=0.026260461700401_prec
	   P(1,2)=0.495048184939705_prec
	   P(2,2)=0.495048184939705_prec
	   P(3,2)=1.0_prec-P(1,2)-P(2,2)
	   P(4,2)=0.005640072604665_prec
	   P(1,3)=0.468716635109574_prec
	   P(2,3)=0.468716635109574_prec
	   P(3,3)=1.0_prec-P(1,3)-P(2,3)
	   P(4,3)=0.015711759181227_prec
	   P(1,4)=0.414521336801277_prec
	   P(2,4)=0.414521336801277_prec
	   P(3,4)=1.0_prec-P(1,4)-P(2,4)
	   P(4,4)=0.023536251252097_prec
	   P(1,5)=0.229399572042831_prec
	   P(2,5)=0.229399572042831_prec
	   P(3,5)=1.0_prec-P(1,5)-P(2,5)
	   P(4,5)=0.023681793268178_prec
	   P(1,6)=0.114424495196330_prec
	   P(2,6)=0.114424495196330_prec
	   P(3,6)=1.0_prec-P(1,6)-P(2,6)
	   P(4,6)=0.015583764522897_prec
	   P(1,7)=0.024811391363459_prec
	   P(2,7)=0.024811391363459_prec
	   P(3,7)=1.0_prec-P(1,7)-P(2,7)
	   P(4,7)=0.003987885732537_prec
	   P(1,8)=0.094853828379579_prec
	   P(2,8)=0.268794997058761_prec
	   P(3,8)=1.0_prec-P(1,8)-P(2,8)
	   P(4,8)=0.018424201364366_prec
	   P(1,9)=0.018100773278807_prec
	   P(2,9)=0.291730066734288_prec
	   P(3,9)=1.0_prec-P(1,9)-P(2,9)
	   P(4,9)=0.008700731651911_prec
	   P(1,10)=0.022233076674090_prec
	   P(2,10)=0.126357385491669_prec
	   P(3,10)=1.0_prec-P(1,10)-P(2,10)
	   P(4,10)=0.007760893419522_prec
	   P(1,11)=P(2,2)
	   P(2,11)=P(3,2)
	   P(3,11)=P(1,2)
	   P(4,11)=P(4,2)
	   P(1,12)=P(2,11)
	   P(2,12)=P(3,11)
	   P(3,12)=P(1,11)
	   P(4,12)=P(4,11)
	   P(1,13)=P(2,3)
	   P(2,13)=P(3,3)
	   P(3,13)=P(1,3)
	   P(4,13)=P(4,3)
	   P(1,14)=P(2,13)
	   P(2,14)=P(3,13)
	   P(3,14)=P(1,13)
	   P(4,14)=P(4,13)
	   P(1,15)=P(2,4)
	   P(2,15)=P(3,4)
	   P(3,15)=P(1,4)
	   P(4,15)=P(4,4)
	   P(1,16)=P(2,15)
	   P(2,16)=P(3,15)
	   P(3,16)=P(1,15)
	   P(4,16)=P(4,15)
	   P(1,17)=P(2,5)
	   P(2,17)=P(3,5)
	   P(3,17)=P(1,5)
	   P(4,17)=P(4,5)
	   P(1,18)=P(2,17)
	   P(2,18)=P(3,17)
	   P(3,18)=P(1,17)
	   P(4,18)=P(4,17)
	   P(1,19)=P(2,6)
	   P(2,19)=P(3,6)
	   P(3,19)=P(1,6)
	   P(4,19)=P(4,6)
	   P(1,20)=P(2,19)
	   P(2,20)=P(3,19)
	   P(3,20)=P(1,19)
	   P(4,20)=P(4,19)
	   P(1,21)=P(2,7)
	   P(2,21)=P(3,7)
	   P(3,21)=P(1,7)
	   P(4,21)=P(4,7)
	   P(1,22)=P(2,21)
	   P(2,22)=P(3,21)
	   P(3,22)=P(1,21)
	   P(4,22)=P(4,21)
	   P(1,23)=P(2,8)
	   P(2,23)=P(3,8)
	   P(3,23)=P(1,8)
	   P(4,23)=P(4,8)
	   P(1,24)=P(3,8)
	   P(2,24)=P(1,8)
	   P(3,24)=P(2,8)
	   P(4,24)=P(4,8)
	   P(1,25)=P(3,8)
	   P(2,25)=P(2,8)
	   P(3,25)=P(1,8)
	   P(4,25)=P(4,8)
	   P(1,26)=P(2,8)
	   P(2,26)=P(1,8)
	   P(3,26)=P(3,8)
	   P(4,26)=P(4,8)
	   P(1,27)=P(1,8)
	   P(2,27)=P(3,8)
	   P(3,27)=P(2,8)
	   P(4,27)=P(4,8)
	   P(1,28)=P(2,9)
	   P(2,28)=P(3,9)
	   P(3,28)=P(1,9)
	   P(4,28)=P(4,9)
	   P(1,29)=P(3,9)
	   P(2,29)=P(1,9)
	   P(3,29)=P(2,9)
	   P(4,29)=P(4,9)
	   P(1,30)=P(3,9)
	   P(2,30)=P(2,9)
	   P(3,30)=P(1,9)
	   P(4,30)=P(4,9)
	   P(1,31)=P(2,9)
	   P(2,31)=P(1,9)
	   P(3,31)=P(3,9)
	   P(4,31)=P(4,9)
	   P(1,32)=P(1,9)
	   P(2,32)=P(3,9)
	   P(3,32)=P(2,9)
	   P(4,32)=P(4,9)
	   P(1,33)=P(2,10)
	   P(2,33)=P(3,10)
	   P(3,33)=P(1,10)
	   P(4,33)=P(4,10)
	   P(1,34)=P(3,10)
	   P(2,34)=P(1,10)
	   P(3,34)=P(2,10)
	   P(4,34)=P(4,10)
	   P(1,35)=P(3,10)
	   P(2,35)=P(2,10)
	   P(3,35)=P(1,10)
	   P(4,35)=P(4,10)
	   P(1,36)=P(2,10)
	   P(2,36)=P(1,10)
	   P(3,36)=P(3,10)
	   P(4,36)=P(4,10)
	   P(1,37)=P(1,10)
	   P(2,37)=P(3,10)
	   P(3,37)=P(2,10)
	   P(4,37)=P(4,10)
	   PtG(1,:)=P(1,:)
	   PtG(2,:)=P(2,:)
	   PtG(3,:)=P(4,:)

	elseif (IS.eq.48) then
!	   %48 points de Gauss, regle de Dunavant degre 15 (source : Burkardt, dunavant.f90,
!	   %rule 15 ; verifiee independamment par exactitude des monomes jusqu'au degre 15)
	   P(1,1)=0.506972916858243_prec
	   P(2,1)=0.506972916858243_prec
	   P(3,1)=1.0_prec-P(1,1)-P(2,1)
	   P(4,1)=0.000958437821424_prec
	   P(1,2)=0.431406354283023_prec
	   P(2,2)=0.431406354283023_prec
	   P(3,2)=1.0_prec-P(1,2)-P(2,2)
	   P(4,2)=0.022124513635572_prec
	   P(1,3)=0.277693644847144_prec
	   P(2,3)=0.277693644847144_prec
	   P(3,3)=1.0_prec-P(1,3)-P(2,3)
	   P(4,3)=0.025593274359426_prec
	   P(1,4)=0.126464891041254_prec
	   P(2,4)=0.126464891041254_prec
	   P(3,4)=1.0_prec-P(1,4)-P(2,4)
	   P(4,4)=0.011843867935344_prec
	   P(1,5)=0.070808385974686_prec
	   P(2,5)=0.070808385974686_prec
	   P(3,5)=1.0_prec-P(1,5)-P(2,5)
	   P(4,5)=0.006644887845011_prec
	   P(1,6)=0.018965170241073_prec
	   P(2,6)=0.018965170241073_prec
	   P(3,6)=1.0_prec-P(1,6)-P(2,6)
	   P(4,6)=0.002374458304096_prec
	   P(1,7)=0.133734161966621_prec
	   P(2,7)=0.261311371140087_prec
	   P(3,7)=1.0_prec-P(1,7)-P(2,7)
	   P(4,7)=0.019275036299796_prec
	   P(1,8)=0.036366677396917_prec
	   P(2,8)=0.388046767090269_prec
	   P(3,8)=1.0_prec-P(1,8)-P(2,8)
	   P(4,8)=0.013607907160312_prec
	   P(1,9)=-0.010174883126571_prec
	   P(2,9)=0.285712220049916_prec
	   P(3,9)=1.0_prec-P(1,9)-P(2,9)
	   P(4,9)=0.001091038683399_prec
	   P(1,10)=0.036843869875878_prec
	   P(2,10)=0.215599664072284_prec
	   P(3,10)=1.0_prec-P(1,10)-P(2,10)
	   P(4,10)=0.010752659923865_prec
	   P(1,11)=0.012459809331199_prec
	   P(2,11)=0.103575616576386_prec
	   P(3,11)=1.0_prec-P(1,11)-P(2,11)
	   P(4,11)=0.003836971315524_prec
	   P(1,12)=P(2,1)
	   P(2,12)=P(3,1)
	   P(3,12)=P(1,1)
	   P(4,12)=P(4,1)
	   P(1,13)=P(2,12)
	   P(2,13)=P(3,12)
	   P(3,13)=P(1,12)
	   P(4,13)=P(4,12)
	   P(1,14)=P(2,2)
	   P(2,14)=P(3,2)
	   P(3,14)=P(1,2)
	   P(4,14)=P(4,2)
	   P(1,15)=P(2,14)
	   P(2,15)=P(3,14)
	   P(3,15)=P(1,14)
	   P(4,15)=P(4,14)
	   P(1,16)=P(2,3)
	   P(2,16)=P(3,3)
	   P(3,16)=P(1,3)
	   P(4,16)=P(4,3)
	   P(1,17)=P(2,16)
	   P(2,17)=P(3,16)
	   P(3,17)=P(1,16)
	   P(4,17)=P(4,16)
	   P(1,18)=P(2,4)
	   P(2,18)=P(3,4)
	   P(3,18)=P(1,4)
	   P(4,18)=P(4,4)
	   P(1,19)=P(2,18)
	   P(2,19)=P(3,18)
	   P(3,19)=P(1,18)
	   P(4,19)=P(4,18)
	   P(1,20)=P(2,5)
	   P(2,20)=P(3,5)
	   P(3,20)=P(1,5)
	   P(4,20)=P(4,5)
	   P(1,21)=P(2,20)
	   P(2,21)=P(3,20)
	   P(3,21)=P(1,20)
	   P(4,21)=P(4,20)
	   P(1,22)=P(2,6)
	   P(2,22)=P(3,6)
	   P(3,22)=P(1,6)
	   P(4,22)=P(4,6)
	   P(1,23)=P(2,22)
	   P(2,23)=P(3,22)
	   P(3,23)=P(1,22)
	   P(4,23)=P(4,22)
	   P(1,24)=P(2,7)
	   P(2,24)=P(3,7)
	   P(3,24)=P(1,7)
	   P(4,24)=P(4,7)
	   P(1,25)=P(3,7)
	   P(2,25)=P(1,7)
	   P(3,25)=P(2,7)
	   P(4,25)=P(4,7)
	   P(1,26)=P(3,7)
	   P(2,26)=P(2,7)
	   P(3,26)=P(1,7)
	   P(4,26)=P(4,7)
	   P(1,27)=P(2,7)
	   P(2,27)=P(1,7)
	   P(3,27)=P(3,7)
	   P(4,27)=P(4,7)
	   P(1,28)=P(1,7)
	   P(2,28)=P(3,7)
	   P(3,28)=P(2,7)
	   P(4,28)=P(4,7)
	   P(1,29)=P(2,8)
	   P(2,29)=P(3,8)
	   P(3,29)=P(1,8)
	   P(4,29)=P(4,8)
	   P(1,30)=P(3,8)
	   P(2,30)=P(1,8)
	   P(3,30)=P(2,8)
	   P(4,30)=P(4,8)
	   P(1,31)=P(3,8)
	   P(2,31)=P(2,8)
	   P(3,31)=P(1,8)
	   P(4,31)=P(4,8)
	   P(1,32)=P(2,8)
	   P(2,32)=P(1,8)
	   P(3,32)=P(3,8)
	   P(4,32)=P(4,8)
	   P(1,33)=P(1,8)
	   P(2,33)=P(3,8)
	   P(3,33)=P(2,8)
	   P(4,33)=P(4,8)
	   P(1,34)=P(2,9)
	   P(2,34)=P(3,9)
	   P(3,34)=P(1,9)
	   P(4,34)=P(4,9)
	   P(1,35)=P(3,9)
	   P(2,35)=P(1,9)
	   P(3,35)=P(2,9)
	   P(4,35)=P(4,9)
	   P(1,36)=P(3,9)
	   P(2,36)=P(2,9)
	   P(3,36)=P(1,9)
	   P(4,36)=P(4,9)
	   P(1,37)=P(2,9)
	   P(2,37)=P(1,9)
	   P(3,37)=P(3,9)
	   P(4,37)=P(4,9)
	   P(1,38)=P(1,9)
	   P(2,38)=P(3,9)
	   P(3,38)=P(2,9)
	   P(4,38)=P(4,9)
	   P(1,39)=P(2,10)
	   P(2,39)=P(3,10)
	   P(3,39)=P(1,10)
	   P(4,39)=P(4,10)
	   P(1,40)=P(3,10)
	   P(2,40)=P(1,10)
	   P(3,40)=P(2,10)
	   P(4,40)=P(4,10)
	   P(1,41)=P(3,10)
	   P(2,41)=P(2,10)
	   P(3,41)=P(1,10)
	   P(4,41)=P(4,10)
	   P(1,42)=P(2,10)
	   P(2,42)=P(1,10)
	   P(3,42)=P(3,10)
	   P(4,42)=P(4,10)
	   P(1,43)=P(1,10)
	   P(2,43)=P(3,10)
	   P(3,43)=P(2,10)
	   P(4,43)=P(4,10)
	   P(1,44)=P(2,11)
	   P(2,44)=P(3,11)
	   P(3,44)=P(1,11)
	   P(4,44)=P(4,11)
	   P(1,45)=P(3,11)
	   P(2,45)=P(1,11)
	   P(3,45)=P(2,11)
	   P(4,45)=P(4,11)
	   P(1,46)=P(3,11)
	   P(2,46)=P(2,11)
	   P(3,46)=P(1,11)
	   P(4,46)=P(4,11)
	   P(1,47)=P(2,11)
	   P(2,47)=P(1,11)
	   P(3,47)=P(3,11)
	   P(4,47)=P(4,11)
	   P(1,48)=P(1,11)
	   P(2,48)=P(3,11)
	   P(3,48)=P(2,11)
	   P(4,48)=P(4,11)
	   PtG(1,:)=P(1,:)
	   PtG(2,:)=P(2,:)
	   PtG(3,:)=P(4,:)

	elseif (IS.eq.64) then
!	   %64 points de gauss calcul�s � partir de 16, ramen�s sur un triangle
!	   %deux fois plus petit que l'on reporte 4 fois sur le triangle unit�
	  !write(*,*) 'dans pgt6'
	   P(1,1)=0.081414823414554_prec
	   P(2,1)=0.459292588292722_prec
	   P(3,1)=1.0_prec-P(1,1)-P(2,1)
	   P(4,1)=0.047545817133642_prec
	   P(1,2)=P(2,1)
	   P(2,2)=P(3,1)
	   P(3,2)=P(1,1)
	   P(4,2)=P(4,1)
	   P(1,3)=P(2,2)
	   P(2,3)=P(3,2)
	   P(3,3)=P(1,2)
	   P(4,3)=P(4,1)
	   P(1,4)=0.898905543365937_prec
	   P(2,4)=0.050547228317031_prec
	   P(3,4)=1.0_prec-P(1,4)-P(2,4)
	   P(4,4)=0.016229248811599_prec
	   P(1,5)=P(2,4)
	   P(2,5)=P(3,4)
	   P(3,5)=P(1,4)
	   P(4,5)=P(4,4)
	   P(1,6)=P(2,5)
	   P(2,6)=P(3,5)
	   P(3,6)=P(1,5)
	   P(4,6)=P(4,4)
	   P(1,7)=1.0_prec/3.0_prec
	   P(2,7)=1.0_prec/3.0_prec
	   P(3,7)=1.0_prec/3.0_prec
	   P(4,7)=0.072157803838893_prec
	   P(1,8)=0.658861384496479_prec
	   P(2,8)=0.170569307751760_prec
	   P(3,8)=1.0_prec-P(1,8)-P(2,8)
	   P(4,8)=0.051608685267359_prec
	   P(1,9)=P(2,8)
	   P(2,9)=P(3,8)
	   P(3,9)=P(1,8)
	   P(4,9)=P(4,8)
	   P(1,10)=P(2,9)
	   P(2,10)=P(3,9)
	   P(3,10)=P(1,9)
	   P(4,10)=P(4,8)
	   P(1,11)=0.008394777409957_prec  !%multiplicit� 6 � discuter
	   P(2,11)=0.728492392955404_prec
	   P(3,11)=1.0_prec-P(1,11)-P(2,11)
	   P(4,11)=0.013615157087217_prec
	   P(1,12)=P(2,11)
	   P(2,12)=P(3,11)
	   P(3,12)=P(1,11)
	   P(4,12)=P(4,11)
	   P(1,13)=P(2,12)
	   P(2,13)=P(3,12)
	   P(3,13)=P(1,12)
	   P(4,13)=P(4,11)
	   P(1,14)=P(2,11)
	   P(2,14)=P(1,11)
	   P(3,14)=P(3,11)
	   P(4,14)=P(4,11)
	   P(1,15)=P(1,11)
	   P(2,15)=P(3,11)
	   P(3,15)=P(2,11)
	   P(4,15)=P(4,11)
	   P(1,16)=P(3,11)
	   P(2,16)=P(2,11)
	   P(3,16)=P(1,11)
	   P(4,16)=P(4,11)

	   P(1,:)=P(1,:)/2.0_prec
	   P(2,:)=P(2,:)/2.0_prec
	   P(4,:)=P(4,:)/4.0_prec
	   !write(*,*) 'après 16'
	   !write(*,*) P(1,:)
	   PtG(1,:)=(/ P(1,1:16),(P(1,1:16)+0.5_prec), P(1,1:16), (0.5_prec-P(1,1:16)) /)
	   PtG(2,:)=(/ P(2,1:16), P(2,1:16),  (P(2,1:16)+0.5_prec), (0.5_prec-P(2,1:16)) /)
	   PtG(3,:)=(/ P(4,1:16), P(4,1:16),  P(4,1:16),  P(4,1:16) /)
	   !write(*,*) 'fin permutation'

	elseif (IS.eq.76) then
!	   %76 points de gauss calcul�s � partir de 19, ramen�s sur un triangle
!	   %deux fois plus petit que l'on reporte 4 fois sur le triangle unit�
	   P(1,1)=0.020634961602525_prec
	   P(2,1)=0.489682519198737_prec
	   P(3,1)=1.0_prec-P(1,1)-P(2,1)
	   P(4,1)=0.015667350113570_prec
	   P(1,2)=P(2,1)
	   P(2,2)=P(3,1)
	   P(3,2)=P(1,1)
	   P(4,2)=P(4,1)
	   P(1,3)=P(2,2)
	   P(2,3)=P(3,2)
	   P(3,3)=P(1,2)
	   P(4,3)=P(4,1)
	   P(1,4)=0.910540973211094_prec
	   P(2,4)=0.044729513394452_prec
	   P(3,4)=1.0_prec-P(1,4)-P(2,4)
	   P(4,4)=0.012788837829349_prec
	   P(1,5)=P(2,4)
	   P(2,5)=P(3,4)
	   P(3,5)=P(1,4)
	   P(4,5)=P(4,4)
	   P(1,6)=P(2,5)
	   P(2,6)=P(3,5)
	   P(3,6)=P(1,5)
	   P(4,6)=P(4,4)
	   P(1,7)=1.0_prec/3.0_prec
	   P(2,7)=1.0_prec/3.0_prec
	   P(3,7)=1.0_prec/3.0_prec
	   P(4,7)=0.048567898141398_prec
	   P(1,8)=0.125820817014129_prec
	   P(2,8)=0.437089591492935_prec
	   P(3,8)=1.0_prec-P(1,8)-P(2,8)
	   P(4,8)=0.038913770502388_prec
	   P(1,9)=P(2,8)
	   P(2,9)=P(3,8)
	   P(3,9)=P(1,8)
	   P(4,9)=P(4,8)
	   P(1,10)=P(2,9)
	   P(2,10)=P(3,9)
	   P(3,10)=P(1,9)
	   P(4,10)=P(4,8)
	   P(1,11)=0.036838412054736_prec
	   P(2,11)=0.741198598784498_prec
	   P(3,11)=1.0_prec-P(1,11)-P(2,11)
	   P(4,11)=0.021641769688645_prec
	   P(1,12)=P(2,11)
	   P(2,12)=P(3,11)
	   P(3,12)=P(1,11)
	   P(4,12)=P(4,11)
	   P(1,13)=P(2,12)
	   P(2,13)=P(3,12)
	   P(3,13)=P(1,12)
	   P(4,13)=P(4,11)
	   P(1,14)=P(3,11)
	   P(2,14)=P(2,11)
	   P(3,14)=P(1,11)
	   P(4,14)=P(4,11)
	   P(1,15)=P(2,11)
	   P(2,15)=P(1,11)
	   P(3,15)=P(3,11)
	   P(4,15)=P(4,11)
	   P(1,16)=P(1,11)
	   P(2,16)=P(3,11)
	   P(3,16)=P(2,11)
	   P(4,16)=P(4,11)
	   P(1,17)=0.623592928761935_prec
	   P(2,17)=0.188203535619032_prec
	   P(3,17)=1.0_prec-P(1,11)-P(2,11)
	   P(4,17)=0.039823869463604_prec
	   P(1,18)=P(2,17)
	   P(2,18)=P(3,17)
	   P(3,18)=P(1,17)
	   P(4,18)=P(4,17)
	   P(1,19)=P(2,18)
	   P(2,19)=P(3,18)
	   P(3,19)=P(1,18)
	   P(4,19)=P(4,18)

	   P(1,:)=P(1,:)/2.0_prec
	   P(2,:)=P(2,:)/2.0_prec
	   P(4,:)=P(4,:)/4.0_prec

	   PtG(1,:)=(/ P(1,1:19),(P(1,1:19)+0.5_prec),P(1,1:19),(0.5_prec-P(1,1:19)) /)
	   PtG(2,:)=(/ P(2,1:19), P(2,1:19), (P(2,1:19)+0.5_prec), (0.5_prec-P(2,1:19)) /)
	   PtG(3,:)=(/ P(4,1:19), P(4,1:19), P(4,1:19), P(4,1:19) /)

	elseif (IS.eq.112) then
!   %112 points de gauss calcul�s � partir de 28, ramen�s sur un triangle
!   %deux fois plus petit que l'on reporte 4 fois sur le triangle unit�
	   P(1,1)=1.0_prec/3.0_prec
	   P(2,1)=1.0_prec/3.0_prec
	   P(3,1)=1.0_prec/3.0_prec
	   P(4,1)=0.08797730116222190_prec/2.0_prec
	   P(1,2)=0.9480217181434233_prec
	   P(2,2)=0.0258914092828833_prec
	   P(3,2)=1.0-P(1,2)-P(2,2)
	   P(4,2)=0.02623293466120857_prec/6.0_prec
	   P(1,3)=0.8114249947041546_prec
	   P(2,3)=0.09428750264792270_prec
	   P(3,3)=1.0-P(1,3)-P(2,3)
	   P(4,3)=0.1142447159818060_prec/6.0_prec
	   P(1,4)=0.01072644996557060_prec
	   P(2,4)=0.4946367750172147_prec
	   P(3,4)=1.0-P(1,4)-P(2,4)
	   P(4,4)=0.05656634416839376_prec/6.0_prec
	   P(1,5)=0.5853132347709715_prec
	   P(2,5)=0.2073433826145142_prec
	   P(3,5)=1.0-P(1,5)-P(2,5)
	   P(4,5)=0.2164790926342230_prec/6.0_prec
	   P(1,6)=0.1221843885990187_prec
	   P(2,6)=0.4389078057004907_prec
	   P(3,6)=1.0-P(1,6)-P(2,6)
	   P(4,6)=0.2079874161166116_prec/6.0_prec
	   P(1,7)=0.0
	   P(2,7)=0.8588702812826364_prec
	   P(3,7)=1.0-P(1,7)-P(2,7)
	   P(4,7)=0.04417430269980344_prec/12.0_prec
	   P(1,8)=0.04484167758913055_prec
	   P(2,8)=0.6779376548825902_prec
	   P(3,8)=1.0-P(1,8)-P(2,8)
	   P(4,8)=0.2463378925757316_prec/12.0_prec
	   P(1,9)=P(2,2)
	   P(2,9)=P(3,2)
	   P(3,9)=P(1,2)
	   P(4,9)=P(4,2)
	   P(1,10)=P(2,9)
	   P(2,10)=P(3,9)
	   P(3,10)=P(1,9)
	   P(4,10)=P(4,9)
	   P(1,11)=P(2,3)
	   P(2,11)=P(3,3)
	   P(3,11)=P(1,3)
	   P(4,11)=P(4,3)
	   P(1,12)=P(2,11)
	   P(2,12)=P(3,11)
	   P(3,12)=P(1,11)
	   P(4,12)=P(4,11)
	   P(1,13)=P(2,4)
	   P(2,13)=P(3,4)
	   P(3,13)=P(1,4)
	   P(4,13)=P(4,4)
	   P(1,14)=P(2,13)
	   P(2,14)=P(3,13)
	   P(3,14)=P(1,13)
	   P(4,14)=P(4,13)
	   P(1,15)=P(2,5)
	   P(2,15)=P(3,5)
	   P(3,15)=P(1,5)
	   P(4,15)=P(4,5)
	   P(1,16)=P(2,15)
	   P(2,16)=P(3,15)
	   P(3,16)=P(1,15)
	   P(4,16)=P(4,15)
	   P(1,17)=P(2,6)
	   P(2,17)=P(3,6)
	   P(3,17)=P(1,6)
	   P(4,17)=P(4,6)
	   P(1,18)=P(2,17)
	   P(2,18)=P(3,17)
	   P(3,18)=P(1,17)
	   P(4,18)=P(4,17)
	   P(1,19)=P(2,7)
	   P(2,19)=P(3,7)
	   P(3,19)=P(1,7)
	   P(4,19)=P(4,7)
	   P(1,20)=P(3,7)
	   P(2,20)=P(1,7)
	   P(3,20)=P(2,7)
	   P(4,20)=P(4,7)
	   P(1,21)=P(3,7)
	   P(2,21)=P(2,7)
	   P(3,21)=P(1,7)
	   P(4,21)=P(4,7)
	   P(1,22)=P(2,7)
	   P(2,22)=P(1,7)
	   P(3,22)=P(3,7)
	   P(4,22)=P(4,7)
	   P(1,23)=P(1,7)
	   P(2,23)=P(3,7)
	   P(3,23)=P(2,7)
	   P(4,23)=P(4,7)
	   P(1,24)=P(2,8)
	   P(2,24)=P(3,8)
	   P(3,24)=P(1,8)
	   P(4,24)=P(4,8)
	   P(1,25)=P(3,8)
	   P(2,25)=P(1,8)
	   P(3,25)=P(2,8)
	   P(4,25)=P(4,8)
	   P(1,26)=P(3,8)
	   P(2,26)=P(2,8)
	   P(3,26)=P(1,8)
	   P(4,26)=P(4,8)
	   P(1,27)=P(2,8)
	   P(2,27)=P(1,8)
	   P(3,27)=P(3,8)
	   P(4,27)=P(4,8)
	   P(1,28)=P(1,8)
	   P(2,28)=P(3,8)
	   P(3,28)=P(2,8)
	   P(4,28)=P(4,8)

	   P(1,:)=P(1,:)/2.0_prec
	   P(2,:)=P(2,:)/2.0_prec
	   P(4,:)=P(4,:)/4.0_prec

	   PtG(1,:)=(/ P(1,1:28),(P(1,1:28)+0.5_prec),P(1,1:28),(0.5_prec-P(1,1:28)) /)
	   PtG(2,:)=(/ P(2,1:28), P(2,1:28), (P(2,1:28)+0.5_prec), (0.5_prec-P(2,1:28)) /)
	   PtG(3,:)=(/ P(4,1:28), P(4,1:28), P(4,1:28), P(4,1:28) /)

!	elseif ((IS<4) & (IS>=3.5))
!   % test d'autres points sur un triangle : 7 points
!
!	   P(1,1)=0.470142064105115
!	   P(2,1)=0.470142064105115
!	   P(3,1)=1-P(1,1)-P(2,1)
!	   P(4,1)=-0.0083
!	   P(1,2)=P(2,1)
!	   P(2,2)=P(3,1)
!	   P(3,2)=P(1,1)
!	   P(4,2)=P(4,1)
!	   P(1,3)=P(2,2)
!	   P(2,3)=P(3,2)
!	   P(3,3)=P(1,2)
!	   P(4,3)=P(4,1)
!	   P(1,4)=0.101286507323456
!	   P(2,4)=0.101286507323456
!	   P(3,4)=1-P(1,4)-P(2,4)
!	   P(4,4)=0.2178
!	   P(1,5)=P(2,4)
!	   P(2,5)=P(3,4)
!	   P(3,5)=P(1,4)
!	   P(4,5)=P(4,4)
!	   P(1,6)=P(2,5)
!	   P(2,6)=P(3,5)
!	   P(3,6)=P(1,5)
!	   P(4,6)=P(4,4)
!	   P(1,7)=1/3
!	   P(2,7)=1/3
!	   P(3,7)=1/3
!	   P(4,7)=-0.1287
!	   PtG=[P(1,:)P(2,:)P(4,:)]
	else
	   write(*,*) 'pgt6 : pas de regle a ',IS,' points.'
	   write(*,*) '       Regles disponibles : 3, 4, 7, 16, 19, 27, 33, 37, 48.'
	   stop 1
	endif

	end subroutine pgt6



    end module module_BEM
