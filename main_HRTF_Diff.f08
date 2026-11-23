! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: MIT
	PROGRAM main_HRTF_Diff
	! Programme MLFMM pour le calcul des HRTFs, en T0 (constant sur élément)
	! Calcul selon Gumerov (JASA)
	use omp_lib
	use parametres
	use module_Acoustiques
	use module_BEM
	use module_Octree
	use module_GMSH
	use module_MLFMM
	IMPLICIT NONE

	character(len=16) :: env_ob          ! OPENBLAS_NUM_THREADS
	integer :: st_ob
	! lecture des fichiers
	integer :: fileunit,arglen,stat
	character(len=15) :: varBid ! variable bidon pour lire les en-tête du fichier texte
	character(:), allocatable :: filename
	character(len=500) :: file_msh,file_save,file_Pext,file_log

	! octree
	type(tree),dimension(:),allocatable :: octree1
	type(FMM0) :: SPH0
	! parametres physiques
	double precision, dimension(:), allocatable :: freq,k,Refl
	double complex, allocatable :: Ared(:)
	character(len=500) :: ligne_refl        ! « Rmax_BF Rmax_HF » ou nom d'un fichier d'admittance (2026-10-07)
	character(len=500) :: ligne9          ! 9e valeur de la config : verbose (ordre Ray) ou impedance (ancien ordre)
	logical :: ancien_ordre
	! champs de pressions
	double precision :: V,r0(1,3),r(1,3)
	double precision, allocatable :: S(:,:)
	double complex, allocatable :: dp0(:,:),dp1(:,:),pref(:,:),pref0(:,:)
	double complex, allocatable :: p1(:,:),p2(:,:),dp2(:,:),p0(:,:)

	! variables locales
	integer :: i1,i2,ind0,iter,deb,fin
	integer :: NL,N,num_core,HRTF_norm   ! NLF vient de parametres (=2)
	! octree ADAPTATIF (ligne NL de la config a -1) : plus grande arete de chaque
	! element, niveaux (voir NL_adaptatif_T0)
	double precision, allocatable :: h_el(:)
	integer :: nhors, ilv
	double precision, allocatable :: Ptous(:,:)   ! CHI adaptatif : [Pc ; Pt] (sources puis cibles)
	double precision :: seuil_adapt           ! cellule >= seuil_adapt x arete (fichier technique)
	character(len=500) :: file_param          ! fichier des parametres techniques MLFMM
	character(len=20000) :: ligne_freq        ! ligne des frequences (une valeur ou une liste ; 512 valeurs tiennent)
	double precision, allocatable :: liste_freq(:)
	integer :: nliste
	double precision :: f_boost
	integer :: dN_BF_lu
	real(prec) :: c_son                 ! celerite du son (m/s), lue dans la config (derniere ligne)
	integer :: ios_c
	logical :: ob_ok                     ! OPENBLAS_NUM_THREADS = 1
	integer :: dt(8)                     ! date du journal
	double precision :: t_f              ! duree d'une frequence (journal)
	double precision :: hmin_m, hmax_m   ! plus petite / plus grande arete du maillage
	double precision :: t1,t2,PT0(3)
	double complex, allocatable :: sol2(:),sol0(:,:),sol1(:,:),sol3(:)
	integer, allocatable :: iter0(:),ordre0(:,:),ordre1(:),opt_mem(:)
	! estimation initiale GMRES par blocs de frequences contigues (guess, 2026-10-01)
	integer :: nblk, ib, jb, ig, ifin, nval, ch_a, ch_z, ch_m, lmax
	integer, allocatable :: b_deb(:), b_fin(:), ordre_blk(:)
	! ordonnancement guide par la RAM (2026-10-02)
	real(prec), allocatable :: mth(:), r_blk(:)
	integer, allocatable :: tri_ram(:)
	logical, allocatable :: pris(:)
	real(prec) :: reserve, budget, rss0, dispo, det(4), rmin_t, rmax_t
	integer :: actifs, i3, imin, imax
	logical :: fini
	logical, allocatable :: fait(:), est_tete(:)
	real(prec), parameter :: Go = 1024.0_prec**3
	complex(kind=8), allocatable :: xg(:), s1(:), s2(:), surf(:)
	complex(kind=8) :: x0_vide(1), surf_vide(1)

	integer :: nm
	double precision :: sj(0:80),dj(0:80),x,x0(3),jn(0:80,3)	
	character(len=1) :: answer

!    Nombre de CPU disponibles, export OMP_NUM_THREADS=36
	num_core=omp_get_max_threads()
!
! ------------------------------------------------------------------------------
! lecture fichier de config
! ------------------------------------------------------------------------------
!
	! test si il y a un fichier en argument
	! --- OpenBLAS (2026-09-28) : le parallelisme est sur les FREQUENCES (OpenMP).
	!     Sans OPENBLAS_NUM_THREADS=1, OpenBLAS (pthread) ajoute son propre pool
	!     de threads -> sursouscription (tuyau : 47 threads au lieu de 24, -29 %
	!     de temps une fois corrige). La variable est lue au chargement de la
	!     bibliotheque : on ne peut que prevenir.
	call get_environment_variable('OPENBLAS_NUM_THREADS', env_ob, status=st_ob)
	ob_ok = (st_ob == 0 .and. trim(env_ob) == '1')
	if (.not. ob_ok) then
		write(*,*) 'WARNING: OPENBLAS_NUM_THREADS is not 1 (value: "'//trim(env_ob)//'")'
		write(*,*) '         run with:  export OPENBLAS_NUM_THREADS=1'
	endif
	if (command_argument_count() .eq. 1) then
		call get_command_argument(number=1, length=arglen)  ! Assume for simplicity success
		allocate (character(arglen) :: filename)
		call get_command_argument(number=1, value=filename, status=stat)
	else
		allocate(character(10) :: filename)
		filename='config.txt'
	endif
	fileunit=10
!	filename='config_G.txt'
	open(unit=fileunit,file=filename)

	read(fileunit,*) varBid ! Nombre de Fréquences
	read(fileunit,*) NF
	read(fileunit,*) varBid ! Fréquence d'échantillonnage (Hz)
	read(fileunit,*) fe
	read(fileunit,*) varBid ! Fréquence initiale
	read(fileunit,'(A)') ligne_freq      ! une valeur, ou la LISTE des frequences (2026-09-30)
	call lire_liste_reels(ligne_freq, liste_freq, nliste)
	if (nliste < 1) then
		write(*,*) 'Case file: empty frequency line (value 3)'
		stop 1
	endif
	inF = liste_freq(1)
	if (nliste > 1) then                  ! liste explicite : NF = nombre de valeurs (NF = 0 accepte)
		if (NF == 0) NF = nliste
		if (NF /= nliste) then
			write(*,*) 'Case file: NF =', NF, ' but', nliste, ' frequencies in the list'
			stop 1
		endif
	endif
	read(fileunit,*) varBid ! Fichier du maillage GMSH, .msh version 4 ascii
	read(fileunit,'(A)') file_msh ; call nom_fichier(file_msh)   ! ligne entiere : chemins avec / acceptes (2026-11-20)
	read(fileunit,*) varBid ! Ficher contenant les points extérieurs
	read(fileunit,'(A)') file_Pext ; call nom_fichier(file_Pext)   ! ligne entiere : chemins avec / acceptes (2026-11-20)
	read(fileunit,*) varBid ! Fichier de sauvegarde
	read(fileunit,'(A)') file_save ; call nom_fichier(file_save)   ! ligne entiere : chemins avec / acceptes (2026-11-20)
	read(fileunit,*) varBid ! Position du point source dans l'oreille > à 2mm de la surface
	read(fileunit,*) PT0    ! x,y,z
	read(fileunit,*) varBid ! Normalisation de la HRTF par pression au centre, oui(1), non(0)
	read(fileunit,*) HRTF_norm
	! --- 2026-11-20 : MEME ORDRE QUE main_HRTF_Ray : verbose, log, c, fichier technique, puis
	!     impedance de peau FACULTATIVE en derniere ligne (absente : rigide). Ancien ordre (avant le
	!     2026-11-20 : impedance en 9e valeur, avant verbose) reconnu : la 9e valeur est alors deux
	!     nombres ou un nom de fichier, et non un entier seul.
	read(fileunit,*) varBid
	read(fileunit,'(A)') ligne9
	ancien_ordre = .not. entier_seul(ligne9)
	if (ancien_ordre) then
		ligne_refl = ligne9                  ! « Rmax_BF Rmax_HF » ou fichier
		read(fileunit,*) varBid ! Verbose
		read(fileunit,*) verbose
	else
		read(ligne9,*) verbose
	endif
	read(fileunit,*) varBid ! Fichier log (stdout si vide)
	read(fileunit,'(A)') file_log ; call nom_fichier(file_log)   ! ligne entiere : chemins avec / acceptes (2026-11-20)
	! --- celerite du son (2026-09-28) : OBLIGATOIRE. Elle n'est plus prise dans parametres.f08
	!     (c y sert au tuyau) : une valeur par defaut changerait les HRTF sans prevenir.
	read(fileunit,*,iostat=ios_c) varBid
	if (ios_c == 0) read(fileunit,*,iostat=ios_c) c_son
	if (ios_c /= 0 .or. c_son <= 0.0_prec) then
		write(*,*) 'Case file: value 11 must be the speed of sound (m/s), for example:'
		write(*,*) '  ! Speed of sound c (m/s)'
		write(*,*) '  340.0'
		stop 1
	endif
	! --- fichier des parametres TECHNIQUES MLFMM (2026-09-30), format « cle = valeur » : voir
	!     parametres.f08, lire_parametres_MLFMM. Ligne absente, vide ou '-' : valeurs par defaut
	!     (reglage rapide du preprint depuis le 2026-11-20).
	file_param = ''
	read(fileunit,*,iostat=ios_c) varBid
	if (ios_c == 0) read(fileunit,'(A)',iostat=ios_c) file_param
	if (ios_c /= 0) file_param = ''
	file_param = adjustl(file_param)
	if (.not. ancien_ordre) then
		! impedance de peau, ligne FACULTATIVE : « Rmax_BF Rmax_HF » ou fichier ; absente : rigide
		ligne_refl = '1.0 1.0'
		read(fileunit,*,iostat=ios_c) varBid
		if (ios_c == 0) read(fileunit,'(A)',iostat=ios_c) ligne_refl
		if (ios_c /= 0 .or. len_trim(ligne_refl) == 0) ligne_refl = '1.0 1.0'
	endif
	close(fileunit)
	call lire_parametres_MLFMM(trim(file_param), NL, seuil_adapt, dN_BF_lu, f_boost)
	fusion_Vc = (liberation_V1 == 1)   ! Diff : V1 fusionne dans V2 (Vc), module_MLFMM, 2026-11-06

	! Ouverture fichier log (tous niveaux de verbose ; nom vide : ecran)
	if (len_trim(file_log) > 0) then
		log_unit = 20
		open(unit=log_unit, file=trim(file_log), status='replace', action='write')
	else
		log_unit = 6
	end if
	call date_and_time(values=dt)
	write(log_unit,'(A,I4.4,"-",I2.2,"-",I2.2,1X,I2.2,":",I2.2)') ' MLFMM2HRTF  main_HRTF_Diff  ', dt(1), dt(2), dt(3), dt(5), dt(6)
	write(log_unit,'(A)') ''
	call log_cle('case file', filename)
	if (ancien_ordre) call log_cle('', 'legacy order (skin impedance as value 9, before the verbose): accepted')
	if (ob_ok) then
		call log_cle('threads', txt_i(num_core)//' (OpenMP)')
	else
		call log_cle('threads', txt_i(num_core)//' (OpenMP)  WARNING: OPENBLAS_NUM_THREADS is not 1')
	endif
	write(*,'(A)') ' main_HRTF_Diff: case file '//filename//', log '//trim(file_log)

!
! ------------------------------------------------------------------------------
!
! Paramètres
!
! ------------------------------------------------------------------------------
!

	! init frequence
	t1=omp_get_wtime()
	allocate(freq(NF))
	allocate(k(NF))
	allocate(Ared(NF))
	allocate(Refl(NF))
	if (nliste > 1) then
		freq = liste_freq            ! liste explicite de la config
	elseif (NF == 1) then
		freq = [inF]
	else
		freq = [(i1*fe/2.0d0/NF,i1=1,NF)]
	end if
	k=2.0d0*pi*freq/c_son      ! celerite lue dans la config (pas parametres%c)
	if (nliste > 1) then
		call log_cle('frequencies', txt_i(NF)//' (list), '//txt_r(freq(1),'(F12.3)')//' to '//txt_r(freq(NF),'(F12.3)')//' Hz')
	else
		call log_cle('frequencies', txt_i(NF)//' (regular), '//txt_r(freq(1),'(F12.3)')//' to '//txt_r(freq(NF),'(F12.3)')//' Hz')
	endif
	call log_cle('c', txt_r(c_son,'(F10.2)')//' m/s')
	if (verbose >= 2) write(log_unit,*) 'k(1) =', k(1), ' k(NF) =', k(NF)
	! magnitude du coefficient de reflexion, lineaire en frequence entre deux
	! bornes FIXES (2026-09-30) : Rmax_BF (paroi rigide) a 0 Hz et Rmax_HF (~0,9)
	! a 24 kHz, quels que soient fe et les frequences calculees (avant :
	! freq(NF), qui dependait de la derniere frequence calculee).
	call admittance_paroi(ligne_refl, freq, NF, Ared)   ! module parametres (commun avec main_HRTF_Ray, 2026-11-06)
	if (any(Ared /= dcmplx(0.0d0,0.0d0))) then
		call log_cle('skin', ligne_refl)
	else
		call log_cle('skin', 'rigid')
	endif
	call log_cle('source', 'point ('//txt_r(PT0(1),'(F10.5)')//', '//txt_r(PT0(2),'(F10.5)')//', '// &
		txt_r(PT0(3),'(F10.5)')//') m (monopole)')



	! Nlevel dans parametres.mod
	Ni=(2**(NLF+1)-1)**3-27 ! nombre d'intéractions dernier niveau quand NLF != 2 (pour stockage)
	N=NL-NLF+1  ! nombre de niveaux pour allocations

!test SPHJ
!	x0(1)=100.0d0
!	x0(2)=2000.0d0
!	x0(3)=3500.0d0
!	x=100.0d0
!	call SPHJ_test(80,x,nm,sj,dj)
!	write(*,*) 'sj1=',sj(80)
!	x=2000.0d0
!	call SPHJ_test(80,x,nm,sj,dj)
!	write(*,*) 'sj2=',sj(80)
!	x=3500.0d0
!	call SPHJ_test(80,x,nm,sj,dj)
!	write(*,*) 'sj3=',sj(80)
!	call sph_j(80,3,x0,jn)
!	write(*,*) 'sph_j=',jn(0,:)
!	write(*,*) 'sph_j=',jn(80,:)
!
! ------------------------------------------------------------------------------
!
! Lecture maillage
! GMSH, .msh, ascii, version 4
!
! ------------------------------------------------------------------------------
!
!	call lecture_gmsh_v4(trim(file_msh),len(trim(file_msh)),3)   ! pour des triangles T3

!
! ------------------------------------------------------------------------------
!
! Lecture maillage
! GMSH, .msh, ascii, version 2
!
! ------------------------------------------------------------------------------
!
	call lecture_maillage_hrtf(trim(file_msh))	  ! T3 : gmsh 2.2 ASCII (.msh) ou STL (.stl, binaire ou ASCII)
!
! ------------------------------------------------------------------------------
!
! pré-calculs BEM (base tangentiel, points centraux, normales, jacobien, init matrice creuse Anear)
!
! ------------------------------------------------------------------------------
!
! à partir du maillage, base tangentielle, Pc, norc et Jac0
	call pre_BEM_T0()
	call aretes_extremes(hmin_m, hmax_m)
	call log_cle('mesh', trim(file_msh)//': '//txt_i(NE)//' triangles, '//txt_i(NP)//' nodes, edges '// &
		txt_r(hmin_m*1d3,'(F10.2)')//' to '//txt_r(hmax_m*1d3,'(F10.2)')//' mm')
	k_chi_max = maxval(k)        ! CHI : profondeur limitee selon la longueur d onde a f max (2026-10-29)
	if (verbose >= 2) write(log_unit,*) 'interp_solveur =', interp_solveur, ' (1: HF translations interpolate then shift)'
	call init_points_cibles()     ! points cibles des lignes : Pc, ou interieurs si chi > 0 (essai CHI)
!
! ------------------------------------------------------------------------------
!
! Création Octree et init Anear pour BEM
!
! ------------------------------------------------------------------------------
!
	if (len_trim(file_param) == 0 .or. trim(file_param) == '-') then
		call log_cle('technical', '- (none: default settings)')
	else
		call log_cle('technical', file_param)
	endif
	call log_cle('setting', merge('precise', 'fast   ', nom_reglage == 'precis'))
	if (verbose >= 2) write(log_unit,'(A,I3,A,F6.2,A,I2,A,I2,A,I2)') ' NL =', NL, '  seuil_adapt =', seuil_adapt, &
		'  FMM_meth =', FMM_meth, '  FMM_pre =', FMM_pre, '  propagation =', methode_propagation
	! boost d'ordre BF (parametres : dN_BF, k_boost) : defaut +4 en adaptatif
	if (dN_BF_lu >= 0) then
		dN_BF = dN_BF_lu
	else
		dN_BF = merge(4, 0, NL < 0)
	endif
	k_boost = 2.0d0*pi*f_boost/c_son
	! CHIEF (cle interne, gardee pour les figures du preprint) : remplace Burton/Miller entre
	! f_min_chief et f_max_chief (ex-cle f_sans_BM, 2026-11-20)
	call lire_points_chief()
	if (f_max_chief > 0.0d0 .and. NC_chief == 0) stop 'Technical file: f_max_chief without CHIEF points (key chief)'
	k_max_chief = 2.0d0*pi*f_max_chief/c_son
	k_min_chief = 2.0d0*pi*f_min_chief/c_son
	if (chi_d > 0.0d0 .and. verbose >= 2) write(log_unit,'(A,F6.3,A)') ' CHI (Cunefare 1989): target points', chi_d, &
		' x element size below the surface'
	call log_cles_internes(f_boost, seuil_adapt, NL)
	call bandes_en_k(c_son)          ! limites des bandes epsilon0/epsilon2 en nombre d'onde
	if (nb_bandes > 0 .and. verbose >= 2) then
		write(log_unit,'(A,10F9.0)') ' bands (Hz) :', f_bandes(1:nb_bandes)
		write(log_unit,'(A,11ES9.1)') '   epsilon0 :', eps0_bandes(1:nb_bandes+1)
		write(log_unit,'(A,11ES9.1)') '   epsilon2 :', eps2_bandes(1:nb_bandes+1)
		write(log_unit,'(A,11F9.2)') '   seuil_HF :', seuilHF_bandes(1:nb_bandes+1)
	endif
	if (verbose >= 2) then
		write(log_unit,'(A,I2,A,F8.1,A,F5.2,A,I3)') ' LF boost: dN_BF =', dN_BF, '  below', f_boost, &
		     ' Hz   seuil_HF =', seuil_HF, '   Nmin =', Nmin_ord
		write(log_unit,'(A,I4,A,ES9.2,A,ES9.2,A,ES9.2)') ' ILUT lfil =', lfil_ilut, '  epsilon1 =', epsilon1, &
			'  epsilon0 =', epsilon0, '  epsilon2 =', epsilon2
	endif
	if (NL < 0) then                ! ligne NL a -1 : octree ADAPTATIF
		call NL_adaptatif_T0(NL)
		fmm_adaptatif = .true.
		N=NL-NLF+1
	endif
	allocate(octree1(N))
	if (fmm_adaptatif .and. chi_d > 0.0d0) then
		! CHI : octree sur 2 NE points, barycentres (sources) puis points interieurs (cibles),
		! chaque cible dans SA cellule (pas d'evaluation hors domaine de convergence)
		cibles_sep = .true.
		allocate(Ptous(2*NE,3))
		Ptous(1:NE,:) = Pc
		Ptous(NE+1:2*NE,:) = Pt
		call make_octree_all(2*NE,Ptous,octree1,NL,NLF)
		deallocate(Ptous)
		if (verbose >= 2) write(log_unit,'(A)') ' CHI: octree on 2 NE points (sources Pc, interior targets Pt)'
	else
		call make_octree_all(NE,Pc,octree1,NL,NLF)
	endif
	if (verbose >= 2) write(log_unit,*) 'cells per level:', octree1(:)%NL
	if (.not. fmm_adaptatif) call log_cle('octree', 'uniform, '//txt_i(NL)//' levels')
! à partir de l'octree, initialisation matrice creuse Anear pour BEM
	if (fmm_adaptatif) then
		! chaque element a SON niveau (cellule >= seuil_adapt x sa plus grande arete),
		! barycentres (cibles) evalues dans leur premiere cellule terminale
		call completer_octree(octree1, N, merge(2*NE, NE, cibles_sep))
		allocate(niv_E(NE))
		call niveaux_elements(octree1, N, NE, h_el, seuil_adapt, niv_E, nhors)
		if (nhors > 0) write(log_unit,'(A,I8,A)') ' WARNING:', nhors, &
			' elements larger than the cell of level NLF'
		call marquer_sources_cibles(octree1, N, NE, niv_E, .true.)
		if (verbose >= 2) then
			write(log_unit,'(A)') ' adaptive octree: level, cell (mm), cells, elements, leaf cells'
			do ilv = 1, N
				write(log_unit,'(I8,F10.3,3I10)') ilv+NLF-1, octree1(ilv)%cell(1)%length*1d3, octree1(ilv)%NL, &
					count(niv_E == ilv), count(octree1(ilv)%cell(:)%feuille)
			end do
		endif
		call init_Anear_T0_adapt(octree1, N)   ! paires au niveau de chaque element (iA, jA)
	else
		call init_Anear(num_K0,octree1,N)
	endif
	if (verbose >= 2) write(log_unit,*) 'near-field pairs num_K0 =', num_K0
!
! ------------------------------------------------------------------------------
!
! Calculs des champs et gradients pour méthode par diffraction
!
! ------------------------------------------------------------------------------
!
! calcul du champ de gradient en inc0 (oreille) pour calcul du gradient sur la tête (cf. Gumerov)
	V=1.0d0
	allocate(dp1(1:NE,1:NF))
	allocate(p1(1:NE,1:NF))
	call champ_monopolaire_derive(dp1,k,NF,Pt,NE,PT0,V,norc)   ! aux points cibles (= Pc sauf CHI)
	call champ_monopolaire(p1,k,NF,Pt,NE,PT0,V)
	if (NC_chief > 0) then      ! CHIEF : champ incident aux points interieurs (second membre -p_inc)
		allocate(pinc_chief(NC_chief,NF))
		call champ_monopolaire(pinc_chief,k,NF,Q_chief,NC_chief,PT0,V)
	endif

! Points de réception (ou source pour HRTF par réciprocité)
! lecture fichier maillage
	fileunit=10
	open(unit=fileunit,file=file_Pext)
	read(fileunit,*) Next
	allocate(Pext(1:Next,1:3))
	do i1=1,Next
		!write(*,*) ind
		read(fileunit,*) Pext(i1,1), Pext(i1,2), Pext(i1,3)
	end do
	close(fileunit)
	call log_cle('directions', trim(file_Pext)//': '//txt_i(Next)//' points')
	call log_cle('normalisation', merge('yes (free field at the head centre)', 'no (raw pressure)                  ', &
		HRTF_norm == 1))

! champs au centre pour normalisation (division, cf Gumerov)
	r(1,1)=0.0d0
	r(1,2)=0.0d0
	r(1,3)=0.0d0
	allocate(pref0(1:Next,1:NF))
	if (HRTF_norm==1) then
		call champ_monopolaire(pref0,k,NF,Pext,Next,r,V)
	else
		pref0=1.0d0
	endif
! champs en Pext pour calcul diffraction
	allocate(pref(Next,NF))
	call champ_monopolaire(pref,k,NF,Pext,Next,PT0,V)
!	pref=dcmplx(0.0d0,0.0d0)
!	write(*,*) 'in main, pref=',pref(1:3,270)
!
!-------------------------------------------------------------
!
! Initialisation BF_MLFMM : rotation et distance pour translation
! non dépendent de la fréquence
!
!-------------------------------------------------------------
!
!	t01=omp_get_wtime()
!	allocate(SPH(N))
	call init_BF_MLFMM_v1(octree1,N,SPH0,19)
!	t02=omp_get_wtime()
!	write(*,*) 'fin Init_BF_MLFMM en ',t02-t01,' sec'

!
!-----------------------------------------------------------------------------------------------
!
! allocations avant boucle
!
!-----------------------------------------------------------------------------------------------
!

	allocate(sol0(Next,NF))
	allocate(sol2(Next))
!	allocate(sol1(Next,NF))
!	allocate(sol3(Next))
	allocate(iter0(NF))
!	allocate(ordre0(N,NF))
	allocate(ordre1(N))
	sol2=dcmplx(0.0d0,0.0d0)
	sol0=dcmplx(0.0d0,0.0d0)
	
	t2=omp_get_wtime()
	call log_cle('preparation', txt_duree(t2-t1))

    ! pour fréquence
!	fileunit=10
!	filename='config_freq_D.txt'
!	open(unit=fileunit,file=filename)
!	read(fileunit,*) i3
!	close(fileunit)

	t1=omp_get_wtime()
	
	!deb=(i3-1)*num_core+1
	!fin=i3*num_core
	!deb=270
	!fin=270
	deb=1
	fin=NF
	! Modification de l'ordre pour optimiser la mémoire
	allocate(opt_mem(NF))
	do i1=1,NF/2
		opt_mem(2*i1-1)=i1
		opt_mem(2*i1)  =NF-(i1-1)
	enddo
	! cas impair : élément du milieu
	if (mod(NF,2) == 1) opt_mem(NF) = (NF+1)/2


	if (guess_mode == 0) then
	call log_entete_frequences(NF)
	!$OMP PARALLEL DO DEFAULT(PRIVATE) shared(octree1,E,Ps) &
	!$OMP SHARED(Pc,norc,Jac0,base_a1,base_a2,num_K0,sol0,Pext,Next,iter0) &
	!$OMP SHARED(freq,k,p1,dp1,pref,pref0,N,coef_ilut,SPH0) &
	!$OMP SHARED(epsilon0,epsilon1,FMM_pre,FMM_meth,Ared,opt_mem,NF,x0_vide,surf_vide,log_unit,verbose) &
	!$OMP SCHEDULE(DYNAMIC,1)
!	!$OMP SCHEDULE(GUIDED,1)
	do i1=1,NF
!	do i1=270,270
		if (verbose >= 2) write(log_unit,'(A,I5,A,F12.3,A)') ' #', opt_mem(i1), ' start: f =', freq(opt_mem(i1)), ' Hz'
		t_f = omp_get_wtime()
		call fait_tout_HRTF_DIFF_T0(opt_mem(i1),N,octree1,SPH0,k(opt_mem(i1)),p1(:,opt_mem(i1)),dp1(:,opt_mem(i1)), &
			pref(:,opt_mem(i1)),pref0(:,opt_mem(i1)),Ared(opt_mem(i1)),iter,sol2,-1,x0_vide,surf_vide)
		call log_frequence(opt_mem(i1), freq(opt_mem(i1)), iter, omp_get_wtime()-t_f, NF)
!		call fait_tout_HRTF_T0(i1,N,octree1,k(i1),dp1(:,i1),pref0(:,i1),pref(:,i1),iter,sol2)
!		call fait_tout_HRTF_Test(i1,N,octree1,k(i1),dp1(:,i1),pref0(:,i1),pref(:,i1),iter,sol2,sol3)
!		call fait_tout_HRTF_Test(i1,N,octree1,k(i1),dp1(:,i1),pref0(:,i1),pref(:,i1),iter,sol2)
!		call fait_tout_FMM_Test(i1,N,octree1,k(i1),dp1(:,i1),pref0(:,i1),pref(:,i1),iter,ordre1)
!		write(*,*) 'sol2=',sol2(1:3),' en ',iter,' itérations'
		sol0(:,opt_mem(i1))=sol2
!		sol1(:,i1)=sol3
!		ordre0(:,i1)=ordre1
		iter0(opt_mem(i1))=iter
		
	end do
	!$OMP END PARALLEL DO
	else
	! ESTIMATION INITIALE GMRES (guess, 2026-10-01) : blocs de frequences CONTIGUES,
	! parcourus en montant par un meme thread ; depart a zero en tete de bloc, puis
	! solution precedente (guess 1) ou extrapolation lineaire en k (guess 2). Ordre
	! des blocs : voir l'ordonnancement guide par la RAM ci-dessous.
	nblk = max(1, min(NF, blocs_par_thread*num_core))
	allocate(b_deb(nblk), b_fin(nblk), ordre_blk(nblk))
	do ib = 1, nblk
		b_deb(ib) = (ib-1)*NF/nblk + 1
		b_fin(ib) = ib*NF/nblk
	end do
	! ORDONNANCEMENT GUIDE PAR LA RAM (2026-10-02) : RAM estimee par frequence
	! (memoire_thread_T0), bloc = son maximum ; blocs pris les plus gourmands
	! d'abord (HF), un bloc ne demarre que si la somme des blocs en cours tient
	! dans le budget (ram_max, sinon 85 % de la RAM disponible), sauf s'il n'y a
	! aucun bloc en cours. OMP_NUM_THREADS = nombre MAXIMAL de threads.
	! Ordre de depart : entrelace gros / petit (voir ENTRELACEMENT).
	allocate(mth(NF), r_blk(nblk), pris(nblk))
	do i1 = 1, NF
		mth(i1) = memoire_thread_T0(octree1, N, k(i1), FMM_meth)
	end do
	do ib = 1, nblk
		r_blk(ib) = maxval(mth(b_deb(ib):b_fin(ib)))
		ordre_blk(ib) = ib
	end do
	do ib = 1, nblk-1                    ! tri decroissant de la RAM (selection, nblk petit)
		i2 = ib
		do jb = ib+1, nblk
			if (r_blk(ordre_blk(jb)) > r_blk(ordre_blk(i2))) i2 = jb
		end do
		i3 = ordre_blk(ib) ; ordre_blk(ib) = ordre_blk(i2) ; ordre_blk(i2) = i3
	end do
	! ENTRELACEMENT (2026-10-02) : alterner le plus gourmand et le moins gourmand des
	! blocs restants (HF, BF, HF, BF...) au lieu de tous les HF d'abord : moins de blocs
	! HF simultanes (GMRES borne par la bande passante memoire), RAM moyenne plus faible.
	! Le guess n'est pas touche (chaque bloc reste calcule en entier par un thread).
	allocate(tri_ram(nblk)) ; tri_ram = ordre_blk
	i2 = 1 ; i3 = nblk
	do ib = 1, nblk
		if (mod(ib,2) == 1) then
			ordre_blk(ib) = tri_ram(i2) ; i2 = i2 + 1
		else
			ordre_blk(ib) = tri_ram(i3) ; i3 = i3 - 1
		endif
	end do
	deallocate(tri_ram)
	rss0  = lit_kb('/proc/self/status', 'VmRSS:')*1024.0_prec
	dispo = lit_kb('/proc/meminfo', 'MemAvailable:')*1024.0_prec
	if (ram_max > 0.0_prec) then
		budget = ram_max*Go - rss0
	else
		budget = 0.85_prec*dispo
	endif
	imin = minloc(mth, 1) ; imax = maxloc(mth, 1)
	rmin_t = mth(imin) ; rmax_t = mth(imax)
	call log_cle('RAM/thread', txt_r(rmin_t/Go,'(F10.2)')//' GB ('//txt_r(freq(imin),'(F10.1)')//' Hz) to '// &
		txt_r(rmax_t/Go,'(F10.2)')//' GB ('//txt_r(freq(imax),'(F10.1)')//' Hz)')
	det = 0.0_prec
	rmax_t = memoire_thread_T0(octree1, N, k(imax), FMM_meth, det)
	if (verbose >= 2) write(log_unit,'(A,4F7.2,A)') '   at the maximum: FMM, near field, ILUT, GMRES =', det/Go, ' GB'
	call log_cle('memory', 'shared '//txt_r(rss0/Go,'(F10.2)')//' GB, available '//txt_r(dispo/Go,'(F10.1)')// &
		' GB, budget of the threads '//txt_r(budget/Go,'(F10.1)')//' GB')
	call log_cle('parallel', txt_i(num_core)//' threads at most: '// &
		txt_i(min(num_core, max(1, int(budget/rmax_t))))//' at the largest frequency, '// &
		txt_i(min(num_core, max(1, int(budget/rmin_t))))//' at the smallest')
	if (rmax_t > budget) write(log_unit,'(A)') ' WARNING: one frequency alone exceeds the RAM budget; it will run alone'
	if (estimer == 1) then             ! estimation seule : table par frequence (log et ecran), puis arret
		write(*,'(A)') '  freq (Hz)   RAM/thread (GB)   FMM  near field  ILUT  GMRES (GB)'
		do i1 = 1, NF
			mth(i1) = memoire_thread_T0(octree1, N, k(i1), FMM_meth, det)
			write(*,'(F11.1,F12.2,6X,4F8.2)') freq(i1), mth(i1)/Go, det/Go
			write(log_unit,'(A,F11.1,F9.2)') ' estimate: freq (Hz), GB/thread ', freq(i1), mth(i1)/Go
		end do
		write(*,'(A,F7.2,A,F7.1,A)') ' shared RSS ', rss0/Go, ' GB; budget of the threads ', budget/Go, ' GB'
		write(*,'(A)') ' estimer = 1: stop without computing'
		stop
	endif
	if (verbose >= 2) write(log_unit,'(A,I2,A,I5,A)') ' guess =', guess_mode, ':', nblk, ' blocks of contiguous frequencies'
	call log_entete_frequences(NF)
	pris = .false. ; reserve = 0.0_prec ; actifs = 0
	! tetes de chaine (2026-11-04) : num_core blocs regulierement espaces (1, 1+nblk/num_core, ...)
	allocate(fait(nblk), est_tete(nblk)) ; fait = .false. ; est_tete = .false.
	do jb = 1, min(num_core, nblk)
		est_tete(1 + (jb-1)*nblk/min(num_core, nblk)) = .true.
	end do
	!$OMP PARALLEL DEFAULT(PRIVATE) shared(octree1,E,Ps) &
	!$OMP SHARED(Pc,norc,Jac0,base_a1,base_a2,num_K0,sol0,Pext,Next,iter0) &
	!$OMP SHARED(freq,k,p1,dp1,pref0,pref,Ared,N,coef_ilut,SPH0,NE,log_unit,guess_mode) &
	!$OMP SHARED(epsilon0,epsilon1,FMM_pre,FMM_meth,nblk,b_deb,b_fin,ordre_blk) &
	!$OMP SHARED(r_blk,pris,reserve,budget,actifs,fait,est_tete,t1,verbose,NF)
	allocate(xg(NE), s1(NE), s2(NE), surf(NE))
	! CHAINAGE (2026-11-04) : chaque thread parcourt une CHAINE de blocs contigus en montant, sans
	! repartir de zero d'un bloc a l'autre (extrapolation continue). Choix du bloc suivant :
	!  1. le bloc ib+1 du bloc qu'il vient de finir, s'il est libre et tient dans le budget ;
	!  2. sinon une tete de chaine libre (num_core tetes espacees), dans l'ordre entrelace par la RAM ;
	!  3. sinon (fin de calcul) il coupe en deux la plus longue suite de blocs libres : il en prend la
	!     seconde moitie si la chaine qui la precede tourne encore, le debut sinon.
	! Departs de zero : ~1 par thread + les coupures de fin, au lieu de 1 par bloc.
	! nval : solutions precedentes contigues valides (0, 1, 2) ; ifin : dernier bloc fini par le thread.
	ib = 0 ; ifin = 0 ; nval = 0
	do
		!$OMP CRITICAL (ordonnancement)
		if (ib > 0) then                   ! fin du bloc ib : liberer sa reservation
			reserve = reserve - r_blk(ib) ; actifs = actifs - 1 ; ifin = ib ; fait(ib) = .true.
			if (verbose >= 2) write(log_unit,'(A,I5,A,F9.1,A,I3,A,I3)') ' end of block', ib, ' at t =', omp_get_wtime()-t1, &
			     ' s, thread', omp_get_thread_num(), ', blocks running', actifs
		endif
		ib = 0
		fini = all(pris)
		if (ifin > 0 .and. ifin < nblk) then                                   ! 1. suite de la chaine
			if (.not. pris(ifin+1) .and. (reserve + r_blk(ifin+1) <= budget .or. actifs == 0)) ib = ifin + 1
		endif
		if (ib == 0) then                                                      ! 2. tete de chaine
			do jb = 1, nblk
				if (pris(ordre_blk(jb)) .or. .not. est_tete(ordre_blk(jb))) cycle
				if (reserve + r_blk(ordre_blk(jb)) <= budget .or. actifs == 0) then
					ib = ordre_blk(jb)
					exit
				endif
			end do
		endif
		if (ib == 0) then                                                      ! 3. coupure de la plus longue suite libre
			lmax = 0 ; jb = 1
			do while (jb <= nblk)
				if (pris(jb)) then
					jb = jb + 1 ; cycle
				endif
				ch_a = jb
				do while (jb <= nblk)
					if (pris(jb)) exit
					jb = jb + 1
				end do
				ch_z = jb - 1 ; ch_m = ch_a
				if (ch_a > 1) then
					if (.not. fait(ch_a-1)) ch_m = ch_a + (ch_z-ch_a+1)/2      ! la chaine precedente tourne : seconde moitie
				endif
				if (ch_z-ch_a+1 > lmax .and. (reserve + r_blk(ch_m) <= budget .or. actifs == 0)) then
					lmax = ch_z - ch_a + 1 ; ib = ch_m
				endif
			end do
		endif
		if (ib > 0) then
			pris(ib) = .true. ; reserve = reserve + r_blk(ib) ; actifs = actifs + 1
			if (ifin == 0 .or. ib /= ifin + 1) nval = 0
			if (verbose >= 2) write(log_unit,'(A,I5,A,F9.1,A,F9.1,A,F6.2,A,I3,A,F7.2,A,A,A,F9.1,A,I3)') ' block', ib, &
			     ' (', freq(b_deb(ib)), ' -', freq(b_fin(ib)), ' Hz,', r_blk(ib)/Go, ' GB):', actifs, ' blocks running, ', &
			     reserve/Go, ' GB reserved', merge(' (chain)', '        ', nval > 0), ' t =', omp_get_wtime()-t1, &
			     ' s, thread', omp_get_thread_num()
		endif
		!$OMP END CRITICAL (ordonnancement)
		if (ib == 0) then
			if (fini) exit
			call sleep(1)              ! plus de place : attendre qu'un bloc se termine
			cycle
		endif
		do i1 = b_deb(ib), b_fin(ib)
			if (nval == 0) then
				ig = 0 ; xg = dcmplx(0.0d0,0.0d0)
			else if (guess_mode == 1 .or. nval == 1) then
				ig = 1 ; xg = s1
			else
				ig = 1 ; xg = s1 + (s1-s2)*((k(i1)-k(i1-1))/(k(i1-1)-k(i1-2)))
			endif
			if (verbose >= 2) write(log_unit,'(A,I5,A,F12.3,A,F10.4,A,I2)') ' #', i1, ' start: f =', freq(i1), &
				' Hz, k =', k(i1), ', initial guess', ig
			t_f = omp_get_wtime()
			call fait_tout_HRTF_DIFF_T0(i1,N,octree1,SPH0,k(i1),p1(:,i1),dp1(:,i1),pref(:,i1),pref0(:,i1), &
				Ared(i1),iter,sol2,ig,xg,surf)
			call log_frequence(i1, freq(i1), iter, omp_get_wtime()-t_f, NF)
			sol0(:,i1)=sol2
			iter0(i1)=iter
			s2 = s1 ; s1 = surf ; nval = min(nval + 1, 2)
		end do
	end do
	deallocate(xg, s1, s2, surf)
	!$OMP END PARALLEL
	deallocate(mth, r_blk, pris, fait, est_tete)
	endif
	t2=omp_get_wtime()
	if (verbose >= 2) write(log_unit,*) 'frequency loop:', t2-t1, ' s'
!
!-----------------------------------------------------------------------------------------------
!
! ecriture resultats
!
!-----------------------------------------------------------------------------------------------
!
	filename=trim(file_save)//'.txt'
	fileunit=10
	open(unit=fileunit,file=filename)
	do i1=deb,fin
!	    write(fileunit,*) itera(i1)
	    do i2=1,Next
!	    	write(fileunit,*) dreal(sol0(i2,i1)), dimag(sol0(i2,i1)), iter0(i1), t2-t1
!	    	write(fileunit,*) dreal(sol0(i2,i1)), dimag(sol0(i2,i1)),dreal(sol1(i2,i1)), dimag(sol1(i2,i1))
	    	write(fileunit,*) dreal(sol0(i2,i1)), dimag(sol0(i2,i1)), iter0(i1), t2-t1
!	    	write(fileunit,*) ordre0(:,i1)
	    enddo
	enddo
	close(fileunit)
	call log_bilan(NF, freq, iter0, t2-t1, trim(file_save)//'.txt')
	write(*,'(A)') ' main_HRTF_Diff: done, '//txt_i(NF)//' frequencies in '//txt_duree(t2-t1)//', results '// &
		trim(file_save)//'.txt'

    ! pour fréquence
!	fileunit=10
!	filename='config_freq_D.txt'
!	open(unit=fileunit,file=filename)
!	write(fileunit,*) i3+1
!	close(fileunit)
	deallocate(sol0)
	deallocate(sol2)
	deallocate(E)
	deallocate(Pc)
	deallocate(norc)
	deallocate(Ps)
	deallocate(octree1)
	deallocate(p1)
	deallocate(dp1)
	deallocate(Pext)
	deallocate(pref)
	deallocate(pref0)
	deallocate(iter0)
	if (log_unit /= 6) close(log_unit)

	contains

	! nom de fichier lu sur une ligne entiere (format '(A)') : la lecture en format libre (*)
	! s'arretait au premier / (fin d'enregistrement en list-directed) et aux espaces.
	! Tabulations -> espaces, espaces en tete retires, guillemets englobants retires.
	! vrai si la ligne ne contient qu'un entier (le verbose), comme dans l'ordre de main_HRTF_Ray
	logical function entier_seul(l)
	character(len=*), intent(in) :: l
	character(len=len(l)) :: t
	integer :: v, ios, k
	t = adjustl(l)
	k = scan(t, '!#')
	if (k > 0) t(k:) = ' '
	entier_seul = .false.
	if (len_trim(t) == 0 .or. index(trim(t), ' ') > 0 .or. scan(trim(t), '.eEdD/') > 0) return
	read(t,*,iostat=ios) v
	entier_seul = (ios == 0)
	end function entier_seul

	subroutine nom_fichier(nom)
	character(len=*), intent(inout) :: nom
	integer :: i, n
	do i = 1, len_trim(nom)
		if (nom(i:i) == achar(9)) nom(i:i) = ' '
	end do
	nom = adjustl(nom)
	n = len_trim(nom)
	if (n >= 2) then
		if ((nom(1:1) == '"' .and. nom(n:n) == '"') .or. (nom(1:1) == "'" .and. nom(n:n) == "'")) nom = nom(2:n-1)
	endif
	end subroutine nom_fichier


!-----------------------------------------------------------------------
! NL ADAPTATIF T0 (ligne NL de la config a -1), 2026-09-29
! Octree jusqu'a NLmax (NL_adaptatif_auto : cellule >= seuil_adapt x la taille
! qui convient a 99,9 % des elements, plafond 10 ; cle NLmax du fichier
! technique pour l'imposer) ; chaque element est ensuite place au niveau le plus
! fin dont la cellule vaut au moins seuil_adapt fois SA plus grande arete.
! Boite comme make_octree_all sur Pc. Remplit h_el.
!-----------------------------------------------------------------------
	subroutine NL_adaptatif_T0(NL)
	integer, intent(out) :: NL
	double precision :: Lc, hmin, cel, h_q
	integer :: iel, a, b
	Lc = maxval(maxval(Pc, dim=1) - minval(Pc, dim=1))*1.0005d0
	allocate(h_el(NE))
	do iel = 1, NE
		h_el(iel) = 0.0d0
		do a = 1, 3
			b = mod(a,3) + 1
			h_el(iel) = max(h_el(iel), norm2(Ps(E(iel,a),:) - Ps(E(iel,b),:)))
		end do
	end do
	hmin = minval(h_el)
	NL = NL_adaptatif_auto(NE, h_el, Lc, seuil_adapt, NLF, NLmax_adapt, h_q)
	cel = Lc/2.0d0**NL
	call log_cle('octree', 'adaptive, '//txt_i(NL)//' levels, box '//txt_r(Lc*1d3,'(F10.1)')//' mm, leaf cell '// &
		txt_r(cel*1d3,'(F10.3)')//' mm')
	if (verbose >= 2) write(log_unit,'(A,F7.3,A)') ' octree: 0.1 % edge quantile', h_q*1d3, ' mm'
	end subroutine NL_adaptatif_T0

	! plus petite et plus grande arete du maillage (journal)
	subroutine aretes_extremes(hmin, hmax)
	double precision, intent(out) :: hmin, hmax
	integer :: iel, a
	double precision :: h
	hmin = huge(1.0d0) ; hmax = 0.0d0
	do iel = 1, NE
		do a = 1, 3
			h = norm2(Ps(E(iel,a),:) - Ps(E(iel,mod(a,3)+1),:))
			hmin = min(hmin, h) ; hmax = max(hmax, h)
		end do
	end do
	end subroutine aretes_extremes

	! valeur (en ko) d'une ligne "cle:  valeur kB" d'un fichier de /proc ; 0 si absente
	! (meme routine que main_beltrami_tube)
	real(prec) function lit_kb(fichier, cle)
	character(len=*), intent(in) :: fichier, cle
	character(len=256) :: l
	integer :: u, ios
	lit_kb = 0.0_prec
	open(newunit=u, file=fichier, status='old', action='read', iostat=ios)
	if (ios /= 0) return
	do
		read(u,'(A)',iostat=ios) l
		if (ios /= 0) exit
		if (index(l, cle) == 1) then
			read(l(len(cle)+1:),*) lit_kb
			exit
		endif
	end do
	close(u)
	end function lit_kb

	END


	
