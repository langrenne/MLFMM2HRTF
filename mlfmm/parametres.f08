! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: MIT
	module parametres
	use module__Precision
	use module__Constants, only: pi, ji => jcomplex
	implicit none
! frequence
	integer :: NF
	real(prec) :: fe,inF,dF     ! dF = pas de frequence (balayage inF + (i-1)*dF)
! FMM
	integer :: Ni,FMM_pre,FMM_meth
	real(prec) :: epsilon0
	integer :: Nmin_ord = 8              ! plancher d'ordre FMM par niveau (accuracy basses freq) ; HRTF : cle Nmin
	integer :: dN_ord = 0                ! increment d'ordre FMM applique apres le plancher (test, defaut 0)
	integer :: gum_debit = 1           ! 13/11 : 1 (defaut) = source Gumerov, debit ramene a 1 (flux manquant sur le triangle source) ; 0 = ancien
	integer :: bascule_reelle = 0        ! essai 10/11 : 1 = test de bascule BF/HF sur N non arrondi (une seule bascule par niveau)
	real(prec) :: plow_coef = 1.0_prec   ! essai 09/11 : plow = 1 - plow_coef log10(eps (1-1/sigma)^1,5) (1 = regle historique)
	integer :: dN_BF = 0                 ! increment d'ordre aux seuls niveaux BF (octree adaptatif, basses freq ; defaut 0)
	real(prec) :: k_boost = huge(1.0_prec) ! dN_BF applique seulement si k < k_boost (HRTF : 2 pi f_boost / c)
	real(prec) :: f_max_chief = 0.0_prec   ! CHIEF (cle interne) : sous f_max_chief (Hz), CHIEF remplace Burton/Miller (ex-f_sans_BM)
	real(prec) :: k_max_chief = -1.0_prec  ! 2 pi f_max_chief / c (main, apres lecture de c)
	real(prec) :: f_min_chief = 0.0_prec ! CHIEF (essai) : inactif sous f_min_chief (Hz), ni BM ni CHIEF (defaut 0)
	real(prec) :: k_min_chief = -1.0_prec ! 2 pi f_min_chief / c (main)
	real(prec) :: C_BM = 3.141592653589793_prec ! Burton/Miller alpha = i k / (k^2 + C_BM) ; C_BM plus grand = alpha plus petit en BF
	! CHI (Cunefare, Koopmann & Brod 1989, essai 2026-10-14) : equations de Burton/Miller ecrites en
	! des points INTERIEURS Pt = Pc - chi_d x sqrt(Jac0) x n (sous chaque element, normale vers le
	! fluide), sans terme 1/2 ni integrale singuliere. chi_d = 0 : collocation sur la surface (BM)
	real(prec) :: chi_d = 0.0_prec
	! CHI, limite par la courbure (2026-10-15) : profondeur <= chi_courbure x rayon de courbure convexe
	! local (sinon les points des bords tres courbes convergent vers le centre de courbure et deux
	! lignes deviennent presque egales). 0 : pas de limite
	real(prec) :: chi_courbure = 0.5_prec
	! CHI, limite par la longueur d'onde (2026-10-29) : a la plus haute frequence du calcul (k_chi_max,
	! fixe par le main), profondeur <= 1 x h si k h <= 2 pi/6 (element valide en lambda/6), 0,75 x h
	! si k h <= 2 pi/3, 0,5 x h au-dela. Les iterations croissent avec k h et chi (sphere, chi_kh.py).
	! k_chi_max <= 0 : pas de limite
	real(prec) :: k_chi_max = -1.0_prec
	! CHIEF (essai 2026-10-09, branche chief) : NC_chief points interieurs ou le champ
	! rayonne par la surface doit etre nul (Ray) ou egal a -p_inc (Diff). Actifs la ou
	! Burton/Miller est coupe (k < k_max_chief) : systeme borde [A Y ; C 0] de taille NE+NC,
	! Y = C^H. Lignes C en dense (quadrature IS_v, comme la propagation).
	character(len=500) :: fichier_chief = ''
	! Ray : modele de la source reciproque (cle source ; ancienne cle source_gum encore lue).
	! 'kreuzer' (ex-piston, defaut depuis 2026-10-15 ; Kreuzer 2009, Mesh2HRTF) : seul le
	! triangle source vibre, dn p = -1/A0 (debit equivalent au monopole unite de la reference
	! au centre), pas de champ incident. 'gumerov' (ex-point, Gumerov 2010, ancien defaut) :
	! champ incident 2G au barycentre du triangle source et -2 dn G ailleurs ; ecart de gain
	! +0,02 dB par mm d'arete du triangle source et erreur HF 2,5 a 4,5 fois plus grande
	! (sphere, CHI, 2026-10-15)
	character(len=20) :: source_ray = 'kreuzer'
	integer :: NC_chief = 0
	! CHIEF-B (essai 2026-10-13, Wu & Seybert 1991) : fichier facultatif chief_comb, matrice
	! creuse G (NC_lignes x NC_chief) : chaque contrainte = combinaison lineaire des points
	! (moyenne ponderee sur un bloc, derivees par differences finies). Sans fichier : G = I.
	character(len=500) :: fichier_chief_comb = ''
	real(prec) :: chief_orth = 0.0_prec   ! > 0 : contraintes orthonormalisees, valeurs propres de C C^H < chief_orth x max ecartees
	integer :: NC_lignes = 0, NZ_comb = 0
	integer, allocatable :: iG_comb(:), jG_comb(:)
	real(prec), allocatable :: vG_comb(:)
	real(prec), allocatable :: Q_chief(:,:)          ! (NC_chief, 3)
	complex(prec), allocatable :: pinc_chief(:,:)    ! Diff : champ incident aux points CHIEF (NC_chief, NF)
	! Bandes de frequence (HRTF, 2026-10-01) : epsilon0 et epsilon2 peuvent prendre une
	! valeur par bande ; limites f_bandes (Hz) croissantes, converties en k_bandes par
	! bandes_en_k (main, apres lecture de c). Une seule valeur = toutes les frequences.
	integer, parameter :: NB_BANDES_MAX = 10
	integer :: nb_bandes = 0
	real(prec) :: f_bandes(NB_BANDES_MAX) = 0.0_prec, k_bandes(NB_BANDES_MAX) = 0.0_prec
	real(prec) :: eps0_bandes(NB_BANDES_MAX+1) = 1.0e-6_prec, eps2_bandes(NB_BANDES_MAX+1) = 1.0e-4_prec
	real(prec) :: seuilHF_bandes(NB_BANDES_MAX+1) = 2.5_prec   ! seuil_HF par bande (HRTF, cle seuil_HF)
	! Estimation initiale GMRES (HRTF, 2026-10-01) : 0 depart a zero ; 1 solution de la
	! frequence precedente ; 2 extrapolation lineaire (en k) des deux precedentes.
	! Si guess > 0, les frequences sont parcourues par BLOCS CONTIGUS (blocs_par_thread
	! blocs par thread, ordre des blocs entrelace pour la RAM), depart a zero en tete de bloc.
	integer :: guess_mode = 2
	integer :: blocs_par_thread = 4
	real(prec) :: ram_max = 0.0_prec   ! HRTF : RAM totale permise (Go) pour l'ordonnancement ; 0 = 85 % de la RAM disponible
	integer :: estimer = 0             ! HRTF : 1 = ecrire la RAM estimee par frequence puis s'arreter (aucun calcul)
	real(prec) :: seuil_HF = 2.5_prec    ! bascule BF -> HF d'un niveau si N - phigh < seuil_HF (init_FMM_T0) ; plus grand = HF plus tot
	integer, parameter :: NLF = 2        ! niveau final FMM (toujours 2)
! Physique
	real(prec), parameter :: c=346.28591544910125_prec     ! célérité (valeur precise ref)
!	real(prec), parameter :: c=340.0_prec     ! célérité
	real(prec), parameter :: rho=1.1844898037900384_prec		! densité du milieu
	real(prec), parameter :: lv=4.466003024146257e-08_prec
	real(prec), parameter :: lh=6.344922504383850e-8_prec
	real(prec), parameter :: Gam=1.402_prec
	real(prec), parameter :: mu=1.83183e-05_prec				! viscosité du flux
	real(prec), parameter :: lvp=mu/(rho*c)

	! ILUT
	real(prec) :: epsilon1,coef_ilut
	integer :: lfil_ilut       ! ILUT : nb max de remplissage par ligne (qualite du preconditionneur, a tatonner)

!	! paramètres gmres
!	integer :: lwork
	integer :: restart = 100            ! GMRES : taille de la base (memoire restart x NE x 16 o) ; HRTF : cle restart
	integer :: iter_max = 150           ! GMRES : iterations max (HRTF : cle iter_max ; tuyau : 100 fixe)
	real(prec) :: epsilon2

! Log
	integer :: verbose  = 1   ! journal : 0 en-tete + bilan, 1 + une ligne par frequence, 2 + details et CHECK
	integer :: log_unit = 6   ! 6=stdout par défaut
	integer :: n_faites = 0    ! journal HRTF : frequences terminees (log_frequence, section critique)
	character(len=10) :: nom_reglage = 'rapide'   ! HRTF : cle reglage du fichier technique (rapide / precis, 2026-11-20)

! Propagation vers les points exterieurs du HRTF (fait_tout_HRTF_T0 et
! fait_tout_HRTF_DIFF_T0) : 0 propagation_FMM_only (ancienne, niveau fixe,
! T_LR a chaque couple) ; 1 propagation_hierarchique (BF, depart au niveau
! NLF, critere 2 x cote)
	integer :: methode_propagation = 1
! Ordre des translations HF entre niveaux d'ordres differents DANS LE SOLVEUR
! (2026-10-30) : 1 (defaut) interpoler PUIS dephaser (montee) et dephaser PUIS
! filtrer (descente), ordre MLFMA ; 0 ancien ordre (repliement). Essai
! interp_20261030 : meme solution a 1e-3 dB pres, +2 % de temps. Voir
! interp_puis_dephase (module_MLFMM), toujours actif pendant la propagation.
	integer :: interp_solveur = 1
! Quadrature du champ proche T0 (ABnear_FMM_T0_adapt, 2026-10-30) : 1 (defaut) nombre
! de points selon S = Hmax/D et k Hmax (7 a 48, erreur <= 1e-7, calibre.py) ; 0 : 33 partout
	integer :: quad_proche = 1
! Quadrature quasi singuliere des points CHI (DMKH_proche_T0, essai 2026-11-03) : 0 (defaut) sinh 16 x 16
! partout ; 1 : selon d/h (Gauss 48 si d/h >= 0,9, sinh 8 / 10 / 12 / 16 points par direction aux seuils
! 0,75 / 0,5 / 0,25), erreur <= 6e-7 sur m3 (quadrature_20261030/proche3.py, proche4.py)
	integer :: quad_chi = 0
	integer :: liberation_V1 = 1     ! Ray T0 adaptatif : 1 libere V1 pendant GMRES (2026-11-06), 0 : ancien comportement
! Profondeur de l'octree adaptatif imposee (0 = automatique, voir
! NL_adaptatif_auto dans module_Octree ; plafonnee a 10)
	integer :: NLmax_adapt = 0

	contains

!------------------------------------------------------------------------------
! lire_parametres_MLFMM (2026-09-30) : parametres TECHNIQUES du HRTF, lus dans
! un fichier « cle = valeur » separe du config.txt (qui ne garde que le cas :
! maillage, points, frequences, source, sauvegarde...). Une ligne par cle,
! dans n'importe quel ordre ; commentaires apres # ou ! ; majuscules
! indifferentes. Une cle ABSENTE prend sa valeur par defaut = REGLAGE RAPIDE du
! preprint (2026-11-20 ; avant : reglage « 60k » du 2026-09-30, sans CHI) ;
! une cle INCONNUE arrete le programme (faute de frappe). fichier vide ou '-'
! : toutes les valeurs par defaut.
!
! CLES PUBLIQUES (MANUAL.md, config.md ; toutes les autres sont internes) :
!   reglage      rapide (defaut) ou precis. Ne pose que bandes / epsilon0 / epsilon2 / Nmin :
!                rapide 5000 8000 / 1e-5 1e-6 1e-7 / 1e-4 / 4 ; precis : sans bandes / 1e-7 / 1e-5 / 8.
!                Lue avant les autres cles, qui la corrigent.
!   estimer      1 : estimation de la RAM puis arret (defaut 0)
!   ram_max      RAM totale permise (Go) pour l'ordonnancement des threads ; 0 (defaut) = 85 % de la RAM disponible
!   source       Ray : kreuzer (defaut, Kreuzer 2009 : seul le triangle source vibre,
!                dn p = -1/A0, pas de champ incident) ou gumerov (Gumerov : 2G + -2 dn G).
!                Ancienne cle source_gum = piston / point encore acceptee (2026-11-20)
!   iter_max     GMRES : nombre maximal d'iterations (defaut 150)
!
! CLES INTERNES (developpement, non documentees ; une ligne « cle interne » dans le log
! quand la valeur differe du defaut, log_cles_internes) :
!   -- reglages figes, caches a l'utilisateur le 2026-11-20 :
!   bandes       limites des bandes en Hz, croissantes ; epsilon0 / epsilon2 / seuil_HF :
!                1 valeur ou nb_limites+1 valeurs
!   epsilon0     precision FMM ; epsilon2 tolerance GMRES
!   FMM_pre      preconditionnement : 0 aucun, 1 ILUT gauche, 2 ILUT droite (defaut)
!   epsilon1     ILUT : seuil d'elimination (1e-4) ; lfil termes max par ligne (25) ;
!   restart      GMRES : taille de la base de Krylov (100), memoire restart x NE x 16 o par thread
!   guess        estimation initiale GMRES : 0 zero, 1 frequence precedente, 2 (defaut)
!                extrapolation lineaire des deux precedentes
!   NL           -1 (defaut) octree adaptatif ; > 2 uniforme de profondeur NL (sans controle
!                cellule / element : ni NLmax ni seuil_adapt n'y interviennent)
!   NLmax        adaptatif seulement : profondeur imposee (0 = automatique, <= 10)
!   coef_ilut    ILUT : facteur de place memoire (1 ; monter par 0,1 si erreur ILUT -2 / -3)
!   Nmin         plancher d'ordre FMM par niveau (pose par reglage : rapide 4, precis 8)
!   seuil_HF     bascule BF -> HF d'un niveau si N - phigh < seuil_HF (defaut 2, une valeur par bande possible)
!   seuil_adapt  adaptatif : cellule >= seuil x plus grande arete de l'element (defaut 1,6)
!   FMM_meth     0 BF, 1 HF, 2 hybride (defaut)
!   blocs        guess > 0 : blocs de frequences par thread (defaut 1000 = frequence par frequence)
!   chi          profondeur CHI (Cunefare 1989) en tailles d'element (defaut 1) ; 0 : Burton/Miller
!                colloque sur la surface
!   chi_courbure CHI : profondeur <= chi_courbure x rayon de courbure convexe local (0,5 ; 0 : sans)
!   gum_debit    source gumerov : 1 (defaut) debit ramene a 1 ; 0 ancien
!   propagation  1 (defaut) propagation_hierarchique ; 0 propagation_FMM_only
!   interp_solveur 1 (defaut) ordre MLFMA (interpoler puis dephaser) ; 0 ancien ordre
!   quad_proche  1 (defaut) points de Gauss du champ proche selon la distance ; 0 : 33 partout
!   quad_chi     1 (defaut) quadrature quasi singuliere CHI selon d/h ; 0 : sinh 16 x 16 partout
!   liberer_V1   1 (defaut) Ray : moments V1 liberes pendant GMRES (-16 % de RAM FMM) ; 0 : gardes
!   -- essais et options de developpement :
!   C_BM         Burton/Miller alpha = i k / (k^2 + C_BM) (defaut 1000)
!   plow_coef    regle de p_low (1) ; bascule_reelle : bascule BF/HF sur N non arrondi (1)
!   dN_BF        increment d'ordre aux niveaux BF sous f_boost Hz (defaut 0 ; -1 : 4 en adaptatif)
!   chief        fichier des points CHIEF (1re ligne NC, puis x y z en m, a l'INTERIEUR du
!                corps, hors des noeuds des modes interieurs) : gardee pour les figures du
!                preprint. CHIEF remplace Burton/Miller entre f_min_chief et f_max_chief (Hz,
!                ex-cle f_sans_BM, encore acceptee) ; sous f_min_chief ni BM ni CHIEF
!   chief_comb   CHIEF-B : matrice de combinaison des points (1re ligne NL NZ, puis NZ lignes
!                i j v : contrainte i += v x point j), blocs de Wu & Seybert 1991
!   chief_orth   > 0 : lignes CHIEF orthonormalisees, valeurs propres < chief_orth x max ecartees
!
! CLES RETIREES le 2026-11-20 (essais sans gain), acceptees seulement a 0 : kr_hf_min,
! bm_partiel, chi_hybride. f_sans_BM (BM coupe en BF) n'existe plus que pour CHIEF.
!------------------------------------------------------------------------------
	subroutine lire_parametres_MLFMM(fichier, NL, seuil_adapt, dN_BF_lu, f_boost)
	character(len=*), intent(in) :: fichier
	integer, intent(out) :: NL, dN_BF_lu
	real(prec), intent(out) :: seuil_adapt, f_boost
	character(len=500) :: ligne, cle, val
	integer :: u, ios, n0, n2, ns, nk
	real(prec), allocatable :: v(:)
	real(prec) :: e0(NB_BANDES_MAX+1), e2(NB_BANDES_MAX+1), sh(NB_BANDES_MAX+1)

	! valeurs par defaut = REGLAGE RAPIDE du preprint (fige le 2026-11-12, ex-parametres_MLFMM_chi_dense.txt) ;
	! la cle « reglage = precis » (lue en premier, quelle que soit sa place) ne change que
	! bandes, epsilon0, epsilon2 et Nmin. Les autres cles du fichier s'appliquent ensuite.
	! (Avant le 2026-11-20 : defauts du reglage « 60k » du 2026-09-30, sans CHI.)
	NL = -1 ; seuil_adapt = 1.6_prec ; FMM_meth = 2 ; FMM_pre = 2
	epsilon1 = 1.0e-4_prec
	coef_ilut = 1.0_prec ; lfil_ilut = 25 ; dN_BF_lu = 0 ; f_boost = 2000.0_prec
	interp_solveur = 1 ; quad_proche = 1 ; quad_chi = 1
	chi_d = 1.0_prec ; chi_courbure = 0.5_prec ; C_BM = 1000.0_prec
	methode_propagation = 1 ; NLmax_adapt = 0 ; guess_mode = 2 ; blocs_par_thread = 1000
	plow_coef = 1.0_prec ; bascule_reelle = 1 ; gum_debit = 1
	iter_max = 150 ; restart = 100 ; ram_max = 0.0_prec ; estimer = 0
	ns = 1 ; sh(1) = 2.0_prec
	nom_reglage = 'rapide'
	if (len_trim(fichier) > 0 .and. trim(fichier) /= '-') then
		open(newunit=u, file=trim(fichier), status='old', action='read', iostat=ios)
		if (ios /= 0) then
			write(*,*) 'Technical file not found: ', trim(fichier)
			stop 1
		endif
		do      ! 1re passe : la cle reglage seule
			read(u,'(A)',iostat=ios) ligne
			if (ios /= 0) exit
			call decouper(ligne, cle, val)
			if (trim(cle) /= 'reglage') cycle
			select case (trim(val))
			case ('rapide', 'fast')    ; nom_reglage = 'rapide'
			case ('precis', 'precise') ; nom_reglage = 'precis'
			case default
				write(*,*) 'Technical file: reglage must be rapide (fast) or precis (precise), not: ', trim(val)
				stop 1
			end select
		end do
		rewind(u)
	endif
	if (nom_reglage == 'precis') then
		nb_bandes = 0 ; n0 = 1 ; e0(1) = 1.0e-7_prec ; n2 = 1 ; e2(1) = 1.0e-5_prec ; Nmin_ord = 8
	else
		nb_bandes = 2 ; f_bandes(1:2) = [5000.0_prec, 8000.0_prec]
		n0 = 3 ; e0(1:3) = [1.0e-5_prec, 1.0e-6_prec, 1.0e-7_prec] ; n2 = 1 ; e2(1) = 1.0e-4_prec ; Nmin_ord = 4
	endif
	if (len_trim(fichier) == 0 .or. trim(fichier) == '-') then
		call fixer_bandes()
		return
	endif

	do
		read(u,'(A)',iostat=ios) ligne
		if (ios /= 0) exit
		call decouper(ligne, cle, val)
		if (len_trim(cle) == 0) cycle
		select case (trim(cle))
		case ('reglage')     ! lue a la 1re passe
		case ('nl')          ; read(val,*,iostat=ios) NL
		case ('seuil_adapt') ; read(val,*,iostat=ios) seuil_adapt
		case ('fmm_meth')    ; read(val,*,iostat=ios) FMM_meth
		case ('fmm_pre')     ; read(val,*,iostat=ios) FMM_pre
		case ('epsilon0')
			call lire_liste_reels(val, v, n0)
			if (n0 < 1 .or. n0 > NB_BANDES_MAX+1) ios = 1
			if (ios == 0) e0(1:n0) = v
		case ('epsilon1')    ; read(val,*,iostat=ios) epsilon1
		case ('epsilon2')
			call lire_liste_reels(val, v, n2)
			if (n2 < 1 .or. n2 > NB_BANDES_MAX+1) ios = 1
			if (ios == 0) e2(1:n2) = v
		case ('bandes')
			call lire_liste_reels(val, v, nb_bandes)
			if (nb_bandes < 1 .or. nb_bandes > NB_BANDES_MAX) ios = 1
			if (ios == 0) f_bandes(1:nb_bandes) = v
			if (ios == 0 .and. nb_bandes > 1) then
				if (any(v(2:) <= v(:nb_bandes-1))) ios = 1
			endif
		case ('coef_ilut')   ; read(val,*,iostat=ios) coef_ilut
		case ('lfil')        ; read(val,*,iostat=ios) lfil_ilut
		case ('dn_bf')       ; read(val,*,iostat=ios) dN_BF_lu
		case ('f_boost')     ; read(val,*,iostat=ios) f_boost
		case ('f_max_chief', 'f_sans_bm') ; read(val,*,iostat=ios) f_max_chief   ! f_sans_bm : ancien nom
		case ('c_bm')        ; read(val,*,iostat=ios) C_BM
		case ('chi')         ; read(val,*,iostat=ios) chi_d
		case ('chi_courbure'); read(val,*,iostat=ios) chi_courbure
		case ('chief')       ; fichier_chief = trim(val)
		case ('chief_comb')  ; fichier_chief_comb = trim(val)
		case ('chief_orth')  ; read(val,*,iostat=ios) chief_orth
		case ('f_min_chief') ; read(val,*,iostat=ios) f_min_chief
		case ('source', 'source_gum')   ! source_gum : ancienne cle (piston / point), 2026-11-20
			select case (trim(val))
			case ('kreuzer', 'Kreuzer', 'piston') ; source_ray = 'kreuzer'
			case ('gumerov', 'Gumerov', 'point')  ; source_ray = 'gumerov'
			case default ; ios = 1
			end select
		case ('propagation') ; read(val,*,iostat=ios) methode_propagation
		case ('interp_solveur'); read(val,*,iostat=ios) interp_solveur
		case ('quad_proche') ; read(val,*,iostat=ios) quad_proche
		case ('quad_chi')    ; read(val,*,iostat=ios) quad_chi
		case ('liberer_v1')  ; read(val,*,iostat=ios) liberation_V1
		case ('kr_hf_min', 'bm_partiel', 'chi_hybride')
			! essais retires le 2026-11-20 : acceptes a 0 (fichiers du preprint), sinon arret
			call lire_liste_reels(val, v, nk)
			if (nk < 1 .or. any(v /= 0.0_prec)) then
				write(*,*) 'Technical file: key removed on 2026-11-20 (trial without gain): ', trim(cle)
				stop 1
			endif
		case ('nlmax')       ; read(val,*,iostat=ios) NLmax_adapt
		case ('seuil_hf')
			call lire_liste_reels(val, v, ns)
			if (ns < 1 .or. ns > NB_BANDES_MAX+1) ios = 1
			if (ios == 0) sh(1:ns) = v
		case ('nmin')        ; read(val,*,iostat=ios) Nmin_ord
		case ('plow_coef')   ; read(val,*,iostat=ios) plow_coef
		case ('bascule_reelle') ; read(val,*,iostat=ios) bascule_reelle
		case ('gum_debit')   ; read(val,*,iostat=ios) gum_debit
		case ('guess')       ; read(val,*,iostat=ios) guess_mode
		case ('iter_max')    ; read(val,*,iostat=ios) iter_max
		case ('restart')     ; read(val,*,iostat=ios) restart
		case ('blocs')       ; read(val,*,iostat=ios) blocs_par_thread
		case ('ram_max')     ; read(val,*,iostat=ios) ram_max
		case ('estimer')     ; read(val,*,iostat=ios) estimer
		case default
			write(*,*) 'Technical file ', trim(fichier), ': unknown key "', trim(cle), '"'
			stop 1
		end select
		if (ios /= 0) then
			write(*,*) 'Technical file: unreadable value for ', trim(cle), ': ', trim(val)
			stop 1
		endif
	end do
	close(u)
	call fixer_bandes()

	contains

	subroutine fixer_bandes()
	! une valeur -> toutes les bandes ; sinon exactement nb_bandes+1 valeurs
	if (n0 /= 1 .and. n0 /= nb_bandes+1) then
		write(*,*) 'Technical file: epsilon0 needs 1 or', nb_bandes+1, ' values (bands)'
		stop 1
	endif
	if (n2 /= 1 .and. n2 /= nb_bandes+1) then
		write(*,*) 'Technical file: epsilon2 needs 1 or', nb_bandes+1, ' values (bands)'
		stop 1
	endif
	if (ns /= 1 .and. ns /= nb_bandes+1) then
		write(*,*) 'Technical file: seuil_HF needs 1 or', nb_bandes+1, ' values (bands)'
		stop 1
	endif
	if (n0 == 1) e0(2:nb_bandes+1) = e0(1)
	if (n2 == 1) e2(2:nb_bandes+1) = e2(1)
	if (ns == 1) sh(2:nb_bandes+1) = sh(1)
	eps0_bandes(1:nb_bandes+1) = e0(1:nb_bandes+1)
	eps2_bandes(1:nb_bandes+1) = e2(1:nb_bandes+1)
	seuilHF_bandes(1:nb_bandes+1) = sh(1:nb_bandes+1)
	epsilon0 = eps0_bandes(1) ; epsilon2 = eps2_bandes(1)   ! valeurs « scalaires » : 1re bande
	seuil_HF = seuilHF_bandes(1)
	end subroutine fixer_bandes

	subroutine decouper(ligne, cle, val)
	! « cle = valeur » : commentaires apres # ou ! retires, cle en minuscules ; ligne vide : cle = ''
	character(len=*), intent(inout) :: ligne
	character(len=*), intent(out) :: cle, val
	integer :: k, j
	cle = '' ; val = ''
	k = scan(ligne,'#!')
	if (k > 0) ligne = ligne(1:k-1)
	if (len_trim(ligne) == 0) return
	k = index(ligne,'=')
	if (k == 0) then
		write(*,*) 'Technical file: line without "=": ', trim(ligne)
		stop 1
	endif
	cle = adjustl(ligne(1:k-1)) ; val = adjustl(ligne(k+1:))
	do j = 1, len_trim(cle)
		if (cle(j:j) >= 'A' .and. cle(j:j) <= 'Z') cle(j:j) = achar(iachar(cle(j:j)) + 32)
	end do
	end subroutine decouper
	end subroutine lire_parametres_MLFMM

!------------------------------------------------------------------------------
! Bandes de frequence : conversion des limites en nombre d'onde (c du config),
! et epsilon0 / epsilon2 / seuil_HF de la bande contenant k (limite incluse dans
! la bande superieure). Sans bandes : les valeurs uniques (le tuyau, qui ne lit
! pas le fichier technique, garde seuil_HF = 2,5).
!------------------------------------------------------------------------------
	subroutine bandes_en_k(c_son)
	real(prec), intent(in) :: c_son
	k_bandes(1:nb_bandes) = 2.0_prec*pi*f_bandes(1:nb_bandes)/c_son
	end subroutine bandes_en_k

	pure integer function bande_de_k(k)
	real(prec), intent(in) :: k
	bande_de_k = 1 + count(k_bandes(1:nb_bandes) <= k)
	end function bande_de_k

	real(prec) function eps0_de_k(k)
	real(prec), intent(in) :: k
	if (nb_bandes == 0) then
		eps0_de_k = epsilon0
	else
		eps0_de_k = eps0_bandes(bande_de_k(k))
	endif
	end function eps0_de_k

	real(prec) function eps2_de_k(k)
	real(prec), intent(in) :: k
	if (nb_bandes == 0) then
		eps2_de_k = epsilon2
	else
		eps2_de_k = eps2_bandes(bande_de_k(k))
	endif
	end function eps2_de_k

	real(prec) function seuil_de_k(k)
	real(prec), intent(in) :: k
	if (nb_bandes == 0) then
		seuil_de_k = seuil_HF
	else
		seuil_de_k = seuilHF_bandes(bande_de_k(k))
	endif
	end function seuil_de_k



!------------------------------------------------------------------------------
! lire_liste_reels (2026-09-30) : tous les nombres d'une ligne de texte
! (separes par des espaces, tabulations ou virgules), dans v(1:n).
! Sert a la ligne « Frequence initiale » du config HRTF : une valeur = mode
! historique (frequence unique si NF = 1, grille i fe/(2 NF) sinon) ;
! plusieurs valeurs = LISTE des frequences a calculer.
!------------------------------------------------------------------------------
	subroutine lire_liste_reels(ligne, v, n)
	character(len=*), intent(in) :: ligne
	real(prec), allocatable, intent(out) :: v(:)
	integer, intent(out) :: n
	integer :: i, ios
	logical :: dans
	character(len=len(ligne)) :: l
	l = ligne
	i = scan(l, '!#')
	if (i > 0) l(i:) = ' '
	do i = 1, len(l)
		if (l(i:i) == ',' .or. l(i:i) == achar(9)) l(i:i) = ' '
	end do
	n = 0 ; dans = .false.
	do i = 1, len_trim(l)
		if (l(i:i) /= ' ' .and. .not. dans) n = n + 1
		dans = (l(i:i) /= ' ')
	end do
	allocate(v(n))
	if (n > 0) then
		read(l,*,iostat=ios) v
		if (ios /= 0) then
			write(*,*) 'Case file: unreadable frequency line: ', trim(ligne)
			stop 1
		endif
	endif
	end subroutine lire_liste_reels

!------------------------------------------------------------------------------
! log_cles_internes (2026-11-20) : une ligne par cle INTERNE (non documentee pour
! l'utilisateur) dont la valeur differe du defaut, pour qu'un resultat reste tracable.
!------------------------------------------------------------------------------
	subroutine log_cles_internes(f_boost, seuil_adapt, NL)
	real(prec), intent(in) :: f_boost, seuil_adapt
	integer, intent(in) :: NL
	if (NL /= -1) write(log_unit,'(A,I0)') ' internal key: NL = ', NL
	if (nom_reglage == 'precis') then
		if (nb_bandes /= 0) write(log_unit,'(A,10F9.0)') ' internal key: bandes = ', f_bandes(1:nb_bandes)
		if (any(abs(eps0_bandes(1:nb_bandes+1) - 1.0e-7_prec) > 1.0e-20_prec)) &
			write(log_unit,'(A,11ES9.1)') ' internal key: epsilon0 = ', eps0_bandes(1:nb_bandes+1)
		if (any(abs(eps2_bandes(1:nb_bandes+1) - 1.0e-5_prec) > 1.0e-20_prec)) &
			write(log_unit,'(A,11ES9.1)') ' internal key: epsilon2 = ', eps2_bandes(1:nb_bandes+1)
	else
		if (nb_bandes /= 2) then
			write(log_unit,'(A,10F9.0)') ' internal key: bandes = ', f_bandes(1:nb_bandes)
		else if (any(abs(f_bandes(1:2) - [5000.0_prec, 8000.0_prec]) > 1.0e-9_prec)) then
			write(log_unit,'(A,10F9.0)') ' internal key: bandes = ', f_bandes(1:nb_bandes)
		endif
		if (nb_bandes /= 2) then
			write(log_unit,'(A,11ES9.1)') ' internal key: epsilon0 = ', eps0_bandes(1:nb_bandes+1)
		else if (any(abs(eps0_bandes(1:3) - [1.0e-5_prec, 1.0e-6_prec, 1.0e-7_prec]) > 1.0e-20_prec)) then
			write(log_unit,'(A,11ES9.1)') ' internal key: epsilon0 = ', eps0_bandes(1:nb_bandes+1)
		endif
		if (any(abs(eps2_bandes(1:nb_bandes+1) - 1.0e-4_prec) > 1.0e-20_prec)) &
			write(log_unit,'(A,11ES9.1)') ' internal key: epsilon2 = ', eps2_bandes(1:nb_bandes+1)
	endif
	if (FMM_pre /= 2) write(log_unit,'(A,I0)') ' internal key: FMM_pre = ', FMM_pre
	if (abs(epsilon1 - 1.0e-4_prec) > 1.0e-20_prec) write(log_unit,'(A,ES9.1)') ' internal key: epsilon1 = ', epsilon1
	if (lfil_ilut /= 25) write(log_unit,'(A,I0)') ' internal key: lfil = ', lfil_ilut
	if (restart /= 100) write(log_unit,'(A,I0)') ' internal key: restart = ', restart
	if (guess_mode /= 2) write(log_unit,'(A,I0)') ' internal key: guess = ', guess_mode
	if (NLmax_adapt /= 0) write(log_unit,'(A,I0)') ' internal key: NLmax = ', NLmax_adapt
	if (abs(coef_ilut - 1.0_prec) > 1.0e-12_prec) write(log_unit,'(A,G0.6)') ' internal key: coef_ilut = ', coef_ilut
	if (Nmin_ord /= merge(8, 4, nom_reglage == 'precis')) write(log_unit,'(A,I0)') ' internal key: Nmin = ', Nmin_ord
	if (any(abs(seuilHF_bandes(1:nb_bandes+1) - 2.0_prec) > 1.0e-12_prec)) &
		write(log_unit,'(A,11F6.2)') ' internal key: seuil_HF = ', seuilHF_bandes(1:nb_bandes+1)
	if (abs(seuil_adapt - 1.6_prec) > 1.0e-12_prec) write(log_unit,'(A,G0.6)') ' internal key: seuil_adapt = ', seuil_adapt
	if (FMM_meth /= 2) write(log_unit,'(A,I0)') ' internal key: FMM_meth = ', FMM_meth
	if (blocs_par_thread /= 1000) write(log_unit,'(A,I0)') ' internal key: blocs = ', blocs_par_thread
	if (abs(chi_d - 1.0_prec) > 1.0e-12_prec) write(log_unit,'(A,G0.6)') ' internal key: chi = ', chi_d
	if (abs(chi_courbure - 0.5_prec) > 1.0e-12_prec) write(log_unit,'(A,G0.6)') ' internal key: chi_courbure = ', chi_courbure
	if (gum_debit /= 1) write(log_unit,'(A,I0)') ' internal key: gum_debit = ', gum_debit
	if (methode_propagation /= 1) write(log_unit,'(A,I0)') ' internal key: propagation = ', methode_propagation
	if (interp_solveur /= 1) write(log_unit,'(A,I0)') ' internal key: interp_solveur = ', interp_solveur
	if (quad_proche /= 1) write(log_unit,'(A,I0)') ' internal key: quad_proche = ', quad_proche
	if (quad_chi /= 1) write(log_unit,'(A,I0)') ' internal key: quad_chi = ', quad_chi
	if (liberation_V1 /= 1) write(log_unit,'(A,I0)') ' internal key: liberer_V1 = ', liberation_V1
	if (abs(C_BM - 1000.0_prec) > 1.0e-12_prec) write(log_unit,'(A,G0.6)') ' internal key: C_BM = ', C_BM
	if (abs(plow_coef - 1.0_prec) > 1.0e-12_prec) write(log_unit,'(A,G0.6)') ' internal key: plow_coef = ', plow_coef
	if (bascule_reelle /= 1) write(log_unit,'(A,I0)') ' internal key: bascule_reelle = ', bascule_reelle
	if (dN_BF /= 0) write(log_unit,'(A,I0,A,G0.6,A)') ' internal key: dN_BF = ', dN_BF, ' below f_boost = ', f_boost, ' Hz'
	if (NC_chief > 0) write(log_unit,'(A,I0,A,G0.6,A,G0.6,A)') ' internal key: chief, ', NC_chief, &
		' points, f_min_chief = ', f_min_chief, ' Hz, f_max_chief = ', f_max_chief, ' Hz'
	end subroutine log_cles_internes

!------------------------------------------------------------------------------
! CHIEF (essai) : lecture des points interieurs (fichier_chief, cle chief du
! fichier technique). 1re ligne : NC ; puis NC lignes x y z (m). Sans fichier :
! NC_chief = 0, rien ne change.
!------------------------------------------------------------------------------
	subroutine lire_points_chief()
	integer :: u, ios, i
	NC_chief = 0
	if (len_trim(fichier_chief) == 0) return
	open(newunit=u, file=trim(fichier_chief), status='old', action='read', iostat=ios)
	if (ios /= 0) then
		write(*,*) 'CHIEF: file not found: ', trim(fichier_chief)
		stop 1
	endif
	read(u,*,iostat=ios) NC_chief
	if (ios /= 0 .or. NC_chief < 1) then
		write(*,*) 'CHIEF: unreadable number of points in ', trim(fichier_chief)
		stop 1
	endif
	allocate(Q_chief(NC_chief,3))
	do i = 1, NC_chief
		read(u,*,iostat=ios) Q_chief(i,:)
		if (ios /= 0) then
			write(*,*) 'CHIEF: unreadable point ', i, ' in ', trim(fichier_chief)
			stop 1
		endif
	end do
	close(u)
	write(log_unit,'(A,I4,A,A)') ' CHIEF (trial):', NC_chief, ' interior points read from ', trim(fichier_chief)
	NC_lignes = NC_chief
	if (len_trim(fichier_chief_comb) == 0) return
	open(newunit=u, file=trim(fichier_chief_comb), status='old', action='read', iostat=ios)
	if (ios /= 0) then
		write(*,*) 'CHIEF: combination file not found: ', trim(fichier_chief_comb)
		stop 1
	endif
	read(u,*,iostat=ios) NC_lignes, NZ_comb
	if (ios /= 0 .or. NC_lignes < 1 .or. NZ_comb < 1) then
		write(*,*) 'CHIEF: unreadable header in ', trim(fichier_chief_comb)
		stop 1
	endif
	allocate(iG_comb(NZ_comb), jG_comb(NZ_comb), vG_comb(NZ_comb))
	do i = 1, NZ_comb
		read(u,*,iostat=ios) iG_comb(i), jG_comb(i), vG_comb(i)
		if (ios /= 0 .or. iG_comb(i) < 1 .or. iG_comb(i) > NC_lignes .or. jG_comb(i) < 1 .or. jG_comb(i) > NC_chief) then
			write(*,*) 'CHIEF: term ', i, ' unreadable or out of range in ', trim(fichier_chief_comb)
			stop 1
		endif
	end do
	close(u)
	write(log_unit,'(A,I5,A,A)') ' CHIEF-B (trial):', NC_lignes, ' combined constraints read from ', trim(fichier_chief_comb)
	end subroutine lire_points_chief



!-----------------------------------------------------------------------
! ADMITTANCE DE PAROI (commune a main_HRTF_Diff et main_HRTF_Ray/Kreuzer, 2026-11-06)
! ligne = « Rmax_BF Rmax_HF » : |R| lineaire de Rmax_BF a 0 Hz a Rmax_HF a 24 kHz, phase 0,17 rad
!   (Ghorbal et al., AES 2020), R complexe, beta = (1 - R)/(1 + R), |R| >= 1 : paroi rigide ;
! sinon : nom d'un fichier (lire_admittance). Formules reprises de main_HRTF_Diff (bit pres).
!-----------------------------------------------------------------------
	subroutine admittance_paroi(ligne, f, nf, beta)
	character(len=*), intent(in) :: ligne
	integer, intent(in) :: nf
	double precision, intent(in) :: f(nf)
	double complex, intent(out) :: beta(nf)
	double precision :: Rmax(2), Refl(nf)
	integer :: ios
	read(ligne,*,iostat=ios) Rmax(2), Rmax(1)          ! BF d'abord, HF ensuite
	if (ios /= 0) then
		call lire_admittance(trim(adjustl(ligne)), f, nf, beta)
		return
	endif
	Refl = (Rmax(1)-Rmax(2)) / 24000.0d0 * f + Rmax(2)
	where (Refl < 1.0d0)
		beta = (1.0d0 - Refl*exp(ji*0.17)) / (Refl*exp(ji*0.17) + 1.0d0) ! 0.17 rad selon Ghorbal, aes 2020
	elsewhere
		beta = dcmplx(0.0d0,0.0d0)                      ! paroi rigide
	end where
	end subroutine admittance_paroi

!-----------------------------------------------------------------------
! ADMITTANCE DE PAROI LUE DANS UN FICHIER (2026-10-07)
! Convention (= NumCalc/Mesh2HRTF) : e^{-iwt}, normale vers le fluide,
! Z = p / vitesse ENTRANT dans la paroi, beta = rho c / Z, dnp = -i k beta p ;
! paroi absorbante : Re(beta) > 0.
! Format : lignes vides et commentaires (# ou !) ignores, sauf une ligne
!   type = reflexion | admittance | impedance      (avec ou sans # devant)
! puis une ligne par frequence, f croissantes :
!   reflexion  : f(Hz)  |R|      phase(rad)   -> beta = (1 - R)/(1 + R), |R| >= 1 = rigide
!   admittance : f(Hz)  Re(beta) Im(beta)     (admittance reduite)
!   impedance  : f(Hz)  Re(z)    Im(z)        (impedance reduite z = Z/(rho c), beta = 1/z)
! Interpolation LINEAIRE en f des deux colonnes (|R| et phase, ou Re et Im) ;
! hors de la plage du fichier : valeur de l'extremite (message dans le log).
!-----------------------------------------------------------------------
	subroutine lire_admittance(fichier, f, nf, beta)
	character(len=*), intent(in) :: fichier
	integer, intent(in) :: nf
	double precision, intent(in) :: f(nf)
	double complex, intent(out) :: beta(nf)
	character(len=500) :: l, t
	character(len=12) :: typ
	double precision :: v(3), w, a, b
	double precision, allocatable :: tf(:), ta(:), tb(:)
	double complex :: Rc
	integer :: u, ios, n, i, j, nhors_bas, nhors_haut
	typ = ' '
	allocate(tf(100000), ta(100000), tb(100000))
	open(newunit=u, file=fichier, status='old', iostat=ios)
	if (ios /= 0) then
		write(*,*) 'Skin impedance file not found: ', trim(fichier)
		stop 1
	endif
	n = 0
	do
		read(u,'(A)',iostat=ios) l
		if (ios /= 0) exit
		t = adjustl(l)
		if (t(1:1) == '#' .or. t(1:1) == '!') t = adjustl(t(2:))
		i = index(t, '=')
		if (i > 0 .and. index(t(1:max(i-1,1)), 'type') > 0) then
			typ = adjustl(t(i+1:))
			cycle
		endif
		t = adjustl(l)
		if (len_trim(t) == 0 .or. t(1:1) == '#' .or. t(1:1) == '!') cycle
		read(t,*,iostat=ios) v
		if (ios /= 0) then
			write(*,*) 'Skin impedance file ', trim(fichier), ': unreadable line: ', trim(l)
			stop 1
		endif
		n = n + 1
		if (n > size(tf)) stop 'Skin impedance file: more than 100000 lines'
		tf(n) = v(1) ; ta(n) = v(2) ; tb(n) = v(3)
		if (n > 1) then
			if (tf(n) <= tf(n-1)) then
				write(*,*) 'Skin impedance file ', trim(fichier), ': frequencies not increasing at ', tf(n)
				stop 1
			endif
		endif
	end do
	close(u)
	if (typ == 'reflection') typ = 'reflexion'          ! orthographe anglaise acceptee
	if (typ /= 'reflexion' .and. typ /= 'admittance' .and. typ /= 'impedance') then
		write(*,*) 'Skin impedance file ', trim(fichier), &
		     ': line "type = reflexion | admittance | impedance" missing or unknown: ', trim(typ)
		stop 1
	endif
	if (n < 1) then
		write(*,*) 'Skin impedance file ', trim(fichier), ': no values'
		stop 1
	endif
	nhors_bas = 0 ; nhors_haut = 0
	do i = 1, nf
		if (f(i) <= tf(1)) then
			a = ta(1) ; b = tb(1)
			if (f(i) < tf(1)) nhors_bas = nhors_bas + 1
		elseif (f(i) >= tf(n)) then
			a = ta(n) ; b = tb(n)
			if (f(i) > tf(n)) nhors_haut = nhors_haut + 1
		else
			j = 1
			do while (tf(j+1) < f(i))
				j = j + 1
			end do
			w = (f(i) - tf(j)) / (tf(j+1) - tf(j))
			a = (1.0d0 - w)*ta(j) + w*ta(j+1)
			b = (1.0d0 - w)*tb(j) + w*tb(j+1)
		endif
		select case (trim(typ))
		case ('reflexion')
			if (a >= 1.0d0) then
				beta(i) = dcmplx(0.0d0, 0.0d0)          ! |R| >= 1 : paroi rigide
			else
				Rc = a*exp(ji*b)
				beta(i) = (1.0d0 - Rc) / (1.0d0 + Rc)
			endif
		case ('admittance')
			beta(i) = dcmplx(a, b)
		case ('impedance')
			if (a == 0.0d0 .and. b == 0.0d0) then
				write(*,*) 'Skin impedance file: zero impedance at ', f(i), ' Hz'
				stop 1
			endif
			beta(i) = 1.0d0 / dcmplx(a, b)
		end select
	end do
	write(log_unit,'(A,A,A,A,A,I6,A,F10.1,A,F10.1,A)') ' skin impedance file ', trim(fichier), &
	     ', type ', trim(typ), ',', n, ' values from ', tf(1), ' to ', tf(n), ' Hz'
	write(log_unit,'(A,2ES12.4,A,2ES12.4)') '   beta(f1) =', beta(1), '   beta(fN) =', beta(nf)
	if (nhors_bas + nhors_haut > 0) write(log_unit,'(A,I6,A,I6,A)') '   WARNING:', nhors_bas, &
	     ' frequencies below the file range,', nhors_haut, ' above: end value used'
	deallocate(tf, ta, tb)
	end subroutine lire_admittance

!------------------------------------------------------------------------------
! JOURNAL HRTF EN ANGLAIS (2026-11-23). Niveaux (valeur verbose du fichier de cas) :
!   0 : en-tete (cas, maillage, reglage, octree, RAM) et bilan ;
!   1 : + une ligne par frequence terminee (production) ;
!   2 : + details du solveur par frequence et lignes CHECK des tests ;
!   3, 4 : + erreur arriere vraie, historique de convergence GMRES.
!------------------------------------------------------------------------------
	subroutine log_cle(cle, texte)          ! « key          : texte »
	character(len=*), intent(in) :: cle, texte
	character(len=13) :: c
	c = cle                                  ! cle alignee a gauche sur 13 caracteres
	write(log_unit,'(1X,A,A,A)') c, ': ', trim(texte)
	end subroutine log_cle

	function txt_r(x, f) result(t)          ! reel -> texte (format f, ex. '(F10.3)')
	real(prec), intent(in) :: x
	character(len=*), intent(in) :: f
	character(len=:), allocatable :: t
	character(len=64) :: b
	write(b,f) x
	t = trim(adjustl(b))
	end function txt_r

	function txt_i(n) result(t)             ! entier -> texte
	integer, intent(in) :: n
	character(len=:), allocatable :: t
	character(len=24) :: b
	write(b,'(I0)') n
	t = trim(b)
	end function txt_i

	function txt_duree(t) result(s)          ! secondes -> « 1 h 05 min », « 12 min 30 s », « 4.2 s »
	real(prec), intent(in) :: t
	character(len=:), allocatable :: s
	character(len=40) :: b
	integer :: m
	if (t >= 3600.0_prec) then
		m = nint(t/60.0_prec)
		write(b,'(I0," h ",I2.2," min")') m/60, mod(m,60)
	elseif (t >= 60.0_prec) then
		m = nint(t)
		write(b,'(I0," min ",I2.2," s")') m/60, mod(m,60)
	else
		write(b,'(F6.1," s")') t
		b = adjustl(b)
	endif
	s = trim(b)
	end function txt_duree

	! en-tete du tableau des frequences (verbose >= 1)
	subroutine log_entete_frequences(nf)
	integer, intent(in) :: nf
	n_faites = 0
	if (verbose < 1) return
	write(log_unit,'(A)') ''
	write(log_unit,'(A,I0,A)') ' frequencies as they finish (', nf, '):'
	write(log_unit,'(A)') '       #   freq (Hz)   iter   time (s)   done'
	end subroutine log_entete_frequences

	! une ligne par frequence terminee (appelee dans la boucle parallele)
	subroutine log_frequence(i, f, iter, t, nf)
	integer, intent(in) :: i, iter, nf
	real(prec), intent(in) :: f, t
	!$omp critical (journal_frequence)
	n_faites = n_faites + 1
	if (verbose >= 1) write(log_unit,'(I8,F12.3,I7,F11.1,3X,I0,"/",I0,A)') i, f, iter, t, n_faites, nf, &
		trim(merge('   iter_max reached', '                   ', iter >= iter_max))
	!$omp end critical (journal_frequence)
	end subroutine log_frequence

	! bilan : duree, iterations, frequences a iter_max, fichier des resultats
	subroutine log_bilan(nf, f, iters, t_total, fichier)
	integer, intent(in) :: nf, iters(nf)
	real(prec), intent(in) :: f(nf), t_total
	character(len=*), intent(in) :: fichier
	integer :: i, n
	character(len=:), allocatable :: l
	write(log_unit,'(A)') ''
	call log_cle('done', txt_i(nf)//' frequencies in '//txt_duree(t_total))
	call log_cle('iterations', txt_i(minval(iters))//' to '//txt_i(maxval(iters))//' (iter_max '//txt_i(iter_max)//')')
	n = count(iters >= iter_max)
	if (n == 0) then
		call log_cle('not converged', 'none')
	else
		l = txt_i(n)//' frequencies reached iter_max (see WARNING lines):'
		do i = 1, nf
			if (iters(i) >= iter_max) l = l//' '//txt_r(f(i), '(F10.1)')
			if (len(l) > 400) exit
		end do
		call log_cle('not converged', l)
	endif
	call log_cle('results', fichier)
	end subroutine log_bilan

	end module parametres
