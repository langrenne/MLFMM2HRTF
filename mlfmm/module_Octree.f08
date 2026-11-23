! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: MIT
    MODULE module_Octree
    USE module__Precision
    IMPLICIT none

	! Octree rangé par niveau : octree(il), il = 1..Nlev, niveau réel il+NLF-1.
	! octree(il)%NL = NOMBRE DE CELLULES du niveau (pas un nombre de niveaux).
	! Toutes les cellules d'un niveau ont le même côté (cell(1)%length).
	! Points (Ptot) : 1..NE barycentres des éléments, NE+1..NE+NP noeuds.
	!
	! Rempli par make_octree_all : child, voisin, interac, center, length.
	!   child : au niveau le plus fin, numéros de POINTS ; ailleurs, numéros
	!           de CELLULES du niveau il+1.
	! Rempli par marquer_sources_cibles :
	!   a_source : la cellule contient un élément de niveau >= il (son
	!              multipôle au niveau il est non nul)
	!   a_cible  : la cellule contient au moins un noeud (point > NE) évalué
	!              à ce niveau ou plus bas
	!   feuille  : (arrêt des cibles) ses noeuds sont évalués à ce niveau
	! Rempli par completer_octree (octree adaptatif, 2026-09-29) :
	!   parent : cellule mère au niveau il-1 (0 au niveau il=1)
	!   pts    : TOUS les points contenus dans la cellule, quel que soit le
	!            niveau (au niveau le plus fin : child, dans le même ordre)
	!   ipt(p) : cellule du point p à ce niveau
	Type cells
		integer(iprec), dimension(:), allocatable :: child, voisin, interac
		integer(iprec), dimension(:), allocatable :: pts
		integer(iprec) :: Nb_c, Nb_v, Nb_i
		integer(iprec) :: Nb_p = 0, parent = 0
		! adaptatif (marquer_sources_cibles) ; VRAI par défaut : le chemin
		! uniforme ne saute rien
		logical :: a_source = .true.   ! multipôle non nul à ce niveau
		logical :: a_cible  = .true.   ! contient au moins un noeud à évaluer ici ou plus bas
		logical :: feuille  = .false.  ! ses noeuds sont évalués à CE niveau (< Nlev)
		real(prec), dimension(1:3) :: center
		real(prec) :: length
	end type cells

	Type tree
		type(cells), allocatable, dimension(:) :: cell
		integer(iprec), allocatable, dimension(:) :: ipt
		integer(iprec) :: NL
	end type tree


	! Profondeur maximale de l'octree. L'ancienne construction (avant le 2026-11-23)
	! numerotait les cellules par i + 2^NL j + 4^NL k en entiers 32 bits : au-dela
	! de NL = 10 les cellules comme les voisinages devenaient FAUX sans erreur
	! (constate le 2026-09-30, id13 et christophe_final a NL 11 et 13).
	! make_octree_all calcule maintenant en 64 bits ; la limite reste a 10 (aucun
	! calcul valide au-dela).
	integer(iprec), parameter :: NL_MAX_OCTREE = 10
	! Facteur sur le cote de la boite racine (2026-10-01, essai) : > 1 agrandit les
	! cellules de tous les niveaux pour caler la cellule fine sur la taille des
	! elements (ex. m2 : elements ~1,95 mm, cellule NL8 1,935 mm -> tous au niveau 7).
	! Defaut 1 : inchange. L'appelant doit multiplier son Lc (NL_adaptatif_auto) de meme.
	real(prec) :: facteur_boite = 1.0_prec

    PRIVATE
    PUBLIC  :: make_octree_all,tree,facteur_boite
    PUBLIC  :: completer_octree,niveaux_elements,compte_paires_near,marquer_sources_cibles
    PUBLIC  :: NL_adaptatif_auto, NL_MAX_OCTREE

    CONTAINS

!
!------------------------------------------------------------------------------------
!
! make_octree_all : octree des Npoints points, niveaux NLF (le plus gros) a NL
! (le plus fin), range dans Octfin(il), il = l-NLF+1 pour le niveau l.
! Reecrit le 2026-11-23 (C. Langrenne, avec Claude) ; principe de la construction
! par identifiants de cellules decrit par S. D. Gedney (methodes integrales
! pour l'electromagnetisme) :
!   - boite cubique de cote Lc englobant les points, cellules de cote Lc/2^l ;
!   - au niveau l, la cellule (i, j, k) a l'identifiant i + 2^l j + 4^l k ;
!   - on trie les objets (points au niveau le plus fin, cellules du niveau
!     inferieur ailleurs) par identifiant de leur cellule : les cellules non
!     vides d'un niveau sont numerotees dans l'ordre croissant des identifiants,
!     leurs enfants (child) dans l'ordre du tri ;
!   - voisins (voisin) : cellules non vides qui touchent (decalages dans
!     {-1,0,1}^3), interactions (interac) : enfants des voisins du parent qui
!     ne touchent pas (decalages de -2 a 3, ou de -3 a 2 selon la parite) ;
!     au niveau NLF, interac = toutes les cellules sauf la cellule et ses voisins.
!     Decalages parcourus dans l'ordre lexicographique (i, j, k).
! Le tri est un tri par tas (Williams 1964, Floyd 1964), non stable : l'ordre
! des enfants d'une cellule est celui du tas, ce qui reproduit au bit l'octree
! d'origine (verifie le 2026-11-23).
!
!------------------------------------------------------------------------------------
!
	subroutine make_octree_all(NPoints,Points,Octfin,NL,NLF)
	implicit none
	integer(iprec) :: Npoints, NL, NLF
	real(prec) :: Points(Npoints,3)
	type(tree) :: Octfin(NL-NLF+1)

	integer, parameter :: i8 = selected_int_kind(18)
	real(prec) :: pmin(3), pmax(3), dx(3), bmin(3), bmax(3), centre(3), coin(3), Lc, dl
	integer(i8), allocatable :: cle(:), ids(:)
	integer(iprec), allocatable :: obj(:), debut(:), ijk(:,:)
	integer(iprec) :: proche(3,26), loin(3,189,0:7)
	integer(i8) :: n1, n2
	integer(iprec) :: l, il, n, ng, g, i, m, nv, ni, par
	integer(iprec) :: liste(189)
	logical, allocatable :: exclu(:)

	if (NL > NL_MAX_OCTREE) then
		write(*,*) 'make_octree_all : NL =', NL, ' > ', NL_MAX_OCTREE, &
			': cell numbers exceed 32-bit integers (octree would be wrong). Stop.'
		stop 1
	endif
	call decalages(proche, loin)

	! --- boite : cube de cote Lc, centre sur la boite englobante elargie de 0,05 %
	pmin = minval(Points, dim=1)
	pmax = maxval(Points, dim=1)
	Lc = maxval(pmax - pmin)*(1.0_prec + 0.0005_prec)*facteur_boite
	dx = (pmax - pmin)/2.0_prec*0.0005_prec
	bmax = pmax + dx
	bmin = pmin - dx
	centre = (bmax + bmin)/2.0_prec
	coin = bmin - ((bmin + Lc/2.0_prec) - centre)       ! coin (xmin, ymin, zmin) du cube
	n1 = 2_i8**NL ; n2 = n1*n1
	dl = Lc/real(n1, prec)

	! --- niveau le plus fin : objets = points
	n = Npoints
	allocate(cle(n), obj(n))
	do i = 1, n
		cle(i) = int(floor((Points(i,1) - coin(1))/dl), i8) + n1*int(floor((Points(i,2) - coin(2))/dl), i8) &
			+ n2*int(floor((Points(i,3) - coin(3))/dl), i8)
		obj(i) = i
	end do

	do l = NL, NLF, -1
		il = l - NLF + 1
		call tri_tas(n, cle, obj)
		! cellules non vides : identifiants distincts, dans l'ordre croissant
		allocate(debut(n+1), ids(n))
		ng = 1 ; debut(1) = 1 ; ids(1) = cle(1)
		do i = 2, n
			if (cle(i) > cle(i-1)) then
				ng = ng + 1 ; debut(ng) = i ; ids(ng) = cle(i)
			endif
		end do
		debut(ng+1) = n + 1
		allocate(ijk(3,ng))
		do g = 1, ng
			ijk(3,g) = int(ids(g)/n2, iprec)
			ijk(2,g) = int((ids(g) - n2*ijk(3,g))/n1, iprec)
			ijk(1,g) = int(ids(g) - n2*ijk(3,g) - n1*ijk(2,g), iprec)
		end do

		Octfin(il)%NL = ng
		allocate(Octfin(il)%cell(ng))
		do g = 1, ng
			associate (c => Octfin(il)%cell(g))
			c%Nb_c = debut(g+1) - debut(g)
			c%child = obj(debut(g):debut(g+1)-1)
			c%center = coin + (real(ijk(:,g), prec) + 0.5_prec)*dl
			c%length = dl
			! voisins qui touchent
			call cellules(proche, 26, ijk(:,g), liste, nv)
			c%Nb_v = nv
			c%voisin = liste(1:nv)
			end associate
		end do
		! interactions
		if (l == NLF) then
			allocate(exclu(ng))
			do g = 1, ng
				associate (c => Octfin(il)%cell(g))
				exclu = .false. ; exclu(g) = .true. ; exclu(c%voisin) = .true.
				c%interac = pack([(m, m = 1, ng)], .not. exclu)
				c%Nb_i = size(c%interac)
				end associate
			end do
			deallocate(exclu)
		else
			do g = 1, ng
				par = merge(0, 4, mod(ijk(1,g), 2) == 0) + merge(0, 2, mod(ijk(2,g), 2) == 0) &
					+ merge(0, 1, mod(ijk(3,g), 2) == 0)
				call cellules(loin(:,:,par), 189, ijk(:,g), liste, ni)
				Octfin(il)%cell(g)%Nb_i = ni
				Octfin(il)%cell(g)%interac = liste(1:ni)
			end do
		endif

		! --- niveau superieur : objets = cellules de ce niveau, cle = cellule mere
		if (l > NLF) then
			deallocate(cle, obj)
			allocate(cle(ng), obj(ng))
			n1 = n1/2 ; n2 = n2/4
			do g = 1, ng
				cle(g) = int(ijk(1,g)/2, i8) + n1*int(ijk(2,g)/2, i8) + n2*int(ijk(3,g)/2, i8)
				obj(g) = g
			end do
			n = ng
			dl = dl*2.0_prec
		endif
		deallocate(debut, ids, ijk)
	end do
	deallocate(cle, obj)

	contains

		! cellules non vides du niveau courant aux positions ijk + dec(:,1:nd), dans
		! l'ordre des decalages ; numeros dans liste(1:nt)
		subroutine cellules(dec, nd, c, liste, nt)
		integer(iprec), intent(in) :: nd, dec(3,nd), c(3)
		integer(iprec), intent(out) :: liste(:), nt
		integer(iprec) :: q, p(3), m
		nt = 0
		do q = 1, nd
			p = c + dec(:,q)
			if (any(p < 0) .or. any(p >= n1)) cycle
			m = cherche(int(p(1), i8) + n1*int(p(2), i8) + n2*int(p(3), i8))
			if (m > 0) then
				nt = nt + 1 ; liste(nt) = m
			endif
		end do
		end subroutine cellules

		! numero de la cellule d'identifiant id (0 si vide) : dichotomie dans ids(1:ng)
		integer(iprec) function cherche(id)
		integer(i8), intent(in) :: id
		integer(iprec) :: a, b, m
		cherche = 0 ; a = 1 ; b = ng
		do while (a <= b)
			m = (a + b)/2
			if (ids(m) == id) then
				cherche = m ; return
			elseif (ids(m) < id) then
				a = m + 1
			else
				b = m - 1
			endif
		end do
		end function cherche

	end subroutine make_octree_all

!
! Decalages (di, dj, dk) des voisins, dans l'ordre lexicographique :
!   proche(:,1:26)        : {-1,0,1}^3 sauf (0,0,0)
!   loin(:,1:189,par)     : par = 4 bi + 2 bj + bk, b = 0 si l'indice de la cellule
!                           est pair (decalages de -2 a 3 sur cet axe), 1 s'il est
!                           impair (de -3 a 2), sauf les 27 de {-1,0,1}^3.
!
	subroutine decalages(proche, loin)
	integer(iprec), intent(out) :: proche(3,26), loin(3,189,0:7)
	integer(iprec) :: a, b, c, n, par, lo(3)
	n = 0
	do a = -1, 1
		do b = -1, 1
			do c = -1, 1
				if (a == 0 .and. b == 0 .and. c == 0) cycle
				n = n + 1 ; proche(:,n) = [a, b, c]
			end do
		end do
	end do
	do par = 0, 7
		lo = merge(-3, -2, [btest(par,2), btest(par,1), btest(par,0)])
		n = 0
		do a = lo(1), lo(1) + 5
			do b = lo(2), lo(2) + 5
				do c = lo(3), lo(3) + 5
					if (abs(a) <= 1 .and. abs(b) <= 1 .and. abs(c) <= 1) cycle
					n = n + 1 ; loin(:,n,par) = [a, b, c]
				end do
			end do
		end do
	end do
	end subroutine decalages

!
! Tri par tas (Williams 1964, Floyd 1964), croissant sur cle ; val suit cle.
!
	subroutine tri_tas(n, cle, val)
	integer, parameter :: i8 = selected_int_kind(18)
	integer(iprec), intent(in) :: n
	integer(i8), intent(inout) :: cle(n)
	integer(iprec), intent(inout) :: val(n)
	integer(iprec) :: i, v
	integer(i8) :: c
	do i = n/2, 1, -1
		call descendre(i, n)
	end do
	do i = n, 2, -1
		c = cle(1) ; cle(1) = cle(i) ; cle(i) = c
		v = val(1) ; val(1) = val(i) ; val(i) = v
		call descendre(1, i-1)
	end do
	contains
		subroutine descendre(haut, bas)      ! retablit le tas de racine haut dans 1..bas
		integer(iprec), intent(in) :: haut, bas
		integer(iprec) :: p, f, v0
		integer(i8) :: c0
		c0 = cle(haut) ; v0 = val(haut)
		p = haut ; f = 2*p
		do while (f <= bas)
			if (f < bas) then
				if (cle(f+1) > cle(f)) f = f + 1
			endif
			if (c0 >= cle(f)) exit
			cle(p) = cle(f) ; val(p) = val(f)
			p = f ; f = 2*f
		end do
		cle(p) = c0 ; val(p) = v0
		end subroutine descendre
	end subroutine tri_tas

!
!------------------------------------------------------------------------------------
!
! completer_octree : ajoute à l'octree de make_octree_all ce qu'il faut pour
! l'octree ADAPTATIF (2026-09-29), sans rien modifier de l'existant :
!   octree(il)%cell(i)%parent  : cellule mère (niveau il-1), 0 au niveau il=1
!   octree(il)%ipt(p)          : cellule du point p au niveau il
!   octree(il)%cell(i)%pts     : tous les points de la cellule (Nb_p) ; au
!                                niveau le plus fin, pts = child (même ordre)
!
!------------------------------------------------------------------------------------
!
	subroutine completer_octree(octree,Nlev,Npoints)
	implicit none
	integer(iprec), intent(in) :: Nlev, Npoints
	type(tree), intent(inout) :: octree(Nlev)
	integer(iprec) :: il, i, j, p, fils

	! parents
	do il=1,Nlev
		octree(il)%cell(:)%parent=0
	enddo
	do il=1,Nlev-1
		do i=1,octree(il)%NL
			do j=1,octree(il)%cell(i)%Nb_c
				octree(il+1)%cell(octree(il)%cell(i)%child(j))%parent=i
			enddo
		enddo
	enddo

	! cellule de chaque point, du niveau le plus fin vers le haut
	do il=1,Nlev
		if (allocated(octree(il)%ipt)) deallocate(octree(il)%ipt)
		allocate(octree(il)%ipt(Npoints))
		octree(il)%ipt=0
	enddo
	do i=1,octree(Nlev)%NL
		do j=1,octree(Nlev)%cell(i)%Nb_c
			octree(Nlev)%ipt(octree(Nlev)%cell(i)%child(j))=i
		enddo
	enddo
	do il=Nlev-1,1,-1
		do p=1,Npoints
			octree(il)%ipt(p)=octree(il+1)%cell(octree(il+1)%ipt(p))%parent
		enddo
	enddo

	! points de chaque cellule, de bas en haut : au niveau le plus fin pts = child
	! (MEME ORDRE), au-dessus pts = concaténation des pts des enfants dans l'ordre
	! de child. Le parcours des pts reproduit donc exactement celui des child du
	! niveau fin (init_Anear_adapt redonne les mêmes Ino/Ino1 en uniforme).
	do i=1,octree(Nlev)%NL
		octree(Nlev)%cell(i)%Nb_p=octree(Nlev)%cell(i)%Nb_c
		if (allocated(octree(Nlev)%cell(i)%pts)) deallocate(octree(Nlev)%cell(i)%pts)
		allocate(octree(Nlev)%cell(i)%pts(octree(Nlev)%cell(i)%Nb_c))
		octree(Nlev)%cell(i)%pts=octree(Nlev)%cell(i)%child(1:octree(Nlev)%cell(i)%Nb_c)
	enddo
	do il=Nlev-1,1,-1
		do i=1,octree(il)%NL
			p=0
			do j=1,octree(il)%cell(i)%Nb_c
				p=p+octree(il+1)%cell(octree(il)%cell(i)%child(j))%Nb_p
			enddo
			octree(il)%cell(i)%Nb_p=p
			if (allocated(octree(il)%cell(i)%pts)) deallocate(octree(il)%cell(i)%pts)
			allocate(octree(il)%cell(i)%pts(p))
			p=0
			do j=1,octree(il)%cell(i)%Nb_c
				fils=octree(il)%cell(i)%child(j)
				octree(il)%cell(i)%pts(p+1:p+octree(il+1)%cell(fils)%Nb_p)=octree(il+1)%cell(fils)%pts
				p=p+octree(il+1)%cell(fils)%Nb_p
			enddo
		enddo
	enddo

	end subroutine completer_octree
!
!------------------------------------------------------------------------------------
!
! niveaux_elements : niveau de chaque élément dans l'octree adaptatif.
! niv(e) = le niveau le PLUS FIN (il le plus grand) dont la cellule vaut au
! moins seuil x h(e), h(e) = taille de l'élément (plus grande arête).
! C'est la règle du NL automatique du tuyau (seuil 1,2) appliquée élément par
! élément. Si même le niveau il=1 est trop petit : niv(e)=1, compté dans nhors.
!
!------------------------------------------------------------------------------------
!
	subroutine niveaux_elements(octree,Nlev,NE,h,seuil,niv,nhors)
	implicit none
	integer(iprec), intent(in) :: Nlev, NE
	type(tree), intent(in) :: octree(Nlev)
	real(prec), intent(in) :: h(NE), seuil
	integer(iprec), intent(out) :: niv(NE), nhors
	integer(iprec) :: e, il

	nhors=0
	do e=1,NE
		niv(e)=0
		do il=Nlev,1,-1
			if (octree(il)%cell(1)%length >= seuil*h(e)) then
				niv(e)=il
				exit
			endif
		enddo
		if (niv(e)==0) then
			niv(e)=1
			nhors=nhors+1
		endif
	enddo

	end subroutine niveaux_elements
!
!------------------------------------------------------------------------------------
!
! NL_adaptatif_auto : profondeur NLmax de l'octree ADAPTATIF (2026-09-30).
! Elle n'est PAS fixee par le plus petit element : un element plus petit que
! sa cellule ne pose aucun probleme de precision (c'est le cas normal en
! uniforme), seul un element plus GRAND que sa cellule en pose. Quelques
! triangles minuscules (christophe_final : 0,03 mm dans les plis du pavillon)
! imposeraient sinon 13 niveaux, des developpements alloues pour des centaines
! de milliers de cellules inutiles et un octree faux (NL > NL_MAX_OCTREE).
! NLmax = le plus grand NL dont la cellule vaut au moins seuil x h_q, h_q =
! taille en dessous de laquelle se trouvent au plus 0,1 % des elements
! (quantile par bissection) ; borne dans [NLF+1, NL_MAX_OCTREE]. NL_force > 0
! impose la profondeur (bornee a NL_MAX_OCTREE). Les elements plus petits que
! seuil x cellule(NLmax) sont poses au niveau NLmax (niveaux_elements).
!
!------------------------------------------------------------------------------------
!
	integer function NL_adaptatif_auto(NE, h, Lc, seuil, NLF, NL_force, h_q)
	implicit none
	integer(iprec), intent(in) :: NE, NLF, NL_force
	real(prec), intent(in) :: h(NE), Lc, seuil
	real(prec), intent(out) :: h_q
	real(prec) :: a, b, m
	integer(iprec) :: it, nq

	nq = max(1, int(0.001_prec*real(NE,prec)))
	a = minval(h) ; b = maxval(h)
	do it = 1, 60
		m = sqrt(a*b)
		if (count(h < m) <= nq) then
			a = m
		else
			b = m
		endif
	end do
	h_q = a
	if (NL_force > 0) then
		NL_adaptatif_auto = min(NL_force, NL_MAX_OCTREE)
	else
		NL_adaptatif_auto = floor(log(Lc/(seuil*h_q))/log(2.0_prec))
		NL_adaptatif_auto = max(NLF+1, min(NL_adaptatif_auto, NL_MAX_OCTREE))
	endif
	end function NL_adaptatif_auto
!
!------------------------------------------------------------------------------------
!
! compte_paires_near : nombre de paires (élément, noeud) du champ proche de
! l'octree adaptatif. L'élément e est vu à SON niveau il=niv(e), dans sa
! cellule C=ipt(e) :
!   n_child  : noeuds (points > NE) de C
!   n_voisin : noeuds des voisines de C au niveau il
! Avec niv(:)=Nlev on retrouve num_K_child et num_K_voisin de
! init_Anear_part_T3/T6_CPU_v3. Nécessite completer_octree.
!
!------------------------------------------------------------------------------------
!
	subroutine compte_paires_near(octree,Nlev,NE,niv,n_child,n_voisin)
	implicit none
	integer(iprec), intent(in) :: Nlev, NE
	type(tree), intent(in) :: octree(Nlev)
	integer(iprec), intent(in) :: niv(NE)
	integer(iprec8), intent(out) :: n_child, n_voisin
	integer(iprec) :: e, il, ic, iv, j, cv

	n_child=0
	n_voisin=0
	do e=1,NE
		il=niv(e)
		ic=octree(il)%ipt(e)
		do j=1,octree(il)%cell(ic)%Nb_p
			if (octree(il)%cell(ic)%pts(j) > NE) n_child=n_child+1
		enddo
		do iv=1,octree(il)%cell(ic)%Nb_v
			cv=octree(il)%cell(ic)%voisin(iv)
			do j=1,octree(il)%cell(cv)%Nb_p
				if (octree(il)%cell(cv)%pts(j) > NE) n_voisin=n_voisin+1
			enddo
		enddo
	enddo

	end subroutine compte_paires_near
!
!------------------------------------------------------------------------------------
!
! marquer_sources_cibles : drapeaux de la FMM adaptative.
!   a_source(il,i) : un élément de niveau niv >= il est dans la cellule (le
!                    multipôle du niveau il ne reçoit que ceux-là)
!   a_cible(il,i)  : un noeud (point > NE) de la cellule est évalué à ce
!                    niveau ou plus bas ; en T0 (octree sur les NE barycentres
!                    seuls) les éléments sont eux-mêmes les cibles
!   feuille(il,i)  : les noeuds de la cellule sont évalués à ce niveau
!
! ARRÊT DES CIBLES (arret = .true.). La cellule C du niveau il est
! TERMINALE si ni elle ni ses voisines ne contiennent d'élément de niveau
! > il. Alors tout élément de niveau > il est séparé de C dès le niveau il
! (M2L à un niveau <= il, déjà dans le développement local de C) et n'est
! jamais en champ proche de ses noeuds : le développement local de C est
! complet, ses noeuds y sont évalués (feuille) et on ne descend plus
! (a_cible = .false. en dessous). Les M2L sautés en dessous sont nuls (les
! listes d'interaction des sous-cellules restent dans le voisinage de C).
! La propriété se transmet aux descendants : chaque noeud est évalué dans
! la PREMIÈRE cellule terminale. Au niveau Nlev toute cellule est terminale.
! arret = .false. : tous les noeuds sont évalués au niveau Nlev.
!
! Les sauts ne changent pas le résultat sans arrêt (ils évitent d'ajouter
! des zéros) ; avec arrêt, l'évaluation se fait plus haut (même précision,
! à la troncature près). Nécessite completer_octree.
!
!------------------------------------------------------------------------------------
!
	subroutine marquer_sources_cibles(octree,Nlev,NE,niv,arret)
	implicit none
	integer(iprec), intent(in) :: Nlev, NE
	type(tree), intent(inout) :: octree(Nlev)
	integer(iprec), intent(in) :: niv(NE)
	logical, intent(in) :: arret
	integer(iprec) :: il, p, i, iv, mx
	integer(iprec), allocatable :: nivmax(:)
	logical, allocatable :: fini(:), fini_haut(:)

	do il=1,Nlev
		octree(il)%cell(:)%a_source=.false.
		octree(il)%cell(:)%a_cible=.false.
		octree(il)%cell(:)%feuille=.false.
		do p=1,size(octree(il)%ipt)
			if (p <= NE) then
				if (niv(p) >= il) octree(il)%cell(octree(il)%ipt(p))%a_source=.true.
				! T0 (octree sur les seuls barycentres, NE points) : chaque
				! element est AUSSI une cible (collocation au barycentre)
				if (size(octree(il)%ipt) == NE) octree(il)%cell(octree(il)%ipt(p))%a_cible=.true.
			else
				octree(il)%cell(octree(il)%ipt(p))%a_cible=.true.
			endif
		enddo
	enddo
	if (.not. arret) return

	! fini(i) : la cellule i (niveau il) est une feuille ou sous une feuille
	allocate(fini_haut(0))
	do il=1,Nlev-1
		allocate(nivmax(octree(il)%NL), fini(octree(il)%NL))
		nivmax=0
		do p=1,NE
			i=octree(il)%ipt(p)
			nivmax(i)=max(nivmax(i),niv(p))
		enddo
		do i=1,octree(il)%NL
			fini(i)=.false.
			if (il > 1) fini(i)=fini_haut(octree(il)%cell(i)%parent)
			if (fini(i)) then
				octree(il)%cell(i)%a_cible=.false.      ! sous une feuille
				cycle
			endif
			if (.not. octree(il)%cell(i)%a_cible) cycle
			mx=nivmax(i)
			do iv=1,octree(il)%cell(i)%Nb_v
				mx=max(mx,nivmax(octree(il)%cell(i)%voisin(iv)))
			enddo
			if (mx <= il) then                          ! terminale
				octree(il)%cell(i)%feuille=.true.
				fini(i)=.true.
			endif
		enddo
		deallocate(fini_haut)
		call move_alloc(fini, fini_haut)
		deallocate(nivmax)
	enddo
	do i=1,octree(Nlev)%NL
		if (fini_haut(octree(Nlev)%cell(i)%parent)) octree(Nlev)%cell(i)%a_cible=.false.
	enddo
	deallocate(fini_haut)

	end subroutine marquer_sources_cibles

END MODULE module_Octree
