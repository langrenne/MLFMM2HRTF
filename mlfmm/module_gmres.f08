! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: MIT
! =============================================================================
! module_gmres : GMRES(m) preconditionne, en communication inverse
!
! Ecrit le 2026-11-23 (C. Langrenne, avec Claude) d'apres :
!   Y. Saad, M. H. Schultz, GMRES: a generalized minimal residual algorithm for
!   solving nonsymmetric linear systems, SIAM J. Sci. Stat. Comput. 7 (1986) 856-869 ;
!   Y. Saad, Iterative Methods for Sparse Linear Systems, 2e ed., SIAM 2003, ch. 6 et 9.
!
! Algorithme : Arnoldi par Gram-Schmidt modifie, rotations de Givens complexes,
! redemarrage tous les m pas. Preconditionnement a gauche (M1), a droite (M2)
! ou des deux cotes : on resout M1^-1 A M2^-1 u = M1^-1 b, x = M2^-1 u.
!
! Communication inverse : gmres_rc rend la main a l'appelant chaque fois qu'il
! faut un produit A*v ou une application de M1^-1 ou M2^-1 ; l'appelant fait
! l'operation dans work et rappelle gmres_rc, jusqu'a irc(1) = 0.
!   irc(1) = 0 : fini (voir info)
!            1 : work(irc(4):irc(4)+n-1) <- A     * work(irc(2):irc(2)+n-1)
!            2 : work(irc(4):irc(4)+n-1) <- M1^-1 * work(irc(2):irc(2)+n-1)
!            3 : work(irc(4):irc(4)+n-1) <- M2^-1 * work(irc(2):irc(2)+n-1)
! Les produits scalaires sont faits ici (pas de demande a l'appelant).
!
! Fil d'execution (OpenMP) : TOUT l'etat entre deux appels est dans la variable
! etat (type etat_gmres) fournie par l'appelant, locale a son sous-programme.
! Aucune variable de module : plusieurs GMRES tournent en parallele sans
! THREADPRIVATE. Une variable etat neuve (ou etat%etape = 0) demarre un calcul.
!
! work (double complex, longueur lwork), a remplir par l'appelant :
!   work(1:n)      x : estimation initiale si icntl(6) = 1 ; solution en sortie
!   work(n+1:2n)   b : second membre
!   le reste est l'espace de travail : r, w, xc (n chacun), base V (n x (m+1),
!   ou n x m si icntl(8) = 1), Hessenberg + second membre reduit ((m+1)^2),
!   rotations et solution reduite (3m). Il faut
!   lwork >= m*m + m*(n+5) + 6*n + 1   (icntl(8) = 0)
!   lwork >= m*m + m*(n+5) + 5*n + 1   (icntl(8) = 1)
!
! icntl (entier, 8) : (1) unite des erreurs ; (2) unite des avertissements
!   (0 : aucun) ; (3) unite de l'historique de convergence (0 : aucun) ;
!   (4) preconditionnement 0 aucun, 1 gauche, 2 droite, 3 les deux ;
!   (5) orthogonalisation : 0 seulement (Gram-Schmidt modifie) ;
!   (6) 0 : x0 = 0, 1 : x0 fourni dans work(1:n) ; (7) iterations maximales ;
!   (8) residu au redemarrage : 0 par recurrence (sans produit A*v), 1 recalcule.
! cntl (reel, 5) : (1) tolerance ; (2:5) inutilises (doivent rester a 0).
!
! Critere d'arret : erreur arriere eta = ||M1^-1 (b - A x)|| / ||M1^-1 b||
! (M1 = I sans preconditionnement a gauche). Pendant Arnoldi on suit
! l'estimation |g(j+1)| ; quand elle passe sous la tolerance (ou au nombre
! maximal d'iterations), le residu est recalcule avec x : on s'arrete s'il
! confirme, sinon on poursuit (et le residu sera recalcule au redemarrage).
!
! info (entier, 3) : (1) 0 converge, -1 n < 1, -2 m < 1, -3 lwork trop petit,
!   -4 pas de convergence en icntl(7) iterations, -5 preconditionnement inconnu ;
!   (2) iterations effectuees ; (3) lwork necessaire.
! rinfo (reel, 2) : (1) eta (systeme preconditionne), (2) ||b - A x|| / ||b||.
! =============================================================================
module module_gmres
	implicit none
	private
	public :: etat_gmres, gmres_rc, gmres_defauts

	integer, parameter :: dp = kind(1.0d0)

	! etats de l'automate : chacun est un point de reprise (ce qu'il faut
	! faire en rentrant dans gmres_rc, en general traiter le resultat demande)
	integer, parameter :: &
		E_DEBUT  = 0,  &   ! nouveau calcul
		E_SPB    = 1,  &   ! r = M1^-1 b recu (si gauche) ; residu initial
		E_R_AV   = 2,  &   ! r = A x recu : r = b - r
		E_R_P    = 3,  &   ! r = residu : w = M1^-1 r
		E_R_FIN  = 4,  &   ! w = residu preconditionne : nouveau cycle
		E_A_PD   = 5,  &   ! debut d'un pas d'Arnoldi : w = M2^-1 V(:,j)
		E_A_AV   = 6,  &   ! w pret : r = A w
		E_A_PG   = 7,  &   ! r = A w recu : w = M1^-1 r
		E_A_GS   = 8,  &   ! w pret : orthogonalisation, rotations, test
		E_C_PD   = 9,  &   ! r = M2^-1 V y recu : xc = x + r, w = A xc
		E_C_AV   = 10, &   ! w = A xc recu : residu vrai
		E_C_PG   = 11, &   ! r = M1^-1 (b - A xc) recu : erreur arriere vraie
		E_M_DEB  = 12, &   ! base pleine : xc = V y, r = M2^-1 xc
		E_M_PD   = 13      ! r recu : x = x + r, residu de redemarrage

	type :: etat_gmres
		integer :: etape = E_DEBUT
		integer :: m = 0, nv = 0           ! taille de la base effective, colonnes de V
		integer :: iV = 0, iH = 0, ics = 0, isn = 0, iy = 0   ! adresses dans work
		integer :: j = 0                   ! pas d'Arnoldi courant (1..m)
		integer :: cycles = 0              ! redemarrages effectues
		integer :: it = 0                  ! iterations totales
		logical :: residu_vrai = .false.   ! recalculer b - A x au redemarrage
		real(dp) :: nb = 0, sPb = 0        ! ||b||, ||M1^-1 b||
		real(dp) :: eta_a = 0, eta = 0     ! erreur arriere estimee (Arnoldi), vraie
		real(dp) :: res_vrai = 0           ! ||b - A x||
	end type etat_gmres

contains

	! -------------------------------------------------------------------------
	subroutine gmres_defauts(icntl, cntl)
		integer, intent(out) :: icntl(8)
		real(dp), intent(out) :: cntl(5)
		icntl = [6, 6, 0, 0, 0, 0, 100, 0]
		cntl  = [1.0e-5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
	end subroutine gmres_defauts

	! -------------------------------------------------------------------------
	subroutine gmres_rc(n, m, lwork, work, irc, icntl, cntl, info, rinfo, etat)
		integer, intent(in) :: n, m, lwork
		complex(dp), intent(inout) :: work(*)
		integer, intent(inout) :: irc(5)
		integer, intent(in) :: icntl(8)
		real(dp), intent(in) :: cntl(5)
		integer, intent(inout) :: info(3)
		real(dp), intent(inout) :: rinfo(2)
		type(etat_gmres), intent(inout) :: etat
		integer :: besoin, ierr

		if (etat%etape == E_DEBUT) then
			ierr = icntl(1)
			info = 0
			rinfo = 0.0_dp
			irc = 0
			if (n < 1) then
				write(ierr,*) 'GMRES : n < 1'
				info(1) = -1 ; return
			end if
			if (m < 1) then
				write(ierr,*) 'GMRES : restart < 1'
				info(1) = -2 ; return
			end if
			if (icntl(4) < 0 .or. icntl(4) > 3) then
				write(ierr,*) 'GMRES : preconditionnement inconnu, icntl(4) =', icntl(4)
				info(1) = -5 ; return
			end if
			if (icntl(2) /= 0) then
				if (icntl(5) /= 0) write(icntl(2),*) 'GMRES : seul Gram-Schmidt modifie est code (icntl(5) ignore)'
				if (any(cntl(2:5) /= 0.0_dp)) write(icntl(2),*) 'GMRES : cntl(2:5) ignores'
			end if
			etat%m = min(m, n)
			etat%nv = merge(etat%m, etat%m + 1, icntl(8) == 1)
			etat%iV  = 5*n + 1
			etat%iH  = etat%iV + n*etat%nv
			etat%ics = etat%iH + (etat%m + 1)**2
			etat%isn = etat%ics + etat%m
			etat%iy  = etat%isn + etat%m
			besoin = etat%iy + etat%m - 1
			info(3) = besoin
			if (lwork < besoin) then
				write(ierr,*) 'GMRES : lwork =', lwork, ' trop petit, il faut', besoin
				info(1) = -3 ; return
			end if
		end if

		call gmres_pas(n, etat%m, etat%nv, work(1), work(n+1), work(2*n+1), work(3*n+1), &
			work(4*n+1), work(etat%iV), work(etat%iH), work(etat%ics), work(etat%isn), &
			work(etat%iy), irc, icntl, cntl, info, rinfo, etat)
	end subroutine gmres_rc

	! -------------------------------------------------------------------------
	! Automate : reprend a etat%etape, avance jusqu'a la prochaine demande a
	! l'appelant (return avec irc(1) > 0) ou jusqu'a la fin (irc(1) = 0).
	! H(1:m+1,1:m) : Hessenberg (triangulaire apres rotations) ; H(:,m+1) = g,
	! second membre reduit du probleme aux moindres carres.
	subroutine gmres_pas(n, m, nv, x, b, r, w, xc, V, H, cs, sn, y, &
			irc, icntl, cntl, info, rinfo, etat)
		integer, intent(in) :: n, m, nv
		complex(dp), intent(inout) :: x(n), b(n), r(n), w(n), xc(n)
		complex(dp), intent(inout) :: V(n,nv), H(m+1,m+1), cs(m), sn(m), y(m)
		integer, intent(inout) :: irc(5)
		integer, intent(in) :: icntl(8)
		real(dp), intent(in) :: cntl(5)
		integer, intent(inout) :: info(3)
		real(dp), intent(inout) :: rinfo(2)
		type(etat_gmres), intent(inout) :: etat

		integer, parameter :: PROD = 1, PREC_G = 2, PREC_D = 3
		integer :: AX, AR, AW, AXC          ! adresses des vecteurs dans work
		logical :: gauche, droite
		integer :: i, j, ihist, iwarn
		real(dp) :: tol, c
		complex(dp) :: h_ij, s

		AX = 1 ; AR = 2*n + 1 ; AW = 3*n + 1 ; AXC = 4*n + 1
		gauche = icntl(4) == 1 .or. icntl(4) == 3
		droite = icntl(4) == 2 .or. icntl(4) == 3
		ihist = icntl(3) ; iwarn = icntl(2)
		tol = cntl(1)
		j = etat%j

		do
			select case (etat%etape)

			! ---------------- demarrage
			case (E_DEBUT)
				etat%cycles = 0 ; etat%it = 0 ; j = 0
				etat%residu_vrai = icntl(8) == 1
				if (icntl(6) == 0) x = (0.0_dp, 0.0_dp)
				etat%nb = norme(b)
				if (ihist /= 0) then
					write(ihist,'(A,I9,A,I5,A,I2,A,ES9.2)') ' GMRES : n =', n, ', restart =', m, &
						', preconditionnement', icntl(4), ', tolerance', tol
					write(ihist,'(A)') ' iteration   eta Arnoldi    eta vraie'
				end if
				if (etat%nb == 0.0_dp) then
					if (iwarn /= 0) write(iwarn,*) 'GMRES : second membre nul, x = 0'
					x = (0.0_dp, 0.0_dp)
					call fin(0, 0.0_dp, 0.0_dp) ; return
				end if
				if (gauche) then
					call demande(PREC_G, n + 1, AR, E_SPB) ; return
				end if
				etat%etape = E_SPB

			case (E_SPB)
				etat%sPb = merge(norme(r), etat%nb, gauche)
				if (icntl(6) /= 0) then
					call demande(PROD, AX, AR, E_R_AV) ; return
				end if
				r = b                                   ! x0 = 0 : r = b sans produit
				etat%etape = E_R_P

			! ---------------- residu r = b - A x, puis w = M1^-1 r
			case (E_R_AV)
				r = b - r
				etat%etape = E_R_P

			case (E_R_P)
				if (gauche) then
					call demande(PREC_G, AR, AW, E_R_FIN) ; return
				end if
				w = r
				etat%etape = E_R_FIN

			case (E_R_FIN)
				call nouveau_cycle(norme(w))
				if (etat%etape == E_DEBUT) return       ! residu nul : fini

			! ---------------- un pas d'Arnoldi : w = M1^-1 A M2^-1 V(:,j)
			case (E_A_PD)
				if (droite) then
					call demande(PREC_D, etat%iV + (j-1)*n, AW, E_A_AV) ; return
				end if
				w = V(:,j)
				etat%etape = E_A_AV

			case (E_A_AV)
				call demande(PROD, AW, AR, E_A_PG) ; return

			case (E_A_PG)
				if (gauche) then
					call demande(PREC_G, AR, AW, E_A_GS) ; return
				end if
				w = r
				etat%etape = E_A_GS

			case (E_A_GS)
				! Gram-Schmidt modifie
				do i = 1, j
					h_ij = dot_product(V(:,i), w)
					H(i,j) = h_ij
					w = w - h_ij*V(:,i)
				end do
				H(j+1,j) = norme(w)
				if (j < nv) then
					if (real(H(j+1,j), dp) > 0.0_dp) then
						V(:,j+1) = w * (1.0_dp/real(H(j+1,j), dp))
					else
						V(:,j+1) = (0.0_dp, 0.0_dp)    ! arret heureux : g(j+1) sera nul
					end if
				end if
				! rotations precedentes, puis celle qui annule H(j+1,j)
				do i = 1, j-1
					call tourne(real(cs(i), dp), sn(i), H(i,j), H(i+1,j))
				end do
				call givens(H(j,j), H(j+1,j), c, s)
				cs(j) = cmplx(c, 0.0_dp, dp) ; sn(j) = s
				call tourne(c, s, H(j,j), H(j+1,j))
				H(j+1,j) = (0.0_dp, 0.0_dp)
				call tourne(c, s, H(j,m+1), H(j+1,m+1))
				etat%it = etat%it + 1
				etat%eta_a = abs(H(j+1,m+1)) / etat%sPb
				if (etat%eta_a <= tol .or. etat%it >= icntl(7)) then
					! controle : x candidat xc = x + M2^-1 V y, puis residu vrai
					call combinaison(j)
					if (droite) then
						call demande(PREC_D, AXC, AR, E_C_PD) ; return
					end if
					r = xc
					etat%etape = E_C_PD
				else
					if (ihist /= 0) write(ihist,'(I8,5X,ES9.2,7X,A)') etat%it, etat%eta_a, '--'
					call pas_suivant()
				end if

			! ---------------- controle du residu vrai
			case (E_C_PD)
				xc = x + r
				call demande(PROD, AXC, AW, E_C_AV) ; return

			case (E_C_AV)
				w = b - w
				etat%res_vrai = norme(w)
				if (gauche) then
					call demande(PREC_G, AW, AR, E_C_PG) ; return
				end if
				r = w
				etat%etape = E_C_PG

			case (E_C_PG)
				etat%eta = norme(r) / etat%sPb
				if (ihist /= 0) write(ihist,'(I8,5X,ES9.2,5X,ES9.2)') etat%it, etat%eta_a, etat%eta
				if (etat%eta <= tol .or. etat%it >= icntl(7)) then
					x = xc
					if (etat%eta <= tol) then
						call fin(0, etat%eta, etat%res_vrai/etat%nb)
					else
						if (iwarn /= 0) write(iwarn,*) 'GMRES : pas de convergence en', etat%it, ' iterations'
						call fin(-4, etat%eta, etat%res_vrai/etat%nb)
					end if
					return
				end if
				! l'estimation d'Arnoldi etait trop optimiste : on continue,
				! et le residu sera recalcule au prochain redemarrage
				etat%residu_vrai = .true.
				call pas_suivant()

			! ---------------- redemarrage : x = x + M2^-1 V y
			case (E_M_DEB)
				call combinaison(m)
				if (droite) then
					call demande(PREC_D, AXC, AR, E_M_PD) ; return
				end if
				r = xc
				etat%etape = E_M_PD

			case (E_M_PD)
				x = x + r
				etat%cycles = etat%cycles + 1
				if (etat%residu_vrai) then
					call demande(PROD, AX, AR, E_R_AV) ; return
				end if
				! residu par recurrence (sans produit) : w = V(:,1:m+1) Q^H (0,...,0,g(m+1))^T,
				! de norme |g(m+1)|
				c = abs(H(m+1,m+1))
				do i = m, 1, -1
					H(i,m+1) = (0.0_dp, 0.0_dp)
					call tourne(real(cs(i), dp), -sn(i), H(i,m+1), H(i+1,m+1))
				end do
				w = matmul(V(:,1:m+1), H(1:m+1,m+1))
				call nouveau_cycle(c)
				if (etat%etape == E_DEBUT) return

			case default
				write(icntl(1),*) 'GMRES : etat inconnu', etat%etape
				stop 1
			end select
		end do

	contains

		! demande a l'appelant : work(col_z) <- operation(work(col_x)), reprise en suite
		subroutine demande(code, col_x, col_z, suite)
			integer, intent(in) :: code, col_x, col_z, suite
			irc(1) = code ; irc(2) = col_x ; irc(3) = 0 ; irc(4) = col_z ; irc(5) = 0
			etat%etape = suite
			etat%j = j
		end subroutine demande

		! nouveau cycle d'Arnoldi a partir de w (residu preconditionne), ||w|| = beta
		subroutine nouveau_cycle(beta)
			real(dp), intent(in) :: beta
			if (beta == 0.0_dp) then             ! x est la solution exacte
				call fin(0, 0.0_dp, 0.0_dp)
				return
			end if
			V(:,1) = w * (1.0_dp/beta)
			H(:,m+1) = (0.0_dp, 0.0_dp)
			H(1,m+1) = cmplx(beta, 0.0_dp, dp)
			j = 1
			etat%j = 1
			etat%etape = E_A_PD
		end subroutine nouveau_cycle

		subroutine pas_suivant()
			if (j < m) then
				j = j + 1
				etat%j = j
				etat%etape = E_A_PD
			else
				etat%etape = E_M_DEB             ! base pleine : redemarrage
			end if
		end subroutine pas_suivant

		! y = H(1:k,1:k)^-1 g(1:k) (H triangulaire superieure), xc = V(:,1:k) y
		subroutine combinaison(k)
			integer, intent(in) :: k
			integer :: l
			do l = k, 1, -1
				y(l) = (H(l,m+1) - sum(H(l,l+1:k)*y(l+1:k))) / H(l,l)
			end do
			xc = (0.0_dp, 0.0_dp)
			do l = 1, k
				xc = xc + y(l)*V(:,l)
			end do
		end subroutine combinaison

		subroutine fin(code, eta, eta_brut)
			integer, intent(in) :: code
			real(dp), intent(in) :: eta, eta_brut
			info(1) = code
			info(2) = etat%it
			rinfo(1) = eta
			rinfo(2) = eta_brut
			irc = 0
			if (ihist /= 0) write(ihist,'(A,I3,A,I5,A,ES9.2)') ' GMRES : info =', code, &
				', iterations', etat%it, ', eta', eta
			etat%etape = E_DEBUT
			etat%j = 0
		end subroutine fin

	end subroutine gmres_pas

	! -------------------------------------------------------------------------
	real(dp) function norme(v)
		complex(dp), intent(in) :: v(:)
		norme = sqrt(real(dot_product(v, v), dp))
	end function norme

	! rotation de Givens complexe : c reel, s complexe tels que
	! [ c  s ; -conj(s)  c ] (a, b)^T = (r, 0)^T
	subroutine givens(a, b, c, s)
		complex(dp), intent(in) :: a, b
		real(dp), intent(out) :: c
		complex(dp), intent(out) :: s
		real(dp) :: na, nr
		na = abs(a)
		if (abs(b) == 0.0_dp) then
			c = 1.0_dp ; s = (0.0_dp, 0.0_dp)
		else if (na == 0.0_dp) then
			c = 0.0_dp ; s = (1.0_dp, 0.0_dp)
		else
			nr = hypot(na, abs(b))
			c = na/nr
			s = (a/na) * conjg(b) / nr
		end if
	end subroutine givens

	! (u, v) <- (c u + s v, c v - conj(s) u)
	subroutine tourne(c, s, u, v)
		real(dp), intent(in) :: c
		complex(dp), intent(in) :: s
		complex(dp), intent(inout) :: u, v
		complex(dp) :: t
		t = c*u + s*v
		v = c*v - conjg(s)*u
		u = t
	end subroutine tourne

end module module_gmres
