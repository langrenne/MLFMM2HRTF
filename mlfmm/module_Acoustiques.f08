! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: MIT
! Contains routines copyrighted by Shanjie Zhang and Jianming Jin (Computation of
! Special Functions, Wiley, 1996), used with their permission: see THIRD_PARTY.md.
module module_Acoustiques
!
! fonctions spéciales "acoustiques" : monopole, gradient, harmoniques sphériques, Bessel, Hankel...
!
use module__Precision
use module__Constants, only: pi, ji => jcomplex
implicit none


! Interfaces génériques
interface Ynm_gum
    module procedure Ynm_gum_tot
    module procedure Ynm_gum_tot_1
end interface Ynm_gum

interface funct_leg
    module procedure funct_leg_tot
    module procedure funct_leg_tot_1
end interface funct_leg

contains

    subroutine SPHJ( n, x, nm, sj, dj )

!*****************************************************************************80
!
!! SPHJ computes spherical Bessel functions jn(x) and their derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However,
!    they give permission to incorporate this routine into a user program
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    15 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, integer N, the order.
!
!    Input, real(prec) X, the argument.
!
!    Output, integer NM, the highest order computed.
!
!    Output, real(prec) SJ(0:N), the values of jn(x).
!
!    Output, real(prec) DJ(0:N), the values of jn'(x).
!
  implicit none

  integer n

  real(prec) cs
  real(prec) dj(0:n)
  real(prec) f
  real(prec) f0
  real(prec) f1
  integer k
  integer m
  integer nm
  real(prec) sa
  real(prec) sb
  real(prec) sj(0:n)
  real(prec) x

  nm = n
!  write(*,*) 'debut nm=',nm
  if ( abs ( x ) <= 1.0e-100_prec ) then
    do k = 0, n
      sj(k) = 0.0_prec
      dj(k) = 0.0_prec
    end do
    sj(0) = 1.0_prec
    dj(1) = 0.3333333333333333_prec
    return
  end if

  sj(0) = sin ( x ) / x
  sj(1) = ( sj(0) - cos ( x ) ) / x

  if ( 2 <= n ) then

    sa = sj(0)
    sb = sj(1)
    m = msta1 ( x, 200 )
!    write(*,*) 'm1=',m
    if ( m < n ) then
!        write(*,*) 'n=',n,' m=',m
      nm = m
      write(*,*) 'attention ligne 84, sph_bessel if m<n nm=',nm
    else
      m = msta2 ( x, n, 20 )
!      write(*,*) 'm2=',m
    end if

    f0 = 0.0_prec
    f1 = 1.0_prec-100.0_prec
    do k = m, 0, -1
      f = ( 2.0_prec * k + 3.0_prec ) * f1 / x - f0
      if ( k <= nm ) then
        sj(k) = f
      end if
      f0 = f1
      f1 = f
    end do

    if ( abs ( sa ) <= abs ( sb ) ) then
      cs = sb / f0
    else
      cs = sa / f
    end if

    do k = 0, nm
!        write(*,*) 'cs(',k,')=',cs,' sj=',sj(k)
      sj(k) = cs * sj(k)
    end do
    do k = nm+1,n
!        write(*,*) 'cs(',k,')=',cs,' sj=',sj(k)
      sj(k) = 0.0_prec
!      write(*,*) 'sj(-2)=',sj(k-2),'sj(-1)=',sj(k-1)
    end do

  end if
!  write(*,*) 'sj(',n,')=',sj(n)

  dj(0) = ( cos(x) - sin(x) / x ) / x
  do k = 1, nm
    dj(k) = sj(k-1) - ( k + 1.0_prec ) * sj(k) / x
  end do

  return
    end


    INTEGER FUNCTION MSTA1(X,MP)

!      ===================================================
!      Purpose: Determine the starting point for backward
!               recurrence such that the magnitude of
!               Jn(x) at that point is about 10^(-MP)
!      Input :  x     --- Argument of Jn(x)
!               MP    --- Value of magnitude
!      Output:  MSTA1 --- Starting point
!      ===================================================

        IMPLICIT real(prec) (A-H,O-Z)
        integer :: MP,IT,N0,NN,N1
        A0=abs(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO IT=1,20
           NN=N1-(N1-N0)/(1.0_prec-F0/F1)
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
           F1=F
        enddo
 20     MSTA1=NN
        RETURN
    END


    INTEGER FUNCTION MSTA2(X,N,MP)

!      ===================================================
!      Purpose: Determine the starting point for backward
!               recurrence such that all Jn(x) has MP
!               significant digits
!      Input :  x  --- Argument of Jn(x)
!               n  --- Order of Jn(x)
!               MP --- Significant digit
!      Output:  MSTA2 --- Starting point
!      ===================================================

         IMPLICIT real(prec) (A-H,O-Z)
         integer :: N,MP,IT,N0,N1,NN

        A0=abs(X)
        HMP=0.5_prec*MP
        EJN=ENVJ(N,A0)
      !  write(*,*) 'EJN=',EJN,'HMP=',HMP
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.2*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        IF(N0.eq.0) then
	  N0=1
	endif
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
       ! write(*,*) 'N1=',N1,'N0=',N0
        F1=ENVJ(N1,A0)-OBJ
        !write(*,*) 'F1=',F1,'F0=',F0,'F0/F1=',F0/F1
        DO IT=1,20
	 !  write(*,*) 'IT=',IT
           NN=N1-(N1-N0)/(1.0_prec-F0/F1)
          ! write(*,*) 'NN =',NN
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
           F1=F
        enddo
20      MSTA2=NN+10
	!write(*,*) 'MSTA2=',MSTA2
        RETURN
    END

    real(prec) FUNCTION ENVJ(N,X)
        integer :: N
        real(prec) X
        ENVJ=0.5_prec*log10(6.28_prec*N)-N*log10(1.36_prec*X/N)
        RETURN
    END

    SUBROUTINE SPHY(N,X,NM,SY,DY)

!       ======================================================
!       Purpose: Compute spherical Bessel functions yn(x) and
!                their derivatives
!       Input :  x --- Argument of yn(x) ( x ò 0 )
!                n --- Order of yn(x) ( n = 0,1,úúú )
!       Output:  SY(n) --- yn(x)
!                DY(n) --- yn'(x)
!                NM --- Highest order computed
!       ======================================================

        IMPLICIT real(prec) (A-H,O-Z)
        integer :: N,NM,K
        !DIMENSION SY(0:250),DY(0:250)
        DIMENSION SY(0:N),DY(0:N)
        NM=N
        IF (X.LT.1.0e-60_prec) THEN
           DO K=0,N
              SY(K)=-1.0e+300_prec
              DY(K)=1.0e+300_prec
           enddo
           RETURN
        ENDIF
        SY(0)=-cos(X)/X
        SY(1)=(SY(0)-sin(X))/X
        F0=SY(0)
        F1=SY(1)
        DO K=2,N
           F=(2.0_prec*K-1.0_prec)*F1/X-F0
           SY(K)=F
           IF (abs(F).GE.1.0e+300_prec) then
                write(*,*) 'attention yn très grand'
                GO TO 20
           endif
           F0=F1
           F1=F
        enddo
20      NM=K-1
           DY(0)=(sin(X)+cos(X)/X)/X
           DO K=1,NM
              DY(K)=SY(K-1)-(K+1.0_prec)*SY(K)/X
           enddo
        RETURN
    END

    Subroutine CBESSJ(z, nu, z1)
!---------------------------------------------------
!                       inf.     (-z^2/4)^k
!   Jnu(z) = (z/2)^nu x Sum  ------------------
!                       k=0  k! x Gamma(nu+k+1)
!  (nu must be >= 0).
!---------------------------------------------------
        integer,parameter :: MAXK=20
        real(prec), parameter :: ZERO=0.0_prec
        !Parameter(MAXK=20,ZERO=0.e0_prec)
        complex(prec) :: z,z1
        Integer :: k,nu
        complex(prec) :: suma,tmp
        suma = cmplx(ZERO, ZERO, kind=prec)
        do k=0, MAXK
                !calculate (-z**2/4)**k
                tmp = (-z*z/4.e0_prec)**k
                !divide by k!
                tmp = tmp / Fact(k)
                !divide by Gamma(nu+k+1)
                tmp = tmp / gamma(real(nu+k+1, prec))
                !actualize sum
                suma = suma + tmp
        end do
        !calculate (z/2)**nu
        tmp = (z/2)**nu
        !multiply (z/2)**nu by sum
        z1 = tmp*suma
        return
    End

    real(prec) Function Fact(K)
        Integer :: i,K
        real(prec)  f
        F=1.e0_prec
        do i=2, K
                F=F*real(i, prec)
        end do
        Fact=F
        return
    End

! Gamma0 supprimée : on utilise la fonction intrinsèque gamma()

!
! -------------------------------------------------------------
!
! -------------------------------------------------------------
!
! Fonctions de bessel sphérique jn par récurrence
! attention, pas super pour les valeurs x petites et non vérifié pour les autres
! -------------------------------------------------------------
!
     subroutine sph_j1(n,x,jn)
     use module_utils
     implicit none
     ! variables d'entrée
     integer :: n
     real(prec) :: x
     ! variable de sortie
     real(prec) :: jn(0:250)
     ! variables locales
     integer :: ind


     jn(0)= sin(x) / x
     jn(1)= (jn(0) - cos(x))/x
     ! récurrence
     do ind=2,n
        jn(ind)=real(2*ind-1, prec)/x * jn(ind-1) - jn(ind-2)
     end do

     end subroutine sph_j1
!
! -------------------------------------------------------------
!
! polygones de Legendre version 1
!
! -------------------------------------------------------------
!
	subroutine poly_leg(N,x,Nx,y)
	implicit none
	integer, intent(in) :: N,Nx
	real(prec), intent(in) :: x(Nx)
	real(prec), intent(out) :: y(Nx,N+1)
	integer :: i1
!	character :: answer

	y(:,1)=1
	y(:,2)=x
!	write(*,*) N
! % récurrence

	do i1=3,N+1
	    y(:,i1)=(2.0_prec*real(i1, prec)-3.0_prec)/(real(i1, prec)-1.0_prec)*x(:)*y(:,i1-1) &
	        -(real(i1, prec)-2.0_prec)/(real(i1, prec)-1.0_prec)*y(:,i1-2)
	enddo

!	write(*,*) N
!	write(*,*) y(1:2,2)
!	write(*,*) y(1:2,3)
!	read(*,*) answer

	end subroutine poly_leg
! ------------------------------------------------------------------------
!
! Harmoniques sphériques conv Gumerov
! Tous les ordres jusqu'à Nmax, soit (Nmax**2)+1, classés en ACN
! pour NB valeurs de theta, phi
!
! ------------------------------------------------------------------------
    subroutine Ynm_gum_tot(Nmax,theta,phi,NB,Ynm)
        implicit none
!
! calcule les harmoniques sphériques de nm (-n<m<n) pour des points theta phi
!
! classement des indices : ACN
!
        integer, intent(in) :: Nmax,NB
        real(prec), intent(in), dimension(NB) :: theta,phi
        complex(prec), dimension((Nmax+1)**2,NB) :: Ynm

        integer :: i1,indn,indm,ind1,ind2
        real(prec) :: costheta(NB)
        complex(prec) :: i
        real(prec), dimension((Nmax+1)**2,NB) :: FL,phi0
        integer, dimension((Nmax+1)**2) :: ACN_n,ACN_m
!        integer, dimension((Nmax+1)**2,NB) :: ACN_n0,ACN_m0
        character(len=1) :: answer
 !       real :: t1,t2


        costheta=cos(theta)
        i=ji
        ! P=zeros(n+1)  % les indices commencent à 0
        ! Ymn=zeros(n+1)
!        write(*,*) 'costheta=',costheta
!        FL=0.0_prec
!        ACN_n=0
!        ACN_m=0
!        call cpu_time(t1)
        call funct_leg_tot(Nmax,costheta,FL,NB,ACN_n,ACN_m)
!        call cpu_time(t2)
!       	write(*,*) 'passe legendre en ',t2-t1,' sec'

!        write(*,*) 'ACN_m=',ACN_m
!        write(*,*) 'exp=',exp(i*real(ACN_m, prec)*phi(2))

!        do i1=1,NB
!                ACN_n0(:,i1)=ACN_n
!                ACN_m0(:,i1)=ACN_m
!        enddo
!        do i1=1,(Nmax+1)**2
!                phi0(i1,:)=phi
!        enddo
!        ACN_n0=spread(ACN_n,2,NB)
!        ACN_m0=spread(ACN_m,2,NB)
!        phi0=spread(phi,1,(Nmax+1)**2)
!        test=exp(i*real(ACN_m0, prec)*phi0)
!        write(*,*) 'phi=',phi
!        write(*,*) 'phi0=',phi0(1,:)
!        write(*,*) 'exp=',test(:,2)
!        read(*,*) answer
!        ind=0
        do indn=0,Nmax
           do indm= 0,indn
                ind1=indn**2+indn+indm+1  ! m>=0
                ind2=indn**2+indn-indm+1  ! m < 0
                Ynm(ind1,:)=(-1.0_prec)**indm*sqrt((2.0_prec*indn+1.0_prec)*gamma(indn-abs(indm)+1.0_prec)/ &
                        (4.0_prec*pi*gamma(indn+abs(indm)+1.0_prec)))*FL(ind1,:)*exp(i*indm*phi)
                Ynm(ind2,:)=conjg(Ynm(ind1,:))
           enddo
        enddo
!        Ynm=cmplx(0.0_prec,0.0_prec, kind=prec)
!        call cpu_time(t1)
!        Ynm=(-1.0_prec)**ACN_m0*sqrt(real(2*ACN_n0+1, prec)*gamma(real(ACN_n0-abs(ACN_m0)+1, prec))) / &
!                (4.0_prec*pi*gamma(real(ACN_n0+abs(ACN_m0)+1, prec)))))*FL*exp(i*real(ACN_m0, prec)*phi0)
!        call cpu_time(t2)
!       	write(*,*) 'passe Ynm en ',t2-t1,' sec'

        ! size(Ynm)
    end subroutine Ynm_gum_tot
!
! ------------------------------------------------------------------------
!
! Harmoniques sphériques conv Gumerov
! Tous les ordres jusqu'à Nmax, soit (Nmax**2)+1, classés en ACN
! pour 1 valeur de theta, phi
!
! ------------------------------------------------------------------------
    subroutine Ynm_gum_tot_1(Nmax,theta,phi,Ynm)
        implicit none
!
! calcule les harmoniques sphériques de nm (-n<m<n) pour des points theta phi
!
! classement des indices : ACN
!
        integer, intent(in) :: Nmax
        real(prec), intent(in) :: theta,phi
        complex(prec), dimension((Nmax+1)**2) :: Ynm

        real(prec) :: costheta
        complex(prec) :: i
        real(prec), dimension((Nmax+1)**2) :: FL
        integer, dimension((Nmax+1)**2) :: ACN_n,ACN_m


        costheta=cos(theta)
        i=ji
        ! P=zeros(n+1)  % les indices commencent à 0
        ! Ymn=zeros(n+1)
!        write(*,*) 'costheta=',costheta
        call funct_leg_tot_1(Nmax,costheta,FL,ACN_n,ACN_m)

        Ynm=(-1.0_prec)**ACN_m*sqrt(real(2*ACN_n+1, prec)*gamma(real(ACN_n-abs(ACN_m)+1, prec)) / &
                (4.0_prec*pi*gamma(real(ACN_n+abs(ACN_m)+1, prec))))*FL*exp(i*real(ACN_m, prec)*phi)

    end subroutine Ynm_gum_tot_1
! ------------------------------------------------------------------------
!
! Fonctions de Legendre, toutes les valeurs classsées selon ACN format
! pour NB valeurs de x
! ------------------------------------------------------------------------
    subroutine funct_leg_tot(Nmax,x,FL,NB,ACN_n,ACN_m)
        implicit none
!
! Calcule la valeur de la fonction de Legendre associ�e
! attention m doit être positif, donc on calcul pour indice m négatif avec abs(m)
!
! pour un vecteur position x(n x 1)
!
        integer, intent(in) :: Nmax,NB
        real(prec), intent(in) :: x(NB)
        integer, dimension((Nmax+1)**2) :: ACN_n,ACN_m

        real(prec) :: FL((Nmax+1)**2,NB)

        integer :: ACN1,ACN2,m,n
        real(prec) :: Pmm(NB),Pmm1(NB),Pmm2(NB)

! Valeur initiale

        do m=0,Nmax
           do n=m,Nmax
                ACN1=n**2+n+m+1
                ACN2=n**2+n-m+1
                ACN_n(ACN1)=n
                ACN_n(ACN2)=n
                ACN_m(ACN1)=m
                ACN_m(ACN2)=-m
                if (m==n) then
                   Pmm=(-1.0_prec)**m * (1.0_prec-x**2)**(real(m, prec)/2.0_prec) * gamma(real(2*m+1, prec)) &
                        / 2.0_prec**m / gamma(real(m+1, prec))  ! modif pour m positif
                   FL(ACN1,:)=Pmm
                   FL(ACN2,:)=Pmm
                elseif (m==n-1) then
                   !Pmm=(-1.0_prec)**m*(1.0_prec-x**2)**(real(m, prec)/2.0_prec)*gamma(real(2*m+1, prec))) &
!                       /2.0_prec**m/gamma(real(m+1, prec)))
                   Pmm1=x*real(2*m+1, prec)*Pmm  ! Pmm déja calculé
                   FL(ACN1,:)=Pmm1
                   FL(ACN2,:)=Pmm1
                else
                   Pmm2=1.0_prec/real(n-m, prec) * (x*real(2*n-1, prec)*Pmm1 - real(n+m-1, prec)*Pmm)  ! Pmm et Pmm1 déjà calculés
                   Pmm=Pmm1
                   Pmm1=Pmm2
                   FL(ACN1,:)=Pmm2
                   FL(ACN2,:)=Pmm2
                endif
           enddo
        enddo

    end subroutine funct_leg_tot
!
! ------------------------------------------------------------------------
!
! Fonctions de Legendre, toutes les valeurs classsées selon ACN format
! pour 1 valeur de x
! ------------------------------------------------------------------------
    subroutine funct_leg_tot_1(Nmax,x,FL,ACN_n,ACN_m)
        implicit none
!
! Calcule la valeur de la fonction de Legendre associ�e
! attention m doit être positif, donc on calcul pour indice m négatif avec abs(m)
!
! pour un vecteur position x(n x 1)
!
        integer, intent(in) :: Nmax
        real(prec), intent(in) :: x

        real(prec) :: FL((Nmax+1)**2)
        integer, dimension((Nmax+1)**2) :: ACN_n,ACN_m

        integer :: ACN1,ACN2,m,n
        real(prec) :: Pmm,Pmm1,Pmm2

! Valeur initiale

        do m=0,Nmax
           do n=m,Nmax
                ACN1=n**2+n+m+1
                ACN2=n**2+n-m+1
                ACN_n(ACN1)=n
                ACN_n(ACN2)=n
                ACN_m(ACN1)=m
                ACN_m(ACN2)=-m
                if (m==n) then
                   Pmm=(-1.0_prec)**m * (1.0_prec-x**2)**(real(m, prec)/2.0_prec) * gamma(real(2*m+1, prec)) &
                        / 2.0_prec**m / gamma(real(m+1, prec))  ! modif pour m positif
                   FL(ACN1)=Pmm
                   FL(ACN2)=Pmm
                elseif (m==n-1) then
                   !Pmm=(-1.0_prec)**m*(1.0_prec-x**2)**(real(m, prec)/2.0_prec)*gamma(real(2*m+1, prec))) &
!                       /2.0_prec**m/gamma(real(m+1, prec)))
                   Pmm1=x*real(2*m+1, prec)*Pmm  ! Pmm déja calculé
                   FL(ACN1)=Pmm1
                   FL(ACN2)=Pmm1
                else
                   Pmm2=1.0_prec/real(n-m, prec) * (x*real(2*n-1, prec)*Pmm1 - real(n+m-1, prec)*Pmm)  ! Pmm et Pmm1 déjà calculés
                   Pmm=Pmm1
                   Pmm1=Pmm2
                   FL(ACN1)=Pmm2
                   FL(ACN2)=Pmm2
                endif
           enddo
        enddo

    end subroutine funct_leg_tot_1
!
!
! ------------------------------------------------------------------------
!
! transformation cartésiens to sphériques
!
! ------------------------------------------------------------------------
!
	subroutine sph2car(Pc,Ps,NP)
!	changement de coordonnées sphériques en coordonnées cartésiennes
	implicit none
	integer, intent(in) :: NP
	real(prec), dimension(1:NP,1:3), intent(in) :: Ps
	real(prec), dimension(1:NP,1:3), intent(out) :: Pc

	Pc(:,1)=Ps(:,1)*sin(Ps(:,2))*cos(Ps(:,3))
	Pc(:,2)=Ps(:,1)*sin(Ps(:,2))*sin(Ps(:,3))
	Pc(:,3)=Ps(:,1)*cos(Ps(:,2))

	end subroutine sph2car
!
!
! ------------------------------------------------------------------------
!
! transformation sphériques to cartésiens
!
! ------------------------------------------------------------------------
!
	subroutine car2sph(Pc,Ps,NP)
!	changement de coordonnées sphériques en coordonnées cartésiennes
	implicit none
	integer, intent(in) :: NP
	real(prec), dimension(1:NP,1:3), intent(in) :: Pc
	real(prec), dimension(1:NP,1:3), intent(out) :: Ps

        Ps(:,1) = sqrt(Pc(:,1)**2 + Pc(:,2)**2 + Pc(:,3)**2)
!        Ps(:,2) = atan(sqrt(Pc(:,1)**2 + Pc(:,2)**2) / Pc(:,3))
!        Ps(:,3) = atan(Pc(:,2)/Pc(:,1))
        Ps(:,2) = acos(Pc(:,3)/Ps(:,1))
        Ps(:,3) = atan2(Pc(:,2),Pc(:,1))
	end subroutine car2sph
!
! ------------------------------------------------------------------------
!
! champ d'un gradient de pression
!
! ------------------------------------------------------------------------
!

	subroutine champ_monopolaire_derive(dpout,k,nf,Ps,np,r0,V,nor)
!
!	Calcul le gradient d'un monopole en r0 au points Ps par rapport au vecteur nor
!
	implicit none
! Variables d'entrée
	integer, intent(in) :: nf,np
	real(prec), dimension(1:nf), intent(in) :: k
	real(prec), dimension(1:np,1:3), intent(in) :: Ps,nor
	real(prec), dimension(1,1:3), intent(in) :: r0
	real(prec) :: V

! Variable de sortie
      	complex(prec), dimension(1:np,1:nf), intent(out) :: dpout

! Variables locales
	real(prec), dimension(1:np,1) :: one
	real(prec), dimension(1:np) :: nr,norR
	real(prec), dimension(1:np,1:3) :: r
	integer :: ind
	complex(prec) :: I

	I=ji
	one(:,1)=1.0_prec
!	write(*,*) k(1)
!	write(*,*) Ps(1,:)
	r=Ps-matmul(one,r0)
!	write(*,*) r0(1,:)
!	write(*,*) Ps(1,:)
!	write(*,*) test
!	write(*,*) r(1,:)
	nr=sqrt(sum(r*r,2))
	norR=r(:,1)*nor(:,1)+r(:,2)*nor(:,2)+r(:,3)*nor(:,3);
!	write(*,*) nr(1:5)
	do ind=1,nf
		dpout(:,ind) = 1 /(4*pi)* V*(I*k(ind)*nr-1)*exp(I*k(ind)*nr) / nr**3*norR;
	enddo
	end subroutine champ_monopolaire_derive

!
! ------------------------------------------------------------------------
!
! champ de pression
!
! ------------------------------------------------------------------------
!
	subroutine champ_monopolaire(pout,k,nf,Ps,np,r0,V)
!
!	Calcul du champ d'un monopole en r0 au points Ps
!
	implicit none
! Variables d'entrée
	integer, intent(in) :: nf,np
	real(prec), dimension(1:nf), intent(in) :: k
	real(prec), dimension(1:np,1:3), intent(in) :: Ps
	real(prec), dimension(1,1:3), intent(in) :: r0
	real(prec) :: V

! Variable de sortie
      	complex(prec), dimension(1:np,1:nf), intent(out) :: pout

! Variables locales
	real(prec), dimension(1:np,1) :: one
	real(prec), dimension(1:np) :: nr
	real(prec), dimension(1:np,1:3) :: r
	integer :: ind
	complex(prec) :: I

	I=ji
	one(:,1)=1.0_prec
!	write(*,*) k(1)
!	write(*,*) Ps(1,:)
	r=Ps-matmul(one,r0)
!	write(*,*) r0(1,:)
!	write(*,*) Ps(1,:)
!	write(*,*) test
!	write(*,*) r(1,:)
	nr=sqrt(sum(r*r,2))
!	write(*,*) nr(1:5)
	do ind=1,nf
		pout(:,ind) = 1.0_prec /(4.0_prec*pi)* V*exp(I*k(ind)*nr) / nr;
	enddo
	end subroutine champ_monopolaire




end module module_Acoustiques
