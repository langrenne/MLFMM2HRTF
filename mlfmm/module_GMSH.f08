! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: MIT
	module module_GMSH
	!
	! module GMSH pour lire les fichiers de maillages GMSH
	!
    use module__Precision
    use module_BEM
    use module_utils
    use unique_sort_mod
    implicit none

    ! variables internes
   	integer, allocatable :: E0(:,:)
	real(prec), allocatable :: Ps0(:,:)

	Type GMSH
		integer, dimension(:,:),allocatable :: E
		integer :: NE,Physical
	end type GMSH
	type(GMSH),allocatable :: mesh_gmsh(:)

	PRIVATE :: mesh_gmsh,E0,Ps0

	! ---------------------------------------------------------------
	! tampons prives du lecteur unifie lecture_gmsh (voir plus bas)
	! ---------------------------------------------------------------
	Type ENT_S                              ! une entite surface GMSH v4
		integer :: tag = 0
		integer, allocatable :: phys(:)     ! tags physiques portes
		integer, allocatable :: bnd(:)      ! courbes de bord (signees)
	end type ENT_S
	type(ENT_S), allocatable :: g_surf(:)   ! entites surfaces (v4 ; vide en v2)
	integer :: g_nsurf = 0
	real(prec), allocatable :: g_Ps(:,:)    ! points, numerotation GMSH brute
	integer, allocatable :: g_E(:,:)        ! elements surface, noeuds GMSH
	integer, allocatable :: g_phys(:)       ! tag physique de chaque element surface
	integer, allocatable :: g_etag(:)       ! tag d'ENTITE de chaque element surface (repli)
	integer, allocatable :: g_lc(:)         ! elements-lignes : tag de la courbe
	integer, allocatable :: g_l1(:), g_l2(:)!                  : les 2 noeuds
	integer :: g_NP0 = 0, g_NE0 = 0, g_nlin = 0

	PRIVATE :: ENT_S, g_surf, g_nsurf, g_Ps, g_E, g_phys, g_etag, g_lc, g_l1, g_l2
	PRIVATE :: g_NP0, g_NE0, g_nlin
	PRIVATE :: lire_v2_brut, lire_v4_brut, construit_Lim, construit_Bord
	PUBLIC  :: detecte_nb_elem
	PRIVATE :: construit_Ps_partition, nb_mots, lire_entiers, libere_tampons
	PRIVATE :: phys_de_entite, repli_par_entite, tri_croissant, type_gmsh

	contains

!--------------------------------------------------------------------
!
! call lecture_simple_v2(file,nb_file,nb_elem)
!
! Lecture simple des fichiers .msh de GMSH, version 2 (v2)
! un seul maillage global
!--------------------------------------------------------------------
!
	subroutine lecture_simple_v2(file,nb_file,nb_elem)
	implicit none

	! variables d'entrées
	integer :: nb_file,nb_elem
	character(len=nb_file) :: file

	! variables locales
	integer :: fileunit,i1,i2,nb,Nelem,ierr,NP0
	integer :: ind,tag1,tag2,tag3,tag4
	integer :: npos1,npos2,npos3,point,Npoints,reduc
	integer, allocatable :: pos1(:),pos2(:),pos3(:)
!	integer, allocatable :: E0(:,:)
!	real(prec), allocatable :: Ps0(:,:)
	character(len=10000) :: string ! chaine a lire
	character(len=15) :: Endvar

	! lecture du maillage .msh de GMSH
	fileunit=10
	open(unit=fileunit,file=file)
	! lecture en-tête
	do i1=1,4
	  read(fileunit,*) Endvar
	enddo

	! lecture points
	read(fileunit,*) NP0
	write(*,*) 'mesh: nodes read:     ',NP0
	allocate(Ps0(1:NP0,1:3))
	do i1=1,NP0
	  read(fileunit,*) ind,Ps0(i1,:)
	enddo
	read(fileunit,*) Endvar

	! lecture matrice de connectivité
	read(fileunit,*) Endvar
	read(fileunit,*) Nelem
	write(*,*) 'mesh: elements read:  ',Nelem
	allocate(E0(Nelem,nb_elem))
	NE=0
	do i1=1,Nelem
		read(10,fmt="(a)",iostat=ierr) string
		nb = 1
		do i2 = 2, len_trim(string)
			if ( (string(i2:i2)==' ') .and. (string(i2-1:i2-1) /= ' ') ) then
				nb=nb+1
			end if
		end do
		backspace(10)  ! reviens à la même ligne
		read(fileunit,*) ind,tag1,tag2,tag3,tag4,E0(i1,1:nb-5)
!	  read(fileunit,*) ind,tag1
	  if (tag1==2 .or. tag1==3) then   !T3 ou Q4
		NE=NE+1
		E0(NE,:)=E0(i1,:)
!		read(fileunit,*) tag2,tag3,tag4,E0(NE,:)
	  endif

	enddo
	read(fileunit,*) Endvar
	close(fileunit)
!
	allocate(E(NE,nb_elem))
	E=E0(1:NE,:)

	! verifie les premiers points si ils sont dans la matrice de connectivité
	allocate(pos1(NE))
	allocate(pos2(NE))
	allocate(pos3(NE))
	Npoints=0
	point=1
	reduc=0
	do while (Npoints==0)
		call FindInVector(NE,E(:,1)==point,npos1,pos1)
		call FindInVector(NE,E(:,2)==point,npos2,pos2)
		call FindInVector(NE,E(:,3)==point,npos3,pos3)
		Npoints=npos1+npos2+npos3
		if (Npoints == 0) then
			reduc=reduc+1
			point=point+1
		endif
	enddo

	NP=NP0-reduc
	allocate(Ps(NP,3))
	Ps=Ps0(reduc+1:NP0,:)
	E=E-reduc

	write(*,*) 'mesh: nodes kept:     ',NP
	write(*,*) 'mesh: elements kept:  ',NE

	end subroutine lecture_simple_v2

!--------------------------------------------------------------------
!
! call lecture_Z_v2(file,nb_file,nb_elem,nb_lim)
!
! Lecture des fichiers .msh de GMSH, version 2 (v2)
! Pour maillage avec impédance
! surface physique 1 : Piston
! surface physique 2 : extrémité ou structure rigide
! surface physique 3 : élément avec impédance
!
!--------------------------------------------------------------------
!
	subroutine lecture_Z_v2(file,nb_file,nb_elem,nb_lim)
! file : character : fichier de lecture avec extension .msh
! nb_file : nombre de caractère pour file
! nb_elem : 3 ou 4 pour T3 ou Q4
! nb_lim : nombre de partitions dans le maillage

	implicit none

	! variables d'entrées
	integer :: nb_file,nb_elem,nb_lim
	character(len=nb_file) :: file

	! variables locales
	integer :: fileunit,i1,i2,nb,Nelem,ierr,NP0
	integer :: ind,tag1,tag2,tag3,tag4
	integer :: npos1,npos2,npos3,point,Npoints,reduc
	integer, allocatable :: pos1(:),pos2(:),pos3(:)
!	integer, allocatable :: E0(:,:)
!	real(prec), allocatable :: Ps0(:,:)
	character(len=10000) :: string ! chaine a lire
	character(len=15) :: Endvar

	! lecture du maillage .msh de GMSH
	fileunit=10
	open(unit=fileunit,file=file)
	! lecture en-tête
	do i1=1,4
	  read(fileunit,*) Endvar
	enddo

	! lecture points
	read(fileunit,*) NP0
	write(*,*) 'mesh: nodes read:     ',NP0
	allocate(Ps0(1:NP0,1:3))
	do i1=1,NP0
	  read(fileunit,*) ind,Ps0(i1,:)
	enddo
	read(fileunit,*) Endvar

	! lecture matrice de connectivité
	read(fileunit,*) Endvar
	read(fileunit,*) Nelem
	write(*,*) 'mesh: elements read:  ',Nelem
	allocate(E0(Nelem,nb_elem))
	allocate(Part(Nelem))
	NE=0
	do i1=1,Nelem
		read(10,fmt="(a)",iostat=ierr) string
!		write(*,*) 'string=',string(1:len_trim(string))
		nb = 1
		do i2 = 2, len_trim(string)
			if ( (string(i2:i2)==' ') .and. (string(i2-1:i2-1) /= ' ') ) then
				nb=nb+1
			end if
		end do
		backspace(10)  ! reviens à la même ligne
		read(fileunit,*) ind,tag1,tag2,tag3,tag4,E0(i1,1:nb-5)
!	  read(fileunit,*) ind,tag1
	  if (tag1==2 .or. tag1==3) then   !T3 ou Q4
		NE=NE+1
		E0(NE,:)=E0(i1,:)
		Part(NE)=tag3   ! surface physique
!		write(*,*) 'E0(',NE,')=',E0(NE,:)
!		write(*,*) 'Part(',NE,')=',Part(NE)
!		read(fileunit,*) tag2,tag3,tag4,E0(NE,:)
	  endif

	enddo
	read(fileunit,*) Endvar
	close(fileunit)
!
	allocate(E(NE,nb_elem))
	E=E0(1:NE,:)

	! tables indices de piston et d'impédance
	allocate(pos1(NE))
	allocate(mesh(nb_lim))
	do i1=1,nb_lim
		call FindInVector(NE,Part(1:NE)==i1,npos1,pos1)
		allocate(mesh(i1)%Lim(npos1))
		mesh(i1)%NL=npos1
		mesh(i1)%Lim(:)=pos1(1:npos1)
	enddo
	! verifie les premiers points si ils sont dans la matrice de connectivité
	allocate(pos2(NE))
	allocate(pos3(NE))
	Npoints=0
	point=1
	reduc=0
	do while (Npoints==0)
		call FindInVector(NE,E(:,1)==point,npos1,pos1)
		call FindInVector(NE,E(:,2)==point,npos2,pos2)
		call FindInVector(NE,E(:,3)==point,npos3,pos3)
		Npoints=npos1+npos2+npos3
		if (Npoints == 0) then
			reduc=reduc+1
			point=point+1
		endif
	enddo

	NP=NP0-reduc
	allocate(Ps(NP,3))
	Ps=Ps0(reduc+1:NP0,:)
	E=E-reduc

	write(*,*) 'mesh: nodes kept:     ',NP
	write(*,*) 'mesh: elements kept:  ',NE


	end subroutine lecture_Z_v2
!--------------------------------------------------------------------
!
! call lecture_gmsh_v4(file,nb_file,nb_elem)
!
! Lecture avec partitions des fichiers .msh de GMSH version 4
! file : nom de fichier avec extension .msh
! nb_file : nombre de caractères du fichier
! nb_elem : type d'élément : 3 == T3, 4 == Q4
! nb_surface : nbre de surface physique =/ nbre de surface élémentaire
!
! lecture sans les limites, donc pas pour laplacien tangentiel mais Cremer
!--------------------------------------------------------------------
!
	subroutine lecture_gmsh_v4(file,nb_file,nb_elem)
	implicit none

	! variables d'entrées
	integer :: nb_file,nb_elem
	character(len=nb_file) :: file

	! variables locales
	integer :: fileunit,i1,nb,NP0,i2,Nelem,numL,ierr ! i3,lim
	integer :: tag1,tag2,tag3,tag4,tag5,tag6,tag7,tag8,boucle ! ind,nPe(nb_elem)
	integer :: blocks,npos1,npos2,npos3,point,Npoints,reduc !,npos4,npos5
	integer, allocatable :: pos1(:),pos2(:),pos3(:)!,pos4(:),pos5(:)
!	integer, allocatable :: E0(:,:) ! E1(:)
!	real(prec), allocatable :: Ps0(:,:)
	!real(prec) :: a1(1:3), a2(1:3), Nor(1:3),PQ(3,nb_elem)
	character(len=500) :: string ! chaine a lire
	character(len=15) :: Endvar
	! DEBUG : character(len=1) :: answer

	integer :: N_P,N_C,N_S,N_V ! no,no1,no2,num,init_E,numE
	integer, allocatable :: num_Pt(:) !,indices(:)
	!real(prec) :: x1,y1,z1,x2,y2,z2
	real(prec), allocatable :: ligne(:)

	!
	!type(BEM), allocatable :: mesh(:)

	! lecture du maillage .msh de GMSH
	fileunit=10
	open(unit=fileunit,file=file)
	! lecture en-tête
	do i1=1,4
	  read(fileunit,*) Endvar
	enddo
	!-----------------------------------------------------------------------------------------------
	!
	! lecture entités
    !
    !-----------------------------------------------------------------------------------------------
	read(fileunit,*) N_P,N_C,N_S,N_V
	! entités points : on ne garde rien
	do i1=1,N_P
		read(fileunit,fmt="(a)",iostat=ierr) string
	enddo
	! entités courbes : on ne garde rien
	do i1=1,N_C
		read(fileunit,fmt="(a)",iostat=ierr) string
	enddo
	! entités surfaces : on ne garde pas les limites
	allocate(mesh_gmsh(N_S))
!	allocate(mesh(1)%NB0(N_C))   ! sauvegadre dans 1 le nombre de courbes, ne sera utilisé que dans 1
!	mesh(1)%NB0(:)=0
	write(*,*) 'nbre de points = ',N_P
	write(*,*) 'nbre de courbes = ',N_C
	write(*,*) 'nbre de surfaces = ',N_S
	write(*,*) 'nbre de volumes = ',N_V
	do i1=1,N_S
		read(fileunit,fmt="(a)",iostat=ierr) string
		nb = 1
		do i2 = 2, len_trim(string)
			if ( (string(i2:i2)==' ') .and. (string(i2-1:i2-1) /= ' ') ) then
				nb=nb+1
			end if
		end do
		backspace(fileunit)  ! reviens à la même ligne
		if (allocated(ligne)) deallocate(ligne)
		allocate(ligne(nb))
		read(fileunit,*) ligne

		mesh_gmsh(i1)%Physical=ligne(9)  ! Surface physique de la surface entité i1

	enddo
!	read(*,*) answer
	! entités Volumes : on ne garde rien
	do i1=1,N_V
		read(fileunit,fmt="(a)",iostat=ierr) string
	enddo

	read(fileunit,*) Endvar	! end entités
	!
	!-----------------------------------------------------------------------------------------------
	!
	! lecture points
	!
	!-----------------------------------------------------------------------------------------------
	!
	read(fileunit,*) Endvar ! debut nodes
	read(fileunit,*) blocks,NP0,tag1,tag2   ! nombre de blocks, nbre de points, min point, max point
	write(*,*) 'mesh: nodes read:     ',NP0
!	write(*,*) 'avant lecture points NB0=',mesh(1)%NB0

	allocate(Ps0(NP0,3))

	do i1=1,blocks
		read(fileunit,*) tag1,tag2,tag3,boucle
		if (allocated(num_Pt)) deallocate(num_Pt)
		allocate(num_Pt(boucle))
		do i2=1,boucle
			read(fileunit,*) num_Pt(i2)
		enddo
		do i2=1,boucle
			read(fileunit,*) Ps0(num_Pt(i2),1),Ps0(num_Pt(i2),2),Ps0(num_Pt(i2),3)
		enddo
	enddo

	read(fileunit,*) Endvar  ! end nodes

	! lecture matrice de connectivité
	read(fileunit,*) Endvar  ! debut éléments
!	write(*,*) 'avant lecture éléments NB0=',mesh(1)%NB0

	read(fileunit,*) blocks,Nelem,tag1,tag2   ! nombre de blocks, nbre d'éléments, min élément, max élément
	write(*,*) 'mesh: elements read:  ',Nelem

	numL=0
	allocate(E0(Nelem,nb_elem+3))  ! on ajoute la dimension, tag d'entité et tag surface physique

	! boucle sur les points
	if (blocks==N_P+N_C+N_S) then

	write(*,*) 'N_P=',N_P
	do i1=1,N_P
		read(fileunit,*) tag1,tag2,tag3,boucle    ! dimension entité, tag entité, type d'éléments, nombre
		do i2=1,boucle
				read(fileunit,*) tag4,tag5
				E0(tag4,1)=tag1
				E0(tag4,2)=tag5
		enddo
	enddo

	! boucle sur les lines
	write(*,*) 'N_C=',N_C
	do i1=1,N_C
		read(fileunit,*) tag1,tag2,tag3,boucle    ! dimension entité, tag entité, type d'éléments, nombre
		do i2=1,boucle
			read(fileunit,*) tag4,tag5,tag6
			E0(tag4,1)=tag1
			E0(tag4,2)=tag2
			E0(tag4,3:4)=(/tag5,tag6/)
		enddo
	enddo

	! boucle sur les surfaces
	write(*,*) 'N_S=',N_S
	do i1=1,N_S
		read(fileunit,*) tag1,tag2,tag3,boucle    ! dimension entité, tag entité, type d'éléments, nombre
		do i2=1,boucle
			if (tag3==3) then   ! Q4
				read(fileunit,*) tag4,tag5,tag6,tag7,tag8
				E0(tag4,1)=tag1
				E0(tag4,2)=tag2
				E0(tag4,3)=mesh_gmsh(i1)%Physical
				E0(tag4,4:7)=(/tag5,tag6,tag7,tag8/)
			elseif (tag3==2) then ! T3
				read(fileunit,*) tag4,tag5,tag6,tag7
				E0(tag4,1)=tag1
				E0(tag4,2)=tag2
				E0(tag4,3)=mesh_gmsh(i1)%Physical
				!write(*,*) 'PHY_ELE(',tag4,')=',mesh_gmsh(i1)%Physical
				E0(tag4,4:6)=(/tag5,tag6,tag7/)
			endif
		enddo
	enddo

	else

	do i1=1,blocks
		read(fileunit,*) tag1,tag2,tag3,boucle    ! dimension entité, tag entité, type d'éléments, nombre
!		write(*,*) 'boucle=',boucle
		do i2=1,boucle
			if (tag1==0) then       ! Dim==0 : points
				read(fileunit,*) tag4,tag5
				E0(tag4,1)=tag1
				E0(tag4,2)=tag5
			elseif (tag1==1) then	! Dim==1 : lines
				read(fileunit,*) tag4,tag5,tag6
				E0(tag4,1)=tag1
				E0(tag4,2)=tag2
				E0(tag4,3:4)=(/tag5,tag6/)
			elseif (tag1==2) then	! Dim==2 : Surfaces
				if (tag3==3) then   ! Q4
					read(fileunit,*) tag4,tag5,tag6,tag7,tag8
					E0(tag4,1)=tag1
					E0(tag4,2)=tag2
					E0(tag4,3)=mesh_gmsh(tag2)%Physical
					E0(tag4,4:7)=(/tag5,tag6,tag7,tag8/)
				elseif (tag3==2) then ! T3
					read(fileunit,*) tag4,tag5,tag6,tag7
					E0(tag4,1)=tag1
					E0(tag4,2)=tag2
					E0(tag4,3)=mesh_gmsh(tag2)%Physical
					E0(tag4,4:6)=(/tag5,tag6,tag7/)
				endif
			endif
		enddo
	enddo
	endif
	read(fileunit,*) Endvar  ! end éléments
	close(fileunit)  ! fermeeture fichier
!	write(*,*) 'NB0=',mesh(1)%NB0
    !
	!-----------------------------------------------------------------------------------------------
	!
	! maillage total
	!
	!-----------------------------------------------------------------------------------------------
	!
	allocate(pos1(Nelem))
	call FindInVector(Nelem,E0(:,1)==2,npos1,pos1)
	allocate(E(npos1,nb_elem))
	E=E0(pos1(1:npos1),4:3+nb_elem)
	NE=npos1

	! verifie les premiers points si ils sont dans la matrice de connectivité
	deallocate(pos1)
	allocate(pos1(NE))
	allocate(pos2(NE))
	allocate(pos3(NE))
	Npoints=0
	point=1
	reduc=0
	do while (Npoints==0)
		call FindInVector(NE,E(:,1)==point,npos1,pos1)
		call FindInVector(NE,E(:,2)==point,npos2,pos2)
		call FindInVector(NE,E(:,3)==point,npos3,pos3)
		Npoints=npos1+npos2+npos3
		if (Npoints == 0) then
			reduc=reduc+1
			point=point+1
		endif
	enddo

	NP=NP0-reduc
	allocate(Ps(NP,3))
	Ps=Ps0(reduc+1:NP0,:)
	E=E-reduc

	write(*,*) 'mesh: nodes kept:     ',NP
	write(*,*) 'mesh: elements kept:  ',NE


	deallocate(pos1,pos2,pos3)
	!inversenormale
	!E=E(:,(/1,3,2/))
	end subroutine lecture_gmsh_v4


!====================================================================
! Lecteur maillage partitionne GMSH v4 AVEC limites (porte de Tuyaux)
! -> T3 nodal + Laplacien tangentiel (remplit mesh%Bord, mesh%NB).
! Geometrie (Pc, base_a1/a2, Jac0, norc) calculee inline => PAS de pre_BEM_T0.
!====================================================================
	subroutine lecture_partitions_ver4(file,nb_file,nb_elem,lim)
	implicit none

	! variables d'entrées
	integer :: nb_file,nb_elem,lim
	character(len=nb_file) :: file

	! variables locales
	integer :: fileunit,i1,i2,i3,nb,Nelem,ierr,NP0
	integer :: ind,tag1,tag2,tag3,tag4,tag5,tag6,tag7,tag8,nPe(nb_elem)
	integer :: npos1,npos2,npos3,npos4,npos5,point,Npoints,reduc
	integer, allocatable :: pos1(:),pos2(:),pos3(:),pos4(:),pos5(:)
	integer, allocatable :: E0(:,:),E1(:)
	double precision, allocatable :: Ps0(:,:)
	double precision :: a1(1:3), a2(1:3), Nor(1:3),PQ(3,nb_elem)
	character(len=500) :: string ! chaine a lire
	character(len=15) :: Endvar
	character(len=1) :: answer

	integer :: N_P,N_C,N_S,N_V,blocks,boucle,no,no1,no2,num,init_E
	integer, allocatable :: num_Pt(:),indices(:)
	double precision :: x1,y1,z1,x2,y2,z2
	double precision, allocatable :: ligne(:)

	!
	!type(BEM), allocatable :: mesh(:)

	! lecture du maillage .msh de GMSH
	fileunit=10
	open(unit=fileunit,file=file)
	! lecture en-tête
	do i1=1,4
	  read(fileunit,*) Endvar
	enddo
	! lecture entités
	read(fileunit,*) N_P,N_C,N_S,N_V
	! entités points : on ne garde rien
	do i1=1,N_P
		read(fileunit,fmt="(a)",iostat=ierr) string
	enddo
	! entités courbes : on ne garde rien
	do i1=1,N_C
		read(fileunit,fmt="(a)",iostat=ierr) string
	enddo
	! entités surfaces : on garde les limites
	allocate(mesh(N_S))
	allocate(mesh(1)%NB0(N_C))   ! sauvegadre dans 1 le nombre de courbes, ne sera utilisé que dans 1
	mesh(1)%NB0(:)=0
	write(*,*) 'nbre de maillage = ',N_S
	write(*,*) 'nbre de courbes = ',N_C
	do i1=1,N_S
		read(fileunit,fmt="(a)",iostat=ierr) string
		nb = 1
		do i2 = 2, len_trim(string)
			if ( (string(i2:i2)==' ') .and. (string(i2-1:i2-1) /= ' ') ) then
				nb=nb+1
			end if
		end do
		backspace(fileunit)  ! reviens à la même ligne
		if (allocated(ligne)) deallocate(ligne)
		allocate(ligne(nb))
		read(fileunit,*) ligne
		write(*,*) 'ligne=',ligne
		write(*,*) 'ligne(1)=',ligne(1)
		write(*,*) 'ligne(10)=',ligne(10)
		! définition des numéros de courbes limites de la surface
		allocate(mesh(int(ligne(1)))%Lim(int(ligne(10))))
		mesh(int(ligne(1)))%NL=int(ligne(10))
		do i2=1,int(ligne(10))
			write(*,*) '10+i2=',10+i2
			mesh(int(ligne(1)))%Lim(i2)=int(ligne(10+i2))
		enddo
	enddo
!	read(*,*) answer
	! entités Volumes : on ne garde rien
	do i1=1,N_V
		read(fileunit,fmt="(a)",iostat=ierr) string
	enddo

	! nombre de Surfaces
	NS=N_S
	read(fileunit,*) Endvar	! end entités

	! lecture points
	read(fileunit,*) Endvar ! debut nodes
	read(fileunit,*) blocks,NP0,tag1,tag2   ! nombre de blocks, nbre de points, min point, max point
	write(*,*) 'mesh: nodes read:     ',NP0
!	write(*,*) 'avant lecture points NB0=',mesh(1)%NB0

	allocate(Ps0(NP0,3))

	do i1=1,blocks
		read(fileunit,*) tag1,tag2,tag3,boucle
		if (allocated(num_Pt)) deallocate(num_Pt)
		allocate(num_Pt(boucle))
		do i2=1,boucle
			read(fileunit,*) num_Pt(i2)
		enddo
		do i2=1,boucle
			read(fileunit,*) Ps0(num_Pt(i2),1),Ps0(num_Pt(i2),2),Ps0(num_Pt(i2),3)
		enddo
	enddo

	read(fileunit,*) Endvar  ! end nodes

	! lecture matrice de connectivité
	read(fileunit,*) Endvar  ! debut éléments
!	write(*,*) 'avant lecture éléments NB0=',mesh(1)%NB0

	read(fileunit,*) blocks,Nelem,tag1,tag2   ! nombre de blocks, nbre d'éléments, min élément, max élément
	write(*,*) 'mesh: elements read:  ',Nelem


	allocate(E0(Nelem,nb_elem+2))  ! on ajoute la dimension et tag d'entité
	do i1=1,blocks
		read(fileunit,*) tag1,tag2,tag3,boucle    ! dimension entité, tag entité, type d'éléments, nombre
		if (tag1==1) then  ! sauvegarde les tailles des limites
			mesh(1)%NB0(tag2)=boucle   ! nombre de lines dans tag2
!			write(*,*) 'init NB0(',tag2,') = ',boucle
		endif
		do i2=1,boucle
			if (tag1==0) then       ! Dim==0 : points
				read(fileunit,*) tag4,tag5
				E0(tag4,1)=tag1
				E0(tag4,2)=tag5
			elseif (tag1==1) then	! Dim==1 : lines
				read(fileunit,*) tag4,tag5,tag6
				E0(tag4,1)=tag1
				E0(tag4,2)=tag2
				E0(tag4,3:4)=(/tag5,tag6/)
			elseif (tag1==2) then	! Dim==2 : Surfaces
				if (tag3==3) then   ! Q4
					read(fileunit,*) tag4,tag5,tag6,tag7,tag8
					E0(tag4,1)=tag1
					E0(tag4,2)=tag2
					E0(tag4,3:6)=(/tag5,tag6,tag7,tag8/)
				elseif (tag3==2) then ! T3
					read(fileunit,*) tag4,tag5,tag6,tag7
					E0(tag4,1)=tag1
					E0(tag4,2)=tag2
					E0(tag4,3:5)=(/tag5,tag6,tag7/)
				endif
			endif
		enddo
	enddo
	read(fileunit,*) Endvar  ! end éléments
	close(fileunit)  ! fermeeture fichier
!	write(*,*) 'NB0=',mesh(1)%NB0
	! crée le maillage total et les partitions
	! maillage total
	allocate(pos1(Nelem))
	call FindInVector(Nelem,E0(:,1)==2,npos1,pos1)
	allocate(E(npos1,nb_elem))
	allocate(Part(npos1))
	E=E0(pos1(1:npos1),3:2+nb_elem)
	NE=npos1

	! verifie les premiers points si ils sont dans la matrice de connectivité
	deallocate(pos1)
	allocate(pos1(NE))
	allocate(pos2(NE))
	allocate(pos3(NE))
	Npoints=0
	point=1
	reduc=0
	do while (Npoints==0)
		call FindInVector(NE,E(:,1)==point,npos1,pos1)
		call FindInVector(NE,E(:,2)==point,npos2,pos2)
		call FindInVector(NE,E(:,3)==point,npos3,pos3)
		Npoints=npos1+npos2+npos3
		if (Npoints == 0) then
			reduc=reduc+1
			point=point+1
		endif
	enddo

	NP=NP0-reduc
	allocate(Ps(NP,3))
	Ps=Ps0(reduc+1:NP0,:)
	E=E-reduc

	write(*,*) 'mesh: nodes kept:     ',NP
	write(*,*) 'mesh: elements kept:  ',NE


	deallocate(pos1,pos2,pos3,Ps0)
	!inversenormale
	!E=E(:,(/1,3,2/))

	! calcul les points centraux et normales correspondantes
	allocate(Pc(NE,3))
	if (nb_elem==3) then
		allocate(Jac0(1:NE))
		allocate(base_a1(NE,3))  !
		allocate(base_a2(NE,3))  !
		allocate(norc(NE,3))
	endif

	do i1=1,NE
	  	  nPe=E(i1,:) ! Points de l'élément
		  PQ=transpose(Ps(nPe,:)) ! Coordonnées des points de l'éléments
		  if (nb_elem==3) then
			Pc(i1,:)=matmul(PQ,(/ 1.0d0/3.0d0,1.0d0/3.0d0,1.0d0/3.0d0 /))
			a1=matmul(PQ,(/-1,1,0/))
			a2=matmul(PQ,(/-1,0,1/))
			base_a1(i1,:)=a1
			base_a2(i1,:)=a2
			Nor(1)=a1(2)*a2(3)-a1(3)*a2(2)
			Nor(2)=a1(3)*a2(1)-a1(1)*a2(3)
			Nor(3)=a1(1)*a2(2)-a1(2)*a2(1)
			Jac0(i1)=sqrt(Nor(1)**2+Nor(2)**2+Nor(3)**2)
			norc(i1,:)=Nor/Jac0(i1)
		  else  ! pour Q4, on ne calcule que le point central
			Pc(i1,:)=matmul(PQ,(/ 1.0d0/4.0d0,1.0d0/4.0d0,1.0d0/4.0d0,1.0d0/4.0d0 /))
		  endif
		  !% calcul de la base tangente naturelle en tous les points de Gauss

	enddo

	! Partitionnement du maillage pour calcul du Laplacien tangentiel
	! position du premier élément dans E0
	allocate(pos2(Nelem))
	call FindInVector(Nelem,(E0(:,1)==2),npos2,pos2)
	init_E=pos2(1)

	allocate(pos1(Nelem))
	do i1=1,NS
		call FindInVector(Nelem,(E0(:,1)==2 .and. E0(:,2)==i1),npos1,pos1)
		allocate(mesh(i1)%E(npos1,nb_elem))
		allocate(mesh(i1)%table_E(npos1))
		mesh(i1)%NE=npos1
		mesh(i1)%E=E0(pos1(1:npos1),3:2+nb_elem)-reduc
		mesh(i1)%table_E=pos1(1:npos1)-init_E+1
		Part(pos1(1:npos1)-init_E+1)=i1    ! pour connaitre la partition de E
		! tri pour Ps et revenir à une numérotation de 1 à npos1 dans E
		if (allocated(E1)) deallocate(E1)
		if (nb_elem==3) then
			allocate(E1(3*npos1))
			E1=reshape(mesh(i1)%E,(/3*npos1/))
			call unique_sort_allocatable_arg(E1, 3*npos1, indices)
		else
			allocate(E1(4*npos1))
			E1=reshape(mesh(i1)%E,(/4*npos1/))
			call unique_sort_allocatable_arg(E1, 4*npos1, indices)
		endif
!		if (allocated(indices)) deallocate(indices)
!		call unique_sort_allocatable_arg(E1, 4*npos1, indices)
!		write(*,*) 'i1=',i1,'3*pos1=',3*npos1,' indices=',size(indices)
!		write(*,*) 'mesh%E(4)=',mesh(i1)%E(4,:)
!		write(*,*) 'indices=',indices(1:5)
		mesh(i1)%NP=size(indices)
		allocate(mesh(i1)%Ps(mesh(i1)%NP,3))
		allocate(mesh(i1)%table_Ps(mesh(i1)%NP))
		mesh(i1)%Ps=Ps(indices,:)
		mesh(i1)%table_Ps=indices
		write(*,*) 'i1=',i1,'mesh%NE=',mesh(i1)%NE
!		write(*,*) 'i1=',i1,' Ps(9)=',mesh(i1)%Ps(9,:)
!		write(*,*) 'indices(9)=',indices(9),' Ps(I9)=',Ps(indices(9),:)
	end do

	if (lim==1) then
	! limites aux surfaces
	! trouve l'élément aux limites et indique le côté à intégrer
	! initialise le nombre d'éléments aux limites par surface
	!allocate(pos2(Nelem))
	allocate(pos3(Nelem))
	allocate(pos4(Nelem))
	allocate(pos5(Nelem))


	do i1=1,NS
		mesh(i1)%NB=0
		do i2=1,mesh(i1)%NL
			mesh(i1)%NB=mesh(i1)%NB+mesh(1)%NB0(abs(mesh(i1)%Lim(i2)))
!			write(*,*) 'Surface ',i1,' NB0=',mesh(1)%NB0(abs(mesh(i1)%Lim(i2)))
		end do
	end do
	write(*,*) 'Eléments limites de 1 à ',NS,' : ',mesh(:)%NB

	do i1=1,NS
		num=0
		allocate(mesh(i1)%Bord(mesh(i1)%NB,2))    ! no d'élément au bord et no du côté
		if (allocated(pos2)) deallocate(pos2)
		allocate(pos2(mesh(i1)%NE))
		do i2=1,mesh(i1)%NL
			!write(*,*) 'i1=',i1,' mesh%Lim(',i2,')=',mesh(i1)%Lim(i2)
			no=abs(mesh(i1)%Lim(i2))
			call FindInVector(Nelem,(E0(:,1)==1 .and. E0(:,2)==no),npos1,pos1)
			do i3=1,npos1
				!if (mesh(i1)%Lim(i2) >0) then
				no1=E0(pos1(i3),3)
				no2=E0(pos1(i3),4)
				!else
				!	no2=E0(pos1(i3),3)
				!	no1=E0(pos1(i3),4)
				!endif

				call FindInVector(mesh(i1)%NE,((mesh(i1)%E(:,1)==no1 .and. mesh(i1)%E(:,2)==no2) .or. &
					(mesh(i1)%E(:,1)==no2 .and. mesh(i1)%E(:,2)==no1)),npos2,pos2)
				if (npos2==0) then
					call FindInVector(mesh(i1)%NE,((mesh(i1)%E(:,2)==no1 .and. mesh(i1)%E(:,3)==no2) .or. &
					(mesh(i1)%E(:,2)==no2 .and. mesh(i1)%E(:,3)==no1)),npos3,pos3)
					if (npos3==0) then
						if (nb_elem==3) then  ! pour dernier côté T3
							call FindInVector(mesh(i1)%NE,((mesh(i1)%E(:,3)==no1 .and. mesh(i1)%E(:,1)==no2) .or. &
							(mesh(i1)%E(:,3)==no2 .and. mesh(i1)%E(:,1)==no1)),npos4,pos4)
							num=num+1
							mesh(i1)%Bord(num,:)=(/pos4(1),3/)  ! côté 3 de T3
						else
							call FindInVector(mesh(i1)%NE,((mesh(i1)%E(:,3)==no1 .and. mesh(i1)%E(:,4)==no2) .or. &
							(mesh(i1)%E(:,3)==no2 .and. mesh(i1)%E(:,4)==no1)),npos4,pos4)
							if (npos4==0) then
								call FindInVector(mesh(i1)%NE,((mesh(i1)%E(:,4)==no1 .and. mesh(i1)%E(:,1)==no2) .or. &
								(mesh(i1)%E(:,4)==no2 .and. mesh(i1)%E(:,1)==no1)),npos5,pos5)
								num=num+1
								mesh(i1)%Bord(num,:)=(/pos5(1),4/)  ! côté 4 de Q4
							else
								num=num+1
								mesh(i1)%Bord(num,:)=(/pos4(1),3/)  ! côté 3 de Q4
							endif
						endif
					else
						num=num+1
						mesh(i1)%Bord(num,:)=(/pos3(1),2/)  ! côté 2 de T3 et Q4
					end if
				else
					num=num+1
					mesh(i1)%Bord(num,:)=(/pos2(1),1/)		! coté 1 de T3 ou Q4
				endif
			end do
		end do
!	write(*,*) 'i1=',i1,'num=',num,' mesh%NB=',mesh(i1)%NB
	end do
	deallocate(pos1,pos2,pos3,pos4,pos5)
	! renumérote la matrice de connectivité avec indices de 1 à NP et aussi Bord
	do i1=1,NS
		deallocate(E0)
		allocate(E0(mesh(i1)%NE,nb_elem))
		E0=mesh(i1)%E
		if (allocated(pos1)) deallocate(pos1)
		if (allocated(pos2)) deallocate(pos2)
		allocate(pos1(mesh(i1)%NB))
		allocate(pos2(mesh(i1)%NE))
		do i2=1,mesh(i1)%NP
			call FindInVector(mesh(i1)%NE,(mesh(i1)%E(:,1)==mesh(i1)%table_Ps(i2)),npos2,pos2)
			if (npos2>0)	mesh(i1)%E(pos2(1:npos2),1)=i2
			call FindInVector(mesh(i1)%NE,(mesh(i1)%E(:,2)==mesh(i1)%table_Ps(i2)),npos2,pos2)
			if (npos2>0)	mesh(i1)%E(pos2(1:npos2),2)=i2
			call FindInVector(mesh(i1)%NE,(mesh(i1)%E(:,3)==mesh(i1)%table_Ps(i2)),npos2,pos2)
			if (npos2>0)	mesh(i1)%E(pos2(1:npos2),3)=i2
			if (nb_elem==4) then
				call FindInVector(mesh(i1)%NE,(mesh(i1)%E(:,4)==mesh(i1)%table_Ps(i2)),npos2,pos2)
				if (npos2>0)	mesh(i1)%E(pos2(1:npos2),4)=i2
			endif
		end do
		do i2=1,mesh(i1)%NE
				call FindInVector(mesh(i1)%NB,(mesh(i1)%Bord(:,1)==mesh(i1)%table_E(i2)),npos1,pos1)
				if (npos1>0) mesh(i1)%Bord(pos1(1:npos1),1)=i2
		enddo
	enddo
	if (allocated(pos2)) deallocate(pos2)
	deallocate(pos1,ligne,E0,E1,indices,num_Pt)

	else
		mesh(1)%NB=0   ! si il n'y a qu'un maillage sans bord
	endif


	end subroutine lecture_partitions_ver4


!====================================================================
!
! call lecture_gmsh(file,nb_file,nb_elem,lim)
!
! Lecteur GMSH UNIFIE v2 / v4.1 : la version est lue dans $MeshFormat
! et la lecture est aiguillee en interne (v2 si < 3.0, sinon v4).
!
! Successeur de lecture_partitions_ver4 (lecteur v4 partitionne present
! sur la branche beltrami, ou il est conserve a cote le temps de la
! re-validation du cas tuyau). Deux corrections de fond :
!
!  1) le regroupement se fait par TAG PHYSIQUE et non par entite
!     geometrique : NS = nombre de surfaces physiques. Plusieurs
!     entites GMSH peuvent porter le meme tag physique (maillage
!     partitionne) ; elles sont fusionnees dans la meme mesh(k).
!  2) les tags d'entites GMSH ne sont plus supposes contigus 1..N_S
!     (ils ne le sont pas sur un maillage partitionne : 5,9,13,20,...).
!
! En prime, les courbes limites d'un groupe physique sont fusionnees
! puis FILTREES : une courbe qui apparait deux fois dans le meme
! groupe est une frontiere interne (bord de partition, ou couture
! d'une surface de revolution), pas une limite physique -> exclue.
!
! Sorties (variables globales de module_BEM) :
!   NP, NE, NS, Ps(NP,3), E(NE,nb_elem), Part(NE)
!   mesh(1:NS) : %NE %E(indices LOCAUX) %table_E %NP %Ps %table_Ps
!                %NL %Lim   %NB %Bord          (Lim/Bord : v4 seulement)
!   geometrie calculee ici : Pc, et pour T3 Jac0, base_a1, base_a2, norc
!   => PAS d'appel a pre_BEM_T0 apres (comme lecture_partitions_ver4).
!
! lim = 1 : construit mesh%Bord = (no element local, no cote 1..nb_elem)
!           pour chaque arete de limite ; necessaire au Laplacien
!           tangentiel. lim = 0 : Bord vide, NB = 0.
!
! Les noeuds non references par un element surface retenu (points de
! construction GMSH, noeuds d'autres blocs) sont ELIMINES et la
! numerotation est compactee : pas de ligne nulle dans la matrice de
! masse du Laplacien tangentiel.
!====================================================================
	subroutine lecture_gmsh(file,nb_file,nb_elem,lim)
	implicit none

	! variables d'entrees
	integer :: nb_file,nb_elem,lim
	character(len=nb_file) :: file

	! variables locales
	integer :: fileunit,ios,i1,i2,k,cnt,nkeep,nzero,reduc
	real(prec) :: version
	character(len=500) :: string
	integer, allocatable :: map(:)
	integer :: nPe(6)
	real(prec) :: PQ(3,6), a1(3), a2(3), Nor(3)
	real(prec) :: F6(6), dF61(6), dF62(6)

	fileunit = 10
	open(unit=fileunit,file=file,status='old',action='read',iostat=ios)
	if (ios /= 0) then
		write(*,*) 'lecture_gmsh : ERREUR, impossible d ouvrir ',file
		stop 1
	endif

	read(fileunit,fmt="(a)",iostat=ios) string   ! $MeshFormat
	read(fileunit,*) version
	read(fileunit,fmt="(a)",iostat=ios) string   ! $EndMeshFormat
	write(*,*) 'lecture_gmsh : version GMSH detectee = ',version

	if (version < 3.0d0) then
		call lire_v2_brut(fileunit,nb_elem)
	else
		if (version < 4.1d0) then
			write(*,*) 'lecture_gmsh : ERREUR, seul le format v4.1 est gere pour les .msh v4'
			write(*,*) '  (les blocs $Nodes de v4.0 ont une autre disposition) ; re-exporter en 4.1'
			stop 1
		endif
		call lire_v4_brut(fileunit,nb_elem)
	endif
	close(fileunit)

	! ---------------------------------------------------------------
	! REPLI : fichier sans AUCUN groupe physique (frequent quand le .msh
	! est exporte sans "Physical Surface"). On regroupe alors par ENTITE
	! geometrique, en renumerotant les tags d'entites (non contigus en
	! general) en 1..N. C'est ce que faisait lecture_partitions_ver4.
	! On propage ces tags synthetiques sur g_surf pour que construit_Lim
	! et construit_Bord fonctionnent sans changement.
	! ---------------------------------------------------------------
	if (g_NE0 > 0) then
	if (all(g_phys(1:g_NE0) == 0)) then
		write(*,*) 'lecture_gmsh : ATTENTION, aucun groupe physique dans le fichier'
		write(*,*) '  -> repli : regroupement par ENTITE geometrique'
		call repli_par_entite()
	endif
	endif

	! ---------------------------------------------------------------
	! elements sans tag physique : ils ne peuvent appartenir a aucune
	! partition, on les ecarte du maillage global (sinon Part serait
	! indefini pour eux).
	! ---------------------------------------------------------------
	nzero = count(g_phys(1:g_NE0)==0)
	if (nzero > 0) then
		write(*,*) 'lecture_gmsh : ATTENTION, ',nzero,' element(s) surface sans tag physique -> ignore(s)'
		cnt = 0
		do i1=1,g_NE0
			if (g_phys(i1) /= 0) then
				cnt = cnt + 1
				g_E(cnt,:) = g_E(i1,:)
				g_phys(cnt) = g_phys(i1)
			endif
		enddo
		g_NE0 = cnt
	endif
	if (g_NE0 == 0) then
		write(*,*) 'lecture_gmsh : ERREUR, aucun element surface du type demande (nb_elem=',nb_elem,')'
		stop 1
	endif

	! ---------------------------------------------------------------
	! compaction des noeuds : on ne garde que ceux reellement utilises
	! par un element surface, et on renumerote 1..NP
	! ---------------------------------------------------------------
	allocate(map(g_NP0))
	map = 0
	do i1=1,g_NE0
		do i2=1,nb_elem
			map(g_E(i1,i2)) = 1
		enddo
	enddo
	nkeep = count(map==1)
	NP = nkeep
	if (allocated(Ps)) deallocate(Ps)
	allocate(Ps(NP,3))
	cnt = 0
	do i1=1,g_NP0
		if (map(i1)==1) then
			cnt = cnt + 1
			Ps(cnt,:) = g_Ps(i1,:)
			map(i1) = cnt            ! map devient la table GMSH -> local
		endif
	enddo
	! l'ancien lecteur supposait que les noeuds inutilises etaient les
	! tout premiers (Ps0(reduc+1:NP0)) ; on signale si ce n'est pas le cas
	reduc = g_NP0 - nkeep
	if (reduc > 0) then
		if (any(map(reduc+1:g_NP0)==0)) then
			write(*,*) 'lecture_gmsh : ',reduc,' noeud(s) inutilise(s) elimine(s), pas tous en tete'
			write(*,*) '  -> la numerotation globale differe de celle de lecture_partitions_ver4'
		else
			write(*,*) 'lecture_gmsh : ',reduc,' noeud(s) de tete inutilise(s) elimine(s)'
		endif
	endif

	NE = g_NE0
	if (allocated(E)) deallocate(E)
	allocate(E(NE,nb_elem))
	do i1=1,NE
		do i2=1,nb_elem
			E(i1,i2) = map(g_E(i1,i2))
		enddo
	enddo
	! memes indices pour les elements-lignes (0 = noeud non retenu)
	do i1=1,g_nlin
		g_l1(i1) = map(g_l1(i1))
		g_l2(i1) = map(g_l2(i1))
	enddo
	deallocate(map)

	write(*,*) 'Nombre de points en sortie   : ',NP
	write(*,*) 'Nombre d elements en sortie  : ',NE

	! ---------------------------------------------------------------
	! partitions par tag physique
	! ---------------------------------------------------------------
	NS = maxval(g_phys(1:NE))
	if (allocated(mesh)) deallocate(mesh)
	allocate(mesh(NS))
	if (allocated(Part)) deallocate(Part)
	allocate(Part(NE))
	Part = 0
	do k=1,NS
		cnt = count(g_phys(1:NE)==k)
		mesh(k)%NE = cnt
		allocate(mesh(k)%E(cnt,nb_elem))
		allocate(mesh(k)%table_E(cnt))
		allocate(mesh(k)%Lim(0))    ; mesh(k)%NL = 0   ! vides par defaut (cas v2)
		allocate(mesh(k)%Bord(0,2)) ; mesh(k)%NB = 0
		cnt = 0
		do i1=1,NE
			if (g_phys(i1)==k) then
				cnt = cnt + 1
				mesh(k)%E(cnt,:)    = E(i1,:)   ! indices GLOBAUX pour l'instant
				mesh(k)%table_E(cnt)= i1
				Part(i1) = k
			endif
		enddo
		if (mesh(k)%NE == 0) write(*,*) 'lecture_gmsh : ATTENTION, partition ',k,' vide'
	enddo
	write(*,*) 'Nombre de surfaces physiques : ',NS

	! ---------------------------------------------------------------
	! limites (v4 seulement : pas de $Entities en v2)
	! l'ordre compte : Bord AVANT la renumerotation locale, car il
	! matche les elements-lignes sur les indices de noeuds GLOBAUX.
	! ---------------------------------------------------------------
	if (g_nsurf > 0) then
		call construit_Lim()
		if (lim == 1) call construit_Bord(nb_elem)
	else if (lim == 1) then
		write(*,*) 'lecture_gmsh : ATTENTION, lim=1 demande mais pas de $Entities (v2) -> Bord vide'
	endif

	call construit_Ps_partition(nb_elem)

	do k=1,NS
		write(*,'(A,I3,A,I8,A,I8,A,I6,A,I4)') ' partition phys=',k, &
			' : NE=',mesh(k)%NE,' NP=',mesh(k)%NP,' NB=',mesh(k)%NB,' NL=',mesh(k)%NL
	enddo

	! ---------------------------------------------------------------
	! geometrie des elements (comme lecture_partitions_ver4)
	! ---------------------------------------------------------------
	if (allocated(Pc)) deallocate(Pc)
	allocate(Pc(NE,3))
	if (nb_elem==3 .or. nb_elem==6) then
		if (allocated(Jac0))    deallocate(Jac0)
		if (allocated(base_a1)) deallocate(base_a1)
		if (allocated(base_a2)) deallocate(base_a2)
		if (allocated(norc))    deallocate(norc)
		allocate(Jac0(NE),base_a1(NE,3),base_a2(NE,3),norc(NE,3))
	endif
	do i1=1,NE
		nPe(1:nb_elem) = E(i1,:)
		PQ(:,1:nb_elem) = transpose(Ps(nPe(1:nb_elem),:))
		if (nb_elem==3) then
			Pc(i1,:) = matmul(PQ(:,1:3),(/ 1.0d0/3.0d0,1.0d0/3.0d0,1.0d0/3.0d0 /))
			a1 = matmul(PQ(:,1:3),(/-1.0d0,1.0d0,0.0d0/))
			a2 = matmul(PQ(:,1:3),(/-1.0d0,0.0d0,1.0d0/))
			base_a1(i1,:) = a1
			base_a2(i1,:) = a2
			Nor(1) = a1(2)*a2(3)-a1(3)*a2(2)
			Nor(2) = a1(3)*a2(1)-a1(1)*a2(3)
			Nor(3) = a1(1)*a2(2)-a1(2)*a2(1)
			Jac0(i1) = sqrt(Nor(1)**2+Nor(2)**2+Nor(3)**2)
			norc(i1,:) = Nor/Jac0(i1)
		elseif (nb_elem==6) then
			! T6 : la geometrie n'est PAS constante sur l'element. Ce qui est
			! stocke ici est la valeur au CENTRE parametrique (1/3,1/3), utile
			! a l'octree et aux diagnostics. Toute integration T6 doit
			! recalculer base_a1/base_a2/Jac au point de Gauss (module_BEM_T6).
			call FORM_T6((/1.0d0/3.0d0, 1.0d0/3.0d0/), 1, F6, dF61, dF62)
			Pc(i1,:) = matmul(PQ(:,1:6), F6)
			a1 = matmul(PQ(:,1:6), dF61)
			a2 = matmul(PQ(:,1:6), dF62)
			base_a1(i1,:) = a1
			base_a2(i1,:) = a2
			Nor(1) = a1(2)*a2(3)-a1(3)*a2(2)
			Nor(2) = a1(3)*a2(1)-a1(1)*a2(3)
			Nor(3) = a1(1)*a2(2)-a1(2)*a2(1)
			Jac0(i1) = sqrt(Nor(1)**2+Nor(2)**2+Nor(3)**2)
			norc(i1,:) = Nor/Jac0(i1)
		else    ! Q4 : seulement le point central
			Pc(i1,:) = matmul(PQ(:,1:4),(/ 0.25d0,0.25d0,0.25d0,0.25d0 /))
		endif
	enddo

	call libere_tampons()

	end subroutine lecture_gmsh

!--------------------------------------------------------------------
! GMSH v2 : le tag physique est porte directement par chaque element,
! il n'y a pas de $Entities (donc pas de Lim/Bord possible).
!--------------------------------------------------------------------
	subroutine lire_v2_brut(fileunit,nb_elem)
	implicit none
	integer, intent(in) :: fileunit,nb_elem

	character(len=2000) :: string
	integer :: ios,i1,n,ntags,ty,nn,NPlus,Nelem0
	integer :: v(64)
	integer, allocatable :: tagn(:)
	real(prec), allocatable :: xyz(:,:)

	! $PhysicalNames optionnel
	read(fileunit,fmt="(a)",iostat=ios) string
	if (index(string,'$PhysicalNames') > 0) then
		read(fileunit,*) n
		do i1=1,n
			read(fileunit,fmt="(a)",iostat=ios) string
		enddo
		read(fileunit,fmt="(a)",iostat=ios) string   ! $EndPhysicalNames
		read(fileunit,fmt="(a)",iostat=ios) string   ! $Nodes
	endif
	! ici string = $Nodes

	read(fileunit,*) NPlus
	write(*,*) 'Nombre de points en entree   : ',NPlus
	allocate(tagn(NPlus),xyz(NPlus,3))
	do i1=1,NPlus
		read(fileunit,*) tagn(i1),xyz(i1,1),xyz(i1,2),xyz(i1,3)
	enddo
	read(fileunit,fmt="(a)",iostat=ios) string   ! $EndNodes

	g_NP0 = maxval(tagn)
	allocate(g_Ps(g_NP0,3))
	g_Ps = 0.0d0
	do i1=1,NPlus
		g_Ps(tagn(i1),:) = xyz(i1,:)
	enddo
	deallocate(tagn,xyz)

	read(fileunit,fmt="(a)",iostat=ios) string   ! $Elements
	read(fileunit,*) Nelem0
	write(*,*) 'Nombre d elements en entree  : ',Nelem0
	allocate(g_E(Nelem0,nb_elem),g_phys(Nelem0),g_etag(Nelem0))
	allocate(g_lc(Nelem0),g_l1(Nelem0),g_l2(Nelem0))
	g_NE0 = 0
	g_nlin = 0
	do i1=1,Nelem0
		read(fileunit,fmt="(a)",iostat=ios) string
		call lire_entiers(string,64,v,n)
		! format v2 : tag type ntags <ntags tags> <noeuds>
		! les 2 premiers tags sont : tag physique, tag entite geometrique
		ty    = v(2)
		ntags = v(3)
		nn    = n - 3 - ntags
		if (ty == type_gmsh(nb_elem)) then   ! 2=T3, 3=Q4, 9=T6
			g_NE0 = g_NE0 + 1
			if (nb_elem == 6) then           ! GMSH A B C mAB mBC mCA -> Matlab A mAB B mBC C mCA
				g_E(g_NE0,1:6) = v((/3+ntags+1,3+ntags+4,3+ntags+2,3+ntags+5,3+ntags+3,3+ntags+6/))
			else
				g_E(g_NE0,1:nb_elem) = v(3+ntags+1 : 3+ntags+nb_elem)
			endif
			if (ntags >= 1) then
				g_phys(g_NE0) = v(4)
			else
				g_phys(g_NE0) = 0
			endif
			if (ntags >= 2) then
				g_etag(g_NE0) = v(5)
			else
				g_etag(g_NE0) = 0
			endif
		elseif (ty == 1 .and. nn == 2) then  ! element-ligne
			g_nlin = g_nlin + 1
			if (ntags >= 2) then
				g_lc(g_nlin) = v(5)          ! tag de l'entite courbe
			else
				g_lc(g_nlin) = 0
			endif
			g_l1(g_nlin) = v(3+ntags+1)
			g_l2(g_nlin) = v(3+ntags+2)
		endif
	enddo
	read(fileunit,fmt="(a)",iostat=ios) string   ! $EndElements

	g_nsurf = 0    ! pas de $Entities en v2 -> pas de Lim/Bord

	end subroutine lire_v2_brut

!--------------------------------------------------------------------
! GMSH v4.1 : les elements referencent une entite geometrique ; il
! faut la table entite -> tag(s) physique(s) lue dans $Entities.
!--------------------------------------------------------------------
	subroutine lire_v4_brut(fileunit,nb_elem)
	implicit none
	integer, intent(in) :: fileunit,nb_elem

	character(len=4000) :: string
	integer :: ios,i1,i2,nb,N_P,N_C,N_S,N_V
	integer :: blocks,boucle,tag1,tag2,tag3,tag4,tag5,tag6,tag7,tag8,tag9,tag10
	integer :: numPhys,numBnd,Nelem0,mintag,maxtag,nnodes
	real(prec), allocatable :: ligne(:)
	integer, allocatable :: num_Pt(:)

	! $PhysicalNames optionnel
	read(fileunit,fmt="(a)",iostat=ios) string
	if (index(string,'$PhysicalNames') > 0) then
		read(fileunit,*) nb
		do i1=1,nb
			read(fileunit,fmt="(a)",iostat=ios) string
		enddo
		read(fileunit,fmt="(a)",iostat=ios) string   ! $EndPhysicalNames
		read(fileunit,fmt="(a)",iostat=ios) string   ! $Entities
	endif
	! ici string = $Entities

	read(fileunit,*) N_P,N_C,N_S,N_V
	write(*,*) 'entites GMSH : ',N_P,' points ',N_C,' courbes ',N_S,' surfaces ',N_V,' volumes'

	do i1=1,N_P          ! entites points : rien a garder
		read(fileunit,fmt="(a)",iostat=ios) string
	enddo
	do i1=1,N_C          ! entites courbes : rien a garder (vues via les surfaces)
		read(fileunit,fmt="(a)",iostat=ios) string
	enddo

	g_nsurf = N_S
	allocate(g_surf(N_S))
	do i1=1,N_S
		! format : tag  xmin ymin zmin xmax ymax zmax
		!          numPhysicalTags <tags>  numBoundingCurves <courbes signees>
		read(fileunit,fmt="(a)",iostat=ios) string
		nb = nb_mots(string)
		backspace(fileunit)
		if (allocated(ligne)) deallocate(ligne)
		allocate(ligne(nb))
		read(fileunit,*) ligne

		g_surf(i1)%tag = int(ligne(1))
		numPhys = int(ligne(8))
		allocate(g_surf(i1)%phys(numPhys))
		if (numPhys > 0) g_surf(i1)%phys = int(ligne(9:8+numPhys))
		numBnd = int(ligne(8+numPhys+1))
		allocate(g_surf(i1)%bnd(numBnd))
		if (numBnd > 0) g_surf(i1)%bnd = int(ligne(8+numPhys+2 : 8+numPhys+1+numBnd))
		if (numPhys > 1) write(*,*) 'lecture_gmsh : ATTENTION, entite surface tag=',g_surf(i1)%tag, &
			' porte ',numPhys,' tags physiques -> seul le 1er est utilise'
	enddo

	do i1=1,N_V          ! entites volumes : rien a garder
		read(fileunit,fmt="(a)",iostat=ios) string
	enddo
	read(fileunit,fmt="(a)",iostat=ios) string   ! $EndEntities

	! ----- $Nodes -----
	read(fileunit,fmt="(a)",iostat=ios) string   ! $Nodes
	read(fileunit,*) blocks,nnodes,mintag,maxtag
	write(*,*) 'Nombre de points en entree   : ',nnodes
	g_NP0 = maxtag
	allocate(g_Ps(g_NP0,3))
	g_Ps = 0.0d0
	do i1=1,blocks
		read(fileunit,*) tag1,tag2,tag3,boucle
		if (allocated(num_Pt)) deallocate(num_Pt)
		allocate(num_Pt(boucle))
		do i2=1,boucle
			read(fileunit,*) num_Pt(i2)
		enddo
		do i2=1,boucle
			read(fileunit,*) g_Ps(num_Pt(i2),1),g_Ps(num_Pt(i2),2),g_Ps(num_Pt(i2),3)
		enddo
	enddo
	read(fileunit,fmt="(a)",iostat=ios) string   ! $EndNodes

	! ----- $Elements -----
	read(fileunit,fmt="(a)",iostat=ios) string   ! $Elements
	read(fileunit,*) blocks,Nelem0,mintag,maxtag
	write(*,*) 'Nombre d elements en entree  : ',Nelem0
	allocate(g_E(Nelem0,nb_elem),g_phys(Nelem0),g_etag(Nelem0))
	allocate(g_lc(Nelem0),g_l1(Nelem0),g_l2(Nelem0))
	g_NE0 = 0
	g_nlin = 0
	do i1=1,blocks
		read(fileunit,*) tag1,tag2,tag3,boucle   ! entityDim entityTag elementType numInBlock
		do i2=1,boucle
			if (tag3==15) then                          ! point : ignore
				read(fileunit,*) tag4,tag5
			elseif (tag3==1) then                       ! ligne
				read(fileunit,*) tag4,tag5,tag6
				g_nlin = g_nlin + 1
				g_lc(g_nlin) = tag2
				g_l1(g_nlin) = tag5
				g_l2(g_nlin) = tag6
			elseif (tag3==2 .and. nb_elem==3) then      ! T3
				read(fileunit,*) tag4,tag5,tag6,tag7
				g_NE0 = g_NE0 + 1
				g_E(g_NE0,1:3) = (/tag5,tag6,tag7/)
				g_phys(g_NE0)  = phys_de_entite(tag2)
				g_etag(g_NE0)  = tag2
			elseif (tag3==3 .and. nb_elem==4) then      ! Q4
				read(fileunit,*) tag4,tag5,tag6,tag7,tag8
				g_NE0 = g_NE0 + 1
				g_E(g_NE0,1:4) = (/tag5,tag6,tag7,tag8/)
				g_phys(g_NE0)  = phys_de_entite(tag2)
				g_etag(g_NE0)  = tag2
			elseif (tag3==9 .and. nb_elem==6) then      ! T6 quadratique
				! GMSH type 9 ordonne : A B C mAB mBC mCA
				! On stocke dans la convention MATLAB (Forme6/d1Forme), qui
				! tourne dans le sens trigo : A mAB B mBC C mCA
				!   -> E_matlab = E_gmsh(/1,4,2,5,3,6/)
				read(fileunit,*) tag4,tag5,tag6,tag7,tag8,tag9,tag10
				g_NE0 = g_NE0 + 1
				g_E(g_NE0,1:6) = (/tag5,tag8,tag6,tag9,tag7,tag10/)
				g_phys(g_NE0)  = phys_de_entite(tag2)
				g_etag(g_NE0)  = tag2
			elseif (tag3==8) then                       ! ligne a 3 noeuds
				! on ne garde que les 2 extremites : Lim/Bord raisonnent sur les aretes
				read(fileunit,*) tag4,tag5,tag6,tag7
				g_nlin = g_nlin + 1
				g_lc(g_nlin) = tag2
				g_l1(g_nlin) = tag5
				g_l2(g_nlin) = tag6
			else                                        ! type non gere
				read(fileunit,fmt="(a)",iostat=ios) string
			endif
		enddo
	enddo
	read(fileunit,fmt="(a)",iostat=ios) string   ! $EndElements

	if (allocated(ligne))  deallocate(ligne)
	if (allocated(num_Pt)) deallocate(num_Pt)

	end subroutine lire_v4_brut

!--------------------------------------------------------------------
! tag physique porte par l'entite surface de tag etag (0 si aucun)
!--------------------------------------------------------------------
	integer function phys_de_entite(etag) result(p)
	implicit none
	integer, intent(in) :: etag
	integer :: i1
	p = 0
	do i1=1,g_nsurf
		if (g_surf(i1)%tag == etag) then
			if (size(g_surf(i1)%phys) > 0) p = g_surf(i1)%phys(1)
			return
		endif
	enddo
	write(*,*) 'lecture_gmsh : ATTENTION, aucune entite surface de tag ',etag
	end function phys_de_entite

!--------------------------------------------------------------------
! Repli quand le fichier ne porte aucun groupe physique : les tags
! d'entites geometriques presents sont tries, dedoublonnes et
! renumerotes en 1..N, puis utilises comme tags physiques -- a la fois
! sur les elements (g_phys) et sur les entites surfaces (g_surf%phys),
! pour que la suite de la chaine soit inchangee.
!--------------------------------------------------------------------
	subroutine repli_par_entite()
	implicit none
	integer :: i1, i2, nent
	integer, allocatable :: etab(:)

	! liste triee sans doublon des tags d'entites portant des elements
	allocate(etab(g_NE0))
	nent = 0
	do i1 = 1, g_NE0
		if (.not. any(etab(1:nent) == g_etag(i1))) then
			nent = nent + 1
			etab(nent) = g_etag(i1)
		endif
	enddo
	call tri_croissant(etab, nent)

	do i1 = 1, g_NE0
		do i2 = 1, nent
			if (etab(i2) == g_etag(i1)) then
				g_phys(i1) = i2
				exit
			endif
		enddo
	enddo

	! memes tags synthetiques sur les entites surfaces (pour Lim/Bord)
	do i1 = 1, g_nsurf
		if (allocated(g_surf(i1)%phys)) deallocate(g_surf(i1)%phys)
		allocate(g_surf(i1)%phys(1))
		g_surf(i1)%phys(1) = 0                 ! entite sans element : ignoree
		do i2 = 1, nent
			if (etab(i2) == g_surf(i1)%tag) then
				g_surf(i1)%phys(1) = i2
				exit
			endif
		enddo
	enddo

	write(*,*) '  -> ', nent, ' entite(s) geometrique(s) -> autant de partitions'
	deallocate(etab)
	end subroutine repli_par_entite

!--------------------------------------------------------------------
	subroutine tri_croissant(v, n)
	implicit none
	integer, intent(inout) :: v(:)
	integer, intent(in) :: n
	integer :: i1, i2, t
	do i1 = 2, n
		t = v(i1)
		i2 = i1 - 1
		do while (i2 >= 1)
			if (v(i2) <= t) exit
			v(i2+1) = v(i2)
			i2 = i2 - 1
		enddo
		v(i2+1) = t
	enddo
	end subroutine tri_croissant

!--------------------------------------------------------------------
! mesh(k)%Lim : reunion des courbes de bord de toutes les entites du
! groupe physique k, SANS les courbes qui y apparaissent deux fois
! (frontiere interne de partition ou couture de revolution).
!--------------------------------------------------------------------
	subroutine construit_Lim()
	implicit none
	integer :: k,i1,i2,nmax,ntot,c
	integer, allocatable :: pool(:),keep(:)

	do k=1,NS
		nmax = 0
		do i1=1,g_nsurf
			if (size(g_surf(i1)%phys) > 0) then
				if (g_surf(i1)%phys(1)==k) nmax = nmax + size(g_surf(i1)%bnd)
			endif
		enddo
		if (allocated(pool)) deallocate(pool)
		allocate(pool(nmax))
		ntot = 0
		do i1=1,g_nsurf
			if (size(g_surf(i1)%phys) > 0) then
				if (g_surf(i1)%phys(1)==k) then
					do i2=1,size(g_surf(i1)%bnd)
						ntot = ntot + 1
						pool(ntot) = abs(g_surf(i1)%bnd(i2))
					enddo
				endif
			endif
		enddo

		if (allocated(keep)) deallocate(keep)
		allocate(keep(ntot))
		mesh(k)%NL = 0
		do i1=1,ntot
			c = count(pool(1:ntot)==pool(i1))
			if (c==1) then
				mesh(k)%NL = mesh(k)%NL + 1
				keep(mesh(k)%NL) = pool(i1)
			endif
		enddo
		if (allocated(mesh(k)%Lim)) deallocate(mesh(k)%Lim)
		allocate(mesh(k)%Lim(mesh(k)%NL))
		mesh(k)%Lim = keep(1:mesh(k)%NL)
	enddo
	if (allocated(pool)) deallocate(pool)
	if (allocated(keep)) deallocate(keep)

	end subroutine construit_Lim

!--------------------------------------------------------------------
! mesh(k)%Bord(:,1:2) = (no element LOCAL, no cote 1..nb_elem) pour
! chaque element-ligne porte par une courbe de mesh(k)%Lim.
! A appeler AVANT construit_Ps_partition : mesh%E est encore en
! indices de noeuds globaux, comme g_l1/g_l2.
!--------------------------------------------------------------------
!----------------------------------------------------------------------
! Detecte le nombre de noeuds par element de SURFACE d'un fichier .msh :
! 3 (type GMSH 2, triangle lineaire) ou 6 (type GMSH 9, triangle
! quadratique). Sert au solveur pour choisir la chaine T3 ou T6 sans
! ligne de configuration supplementaire -- une ligne de config en plus,
! c'est un decalage de lecture en puissance dans les 7 autres fichiers.
!
! Renvoie 0 si aucun element de surface n'est trouve.
!----------------------------------------------------------------------
	integer function detecte_nb_elem(file,nb_file) result(nne)
	implicit none
	integer, intent(in) :: nb_file
	character(len=nb_file), intent(in) :: file
	integer :: fileunit,ios,i1,i2,nel,nblocs,typ,nb,dim0,tag0
	real(prec) :: version
	character(len=500) :: string
	integer :: v(10)

	nne = 0
	fileunit = 11
	open(unit=fileunit,file=file,status='old',action='read',iostat=ios)
	if (ios /= 0) then
		write(*,*) 'detecte_nb_elem : ERREUR, impossible d ouvrir ',file
		stop 1
	endif
	read(fileunit,fmt="(a)",iostat=ios) string   ! $MeshFormat
	read(fileunit,*) version
	read(fileunit,fmt="(a)",iostat=ios) string   ! $EndMeshFormat

	do
		read(fileunit,fmt="(a)",iostat=ios) string
		if (ios /= 0) exit
		if (trim(string) == '$Elements') exit
	end do
	if (ios /= 0) then
		close(fileunit) ; return
	endif

	if (version < 3.0_prec) then
		! v2 : une ligne d'entete (nombre d'elements) puis
		!      "id type ntags ..." -> le TYPE est le 2e entier
		read(fileunit,*) nel
		do i1 = 1, nel
			read(fileunit,*,iostat=ios) v(1), v(2)
			if (ios /= 0) exit
			if (v(2) == 2) then ; nne = 3 ; exit ; endif
			if (v(2) == 9) then ; nne = 6 ; exit ; endif
		end do
	else
		! v4 : "nblocs nelem min max" puis, par bloc,
		!      "dim tag type nb" suivi de nb lignes
		read(fileunit,*) nblocs
		do i1 = 1, nblocs
			read(fileunit,*,iostat=ios) dim0, tag0, typ, nb
			if (ios /= 0) exit
			if (typ == 2) then ; nne = 3 ; exit ; endif
			if (typ == 9) then ; nne = 6 ; exit ; endif
			do i2 = 1, nb
				read(fileunit,fmt="(a)",iostat=ios) string
			end do
		end do
	endif
	close(fileunit)

	end function detecte_nb_elem

	subroutine construit_Bord(nb_elem)
	implicit none
	integer, intent(in) :: nb_elem
	integer :: k,i1,i2,j,n1,n2,num,ielem,icote,s1,s2,nmiss,ncote

	! nombre de cotes de l'element : 3 pour un T6 (6 noeuds, 3 aretes)
	ncote = nb_elem
	if (nb_elem == 6) ncote = 3

	do k=1,NS
		num = 0
		do i1=1,g_nlin
			if (any(mesh(k)%Lim == g_lc(i1))) num = num + 1
		enddo
		mesh(k)%NB = num
		if (allocated(mesh(k)%Bord)) deallocate(mesh(k)%Bord)
		allocate(mesh(k)%Bord(num,2))

		num = 0
		nmiss = 0
		do i1=1,g_nlin
			if (.not. any(mesh(k)%Lim == g_lc(i1))) cycle
			n1 = g_l1(i1)
			n2 = g_l2(i1)
			ielem = 0
			icote = 0
			do i2=1,mesh(k)%NE
				do j=1,ncote
					if (nb_elem == 6) then
						! T6 : les SOMMETS sont aux positions 1, 3, 5 (convention
						! Matlab, sommet/milieu alternes). Le cote j relie les
						! sommets 2j-1 et 2j+1 (mod 6), et l'element-ligne GMSH
						! de type 8 ne donne que ses deux extremites (le noeud
						! milieu est ignore a la lecture). Comparer des paires
						! consecutives comme en T3 ne matcherait jamais.
						s1 = mesh(k)%E(i2,2*j-1)
						s2 = mesh(k)%E(i2,mod(2*j,6)+1)
					else
						s1 = mesh(k)%E(i2,j)
						if (j < nb_elem) then
							s2 = mesh(k)%E(i2,j+1)
						else
							s2 = mesh(k)%E(i2,1)
						endif
					endif
					if ((s1==n1 .and. s2==n2) .or. (s1==n2 .and. s2==n1)) then
						ielem = i2
						icote = j
						exit
					endif
				enddo
				if (ielem > 0) exit
			enddo
			num = num + 1
			if (ielem == 0) nmiss = nmiss + 1
			mesh(k)%Bord(num,:) = (/ielem,icote/)
		enddo
		if (nmiss > 0) write(*,*) 'lecture_gmsh : ATTENTION, partition ',k,' : ',nmiss, &
			' arete(s) de limite non retrouvee(s) dans un element'
	enddo

	end subroutine construit_Bord

!--------------------------------------------------------------------
! mesh(k)%Ps / %table_Ps (sous-ensemble de points du groupe) et
! renumerotation de mesh(k)%E des indices GLOBAUX vers 1..mesh(k)%NP.
! A appeler EN DERNIER (construit_Bord a besoin des indices globaux ;
! Bord reference des numeros d'elements, pas de noeuds, donc il n'est
! pas affecte par cette renumerotation).
!--------------------------------------------------------------------
	subroutine construit_Ps_partition(nb_elem)
	implicit none
	integer, intent(in) :: nb_elem
	integer :: k,i1,i2,cnt
	integer, allocatable :: lookup(:)

	allocate(lookup(NP))
	do k=1,NS
		lookup = 0
		do i1=1,mesh(k)%NE
			do i2=1,nb_elem
				lookup(mesh(k)%E(i1,i2)) = 1
			enddo
		enddo
		cnt = count(lookup==1)
		mesh(k)%NP = cnt
		if (allocated(mesh(k)%table_Ps)) deallocate(mesh(k)%table_Ps)
		if (allocated(mesh(k)%Ps))       deallocate(mesh(k)%Ps)
		allocate(mesh(k)%table_Ps(cnt))
		allocate(mesh(k)%Ps(cnt,3))
		cnt = 0
		do i1=1,NP
			if (lookup(i1)==1) then
				cnt = cnt + 1
				mesh(k)%table_Ps(cnt) = i1
				mesh(k)%Ps(cnt,:) = Ps(i1,:)
				lookup(i1) = cnt      ! lookup devient la table global -> local
			endif
		enddo
		do i1=1,mesh(k)%NE
			do i2=1,nb_elem
				mesh(k)%E(i1,i2) = lookup(mesh(k)%E(i1,i2))
			enddo
		enddo
	enddo
	deallocate(lookup)

	end subroutine construit_Ps_partition

!--------------------------------------------------------------------
! type d'element GMSH correspondant au nombre de noeuds demande
!--------------------------------------------------------------------
	integer function type_gmsh(nb_elem) result(ty)
	implicit none
	integer, intent(in) :: nb_elem
	select case (nb_elem)
		case (3) ; ty = 2       ! triangle 3 noeuds
		case (4) ; ty = 3       ! quadrangle 4 noeuds
		case (6) ; ty = 9       ! triangle 6 noeuds (quadratique)
		case default
			write(*,*) 'lecture_gmsh : ERREUR, nb_elem=',nb_elem,' non gere (3, 4 ou 6)'
			stop 1
	end select
	end function type_gmsh

!--------------------------------------------------------------------
	integer function nb_mots(s) result(nb)
	implicit none
	character(len=*), intent(in) :: s
	integer :: i2
	nb = 1
	do i2 = 2, len_trim(s)
		if ( (s(i2:i2)==' ') .and. (s(i2-1:i2-1) /= ' ') ) nb = nb + 1
	enddo
	end function nb_mots

!--------------------------------------------------------------------
! extrait tous les entiers d'une ligne (separateurs espace/tabulation)
!--------------------------------------------------------------------
	subroutine lire_entiers(s,nmax,v,n)
	implicit none
	character(len=*), intent(in) :: s
	integer, intent(in) :: nmax
	integer, intent(out) :: v(nmax),n
	integer :: i1,L,deb,ios
	logical :: dans
	character(len=1) :: c

	n = 0
	deb = 1
	dans = .false.
	L = len_trim(s)
	do i1=1,L
		c = s(i1:i1)
		if (c /= ' ' .and. c /= char(9)) then
			if (.not.dans) then
				dans = .true.
				deb = i1
			endif
		else
			if (dans) then
				n = n + 1
				if (n > nmax) then
					write(*,*) 'lecture_gmsh : ERREUR, ligne trop longue (>',nmax,' entiers)'
					stop 1
				endif
				read(s(deb:i1-1),*,iostat=ios) v(n)
				dans = .false.
			endif
		endif
	enddo
	if (dans) then
		n = n + 1
		if (n > nmax) then
			write(*,*) 'lecture_gmsh : ERREUR, ligne trop longue (>',nmax,' entiers)'
			stop 1
		endif
		read(s(deb:L),*,iostat=ios) v(n)
	endif

	end subroutine lire_entiers

!--------------------------------------------------------------------
	subroutine libere_tampons()
	implicit none
	integer :: i1
	if (allocated(g_surf)) then
		do i1=1,size(g_surf)
			if (allocated(g_surf(i1)%phys)) deallocate(g_surf(i1)%phys)
			if (allocated(g_surf(i1)%bnd))  deallocate(g_surf(i1)%bnd)
		enddo
		deallocate(g_surf)
	endif
	if (allocated(g_Ps))   deallocate(g_Ps)
	if (allocated(g_E))    deallocate(g_E)
	if (allocated(g_phys)) deallocate(g_phys)
	if (allocated(g_etag)) deallocate(g_etag)
	if (allocated(g_lc))   deallocate(g_lc)
	if (allocated(g_l1))   deallocate(g_l1)
	if (allocated(g_l2))   deallocate(g_l2)
	g_nsurf = 0
	g_NP0 = 0
	g_NE0 = 0
	g_nlin = 0
	end subroutine libere_tampons


!--------------------------------------------------------------------
!
! call lecture_stl(file,nb_file)   (2026-10-02)
!
! Lecture d'un maillage STL (binaire ou ASCII, detecte par la taille du
! fichier : binaire si taille = 84 + 50 x nombre de facettes), triangles T3,
! pour les utilisateurs sans gmsh. Remplit NE, NP, Ps(NP,3), E(NE,3) comme
! lecture_simple_v2 :
!   - les sommets, repetes a chaque facette dans un STL, sont FUSIONNES
!     (coordonnees identiques au bit pres, tri lexicographique x, y, z) ;
!   - les triangles gardent l'ordre des facettes du fichier : le numero du
!     triangle source de Ray est le rang de la facette (a partir de 1), comme
!     dans python/source_oreille.py ;
!   - unites : le solveur travaille en metres ; si la plus grande etendue du
!     maillage depasse 2, le fichier est suppose en MILLIMETRES et converti
!     (x 1e-3), avec un message ;
!   - l'orientation des facettes (ordre des sommets) est gardee : normales
!     vers l'air, a verifier comme pour un .msh.
!--------------------------------------------------------------------
!
	subroutine lecture_stl(file,nb_file)
	use, intrinsic :: iso_fortran_env, only: int16, int32, int64, real32
	implicit none
	integer, intent(in) :: nb_file
	character(len=nb_file), intent(in) :: file
	integer :: u, ios, nt, nv, i, j, k
	integer(int32) :: n32
	integer(int16) :: attr
	integer(int64) :: taille
	real(real32) :: f(12)
	character(len=80) :: entete
	character(len=512) :: ligne
	real(prec), allocatable :: V(:,:)
	integer, allocatable :: idx(:), tmp(:), id(:)
	real(prec) :: ext
	logical :: binaire

	! binaire ou ASCII
	open(newunit=u, file=file, access='stream', form='unformatted', status='old', action='read')
	inquire(unit=u, size=taille)
	read(u) entete, n32
	binaire = (taille == 84_int64 + 50_int64*int(n32,int64))
	if (binaire) then
		nt = int(n32)
		nv = 3*nt
		allocate(V(nv,3))
		do i = 1, nt
			read(u) f, attr                         ! normale (3), 3 sommets (9), attribut
			do j = 1, 3
				V(3*(i-1)+j,:) = real(f(3*j+1:3*j+3), prec)
			end do
		end do
		close(u)
	else
		close(u)
		open(newunit=u, file=file, status='old', action='read')
		nv = 0
		do                                         ! 1re passe : nombre de sommets
			read(u,'(A)',iostat=ios) ligne
			if (ios /= 0) exit
			if (index(adjustl(ligne), 'vertex') == 1) nv = nv + 1
		end do
		rewind(u)
		allocate(V(nv,3))
		k = 0
		do
			read(u,'(A)',iostat=ios) ligne
			if (ios /= 0) exit
			ligne = adjustl(ligne)
			if (index(ligne, 'vertex') == 1) then
				k = k + 1
				read(ligne(7:),*) V(k,:)
			endif
		end do
		close(u)
		if (mod(nv,3) /= 0) then
			write(*,*) 'STL: number of vertices not a multiple of 3: ', nv
			stop 1
		endif
		nt = nv/3
	endif
	write(*,*) 'STL ', merge('binary', 'ASCII ', binaire), ': facets read: ', nt

	! fusion des sommets identiques : tri lexicographique (x, y, z) des indices
	allocate(idx(nv), tmp(nv), id(nv))
	idx = (/ (i, i=1,nv) /)
	call tri_fusion_xyz(1, nv)
	NP = 0
	do i = 1, nv
		if (i == 1) then
			NP = NP + 1
		elseif (any(V(idx(i),:) /= V(idx(i-1),:))) then
			NP = NP + 1
		endif
		id(idx(i)) = NP
	end do
	if (allocated(Ps)) deallocate(Ps)
	if (allocated(E)) deallocate(E)
	allocate(Ps(NP,3), E(nt,3))
	do i = 1, nv
		Ps(id(i),:) = V(i,:)
	end do
	NE = nt
	do i = 1, nt
		E(i,:) = id(3*(i-1)+1:3*i)
	end do
	ext = maxval(maxval(Ps,dim=1) - minval(Ps,dim=1))
	if (ext > 2.0_prec) then
		Ps = Ps*1.0e-3_prec
		write(*,'(A,F10.2,A)') ' STL: extent ', ext, ' > 2: file taken as millimetres, converted to metres'
	endif
	write(*,*) 'mesh: nodes kept:     ',NP
	write(*,*) 'mesh: elements kept:  ',NE
	deallocate(V, idx, tmp, id)

	contains

	! tri fusion recursif de idx(a:b) selon (x, y, z) de V
	recursive subroutine tri_fusion_xyz(a, b)
	integer, intent(in) :: a, b
	integer :: m, i1, i2, kk
	if (b <= a) return
	m = (a + b)/2
	call tri_fusion_xyz(a, m)
	call tri_fusion_xyz(m+1, b)
	i1 = a ; i2 = m+1 ; kk = a
	do while (i1 <= m .and. i2 <= b)
		if (avant(idx(i2), idx(i1))) then
			tmp(kk) = idx(i2) ; i2 = i2 + 1
		else
			tmp(kk) = idx(i1) ; i1 = i1 + 1
		endif
		kk = kk + 1
	end do
	do while (i1 <= m)
		tmp(kk) = idx(i1) ; i1 = i1 + 1 ; kk = kk + 1
	end do
	do while (i2 <= b)
		tmp(kk) = idx(i2) ; i2 = i2 + 1 ; kk = kk + 1
	end do
	idx(a:b) = tmp(a:b)
	end subroutine tri_fusion_xyz

	! vrai si le sommet p precede strictement le sommet q (ordre x, puis y, puis z)
	logical function avant(p, q)
	integer, intent(in) :: p, q
	if (V(p,1) /= V(q,1)) then
		avant = V(p,1) < V(q,1)
	elseif (V(p,2) /= V(q,2)) then
		avant = V(p,2) < V(q,2)
	else
		avant = V(p,3) < V(q,3)
	endif
	end function avant
	end subroutine lecture_stl


!--------------------------------------------------------------------
!
! call lecture_msh2_binaire(file,nb_file)   (2026-10-02)
!
! Lecture d'un maillage gmsh 2.2 BINAIRE (en-tete "2.2 1 8" : entiers 4
! octets, reels 8 octets, little-endian), triangles T3 seulement. Remplit NE,
! NP, Ps, E comme lecture_simple_v2. Les noeuds sont reperes par leur NUMERO
! (pas par leur rang dans le fichier) ; ceux qui ne servent a aucun triangle
! sont retires. Les triangles gardent leur ordre dans le fichier (numero du
! triangle source de Ray, comme python/source_oreille.py).
!--------------------------------------------------------------------
!
	subroutine lecture_msh2_binaire(file,nb_file)
	use, intrinsic :: iso_fortran_env, only: int32, real64
	implicit none
	integer, intent(in) :: nb_file
	character(len=nb_file), intent(in) :: file
	integer :: u, nn, nel, i, j, k, typ, nblk, ntag, nnod, idmax, ntri
	integer(int32) :: un, i4, h(3)
	integer(int32), allocatable :: rec(:), ids(:), tri(:,:)
	real(real64), allocatable :: X(:,:)
	real(real64) :: xyz(3)
	integer, allocatable :: rang(:), neuf(:)
	character(len=256) :: l

	open(newunit=u, file=file, access='stream', form='unformatted', status='old', action='read')
	call ligne(l)                                    ! $MeshFormat
	call ligne(l)                                    ! 2.2 1 8
	read(u) un
	if (un /= 1) then
		write(*,*) 'msh 2.2 binary: byte order not supported (big-endian)'
		stop 1
	endif
	call ligne(l)                                    ! fin de la ligne binaire
	call ligne(l)                                    ! $EndMeshFormat
	do                                               ! sections jusqu'a $Nodes
		call ligne(l)
		if (trim(l) == '$Nodes') exit
	end do
	call ligne(l) ; read(l,*) nn
	allocate(ids(nn), X(nn,3))
	do i = 1, nn
		read(u) i4, xyz
		ids(i) = i4 ; X(i,:) = xyz
	end do
	call ligne(l)                                    ! fin de la derniere ligne binaire
	call ligne(l)                                    ! $EndNodes
	call ligne(l)                                    ! $Elements
	call ligne(l) ; read(l,*) nel
	allocate(tri(nel,3))
	ntri = 0 ; k = 0
	do while (k < nel)                               ! blocs : type, nombre, nombre de tags
		read(u) h
		typ = h(1) ; nblk = h(2) ; ntag = h(3)
		select case (typ)
		case (1) ; nnod = 2
		case (2) ; nnod = 3
		case (3) ; nnod = 4
		case (4) ; nnod = 4
		case (8) ; nnod = 3
		case (9) ; nnod = 6
		case (15) ; nnod = 1
		case default
			write(*,*) 'msh 2.2 binary: element type not supported: ', typ
			stop 1
		end select
		allocate(rec(1+ntag+nnod))
		do j = 1, nblk
			read(u) rec
			if (typ == 2) then
				ntri = ntri + 1
				tri(ntri,:) = rec(2+ntag:4+ntag)
			endif
		end do
		deallocate(rec)
		k = k + nblk
	end do
	close(u)
	write(*,*) 'msh 2.2 binary: nodes ', nn, ', elements ', nel, ', triangles ', ntri

	! numero gmsh -> rang ; seuls les noeuds des triangles sont gardes
	idmax = maxval(ids)
	allocate(rang(idmax), neuf(nn))
	rang = 0
	do i = 1, nn
		rang(ids(i)) = i
	end do
	neuf = 0
	do i = 1, ntri
		do j = 1, 3
			neuf(rang(tri(i,j))) = 1
		end do
	end do
	NP = 0
	do i = 1, nn
		if (neuf(i) > 0) then
			NP = NP + 1 ; neuf(i) = NP
		endif
	end do
	if (allocated(Ps)) deallocate(Ps)
	if (allocated(E)) deallocate(E)
	allocate(Ps(NP,3), E(ntri,3))
	do i = 1, nn
		if (neuf(i) > 0) Ps(neuf(i),:) = real(X(i,:), prec)
	end do
	NE = ntri
	do i = 1, ntri
		do j = 1, 3
			E(i,j) = neuf(rang(tri(i,j)))
		end do
	end do
	write(*,*) 'mesh: nodes kept:     ',NP
	write(*,*) 'mesh: elements kept:  ',NE

	contains

	! lit une ligne de texte (jusqu'a LF) dans le fichier ouvert en flux
	subroutine ligne(t)
	character(len=*), intent(out) :: t
	character :: c
	integer :: n
	t = '' ; n = 0
	do
		read(u) c
		if (c == achar(10)) exit
		if (c == achar(13)) cycle
		n = n + 1
		if (n <= len(t)) t(n:n) = c
	end do
	end subroutine ligne
	end subroutine lecture_msh2_binaire
!--------------------------------------------------------------------
! lecture du maillage HRTF (2026-10-02), triangles T3 :
!   .stl / .STL                      -> lecture_stl (binaire ou ASCII)
!   gmsh 2.2 binaire (en-tete 2.2 1) -> lecture_msh2_binaire
!   sinon gmsh 2.2 ASCII             -> lecture_simple_v2
!--------------------------------------------------------------------
	subroutine lecture_maillage_hrtf(file)
	character(len=*), intent(in) :: file
	integer :: n, u, ios, ftype
	character(len=256) :: l1, l2
	real(prec) :: ver
	n = len_trim(file)
	if (n > 4) then
		if (file(n-3:n) == '.stl' .or. file(n-3:n) == '.STL') then
			call lecture_stl(trim(file), n)
			return
		endif
	endif
	ftype = 0
	open(newunit=u, file=trim(file), status='old', action='read', iostat=ios)
	if (ios == 0) then
		read(u,'(A)',iostat=ios) l1
		read(u,'(A)',iostat=ios) l2
		if (ios == 0 .and. index(l1, '$MeshFormat') > 0) read(l2,*,iostat=ios) ver, ftype
		close(u)
	endif
	if (ftype == 1) then
		call lecture_msh2_binaire(trim(file), n)
	else
		call lecture_simple_v2(trim(file), n, 3)
	endif
	end subroutine lecture_maillage_hrtf
	end module module_GMSH
