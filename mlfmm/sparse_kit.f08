! Copyright (c) 2024-2026 Conservatoire national des arts et metiers (Cnam), Christophe Langrenne
! SPDX-License-Identifier: LGPL-2.1-or-later
! Contains routines of SPARSKIT (Y. Saad), Copyright (C) 2005 the Regents of the
! University of Minnesota, GNU LGPL 2.1 or later: this file is LGPL-2.1-or-later.
      subroutine coo_zmatvec(nla,nca,nnz,A,ia,ja,x,y)
! multiplication Ax=y avec A une matrice creuse
! en coordonnées coo
      implicit none
! variables d'entrée      
      integer, intent(in) :: nla,nca,nnz
      integer, intent(in),pointer :: ia(:),ja(:)
      double complex, intent(in),pointer :: A(:)
      double complex, intent(in) :: x(nca)
! variables de sortie
      double complex, intent(inout) :: y(nla)
! variables locales
      integer :: i1
      
      y(:)=dcmplx(0.0d0,0.0d0)
      do i1=1,nnz
	  !write(*,*) ia(i1),ja(i1)
	 y(ia(i1))=y(ia(i1))+A(i1)*x(ja(i1))
      end do
      
      end subroutine coo_zmatvec
!------------- end of coo_zmatvec ------------------------------------------- 
!----------------------------------------------------------------------- 
      
      subroutine zcoocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
!----------------------------------------------------------------------- 
      double complex a(*),ao(*),x
      integer ir(*),jc(*),jao(*),iao(*)
!      integer, intent(inout) :: iao(nrow+1)
!-----------------------------------------------------------------------
!  Coordinate     to   Compressed Sparse Row 
!----------------------------------------------------------------------- 
! converts a matrix that is stored in coordinate format
!  a, ir, jc into a row general sparse ao, jao, iao format.
!
! on entry:
!--------- 
! nrow	= dimension of the matrix 
! nnz	= number of nonzero elements in matrix
! a,
! ir, 
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
! 	  the elements, ir(k) = its row number and jc(k) = its column 
!	  number. The order of the elements is arbitrary. 
!
! on return:
!----------- 
! ir 	is destroyed
!
! ao, jao, iao = matrix in general sparse matrix format with ao 
! 	continung the real values, jao containing the column indices, 
!	and iao being the pointer to the beginning of the row, 
!	in arrays ao, jao.
!
! Notes:
!------ This routine is NOT in place.  See coicsr
!
!------------------------------------------------------------------------
      integer :: k
      
!      write(*,*) 'in zcoocsr'
!      iao=[(0,k=1,nrow+1)]
      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
!      write(*,*) 'raz ia0'
! determine row-lengths.
      do 2 k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
 2    continue
!      write(*,*) 'fin allocation iao'
! starting position of each row..
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
! go through the structure  once more. Fill in output matrix.
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
! shift back iao
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
!------------- end of coocsr ------------------------------------------- 
!----------------------------------------------------------------------- 
      end
!----------------------------------------------------------------------- 
      subroutine zcoicsr (n,nnz,job,a,ja,ia,iwk)
      integer ia(nnz),ja(nnz),iwk(n+1) 
      double complex a(*)
!------------------------------------------------------------------------
! IN-PLACE coo-csr conversion routine.
!------------------------------------------------------------------------
! this subroutine converts a matrix stored in coordinate format into 
! the csr format. The conversion is done in place in that the arrays 
! a,ja,ia of the result are overwritten onto the original arrays.
!------------------------------------------------------------------------
! on entry:
!--------- 
! n	= integer. row dimension of A.
! nnz	= integer. number of nonzero elements in A.
! job   = integer. Job indicator. when job=1, the real values in a are
!         filled. Otherwise a is not touched and the structure of the
!         array only (i.e. ja, ia)  is obtained.
! a	= real array of size nnz (number of nonzero elements in A)
!         containing the nonzero elements 
! ja	= integer array of length nnz containing the column positions
! 	  of the corresponding elements in a.
! ia	= integer array of length nnz containing the row positions
! 	  of the corresponding elements in a.
! iwk	= integer work array of length n+1 
! on return:
!----------
! a
! ja 
! ia	= contains the compressed sparse row data structure for the 
!         resulting matrix.
! Note: 
!-------
!         the entries of the output matrix are not sorted (the column
!         indices in each are not in increasing order) use coocsr
!         if you want them sorted.
!----------------------------------------------------------------------c
!  Coded by Y. Saad, Sep. 26 1989                                      c
!----------------------------------------------------------------------c
      real*8 t,tnext
      logical values
!----------------------------------------------------------------------- 
      values = (job .eq. 1) 
! find pointer array for resulting matrix. 
      do 35 i=1,n+1
         iwk(i) = 0
 35   continue
      do 4 k=1,nnz
         i = ia(k)
         iwk(i+1) = iwk(i+1)+1
 4    continue 
!----------------------------------------------------------------------
      iwk(1) = 1 
      do 44 i=2,n
         iwk(i) = iwk(i-1) + iwk(i)
 44   continue 
!
!     loop for a cycle in chasing process. 
!
      init = 1
      k = 0
 5    if (values) t = a(init)
      i = ia(init)
      j = ja(init)
      ia(init) = -1
!----------------------------------------------------------------------
 6    k = k+1 		
!     current row number is i.  determine  where to go. 
      ipos = iwk(i)
!     save the chased element. 
      if (values) tnext = a(ipos)
      inext = ia(ipos)
      jnext = ja(ipos)
!     then occupy its location.
      if (values) a(ipos)  = t
      ja(ipos) = j
!     update pointer information for next element to come in row i. 
      iwk(i) = ipos+1
!     determine  next element to be chased,
      if (ia(ipos) .lt. 0) goto 65
      t = tnext
      i = inext
      j = jnext 
      ia(ipos) = -1
      if (k .lt. nnz) goto 6
      goto 70
 65   init = init+1
      if (init .gt. nnz) goto 70
      if (ia(init) .lt. 0) goto 65
!     restart chasing --	
      goto 5
 70   do 80 i=1,n 
         ia(i+1) = iwk(i)
 80   continue
      ia(1) = 1
      return
!----------------- end of coicsr ----------------------------------------
!------------------------------------------------------------------------
      end
      
      subroutine ilut(n,nnz,a,ja,ia,lfil,droptol,alu,jlu,ju,iwk,w,jw,ierr)
!-----------------------------------------------------------------------
      implicit none 
      integer n,nnz
      double precision droptol
      double complex a(nnz),alu(iwk),w(n+1)
      integer ja(nnz),ia(n+1),jlu(iwk),ju(n),jw(2*n),lfil,iwk,ierr
!----------------------------------------------------------------------*
!                      *** ILUT preconditioner ***                     *
!      incomplete LU factorization with dual truncation mechanism      *
!----------------------------------------------------------------------*
!     Author: Yousef Saad *May, 5, 1990, Latest revision, August 1996  *
!----------------------------------------------------------------------*
! PARAMETERS                                                           
!-----------                                                           
!
! on entry:
!========== 
! n       = integer. The row dimension of the matrix A. The matrix 
!
! a,ja,ia = matrix stored in Compressed Sparse Row format.              
!
! lfil    = integer. The fill-in parameter. Each row of L and each row
!           of U will have a maximum of lfil elements (excluding the 
!           diagonal element). lfil must be .ge. 0.
!           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
!           EARLIER VERSIONS. 
!
! droptol = real*8. Sets the threshold for dropping small terms in the
!           factorization. See below for details on dropping strategy.
!
!  
! iwk     = integer. The lengths of arrays alu and jlu. If the arrays
!           are not big enough to store the ILU factorizations, ilut
!           will stop with an error message. 
!
! On return:
!===========
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju      = integer array of length n containing the pointers to
!           the beginning of each row of U in the matrix alu,jlu.
!
! ierr    = integer. Error message with the following meaning.
!           ierr  = 0    --> successful return.
!           ierr .gt. 0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.
!                            (The elimination process has generated a
!                            row in L or U whose length is .gt.  n.)
!           ierr  = -2   --> The matrix L overflows the array al.
!           ierr  = -3   --> The matrix U overflows the array alu.
!           ierr  = -4   --> Illegal value for lfil.
!           ierr  = -5   --> zero row encountered.
!
! work arrays:
!=============
! jw      = integer work array of length 2*n.
! w       = real work array of length n+1.
!  
!----------------------------------------------------------------------
! w, ju (1:n) store the working array [1:ii-1 = L-part, ii:n = u] 
! jw(n+1:2n)  stores nonzero indicators
! 
! Notes:
! ------
! The diagonal elements of the input matrix must be  nonzero (at least
! 'structurally'). 
!
!----------------------------------------------------------------------* 
!---- Dual drop strategy works as follows.                             *
!                                                                      *
!     1) Theresholding in L and U as set by droptol. Any element whose *
!        magnitude is less than some tolerance (relative to the abs    *
!        value of diagonal element in u) is dropped.                   *
!                                                                      *
!     2) Keeping only the largest lfil elements in the i-th row of L   * 
!        and the largest lfil elements in the i-th row of U (excluding *
!        diagonal elements).                                           *
!                                                                      *
! Flexibility: one  can use  droptol=0  to get  a strategy  based on   *
! keeping  the largest  elements in  each row  of L  and U.   Taking   *
! droptol .ne.  0 but lfil=n will give  the usual threshold strategy   *
! (however, fill-in is then mpredictible).                             *
!----------------------------------------------------------------------*
!     locals
      integer ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,len 
      double precision tnorm, abs
      double complex t, s, fact
      
      if (lfil .lt. 0) goto 998
!-----------------------------------------------------------------------
!     initialize ju0 (points to next element to be added to alu,jlu)
!     and pointer array.
!-----------------------------------------------------------------------
!      write(*,*) 'init ju0 et jlu(1)'
      ju0 = n+2
      jlu(1) = ju0
!
!     initialize nonzero indicator array. 
!
!      write(*,*) 'init jw'
      do 1 j=1,n
!	  write(*,*) j
         jw(n+j)  = 0
 1    continue
!-----------------------------------------------------------------------
!     beginning of main loop.
!-----------------------------------------------------------------------
!      write(*,*) 'begin boucle'
      do 500 ii = 1, n
!	  write(*,*) 'boucle ii=',ii
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
         tnorm = 0.0d0
         do 501 k=j1,j2
            tnorm = tnorm+abs(a(k))
 501     continue
         if (tnorm .eq. 0.0) goto 999
         tnorm = tnorm/real(j2-j1+1)
!     
!     unpack L-part and U-part of row of A in arrays w 
!     
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = dcmplx(0.0d0,0.0d0)
         jw(n+ii) = ii
!
         do 170  j = j1, j2
            k = ja(j)
            t = a(j)
!	    if (isnan(abs(t))) then
!		    write(*,*) 't isnan = ',i
!	    end if
            
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n+k) = lenl
            else if (k .eq. ii) then
               w(ii) = t
            else
               lenu = lenu+1
               jpos = ii+lenu-1 
               jw(jpos) = k
               w(jpos) = t
               jw(n+k) = jpos
            endif
 170     continue
         jj = 0
         len = 0 
!     
!     eliminate previous rows
!     
 150     jj = jj+1
         if (jj .gt. lenl) goto 160
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.
!-----------------------------------------------------------------------
         jrow = jw(jj)
         k = jj
!     
!     determine smallest column index
!     
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            endif
 151     continue
!
         if (k .ne. jj) then
!     exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
!     exchange in jr
            jw(n+jrow) = jj
            jw(n+j) = k
!     exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
!
!     zero out element in row by setting jw(n+jrow) to zero.
!     
         jw(n+jrow) = 0
!
!     get the multiplier for row to be eliminated (jrow).
!     
         fact = w(jj)*alu(jrow)
         if (abs(fact) .le. droptol) goto 150
!     
!     combine current row and row jrow
!
         do 203 k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
!	    if (isnan(abs(s))) then
!		write(*,*) 's isnan = ',fact,alu(k)
!	    end if
           
            j = jlu(k)
            jpos = jw(n+j)
            if (j .ge. ii) then
!     
!     dealing with upper part.
!     
               if (jpos .eq. 0) then
!
!     this is a fill-in element
!     
                  lenu = lenu+1
                  if (lenu .gt. n) goto 995
                  i = ii+lenu-1
                  jw(i) = j
                  jw(n+j) = i
                  w(i) = - s
               else
!
!     this is not a fill-in element 
!
                  w(jpos) = w(jpos) - s

               endif
            else
!     
!     dealing  with lower part.
!     
               if (jpos .eq. 0) then
!
!     this is a fill-in element
!     
                  lenl = lenl+1
                  if (lenl .gt. n) goto 995
                  jw(lenl) = j
                  jw(n+j) = lenl
                  w(lenl) = - s
               else
!     
!     this is not a fill-in element 
!     
                  w(jpos) = w(jpos) - s
               endif
            endif
 203     continue
!     
!     store this pivot element -- (from left to right -- no danger of
!     overlap with the working elements in L (pivots). 
!     
!	    if (isnan(abs(fact))) then
!		write(*,*) 'fact isnan = ',fact
!	    end if
         len = len+1 
         w(len) = fact
         jw(len)  = jrow
         goto 150
 160     continue
!     
!     reset double-pointer to zero (U-part)
!     
         do 308 k=1, lenu
            jw(n+jw(ii+k-1)) = 0
 308     continue
!     
!     update L-matrix
!     
         lenl = len 
         len = min0(lenl,lfil)
!     
!     sort by quick-split
!
         call qsplit (w,jw,lenl,len)
!
!     store L-part
! 
         do 204 k=1, len 
            if (ju0 .gt. iwk) goto 996
            alu(ju0) =  w(k)
            jlu(ju0) =  jw(k)
            ju0 = ju0+1
 204     continue
!     
!     save pointer to beginning of row ii of U
!     
         ju(ii) = ju0
!
!     update U-matrix -- first apply dropping strategy 
!
         len = 0
         do k=1, lenu-1
            if (abs(w(ii+k)) .gt. droptol*tnorm) then 
               len = len+1
               w(ii+len) = w(ii+k) 
               jw(ii+len) = jw(ii+k) 
            endif
         enddo
         lenu = len+1
         len = min0(lenu,lfil)
!
         call qsplit (w(ii+1), jw(ii+1), lenu-1,len)
!
!     copy
! 
         t = abs(w(ii))
         if (len + ju0 .gt. iwk) goto 997
         do 302 k=ii+1,ii+len-1 
            jlu(ju0) = jw(k)
            alu(ju0) = w(k)
            t = t + abs(w(k) )
            ju0 = ju0+1
 302     continue
!     
!     store inverse of diagonal element of u
!     
         if (abs(w(ii)) .eq. 0.0d0) w(ii) = (dcmplx(0.0001d0,0.0d0) + droptol)*tnorm
!     	
!	 if (isnan(abs(w(ii)))) then
!	    write(*,*) 'w(ii) =',ii
!	 endif
         alu(ii) = 1.0d0/ w(ii) 
!     
!     update pointer to beginning of next row of U.
!     
         jlu(ii+1) = ju0
!-----------------------------------------------------------------------
!     end main loop
!-----------------------------------------------------------------------
 500  continue
      ierr = 0
      return
!
!     incomprehensible error. Matrix must be wrong.
!     
 995  ierr = -1
      return
!     
!     insufficient storage in L.
!     
 996  ierr = -2
      return
!     
!     insufficient storage in U.
!     
 997  ierr = -3
      return
!     
!     illegal lfil entered.
!     
 998  ierr = -4
      return
!     
!     zero row encountered
!     
 999  ierr = -5
      return
!----------------end-of-ilut--------------------------------------------
!-----------------------------------------------------------------------
      end
      
! fonction isnan      
!      logical function isnan(a)
!      double precision a
!      if (a.ne.a) then
!	  isnan = .true.
!	else
!	  isnan = .false.
!	end if
!	return
!	end
     
	subroutine lusol(n, y, x, alu, jlu, ju,iwk)
	implicit none
	integer n,iwk
	integer jlu(iwk),ju(n)
        double complex x(n), y(n), alu(iwk)
!-----------------------------------------------------------------------
!
! This routine solves the system (LU) x = y, 
! given an LU decomposition of a matrix stored in (alu, jlu, ju) 
! modified sparse row format 
!
!-----------------------------------------------------------------------
! on entry:
! n   = dimension of system 
! y   = the right-hand-side vector
! alu, jlu, ju 
!     = the LU matrix as provided from the ILU routines. 
!
! on return
! x   = solution of LU x = y.     
!-----------------------------------------------------------------------
! 
! Note: routine is in place: call lusol (n, x, x, alu, jlu, ju) 
!       will solve the system with rhs x and overwrite the result on x . 
!
!-----------------------------------------------------------------------
! local variables
!
        integer i,k
!
! forward solve
!
!	write(*,*) 'debut forward'
	
        do 40 i = 1, n
           x(i) = y(i)
           do 41 k=jlu(i),ju(i)-1
!	      write(*,*) 'k =',k
              x(i) = x(i) - alu(k)* x(jlu(k))
 41        continue
 40     continue
!
!     backward solve.
!
!	write(*,*) 'debut backward'
	do 90 i = n, 1, -1
	   do 91 k=ju(i),jlu(i+1)-1
              x(i) = x(i) - alu(k)*x(jlu(k))
 91	   continue
           x(i) = alu(i)*x(i)
 90     continue
!
  	return
!----------------end of lusol ------------------------------------------
!-----------------------------------------------------------------------
	end
	
      subroutine ilutp(n,nnz,a,ja,ia,lfil,droptol,permtol,mbloc,alu, &
          jlu,ju,iwk,w,jw,iperm,ierr)
!-----------------------------------------------------------------------
!     implicit none
      integer n,ja(nnz),ia(n+1),lfil,jlu(iwk),ju(n),jw(2*n),iwk, &
          iperm(2*n),ierr
      double complex a(nnz), alu(iwk), w(n+1)
      double precision droptol
!----------------------------------------------------------------------*
!       *** ILUTP preconditioner -- ILUT with pivoting  ***            *
!      incomplete LU factorization with dual truncation mechanism      *
!----------------------------------------------------------------------*
! author Yousef Saad *Sep 8, 1993 -- Latest revision, August 1996.     *
!----------------------------------------------------------------------*
! on entry:
!==========
! n       = integer. The dimension of the matrix A.
!
! a,ja,ia = matrix stored in Compressed Sparse Row format.
!           ON RETURN THE COLUMNS OF A ARE PERMUTED. SEE BELOW FOR 
!           DETAILS. 
!
! lfil    = integer. The fill-in parameter. Each row of L and each row
!           of U will have a maximum of lfil elements (excluding the 
!           diagonal element). lfil must be .ge. 0.
!           ** WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
!           EARLIER VERSIONS. 
!
! droptol = real*8. Sets the threshold for dropping small terms in the
!           factorization. See below for details on dropping strategy.
!
! lfil    = integer. The fill-in parameter. Each row of L and
!           each row of U will have a maximum of lfil elements.
!           WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
!           EARLIER VERSIONS. 
!           lfil must be .ge. 0.
!
! permtol = tolerance ratio used to  determne whether or not to permute
!           two columns.  At step i columns i and j are permuted when 
!
!                     abs(a(i,j))*permtol .gt. abs(a(i,i))
!
!           [0 --> never permute; good values 0.1 to 0.01]
!
! mbloc   = if desired, permuting can be done only within the diagonal
!           blocks of size mbloc. Useful for PDE problems with several
!           degrees of freedom.. If feature not wanted take mbloc=n.
!
!  
! iwk     = integer. The lengths of arrays alu and jlu. If the arrays
!           are not big enough to store the ILU factorizations, ilut
!           will stop with an error message. 
!
! On return:
!===========
!
! alu,jlu = matrix stored in Modified Sparse Row (MSR) format containing
!           the L and U factors together. The diagonal (stored in
!           alu(1:n) ) is inverted. Each i-th row of the alu,jlu matrix
!           contains the i-th row of L (excluding the diagonal entry=1)
!           followed by the i-th row of U.
!
! ju      = integer array of length n containing the pointers to
!           the beginning of each row of U in the matrix alu,jlu.
!
! iperm   = contains the permutation arrays. 
!           iperm(1:n) = old numbers of unknowns
!           iperm(n+1:2*n) = reverse permutation = new unknowns.
!
! ierr    = integer. Error message with the following meaning.
!           ierr  = 0    --> successful return.
!           ierr .gt. 0  --> zero pivot encountered at step number ierr.
!           ierr  = -1   --> Error. input matrix may be wrong.
!                            (The elimination process has generated a
!                            row in L or U whose length is .gt.  n.)
!           ierr  = -2   --> The matrix L overflows the array al.
!           ierr  = -3   --> The matrix U overflows the array alu.
!           ierr  = -4   --> Illegal value for lfil.
!           ierr  = -5   --> zero row encountered.
!
! work arrays:
!=============
! jw      = integer work array of length 2*n.
! w       = real work array of length n 
!
! IMPORTANR NOTE:
! --------------
! TO AVOID PERMUTING THE SOLUTION VECTORS ARRAYS FOR EACH LU-SOLVE, 
! THE MATRIX A IS PERMUTED ON RETURN. [all column indices are
! changed]. SIMILARLY FOR THE U MATRIX. 
! To permute the matrix back to its original state use the loop:
!
!      do k=ia(1), ia(n+1)-1
!         ja(k) = iperm(ja(k)) 
!      enddo
! 
!-----------------------------------------------------------------------
!     local variables
!
      integer k,i,j,jrow,ju0,ii,j1,j2,jpos,len,imax,lenu,lenl,jj,mbloc,&
          icut
      double precision tnorm, abs, permtol
      double complex s,t,fact,tmp,xmax,xmax0
      character(len=1) :: answer
!     
      write(*,*) 'in ilutp'
      if (lfil .lt. 0) goto 998
!----------------------------------------------------------------------- 
!     initialize ju0 (points to next element to be added to alu,jlu)
!     and pointer array.
!-----------------------------------------------------------------------
      ju0 = n+2
      jlu(1) = ju0
!
!  integer double pointer array.
!
      do 1 j=1, n
         jw(n+j)  = 0
         iperm(j) = j
         iperm(n+j) = j
 1    continue
     write(*,*) 'fin initialisation'
!      write(*,*) 'shape iperm =',shape(iperm)
!      write(*,*) iperm(22126:22127)
!      read(*,*) answer
!-----------------------------------------------------------------------
!     beginning of main loop.
!-----------------------------------------------------------------------
      do 500 ii = 1, n
	 write(*,*) 'main boucle =',ii,' sur ',n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
!         write(*,*) 'j1,j2 =',j1,j2
         tnorm = 0.0d0
         do 501 k=j1,j2
            tnorm = tnorm+abs(a(k))
 501     continue
         if (tnorm .eq. 0.0) goto 999
         tnorm = tnorm/(j2-j1+1)
!
!     unpack L-part and U-part of row of A in arrays  w  --
!
         lenu = 1
         lenl = 0
         jw(ii) = ii
         w(ii) = dcmplx(0.0d0,0.0d0)
         jw(n+ii) = ii
!
         do 170  j = j1, j2
!            if (iperm(22127).ne.22127) then
!	      write(*,*) 'problem debut 170, j1 et j2',j1,j2
!	      read(*,*) answer
!	    endif

	    !write(*,*) 'n+ja =',n+ja(j)
            k = iperm(n+ja(j))
            t = a(j)
            !write(*,*) 'k = ',k
!            if (iperm(22127).ne.22127) then
!	      write(*,*) 'problem 170 avant if',j1,j2
!	      read(*,*) answer
!	    endif
            
            if (k .lt. ii) then
               lenl = lenl+1
               jw(lenl) = k
               w(lenl) = t
               jw(n+k) = lenl
!	    if (iperm(22127).ne.22127) then
!	      write(*,*) 'problem iperm k<ii',lenl,n+k
!	      read(*,*) answer
!	    endif
 !           if (n+k.gt.2*n) then
!		write(*,*) 'problem k<ii',n+k
!		read(*,*) answer
!	    endif
!	    if (lenl.gt.n+1) then
!		write(*,*) 'problem k<ii',lenl
!		read(*,*) answer
!	    endif
               
            else if (k .eq. ii) then
          
               w(ii) = t
            if (iperm(22127).ne.22127) then
	      write(*,*) 'problem iperm k=ii',ii
	      read(*,*) answer
	    endif
               
            else

               lenu = lenu+1
               jpos = ii+lenu-1 
               jw(jpos) = k   
               w(jpos) = t
               jw(n+k) = jpos
            
               
            endif

 170     continue
!	 write(*,*) 'boucle 170'
         jj = 0
         len = 0 

         !
!     eliminate previous rows
!
 150     jj = jj+1
         if (jj .gt. lenl) goto 160
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.
!-----------------------------------------------------------------------
         jrow = jw(jj)
         k = jj
!
!     determine smallest column index
!
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            endif
 151     continue
!	 write(*,*) 'boucle 151'

         if (k .ne. jj) then
!     exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
!     exchange in jr
            jw(n+jrow) = jj
            jw(n+j) = k
!     exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
!
!     zero out element in row by resetting jw(n+jrow) to zero.
!     
         jw(n+jrow) = 0
!
!     get the multiplier for row to be eliminated: jrow
!
         fact = w(jj)*alu(jrow)
!
!     drop term if small
!     
         if (abs(fact) .le. droptol) goto 150
!
!     combine current row and row jrow
!
         do 203 k = ju(jrow), jlu(jrow+1)-1
            s = fact*alu(k)
!     new column number
            j = iperm(n+jlu(k))
            jpos = jw(n+j)
            if (j .ge. ii) then
!
!     dealing with upper part.
!
               if (jpos .eq. 0) then
!
!     this is a fill-in element
!
                  lenu = lenu+1
                  i = ii+lenu-1 
                  if (lenu .gt. n) goto 995
                  jw(i) = j
                  jw(n+j) = i 
                  w(i) = - s
               else
!     no fill-in element --
                  w(jpos) = w(jpos) - s
               endif
            else
!
!     dealing with lower part.
!
               if (jpos .eq. 0) then
!
!     this is a fill-in element
!
                 lenl = lenl+1
                 if (lenl .gt. n) goto 995
                 jw(lenl) = j
                 jw(n+j) = lenl
                 w(lenl) = - s
              else
!
!     this is not a fill-in element
!
                 w(jpos) = w(jpos) - s
              endif
           endif
 203	continue
!      write(*,*) 'boucle 203'
     
!     
!     store this pivot element -- (from left to right -- no danger of
!     overlap with the working elements in L (pivots). 
!     
        len = len+1 
        w(len) = fact
        jw(len)  = jrow
	goto 150
 160    continue
!
!     reset double-pointer to zero (U-part)
!     
        do 308 k=1, lenu
           jw(n+jw(ii+k-1)) = 0
 308	continue
!	write(*,*) 'boucle 308'
	
!
!     update L-matrix
!
        lenl = len 
        len = min0(lenl,lfil)
!     
!     sort by quick-split
!
        call qsplit (w,jw,lenl,len)
!
!     store L-part -- in original coordinates ..
!
        do 204 k=1, len
           if (ju0 .gt. iwk) goto 996
           alu(ju0) =  w(k)  
           jlu(ju0) = iperm(jw(k))
           ju0 = ju0+1
 204    continue
!	write(*,*) 'boucle 204'
	
!
!     save pointer to beginning of row ii of U
!
        ju(ii) = ju0
!
!     update U-matrix -- first apply dropping strategy 
!
         len = 0
         do k=1, lenu-1
            if (abs(w(ii+k)) .gt. droptol*tnorm) then 
               len = len+1
               w(ii+len) = w(ii+k) 
               jw(ii+len) = jw(ii+k) 
            endif
         enddo
         lenu = len+1
         len = min0(lenu,lfil)
         call qsplit (w(ii+1), jw(ii+1), lenu-1,len)
!         write(*,*) 'update U'
        
!
!     determine next pivot -- 
!
        imax = ii
        xmax = abs(w(imax))
        xmax0 = xmax
        icut = ii - 1 + mbloc - mod(ii-1,mbloc)
        do k=ii+1,ii+len-1
           t = abs(w(k))
           if (abs(t) .gt. abs(xmax) .and. abs(t*permtol) .gt. abs(xmax0) .and. &
               jw(k) .le. icut) then
              imax = k
              xmax = t
           endif
        enddo
!        write(*,*) ' determine next pivot'
       
!
!     exchange w's
!
!	write(*,*) 'exchange'
        tmp = w(ii)
        w(ii) = w(imax)
        w(imax) = tmp
 
!
!     update iperm and reverse iperm
!
!	write(*,*) 'update iperm'
       j = jw(imax)
        i = iperm(ii)
!        write(*,*) 'indice i et j',i,j
        iperm(ii) = iperm(j)
        iperm(j) = i
!
!     reverse iperm
!
!	write(*,*) 'reverse iperm'
!       write(*,*) 'reverse ii',n,ii,shape(iperm)
!        write(*,*) iperm(ii)
        iperm(n+iperm(ii)) = ii
!        write(*,*) 'reverse j',j,iperm(j)
        iperm(n+iperm(j)) = j
!        write(*,*) 'reverse j'
!----------------------------------------------------------------------- 
!
        if (len + ju0 .gt. iwk) goto 997
!
!     copy U-part in original coordinates
!     
!	write(*,*) 'avant boucle 302'
        do 302 k=ii+1,ii+len-1 
           jlu(ju0) = iperm(jw(k))
           alu(ju0) = w(k)
           ju0 = ju0+1
 302	continue
!	write(*,*) 'boucle 302'
!
!     store inverse of diagonal element of u
!
        if (w(ii) .eq. 0.0) w(ii) = (1.0D-4 + droptol)*tnorm
        alu(ii) = 1.0d0/ w(ii) 
!
!     update pointer to beginning of next row of U.
!
!	write(*,*) 'update pointer'
	jlu(ii+1) = ju0
!	write(*,*) 'bas de boucle'
!-----------------------------------------------------------------------
!     end main loop
!-----------------------------------------------------------------------
 500  continue
!      write(*,*) 'fin boucle principale'
!
!     permute all column indices of LU ...
!
      do k = jlu(1),jlu(n+1)-1
         jlu(k) = iperm(n+jlu(k))
      enddo

!     ...and of A
!
      do k=ia(1), ia(n+1)-1
         ja(k) = iperm(n+ja(k))
      enddo
!
      ierr = 0
      return
!
!     incomprehensible error. Matrix must be wrong.
!
 995  ierr = -1
      return
!
!     insufficient storage in L.
!
 996  ierr = -2
      return
!
!     insufficient storage in U.
!
 997  ierr = -3
      return
!
!     illegal lfil entered.
!
 998  ierr = -4
      return
!
!     zero row encountered
!
 999  ierr = -5
      return
!----------------end-of-ilutp-------------------------------------------
!-----------------------------------------------------------------------
      end
     
      subroutine qsplit(a,ind,n,ncut)
      double complex a(n)
        integer ind(n), n, ncut
!-----------------------------------------------------------------------
!     does a quick-sort split of a real array.
!     on input a(1:n). is a real array
!     on output a(1:n) is permuted such that its elements satisfy:
!
!     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
!     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
!
!     ind(1:n) is an integer array which permuted in the same way as a(*).
!-----------------------------------------------------------------------
        double complex tmp
        double precision abskey
        integer itmp, first, last
!-----
        first = 1
        last = n
        if (ncut .lt. first .or. ncut .gt. last) return
!
!     outer loop -- while mid .ne. ncut do
!
 1      mid = first
        abskey = abs(a(mid))
        do 2 j=first+1, last
           if (abs(a(j)) .gt. abskey) then
              mid = mid+1
!     interchange
              tmp = a(mid)
              itmp = ind(mid)
              a(mid) = a(j)
              ind(mid) = ind(j)
              a(j)  = tmp
              ind(j) = itmp
           endif
 2      continue
!
!     interchange
!
        tmp = a(mid)
        a(mid) = a(first)
        a(first)  = tmp
!
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
!
!     test for while loop
!
        if (mid .eq. ncut) return
        if (mid .gt. ncut) then
           last = mid-1
        else
           first = mid+1
        endif
        goto 1
!----------------end-of-qsplit------------------------------------------
!-----------------------------------------------------------------------
        end

! ==================================================================
! Helpers CSR/COO reels ajoutes pour le Laplacien Beltrami (MUMPS)
! Portes verbatim depuis l'ancien sparse_kit (SPARSKIT + utilitaires)
! ==================================================================
      subroutine coo_zmatvec_sym(nla,nca,nnz,A,ia,ja,x,y)
! multiplication Ax=y avec A une matrice creuse symétrique
! en coordonnées coo, stockage A matrice supérieure ou inférieure
      implicit none
! variables d'entrée
      integer, intent(in) :: nla,nca,nnz
      integer, intent(in) :: ia(nnz),ja(nnz)
      double complex, intent(in) :: A(nnz)
      double complex, intent(in) :: x(nca)
! variables de sortie
      double complex, intent(inout) :: y(nla)
! variables locales
      integer :: i1

      y(:)=dcmplx(0.0d0,0.0d0)
      do i1=1,nnz
	  !write(*,*) ia(i1),ja(i1)
      y(ia(i1))=y(ia(i1))+A(i1)*x(ja(i1))
      if (ia(i1) /= ja(i1)) then
            y(ja(i1))=y(ja(i1))+A(i1)*x(ia(i1))
      endif
      end do

      end subroutine coo_zmatvec_sym
      subroutine coo_dmatvec(nla,nca,nnz,A,ia,ja,x,y)
! multiplication Ax=y avec A une matrice creuse
! en coordonnées coo
      implicit none
! variables d'entrée
      integer, intent(in) :: nla,nca,nnz
      integer, intent(in) :: ia(:),ja(:)
      double precision, intent(in) :: A(:)
      double precision, intent(in) :: x(nca)
! variables de sortie
      double precision, intent(inout) :: y(nla)
! variables locales
      integer :: i1

      y(:)=dcmplx(0.0d0,0.0d0)
      do i1=1,nnz
	  !write(*,*) ia(i1),ja(i1)
	 y(ia(i1))=y(ia(i1))+A(i1)*x(ja(i1))
      end do

      end subroutine coo_dmatvec

      subroutine coo_dmatvec_sym(nla,nca,nnz,A,ia,ja,x,y)
! multiplication Ax=y avec A une matrice creuse symétrique
! en coordonnées coo, stockage A matrice supérieure ou inférieure
      implicit none
! variables d'entrée
      integer, intent(in) :: nla,nca,nnz
      integer, intent(in) :: ia(nnz),ja(nnz)
      double precision, intent(in) :: A(nnz)
      double precision, intent(in) :: x(nca)
! variables de sortie
      double precision, intent(inout) :: y(nla)
! variables locales
      integer :: i1

      y(:)=dcmplx(0.0d0,0.0d0)
      do i1=1,nnz
	  !write(*,*) ia(i1),ja(i1)
      y(ia(i1))=y(ia(i1))+A(i1)*x(ja(i1))
      if (ia(i1) /= ja(i1)) then
            y(ja(i1))=y(ja(i1))+A(i1)*x(ia(i1))
      endif
      end do

      end subroutine coo_dmatvec_sym
!------------- end of coo_zmatvec -------------------------------------------
!-----------------------------------------------------------------------

      subroutine dcoocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
!-----------------------------------------------------------------------
      integer nrow,nnz
      double precision a(*),ao(*),x
      integer ir(*),jc(*),jao(*),iao(*)
      integer k0,i,j,iad
!      integer, intent(inout) :: iao(nrow+1)
!-----------------------------------------------------------------------
!  Coordinate     to   Compressed Sparse Row
!-----------------------------------------------------------------------
! converts a matrix that is stored in coordinate format
!  a, ir, jc into a row general sparse ao, jao, iao format.
!
! on entry:
!---------
! nrow	= dimension of the matrix
! nnz	= number of nonzero elements in matrix
! a,
! ir,
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
! 	  the elements, ir(k) = its row number and jc(k) = its column
!	  number. The order of the elements is arbitrary.
!
! on return:
!-----------
! ir 	is destroyed
!
! ao, jao, iao = matrix in general sparse matrix format with ao
! 	continung the real values, jao containing the column indices,
!	and iao being the pointer to the beginning of the row,
!	in arrays ao, jao.
!
! Notes:
!------ This routine is NOT in place.  See coicsr
!
!------------------------------------------------------------------------
      integer :: k

!      write(*,*) 'in zcoocsr'
!      iao=[(0,k=1,nrow+1)]
      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
!      write(*,*) 'raz ia0'
! determine row-lengths.
      do 2 k=1, nnz
            !write(*,*) 'k=',k,' ir(k)=',ir(k)
         iao(ir(k)) = iao(ir(k))+1
 2    continue
!      write(*,*) 'fin allocation iao'
! starting position of each row..
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
! go through the structure  once more. Fill in output matrix.
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
! shift back iao
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
!------------- end of dcoocsr -------------------------------------------
!-----------------------------------------------------------------------
      end
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      subroutine dcsrcoo (nrow,job,nzmax,a,ja,ia,nnz,ao,ir,jc,ierr)
!-----------------------------------------------------------------------
      integer nrow,nnz,job,nzmax

      real*8 a(*),ao(*)
      integer ir(*),jc(*),ja(*),ia(nrow+1),ierr
      integer k,i,k1,k2
!-----------------------------------------------------------------------
!  Compressed Sparse Row      to      Coordinate
!-----------------------------------------------------------------------
! converts a matrix that is stored in coordinate format
!  a, ir, jc into a row general sparse ao, jao, iao format.
!
! on entry:
!---------
! nrow	= dimension of the matrix.
! job   = integer serving as a job indicator.
!         if job = 1 fill in only the array ir, ignore jc, and ao.
!         if job = 2 fill in ir, and jc but not ao
!         if job = 3 fill in everything.
!         The reason why these options are provided is that on return
!         ao and jc are the same as a, ja. So when job = 3, a and ja are
!         simply copied into ao, jc.  When job=2, only jc and ir are
!         returned. With job=1 only the array ir is returned. Moreover,
!         the algorithm is in place:
!	     call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr)
!         will write the output matrix in coordinate format on a, ja,ia.
!
! a,
! ja,
! ia    = matrix in compressed sparse row format.
! nzmax = length of space available in a, ir, jc.
!         the code will stop immediatly if the number of
!         nonzero elements found in input matrix exceeds nzmax.
!
! on return:
!-----------
! ao, ir, jc = matrix in coordinate format.
!
! nnz        = number of nonzero elements in matrix.
! ierr       = integer error indicator.
!         ierr .eq. 0 means normal retur
!         ierr .eq. 1 means that the the code stopped
!         because there was no space in ao, ir, jc
!         (according to the value of  nzmax).
!
! NOTES: 1)This routine is PARTIALLY in place: csrcoo can be called with
!         ao being the same array as as a, and jc the same array as ja.
!         but ir CANNOT be the same as ia.
!         2) note the order in the output arrays,
!------------------------------------------------------------------------
      ierr = 0
      nnz = ia(nrow+1)-1
      if (nnz .gt. nzmax) then
         ierr = 1
         return
      endif
!------------------------------------------------------------------------
      goto (3,2,1) job
 1    do 10 k=1,nnz
         ao(k) = a(k)
 10   continue
 2    do 11 k=1,nnz
         jc(k) = ja(k)
 11   continue
!
!     copy backward to allow for in-place processing.
!
 3    do 13 i=nrow,1,-1
         k1 = ia(i+1)-1
         k2 = ia(i)
         do 12 k=k1,k2,-1
            ir(k) = i
 12      continue
 13   continue
      return
!------------- end-of-csrcoo -------------------------------------------
!-----------------------------------------------------------------------
      end
      subroutine dcsr_sum_duplicates(n, nnzin, Ap, Aj, Ax)
! Sum together duplicate column entries in each row of CSR matrix A
! The column indicies within each row must be in sorted order.
! Explicit zeros are retained.
! Ap, Aj, and Ax will be modified *inplace* ; Ap(n+1)-1 = nouveau nnz.
! Shape explicite (tailles en argument) => appelable sans interface.
      implicit none
      integer, intent(in) :: n, nnzin
      integer, intent(inout) :: Ap(n+1), Aj(nnzin)
      double precision, intent(inout) :: Ax(nnzin)
      integer :: nnz, r1, r2, i, j, jj
      double precision :: x
      nnz = 1
      r2 = 1
      do i = 1, n
            r1 = r2
            r2 = Ap(i+1)
            jj = r1
            do while (jj < r2)
                  j = Aj(jj)
                  x = Ax(jj)
                  jj = jj + 1
                  do while (jj < r2)
                        if (Aj(jj) == j) then
                              x = x + Ax(jj)
                              jj = jj + 1
                        else
                              exit
                        end if
                  end do
                  Aj(nnz) = j
                  Ax(nnz) = x
                  nnz = nnz + 1
            end do
            Ap(i+1) = nnz
      end do
      end subroutine dcsr_sum_duplicates
      subroutine dcsr_sort(n, nnzin, Ap, Aj, Ax)
! Sort CSR column indices inplace (tri par insertion, par ligne — autonome)
! Shape explicite (tailles en argument) => appelable sans interface.
      implicit none
      integer, intent(in) :: n, nnzin
      integer, intent(inout) :: Ap(n+1), Aj(nnzin)
      double precision, intent(inout) :: Ax(nnzin)
      integer :: i, r1, r2, j, k, ctmp
      double precision :: xtmp
      do i = 1, n
            r1 = Ap(i)
            r2 = Ap(i+1)-1
            do j = r1+1, r2
                  ctmp = Aj(j)
                  xtmp = Ax(j)
                  k = j-1
                  do while (k >= r1)
                        if (Aj(k) <= ctmp) exit
                        Aj(k+1) = Aj(k)
                        Ax(k+1) = Ax(k)
                        k = k-1
                  end do
                  Aj(k+1) = ctmp
                  Ax(k+1) = xtmp
            end do
      end do
      end subroutine dcsr_sort


! ==================================================================
! Versions COMPLEXES de tri/somme-doublons CSR (pour near-field T3 avant ILUT)
! Shape explicite (tailles en argument) => appelables sans interface.
! ==================================================================
      subroutine zcsr_sort(n, nnzin, Ap, Aj, Ax)
! Tri des indices colonnes par ligne (insertion), valeurs complexes.
      implicit none
      integer, intent(in) :: n, nnzin
      integer, intent(inout) :: Ap(n+1), Aj(nnzin)
      double complex, intent(inout) :: Ax(nnzin)
      integer :: i, r1, r2, j, k, ctmp
      double complex :: xtmp
      do i = 1, n
            r1 = Ap(i)
            r2 = Ap(i+1)-1
            do j = r1+1, r2
                  ctmp = Aj(j)
                  xtmp = Ax(j)
                  k = j-1
                  do while (k >= r1)
                        if (Aj(k) <= ctmp) exit
                        Aj(k+1) = Aj(k)
                        Ax(k+1) = Ax(k)
                        k = k-1
                  end do
                  Aj(k+1) = ctmp
                  Ax(k+1) = xtmp
            end do
      end do
      end subroutine zcsr_sort

      subroutine zcsr_sum_duplicates(n, nnzin, Ap, Aj, Ax)
! Somme les entrees dupliquees (memes colonnes triees) par ligne, valeurs complexes.
! Ap(n+1)-1 donne le nouveau nnz apres compactage.
      implicit none
      integer, intent(in) :: n, nnzin
      integer, intent(inout) :: Ap(n+1), Aj(nnzin)
      double complex, intent(inout) :: Ax(nnzin)
      integer :: nnz, r1, r2, i, j, jj
      double complex :: x
      nnz = 1
      r2 = 1
      do i = 1, n
            r1 = r2
            r2 = Ap(i+1)
            jj = r1
            do while (jj < r2)
                  j = Aj(jj)
                  x = Ax(jj)
                  jj = jj + 1
                  do while (jj < r2)
                        if (Aj(jj) == j) then
                              x = x + Ax(jj)
                              jj = jj + 1
                        else
                              exit
                        end if
                  end do
                  Aj(nnz) = j
                  Ax(nnz) = x
                  nnz = nnz + 1
            end do
            Ap(i+1) = nnz
      end do
      end subroutine zcsr_sum_duplicates
