      SUBROUTINE LSQ(m,meq,n,nl,la,l,g,a,b,xl,xu,x,y,w,jw,mode)

C   MINIMIZE with respect to X

C             ||E*X - F||
C                                      1/2  T
C   WITH UPPER TRIANGULAR MATRIX E = +D   *L ,

C                                      -1/2  -1
C                     AND VECTOR F = -D    *L  *G,

C  WHERE THE UNIT LOWER TRIDIANGULAR MATRIX L IS STORED COLUMNWISE
C  DENSE IN THE N*(N+1)/2 ARRAY L WITH VECTOR D STORED IN ITS
C 'DIAGONAL' THUS SUBSTITUTING THE ONE-ELEMENTS OF L

C   SUBJECT TO

C             A(J)*X - B(J) = 0 ,         J=1,...,MEQ,
C             A(J)*X - B(J) >=0,          J=MEQ+1,...,M,
C             XL(I) <= X(I) <= XU(I),     I=1,...,N,
C     ON ENTRY, THE USER HAS TO PROVIDE THE ARRAYS L, G, A, B, XL, XU.
C     WITH DIMENSIONS: L(N*(N+1)/2), G(N), A(LA,N), B(M), XL(N), XU(N)
C     THE WORKING ARRAY W MUST HAVE AT LEAST THE FOLLOWING DIMENSION:
c     DIM(W) =        (3*N+M)*(N+1)                        for LSQ
c                    +(N-MEQ+1)*(MINEQ+2) + 2*MINEQ        for LSI
c                    +(N+MINEQ)*(N-MEQ) + 2*MEQ + N        for LSEI
c                      with MINEQ = M - MEQ + 2*N
C     ON RETURN, NO ARRAY WILL BE CHANGED BY THE SUBROUTINE.
C     X     STORES THE N-DIMENSIONAL SOLUTION VECTOR
C     Y     STORES THE VECTOR OF LAGRANGE MULTIPLIERS OF DIMENSION
C           M+N+N (CONSTRAINTS+LOWER+UPPER BOUNDS)
C     MODE  IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING MEANINGS:
C          MODE=1: SUCCESSFUL COMPUTATION
C               2: ERROR RETURN BECAUSE OF WRONG DIMENSIONS (N<1)
C               3: ITERATION COUNT EXCEEDED BY NNLS
C               4: INEQUALITY CONSTRAINTS INCOMPATIBLE
C               5: MATRIX E IS NOT OF FULL RANK
C               6: MATRIX C IS NOT OF FULL RANK
C               7: RANK DEFECT IN HFTI

c     coded            Dieter Kraft, april 1987
c     revised                        march 1989

      DOUBLE PRECISION l,g,a,b,w,xl,xu,x,y,
     .                 diag,ZERO,one,DDOT,xnorm

      INTEGER          jw(*),i,ic,id,ie,IF,ig,ih,il,im,ip,iu,iw,
     .     i1,i2,i3,i4,la,m,meq,mineq,mode,m1,n,nl,n1,n2,n3,
     .     nancnt,j

      DIMENSION        a(la,n), b(la), g(n), l(nl),
     .                 w(*), x(n), xl(n), xu(n), y(m+n+n)

      DATA             ZERO/0.0d0/, one/1.0d0/

      n1 = n + 1
      mineq = m - meq
      m1 = mineq + n + n

c  determine whether to solve problem
c  with inconsistent linerarization (n2=1)
c  or not (n2=0)

      n2 = n1*n/2 + 1
      IF (n2.EQ.nl) THEN
          n2 = 0
      ELSE
          n2 = 1
      ENDIF
      n3 = n-n2

C  RECOVER MATRIX E AND VECTOR F FROM L AND G

      i2 = 1
      i3 = 1
      i4 = 1
      ie = 1
      IF = n*n+1
      DO 10 i=1,n3
         i1 = n1-i
         diag = SQRT (l(i2))
         w(i3) = ZERO
         CALL DCOPY (i1  ,  w(i3), 0, w(i3), 1)
         CALL DCOPY (i1-n2, l(i2), 1, w(i3), n)
         CALL DSCAL (i1-n2,     diag, w(i3), n)
         w(i3) = diag
         w(IF-1+i) = (g(i) - DDOT (i-1, w(i4), 1, w(IF), 1))/diag
         i2 = i2 + i1 - n2
         i3 = i3 + n1
         i4 = i4 + n
   10 CONTINUE
      IF (n2.EQ.1) THEN
          w(i3) = l(nl)
          w(i4) = ZERO
          CALL DCOPY (n3, w(i4), 0, w(i4), 1)
          w(IF-1+n) = ZERO
      ENDIF
      CALL DSCAL (n, - one, w(IF), 1)

      ic = IF + n
      id = ic + meq*n

      IF (meq .GT. 0) THEN

C  RECOVER MATRIX C FROM UPPER PART OF A

          DO 20 i=1,meq
              CALL DCOPY (n, a(i,1), la, w(ic-1+i), meq)
   20     CONTINUE

C  RECOVER VECTOR D FROM UPPER PART OF B

          CALL DCOPY (meq, b(1), 1, w(id), 1)
          CALL DSCAL (meq,   - one, w(id), 1)

      ENDIF

      ig = id + meq

C  RECOVER MATRIX G FROM LOWER PART OF A
C  The matrix G(mineq+2*n,m1) is stored at w(ig)
C  Not all rows will be filled if some of the upper/lower
C  bounds are unbounded.

      IF (mineq .GT. 0) THEN

          DO 30 i=1,mineq
              CALL DCOPY (n, a(meq+i,1), la, w(ig-1+i), m1)
   30     CONTINUE

      ENDIF

      ih = ig + m1*n
      iw = ih + mineq + 2*n

      IF (mineq .GT. 0) THEN

C  RECOVER H FROM LOWER PART OF B
C  The vector H(mineq+2*n) is stored at w(ih)

          CALL DCOPY (mineq, b(meq+1), 1, w(ih), 1)
          CALL DSCAL (mineq,       - one, w(ih), 1)

      ENDIF

C  AUGMENT MATRIX G BY +I AND -I, AND,
C  AUGMENT VECTOR H BY XL AND XU
C  NaN value indicates no bound

      ip = ig + mineq
      il = ih + mineq
      nancnt = 0

      DO 40 i=1,n
         if (xl(i).eq.xl(i)) then
            w(il) = xl(i)
            do 41 j=1,n
               w(ip + m1*(j-1)) = 0
 41         continue
            w(ip + m1*(i-1)) = 1
            ip = ip + 1
            il = il + 1
         else
            nancnt = nancnt + 1
         end if
   40 CONTINUE

      DO 50 i=1,n
         if (xu(i).eq.xu(i)) then
            w(il) = -xu(i)
            do 51 j=1,n
               w(ip + m1*(j-1)) = 0
 51         continue
            w(ip + m1*(i-1)) = -1
            ip = ip + 1
            il = il + 1
         else
            nancnt = nancnt + 1
         end if
 50   CONTINUE

      CALL LSEI (w(ic), w(id), w(ie), w(IF), w(ig), w(ih), MAX(1,meq),
     .           meq, n, n, m1, m1-nancnt, n, x, xnorm, w(iw), jw, mode)

      IF (mode .EQ. 1) THEN

c   restore Lagrange multipliers (only for user-defined variables)

          CALL DCOPY (m,  w(iw),     1, y(1),      1)

c   set rest of the multipliers to nan (they are not used)

          IF (n3 .GT. 0) THEN
             y(m+1) = 0
             y(m+1) = 0 / y(m+1)
             do 60 i=m+2,m+n3+n3
                y(i) = y(m+1)
 60          continue
          ENDIF

      ENDIF
      call bound(n, x, xl, xu)

C   END OF SUBROUTINE LSQ

      END

      subroutine bound(n, x, xl, xu)
      integer n, i
      double precision x(n), xl(n), xu(n)
      do i = 1, n
C        Note that xl(i) and xu(i) may be NaN to indicate no bound
         if(xl(i).eq.xl(i).and.x(i) < xl(i))then
            x(i) = xl(i)
         else if(xu(i).eq.xu(i).and.x(i) > xu(i))then
            x(i) = xu(i)
            end if
         end do
         end subroutine bound