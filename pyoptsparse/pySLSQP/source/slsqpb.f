      SUBROUTINE SLSQPB (m, meq, la, n, x, xl, xu, f, c, g, a, acc,
     *                   iter, mode, r, l, x0, mu, s, u, v, w, iw,
     *                   alpha, f0, gs, h1, h2, h3, h4, t, t0, tol,
     *                   iexact, incons, ireset, itermx, line, 
     *                   n1, n2, n3)

C   NONLINEAR PROGRAMMING BY SOLVING SEQUENTIALLY QUADRATIC PROGRAMS

C        -  L1 - LINE SEARCH,  POSITIVE DEFINITE  BFGS UPDATE  -

C                      BODY SUBROUTINE FOR SLSQP

      INTEGER          iw(*), i, iexact, incons, ireset, iter, itermx,
     *                 k, j, la, line, m, meq, mode, n, n1, n2, n3
      LOGICAL          badlin

      DOUBLE PRECISION a(la,n+1), c(la), g(n+1), l((n+1)*(n+2)/2),
     *                 mu(la), r(m+n+n+2), s(n+1), u(n+1), v(n+1), w(*),
     *                 x(n), xl(n), xu(n), x0(n),
     *                 DDOT, DNRM2, LINMIN,
     *                 acc, alfmin, alpha, f, f0, gs, h1, h2, h3, h4,
     *                 hun, one, t, t0, ten, tol, two, ZERO

c     dim(W) =         N1*(N1+1) + MEQ*(N1+1) + MINEQ*(N1+1)  for LSQ
c                     +(N1-MEQ+1)*(MINEQ+2) + 2*MINEQ
c                     +(N1+MINEQ)*(N1-MEQ) + 2*MEQ + N1       for LSEI
c                      with MINEQ = M - MEQ + 2*N1  &  N1 = N+1

      DATA             ZERO /0.0d0/, one /1.0d0/, alfmin /1.0d-1/,
     *                 hun /1.0d+2/, ten /1.0d+1/, two /2.0d0/

C     The badlin flag keeps track whether the SQP problem on the current
C     iteration was inconsistent or not.
      badlin = .false.

      IF (mode) 260, 100, 220

  100 itermx = iter
      IF (acc.GE.ZERO) THEN
          iexact = 0
      ELSE
          iexact = 1
      ENDIF
      acc = ABS(acc)
      tol = ten*acc
      iter = 0
      ireset = 0
      n1 = n + 1
      n2 = n1*n/2
      n3 = n2 + 1
      s(1) = ZERO
      mu(1) = ZERO
      CALL DCOPY(n, s(1),  0, s,  1)
      CALL DCOPY(m, mu(1), 0, mu, 1)

C   RESET BFGS MATRIX

  110 ireset = ireset + 1
      IF (ireset.GT.5) GO TO 255
      l(1) = ZERO
      CALL DCOPY(n2, l(1), 0, l, 1)
      j = 1
      DO 120 i=1,n
         l(j) = one
         j = j + n1 - i
  120 CONTINUE

C   MAIN ITERATION : SEARCH DIRECTION, STEPLENGTH, LDL'-UPDATE

  130 iter = iter + 1
      mode = 9
      IF (iter.GT.itermx) GO TO 330

C   SEARCH DIRECTION AS SOLUTION OF QP - SUBPROBLEM

      CALL DCOPY(n, xl, 1, u, 1)
      CALL DCOPY(n, xu, 1, v, 1)
      CALL DAXPY(n, -one, x, 1, u, 1)
      CALL DAXPY(n, -one, x, 1, v, 1)
      h4 = one
      CALL LSQ (m, meq, n , n3, la, l, g, a, c, u, v, s, r, w, iw, mode)

C   AUGMENTED PROBLEM FOR INCONSISTENT LINEARIZATION
C
C   If it turns out that the original SQP problem is inconsistent,
C   disallow termination with convergence on this iteration,
C   even if the augmented problem was solved.

      badlin = .false.
      IF (mode.EQ.6) THEN
          IF (n.EQ.meq) THEN
              mode = 4
          ENDIF
      ENDIF
      IF (mode.EQ.4) THEN
          badlin = .true.
          DO 140 j=1,m
             IF (j.LE.meq) THEN
                 a(j,n1) = -c(j)
             ELSE
                 a(j,n1) = MAX(-c(j),ZERO)
             ENDIF
  140     CONTINUE
          s(1) = ZERO
          CALL DCOPY(n, s(1), 0, s, 1)
          h3 = ZERO
          g(n1) = ZERO
          l(n3) = hun
          s(n1) = one
          u(n1) = ZERO
          v(n1) = one
          incons = 0
  150     CALL LSQ (m, meq, n1, n3, la, l, g, a, c, u, v, s, r,
     *              w, iw, mode)
          h4 = one - s(n1)
          IF (mode.EQ.4) THEN
              l(n3) = ten*l(n3)
              incons = incons + 1
              IF (incons.GT.5) GO TO 330
              GOTO 150
          ELSE IF (mode.NE.1) THEN
              GOTO 330
          ENDIF
      ELSE IF (mode.NE.1) THEN
          GOTO 330
      ENDIF

C   UPDATE MULTIPLIERS FOR L1-TEST

      DO 160 i=1,n
         v(i) = g(i) - DDOT(m,a(1,i),1,r,1)
  160 CONTINUE
      f0 = f
      CALL DCOPY(n, x, 1, x0, 1)
      gs = DDOT(n, g, 1, s, 1)
      h1 = ABS(gs)
      h2 = ZERO
      DO 170 j=1,m
         IF (j.LE.meq) THEN
             h3 = c(j)
         ELSE
             h3 = ZERO
         ENDIF
         h2 = h2 + MAX(-c(j),h3)
         h3 = ABS(r(j))
         mu(j) = MAX(h3,(mu(j)+h3)/two)
         h1 = h1 + h3*ABS(c(j))
  170 CONTINUE

C   CHECK CONVERGENCE

      mode = 0
      IF (h1.LT.acc .AND. h2.LT.acc .AND. .NOT. badlin
     *     .AND. f .EQ. f) GO TO 330
      h1 = ZERO
      DO 180 j=1,m
         IF (j.LE.meq) THEN
             h3 = c(j)
         ELSE
             h3 = ZERO
         ENDIF
         h1 = h1 + mu(j)*MAX(-c(j),h3)
  180 CONTINUE
      t0 = f + h1
      h3 = gs - h1*h4
      mode = 8
      IF (h3.GE.ZERO) GO TO 110

C   LINE SEARCH WITH AN L1-TESTFUNCTION

      line = 0
      alpha = one
      IF (iexact.EQ.1) GOTO 210

C   INEXACT LINESEARCH

  190     line = line + 1
          h3 = alpha*h3
          CALL DSCAL(n, alpha, s, 1)
          CALL DCOPY(n, x0, 1, x, 1)
          CALL DAXPY(n, one, s, 1, x, 1)
          mode = 1
          GO TO 330
  200         IF (h1.LE.h3/ten .OR. line.GT.10) GO TO 240
              alpha = MAX(h3/(two*(h3-h1)),alfmin)
              GO TO 190

C   EXACT LINESEARCH

  210 IF (line.NE.3) THEN
          alpha = LINMIN(line,alfmin,one,t,tol)
          CALL DCOPY(n, x0, 1, x, 1)
          CALL DAXPY(n, alpha, s, 1, x, 1)
          mode = 1
          GOTO 330
      ENDIF
      CALL DSCAL(n, alpha, s, 1)
      GOTO 240

C   CALL FUNCTIONS AT CURRENT X

  220     t = f
          DO 230 j=1,m
             IF (j.LE.meq) THEN
                 h1 = c(j)
             ELSE
                 h1 = ZERO
             ENDIF
             t = t + mu(j)*MAX(-c(j),h1)
  230     CONTINUE
          h1 = t - t0
          GOTO (200, 210) iexact+1

C   CHECK CONVERGENCE

  240 h3 = ZERO
      DO 250 j=1,m
         IF (j.LE.meq) THEN
             h1 = c(j)
         ELSE
             h1 = ZERO
         ENDIF
         h3 = h3 + MAX(-c(j),h1)
  250 CONTINUE
      IF ((ABS(f-f0).LT.acc .OR. DNRM2(n,s,1).LT.acc) .AND. h3.LT.acc
     *     .AND. .NOT. badlin .AND. f .EQ. f)
     *   THEN
            mode = 0
         ELSE
            mode = -1
         ENDIF
      GO TO 330

C   CHECK relaxed CONVERGENCE in case of positive directional derivative

  255 CONTINUE
      h3 = ZERO
      DO 256 j=1,m
         IF (j.LE.meq) THEN
             h1 = c(j)
         ELSE
             h1 = ZERO
         ENDIF
         h3 = h3 + MAX(-c(j),h1)
  256 CONTINUE
      IF ((ABS(f-f0).LT.tol .OR. DNRM2(n,s,1).LT.tol) .AND. h3.LT.tol
     *     .AND. .NOT. badlin .AND. f .EQ. f)
     *   THEN
            mode = 0
         ELSE
            mode = 8
         ENDIF
      GO TO 330

C   CALL JACOBIAN AT CURRENT X

C   UPDATE CHOLESKY-FACTORS OF HESSIAN MATRIX BY MODIFIED BFGS FORMULA

  260 DO 270 i=1,n
         u(i) = g(i) - DDOT(m,a(1,i),1,r,1) - v(i)
  270 CONTINUE

C   L'*S

      k = 0
      DO 290 i=1,n
         h1 = ZERO
         k = k + 1
         DO 280 j=i+1,n
            k = k + 1
            h1 = h1 + l(k)*s(j)
  280    CONTINUE
         v(i) = s(i) + h1
  290 CONTINUE

C   D*L'*S

      k = 1
      DO 300 i=1,n
         v(i) = l(k)*v(i)
         k = k + n1 - i
  300 CONTINUE

C   L*D*L'*S

      DO 320 i=n,1,-1
         h1 = ZERO
         k = i
         DO 310 j=1,i - 1
            h1 = h1 + l(k)*v(j)
            k = k + n - j
  310    CONTINUE
         v(i) = v(i) + h1
  320 CONTINUE

      h1 = DDOT(n,s,1,u,1)
      h2 = DDOT(n,s,1,v,1)
      h3 = 0.2d0*h2
      IF (h1.LT.h3) THEN
          h4 = (h2-h3)/(h2-h1)
          h1 = h3
          CALL DSCAL(n, h4, u, 1)
          CALL DAXPY(n, one-h4, v, 1, u, 1)
      ENDIF
      IF (h1.EQ.0 .or. h2.EQ.0) THEN
C         Singular update: reset hessian.
          GO TO 110
      end if
      CALL ldl(n, l, u, +one/h1, v)
      CALL ldl(n, l, v, -one/h2, u)

C   END OF MAIN ITERATION

      GO TO 130

C   END OF SLSQPB

  330 END
      