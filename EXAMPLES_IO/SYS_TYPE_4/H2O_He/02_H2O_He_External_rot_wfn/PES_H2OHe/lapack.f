C    -------------   BELOW IS DSYEVX  --------------------
CAT > DSYEVX.F <<'CUT HERE............'
      SUBROUTINE DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
     $                   ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,
     $                   IFAIL, INFO )
*
*  -- LAPACK DRIVER ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IFAIL( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DSYEVX COMPUTES SELECTED EIGENVALUES AND, OPTIONALLY, EIGENVECTORS
*  OF A REAL SYMMETRIC MATRIX A.  EIGENVALUES AND EIGENVECTORS CAN BE
*  SELECTED BY SPECIFYING EITHER A RANGE OF VALUES OR A RANGE OF INDICES
*  FOR THE DESIRED EIGENVALUES.
*
*  ARGUMENTS
*  =========
*
*  JOBZ    (INPUT) CHARACTER*1
*          = 'N':  COMPUTE EIGENVALUES ONLY;
*          = 'V':  COMPUTE EIGENVALUES AND EIGENVECTORS.
*
*  RANGE   (INPUT) CHARACTER*1
*          = 'A': ALL EIGENVALUES WILL BE FOUND.
*          = 'V': ALL EIGENVALUES IN THE HALF-OPEN INTERVAL (VL,VU]
*                 WILL BE FOUND.
*          = 'I': THE IL-TH THROUGH IU-TH EIGENVALUES WILL BE FOUND.
*
*  UPLO    (INPUT) CHARACTER*1
*          = 'U':  UPPER TRIANGLE OF A IS STORED;
*          = 'L':  LOWER TRIANGLE OF A IS STORED.
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (LDA, N)
*          ON ENTRY, THE SYMMETRIC MATRIX A.  IF UPLO = 'U', THE
*          LEADING N-BY-N UPPER TRIANGULAR PART OF A CONTAINS THE
*          UPPER TRIANGULAR PART OF THE MATRIX A.  IF UPLO = 'L',
*          THE LEADING N-BY-N LOWER TRIANGULAR PART OF A CONTAINS
*          THE LOWER TRIANGULAR PART OF THE MATRIX A.
*          ON EXIT, THE LOWER TRIANGLE (IF UPLO='L') OR THE UPPER
*          TRIANGLE (IF UPLO='U') OF A, INCLUDING THE DIAGONAL, IS
*          DESTROYED.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,N).
*
*  VL      (INPUT) DOUBLE PRECISION
*          IF RANGE='V', THE LOWER BOUND OF THE INTERVAL TO BE SEARCHED
*          FOR EIGENVALUES.  NOT REFERENCED IF RANGE = 'A' OR 'I'.
*
*  VU      (INPUT) DOUBLE PRECISION
*          IF RANGE='V', THE UPPER BOUND OF THE INTERVAL TO BE SEARCHED
*          FOR EIGENVALUES.  NOT REFERENCED IF RANGE = 'A' OR 'I'.
*
*  IL      (INPUT) INTEGER
*          IF RANGE='I', THE INDEX (FROM SMALLEST TO LARGEST) OF THE
*          SMALLEST EIGENVALUE TO BE RETURNED.  IL >= 1.
*          NOT REFERENCED IF RANGE = 'A' OR 'V'.
*
*  IU      (INPUT) INTEGER
*          IF RANGE='I', THE INDEX (FROM SMALLEST TO LARGEST) OF THE
*          LARGEST EIGENVALUE TO BE RETURNED.  MIN(IL,N) <= IU <= N.
*          NOT REFERENCED IF RANGE = 'A' OR 'V'.
*
*  ABSTOL  (INPUT) DOUBLE PRECISION
*          THE ABSOLUTE ERROR TOLERANCE FOR THE EIGENVALUES.
*          AN APPROXIMATE EIGENVALUE IS ACCEPTED AS CONVERGED
*          WHEN IT IS DETERMINED TO LIE IN AN INTERVAL [A,B]
*          OF WIDTH LESS THAN OR EQUAL TO
*
*                  ABSTOL + EPS *   MAX( |A|,|B| ) ,
*
*          WHERE EPS IS THE MACHINE PRECISION.  IF ABSTOL IS LESS THAN
*          OR EQUAL TO ZERO, THEN  EPS*|T|  WILL BE USED IN ITS PLACE,
*          WHERE |T| IS THE 1-NORM OF THE TRIDIAGONAL MATRIX OBTAINED
*          BY REDUCING A TO TRIDIAGONAL FORM.
*
*          SEE "COMPUTING SMALL SINGULAR VALUES OF BIDIAGONAL MATRICES
*          WITH GUARANTEED HIGH RELATIVE ACCURACY," BY DEMMEL AND
*          KAHAN, LAPACK WORKING NOTE #3.
*
*  M       (OUTPUT) INTEGER
*          THE TOTAL NUMBER OF EIGENVALUES FOUND.  0 <= M <= N.
*          IF RANGE = 'A', M = N, AND IF RANGE = 'I', M = IU-IL+1.
*
*  W       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          ON NORMAL EXIT, THE FIRST M ENTRIES CONTAIN THE SELECTED
*          EIGENVALUES IN ASCENDING ORDER.
*
*  Z       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDZ, MAX(1,M))
*          IF JOBZ = 'V', THEN IF INFO = 0, THE FIRST M COLUMNS OF Z
*          CONTAIN THE ORTHONORMAL EIGENVECTORS OF THE MATRIX
*          CORRESPONDING TO THE SELECTED EIGENVALUES.  IF AN EIGENVECTOR
*          FAILS TO CONVERGE, THEN THAT COLUMN OF Z CONTAINS THE LATEST
*          APPROXIMATION TO THE EIGENVECTOR, AND THE INDEX OF THE
*          EIGENVECTOR IS RETURNED IN IFAIL.
*          IF JOBZ = 'N', THEN Z IS NOT REFERENCED.
*          NOTE: THE USER MUST ENSURE THAT AT LEAST MAX(1,M) COLUMNS ARE
*          SUPPLIED IN THE ARRAY Z; IF RANGE = 'V', THE EXACT VALUE OF M
*          IS NOT KNOWN IN ADVANCE AND AN UPPER BOUND MUST BE USED.
*
*  LDZ     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY Z.  LDZ >= 1, AND IF
*          JOBZ = 'V', LDZ >= MAX(1,N).
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (LWORK)
*          ON EXIT, IF INFO = 0, WORK(1) RETURNS THE OPTIMAL LWORK.
*
*  LWORK   (INPUT) INTEGER
*          THE LENGTH OF THE ARRAY WORK.  LWORK >= MAX(1,8*N).
*          FOR OPTIMAL EFFICIENCY, LWORK >= (NB+3)*N,
*          WHERE NB IS THE BLOCKSIZE FOR DSYTRD RETURNED BY ILAENV.
*
*  IWORK   (WORKSPACE) INTEGER ARRAY, DIMENSION (5*N)
*
*  IFAIL   (OUTPUT) INTEGER ARRAY, DIMENSION (N)
*          IF JOBZ = 'V', THEN IF INFO = 0, THE FIRST M ELEMENTS OF
*          IFAIL ARE ZERO.  IF INFO > 0, THEN IFAIL CONTAINS THE
*          INDICES OF THE EIGENVECTORS THAT FAILED TO CONVERGE.
*          IF JOBZ = 'N', THEN IFAIL IS NOT REFERENCED.
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0:  IF INFO = I, THEN I EIGENVECTORS FAILED TO CONVERGE.
*                THEIR INDICES ARE STORED IN ARRAY IFAIL.
*
* =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            ALLEIG, INDEIG, LOWER, VALEIG, WANTZ
      CHARACTER          ORDER
      INTEGER            I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL,
     $                   INDISP, INDIWO, INDTAU, INDWKN, INDWRK, ISCALE,
     $                   ITMP1, J, JJ, LLWORK, LLWRKN, LOPT, NSPLIT
      DOUBLE PRECISION   ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN,
     $                   SIGMA, SMLNUM, TMP1, VLL, VUU
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           LSAME, DLAMCH, DLANSY
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DCOPY, DLACPY, DORGTR, DORMTR, DSCAL, DSTEBZ,
     $                   DSTEIN, DSTEQR, DSTERF, DSWAP, DSYTRD, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      LOWER = LSAME( UPLO, 'L' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( VALEIG .AND. N.GT.0 .AND. VU.LE.VL ) THEN
         INFO = -8
      ELSE IF( INDEIG .AND. IL.LT.1 ) THEN
         INFO = -9
      ELSE IF( INDEIG .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) ) THEN
         INFO = -10
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -15
      ELSE IF( LWORK.LT.MAX( 1, 8*N ) ) THEN
         INFO = -17
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYEVX', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      M = 0
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( N.EQ.1 ) THEN
         WORK( 1 ) = 7
         IF( ALLEIG .OR. INDEIG ) THEN
            M = 1
            W( 1 ) = A( 1, 1 )
         ELSE
            IF( VL.LT.A( 1, 1 ) .AND. VU.GE.A( 1, 1 ) ) THEN
               M = 1
               W( 1 ) = A( 1, 1 )
            END IF
         END IF
         IF( WANTZ )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     GET MACHINE CONSTANTS.
*
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      EPS = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
*     SCALE MATRIX TO ALLOWABLE RANGE, IF NECESSARY.
*
      ISCALE = 0
      ABSTLL = ABSTOL
      IF( VALEIG ) THEN
         VLL = VL
         VUU = VU
      END IF
      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         IF( LOWER ) THEN
            DO 10 J = 1, N
               CALL DSCAL( N-J+1, SIGMA, A( J, J ), 1 )
   10       CONTINUE
         ELSE
            DO 20 J = 1, N
               CALL DSCAL( J, SIGMA, A( 1, J ), 1 )
   20       CONTINUE
         END IF
         IF( ABSTOL.GT.0 )
     $      ABSTLL = ABSTOL*SIGMA
         IF( VALEIG ) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         END IF
      END IF
*
*     CALL DSYTRD TO REDUCE SYMMETRIC MATRIX TO TRIDIAGONAL FORM.
*
      INDTAU = 1
      INDE = INDTAU + N
      INDD = INDE + N
      INDWRK = INDD + N
      LLWORK = LWORK - INDWRK + 1
      CALL DSYTRD( UPLO, N, A, LDA, WORK( INDD ), WORK( INDE ),
     $             WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO )
      LOPT = 3*N + WORK( INDWRK )
*
*     IF ALL EIGENVALUES ARE DESIRED AND ABSTOL IS LESS THAN OR EQUAL TO
*     ZERO, THEN CALL DSTERF OR DORGTR AND SSTEQR.  IF THIS FAILS FOR
*     SOME EIGENVALUE, THEN TRY DSTEBZ.
*
      IF( ( ALLEIG .OR. ( INDEIG .AND. IL.EQ.1 .AND. IU.EQ.N ) ) .AND.
     $    ( ABSTOL.LE.ZERO ) ) THEN
         CALL DCOPY( N, WORK( INDD ), 1, W, 1 )
         INDEE = INDWRK + 2*N
         IF( .NOT.WANTZ ) THEN
            CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL DSTERF( N, W, WORK( INDEE ), INFO )
         ELSE
            CALL DLACPY( 'A', N, N, A, LDA, Z, LDZ )
            CALL DORGTR( UPLO, N, Z, LDZ, WORK( INDTAU ),
     $                   WORK( INDWRK ), LLWORK, IINFO )
            CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL DSTEQR( JOBZ, N, W, WORK( INDEE ), Z, LDZ,
     $                   WORK( INDWRK ), INFO )
            IF( INFO.EQ.0 ) THEN
               DO 30 I = 1, N
                  IFAIL( I ) = 0
   30          CONTINUE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            M = N
            GO TO 40
         END IF
         INFO = 0
      END IF
*
*     OTHERWISE, CALL DSTEBZ AND, IF EIGENVECTORS ARE DESIRED, SSTEIN.
*
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWO = INDISP + N
      CALL DSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL,
     $             WORK( INDD ), WORK( INDE ), M, NSPLIT, W,
     $             IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWRK ),
     $             IWORK( INDIWO ), INFO )
*
      IF( WANTZ ) THEN
         CALL DSTEIN( N, WORK( INDD ), WORK( INDE ), M, W,
     $                IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ,
     $                WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO )
*
*        APPLY ORTHOGONAL MATRIX USED IN REDUCTION TO TRIDIAGONAL
*        FORM TO EIGENVECTORS RETURNED BY DSTEIN.
*
         INDWKN = INDE
         LLWRKN = LWORK - INDWKN + 1
         CALL DORMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z,
     $                LDZ, WORK( INDWKN ), LLWRKN, IINFO )
      END IF
*
*     IF MATRIX WAS SCALED, THEN RESCALE EIGENVALUES APPROPRIATELY.
*
   40 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = M
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
*     IF EIGENVALUES ARE NOT IN ORDER, THEN SORT THEM, ALONG WITH
*     EIGENVECTORS.
*
      IF( WANTZ ) THEN
         DO 60 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 50 JJ = J + 1, M
               IF( W( JJ ).LT.TMP1 ) THEN
                  I = JJ
                  TMP1 = W( JJ )
               END IF
   50       CONTINUE
*
            IF( I.NE.0 ) THEN
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
               IF( INFO.NE.0 ) THEN
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               END IF
            END IF
   60    CONTINUE
      END IF
*
*     SET WORK(1) TO OPTIMAL WORKSPACE SIZE.
*
      WORK( 1 ) = MAX( 7*N, LOPT )
*
      RETURN
*
*     END OF DSYEVX
*
      END
CUT HERE............
CAT > DORMTR.F <<'CUT HERE............'
      SUBROUTINE DORMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          SIDE, TRANS, UPLO
      INTEGER            INFO, LDA, LDC, LWORK, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ),
     $                   WORK( LWORK )
*     ..
*
*  PURPOSE
*  =======
*
*  DORMTR OVERWRITES THE GENERAL REAL M-BY-N MATRIX C WITH
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'T':      Q**T * C       C * Q**T
*
*  WHERE Q IS A REAL ORTHOGONAL MATRIX OF ORDER NQ, WITH NQ = M IF
*  SIDE = 'L' AND NQ = N IF SIDE = 'R'. Q IS DEFINED AS THE PRODUCT OF
*  NQ-1 ELEMENTARY REFLECTORS, AS RETURNED BY DSYTRD:
*
*  IF UPLO = 'U', Q = H(NQ-1) . . . H(2) H(1);
*
*  IF UPLO = 'L', Q = H(1) H(2) . . . H(NQ-1).
*
*  ARGUMENTS
*  =========
*
*  SIDE    (INPUT) CHARACTER*1
*          = 'L': APPLY Q OR Q**T FROM THE LEFT;
*          = 'R': APPLY Q OR Q**T FROM THE RIGHT.
*
*  UPLO    (INPUT) CHARACTER*1
*          = 'U': UPPER TRIANGLE OF A CONTAINS ELEMENTARY REFLECTORS
*                 FROM DSYTRD;
*          = 'L': LOWER TRIANGLE OF A CONTAINS ELEMENTARY REFLECTORS
*                 FROM DSYTRD.
*
*  TRANS   (INPUT) CHARACTER*1
*          = 'N':  NO TRANSPOSE, APPLY Q;
*          = 'T':  TRANSPOSE, APPLY Q**T.
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX C. M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX C. N >= 0.
*
*  A       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION
*                               (LDA,M) IF SIDE = 'L'
*                               (LDA,N) IF SIDE = 'R'
*          THE VECTORS WHICH DEFINE THE ELEMENTARY REFLECTORS, AS
*          RETURNED BY DSYTRD.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.
*          LDA >= MAX(1,M) IF SIDE = 'L'; LDA >= MAX(1,N) IF SIDE = 'R'.
*
*  TAU     (INPUT) DOUBLE PRECISION ARRAY, DIMENSION
*                               (M-1) IF SIDE = 'L'
*                               (N-1) IF SIDE = 'R'
*          TAU(I) MUST CONTAIN THE SCALAR FACTOR OF THE ELEMENTARY
*          REFLECTOR H(I), AS RETURNED BY DSYTRD.
*
*  C       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDC,N)
*          ON ENTRY, THE M-BY-N MATRIX C.
*          ON EXIT, C IS OVERWRITTEN BY Q*C OR Q**T*C OR C*Q**T OR C*Q.
*
*  LDC     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY C. LDC >= MAX(1,M).
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (LWORK)
*          ON EXIT, IF INFO = 0, WORK(1) RETURNS THE OPTIMAL LWORK.
*
*  LWORK   (INPUT) INTEGER
*          THE DIMENSION OF THE ARRAY WORK.
*          IF SIDE = 'L', LWORK >= MAX(1,N);
*          IF SIDE = 'R', LWORK >= MAX(1,M).
*          FOR OPTIMUM PERFORMANCE LWORK >= N*NB IF SIDE = 'L', AND
*          LWORK >= M*NB IF SIDE = 'R', WHERE NB IS THE OPTIMAL
*          BLOCKSIZE.
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*
*  =====================================================================
*
*     .. LOCAL SCALARS ..
      LOGICAL            LEFT, UPPER
      INTEGER            I1, I2, IINFO, MI, NI, NQ, NW
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DORMQL, DORMQR, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT ARGUMENTS
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      UPPER = LSAME( UPLO, 'U' )
*
*     NQ IS THE ORDER OF Q AND NW IS THE MINIMUM DIMENSION OF WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'T' ) )
     $          THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMTR', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. NQ.EQ.1 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( LEFT ) THEN
         MI = M - 1
         NI = N
      ELSE
         MI = M
         NI = N - 1
      END IF
*
      IF( UPPER ) THEN
*
*        Q WAS DETERMINED BY A CALL TO DSYTRD WITH UPLO = 'U'
*
         CALL DORMQL( SIDE, TRANS, MI, NI, NQ-1, A( 1, 2 ), LDA, TAU, C,
     $                LDC, WORK, LWORK, IINFO )
      ELSE
*
*        Q WAS DETERMINED BY A CALL TO DSYTRD WITH UPLO = 'L'
*
         IF( LEFT ) THEN
            I1 = 2
            I2 = 1
         ELSE
            I1 = 1
            I2 = 2
         END IF
         CALL DORMQR( SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU,
     $                C( I1, I2 ), LDC, WORK, LWORK, IINFO )
      END IF
      RETURN
*
*     END OF DORMTR
*
      END
CUT HERE............
CAT > DORMQR.F <<'CUT HERE............'
      SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ),
     $                   WORK( LWORK )
*     ..
*
*  PURPOSE
*  =======
*
*  DORMQR OVERWRITES THE GENERAL REAL M-BY-N MATRIX C WITH
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'T':      Q**T * C       C * Q**T
*
*  WHERE Q IS A REAL ORTHOGONAL MATRIX DEFINED AS THE PRODUCT OF K
*  ELEMENTARY REFLECTORS
*
*        Q = H(1) H(2) . . . H(K)
*
*  AS RETURNED BY DGEQRF. Q IS OF ORDER M IF SIDE = 'L' AND OF ORDER N
*  IF SIDE = 'R'.
*
*  ARGUMENTS
*  =========
*
*  SIDE    (INPUT) CHARACTER*1
*          = 'L': APPLY Q OR Q**T FROM THE LEFT;
*          = 'R': APPLY Q OR Q**T FROM THE RIGHT.
*
*  TRANS   (INPUT) CHARACTER*1
*          = 'N':  NO TRANSPOSE, APPLY Q;
*          = 'T':  TRANSPOSE, APPLY Q**T.
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX C. M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX C. N >= 0.
*
*  K       (INPUT) INTEGER
*          THE NUMBER OF ELEMENTARY REFLECTORS WHOSE PRODUCT DEFINES
*          THE MATRIX Q.
*          IF SIDE = 'L', M >= K >= 0;
*          IF SIDE = 'R', N >= K >= 0.
*
*  A       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,K)
*          THE I-TH COLUMN MUST CONTAIN THE VECTOR WHICH DEFINES THE
*          ELEMENTARY REFLECTOR H(I), FOR I = 1,2,...,K, AS RETURNED BY
*          DGEQRF IN THE FIRST K COLUMNS OF ITS ARRAY ARGUMENT A.
*          A IS MODIFIED BY THE ROUTINE BUT RESTORED ON EXIT.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.
*          IF SIDE = 'L', LDA >= MAX(1,M);
*          IF SIDE = 'R', LDA >= MAX(1,N).
*
*  TAU     (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (K)
*          TAU(I) MUST CONTAIN THE SCALAR FACTOR OF THE ELEMENTARY
*          REFLECTOR H(I), AS RETURNED BY DGEQRF.
*
*  C       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDC,N)
*          ON ENTRY, THE M-BY-N MATRIX C.
*          ON EXIT, C IS OVERWRITTEN BY Q*C OR Q**T*C OR C*Q**T OR C*Q.
*
*  LDC     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY C. LDC >= MAX(1,M).
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (LWORK)
*          ON EXIT, IF INFO = 0, WORK(1) RETURNS THE OPTIMAL LWORK.
*
*  LWORK   (INPUT) INTEGER
*          THE DIMENSION OF THE ARRAY WORK.
*          IF SIDE = 'L', LWORK >= MAX(1,N);
*          IF SIDE = 'R', LWORK >= MAX(1,M).
*          FOR OPTIMUM PERFORMANCE LWORK >= N*NB IF SIDE = 'L', AND
*          LWORK >= M*NB IF SIDE = 'R', WHERE NB IS THE OPTIMAL
*          BLOCKSIZE.
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*
*  =====================================================================
*
*     .. PARAMETERS ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK,
     $                   MI, NB, NBMIN, NI, NQ, NW
*     ..
*     .. LOCAL ARRAYS ..
      DOUBLE PRECISION   T( LDT, NBMAX )
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLARFB, DLARFT, DORM2R, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT ARGUMENTS
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ IS THE ORDER OF Q AND NW IS THE MINIMUM DIMENSION OF WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMQR', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     DETERMINE THE BLOCK SIZE.  NB MAY BE AT MOST NBMAX, WHERE NBMAX
*     IS USED TO DEFINE THE LOCAL ARRAY T.
*
      NB = MIN( NBMAX, ILAENV( 1, 'DORMQR', SIDE // TRANS, M, N, K,
     $     -1 ) )
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMQR', SIDE // TRANS, M, N, K,
     $              -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
*
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
*
*        USE UNBLOCKED CODE
*
         CALL DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,
     $                IINFO )
      ELSE
*
*        USE BLOCKED CODE
*
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR.
     $       ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
*
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
*
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
*
*           FORM THE TRIANGULAR FACTOR OF THE BLOCK REFLECTOR
*           H = H(I) H(I+1) . . . H(I+IB-1)
*
            CALL DLARFT( 'FORWARD', 'COLUMNWISE', NQ-I+1, IB, A( I, I ),
     $                   LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
*
*              H OR H' IS APPLIED TO C(I:M,1:N)
*
               MI = M - I + 1
               IC = I
            ELSE
*
*              H OR H' IS APPLIED TO C(1:M,I:N)
*
               NI = N - I + 1
               JC = I
            END IF
*
*           APPLY H OR H'
*
            CALL DLARFB( SIDE, TRANS, 'FORWARD', 'COLUMNWISE', MI, NI,
     $                   IB, A( I, I ), LDA, T, LDT, C( IC, JC ), LDC,
     $                   WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = IWS
      RETURN
*
*     END OF DORMQR
*
      END
CUT HERE............
CAT > DORM2R.F <<'CUT HERE............'
      SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 29, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DORM2R OVERWRITES THE GENERAL REAL M BY N MATRIX C WITH
*
*        Q * C  IF SIDE = 'L' AND TRANS = 'N', OR
*
*        Q'* C  IF SIDE = 'L' AND TRANS = 'T', OR
*
*        C * Q  IF SIDE = 'R' AND TRANS = 'N', OR
*
*        C * Q' IF SIDE = 'R' AND TRANS = 'T',
*
*  WHERE Q IS A REAL ORTHOGONAL MATRIX DEFINED AS THE PRODUCT OF K
*  ELEMENTARY REFLECTORS
*
*        Q = H(1) H(2) . . . H(K)
*
*  AS RETURNED BY DGEQRF. Q IS OF ORDER M IF SIDE = 'L' AND OF ORDER N
*  IF SIDE = 'R'.
*
*  ARGUMENTS
*  =========
*
*  SIDE    (INPUT) CHARACTER*1
*          = 'L': APPLY Q OR Q' FROM THE LEFT
*          = 'R': APPLY Q OR Q' FROM THE RIGHT
*
*  TRANS   (INPUT) CHARACTER*1
*          = 'N': APPLY Q  (NO TRANSPOSE)
*          = 'T': APPLY Q' (TRANSPOSE)
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX C. M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX C. N >= 0.
*
*  K       (INPUT) INTEGER
*          THE NUMBER OF ELEMENTARY REFLECTORS WHOSE PRODUCT DEFINES
*          THE MATRIX Q.
*          IF SIDE = 'L', M >= K >= 0;
*          IF SIDE = 'R', N >= K >= 0.
*
*  A       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,K)
*          THE I-TH COLUMN MUST CONTAIN THE VECTOR WHICH DEFINES THE
*          ELEMENTARY REFLECTOR H(I), FOR I = 1,2,...,K, AS RETURNED BY
*          DGEQRF IN THE FIRST K COLUMNS OF ITS ARRAY ARGUMENT A.
*          A IS MODIFIED BY THE ROUTINE BUT RESTORED ON EXIT.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.
*          IF SIDE = 'L', LDA >= MAX(1,M);
*          IF SIDE = 'R', LDA >= MAX(1,N).
*
*  TAU     (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (K)
*          TAU(I) MUST CONTAIN THE SCALAR FACTOR OF THE ELEMENTARY
*          REFLECTOR H(I), AS RETURNED BY DGEQRF.
*
*  C       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDC,N)
*          ON ENTRY, THE M BY N MATRIX C.
*          ON EXIT, C IS OVERWRITTEN BY Q*C OR Q'*C OR C*Q' OR C*Q.
*
*  LDC     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY C. LDC >= MAX(1,M).
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION
*                                   (N) IF SIDE = 'L',
*                                   (M) IF SIDE = 'R'
*
*  INFO    (OUTPUT) INTEGER
*          = 0: SUCCESSFUL EXIT
*          < 0: IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      DOUBLE PRECISION   AII
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLARF, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT ARGUMENTS
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ IS THE ORDER OF Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2R', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) )
     $     THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(I) IS APPLIED TO C(I:M,1:N)
*
            MI = M - I + 1
            IC = I
         ELSE
*
*           H(I) IS APPLIED TO C(1:M,I:N)
*
            NI = N - I + 1
            JC = I
         END IF
*
*        APPLY H(I)
*
         AII = A( I, I )
         A( I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ),
     $               LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
*
*     END OF DORM2R
*
      END
CUT HERE............
CAT > DORMQL.F <<'CUT HERE............'
      SUBROUTINE DORMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ),
     $                   WORK( LWORK )
*     ..
*
*  PURPOSE
*  =======
*
*  DORMQL OVERWRITES THE GENERAL REAL M-BY-N MATRIX C WITH
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'T':      Q**T * C       C * Q**T
*
*  WHERE Q IS A REAL ORTHOGONAL MATRIX DEFINED AS THE PRODUCT OF K
*  ELEMENTARY REFLECTORS
*
*        Q = H(K) . . . H(2) H(1)
*
*  AS RETURNED BY DGEQLF. Q IS OF ORDER M IF SIDE = 'L' AND OF ORDER N
*  IF SIDE = 'R'.
*
*  ARGUMENTS
*  =========
*
*  SIDE    (INPUT) CHARACTER*1
*          = 'L': APPLY Q OR Q**T FROM THE LEFT;
*          = 'R': APPLY Q OR Q**T FROM THE RIGHT.
*
*  TRANS   (INPUT) CHARACTER*1
*          = 'N':  NO TRANSPOSE, APPLY Q;
*          = 'T':  TRANSPOSE, APPLY Q**T.
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX C. M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX C. N >= 0.
*
*  K       (INPUT) INTEGER
*          THE NUMBER OF ELEMENTARY REFLECTORS WHOSE PRODUCT DEFINES
*          THE MATRIX Q.
*          IF SIDE = 'L', M >= K >= 0;
*          IF SIDE = 'R', N >= K >= 0.
*
*  A       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,K)
*          THE I-TH COLUMN MUST CONTAIN THE VECTOR WHICH DEFINES THE
*          ELEMENTARY REFLECTOR H(I), FOR I = 1,2,...,K, AS RETURNED BY
*          DGEQLF IN THE LAST K COLUMNS OF ITS ARRAY ARGUMENT A.
*          A IS MODIFIED BY THE ROUTINE BUT RESTORED ON EXIT.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.
*          IF SIDE = 'L', LDA >= MAX(1,M);
*          IF SIDE = 'R', LDA >= MAX(1,N).
*
*  TAU     (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (K)
*          TAU(I) MUST CONTAIN THE SCALAR FACTOR OF THE ELEMENTARY
*          REFLECTOR H(I), AS RETURNED BY DGEQLF.
*
*  C       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDC,N)
*          ON ENTRY, THE M-BY-N MATRIX C.
*          ON EXIT, C IS OVERWRITTEN BY Q*C OR Q**T*C OR C*Q**T OR C*Q.
*
*  LDC     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY C. LDC >= MAX(1,M).
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (LWORK)
*          ON EXIT, IF INFO = 0, WORK(1) RETURNS THE OPTIMAL LWORK.
*
*  LWORK   (INPUT) INTEGER
*          THE DIMENSION OF THE ARRAY WORK.
*          IF SIDE = 'L', LWORK >= MAX(1,N);
*          IF SIDE = 'R', LWORK >= MAX(1,M).
*          FOR OPTIMUM PERFORMANCE LWORK >= N*NB IF SIDE = 'L', AND
*          LWORK >= M*NB IF SIDE = 'R', WHERE NB IS THE OPTIMAL
*          BLOCKSIZE.
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*
*  =====================================================================
*
*     .. PARAMETERS ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IINFO, IWS, LDWORK, MI, NB,
     $                   NBMIN, NI, NQ, NW
*     ..
*     .. LOCAL ARRAYS ..
      DOUBLE PRECISION   T( LDT, NBMAX )
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLARFB, DLARFT, DORM2L, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT ARGUMENTS
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ IS THE ORDER OF Q AND NW IS THE MINIMUM DIMENSION OF WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMQL', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     DETERMINE THE BLOCK SIZE.  NB MAY BE AT MOST NBMAX, WHERE NBMAX
*     IS USED TO DEFINE THE LOCAL ARRAY T.
*
      NB = MIN( NBMAX, ILAENV( 1, 'DORMQL', SIDE // TRANS, M, N, K,
     $     -1 ) )
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DORMQL', SIDE // TRANS, M, N, K,
     $              -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
*
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
*
*        USE UNBLOCKED CODE
*
         CALL DORM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,
     $                IINFO )
      ELSE
*
*        USE BLOCKED CODE
*
         IF( ( LEFT .AND. NOTRAN ) .OR.
     $       ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
*
         IF( LEFT ) THEN
            NI = N
         ELSE
            MI = M
         END IF
*
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
*
*           FORM THE TRIANGULAR FACTOR OF THE BLOCK REFLECTOR
*           H = H(I+IB-1) . . . H(I+1) H(I)
*
            CALL DLARFT( 'BACKWARD', 'COLUMNWISE', NQ-K+I+IB-1, IB,
     $                   A( 1, I ), LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
*
*              H OR H' IS APPLIED TO C(1:M-K+I+IB-1,1:N)
*
               MI = M - K + I + IB - 1
            ELSE
*
*              H OR H' IS APPLIED TO C(1:M,1:N-K+I+IB-1)
*
               NI = N - K + I + IB - 1
            END IF
*
*           APPLY H OR H'
*
            CALL DLARFB( SIDE, TRANS, 'BACKWARD', 'COLUMNWISE', MI, NI,
     $                   IB, A( 1, I ), LDA, T, LDT, C, LDC, WORK,
     $                   LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = IWS
      RETURN
*
*     END OF DORMQL
*
      END
CUT HERE............
CAT > DORM2L.F <<'CUT HERE............'
      SUBROUTINE DORM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 29, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DORM2L OVERWRITES THE GENERAL REAL M BY N MATRIX C WITH
*
*        Q * C  IF SIDE = 'L' AND TRANS = 'N', OR
*
*        Q'* C  IF SIDE = 'L' AND TRANS = 'T', OR
*
*        C * Q  IF SIDE = 'R' AND TRANS = 'N', OR
*
*        C * Q' IF SIDE = 'R' AND TRANS = 'T',
*
*  WHERE Q IS A REAL ORTHOGONAL MATRIX DEFINED AS THE PRODUCT OF K
*  ELEMENTARY REFLECTORS
*
*        Q = H(K) . . . H(2) H(1)
*
*  AS RETURNED BY DGEQLF. Q IS OF ORDER M IF SIDE = 'L' AND OF ORDER N
*  IF SIDE = 'R'.
*
*  ARGUMENTS
*  =========
*
*  SIDE    (INPUT) CHARACTER*1
*          = 'L': APPLY Q OR Q' FROM THE LEFT
*          = 'R': APPLY Q OR Q' FROM THE RIGHT
*
*  TRANS   (INPUT) CHARACTER*1
*          = 'N': APPLY Q  (NO TRANSPOSE)
*          = 'T': APPLY Q' (TRANSPOSE)
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX C. M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX C. N >= 0.
*
*  K       (INPUT) INTEGER
*          THE NUMBER OF ELEMENTARY REFLECTORS WHOSE PRODUCT DEFINES
*          THE MATRIX Q.
*          IF SIDE = 'L', M >= K >= 0;
*          IF SIDE = 'R', N >= K >= 0.
*
*  A       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,K)
*          THE I-TH COLUMN MUST CONTAIN THE VECTOR WHICH DEFINES THE
*          ELEMENTARY REFLECTOR H(I), FOR I = 1,2,...,K, AS RETURNED BY
*          DGEQLF IN THE LAST K COLUMNS OF ITS ARRAY ARGUMENT A.
*          A IS MODIFIED BY THE ROUTINE BUT RESTORED ON EXIT.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.
*          IF SIDE = 'L', LDA >= MAX(1,M);
*          IF SIDE = 'R', LDA >= MAX(1,N).
*
*  TAU     (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (K)
*          TAU(I) MUST CONTAIN THE SCALAR FACTOR OF THE ELEMENTARY
*          REFLECTOR H(I), AS RETURNED BY DGEQLF.
*
*  C       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDC,N)
*          ON ENTRY, THE M BY N MATRIX C.
*          ON EXIT, C IS OVERWRITTEN BY Q*C OR Q'*C OR C*Q' OR C*Q.
*
*  LDC     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY C. LDC >= MAX(1,M).
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION
*                                   (N) IF SIDE = 'L',
*                                   (M) IF SIDE = 'R'
*
*  INFO    (OUTPUT) INTEGER
*          = 0: SUCCESSFUL EXIT
*          < 0: IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, MI, NI, NQ
      DOUBLE PRECISION   AII
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLARF, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT ARGUMENTS
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ IS THE ORDER OF Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2L', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      IF( ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN ) )
     $     THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
      ELSE
         MI = M
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(I) IS APPLIED TO C(1:M-K+I,1:N)
*
            MI = M - K + I
         ELSE
*
*           H(I) IS APPLIED TO C(1:M,1:N-K+I)
*
            NI = N - K + I
         END IF
*
*        APPLY H(I)
*
         AII = A( NQ-K+I, I )
         A( NQ-K+I, I ) = ONE
         CALL DLARF( SIDE, MI, NI, A( 1, I ), 1, TAU( I ), C, LDC,
     $               WORK )
         A( NQ-K+I, I ) = AII
   10 CONTINUE
      RETURN
*
*     END OF DORM2L
*
      END
CUT HERE............
CAT > DSTEIN.F <<'CUT HERE............'
      SUBROUTINE DSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK,
     $                   IWORK, IFAIL, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, LDZ, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ),
     $                   IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DSTEIN COMPUTES THE EIGENVECTORS OF A REAL SYMMETRIC TRIDIAGONAL
*  MATRIX T CORRESPONDING TO SPECIFIED EIGENVALUES, USING INVERSE
*  ITERATION.
*
*  THE MAXIMUM NUMBER OF ITERATIONS ALLOWED FOR EACH EIGENVECTOR IS
*  SPECIFIED BY AN INTERNAL PARAMETER MAXITS (CURRENTLY SET TO 5).
*
*  ARGUMENTS
*  =========
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX.  N >= 0.
*
*  D       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          THE N DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX T.
*
*  E       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          THE (N-1) SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX
*          T, IN ELEMENTS 1 TO N-1.  E(N) NEED NOT BE SET.
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF EIGENVECTORS TO BE FOUND.  0 <= M <= N.
*
*  W       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          THE FIRST M ELEMENTS OF W CONTAIN THE EIGENVALUES FOR
*          WHICH EIGENVECTORS ARE TO BE COMPUTED.  THE EIGENVALUES
*          SHOULD BE GROUPED BY SPLIT-OFF BLOCK AND ORDERED FROM
*          SMALLEST TO LARGEST WITHIN THE BLOCK.  ( THE OUTPUT ARRAY
*          W FROM DSTEBZ WITH ORDER = 'B' IS EXPECTED HERE. )
*
*  IBLOCK  (INPUT) INTEGER ARRAY, DIMENSION (N)
*          THE SUBMATRIX INDICES ASSOCIATED WITH THE CORRESPONDING
*          EIGENVALUES IN W; IBLOCK(I)=1 IF EIGENVALUE W(I) BELONGS TO
*          THE FIRST SUBMATRIX FROM THE TOP, =2 IF W(I) BELONGS TO
*          THE SECOND SUBMATRIX, ETC.  ( THE OUTPUT ARRAY IBLOCK
*          FROM DSTEBZ IS EXPECTED HERE. )
*
*  ISPLIT  (INPUT) INTEGER ARRAY, DIMENSION (N)
*          THE SPLITTING POINTS, AT WHICH T BREAKS UP INTO SUBMATRICES.
*          THE FIRST SUBMATRIX CONSISTS OF ROWS/COLUMNS 1 TO
*          ISPLIT( 1 ), THE SECOND OF ROWS/COLUMNS ISPLIT( 1 )+1
*          THROUGH ISPLIT( 2 ), ETC.
*          ( THE OUTPUT ARRAY ISPLIT FROM DSTEBZ IS EXPECTED HERE. )
*
*  Z       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDZ, M)
*          THE COMPUTED EIGENVECTORS.  THE EIGENVECTOR ASSOCIATED
*          WITH THE EIGENVALUE W(I) IS STORED IN THE I-TH COLUMN OF
*          Z.  ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ITS CURRENT
*          ITERATE AFTER MAXITS ITERATIONS.
*
*  LDZ     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY Z.  LDZ >= MAX(1,N).
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (5*N)
*
*  IWORK   (WORKSPACE) INTEGER ARRAY, DIMENSION (N)
*
*  IFAIL   (OUTPUT) INTEGER ARRAY, DIMENSION (M)
*          ON NORMAL EXIT, ALL ELEMENTS OF IFAIL ARE ZERO.
*          IF ONE OR MORE EIGENVECTORS FAIL TO CONVERGE AFTER
*          MAXITS ITERATIONS, THEN THEIR INDICES ARE STORED IN
*          ARRAY IFAIL.
*
*  INFO    (OUTPUT) INTEGER
*          = 0: SUCCESSFUL EXIT.
*          < 0: IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0: IF INFO = I, THEN I EIGENVECTORS FAILED TO CONVERGE
*               IN MAXITS ITERATIONS.  THEIR INDICES ARE STORED IN
*               ARRAY IFAIL.
*
*  INTERNAL PARAMETERS
*  ===================
*
*  MAXITS  INTEGER, DEFAULT = 5
*          THE MAXIMUM NUMBER OF ITERATIONS PERFORMED.
*
*  EXTRA   INTEGER, DEFAULT = 2
*          THE NUMBER OF ITERATIONS PERFORMED AFTER NORM GROWTH
*          CRITERION IS SATISFIED, SHOULD BE AT LEAST 1.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO, ONE, TEN, ODM3, ODM1
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TEN = 1.0D1,
     $                   ODM3 = 1.0D-3, ODM1 = 1.0D-1 )
      INTEGER            MAXITS, EXTRA
      PARAMETER          ( MAXITS = 5, EXTRA = 2 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            B1, BLKSIZ, BN, GPIND, I, IINFO, INDRV1,
     $                   INDRV2, INDRV3, INDRV4, INDRV5, ITS, J, J1,
     $                   JBLK, JMAX, NBLK, NRMCHK
      DOUBLE PRECISION   DTPCRT, EPS, EPS1, NRM, ONENRM, ORTOL, PERTOL,
     $                   SCL, SEP, TOL, XJ, XJM, ZTR
*     ..
*     .. LOCAL ARRAYS ..
      INTEGER            ISEED( 4 )
*     ..
*     .. EXTERNAL FUNCTIONS ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM, DDOT, DLAMCH, DNRM2
      EXTERNAL           IDAMAX, DASUM, DDOT, DLAMCH, DNRM2
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DAXPY, DCOPY, DLAGTF, DLAGTS, DLARNV, DSCAL,
     $                   XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      DO 10 I = 1, M
         IFAIL( I ) = 0
   10 CONTINUE
*
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 .OR. M.GT.N ) THEN
         INFO = -4
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE
         DO 20 J = 2, M
            IF( IBLOCK( J ).LT.IBLOCK( J-1 ) ) THEN
               INFO = -6
               GO TO 30
            END IF
            IF( IBLOCK( J ).EQ.IBLOCK( J-1 ) .AND. W( J ).LT.W( J-1 ) )
     $           THEN
               INFO = -5
               GO TO 30
            END IF
   20    CONTINUE
   30    CONTINUE
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEIN', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.EQ.0 .OR. M.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     GET MACHINE CONSTANTS.
*
      EPS = DLAMCH( 'PRECISION' )
*
*     INITIALIZE SEED FOR RANDOM NUMBER GENERATOR DLARNV.
*
      DO 40 I = 1, 4
         ISEED( I ) = 1
   40 CONTINUE
*
*     INITIALIZE POINTERS.
*
      INDRV1 = 0
      INDRV2 = INDRV1 + N
      INDRV3 = INDRV2 + N
      INDRV4 = INDRV3 + N
      INDRV5 = INDRV4 + N
*
*     COMPUTE EIGENVECTORS OF MATRIX BLOCKS.
*
      GPIND = 1
      J1 = 1
      DO 160 NBLK = 1, IBLOCK( M )
*
*        FIND STARTING AND ENDING INDICES OF BLOCK NBLK.
*
         IF( NBLK.EQ.1 ) THEN
            B1 = 1
         ELSE
            B1 = ISPLIT( NBLK-1 ) + 1
         END IF
         BN = ISPLIT( NBLK )
         BLKSIZ = BN - B1 + 1
         IF( BLKSIZ.EQ.1 )
     $      GO TO 60
*
*        COMPUTE REORTHOGONALIZATION CRITERION AND STOPPING CRITERION.
*
         ONENRM = ABS( D( B1 ) ) + ABS( E( B1 ) )
         ONENRM = MAX( ONENRM, ABS( D( BN ) )+ABS( E( BN-1 ) ) )
         DO 50 I = B1 + 1, BN - 1
            ONENRM = MAX( ONENRM, ABS( D( I ) )+ABS( E( I-1 ) )+
     $               ABS( E( I ) ) )
   50    CONTINUE
         ORTOL = ODM3*ONENRM
*
         DTPCRT = SQRT( ODM1 / BLKSIZ )
*
*        LOOP THROUGH EIGENVALUES OF BLOCK NBLK.
*
   60    CONTINUE
         JBLK = 0
         DO 150 J = J1, M
            IF( IBLOCK( J ).NE.NBLK ) THEN
               J1 = J
               GO TO 160
            END IF
            JBLK = JBLK + 1
            XJ = W( J )
*
*           SKIP ALL THE WORK IF THE BLOCK SIZE IS ONE.
*
            IF( BLKSIZ.EQ.1 ) THEN
               WORK( INDRV1+1 ) = ONE
               GO TO 120
            END IF
*
*           IF EIGENVALUES J AND J-1 ARE TOO CLOSE, ADD A RELATIVELY
*           SMALL PERTURBATION.
*
            IF( JBLK.GT.1 ) THEN
               EPS1 = ABS( EPS*XJ )
               PERTOL = TEN*EPS1
               SEP = XJ - XJM
               IF( SEP.LT.PERTOL )
     $            XJ = XJM + PERTOL
            END IF
*
            ITS = 0
            NRMCHK = 0
*
*           GET RANDOM STARTING VECTOR.
*
            CALL DLARNV( 3, ISEED, BLKSIZ, WORK( INDRV1+1 ) )
*
*           COPY THE MATRIX T SO IT WON'T BE DESTROYED IN FACTORIZATION.
*
            CALL DCOPY( BLKSIZ, D( B1 ), 1, WORK( INDRV4+1 ), 1 )
            CALL DCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV2+2 ), 1 )
            CALL DCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV3+1 ), 1 )
*
*           COMPUTE LU FACTORS WITH PARTIAL PIVOTING  ( PT = LU )
*
            TOL = ZERO
            CALL DLAGTF( BLKSIZ, WORK( INDRV4+1 ), XJ, WORK( INDRV2+2 ),
     $                   WORK( INDRV3+1 ), TOL, WORK( INDRV5+1 ), IWORK,
     $                   IINFO )
*
*           UPDATE ITERATION COUNT.
*
   70       CONTINUE
            ITS = ITS + 1
            IF( ITS.GT.MAXITS )
     $         GO TO 100
*
*           NORMALIZE AND SCALE THE RIGHTHAND SIDE VECTOR PB.
*
            SCL = BLKSIZ*ONENRM*MAX( EPS,
     $            ABS( WORK( INDRV4+BLKSIZ ) ) ) /
     $            DASUM( BLKSIZ, WORK( INDRV1+1 ), 1 )
            CALL DSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
*
*           SOLVE THE SYSTEM LU = PB.
*
            CALL DLAGTS( -1, BLKSIZ, WORK( INDRV4+1 ), WORK( INDRV2+2 ),
     $                   WORK( INDRV3+1 ), WORK( INDRV5+1 ), IWORK,
     $                   WORK( INDRV1+1 ), TOL, IINFO )
*
*           REORTHOGONALIZE BY MODIFIED GRAM-SCHMIDT IF EIGENVALUES ARE
*           CLOSE ENOUGH.
*
            IF( JBLK.EQ.1 )
     $         GO TO 90
            IF( ABS( XJ-XJM ).GT.ORTOL )
     $         GPIND = J
            IF( GPIND.NE.J ) THEN
               DO 80 I = GPIND, J - 1
                  ZTR = -DDOT( BLKSIZ, WORK( INDRV1+1 ), 1, Z( B1, I ),
     $                  1 )
                  CALL DAXPY( BLKSIZ, ZTR, Z( B1, I ), 1,
     $                        WORK( INDRV1+1 ), 1 )
   80          CONTINUE
            END IF
*
*           CHECK THE INFINITY NORM OF THE ITERATE.
*
   90       CONTINUE
            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            NRM = ABS( WORK( INDRV1+JMAX ) )
*
*           CONTINUE FOR ADDITIONAL ITERATIONS AFTER NORM REACHES
*           STOPPING CRITERION.
*
            IF( NRM.LT.DTPCRT )
     $         GO TO 70
            NRMCHK = NRMCHK + 1
            IF( NRMCHK.LT.EXTRA+1 )
     $         GO TO 70
*
            GO TO 110
*
*           IF STOPPING CRITERION WAS NOT SATISFIED, UPDATE INFO AND
*           STORE EIGENVECTOR NUMBER IN ARRAY IFAIL.
*
  100       CONTINUE
            INFO = INFO + 1
            IFAIL( INFO ) = J
*
*           ACCEPT ITERATE AS JTH EIGENVECTOR.
*
  110       CONTINUE
            SCL = ONE / DNRM2( BLKSIZ, WORK( INDRV1+1 ), 1 )
            JMAX = IDAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
            IF( WORK( INDRV1+JMAX ).LT.ZERO )
     $         SCL = -SCL
            CALL DSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
  120       CONTINUE
            DO 130 I = 1, N
               Z( I, J ) = ZERO
  130       CONTINUE
            DO 140 I = 1, BLKSIZ
               Z( B1+I-1, J ) = WORK( INDRV1+I )
  140       CONTINUE
*
*           SAVE THE SHIFT TO CHECK EIGENVALUE SPACING AT NEXT
*           ITERATION.
*
            XJM = XJ
*
  150    CONTINUE
  160 CONTINUE
*
      RETURN
*
*     END OF DSTEIN
*
      END
CUT HERE............
CAT > DLAGTS.F <<'CUT HERE............'
      SUBROUTINE DLAGTS( JOB, N, A, B, C, D, IN, Y, TOL, INFO )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, JOB, N
      DOUBLE PRECISION   TOL
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IN( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * ), Y( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLAGTS MAY BE USED TO SOLVE ONE OF THE SYSTEMS OF EQUATIONS
*
*     (T - LAMBDA*I)*X = Y   OR   (T - LAMBDA*I)'*X = Y,
*
*  WHERE T IS AN N BY N TRIDIAGONAL MATRIX, FOR X, FOLLOWING THE
*  FACTORIZATION OF (T - LAMBDA*I) AS
*
*     (T - LAMBDA*I) = P*L*U ,
*
*  BY ROUTINE DLAGTF. THE CHOICE OF EQUATION TO BE SOLVED IS
*  CONTROLLED BY THE ARGUMENT JOB, AND IN EACH CASE THERE IS AN OPTION
*  TO PERTURB ZERO OR VERY SMALL DIAGONAL ELEMENTS OF U, THIS OPTION
*  BEING INTENDED FOR USE IN APPLICATIONS SUCH AS INVERSE ITERATION.
*
*  ARGUMENTS
*  =========
*
*  JOB     (INPUT) INTEGER
*          SPECIFIES THE JOB TO BE PERFORMED BY DLAGTS AS FOLLOWS:
*          =  1: THE EQUATIONS  (T - LAMBDA*I)X = Y  ARE TO BE SOLVED,
*                BUT DIAGONAL ELEMENTS OF U ARE NOT TO BE PERTURBED.
*          = -1: THE EQUATIONS  (T - LAMBDA*I)X = Y  ARE TO BE SOLVED
*                AND, IF OVERFLOW WOULD OTHERWISE OCCUR, THE DIAGONAL
*                ELEMENTS OF U ARE TO BE PERTURBED. SEE ARGUMENT TOL
*                BELOW.
*          =  2: THE EQUATIONS  (T - LAMBDA*I)'X = Y  ARE TO BE SOLVED,
*                BUT DIAGONAL ELEMENTS OF U ARE NOT TO BE PERTURBED.
*          = -2: THE EQUATIONS  (T - LAMBDA*I)'X = Y  ARE TO BE SOLVED
*                AND, IF OVERFLOW WOULD OTHERWISE OCCUR, THE DIAGONAL
*                ELEMENTS OF U ARE TO BE PERTURBED. SEE ARGUMENT TOL
*                BELOW.
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX T.
*
*  A       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          ON ENTRY, A MUST CONTAIN THE DIAGONAL ELEMENTS OF U AS
*          RETURNED FROM DLAGTF.
*
*  B       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          ON ENTRY, B MUST CONTAIN THE FIRST SUPER-DIAGONAL ELEMENTS OF
*          U AS RETURNED FROM DLAGTF.
*
*  C       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          ON ENTRY, C MUST CONTAIN THE SUB-DIAGONAL ELEMENTS OF L AS
*          RETURNED FROM DLAGTF.
*
*  D       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-2)
*          ON ENTRY, D MUST CONTAIN THE SECOND SUPER-DIAGONAL ELEMENTS
*          OF U AS RETURNED FROM DLAGTF.
*
*  IN      (INPUT) INTEGER ARRAY, DIMENSION (N)
*          ON ENTRY, IN MUST CONTAIN DETAILS OF THE MATRIX P AS RETURNED
*          FROM DLAGTF.
*
*  Y       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          ON ENTRY, THE RIGHT HAND SIDE VECTOR Y.
*
*          ON EXIT, Y IS OVERWRITTEN BY THE SOLUTION VECTOR X.
*
*  TOL     (INPUT/OUTPUT) DOUBLE PRECISION
*          ON ENTRY WITH  JOB .LT. 0, TOL SHOULD BE THE MINIMUM
*          PERTURBATION TO BE MADE TO VERY SMALL DIAGONAL ELEMENTS OF U.
*          TOL SHOULD NORMALLY BE CHOSEN AS ABOUT EPS*NORM(U), WHERE EPS
*          IS THE RELATIVE MACHINE PRECISION, BUT IF TOL IS SUPPLIED AS
*          NON-POSITIVE, THEN IT IS RESET TO EPS*MAX( ABS( U(I,J) ) ).
*          IF  JOB .GT. 0  THEN TOL IS NOT REFERENCED.
*
*          ON EXIT, TOL IS CHANGED AS DESCRIBED ABOVE, ONLY IF TOL IS
*          NON-POSITIVE ON ENTRY. OTHERWISE TOL IS UNCHANGED.
*
*  INFO    (OUTPUT)
*          = 0   : SUCCESSFUL EXIT
*          .LT. 0: IF INFO = -K, THE KTH ARGUMENT HAD AN ILLEGAL VALUE
*          .GT. 0: OVERFLOW WOULD OCCUR WHEN COMPUTING THE INFO(TH)
*                  ELEMENT OF THE SOLUTION VECTOR X. THIS CAN ONLY OCCUR
*                  WHEN JOB IS SUPPLIED AS POSITIVE AND EITHER MEANS
*                  THAT A DIAGONAL ELEMENT OF U IS VERY SMALL, OR THAT
*                  THE ELEMENTS OF THE RIGHT-HAND SIDE VECTOR Y ARE VERY
*                  LARGE.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            K
      DOUBLE PRECISION   ABSAK, AK, BIGNUM, EPS, PERT, SFMIN, TEMP
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, MAX, SIGN
*     ..
*     .. EXTERNAL FUNCTIONS ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           XERBLA
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      INFO = 0
      IF( ( ABS( JOB ).GT.2 ) .OR. ( JOB.EQ.0 ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAGTS', -INFO )
         RETURN
      END IF
*
      IF( N.EQ.0 )
     $   RETURN
*
      EPS = DLAMCH( 'EPSILON' )
      SFMIN = DLAMCH( 'SAFE MINIMUM' )
      BIGNUM = ONE / SFMIN
*
      IF( JOB.LT.0 ) THEN
         IF( TOL.LE.ZERO ) THEN
            TOL = ABS( A( 1 ) )
            IF( N.GT.1 )
     $         TOL = MAX( TOL, ABS( A( 2 ) ), ABS( B( 1 ) ) )
            DO 10 K = 3, N
               TOL = MAX( TOL, ABS( A( K ) ), ABS( B( K-1 ) ),
     $               ABS( D( K-2 ) ) )
   10       CONTINUE
            TOL = TOL*EPS
            IF( TOL.EQ.ZERO )
     $         TOL = EPS
         END IF
      END IF
*
      IF( ABS( JOB ).EQ.1 ) THEN
         DO 20 K = 2, N
            IF( IN( K-1 ).EQ.0 ) THEN
               Y( K ) = Y( K ) - C( K-1 )*Y( K-1 )
            ELSE
               TEMP = Y( K-1 )
               Y( K-1 ) = Y( K )
               Y( K ) = TEMP - C( K-1 )*Y( K )
            END IF
   20    CONTINUE
         IF( JOB.EQ.1 ) THEN
            DO 30 K = N, 1, -1
               IF( K.LE.N-2 ) THEN
                  TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 )
               ELSE IF( K.EQ.N-1 ) THEN
                  TEMP = Y( K ) - B( K )*Y( K+1 )
               ELSE
                  TEMP = Y( K )
               END IF
               AK = A( K )
               ABSAK = ABS( AK )
               IF( ABSAK.LT.ONE ) THEN
                  IF( ABSAK.LT.SFMIN ) THEN
                     IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK )
     $                    THEN
                        INFO = K
                        RETURN
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
                     INFO = K
                     RETURN
                  END IF
               END IF
               Y( K ) = TEMP / AK
   30       CONTINUE
         ELSE
            DO 50 K = N, 1, -1
               IF( K.LE.N-2 ) THEN
                  TEMP = Y( K ) - B( K )*Y( K+1 ) - D( K )*Y( K+2 )
               ELSE IF( K.EQ.N-1 ) THEN
                  TEMP = Y( K ) - B( K )*Y( K+1 )
               ELSE
                  TEMP = Y( K )
               END IF
               AK = A( K )
               PERT = SIGN( TOL, AK )
   40          CONTINUE
               ABSAK = ABS( AK )
               IF( ABSAK.LT.ONE ) THEN
                  IF( ABSAK.LT.SFMIN ) THEN
                     IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK )
     $                    THEN
                        AK = AK + PERT
                        PERT = 2*PERT
                        GO TO 40
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
                     AK = AK + PERT
                     PERT = 2*PERT
                     GO TO 40
                  END IF
               END IF
               Y( K ) = TEMP / AK
   50       CONTINUE
         END IF
      ELSE
*
*        COME TO HERE IF  JOB = 2 OR -2
*
         IF( JOB.EQ.2 ) THEN
            DO 60 K = 1, N
               IF( K.GE.3 ) THEN
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 )
               ELSE IF( K.EQ.2 ) THEN
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 )
               ELSE
                  TEMP = Y( K )
               END IF
               AK = A( K )
               ABSAK = ABS( AK )
               IF( ABSAK.LT.ONE ) THEN
                  IF( ABSAK.LT.SFMIN ) THEN
                     IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK )
     $                    THEN
                        INFO = K
                        RETURN
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
                     INFO = K
                     RETURN
                  END IF
               END IF
               Y( K ) = TEMP / AK
   60       CONTINUE
         ELSE
            DO 80 K = 1, N
               IF( K.GE.3 ) THEN
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 ) - D( K-2 )*Y( K-2 )
               ELSE IF( K.EQ.2 ) THEN
                  TEMP = Y( K ) - B( K-1 )*Y( K-1 )
               ELSE
                  TEMP = Y( K )
               END IF
               AK = A( K )
               PERT = SIGN( TOL, AK )
   70          CONTINUE
               ABSAK = ABS( AK )
               IF( ABSAK.LT.ONE ) THEN
                  IF( ABSAK.LT.SFMIN ) THEN
                     IF( ABSAK.EQ.ZERO .OR. ABS( TEMP )*SFMIN.GT.ABSAK )
     $                    THEN
                        AK = AK + PERT
                        PERT = 2*PERT
                        GO TO 70
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF( ABS( TEMP ).GT.ABSAK*BIGNUM ) THEN
                     AK = AK + PERT
                     PERT = 2*PERT
                     GO TO 70
                  END IF
               END IF
               Y( K ) = TEMP / AK
   80       CONTINUE
         END IF
*
         DO 90 K = N, 2, -1
            IF( IN( K-1 ).EQ.0 ) THEN
               Y( K-1 ) = Y( K-1 ) - C( K-1 )*Y( K )
            ELSE
               TEMP = Y( K-1 )
               Y( K-1 ) = Y( K )
               Y( K ) = TEMP - C( K-1 )*Y( K )
            END IF
   90    CONTINUE
      END IF
*
*     END OF DLAGTS
*
      END
CUT HERE............
CAT > DLAGTF.F <<'CUT HERE............'
      SUBROUTINE DLAGTF( N, A, LAMBDA, B, C, TOL, D, IN, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, N
      DOUBLE PRECISION   LAMBDA, TOL
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IN( * )
      DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLAGTF FACTORIZES THE MATRIX (T - LAMBDA*I), WHERE T IS AN N BY N
*  TRIDIAGONAL MATRIX AND LAMBDA IS A SCALAR, AS
*
*     T - LAMBDA*I = PLU,
*
*  WHERE P IS A PERMUTATION MATRIX, L IS A UNIT LOWER TRIDIAGONAL MATRIX
*  WITH AT MOST ONE NON-ZERO SUB-DIAGONAL ELEMENTS PER COLUMN AND U IS
*  AN UPPER TRIANGULAR MATRIX WITH AT MOST TWO NON-ZERO SUPER-DIAGONAL
*  ELEMENTS PER COLUMN.
*
*  THE FACTORIZATION IS OBTAINED BY GAUSSIAN ELIMINATION WITH PARTIAL
*  PIVOTING AND IMPLICIT ROW SCALING.
*
*  THE PARAMETER LAMBDA IS INCLUDED IN THE ROUTINE SO THAT DLAGTF MAY
*  BE USED, IN CONJUNCTION WITH DLAGTS, TO OBTAIN EIGENVECTORS OF T BY
*  INVERSE ITERATION.
*
*  ARGUMENTS
*  =========
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX T.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          ON ENTRY, A MUST CONTAIN THE DIAGONAL ELEMENTS OF T.
*
*          ON EXIT, A IS OVERWRITTEN BY THE N DIAGONAL ELEMENTS OF THE
*          UPPER TRIANGULAR MATRIX U OF THE FACTORIZATION OF T.
*
*  LAMBDA  (INPUT) DOUBLE PRECISION
*          ON ENTRY, THE SCALAR LAMBDA.
*
*  B       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          ON ENTRY, B MUST CONTAIN THE (N-1) SUPER-DIAGONAL ELEMENTS OF
*          T.
*
*          ON EXIT, B IS OVERWRITTEN BY THE (N-1) SUPER-DIAGONAL
*          ELEMENTS OF THE MATRIX U OF THE FACTORIZATION OF T.
*
*  C       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          ON ENTRY, C MUST CONTAIN THE (N-1) SUB-DIAGONAL ELEMENTS OF
*          T.
*
*          ON EXIT, C IS OVERWRITTEN BY THE (N-1) SUB-DIAGONAL ELEMENTS
*          OF THE MATRIX L OF THE FACTORIZATION OF T.
*
*  TOL     (INPUT) DOUBLE PRECISION
*          ON ENTRY, A RELATIVE TOLERANCE USED TO INDICATE WHETHER OR
*          NOT THE MATRIX (T - LAMBDA*I) IS NEARLY SINGULAR. TOL SHOULD
*          NORMALLY BE CHOSE AS APPROXIMATELY THE LARGEST RELATIVE ERROR
*          IN THE ELEMENTS OF T. FOR EXAMPLE, IF THE ELEMENTS OF T ARE
*          CORRECT TO ABOUT 4 SIGNIFICANT FIGURES, THEN TOL SHOULD BE
*          SET TO ABOUT 5*10**(-4). IF TOL IS SUPPLIED AS LESS THAN EPS,
*          WHERE EPS IS THE RELATIVE MACHINE PRECISION, THEN THE VALUE
*          EPS IS USED IN PLACE OF TOL.
*
*  D       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-2)
*          ON EXIT, D IS OVERWRITTEN BY THE (N-2) SECOND SUPER-DIAGONAL
*          ELEMENTS OF THE MATRIX U OF THE FACTORIZATION OF T.
*
*  IN      (OUTPUT) INTEGER ARRAY, DIMENSION (N)
*          ON EXIT, IN CONTAINS DETAILS OF THE PERMUTATION MATRIX P. IF
*          AN INTERCHANGE OCCURRED AT THE KTH STEP OF THE ELIMINATION,
*          THEN IN(K) = 1, OTHERWISE IN(K) = 0. THE ELEMENT IN(N)
*          RETURNS THE SMALLEST POSITIVE INTEGER J SUCH THAT
*
*             ABS( U(J,J) ).LE. NORM( (T - LAMBDA*I)(J) )*TOL,
*
*          WHERE NORM( A(J) ) DENOTES THE SUM OF THE ABSOLUTE VALUES OF
*          THE JTH ROW OF THE MATRIX A. IF NO SUCH J EXISTS THEN IN(N)
*          IS RETURNED AS ZERO. IF IN(N) IS RETURNED AS POSITIVE, THEN A
*          DIAGONAL ELEMENT OF U IS SMALL, INDICATING THAT
*          (T - LAMBDA*I) IS SINGULAR OR NEARLY SINGULAR,
*
*  INFO    (OUTPUT)
*          = 0   : SUCCESSFUL EXIT
*          .LT. 0: IF INFO = -K, THE KTH ARGUMENT HAD AN ILLEGAL VALUE
*
* =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            K
      DOUBLE PRECISION   EPS, MULT, PIV1, PIV2, SCALE1, SCALE2, TEMP, TL
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, MAX
*     ..
*     .. EXTERNAL FUNCTIONS ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           XERBLA
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DLAGTF', -INFO )
         RETURN
      END IF
*
      IF( N.EQ.0 )
     $   RETURN
*
      A( 1 ) = A( 1 ) - LAMBDA
      IN( N ) = 0
      IF( N.EQ.1 ) THEN
         IF( A( 1 ).EQ.ZERO )
     $      IN( 1 ) = 1
         RETURN
      END IF
*
      EPS = DLAMCH( 'EPSILON' )
*
      TL = MAX( TOL, EPS )
      SCALE1 = ABS( A( 1 ) ) + ABS( B( 1 ) )
      DO 10 K = 1, N - 1
         A( K+1 ) = A( K+1 ) - LAMBDA
         SCALE2 = ABS( C( K ) ) + ABS( A( K+1 ) )
         IF( K.LT.( N-1 ) )
     $      SCALE2 = SCALE2 + ABS( B( K+1 ) )
         IF( A( K ).EQ.ZERO ) THEN
            PIV1 = ZERO
         ELSE
            PIV1 = ABS( A( K ) ) / SCALE1
         END IF
         IF( C( K ).EQ.ZERO ) THEN
            IN( K ) = 0
            PIV2 = ZERO
            SCALE1 = SCALE2
            IF( K.LT.( N-1 ) )
     $         D( K ) = ZERO
         ELSE
            PIV2 = ABS( C( K ) ) / SCALE2
            IF( PIV2.LE.PIV1 ) THEN
               IN( K ) = 0
               SCALE1 = SCALE2
               C( K ) = C( K ) / A( K )
               A( K+1 ) = A( K+1 ) - C( K )*B( K )
               IF( K.LT.( N-1 ) )
     $            D( K ) = ZERO
            ELSE
               IN( K ) = 1
               MULT = A( K ) / C( K )
               A( K ) = C( K )
               TEMP = A( K+1 )
               A( K+1 ) = B( K ) - MULT*TEMP
               IF( K.LT.( N-1 ) ) THEN
                  D( K ) = B( K+1 )
                  B( K+1 ) = -MULT*D( K )
               END IF
               B( K ) = TEMP
               C( K ) = MULT
            END IF
         END IF
         IF( ( MAX( PIV1, PIV2 ).LE.TL ) .AND. ( IN( N ).EQ.0 ) )
     $      IN( N ) = K
   10 CONTINUE
      IF( ( ABS( A( N ) ).LE.SCALE1*TL ) .AND. ( IN( N ).EQ.0 ) )
     $   IN( N ) = N
*
      RETURN
*
*     END OF DLAGTF
*
      END
CUT HERE............
CAT > DLARNV.F <<'CUT HERE............'
      SUBROUTINE DLARNV( IDIST, ISEED, N, X )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            IDIST, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   X( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLARNV RETURNS A VECTOR OF N RANDOM REAL NUMBERS FROM A UNIFORM OR
*  NORMAL DISTRIBUTION.
*
*  ARGUMENTS
*  =========
*
*  IDIST   (INPUT) INTEGER
*          SPECIFIES THE DISTRIBUTION OF THE RANDOM NUMBERS:
*          = 1:  UNIFORM (0,1)
*          = 2:  UNIFORM (-1,1)
*          = 3:  NORMAL (0,1)
*
*  ISEED   (INPUT/OUTPUT) INTEGER ARRAY, DIMENSION (4)
*          ON ENTRY, THE SEED OF THE RANDOM NUMBER GENERATOR; THE ARRAY
*          ELEMENTS MUST BE BETWEEN 0 AND 4095, AND ISEED(4) MUST BE
*          ODD.
*          ON EXIT, THE SEED IS UPDATED.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF RANDOM NUMBERS TO BE GENERATED.
*
*  X       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          THE GENERATED RANDOM NUMBERS.
*
*  FURTHER DETAILS
*  ===============
*
*  THIS ROUTINE CALLS THE AUXILIARY ROUTINE DLARUV TO GENERATE RANDOM
*  REAL NUMBERS FROM A UNIFORM (0,1) DISTRIBUTION, IN BATCHES OF UP TO
*  128 USING VECTORISABLE CODE. THE BOX-MULLER METHOD IS USED TO
*  TRANSFORM NUMBERS FROM A UNIFORM TO A NORMAL DISTRIBUTION.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, TWO
      PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 )
      INTEGER            LV
      PARAMETER          ( LV = 128 )
      DOUBLE PRECISION   TWOPI
      PARAMETER          ( TWOPI = 6.28318530717958623199592D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, IL, IL2, IV
*     ..
*     .. LOCAL ARRAYS ..
      DOUBLE PRECISION   U( LV )
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          COS, LOG, MIN, SQRT
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLARUV
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      DO 40 IV = 1, N, LV / 2
         IL = MIN( LV / 2, N-IV+1 )
         IF( IDIST.EQ.3 ) THEN
            IL2 = 2*IL
         ELSE
            IL2 = IL
         END IF
*
*        CALL DLARUV TO GENERATE IL2 NUMBERS FROM A UNIFORM (0,1)
*        DISTRIBUTION (IL2 <= LV)
*
         CALL DLARUV( ISEED, IL2, U )
*
         IF( IDIST.EQ.1 ) THEN
*
*           COPY GENERATED NUMBERS
*
            DO 10 I = 1, IL
               X( IV+I-1 ) = U( I )
   10       CONTINUE
         ELSE IF( IDIST.EQ.2 ) THEN
*
*           CONVERT GENERATED NUMBERS TO UNIFORM (-1,1) DISTRIBUTION
*
            DO 20 I = 1, IL
               X( IV+I-1 ) = TWO*U( I ) - ONE
   20       CONTINUE
         ELSE IF( IDIST.EQ.3 ) THEN
*
*           CONVERT GENERATED NUMBERS TO NORMAL (0,1) DISTRIBUTION
*
            DO 30 I = 1, IL
               X( IV+I-1 ) = SQRT( -TWO*LOG( U( 2*I-1 ) ) )*
     $                       COS( TWOPI*U( 2*I ) )
   30       CONTINUE
         END IF
   40 CONTINUE
      RETURN
*
*     END OF DLARNV
*
      END
CUT HERE............
CAT > DLARUV.F <<'CUT HERE............'
      SUBROUTINE DLARUV( ISEED, N, X )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            ISEED( 4 )
      DOUBLE PRECISION   X( N )
*     ..
*
*  PURPOSE
*  =======
*
*  DLARUV RETURNS A VECTOR OF N RANDOM REAL NUMBERS FROM A UNIFORM (0,1)
*  DISTRIBUTION (N <= 128).
*
*  THIS IS AN AUXILIARY ROUTINE CALLED BY DLARNV AND ZLARNV.
*
*  ARGUMENTS
*  =========
*
*  ISEED   (INPUT/OUTPUT) INTEGER ARRAY, DIMENSION (4)
*          ON ENTRY, THE SEED OF THE RANDOM NUMBER GENERATOR; THE ARRAY
*          ELEMENTS MUST BE BETWEEN 0 AND 4095, AND ISEED(4) MUST BE
*          ODD.
*          ON EXIT, THE SEED IS UPDATED.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF RANDOM NUMBERS TO BE GENERATED. N <= 128.
*
*  X       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          THE GENERATED RANDOM NUMBERS.
*
*  FURTHER DETAILS
*  ===============
*
*  THIS ROUTINE USES A MULTIPLICATIVE CONGRUENTIAL METHOD WITH MODULUS
*  2**48 AND MULTIPLIER 33952834046453 (SEE G.S.FISHMAN,
*  'MULTIPLICATIVE CONGRUENTIAL RANDOM NUMBER GENERATORS WITH MODULUS
*  2**B: AN EXHAUSTIVE ANALYSIS FOR B = 32 AND A PARTIAL ANALYSIS FOR
*  B = 48', MATH. COMP. 189, PP 331-344, 1990).
*
*  48-BIT INTEGERS ARE STORED IN 4 INTEGER ARRAY ELEMENTS WITH 12 BITS
*  PER ELEMENT. HENCE THE ROUTINE IS PORTABLE ACROSS MACHINES WITH
*  INTEGERS OF 32 BITS OR MORE.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      INTEGER            LV, IPW2
      DOUBLE PRECISION   R
      PARAMETER          ( LV = 128, IPW2 = 4096, R = ONE / IPW2 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, I1, I2, I3, I4, IT1, IT2, IT3, IT4, J
*     ..
*     .. LOCAL ARRAYS ..
      INTEGER            MM( LV, 4 )
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          DBLE, MIN, MOD
*     ..
*     .. DATA STATEMENTS ..
      DATA               ( MM( 1, J ), J = 1, 4 ) / 494, 322, 2508,
     $                   2549 /
      DATA               ( MM( 2, J ), J = 1, 4 ) / 2637, 789, 3754,
     $                   1145 /
      DATA               ( MM( 3, J ), J = 1, 4 ) / 255, 1440, 1766,
     $                   2253 /
      DATA               ( MM( 4, J ), J = 1, 4 ) / 2008, 752, 3572,
     $                   305 /
      DATA               ( MM( 5, J ), J = 1, 4 ) / 1253, 2859, 2893,
     $                   3301 /
      DATA               ( MM( 6, J ), J = 1, 4 ) / 3344, 123, 307,
     $                   1065 /
      DATA               ( MM( 7, J ), J = 1, 4 ) / 4084, 1848, 1297,
     $                   3133 /
      DATA               ( MM( 8, J ), J = 1, 4 ) / 1739, 643, 3966,
     $                   2913 /
      DATA               ( MM( 9, J ), J = 1, 4 ) / 3143, 2405, 758,
     $                   3285 /
      DATA               ( MM( 10, J ), J = 1, 4 ) / 3468, 2638, 2598,
     $                   1241 /
      DATA               ( MM( 11, J ), J = 1, 4 ) / 688, 2344, 3406,
     $                   1197 /
      DATA               ( MM( 12, J ), J = 1, 4 ) / 1657, 46, 2922,
     $                   3729 /
      DATA               ( MM( 13, J ), J = 1, 4 ) / 1238, 3814, 1038,
     $                   2501 /
      DATA               ( MM( 14, J ), J = 1, 4 ) / 3166, 913, 2934,
     $                   1673 /
      DATA               ( MM( 15, J ), J = 1, 4 ) / 1292, 3649, 2091,
     $                   541 /
      DATA               ( MM( 16, J ), J = 1, 4 ) / 3422, 339, 2451,
     $                   2753 /
      DATA               ( MM( 17, J ), J = 1, 4 ) / 1270, 3808, 1580,
     $                   949 /
      DATA               ( MM( 18, J ), J = 1, 4 ) / 2016, 822, 1958,
     $                   2361 /
      DATA               ( MM( 19, J ), J = 1, 4 ) / 154, 2832, 2055,
     $                   1165 /
      DATA               ( MM( 20, J ), J = 1, 4 ) / 2862, 3078, 1507,
     $                   4081 /
      DATA               ( MM( 21, J ), J = 1, 4 ) / 697, 3633, 1078,
     $                   2725 /
      DATA               ( MM( 22, J ), J = 1, 4 ) / 1706, 2970, 3273,
     $                   3305 /
      DATA               ( MM( 23, J ), J = 1, 4 ) / 491, 637, 17,
     $                   3069 /
      DATA               ( MM( 24, J ), J = 1, 4 ) / 931, 2249, 854,
     $                   3617 /
      DATA               ( MM( 25, J ), J = 1, 4 ) / 1444, 2081, 2916,
     $                   3733 /
      DATA               ( MM( 26, J ), J = 1, 4 ) / 444, 4019, 3971,
     $                   409 /
      DATA               ( MM( 27, J ), J = 1, 4 ) / 3577, 1478, 2889,
     $                   2157 /
      DATA               ( MM( 28, J ), J = 1, 4 ) / 3944, 242, 3831,
     $                   1361 /
      DATA               ( MM( 29, J ), J = 1, 4 ) / 2184, 481, 2621,
     $                   3973 /
      DATA               ( MM( 30, J ), J = 1, 4 ) / 1661, 2075, 1541,
     $                   1865 /
      DATA               ( MM( 31, J ), J = 1, 4 ) / 3482, 4058, 893,
     $                   2525 /
      DATA               ( MM( 32, J ), J = 1, 4 ) / 657, 622, 736,
     $                   1409 /
      DATA               ( MM( 33, J ), J = 1, 4 ) / 3023, 3376, 3992,
     $                   3445 /
      DATA               ( MM( 34, J ), J = 1, 4 ) / 3618, 812, 787,
     $                   3577 /
      DATA               ( MM( 35, J ), J = 1, 4 ) / 1267, 234, 2125,
     $                   77 /
      DATA               ( MM( 36, J ), J = 1, 4 ) / 1828, 641, 2364,
     $                   3761 /
      DATA               ( MM( 37, J ), J = 1, 4 ) / 164, 4005, 2460,
     $                   2149 /
      DATA               ( MM( 38, J ), J = 1, 4 ) / 3798, 1122, 257,
     $                   1449 /
      DATA               ( MM( 39, J ), J = 1, 4 ) / 3087, 3135, 1574,
     $                   3005 /
      DATA               ( MM( 40, J ), J = 1, 4 ) / 2400, 2640, 3912,
     $                   225 /
      DATA               ( MM( 41, J ), J = 1, 4 ) / 2870, 2302, 1216,
     $                   85 /
      DATA               ( MM( 42, J ), J = 1, 4 ) / 3876, 40, 3248,
     $                   3673 /
      DATA               ( MM( 43, J ), J = 1, 4 ) / 1905, 1832, 3401,
     $                   3117 /
      DATA               ( MM( 44, J ), J = 1, 4 ) / 1593, 2247, 2124,
     $                   3089 /
      DATA               ( MM( 45, J ), J = 1, 4 ) / 1797, 2034, 2762,
     $                   1349 /
      DATA               ( MM( 46, J ), J = 1, 4 ) / 1234, 2637, 149,
     $                   2057 /
      DATA               ( MM( 47, J ), J = 1, 4 ) / 3460, 1287, 2245,
     $                   413 /
      DATA               ( MM( 48, J ), J = 1, 4 ) / 328, 1691, 166,
     $                   65 /
      DATA               ( MM( 49, J ), J = 1, 4 ) / 2861, 496, 466,
     $                   1845 /
      DATA               ( MM( 50, J ), J = 1, 4 ) / 1950, 1597, 4018,
     $                   697 /
      DATA               ( MM( 51, J ), J = 1, 4 ) / 617, 2394, 1399,
     $                   3085 /
      DATA               ( MM( 52, J ), J = 1, 4 ) / 2070, 2584, 190,
     $                   3441 /
      DATA               ( MM( 53, J ), J = 1, 4 ) / 3331, 1843, 2879,
     $                   1573 /
      DATA               ( MM( 54, J ), J = 1, 4 ) / 769, 336, 153,
     $                   3689 /
      DATA               ( MM( 55, J ), J = 1, 4 ) / 1558, 1472, 2320,
     $                   2941 /
      DATA               ( MM( 56, J ), J = 1, 4 ) / 2412, 2407, 18,
     $                   929 /
      DATA               ( MM( 57, J ), J = 1, 4 ) / 2800, 433, 712,
     $                   533 /
      DATA               ( MM( 58, J ), J = 1, 4 ) / 189, 2096, 2159,
     $                   2841 /
      DATA               ( MM( 59, J ), J = 1, 4 ) / 287, 1761, 2318,
     $                   4077 /
      DATA               ( MM( 60, J ), J = 1, 4 ) / 2045, 2810, 2091,
     $                   721 /
      DATA               ( MM( 61, J ), J = 1, 4 ) / 1227, 566, 3443,
     $                   2821 /
      DATA               ( MM( 62, J ), J = 1, 4 ) / 2838, 442, 1510,
     $                   2249 /
      DATA               ( MM( 63, J ), J = 1, 4 ) / 209, 41, 449,
     $                   2397 /
      DATA               ( MM( 64, J ), J = 1, 4 ) / 2770, 1238, 1956,
     $                   2817 /
      DATA               ( MM( 65, J ), J = 1, 4 ) / 3654, 1086, 2201,
     $                   245 /
      DATA               ( MM( 66, J ), J = 1, 4 ) / 3993, 603, 3137,
     $                   1913 /
      DATA               ( MM( 67, J ), J = 1, 4 ) / 192, 840, 3399,
     $                   1997 /
      DATA               ( MM( 68, J ), J = 1, 4 ) / 2253, 3168, 1321,
     $                   3121 /
      DATA               ( MM( 69, J ), J = 1, 4 ) / 3491, 1499, 2271,
     $                   997 /
      DATA               ( MM( 70, J ), J = 1, 4 ) / 2889, 1084, 3667,
     $                   1833 /
      DATA               ( MM( 71, J ), J = 1, 4 ) / 2857, 3438, 2703,
     $                   2877 /
      DATA               ( MM( 72, J ), J = 1, 4 ) / 2094, 2408, 629,
     $                   1633 /
      DATA               ( MM( 73, J ), J = 1, 4 ) / 1818, 1589, 2365,
     $                   981 /
      DATA               ( MM( 74, J ), J = 1, 4 ) / 688, 2391, 2431,
     $                   2009 /
      DATA               ( MM( 75, J ), J = 1, 4 ) / 1407, 288, 1113,
     $                   941 /
      DATA               ( MM( 76, J ), J = 1, 4 ) / 634, 26, 3922,
     $                   2449 /
      DATA               ( MM( 77, J ), J = 1, 4 ) / 3231, 512, 2554,
     $                   197 /
      DATA               ( MM( 78, J ), J = 1, 4 ) / 815, 1456, 184,
     $                   2441 /
      DATA               ( MM( 79, J ), J = 1, 4 ) / 3524, 171, 2099,
     $                   285 /
      DATA               ( MM( 80, J ), J = 1, 4 ) / 1914, 1677, 3228,
     $                   1473 /
      DATA               ( MM( 81, J ), J = 1, 4 ) / 516, 2657, 4012,
     $                   2741 /
      DATA               ( MM( 82, J ), J = 1, 4 ) / 164, 2270, 1921,
     $                   3129 /
      DATA               ( MM( 83, J ), J = 1, 4 ) / 303, 2587, 3452,
     $                   909 /
      DATA               ( MM( 84, J ), J = 1, 4 ) / 2144, 2961, 3901,
     $                   2801 /
      DATA               ( MM( 85, J ), J = 1, 4 ) / 3480, 1970, 572,
     $                   421 /
      DATA               ( MM( 86, J ), J = 1, 4 ) / 119, 1817, 3309,
     $                   4073 /
      DATA               ( MM( 87, J ), J = 1, 4 ) / 3357, 676, 3171,
     $                   2813 /
      DATA               ( MM( 88, J ), J = 1, 4 ) / 837, 1410, 817,
     $                   2337 /
      DATA               ( MM( 89, J ), J = 1, 4 ) / 2826, 3723, 3039,
     $                   1429 /
      DATA               ( MM( 90, J ), J = 1, 4 ) / 2332, 2803, 1696,
     $                   1177 /
      DATA               ( MM( 91, J ), J = 1, 4 ) / 2089, 3185, 1256,
     $                   1901 /
      DATA               ( MM( 92, J ), J = 1, 4 ) / 3780, 184, 3715,
     $                   81 /
      DATA               ( MM( 93, J ), J = 1, 4 ) / 1700, 663, 2077,
     $                   1669 /
      DATA               ( MM( 94, J ), J = 1, 4 ) / 3712, 499, 3019,
     $                   2633 /
      DATA               ( MM( 95, J ), J = 1, 4 ) / 150, 3784, 1497,
     $                   2269 /
      DATA               ( MM( 96, J ), J = 1, 4 ) / 2000, 1631, 1101,
     $                   129 /
      DATA               ( MM( 97, J ), J = 1, 4 ) / 3375, 1925, 717,
     $                   1141 /
      DATA               ( MM( 98, J ), J = 1, 4 ) / 1621, 3912, 51,
     $                   249 /
      DATA               ( MM( 99, J ), J = 1, 4 ) / 3090, 1398, 981,
     $                   3917 /
      DATA               ( MM( 100, J ), J = 1, 4 ) / 3765, 1349, 1978,
     $                   2481 /
      DATA               ( MM( 101, J ), J = 1, 4 ) / 1149, 1441, 1813,
     $                   3941 /
      DATA               ( MM( 102, J ), J = 1, 4 ) / 3146, 2224, 3881,
     $                   2217 /
      DATA               ( MM( 103, J ), J = 1, 4 ) / 33, 2411, 76,
     $                   2749 /
      DATA               ( MM( 104, J ), J = 1, 4 ) / 3082, 1907, 3846,
     $                   3041 /
      DATA               ( MM( 105, J ), J = 1, 4 ) / 2741, 3192, 3694,
     $                   1877 /
      DATA               ( MM( 106, J ), J = 1, 4 ) / 359, 2786, 1682,
     $                   345 /
      DATA               ( MM( 107, J ), J = 1, 4 ) / 3316, 382, 124,
     $                   2861 /
      DATA               ( MM( 108, J ), J = 1, 4 ) / 1749, 37, 1660,
     $                   1809 /
      DATA               ( MM( 109, J ), J = 1, 4 ) / 185, 759, 3997,
     $                   3141 /
      DATA               ( MM( 110, J ), J = 1, 4 ) / 2784, 2948, 479,
     $                   2825 /
      DATA               ( MM( 111, J ), J = 1, 4 ) / 2202, 1862, 1141,
     $                   157 /
      DATA               ( MM( 112, J ), J = 1, 4 ) / 2199, 3802, 886,
     $                   2881 /
      DATA               ( MM( 113, J ), J = 1, 4 ) / 1364, 2423, 3514,
     $                   3637 /
      DATA               ( MM( 114, J ), J = 1, 4 ) / 1244, 2051, 1301,
     $                   1465 /
      DATA               ( MM( 115, J ), J = 1, 4 ) / 2020, 2295, 3604,
     $                   2829 /
      DATA               ( MM( 116, J ), J = 1, 4 ) / 3160, 1332, 1888,
     $                   2161 /
      DATA               ( MM( 117, J ), J = 1, 4 ) / 2785, 1832, 1836,
     $                   3365 /
      DATA               ( MM( 118, J ), J = 1, 4 ) / 2772, 2405, 1990,
     $                   361 /
      DATA               ( MM( 119, J ), J = 1, 4 ) / 1217, 3638, 2058,
     $                   2685 /
      DATA               ( MM( 120, J ), J = 1, 4 ) / 1822, 3661, 692,
     $                   3745 /
      DATA               ( MM( 121, J ), J = 1, 4 ) / 1245, 327, 1194,
     $                   2325 /
      DATA               ( MM( 122, J ), J = 1, 4 ) / 2252, 3660, 20,
     $                   3609 /
      DATA               ( MM( 123, J ), J = 1, 4 ) / 3904, 716, 3285,
     $                   3821 /
      DATA               ( MM( 124, J ), J = 1, 4 ) / 2774, 1842, 2046,
     $                   3537 /
      DATA               ( MM( 125, J ), J = 1, 4 ) / 997, 3987, 2107,
     $                   517 /
      DATA               ( MM( 126, J ), J = 1, 4 ) / 2573, 1368, 3508,
     $                   3017 /
      DATA               ( MM( 127, J ), J = 1, 4 ) / 1148, 1848, 3525,
     $                   2141 /
      DATA               ( MM( 128, J ), J = 1, 4 ) / 545, 2366, 3801,
     $                   1537 /
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      I1 = ISEED( 1 )
      I2 = ISEED( 2 )
      I3 = ISEED( 3 )
      I4 = ISEED( 4 )
*
      DO 10 I = 1, MIN( N, LV )
*
*        MULTIPLY THE SEED BY I-TH POWER OF THE MULTIPLIER MODULO 2**48
*
         IT4 = I4*MM( I, 4 )
         IT3 = IT4 / IPW2
         IT4 = IT4 - IPW2*IT3
         IT3 = IT3 + I3*MM( I, 4 ) + I4*MM( I, 3 )
         IT2 = IT3 / IPW2
         IT3 = IT3 - IPW2*IT2
         IT2 = IT2 + I2*MM( I, 4 ) + I3*MM( I, 3 ) + I4*MM( I, 2 )
         IT1 = IT2 / IPW2
         IT2 = IT2 - IPW2*IT1
         IT1 = IT1 + I1*MM( I, 4 ) + I2*MM( I, 3 ) + I3*MM( I, 2 ) +
     $         I4*MM( I, 1 )
         IT1 = MOD( IT1, IPW2 )
*
*        CONVERT 48-BIT INTEGER TO A REAL NUMBER IN THE INTERVAL (0,1)
*
         X( I ) = R*( DBLE( IT1 )+R*( DBLE( IT2 )+R*( DBLE( IT3 )+R*
     $            DBLE( IT4 ) ) ) )
   10 CONTINUE
*
*     RETURN FINAL VALUE OF SEED
*
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
      RETURN
*
*     END OF DLARUV
*
      END
CUT HERE............
CAT > DSTEBZ.F <<'CUT HERE............'
      SUBROUTINE DSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E,
     $                   M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK,
     $                   INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          ORDER, RANGE
      INTEGER            IL, INFO, IU, M, N, NSPLIT
      DOUBLE PRECISION   ABSTOL, VL, VU
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), W( * ), WORK( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DSTEBZ COMPUTES THE EIGENVALUES OF A SYMMETRIC TRIDIAGONAL
*  MATRIX T.  THE USER MAY ASK FOR ALL EIGENVALUES, ALL EIGENVALUES
*  IN THE HALF-OPEN INTERVAL (VL, VU], OR THE IL-TH THROUGH IU-TH
*  EIGENVALUES.
*
*  SEE W. KAHAN "ACCURATE EIGENVALUES OF A SYMMETRIC TRIDIAGONAL
*  MATRIX", REPORT CS41, COMPUTER SCIENCE DEPT., STANFORD
*  UNIVERSITY, JULY 21, 1966.
*
*  ARGUMENTS
*  =========
*
*  RANGE   (INPUT) CHARACTER
*          = 'A': ("ALL")   ALL EIGENVALUES WILL BE FOUND.
*          = 'V': ("VALUE") ALL EIGENVALUES IN THE HALF-OPEN INTERVAL
*                           (VL, VU] WILL BE FOUND.
*          = 'I': ("INDEX") THE IL-TH THROUGH IU-TH EIGENVALUES (OF THE
*                           ENTIRE MATRIX) WILL BE FOUND.
*
*  ORDER   (INPUT) CHARACTER
*          = 'B': ("BY BLOCK") THE EIGENVALUES WILL BE GROUPED BY
*                              SPLIT-OFF BLOCK (SEE IBLOCK, ISPLIT) AND
*                              ORDERED FROM SMALLEST TO LARGEST WITHIN
*                              THE BLOCK.
*          = 'E': ("ENTIRE MATRIX")
*                              THE EIGENVALUES FOR THE ENTIRE MATRIX
*                              WILL BE ORDERED FROM SMALLEST TO
*                              LARGEST.
*
*  N       (INPUT) INTEGER
*          THE DIMENSION OF THE TRIDIAGONAL MATRIX T.  N >= 0.
*
*  VL      (INPUT) DOUBLE PRECISION
*          IF RANGE='V', THE LOWER BOUND OF THE INTERVAL TO BE SEARCHED
*          FOR EIGENVALUES.  EIGENVALUES LESS THAN OR EQUAL TO VL WILL
*          NOT BE RETURNED.  NOT REFERENCED IF RANGE='A' OR 'I'.
*
*  VU      (INPUT) DOUBLE PRECISION
*          IF RANGE='V', THE UPPER BOUND OF THE INTERVAL TO BE SEARCHED
*          FOR EIGENVALUES.  EIGENVALUES GREATER THAN VU WILL NOT BE
*          RETURNED.  VU MUST BE GREATER THAN VL.  NOT REFERENCED IF
*          RANGE='A' OR 'I'.
*
*  IL      (INPUT) INTEGER
*          IF RANGE='I', THE INDEX (FROM SMALLEST TO LARGEST) OF THE
*          SMALLEST EIGENVALUE TO BE RETURNED.  IL MUST BE AT LEAST 1.
*          NOT REFERENCED IF RANGE='A' OR 'V'.
*
*  IU      (INPUT) INTEGER
*          IF RANGE='I', THE INDEX (FROM SMALLEST TO LARGEST) OF THE
*          LARGEST EIGENVALUE TO BE RETURNED.  IU MUST BE AT LEAST IL
*          AND NO GREATER THAN N.  NOT REFERENCED IF RANGE='A' OR 'V'.
*
*  ABSTOL  (INPUT) DOUBLE PRECISION
*          THE ABSOLUTE TOLERANCE FOR THE EIGENVALUES.  AN EIGENVALUE
*          (OR CLUSTER) IS CONSIDERED TO BE LOCATED IF IT HAS BEEN
*          DETERMINED TO LIE IN AN INTERVAL WHOSE WIDTH IS ABSTOL OR
*          LESS.  IF ABSTOL IS LESS THAN OR EQUAL TO ZERO, THEN ULP*|T|
*          WILL BE USED, WHERE |T| MEANS THE 1-NORM OF T.
*
*  D       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          THE N DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX T.  TO
*          AVOID OVERFLOW, THE MATRIX MUST BE SCALED SO THAT ITS LARGEST
*          ENTRY IS NO GREATER THAN OVERFLOW**(1/2) * UNDERFLOW**(1/4)
*          IN ABSOLUTE VALUE, AND FOR GREATEST ACCURACY, IT SHOULD NOT
*          BE MUCH SMALLER THAN THAT.
*
*  E       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          THE (N-1) OFF-DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX T.
*          TO AVOID OVERFLOW, THE MATRIX MUST BE SCALED SO THAT ITS
*          LARGEST ENTRY IS NO GREATER THAN OVERFLOW**(1/2) *
*          UNDERFLOW**(1/4) IN ABSOLUTE VALUE, AND FOR GREATEST
*          ACCURACY, IT SHOULD NOT BE MUCH SMALLER THAN THAT.
*
*  M       (OUTPUT) INTEGER
*          THE ACTUAL NUMBER OF EIGENVALUES FOUND. 0 <= M <= N.
*          (SEE ALSO THE DESCRIPTION OF INFO=2,3.)
*
*  NSPLIT  (OUTPUT) INTEGER
*          THE NUMBER OF DIAGONAL BLOCKS IN THE MATRIX T.
*          1 <= NSPLIT <= N.
*
*  W       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          ON EXIT, THE FIRST M ELEMENTS OF W WILL CONTAIN THE
*          EIGENVALUES.  (DSTEBZ MAY USE THE REMAINING N-M ELEMENTS AS
*          WORKSPACE.)
*
*  IBLOCK  (OUTPUT) INTEGER ARRAY, DIMENSION (N)
*          AT EACH ROW/COLUMN J WHERE E(J) IS ZERO OR SMALL, THE
*          MATRIX T IS CONSIDERED TO SPLIT INTO A BLOCK DIAGONAL
*          MATRIX.  ON EXIT, IBLOCK(I) SPECIFIES WHICH BLOCK (FROM 1 TO
*          THE NUMBER OF BLOCKS) THE EIGENVALUE W(I) BELONGS TO.
*          (DSTEBZ MAY USE THE REMAINING N-M ELEMENTS AS WORKSPACE.)
*
*  ISPLIT  (OUTPUT) INTEGER ARRAY, DIMENSION (N)
*          THE SPLITTING POINTS, AT WHICH T BREAKS UP INTO SUBMATRICES.
*          THE FIRST SUBMATRIX CONSISTS OF ROWS/COLUMNS 1 TO ISPLIT(1),
*          THE SECOND OF ROWS/COLUMNS ISPLIT(1)+1 THROUGH ISPLIT(2),
*          ETC., AND THE NSPLIT-TH CONSISTS OF ROWS/COLUMNS
*          ISPLIT(NSPLIT-1)+1 THROUGH ISPLIT(NSPLIT)=N.
*          (ONLY THE FIRST NSPLIT ELEMENTS WILL ACTUALLY BE USED, BUT
*          SINCE THE USER CANNOT KNOW A PRIORI WHAT VALUE NSPLIT WILL
*          HAVE, N WORDS MUST BE RESERVED FOR ISPLIT.)
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (4*N)
*
*  IWORK   (WORKSPACE) INTEGER ARRAY, DIMENSION (3*N)
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0:  SOME OR ALL OF THE EIGENVALUES FAILED TO CONVERGE OR
*                WERE NOT COMPUTED:
*                =1 OR 3: BISECTION FAILED TO CONVERGE FOR SOME
*                        EIGENVALUES; THESE EIGENVALUES ARE FLAGGED BY A
*                        NEGATIVE BLOCK NUMBER.  THE EFFECT IS THAT THE
*                        EIGENVALUES MAY NOT BE AS ACCURATE AS THE
*                        ABSOLUTE AND RELATIVE TOLERANCES.  THIS IS
*                        GENERALLY CAUSED BY UNEXPECTEDLY INACCURATE
*                        ARITHMETIC.
*                =2 OR 3: RANGE='I' ONLY: NOT ALL OF THE EIGENVALUES
*                        IL:IU WERE FOUND.
*                        EFFECT: M < IU+1-IL
*                        CAUSE:  NON-MONOTONIC ARITHMETIC, CAUSING THE
*                                STURM SEQUENCE TO BE NON-MONOTONIC.
*                        CURE:   RECALCULATE, USING RANGE='A', AND PICK
*                                OUT EIGENVALUES IL:IU.  IN SOME CASES,
*                                INCREASING THE PARAMETER "FUDGE" MAY
*                                MAKE THINGS WORK.
*                = 4:    RANGE='I', AND THE GERSHGORIN INTERVAL
*                        INITIALLY USED WAS TOO SMALL.  NO EIGENVALUES
*                        WERE COMPUTED.
*                        PROBABLE CAUSE: YOUR MACHINE HAS SLOPPY
*                                        FLOATING-POINT ARITHMETIC.
*                        CURE: INCREASE THE PARAMETER "FUDGE",
*                              RECOMPILE, AND TRY AGAIN.
*
*  INTERNAL PARAMETERS
*  ===================
*
*  RELFAC  DOUBLE PRECISION, DEFAULT = 2.0E0
*          THE RELATIVE TOLERANCE.  AN INTERVAL (A,B] LIES WITHIN
*          "RELATIVE TOLERANCE" IF  B-A < RELFAC*ULP*MAX(|A|,|B|),
*          WHERE "ULP" IS THE MACHINE PRECISION (DISTANCE FROM 1 TO
*          THE NEXT LARGER FLOATING POINT NUMBER.)
*
*  FUDGE   DOUBLE PRECISION, DEFAULT = 2
*          A "FUDGE FACTOR" TO WIDEN THE GERSHGORIN INTERVALS.  IDEALLY,
*          A VALUE OF 1 SHOULD WORK, BUT ON MACHINES WITH SLOPPY
*          ARITHMETIC, THIS NEEDS TO BE LARGER.  THE DEFAULT FOR
*          PUBLICLY RELEASED VERSIONS SHOULD BE LARGE ENOUGH TO HANDLE
*          THE WORST MACHINE AROUND.  NOTE THAT THIS HAS NO EFFECT
*          ON ACCURACY OF THE SOLUTION.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO, ONE, TWO, HALF
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   HALF = 1.0D0 / TWO )
      DOUBLE PRECISION   FUDGE, RELFAC
      PARAMETER          ( FUDGE = 2.0D0, RELFAC = 2.0D0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            NCNVRG, TOOFEW
      INTEGER            IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO,
     $                   IM, IN, IOFF, IORDER, IOUT, IRANGE, ITMAX,
     $                   ITMP1, IW, IWOFF, J, JB, JDISC, JE, NB, NWL,
     $                   NWU
      DOUBLE PRECISION   ATOLI, BNORM, GL, GU, PIVMIN, RTOLI, SAFEMN,
     $                   TMP1, TMP2, TNORM, ULP, WKILL, WL, WLU, WU, WUL
*     ..
*     .. LOCAL ARRAYS ..
      INTEGER            IDUMMA( 1 )
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, ILAENV, DLAMCH
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLAEBZ, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, INT, LOG, MAX, MIN, SQRT
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      INFO = 0
*
*     DECODE RANGE
*
      IF( LSAME( RANGE, 'A' ) ) THEN
         IRANGE = 1
      ELSE IF( LSAME( RANGE, 'V' ) ) THEN
         IRANGE = 2
      ELSE IF( LSAME( RANGE, 'I' ) ) THEN
         IRANGE = 3
      ELSE
         IRANGE = 0
      END IF
*
*     DECODE ORDER
*
      IF( LSAME( ORDER, 'B' ) ) THEN
         IORDER = 2
      ELSE IF( LSAME( ORDER, 'E' ) .OR. LSAME( ORDER, 'A' ) ) THEN
         IORDER = 1
      ELSE
         IORDER = 0
      END IF
*
*     CHECK FOR ERRORS
*
      IF( IRANGE.LE.0 ) THEN
         INFO = -1
      ELSE IF( IORDER.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( IRANGE.EQ.2 .AND. VL.GE.VU ) THEN
         INFO = -5
      ELSE IF( IRANGE.EQ.3 .AND. ( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) )
     $          THEN
         INFO = -6
      ELSE IF( IRANGE.EQ.3 .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) )
     $          THEN
         INFO = -7
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEBZ', -INFO )
         RETURN
      END IF
*
*     INITIALIZE ERROR FLAGS
*
      INFO = 0
      NCNVRG = .FALSE.
      TOOFEW = .FALSE.
*
*     QUICK RETURN IF POSSIBLE
*
      M = 0
      IF( N.EQ.0 )
     $   RETURN
*
*     SIMPLIFICATIONS:
*
      IF( IRANGE.EQ.3 .AND. IL.EQ.1 .AND. IU.EQ.N )
     $   IRANGE = 1
*
*     GET MACHINE CONSTANTS
*     NB IS THE MINIMUM VECTOR LENGTH FOR VECTOR BISECTION, OR 0
*     IF ONLY SCALAR IS TO BE DONE.
*
      SAFEMN = DLAMCH( 'S' )
      ULP = DLAMCH( 'P' )
      RTOLI = ULP*RELFAC
      NB = ILAENV( 1, 'DSTEBZ', ' ', N, -1, -1, -1 )
      IF( NB.LE.1 )
     $   NB = 0
*
*     SPECIAL CASE WHEN N=1
*
      IF( N.EQ.1 ) THEN
         NSPLIT = 1
         ISPLIT( 1 ) = 1
         IF( IRANGE.EQ.2 .AND. ( VL.GE.D( 1 ) .OR. VU.LT.D( 1 ) ) ) THEN
            M = 0
         ELSE
            W( 1 ) = D( 1 )
            IBLOCK( 1 ) = 1
            M = 1
         END IF
         RETURN
      END IF
*
*     COMPUTE SPLITTING POINTS
*
      NSPLIT = 1
      WORK( N ) = ZERO
      PIVMIN = ONE
*
      DO 10 J = 2, N
         TMP1 = E( J-1 )**2
         IF( ABS( D( J )*D( J-1 ) )*ULP**2+SAFEMN.GT.TMP1 ) THEN
            ISPLIT( NSPLIT ) = J - 1
            NSPLIT = NSPLIT + 1
            WORK( J-1 ) = ZERO
         ELSE
            WORK( J-1 ) = TMP1
            PIVMIN = MAX( PIVMIN, TMP1 )
         END IF
   10 CONTINUE
      ISPLIT( NSPLIT ) = N
      PIVMIN = PIVMIN*SAFEMN
*
*     COMPUTE INTERVAL AND ATOLI
*
      IF( IRANGE.EQ.3 ) THEN
*
*        RANGE='I': COMPUTE THE INTERVAL CONTAINING EIGENVALUES
*                   IL THROUGH IU.
*
*        COMPUTE GERSHGORIN INTERVAL FOR ENTIRE (SPLIT) MATRIX
*        AND USE IT AS THE INITIAL INTERVAL
*
         GU = D( 1 )
         GL = D( 1 )
         TMP1 = ZERO
*
         DO 20 J = 1, N - 1
            TMP2 = SQRT( WORK( J ) )
            GU = MAX( GU, D( J )+TMP1+TMP2 )
            GL = MIN( GL, D( J )-TMP1-TMP2 )
            TMP1 = TMP2
   20    CONTINUE
*
         GU = MAX( GU, D( N )+TMP1 )
         GL = MIN( GL, D( N )-TMP1 )
         TNORM = MAX( ABS( GL ), ABS( GU ) )
         GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN
         GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN
*
*        COMPUTE ITERATION PARAMETERS
*
         ITMAX = INT( ( LOG( TNORM+PIVMIN )-LOG( PIVMIN ) ) /
     $           LOG( TWO ) ) + 2
         IF( ABSTOL.LE.ZERO ) THEN
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         END IF
*
         WORK( N+1 ) = GL
         WORK( N+2 ) = GL
         WORK( N+3 ) = GU
         WORK( N+4 ) = GU
         WORK( N+5 ) = GL
         WORK( N+6 ) = GU
         IWORK( 1 ) = -1
         IWORK( 2 ) = -1
         IWORK( 3 ) = N + 1
         IWORK( 4 ) = N + 1
         IWORK( 5 ) = IL - 1
         IWORK( 6 ) = IU
*
         CALL DLAEBZ( 3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E,
     $                WORK, IWORK( 5 ), WORK( N+1 ), WORK( N+5 ), IOUT,
     $                IWORK, W, IBLOCK, IINFO )
*
         IF( IWORK( 6 ).EQ.IU ) THEN
            WL = WORK( N+1 )
            WLU = WORK( N+3 )
            NWL = IWORK( 1 )
            WU = WORK( N+4 )
            WUL = WORK( N+2 )
            NWU = IWORK( 4 )
         ELSE
            WL = WORK( N+2 )
            WLU = WORK( N+4 )
            NWL = IWORK( 2 )
            WU = WORK( N+3 )
            WUL = WORK( N+1 )
            NWU = IWORK( 3 )
         END IF
*
         IF( NWL.LT.0 .OR. NWL.GE.N .OR. NWU.LT.1 .OR. NWU.GT.N ) THEN
            INFO = 4
            RETURN
         END IF
      ELSE
*
*        RANGE='A' OR 'V' -- SET ATOLI
*
         TNORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ),
     $           ABS( D( N ) )+ABS( E( N-1 ) ) )
*
         DO 30 J = 2, N - 1
            TNORM = MAX( TNORM, ABS( D( J ) )+ABS( E( J-1 ) )+
     $              ABS( E( J ) ) )
   30    CONTINUE
*
         IF( ABSTOL.LE.ZERO ) THEN
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         END IF
*
         IF( IRANGE.EQ.2 ) THEN
            WL = VL
            WU = VU
         END IF
      END IF
*
*     FIND EIGENVALUES -- LOOP OVER BLOCKS AND RECOMPUTE NWL AND NWU.
*     NWL ACCUMULATES THE NUMBER OF EIGENVALUES .LE. WL,
*     NWU ACCUMULATES THE NUMBER OF EIGENVALUES .LE. WU
*
      M = 0
      IEND = 0
      INFO = 0
      NWL = 0
      NWU = 0
*
      DO 70 JB = 1, NSPLIT
         IOFF = IEND
         IBEGIN = IOFF + 1
         IEND = ISPLIT( JB )
         IN = IEND - IOFF
*
         IF( IN.EQ.1 ) THEN
*
*           SPECIAL CASE -- IN=1
*
            IF( IRANGE.EQ.1 .OR. WL.GE.D( IBEGIN )-PIVMIN )
     $         NWL = NWL + 1
            IF( IRANGE.EQ.1 .OR. WU.GE.D( IBEGIN )-PIVMIN )
     $         NWU = NWU + 1
            IF( IRANGE.EQ.1 .OR. ( WL.LT.D( IBEGIN )-PIVMIN .AND. WU.GE.
     $          D( IBEGIN )-PIVMIN ) ) THEN
               M = M + 1
               W( M ) = D( IBEGIN )
               IBLOCK( M ) = JB
            END IF
         ELSE
*
*           GENERAL CASE -- IN > 1
*
*           COMPUTE GERSHGORIN INTERVAL
*           AND USE IT AS THE INITIAL INTERVAL
*
            GU = D( IBEGIN )
            GL = D( IBEGIN )
            TMP1 = ZERO
*
            DO 40 J = IBEGIN, IEND - 1
               TMP2 = ABS( E( J ) )
               GU = MAX( GU, D( J )+TMP1+TMP2 )
               GL = MIN( GL, D( J )-TMP1-TMP2 )
               TMP1 = TMP2
   40       CONTINUE
*
            GU = MAX( GU, D( IEND )+TMP1 )
            GL = MIN( GL, D( IEND )-TMP1 )
            BNORM = MAX( ABS( GL ), ABS( GU ) )
            GL = GL - FUDGE*BNORM*ULP*IN - FUDGE*PIVMIN
            GU = GU + FUDGE*BNORM*ULP*IN + FUDGE*PIVMIN
*
            IF( IRANGE.GT.1 ) THEN
               GL = MAX( GL, WL )
               GU = MIN( GU, WU )
               IF( GL.GE.GU )
     $            GO TO 70
            END IF
*
*           SET UP INITIAL INTERVAL
*
            WORK( N+1 ) = GL
            WORK( N+IN+1 ) = GU
            CALL DLAEBZ( 1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
     $                   D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),
     $                   IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IM,
     $                   IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
*
            NWL = NWL + IWORK( 1 )
            NWU = NWU + IWORK( IN+1 )
            IWOFF = M - IWORK( 1 )
*
*           COMPUTE EIGENVALUES
*
            ITMAX = INT( ( LOG( GU-GL+PIVMIN )-LOG( PIVMIN ) ) /
     $              LOG( TWO ) ) + 2
            CALL DLAEBZ( 2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
     $                   D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),
     $                   IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IOUT,
     $                   IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
*
*           COPY EIGENVALUES INTO W AND IBLOCK
*           USE -JB FOR BLOCK NUMBER FOR UNCONVERGED EIGENVALUES.
*
            DO 60 J = 1, IOUT
               TMP1 = HALF*( WORK( J+N )+WORK( J+IN+N ) )
*
*              FLAG NON-CONVERGENCE.
*
               IF( J.GT.IOUT-IINFO ) THEN
                  NCNVRG = .TRUE.
                  IB = -JB
               ELSE
                  IB = JB
               END IF
               DO 50 JE = IWORK( J ) + 1 + IWOFF,
     $                 IWORK( J+IN ) + IWOFF
                  W( JE ) = TMP1
                  IBLOCK( JE ) = IB
   50          CONTINUE
   60       CONTINUE
*
            M = M + IM
         END IF
   70 CONTINUE
*
*     IF RANGE='I', THEN (WL,WU) CONTAINS EIGENVALUES NWL+1,...,NWU
*     IF NWL+1 < IL OR NWU > IU, DISCARD EXTRA EIGENVALUES.
*
      IF( IRANGE.EQ.3 ) THEN
         IM = 0
         IDISCL = IL - 1 - NWL
         IDISCU = NWU - IU
*
         IF( IDISCL.GT.0 .OR. IDISCU.GT.0 ) THEN
            DO 80 JE = 1, M
               IF( W( JE ).LE.WLU .AND. IDISCL.GT.0 ) THEN
                  IDISCL = IDISCL - 1
               ELSE IF( W( JE ).GE.WUL .AND. IDISCU.GT.0 ) THEN
                  IDISCU = IDISCU - 1
               ELSE
                  IM = IM + 1
                  W( IM ) = W( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               END IF
   80       CONTINUE
            M = IM
         END IF
         IF( IDISCL.GT.0 .OR. IDISCU.GT.0 ) THEN
*
*           CODE TO DEAL WITH EFFECTS OF BAD ARITHMETIC:
*           SOME LOW EIGENVALUES TO BE DISCARDED ARE NOT IN (WL,WLU],
*           OR HIGH EIGENVALUES TO BE DISCARDED ARE NOT IN (WUL,WU]
*           SO JUST KILL OFF THE SMALLEST IDISCL/LARGEST IDISCU
*           EIGENVALUES, BY SIMPLY FINDING THE SMALLEST/LARGEST
*           EIGENVALUE(S).
*
*           (IF N(W) IS MONOTONE NON-DECREASING, THIS SHOULD NEVER
*               HAPPEN.)
*
            IF( IDISCL.GT.0 ) THEN
               WKILL = WU
               DO 100 JDISC = 1, IDISCL
                  IW = 0
                  DO 90 JE = 1, M
                     IF( IBLOCK( JE ).NE.0 .AND.
     $                   ( W( JE ).LT.WKILL .OR. IW.EQ.0 ) ) THEN
                        IW = JE
                        WKILL = W( JE )
                     END IF
   90             CONTINUE
                  IBLOCK( IW ) = 0
  100          CONTINUE
            END IF
            IF( IDISCU.GT.0 ) THEN
*
               WKILL = WL
               DO 120 JDISC = 1, IDISCU
                  IW = 0
                  DO 110 JE = 1, M
                     IF( IBLOCK( JE ).NE.0 .AND.
     $                   ( W( JE ).GT.WKILL .OR. IW.EQ.0 ) ) THEN
                        IW = JE
                        WKILL = W( JE )
                     END IF
  110             CONTINUE
                  IBLOCK( IW ) = 0
  120          CONTINUE
            END IF
            IM = 0
            DO 130 JE = 1, M
               IF( IBLOCK( JE ).NE.0 ) THEN
                  IM = IM + 1
                  W( IM ) = W( JE )
                  IBLOCK( IM ) = IBLOCK( JE )
               END IF
  130       CONTINUE
            M = IM
         END IF
         IF( IDISCL.LT.0 .OR. IDISCU.LT.0 ) THEN
            TOOFEW = .TRUE.
         END IF
      END IF
*
*     IF ORDER='B', DO NOTHING -- THE EIGENVALUES ARE ALREADY SORTED
*        BY BLOCK.
*     IF ORDER='E' OR 'A', SORT THE EIGENVALUES FROM SMALLEST TO LARGEST
*
      IF( IORDER.EQ.1 .AND. NSPLIT.GT.1 ) THEN
         DO 150 JE = 1, M - 1
            IE = 0
            TMP1 = W( JE )
            DO 140 J = JE + 1, M
               IF( W( J ).LT.TMP1 ) THEN
                  IE = J
                  TMP1 = W( J )
               END IF
  140       CONTINUE
*
            IF( IE.NE.0 ) THEN
               ITMP1 = IBLOCK( IE )
               W( IE ) = W( JE )
               IBLOCK( IE ) = IBLOCK( JE )
               W( JE ) = TMP1
               IBLOCK( JE ) = ITMP1
            END IF
  150    CONTINUE
      END IF
*
      INFO = 0
      IF( NCNVRG )
     $   INFO = INFO + 1
      IF( TOOFEW )
     $   INFO = INFO + 2
      RETURN
*
*     END OF DSTEBZ
*
      END
CUT HERE............
CAT > DLAEBZ.F <<'CUT HERE............'
      SUBROUTINE DLAEBZ( IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL,
     $                   RELTOL, PIVMIN, D, E, E2, NVAL, AB, C, MOUT,
     $                   NAB, WORK, IWORK, INFO )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            IJOB, INFO, MINP, MMAX, MOUT, N, NBMIN, NITMAX
      DOUBLE PRECISION   ABSTOL, PIVMIN, RELTOL
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IWORK( * ), NAB( MMAX, * ), NVAL( * )
      DOUBLE PRECISION   AB( MMAX, * ), C( * ), D( * ), E( * ), E2( * ),
     $                   WORK( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLAEBZ CONTAINS THE ITERATION LOOPS WHICH COMPUTE AND USE THE
*  FUNCTION N(W), WHICH IS THE COUNT OF EIGENVALUES OF A SYMMETRIC
*  TRIDIAGONAL MATRIX T LESS THAN OR EQUAL TO ITS ARGUMENT  W.  IT
*  PERFORMS A CHOICE OF TWO TYPES OF LOOPS:
*
*  IJOB=1, FOLLOWED BY
*  IJOB=2: IT TAKES AS INPUT A LIST OF INTERVALS AND RETURNS A LIST OF
*          SUFFICIENTLY SMALL INTERVALS WHOSE UNION CONTAINS THE SAME
*          EIGENVALUES AS THE UNION OF THE ORIGINAL INTERVALS.
*          THE INPUT INTERVALS ARE (AB(J,1),AB(J,2)], J=1,...,MINP.
*          THE OUTPUT INTERVAL (AB(J,1),AB(J,2)] WILL CONTAIN
*          EIGENVALUES NAB(J,1)+1,...,NAB(J,2), WHERE 1 <= J <= MOUT.
*
*  IJOB=3: IT PERFORMS A BINARY SEARCH IN EACH INPUT INTERVAL
*          (AB(J,1),AB(J,2)] FOR A POINT  W(J)  SUCH THAT
*          N(W(J))=NVAL(J), AND USES  C(J)  AS THE STARTING POINT OF
*          THE SEARCH.  IF SUCH A W(J) IS FOUND, THEN ON OUTPUT
*          AB(J,1)=AB(J,2)=W.  IF NO SUCH W(J) IS FOUND, THEN ON OUTPUT
*          (AB(J,1),AB(J,2)] WILL BE A SMALL INTERVAL CONTAINING THE
*          POINT WHERE N(W) JUMPS THROUGH NVAL(J), UNLESS THAT POINT
*          LIES OUTSIDE THE INITIAL INTERVAL.
*
*  NOTE THAT THE INTERVALS ARE IN ALL CASES HALF-OPEN INTERVALS,
*  I.E., OF THE FORM  (A,B] , WHICH INCLUDES  B  BUT NOT  A .
*
*  SEE W. KAHAN "ACCURATE EIGENVALUES OF A SYMMETRIC TRIDIAGONAL
*  MATRIX", REPORT CS41, COMPUTER SCIENCE DEPT., STANFORD
*  UNIVERSITY, JULY 21, 1966
*
*  NOTE: THE ARGUMENTS ARE, IN GENERAL, *NOT* CHECKED FOR UNREASONABLE
*  VALUES.
*
*  ARGUMENTS
*  =========
*
*  IJOB    (INPUT) INTEGER
*          SPECIFIES WHAT IS TO BE DONE:
*          = 1:  COMPUTE NAB FOR THE INITIAL INTERVALS.
*          = 2:  PERFORM BISECTION ITERATION TO FIND EIGENVALUES OF T.
*          = 3:  PERFORM BISECTION ITERATION TO INVERT N(W), I.E.,
*                TO FIND A POINT WHICH HAS A SPECIFIED NUMBER OF
*                EIGENVALUES OF T TO ITS LEFT.
*          OTHER VALUES WILL CAUSE DLAEBZ TO RETURN WITH INFO=-1.
*
*  NITMAX  (INPUT) INTEGER
*          THE MAXIMUM NUMBER OF "LEVELS" OF BISECTION TO BE
*          PERFORMED, I.E., AN INTERVAL OF WIDTH W WILL NOT BE MADE
*          SMALLER THAN 2^(-NITMAX) * W.  IF NOT ALL INTERVALS
*          HAVE CONVERGED AFTER NITMAX ITERATIONS, THEN INFO IS SET
*          TO THE NUMBER OF NON-CONVERGED INTERVALS.
*
*  N       (INPUT) INTEGER
*          THE DIMENSION N OF THE TRIDIAGONAL MATRIX T.  IT MUST BE AT
*          LEAST 1.
*
*  MMAX    (INPUT) INTEGER
*          THE MAXIMUM NUMBER OF INTERVALS.  IF MORE THAN MMAX INTERVALS
*          ARE GENERATED, THEN DLAEBZ WILL QUIT WITH INFO=MMAX+1.
*
*  MINP    (INPUT) INTEGER
*          THE INITIAL NUMBER OF INTERVALS.  IT MAY NOT BE GREATER THAN
*          MMAX.
*
*  NBMIN   (INPUT) INTEGER
*          THE SMALLEST NUMBER OF INTERVALS THAT SHOULD BE PROCESSED
*          USING A VECTOR LOOP.  IF ZERO, THEN ONLY THE SCALAR LOOP
*          WILL BE USED.
*
*  ABSTOL  (INPUT) DOUBLE PRECISION
*          THE MINIMUM (ABSOLUTE) WIDTH OF AN INTERVAL.  WHEN AN
*          INTERVAL IS NARROWER THAN ABSTOL, OR THAN RELTOL TIMES THE
*          LARGER (IN MAGNITUDE) ENDPOINT, THEN IT IS CONSIDERED TO BE
*          SUFFICIENTLY SMALL, I.E., CONVERGED.  THIS MUST BE AT LEAST
*          ZERO.
*
*  RELTOL  (INPUT) DOUBLE PRECISION
*          THE MINIMUM RELATIVE WIDTH OF AN INTERVAL.  WHEN AN INTERVAL
*          IS NARROWER THAN ABSTOL, OR THAN RELTOL TIMES THE LARGER (IN
*          MAGNITUDE) ENDPOINT, THEN IT IS CONSIDERED TO BE
*          SUFFICIENTLY SMALL, I.E., CONVERGED.  NOTE: THIS SHOULD
*          ALWAYS BE AT LEAST RADIX*MACHINE EPSILON.
*
*  PIVMIN  (INPUT) DOUBLE PRECISION
*          THE MINIMUM ABSOLUTE VALUE OF A "PIVOT" IN THE STURM
*          SEQUENCE LOOP.  THIS *MUST* BE AT LEAST  MAX |E(J)**2| *
*          SAFE_MIN  AND AT LEAST SAFE_MIN, WHERE SAFE_MIN IS AT LEAST
*          THE SMALLEST NUMBER THAT CAN DIVIDE ONE WITHOUT OVERFLOW.
*
*  D       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX T.  TO AVOID
*          UNDERFLOW, THE MATRIX SHOULD BE SCALED SO THAT ITS LARGEST
*          ENTRY IS NO GREATER THAN  OVERFLOW**(1/2) * UNDERFLOW**(1/4)
*          IN ABSOLUTE VALUE.  TO ASSURE THE MOST ACCURATE COMPUTATION
*          OF SMALL EIGENVALUES, THE MATRIX SHOULD BE SCALED TO BE
*          NOT MUCH SMALLER THAN THAT, EITHER.
*
*  E       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          THE OFFDIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX T IN
*          POSITIONS 1 THROUGH N-1.  E(N) IS ARBITRARY.
*          TO AVOID UNDERFLOW, THE
*          MATRIX SHOULD BE SCALED SO THAT ITS LARGEST ENTRY IS NO
*          GREATER THAN  OVERFLOW**(1/2) * UNDERFLOW**(1/4) IN ABSOLUTE
*          VALUE.  TO ASSURE THE MOST ACCURATE COMPUTATION OF SMALL
*          EIGENVALUES, THE MATRIX SHOULD BE SCALED TO BE NOT MUCH
*          SMALLER THAN THAT, EITHER.
*
*  E2      (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          THE SQUARES OF THE OFFDIAGONAL ELEMENTS OF THE TRIDIAGONAL
*          MATRIX T.  E2(N) IS IGNORED.
*
*  NVAL    (INPUT/OUTPUT) INTEGER ARRAY, DIMENSION (MINP)
*          IF IJOB=1 OR 2, NOT REFERENCED.
*          IF IJOB=3, THE DESIRED VALUES OF N(W).  THE ELEMENTS OF NVAL
*          WILL BE REORDERED TO CORRESPOND WITH THE INTERVALS IN AB.
*          THUS, NVAL(J) ON OUTPUT WILL NOT, IN GENERAL BE THE SAME AS
*          NVAL(J) ON INPUT, BUT IT WILL CORRESPOND WITH THE INTERVAL
*          (AB(J,1),AB(J,2)] ON OUTPUT.
*
*  AB      (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (MMAX,2)
*          THE ENDPOINTS OF THE INTERVALS.  AB(J,1) IS  A(J), THE LEFT
*          ENDPOINT OF THE J-TH INTERVAL, AND AB(J,2) IS B(J), THE
*          RIGHT ENDPOINT OF THE J-TH INTERVAL.  THE INPUT INTERVALS
*          WILL, IN GENERAL, BE MODIFIED, SPLIT, AND REORDERED BY THE
*          CALCULATION.
*
*  C       (INPUT/WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (MMAX)
*          IF IJOB=1, IGNORED.
*          IF IJOB=2, WORKSPACE.
*          IF IJOB=3, THEN ON INPUT C(J) SHOULD BE INITIALIZED TO THE
*          FIRST SEARCH POINT IN THE BINARY SEARCH.
*
*  MOUT    (OUTPUT) INTEGER
*          IF IJOB=1, THE NUMBER OF EIGENVALUES IN THE INTERVALS.
*          IF IJOB=2 OR 3, THE NUMBER OF INTERVALS OUTPUT.
*          IF IJOB=3, MOUT WILL EQUAL MINP.
*
*  NAB     (INPUT/OUTPUT) INTEGER ARRAY, DIMENSION (MMAX,2)
*          IF IJOB=1, THEN ON OUTPUT NAB(I,J) WILL BE SET TO N(AB(I,J)).
*          IF IJOB=2, THEN ON INPUT, NAB(I,J) SHOULD BE SET.  IT MUST
*             SATISFY THE CONDITION:
*             N(AB(I,1)) <= NAB(I,1) <= NAB(I,2) <= N(AB(I,2)),
*             WHICH MEANS THAT IN INTERVAL I ONLY EIGENVALUES
*             NAB(I,1)+1,...,NAB(I,2) WILL BE CONSIDERED.  USUALLY,
*             NAB(I,J)=N(AB(I,J)), FROM A PREVIOUS CALL TO DLAEBZ WITH
*             IJOB=1.
*             ON OUTPUT, NAB(I,J) WILL CONTAIN
*             MAX(NA(K),MIN(NB(K),N(AB(I,J)))), WHERE K IS THE INDEX OF
*             THE INPUT INTERVAL THAT THE OUTPUT INTERVAL
*             (AB(J,1),AB(J,2)] CAME FROM, AND NA(K) AND NB(K) ARE THE
*             THE INPUT VALUES OF NAB(K,1) AND NAB(K,2).
*          IF IJOB=3, THEN ON OUTPUT, NAB(I,J) CONTAINS N(AB(I,J)),
*             UNLESS N(W) > NVAL(I) FOR ALL SEARCH POINTS  W , IN WHICH
*             CASE NAB(I,1) WILL NOT BE MODIFIED, I.E., THE OUTPUT
*             VALUE WILL BE THE SAME AS THE INPUT VALUE (MODULO
*             REORDERINGS -- SEE NVAL AND AB), OR UNLESS N(W) < NVAL(I)
*             FOR ALL SEARCH POINTS  W , IN WHICH CASE NAB(I,2) WILL
*             NOT BE MODIFIED.  NORMALLY, NAB SHOULD BE SET TO SOME
*             DISTINCTIVE VALUE(S) BEFORE DLAEBZ IS CALLED.
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (MMAX)
*          WORKSPACE.
*
*  IWORK   (WORKSPACE) INTEGER ARRAY, DIMENSION (MMAX)
*          WORKSPACE.
*
*  INFO    (OUTPUT) INTEGER
*          = 0:       ALL INTERVALS CONVERGED.
*          = 1--MMAX: THE LAST INFO INTERVALS DID NOT CONVERGE.
*          = MMAX+1:  MORE THAN MMAX INTERVALS WERE GENERATED.
*
*  FURTHER DETAILS
*  ===============
*
*      THIS ROUTINE IS INTENDED TO BE CALLED ONLY BY OTHER LAPACK
*  ROUTINES, THUS THE INTERFACE IS LESS USER-FRIENDLY.  IT IS INTENDED
*  FOR TWO PURPOSES:
*
*  (A) FINDING EIGENVALUES.  IN THIS CASE, DLAEBZ SHOULD HAVE ONE OR
*      MORE INITIAL INTERVALS SET UP IN AB, AND DLAEBZ SHOULD BE CALLED
*      WITH IJOB=1.  THIS SETS UP NAB, AND ALSO COUNTS THE EIGENVALUES.
*      INTERVALS WITH NO EIGENVALUES WOULD USUALLY BE THROWN OUT AT
*      THIS POINT.  ALSO, IF NOT ALL THE EIGENVALUES IN AN INTERVAL I
*      ARE DESIRED, NAB(I,1) CAN BE INCREASED OR NAB(I,2) DECREASED.
*      FOR EXAMPLE, SET NAB(I,1)=NAB(I,2)-1 TO GET THE LARGEST
*      EIGENVALUE.  DLAEBZ IS THEN CALLED WITH IJOB=2 AND MMAX
*      NO SMALLER THAN THE VALUE OF MOUT RETURNED BY THE CALL WITH
*      IJOB=1.  AFTER THIS (IJOB=2) CALL, EIGENVALUES NAB(I,1)+1
*      THROUGH NAB(I,2) ARE APPROXIMATELY AB(I,1) (OR AB(I,2)) TO THE
*      TOLERANCE SPECIFIED BY ABSTOL AND RELTOL.
*
*  (B) FINDING AN INTERVAL (A',B'] CONTAINING EIGENVALUES W(F),...,W(L).
*      IN THIS CASE, START WITH A GERSHGORIN INTERVAL  (A,B).  SET UP
*      AB TO CONTAIN 2 SEARCH INTERVALS, BOTH INITIALLY (A,B).  ONE
*      NVAL ENTRY SHOULD CONTAIN  F-1  AND THE OTHER SHOULD CONTAIN  L
*      , WHILE C SHOULD CONTAIN A AND B, RESP.  NAB(I,1) SHOULD BE -1
*      AND NAB(I,2) SHOULD BE N+1, TO FLAG AN ERROR IF THE DESIRED
*      INTERVAL DOES NOT LIE IN (A,B).  DLAEBZ IS THEN CALLED WITH
*      IJOB=3.  ON EXIT, IF W(F-1) < W(F), THEN ONE OF THE INTERVALS --
*      J -- WILL HAVE AB(J,1)=AB(J,2) AND NAB(J,1)=NAB(J,2)=F-1, WHILE
*      IF, TO THE SPECIFIED TOLERANCE, W(F-K)=...=W(F+R), K > 0 AND R
*      >= 0, THEN THE INTERVAL WILL HAVE  N(AB(J,1))=NAB(J,1)=F-K AND
*      N(AB(J,2))=NAB(J,2)=F+R.  THE CASES W(L) < W(L+1) AND
*      W(L-R)=...=W(L+K) ARE HANDLED SIMILARLY.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO, TWO, HALF
      PARAMETER          ( ZERO = 0.0D0, TWO = 2.0D0,
     $                   HALF = 1.0D0 / TWO )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            ITMP1, ITMP2, J, JI, JIT, JP, KF, KFNEW, KL,
     $                   KLNEW
      DOUBLE PRECISION   TMP1, TMP2
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     CHECK FOR ERRORS
*
      INFO = 0
      IF( IJOB.LT.1 .OR. IJOB.GT.3 ) THEN
         INFO = -1
         RETURN
      END IF
*
*     INITIALIZE NAB
*
      IF( IJOB.EQ.1 ) THEN
*
*        COMPUTE THE NUMBER OF EIGENVALUES IN THE INITIAL INTERVALS.
*
         MOUT = 0
         DO 30 JI = 1, MINP
            DO 20 JP = 1, 2
               TMP1 = D( 1 ) - AB( JI, JP )
               IF( ABS( TMP1 ).LT.PIVMIN )
     $            TMP1 = -PIVMIN
               NAB( JI, JP ) = 0
               IF( TMP1.LE.ZERO )
     $            NAB( JI, JP ) = 1
*
               DO 10 J = 2, N
                  TMP1 = D( J ) - E2( J-1 ) / TMP1 - AB( JI, JP )
                  IF( ABS( TMP1 ).LT.PIVMIN )
     $               TMP1 = -PIVMIN
                  IF( TMP1.LE.ZERO )
     $               NAB( JI, JP ) = NAB( JI, JP ) + 1
   10          CONTINUE
   20       CONTINUE
            MOUT = MOUT + NAB( JI, 2 ) - NAB( JI, 1 )
   30    CONTINUE
         RETURN
      END IF
*
*     INITIALIZE FOR LOOP
*
*     KF AND KL HAVE THE FOLLOWING MEANING:
*        INTERVALS 1,...,KF-1 HAVE CONVERGED.
*        INTERVALS KF,...,KL  STILL NEED TO BE REFINED.
*
      KF = 1
      KL = MINP
*
*     IF IJOB=2, INITIALIZE C.
*     IF IJOB=3, USE THE USER-SUPPLIED STARTING POINT.
*
      IF( IJOB.EQ.2 ) THEN
         DO 40 JI = 1, MINP
            C( JI ) = HALF*( AB( JI, 1 )+AB( JI, 2 ) )
   40    CONTINUE
      END IF
*
*     ITERATION LOOP
*
      DO 130 JIT = 1, NITMAX
*
*        LOOP OVER INTERVALS
*
         IF( KL-KF+1.GE.NBMIN .AND. NBMIN.GT.0 ) THEN
*
*           BEGIN OF PARALLEL VERSION OF THE LOOP
*
            DO 60 JI = KF, KL
*
*              COMPUTE N(C), THE NUMBER OF EIGENVALUES LESS THAN C
*
               WORK( JI ) = D( 1 ) - C( JI )
               IWORK( JI ) = 0
               IF( WORK( JI ).LE.PIVMIN ) THEN
                  IWORK( JI ) = 1
                  WORK( JI ) = MIN( WORK( JI ), -PIVMIN )
               END IF
*
               DO 50 J = 2, N
                  WORK( JI ) = D( J ) - E2( J-1 ) / WORK( JI ) - C( JI )
                  IF( WORK( JI ).LE.PIVMIN ) THEN
                     IWORK( JI ) = IWORK( JI ) + 1
                     WORK( JI ) = MIN( WORK( JI ), -PIVMIN )
                  END IF
   50          CONTINUE
   60       CONTINUE
*
            IF( IJOB.LE.2 ) THEN
*
*              IJOB=2: CHOOSE ALL INTERVALS CONTAINING EIGENVALUES.
*
               KLNEW = KL
               DO 70 JI = KF, KL
*
*                 INSURE THAT N(W) IS MONOTONE
*
                  IWORK( JI ) = MIN( NAB( JI, 2 ),
     $                          MAX( NAB( JI, 1 ), IWORK( JI ) ) )
*
*                 UPDATE THE QUEUE -- ADD INTERVALS IF BOTH HALVES
*                 CONTAIN EIGENVALUES.
*
                  IF( IWORK( JI ).EQ.NAB( JI, 2 ) ) THEN
*
*                    NO EIGENVALUE IN THE UPPER INTERVAL:
*                    JUST USE THE LOWER INTERVAL.
*
                     AB( JI, 2 ) = C( JI )
*
                  ELSE IF( IWORK( JI ).EQ.NAB( JI, 1 ) ) THEN
*
*                    NO EIGENVALUE IN THE LOWER INTERVAL:
*                    JUST USE THE UPPER INTERVAL.
*
                     AB( JI, 1 ) = C( JI )
                  ELSE
                     KLNEW = KLNEW + 1
                     IF( KLNEW.LE.MMAX ) THEN
*
*                       EIGENVALUE IN BOTH INTERVALS -- ADD UPPER TO
*                       QUEUE.
*
                        AB( KLNEW, 2 ) = AB( JI, 2 )
                        NAB( KLNEW, 2 ) = NAB( JI, 2 )
                        AB( KLNEW, 1 ) = C( JI )
                        NAB( KLNEW, 1 ) = IWORK( JI )
                        AB( JI, 2 ) = C( JI )
                        NAB( JI, 2 ) = IWORK( JI )
                     ELSE
                        INFO = MMAX + 1
                     END IF
                  END IF
   70          CONTINUE
               IF( INFO.NE.0 )
     $            RETURN
               KL = KLNEW
            ELSE
*
*              IJOB=3: BINARY SEARCH.  KEEP ONLY THE INTERVAL CONTAINING
*                      W   S.T. N(W) = NVAL
*
               DO 80 JI = KF, KL
                  IF( IWORK( JI ).LE.NVAL( JI ) ) THEN
                     AB( JI, 1 ) = C( JI )
                     NAB( JI, 1 ) = IWORK( JI )
                  END IF
                  IF( IWORK( JI ).GE.NVAL( JI ) ) THEN
                     AB( JI, 2 ) = C( JI )
                     NAB( JI, 2 ) = IWORK( JI )
                  END IF
   80          CONTINUE
            END IF
*
         ELSE
*
*           END OF PARALLEL VERSION OF THE LOOP
*
*           BEGIN OF SERIAL VERSION OF THE LOOP
*
            KLNEW = KL
            DO 100 JI = KF, KL
*
*              COMPUTE N(W), THE NUMBER OF EIGENVALUES LESS THAN W
*
               TMP1 = C( JI )
               TMP2 = D( 1 ) - TMP1
               ITMP1 = 0
               IF( TMP2.LE.PIVMIN ) THEN
                  ITMP1 = 1
                  TMP2 = MIN( TMP2, -PIVMIN )
               END IF
*
*              A SERIES OF COMPILER DIRECTIVES TO DEFEAT VECTORIZATION
*              FOR THE NEXT LOOP
*
*$PL$ CMCHAR=' '
CDIR$          NEXTSCALAR
C$DIR          SCALAR
CDIR$          NEXT SCALAR
CVD$L          NOVECTOR
CDEC$          NOVECTOR
CVD$           NOVECTOR
*VDIR          NOVECTOR
*VOCL          LOOP,SCALAR
CIBM           PREFER SCALAR
*$PL$ CMCHAR='*'
*
               DO 90 J = 2, N
                  TMP2 = D( J ) - E2( J-1 ) / TMP2 - TMP1
                  IF( TMP2.LE.PIVMIN ) THEN
                     ITMP1 = ITMP1 + 1
                     TMP2 = MIN( TMP2, -PIVMIN )
                  END IF
   90          CONTINUE
*
               IF( IJOB.LE.2 ) THEN
*
*                 IJOB=2: CHOOSE ALL INTERVALS CONTAINING EIGENVALUES.
*
*                 INSURE THAT N(W) IS MONOTONE
*
                  ITMP1 = MIN( NAB( JI, 2 ),
     $                    MAX( NAB( JI, 1 ), ITMP1 ) )
*
*                 UPDATE THE QUEUE -- ADD INTERVALS IF BOTH HALVES
*                 CONTAIN EIGENVALUES.
*
                  IF( ITMP1.EQ.NAB( JI, 2 ) ) THEN
*
*                    NO EIGENVALUE IN THE UPPER INTERVAL:
*                    JUST USE THE LOWER INTERVAL.
*
                     AB( JI, 2 ) = TMP1
*
                  ELSE IF( ITMP1.EQ.NAB( JI, 1 ) ) THEN
*
*                    NO EIGENVALUE IN THE LOWER INTERVAL:
*                    JUST USE THE UPPER INTERVAL.
*
                     AB( JI, 1 ) = TMP1
                  ELSE IF( KLNEW.LT.MMAX ) THEN
*
*                    EIGENVALUE IN BOTH INTERVALS -- ADD UPPER TO QUEUE.
*
                     KLNEW = KLNEW + 1
                     AB( KLNEW, 2 ) = AB( JI, 2 )
                     NAB( KLNEW, 2 ) = NAB( JI, 2 )
                     AB( KLNEW, 1 ) = TMP1
                     NAB( KLNEW, 1 ) = ITMP1
                     AB( JI, 2 ) = TMP1
                     NAB( JI, 2 ) = ITMP1
                  ELSE
                     INFO = MMAX + 1
                     RETURN
                  END IF
               ELSE
*
*                 IJOB=3: BINARY SEARCH.  KEEP ONLY THE INTERVAL
*                         CONTAINING  W  S.T. N(W) = NVAL
*
                  IF( ITMP1.LE.NVAL( JI ) ) THEN
                     AB( JI, 1 ) = TMP1
                     NAB( JI, 1 ) = ITMP1
                  END IF
                  IF( ITMP1.GE.NVAL( JI ) ) THEN
                     AB( JI, 2 ) = TMP1
                     NAB( JI, 2 ) = ITMP1
                  END IF
               END IF
  100       CONTINUE
            KL = KLNEW
*
*           END OF SERIAL VERSION OF THE LOOP
*
         END IF
*
*        CHECK FOR CONVERGENCE
*
         KFNEW = KF
         DO 110 JI = KF, KL
            TMP1 = ABS( AB( JI, 2 )-AB( JI, 1 ) )
            TMP2 = MAX( ABS( AB( JI, 2 ) ), ABS( AB( JI, 1 ) ) )
            IF( TMP1.LT.MAX( ABSTOL, PIVMIN, RELTOL*TMP2 ) .OR.
     $          NAB( JI, 1 ).GE.NAB( JI, 2 ) ) THEN
*
*              CONVERGED -- SWAP WITH POSITION KFNEW,
*                           THEN INCREMENT KFNEW
*
               IF( JI.GT.KFNEW ) THEN
                  TMP1 = AB( JI, 1 )
                  TMP2 = AB( JI, 2 )
                  ITMP1 = NAB( JI, 1 )
                  ITMP2 = NAB( JI, 2 )
                  AB( JI, 1 ) = AB( KFNEW, 1 )
                  AB( JI, 2 ) = AB( KFNEW, 2 )
                  NAB( JI, 1 ) = NAB( KFNEW, 1 )
                  NAB( JI, 2 ) = NAB( KFNEW, 2 )
                  AB( KFNEW, 1 ) = TMP1
                  AB( KFNEW, 2 ) = TMP2
                  NAB( KFNEW, 1 ) = ITMP1
                  NAB( KFNEW, 2 ) = ITMP2
                  IF( IJOB.EQ.3 ) THEN
                     ITMP1 = NVAL( JI )
                     NVAL( JI ) = NVAL( KFNEW )
                     NVAL( KFNEW ) = ITMP1
                  END IF
               END IF
               KFNEW = KFNEW + 1
            END IF
  110    CONTINUE
         KF = KFNEW
*
*        CHOOSE MIDPOINTS
*
         DO 120 JI = KF, KL
            C( JI ) = HALF*( AB( JI, 1 )+AB( JI, 2 ) )
  120    CONTINUE
*
*        IF NO MORE INTERVALS TO REFINE, QUIT.
*
         IF( KF.GT.KL )
     $      GO TO 140
  130 CONTINUE
*
*     CONVERGED
*
  140 CONTINUE
      INFO = MAX( KL+1-KF, 0 )
      MOUT = KL
*
      RETURN
*
*     END OF DLAEBZ
*
      END
CUT HERE............
CAT > DSTEQR.F <<'CUT HERE............'
      SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DSTEQR COMPUTES ALL EIGENVALUES AND, OPTIONALLY, EIGENVECTORS OF A
*  SYMMETRIC TRIDIAGONAL MATRIX USING THE IMPLICIT QL OR QR METHOD.
*  THE EIGENVECTORS OF A FULL OR BAND SYMMETRIC MATRIX CAN ALSO BE FOUND
*  IF DSYTRD OR DSPTRD OR DSBTRD HAS BEEN USED TO REDUCE THIS MATRIX TO
*  TRIDIAGONAL FORM.
*
*  ARGUMENTS
*  =========
*
*  COMPZ   (INPUT) CHARACTER*1
*          = 'N':  COMPUTE EIGENVALUES ONLY.
*          = 'V':  COMPUTE EIGENVALUES AND EIGENVECTORS OF THE ORIGINAL
*                  SYMMETRIC MATRIX.  ON ENTRY, Z MUST CONTAIN THE
*                  ORTHOGONAL MATRIX USED TO REDUCE THE ORIGINAL MATRIX
*                  TO TRIDIAGONAL FORM.
*          = 'I':  COMPUTE EIGENVALUES AND EIGENVECTORS OF THE
*                  TRIDIAGONAL MATRIX.  Z IS INITIALIZED TO THE IDENTITY
*                  MATRIX.
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX.  N >= 0.
*
*  D       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          ON ENTRY, THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
*          ON EXIT, IF INFO = 0, THE EIGENVALUES IN ASCENDING ORDER.
*
*  E       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          ON ENTRY, THE (N-1) SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
*          MATRIX.
*          ON EXIT, E HAS BEEN DESTROYED.
*
*  Z       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDZ, N)
*          ON ENTRY, IF  COMPZ = 'V', THEN Z CONTAINS THE ORTHOGONAL
*          MATRIX USED IN THE REDUCTION TO TRIDIAGONAL FORM.
*          ON EXIT, IF  COMPZ = 'V', Z CONTAINS THE ORTHONORMAL
*          EIGENVECTORS OF THE ORIGINAL SYMMETRIC MATRIX, AND IF
*          COMPZ = 'I', Z CONTAINS THE ORTHONORMAL EIGENVECTORS OF
*          THE SYMMETRIC TRIDIAGONAL MATRIX.  IF AN ERROR EXIT IS
*          MADE, Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE
*          STORED EIGENVALUES.
*          IF COMPZ = 'N', THEN Z IS NOT REFERENCED.
*
*  LDZ     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY Z.  LDZ >= 1, AND IF
*          EIGENVECTORS ARE DESIRED, THEN  LDZ >= MAX(1,N).
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (MAX(1,2*N-2))
*          IF COMPZ = 'N', THEN WORK IS NOT REFERENCED.
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0:  THE ALGORITHM HAS FAILED TO FIND ALL THE EIGENVALUES IN
*                A TOTAL OF 30*N ITERATIONS; IF INFO = I, THEN I
*                ELEMENTS OF E HAVE NOT CONVERGED TO ZERO; ON EXIT, D
*                AND E CONTAIN THE ELEMENTS OF A SYMMETRIC TRIDIAGONAL
*                MATRIX WHICH IS ORTHOGONALLY SIMILAR TO THE ORIGINAL
*                MATRIX.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, ICOMPZ, II, J, JTOT, K, L, L1, LEND, LENDM1,
     $                   LENDP1, LM1, M, MM, MM1, NM1, NMAXIT
      DOUBLE PRECISION   B, C, EPS, F, G, P, R, RT1, RT2, S, TST
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLAPY2
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASR, DLAZRO, DSWAP,
     $                   XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, MAX, SIGN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEQR', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.GT.0 )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     DETERMINE THE UNIT ROUNDOFF FOR THIS ENVIRONMENT.
*
      EPS = DLAMCH( 'E' )
*
*     COMPUTE THE EIGENVALUES AND EIGENVECTORS OF THE TRIDIAGONAL
*     MATRIX.
*
      IF( ICOMPZ.EQ.2 )
     $   CALL DLAZRO( N, N, ZERO, ONE, Z, LDZ )
*
      NMAXIT = N*MAXIT
      JTOT = 0
*
*     DETERMINE WHERE THE MATRIX SPLITS AND CHOOSE QL OR QR ITERATION
*     FOR EACH BLOCK, ACCORDING TO WHETHER TOP OR BOTTOM DIAGONAL
*     ELEMENT IS SMALLER.
*
      L1 = 1
      NM1 = N - 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.LE.EPS*( ABS( D( M ) )+ABS( D( M+1 ) ) ) )
     $         GO TO 30
   20    CONTINUE
      END IF
      M = N
*
   30 CONTINUE
      L = L1
      LEND = M
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         L = LEND
         LEND = L1
      END IF
      L1 = M + 1
*
      IF( LEND.GE.L ) THEN
*
*        QL ITERATION
*
*        LOOK FOR SMALL SUBDIAGONAL ELEMENT.
*
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )
               IF( TST.LE.EPS*( ABS( D( M ) )+ABS( D( M+1 ) ) ) )
     $            GO TO 60
   50       CONTINUE
         END IF
*
         M = LEND
*
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 80
*
*        IF REMAINING MATRIX IS 2-BY-2, USE DLAE2 OR DLAEV2
*        TO COMPUTE ITS EIGENSYSTEM.
*
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL DLASR( 'R', 'V', 'B', N, 2, WORK( L ),
     $                     WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 10
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        FORM SHIFT.
*
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        INNER LOOP
*
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
*
*           IF EIGENVECTORS ARE DESIRED, THEN SAVE ROTATIONS.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
*
   70    CONTINUE
*
*        IF EIGENVECTORS ARE DESIRED, THEN APPLY SAVED ROTATIONS.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL DLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
*
*        EIGENVALUE FOUND.
*
   80    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 10
*
      ELSE
*
*        QR ITERATION
*
*        LOOK FOR SMALL SUPERDIAGONAL ELEMENT.
*
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )
               IF( TST.LE.EPS*( ABS( D( M ) )+ABS( D( M-1 ) ) ) )
     $            GO TO 110
  100       CONTINUE
         END IF
*
         M = LEND
*
  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 130
*
*        IF REMAINING MATRIX IS 2-BY-2, USE DLAE2 OR DLAEV2
*        TO COMPUTE ITS EIGENSYSTEM.
*
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL DLASR( 'R', 'V', 'F', N, 2, WORK( M ),
     $                     WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 10
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        FORM SHIFT.
*
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        INNER LOOP
*
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
*
*           IF EIGENVECTORS ARE DESIRED, THEN SAVE ROTATIONS.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
*
  120    CONTINUE
*
*        IF EIGENVECTORS ARE DESIRED, THEN APPLY SAVED ROTATIONS.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL DLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, M ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
*
*        EIGENVALUE FOUND.
*
  130    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 10
*
      END IF
*
*     SET ERROR -- NO CONVERGENCE TO AN EIGENVALUE AFTER A TOTAL
*     OF N*MAXIT ITERATIONS.
*
  140 CONTINUE
      DO 150 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  150 CONTINUE
      RETURN
*
*     ORDER EIGENVALUES AND EIGENVECTORS.
*
  160 CONTINUE
      DO 180 II = 2, N
         I = II - 1
         K = I
         P = D( I )
         DO 170 J = II, N
            IF( D( J ).LT.P ) THEN
               K = J
               P = D( J )
            END IF
  170    CONTINUE
         IF( K.NE.I ) THEN
            D( K ) = D( I )
            D( I ) = P
            IF( ICOMPZ.GT.0 )
     $         CALL DSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
         END IF
  180 CONTINUE
*
      RETURN
*
*     END OF DSTEQR
*
      END
CUT HERE............
CAT > DLARTG.F <<'CUT HERE............'
      SUBROUTINE DLARTG( F, G, CS, SN, R )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   CS, F, G, R, SN
*     ..
*
*  PURPOSE
*  =======
*
*  DLARTG GENERATE A PLANE ROTATION SO THAT
*
*     [  CS  SN  ]  .  [ F ]  =  [ R ]   WHERE CS**2 + SN**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  THIS IS A FASTER VERSION OF THE BLAS1 ROUTINE DROTG, EXCEPT FOR
*  THE FOLLOWING DIFFERENCES:
*     F AND G ARE UNCHANGED ON RETURN.
*     IF G=0, THEN CS=1 AND SN=0.
*     IF F=0 AND (G .NE. 0), THEN CS=0 AND SN=1 WITHOUT DOING ANY
*        FLOATING POINT OPERATIONS (SAVES WORK IN DBDSQR WHEN
*        THERE ARE ZEROS ON THE DIAGONAL).
*
*  ARGUMENTS
*  =========
*
*  F       (INPUT) DOUBLE PRECISION
*          THE FIRST COMPONENT OF VECTOR TO BE ROTATED.
*
*  G       (INPUT) DOUBLE PRECISION
*          THE SECOND COMPONENT OF VECTOR TO BE ROTATED.
*
*  CS      (OUTPUT) DOUBLE PRECISION
*          THE COSINE OF THE ROTATION.
*
*  SN      (OUTPUT) DOUBLE PRECISION
*          THE SINE OF THE ROTATION.
*
*  R       (OUTPUT) DOUBLE PRECISION
*          THE NONZERO COMPONENT OF THE ROTATED VECTOR.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. LOCAL SCALARS ..
      DOUBLE PRECISION   T, TT
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      IF( G.EQ.ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F.EQ.ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         IF( ABS( F ).GT.ABS( G ) ) THEN
            T = G / F
            TT = SQRT( ONE+T*T )
            CS = ONE / TT
            SN = T*CS
            R = F*TT
         ELSE
            T = F / G
            TT = SQRT( ONE+T*T )
            SN = ONE / TT
            CS = T*SN
            R = G*TT
         END IF
      END IF
      RETURN
*
*     END OF DLARTG
*
      END
CUT HERE............
CAT > DLASR.F <<'CUT HERE............'
      SUBROUTINE DLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), S( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLASR   PERFORMS THE TRANSFORMATION
*
*     A := P*A,   WHEN SIDE = 'L' OR 'L'  (  LEFT-HAND SIDE )
*
*     A := A*P',  WHEN SIDE = 'R' OR 'R'  ( RIGHT-HAND SIDE )
*
*  WHERE A IS AN M BY N REAL MATRIX AND P IS AN ORTHOGONAL MATRIX,
*  CONSISTING OF A SEQUENCE OF PLANE ROTATIONS DETERMINED BY THE
*  PARAMETERS PIVOT AND DIRECT AS FOLLOWS ( Z = M WHEN SIDE = 'L' OR 'L'
*  AND Z = N WHEN SIDE = 'R' OR 'R' ):
*
*  WHEN  DIRECT = 'F' OR 'F'  ( FORWARD SEQUENCE ) THEN
*
*     P = P( Z - 1 )*...*P( 2 )*P( 1 ),
*
*  AND WHEN DIRECT = 'B' OR 'B'  ( BACKWARD SEQUENCE ) THEN
*
*     P = P( 1 )*P( 2 )*...*P( Z - 1 ),
*
*  WHERE  P( K ) IS A PLANE ROTATION MATRIX FOR THE FOLLOWING PLANES:
*
*     WHEN  PIVOT = 'V' OR 'V'  ( VARIABLE PIVOT ),
*        THE PLANE ( K, K + 1 )
*
*     WHEN  PIVOT = 'T' OR 'T'  ( TOP PIVOT ),
*        THE PLANE ( 1, K + 1 )
*
*     WHEN  PIVOT = 'B' OR 'B'  ( BOTTOM PIVOT ),
*        THE PLANE ( K, Z )
*
*  C( K ) AND S( K )  MUST CONTAIN THE  COSINE AND SINE THAT DEFINE THE
*  MATRIX  P( K ).  THE TWO BY TWO PLANE ROTATION PART OF THE MATRIX
*  P( K ), R( K ), IS ASSUMED TO BE OF THE FORM
*
*     R( K ) = (  C( K )  S( K ) ).
*              ( -S( K )  C( K ) )
*
*  THIS VERSION VECTORISES ACROSS ROWS OF THE ARRAY A WHEN SIDE = 'L'.
*
*  ARGUMENTS
*  =========
*
*  SIDE    (INPUT) CHARACTER*1
*          SPECIFIES WHETHER THE PLANE ROTATION MATRIX P IS APPLIED TO
*          A ON THE LEFT OR THE RIGHT.
*          = 'L':  LEFT, COMPUTE A := P*A
*          = 'R':  RIGHT, COMPUTE A:= A*P'
*
*  DIRECT  (INPUT) CHARACTER*1
*          SPECIFIES WHETHER P IS A FORWARD OR BACKWARD SEQUENCE OF
*          PLANE ROTATIONS.
*          = 'F':  FORWARD, P = P( Z - 1 )*...*P( 2 )*P( 1 )
*          = 'B':  BACKWARD, P = P( 1 )*P( 2 )*...*P( Z - 1 )
*
*  PIVOT   (INPUT) CHARACTER*1
*          SPECIFIES THE PLANE FOR WHICH P(K) IS A PLANE ROTATION
*          MATRIX.
*          = 'V':  VARIABLE PIVOT, THE PLANE (K,K+1)
*          = 'T':  TOP PIVOT, THE PLANE (1,K+1)
*          = 'B':  BOTTOM PIVOT, THE PLANE (K,Z)
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX A.  IF M <= 1, AN IMMEDIATE
*          RETURN IS EFFECTED.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX A.  IF N <= 1, AN
*          IMMEDIATE RETURN IS EFFECTED.
*
*  C, S    (INPUT) DOUBLE PRECISION ARRAYS, DIMENSION
*                  (M-1) IF SIDE = 'L'
*                  (N-1) IF SIDE = 'R'
*          C(K) AND S(K) CONTAIN THE COSINE AND SINE THAT DEFINE THE
*          MATRIX P(K).  THE TWO BY TWO PLANE ROTATION PART OF THE
*          MATRIX P(K), R(K), IS ASSUMED TO BE OF THE FORM
*          R( K ) = (  C( K )  S( K ) ).
*                   ( -S( K )  C( K ) )
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          THE M BY N MATRIX A.  ON EXIT, A IS OVERWRITTEN BY P*A IF
*          SIDE = 'R' OR BY A*P' IF SIDE = 'L'.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,M).
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP, TEMP
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS
*
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASR ', INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        FORM  P * A
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        FORM A * P'
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     END OF DLASR
*
      END
CUT HERE............
CAT > DLAEV2.F <<'CUT HERE............'
      SUBROUTINE DLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
*     ..
*
*  PURPOSE
*  =======
*
*  DLAEV2 COMPUTES THE EIGENDECOMPOSITION OF A 2-BY-2 SYMMETRIC MATRIX
*     [  A   B  ]
*     [  B   C  ].
*  ON RETURN, RT1 IS THE EIGENVALUE OF LARGER ABSOLUTE VALUE, RT2 IS THE
*  EIGENVALUE OF SMALLER ABSOLUTE VALUE, AND (CS1,SN1) IS THE UNIT RIGHT
*  EIGENVECTOR FOR RT1, GIVING THE DECOMPOSITION
*
*     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
*     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
*
*  ARGUMENTS
*  =========
*
*  A       (INPUT) DOUBLE PRECISION
*          THE (1,1) ENTRY OF THE 2-BY-2 MATRIX.
*
*  B       (INPUT) DOUBLE PRECISION
*          THE (1,2) ENTRY AND THE CONJUGATE OF THE (2,1) ENTRY OF THE
*          2-BY-2 MATRIX.
*
*  C       (INPUT) DOUBLE PRECISION
*          THE (2,2) ENTRY OF THE 2-BY-2 MATRIX.
*
*  RT1     (OUTPUT) DOUBLE PRECISION
*          THE EIGENVALUE OF LARGER ABSOLUTE VALUE.
*
*  RT2     (OUTPUT) DOUBLE PRECISION
*          THE EIGENVALUE OF SMALLER ABSOLUTE VALUE.
*
*  CS1     (OUTPUT) DOUBLE PRECISION
*  SN1     (OUTPUT) DOUBLE PRECISION
*          THE VECTOR (CS1, SN1) IS A UNIT RIGHT EIGENVECTOR FOR RT1.
*
*  FURTHER DETAILS
*  ===============
*
*  RT1 IS ACCURATE TO A FEW ULPS BARRING OVER/UNDERFLOW.
*
*  RT2 MAY BE INACCURATE IF THERE IS MASSIVE CANCELLATION IN THE
*  DETERMINANT A*C-B*B; HIGHER PRECISION OR CORRECTLY ROUNDED OR
*  CORRECTLY TRUNCATED ARITHMETIC WOULD BE NEEDED TO COMPUTE RT2
*  ACCURATELY IN ALL CASES.
*
*  CS1 AND SN1 ARE ACCURATE TO A FEW ULPS BARRING OVER/UNDERFLOW.
*
*  OVERFLOW IS POSSIBLE ONLY IF RT1 IS WITHIN A FACTOR OF 5 OF OVERFLOW.
*  UNDERFLOW IS HARMLESS IF THE INPUT DATA IS 0 OR EXCEEDS
*     UNDERFLOW_THRESHOLD / MACHEPS.
*
* =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            SGN1, SGN2
      DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM,
     $                   TB, TN
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     COMPUTE THE EIGENVALUES
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        INCLUDES CASE AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
*
*        ORDER OF EXECUTION IMPORTANT.
*        TO GET FULLY ACCURATE SMALLER EIGENVALUE,
*        NEXT LINE NEEDS TO BE EXECUTED IN HIGHER PRECISION.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
*
*        ORDER OF EXECUTION IMPORTANT.
*        TO GET FULLY ACCURATE SMALLER EIGENVALUE,
*        NEXT LINE NEEDS TO BE EXECUTED IN HIGHER PRECISION.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        INCLUDES CASE RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
*
*     COMPUTE THE EIGENVECTOR
*
      IF( DF.GE.ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS.GT.AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB.EQ.ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1.EQ.SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
*
*     END OF DLAEV2
*
      END
CUT HERE............
CAT > DLAZRO.F <<'CUT HERE............'
      SUBROUTINE DLAZRO( M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            LDA, M, N
      DOUBLE PRECISION   ALPHA, BETA
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLAZRO INITIALIZES A 2-D ARRAY A TO BETA ON THE DIAGONAL AND
*  ALPHA ON THE OFFDIAGONALS.
*
*  ARGUMENTS
*  =========
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX A.  M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX A.  N >= 0.
*
*  ALPHA   (INPUT) DOUBLE PRECISION
*          THE CONSTANT TO WHICH THE OFFDIAGONAL ELEMENTS ARE TO BE SET.
*
*  BETA    (INPUT) DOUBLE PRECISION
*          THE CONSTANT TO WHICH THE DIAGONAL ELEMENTS ARE TO BE SET.
*
*  A       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON EXIT, THE LEADING M BY N SUBMATRIX OF A IS SET SUCH THAT
*             A(I,J) = ALPHA,  1 <= I <= M, 1 <= J <= N, I <> J
*             A(I,I) = BETA,   1 <= I <= MIN(M,N).
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,M).
*
*  =====================================================================
*
*     .. LOCAL SCALARS ..
      INTEGER            I, J
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      DO 20 J = 1, N
         DO 10 I = 1, M
            A( I, J ) = ALPHA
   10    CONTINUE
   20 CONTINUE
*
      DO 30 I = 1, MIN( M, N )
         A( I, I ) = BETA
   30 CONTINUE
*
      RETURN
*
*     END OF DLAZRO
*
      END
CUT HERE............
CAT > DORGTR.F <<'CUT HERE............'
      SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  PURPOSE
*  =======
*
*  DORGTR GENERATES A REAL ORTHOGONAL MATRIX Q WHICH IS DEFINED AS THE
*  PRODUCT OF N-1 ELEMENTARY REFLECTORS OF ORDER N, AS RETURNED BY
*  DSYTRD:
*
*  IF UPLO = 'U', Q = H(N-1) . . . H(2) H(1),
*
*  IF UPLO = 'L', Q = H(1) H(2) . . . H(N-1).
*
*  ARGUMENTS
*  =========
*
*  UPLO    (INPUT) CHARACTER*1
*          = 'U': UPPER TRIANGLE OF A CONTAINS ELEMENTARY REFLECTORS
*                 FROM DSYTRD;
*          = 'L': LOWER TRIANGLE OF A CONTAINS ELEMENTARY REFLECTORS
*                 FROM DSYTRD.
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX Q. N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE VECTORS WHICH DEFINE THE ELEMENTARY REFLECTORS,
*          AS RETURNED BY DSYTRD.
*          ON EXIT, THE N-BY-N ORTHOGONAL MATRIX Q.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A. LDA >= MAX(1,N).
*
*  TAU     (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          TAU(I) MUST CONTAIN THE SCALAR FACTOR OF THE ELEMENTARY
*          REFLECTOR H(I), AS RETURNED BY DSYTRD.
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (LWORK)
*          ON EXIT, IF INFO = 0, WORK(1) RETURNS THE OPTIMAL LWORK.
*
*  LWORK   (INPUT) INTEGER
*          THE DIMENSION OF THE ARRAY WORK. LWORK >= MAX(1,N-1).
*          FOR OPTIMUM PERFORMANCE LWORK >= (N-1)*NB, WHERE NB IS
*          THE OPTIMAL BLOCKSIZE.
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            UPPER
      INTEGER            I, IINFO, J
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DORGQL, DORGQR, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT ARGUMENTS
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N-1 ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGTR', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Q WAS DETERMINED BY A CALL TO DSYTRD WITH UPLO = 'U'
*
*        SHIFT THE VECTORS WHICH DEFINE THE ELEMENTARY REFLECTORS ONE
*        COLUMN TO THE LEFT, AND SET THE LAST ROW AND COLUMN OF Q TO
*        THOSE OF THE UNIT MATRIX
*
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               A( I, J ) = A( I, J+1 )
   10       CONTINUE
            A( N, J ) = ZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            A( I, N ) = ZERO
   30    CONTINUE
         A( N, N ) = ONE
*
*        GENERATE Q(1:N-1,1:N-1)
*
         CALL DORGQL( N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
*
      ELSE
*
*        Q WAS DETERMINED BY A CALL TO DSYTRD WITH UPLO = 'L'.
*
*        SHIFT THE VECTORS WHICH DEFINE THE ELEMENTARY REFLECTORS ONE
*        COLUMN TO THE RIGHT, AND SET THE FIRST ROW AND COLUMN OF Q TO
*        THOSE OF THE UNIT MATRIX
*
         DO 50 J = N, 2, -1
            A( 1, J ) = ZERO
            DO 40 I = J + 1, N
               A( I, J ) = A( I, J-1 )
   40       CONTINUE
   50    CONTINUE
         A( 1, 1 ) = ONE
         DO 60 I = 2, N
            A( I, 1 ) = ZERO
   60    CONTINUE
         IF( N.GT.1 ) THEN
*
*           GENERATE Q(2:N,2:N)
*
            CALL DORGQR( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK,
     $                   LWORK, IINFO )
         END IF
      END IF
      RETURN
*
*     END OF DORGTR
*
      END
CUT HERE............
CAT > DORGQR.F <<'CUT HERE............'
      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  PURPOSE
*  =======
*
*  DORGQR GENERATES AN M-BY-N REAL MATRIX Q WITH ORTHONORMAL COLUMNS,
*  WHICH IS DEFINED AS THE FIRST N COLUMNS OF A PRODUCT OF K ELEMENTARY
*  REFLECTORS OF ORDER M
*
*        Q  =  H(1) H(2) . . . H(K)
*
*  AS RETURNED BY DGEQRF.
*
*  ARGUMENTS
*  =========
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX Q. M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX Q. M >= N >= 0.
*
*  K       (INPUT) INTEGER
*          THE NUMBER OF ELEMENTARY REFLECTORS WHOSE PRODUCT DEFINES THE
*          MATRIX Q. N >= K >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE I-TH COLUMN MUST CONTAIN THE VECTOR WHICH
*          DEFINES THE ELEMENTARY REFLECTOR H(I), FOR I = 1,2,...,K, AS
*          RETURNED BY DGEQRF IN THE FIRST K COLUMNS OF ITS ARRAY
*          ARGUMENT A.
*          ON EXIT, THE M-BY-N MATRIX Q.
*
*  LDA     (INPUT) INTEGER
*          THE FIRST DIMENSION OF THE ARRAY A. LDA >= MAX(1,M).
*
*  TAU     (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (K)
*          TAU(I) MUST CONTAIN THE SCALAR FACTOR OF THE ELEMENTARY
*          REFLECTOR H(I), AS RETURNED BY DGEQRF.
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (LWORK)
*          ON EXIT, IF INFO = 0, WORK(1) RETURNS THE OPTIMAL LWORK.
*
*  LWORK   (INPUT) INTEGER
*          THE DIMENSION OF THE ARRAY WORK. LWORK >= MAX(1,N).
*          FOR OPTIMUM PERFORMANCE LWORK >= N*NB, WHERE NB IS THE
*          OPTIMAL BLOCKSIZE.
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAS AN ILLEGAL VALUE
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, NB,
     $                   NBMIN, NX
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN
*     ..
*     .. EXTERNAL FUNCTIONS ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT ARGUMENTS
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     DETERMINE THE BLOCK SIZE.
*
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        DETERMINE WHEN TO CROSS OVER FROM BLOCKED TO UNBLOCKED CODE.
*
         NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           DETERMINE IF WORKSPACE IS LARGE ENOUGH FOR BLOCKED CODE.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              NOT ENOUGH WORKSPACE TO USE OPTIMAL NB:  REDUCE NB AND
*              DETERMINE THE MINIMUM VALUE OF NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        USE BLOCKED CODE AFTER THE LAST BLOCK.
*        THE FIRST KK COLUMNS ARE HANDLED BY THE BLOCK METHOD.
*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
*
*        SET A(1:KK,KK+1:N) TO ZERO.
*
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     USE UNBLOCKED CODE FOR THE LAST OR ONLY BLOCK.
*
      IF( KK.LT.N )
     $   CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        USE BLOCKED CODE
*
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
*
*              FORM THE TRIANGULAR FACTOR OF THE BLOCK REFLECTOR
*              H = H(I) H(I+1) . . . H(I+IB-1)
*
               CALL DLARFT( 'FORWARD', 'COLUMNWISE', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              APPLY H TO A(I:M,I+IB:N) FROM THE LEFT
*
               CALL DLARFB( 'LEFT', 'NO TRANSPOSE', 'FORWARD',
     $                      'COLUMNWISE', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
*
*           APPLY H TO ROWS I:M OF CURRENT BLOCK
*
            CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
*
*           SET ROWS 1:I-1 OF CURRENT BLOCK TO ZERO
*
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     END OF DORGQR
*
      END
CUT HERE............
CAT > DORG2R.F <<'CUT HERE............'
      SUBROUTINE DORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 29, 1992
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DORG2R GENERATES AN M BY N REAL MATRIX Q WITH ORTHONORMAL COLUMNS,
*  WHICH IS DEFINED AS THE FIRST N COLUMNS OF A PRODUCT OF K ELEMENTARY
*  REFLECTORS OF ORDER M
*
*        Q  =  H(1) H(2) . . . H(K)
*
*  AS RETURNED BY DGEQRF.
*
*  ARGUMENTS
*  =========
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX Q. M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX Q. M >= N >= 0.
*
*  K       (INPUT) INTEGER
*          THE NUMBER OF ELEMENTARY REFLECTORS WHOSE PRODUCT DEFINES THE
*          MATRIX Q. N >= K >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE I-TH COLUMN MUST CONTAIN THE VECTOR WHICH
*          DEFINES THE ELEMENTARY REFLECTOR H(I), FOR I = 1,2,...,K, AS
*          RETURNED BY DGEQRF IN THE FIRST K COLUMNS OF ITS ARRAY
*          ARGUMENT A.
*          ON EXIT, THE M-BY-N MATRIX Q.
*
*  LDA     (INPUT) INTEGER
*          THE FIRST DIMENSION OF THE ARRAY A. LDA >= MAX(1,M).
*
*  TAU     (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (K)
*          TAU(I) MUST CONTAIN THE SCALAR FACTOR OF THE ELEMENTARY
*          REFLECTOR H(I), AS RETURNED BY DGEQRF.
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (N)
*
*  INFO    (OUTPUT) INTEGER
*          = 0: SUCCESSFUL EXIT
*          < 0: IF INFO = -I, THE I-TH ARGUMENT HAS AN ILLEGAL VALUE
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, J, L
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLARF, DSCAL, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT ARGUMENTS
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2R', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.LE.0 )
     $   RETURN
*
*     INITIALISE COLUMNS K+1:N TO COLUMNS OF THE UNIT MATRIX
*
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = K, 1, -1
*
*        APPLY H(I) TO A(I:M,I:N) FROM THE LEFT
*
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL DLARF( 'LEFT', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M )
     $      CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
*
*        SET A(1:I-1,I) TO ZERO
*
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     END OF DORG2R
*
      END
CUT HERE............
CAT > DORGQL.F <<'CUT HERE............'
      SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( LWORK )
*     ..
*
*  PURPOSE
*  =======
*
*  DORGQL GENERATES AN M-BY-N REAL MATRIX Q WITH ORTHONORMAL COLUMNS,
*  WHICH IS DEFINED AS THE LAST N COLUMNS OF A PRODUCT OF K ELEMENTARY
*  REFLECTORS OF ORDER M
*
*        Q  =  H(K) . . . H(2) H(1)
*
*  AS RETURNED BY DGEQLF.
*
*  ARGUMENTS
*  =========
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX Q. M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX Q. M >= N >= 0.
*
*  K       (INPUT) INTEGER
*          THE NUMBER OF ELEMENTARY REFLECTORS WHOSE PRODUCT DEFINES THE
*          MATRIX Q. N >= K >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE (N-K+I)-TH COLUMN MUST CONTAIN THE VECTOR WHICH
*          DEFINES THE ELEMENTARY REFLECTOR H(I), FOR I = 1,2,...,K, AS
*          RETURNED BY DGEQLF IN THE LAST K COLUMNS OF ITS ARRAY
*          ARGUMENT A.
*          ON EXIT, THE M-BY-N MATRIX Q.
*
*  LDA     (INPUT) INTEGER
*          THE FIRST DIMENSION OF THE ARRAY A. LDA >= MAX(1,M).
*
*  TAU     (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (K)
*          TAU(I) MUST CONTAIN THE SCALAR FACTOR OF THE ELEMENTARY
*          REFLECTOR H(I), AS RETURNED BY DGEQLF.
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (LWORK)
*          ON EXIT, IF INFO = 0, WORK(1) RETURNS THE OPTIMAL LWORK.
*
*  LWORK   (INPUT) INTEGER
*          THE DIMENSION OF THE ARRAY WORK. LWORK >= MAX(1,N).
*          FOR OPTIMUM PERFORMANCE LWORK >= N*NB, WHERE NB IS THE
*          OPTIMAL BLOCKSIZE.
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAS AN ILLEGAL VALUE
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, NB, NBMIN,
     $                   NX
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLARFB, DLARFT, DORG2L, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN
*     ..
*     .. EXTERNAL FUNCTIONS ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT ARGUMENTS
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQL', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     DETERMINE THE BLOCK SIZE.
*
      NB = ILAENV( 1, 'DORGQL', ' ', M, N, K, -1 )
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        DETERMINE WHEN TO CROSS OVER FROM BLOCKED TO UNBLOCKED CODE.
*
         NX = MAX( 0, ILAENV( 3, 'DORGQL', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           DETERMINE IF WORKSPACE IS LARGE ENOUGH FOR BLOCKED CODE.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              NOT ENOUGH WORKSPACE TO USE OPTIMAL NB:  REDUCE NB AND
*              DETERMINE THE MINIMUM VALUE OF NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQL', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        USE BLOCKED CODE AFTER THE FIRST BLOCK.
*        THE LAST KK COLUMNS ARE HANDLED BY THE BLOCK METHOD.
*
         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
*
*        SET A(M-KK+1:M,1:N-KK) TO ZERO.
*
         DO 20 J = 1, N - KK
            DO 10 I = M - KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     USE UNBLOCKED CODE FOR THE FIRST OR ONLY BLOCK.
*
      CALL DORG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        USE BLOCKED CODE
*
         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            IF( N-K+I.GT.1 ) THEN
*
*              FORM THE TRIANGULAR FACTOR OF THE BLOCK REFLECTOR
*              H = H(I+IB-1) . . . H(I+1) H(I)
*
               CALL DLARFT( 'BACKWARD', 'COLUMNWISE', M-K+I+IB-1, IB,
     $                      A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
*
*              APPLY H TO A(1:M-K+I+IB-1,1:N-K+I-1) FROM THE LEFT
*
               CALL DLARFB( 'LEFT', 'NO TRANSPOSE', 'BACKWARD',
     $                      'COLUMNWISE', M-K+I+IB-1, N-K+I-1, IB,
     $                      A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA,
     $                      WORK( IB+1 ), LDWORK )
            END IF
*
*           APPLY H TO ROWS 1:M-K+I+IB-1 OF CURRENT BLOCK
*
            CALL DORG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA,
     $                   TAU( I ), WORK, IINFO )
*
*           SET ROWS M-K+I+IB:M OF CURRENT BLOCK TO ZERO
*
            DO 40 J = N - K + I, N - K + I + IB - 1
               DO 30 L = M - K + I + IB, M
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     END OF DORGQL
*
      END
CUT HERE............
CAT > DLARFB.F <<'CUT HERE............'
      SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 29, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLARFB APPLIES A REAL BLOCK REFLECTOR H OR ITS TRANSPOSE H' TO A
*  REAL M BY N MATRIX C, FROM EITHER THE LEFT OR THE RIGHT.
*
*  ARGUMENTS
*  =========
*
*  SIDE    (INPUT) CHARACTER*1
*          = 'L': APPLY H OR H' FROM THE LEFT
*          = 'R': APPLY H OR H' FROM THE RIGHT
*
*  TRANS   (INPUT) CHARACTER*1
*          = 'N': APPLY H (NO TRANSPOSE)
*          = 'T': APPLY H' (TRANSPOSE)
*
*  DIRECT  (INPUT) CHARACTER*1
*          INDICATES HOW H IS FORMED FROM A PRODUCT OF ELEMENTARY
*          REFLECTORS
*          = 'F': H = H(1) H(2) . . . H(K) (FORWARD)
*          = 'B': H = H(K) . . . H(2) H(1) (BACKWARD)
*
*  STOREV  (INPUT) CHARACTER*1
*          INDICATES HOW THE VECTORS WHICH DEFINE THE ELEMENTARY
*          REFLECTORS ARE STORED:
*          = 'C': COLUMNWISE
*          = 'R': ROWWISE
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX C.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX C.
*
*  K       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX T (= THE NUMBER OF ELEMENTARY
*          REFLECTORS WHOSE PRODUCT DEFINES THE BLOCK REFLECTOR).
*
*  V       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION
*                                (LDV,K) IF STOREV = 'C'
*                                (LDV,M) IF STOREV = 'R' AND SIDE = 'L'
*                                (LDV,N) IF STOREV = 'R' AND SIDE = 'R'
*          THE MATRIX V. SEE FURTHER DETAILS.
*
*  LDV     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY V.
*          IF STOREV = 'C' AND SIDE = 'L', LDV >= MAX(1,M);
*          IF STOREV = 'C' AND SIDE = 'R', LDV >= MAX(1,N);
*          IF STOREV = 'R', LDV >= K.
*
*  T       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDT,K)
*          THE TRIANGULAR K BY K MATRIX T IN THE REPRESENTATION OF THE
*          BLOCK REFLECTOR.
*
*  LDT     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY T. LDT >= K.
*
*  C       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDC,N)
*          ON ENTRY, THE M BY N MATRIX C.
*          ON EXIT, C IS OVERWRITTEN BY H*C OR H'*C OR C*H OR C*H'.
*
*  LDC     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY C. LDA >= MAX(1,M).
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (LDWORK,K)
*
*  LDWORK  (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY WORK.
*          IF SIDE = 'L', LDWORK >= MAX(1,N);
*          IF SIDE = 'R', LDWORK >= MAX(1,M).
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      CHARACTER          TRANST
      INTEGER            I, J
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DCOPY, DGEMM, DTRMM
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     QUICK RETURN IF POSSIBLE
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
*
      IF( LSAME( STOREV, 'C' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           LET  V =  ( V1 )    (FIRST K ROWS)
*                     ( V2 )
*           WHERE  V1  IS UNIT LOWER TRIANGULAR.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              FORM  H * C  OR  H' * C  WHERE  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (STORED IN WORK)
*
*              W := C1'
*
               DO 10 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
   10          CONTINUE
*
*              W := W * V1
*
               CALL DTRMM( 'RIGHT', 'LOWER', 'NO TRANSPOSE', 'UNIT', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2'*V2
*
                  CALL DGEMM( 'TRANSPOSE', 'NO TRANSPOSE', N, K, M-K,
     $                        ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  OR  W * T
*
               CALL DTRMM( 'RIGHT', 'UPPER', TRANST, 'NON-UNIT', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2 * W'
*
                  CALL DGEMM( 'NO TRANSPOSE', 'TRANSPOSE', M-K, N, K,
     $                        -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL DTRMM( 'RIGHT', 'LOWER', 'TRANSPOSE', 'UNIT', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
   20             CONTINUE
   30          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              FORM  C * H  OR  C * H'  WHERE  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (STORED IN WORK)
*
*              W := C1
*
               DO 40 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
*
*              W := W * V1
*
               CALL DTRMM( 'RIGHT', 'LOWER', 'NO TRANSPOSE', 'UNIT', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2
*
                  CALL DGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  OR  W * T'
*
               CALL DTRMM( 'RIGHT', 'UPPER', TRANS, 'NON-UNIT', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2'
*
                  CALL DGEMM( 'NO TRANSPOSE', 'TRANSPOSE', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL DTRMM( 'RIGHT', 'LOWER', 'TRANSPOSE', 'UNIT', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
*
         ELSE
*
*           LET  V =  ( V1 )
*                     ( V2 )    (LAST K ROWS)
*           WHERE  V2  IS UNIT UPPER TRIANGULAR.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              FORM  H * C  OR  H' * C  WHERE  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (STORED IN WORK)
*
*              W := C2'
*
               DO 70 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
   70          CONTINUE
*
*              W := W * V2
*
               CALL DTRMM( 'RIGHT', 'UPPER', 'NO TRANSPOSE', 'UNIT', N,
     $                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1'*V1
*
                  CALL DGEMM( 'TRANSPOSE', 'NO TRANSPOSE', N, K, M-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  OR  W * T
*
               CALL DTRMM( 'RIGHT', 'LOWER', TRANST, 'NON-UNIT', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1 * W'
*
                  CALL DGEMM( 'NO TRANSPOSE', 'TRANSPOSE', M-K, N, K,
     $                        -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL DTRMM( 'RIGHT', 'UPPER', 'TRANSPOSE', 'UNIT', N, K,
     $                     ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W'
*
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             CONTINUE
   90          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              FORM  C * H  OR  C * H'  WHERE  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (STORED IN WORK)
*
*              W := C2
*
               DO 100 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
*
*              W := W * V2
*
               CALL DTRMM( 'RIGHT', 'UPPER', 'NO TRANSPOSE', 'UNIT', M,
     $                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1
*
                  CALL DGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  OR  W * T'
*
               CALL DTRMM( 'RIGHT', 'LOWER', TRANS, 'NON-UNIT', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1'
*
                  CALL DGEMM( 'NO TRANSPOSE', 'TRANSPOSE', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL DTRMM( 'RIGHT', 'UPPER', 'TRANSPOSE', 'UNIT', M, K,
     $                     ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W
*
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
*
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           LET  V =  ( V1  V2 )    (V1: FIRST K COLUMNS)
*           WHERE  V1  IS UNIT UPPER TRIANGULAR.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              FORM  H * C  OR  H' * C  WHERE  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (STORED IN WORK)
*
*              W := C1'
*
               DO 130 J = 1, K
                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          CONTINUE
*
*              W := W * V1'
*
               CALL DTRMM( 'RIGHT', 'UPPER', 'TRANSPOSE', 'UNIT', N, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2'*V2'
*
                  CALL DGEMM( 'TRANSPOSE', 'TRANSPOSE', N, K, M-K, ONE,
     $                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,
     $                        WORK, LDWORK )
               END IF
*
*              W := W * T'  OR  W * T
*
               CALL DTRMM( 'RIGHT', 'UPPER', TRANST, 'NON-UNIT', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2' * W'
*
                  CALL DGEMM( 'TRANSPOSE', 'TRANSPOSE', M-K, N, K, -ONE,
     $                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL DTRMM( 'RIGHT', 'UPPER', 'NO TRANSPOSE', 'UNIT', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - WORK( I, J )
  140             CONTINUE
  150          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              FORM  C * H  OR  C * H'  WHERE  C = ( C1  C2 )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (STORED IN WORK)
*
*              W := C1
*
               DO 160 J = 1, K
                  CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
*
*              W := W * V1'
*
               CALL DTRMM( 'RIGHT', 'UPPER', 'TRANSPOSE', 'UNIT', M, K,
     $                     ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2'
*
                  CALL DGEMM( 'NO TRANSPOSE', 'TRANSPOSE', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  OR  W * T'
*
               CALL DTRMM( 'RIGHT', 'UPPER', TRANS, 'NON-UNIT', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2
*
                  CALL DGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL DTRMM( 'RIGHT', 'UPPER', 'NO TRANSPOSE', 'UNIT', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
*
            END IF
*
         ELSE
*
*           LET  V =  ( V1  V2 )    (V2: LAST K COLUMNS)
*           WHERE  V2  IS UNIT LOWER TRIANGULAR.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              FORM  H * C  OR  H' * C  WHERE  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (STORED IN WORK)
*
*              W := C2'
*
               DO 190 J = 1, K
                  CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          CONTINUE
*
*              W := W * V2'
*
               CALL DTRMM( 'RIGHT', 'LOWER', 'TRANSPOSE', 'UNIT', N, K,
     $                     ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1'*V1'
*
                  CALL DGEMM( 'TRANSPOSE', 'TRANSPOSE', N, K, M-K, ONE,
     $                        C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  OR  W * T
*
               CALL DTRMM( 'RIGHT', 'LOWER', TRANST, 'NON-UNIT', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1' * W'
*
                  CALL DGEMM( 'TRANSPOSE', 'TRANSPOSE', M-K, N, K, -ONE,
     $                        V, LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL DTRMM( 'RIGHT', 'LOWER', 'NO TRANSPOSE', 'UNIT', N,
     $                     K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W'
*
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
  200             CONTINUE
  210          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              FORM  C * H  OR  C * H'  WHERE  C = ( C1  C2 )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (STORED IN WORK)
*
*              W := C2
*
               DO 220 J = 1, K
                  CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
*
*              W := W * V2'
*
               CALL DTRMM( 'RIGHT', 'LOWER', 'TRANSPOSE', 'UNIT', M, K,
     $                     ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1'
*
                  CALL DGEMM( 'NO TRANSPOSE', 'TRANSPOSE', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  OR  W * T'
*
               CALL DTRMM( 'RIGHT', 'LOWER', TRANS, 'NON-UNIT', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1
*
                  CALL DGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL DTRMM( 'RIGHT', 'LOWER', 'NO TRANSPOSE', 'UNIT', M,
     $                     K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
*
            END IF
*
         END IF
      END IF
*
      RETURN
*
*     END OF DLARFB
*
      END
CUT HERE............
CAT > DLARFT.F <<'CUT HERE............'
      SUBROUTINE DLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 29, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLARFT FORMS THE TRIANGULAR FACTOR T OF A REAL BLOCK REFLECTOR H
*  OF ORDER N, WHICH IS DEFINED AS A PRODUCT OF K ELEMENTARY REFLECTORS.
*
*  IF DIRECT = 'F', H = H(1) H(2) . . . H(K) AND T IS UPPER TRIANGULAR;
*
*  IF DIRECT = 'B', H = H(K) . . . H(2) H(1) AND T IS LOWER TRIANGULAR.
*
*  IF STOREV = 'C', THE VECTOR WHICH DEFINES THE ELEMENTARY REFLECTOR
*  H(I) IS STORED IN THE I-TH COLUMN OF THE ARRAY V, AND
*
*     H  =  I - V * T * V'
*
*  IF STOREV = 'R', THE VECTOR WHICH DEFINES THE ELEMENTARY REFLECTOR
*  H(I) IS STORED IN THE I-TH ROW OF THE ARRAY V, AND
*
*     H  =  I - V' * T * V
*
*  ARGUMENTS
*  =========
*
*  DIRECT  (INPUT) CHARACTER*1
*          SPECIFIES THE ORDER IN WHICH THE ELEMENTARY REFLECTORS ARE
*          MULTIPLIED TO FORM THE BLOCK REFLECTOR:
*          = 'F': H = H(1) H(2) . . . H(K) (FORWARD)
*          = 'B': H = H(K) . . . H(2) H(1) (BACKWARD)
*
*  STOREV  (INPUT) CHARACTER*1
*          SPECIFIES HOW THE VECTORS WHICH DEFINE THE ELEMENTARY
*          REFLECTORS ARE STORED (SEE ALSO FURTHER DETAILS):
*          = 'C': COLUMNWISE
*          = 'R': ROWWISE
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE BLOCK REFLECTOR H. N >= 0.
*
*  K       (INPUT) INTEGER
*          THE ORDER OF THE TRIANGULAR FACTOR T (= THE NUMBER OF
*          ELEMENTARY REFLECTORS). K >= 1.
*
*  V       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION
*                               (LDV,K) IF STOREV = 'C'
*                               (LDV,N) IF STOREV = 'R'
*          THE MATRIX V. SEE FURTHER DETAILS.
*
*  LDV     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY V.
*          IF STOREV = 'C', LDV >= MAX(1,N); IF STOREV = 'R', LDV >= K.
*
*  TAU     (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (K)
*          TAU(I) MUST CONTAIN THE SCALAR FACTOR OF THE ELEMENTARY
*          REFLECTOR H(I).
*
*  T       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDT,K)
*          THE K BY K TRIANGULAR FACTOR T OF THE BLOCK REFLECTOR.
*          IF DIRECT = 'F', T IS UPPER TRIANGULAR; IF DIRECT = 'B', T IS
*          LOWER TRIANGULAR. THE REST OF THE ARRAY IS NOT USED.
*
*  LDT     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY T. LDT >= K.
*
*  FURTHER DETAILS
*  ===============
*
*  THE SHAPE OF THE MATRIX V AND THE STORAGE OF THE VECTORS WHICH DEFINE
*  THE H(I) IS BEST ILLUSTRATED BY THE FOLLOWING EXAMPLE WITH N = 5 AND
*  K = 3. THE ELEMENTS EQUAL TO 1 ARE NOT STORED; THE CORRESPONDING
*  ARRAY ELEMENTS ARE MODIFIED BUT RESTORED ON EXIT. THE REST OF THE
*  ARRAY IS NOT USED.
*
*  DIRECT = 'F' AND STOREV = 'C':         DIRECT = 'F' AND STOREV = 'R':
*
*               V = (  1       )                 V = (  1 V1 V1 V1 V1 )
*                   ( V1  1    )                     (     1 V2 V2 V2 )
*                   ( V1 V2  1 )                     (        1 V3 V3 )
*                   ( V1 V2 V3 )
*                   ( V1 V2 V3 )
*
*  DIRECT = 'B' AND STOREV = 'C':         DIRECT = 'B' AND STOREV = 'R':
*
*               V = ( V1 V2 V3 )                 V = ( V1 V1  1       )
*                   ( V1 V2 V3 )                     ( V2 V2 V2  1    )
*                   (  1 V2 V3 )                     ( V3 V3 V3 V3  1 )
*                   (     1 V3 )
*                   (        1 )
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, J
      DOUBLE PRECISION   VII
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DGEMV, DTRMV
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( LSAME( DIRECT, 'F' ) ) THEN
         DO 20 I = 1, K
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(I)  =  I
*
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
*
*              GENERAL CASE
*
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
*
*                 T(1:I-1,I) := - TAU(I) * V(I:N,1:I-1)' * V(I:N,I)
*
                  CALL DGEMV( 'TRANSPOSE', N-I+1, I-1, -TAU( I ),
     $                        V( I, 1 ), LDV, V( I, I ), 1, ZERO,
     $                        T( 1, I ), 1 )
               ELSE
*
*                 T(1:I-1,I) := - TAU(I) * V(1:I-1,I:N) * V(I,I:N)'
*
                  CALL DGEMV( 'NO TRANSPOSE', I-1, N-I+1, -TAU( I ),
     $                        V( 1, I ), LDV, V( I, I ), LDV, ZERO,
     $                        T( 1, I ), 1 )
               END IF
               V( I, I ) = VII
*
*              T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)
*
               CALL DTRMV( 'UPPER', 'NO TRANSPOSE', 'NON-UNIT', I-1, T,
     $                     LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
            END IF
   20    CONTINUE
      ELSE
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(I)  =  I
*
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
*
*              GENERAL CASE
*
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
*
*                    T(I+1:K,I) :=
*                            - TAU(I) * V(1:N-K+I,I+1:K)' * V(1:N-K+I,I)
*
                     CALL DGEMV( 'TRANSPOSE', N-K+I, K-I, -TAU( I ),
     $                           V( 1, I+1 ), LDV, V( 1, I ), 1, ZERO,
     $                           T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
*
*                    T(I+1:K,I) :=
*                            - TAU(I) * V(I+1:K,1:N-K+I) * V(I,1:N-K+I)'
*
                     CALL DGEMV( 'NO TRANSPOSE', K-I, N-K+I, -TAU( I ),
     $                           V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO,
     $                           T( I+1, I ), 1 )
                     V( I, N-K+I ) = VII
                  END IF
*
*                 T(I+1:K,I) := T(I+1:K,I+1:K) * T(I+1:K,I)
*
                  CALL DTRMV( 'LOWER', 'NO TRANSPOSE', 'NON-UNIT', K-I,
     $                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
*
*     END OF DLARFT
*
      END
CUT HERE............
CAT > DORG2L.F <<'CUT HERE............'
      SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 29, 1992
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DORG2L GENERATES AN M BY N REAL MATRIX Q WITH ORTHONORMAL COLUMNS,
*  WHICH IS DEFINED AS THE LAST N COLUMNS OF A PRODUCT OF K ELEMENTARY
*  REFLECTORS OF ORDER M
*
*        Q  =  H(K) . . . H(2) H(1)
*
*  AS RETURNED BY DGEQLF.
*
*  ARGUMENTS
*  =========
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX Q. M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX Q. M >= N >= 0.
*
*  K       (INPUT) INTEGER
*          THE NUMBER OF ELEMENTARY REFLECTORS WHOSE PRODUCT DEFINES THE
*          MATRIX Q. N >= K >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE (N-K+I)-TH COLUMN MUST CONTAIN THE VECTOR WHICH
*          DEFINES THE ELEMENTARY REFLECTOR H(I), FOR I = 1,2,...,K, AS
*          RETURNED BY DGEQLF IN THE LAST K COLUMNS OF ITS ARRAY
*          ARGUMENT A.
*          ON EXIT, THE M BY N MATRIX Q.
*
*  LDA     (INPUT) INTEGER
*          THE FIRST DIMENSION OF THE ARRAY A. LDA >= MAX(1,M).
*
*  TAU     (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (K)
*          TAU(I) MUST CONTAIN THE SCALAR FACTOR OF THE ELEMENTARY
*          REFLECTOR H(I), AS RETURNED BY DGEQLF.
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (N)
*
*  INFO    (OUTPUT) INTEGER
*          = 0: SUCCESSFUL EXIT
*          < 0: IF INFO = -I, THE I-TH ARGUMENT HAS AN ILLEGAL VALUE
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, II, J, L
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLARF, DSCAL, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT ARGUMENTS
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORG2L', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.LE.0 )
     $   RETURN
*
*     INITIALISE COLUMNS 1:N-K TO COLUMNS OF THE UNIT MATRIX
*
      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = 1, K
         II = N - K + I
*
*        APPLY H(I) TO A(1:M-K+I,1:N-K+I) FROM THE LEFT
*
         A( M-N+II, II ) = ONE
         CALL DLARF( 'LEFT', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A,
     $               LDA, WORK )
         CALL DSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
         A( M-N+II, II ) = ONE - TAU( I )
*
*        SET A(M-K+I+1:M,N-K+I) TO ZERO
*
         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     END OF DORG2L
*
      END
CUT HERE............
CAT > DLARF.F <<'CUT HERE............'
      SUBROUTINE DLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 29, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      DOUBLE PRECISION   TAU
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLARF APPLIES A REAL ELEMENTARY REFLECTOR H TO A REAL M BY N MATRIX
*  C, FROM EITHER THE LEFT OR THE RIGHT. H IS REPRESENTED IN THE FORM
*
*        H = I - TAU * V * V'
*
*  WHERE TAU IS A REAL SCALAR AND V IS A REAL VECTOR.
*
*  IF TAU = 0, THEN H IS TAKEN TO BE THE UNIT MATRIX.
*
*  ARGUMENTS
*  =========
*
*  SIDE    (INPUT) CHARACTER*1
*          = 'L': FORM  H * C
*          = 'R': FORM  C * H
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX C.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX C.
*
*  V       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION
*                     (1 + (M-1)*ABS(INCV)) IF SIDE = 'L'
*                  OR (1 + (N-1)*ABS(INCV)) IF SIDE = 'R'
*          THE VECTOR V IN THE REPRESENTATION OF H. V IS NOT USED IF
*          TAU = 0.
*
*  INCV    (INPUT) INTEGER
*          THE INCREMENT BETWEEN ELEMENTS OF V. INCV <> 0.
*
*  TAU     (INPUT) DOUBLE PRECISION
*          THE VALUE TAU IN THE REPRESENTATION OF H.
*
*  C       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDC,N)
*          ON ENTRY, THE M BY N MATRIX C.
*          ON EXIT, C IS OVERWRITTEN BY THE MATRIX H * C IF SIDE = 'L',
*          OR C * H IF SIDE = 'R'.
*
*  LDC     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY C. LDC >= MAX(1,M).
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION
*                         (N) IF SIDE = 'L'
*                      OR (M) IF SIDE = 'R'
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DGEMV, DGER
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        FORM  H * C
*
         IF( TAU.NE.ZERO ) THEN
*
*           W := C' * V
*
            CALL DGEMV( 'TRANSPOSE', M, N, ONE, C, LDC, V, INCV, ZERO,
     $                  WORK, 1 )
*
*           C := C - V * W'
*
            CALL DGER( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
*
*        FORM  C * H
*
         IF( TAU.NE.ZERO ) THEN
*
*           W := C * V
*
            CALL DGEMV( 'NO TRANSPOSE', M, N, ONE, C, LDC, V, INCV,
     $                  ZERO, WORK, 1 )
*
*           C := C - W * V'
*
            CALL DGER( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
*
*     END OF DLARF
*
      END
CUT HERE............
CAT > DLACPY.F <<'CUT HERE............'
      SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 29, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLACPY COPIES ALL OR PART OF A TWO-DIMENSIONAL MATRIX A TO ANOTHER
*  MATRIX B.
*
*  ARGUMENTS
*  =========
*
*  UPLO    (INPUT) CHARACTER*1
*          SPECIFIES THE PART OF THE MATRIX A TO BE COPIED TO B.
*          = 'U':      UPPER TRIANGULAR PART
*          = 'L':      LOWER TRIANGULAR PART
*          OTHERWISE:  ALL OF THE MATRIX A
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX A.  M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          THE M BY N MATRIX A.  IF UPLO = 'U', ONLY THE UPPER TRIANGLE
*          OR TRAPEZOID IS ACCESSED; IF UPLO = 'L', ONLY THE LOWER
*          TRIANGLE OR TRAPEZOID IS ACCESSED.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,M).
*
*  B       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDB,N)
*          ON EXIT, B = A IN THE LOCATIONS SPECIFIED BY UPLO.
*
*  LDB     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY B.  LDB >= MAX(1,M).
*
*  =====================================================================
*
*     .. LOCAL SCALARS ..
      INTEGER            I, J
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
*
*     END OF DLACPY
*
      END
CUT HERE............
CAT > DSTERF.F <<'CUT HERE............'
      SUBROUTINE DSTERF( N, D, E, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   D( * ), E( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DSTERF COMPUTES ALL EIGENVALUES OF A SYMMETRIC TRIDIAGONAL MATRIX
*  USING THE PAL-WALKER-KAHAN VARIANT OF THE QL OR QR ALGORITHM.
*
*  ARGUMENTS
*  =========
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX.  N >= 0.
*
*  D       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          ON ENTRY, THE N DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.
*          ON EXIT, IF INFO = 0, THE EIGENVALUES IN ASCENDING ORDER.
*
*  E       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          ON ENTRY, THE (N-1) SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
*          MATRIX.
*          ON EXIT, E HAS BEEN DESTROYED.
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0:  THE ALGORITHM FAILED TO FIND ALL OF THE EIGENVALUES IN
*                A TOTAL OF 30*N ITERATIONS; IF INFO = I, THEN I
*                ELEMENTS OF E HAVE NOT CONVERGED TO ZERO.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, II, J, JTOT, K, L, L1, LEND, LENDM1, LENDP1,
     $                   LM1, M, MM1, NM1, NMAXIT
      DOUBLE PRECISION   ALPHA, BB, C, EPS, GAMMA, OLDC, OLDGAM, P, R,
     $                   RT1, RT2, RTE, S, SIGMA, TST
*     ..
*     .. EXTERNAL FUNCTIONS ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           DLAMCH, DLAPY2
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLAE2, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, SIGN, SQRT
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DSTERF', -INFO )
         RETURN
      END IF
      IF( N.LE.1 )
     $   RETURN
*
*     DETERMINE THE UNIT ROUNDOFF FOR THIS ENVIRONMENT.
*
      EPS = DLAMCH( 'E' )
*
*     COMPUTE THE EIGENVALUES OF THE TRIDIAGONAL MATRIX.
*
      DO 10 I = 1, N - 1
         E( I ) = E( I )**2
   10 CONTINUE
*
      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0
*
*     DETERMINE WHERE THE MATRIX SPLITS AND CHOOSE QL OR QR ITERATION
*     FOR EACH BLOCK, ACCORDING TO WHETHER TOP OR BOTTOM DIAGONAL
*     ELEMENT IS SMALLER.
*
      L1 = 1
      NM1 = N - 1
*
   20 CONTINUE
      IF( L1.GT.N )
     $   GO TO 170
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 30 M = L1, NM1
            TST = SQRT( ABS( E( M ) ) )
            IF( TST.LE.EPS*( ABS( D( M ) )+ABS( D( M+1 ) ) ) )
     $         GO TO 40
   30    CONTINUE
      END IF
      M = N
*
   40 CONTINUE
      L = L1
      LEND = M
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         L = LEND
         LEND = L1
      END IF
      L1 = M + 1
*
      IF( LEND.GE.L ) THEN
*
*        QL ITERATION
*
*        LOOK FOR SMALL SUBDIAGONAL ELEMENT.
*
   50    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 60 M = L, LENDM1
               TST = SQRT( ABS( E( M ) ) )
               IF( TST.LE.EPS*( ABS( D( M ) )+ABS( D( M+1 ) ) ) )
     $            GO TO 70
   60       CONTINUE
         END IF
*
         M = LEND
*
   70    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 90
*
*        IF REMAINING MATRIX IS 2 BY 2, USE DLAE2 TO COMPUTE ITS
*        EIGENVALUES.
*
         IF( M.EQ.L+1 ) THEN
            RTE = SQRT( E( L ) )
            CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 50
            GO TO 20
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        FORM SHIFT.
*
         RTE = SQRT( E( L ) )
         SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        INNER LOOP
*
         MM1 = M - 1
         DO 80 I = MM1, L, -1
            BB = E( I )
            R = P + BB
            IF( I.NE.M-1 )
     $         E( I+1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
   80    CONTINUE
*
         E( L ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 50
*
*        EIGENVALUE FOUND.
*
   90    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 50
         GO TO 20
*
      ELSE
*
*        QR ITERATION
*
*        LOOK FOR SMALL SUPERDIAGONAL ELEMENT.
*
  100    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 110 M = L, LENDP1, -1
               TST = SQRT( ABS( E( M-1 ) ) )
               IF( TST.LE.EPS*( ABS( D( M ) )+ABS( D( M-1 ) ) ) )
     $            GO TO 120
  110       CONTINUE
         END IF
*
         M = LEND
*
  120    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 140
*
*        IF REMAINING MATRIX IS 2 BY 2, USE DLAE2 TO COMPUTE ITS
*        EIGENVALUES.
*
         IF( M.EQ.L-1 ) THEN
            RTE = SQRT( E( L-1 ) )
            CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
            D( L ) = RT1
            D( L-1 ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 100
            GO TO 20
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 150
         JTOT = JTOT + 1
*
*        FORM SHIFT.
*
         RTE = SQRT( E( L-1 ) )
         SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
*
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
*
*        INNER LOOP
*
         LM1 = L - 1
         DO 130 I = M, LM1
            BB = E( I )
            R = P + BB
            IF( I.NE.M )
     $         E( I-1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I+1 )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
  130    CONTINUE
*
         E( LM1 ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 100
*
*        EIGENVALUE FOUND.
*
  140    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 100
         GO TO 20
*
      END IF
*
*     SET ERROR -- NO CONVERGENCE TO AN EIGENVALUE AFTER A TOTAL
*     OF N*MAXIT ITERATIONS.
*
  150 CONTINUE
      DO 160 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  160 CONTINUE
      RETURN
*
*     SORT EIGENVALUES IN INCREASING ORDER.
*
  170 CONTINUE
      DO 190 II = 2, N
         I = II - 1
         K = I
         P = D( I )
         DO 180 J = II, N
            IF( D( J ).LT.P ) THEN
               K = J
               P = D( J )
            END IF
  180    CONTINUE
         IF( K.NE.I ) THEN
            D( K ) = D( I )
            D( I ) = P
         END IF
  190 CONTINUE
*
      RETURN
*
*     END OF DSTERF
*
      END
CUT HERE............
CAT > DLAE2.F <<'CUT HERE............'
      SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   A, B, C, RT1, RT2
*     ..
*
*  PURPOSE
*  =======
*
*  DLAE2  COMPUTES THE EIGENVALUES OF A 2-BY-2 SYMMETRIC MATRIX
*     [  A   B  ]
*     [  B   C  ].
*  ON RETURN, RT1 IS THE EIGENVALUE OF LARGER ABSOLUTE VALUE, AND RT2
*  IS THE EIGENVALUE OF SMALLER ABSOLUTE VALUE.
*
*  ARGUMENTS
*  =========
*
*  A       (INPUT) DOUBLE PRECISION
*          THE (1,1) ENTRY OF THE 2-BY-2 MATRIX.
*
*  B       (INPUT) DOUBLE PRECISION
*          THE (1,2) AND (2,1) ENTRIES OF THE 2-BY-2 MATRIX.
*
*  C       (INPUT) DOUBLE PRECISION
*          THE (2,2) ENTRY OF THE 2-BY-2 MATRIX.
*
*  RT1     (OUTPUT) DOUBLE PRECISION
*          THE EIGENVALUE OF LARGER ABSOLUTE VALUE.
*
*  RT2     (OUTPUT) DOUBLE PRECISION
*          THE EIGENVALUE OF SMALLER ABSOLUTE VALUE.
*
*  FURTHER DETAILS
*  ===============
*
*  RT1 IS ACCURATE TO A FEW ULPS BARRING OVER/UNDERFLOW.
*
*  RT2 MAY BE INACCURATE IF THERE IS MASSIVE CANCELLATION IN THE
*  DETERMINANT A*C-B*B; HIGHER PRECISION OR CORRECTLY ROUNDED OR
*  CORRECTLY TRUNCATED ARITHMETIC WOULD BE NEEDED TO COMPUTE RT2
*  ACCURATELY IN ALL CASES.
*
*  OVERFLOW IS POSSIBLE ONLY IF RT1 IS WITHIN A FACTOR OF 5 OF OVERFLOW.
*  UNDERFLOW IS HARMLESS IF THE INPUT DATA IS 0 OR EXCEEDS
*     UNDERFLOW_THRESHOLD / MACHEPS.
*
* =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D0 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D0 )
*     ..
*     .. LOCAL SCALARS ..
      DOUBLE PRECISION   AB, ACMN, ACMX, ADF, DF, RT, SM, TB
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     COMPUTE THE EIGENVALUES
*
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A ).GT.ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF.GT.AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF.LT.AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
*
*        INCLUDES CASE AB=ADF=0
*
         RT = AB*SQRT( TWO )
      END IF
      IF( SM.LT.ZERO ) THEN
         RT1 = HALF*( SM-RT )
*
*        ORDER OF EXECUTION IMPORTANT.
*        TO GET FULLY ACCURATE SMALLER EIGENVALUE,
*        NEXT LINE NEEDS TO BE EXECUTED IN HIGHER PRECISION.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM.GT.ZERO ) THEN
         RT1 = HALF*( SM+RT )
*
*        ORDER OF EXECUTION IMPORTANT.
*        TO GET FULLY ACCURATE SMALLER EIGENVALUE,
*        NEXT LINE NEEDS TO BE EXECUTED IN HIGHER PRECISION.
*
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
*
*        INCLUDES CASE RT1 = RT2 = 0
*
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
*
*     END OF DLAE2
*
      END
CUT HERE............
CAT > DSYTRD.F <<'CUT HERE............'
      SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ),
     $                   WORK( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DSYTRD REDUCES A REAL SYMMETRIC MATRIX A TO REAL SYMMETRIC
*  TRIDIAGONAL FORM T BY AN ORTHOGONAL SIMILARITY TRANSFORMATION:
*  Q**T * A * Q = T.
*
*  ARGUMENTS
*  =========
*
*  UPLO    (INPUT) CHARACTER*1
*          = 'U':  UPPER TRIANGLE OF A IS STORED;
*          = 'L':  LOWER TRIANGLE OF A IS STORED.
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE SYMMETRIC MATRIX A.  IF UPLO = 'U', THE LEADING
*          N-BY-N UPPER TRIANGULAR PART OF A CONTAINS THE UPPER
*          TRIANGULAR PART OF THE MATRIX A, AND THE STRICTLY LOWER
*          TRIANGULAR PART OF A IS NOT REFERENCED.  IF UPLO = 'L', THE
*          LEADING N-BY-N LOWER TRIANGULAR PART OF A CONTAINS THE LOWER
*          TRIANGULAR PART OF THE MATRIX A, AND THE STRICTLY UPPER
*          TRIANGULAR PART OF A IS NOT REFERENCED.
*          ON EXIT, IF UPLO = 'U', THE DIAGONAL AND FIRST SUPERDIAGONAL
*          OF A ARE OVERWRITTEN BY THE CORRESPONDING ELEMENTS OF THE
*          TRIDIAGONAL MATRIX T, AND THE ELEMENTS ABOVE THE FIRST
*          SUPERDIAGONAL, WITH THE ARRAY TAU, REPRESENT THE ORTHOGONAL
*          MATRIX Q AS A PRODUCT OF ELEMENTARY REFLECTORS; IF UPLO
*          = 'L', THE DIAGONAL AND FIRST SUBDIAGONAL OF A ARE OVER-
*          WRITTEN BY THE CORRESPONDING ELEMENTS OF THE TRIDIAGONAL
*          MATRIX T, AND THE ELEMENTS BELOW THE FIRST SUBDIAGONAL, WITH
*          THE ARRAY TAU, REPRESENT THE ORTHOGONAL MATRIX Q AS A PRODUCT
*          OF ELEMENTARY REFLECTORS. SEE FURTHER DETAILS.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,N).
*
*  D       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX T:
*          D(I) = A(I,I).
*
*  E       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          THE OFF-DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX T:
*          E(I) = A(I,I+1) IF UPLO = 'U', E(I) = A(I+1,I) IF UPLO = 'L'.
*
*  TAU     (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          THE SCALAR FACTORS OF THE ELEMENTARY REFLECTORS (SEE FURTHER
*          DETAILS).
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (LWORK)
*          ON EXIT, IF INFO = 0, WORK(1) RETURNS THE OPTIMAL LWORK.
*
*  LWORK   (INPUT) INTEGER
*          THE DIMENSION OF THE ARRAY WORK.  LWORK >= 1.
*          FOR OPTIMUM PERFORMANCE LWORK >= N*NB, WHERE NB IS THE
*          OPTIMAL BLOCKSIZE.
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*
*  FURTHER DETAILS
*  ===============
*
*  IF UPLO = 'U', THE MATRIX Q IS REPRESENTED AS A PRODUCT OF ELEMENTARY
*  REFLECTORS
*
*     Q = H(N-1) . . . H(2) H(1).
*
*  EACH H(I) HAS THE FORM
*
*     H(I) = I - TAU * V * V'
*
*  WHERE TAU IS A REAL SCALAR, AND V IS A REAL VECTOR WITH
*  V(I+1:N) = 0 AND V(I) = 1; V(1:I-1) IS STORED ON EXIT IN
*  A(1:I-1,I+1), AND TAU IN TAU(I).
*
*  IF UPLO = 'L', THE MATRIX Q IS REPRESENTED AS A PRODUCT OF ELEMENTARY
*  REFLECTORS
*
*     Q = H(1) H(2) . . . H(N-1).
*
*  EACH H(I) HAS THE FORM
*
*     H(I) = I - TAU * V * V'
*
*  WHERE TAU IS A REAL SCALAR, AND V IS A REAL VECTOR WITH
*  V(1:I) = 0 AND V(I+1) = 1; V(I+2:N) IS STORED ON EXIT IN A(I+2:N,I),
*  AND TAU IN TAU(I).
*
*  THE CONTENTS OF A ON EXIT ARE ILLUSTRATED BY THE FOLLOWING EXAMPLES
*  WITH N = 5:
*
*  IF UPLO = 'U':                       IF UPLO = 'L':
*
*    (  D   E   V2  V3  V4 )              (  D                  )
*    (      D   E   V3  V4 )              (  E   D              )
*    (          D   E   V4 )              (  V1  E   D          )
*    (              D   E  )              (  V1  V2  E   D      )
*    (                  D  )              (  V1  V2  V3  E   D  )
*
*  WHERE D AND E DENOTE DIAGONAL AND OFF-DIAGONAL ELEMENTS OF T, AND VI
*  DENOTES AN ELEMENT OF THE VECTOR DEFINING H(I).
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            UPPER
      INTEGER            I, IINFO, IWS, J, KK, LDWORK, NB, NBMIN, NX
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLATRD, DSYR2K, DSYTD2, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRD', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     DETERMINE THE BLOCK SIZE.
*
      NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
*
*        DETERMINE WHEN TO CROSS OVER FROM BLOCKED TO UNBLOCKED CODE
*        (LAST BLOCK IS ALWAYS HANDLED BY UNBLOCKED CODE).
*
         NX = MAX( NB, ILAENV( 3, 'DSYTRD', UPLO, N, -1, -1, -1 ) )
         IF( NX.LT.N ) THEN
*
*           DETERMINE IF WORKSPACE IS LARGE ENOUGH FOR BLOCKED CODE.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              NOT ENOUGH WORKSPACE TO USE OPTIMAL NB:  DETERMINE THE
*              MINIMUM VALUE OF NB, AND REDUCE NB OR FORCE USE OF
*              UNBLOCKED CODE BY SETTING NX = N.
*
               NB = LWORK / LDWORK
               NBMIN = ILAENV( 2, 'DSYTRD', UPLO, N, -1, -1, -1 )
               IF( NB.LT.NBMIN )
     $            NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
*
      IF( UPPER ) THEN
*
*        REDUCE THE UPPER TRIANGLE OF A.
*        COLUMNS 1:KK ARE HANDLED BY THE UNBLOCKED METHOD.
*
         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB
*
*           REDUCE COLUMNS I:I+NB-1 TO TRIDIAGONAL FORM AND FORM THE
*           MATRIX W WHICH IS NEEDED TO UPDATE THE UNREDUCED PART OF
*           THE MATRIX
*
            CALL DLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK,
     $                   LDWORK )
*
*           UPDATE THE UNREDUCED SUBMATRIX A(1:I-1,1:I-1), USING AN
*           UPDATE OF THE FORM:  A := A - V*W' - W*V'
*
            CALL DSYR2K( UPLO, 'NO TRANSPOSE', I-1, NB, -ONE, A( 1, I ),
     $                   LDA, WORK, LDWORK, ONE, A, LDA )
*
*           COPY SUPERDIAGONAL ELEMENTS BACK INTO A, AND DIAGONAL
*           ELEMENTS INTO D
*
            DO 10 J = I, I + NB - 1
               A( J-1, J ) = E( J-1 )
               D( J ) = A( J, J )
   10       CONTINUE
   20    CONTINUE
*
*        USE UNBLOCKED CODE TO REDUCE THE LAST OR ONLY BLOCK
*
         CALL DSYTD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
      ELSE
*
*        REDUCE THE LOWER TRIANGLE OF A
*
         DO 40 I = 1, N - NX, NB
*
*           REDUCE COLUMNS I:I+NB-1 TO TRIDIAGONAL FORM AND FORM THE
*           MATRIX W WHICH IS NEEDED TO UPDATE THE UNREDUCED PART OF
*           THE MATRIX
*
            CALL DLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ),
     $                   TAU( I ), WORK, LDWORK )
*
*           UPDATE THE UNREDUCED SUBMATRIX A(I+IB:N,I+IB:N), USING
*           AN UPDATE OF THE FORM:  A := A - V*W' - W*V'
*
            CALL DSYR2K( UPLO, 'NO TRANSPOSE', N-I-NB+1, NB, -ONE,
     $                   A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE,
     $                   A( I+NB, I+NB ), LDA )
*
*           COPY SUBDIAGONAL ELEMENTS BACK INTO A, AND DIAGONAL
*           ELEMENTS INTO D
*
            DO 30 J = I, I + NB - 1
               A( J+1, J ) = E( J )
               D( J ) = A( J, J )
   30       CONTINUE
   40    CONTINUE
*
*        USE UNBLOCKED CODE TO REDUCE THE LAST OR ONLY BLOCK
*
         CALL DSYTD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ),
     $                TAU( I ), IINFO )
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     END OF DSYTRD
*
      END
CUT HERE............
CAT > DSYTD2.F <<'CUT HERE............'
      SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DSYTD2 REDUCES A REAL SYMMETRIC MATRIX A TO SYMMETRIC TRIDIAGONAL
*  FORM T BY AN ORTHOGONAL SIMILARITY TRANSFORMATION: Q' * A * Q = T.
*
*  ARGUMENTS
*  =========
*
*  UPLO    (INPUT) CHARACTER*1
*          SPECIFIES WHETHER THE UPPER OR LOWER TRIANGULAR PART OF THE
*          SYMMETRIC MATRIX A IS STORED:
*          = 'U':  UPPER TRIANGULAR
*          = 'L':  LOWER TRIANGULAR
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE SYMMETRIC MATRIX A.  IF UPLO = 'U', THE LEADING
*          N-BY-N UPPER TRIANGULAR PART OF A CONTAINS THE UPPER
*          TRIANGULAR PART OF THE MATRIX A, AND THE STRICTLY LOWER
*          TRIANGULAR PART OF A IS NOT REFERENCED.  IF UPLO = 'L', THE
*          LEADING N-BY-N LOWER TRIANGULAR PART OF A CONTAINS THE LOWER
*          TRIANGULAR PART OF THE MATRIX A, AND THE STRICTLY UPPER
*          TRIANGULAR PART OF A IS NOT REFERENCED.
*          ON EXIT, IF UPLO = 'U', THE DIAGONAL AND FIRST SUPERDIAGONAL
*          OF A ARE OVERWRITTEN BY THE CORRESPONDING ELEMENTS OF THE
*          TRIDIAGONAL MATRIX T, AND THE ELEMENTS ABOVE THE FIRST
*          SUPERDIAGONAL, WITH THE ARRAY TAU, REPRESENT THE ORTHOGONAL
*          MATRIX Q AS A PRODUCT OF ELEMENTARY REFLECTORS; IF UPLO
*          = 'L', THE DIAGONAL AND FIRST SUBDIAGONAL OF A ARE OVER-
*          WRITTEN BY THE CORRESPONDING ELEMENTS OF THE TRIDIAGONAL
*          MATRIX T, AND THE ELEMENTS BELOW THE FIRST SUBDIAGONAL, WITH
*          THE ARRAY TAU, REPRESENT THE ORTHOGONAL MATRIX Q AS A PRODUCT
*          OF ELEMENTARY REFLECTORS. SEE FURTHER DETAILS.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,N).
*
*  D       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX T:
*          D(I) = A(I,I).
*
*  E       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          THE OFF-DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX T:
*          E(I) = A(I,I+1) IF UPLO = 'U', E(I) = A(I+1,I) IF UPLO = 'L'.
*
*  TAU     (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          THE SCALAR FACTORS OF THE ELEMENTARY REFLECTORS (SEE FURTHER
*          DETAILS).
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE.
*
*  FURTHER DETAILS
*  ===============
*
*  IF UPLO = 'U', THE MATRIX Q IS REPRESENTED AS A PRODUCT OF ELEMENTARY
*  REFLECTORS
*
*     Q = H(N-1) . . . H(2) H(1).
*
*  EACH H(I) HAS THE FORM
*
*     H(I) = I - TAU * V * V'
*
*  WHERE TAU IS A REAL SCALAR, AND V IS A REAL VECTOR WITH
*  V(I+1:N) = 0 AND V(I) = 1; V(1:I-1) IS STORED ON EXIT IN
*  A(1:I-1,I+1), AND TAU IN TAU(I).
*
*  IF UPLO = 'L', THE MATRIX Q IS REPRESENTED AS A PRODUCT OF ELEMENTARY
*  REFLECTORS
*
*     Q = H(1) H(2) . . . H(N-1).
*
*  EACH H(I) HAS THE FORM
*
*     H(I) = I - TAU * V * V'
*
*  WHERE TAU IS A REAL SCALAR, AND V IS A REAL VECTOR WITH
*  V(1:I) = 0 AND V(I+1) = 1; V(I+2:N) IS STORED ON EXIT IN A(I+2:N,I),
*  AND TAU IN TAU(I).
*
*  THE CONTENTS OF A ON EXIT ARE ILLUSTRATED BY THE FOLLOWING EXAMPLES
*  WITH N = 5:
*
*  IF UPLO = 'U':                       IF UPLO = 'L':
*
*    (  D   E   V2  V3  V4 )              (  D                  )
*    (      D   E   V3  V4 )              (  E   D              )
*    (          D   E   V4 )              (  V1  E   D          )
*    (              D   E  )              (  V1  V2  E   D      )
*    (                  D  )              (  V1  V2  V3  E   D  )
*
*  WHERE D AND E DENOTE DIAGONAL AND OFF-DIAGONAL ELEMENTS OF T, AND VI
*  DENOTES AN ELEMENT OF THE VECTOR DEFINING H(I).
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO, HALF
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0,
     $                   HALF = 1.0D0 / 2.0D0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            UPPER
      INTEGER            I
      DOUBLE PRECISION   ALPHA, TAUI
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DAXPY, DLARFG, DSYMV, DSYR2, XERBLA
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTD2', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        REDUCE THE UPPER TRIANGLE OF A
*
         DO 10 I = N - 1, 1, -1
*
*           GENERATE ELEMENTARY REFLECTOR H(I) = I - TAU * V * V'
*           TO ANNIHILATE A(1:I-1,I+1)
*
            CALL DLARFG( I, A( I, I+1 ), A( 1, I+1 ), 1, TAUI )
            E( I ) = A( I, I+1 )
*
            IF( TAUI.NE.ZERO ) THEN
*
*              APPLY H(I) FROM BOTH SIDES TO A(1:I,1:I)
*
               A( I, I+1 ) = ONE
*
*              COMPUTE  X := TAU * A * V  STORING X IN TAU(1:I)
*
               CALL DSYMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO,
     $                     TAU, 1 )
*
*              COMPUTE  W := X - 1/2 * TAU * (X'*V) * V
*
               ALPHA = -HALF*TAUI*DDOT( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL DAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
*
*              APPLY THE TRANSFORMATION AS A RANK-2 UPDATE:
*                 A := A - V * W' - W * V'
*
               CALL DSYR2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A,
     $                     LDA )
*
               A( I, I+1 ) = E( I )
            END IF
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE
*
*        REDUCE THE LOWER TRIANGLE OF A
*
         DO 20 I = 1, N - 1
*
*           GENERATE ELEMENTARY REFLECTOR H(I) = I - TAU * V * V'
*           TO ANNIHILATE A(I+2:N,I)
*
            CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                   TAUI )
            E( I ) = A( I+1, I )
*
            IF( TAUI.NE.ZERO ) THEN
*
*              APPLY H(I) FROM BOTH SIDES TO A(I+1:N,I+1:N)
*
               A( I+1, I ) = ONE
*
*              COMPUTE  X := TAU * A * V  STORING Y IN TAU(I:N-1)
*
               CALL DSYMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, TAU( I ), 1 )
*
*              COMPUTE  W := X - 1/2 * TAU * (X'*V) * V
*
               ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), 1, A( I+1, I ),
     $                 1 )
               CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
*
*              APPLY THE TRANSFORMATION AS A RANK-2 UPDATE:
*                 A := A - V * W' - W * V'
*
               CALL DSYR2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1,
     $                     A( I+1, I+1 ), LDA )
*
               A( I+1, I ) = E( I )
            END IF
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF
*
      RETURN
*
*     END OF DSYTD2
*
      END
CUT HERE............
CAT > DLATRD.F <<'CUT HERE............'
      SUBROUTINE DLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), E( * ), TAU( * ), W( LDW, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLATRD REDUCES NB ROWS AND COLUMNS OF A REAL SYMMETRIC MATRIX A TO
*  SYMMETRIC TRIDIAGONAL FORM BY AN ORTHOGONAL SIMILARITY
*  TRANSFORMATION Q' * A * Q, AND RETURNS THE MATRICES V AND W WHICH ARE
*  NEEDED TO APPLY THE TRANSFORMATION TO THE UNREDUCED PART OF A.
*
*  IF UPLO = 'U', DLATRD REDUCES THE LAST NB ROWS AND COLUMNS OF A
*  MATRIX, OF WHICH THE UPPER TRIANGLE IS SUPPLIED;
*  IF UPLO = 'L', DLATRD REDUCES THE FIRST NB ROWS AND COLUMNS OF A
*  MATRIX, OF WHICH THE LOWER TRIANGLE IS SUPPLIED.
*
*  THIS IS AN AUXILIARY ROUTINE CALLED BY DSYTRD.
*
*  ARGUMENTS
*  =========
*
*  UPLO    (INPUT) CHARACTER
*          SPECIFIES WHETHER THE UPPER OR LOWER TRIANGULAR PART OF THE
*          SYMMETRIC MATRIX A IS STORED:
*          = 'U': UPPER TRIANGULAR
*          = 'L': LOWER TRIANGULAR
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX A.
*
*  NB      (INPUT) INTEGER
*          THE NUMBER OF ROWS AND COLUMNS TO BE REDUCED.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE SYMMETRIC MATRIX A.  IF UPLO = 'U', THE LEADING
*          N-BY-N UPPER TRIANGULAR PART OF A CONTAINS THE UPPER
*          TRIANGULAR PART OF THE MATRIX A, AND THE STRICTLY LOWER
*          TRIANGULAR PART OF A IS NOT REFERENCED.  IF UPLO = 'L', THE
*          LEADING N-BY-N LOWER TRIANGULAR PART OF A CONTAINS THE LOWER
*          TRIANGULAR PART OF THE MATRIX A, AND THE STRICTLY UPPER
*          TRIANGULAR PART OF A IS NOT REFERENCED.
*          ON EXIT:
*          IF UPLO = 'U', THE LAST NB COLUMNS HAVE BEEN REDUCED TO
*            TRIDIAGONAL FORM, WITH THE DIAGONAL ELEMENTS OVERWRITING
*            THE DIAGONAL ELEMENTS OF A; THE ELEMENTS ABOVE THE DIAGONAL
*            WITH THE ARRAY TAU, REPRESENT THE ORTHOGONAL MATRIX Q AS A
*            PRODUCT OF ELEMENTARY REFLECTORS;
*          IF UPLO = 'L', THE FIRST NB COLUMNS HAVE BEEN REDUCED TO
*            TRIDIAGONAL FORM, WITH THE DIAGONAL ELEMENTS OVERWRITING
*            THE DIAGONAL ELEMENTS OF A; THE ELEMENTS BELOW THE DIAGONAL
*            WITH THE ARRAY TAU, REPRESENT THE  ORTHOGONAL MATRIX Q AS A
*            PRODUCT OF ELEMENTARY REFLECTORS.
*          SEE FURTHER DETAILS.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= (1,N).
*
*  E       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          IF UPLO = 'U', E(N-NB:N-1) CONTAINS THE SUPERDIAGONAL
*          ELEMENTS OF THE LAST NB COLUMNS OF THE REDUCED MATRIX;
*          IF UPLO = 'L', E(1:NB) CONTAINS THE SUBDIAGONAL ELEMENTS OF
*          THE FIRST NB COLUMNS OF THE REDUCED MATRIX.
*
*  TAU     (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N-1)
*          THE SCALAR FACTORS OF THE ELEMENTARY REFLECTORS, STORED IN
*          TAU(N-NB:N-1) IF UPLO = 'U', AND IN TAU(1:NB) IF UPLO = 'L'.
*          SEE FURTHER DETAILS.
*
*  W       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDW,NB)
*          THE N-BY-NB MATRIX W REQUIRED TO UPDATE THE UNREDUCED PART
*          OF A.
*
*  LDW     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY W. LDW >= MAX(1,N).
*
*  FURTHER DETAILS
*  ===============
*
*  IF UPLO = 'U', THE MATRIX Q IS REPRESENTED AS A PRODUCT OF ELEMENTARY
*  REFLECTORS
*
*     Q = H(N) H(N-1) . . . H(N-NB+1).
*
*  EACH H(I) HAS THE FORM
*
*     H(I) = I - TAU * V * V'
*
*  WHERE TAU IS A REAL SCALAR, AND V IS A REAL VECTOR WITH
*  V(I:N) = 0 AND V(I-1) = 1; V(1:I-1) IS STORED ON EXIT IN A(1:I-1,I),
*  AND TAU IN TAU(I-1).
*
*  IF UPLO = 'L', THE MATRIX Q IS REPRESENTED AS A PRODUCT OF ELEMENTARY
*  REFLECTORS
*
*     Q = H(1) H(2) . . . H(NB).
*
*  EACH H(I) HAS THE FORM
*
*     H(I) = I - TAU * V * V'
*
*  WHERE TAU IS A REAL SCALAR, AND V IS A REAL VECTOR WITH
*  V(1:I) = 0 AND V(I+1) = 1; V(I+1:N) IS STORED ON EXIT IN A(I+1:N,I),
*  AND TAU IN TAU(I).
*
*  THE ELEMENTS OF THE VECTORS V TOGETHER FORM THE N-BY-NB MATRIX V
*  WHICH IS NEEDED, WITH W, TO APPLY THE TRANSFORMATION TO THE UNREDUCED
*  PART OF THE MATRIX, USING A SYMMETRIC RANK-2K UPDATE OF THE FORM:
*  A := A - V*W' - W*V'.
*
*  THE CONTENTS OF A ON EXIT ARE ILLUSTRATED BY THE FOLLOWING EXAMPLES
*  WITH N = 5 AND NB = 2:
*
*  IF UPLO = 'U':                       IF UPLO = 'L':
*
*    (  A   A   A   V4  V5 )              (  D                  )
*    (      A   A   V4  V5 )              (  1   D              )
*    (          A   1   V5 )              (  V1  1   A          )
*    (              D   1  )              (  V1  V2  A   A      )
*    (                  D  )              (  V1  V2  A   A   A  )
*
*  WHERE D DENOTES A DIAGONAL ELEMENT OF THE REDUCED MATRIX, A DENOTES
*  AN ELEMENT OF THE ORIGINAL MATRIX THAT IS UNCHANGED, AND VI DENOTES
*  AN ELEMENT OF THE VECTOR DEFINING H(I).
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, IW
      DOUBLE PRECISION   ALPHA
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DAXPY, DGEMV, DLARFG, DSCAL, DSYMV
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        REDUCE LAST NB COLUMNS OF UPPER TRIANGLE
*
         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN
*
*              UPDATE A(1:I,I)
*
               CALL DGEMV( 'NO TRANSPOSE', I, N-I, -ONE, A( 1, I+1 ),
     $                     LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL DGEMV( 'NO TRANSPOSE', I, N-I, -ONE, W( 1, IW+1 ),
     $                     LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
            END IF
            IF( I.GT.1 ) THEN
*
*              GENERATE ELEMENTARY REFLECTOR H(I) TO ANNIHILATE
*              A(1:I-2,I)
*
               CALL DLARFG( I-1, A( I-1, I ), A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = A( I-1, I )
               A( I-1, I ) = ONE
*
*              COMPUTE W(1:I-1,I)
*
               CALL DSYMV( 'UPPER', I-1, ONE, A, LDA, A( 1, I ), 1,
     $                     ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL DGEMV( 'TRANSPOSE', I-1, N-I, ONE, W( 1, IW+1 ),
     $                        LDW, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'NO TRANSPOSE', I-1, N-I, -ONE,
     $                        A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
                  CALL DGEMV( 'TRANSPOSE', I-1, N-I, ONE, A( 1, I+1 ),
     $                        LDA, A( 1, I ), 1, ZERO, W( I+1, IW ), 1 )
                  CALL DGEMV( 'NO TRANSPOSE', I-1, N-I, -ONE,
     $                        W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
               END IF
               CALL DSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*DDOT( I-1, W( 1, IW ), 1,
     $                 A( 1, I ), 1 )
               CALL DAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF
*
   10    CONTINUE
      ELSE
*
*        REDUCE FIRST NB COLUMNS OF LOWER TRIANGLE
*
         DO 20 I = 1, NB
*
*           UPDATE A(I:N,I)
*
            CALL DGEMV( 'NO TRANSPOSE', N-I+1, I-1, -ONE, A( I, 1 ),
     $                  LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL DGEMV( 'NO TRANSPOSE', N-I+1, I-1, -ONE, W( I, 1 ),
     $                  LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            IF( I.LT.N ) THEN
*
*              GENERATE ELEMENTARY REFLECTOR H(I) TO ANNIHILATE
*              A(I+2:N,I)
*
               CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1,
     $                      TAU( I ) )
               E( I ) = A( I+1, I )
               A( I+1, I ) = ONE
*
*              COMPUTE W(I+1:N,I)
*
               CALL DSYMV( 'LOWER', N-I, ONE, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
               CALL DGEMV( 'TRANSPOSE', N-I, I-1, ONE, W( I+1, 1 ), LDW,
     $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'NO TRANSPOSE', N-I, I-1, -ONE, A( I+1, 1 ),
     $                     LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DGEMV( 'TRANSPOSE', N-I, I-1, ONE, A( I+1, 1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( 1, I ), 1 )
               CALL DGEMV( 'NO TRANSPOSE', N-I, I-1, -ONE, W( I+1, 1 ),
     $                     LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL DSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*DDOT( N-I, W( I+1, I ), 1,
     $                 A( I+1, I ), 1 )
               CALL DAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF
*
   20    CONTINUE
      END IF
*
      RETURN
*
*     END OF DLATRD
*
      END
CUT HERE............
CAT > DLARFG.F <<'CUT HERE............'
      SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 29, 1992
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INCX, N
      DOUBLE PRECISION   ALPHA, TAU
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLARFG GENERATES A REAL ELEMENTARY REFLECTOR H OF ORDER N, SUCH
*  THAT
*
*        H * ( ALPHA ) = ( BETA ),   H' * H = I.
*            (   X   )   (   0  )
*
*  WHERE ALPHA AND BETA ARE SCALARS, AND X IS AN (N-1)-ELEMENT REAL
*  VECTOR. H IS REPRESENTED IN THE FORM
*
*        H = I - TAU * ( 1 ) * ( 1 V' ) ,
*                      ( V )
*
*  WHERE TAU IS A REAL SCALAR AND V IS A REAL (N-1)-ELEMENT
*  VECTOR.
*
*  IF THE ELEMENTS OF X ARE ALL ZERO, THEN TAU = 0 AND H IS TAKEN TO BE
*  THE UNIT MATRIX.
*
*  OTHERWISE  1 <= TAU <= 2.
*
*  ARGUMENTS
*  =========
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE ELEMENTARY REFLECTOR.
*
*  ALPHA   (INPUT/OUTPUT) DOUBLE PRECISION
*          ON ENTRY, THE VALUE ALPHA.
*          ON EXIT, IT IS OVERWRITTEN WITH THE VALUE BETA.
*
*  X       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION
*                         (1+(N-2)*ABS(INCX))
*          ON ENTRY, THE VECTOR X.
*          ON EXIT, IT IS OVERWRITTEN WITH THE VECTOR V.
*
*  INCX    (INPUT) INTEGER
*          THE INCREMENT BETWEEN ELEMENTS OF X. INCX <> 0.
*
*  TAU     (OUTPUT) DOUBLE PRECISION
*          THE VALUE TAU.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            J, KNT
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. EXTERNAL FUNCTIONS ..
      DOUBLE PRECISION   DLAMCH, DLAPY2, DNRM2
      EXTERNAL           DLAMCH, DLAPY2, DNRM2
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, SIGN
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DSCAL
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      IF( N.LE.1 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = DNRM2( N-1, X, INCX )
*
      IF( XNORM.EQ.ZERO ) THEN
*
*        H  =  I
*
         TAU = ZERO
      ELSE
*
*        GENERAL CASE
*
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = DLAMCH( 'S' )
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA MAY BE INACCURATE; SCALE X AND RECOMPUTE THEM
*
            RSAFMN = ONE / SAFMIN
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL DSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
*           NEW BETA IS AT MOST 1, AT LEAST SAFMIN
*
            XNORM = DNRM2( N-1, X, INCX )
            BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
            TAU = ( BETA-ALPHA ) / BETA
            CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
*
*           IF ALPHA IS SUBNORMAL, IT MAY LOSE RELATIVE ACCURACY
*
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = ( BETA-ALPHA ) / BETA
            CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
            ALPHA = BETA
         END IF
      END IF
*
      RETURN
*
*     END OF DLARFG
*
      END
CUT HERE............
CAT > DLAMCH.F <<'CUT HERE............'
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          CMACH
*     ..
*
*  PURPOSE
*  =======
*
*  DLAMCH DETERMINES DOUBLE PRECISION MACHINE PARAMETERS.
*
*  ARGUMENTS
*  =========
*
*  CMACH   (INPUT) CHARACTER*1
*          SPECIFIES THE VALUE TO BE RETURNED BY DLAMCH:
*          = 'E' OR 'E',   DLAMCH := EPS
*          = 'S' OR 'S ,   DLAMCH := SFMIN
*          = 'B' OR 'B',   DLAMCH := BASE
*          = 'P' OR 'P',   DLAMCH := EPS*BASE
*          = 'N' OR 'N',   DLAMCH := T
*          = 'R' OR 'R',   DLAMCH := RND
*          = 'M' OR 'M',   DLAMCH := EMIN
*          = 'U' OR 'U',   DLAMCH := RMIN
*          = 'L' OR 'L',   DLAMCH := EMAX
*          = 'O' OR 'O',   DLAMCH := RMAX
*
*          WHERE
*
*          EPS   = RELATIVE MACHINE PRECISION
*          SFMIN = SAFE MINIMUM, SUCH THAT 1/SFMIN DOES NOT OVERFLOW
*          BASE  = BASE OF THE MACHINE
*          PREC  = EPS*BASE
*          T     = NUMBER OF (BASE) DIGITS IN THE MANTISSA
*          RND   = 1.0 WHEN ROUNDING OCCURS IN ADDITION, 0.0 OTHERWISE
*          EMIN  = MINIMUM EXPONENT BEFORE (GRADUAL) UNDERFLOW
*          RMIN  = UNDERFLOW THRESHOLD - BASE**(EMIN-1)
*          EMAX  = LARGEST EXPONENT BEFORE OVERFLOW
*          RMAX  = OVERFLOW THRESHOLD  - (BASE**EMAX)*(1-EPS)
*
* =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLAMC2
*     ..
*     .. SAVE STATEMENT ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. DATA STATEMENTS ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
*
*           USE SMALL PLUS A BIT, TO AVOID THE POSSIBILITY OF ROUNDING
*           CAUSING OVERFLOW WHEN COMPUTING  1/SFMIN.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
*
      DLAMCH = RMACH
      RETURN
*
*     END OF DLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  PURPOSE
*  =======
*
*  DLAMC1 DETERMINES THE MACHINE PARAMETERS GIVEN BY BETA, T, RND, AND
*  IEEE1.
*
*  ARGUMENTS
*  =========
*
*  BETA    (OUTPUT) INTEGER
*          THE BASE OF THE MACHINE.
*
*  T       (OUTPUT) INTEGER
*          THE NUMBER OF ( BETA ) DIGITS IN THE MANTISSA.
*
*  RND     (OUTPUT) LOGICAL
*          SPECIFIES WHETHER PROPER ROUNDING  ( RND = .TRUE. )  OR
*          CHOPPING  ( RND = .FALSE. )  OCCURS IN ADDITION. THIS MAY NOT
*          BE A RELIABLE GUIDE TO THE WAY IN WHICH THE MACHINE PERFORMS
*          ITS ARITHMETIC.
*
*  IEEE1   (OUTPUT) LOGICAL
*          SPECIFIES WHETHER ROUNDING APPEARS TO BE DONE IN THE IEEE
*          'ROUND TO NEAREST' STYLE.
*
*  FURTHER DETAILS
*  ===============
*
*  THE ROUTINE IS BASED ON THE ROUTINE  ENVRON  BY MALCOLM AND
*  INCORPORATES SUGGESTIONS BY GENTLEMAN AND MAROVICH. SEE
*
*     MALCOLM M. A. (1972) ALGORITHMS TO REVEAL PROPERTIES OF
*        FLOATING-POINT ARITHMETIC. COMMS. OF THE ACM, 15, 949-951.
*
*     GENTLEMAN W. M. AND MAROVICH S. B. (1974) MORE ON ALGORITHMS
*        THAT REVEAL PROPERTIES OF FLOATING POINT ARITHMETIC UNITS.
*        COMMS. OF THE ACM, 17, 276-277.
*
* =====================================================================
*
*     .. LOCAL SCALARS ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. EXTERNAL FUNCTIONS ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. SAVE STATEMENT ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. DATA STATEMENTS ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
*
*        LBETA,  LIEEE1,  LT AND  LRND  ARE THE  LOCAL VALUES  OF  BETA,
*        IEEE1, T AND RND.
*
*        THROUGHOUT THIS ROUTINE  WE USE THE FUNCTION  DLAMC3  TO ENSURE
*        THAT RELEVANT VALUES ARE  STORED AND NOT HELD IN REGISTERS,  OR
*        ARE NOT AFFECTED BY OPTIMIZERS.
*
*        COMPUTE  A = 2.0**M  WITH THE  SMALLEST POSITIVE INTEGER M SUCH
*        THAT
*
*           FL( A + 1.0 ) = A.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        NOW COMPUTE  B = 2.0**M  WITH THE SMALLEST POSITIVE INTEGER M
*        SUCH THAT
*
*           FL( A + B ) .GT. A.
*
         B = 1
         C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        NOW COMPUTE THE BASE.  A AND C  ARE NEIGHBOURING FLOATING POINT
*        NUMBERS  IN THE  INTERVAL  ( BETA**T, BETA**( T + 1 ) )  AND SO
*        THEIR DIFFERENCE IS BETA. ADDING 0.25 TO C IS TO ENSURE THAT IT
*        IS TRUNCATED TO BETA AND NOT ( BETA - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
*
*        NOW DETERMINE WHETHER ROUNDING OR CHOPPING OCCURS,  BY ADDING A
*        BIT  LESS  THAN  BETA/2  AND A  BIT  MORE  THAN  BETA/2  TO  A.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        TRY AND DECIDE WHETHER ROUNDING IS DONE IN THE  IEEE  'ROUND TO
*        NEAREST' STYLE. B/2 IS HALF A UNIT IN THE LAST PLACE OF THE TWO
*        NUMBERS A AND SAVEC. FURTHERMORE, A IS EVEN, I.E. HAS LAST  BIT
*        ZERO, AND SAVEC IS ODD. THUS ADDING B/2 TO A SHOULD NOT  CHANGE
*        A, BUT ADDING B/2 TO SAVEC SHOULD CHANGE SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        NOW FIND  THE  MANTISSA, T.  IT SHOULD  BE THE  INTEGER PART OF
*        LOG TO THE BASE BETA OF A,  HOWEVER IT IS SAFER TO DETERMINE  T
*        BY POWERING.  SO WE FIND T AS THE SMALLEST POSITIVE INTEGER FOR
*        WHICH
*
*           FL( BETA**T + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
*
*     END OF DLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
*
*  PURPOSE
*  =======
*
*  DLAMC2 DETERMINES THE MACHINE PARAMETERS SPECIFIED IN ITS ARGUMENT
*  LIST.
*
*  ARGUMENTS
*  =========
*
*  BETA    (OUTPUT) INTEGER
*          THE BASE OF THE MACHINE.
*
*  T       (OUTPUT) INTEGER
*          THE NUMBER OF ( BETA ) DIGITS IN THE MANTISSA.
*
*  RND     (OUTPUT) LOGICAL
*          SPECIFIES WHETHER PROPER ROUNDING  ( RND = .TRUE. )  OR
*          CHOPPING  ( RND = .FALSE. )  OCCURS IN ADDITION. THIS MAY NOT
*          BE A RELIABLE GUIDE TO THE WAY IN WHICH THE MACHINE PERFORMS
*          ITS ARITHMETIC.
*
*  EPS     (OUTPUT) DOUBLE PRECISION
*          THE SMALLEST POSITIVE NUMBER SUCH THAT
*
*             FL( 1.0 - EPS ) .LT. 1.0,
*
*          WHERE FL DENOTES THE COMPUTED VALUE.
*
*  EMIN    (OUTPUT) INTEGER
*          THE MINIMUM EXPONENT BEFORE (GRADUAL) UNDERFLOW OCCURS.
*
*  RMIN    (OUTPUT) DOUBLE PRECISION
*          THE SMALLEST NORMALIZED NUMBER FOR THE MACHINE, GIVEN BY
*          BASE**( EMIN - 1 ), WHERE  BASE  IS THE FLOATING POINT VALUE
*          OF BETA.
*
*  EMAX    (OUTPUT) INTEGER
*          THE MAXIMUM EXPONENT BEFORE OVERFLOW OCCURS.
*
*  RMAX    (OUTPUT) DOUBLE PRECISION
*          THE LARGEST POSITIVE NUMBER FOR THE MACHINE, GIVEN BY
*          BASE**EMAX * ( 1 - EPS ), WHERE  BASE  IS THE FLOATING POINT
*          VALUE OF BETA.
*
*  FURTHER DETAILS
*  ===============
*
*  THE COMPUTATION OF  EPS  IS BASED ON A ROUTINE PARANOIA BY
*  W. KAHAN OF THE UNIVERSITY OF CALIFORNIA AT BERKELEY.
*
* =====================================================================
*
*     .. LOCAL SCALARS ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. EXTERNAL FUNCTIONS ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. SAVE STATEMENT ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. DATA STATEMENTS ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN AND LRMIN  ARE THE LOCAL VALUES OF
*        BETA, T, RND, EPS, EMIN AND RMIN.
*
*        THROUGHOUT THIS ROUTINE  WE USE THE FUNCTION  DLAMC3  TO ENSURE
*        THAT RELEVANT VALUES ARE STORED  AND NOT HELD IN REGISTERS,  OR
*        ARE NOT AFFECTED BY OPTIMIZERS.
*
*        DLAMC1 RETURNS THE PARAMETERS  LBETA, LT, LRND AND LIEEE1.
*
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        START TO FIND EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        TRY SOME TRICKS TO SEE WHETHER OR NOT THIS IS THE CORRECT  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        COMPUTATION OF EPS COMPLETE.
*
*        NOW FIND  EMIN.  LET A = + OR - 1, AND + OR - (1 + BASE**(-3)).
*        KEEP DIVIDING  A BY BETA UNTIL (GRADUAL) UNDERFLOW OCCURS. THIS
*        IS DETECTED WHEN WE CANNOT RECOVER THE PREVIOUS A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( NON TWOS-COMPLEMENT MACHINES, NO GRADUAL UNDERFLOW;
*              E.G.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( NON TWOS-COMPLEMENT MACHINES, WITH GRADUAL UNDERFLOW;
*              E.G., IEEE STANDARD FOLLOWERS )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A GUESS; NO KNOWN MACHINE )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( TWOS-COMPLEMENT MACHINES, NO GRADUAL UNDERFLOW;
*              E.G., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A GUESS; NO KNOWN MACHINE )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( TWOS-COMPLEMENT MACHINES WITH GRADUAL UNDERFLOW;
*              NO KNOWN MACHINE )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A GUESS; NO KNOWN MACHINE )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A GUESS; NO KNOWN MACHINE )
            IWARN = .TRUE.
         END IF
***
* COMMENT OUT THIS IF BLOCK IF EMIN IS OK
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        ASSUME IEEE ARITHMETIC IF WE FOUND DENORMALISED  NUMBERS ABOVE,
*        OR IF ARITHMETIC SEEMS TO ROUND IN THE  IEEE STYLE,  DETERMINED
*        IN ROUTINE DLAMC1. A TRUE IEEE MACHINE SHOULD HAVE BOTH  THINGS
*        TRUE; HOWEVER, FAULTY MACHINES MAY HAVE ONE OR THE OTHER.
*
         IEEE = IEEE .OR. LIEEE1
*
*        COMPUTE  RMIN BY SUCCESSIVE DIVISION BY  BETA. WE COULD COMPUTE
*        RMIN AS BASE**( EMIN - 1 ),  BUT SOME MACHINES UNDERFLOW DURING
*        THIS COMPUTATION.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        FINALLY, CALL DLAMC5 TO COMPUTE EMAX AND RMAX.
*
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. THE VALUE EMIN MAY BE INCORRECT:-',
     $      '  EMIN = ', I8, /
     $      ' IF, AFTER INSPECTION, THE VALUE EMIN LOOKS',
     $      ' ACCEPTABLE PLEASE COMMENT OUT ',
     $      / ' THE IF BLOCK AS MARKED WITHIN THE CODE OF ROUTINE',
     $      ' DLAMC2,', / ' OTHERWISE SUPPLY EMIN EXPLICITLY.', / )
*
*     END OF DLAMC2
*
      END
*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   A, B
*     ..
*
*  PURPOSE
*  =======
*
*  DLAMC3  IS INTENDED TO FORCE  A  AND  B  TO BE STORED PRIOR TO DOING
*  THE ADDITION OF  A  AND  B ,  FOR USE IN SITUATIONS WHERE OPTIMIZERS
*  MIGHT HOLD ONE OF THESE IN A REGISTER.
*
*  ARGUMENTS
*  =========
*
*  A, B    (INPUT) DOUBLE PRECISION
*          THE VALUES A AND B.
*
* =====================================================================
*
*     .. EXECUTABLE STATEMENTS ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     END OF DLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
*
*  PURPOSE
*  =======
*
*  DLAMC4 IS A SERVICE ROUTINE FOR DLAMC2.
*
*  ARGUMENTS
*  =========
*
*  EMIN    (OUTPUT) EMIN
*          THE MINIMUM EXPONENT BEFORE (GRADUAL) UNDERFLOW, COMPUTED BY
*          SETTING A = START AND DIVIDING BY BASE UNTIL THE PREVIOUS A
*          CAN NOT BE RECOVERED.
*
*  START   (INPUT) DOUBLE PRECISION
*          THE STARTING POINT FOR DETERMINING EMIN.
*
*  BASE    (INPUT) INTEGER
*          THE BASE OF THE MACHINE.
*
* =====================================================================
*
*     .. LOCAL SCALARS ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. EXTERNAL FUNCTIONS ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     END OF DLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
*
*  PURPOSE
*  =======
*
*  DLAMC5 ATTEMPTS TO COMPUTE RMAX, THE LARGEST MACHINE FLOATING-POINT
*  NUMBER, WITHOUT OVERFLOW.  IT ASSUMES THAT EMAX + ABS(EMIN) SUM
*  APPROXIMATELY TO A POWER OF 2.  IT WILL FAIL ON MACHINES WHERE THIS
*  ASSUMPTION DOES NOT HOLD, FOR EXAMPLE, THE CYBER 205 (EMIN = -28625,
*  EMAX = 28718).  IT WILL ALSO FAIL IF THE VALUE SUPPLIED FOR EMIN IS
*  TOO LARGE (I.E. TOO CLOSE TO ZERO), PROBABLY WITH OVERFLOW.
*
*  ARGUMENTS
*  =========
*
*  BETA    (INPUT) INTEGER
*          THE BASE OF FLOATING-POINT ARITHMETIC.
*
*  P       (INPUT) INTEGER
*          THE NUMBER OF BASE BETA DIGITS IN THE MANTISSA OF A
*          FLOATING-POINT VALUE.
*
*  EMIN    (INPUT) INTEGER
*          THE MINIMUM EXPONENT BEFORE (GRADUAL) UNDERFLOW.
*
*  IEEE    (INPUT) LOGICAL
*          A LOGICAL FLAG SPECIFYING WHETHER OR NOT THE ARITHMETIC
*          SYSTEM IS THOUGHT TO COMPLY WITH THE IEEE STANDARD.
*
*  EMAX    (OUTPUT) INTEGER
*          THE LARGEST EXPONENT BEFORE OVERFLOW
*
*  RMAX    (OUTPUT) DOUBLE PRECISION
*          THE LARGEST MACHINE FLOATING-POINT NUMBER.
*
* =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. EXTERNAL FUNCTIONS ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MOD
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     FIRST COMPUTE LEXP AND UEXP, TWO POWERS OF 2 THAT BOUND
*     ABS(EMIN). WE THEN ASSUME THAT EMAX + ABS(EMIN) WILL SUM
*     APPROXIMATELY TO THE BOUND THAT IS CLOSEST TO ABS(EMIN).
*     (EMAX IS THE EXPONENT OF THE REQUIRED NUMBER RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     NOW -LEXP IS LESS THAN OR EQUAL TO EMIN, AND -UEXP IS GREATER
*     THAN OR EQUAL TO EMIN. EXBITS IS THE NUMBER OF BITS NEEDED TO
*     STORE THE EXPONENT.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM IS THE EXPONENT RANGE, APPROXIMATELY EQUAL TO
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS IS THE TOTAL NUMBER OF BITS NEEDED TO STORE A
*     FLOATING-POINT NUMBER.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        EITHER THERE ARE AN ODD NUMBER OF BITS USED TO STORE A
*        FLOATING-POINT NUMBER, WHICH IS UNLIKELY, OR SOME BITS ARE
*        NOT USED IN THE REPRESENTATION OF NUMBERS, WHICH IS POSSIBLE,
*        (E.G. CRAY MACHINES) OR THE MANTISSA HAS AN IMPLICIT BIT,
*        (E.G. IEEE MACHINES, DEC VAX MACHINES), WHICH IS PERHAPS THE
*        MOST LIKELY. WE HAVE TO ASSUME THE LAST ALTERNATIVE.
*        IF THIS IS TRUE, THEN WE NEED TO REDUCE EMAX BY ONE BECAUSE
*        THERE MUST BE SOME WAY OF REPRESENTING ZERO IN AN IMPLICIT-BIT
*        SYSTEM. ON MACHINES LIKE CRAY, WE ARE REDUCING EMAX BY ONE
*        UNNECESSARILY.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        ASSUME WE ARE ON AN IEEE MACHINE WHICH RESERVES ONE EXPONENT
*        FOR INFINITY AND NAN.
*
         EMAX = EMAX - 1
      END IF
*
*     NOW CREATE RMAX, THE LARGEST MACHINE NUMBER, WHICH SHOULD
*     BE EQUAL TO (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     FIRST COMPUTE 1.0 - BETA**(-P), BEING CAREFUL THAT THE
*     RESULT IS LESS THAN 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     NOW MULTIPLY BY BETA**EMAX TO GET RMAX.
*
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     END OF DLAMC5
*
      END
CUT HERE............
CAT > DLAPY2.F <<'CUT HERE............'
      DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION   X, Y
*     ..
*
*  PURPOSE
*  =======
*
*  DLAPY2 RETURNS SQRT(X**2+Y**2), TAKING CARE NOT TO CAUSE UNNECESSARY
*  OVERFLOW.
*
*  ARGUMENTS
*  =========
*
*  X       (INPUT) DOUBLE PRECISION
*  Y       (INPUT) DOUBLE PRECISION
*          X AND Y SPECIFY THE VALUES X AND Y.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
*     ..
*     .. LOCAL SCALARS ..
      DOUBLE PRECISION   W, XABS, YABS, Z
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z.EQ.ZERO ) THEN
         DLAPY2 = W
      ELSE
         DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
*
*     END OF DLAPY2
*
      END
CUT HERE............
CAT > ILAENV.F <<'CUT HERE............'
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK AUXILIARY ROUTINE (PRELIMINARY VERSION) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 20, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  PURPOSE
*  =======
*
*  ILAENV IS CALLED FROM THE LAPACK ROUTINES TO CHOOSE PROBLEM-DEPENDENT
*  PARAMETERS FOR THE LOCAL ENVIRONMENT.  SEE ISPEC FOR A DESCRIPTION OF
*  THE PARAMETERS.
*
*  THIS VERSION PROVIDES A SET OF PARAMETERS WHICH SHOULD GIVE GOOD,
*  BUT NOT OPTIMAL, PERFORMANCE ON MANY OF THE CURRENTLY AVAILABLE
*  COMPUTERS.  USERS ARE ENCOURAGED TO MODIFY THIS SUBROUTINE TO SET
*  THE TUNING PARAMETERS FOR THEIR PARTICULAR MACHINE USING THE OPTION
*  AND PROBLEM SIZE INFORMATION IN THE ARGUMENTS.
*
*  THIS ROUTINE WILL NOT FUNCTION CORRECTLY IF IT IS CONVERTED TO ALL
*  LOWER CASE.  CONVERTING IT TO ALL UPPER CASE IS ALLOWED.
*
*  ARGUMENTS
*  =========
*
*  ISPEC   (INPUT) INTEGER
*          SPECIFIES THE PARAMETER TO BE RETURNED AS THE VALUE OF
*          ILAENV.
*          = 1: THE OPTIMAL BLOCKSIZE; IF THIS VALUE IS 1, AN UNBLOCKED
*               ALGORITHM WILL GIVE THE BEST PERFORMANCE.
*          = 2: THE MINIMUM BLOCK SIZE FOR WHICH THE BLOCK ROUTINE
*               SHOULD BE USED; IF THE USABLE BLOCK SIZE IS LESS THAN
*               THIS VALUE, AN UNBLOCKED ROUTINE SHOULD BE USED.
*          = 3: THE CROSSOVER POINT (IN A BLOCK ROUTINE, FOR N LESS
*               THAN THIS VALUE, AN UNBLOCKED ROUTINE SHOULD BE USED)
*          = 4: THE NUMBER OF SHIFTS, USED IN THE NONSYMMETRIC
*               EIGENVALUE ROUTINES
*          = 5: THE MINIMUM COLUMN DIMENSION FOR BLOCKING TO BE USED;
*               RECTANGULAR BLOCKS MUST HAVE DIMENSION AT LEAST K BY M,
*               WHERE K IS GIVEN BY ILAENV(2,...) AND M BY ILAENV(5,...)
*          = 6: THE CROSSOVER POINT FOR THE SVD (WHEN REDUCING AN M BY N
*               MATRIX TO BIDIAGONAL FORM, IF MAX(M,N)/MIN(M,N) EXCEEDS
*               THIS VALUE, A QR FACTORIZATION IS USED FIRST TO REDUCE
*               THE MATRIX TO A TRIANGULAR FORM.)
*          = 7: THE NUMBER OF PROCESSORS
*          = 8: THE CROSSOVER POINT FOR THE MULTISHIFT QR AND QZ METHODS
*               FOR NONSYMMETRIC EIGENVALUE PROBLEMS.
*
*  NAME    (INPUT) CHARACTER*(*)
*          THE NAME OF THE CALLING SUBROUTINE, IN EITHER UPPER CASE OR
*          LOWER CASE.
*
*  OPTS    (INPUT) CHARACTER*(*)
*          THE CHARACTER OPTIONS TO THE SUBROUTINE NAME, CONCATENATED
*          INTO A SINGLE CHARACTER STRING.  FOR EXAMPLE, UPLO = 'U',
*          TRANS = 'T', AND DIAG = 'N' FOR A TRIANGULAR ROUTINE WOULD
*          BE SPECIFIED AS OPTS = 'UTN'.
*
*  N1      (INPUT) INTEGER
*  N2      (INPUT) INTEGER
*  N3      (INPUT) INTEGER
*  N4      (INPUT) INTEGER
*          PROBLEM DIMENSIONS FOR THE SUBROUTINE NAME; THESE MAY NOT ALL
*          BE REQUIRED.
*
* (ILAENV) (OUTPUT) INTEGER
*          >= 0: THE VALUE OF THE PARAMETER SPECIFIED BY ISPEC
*          < 0:  IF ILAENV = -K, THE K-TH ARGUMENT HAD AN ILLEGAL VALUE.
*
*  FURTHER DETAILS
*  ===============
*
*  THE FOLLOWING CONVENTIONS HAVE BEEN USED WHEN CALLING ILAENV FROM THE
*  LAPACK ROUTINES:
*  1)  OPTS IS A CONCATENATION OF ALL OF THE CHARACTER OPTIONS TO
*      SUBROUTINE NAME, IN THE SAME ORDER THAT THEY APPEAR IN THE
*      ARGUMENT LIST FOR NAME, EVEN IF THEY ARE NOT USED IN DETERMINING
*      THE VALUE OF THE PARAMETER SPECIFIED BY ISPEC.
*  2)  THE PROBLEM DIMENSIONS N1, N2, N3, N4 ARE SPECIFIED IN THE ORDER
*      THAT THEY APPEAR IN THE ARGUMENT LIST FOR NAME.  N1 IS USED
*      FIRST, N2 SECOND, AND SO ON, AND UNUSED PROBLEM DIMENSIONS ARE
*      PASSED A VALUE OF -1.
*  3)  THE PARAMETER VALUE RETURNED BY ILAENV IS CHECKED FOR VALIDITY IN
*      THE CALLING SUBROUTINE.  FOR EXAMPLE, ILAENV IS USED TO RETRIEVE
*      THE OPTIMAL BLOCKSIZE FOR STRTRI AS FOLLOWS:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. LOCAL SCALARS ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
*
*     INVALID VALUE FOR ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     CONVERT NAME TO UPPER CASE IF THE FIRST CHARACTER IS LOWER CASE.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII CHARACTER SET
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC CHARACTER SET
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        PRIME MACHINES:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  BLOCK SIZE
*
*     IN THESE EXAMPLES, SEPARATE CODE IS PROVIDED FOR SETTING NB FOR
*     REAL AND COMPLEX.  WE ASSUME THAT NB WILL TAKE THE SAME VALUE IN
*     SINGLE OR DOUBLE PRECISION.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  MINIMUM BLOCK SIZE
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  CROSSOVER POINT
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  NUMBER OF SHIFTS (USED BY XHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  MINIMUM COLUMN DIMENSION (NOT USED)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE
*
*     ISPEC = 6:  CROSSOVER POINT FOR SVD (USED BY XGELSS AND XGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  NUMBER OF PROCESSORS (NOT USED)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  CROSSOVER POINT FOR MULTISHIFT (USED BY XHSEQR)
*
      ILAENV = 50
      RETURN
*
*     END OF ILAENV
*
      END
CUT HERE............
CAT > DLANSY.F <<'CUT HERE............'
      DOUBLE PRECISION FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLANSY  RETURNS THE VALUE OF THE ONE NORM,  OR THE FROBENIUS NORM, OR
*  THE  INFINITY NORM,  OR THE  ELEMENT OF  LARGEST ABSOLUTE VALUE  OF A
*  REAL SYMMETRIC MATRIX A.
*
*  DESCRIPTION
*  ===========
*
*  DLANSY RETURNS THE VALUE
*
*     DLANSY = ( MAX(ABS(A(I,J))), NORM = 'M' OR 'M'
*              (
*              ( NORM1(A),         NORM = '1', 'O' OR 'O'
*              (
*              ( NORMI(A),         NORM = 'I' OR 'I'
*              (
*              ( NORMF(A),         NORM = 'F', 'F', 'E' OR 'E'
*
*  WHERE  NORM1  DENOTES THE  ONE NORM OF A MATRIX (MAXIMUM COLUMN SUM),
*  NORMI  DENOTES THE  INFINITY NORM  OF A MATRIX  (MAXIMUM ROW SUM) AND
*  NORMF  DENOTES THE  FROBENIUS NORM OF A MATRIX (SQUARE ROOT OF SUM OF
*  SQUARES).  NOTE THAT  MAX(ABS(A(I,J)))  IS NOT A  MATRIX NORM.
*
*  ARGUMENTS
*  =========
*
*  NORM    (INPUT) CHARACTER*1
*          SPECIFIES THE VALUE TO BE RETURNED IN DLANSY AS DESCRIBED
*          ABOVE.
*
*  UPLO    (INPUT) CHARACTER*1
*          SPECIFIES WHETHER THE UPPER OR LOWER TRIANGULAR PART OF THE
*          SYMMETRIC MATRIX A IS TO BE REFERENCED.
*          = 'U':  UPPER TRIANGULAR PART OF A IS REFERENCED
*          = 'L':  LOWER TRIANGULAR PART OF A IS REFERENCED
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX A.  N >= 0.  WHEN N = 0, DLANSY IS
*          SET TO ZERO.
*
*  A       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          THE SYMMETRIC MATRIX A.  IF UPLO = 'U', THE LEADING N BY N
*          UPPER TRIANGULAR PART OF A CONTAINS THE UPPER TRIANGULAR PART
*          OF THE MATRIX A, AND THE STRICTLY LOWER TRIANGULAR PART OF A
*          IS NOT REFERENCED.  IF UPLO = 'L', THE LEADING N BY N LOWER
*          TRIANGULAR PART OF A CONTAINS THE LOWER TRIANGULAR PART OF
*          THE MATRIX A, AND THE STRICTLY UPPER TRIANGULAR PART OF A IS
*          NOT REFERENCED.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(N,1).
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (LWORK),
*          WHERE LWORK >= N WHEN NORM = 'I' OR '1' OR 'O'; OTHERWISE,
*          WORK IS NOT REFERENCED.
*
* =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, J
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLASSQ
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        FIND MAX(ABS(A(I,J))).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, J
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               DO 30 I = J, N
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR.
     $         ( NORM.EQ.'1' ) ) THEN
*
*        FIND NORMI(A) ( = NORM1(A), SINCE A IS SYMMETRIC).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( A( J, J ) )
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( A( J, J ) )
               DO 90 I = J + 1, N
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        FIND NORMF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL DLASSQ( J-1, A( 1, J ), 1, SCALE, SUM )
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL DLASSQ( N-J, A( J+1, J ), 1, SCALE, SUM )
  120       CONTINUE
         END IF
         SUM = 2*SUM
         CALL DLASSQ( N, A, LDA+1, SCALE, SUM )
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      DLANSY = VALUE
      RETURN
*
*     END OF DLANSY
*
      END
CUT HERE............
CAT > DLASSQ.F <<'CUT HERE............'
      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLASSQ  RETURNS THE VALUES  SCL  AND  SMSQ  SUCH THAT
*
*     ( SCL**2 )*SMSQ = X( 1 )**2 +...+ X( N )**2 + ( SCALE**2 )*SUMSQ,
*
*  WHERE  X( I ) = X( 1 + ( I - 1 )*INCX ). THE VALUE OF  SUMSQ  IS
*  ASSUMED TO BE NON-NEGATIVE AND  SCL  RETURNS THE VALUE
*
*     SCL = MAX( SCALE, ABS( X( I ) ) ).
*
*  SCALE AND SUMSQ MUST BE SUPPLIED IN SCALE AND SUMSQ AND
*  SCL AND SMSQ ARE OVERWRITTEN ON SCALE AND SUMSQ RESPECTIVELY.
*
*  THE ROUTINE MAKES ONLY ONE PASS THROUGH THE VECTOR X.
*
*  ARGUMENTS
*  =========
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF ELEMENTS TO BE USED FROM THE VECTOR X.
*
*  X       (INPUT) DOUBLE PRECISION
*          THE VECTOR FOR WHICH A SCALED SUM OF SQUARES IS COMPUTED.
*             X( I )  = X( 1 + ( I - 1 )*INCX ), 1 <= I <= N.
*
*  INCX    (INPUT) INTEGER
*          THE INCREMENT BETWEEN SUCCESSIVE VALUES OF THE VECTOR X.
*          INCX > 0.
*
*  SCALE   (INPUT/OUTPUT) DOUBLE PRECISION
*          ON ENTRY, THE VALUE  SCALE  IN THE EQUATION ABOVE.
*          ON EXIT, SCALE IS OVERWRITTEN WITH  SCL , THE SCALING FACTOR
*          FOR THE SUM OF SQUARES.
*
*  SUMSQ   (INPUT/OUTPUT) DOUBLE PRECISION
*          ON ENTRY, THE VALUE  SUMSQ  IN THE EQUATION ABOVE.
*          ON EXIT, SUMSQ IS OVERWRITTEN WITH  SMSQ , THE BASIC SUM OF
*          SQUARES FROM WHICH  SCL  HAS BEEN FACTORED OUT.
*
* =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     END OF DLASSQ
*
      END
CUT HERE............
C    -----------------      BELOW IS DSYTRF  -------------------
CAT > DSYTRF.F <<'CUT HERE..........'
      SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.0B) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 29, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( LWORK )
*     ..
*
*  PURPOSE
*  =======
*
*  DSYTRF COMPUTES THE FACTORIZATION OF A REAL SYMMETRIC MATRIX A USING
*  THE BUNCH-KAUFMAN DIAGONAL PIVOTING METHOD:
*
*     A = U*D*U'  OR  A = L*D*L'
*
*  WHERE U (OR L) IS A PRODUCT OF PERMUTATION AND UNIT UPPER (LOWER)
*  TRIANGULAR MATRICES, U' IS THE TRANSPOSE OF U, AND D IS SYMMETRIC AND
*  BLOCK DIAGONAL WITH 1-BY-1 AND 2-BY-2 DIAGONAL BLOCKS.
*
*  THIS IS THE BLOCKED VERSION OF THE ALGORITHM, CALLING LEVEL 3 BLAS.
*
*  ARGUMENTS
*  =========
*
*  UPLO    (INPUT) CHARACTER*1
*          SPECIFIES WHETHER THE UPPER OR LOWER TRIANGULAR PART OF THE
*          SYMMETRIC MATRIX A IS STORED:
*          = 'U':  UPPER TRIANGULAR
*          = 'L':  LOWER TRIANGULAR
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE SYMMETRIC MATRIX A.  IF UPLO = 'U', THE LEADING
*          N-BY-N UPPER TRIANGULAR PART OF A CONTAINS THE UPPER
*          TRIANGULAR PART OF THE MATRIX A, AND THE STRICTLY LOWER
*          TRIANGULAR PART OF A IS NOT REFERENCED.  IF UPLO = 'L', THE
*          LEADING N-BY-N LOWER TRIANGULAR PART OF A CONTAINS THE LOWER
*          TRIANGULAR PART OF THE MATRIX A, AND THE STRICTLY UPPER
*          TRIANGULAR PART OF A IS NOT REFERENCED.
*
*          ON EXIT, THE BLOCK DIAGONAL MATRIX D AND THE MULTIPLIERS USED
*          TO OBTAIN THE FACTOR U OR L (SEE BELOW FOR FURTHER DETAILS).
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,N).
*
*  IPIV    (OUTPUT) INTEGER ARRAY, DIMENSION (N)
*          DETAILS OF THE INTERCHANGES AND THE BLOCK STRUCTURE OF D.
*          IF IPIV(K) > 0, THEN ROWS AND COLUMNS K AND IPIV(K) WERE
*          INTERCHANGED AND D(K,K) IS A 1-BY-1 DIAGONAL BLOCK.
*          IF UPLO = 'U' AND IPIV(K) = IPIV(K-1) < 0, THEN ROWS AND
*          COLUMNS K-1 AND -IPIV(K) WERE INTERCHANGED AND D(K-1:K,K-1:K)
*          IS A 2-BY-2 DIAGONAL BLOCK.  IF UPLO = 'L' AND IPIV(K) =
*          IPIV(K+1) < 0, THEN ROWS AND COLUMNS K+1 AND -IPIV(K) WERE
*          INTERCHANGED AND D(K:K+1,K:K+1) IS A 2-BY-2 DIAGONAL BLOCK.
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (LWORK)
*          IF INFO RETURNS 0, THEN WORK(1) RETURNS N*NB, THE MINIMUM
*          VALUE OF LWORK REQUIRED TO USE THE OPTIMAL BLOCKSIZE.
*
*  LWORK   (INPUT) INTEGER
*          THE LENGTH OF WORK.  LWORK SHOULD BE >= N*NB, WHERE NB IS THE
*          BLOCK SIZE RETURNED BY ILAENV.
*
*  INFO    (OUTPUT) INTEGER
*          = 0: SUCCESSFUL EXIT
*          < 0: IF INFO = -K, THE K-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0: IF INFO = K, D(K,K) IS EXACTLY ZERO.  THE FACTORIZATION
*               HAS BEEN COMPLETED, BUT THE BLOCK DIAGONAL MATRIX D IS
*               EXACTLY SINGULAR, AND DIVISION BY ZERO WILL OCCUR IF IT
*               IS USED TO SOLVE A SYSTEM OF EQUATIONS.
*
*  FURTHER DETAILS
*  ===============
*
*  IF UPLO = 'U', THEN A = U*D*U', WHERE
*     U = P(N)*U(N)* ... *P(K)U(K)* ...,
*  I.E., U IS A PRODUCT OF TERMS P(K)*U(K), WHERE K DECREASES FROM N TO
*  1 IN STEPS OF 1 OR 2, AND D IS A BLOCK DIAGONAL MATRIX WITH 1-BY-1
*  AND 2-BY-2 DIAGONAL BLOCKS D(K).  P(K) IS A PERMUTATION MATRIX AS
*  DEFINED BY IPIV(K), AND U(K) IS A UNIT UPPER TRIANGULAR MATRIX, SUCH
*  THAT IF THE DIAGONAL BLOCK D(K) IS OF ORDER S (S = 1 OR 2), THEN
*
*             (   I    V    0   )   K-S
*     U(K) =  (   0    I    0   )   S
*             (   0    0    I   )   N-K
*                K-S   S   N-K
*
*  IF S = 1, D(K) OVERWRITES A(K,K), AND V OVERWRITES A(1:K-1,K).
*  IF S = 2, THE UPPER TRIANGLE OF D(K) OVERWRITES A(K-1,K-1), A(K-1,K),
*  AND A(K,K), AND V OVERWRITES A(1:K-2,K-1:K).
*
*  IF UPLO = 'L', THEN A = L*D*L', WHERE
*     L = P(1)*L(1)* ... *P(K)*L(K)* ...,
*  I.E., L IS A PRODUCT OF TERMS P(K)*L(K), WHERE K INCREASES FROM 1 TO
*  N IN STEPS OF 1 OR 2, AND D IS A BLOCK DIAGONAL MATRIX WITH 1-BY-1
*  AND 2-BY-2 DIAGONAL BLOCKS D(K).  P(K) IS A PERMUTATION MATRIX AS
*  DEFINED BY IPIV(K), AND L(K) IS A UNIT LOWER TRIANGULAR MATRIX, SUCH
*  THAT IF THE DIAGONAL BLOCK D(K) IS OF ORDER S (S = 1 OR 2), THEN
*
*             (   I    0     0   )  K-1
*     L(K) =  (   0    I     0   )  S
*             (   0    V     I   )  N-K-S+1
*                K-1   S  N-K-S+1
*
*  IF S = 1, D(K) OVERWRITES A(K,K), AND V OVERWRITES A(K+1:N,K).
*  IF S = 2, THE LOWER TRIANGLE OF D(K) OVERWRITES A(K,K), A(K+1,K),
*  AND A(K+1,K+1), AND V OVERWRITES A(K+2:N,K:K+1).
*
*  =====================================================================
*
*     .. LOCAL SCALARS ..
      LOGICAL            UPPER
      INTEGER            IINFO, IWS, J, K, KB, LDWORK, NB, NBMIN
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLASYF, DSYTF2, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRF', -INFO )
         RETURN
      END IF
*
*     DETERMINE THE BLOCK SIZE
*
      NB = ILAENV( 1, 'DSYTRF', UPLO, N, -1, -1, -1 )
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = LDWORK*NB
         IF( LWORK.LT.IWS ) THEN
            NB = MAX( LWORK / LDWORK, 1 )
            NBMIN = MAX( 2, ILAENV( 2, 'DSYTRF', UPLO, N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = 1
      END IF
      IF( NB.LT.NBMIN )
     $   NB = N
*
      IF( UPPER ) THEN
*
*        FACTORIZE A AS U*D*U' USING THE UPPER TRIANGLE OF A
*
*        K IS THE MAIN LOOP INDEX, DECREASING FROM N TO 1 IN STEPS OF
*        KB, WHERE KB IS THE NUMBER OF COLUMNS FACTORIZED BY DLASYF;
*        KB IS EITHER NB OR NB-1, OR K FOR THE LAST BLOCK
*
         K = N
   10    CONTINUE
*
*        IF K < 1, EXIT FROM LOOP
*
         IF( K.LT.1 )
     $      GO TO 40
*
         IF( K.GT.NB ) THEN
*
*           FACTORIZE COLUMNS K-KB+1:K OF A AND USE BLOCKED CODE TO
*           UPDATE COLUMNS 1:K-KB
*
            CALL DLASYF( UPLO, K, NB, KB, A, LDA, IPIV, WORK, LDWORK,
     $                   IINFO )
         ELSE
*
*           USE UNBLOCKED CODE TO FACTORIZE COLUMNS 1:K OF A
*
            CALL DSYTF2( UPLO, K, A, LDA, IPIV, IINFO )
            KB = K
         END IF
*
*        SET INFO ON THE FIRST OCCURRENCE OF A ZERO PIVOT
*
         IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO
*
*        DECREASE K AND RETURN TO THE START OF THE MAIN LOOP
*
         K = K - KB
         GO TO 10
*
      ELSE
*
*        FACTORIZE A AS L*D*L' USING THE LOWER TRIANGLE OF A
*
*        K IS THE MAIN LOOP INDEX, INCREASING FROM 1 TO N IN STEPS OF
*        KB, WHERE KB IS THE NUMBER OF COLUMNS FACTORIZED BY DLASYF;
*        KB IS EITHER NB OR NB-1, OR N-K+1 FOR THE LAST BLOCK
*
         K = 1
   20    CONTINUE
*
*        IF K > N, EXIT FROM LOOP
*
         IF( K.GT.N )
     $      GO TO 40
*
         IF( K.LE.N-NB ) THEN
*
*           FACTORIZE COLUMNS K:K+KB-1 OF A AND USE BLOCKED CODE TO
*           UPDATE COLUMNS K+KB:N
*
            CALL DLASYF( UPLO, N-K+1, NB, KB, A( K, K ), LDA, IPIV( K ),
     $                   WORK, LDWORK, IINFO )
         ELSE
*
*           USE UNBLOCKED CODE TO FACTORIZE COLUMNS K:N OF A
*
            CALL DSYTF2( UPLO, N-K+1, A( K, K ), LDA, IPIV( K ), IINFO )
            KB = N - K + 1
         END IF
*
*        SET INFO ON THE FIRST OCCURRENCE OF A ZERO PIVOT
*
         IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $      INFO = IINFO + K - 1
*
*        ADJUST IPIV
*
         DO 30 J = K, K + KB - 1
            IF( IPIV( J ).GT.0 ) THEN
               IPIV( J ) = IPIV( J ) + K - 1
            ELSE
               IPIV( J ) = IPIV( J ) - K + 1
            END IF
   30    CONTINUE
*
*        INCREASE K AND RETURN TO THE START OF THE MAIN LOOP
*
         K = K + KB
         GO TO 20
*
      END IF
*
   40 CONTINUE
      WORK( 1 ) = IWS
      RETURN
*
*     END OF DSYTRF
*
      END
CUT HERE..........
CAT > DSYTF2.F <<'CUT HERE..........'
      SUBROUTINE DSYTF2( UPLO, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.0B) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DSYTF2 COMPUTES THE FACTORIZATION OF A REAL SYMMETRIC MATRIX A USING
*  THE BUNCH-KAUFMAN DIAGONAL PIVOTING METHOD:
*
*     A = U*D*U'  OR  A = L*D*L'
*
*  WHERE U (OR L) IS A PRODUCT OF PERMUTATION AND UNIT UPPER (LOWER)
*  TRIANGULAR MATRICES, U' IS THE TRANSPOSE OF U, AND D IS SYMMETRIC AND
*  BLOCK DIAGONAL WITH 1-BY-1 AND 2-BY-2 DIAGONAL BLOCKS.
*
*  THIS IS THE UNBLOCKED VERSION OF THE ALGORITHM, CALLING LEVEL 2 BLAS.
*
*  ARGUMENTS
*  =========
*
*  UPLO    (INPUT) CHARACTER*1
*          SPECIFIES WHETHER THE UPPER OR LOWER TRIANGULAR PART OF THE
*          SYMMETRIC MATRIX A IS STORED:
*          = 'U':  UPPER TRIANGULAR
*          = 'L':  LOWER TRIANGULAR
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE SYMMETRIC MATRIX A.  IF UPLO = 'U', THE LEADING
*          N-BY-N UPPER TRIANGULAR PART OF A CONTAINS THE UPPER
*          TRIANGULAR PART OF THE MATRIX A, AND THE STRICTLY LOWER
*          TRIANGULAR PART OF A IS NOT REFERENCED.  IF UPLO = 'L', THE
*          LEADING N-BY-N LOWER TRIANGULAR PART OF A CONTAINS THE LOWER
*          TRIANGULAR PART OF THE MATRIX A, AND THE STRICTLY UPPER
*          TRIANGULAR PART OF A IS NOT REFERENCED.
*
*          ON EXIT, THE BLOCK DIAGONAL MATRIX D AND THE MULTIPLIERS USED
*          TO OBTAIN THE FACTOR U OR L (SEE BELOW FOR FURTHER DETAILS).
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,N).
*
*  IPIV    (OUTPUT) INTEGER ARRAY, DIMENSION (N)
*          DETAILS OF THE INTERCHANGES AND THE BLOCK STRUCTURE OF D.
*          IF IPIV(K) > 0, THEN ROWS AND COLUMNS K AND IPIV(K) WERE
*          INTERCHANGED AND D(K,K) IS A 1-BY-1 DIAGONAL BLOCK.
*          IF UPLO = 'U' AND IPIV(K) = IPIV(K-1) < 0, THEN ROWS AND
*          COLUMNS K-1 AND -IPIV(K) WERE INTERCHANGED AND D(K-1:K,K-1:K)
*          IS A 2-BY-2 DIAGONAL BLOCK.  IF UPLO = 'L' AND IPIV(K) =
*          IPIV(K+1) < 0, THEN ROWS AND COLUMNS K+1 AND -IPIV(K) WERE
*          INTERCHANGED AND D(K:K+1,K:K+1) IS A 2-BY-2 DIAGONAL BLOCK.
*
*  INFO    (OUTPUT) INTEGER
*          = 0: SUCCESSFUL EXIT
*          < 0: IF INFO = -K, THE K-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0: IF INFO = K, D(K,K) IS EXACTLY ZERO.  THE FACTORIZATION
*               HAS BEEN COMPLETED, BUT THE BLOCK DIAGONAL MATRIX D IS
*               EXACTLY SINGULAR, AND DIVISION BY ZERO WILL OCCUR IF IT
*               IS USED TO SOLVE A SYSTEM OF EQUATIONS.
*
*  FURTHER DETAILS
*  ===============
*
*  IF UPLO = 'U', THEN A = U*D*U', WHERE
*     U = P(N)*U(N)* ... *P(K)U(K)* ...,
*  I.E., U IS A PRODUCT OF TERMS P(K)*U(K), WHERE K DECREASES FROM N TO
*  1 IN STEPS OF 1 OR 2, AND D IS A BLOCK DIAGONAL MATRIX WITH 1-BY-1
*  AND 2-BY-2 DIAGONAL BLOCKS D(K).  P(K) IS A PERMUTATION MATRIX AS
*  DEFINED BY IPIV(K), AND U(K) IS A UNIT UPPER TRIANGULAR MATRIX, SUCH
*  THAT IF THE DIAGONAL BLOCK D(K) IS OF ORDER S (S = 1 OR 2), THEN
*
*             (   I    V    0   )   K-S
*     U(K) =  (   0    I    0   )   S
*             (   0    0    I   )   N-K
*                K-S   S   N-K
*
*  IF S = 1, D(K) OVERWRITES A(K,K), AND V OVERWRITES A(1:K-1,K).
*  IF S = 2, THE UPPER TRIANGLE OF D(K) OVERWRITES A(K-1,K-1), A(K-1,K),
*  AND A(K,K), AND V OVERWRITES A(1:K-2,K-1:K).
*
*  IF UPLO = 'L', THEN A = L*D*L', WHERE
*     L = P(1)*L(1)* ... *P(K)*L(K)* ...,
*  I.E., L IS A PRODUCT OF TERMS P(K)*L(K), WHERE K INCREASES FROM 1 TO
*  N IN STEPS OF 1 OR 2, AND D IS A BLOCK DIAGONAL MATRIX WITH 1-BY-1
*  AND 2-BY-2 DIAGONAL BLOCKS D(K).  P(K) IS A PERMUTATION MATRIX AS
*  DEFINED BY IPIV(K), AND L(K) IS A UNIT LOWER TRIANGULAR MATRIX, SUCH
*  THAT IF THE DIAGONAL BLOCK D(K) IS OF ORDER S (S = 1 OR 2), THEN
*
*             (   I    0     0   )  K-1
*     L(K) =  (   0    I     0   )  S
*             (   0    V     I   )  N-K-S+1
*                K-1   S  N-K-S+1
*
*  IF S = 1, D(K) OVERWRITES A(K,K), AND V OVERWRITES A(K+1:N,K).
*  IF S = 2, THE LOWER TRIANGLE OF D(K) OVERWRITES A(K,K), A(K+1,K),
*  AND A(K+1,K+1), AND V OVERWRITES A(K+2:N,K:K+1).
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            UPPER
      INTEGER            IMAX, JMAX, K, KK, KP, KSTEP
      DOUBLE PRECISION   ABSAKK, ALPHA, C, COLMAX, R1, R2, ROWMAX, S, T
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      EXTERNAL           LSAME, IDAMAX
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLAEV2, DROT, DSCAL, DSWAP, DSYR, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTF2', -INFO )
         RETURN
      END IF
*
*     INITIALIZE ALPHA FOR USE IN CHOOSING PIVOT BLOCK SIZE.
*
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
*
      IF( UPPER ) THEN
*
*        FACTORIZE A AS U*D*U' USING THE UPPER TRIANGLE OF A
*
*        K IS THE MAIN LOOP INDEX, DECREASING FROM N TO 1 IN STEPS OF
*        1 OR 2
*
         K = N
   10    CONTINUE
*
*        IF K < 1, EXIT FROM LOOP
*
         IF( K.LT.1 )
     $      GO TO 30
         KSTEP = 1
*
*        DETERMINE ROWS AND COLUMNS TO BE INTERCHANGED AND WHETHER
*        A 1-BY-1 OR 2-BY-2 PIVOT BLOCK WILL BE USED
*
         ABSAKK = ABS( A( K, K ) )
*
*        IMAX IS THE ROW-INDEX OF THE LARGEST OFF-DIAGONAL ELEMENT IN
*        COLUMN K, AND COLMAX IS ITS ABSOLUTE VALUE
*
         IF( K.GT.1 ) THEN
            IMAX = IDAMAX( K-1, A( 1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*
*           COLUMN K IS ZERO: SET INFO AND CONTINUE
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              NO INTERCHANGE, USE 1-BY-1 PIVOT BLOCK
*
               KP = K
            ELSE
*
*              JMAX IS THE COLUMN-INDEX OF THE LARGEST OFF-DIAGONAL
*              ELEMENT IN ROW IMAX, AND ROWMAX IS ITS ABSOLUTE VALUE
*
               JMAX = IMAX + IDAMAX( K-IMAX, A( IMAX, IMAX+1 ), LDA )
               ROWMAX = ABS( A( IMAX, JMAX ) )
               IF( IMAX.GT.1 ) THEN
                  JMAX = IDAMAX( IMAX-1, A( 1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( A( JMAX, IMAX ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 NO INTERCHANGE, USE 1-BY-1 PIVOT BLOCK
*
                  KP = K
               ELSE IF( ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 INTERCHANGE ROWS AND COLUMNS K AND IMAX, USE 1-BY-1
*                 PIVOT BLOCK
*
                  KP = IMAX
               ELSE
*
*                 INTERCHANGE ROWS AND COLUMNS K-1 AND IMAX, USE 2-BY-2
*                 PIVOT BLOCK
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
            KK = K - KSTEP + 1
            IF( KP.NE.KK ) THEN
*
*              INTERCHANGE ROWS AND COLUMNS KK AND KP IN THE LEADING
*              SUBMATRIX A(1:K,1:K)
*
               CALL DSWAP( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
               CALL DSWAP( KK-KP-1, A( KP+1, KK ), 1, A( KP, KP+1 ),
     $                     LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = A( K-1, K )
                  A( K-1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF
            END IF
*
*           UPDATE THE LEADING SUBMATRIX
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-BY-1 PIVOT BLOCK D(K): COLUMN K NOW HOLDS
*
*              W(K) = U(K)*D(K)
*
*              WHERE U(K) IS THE K-TH COLUMN OF U
*
*              PERFORM A RANK-1 UPDATE OF A(1:K-1,1:K-1) AS
*
*              A := A - U(K)*D(K)*U(K)' = A - W(K)*1/D(K)*W(K)'
*
               R1 = ONE / A( K, K )
               CALL DSYR( UPLO, K-1, -R1, A( 1, K ), 1, A, LDA )
*
*              STORE U(K) IN COLUMN K
*
               CALL DSCAL( K-1, R1, A( 1, K ), 1 )
            ELSE
*
*              2-BY-2 PIVOT BLOCK D(K): COLUMNS K AND K-1 NOW HOLD
*
*              ( W(K-1) W(K) ) = ( U(K-1) U(K) )*D(K)
*
*              WHERE U(K) AND U(K-1) ARE THE K-TH AND (K-1)-TH COLUMNS
*              OF U
*
*              PERFORM A RANK-2 UPDATE OF A(1:K-2,1:K-2) AS
*
*              A := A - ( U(K-1) U(K) )*D(K)*( U(K-1) U(K) )'
*                 = A - ( W(K-1) W(K) )*INV(D(K))*( W(K-1) W(K) )'
*
*              CONVERT THIS TO TWO RANK-1 UPDATES BY USING THE EIGEN-
*              DECOMPOSITION OF D(K)
*
               CALL DLAEV2( A( K-1, K-1 ), A( K-1, K ), A( K, K ), R1,
     $                      R2, C, S )
               R1 = ONE / R1
               R2 = ONE / R2
               CALL DROT( K-2, A( 1, K-1 ), 1, A( 1, K ), 1, C, S )
               CALL DSYR( UPLO, K-2, -R1, A( 1, K-1 ), 1, A, LDA )
               CALL DSYR( UPLO, K-2, -R2, A( 1, K ), 1, A, LDA )
*
*              STORE U(K) AND U(K-1) IN COLUMNS K AND K-1
*
               CALL DSCAL( K-2, R1, A( 1, K-1 ), 1 )
               CALL DSCAL( K-2, R2, A( 1, K ), 1 )
               CALL DROT( K-2, A( 1, K-1 ), 1, A( 1, K ), 1, C, -S )
            END IF
         END IF
*
*        STORE DETAILS OF THE INTERCHANGES IN IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF
*
*        DECREASE K AND RETURN TO THE START OF THE MAIN LOOP
*
         K = K - KSTEP
         GO TO 10
*
      ELSE
*
*        FACTORIZE A AS L*D*L' USING THE LOWER TRIANGLE OF A
*
*        K IS THE MAIN LOOP INDEX, INCREASING FROM 1 TO N IN STEPS OF
*        1 OR 2
*
         K = 1
   20    CONTINUE
*
*        IF K > N, EXIT FROM LOOP
*
         IF( K.GT.N )
     $      GO TO 30
         KSTEP = 1
*
*        DETERMINE ROWS AND COLUMNS TO BE INTERCHANGED AND WHETHER
*        A 1-BY-1 OR 2-BY-2 PIVOT BLOCK WILL BE USED
*
         ABSAKK = ABS( A( K, K ) )
*
*        IMAX IS THE ROW-INDEX OF THE LARGEST OFF-DIAGONAL ELEMENT IN
*        COLUMN K, AND COLMAX IS ITS ABSOLUTE VALUE
*
         IF( K.LT.N ) THEN
            IMAX = K + IDAMAX( N-K, A( K+1, K ), 1 )
            COLMAX = ABS( A( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*
*           COLUMN K IS ZERO: SET INFO AND CONTINUE
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              NO INTERCHANGE, USE 1-BY-1 PIVOT BLOCK
*
               KP = K
            ELSE
*
*              JMAX IS THE COLUMN-INDEX OF THE LARGEST OFF-DIAGONAL
*              ELEMENT IN ROW IMAX, AND ROWMAX IS ITS ABSOLUTE VALUE
*
               JMAX = K - 1 + IDAMAX( IMAX-K, A( IMAX, K ), LDA )
               ROWMAX = ABS( A( IMAX, JMAX ) )
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + IDAMAX( N-IMAX, A( IMAX+1, IMAX ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( A( JMAX, IMAX ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 NO INTERCHANGE, USE 1-BY-1 PIVOT BLOCK
*
                  KP = K
               ELSE IF( ABS( A( IMAX, IMAX ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 INTERCHANGE ROWS AND COLUMNS K AND IMAX, USE 1-BY-1
*                 PIVOT BLOCK
*
                  KP = IMAX
               ELSE
*
*                 INTERCHANGE ROWS AND COLUMNS K+1 AND IMAX, USE 2-BY-2
*                 PIVOT BLOCK
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
            KK = K + KSTEP - 1
            IF( KP.NE.KK ) THEN
*
*              INTERCHANGE ROWS AND COLUMNS KK AND KP IN THE TRAILING
*              SUBMATRIX A(K:N,K:N)
*
               IF( KP.LT.N )
     $            CALL DSWAP( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
               CALL DSWAP( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ),
     $                     LDA )
               T = A( KK, KK )
               A( KK, KK ) = A( KP, KP )
               A( KP, KP ) = T
               IF( KSTEP.EQ.2 ) THEN
                  T = A( K+1, K )
                  A( K+1, K ) = A( KP, K )
                  A( KP, K ) = T
               END IF
            END IF
*
*           UPDATE THE TRAILING SUBMATRIX
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-BY-1 PIVOT BLOCK D(K): COLUMN K NOW HOLDS
*
*              W(K) = L(K)*D(K)
*
*              WHERE L(K) IS THE K-TH COLUMN OF L
*
               IF( K.LT.N ) THEN
*
*                 PERFORM A RANK-1 UPDATE OF A(K+1:N,K+1:N) AS
*
*                 A := A - L(K)*D(K)*L(K)' = A - W(K)*(1/D(K))*W(K)'
*
                  R1 = ONE / A( K, K )
                  CALL DSYR( UPLO, N-K, -R1, A( K+1, K ), 1,
     $                       A( K+1, K+1 ), LDA )
*
*                 STORE L(K) IN COLUMN K
*
                  CALL DSCAL( N-K, R1, A( K+1, K ), 1 )
               END IF
            ELSE
*
*              2-BY-2 PIVOT BLOCK D(K): COLUMNS K AND K+1 NOW HOLD
*
*              ( W(K) W(K+1) ) = ( L(K) L(K+1) )*D(K)
*
*              WHERE L(K) AND L(K+1) ARE THE K-TH AND (K+1)-TH COLUMNS
*              OF L
*
               IF( K.LT.N-1 ) THEN
*
*                 PERFORM A RANK-2 UPDATE OF A(K+2:N,K+2:N) AS
*
*                 A := A - ( L(K) L(K+1) )*D(K)*( L(K) L(K+1) )'
*                    = A - ( W(K) W(K+1) )*INV(D(K))*( W(K) W(K+1) )'
*
*                 CONVERT THIS TO TWO RANK-1 UPDATES BY USING THE EIGEN-
*                 DECOMPOSITION OF D(K)
*
                  CALL DLAEV2( A( K, K ), A( K+1, K ), A( K+1, K+1 ),
     $                         R1, R2, C, S )
                  R1 = ONE / R1
                  R2 = ONE / R2
                  CALL DROT( N-K-1, A( K+2, K ), 1, A( K+2, K+1 ), 1, C,
     $                       S )
                  CALL DSYR( UPLO, N-K-1, -R1, A( K+2, K ), 1,
     $                       A( K+2, K+2 ), LDA )
                  CALL DSYR( UPLO, N-K-1, -R2, A( K+2, K+1 ), 1,
     $                       A( K+2, K+2 ), LDA )
*
*                 STORE L(K) AND L(K+1) IN COLUMNS K AND K+1
*
                  CALL DSCAL( N-K-1, R1, A( K+2, K ), 1 )
                  CALL DSCAL( N-K-1, R2, A( K+2, K+1 ), 1 )
                  CALL DROT( N-K-1, A( K+2, K ), 1, A( K+2, K+1 ), 1, C,
     $                       -S )
               END IF
            END IF
         END IF
*
*        STORE DETAILS OF THE INTERCHANGES IN IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF
*
*        INCREASE K AND RETURN TO THE START OF THE MAIN LOOP
*
         K = K + KSTEP
         GO TO 20
*
      END IF
*
   30 CONTINUE
      RETURN
*
*     END OF DSYTF2
*
      END
CUT HERE..........
CAT > DLASYF.F <<'CUT HERE..........'
      SUBROUTINE DLASYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.0B) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     FEBRUARY 29, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          UPLO
      INTEGER            INFO, KB, LDA, LDW, N, NB
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), W( LDW, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLASYF COMPUTES A PARTIAL FACTORIZATION OF A REAL SYMMETRIC MATRIX A
*  USING THE BUNCH-KAUFMAN DIAGONAL PIVOTING METHOD. THE PARTIAL
*  FACTORIZATION HAS THE FORM:
*
*  A  =  ( I  U12 ) ( A11  0  ) (  I    0   )  IF UPLO = 'U', OR:
*        ( 0  U22 ) (  0   D  ) ( U12' U22' )
*
*  A  =  ( L11  0 ) (  D   0  ) ( L11' L21' )  IF UPLO = 'L'
*        ( L21  I ) (  0  A22 ) (  0    I   )
*
*  WHERE THE ORDER OF D IS AT MOST NB. THE ACTUAL ORDER IS RETURNED IN
*  THE ARGUMENT KB, AND IS EITHER NB OR NB-1, OR N IF N <= NB.
*
*  DLASYF IS AN AUXILIARY ROUTINE CALLED BY DSYTRF. IT USES BLOCKED CODE
*  (CALLING LEVEL 3 BLAS) TO UPDATE THE SUBMATRIX A11 (IF UPLO = 'U') OR
*  A22 (IF UPLO = 'L').
*
*  ARGUMENTS
*  =========
*
*  UPLO    (INPUT) CHARACTER*1
*          SPECIFIES WHETHER THE UPPER OR LOWER TRIANGULAR PART OF THE
*          SYMMETRIC MATRIX A IS STORED:
*          = 'U':  UPPER TRIANGULAR
*          = 'L':  LOWER TRIANGULAR
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX A.  N >= 0.
*
*  NB      (INPUT) INTEGER
*          THE MAXIMUM NUMBER OF COLUMNS OF THE MATRIX A THAT SHOULD BE
*          FACTORED.  NB SHOULD BE AT LEAST 2 TO ALLOW FOR 2-BY-2 PIVOT
*          BLOCKS.
*
*  KB      (OUTPUT) INTEGER
*          THE NUMBER OF COLUMNS OF A THAT WERE ACTUALLY FACTORED.
*          KB IS EITHER NB-1 OR NB, OR N IF N <= NB.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE SYMMETRIC MATRIX A.  IF UPLO = 'U', THE LEADING
*          N-BY-N UPPER TRIANGULAR PART OF A CONTAINS THE UPPER
*          TRIANGULAR PART OF THE MATRIX A, AND THE STRICTLY LOWER
*          TRIANGULAR PART OF A IS NOT REFERENCED.  IF UPLO = 'L', THE
*          LEADING N-BY-N LOWER TRIANGULAR PART OF A CONTAINS THE LOWER
*          TRIANGULAR PART OF THE MATRIX A, AND THE STRICTLY UPPER
*          TRIANGULAR PART OF A IS NOT REFERENCED.
*          ON EXIT, A CONTAINS DETAILS OF THE PARTIAL FACTORIZATION.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,N).
*
*  IPIV    (OUTPUT) INTEGER ARRAY, DIMENSION (N)
*          DETAILS OF THE INTERCHANGES AND THE BLOCK STRUCTURE OF D.
*          IF UPLO = 'U', ONLY THE LAST KB ELEMENTS OF IPIV ARE SET;
*          IF UPLO = 'L', ONLY THE FIRST KB ELEMENTS ARE SET.
*
*          IF IPIV(K) > 0, THEN ROWS AND COLUMNS K AND IPIV(K) WERE
*          INTERCHANGED AND D(K,K) IS A 1-BY-1 DIAGONAL BLOCK.
*          IF UPLO = 'U' AND IPIV(K) = IPIV(K-1) < 0, THEN ROWS AND
*          COLUMNS K-1 AND -IPIV(K) WERE INTERCHANGED AND D(K-1:K,K-1:K)
*          IS A 2-BY-2 DIAGONAL BLOCK.  IF UPLO = 'L' AND IPIV(K) =
*          IPIV(K+1) < 0, THEN ROWS AND COLUMNS K+1 AND -IPIV(K) WERE
*          INTERCHANGED AND D(K:K+1,K:K+1) IS A 2-BY-2 DIAGONAL BLOCK.
*
*  W       (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (LDW,NB)
*
*  LDW     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY W.  LDW >= MAX(1,N).
*
*  INFO    (OUTPUT) INTEGER
*          = 0: SUCCESSFUL EXIT
*          > 0: IF INFO = K, D(K,K) IS EXACTLY ZERO.  THE FACTORIZATION
*               HAS BEEN COMPLETED, BUT THE BLOCK DIAGONAL MATRIX D IS
*               EXACTLY SINGULAR.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            IMAX, J, JB, JJ, JMAX, JP, K, KK, KKW, KP,
     $                   KSTEP, KW
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, D11, D21, D22, R1,
     $                   ROWMAX, T
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      EXTERNAL           LSAME, IDAMAX
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DCOPY, DGEMM, DGEMV, DSCAL, DSWAP
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
      INFO = 0
*
*     INITIALIZE ALPHA FOR USE IN CHOOSING PIVOT BLOCK SIZE.
*
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        FACTORIZE THE TRAILING COLUMNS OF A USING THE UPPER TRIANGLE
*        OF A AND WORKING BACKWARDS, AND COMPUTE THE MATRIX W = U12*D
*        FOR USE IN UPDATING A11
*
*        K IS THE MAIN LOOP INDEX, DECREASING FROM N IN STEPS OF 1 OR 2
*
*        KW IS THE COLUMN OF W WHICH CORRESPONDS TO COLUMN K OF A
*
         K = N
   10    CONTINUE
         KW = NB + K - N
*
*        EXIT FROM LOOP
*
         IF( ( K.LE.N-NB+1 .AND. NB.LT.N ) .OR. K.LT.1 )
     $      GO TO 30
*
*        COPY COLUMN K OF A TO COLUMN KW OF W AND UPDATE IT
*
         CALL DCOPY( K, A( 1, K ), 1, W( 1, KW ), 1 )
         IF( K.LT.N )
     $      CALL DGEMV( 'NO TRANSPOSE', K, N-K, -ONE, A( 1, K+1 ), LDA,
     $                  W( K, KW+1 ), LDW, ONE, W( 1, KW ), 1 )
*
         KSTEP = 1
*
*        DETERMINE ROWS AND COLUMNS TO BE INTERCHANGED AND WHETHER
*        A 1-BY-1 OR 2-BY-2 PIVOT BLOCK WILL BE USED
*
         ABSAKK = ABS( W( K, KW ) )
*
*        IMAX IS THE ROW-INDEX OF THE LARGEST OFF-DIAGONAL ELEMENT IN
*        COLUMN K, AND COLMAX IS ITS ABSOLUTE VALUE
*
         IF( K.GT.1 ) THEN
            IMAX = IDAMAX( K-1, W( 1, KW ), 1 )
            COLMAX = ABS( W( IMAX, KW ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*
*           COLUMN K IS ZERO: SET INFO AND CONTINUE
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              NO INTERCHANGE, USE 1-BY-1 PIVOT BLOCK
*
               KP = K
            ELSE
*
*              COPY COLUMN IMAX TO COLUMN KW-1 OF W AND UPDATE IT
*
               CALL DCOPY( IMAX, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 )
               CALL DCOPY( K-IMAX, A( IMAX, IMAX+1 ), LDA,
     $                     W( IMAX+1, KW-1 ), 1 )
               IF( K.LT.N )
     $            CALL DGEMV( 'NO TRANSPOSE', K, N-K, -ONE, A( 1, K+1 ),
     $                        LDA, W( IMAX, KW+1 ), LDW, ONE,
     $                        W( 1, KW-1 ), 1 )
*
*              JMAX IS THE COLUMN-INDEX OF THE LARGEST OFF-DIAGONAL
*              ELEMENT IN ROW IMAX, AND ROWMAX IS ITS ABSOLUTE VALUE
*
               JMAX = IMAX + IDAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 )
               ROWMAX = ABS( W( JMAX, KW-1 ) )
               IF( IMAX.GT.1 ) THEN
                  JMAX = IDAMAX( IMAX-1, W( 1, KW-1 ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( W( JMAX, KW-1 ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 NO INTERCHANGE, USE 1-BY-1 PIVOT BLOCK
*
                  KP = K
               ELSE IF( ABS( W( IMAX, KW-1 ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 INTERCHANGE ROWS AND COLUMNS K AND IMAX, USE 1-BY-1
*                 PIVOT BLOCK
*
                  KP = IMAX
*
*                 COPY COLUMN KW-1 OF W TO COLUMN KW
*
                  CALL DCOPY( K, W( 1, KW-1 ), 1, W( 1, KW ), 1 )
               ELSE
*
*                 INTERCHANGE ROWS AND COLUMNS K-1 AND IMAX, USE 2-BY-2
*                 PIVOT BLOCK
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
            KK = K - KSTEP + 1
            KKW = NB + KK - N
*
*           UPDATED COLUMN KP IS ALREADY STORED IN COLUMN KKW OF W
*
            IF( KP.NE.KK ) THEN
*
*              COPY NON-UPDATED COLUMN KK TO COLUMN KP
*
               A( KP, K ) = A( KK, K )
               CALL DCOPY( K-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ),
     $                     LDA )
               CALL DCOPY( KP, A( 1, KK ), 1, A( 1, KP ), 1 )
*
*              INTERCHANGE ROWS KK AND KP IN LAST KK COLUMNS OF A AND W
*
               CALL DSWAP( N-KK+1, A( KK, KK ), LDA, A( KP, KK ), LDA )
               CALL DSWAP( N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ),
     $                     LDW )
            END IF
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-BY-1 PIVOT BLOCK D(K): COLUMN KW OF W NOW HOLDS
*
*              W(K) = U(K)*D(K)
*
*              WHERE U(K) IS THE K-TH COLUMN OF U
*
*              STORE U(K) IN COLUMN K OF A
*
               CALL DCOPY( K, W( 1, KW ), 1, A( 1, K ), 1 )
               R1 = ONE / A( K, K )
               CALL DSCAL( K-1, R1, A( 1, K ), 1 )
            ELSE
*
*              2-BY-2 PIVOT BLOCK D(K): COLUMNS KW AND KW-1 OF W NOW
*              HOLD
*
*              ( W(K-1) W(K) ) = ( U(K-1) U(K) )*D(K)
*
*              WHERE U(K) AND U(K-1) ARE THE K-TH AND (K-1)-TH COLUMNS
*              OF U
*
               IF( K.GT.2 ) THEN
*
*                 STORE U(K) AND U(K-1) IN COLUMNS K AND K-1 OF A
*
                  D21 = W( K-1, KW )
                  D11 = W( K, KW ) / D21
                  D22 = W( K-1, KW-1 ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
                  DO 20 J = 1, K - 2
                     A( J, K-1 ) = D21*( D11*W( J, KW-1 )-W( J, KW ) )
                     A( J, K ) = D21*( D22*W( J, KW )-W( J, KW-1 ) )
   20             CONTINUE
               END IF
*
*              COPY D(K) TO A
*
               A( K-1, K-1 ) = W( K-1, KW-1 )
               A( K-1, K ) = W( K-1, KW )
               A( K, K ) = W( K, KW )
            END IF
         END IF
*
*        STORE DETAILS OF THE INTERCHANGES IN IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF
*
*        DECREASE K AND RETURN TO THE START OF THE MAIN LOOP
*
         K = K - KSTEP
         GO TO 10
*
   30    CONTINUE
*
*        UPDATE THE UPPER TRIANGLE OF A11 (= A(1:K,1:K)) AS
*
*        A11 := A11 - U12*D*U12' = A11 - U12*W'
*
*        COMPUTING BLOCKS OF NB COLUMNS AT A TIME
*
         DO 50 J = ( ( K-1 ) / NB )*NB + 1, 1, -NB
            JB = MIN( NB, K-J+1 )
*
*           UPDATE THE UPPER TRIANGLE OF THE DIAGONAL BLOCK
*
            DO 40 JJ = J, J + JB - 1
               CALL DGEMV( 'NO TRANSPOSE', JJ-J+1, N-K, -ONE,
     $                     A( J, K+1 ), LDA, W( JJ, KW+1 ), LDW, ONE,
     $                     A( J, JJ ), 1 )
   40       CONTINUE
*
*           UPDATE THE RECTANGULAR SUPERDIAGONAL BLOCK
*
            CALL DGEMM( 'NO TRANSPOSE', 'TRANSPOSE', J-1, JB, N-K, -ONE,
     $                  A( 1, K+1 ), LDA, W( J, KW+1 ), LDW, ONE,
     $                  A( 1, J ), LDA )
   50    CONTINUE
*
*        PUT U12 IN STANDARD FORM BY PARTIALLY UNDOING THE INTERCHANGES
*        IN COLUMNS K+1:N
*
         J = K + 1
   60    CONTINUE
         JJ = J
         JP = IPIV( J )
         IF( JP.LT.0 ) THEN
            JP = -JP
            J = J + 1
         END IF
         J = J + 1
         IF( JP.NE.JJ .AND. J.LE.N )
     $      CALL DSWAP( N-J+1, A( JP, J ), LDA, A( JJ, J ), LDA )
         IF( J.LE.N )
     $      GO TO 60
*
*        SET KB TO THE NUMBER OF COLUMNS FACTORIZED
*
         KB = N - K
*
      ELSE
*
*        FACTORIZE THE LEADING COLUMNS OF A USING THE LOWER TRIANGLE
*        OF A AND WORKING FORWARDS, AND COMPUTE THE MATRIX W = L21*D
*        FOR USE IN UPDATING A22
*
*        K IS THE MAIN LOOP INDEX, INCREASING FROM 1 IN STEPS OF 1 OR 2
*
         K = 1
   70    CONTINUE
*
*        EXIT FROM LOOP
*
         IF( ( K.GE.NB .AND. NB.LT.N ) .OR. K.GT.N )
     $      GO TO 90
*
*        COPY COLUMN K OF A TO COLUMN K OF W AND UPDATE IT
*
         CALL DCOPY( N-K+1, A( K, K ), 1, W( K, K ), 1 )
         CALL DGEMV( 'NO TRANSPOSE', N-K+1, K-1, -ONE, A( K, 1 ), LDA,
     $               W( K, 1 ), LDW, ONE, W( K, K ), 1 )
*
         KSTEP = 1
*
*        DETERMINE ROWS AND COLUMNS TO BE INTERCHANGED AND WHETHER
*        A 1-BY-1 OR 2-BY-2 PIVOT BLOCK WILL BE USED
*
         ABSAKK = ABS( W( K, K ) )
*
*        IMAX IS THE ROW-INDEX OF THE LARGEST OFF-DIAGONAL ELEMENT IN
*        COLUMN K, AND COLMAX IS ITS ABSOLUTE VALUE
*
         IF( K.LT.N ) THEN
            IMAX = K + IDAMAX( N-K, W( K+1, K ), 1 )
            COLMAX = ABS( W( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
*
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
*
*           COLUMN K IS ZERO: SET INFO AND CONTINUE
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
*
*              NO INTERCHANGE, USE 1-BY-1 PIVOT BLOCK
*
               KP = K
            ELSE
*
*              COPY COLUMN IMAX TO COLUMN K+1 OF W AND UPDATE IT
*
               CALL DCOPY( IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1 )
               CALL DCOPY( N-IMAX+1, A( IMAX, IMAX ), 1, W( IMAX, K+1 ),
     $                     1 )
               CALL DGEMV( 'NO TRANSPOSE', N-K+1, K-1, -ONE, A( K, 1 ),
     $                     LDA, W( IMAX, 1 ), LDW, ONE, W( K, K+1 ), 1 )
*
*              JMAX IS THE COLUMN-INDEX OF THE LARGEST OFF-DIAGONAL
*              ELEMENT IN ROW IMAX, AND ROWMAX IS ITS ABSOLUTE VALUE
*
               JMAX = K - 1 + IDAMAX( IMAX-K, W( K, K+1 ), 1 )
               ROWMAX = ABS( W( JMAX, K+1 ) )
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + IDAMAX( N-IMAX, W( IMAX+1, K+1 ), 1 )
                  ROWMAX = MAX( ROWMAX, ABS( W( JMAX, K+1 ) ) )
               END IF
*
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
*
*                 NO INTERCHANGE, USE 1-BY-1 PIVOT BLOCK
*
                  KP = K
               ELSE IF( ABS( W( IMAX, K+1 ) ).GE.ALPHA*ROWMAX ) THEN
*
*                 INTERCHANGE ROWS AND COLUMNS K AND IMAX, USE 1-BY-1
*                 PIVOT BLOCK
*
                  KP = IMAX
*
*                 COPY COLUMN K+1 OF W TO COLUMN K
*
                  CALL DCOPY( N-K+1, W( K, K+1 ), 1, W( K, K ), 1 )
               ELSE
*
*                 INTERCHANGE ROWS AND COLUMNS K+1 AND IMAX, USE 2-BY-2
*                 PIVOT BLOCK
*
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
*
            KK = K + KSTEP - 1
*
*           UPDATED COLUMN KP IS ALREADY STORED IN COLUMN KK OF W
*
            IF( KP.NE.KK ) THEN
*
*              COPY NON-UPDATED COLUMN KK TO COLUMN KP
*
               A( KP, K ) = A( KK, K )
               CALL DCOPY( KP-K-1, A( K+1, KK ), 1, A( KP, K+1 ), LDA )
               CALL DCOPY( N-KP+1, A( KP, KK ), 1, A( KP, KP ), 1 )
*
*              INTERCHANGE ROWS KK AND KP IN FIRST KK COLUMNS OF A AND W
*
               CALL DSWAP( KK, A( KK, 1 ), LDA, A( KP, 1 ), LDA )
               CALL DSWAP( KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW )
            END IF
*
            IF( KSTEP.EQ.1 ) THEN
*
*              1-BY-1 PIVOT BLOCK D(K): COLUMN K OF W NOW HOLDS
*
*              W(K) = L(K)*D(K)
*
*              WHERE L(K) IS THE K-TH COLUMN OF L
*
*              STORE L(K) IN COLUMN K OF A
*
               CALL DCOPY( N-K+1, W( K, K ), 1, A( K, K ), 1 )
               IF( K.LT.N ) THEN
                  R1 = ONE / A( K, K )
                  CALL DSCAL( N-K, R1, A( K+1, K ), 1 )
               END IF
            ELSE
*
*              2-BY-2 PIVOT BLOCK D(K): COLUMNS K AND K+1 OF W NOW HOLD
*
*              ( W(K) W(K+1) ) = ( L(K) L(K+1) )*D(K)
*
*              WHERE L(K) AND L(K+1) ARE THE K-TH AND (K+1)-TH COLUMNS
*              OF L
*
               IF( K.LT.N-1 ) THEN
*
*                 STORE L(K) AND L(K+1) IN COLUMNS K AND K+1 OF A
*
                  D21 = W( K+1, K )
                  D11 = W( K+1, K+1 ) / D21
                  D22 = W( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  D21 = T / D21
                  DO 80 J = K + 2, N
                     A( J, K ) = D21*( D11*W( J, K )-W( J, K+1 ) )
                     A( J, K+1 ) = D21*( D22*W( J, K+1 )-W( J, K ) )
   80             CONTINUE
               END IF
*
*              COPY D(K) TO A
*
               A( K, K ) = W( K, K )
               A( K+1, K ) = W( K+1, K )
               A( K+1, K+1 ) = W( K+1, K+1 )
            END IF
         END IF
*
*        STORE DETAILS OF THE INTERCHANGES IN IPIV
*
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF
*
*        INCREASE K AND RETURN TO THE START OF THE MAIN LOOP
*
         K = K + KSTEP
         GO TO 70
*
   90    CONTINUE
*
*        UPDATE THE LOWER TRIANGLE OF A22 (= A(K:N,K:N)) AS
*
*        A22 := A22 - L21*D*L21' = A22 - L21*W'
*
*        COMPUTING BLOCKS OF NB COLUMNS AT A TIME
*
         DO 110 J = K, N, NB
            JB = MIN( NB, N-J+1 )
*
*           UPDATE THE LOWER TRIANGLE OF THE DIAGONAL BLOCK
*
            DO 100 JJ = J, J + JB - 1
               CALL DGEMV( 'NO TRANSPOSE', J+JB-JJ, K-1, -ONE,
     $                     A( JJ, 1 ), LDA, W( JJ, 1 ), LDW, ONE,
     $                     A( JJ, JJ ), 1 )
  100       CONTINUE
*
*           UPDATE THE RECTANGULAR SUBDIAGONAL BLOCK
*
            IF( J+JB.LE.N )
     $         CALL DGEMM( 'NO TRANSPOSE', 'TRANSPOSE', N-J-JB+1, JB,
     $                     K-1, -ONE, A( J+JB, 1 ), LDA, W( J, 1 ), LDW,
     $                     ONE, A( J+JB, J ), LDA )
  110    CONTINUE
*
*        PUT L21 IN STANDARD FORM BY PARTIALLY UNDOING THE INTERCHANGES
*        IN COLUMNS 1:K-1
*
         J = K - 1
  120    CONTINUE
         JJ = J
         JP = IPIV( J )
         IF( JP.LT.0 ) THEN
            JP = -JP
            J = J - 1
         END IF
         J = J - 1
         IF( JP.NE.JJ .AND. J.GE.1 )
     $      CALL DSWAP( J, A( JP, 1 ), LDA, A( JJ, 1 ), LDA )
         IF( J.GE.1 )
     $      GO TO 120
*
*        SET KB TO THE NUMBER OF COLUMNS FACTORIZED
*
         KB = K - 1
*
      END IF
      RETURN
*
*     END OF DLASYF
*
      END
CUT HERE..........
C    ---------------------  BELOW IS DSYTRI  ------------------
      SUBROUTINE DSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.0B) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  PURPOSE
*  =======
*
*  DSYTRI COMPUTES THE INVERSE OF A REAL SYMMETRIC INDEFINITE MATRIX
*  A USING THE FACTORIZATION A = U*D*U' OR A = L*D*L' COMPUTED BY
*  DSYTRF.
*
*  ARGUMENTS
*  =========
*
*  UPLO    (INPUT) CHARACTER*1
*          SPECIFIES WHETHER THE DETAILS OF THE FACTORIZATION ARE STORED
*          AS AN UPPER OR LOWER TRIANGULAR MATRIX.
*          = 'U':  UPPER TRIANGULAR (FORM IS A = U*D*U')
*          = 'L':  LOWER TRIANGULAR (FORM IS A = L*D*L')
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE BLOCK DIAGONAL MATRIX D AND THE MULTIPLIERS
*          USED TO OBTAIN THE FACTOR U OR L AS COMPUTED BY DSYTRF.
*
*          ON EXIT, IF INFO = 0, THE (SYMMETRIC) INVERSE OF THE ORIGINAL
*          MATRIX.  IF UPLO = 'U', THE UPPER TRIANGULAR PART OF THE
*          INVERSE IS FORMED AND THE PART OF A BELOW THE DIAGONAL IS NOT
*          REFERENCED; IF UPLO = 'L' THE LOWER TRIANGULAR PART OF THE
*          INVERSE IS FORMED AND THE PART OF A ABOVE THE DIAGONAL IS
*          NOT REFERENCED.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,N).
*
*  IPIV    (INPUT) INTEGER ARRAY, DIMENSION (N)
*          DETAILS OF THE INTERCHANGES AND THE BLOCK STRUCTURE OF D
*          AS DETERMINED BY DSYTRF.
*
*  WORK    (WORKSPACE) DOUBLE PRECISION ARRAY, DIMENSION (N)
*
*  INFO    (OUTPUT) INTEGER
*          = 0: SUCCESSFUL EXIT
*          < 0: IF INFO = -K, THE K-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0: IF INFO = K, D(K,K) = 0; THE MATRIX IS SINGULAR AND ITS
*               INVERSE COULD NOT BE COMPUTED.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            UPPER
      INTEGER            K, KP, KSTEP
      DOUBLE PRECISION   AK, AKKP1, AKP1, D, T, TEMP
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DCOPY, DSWAP, DSYMV, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          ABS, MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRI', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.EQ.0 )
     $   RETURN
*
*     CHECK THAT THE DIAGONAL MATRIX D IS NONSINGULAR.
*
      IF( UPPER ) THEN
*
*        UPPER TRIANGULAR STORAGE: EXAMINE D FROM BOTTOM TO TOP
*
         DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
      ELSE
*
*        LOWER TRIANGULAR STORAGE: EXAMINE D FROM TOP TO BOTTOM.
*
         DO 20 INFO = 1, N
            IF( IPIV( INFO ).GT.0 .AND. A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   20    CONTINUE
      END IF
      INFO = 0
*
      IF( UPPER ) THEN
*
*        COMPUTE INV(A) FROM THE FACTORIZATION A = U*D*U'.
*
*        K IS THE MAIN LOOP INDEX, INCREASING FROM 1 TO N IN STEPS OF
*        1 OR 2, DEPENDING ON THE SIZE OF THE DIAGONAL BLOCKS.
*
         K = 1
   30    CONTINUE
*
*        IF K > N, EXIT FROM LOOP.
*
         IF( K.GT.N )
     $      GO TO 40
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 X 1 DIAGONAL BLOCK
*
*           INVERT THE DIAGONAL BLOCK.
*
            A( K, K ) = ONE / A( K, K )
*
*           COMPUTE COLUMN K OF THE INVERSE.
*
            IF( K.GT.1 ) THEN
               CALL DCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO,
     $                     A( 1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( K-1, WORK, 1, A( 1, K ),
     $                     1 )
            END IF
            KSTEP = 1
         ELSE
*
*           2 X 2 DIAGONAL BLOCK
*
*           INVERT THE DIAGONAL BLOCK.
*
            T = ABS( A( K, K+1 ) )
            AK = A( K, K ) / T
            AKP1 = A( K+1, K+1 ) / T
            AKKP1 = A( K, K+1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K, K ) = AKP1 / D
            A( K+1, K+1 ) = AK / D
            A( K, K+1 ) = -AKKP1 / D
*
*           COMPUTE COLUMNS K AND K+1 OF THE INVERSE.
*
            IF( K.GT.1 ) THEN
               CALL DCOPY( K-1, A( 1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO,
     $                     A( 1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( K-1, WORK, 1, A( 1, K ),
     $                     1 )
               A( K, K+1 ) = A( K, K+1 ) -
     $                       DDOT( K-1, A( 1, K ), 1, A( 1, K+1 ), 1 )
               CALL DCOPY( K-1, A( 1, K+1 ), 1, WORK, 1 )
               CALL DSYMV( UPLO, K-1, -ONE, A, LDA, WORK, 1, ZERO,
     $                     A( 1, K+1 ), 1 )
               A( K+1, K+1 ) = A( K+1, K+1 ) -
     $                         DDOT( K-1, WORK, 1, A( 1, K+1 ), 1 )
            END IF
            KSTEP = 2
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           INTERCHANGE ROWS AND COLUMNS K AND KP IN THE LEADING
*           SUBMATRIX A(1:K+1,1:K+1)
*
            CALL DSWAP( KP-1, A( 1, K ), 1, A( 1, KP ), 1 )
            CALL DSWAP( K-KP-1, A( KP+1, K ), 1, A( KP, KP+1 ), LDA )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = A( K, K+1 )
               A( K, K+1 ) = A( KP, K+1 )
               A( KP, K+1 ) = TEMP
            END IF
         END IF
*
         K = K + KSTEP
         GO TO 30
   40    CONTINUE
*
      ELSE
*
*        COMPUTE INV(A) FROM THE FACTORIZATION A = L*D*L'.
*
*        K IS THE MAIN LOOP INDEX, INCREASING FROM 1 TO N IN STEPS OF
*        1 OR 2, DEPENDING ON THE SIZE OF THE DIAGONAL BLOCKS.
*
         K = N
   50    CONTINUE
*
*        IF K < 1, EXIT FROM LOOP.
*
         IF( K.LT.1 )
     $      GO TO 60
*
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 X 1 DIAGONAL BLOCK
*
*           INVERT THE DIAGONAL BLOCK.
*
            A( K, K ) = ONE / A( K, K )
*
*           COMPUTE COLUMN K OF THE INVERSE.
*
            IF( K.LT.N ) THEN
               CALL DCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1,
     $                     ZERO, A( K+1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( N-K, WORK, 1, A( K+1, K ),
     $                     1 )
            END IF
            KSTEP = 1
         ELSE
*
*           2 X 2 DIAGONAL BLOCK
*
*           INVERT THE DIAGONAL BLOCK.
*
            T = ABS( A( K, K-1 ) )
            AK = A( K-1, K-1 ) / T
            AKP1 = A( K, K ) / T
            AKKP1 = A( K, K-1 ) / T
            D = T*( AK*AKP1-ONE )
            A( K-1, K-1 ) = AKP1 / D
            A( K, K ) = AK / D
            A( K, K-1 ) = -AKKP1 / D
*
*           COMPUTE COLUMNS K-1 AND K OF THE INVERSE.
*
            IF( K.LT.N ) THEN
               CALL DCOPY( N-K, A( K+1, K ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1,
     $                     ZERO, A( K+1, K ), 1 )
               A( K, K ) = A( K, K ) - DDOT( N-K, WORK, 1, A( K+1, K ),
     $                     1 )
               A( K, K-1 ) = A( K, K-1 ) -
     $                       DDOT( N-K, A( K+1, K ), 1, A( K+1, K-1 ),
     $                       1 )
               CALL DCOPY( N-K, A( K+1, K-1 ), 1, WORK, 1 )
               CALL DSYMV( UPLO, N-K, -ONE, A( K+1, K+1 ), LDA, WORK, 1,
     $                     ZERO, A( K+1, K-1 ), 1 )
               A( K-1, K-1 ) = A( K-1, K-1 ) -
     $                         DDOT( N-K, WORK, 1, A( K+1, K-1 ), 1 )
            END IF
            KSTEP = 2
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           INTERCHANGE ROWS AND COLUMNS K AND KP IN THE TRAILING
*           SUBMATRIX A(K-1:N,K-1:N)
*
            IF( KP.LT.N )
     $         CALL DSWAP( N-KP, A( KP+1, K ), 1, A( KP+1, KP ), 1 )
            CALL DSWAP( KP-K-1, A( K+1, K ), 1, A( KP, K+1 ), LDA )
            TEMP = A( K, K )
            A( K, K ) = A( KP, KP )
            A( KP, KP ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = A( K, K-1 )
               A( K, K-1 ) = A( KP, K-1 )
               A( KP, K-1 ) = TEMP
            END IF
         END IF
*
         K = K - KSTEP
         GO TO 50
   60    CONTINUE
      END IF
*
      RETURN
*
*     END OF DSYTRI
*
      END

C    ------------------    BELOW IS DGESV  ------------------------
      SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK DRIVER ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGESV COMPUTES THE SOLUTION TO A REAL SYSTEM OF LINEAR EQUATIONS
*     A * X = B,
*  WHERE A IS AN N-BY-N MATRIX AND X AND B ARE N-BY-NRHS MATRICES.
*
*  THE LU DECOMPOSITION WITH PARTIAL PIVOTING AND ROW INTERCHANGES IS
*  USED TO FACTOR A AS
*     A = P * L * U,
*  WHERE P IS A PERMUTATION MATRIX, L IS UNIT LOWER TRIANGULAR, AND U IS
*  UPPER TRIANGULAR.  THE FACTORED FORM OF A IS THEN USED TO SOLVE THE
*  SYSTEM OF EQUATIONS A * X = B.
*
*  ARGUMENTS
*  =========
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF LINEAR EQUATIONS, I.E., THE ORDER OF THE
*          MATRIX A.  N >= 0.
*
*  NRHS    (INPUT) INTEGER
*          THE NUMBER OF RIGHT HAND SIDES, I.E., THE NUMBER OF COLUMNS
*          OF THE MATRIX B.  NRHS >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE N-BY-N COEFFICIENT MATRIX A.
*          ON EXIT, THE FACTORS L AND U FROM THE FACTORIZATION
*          A = P*L*U; THE UNIT DIAGONAL ELEMENTS OF L ARE NOT STORED.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,N).
*
*  IPIV    (OUTPUT) INTEGER ARRAY, DIMENSION (N)
*          THE PIVOT INDICES THAT DEFINE THE PERMUTATION MATRIX P;
*          ROW I OF THE MATRIX WAS INTERCHANGED WITH ROW IPIV(I).
*
*  B       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDB,NRHS)
*          ON ENTRY, THE N-BY-NRHS MATRIX OF RIGHT HAND SIDE MATRIX B.
*          ON EXIT, IF INFO = 0, THE N-BY-NRHS SOLUTION MATRIX X.
*
*  LDB     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY B.  LDB >= MAX(1,N).
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0:  IF INFO = I, U(I,I) IS EXACTLY ZERO.  THE FACTORIZATION
*                HAS BEEN COMPLETED, BUT THE FACTOR U IS EXACTLY
*                SINGULAR, SO THE SOLUTION COULD NOT BE COMPUTED.
*
*  =====================================================================
*
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DGETRF, DGETRS, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGESV ', -INFO )
         RETURN
      END IF
*
*     COMPUTE THE LU FACTORIZATION OF A.
*
      CALL DGETRF( N, N, A, LDA, IPIV, INFO )
      IF( INFO.EQ.0 ) THEN
*
*        SOLVE THE SYSTEM A*X = B, OVERWRITING B WITH X.
*
         CALL DGETRS( 'NO TRANSPOSE', N, NRHS, A, LDA, IPIV, B, LDB,
     $                INFO )
      END IF
      RETURN
*
*     END OF DGESV
*
      END
CUT HERE............
CAT > DGETRS.F <<'CUT HERE............'
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGETRS SOLVES A SYSTEM OF LINEAR EQUATIONS
*     A * X = B  OR  A' * X = B
*  WITH A GENERAL N-BY-N MATRIX A USING THE LU FACTORIZATION COMPUTED
*  BY DGETRF.
*
*  ARGUMENTS
*  =========
*
*  TRANS   (INPUT) CHARACTER*1
*          SPECIFIES THE FORM OF THE SYSTEM OF EQUATIONS:
*          = 'N':  A * X = B  (NO TRANSPOSE)
*          = 'T':  A'* X = B  (TRANSPOSE)
*          = 'C':  A'* X = B  (CONJUGATE TRANSPOSE = TRANSPOSE)
*
*  N       (INPUT) INTEGER
*          THE ORDER OF THE MATRIX A.  N >= 0.
*
*  NRHS    (INPUT) INTEGER
*          THE NUMBER OF RIGHT HAND SIDES, I.E., THE NUMBER OF COLUMNS
*          OF THE MATRIX B.  NRHS >= 0.
*
*  A       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          THE FACTORS L AND U FROM THE FACTORIZATION A = P*L*U
*          AS COMPUTED BY DGETRF.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,N).
*
*  IPIV    (INPUT) INTEGER ARRAY, DIMENSION (N)
*          THE PIVOT INDICES FROM DGETRF; FOR 1<=I<=N, ROW I OF THE
*          MATRIX WAS INTERCHANGED WITH ROW IPIV(I).
*
*  B       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDB,NRHS)
*          ON ENTRY, THE RIGHT HAND SIDE MATRIX B.
*          ON EXIT, THE SOLUTION MATRIX X.
*
*  LDB     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY B.  LDB >= MAX(1,N).
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      LOGICAL            NOTRAN
*     ..
*     .. EXTERNAL FUNCTIONS ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        SOLVE A * X = B.
*
*        APPLY ROW INTERCHANGES TO THE RIGHT HAND SIDES.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        SOLVE L*X = B, OVERWRITING B WITH X.
*
         CALL DTRSM( 'LEFT', 'LOWER', 'NO TRANSPOSE', 'UNIT', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        SOLVE U*X = B, OVERWRITING B WITH X.
*
         CALL DTRSM( 'LEFT', 'UPPER', 'NO TRANSPOSE', 'NON-UNIT', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        SOLVE A' * X = B.
*
*        SOLVE U'*X = B, OVERWRITING B WITH X.
*
         CALL DTRSM( 'LEFT', 'UPPER', 'TRANSPOSE', 'NON-UNIT', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        SOLVE L'*X = B, OVERWRITING B WITH X.
*
         CALL DTRSM( 'LEFT', 'LOWER', 'TRANSPOSE', 'UNIT', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        APPLY ROW INTERCHANGES TO THE SOLUTION VECTORS.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     END OF DGETRS
*
      END
CUT HERE............
CAT > DGETRF.F <<'CUT HERE............'
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     MARCH 31, 1993
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGETRF COMPUTES AN LU FACTORIZATION OF A GENERAL M-BY-N MATRIX A
*  USING PARTIAL PIVOTING WITH ROW INTERCHANGES.
*
*  THE FACTORIZATION HAS THE FORM
*     A = P * L * U
*  WHERE P IS A PERMUTATION MATRIX, L IS LOWER TRIANGULAR WITH UNIT
*  DIAGONAL ELEMENTS (LOWER TRAPEZOIDAL IF M > N), AND U IS UPPER
*  TRIANGULAR (UPPER TRAPEZOIDAL IF M < N).
*
*  THIS IS THE RIGHT-LOOKING LEVEL 3 BLAS VERSION OF THE ALGORITHM.
*
*  ARGUMENTS
*  =========
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX A.  M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE M-BY-N MATRIX TO BE FACTORED.
*          ON EXIT, THE FACTORS L AND U FROM THE FACTORIZATION
*          A = P*L*U; THE UNIT DIAGONAL ELEMENTS OF L ARE NOT STORED.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,M).
*
*  IPIV    (OUTPUT) INTEGER ARRAY, DIMENSION (MIN(M,N))
*          THE PIVOT INDICES; FOR 1 <= I <= MIN(M,N), ROW I OF THE
*          MATRIX WAS INTERCHANGED WITH ROW IPIV(I).
*
*  INFO    (OUTPUT) INTEGER
*          = 0:  SUCCESSFUL EXIT
*          < 0:  IF INFO = -I, THE I-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0:  IF INFO = I, U(I,I) IS EXACTLY ZERO. THE FACTORIZATION
*                HAS BEEN COMPLETED, BUT THE FACTOR U IS EXACTLY
*                SINGULAR, AND DIVISION BY ZERO WILL OCCUR IF IT IS USED
*                TO SOLVE A SYSTEM OF EQUATIONS.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. EXTERNAL FUNCTIONS ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     DETERMINE THE BLOCK SIZE FOR THIS ENVIRONMENT.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        USE UNBLOCKED CODE.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        USE BLOCKED CODE.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           FACTOR DIAGONAL AND SUBDIAGONAL BLOCKS AND TEST FOR EXACT
*           SINGULARITY.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           ADJUST INFO AND THE PIVOT INDICES.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           APPLY INTERCHANGES TO COLUMNS 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              APPLY INTERCHANGES TO COLUMNS J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              COMPUTE BLOCK ROW OF U.
*
               CALL DTRSM( 'LEFT', 'LOWER', 'NO TRANSPOSE', 'UNIT', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 UPDATE TRAILING SUBMATRIX.
*
                  CALL DGEMM( 'NO TRANSPOSE', 'NO TRANSPOSE', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     END OF DGETRF
*
      END
CUT HERE............
CAT > DLASWP.F <<'CUT HERE............'
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK AUXILIARY ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     OCTOBER 31, 1992
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DLASWP PERFORMS A SERIES OF ROW INTERCHANGES ON THE MATRIX A.
*  ONE ROW INTERCHANGE IS INITIATED FOR EACH OF ROWS K1 THROUGH K2 OF A.
*
*  ARGUMENTS
*  =========
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX A.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE MATRIX OF COLUMN DIMENSION N TO WHICH THE ROW
*          INTERCHANGES WILL BE APPLIED.
*          ON EXIT, THE PERMUTED MATRIX.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.
*
*  K1      (INPUT) INTEGER
*          THE FIRST ELEMENT OF IPIV FOR WHICH A ROW INTERCHANGE WILL
*          BE DONE.
*
*  K2      (INPUT) INTEGER
*          THE LAST ELEMENT OF IPIV FOR WHICH A ROW INTERCHANGE WILL
*          BE DONE.
*
*  IPIV    (INPUT) INTEGER ARRAY, DIMENSION (M*ABS(INCX))
*          THE VECTOR OF PIVOT INDICES.  ONLY THE ELEMENTS IN POSITIONS
*          K1 THROUGH K2 OF IPIV ARE ACCESSED.
*          IPIV(K) = L IMPLIES ROWS K AND L ARE TO BE INTERCHANGED.
*
*  INCX    (INPUT) INTEGER
*          THE INCREMENT BETWEEN SUCCESSIVE VALUES OF IPIV.  IF IPIV
*          IS NEGATIVE, THE PIVOTS ARE APPLIED IN REVERSE ORDER.
*
* =====================================================================
*
*     .. LOCAL SCALARS ..
      INTEGER            I, IP, IX
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DSWAP
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     INTERCHANGE ROW I WITH ROW IPIV(I) FOR EACH OF ROWS K1 THROUGH K2.
*
      IF( INCX.EQ.0 )
     $   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
*
      RETURN
*
*     END OF DLASWP
*
      END
CUT HERE............
CAT > DGETF2.F <<'CUT HERE............'
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK ROUTINE (VERSION 1.1) --
*     UNIV. OF TENNESSEE, UNIV. OF CALIFORNIA BERKELEY, NAG LTD.,
*     COURANT INSTITUTE, ARGONNE NATIONAL LAB, AND RICE UNIVERSITY
*     JUNE 30, 1992
*
*     .. SCALAR ARGUMENTS ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. ARRAY ARGUMENTS ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  PURPOSE
*  =======
*
*  DGETF2 COMPUTES AN LU FACTORIZATION OF A GENERAL M-BY-N MATRIX A
*  USING PARTIAL PIVOTING WITH ROW INTERCHANGES.
*
*  THE FACTORIZATION HAS THE FORM
*     A = P * L * U
*  WHERE P IS A PERMUTATION MATRIX, L IS LOWER TRIANGULAR WITH UNIT
*  DIAGONAL ELEMENTS (LOWER TRAPEZOIDAL IF M > N), AND U IS UPPER
*  TRIANGULAR (UPPER TRAPEZOIDAL IF M < N).
*
*  THIS IS THE RIGHT-LOOKING LEVEL 2 BLAS VERSION OF THE ALGORITHM.
*
*  ARGUMENTS
*  =========
*
*  M       (INPUT) INTEGER
*          THE NUMBER OF ROWS OF THE MATRIX A.  M >= 0.
*
*  N       (INPUT) INTEGER
*          THE NUMBER OF COLUMNS OF THE MATRIX A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          ON ENTRY, THE M BY N MATRIX TO BE FACTORED.
*          ON EXIT, THE FACTORS L AND U FROM THE FACTORIZATION
*          A = P*L*U; THE UNIT DIAGONAL ELEMENTS OF L ARE NOT STORED.
*
*  LDA     (INPUT) INTEGER
*          THE LEADING DIMENSION OF THE ARRAY A.  LDA >= MAX(1,M).
*
*  IPIV    (OUTPUT) INTEGER ARRAY, DIMENSION (MIN(M,N))
*          THE PIVOT INDICES; FOR 1 <= I <= MIN(M,N), ROW I OF THE
*          MATRIX WAS INTERCHANGED WITH ROW IPIV(I).
*
*  INFO    (OUTPUT) INTEGER
*          = 0: SUCCESSFUL EXIT
*          < 0: IF INFO = -K, THE K-TH ARGUMENT HAD AN ILLEGAL VALUE
*          > 0: IF INFO = K, U(K,K) IS EXACTLY ZERO. THE FACTORIZATION
*               HAS BEEN COMPLETED, BUT THE FACTOR U IS EXACTLY
*               SINGULAR, AND DIVISION BY ZERO WILL OCCUR IF IT IS USED
*               TO SOLVE A SYSTEM OF EQUATIONS.
*
*  =====================================================================
*
*     .. PARAMETERS ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. LOCAL SCALARS ..
      INTEGER            J, JP
*     ..
*     .. EXTERNAL FUNCTIONS ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. EXTERNAL SUBROUTINES ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. INTRINSIC FUNCTIONS ..
      INTRINSIC          MAX, MIN
*     ..
*     .. EXECUTABLE STATEMENTS ..
*
*     TEST THE INPUT PARAMETERS.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     QUICK RETURN IF POSSIBLE
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      DO 10 J = 1, MIN( M, N )
*
*        FIND PIVOT AND TEST FOR SINGULARITY.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           APPLY THE INTERCHANGE TO COLUMNS 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           COMPUTE ELEMENTS J+1:M OF J-TH COLUMN.
*
            IF( J.LT.M )
     $         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           UPDATE TRAILING SUBMATRIX.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     END OF DGETF2
*
      END
CUT HERE............
