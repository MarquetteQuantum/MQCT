      SUBROUTINE  EIGEN_VAL(A,B,C,J_inp,KA_inp,KC_inp,E_ch) !!! COMPUTES ASYMETRIC 
! This subroutine is updated by Bikramaditya Mandal
!----------------------------------------------------------------------------------------------------------------------------------------------------
! Comments by Carolin - 08/26/2024

! Purpose:
! This subroutine computes the eigenvalues (energies) and eigenvectors of an 
! asymmetric Hamiltonian matrix..
!
! INPUT PARAMETERS:
! A, B, C                - Rotational constants of the molecule
! J_inp, KA_inp, KC_inp  - Total angular momentum quantum number & Projection quantum numbers
!
! OUTPUT:
! E_ch            - Computed energy level (eigenvalue) for the given quantum numbers
!----------------------------------------------------------------------------------------------------------------------------------------------------

      IMPLICIT NONE!!! A,B,C ,J, Ka,Kc - INPUT, E - output
      LOGICAL exst	  
      REAL*8 A,B,C,E_ch	  
      INTEGER J1,K1,J2,K2,KA,KC,J_inp,KA_inp,KC_inp
      INTEGER Jmax,J
      INTEGER N, count1,count2,n_buff
      REAL*8 E, norm
      REAL*8,ALLOCATABLE :: hamil(:,:)  ! 2D array to store the Hamiltonian matrix
      REAL*8,ALLOCATABLE :: w(:)        ! array to store eigenvalues (energy)
      REAL*8,ALLOCATABLE :: z(:,:)      ! array to store eigenvectors
      REAL*8,ALLOCATABLE :: work(:)
      INTEGER, allocatable :: iwork(:)
      INTEGER, allocatable :: ifail(:)
      CHARACTER*1 jobz
      CHARACTER*1 range
      CHARACTER*1 uplo
      REAL*8 vl,vu,abstol
      INTEGER il,iu,M,ldz,lwork,info,lda 
      INTEGER i,k
c      OPEN(1234,file = 'energies.dat')

!      OPEN(324, FILE='hamiltonian_matrix.txt') ! Uncomment if one needs to print and see Hamiltonian matrix.																											 
      DO J=J_inp,J_inp
      N = (J+1)**2-J**2  ! matrix size calculation
      jobz = 'V'
      RANGE = 'A'
      UPLO = 'L'
      lda = N
      VL = 0d0
      VU = 0d0	  
	  il = 1
	  iu = 1
      abstol = 0d0
      M = N
      ALLOCATE(w(N))
      ldz = N
      ALLOCATE(z(N,N))
      lwork = 8*N
      ALLOCATE(work(lwork))
      ALLOCATE(iwork(5*N))
      ALLOCATE(ifail(N))
      count1 = 0
      count2 = 0
      ALLOCATE(hamil(N,N))
      hamil = 0d0
! This is a nested loop that builds the Hamiltonian matrix by calling hamilt_kyro subroutine for each matrix element.																													  
      DO J1 = J,J
      DO K1= -J1,J1
      count1 = count1 +1
      DO J2 = J,J
      DO K2= -J2,J2
      count2 = count2 +1
      CALL hamilt_kyro(E,J1,J2,K1,K2,A,B,C)
      hamil(count1,count2) = E
      ENDDO
      ENDDO 
      count2 = 0
      ENDDO
      ENDDO
	
! Print the Hamiltonian matrix
!	  WRITE(324,*) "Hamiltonian Matrix for J, A, B, C =", J, A, B, C
!	  DO i = 1, N
!		 WRITE(324, '(*(F12.6))') (hamil(i,k), k = 1, N)
!	  END DO
!	  WRITE(324,*) "" 
		  
! This calls the LAPACK subroutine DSYEVX to compute eigenvalues and eigenvectors of the Hamiltonian matrix.
      CALL DSYEVX( JOBZ, RANGE, UPLO, N, hamil, LDA, VL, VU, IL, IU,
     $                   ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,
     $                   IFAIL, INFO )
! This loop processes the result, calling subroutine STATE and setting up the E_ch output value based on the corresponding KA & KC
      DO J1 = J,J
      DO K1 = -J1,J1!-(J1 - mod(J1,2)),J1 - mod(J1,2),2 
      CALL STATE(J1,K1,KA,KC,exst)	   
!      WRITE (*,'(i3,1x,i3,1x,i3,1x,F9.3)') J1,KA,KC, w(J1+K1+1)
      IF(KA_inp.eq.KA .and. KC_inp.eq.KC) E_ch = w(J1+K1+1)	  
      ENDDO
      ENDDO 
	   
      DEALLOCATE(w,z,work,iwork,ifail,hamil)
      ENDDO
		close(324)			
      END

      SUBROUTINE hamilt_kyro(E,J1,J2,K1,K2,A,B,C)!!! THE MATRIX ELEMENT OF ASYMETRIC TOP HAMILTONIN IN THE BASIS OF WIGNER D-functions
! This subroutine is updated by Bikramaditya Mandal
!----------------------------------------------------------------------------------------------------------------------------------------------------
! Comments by Carolin - 08/26/2024

! Purpose: This subroutine computes elements of hamiltonian_matrix for the asymmetric top molecule when called repeatedly with different
! J and K values, it builds up the complete matrix needed for energy level calculations.
!----------------------------------------------------------------------------------------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER J1,J2,K1,K2,Kav
      REAL*8 E,ebuff,jxy,jzz,jj,jzz2
      REAL*8 A,B,C
      REAL*8 PKKJ,PK,PPK,GK,sw
      Ebuff =0d0
      jzz = k1**2
      jzz2 = k2**2
	  jj = J1*(J1+1d0)
      IF(J1 .eq. J2) THEN  ! checking if we're calculating a diagonal or off-diagonal element in J
      IF(K1 .ne. K2) THEN  ! checking if we're calculating an off-diagonal element in K.
      IF (K1 .eq. (K2 - 2) .or. K1 .eq. (K2+2)) THEN
      Kav = (K1 + K2)/2
      jxy = sqrt((J1+1d0)*J1 - Kav*K1)*
     & sqrt((J1+1d0)*J1 - Kav*K2)/2d0
      Ebuff = 0.5*(A-B)*jxy
      ELSE
      Ebuff = 0d0
      ENDIF
      ELSE
      Ebuff = 0.5*(A+B)*jj + 
     & (C-0.5d0*(A+B))*jzz	  
      ENDIF
      ENDIF
      E = Ebuff	  
      END

      SUBROUTINE STATE(J,DK,KA,KC,exst) !!! CREATE A STATE OF ASYMETRIC TOP
! This subroutine is updated by Bikramaditya Mandal
!----------------------------------------------------------------------------------------------------------------------------------------------------
! Comments by Carolin - 08/26/2024
! This subroutine creates a state of an asymmetric top molecule
!
! Input Parameters:
!   J  - Integer: Total angular momentum quantum number
!   DK - Integer: -J < DK < J
!   D_K - Integer: Difference between Ka and Kc
! Output Parameters:
!   KA   - Integer: Ka quantum number
!   KC   - Integer: Kc quantum number
!   exst - Logical: Indicates if a valid state was found (.TRUE.) or not (.FALSE.)
!------------------------------------------------------------------------------------------------
      IMPLICIT NONE!!! J, -J<DK<J -INPUT, Ka,Kc - output, exst - check, if succesfully finished
      LOGICAL exst
      INTEGER J,DK,KA,KC,i
      INTEGER K_A(2*J+1),K_C(2*J+1),D_K(2*J+1)
      K_A(1) = 0
      exst = .FALSE.	  
      DO i=2,J*2,2
      K_A(i) = i/2
      K_A(i+1) = i/2 	  
      ENDDO	
      DO i=1,2*J+1
      K_C(i) = K_A(2*J+2-i) 
      D_K(i) = K_A(i) - K_C(i)
      IF(DK.eq.D_K(i)) THEN
      KA = K_A(i)
      KC = K_C(i)
      exst = .TRUE.	  
      ENDIF	  
      ENDDO	  
      END SUBROUTINE STATE

      SUBROUTINE hpsort(n,ra,ja) !! Heapsort !! SEE NUMERICAL RECIPIES
! This subroutine is updated by Bikramaditya Mandal
!----------------------------------------------------------------------------------------------------------------------------------------------------
! Comments by Carolin - 08/26/2024
! This subroutine implements the Heapsort algorithm to sort an array ra of real numbers in ascending order.
! It also provides a corresponding integer array ja that keeps track of the orginal indices of the elements in ra.
! 
! INPUT PARAMETERS:
! n : The number of elements to sort
! ra: An array of real numbers to be sorted
! ja: An integer array for tracking original indices
!
! note: This in-place sorting means that the subroutine modifies the input arrays directly, rather than creating new arrays for the output. 
! This approach is memory-efficient as it doesn't require additional storage proportional to the input size.
!-----------------------------------------------------------------------------------------------------------------------------------------------------
      INTEGER n,ja(n)
      REAL*8 ra(n)
*       Sorts an array ra(1:n) into ascending numerical order using
*       the Heapsort algorithm.
      INTEGER i,ir,j,l,kja
      REAL*8 rra
      DO i=1,n
        ja(i)=i
      ENDDO
      IF ( n<2 ) RETURN
      l=n/2+1
      ir=n
   10 CONTINUE
        IF ( l>1 ) THEN
	     l=l-1
	     rra=ra(l)
	     kja=ja(l)
        ELSE
	     rra=ra(ir)
	     kja=ja(ir)
	     ra(ir)=ra(1)
	     ja(ir)=ja(1)
	     ir=ir-1
	  if ( ir==1 ) then
	      ra(1)=rra
	      ja(1)=kja
        RETURN
       ENDIF
        ENDIF
        i=l
       j=l+l
   20   IF ( j<=ir ) THEN
          IF ( j<ir ) THEN
      IF ( ra(j)<ra(j+1) ) j=j+1
      ENDIF
	  IF ( rra<ra(j) ) THEN
      ra(i)=ra(j)
      ja(i)=ja(j)
      i=j
      j=j+j
      ELSE
      j=ir+1
      ENDIF
      GOTO 20
      ENDIF	  
      ra(i)=rra
      ja(i)=kja
      GOTO 10
      END

      FUNCTION KRONEKER(i,j)  !!! KRONEKER
! This function is updated by Bikramaditya Mandal
      INTEGER KRONEKER,i,j
      KRONEKER = 0
      IF(i.eq.j) KRONEKER = 1	  
      END FUNCTION KRONEKER	
	  
      FUNCTION EXST_STATE(J,KA,KC) !!! CHECK IF STATE EXIST
! This function is updated by Bikramaditya Mandal
!----------------------------------------------------------------------------------------------------------------------------------------------------
! Comments by Carolin - 08/26/2024
!
! Purpose:
! This function checks for the existence of a specific quantum state based on J. 
! INPUT PARAMETERS:
! J, KA, KC 
! Returns:
! Logiacal: .TRUE. if the state exists, .FALSE. otherwise  
!----------------------------------------------------------------------------------------------------------------------------------------------------

      IMPLICIT NONE
      LOGICAL EXST_STATE
      INTEGER J,DK,KA,KC,i
      INTEGER K_A(2*J+1),K_C(2*J+1),D_K(2*J+1)
      K_A(1) = 0
      EXST_STATE = .FALSE.	  
      DO i=2,J*2,2
      K_A(i) = i/2
      K_A(i+1) = i/2 	  
      ENDDO	
      DO i=1,2*J+1
      K_C(i) = K_A(2*J+2-i) 
      D_K(i) = K_A(i) - K_C(i)
      IF(KA.eq.K_A(i) .and. KC.eq.K_C(i)) THEN
      EXST_STATE = .TRUE.
      RETURN	   
      ENDIF	  
      ENDDO	  
      END FUNCTION EXST_STATE
	  
      SUBROUTINE  EIGEN_VEC(A,B,C,J_inp,KA_inp,KC_inp,Z_OUT) !!! COMPUTE eigenvectors OF ASYMETRIC TOP
! This subroutine is updated by Bikramaditya Mandal
!----------------------------------------------------------------------------------------------------------------------------------------------------
! Comments by Carolin - 08/26/2024

! Purpose:

!   This subroutine computes the eigenvectors of an asymmetric top molecule's
!   Hamiltonian for a specific J, Ka, and Kc. It constructs the Hamiltonian matrix, solves for all eigenstates
!   using LAPACK, and extracts the eigenvector for the specified state.
!
! Input Parameters:
!   A, B, C   - Real*8  : Rotational constants of the molecule
!   J_inp     - Integer : Total angular momentum quantum number
!   KA_inp    - Integer : Ka quantum number of the desired state
!   KC_inp    - Integer : Kc quantum number of the desired state
!
! Output Parameters:
!   Z_OUT     - Real*8 Array(2*J_inp+1) : Eigenvector for the 
!                                         specified state
!----------------------------------------------------------------------------------------------------------------------------------------------------

      IMPLICIT NONE
      LOGICAL exst	  
      REAL*8 A,B,C,E_ch	  
      INTEGER J1,K1,J2,K2,KA,KC,J_inp,KA_inp,KC_inp
      INTEGER Jmax,J
      INTEGER N, count1,count2,n_buff
      REAL*8 E, norm
      REAL*8,ALLOCATABLE :: hamil(:,:)
      REAL*8,ALLOCATABLE :: w(:)
      REAL*8,ALLOCATABLE :: z(:,:)
      REAL*8,ALLOCATABLE :: work(:)
      INTEGER, allocatable :: iwork(:)
      INTEGER, allocatable :: ifail(:)
      REAL*8 Z_OUT(2*J_inp+1)	  
! input
      CHARACTER*1 jobz
      CHARACTER*1 range
      CHARACTER*1 uplo
      REAL*8 vl,vu,abstol
      INTEGER il,iu,M,ldz,lwork,info,lda 
      INTEGER i,k
		
!      OPEN(1234,file = 'energies.dat')
 !     OPEN(UNIT=125, FILE='eigenvectors.txt')
!      OPEN(UNIT=121, FILE='Hamiltonian-eigenvectors.txt')

      DO J=J_inp,J_inp
      N = (J+1)**2-J**2
      jobz = 'V'
      RANGE = 'A'
      UPLO = 'L'
      lda = N
      VL = 0d0
      VU = 0d0	  
	  il = 1
	  iu = 1
      abstol = 0d0
      M = N
      ALLOCATE(w(N))
      ldz = N
      ALLOCATE(z(N,N))
      lwork = 8*N
      ALLOCATE(work(lwork))
      ALLOCATE(iwork(5*N))
      ALLOCATE(ifail(N))
      count1 = 0
      count2 = 0
      ALLOCATE(hamil(N,N))
      hamil = 0d0
      DO J1 = J,J
      DO K1= -J1,J1
      count1 = count1 +1
      DO J2 = J,J  !The loop is iterating over all combinations of K1 and K2.
      DO K2= -J2,J2
      count2 = count2 +1
      CALL hamilt_kyro(E,J1,J2,K1,K2,A,B,C)  ! The hamil array is being built up gradually with each iteration. 
      hamil(count1,count2) = E               ! The hamilt_kyro is calculating the matrix elements one by one and they are being placed into the correct positions in the hamil array.
      ENDDO
      ENDDO 
      count2 = 0
      ENDDO
      ENDDO
	
      CALL DSYEVX( JOBZ, RANGE, UPLO, N, hamil, LDA, VL, VU, IL, IU,
     $                   ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,
     $                   IFAIL, INFO )

      DO J1 = J,J
      DO K1 = -J1,J1!-(J1 - mod(J1,2)),J1 - mod(J1,2),2 
      CALL STATE(J1,K1,KA,KC,exst)	   
!      WRITE (*,'(i3,1x,i3,1x,i3,1x,F9.3)') J1,KA,KC, w(J1+K1+1)
      IF(KA_inp.eq.KA .and. KC_inp.eq.KC) THEN 
      E_ch = w(J1+K1+1)
!		WRITE(120,*) J1, K1,K2, KA_inp, KC_inp, w(J1+K1+1)													 
      DO K2=-J1,J1
      Z_OUT(K2+J1+1) = z(J1+K2+1,J1+K1+1)  ! eigenvectors
!!! Write each component of the eigenvector
!		WRITE(121,*) J1, K1, K2, KA_inp, KC_inp, w(J1+K1+1), Z_OUT(K2+J1+1)
      ENDDO
!	   WRITE(121,*) ! Empty line for readability 											   
      ENDIF	  
      ENDDO
      ENDDO 
	   
      DEALLOCATE(w,z,work,iwork,ifail,hamil)
      ENDDO
      END 	  