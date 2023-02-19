      SUBROUTINE  EIGEN_VAL(A,B,C,J_inp,KA_inp,KC_inp,E_ch) !!! COMPUTES ASYMETRIC 
! This subroutine is updated by Bikramaditya Mandal
      IMPLICIT NONE!!! A,B,C ,J, Ka,Kc - INPUT, E - output
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
! input
      CHARACTER*1 jobz
      CHARACTER*1 range
      CHARACTER*1 uplo
      REAL*8 vl,vu,abstol
      INTEGER il,iu,M,ldz,lwork,info,lda 
      INTEGER i,k
c      OPEN(1234,file = 'energies.dat')
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
	
      CALL DSYEVX( JOBZ, RANGE, UPLO, N, hamil, LDA, VL, VU, IL, IU,
     $                   ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,
     $                   IFAIL, INFO )

      DO J1 = J,J
      DO K1 = -J1,J1!-(J1 - mod(J1,2)),J1 - mod(J1,2),2 
      CALL STATE(J1,K1,KA,KC,exst)	   
c      WRITE (*,'(i3,1x,i3,1x,i3,1x,F9.3)') J1,KA,KC, w(J1+K1+1)
      IF(KA_inp.eq.KA .and. KC_inp.eq.KC) E_ch = w(J1+K1+1)	  
      ENDDO
      ENDDO 
	   
      DEALLOCATE(w,z,work,iwork,ifail,hamil)
      ENDDO
      END

      SUBROUTINE hamilt_kyro(E,J1,J2,K1,K2,A,B,C)!!! THE MATRIX ELEMENT OF ASYMETRIC TOP HAMILTONIN IN THE BASIS OF WIGNER D-functions
! This subroutine is updated by Bikramaditya Mandal
      IMPLICIT NONE
      INTEGER J1,J2,K1,K2,Kav
      REAL*8 E,ebuff,jxy,jzz,jj,jzz2
      REAL*8 A,B,C
      REAL*8 PKKJ,PK,PPK,GK,sw
      Ebuff =0d0
      jzz = k1**2
      jzz2 = k2**2
	  jj = J1*(J1+1d0)
      IF(J1 .eq. J2) THEN
      IF(K1 .ne. K2) THEN
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
	  
      SUBROUTINE  EIGEN_VEC(A,B,C,J_inp,KA_inp,KC_inp,Z_OUT) !!! COMPUTE EIGENVALUES OF ASYMETRIC TOP
! This subroutine is updated by Bikramaditya Mandal
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
c      OPEN(1234,file = 'energies.dat')
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
	
      CALL DSYEVX( JOBZ, RANGE, UPLO, N, hamil, LDA, VL, VU, IL, IU,
     $                   ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,
     $                   IFAIL, INFO )

      DO J1 = J,J
      DO K1 = -J1,J1!-(J1 - mod(J1,2)),J1 - mod(J1,2),2 
      CALL STATE(J1,K1,KA,KC,exst)	   
c      WRITE (*,'(i3,1x,i3,1x,i3,1x,F9.3)') J1,KA,KC, w(J1+K1+1)
      IF(KA_inp.eq.KA .and. KC_inp.eq.KC) THEN 
      E_ch = w(J1+K1+1)
      DO K2=-J1,J1
      Z_OUT(K2+J1+1) = z(J1+K2+1,J1+K1+1)
      ENDDO
      ENDIF	  
      ENDDO
      ENDDO 
	   
      DEALLOCATE(w,z,work,iwork,ifail,hamil)
      ENDDO
      END 	  