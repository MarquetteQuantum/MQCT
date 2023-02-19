      SUBROUTINE ASYM_TOP_VECTORS !!! COMPUTES ASSYMETRIC TOP EIGNEVALUES
      USE VARIABLES
      USE EIGEN_VECTORS
      IMPLICIT NONE
      INTEGER j1_t,j2_t,ka1_t,ka2_t,kc1_t,kc2_t,i
      INTEGER j_t,ka_t,kc_t
      INTEGER par_1,par_2,dk_var,dk_chann
      REAL*8 dk_mult
      IF(eig_vec_def) RETURN
      IF(user_defined) THEN
      IF(coll_type.gt.0) THEN
      M_VECTORS	 = M_EIGEN_VECTORS_USER  !!! READING USER INPUT
      ELSE
      ALLOCATE(M1_VECTORS(number_of_channels,jmax_included*2+1))
      ALLOCATE(M2_VECTORS(number_of_channels,jmax_included*2+1))
      DO i=1,number_of_channels
      j1_t = j1_ch(i)
      j2_t = j2_ch(i)	
      M1_VECTORS(i,1:1+2*j1_t) = M_EIGEN_VECTORS_USER(i,1:1+2*j1_t)
      M2_VECTORS(i,1:1+2*j2_t) = M_EIGEN_VECTORS_USER(i,2+2*j1_t:2+2
     & *(j1_t+j2_t))	  !!!! CASE OF ASYM + ASYM TOP
      ENDDO	  
      ENDIF
	  
	  SELECT CASE(coll_type)
      CASE(4)
      A_I = A
      B_I = B
      C_I = C	  
      DO i=1,number_of_channels
      j_t = j_ch(i)
      ka_t = ka_ch(i)
      kc_t = kc_ch(i)
      ENDDO
      CASE(8)
      A_I = A
      B_I = B
      C_I = C	
      DO i=1,number_of_channels
      j_t = j1_ch(i)
      ka_t = ka1_ch(i)
      kc_t = kc1_ch(i)
      ENDDO	
      CASE(9)
      A_I = A1
      B_I = B1
      C_I = C1
      DO i=1,number_of_channels
      j_t = j1_ch(i)
      ka_t = ka1_ch(i)
      kc_t = kc1_ch(i)
      ENDDO
      CASE(0)
      IF(identical_particles_defined) then
      ALLOCATE(parity_inversion(number_of_channels))
      ALLOCATE(p1p2_bk(2,number_of_channels))
	  end if
      A1_I = A1
      B1_I = B1
      C1_I = C1
      A2_I = A2
      B2_I = B2
      C2_I = C2
      DO i=1,number_of_channels
      j1_t = j1_ch(i)
      j2_t = j2_ch(i)
      ka1_t = ka1_ch(i)
      ka2_t = ka2_ch(i)
      kc1_t = kc1_ch(i)
      kc2_t = kc2_ch(i)
      IF(identical_particles_defined) THEN         !!! CASE OF IDENTICAL PARTICLES : PARITY INVERSION
      parity_inversion(i) = 1	  
      dk_mult = 0d0	  
      DO dk_chann=1,j1_t+1
      dk_var = j1_t*2+2 - dk_chann

      dk_mult = M1_VECTORS(i,dk_chann)*M1_VECTORS(i,dk_var) !!!! THE SIGN
!      IF(i.eq.5) PRINT*,i,	dk_mult   
      IF(dk_mult.eq.0) THEN
      IF(ABS(M1_VECTORS(i,dk_chann))+
     & ABS(M1_VECTORS(i,dk_var)).ne.0d0) STOP           !!! CHEKING THE SYMMMETRY
     & "STOP: ASYMETRIC TOP DIAG WRONG 1"
      CYCLE	  
      ELSE
      IF(dk_mult.gt.0d0) par_1=1      !!! IF BJ_-k = BJ_+k
      IF(dk_mult.lt.0d0) par_1=-1       !!! IF BJ_k - = -BJ_+k 
      EXIT	  
      ENDIF 	  
      ENDDO
!	  write(*,'(a, 1x, 3(i0), 1x, i0)')'1st',j1_t,ka1_t,kc1_t,par_1
	  
      dk_mult = 0d0	  
      DO dk_chann=1,j2_t+1
      dk_var = j2_t*2+2 - dk_chann
      dk_mult = M2_VECTORS(i,dk_chann)*M2_VECTORS(i,dk_var)      !!! DEFINING PARITY FOR THE SECOND MOLECULE
      IF(dk_mult.eq.0) THEN
      IF(ABS(M2_VECTORS(i,dk_chann))+
     & ABS(M2_VECTORS(i,dk_var)).ne.0d0) STOP
     & "STOP: ASYMETRIC TOP DIAG WRONG 2"	  
      CYCLE	  
      ELSE
      IF(dk_mult.gt.0d0) par_2=1
      IF(dk_mult.lt.0d0) par_2=-1
      EXIT	  
      ENDIF 	  	  
      ENDDO	
!	  write(*,'(a, 1x, 3(i0), 1x, i0)')'2nd',j2_t,ka2_t,kc2_t,par_2
      parity_inversion(i) = par_1*par_2
	  p1p2_bk(1,i) = par_1
	  p1p2_bk(2,i) = par_2
!      PRINT*,"parity_iversion",i,parity_inversion(i) !!!! TOTAL PARITY INVERSION	  
      ENDIF
	  
      ENDDO	
      END SELECT
	  
      eig_vec_def = .TRUE.	  
	  
	  
      else	  
      SELECT CASE(coll_type)
      CASE(4)
      A_I = A
      B_I = B
      C_I = C	  
      ALLOCATE(M_VECTORS(number_of_channels,jmax_included*2+1))
      M_VECTORS = 0d0
      DO i=1,number_of_channels                   !!! LOOP OVER CHANNELS
      j_t = j_ch(i)

      ka_t = ka_ch(i)

      kc_t = kc_ch(i)

      CALL EIGEN_VEC
     & (A,B,C,j_t,ka_t,kc_t,M_VECTORS(i,1:j_t*2+1))	  !!!! COMPUTING EIGENCVECTORS
	  
      ENDDO	
      eig_vec_def = .TRUE.
      CASE(8)
      A_I = A
      B_I = B
      C_I = C	  
      ALLOCATE(M_VECTORS(number_of_channels,jmax_included*2+1))
      M_VECTORS = 0d0
      DO i=1,number_of_channels
      j_t = j1_ch(i)

      ka_t = ka1_ch(i)

      kc_t = kc1_ch(i)

      CALL EIGEN_VEC
     & (A,B,C,j_t,ka_t,kc_t,M_VECTORS(i,1:j_t*2+1))	  
	  
      ENDDO	
      eig_vec_def = .TRUE.	  
      CASE(9)
      A_I = A1
      B_I = B1
      C_I = C1	  
      ALLOCATE(M_VECTORS(number_of_channels,jmax_included*2+1))
      M_VECTORS = 0d0
      DO i=1,number_of_channels
      j_t = j1_ch(i)

      ka_t = ka1_ch(i)

      kc_t = kc1_ch(i)

      CALL EIGEN_VEC
     & (A,B,C,j_t,ka_t,kc_t,M_VECTORS(i,1:j_t*2+1))	  
	  
      ENDDO	
      eig_vec_def = .TRUE.	  
      CASE(0)
      ALLOCATE(M1_VECTORS(number_of_channels,jmax_included*2+1))
      ALLOCATE(M2_VECTORS(number_of_channels,jmax_included*2+1))
      IF(identical_particles_defined) then
      ALLOCATE(parity_inversion(number_of_channels))
      ALLOCATE(p1p2_bk(2,number_of_channels))
	  end if
      M1_VECTORS = 0d0
      M1_VECTORS = 0d0
      A1_I = A1
      B1_I = B1
      C1_I = C1
      A2_I = A2
      B2_I = B2
      C2_I = C2
      DO i=1,number_of_channels
      j1_t = j1_ch(i)
      j2_t = j2_ch(i)
      ka1_t = ka1_ch(i)
      ka2_t = ka2_ch(i)
      kc1_t = kc1_ch(i)
      kc2_t = kc2_ch(i)
      CALL EIGEN_VEC
     & (A1,B1,C1,j1_t,ka1_t,kc1_t,M1_VECTORS(i,1:j1_t*2+1))	  
      CALL EIGEN_VEC
     & (A2,B2,C2,j2_t,ka2_t,kc2_t,M2_VECTORS(i,1:j2_t*2+1))
      IF(identical_particles_defined) THEN         !!! CASE OF IDENTICAL PARTICLES : PARITY INVERSION
      parity_inversion(i) = 1	  
      dk_mult = 0d0	  
      DO dk_chann=1,j1_t+1
      dk_var = j1_t*2+2 - dk_chann

      dk_mult = M1_VECTORS(i,dk_chann)*M1_VECTORS(i,dk_var) !!!! THE SIGN
!      IF(i.eq.5) PRINT*,i,	dk_mult   
      IF(dk_mult.eq.0) THEN
      IF(ABS(M1_VECTORS(i,dk_chann))+
     & ABS(M1_VECTORS(i,dk_var)).ne.0d0) STOP           !!! CHEKING THE SYMMMETRY
     & "STOP: ASYMETRIC TOP DIAG WRONG 1"
      CYCLE	  
      ELSE
      IF(dk_mult.gt.0d0) par_1=1      !!! IF BJ_-k = BJ_+k
      IF(dk_mult.lt.0d0) par_1=-1       !!! IF BJ_k - = -BJ_+k 
      EXIT	  
      ENDIF 	  
      ENDDO
	  
      dk_mult = 0d0	  
      DO dk_chann=1,j2_t+1
      dk_var = j2_t*2+2 - dk_chann
      dk_mult = M2_VECTORS(i,dk_chann)*M2_VECTORS(i,dk_var)      !!! DEFINING PARITY FOR THE SECOND MOLECULE
      IF(dk_mult.eq.0) THEN
      IF(ABS(M2_VECTORS(i,dk_chann))+
     & ABS(M2_VECTORS(i,dk_var)).ne.0d0) STOP
     & "STOP: ASYMETRIC TOP DIAG WRONG 2"	  
      CYCLE	  
      ELSE
      IF(dk_mult.gt.0d0) par_2=1
      IF(dk_mult.lt.0d0) par_2=-1
      EXIT	  
      ENDIF 	  	  
      ENDDO	

      parity_inversion(i) = par_1*par_2
	  p1p2_bk(1,i) = par_1
	  p1p2_bk(2,i) = par_2
!	  write(*,'(a,1x,3(i0),2(1x,i0))')
!     & '1st',j1_t,ka1_t,kc1_t,par_1,p1p2_bk(1)
!	  write(*,'(a,1x,3(i0),2(1x,i0))')
!     & '2nd',j2_t,ka2_t,kc2_t,par_2,p1p2_bk(2)
!      PRINT*,"parity_iversion",i,parity_inversion(i) !!!! TOTAL PARITY INVERSION	  
      ENDIF
	  
      ENDDO	
	  
      eig_vec_def = .TRUE.	  
      END SELECT
c      PRINT*, "VECTORS_COMPUTED", M2_VECTORS(:,:)     	  
	  end if
      END SUBROUTINE ASYM_TOP_VECTORS	  
	  
      SUBROUTINE Djkm(real_part,imag_part,j,k,m,alpha,beta,gamma)!!! EIGENFUNCTION OF ASYMETRIC TOP/WIGNER D_FUNCTION
      USE FACTORIAL	  
      IMPLICIT NONE
      INTEGER j,k,m,  a,b,n, lambda
      REAL*8 cx(0:2*j),C_N_K
      EXTERNAL C_N_K	  
      REAL*8 real_part,imag_part,alpha,beta,gamma,dkm
      IF(j.lt.max(abs(k),abs(m))) THEN
	  PRINT*,"ERROR:D_WIGNER NOT DEFINED"
      PRINT*,"j=",j
      PRINT*,"k=",k
      STOP	  
      RETURN	  
      ENDIF	  
      cx = 0d0	  
      n = min(j+k,j-k,j+m,j-m)
c      PRINT*," n",n
      IF(n.eq.j+m) THEN
      a = k-m
!      lambda = k - m
      lambda = 0
      ENDIF  
      IF(n.eq.j-m) THEN
      a = m-k
!      lambda = 0
      lambda = k - m
      ENDIF
      IF(n.eq.j+k) THEN
      a = m-k
!      lambda = 0
      lambda = k - m
      ENDIF
      IF(n.eq.j-k) THEN
      a = k-m
!      lambda = k - m
      lambda = 0
      ENDIF	
      b= 2*j -2*n - a
c      PRINT*," a",a	  
      IF(a.lt.0) STOP "ERORR: WIGNER a<0"
      IF(n.lt.0) STOP "ERORR: WIGNER n<0"
      IF(b.lt.0) STOP "ERORR: WIGNER b<0"	  
c      ALLOCATE(cx(0:n))	  
      CALL jacobi_poly ( n, dble(a), dble(b), cos(beta), cx )
      dkm  = cx(n)*(-1)**lambda*dsqrt(dble(C_N_K(2*j-n,n+a)))/
     & dsqrt(dble(C_N_K(n+b,b)))*dsin(beta/2d0)**a*dcos(beta/2d0)**b
c      PRINT*,"dkm", dkm
      real_part = dkm*dcos(k*gamma+m*alpha)
      imag_part = dkm*dsin(k*gamma+m*alpha)	
!	  write(*,'(3(i4,2x),8(f12.5,2x))')J, K, M, alpha, gamma, beta, 
!     & dkm, dcos(k*gamma+m*alpha), dsin(k*gamma+m*alpha), 
!     & real_part, imag_part
c      DEALLOCATE(cx) 	  
      END SUBROUTINE Djkm

      SUBROUTINE INI_FACT !!!! INTILIZING FACTORIAL
      USE FACTORIAL
      USE VARIABLES	  
      IMPLICIT NONE
      INTEGER i, fact_max
      fact_max=max(99,4*jmax_included)	  
      ALLOCATE(LOGFACT(0:fact_max))	  
      LOGFACT(0) = 0d0      
      DO i=1,fact_max
      LOGFACT(i) = LOGFACT(i-1)+dlog(dble(i))  
      ENDDO	  
      END SUBROUTINE INI_FACT

      REAL*8 FUNCTION C_N_K(n,k) !!!!! BINOMINAL COEFFICENT
      USE FACTORIAL
      IMPLICIT NONE
      INTEGER n,k
      IF(n.lt.k) STOP" ERROR: N>K FOR C_N_K"
      IF(k.lt.0) STOP "ERORRL K<0 FOR C_N_K"	  
      C_N_K = exp(LOGFACT(n)-LOGFACT(n-k)-LOGFACT(k))	  
      END FUNCTION C_N_K

      SUBROUTINE INI_POT_STORAGE  !!! INTIALIZING THE GRID ARRAYS FOR NUMERICAL INTEGRATION OF POTENTIAL OVER WAVEFUNCTIONS
      USE POT_STORE
      USE GRID_INTEGR	  
      USE VARIABLES
      USE CONSTANTS
      USE FACTORIAL	
      USE MPI_DATA
      USE MPI	  
      IMPLICIT NONE		  
      REAL*8, ALLOCATABLE :: buffer(:,:,:)
      REAL*8, ALLOCATABLE :: buffer_3_3(:,:,:,:,:)	  
      REAL*8 x1,x2,x,R,V,r_vib1,r_vib2
      REAL*8, ALLOCATABLE :: wpsir(:),wpsir_storage(:,:),
     & wpsir_storage_2(:,:,:)
      REAL*8 effect_length
      LOGICAL, ALLOCATABLE :: v_all_unique(:),v_all_unique1(:),
     & v_all_unique2(:)
      INTEGER i1,i2,i3,i4,i5,i6,i,k,vw,i_vib1,i_vib2
	  REAL*8 morse_r_db, r_vib_DB													!Dulat
      INTEGER task_size,tag2,dest,source	  
      INTEGER coll_type_checker, 
     & nb1_checker,na1_checker,ng1_checker,
     & nb2_checker,na2_checker,ng2_checker,nr_checker
      INTEGER status(MPI_STATUS_SIZE)
      INTEGER CASE_LOOP
      LOGICAL file_exst	  
      conv_unit_r = 1d0
      IF(make_grid_file.and. nproc.lt.n_r_coll .and. coll_type.gt.4)
     & THEN
      IF(myid.eq.0 .and. .not.grid_file_found)
     & PRINT*,"WARNING:NUMBER OF R_G_POINTS>=N_PROCESSORS"	 
      ENDIF	 
      IF(angs_unit) conv_unit_r = a_bohr
      IF(ALLOCATED(xb)) THEN
      IF(MYID.eq.0) PRINT*,
     & "GAUSS_LEGENDER GRID AND/OR VGRID ALREADY ALLOCATED"
      RETURN	 
      ENDIF	  
      IF(MYID.EQ.0) PRINT*, "THE GRID INI"	 	  
      CASE_FORM: DO CASE_LOOP=1,1	  
      SELECT CASE(coll_type)
      CASE(1)
      ALLOCATE(V_2_1(n_beta,n_r_coll))!!! ANGLE BETA GRID
      x2 = pi
      x1 = 0
      xbm = 0.5 * (x2 + x1)
   	  xbl = 0.5 * (x2 - x1)
      ALLOCATE(xb(n_beta),wb(n_beta))	  
      CALL gauleg(n_beta, xb, wb)
      IF(grid_file_found) EXIT CASE_FORM	  
      DO i=1,n_r_coll	  
      DO i2=1,n_beta
      R =  R_COM(i)*conv_unit_r
      beta = xbm + xbl*xb(i2)	  
      CALL V_POT_DIATOM_ATOM(V,R,beta) !!! SAVING POTENTIAL
      V_2_1(i2,i) = V
      IF(V.ne.V) THEN
      WRITE(*,*) i2,i,V
      STOP "ERROR IN POTENTIAL ROUTINE"	  
      ENDIF	  
      ENDDO
      ENDDO
      CASE(2)
      ALLOCATE(V_VIB_2(n_beta,n_r_vib,n_r_coll)) !!! potential storage
      x2 = pi
      x1 = 0
      xbm = 0.5 * (x2 + x1)
   	  xbl = 0.5 * (x2 - x1)
      ALLOCATE(xb(n_beta),wb(n_beta))	  !!!! ANGULAR GRID FOR BETA
      CALL gauleg(n_beta, xb, wb)

c      PRINT*,x2,x1
      IF(.not.user_defined) THEN
      x2 = r_vib_diatom_max
      x1 = r_vib_diatom_min	  
      xgm = 0.5 * (x2 + x1)
   	  xgl = 0.5 * (x2 - x1)
      ELSE
      xgl = 1d0
      xgm = 0d0	  
      ENDIF	  
      ALLOCATE(wpsir(n_r_vib))!!!! VIB WAVEFUNCTION ARRAY - BUFFER
      ALLOCATE(wf_vib_part_diatom(n_r_vib,number_of_channels),
     & r_vb_dt_inegrator(n_r_vib)) 	  !!! ARRAYS CONTAIN GRID AND WAVEFUNCTIONS
      ALLOCATE(xg(n_r_vib),wg(n_r_vib))
      IF(.not.user_defined) THEN  !!!! MORSE DEFINED
	  
! Generation of Gauss-Legendre grid and weights 
! original by Semenov, commented out by Dulat
!      CALL gauleg(n_r_vib, xg, wg)
!      DO i=1,n_r_vib	  
!      r_vb_dt_inegrator(i) = xgm + xgl*xg(i)
!      ENDDO

! Generation of equidistant grid and weight
! Dulat Bostan 4/23/2021 step size
	   DO i=1,n_r_vib	  
        r_vb_dt_inegrator(i)=x1+(DBLE(i)-0.5d0)*2*xgl/n_r_vib
        wg(i)= 2.d0/n_r_vib
       ENDDO
!Dulat Bostan: end

      ELSE
      r_vb_dt_inegrator = r_grid_vib  !!!! WHAT USER  GAVE
      wg = jacobian_vib
      wf_vib_part_diatom = vibrational_wavefunctions
      ENDIF  
!      STOP 
      IF(MYID.EQ.0) PRINT*,"ARRAYS ALLOCATED" 
      DO i=1,n_r_coll	  
      DO i3=1,n_r_vib
      DO i2=1,n_beta
      R =  R_COM(i)*conv_unit_r
      beta = xbm + xbl*xb(i2)	  													!CALLING PES SUBROUTINE
	  r_vib_DB = r_vb_dt_inegrator(i3)*conv_unit_r									!Dulat: converting vib-distance depending on the input
      CALL V_POT_VIB_DIATOM_ATOM(V,R,r_vb_dt_inegrator(i3),beta)
      V_VIB_2(i2,i3,i) = V
      IF(V.eq.0d0) THEN
      WRITE(*,*) i2,i3,i,V
      STOP "ERROR IN POTENTIAL ROUTINE"	  
      ENDIF	  
      ENDDO
      ENDDO
      ENDDO	  

     !!!! TESTING	  
      IF(.not.user_defined) THEN													!!!! FOR MORSE  
!      morse_de = WE**2/4d0/XE/eVtown/autoeV										!Dulat: moved down inside if condition because it now takes care of input argument
!      morse_a = (WE/eVtown/autoeV)*dsqrt(atomic_red_mass/2d0/morse_de
!     & * amutoau)

! Dulat Bostan - Morse parameters check 7/15/2022
	  if (morse_de.le.0) then 														! If Morse depth is not specified in the INPUT file
	  if(myid.eq.0) print*, 
     & "Morse Depth is not specified in the input"
	  morse_de = WE**2/4d0/XE/eVtown/autoeV
	  if(myid.eq.0) print*, "Morse Depth: morse_de =", 
     & morse_de, "a.u."
	  else
	  morse_de = morse_de/eVtown/autoeV 											! Converting Morse depth from cm-1 to a.u.
	  if(myid.eq.0) print*, "Morse Depth: morse_de =", 
     & morse_de, "a.u."
	  endif
	  
	  if (morse_a.le.0) then  														! If Morse width is not specified in the INPUT file
	  if(myid.eq.0) print*, 
     & "Morse Width is not specified in the input"
	  morse_a = (WE/eVtown/autoeV)*dsqrt(atomic_red_mass/2d0/morse_de
     & * amutoau)
	  if(myid.eq.0) print*, "Morse Width: morse_a =", 
     & morse_a, "a.u."
	  else
	  if(myid.eq.0) print*, "Morse Width: morse_a =", 
     & morse_a, "a.u."
	  endif

      IF(MYID.eq.0) PRINT*,"COMPUTING VIBRATIONAL WAVE_FUNCTIONS"
!      PRINT*,"morse_de",morse_de
!      PRINT*,"morse_a",morse_a
!      STOP
      ALLOCATE(v_all_unique(vmax_included+1))
      ALLOCATE(wpsir_storage(n_r_vib,vmax_included+1))	  
      DO k=1,vmax_included+1
      v_all_unique(k) = .FALSE.	  
      ENDDO	  
      DO k=1,number_of_channels
      vw = v_ch(k)
      IF(v_all_unique(vw+1)) THEN
      wf_vib_part_diatom(1:n_r_vib,k)=wpsir_storage(1:n_r_vib,vw+1)
      CYCLE	  
      ENDIF	  
!      PRINT*,"OLD",r_vb_dt_inegrator(1)	  
      CALL VIBRATIONAL_WAVEFUNCTION !!! COMPUTING WAVEFUNCTION
     & (wpsir,r_vb_dt_inegrator,n_r_vib,vw,atomic_red_mass*amutoau,
     &  morse_de,morse_a,morse_re)
      wf_vib_part_diatom(1:n_r_vib,k) = wpsir(1:n_r_vib)
      wpsir_storage(1:n_r_vib,vw+1) = wpsir(1:n_r_vib)
      v_all_unique(vw+1) = .TRUE.	  
!      PRINT*,"NEW",r_vb_dt_inegrator(1)	  
      ENDDO
!      PRINT*,"NEW",r_vb_dt_inegrator(1)
!      STOP	  
!      DEALLOCATE(wpsir)
      IF(MYID.EQ.0) THEN
      OPEN(456,FILE="VIBRATIONAL_WAVE_FUNCTIONS.out",
     & ACTION="WRITE")
!      ! WRITTING A GRID	 
      WRITE(456,*) "R(i),J(i)" ! A HEADER
      DO i=1,n_r_vib! LOOP OVER VIBRATIONAL GRID	  
      WRITE(456,*)r_vb_dt_inegrator(i),wg(i)*xgl ! r_vb_dt_inegrator(i)-					!Dulat: multiplied by xgl because it's already part of the weight
! r-value  for i-point of the grid, wg(i) - weight(jacobian)
      ENDDO
      WRITE(456,*) "WAVEFUNCTIONS ARE LISTED BELOW"	  
      DO k=1,number_of_channels !! k - channel number
c      !! WRITING THE WAVEFUNCTIONS	  
      WRITE(456,*) "CHANNEL=",k !!! A HEADER
      DO i=1,n_r_vib	  
      WRITE(456,*)wf_vib_part_diatom(i,k) ! VALUE OF WAVEFUNCTION FOR k-channel at i-point of the grid
      ENDDO
      WRITE(456,*) !! LEAVE IT BLANK	  
      ENDDO	
      CLOSE(456)
      ENDIF
      ENDIF
!      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      CASE(3)
      ALLOCATE(V_3_1(n_gamma,n_beta,n_r_coll))	 !!! POTENTIAL STORAGE 
      x2 = pi
      x1 = 0
      xbm = 0.5 * (x2 + x1)
   	  xbl = 0.5 * (x2 - x1)
      ALLOCATE(xb(n_beta),wb(n_beta))	  !!! GRID FOR BETA
      CALL gauleg(n_beta, xb, wb)
      x2 = pi*2
      x1 = 0
      xgm = 0.5 * (x2 + x1)
   	  xgl = 0.5 * (x2 - x1)
      ALLOCATE(xg(n_gamma),wg(n_gamma))	  !!!! GRID FOR GAMMA

! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      CALL gauleg(n_gamma, xg, wg)
	  else
	  if(myid.eq.0) 
     & print*, "Using Equidistant Grid for Angle Gamma"
	  do i = 1, n_gamma
	  xg(i) = x1 + (dble(i) - 0.50d0)*2d0*xgl/n_gamma
      wg(i)= 2.0d0/n_gamma
	  end do
	  end if
! Bikram End.
	  
      DO i=1,n_r_coll	  
      DO i2=1,n_beta
      DO i3=1,n_gamma
      beta = xbm + xbl*xb(i2)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i3)
	  else
	  gamma = xg(i3)
	  end if
! Bikram End.
      R =  R_COM(i)*conv_unit_r	  

      CALL V_POT_TOP_ATOM(V,R,beta,gamma)	  !!! CALLING PES SUBROUTINE
c      integrant	=  (- dr1*di2+di1*dr2)*wb(i2)*wg(i3)*dsin(beta)
c      intgrli  = intgrli+integrant*xbl*xgl*2d0*pi
      V_3_1(i3,i2,i) = V	  
      ENDDO
      ENDDO
      ENDDO
      CASE(4)
      ALLOCATE(V_3_1(n_gamma,n_beta,n_r_coll))	  !! THE SAME AS CASE(3)
      x2 = pi
      x1 = 0
      xbm = 0.5 * (x2 + x1)
   	  xbl = 0.5 * (x2 - x1)
      ALLOCATE(xb(n_beta),wb(n_beta))	  
      CALL gauleg(n_beta, xb, wb)
      x2 = pi*2
      x1 = 0
      xgm = 0.5 * (x2 + x1)
   	  xgl = 0.5 * (x2 - x1)
      ALLOCATE(xg(n_gamma),wg(n_gamma))	  
!      CALL gauleg(n_gamma, xg, wg)

! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      CALL gauleg(n_gamma, xg, wg)	
	  else
	  if(myid.eq.0) 
     & print*, "Using Equidistant Grid for Angle Gamma"
	  do i = 1, n_gamma
	  xg(i) = x1 + (dble(i) - 0.50d0)*2d0*xgl/n_gamma
      wg(i)= 2.0d0/n_gamma
	  end do
	  end if
! Bikram End.

      DO i=1,n_r_coll	  
      DO i2=1,n_beta
      DO i3=1,n_gamma
      beta = xbm + xbl*xb(i2)
!      gamma = xgm + xgl*xg(i3)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i3)
	  else
	  gamma = xg(i3)
	  end if
! Bikram End.
      R =  R_COM(i)*conv_unit_r	  

      CALL V_POT_TOP_ATOM(V,R,beta,gamma)	  
c      integrant	=  (- dr1*di2+di1*dr2)*wb(i2)*wg(i3)*dsin(beta)
c      intgrli  = intgrli+integrant*xbl*xgl*2d0*pi
      V_3_1(i3,i2,i) = V	  
      ENDDO
      ENDDO
      ENDDO 
      CASE(5)
      ALLOCATE(V_2_2(n_alpha2,n_beta2,n_beta1,n_r_coll)) !!! PES STORAGE
      x2 = pi
      x1 = 0
      xbm = 0.5 * (x2 + x1)
   	  xbl = 0.5 * (x2 - x1)
      ALLOCATE(xb(n_beta1),wb(n_beta1))	  !!! ALLOCATION OF ARRAYS: FOR BETA1
      CALL gauleg(n_beta1, xb, wb)
      x2 = pi
      x1 = 0
      xbbm = 0.5 * (x2 + x1)
   	  xbbl = 0.5 * (x2 - x1)
      ALLOCATE(xbb(n_beta2),wbb(n_beta2))	  !!! BETA 2
      CALL gauleg(n_beta2, xbb, wbb)
      x2 = pi*2
      x1 = 0
      xaam = 0.5 * (x2 + x1)
   	  xaal = 0.5 * (x2 - x1)
      ALLOCATE(xaa(n_alpha2),waa(n_alpha2))
!      x2 = pi*2
!      x1 = 0
!      xgm = 0.5 * (x2 + x1)
!   	  xgl = 0.5 * (x2 - x1)
!      ALLOCATE(xg(n_gamma1),wg(n_gamma1)) !!! GAMMA1
!      CALL gauleg(n_gamma1, xg, wg)

! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
!      CALL gauleg(n_gamma1, xg, wg)	
      CALL gauleg(n_alpha2, xaa, waa)	
	  else
	  if(myid.eq.0) 
     & print*, "Using Equidistant Grid for Angle Alpha"
!	  do i = 1, n_gamma1
!	  xg(i) = x1 + (dble(i) - 0.50d0)*2d0*xgl/n_gamma1
!      wg(i)= 2.0d0/n_gamma1
!	  end do
	  do i = 1, n_alpha2
	  xaa(i) = x1 + (dble(i) - 0.50d0)*2d0*xaal/n_alpha2
      waa(i)= 2.0d0/n_alpha2
	  end do
	  end if
! Bikram End.

      IF(make_grid_file) THEN	  
      IF(grid_file_found) THEN	  !!!! GRID FILE READING
      OPEN(1,FILE=potential_file_name,ACTION="READ",STATUS="OLD")
      READ(1,*)	coll_type_checker
      READ(1,*) nb1_checker,nb2_checker,na2_checker,nr_checker
      IF(coll_type_checker.ne.coll_type .or. nb1_checker.ne.n_beta1 .or.
     & nb2_checker.ne.n_beta2 .or. na2_checker.ne.n_alpha2 .or.
     & nr_checker.ne.n_r_coll  )
     & STOP "WRONG GRID SIZE(S) OR COLLISION TYPE"	 
      DO i=1,n_r_coll
      DO i1=1,n_beta1
      DO i2=1,n_beta2
      DO i3=1,n_alpha2
      R =  R_COM(i)*conv_unit_r	  
      beta =  xbm + xbl*xb(i1)
      bbeta = xbbm + xbbl*xbb(i2)
!      gamma = xgm + xgl*xg(i3)	  
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
!      gamma = xgm + xgl*xg(i3)
      aalpha = xaam + xaal*xaa(i3)
	  else
!	  gamma = xg(i3)
	  aalpha = xaa(i3)
	  end if
! Bikram End.

!      CALL V_POT_DIAT_DIAT(V,R,alpha,beta,gamma)
!      V_2_2(i3,i2,i1,i) = V 	  
      READ(1,*) V_2_2(i3,i2,i1,i)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      CLOSE(1)
	  
      ELSE
      ALLOCATE(buffer(n_alpha2,n_beta2,n_beta1))
      tag2= 2
      buffer = 0d0
!      source = MYID 	  
      IF(MYID.le.n_r_coll-1) THEN	  !!! MAKING GRID FILE !!! 
      DO i=MYID+1,MYID+1  !! EACH PROCCESOR MAKING ITS OWN R_POINT.
      DO i1=1,n_beta1      !! NEEDS MORE PARALEZATION IN FUTURE
      DO i2=1,n_beta2
      DO i3=1,n_alpha2
      beta =  xbm + xbl*xb(i1)
      bbeta = xbbm + xbbl*xbb(i2)
!      gamma = xgm + xgl*xg(i3)	  
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
!      gamma = xgm + xgl*xg(i3)
      aalpha = xaam + xaal*xaa(i3)
	  else
!	  gamma = xg(i3)
	  aalpha = xaa(i3)
	  end if
! Bikram End.
      R =  R_COM(i)*conv_unit_r	  

      CALL V_POT_DIAT_DIAT(V,R,beta,bbeta,aalpha)
      buffer(i3,i2,i1) = V
     	  
!      IF(MYID.EQ.0) THEN
!      WRITE(2,*) R, buffer(i3,i2,i1)	  
!      ENDIF	  
!      READ(1,*) V_2_2(i3,i2,i1,i)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      task_size = n_beta1*n_beta2*n_alpha2
      IF(MYID.EQ.0) THEN!!! GATHERING ALL BUFFERS TO ROOT PROCESSOR
      V_2_2(:,:,:,1) = buffer 	  
      DO source=1,min(n_r_coll,nproc)-1	  
      CALL MPI_RECV(V_2_2(:,:,:,source+1), task_size, MPI_REAL8, 
     & source, tag2, MPI_COMM_WORLD, status, ierr_mpi)
      ENDDO
      ENDIF	  
      IF(MYID.ne.0 .and. MYID.le.n_r_coll-1 ) THEN
      dest = 0	  
      CALL MPI_SEND(buffer, task_size, MPI_REAL8, dest, 
     &      tag2, MPI_COMM_WORLD, ierr_mpi)
      ENDIF

      IF(MYID.EQ.0 .and.make_grid_file) THEN	 
      OPEN(1,FILE=potential_file_name,ACTION="WRITE",STATUS="NEW")
      WRITE(1,*) coll_type!!! WRITING GRID FILE
      WRITE(1,*) n_beta1,n_beta2,n_alpha2,n_r_coll
      DO i=1,n_r_coll  
      DO i1=1,n_beta1
      DO i2=1,n_beta2
      DO i3=1,n_alpha2
      WRITE(1,*) V_2_2(i3,i2,i1,i)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      CLOSE(1)
      ENDIF	  
      DEALLOCATE(buffer) 
      grid_file_found = .TRUE.		  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      ENDIF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
 	  
!      CALL MPI_BCAST(V_2_2, n_r_coll*task_size, MPI_REAL8,0,
!     &  MPI_COMM_WORLD,ierr_mpi)
!      IF(myid.le.n_r_coll-1) THEN
!      DO i=myid+1,myid+1  
!      DO i1=1,n_beta1
!      DO i2=1,n_beta2
!      DO i3=1,n_gamma1
!      IF(V_2_2(i3,i2,i1,i).ne.buffer(i3,i2,i1)) THEN
!      PRINT*,myid,i3,i2,i1
!      STOP "ERROR IN DATA PASSING"
!      ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO	  
!      ENDIF
      ELSE
      DO i=1,n_r_coll
      DO i1=1,n_beta1
      DO i2=1,n_beta2
      DO i3=1,n_alpha2
      R =  R_COM(i)*conv_unit_r	  
      beta =  xbm + xbl*xb(i1)
      bbeta = xbbm + xbbl*xbb(i2)
!      gamma = xgm + xgl*xg(i3)	  
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
!      gamma = xgm + xgl*xg(i3)
      aalpha = xaam + xaal*xaa(i3)
	  else
!	  gamma = xg(i3)
	  aalpha = xaa(i3)
	  end if
! Bikram End.

      CALL V_POT_DIAT_DIAT(V,R,beta,bbeta,aalpha)
      V_2_2(i3,i2,i1,i) = V
!      READ(1,*) V_2_2(i3,i2,i1,i)
      ENDDO
      ENDDO
      ENDDO
      ENDDO	 
      ENDIF   	  
!      CALL MPI_GATHER(buffer,task_size,MPI_DOUBLE_PRECISION,V_2_2,
!     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
      CASE(6)
!      ALLOCATE(V_VIB_2_2(n_gamma1,n_beta2,n_beta1,n_r_vib2,
!     & n_r_vib1,n_r_coll))
c      ALLOCATE(buffer(n_beta1,n_beta2,n_gamma1))
c      PRINT*,myid
      x2 = pi
      x1 = 0
      xbm = 0.5 * (x2 + x1)
   	  xbl = 0.5 * (x2 - x1)
      ALLOCATE(xb(n_beta1),wb(n_beta1))	  
      CALL gauleg(n_beta1, xb, wb)
      x2 = pi
      x1 = 0
      xbbm = 0.5 * (x2 + x1)
   	  xbbl = 0.5 * (x2 - x1)
      ALLOCATE(xbb(n_beta2),wbb(n_beta2))	  
      CALL gauleg(n_beta2, xbb, wbb)
      x2 = pi*2
      x1 = 0
      xaam = 0.5 * (x2 + x1)
   	  xaal = 0.5 * (x2 - x1)
      ALLOCATE(xaa(n_alpha2),waa(n_alpha2)) 
!      CALL gauleg(n_gamma1, xg, wg)

! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
!      CALL gauleg(n_gamma1, xg, wg)	
      CALL gauleg(n_alpha2, xaa, waa)	
	  else
	  if(myid.eq.0) 
     & print*, "Using Equidistant Grid for Angle Gamma"
	  do i = 1, n_alpha2
	  xaa(i) = x1 + (dble(i) - 0.50d0)*2d0*xaal/n_alpha2
      waa(i)= 2.0d0/n_alpha2
	  end do
	  end if
! Bikram End.

      ALLOCATE(xgg(n_r_vib1),wgg(n_r_vib1))
      ALLOCATE(xg(n_r_vib2),wg(n_r_vib2))	  
!      IF(.not.user_vib_mat_el_defined) THEN											!Dulat: this if condition is corrected
      IF(.not.user_defined) THEN
      x2 = r_vib_diatom_max1
      x1 = r_vib_diatom_min1
      xggm = 0.5 * (x2 + x1)
   	  xggl = 0.5 * (x2 - x1)
      x2 = r_vib_diatom_max2
      x1 = r_vib_diatom_min2
      xgm = 0.5 * (x2 + x1)
   	  xgl = 0.5 * (x2 - x1)	  
      CALL gauleg(n_r_vib1, xgg, wgg)
      CALL gauleg(n_r_vib2, xg, wg)	  
      ELSE
      xgl = 1d0
      xggl = 1d0	  
      ENDIF
      ALLOCATE(r_vb_dt_integrator1(n_r_vib1))
      ALLOCATE(r_vb_dt_integrator2(n_r_vib2))
      ALLOCATE(wf_vib_part_diatom1(n_r_vib1,number_of_channels))
      ALLOCATE(wf_vib_part_diatom2(n_r_vib2,number_of_channels))	  
      IF(.not.user_defined) THEN

! Generation of Gauss-Legendre grid and weights 
! original by Semenov, commented by Bostan 08/08/2022
!      DO i=1,n_r_vib1	  
!      r_vb_dt_integrator1(i) = xggm + xggl*xgg(i)
!      ENDDO	  
!      DO i=1,n_r_vib2	  
!      r_vb_dt_integrator2(i) = xbbm + xbbl*xbb(i)
!      ENDDO
! Generation of equidistant grid and weight

! Dulat Bostan 08/08/2022 step size
	   DO i=1,n_r_vib1	  
        r_vb_dt_integrator1(i)= r_vib_diatom_min1
     & +(DBLE(i)-0.5d0)*2*xggl/n_r_vib1
        wgg(i)= 2.d0/n_r_vib1
       ENDDO
	   
	   DO i=1,n_r_vib2	  
        r_vb_dt_integrator2(i)= r_vib_diatom_min2
     & +(DBLE(i)-0.5d0)*2*xgl/n_r_vib2
        wg(i)= 2.d0/n_r_vib2
       ENDDO
!Dulat Bostan: end
      ENDIF	  
      IF(.not.user_defined) THEN	  
!      morse_de = WE1**2/4d0/XE1/eVtown/autoeV
!      morse_a = (WE1/eVtown/autoeV)*dsqrt(atomic_red_mass1/2d0/morse_de
!     & * amutoau)

! Dulat Bostan 08/08/2022
	  if (morse12_de(1).le.0) then 	! If Morse depth is not specified in the INPUT file
	  if(myid.eq.0) print*, "Molecule #1"
	  if(myid.eq.0) print*, 
     & "Morse Depth is not specified in the input"
	  morse12_de(1) = WE1**2/4d0/XE1/eVtown/autoeV
	  if(myid.eq.0) print*, "Morse Depth: morse_de =", 
     & morse12_de(1), "a.u."
	  else
	  morse12_de(1) = morse12_de(1)/eVtown/autoeV ! Converting Morse depth from cm-1 to a.u.
	  if(myid.eq.0) print*, "Morse Depth: morse_de =", 
     & morse12_de(1), "a.u."
	  endif
	  
	  if (morse12_a(1).le.0) then  ! If Morse width is not specified in the INPUT file
!	  if(myid.eq.0) print*, "Molecule #1"
	  if(myid.eq.0) print*, 
     & "Morse Width is not specified in the input"	  
      morse12_a(1) = (WE1/eVtown/autoeV)*dsqrt(atomic_red_mass1/2d0
     & /morse12_de(1)* amutoau)
	  if(myid.eq.0) print*, "Morse Width: morse_a =", 
     & morse12_a(1), "a.u."
	  else
	  if(myid.eq.0) print*, "Morse Width: morse_a =", 
     & morse12_a(1), "a.u."
	  endif
!	end 08/08/2022

!      PRINT*,"morse_de",morse_de
!      PRINT*,"morse_a",morse_a
!      STOP
      ALLOCATE(wpsir(n_r_vib1))
      DO k=1,number_of_channels
      vw = v1_ch(k)	  
!      PRINT*,"OLD",r_vb_dt_inegrator(1)
      CALL VIBRATIONAL_WAVEFUNCTION											! Dulat: changed the call arguments
     & (wpsir,r_vb_dt_integrator1,n_r_vib1,
     & vw,atomic_red_mass1*amutoau,
     &  morse12_de(1),morse12_a(1),morse12_re(1))
      wf_vib_part_diatom1(1:n_r_vib1,k) = wpsir(1:n_r_vib1)
!      PRINT*,"NEW",r_vb_dt_inegrator(1)	  
      ENDDO
!      morse_de = WE2**2/4d0/XE2/eVtown/autoeV
!      morse_a = (WE2/eVtown/autoeV)*dsqrt(atomic_red_mass2/2d0/morse_de
!     & * amutoau)

! Dulat Bostan 08/08/2022
	  if (morse12_de(2).le.0) then 	! If Morse depth is not specified in the INPUT file
	  if(myid.eq.0) print*, "Molecule #2"
	  if(myid.eq.0) print*, 
     & "Morse Depth is not specified in the input"
	  morse12_de(2) = WE2**2/4d0/XE2/eVtown/autoeV
	  if(myid.eq.0) print*, "Morse Depth: morse_de =", 
     & morse12_de(2), "a.u."
	  else
	  morse12_de(2) = morse12_de(2)/eVtown/autoeV ! Converting Morse depth from cm-1 to a.u.
      print*, "Morse Depth: morse_de =", morse12_de(2), "a.u."
	  endif
	  
	  if (morse12_a(2).le.0) then  ! If Morse width is not specified in the INPUT file
!	  print*, "Molecule #2"
	  if(myid.eq.0) print*, 
     & "Morse Width is not specified in the input"	  
      morse12_a(2) = (WE2/eVtown/autoeV)*dsqrt(atomic_red_mass2/2d0
     & /morse12_de(2)* amutoau)
	  if(myid.eq.0) print*, "Morse Width: morse_a =", 
     & morse12_a(2), "a.u."
	  else
	  if(myid.eq.0) print*, "Morse Width: morse_a =", 
     & morse12_a(2), "a.u."
	  endif      
!	end 08/08/2022

      DEALLOCATE(wpsir)
      ALLOCATE(wpsir(n_r_vib2))      	  
      DO k=1,number_of_channels
      vw = v2_ch(k)	  
!      PRINT*,"OLD",r_vb_dt_inegrator(1)	  
      CALL VIBRATIONAL_WAVEFUNCTION													! Dulat: changed the call arguments
     & (wpsir,r_vb_dt_integrator2,n_r_vib2,
     & vw,atomic_red_mass2*amutoau,
     &  morse12_de(2),morse12_a(2),morse12_re(2))
      wf_vib_part_diatom2(1:n_r_vib2,k) = wpsir(1:n_r_vib2)
!      PRINT*,"NEW",r_vb_dt_inegrator(1)	  
      ENDDO
      DEALLOCATE(wpsir) 	  
!      PRINT*,"NEW",r_vb_dt_inegrator(1)
!      STOP	  
!      DEALLOCATE(wpsir)
      OPEN(456,FILE="WAVEFUNCTIONS.out",
     & ACTION="WRITE")
c      ! WRITTING A GRID	 
      WRITE(456,*) "R1(i),J1(i)" ! A HEADER
      DO i=1,n_r_vib1! LOOP OVER VIBRATIONAL GRID	  
      WRITE(456,*)r_vb_dt_integrator1(i),wgg(i)*xggl								!Dulat: multiplied by xggl
	  ! r_vb_dt_inegrator(i)-
! r-value  for i-point of the grid, wg(i) - weight(jacobian)
      ENDDO
      WRITE(456,*) "R2(i),J2(i)" ! A HEADER
      DO i=1,n_r_vib2! LOOP OVER VIBRATIONAL GRID	  
      WRITE(456,*)r_vb_dt_integrator2(i),wg(i)*xgl									!Dulat: multiplied by xgl
	  ! r_vb_dt_inegrator(i)-
! r-value  for i-point of the grid, wg(i) - weight(jacobian)
      ENDDO	  
      WRITE(456,*) "WAVEFUNCTIONS ARE LISTED BELOW"	  
      DO k=1,number_of_channels !! k - channel number
c      !! WRITING THE WAVEFUNCTIONS	  
      WRITE(456,*) "MOLECULE#1,CHANNEL=",k !!! A HEADER
      DO i=1,n_r_vib1	  
      WRITE(456,*)wf_vib_part_diatom1(i,k) ! VALUE OF WAVEFUNCTION FOR k-channel at i-point of the grid
      ENDDO
      WRITE(456,*) !! LEAVE IT BLANK	  
      ENDDO

      DO k=1,number_of_channels !! k - channel number
c      !! WRITING THE WAVEFUNCTIONS	  
      WRITE(456,*) "MOLECULE#2,CHANNEL=",k !!! A HEADER
      DO i=1,n_r_vib2	  
      WRITE(456,*)wf_vib_part_diatom2(i,k) ! VALUE OF WAVEFUNCTION FOR k-channel at i-point of the grid
      ENDDO
      WRITE(456,*) !! LEAVE IT BLANK	  
      ENDDO	  
      CLOSE(456)
      ELSE
      r_vb_dt_integrator1 = r_grid_vib1
      r_vb_dt_integrator2 = r_grid_vib2
      wgg = jacobian_vib1
      wg = jacobian_vib2
      wf_vib_part_diatom1 = vibrational_wavefunctions_1
      wf_vib_part_diatom2 = vibrational_wavefunctions_2	  
      ENDIF	 	  
	  

      IF(.NOT.make_grid_file) EXIT CASE_FORM	  
      IF(grid_file_found) EXIT	  
      ALLOCATE(V_VIB_2_2_int_buffer(n_alpha2,n_beta2,n_beta1,n_r_vib2,
     & n_r_vib1))
      tag2= 2
      V_VIB_2_2_int_buffer = 0d0
!      source = MYID 	  
      IF(MYID.le.n_r_coll-1) THEN	  
      DO i=MYID+1,MYID+1
      DO i_vib1 = 1,n_r_vib1
      DO i_vib2 = 1,n_r_vib2	  
      DO i1=1,n_beta1
      DO i2=1,n_beta2
      DO i3=1,n_alpha2
      beta =  xbm + xbl*xb(i1)
      bbeta = xbbm + xbbl*xbb(i2)
!      gamma = xgm + xgl*xg(i3)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
!      gamma = xgm + xgl*xg(i3)
      aalpha = xaam + xaal*xaa(i3)
	  else
!	  gamma = xg(i3)
	  aalpha = xaa(i3)
	  end if
! Bikram End.
      R =  R_COM(i)*conv_unit_r	  
      r_vib1 = r_vb_dt_integrator1(i_vib1)*conv_unit_r
      r_vib2 = r_vb_dt_integrator2(i_vib2)*conv_unit_r
	  
      CALL V_POT_VIB_DIAT_DIAT(V,R,r_vib1,r_vib2,
     & beta,bbeta,aalpha)
      V_VIB_2_2_int_buffer(i3,i2,i1,i_vib2,
     & i_vib1) = V
     	  
!      IF(MYID.EQ.0) THEN
!      WRITE(2,*) R, buffer(i3,i2,i1)	  
!      ENDIF	  
!      READ(1,*) V_2_2(i3,i2,i1,i)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO	  
      ENDIF
      task_size = n_beta1*n_beta2*n_alpha2*n_r_vib1*n_r_vib2
      IF(MYID.ne.0 .and. MYID.le.n_r_coll-1 ) THEN
      dest = 0	  
      CALL MPI_SEND(V_VIB_2_2_int_buffer, task_size, MPI_REAL8, dest, 
     &      tag2, MPI_COMM_WORLD, ierr_mpi)
      ENDIF
      IF(MYID.EQ.0 .and. make_grid_file) THEN
      OPEN(1,FILE=potential_file_name,ACTION="WRITE",STATUS="NEW")
      WRITE(1,*) coll_type
      WRITE(1,*) n_beta1,n_beta2,n_alpha2,n_r_vib1,
     & n_r_vib2, n_r_coll	  
      	  
      DO i=1,min(n_r_coll,nproc)
      IF(i.gt.1) THEN	  
      CALL MPI_RECV(V_VIB_2_2_int_buffer, task_size, MPI_REAL8, 
     & i-1, tag2, MPI_COMM_WORLD, status, ierr_mpi)
      ENDIF

      V = V_VIB_2_2_int_buffer(i3,i2,i1,i_vib2,
     & i_vib1)
      WRITE(1,*) V	 

      ENDDO
	  
      CLOSE(1)      	  
      ENDIF	  
      IF(myid.le.n_r_coll-1) DEALLOCATE(V_VIB_2_2_int_buffer) 
      grid_file_found = .TRUE.		  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      CASE(7)
	  !      IF(MYID.EQ.0) ALLOCATE(V_3_3
!     & (n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,n_r_coll))

!      buffer_3_3 = 0d0	 
      x2 = pi*2
      x1 = 0
      xam = 0.5 * (x2 + x1)
   	  xal = 0.5 * (x2 - x1)
      IF(n_alpha1.le.0) n_alpha1 = n_alpha2  	  
      ALLOCATE(xa(n_alpha1),wa(n_alpha1))
!      CALL gauleg(n_alpha1, xa, wa)		  
      x2 = pi*2
      x1 = 0
      xaam = 0.5 * (x2 + x1)
   	  xaal = 0.5 * (x2 - x1)
      IF(n_alpha2.le.0) n_alpha2 = n_alpha1  	  
      ALLOCATE(xaa(n_alpha2),waa(n_alpha2)) 
!      CALL gauleg(n_alpha2, xaa, waa)	  	  

! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      CALL gauleg(n_alpha1, xa, wa)	
      CALL gauleg(n_alpha2, xaa, waa)
	  else
	  if(myid.eq.0) 
     & print*, "Using Equidistant Grid for Angle Alpha"
	  do i = 1, n_alpha1
	  xa(i) = x1 + (dble(i) - 0.50d0)*2d0*xal/n_alpha1
      wa(i)= 2.0d0/n_alpha1
	  end do
	  do i = 1, n_alpha2
	  xaa(i) = x1 + (dble(i) - 0.50d0)*2d0*xaal/n_alpha2
      waa(i)= 2.0d0/n_alpha2
	  end do
	  end if
! Bikram End.
	  
      x2 = pi
      x1 = 0
      xbm = 0.5 * (x2 + x1)
   	  xbl = 0.5 * (x2 - x1)
      xbbm = 0.5 * (x2 + x1)
   	  xbbl = 0.5 * (x2 - x1)		  
      ALLOCATE(xb(n_beta1),wb(n_beta1))
      CALL gauleg(n_beta1, xb, wb)	
      ALLOCATE(xbb(n_beta2),wbb(n_beta2))
      CALL gauleg(n_beta2, xbb, wbb)			  
      x2 = pi*2
      x1 = 0
      xgm = 0.5 * (x2 + x1)
   	  xgl = 0.5 * (x2 - x1)
      xggm = 0.5 * (x2 + x1)
   	  xggl = 0.5 * (x2 - x1)	  
      ALLOCATE(xg(n_gamma1),wg(n_gamma1))  
!      CALL gauleg(n_gamma1, xg, wg)	
!      ALLOCATE(xgg(n_gamma2),wgg(n_gamma2))
!      CALL gauleg(n_gamma2, xgg, wgg)

! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      CALL gauleg(n_gamma1, xg, wg)	
	  else
	  if(myid.eq.0) 
     & print*, "Using Equidistant Grid for Angle Gamma"
	  do i = 1, n_gamma1
	  xg(i) = x1 + (dble(i) - 0.50d0)*2d0*xgl/n_gamma1
      wg(i)= 2.0d0/n_gamma1
	  end do
	  end if
! Bikram End.
	  
      IF(.NOT.make_grid_file) THEN
	  ALLOCATE(V_3_2(n_gamma1,n_beta2,n_beta1,n_alpha2,n_r_coll))
      DO i=1,n_r_coll
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      beta = xbm + xbl*xb(i3)
      bbeta = xbbm + xbbl*xbb(i4)
!      gamma = xgm + xgl*xg(i5)
!      aalpha =  xaam + xaal*xaa(i2)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
	  else
	  gamma = xg(i5)
	  aalpha = xaa(i2)
	  end if
! Bikram End.
      R =  R_COM(i)*conv_unit_r	   
	  
      CALL V_POT_SYM_DIATOM(V,R,0d0,beta,gamma,
     & aalpha,bbeta)
      V_3_2(i5,i4,i3,i2,i) = V	  
!      READ(1,*) V_2_2(i3,i2,i1,i)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO	 	  
      ENDIF      
      IF(.NOT.make_grid_file) EXIT CASE_FORM
      IF(grid_file_found) EXIT CASE_FORM
	 !!! CHECK
!      EXIT CASE_FORM	 
    !!	 
      IF(myid.le.n_r_coll-1)
     & ALLOCATE(V_3_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      IF(myid.le.n_r_coll-1) THEN	  
      DO i=myid+1,myid+1
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      beta = xbm + xbl*xb(i3)
      bbeta = xbbm + xbbl*xbb(i4)
!      gamma = xgm + xgl*xg(i5)
!      aalpha =  xaam + xaal*xaa(i2)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
	  else
	  gamma = xg(i5)
	  aalpha = xaa(i2)
	  end if
! Bikram End.
      R =  R_COM(i)*conv_unit_r	  
	  
      CALL V_POT_SYM_DIATOM(V,R,0d0,beta,gamma,
     & aalpha,bbeta)	 
!      V_3_3(i6,i5,i4,i3,i2,i) = V
      V_3_2_int_buffer(i5,i4,i3,i2) = V
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      task_size = n_alpha2*n_beta1*n_beta2*n_gamma1
      IF(MYID.ne.0 .and. MYID.le.n_r_coll-1 ) THEN
      dest = 0
      tag2 = 2	  
      CALL MPI_SEND(V_3_2_int_buffer, task_size, MPI_REAL8, dest, 
     &      tag2, MPI_COMM_WORLD, ierr_mpi)
	  print *, 'sent', myid
      ENDIF
      IF(MYID.EQ.0 .and.make_grid_file) THEN	 
      OPEN(1,FILE=potential_file_name,ACTION="WRITE",STATUS="NEW"
     & ,FORM="UNFORMATTED")
      WRITE(1) coll_type
      WRITE(1) n_alpha2,n_beta1,n_beta2,n_gamma1,n_r_coll
      DO i=1,min(n_r_coll,nproc)
      IF(i.gt.1) THEN
      tag2 = 2
      V_3_2_int_buffer = 0d0	  
      CALL MPI_RECV(V_3_2_int_buffer, task_size, MPI_REAL8, 
     & i-1, tag2, MPI_COMM_WORLD, status, ierr_mpi)	  
	  print *, 'rcvd', i-1
      ENDIF	  


      WRITE(1)V_3_2_int_buffer
      ENDDO

      CLOSE(1)
      ENDIF 
      grid_file_found = .TRUE.
      IF(myid.le.n_r_coll-1)        DEALLOCATE(V_3_2_int_buffer)		  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )		  
      CASE(8)
	  !      IF(MYID.EQ.0) ALLOCATE(V_3_3
!     & (n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,n_r_coll))

!      buffer_3_3 = 0d0	 
      x2 = pi*2
      x1 = 0
      xam = 0.5 * (x2 + x1)
   	  xal = 0.5 * (x2 - x1)
      IF(n_alpha1.le.0) n_alpha1 = n_alpha2  	  
      ALLOCATE(xa(n_alpha1),wa(n_alpha1))
!      CALL gauleg(n_alpha1, xa, wa)		  
      x2 = pi*2
      x1 = 0
      xaam = 0.5 * (x2 + x1)
   	  xaal = 0.5 * (x2 - x1)
      IF(n_alpha2.le.0) n_alpha2 = n_alpha1  	  
      ALLOCATE(xaa(n_alpha2),waa(n_alpha2))
!      CALL gauleg(n_alpha2, xaa, waa)	  

! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      CALL gauleg(n_alpha1, xa, wa)	
      CALL gauleg(n_alpha2, xaa, waa)
	  else
	  if(myid.eq.0) 
     & print*, "Using Equidistant Grid for Angle Alpha"
	  do i = 1, n_alpha1
	  xa(i) = x1 + (dble(i) - 0.50d0)*2d0*xal/n_alpha1
      wa(i)= 2.0d0/n_alpha1
	  end do
	  do i = 1, n_alpha2
	  xaa(i) = x1 + (dble(i) - 0.50d0)*2d0*xaal/n_alpha2
      waa(i)= 2.0d0/n_alpha2
	  end do
	  end if
! Bikram End.
	  
      x2 = pi
      x1 = 0
      xbm = 0.5 * (x2 + x1)
   	  xbl = 0.5 * (x2 - x1)
      xbbm = 0.5 * (x2 + x1)
   	  xbbl = 0.5 * (x2 - x1)		  
      ALLOCATE(xb(n_beta1),wb(n_beta1))
      CALL gauleg(n_beta1, xb, wb)	
      ALLOCATE(xbb(n_beta2),wbb(n_beta2))
      CALL gauleg(n_beta2, xbb, wbb)			  
      x2 = pi*2
      x1 = 0
      xgm = 0.5 * (x2 + x1)
   	  xgl = 0.5 * (x2 - x1)
      xggm = 0.5 * (x2 + x1)
   	  xggl = 0.5 * (x2 - x1)	  
      ALLOCATE(xg(n_gamma1),wg(n_gamma1))
!      CALL gauleg(n_gamma1, xg, wg)	
!      ALLOCATE(xgg(n_gamma2),wgg(n_gamma2))
!      CALL gauleg(n_gamma2, xgg, wgg)

! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      CALL gauleg(n_gamma1, xg, wg)	
	  else
	  if(myid.eq.0) 
     & print*, "Using Equidistant Grid for Angle Gamma"
	  do i = 1, n_gamma1
	  xg(i) = x1 + (dble(i) - 0.50d0)*2d0*xgl/n_gamma1
      wg(i)= 2.0d0/n_gamma1
	  end do
	  end if
! Bikram End.
	  
      IF(.NOT.make_grid_file) THEN
	  ALLOCATE(V_3_2(n_gamma1,n_beta2,n_beta1,n_alpha2,n_r_coll))
      DO i=1,n_r_coll
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      beta = xbm + xbl*xb(i3)
      bbeta = xbbm + xbbl*xbb(i4)
!      gamma = xgm + xgl*xg(i5)
!      aalpha =  xaam + xaal*xaa(i2)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
	  else
	  gamma = xg(i5)
	  aalpha = xaa(i2)
	  end if
! Bikram End.
      R =  R_COM(i)*conv_unit_r	   
	  
      CALL V_POT_ASYM_DIATOM(V,R,0d0,beta,gamma,
     & aalpha,bbeta)
      V_3_2(i5,i4,i3,i2,i) = V	  
!      READ(1,*) V_2_2(i3,i2,i1,i)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO	 	  
      ENDIF  
      IF(.NOT.make_grid_file) EXIT CASE_FORM
      IF(grid_file_found) EXIT CASE_FORM
	 !!! CHECK
!      EXIT CASE_FORM	 
    !!	 
      IF(myid.le.n_r_coll-1)
     & ALLOCATE(V_3_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      IF(myid.le.n_r_coll-1) THEN	  
      DO i=myid+1,myid+1
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      beta = xbm + xbl*xb(i3)
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
      bbeta = xbbm + xbbl*xbb(i4)
      R =  R_COM(i)*conv_unit_r	  
      CALL V_POT_ASYM_DIATOM(V,R,0d0,beta,gamma,
     & aalpha,bbeta)	 
!      V_3_3(i6,i5,i4,i3,i2,i) = V
      V_3_2_int_buffer(i5,i4,i3,i2) = V
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      task_size = n_alpha2*n_beta1*n_beta2*n_gamma1
      IF(MYID.ne.0 .and. MYID.le.n_r_coll-1 ) THEN
      dest = 0
      tag2 = 2	  
      CALL MPI_SEND(V_3_2_int_buffer, task_size, MPI_REAL8, dest, 
     &      tag2, MPI_COMM_WORLD, ierr_mpi)
      ENDIF
      IF(MYID.EQ.0 .and.make_grid_file) THEN	 
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",ACTION="WRITE",STATUS="NEW")
      WRITE(1) coll_type
      WRITE(1) n_alpha2,n_beta1,n_beta2,n_gamma1,n_r_coll
      DO i=1,min(n_r_coll,nproc)
      IF(i.gt.1) THEN
      tag2 = 2
      V_3_2_int_buffer = 0d0	  
      CALL MPI_RECV(V_3_2_int_buffer, task_size, MPI_REAL8, 
     & i-1, tag2, MPI_COMM_WORLD, status, ierr_mpi)	  
      ENDIF	  
  

      WRITE(1)V_3_2_int_buffer

      ENDDO
      CLOSE(1)
      ENDIF 
      grid_file_found = .TRUE.
      IF(myid.le.n_r_coll-1)  DEALLOCATE(V_3_2_int_buffer)	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )		  
      CASE(9)
	  !      IF(MYID.EQ.0) ALLOCATE(V_3_3
!     & (n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,n_r_coll))

!      buffer_3_3 = 0d0	 
      x2 = pi*2
      x1 = 0
      xam = 0.5 * (x2 + x1)
   	  xal = 0.5 * (x2 - x1)
      IF(n_alpha1.le.0) n_alpha1 = n_alpha2  	  
      ALLOCATE(xa(n_alpha1),wa(n_alpha1))
!      CALL gauleg(n_alpha1, xa, wa)		  
      x2 = pi*2
      x1 = 0
      xaam = 0.5 * (x2 + x1)
   	  xaal = 0.5 * (x2 - x1)
      IF(n_alpha2.le.0) n_alpha2 = n_alpha1  	  
      ALLOCATE(xaa(n_alpha2),waa(n_alpha2))
!      CALL gauleg(n_alpha2, xaa, waa)	   

! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      CALL gauleg(n_alpha1, xa, wa)	
      CALL gauleg(n_alpha2, xaa, waa)
	  else
	  if(myid.eq.0) 
     & print*, "Using Equidistant Grid for Angle Alpha"
	  do i = 1, n_alpha1
	  xa(i) = x1 + (dble(i) - 0.50d0)*2d0*xal/n_alpha1
      wa(i)= 2.0d0/n_alpha1
	  end do
	  do i = 1, n_alpha2
	  xaa(i) = x1 + (dble(i) - 0.50d0)*2d0*xaal/n_alpha2
      waa(i)= 2.0d0/n_alpha2
	  end do
	  end if
! Bikram End.
	  
      x2 = pi
      x1 = 0
      xbm = 0.5 * (x2 + x1)
   	  xbl = 0.5 * (x2 - x1)
      xbbm = 0.5 * (x2 + x1)
   	  xbbl = 0.5 * (x2 - x1)		  
      ALLOCATE(xb(n_beta1),wb(n_beta1))
      CALL gauleg(n_beta1, xb, wb)	
      ALLOCATE(xbb(n_beta2),wbb(n_beta2))
      CALL gauleg(n_beta2, xbb, wbb)			  
      x2 = pi*2
      x1 = 0
      xgm = 0.5 * (x2 + x1)
   	  xgl = 0.5 * (x2 - x1)
      xggm = 0.5 * (x2 + x1)
   	  xggl = 0.5 * (x2 - x1)	  
      ALLOCATE(xg(n_gamma1),wg(n_gamma1))
!      CALL gauleg(n_gamma1, xg, wg)	
      ALLOCATE(xgg(n_gamma2),wgg(n_gamma2))
!      CALL gauleg(n_gamma2, xgg, wgg)

! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      CALL gauleg(n_gamma1, xg, wg)	
      CALL gauleg(n_gamma2, xgg, wgg)
	  else
	  if(myid.eq.0) 
     & print*, "Using Equidistant Grid for Angle Gamma"
	  do i = 1, n_gamma1
	  xg(i) = x1 + (dble(i) - 0.50d0)*2d0*xgl/n_gamma1
      wg(i)= 2.0d0/n_gamma1
	  end do
	  do i = 1, n_gamma2
	  xgg(i) = x1 + (dble(i) - 0.50d0)*2d0*xggl/n_gamma2
      wgg(i)= 2.0d0/n_gamma2
	  end do
	  end if
! Bikram End.
	  
      IF(.NOT.make_grid_file) EXIT CASE_FORM
      IF(grid_file_found) EXIT CASE_FORM
!!   CREATION OF GRID FILE 
      IF(myid.le.n_r_coll-1)
     & ALLOCATE(buffer_3_3(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      IF(myid.le.n_r_coll-1) THEN	  
      DO i=myid+1,myid+1
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      DO i6=1,n_gamma2
      beta = xbm + xbl*xb(i3)
      bbeta = xbbm + xbbl*xbb(i4)
!      gamma = xgm + xgl*xg(i5)
!      aalpha =  xaam + xaal*xaa(i2)
!      ggamma = xggm + xggl*xgg(i6)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
      ggamma = xggm + xggl*xgg(i6)
	  else
	  gamma = xg(i5)
	  aalpha = xaa(i2)
	  ggamma = xgg(i6)
	  end if
! Bikram End.
      R =  R_COM(i)*conv_unit_r	 
	  
      CALL V_POT_ASYM_ASYM(V,R,0d0,beta,gamma,
     & aalpha,bbeta,ggamma)	 
!      V_3_3(i6,i5,i4,i3,i2,i) = V
      buffer_3_3(i6,i5,i4,i3,i2) = V
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      task_size = n_alpha2*n_beta1*n_beta2*n_gamma1*n_gamma2
      IF(MYID.ne.0 .and. MYID.le.n_r_coll-1 ) THEN
      dest = 0
      tag2 = 2	  
      CALL MPI_SEND(buffer_3_3, task_size, MPI_REAL8, dest, 
     &      tag2, MPI_COMM_WORLD, ierr_mpi)
      ENDIF
      IF(MYID.EQ.0 .and.make_grid_file) THEN	 
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",ACTION="WRITE",STATUS="NEW")
      WRITE(1) coll_type
      WRITE(1) n_alpha2,n_beta1,n_beta2,n_gamma1,n_gamma2,n_r_coll
      DO i=1,min(n_r_coll,nproc)
      IF(i.gt.1) THEN
      tag2 = 2
      buffer_3_3 = 0d0	  
      CALL MPI_RECV(buffer_3_3, task_size, MPI_REAL8, 
     & i-1, tag2, MPI_COMM_WORLD, status, ierr_mpi)	  
      ENDIF	  
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      DO i6=1,n_gamma2

      WRITE(1)buffer_3_3(i6,i5,i4,i3,i2)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      CLOSE(1)
      ENDIF	  
      IF(myid.le.n_r_coll-1)        DEALLOCATE(buffer_3_3) 
      grid_file_found = .TRUE.		  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	
      CASE(0)
!      IF(MYID.EQ.0) ALLOCATE(V_3_3
!     & (n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,n_r_coll))

!      buffer_3_3 = 0d0	 
      x2 = pi*2
      x1 = 0
      xam = 0.5 * (x2 + x1)
   	  xal = 0.5 * (x2 - x1)
      IF(n_alpha1.le.0) n_alpha1 = n_alpha2 
      IF(n_alpha2.le.0) STOP "CHECK INPUT,a2" 	  
      ALLOCATE(xa(n_alpha1),wa(n_alpha1))
!      CALL gauleg(n_alpha1, xa, wa)		  
      x2 = pi*2
      x1 = 0
      xaam = 0.5 * (x2 + x1)
   	  xaal = 0.5 * (x2 - x1)
      IF(n_alpha2.le.0) n_alpha2 = n_alpha1
      IF(n_alpha1.le.0) STOP "CHECK INPUT,a1"  	  
      ALLOCATE(xaa(n_alpha2),waa(n_alpha2))
!      CALL gauleg(n_alpha2, xaa, waa)	  

! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      CALL gauleg(n_alpha1, xa, wa)	
      CALL gauleg(n_alpha2, xaa, waa)
	  else
	  if(myid.eq.0) 
     & print*, "Using Equidistant Grid for Angle Alpha"
	  do i = 1, n_alpha1
	  xa(i) = x1 + (dble(i) - 0.50d0)*2d0*xal/n_alpha1
      wa(i)= 2.0d0/n_alpha1
	  end do
	  do i = 1, n_alpha2
	  xaa(i) = x1 + (dble(i) - 0.50d0)*2d0*xaal/n_alpha2
      waa(i)= 2.0d0/n_alpha2
	  end do
	  end if
! Bikram End.
	  
      x2 = pi
      x1 = 0
      xbm = 0.5 * (x2 + x1)
   	  xbl = 0.5 * (x2 - x1)
      xbbm = 0.5 * (x2 + x1)
   	  xbbl = 0.5 * (x2 - x1)
      IF(n_beta1.le.0) STOP "CHECK INPUT,b1"
      IF(n_beta2.le.0) STOP "CHECK INPUT,b2"	  
      ALLOCATE(xb(n_beta1),wb(n_beta1))
      CALL gauleg(n_beta1, xb, wb)	
      ALLOCATE(xbb(n_beta2),wbb(n_beta2))
      CALL gauleg(n_beta2, xbb, wbb)			  
      x2 = pi*2
      x1 = 0
      xgm = 0.5 * (x2 + x1)
   	  xgl = 0.5 * (x2 - x1)
      xggm = 0.5 * (x2 + x1)
   	  xggl = 0.5 * (x2 - x1)
      IF(n_gamma1.le.0) STOP "CHECK INPUT,g1"
      IF(n_gamma2.le.0) STOP "CHECK INPUT,g2"	  
      ALLOCATE(xg(n_gamma1),wg(n_gamma1))
!      CALL gauleg(n_gamma1, xg, wg)	
      ALLOCATE(xgg(n_gamma2),wgg(n_gamma2))
!      CALL gauleg(n_gamma2, xgg, wgg)

! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      CALL gauleg(n_gamma1, xg, wg)	
      CALL gauleg(n_gamma2, xgg, wgg)
	  else
	  if(myid.eq.0) 
     & print*, "Using Equidistant Grid for Angle Gamma"
	  do i = 1, n_gamma1
	  xg(i) = x1 + (dble(i) - 0.50d0)*2d0*xgl/n_gamma1
      wg(i)= 2.0d0/n_gamma1
	  end do
	  do i = 1, n_gamma2
	  xgg(i) = x1 + (dble(i) - 0.50d0)*2d0*xggl/n_gamma2
      wgg(i)= 2.0d0/n_gamma2
	  end do
	  end if
! Bikram End.
	  
      IF(.NOT.make_grid_file) EXIT CASE_FORM
      IF(grid_file_found) EXIT CASE_FORM
!!  CREATION OF GRID FILE
      IF(myid.le.n_r_coll-1) THEN
      ALLOCATE(buffer_3_3(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2))
!      PRINT*,"buffer_has_beenn_allocated for i_r",myid+1	  
      ENDIF	  
      IF(myid.le.n_r_coll-1) THEN	  
      DO i=myid+1,myid+1
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      DO i6=1,n_gamma2
      beta = xbm + xbl*xb(i3)
      bbeta = xbbm + xbbl*xbb(i4)
!      gamma = xgm + xgl*xg(i5)
!      aalpha =  xaam + xaal*xaa(i2)
!      ggamma = xggm + xggl*xgg(i6)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
      ggamma = xggm + xggl*xgg(i6)
	  else
	  gamma = xg(i5)
	  aalpha = xaa(i2)
	  ggamma = xgg(i6)
	  end if
! Bikram End.
      R =  R_COM(i)*conv_unit_r
	  
      CALL V_POT_ASYM_ASYM(V,R,0d0,beta,gamma,
     & aalpha,bbeta,ggamma)
!      CALL POTENTIAL(V,R,0d0,0d0,beta,gamma,
!     & aalpha,bbeta,ggamma)	  
!      V_3_3(i6,i5,i4,i3,i2,i) = V
      buffer_3_3(i6,i5,i4,i3,i2) = V
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      task_size = n_alpha2*n_beta1*n_beta2*n_gamma1*n_gamma2
      IF(MYID.ne.0 .and. MYID.le.n_r_coll-1 ) THEN
      dest = 0
      tag2 = 2	  
      CALL MPI_SEND(buffer_3_3, task_size, MPI_REAL8, dest, 
     &      tag2, MPI_COMM_WORLD, ierr_mpi)
      ENDIF
      IF(MYID.EQ.0 .and.make_grid_file) THEN	 
      OPEN(1,FILE=potential_file_name,FORM="UNFORMATTED",
     & ACTION="WRITE",STATUS="NEW")
      WRITE(1) coll_type
      WRITE(1) n_alpha2,n_beta1,n_beta2,n_gamma1,n_gamma2,n_r_coll
      DO i=1,min(n_r_coll,nproc)
      IF(i.gt.1) THEN
      tag2 = 2
      buffer_3_3 = 0d0	  
      CALL MPI_RECV(buffer_3_3, task_size, MPI_REAL8, 
     & i-1, tag2, MPI_COMM_WORLD, status, ierr_mpi)	  
      ENDIF	  

      WRITE(1)buffer_3_3
      ENDDO
      CLOSE(1)
      ENDIF	  
      IF(myid.le.n_r_coll-1)        DEALLOCATE(buffer_3_3) 
      grid_file_found = .TRUE.		  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	
      END SELECT
      ENDDO CASE_FORM
      GAUSS_LEGENDRE_GRID_DEFINED = .TRUE.
      IF(MYID.EQ.0) PRINT*,"GAUSS-LEGENDRE GRID DONE"
      END SUBROUTINE INI_POT_STORAGE
	  
      SUBROUTINE INTEGRATOR_TEMP(intgrlr,k,i_r_point, tmp1)
      USE POT_STORE
      USE GRID_INTEGR	  
      USE VARIABLES
      USE CONSTANTS
      USE MPI_DATA	  
      IMPLICIT NONE
      INTEGER st_1,st_2,parity_ident,i_r_point,k	  
      INTEGER, PARAMETER :: n_g_p = 30
      REAL*8 x1,x2,x,intgrlr, intgrli,integrant
      INTEGER i1,i2,i3,i4,i5,i6,j1_t,k1_t,ka1_t,eps1_t,kc1_t,v1_t,m1_t
      INTEGER j2_t,k2_t,ka2_t,eps2_t,kc2_t,v2_t,m2_t
      INTEGER j12_1_t,j12_2_t,j1_1_t,j1_2_t,j2_1_t,j2_2_t,
     & p_1,p_2,l_1,l_2	  
      INTEGER chann_1,chann_2,i_vib1,i_vib2
      INTEGER i_spin_count, tmp1
      REAL*8 dr1,dr2,di1,di2,r_vib1,r_vib2
      REAL*8, ALLOCATABLE :: wfr1_fine(:),wfi1_fine(:)
      REAL*8, ALLOCATABLE :: wfr2_fine(:),wfi2_fine(:)	  
      REAL*8 V,R,V1
      intgrlr = 0d0
      intgrli = 0d0
      st_1 = ind_mat(1,k)
      st_2 = ind_mat(2,k)	  
      conv_unit_r = 1d0	  
      IF(angs_unit) conv_unit_r = a_bohr	  
      SELECT CASE(coll_type)
      CASE(1)
      IF(.NOT. expansion_defined) THEN 	  
      IF(.not.GAUSS_LEGENDRE_GRID_DEFINED)
     & CALL INI_POT_STORAGE
      chann_1 = indx_chann(st_1)
      chann_2 = indx_chann(st_2)
      j1_t = j_ch(chann_1)
      j2_t = j_ch(chann_2)
      m1_t = m12(st_1)
      m2_t = m12(st_2)

      IF(.not.fine_structure_defined) THEN
      DO i2=1,n_beta	  
      beta = xbm + xbl*xb(i2)	  
      CALL WF_DIAT_TOP(j1_t,m1_t,beta,0d0,dr1,di1)
      CALL WF_DIAT_TOP(j2_t,m2_t,beta,0d0,dr2,di2)
c      R =  R_COM(i)*conv_unit_r
c      beta = xbm + xbl*xb(i2)	  
c      CALL V_POT_DIATOM_ATOM(V,R,beta)	  
      V = V_2_1(i2,i_r_point)	  
      integrant	 = (dr1*dr2+di1*di2)
     & *wb(i2)*dsin(beta)*V
      intgrlr  = intgrlr+integrant*xbl*2d0*pi
      ENDDO
      ELSE
      IF(LORB_FINE.ne.0) STOP "ERROR: ONLY EXPANSION ACCEPTABLE"
      IF(SPIN_FINE.ne.3) STOP "ERROR: MULTIPLICITY 3 ACCEPTABLE"	  
	  ALLOCATE(wfr1_fine(SPIN_FINE+1))
	  ALLOCATE(wfi1_fine(SPIN_FINE+1))
	  ALLOCATE(wfr2_fine(SPIN_FINE+1))
	  ALLOCATE(wfi2_fine(SPIN_FINE+1))
!	  write(*,*) n_beta,"n_beta"
      DO i2=1,n_beta
      beta = xbm + xbl*xb(i2)	  	  
      CALL WF_DIAT_TOP_FINE(wfr1_fine,wfi1_fine,
     & (SPIN_FINE-1)/2,st_1,beta,0d0,m1_t)

      CALL WF_DIAT_TOP_FINE(wfr2_fine,wfi2_fine,
     & (SPIN_FINE-1)/2,st_2,beta,0d0,m2_t)		 
!      R =  R_COM(i)*conv_unit_r !!!!!!!! CAREFULL HERE
!      beta = xbm + xbl*xb(i2)	  
!      CALL V_POT_DIATOM_ATOM(V,R,beta)	  
      V = V_2_1(i2,i_r_point)
	   integrant = 0d0 
      DO i_spin_count=1,SPIN_FINE*2+1
      dr1 = wfr1_fine(i_spin_count)
      dr2 = wfr2_fine(i_spin_count)
      di1 = wfi1_fine(i_spin_count)
      di2 = wfi2_fine(i_spin_count)
      IF(i_r_point*k.eq.-1) THEN
      PRINT*,beta,i_spin_count
      PRINT*,dr1
      PRINT*,dr2
      PRINT*,di1
      PRINT*,di2	  
      ENDIF	  
      integrant	 =integrant +  (dr1*dr2+di1*di2)
     & *wb(i2)*dsin(beta)*V	  
      ENDDO	
      intgrlr  = intgrlr+integrant*xbl*2d0*pi
      ENDDO	  
	  DEALLOCATE(wfr1_fine)
	  DEALLOCATE(wfi1_fine)
	  DEALLOCATE(wfr2_fine)
	  DEALLOCATE(wfi2_fine)		  
      ENDIF	  
      ELSE
      CALL EXPANSION_MATRIX_ELEMENT(intgrlr,k,i_r_point,tmp1)	  
      ENDIF	  
	 
      CASE(2)
      IF(.NOT. expansion_defined) THEN 	  
      IF(.not.GAUSS_LEGENDRE_GRID_DEFINED)
     & CALL INI_POT_STORAGE	
c      PRINT*,"states",st_1,st_2	 
      chann_1 = indx_chann(st_1)
      chann_2 = indx_chann(st_2)
      j1_t = j_ch(chann_1)
      j2_t = j_ch(chann_2)
      m1_t = m12(st_1)
      m2_t = m12(st_2)
      v1_t = v_ch(chann_1)
      v2_t = v_ch(chann_2)
  
      DO i3=1,n_r_vib
      DO i2=1,n_beta	  
      beta = xbm + xbl*xb(i2)	  
      CALL WF_DIAT_TOP(j1_t,m1_t,beta,0d0,dr1,di1)
      CALL WF_DIAT_TOP(j2_t,m2_t,beta,0d0,dr2,di2)

!      R =  R_COM(i)*conv_unit_r
!      beta = xbm + xbl*xb(i2)	  
!      CALL V_POT_VIB_DIATOM_ATOM(V,R,r_vb_dt_inegrator(i3),beta)
!      V = V_VIB_2(i2,i3,i_r_point)     
      dr1 = dr1*wf_vib_part_diatom(i3,chann_1)
      dr2 = dr2*wf_vib_part_diatom(i3,chann_2)
!      di1 = di1*wf_vib_part_diatom(i3,chann_1)
!      di2 = di2*wf_vib_part_diatom(i3,chann_2)	  
      integrant	 = (dr1*dr2)*wg(i3)
     & *wb(i2)*dsin(beta)*V_VIB_2(i2,i3,i_r_point)
      intgrlr  = intgrlr+integrant*xgl
     & *xbl*2d0*pi
!      IF(V_VIB_2(i2,i3,i_r_point).eq.0d0)PRINT*,i3,i2,i_r_point
!      IF(wg(i3)
!     & *wb(i2).eq.0d0)PRINT*,i3,i2,i_r_point 
c      integrant = wg(i3)*
c     & wf_vib_part_diatom(i3,chann_1)*wf_vib_part_diatom(i3,chann_2)
c      intgrlr  = intgrlr+integrant*xgl
	  
      ENDDO	  
      ENDDO
      ELSE
      CALL EXPANSION_MATRIX_ELEMENT(intgrlr,k,i_r_point,tmp1)	  
      ENDIF	       	  
      CASE(3)
      IF(.NOT. expansion_defined) THEN 		  
      IF(.not.GAUSS_LEGENDRE_GRID_DEFINED)
     & CALL INI_POT_STORAGE	  
      chann_1 = indx_chann(st_1)
      chann_2 = indx_chann(st_2)
      j1_t = j_ch(chann_1)
      j2_t = j_ch(chann_2)
      k1_t = k_ch(chann_1)
      k2_t = k_ch(chann_2)
      eps1_t = eps_ch(chann_1)
      eps2_t = eps_ch(chann_2)
      m1_t = m12(st_1)
      m2_t = m12(st_2)
	  
      DO i3=1,n_gamma
      DO i2=1,n_beta
      beta = xbm + xbl*xb(i2)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i3)
	  else
	  gamma = xg(i3)
	  end if
! Bikram End.
	  
      CALL WF_SYM_TOP(j1_t,k1_t,eps1_t,m1_t,dr1,di1,0,beta,gamma)
      CALL WF_SYM_TOP(j2_t,k2_t,eps2_t,m2_t,dr2,di2,0,beta,gamma)
      V = V_3_1(i3,i2,i_r_point)	  
      integrant	 = (dr1*dr2+di1*di2)*wb(i2)*wg(i3)*dsin(beta)
      intgrlr  = intgrlr+integrant*xbl*xgl*2d0*pi*V
c      integrant	=  (- dr1*di2+di1*dr2)*wb(i2)*wg(i3)*dsin(beta)
c      intgrli  = intgrli+integrant*xbl*xgl*2d0*pi
	  
      ENDDO
      ENDDO	  
	  
! This part is to check the ortho-normality of the wavefunction 
! for the system of symmetric to rotor + atom	  
	  if(myid.eq.-10) then	  
	  print*, 'Starting Ortho-Normality Check 1'
	  do i4 = 1, size(j_ch)
	  do i5 = 1, size(j_ch)
!	  CALL INI_POT_STORAGE
!	  chann_1 = indx_chann(i4)
!      chann_2 = indx_chann(i5)
      j1_t = j_ch(i4)
      j2_t = j_ch(i5)
      k1_t = k_ch(i4)
      k2_t = k_ch(i5)
      eps1_t = eps_ch(i4)
      eps2_t = eps_ch(i5)
      m1_t = 0!m12(i4)
      m2_t = 0!m12(i5)	
      intgrlr = 0.d0
      intgrli = 0.d0
      DO i3=1,n_gamma
      DO i2=1,n_beta
      beta = xbm + xbl*xb(i2)
      gamma = xgm + xgl*xg(i3)
	  CALL WF_SYM_TOP(j1_t,k1_t,eps1_t,m1_t,dr1,di1,0,beta,gamma)
      CALL WF_SYM_TOP(j2_t,k2_t,eps2_t,m2_t,dr2,di2,0,beta,gamma)
	  integrant	 = (dr1*dr2 + di1*di2) *wb(i2)*wg(i3)*dsin(beta)
	  intgrlr  = intgrlr+integrant*xbl*xgl*2d0*pi
      integrant	=  (- dr1*di2+di1*dr2)*wb(i2)*wg(i3)*dsin(beta)
      intgrli  = intgrli+integrant*xbl*xgl*2d0*pi
	  end do
	  end do
	  write(*,'(e19.12,1x)',advance='no')intgrlr
	  end do
	  write(*,*)
	  end do
!	  stop
	  end if
	  
      ELSE
	  
! This part is to check the ortho-normality of the wavefunction 
! for the system of symmetric to rotor + atom	  
	  if(myid.eq.-10) then	  
	  print*, 'Starting Ortho-Normality Check 2'
      IF(.not.GAUSS_LEGENDRE_GRID_DEFINED)
     & CALL INI_POT_STORAGE	  
      chann_1 = indx_chann(st_1)
      chann_2 = indx_chann(st_2)
      j1_t = j_ch(chann_1)
      j2_t = j_ch(chann_2)
      k1_t = k_ch(chann_1)
      k2_t = k_ch(chann_2)
      eps1_t = eps_ch(chann_1)
      eps2_t = eps_ch(chann_2)
      m1_t = m12(st_1)
      m2_t = m12(st_2)	  
	  
	  do i4 = 1, size(j_ch)
	  do i5 = 1, size(j_ch)
!	  CALL INI_POT_STORAGE
!	  chann_1 = indx_chann(i4)
!      chann_2 = indx_chann(i5)
      j1_t = j_ch(i4)
      j2_t = j_ch(i5)
      k1_t = k_ch(i4)
      k2_t = k_ch(i5)
      eps1_t = eps_ch(i4)
      eps2_t = eps_ch(i5)
      m1_t = 0!m12(i4)
      m2_t = 0!m12(i5)	
      intgrlr = 0.d0
      DO i2=1,n_beta
      DO i3=1,n_gamma
      beta = xbm + xbl*xb(i2)
      gamma = xgm + xgl*xg(i3)
	  CALL WF_SYM_TOP(j1_t,k1_t,eps1_t,m1_t,dr1,di1,0,beta,gamma)
      CALL WF_SYM_TOP(j2_t,k2_t,eps2_t,m2_t,dr2,di2,0,beta,gamma)
	  integrant	 = (dr1*dr2 + di1*di2) *wb(i2)*wg(i3)*dsin(beta)
	  intgrlr  = intgrlr+integrant*xbl*xgl*2d0*pi
	  end do
	  end do
	  write(*,'(e19.12,1x)',advance='no')intgrlr
	  end do
	  write(*,*)
	  end do
	  stop
	  end if
	  
      CALL EXPANSION_MATRIX_ELEMENT(intgrlr,k,i_r_point,tmp1)	  
      ENDIF	 	  
      CASE(4)
      IF(.NOT. expansion_defined) THEN 		  
      IF(.not.GAUSS_LEGENDRE_GRID_DEFINED)
     & CALL INI_POT_STORAGE	  
      chann_1 = indx_chann(st_1)
      chann_2 = indx_chann(st_2)
      j1_t = j_ch(chann_1)
      j2_t = j_ch(chann_2)
      ka1_t = ka_ch(chann_1)
      ka2_t = ka_ch(chann_2)
      kc1_t = kc_ch(chann_1)
      kc2_t = kc_ch(chann_2)
      m1_t = m12(st_1)
      m2_t = m12(st_2)
      DO i2=1,n_beta
      DO i3=1,n_gamma
      beta = xbm + xbl*xb(i2)
!      gamma = xgm + xgl*xg(i3)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i3)
	  else
	  gamma = xg(i3)
	  end if
! Bikram End.
	  
      CALL WF_ASYM_TOP(A,B,C,j1_t,1,chann_1,m1_t,dr1,di1,
     & 0d0,beta,gamma)
      CALL WF_ASYM_TOP(A,B,C,j2_t,1,chann_2,m2_t,dr2,di2,
     & 0d0,beta,gamma)
      V = V_3_1(i3,i2,i_r_point)	  
      integrant	 = (dr1*dr2+di1*di2)*wb(i2)*wg(i3)*dsin(beta)
      intgrlr  = intgrlr+integrant*xbl*xgl*2d0*pi*V
c      integrant	=  (- dr1*di2+di1*dr2)*wb(i2)*wg(i3)*dsin(beta)
c      intgrli  = intgrli+integrant*xbl*xgl*2d0*pi
      ENDDO
      ENDDO
      ELSE
      CALL EXPANSION_MATRIX_ELEMENT(intgrlr,k,i_r_point,tmp1)	  
      ENDIF		  
      CASE(5)

      IF(.NOT. expansion_defined) THEN 
      IF(.not.GAUSS_LEGENDRE_GRID_DEFINED) 
     & CALL INI_POT_STORAGE
c      STOP
      chann_1 = indx_chann(st_1)
      chann_2 = indx_chann(st_2)
c      PRINT*,chann_1,chann_2	  
      j1_1_t = j1_ch(chann_1)
      j2_1_t = j2_ch(chann_1)
      j1_2_t = j1_ch(chann_2)
      j2_2_t = j2_ch(chann_2)	  
c      PRINT*,j1_t,j2_t	  
      j12_1_t = j12(st_1)
      j12_2_t = j12(st_2)	  
      m1_t = m12(st_1)
      m2_t = m12(st_2)
	  
      IF(identical_particles_defined) THEN
      p_1 = parity_state(st_1)
      p_2 = parity_state(st_2)
      ENDIF
      DO i1=1,n_beta1
      DO i2=1,n_beta2
      DO i3=1,n_alpha2	  
      beta =  xbm + xbl*xb(i1)
      bbeta = xbbm + xbbl*xbb(i2)
!      gamma = xgm + xgl*xg(i3)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
!      gamma = xgm + xgl*xg(i3)
      aalpha = xaam + xaal*xaa(i3)
	  else
!	  gamma = xg(i3)
	  aalpha = xaa(i3)
	  end if
! Bikram End.

      IF(.not.identical_particles_defined) THEN	  
      CALL WF_DIAT_DIAT(j12_1_t,m1_t,j1_1_t,j2_1_t,
     & beta,bbeta,aalpha,dr1,di1)
      CALL WF_DIAT_DIAT(j12_2_t,m2_t,j1_2_t,j2_2_t,
     & beta,bbeta,aalpha,dr2,di2)
      ELSE
      CALL
     & WF_DIAT_DIAT_IDENT(j12_1_t,m1_t,j1_1_t,j2_1_t,beta,bbeta,
     & aalpha,dr1,di1,p_1)
      CALL
     & WF_DIAT_DIAT_IDENT(j12_2_t,m2_t,j1_2_t,j2_2_t,beta,bbeta,
     & aalpha,dr2,di2,p_2)
	  
      ENDIF	
      V = V_2_2(i3,i2,i1,i_r_point)
!      CALL V_CO_CO(V,R_COM(i_r_point),alpha,beta,gamma)
      integrant	 = (dr1*dr2+di1*di2)*wb(i1)*wbb(i2)*waa(i3)*dsin(beta)
     & *dsin(bbeta)*1d0*V	  
      intgrlr  = intgrlr + integrant*xbl*xbbl*xaal*2d0*pi 
      ENDDO
      ENDDO
      ENDDO	
      ELSE
      CALL EXPANSION_MATRIX_ELEMENT(intgrlr,k,i_r_point,tmp1)
!      IF(k.eq.2) PRINT*,i_r_point,intgrlr	  
      ENDIF
      CASE(6)

      IF(.NOT. expansion_defined) THEN 
      IF(.not.GAUSS_LEGENDRE_GRID_DEFINED) 
     & CALL INI_POT_STORAGE
c      STOP
      chann_1 = indx_chann(st_1)
      chann_2 = indx_chann(st_2)
c      PRINT*,chann_1,chann_2	  
      j1_1_t = j1_ch(chann_1)
      j2_1_t = j2_ch(chann_1)
      j1_2_t = j1_ch(chann_2)
      j2_2_t = j2_ch(chann_2)	  
c      PRINT*,j1_t,j2_t	  
      j12_1_t = j12(st_1)
      j12_2_t = j12(st_2)	  
      m1_t = m12(st_1)
      m2_t = m12(st_2)
	  
      IF(identical_particles_defined) THEN
      p_1 = parity_state(st_1)
      p_2 = parity_state(st_2)
      ENDIF
      DO i_vib1 = 1,n_r_vib1
      DO i_vib2 = 1,n_r_vib2	  
      DO i1=1,n_beta1
      DO i2=1,n_beta2
      DO i3=1,n_alpha2	  
      beta =  xbm + xbl*xb(i1)
      bbeta = xbbm + xbbl*xbb(i2)
!      gamma = xgm + xgl*xg(i3)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
!      gamma = xgm + xgl*xg(i3)
      aalpha = xaam + xaal*xaa(i3)
	  else
!	  gamma = xg(i3)
	  aalpha = xaa(i3)
	  end if
! Bikram End.
      R =  R_COM(i_r_point)*conv_unit_r	  
      r_vib1 = r_vb_dt_integrator1(i_vib1)
      r_vib2 = r_vb_dt_integrator1(i_vib2)
	  
      IF(.not.identical_particles_defined) THEN	  
      CALL
     & WF_VIB_DIAT_DIAT(st_1,
     & i_vib1,i_vib2,beta,bbeta,aalpha,dr1,di1)
      CALL
     & WF_VIB_DIAT_DIAT(st_2,
     & i_vib1,i_vib2,beta,bbeta,aalpha,dr2,di2)
      ELSE
      CALL
     & WF_VIB_DIAT_DIAT_IDENT(st_1,
     & i_vib1,i_vib2,beta,bbeta,aalpha,dr1,di1,p_1)
      CALL
     & WF_VIB_DIAT_DIAT_IDENT(st_2,
     & i_vib1,i_vib2,beta,bbeta,aalpha,dr2,di2,p_2)
	  
      ENDIF
      IF(make_grid_file) THEN	  
      V = V_VIB_2_2_int_buffer(i3,i2,i1,i_vib2,
     & i_vib1)
      ELSE
      CALL V_POT_VIB_DIAT_DIAT(V,R,
     & r_vib1,r_vib2,
     & beta,bbeta,aalpha)	  
      ENDIF	  
!      CALL V_CO_CO(V,R_COM(i_r_point),alpha,beta,gamma)	  
      integrant	 = (dr1*dr2+di1*di2)*wb(i1)*wbb(i2)*waa(i3)*dsin(beta)
     & *dsin(bbeta)*1d0*V*wgg(i_vib1)*wg(i_vib2)	  
      intgrlr  = intgrlr + integrant*xbl*xbbl*xaal*2d0*pi*
     & xggl*xgl	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO	  
      ELSE
      CALL EXPANSION_MATRIX_ELEMENT(intgrlr,k,i_r_point,tmp1)
!      IF(k.eq.2) PRINT*,i_r_point,intgrlr	  
      ENDIF
      CASE(7)
      IF(.NOT. expansion_defined) THEN
      IF(.not.GAUSS_LEGENDRE_GRID_DEFINED) THEN
      IF(MYID.EQ.0)PRINT*,"GAUSS-LEGENDRE GRID INTIALIZATION STARTED"
      CALL INI_POT_STORAGE
      ENDIF
!      DO i1=1,n_alpha1	 
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      beta = xbm + xbl*xb(i3)
      bbeta = xbbm + xbbl*xbb(i4)
!      gamma = xgm + xgl*xg(i5)
!      aalpha =  xaam + xaal*xaa(i2)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
	  else
	  gamma = xg(i5)
	  aalpha = xaa(i2)
	  end if
! Bikram End.
      R =  R_COM(i_r_point)*conv_unit_r
	 
      CALL WF_SYM_TOP_DIAT_TOP(dr1,di1,st_1,0,beta,gamma,
     & aalpha,bbeta,0d0)
      CALL WF_SYM_TOP_DIAT_TOP(dr2,di2,st_2,0,beta,gamma,
     & aalpha,bbeta,0d0)
	 
      IF(make_grid_file) THEN
!      IF(i_r_point*i2*i3*i4*i5*i6
!     & .le.-1) 
!     & PRINT*,MYID,k,V_3_2_int_buffer(i5,i4,i3,i2)	  
      V = V_3_2_int_buffer(i5,i4,i3,i2)	  
      ELSE	  
      CALL V_POT_SYM_DIATOM(V,R,0d0,beta,gamma,
     & aalpha,bbeta)
!      CALL V_H2O_H2O(V,R,beta,gamma,
!     & aalpha,bbeta,ggamma)
      ENDIF	 

!      CALL V_POT_ASYM_ASYM(V,R,aalpha,bbeta,ggamma,
!     & alpha,beta,gamma)
      integrant	 = (dr1*dr2+di1*di2)*wb(i3)*wg(i5)*dsin(beta)*
     & waa(i2)*wbb(i4)*dsin(bbeta)*V!*wa(i1)
      intgrlr  = intgrlr + integrant*xbl*xgl!*xal
     & *xaal*xbbl*2*pi!**2
	 
!      integrant=(- dr1*di2+di1*dr2)*wb(i2)*wg(i3)*dsin(beta)*
!     & waa(i4)*wbb(i5)*wgg(i6)*dsin(bbeta)*V		  
!      intgrli  = intgrli+integrant*xbl*xgl
!     & *xaal*xbbl*xggl*2*pi!**2	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ELSE
      CALL EXPANSION_MATRIX_ELEMENT(intgrlr,k,i_r_point,tmp1)	  
      ENDIF
	  
      CASE(8)
      IF(.NOT. expansion_defined) THEN
      IF(.not.GAUSS_LEGENDRE_GRID_DEFINED) THEN
      IF(MYID.EQ.0)PRINT*,"GAUSS-LEGENDRE GRID INTIALIZATION STARTED"
      CALL INI_POT_STORAGE
      ENDIF
!      DO i1=1,n_alpha1	 
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      beta = xbm + xbl*xb(i3)
      bbeta = xbbm + xbbl*xbb(i4)
!      gamma = xgm + xgl*xg(i5)
!      aalpha =  xaam + xaal*xaa(i2)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
	  else
	  gamma = xg(i5)
	  aalpha = xaa(i2)
	  end if
! Bikram End.
      R =  R_COM(i_r_point)*conv_unit_r

      CALL WF_ASYM_TOP_DIAT_TOP(dr1,di1,st_1,0,beta,gamma,
     & aalpha,bbeta,0d0)
      CALL WF_ASYM_TOP_DIAT_TOP(dr2,di2,st_2,0,beta,gamma,
     & aalpha,bbeta,0d0)
      IF(make_grid_file) THEN
  
      V = V_3_2_int_buffer(i5,i4,i3,i2)	  
      ELSE	  
      CALL V_POT_ASYM_DIATOM(V,R,0d0,beta,gamma,
     & aalpha,bbeta)
      ENDIF	 


      integrant	 = (dr1*dr2+di1*di2)*wb(i3)*wg(i5)*dsin(beta)*
     & waa(i2)*wbb(i4)*dsin(bbeta)*V
      intgrlr  = intgrlr + integrant*xbl*xgl
     & *xaal*xbbl*2*pi
	 
      ENDDO
      ENDDO
      ENDDO
      ENDDO
	  
      ELSE
      CALL EXPANSION_MATRIX_ELEMENT(intgrlr,k,i_r_point,tmp1)	  
      ENDIF
	  
      CASE(9)
      IF(.NOT. expansion_defined) THEN
      IF(.not.GAUSS_LEGENDRE_GRID_DEFINED) THEN
      IF(MYID.EQ.0)PRINT*,"GAUSS-LEGENDRE GRID INTIALIZATION STARTED"
      CALL INI_POT_STORAGE
      ENDIF
!      DO i1=1,n_alpha1	 
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      DO i6=1,n_gamma2
      beta = xbm + xbl*xb(i3)
      bbeta = xbbm + xbbl*xbb(i4)
!      gamma = xgm + xgl*xg(i5)
!      aalpha =  xaam + xaal*xaa(i2)
!      ggamma = xggm + xggl*xgg(i6)	  
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
      ggamma = xggm + xggl*xgg(i6)	  
	  else
	  gamma = xg(i5)
	  aalpha = xaa(i2)
	  ggamma = xgg(i6)
	  end if
! Bikram End.
      R =  R_COM(i_r_point)*conv_unit_r

      CALL WF_ASYM_TOP_SYM_TOP(dr1,di1,st_1,0,beta,gamma,
     & aalpha,bbeta,ggamma)
      CALL WF_ASYM_TOP_SYM_TOP(dr2,di2,st_2,0,beta,gamma,
     & aalpha,bbeta,ggamma)	 
      IF(make_grid_file) THEN
      IF(i_r_point*i2*i3*i4*i5*i6
     & .le.-1) 
     & PRINT*,MYID,k,V_3_3_int_buffer(i6,i5,i4,i3,i2)	  
      V = V_3_3_int_buffer(i6,i5,i4,i3,i2)	  
      ELSE	  
      CALL V_POT_ASYM_SYM(V,R,0d0,beta,gamma,
     & aalpha,bbeta,ggamma)
!      CALL V_H2O_H2O(V,R,beta,gamma,
!     & aalpha,bbeta,ggamma)
      ENDIF	 

!      CALL V_POT_ASYM_ASYM(V,R,aalpha,bbeta,ggamma,
!     & alpha,beta,gamma)
      integrant	 = (dr1*dr2+di1*di2)*wb(i3)*wg(i5)*dsin(beta)*
     & waa(i2)*wbb(i4)*wgg(i6)*dsin(bbeta)*V!*wa(i1)
      intgrlr  = intgrlr + integrant*xbl*xgl!*xal
     & *xaal*xbbl*xggl*2*pi!**2
	 
c      integrant=(- dr1*di2+di1*dr2)*wb(i2)*wg(i3)*dsin(beta)*
c     & waa(i4)*wbb(i5)*wgg(i6)*dsin(bbeta)*V		  
c      intgrli  = intgrli+integrant*xbl*xgl
c     & *xaal*xbbl*xggl*2*pi!**2	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
	  
      ELSE
      CALL EXPANSION_MATRIX_ELEMENT(intgrlr,k,i_r_point,tmp1)	  
      ENDIF	  	  
      CASE(0)
      IF(.NOT. expansion_defined) THEN
      IF(.not.GAUSS_LEGENDRE_GRID_DEFINED) THEN
      IF(MYID.EQ.0)PRINT*,"GAUSS-LEGENDRE GRID INTIALIZATION STARTED"
      CALL INI_POT_STORAGE
      ENDIF
!      DO i1=1,n_alpha1	 
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      DO i6=1,n_gamma2
      beta = xbm + xbl*xb(i3)
      bbeta = xbbm + xbbl*xbb(i4)
!      gamma = xgm + xgl*xg(i5)
!      aalpha =  xaam + xaal*xaa(i2)
!      ggamma = xggm + xggl*xgg(i6)	  
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
      ggamma = xggm + xggl*xgg(i6)	  
	  else
	  gamma = xg(i5)
	  aalpha = xaa(i2)
	  ggamma = xgg(i6)
	  end if
! Bikram End.
      R =  R_COM(i_r_point)*conv_unit_r

      IF(.not.identical_particles_defined) THEN
      CALL WF_ASYM_TOP_ASYM_TOP(dr1,di1,st_1,0,beta,gamma,
     & aalpha,bbeta,ggamma)
      CALL WF_ASYM_TOP_ASYM_TOP(dr2,di2,st_2,0,beta,gamma,
     & aalpha,bbeta,ggamma)
      ELSE
      CALL WF_ASYM_TOP_ASYM_TOP_IDENT(dr1,di1,st_1,0,beta,gamma,
     & aalpha,bbeta,ggamma)
      CALL WF_ASYM_TOP_ASYM_TOP_IDENT(dr2,di2,st_2,0,beta,gamma,
     & aalpha,bbeta,ggamma)	  
      ENDIF	 
	 
      IF(make_grid_file) THEN
      IF(i_r_point*i2*i3*i4*i5*i6
     & .le.-1) 
     & PRINT*,MYID,k,V_3_3_int_buffer(i6,i5,i4,i3,i2)	  
      V = V_3_3_int_buffer(i6,i5,i4,i3,i2)	  
      ELSE	  
      CALL V_POT_ASYM_ASYM(V,R,0d0,beta,gamma,
     & aalpha,bbeta,ggamma)
!      CALL V_H2O_H2O(V,R,beta,gamma,
!     & aalpha,bbeta,ggamma)
      ENDIF	 

!      CALL V_POT_ASYM_ASYM(V,R,aalpha,bbeta,ggamma,
!     & alpha,beta,gamma)
      integrant	 = (dr1*dr2+di1*di2)*wb(i3)*wg(i5)*dsin(beta)*
     & waa(i2)*wbb(i4)*wgg(i6)*dsin(bbeta)*V!*wa(i1)
      intgrlr  = intgrlr + integrant*xbl*xgl!*xal
     & *xaal*xbbl*xggl*2*pi!**2
	 
c      integrant=(- dr1*di2+di1*dr2)*wb(i2)*wg(i3)*dsin(beta)*
c     & waa(i4)*wbb(i5)*wgg(i6)*dsin(bbeta)*V		  
c      intgrli  = intgrli+integrant*xbl*xgl
c     & *xaal*xbbl*xggl*2*pi!**2	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
	  
      ELSE
      CALL EXPANSION_MATRIX_ELEMENT(intgrlr,k,i_r_point,tmp1)	 
      ENDIF	  
      END SELECT	  
      END SUBROUTINE INTEGRATOR_TEMP	  
  
      SUBROUTINE EXPANSION_TERM
     & (l,l1,l2,nju1,nju2,wfr,wfi,
     & beta1,gamma1,alpha2,beta2,gamma2,coll_type,check_fact)
      IMPLICIT NONE
      INTEGER l,l1,l2,nju1,nju2,m,coll_type
      INTEGER m1,m2
      REAL*8 wfr,wfi,alpha2,beta2,gamma2,dwr1,dwi1,dwr2,dwi2,dwr3,dwi3
      REAL*8 CG,beta1,gamma1,clb_grd,delta,plgndr,W3JS,CG_bikram
	  integer kroneker
	  REAL*8 signfactor
      REAL*8, PARAMETER :: pi = dacos(-1d0)
	  logical check_fact
      EXTERNAL CG,delta,plgndr,W3JS, kroneker, CG_bikram
      wfr = 0d0
      wfi = 0d0
!c      PRINT*,coll_type	  
      SELECT CASE(coll_type)
      CASE(1)
	  CALL WF_DIAT_TOP(l,0,beta1,0d0,wfr,wfi)
      CASE(2)
      wfr =plgndr(l,0,dcos(beta1))*(2d0*l+1)/4/pi	  
      CASE(3)
      CALL WF_DIAT_TOP(l,nju1,beta1,gamma1,dwr1,dwi1)					
      wfr = dwr1*(-1)**nju1
      wfi = dwi1*(-1)**nju1
      CASE(4)
      CALL WF_DIAT_TOP(l,nju1,beta1,gamma1,dwr1,dwi1)
      wfr = dwr1*(-1)**nju1
      wfi = dwi1*(-1)**nju1
      CASE(5)
      IF(max(l,l1,l2)*2 .gt.l1+l2+l ) THEN
      PRINT*,l,l1,l2
      STOP"ERROR:TRINAGULAR RULE"
      ENDIF		  
!      CALL WF_DIAT_DIAT(l,0,l1,l2,beta1,beta2,gamma1,wfr,wfi)

      DO m = -l1, l1
      IF(ABS(m).gt.l2) cycle	  
      CALL WF_DIAT_TOP(l1,m,beta1,0.d0,dwr1,dwi1)
      CALL WF_DIAT_TOP(l2,-m,beta2,alpha2,dwr2,dwi2)	  
	 
      clb_grd	= CG_bikram(l1,l2,l,m,-m,0,check_fact)	 
      wfr = wfr +  clb_grd*(dwr1*dwr2 - dwi1*dwi2)
      wfi = wfi +  clb_grd*(dwr1*dwi2 + dwr2*dwi1)
      ENDDO
	  wfr = wfr!*dsqrt(4*pi/(2d0*l+1d0))
	  wfi = wfi!*dsqrt(4*pi/(2d0*l+1d0))

      CASE(6)
      IF(max(l,l1,l2)*2 .gt.l1+l2+l ) THEN
      PRINT*,l,l1,l2
      STOP"ERROR:TRINAGULAR RULE"
      ENDIF		  
      CALL WF_DIAT_DIAT(l,0,l1,l2,beta1,beta2,gamma1,wfr,wfi)	  
      CASE(7)
      IF(max(l,l1,l2)*2 .gt.l1+l2+l ) THEN
      PRINT*,l,l1,l2
      STOP"ERROR:TRINAGULAR RULE"
      ENDIF	  


      DO m=-l2,l2
      IF(abs(m).gt.l1) CYCLE	  
      CALL Djkm(dwr1,dwi1,l1,nju1,m,0d0,beta1,gamma1)
      CALL WF_DIAT_TOP(l2,-m,beta2,alpha2,dwr2,dwi2)
      clb_grd=CG_bikram(l1,l2,l,m,-m,0,check_fact)  
      wfr =wfr +  clb_grd*(dwr1*dwr2 - dwi1*dwi2)
      wfi =wfi +  clb_grd*(dwr1*dwi2 + dwr2*dwi1)
      ENDDO
      wfr = wfr*dsqrt((2d0*dble(l1)+1d0)/4d0/pi)
      wfi = wfi*dsqrt((2d0*dble(l1)+1d0)/4d0/pi)	
      CASE(8)
      IF(max(l,l1,l2)*2 .gt.l1+l2+l ) THEN
      PRINT*,l,l1,l2
      STOP"ERROR:TRINAGULAR RULE"
      ENDIF	  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Bikram making changes to the projection based on the manual of MQCT2022 code
	  DO m = -l2, l2
      IF(abs(m).gt.l1) CYCLE	  
      CALL Djkm(dwr1,dwi1,l1,nju1,m,0d0,beta1,gamma1)
      CALL WF_DIAT_TOP(l2,-m,beta2,alpha2,dwr2,dwi2)
      clb_grd = CG_bikram(l1,l2,l,m,-m,0,check_fact)  
      wfr = wfr +  clb_grd*(dwr1*dwr2 - dwi1*dwi2)
      wfi = wfi +  clb_grd*(dwr1*dwi2 + dwr2*dwi1)
      ENDDO
      wfr = wfr*dsqrt((2d0*dble(l1) + 1d0)/4d0/pi)
      wfi = wfi*dsqrt((2d0*dble(l1) + 1d0)/4d0/pi)
	  return
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! The below piece is old code from Alex.
      DO m=-l2,l2
      IF(abs(m).gt.l1) CYCLE	  
      CALL Djkm(dwr1,dwi1,l1,nju1,m,0d0,beta1,gamma1)
      CALL WF_DIAT_TOP(l2,-m,beta2,alpha2,dwr2,dwi2)
      clb_grd=CG_bikram(l1,l2,l,m,-m,0,check_fact)  
      wfr =wfr +  clb_grd*(dwr1*dwr2 - 
     &	dwi1*dwi2)
      wfi =wfi +  clb_grd*(dwr1*dwi2 + 
     &	dwr2*dwi1)
      ENDDO
	  signfactor = (-1)**(l1+l2+nju1+l)
      DO m=-l2,l2
      IF(abs(m).gt.l1) CYCLE	  
      CALL Djkm(dwr1,dwi1,l1,-nju1,m,0d0,beta1,gamma1)
      CALL WF_DIAT_TOP(l2,-m,beta2,alpha2,dwr2,dwi2)
      clb_grd=CG_bikram(l1,l2,l,m,-m,0,check_fact)
      wfr =wfr +  clb_grd*(dwr1*dwr2 - 
     &	dwi1*dwi2)*signfactor
      wfi =wfi +  clb_grd*(dwr1*dwi2 + 
     &	dwr2*dwi1)*signfactor
      ENDDO	  
      wfr = wfr*dsqrt((2d0*dble(l1)+1d0)/4d0/pi)/sqrt(2d0)
     & *(-1)**(l1+l2)
      wfi = wfi*dsqrt((2d0*dble(l1)+1d0)/4d0/pi)/sqrt(2d0)
     & *(-1)**(l1+l2)	  
      IF(nju1.eq.0) THEN
	  wfr=wfr/sqrt(2d0)
	  wfi=wfi/sqrt(2d0)
      ENDIF	  
      CASE(9)
      IF(max(l,l1,l2)*2 .gt.l1+l2+l ) THEN
      PRINT*,l,l1,l2
      STOP"ERROR:TRINAGULAR RULE"
      ENDIF	  
      DO m=-min(l1,l2),min(l1,l2)
      CALL Djkm(dwr1,dwi1,l1,nju1,m,0d0,beta1,gamma1)
      dwr1 = dwr1*dsqrt((2d0*dble(l1)+1d0)/8d0/pi**2)
      dwi1 = dwi1*dsqrt((2d0*dble(l1)+1d0)/8d0/pi**2)
      CALL Djkm(dwr2,dwi2,l2,nju2,-m,alpha2,beta2,gamma2)
      dwr2 = dwr2*dsqrt((2d0*dble(l2)+1d0)/8d0/pi**2)
      dwi2 = dwi2*dsqrt((2d0*dble(l2)+1d0)/8d0/pi**2)
      clb_grd	= CG_bikram(l1,l2,l,m,-m,0,check_fact)
      wfr =wfr +  clb_grd*(dwr1*dwr2 - 
     &	dwi1*dwi2)
      wfi =wfi +  clb_grd*(dwr1*dwi2 + 
     &	dwr2*dwi1) 	  
      ENDDO	  
      CASE(0)      	  
      IF(max(l,l1,l2)*2 .gt.l1+l2+l ) THEN
      PRINT*,l,l1,l2
      STOP"ERROR:TRINAGULAR RULE"
      ENDIF	  
      DO m=-min(l1,l2),min(l1,l2)
      CALL Djkm(dwr1,dwi1,l1,nju1,m,0d0,beta1,gamma1)
      dwr1 = dwr1*dsqrt((2d0*dble(l1)+1d0)/8d0/pi**2)
      dwi1 = dwi1*dsqrt((2d0*dble(l1)+1d0)/8d0/pi**2)
      CALL Djkm(dwr2,dwi2,l2,nju2,-m,alpha2,beta2,gamma2)
      dwr2 = dwr2*dsqrt((2d0*dble(l2)+1d0)/8d0/pi**2)
      dwi2 = dwi2*dsqrt((2d0*dble(l2)+1d0)/8d0/pi**2)
      clb_grd	= CG_bikram(l1,l2,l,m,-m,0,check_fact)
      wfr =wfr +  clb_grd*(dwr1*dwr2 - 
     &	dwi1*dwi2)
      wfi =wfi +  clb_grd*(dwr1*dwi2 + 
     &	dwr2*dwi1) 	  
      ENDDO
      END SELECT	  
      END SUBROUTINE EXPANSION_TERM	  
	  
      SUBROUTINE CALC_EXPANSION(coll_type,n_r_coll,L_NJU,
     & bk_idnt_trms,bk_idnt_pes, du1, du2, dlamda1, dlamda2)
      USE GRID_INTEGR
      USE MPI_DATA	  
      IMPLICIT NONE
      INTEGER coll_type,n_r_coll,L_NJU(8)
      INTEGER L1,L2,L, NJU1,NJU2
      INTEGER L1_MAX,L2_MAX,L_MAX,NJU_MAX,NJU1_MAX,NJU2_MAX	
      INTEGER L1_MIN,L2_MIN
      INTEGER L2_RANGE_MAX	  
      INTEGER N_TERMS,i_vib_terms, n_tmp, i, n1r, n2r
	  integer, allocatable :: l1s(:), n1s(:), l2s(:), n2s(:), ls(:)
	  integer du1, du2, dlamda1, dlamda2
	  logical same_term, exst, bk_idnt_trms, bk_idnt_pes

      inquire(file = "tau_functions.tmp", exist = exst)
	  if(exst) call system ("rm tau_functions.tmp")
	  
	  N_TERMS = 0
      L_MAX = L_NJU(1)
      L1_MAX=L_NJU(2)
      L2_MAX = L_NJU(3)
      NJU_MAX = L_NJU(4)
      NJU1_MAX = L_NJU(5)
      NJU2_MAX = L_NJU(6)
      L1_MIN =  L_NJU(7)
      L2_MIN = L_NJU(8)	  
      SELECT CASE(coll_type)
      CASE(1)
      IF(L_MAX.lt.0) THEN
      CALL READ_USER_TERMS
      CALL CALC_EXPANSION_TERMS	  
      RETURN	  
      end if
	  
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! starting generation of terms
	  open(3, file = "EXPANSION_TERMS_LIST.DAT")
	  
	  if(.not.bk_idnt_pes) then
	  DO L1 = 0,L_MAX, dlamda1
      N_TERMS = N_TERMS + 1
	  write(3,'(1(i0,a))') L1,','
      ENDDO
	  
	  else
	  print*, "Error in terms generation."
	  stop
	  end if
	  
	  close(3)
	  Write(*,'(a,a)') "The list of expansion terms was generated ", 
     & "and saved in the file named EXPANSION_TERMS_LIST.DAT"
	  print*, "Program will stop now."
	  stop
! Bikram End.
	  
      DO L=0,L_MAX
      N_TERMS = N_TERMS + 1	  
      ENDDO
      N_TERMS_EXP = N_TERMS 	  
      ALLOCATE(A_TOP(5,N_TERMS))
      ALLOCATE(expansion_terms(n_r_coll,N_TERMS))
      ALLOCATE(expansion_terms_im(n_r_coll,N_TERMS))
      ALLOCATE(weig_terms(N_TERMS))
      expansion_terms = 0
      A_TOP = 0
      weig_terms = 0	  
      N_TERMS = 0	  

      DO L=0,L_MAX
      N_TERMS = N_TERMS + 1
      A_TOP(1,N_TERMS) = L

      ENDDO
  

      CASE(2)
      IF(L_MAX.lt.0) THEN
      CALL READ_USER_TERMS
      CALL CALC_EXPANSION_TERMS	  
      RETURN	  
      ELSE
      DO L=0,L_MAX
      N_TERMS = N_TERMS + 1	  
      ENDDO
      N_TERMS_EXP = N_TERMS 	  
      ALLOCATE(A_TOP(5,N_TERMS))
      ALLOCATE(expansion_terms(n_r_coll,N_TERMS))
      ALLOCATE(expansion_terms_im(n_r_coll,N_TERMS))
      ALLOCATE(weig_terms(N_TERMS))
      expansion_terms = 0
      A_TOP = 0
      weig_terms = 0	  
      N_TERMS = 0	  

      DO L=0,L_MAX
      N_TERMS = N_TERMS + 1
      A_TOP(1,N_TERMS) = L

      ENDDO
  
      ENDIF
      CASE(3)
      IF(L_MAX.lt.0) THEN
      CALL READ_USER_TERMS
      CALL CALC_EXPANSION_TERMS	  
      RETURN
	  end if
	  
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! starting generation of terms
	  open(3, file = "EXPANSION_TERMS_LIST.DAT")
	  if(NJU_MAX.lt.0) NJU_MAX = L_MAX
	  
	  if(.not.bk_idnt_pes) then
	  DO L1 = 0,L_MAX, dlamda1
      DO NJU1 = 0, min(NJU_MAX,L1), du1
      N_TERMS = N_TERMS + 1
	  write(3,'(2(i0,a))') L1,',',NJU1
      ENDDO
      ENDDO
	  
	  else
	  print*, "Error in terms generation."
	  stop
	  end if
	  
	  close(3)
	  Write(*,'(a,a)') "The list of expansion terms was generated ", 
     & "and saved in the file named EXPANSION_TERMS_LIST.DAT"
	  print*, "Program will stop now."
	  stop
! Bikram End.

      IF(NJU_MAX.lt.0) NJU_MAX  = 	L_MAX  
      DO L=0,L_MAX
      DO NJU1 = 0,MIN(NJU_MAX,L)	  
      N_TERMS = N_TERMS + 1	  
      ENDDO
      ENDDO	  
      N_TERMS_EXP = N_TERMS 	  
      ALLOCATE(A_TOP(5,N_TERMS))
      ALLOCATE(expansion_terms(n_r_coll,N_TERMS))
      ALLOCATE(expansion_terms_im(n_r_coll,N_TERMS))
      ALLOCATE(weig_terms(N_TERMS))
      expansion_terms = 0
      A_TOP = 0
      weig_terms = 0	  
      N_TERMS = 0	  

      DO L=0,L_MAX
      DO NJU1 = 0,MIN(NJU_MAX,L)
      N_TERMS = N_TERMS + 1
      A_TOP(1,N_TERMS) = L
      A_TOP(2,N_TERMS) = NJU1
      ENDDO
      ENDDO	  

      CASE(4)
      IF(L_MAX.lt.0) THEN
      CALL READ_USER_TERMS
      CALL CALC_EXPANSION_TERMS	  
      RETURN
	  end if
	  	  
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! starting generation of terms
	  open(3, file = "EXPANSION_TERMS_LIST.DAT")
	  if(NJU_MAX.lt.0) NJU_MAX = L_MAX
	  
	  if(.not.bk_idnt_pes) then
	  DO L1 = 0,L_MAX, dlamda1
      DO NJU1 = 0, min(NJU_MAX,L1), du1
      N_TERMS = N_TERMS + 1
	  write(3,'(2(i0,a))') L1,',',NJU1
      ENDDO
      ENDDO
	  
	  else
	  print*, "Error in terms generation."
	  stop
	  end if
	  
	  close(3)
	  Write(*,'(a,a)') "The list of expansion terms was generated ", 
     & "and saved in the file named EXPANSION_TERMS_LIST.DAT"
	  print*, "Program will stop now."
	  stop
! Bikram End.
	  
      IF(NJU_MAX.lt.0) NJU_MAX  = 	L_MAX  
      DO L=0,L_MAX
      DO NJU1 = 0,MIN(NJU_MAX,L)	  
      N_TERMS = N_TERMS + 1	  
      ENDDO
      ENDDO	  
      N_TERMS_EXP = N_TERMS 	  
      ALLOCATE(A_TOP(5,N_TERMS))
      ALLOCATE(expansion_terms(n_r_coll,N_TERMS))
      ALLOCATE(expansion_terms_im(n_r_coll,N_TERMS))
      ALLOCATE(weig_terms(N_TERMS))
      expansion_terms = 0
      A_TOP = 0
      weig_terms = 0	  
      N_TERMS = 0	  

      DO L=0,L_MAX
      DO NJU1 = 0,MIN(NJU_MAX,L)
      N_TERMS = N_TERMS + 1
      A_TOP(1,N_TERMS) = L
      A_TOP(2,N_TERMS) = NJU1
      ENDDO
      ENDDO	  
  
      CASE(5)
      IF(L1_MAX.lt.0 .or. L2_MAX.lt.0) THEN
      CALL READ_USER_TERMS
      CALL CALC_EXPANSION_TERMS	  
      RETURN
	  end if
	  
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! starting generation of terms
	  open(3, file = "EXPANSION_TERMS_LIST.DAT")
	  if(L_MAX.lt.0) L_MAX = L1_MAX+L2_MAX
	  
	  if(.not.bk_idnt_pes) then
	  DO L1 = L1_MIN,L1_MAX, dlamda1
      DO L2 = L2_MIN,L2_MAX, dlamda2
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      IF(.not.bk_idnt_trms .and. mod((l1+l2+l),2).ne.0) CYCLE 	  
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"
      N_TERMS = N_TERMS + 1
	  write(3,'(3(i0,a))') L1, ',', L2, ',', L, ','
      ENDDO
      ENDDO
      ENDDO
	  
	  else if(bk_idnt_pes) then
	  DO L1 = L1_MIN,L1_MAX, dlamda1
      DO L2 = L2_MIN,min(L2_MAX,L1), dlamda2
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      IF(.not.bk_idnt_trms .and. mod((l1+l2+l),2).ne.0) CYCLE 	  
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"
      N_TERMS = N_TERMS + 1
	  write(3,'(3(i0,a))') L1, ',', L2, ',', L, ','
      ENDDO
      ENDDO
      ENDDO	 
	  
	  else
	  print*, "Error in terms generation."
	  stop
	  end if
	  
	  close(3)
	  Write(*,'(a,a)') "The list of expansion terms was generated ", 
     & "and saved in the file named EXPANSION_TERMS_LIST.DAT"
	  print*, "Program will stop now."
	  stop
! Bikram End.
	  
      DO L1 = 0,L1_MAX
      DO L2 = 0,L2_MAX
      DO L=abs(L1-L2),min(L1+L2,L)
      IF(((L1+L2+L)/2)*2 .eq.L1+L2+L) THEN
      N_TERMS = N_TERMS + 1
      ENDIF	  
      ENDDO
      ENDDO
      ENDDO
      N_TERMS_EXP = N_TERMS 	  
      ALLOCATE(A_TOP(5,N_TERMS))
      ALLOCATE(expansion_terms(n_r_coll,N_TERMS))
      ALLOCATE(expansion_terms_im(n_r_coll,N_TERMS))
      ALLOCATE(weig_terms(N_TERMS))
      expansion_terms = 0
      A_TOP = 0
      weig_terms = 0	  
      N_TERMS = 0	  
      DO L1 = 0,L1_MAX
      DO L2 = 0,L2_MAX
      DO L=abs(L1-L2),L1+L2
      IF(((L1+L2+L)/2)*2 .eq.L1+L2+L) THEN	  
      N_TERMS = N_TERMS + 1
      A_TOP(1,N_TERMS) = L1
      A_TOP(2,N_TERMS) = L2
      A_TOP(3,N_TERMS) = L
      ENDIF	  
      ENDDO
      ENDDO
      ENDDO
      CASE(6)
      IF(L1_MAX.lt.0 .or. L2_MAX.lt.0) THEN
      CALL READ_USER_TERMS
      CALL CALC_EXPANSION_TERMS	  
      RETURN	  
      ELSE	  
      DO L1 = 0,L1_MAX
      DO L2 = 0,L2_MAX
      DO L=abs(L1-L2),min(L1+L2,L)
      IF(((L1+L2+L)/2)*2 .eq.L1+L2+L) THEN
      N_TERMS = N_TERMS + 1
      ENDIF	  
      ENDDO
      ENDDO
      ENDDO
      N_TERMS_EXP = N_TERMS 	  
      ALLOCATE(A_TOP(5,N_TERMS))
      ALLOCATE(expansion_terms(n_r_coll,N_TERMS))
      ALLOCATE(expansion_terms_im(n_r_coll,N_TERMS))
      ALLOCATE(weig_terms(N_TERMS))
      expansion_terms = 0
      A_TOP = 0
      weig_terms = 0	  
      N_TERMS = 0	  
      DO L1 = 0,L1_MAX
      DO L2 = 0,L2_MAX
      DO L=abs(L1-L2),L1+L2
      IF(((L1+L2+L)/2)*2 .eq.L1+L2+L) THEN	  
      N_TERMS = N_TERMS + 1
      A_TOP(1,N_TERMS) = L1
      A_TOP(2,N_TERMS) = L2
      A_TOP(3,N_TERMS) = L
      ENDIF	  
      ENDDO
      ENDDO
      ENDDO
      ENDIF	  
      CASE(7)
      IF(L1_MAX.lt.0 .or. L2_MAX.lt.0) THEN
      CALL READ_USER_TERMS
      CALL CALC_EXPANSION_TERMS	  
      RETURN
      ENDIF
	  	  
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! starting generation of terms
	  open(3, file = "EXPANSION_TERMS_LIST.DAT")
	  if(NJU1_MAX.lt.0) NJU1_MAX = L1_MAX
	  if(NJU2_MAX.lt.0) NJU2_MAX = L2_MAX
	  if(L_MAX.lt.0) L_MAX = L1_MAX+L2_MAX
	  
	  if(.not.bk_idnt_pes) then
	  DO L1 = L1_MIN,L1_MAX, dlamda1
      DO NJU1 = 0, min(L1,NJU1_MAX), du1
      DO L2 = L2_MIN,L2_MAX, dlamda2
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      IF(NJU1.eq.0 .and. NJU2.lt.0) CYCLE
      IF(.not.bk_idnt_trms .and. mod((l1+l2+l),2).ne.0) CYCLE 	  
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"
      N_TERMS = N_TERMS + 1
	  write(3,'(4(i0,a))') L1,',',NJU1,',',L2,',',L,','
      ENDDO
      ENDDO
      ENDDO
      ENDDO	 
	  
	  else
	  print*, "Error in terms generation."
	  stop
	  end if
	  
	  close(3)
	  Write(*,'(a,a)') "The list of expansion terms was generated ", 
     & "and saved in the file named EXPANSION_TERMS_LIST.DAT"
	  print*, "Program will stop now."
	  stop
! Bikram End.

      DO L1 = L1_MIN,L1_MAX
      DO L2 = L2_MIN,L2_MAX! min(L2_MAX,L1)
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      DO NJU1 = 0, min(L1,NJU1_MAX)
      N_TERMS = N_TERMS + 1
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ALLOCATE(A_TOP(5,N_TERMS))
      ALLOCATE(expansion_terms(n_r_coll,N_TERMS))
      ALLOCATE(expansion_terms_im(n_r_coll,N_TERMS))	  
	  
      ALLOCATE(weig_terms(N_TERMS))	  
      IF(MYID.EQ.0) PRINT*, "NTERMS",N_TERMS  
      N_TERMS = 0
      A_TOP = 0  	  
      DO L1 = L1_MIN,L1_MAX
      DO L2 = L2_MIN,L2_MAX
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      DO NJU1 = 0, L1
      N_TERMS = N_TERMS + 1
      A_TOP(1,N_TERMS) = L1
      A_TOP(2,N_TERMS) = NJU1
      A_TOP(3,N_TERMS) = L2
      A_TOP(4,N_TERMS) = L
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      N_TERMS_EXP = N_TERMS	  
      CASE(8)
      IF(L1_MAX.lt.0 .or. L2_MAX.lt.0) THEN
      CALL READ_USER_TERMS
      CALL CALC_EXPANSION_TERMS	  
      RETURN
      ENDIF
	  	  
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! starting generation of terms
	  open(3, file = "EXPANSION_TERMS_LIST.DAT")
	  if(NJU1_MAX.lt.0) NJU1_MAX = L1_MAX
	  if(NJU2_MAX.lt.0) NJU2_MAX = L2_MAX
	  if(L_MAX.lt.0) L_MAX = L1_MAX+L2_MAX
	  
	  if(.not.bk_idnt_pes) then
	  DO L1 = L1_MIN,L1_MAX, dlamda1
      DO NJU1 = 0, min(L1,NJU1_MAX), du1
      DO L2 = L2_MIN,L2_MAX, dlamda2
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      IF(NJU1.eq.0 .and. NJU2.lt.0) CYCLE
      IF(.not.bk_idnt_trms .and. mod((l1+l2+l),2).ne.0) CYCLE 	  
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"
      N_TERMS = N_TERMS + 1
	  write(3,'(4(i0,a))') L1,',',NJU1,',',L2,',',L,','
      ENDDO
      ENDDO
      ENDDO
      ENDDO	 
	  
	  else
	  print*, "Error in terms generation."
	  stop
	  end if
	  
	  close(3)
	  Write(*,'(a,a)') "The list of expansion terms was generated ", 
     & "and saved in the file named EXPANSION_TERMS_LIST.DAT"
	  print*, "Program will stop now."
	  stop
! Bikram End.

      DO L1 = L1_MIN,L1_MAX
      DO L2 = L2_MIN,L2_MAX! min(L2_MAX,L1)
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      DO NJU1 = 0, min(L1,NJU1_MAX)
      N_TERMS = N_TERMS + 1
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ALLOCATE(A_TOP(5,N_TERMS))
      ALLOCATE(expansion_terms(n_r_coll,N_TERMS))
      ALLOCATE(expansion_terms_im(n_r_coll,N_TERMS))	  
	  
      ALLOCATE(weig_terms(N_TERMS))	  
      IF(MYID.EQ.0) PRINT*, "NTERMS",N_TERMS  
      N_TERMS = 0
      A_TOP = 0  	  
      DO L1 = L1_MIN,L1_MAX
      DO L2 = L2_MIN,L2_MAX
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      DO NJU1 = 0, L1
      N_TERMS = N_TERMS + 1
      A_TOP(1,N_TERMS) = L1
      A_TOP(2,N_TERMS) = NJU1
      A_TOP(3,N_TERMS) = L2
      A_TOP(4,N_TERMS) = L
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      N_TERMS_EXP = N_TERMS	  
      CASE(9)
      IF(L1_MAX.lt.0 .or. L2_MAX.lt.0) THEN
      CALL READ_USER_TERMS
      CALL CALC_EXPANSION_TERMS	  
      RETURN
      ENDIF
	  
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! starting generation of terms
	  open(3, file = "EXPANSION_TERMS_LIST.DAT")
	  if(NJU1_MAX.lt.0) NJU1_MAX = L1_MAX
	  if(NJU2_MAX.lt.0) NJU2_MAX = L2_MAX
	  if(L_MAX.lt.0) L_MAX = L1_MAX+L2_MAX
	  
	  if(.not.bk_idnt_pes) then
	  DO L1 = L1_MIN,L1_MAX, dlamda1
      DO NJU1 = 0, min(L1,NJU1_MAX), du1
      DO L2 = L2_MIN,L2_MAX, dlamda2! min(L2_MAX,L1)
      DO NJU2 = -min(L2,NJU2_MAX), min(L2,NJU2_MAX), du2
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      IF(NJU1.eq.0 .and. NJU2.lt.0) CYCLE
      IF(.not.bk_idnt_trms .and. mod((l1+l2+l),2).ne.0) CYCLE 	  
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"
      N_TERMS = N_TERMS + 1
	  write(3,'(5(i0,a))') L1,',',NJU1,',',L2,',',NJU2,',',L,','
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO	 
	  	  
	  else
	  print*, "Error in terms generation."
	  stop
	  end if
	  
	  close(3)
	  Write(*,'(a,a)') "The list of expansion terms was generated ", 
     & "and saved in the file named EXPANSION_TERMS_LIST.DAT"
	  print*, "Program will stop now."
	  stop
! Bikram End.
	  
      DO L1 = L1_MIN,L1_MAX
      DO L2 = L2_MIN,L2_MAX! min(L2_MAX,L1)
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      DO NJU1 = 0, min(L1,NJU1_MAX)
      DO NJU2 = -min(L2,NJU2_MAX), min(L2,NJU2_MAX)
      IF(NJU1.eq.0 .and. NJU2.lt.0) CYCLE 	  
      N_TERMS = N_TERMS + 1
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO	  
      ALLOCATE(A_TOP(5,N_TERMS))
      ALLOCATE(expansion_terms(n_r_coll,N_TERMS))
      ALLOCATE(expansion_terms_im(n_r_coll,N_TERMS))	  
	  
      ALLOCATE(weig_terms(N_TERMS))	  
      IF(MYID.EQ.0) PRINT*, "NTERMS",N_TERMS  
      N_TERMS = 0
      A_TOP = 0  	  
      DO L1 = L1_MIN,L1_MAX
      DO L2 = L2_MIN,L2_MAX
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      DO NJU1 = 0, L1
      DO NJU2 = -L2, L2
      IF(NJU1.eq.0 .and. NJU2.lt.0) CYCLE   	  
      N_TERMS = N_TERMS + 1
      A_TOP(1,N_TERMS) = L1
      A_TOP(2,N_TERMS) = NJU1
      A_TOP(3,N_TERMS) = L2
      A_TOP(4,N_TERMS) = NJU2
      A_TOP(5,N_TERMS) = L
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      N_TERMS_EXP = N_TERMS 	  
      CASE(0)
      IF(L1_MAX.lt.0 .or. L2_MAX.lt.0) THEN
      CALL READ_USER_TERMS
      CALL CALC_EXPANSION_TERMS	  
      RETURN
      ENDIF	  
	  
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! starting generation of terms
	  open(3, file = "EXPANSION_TERMS_LIST.DAT")
	  if(NJU1_MAX.lt.0) NJU1_MAX = L1_MAX
	  if(NJU2_MAX.lt.0) NJU2_MAX = L2_MAX
	  if(L_MAX.lt.0) L_MAX = L1_MAX+L2_MAX
	  
	  if(.not.bk_idnt_pes) then
	  DO L1 = L1_MIN,L1_MAX, dlamda1
      DO NJU1 = 0, min(L1,NJU1_MAX), du1
      DO L2 = L2_MIN,L2_MAX, dlamda2! min(L2_MAX,L1)
      DO NJU2 = -min(L2,NJU2_MAX), min(L2,NJU2_MAX), du2
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      IF(NJU1.eq.0 .and. NJU2.lt.0) CYCLE
      IF(.not.bk_idnt_trms .and. mod((l1+l2+l),2).ne.0) CYCLE 	  
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"
      N_TERMS = N_TERMS + 1
	  write(3,'(5(i0,a))') L1,',',NJU1,',',L2,',',NJU2,',',L,','
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO	 
	  
	  else if(bk_idnt_pes) then
	  DO L1 = L1_MIN,L1_MAX, dlamda1
      DO NJU1 = 0, min(L1,NJU1_MAX), du1
      DO L2 = L2_MIN, min(L2_MAX,L1), dlamda2
      DO NJU2 = -min(L2,NJU2_MAX), min(L2,NJU2_MAX), du2
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      IF(NJU1.eq.0 .and. NJU2.lt.0) CYCLE 	  
	  if(L1.eq.L2 .and. abs(NJU1).lt.abs(NJU2)) cycle
      IF(.not.bk_idnt_trms .and. mod((l1+l2+l),2).ne.0) CYCLE 	  
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"

	  if(N_TERMS.gt.0) then
	  if(allocated(l1s)) deallocate(l1s)
	  if(allocated(n1s)) deallocate(n1s)
	  if(allocated(l2s)) deallocate(l2s)
	  if(allocated(n2s)) deallocate(n2s)
	  if(allocated(ls)) deallocate(ls)
	  allocate(l1s(N_TERMS), n1s(N_TERMS), l2s(N_TERMS), n2s(N_TERMS),
     & ls(N_TERMS))
	  open(1,file = "tau_functions.tmp")
	  do i = 1, N_TERMS
	  read(1,*) l1s(i), n1s(i), l2s(i), n2s(i), ls(i)
	  end do
	  close(1)
	  
	  same_term = .false.
	  do i = 1, N_TERMS
	  if(l1s(i).eq.L2 .and. n1s(i).eq.NJU1 .and. l2s(i).eq.L1 .and. 
     & n2s(i).eq.NJU1 .and. ls(i).eq.L) then
	  same_term = .true.
	  exit
	  end if
	  end do
	  end if

	  if(.not.same_term) then
	  open(2, file = "tau_functions.tmp", access = "append")
	  write(2,'(5(i0,a))') L1,',',NJU1,',',L2,',',NJU2,',',L,','
	  write(3,'(5(i0,a))') L1,',',NJU1,',',L2,',',NJU2,',',L,','
	  close(2)
	  N_TERMS = N_TERMS + 1
	  end if
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO	 
	  
	  else
	  print*, "Error in terms generation."
	  stop
	  end if
	  
	  close(3)
	  Write(*,'(a,a)') "The list of expansion terms was generated ", 
     & "and saved in the file named EXPANSION_TERMS_LIST.DAT"
	  print*, "Program will stop now."
	  stop
! Bikram End.
	  
!      DO L1 = L1_MIN,L1_MAX
!      DO L2 = L2_MIN,L2_MAX! min(L2_MAX,L1)
!      DO L = abs(L1-L2),min(L1+L2,L_MAX)
!      DO NJU1 = 0, min(L1,NJU1_MAX)
!      DO NJU2 = -min(L2,NJU2_MAX), min(L2,NJU2_MAX)
!      IF(NJU1.eq.0 .and. NJU2.lt.0) CYCLE 	  
!      IF(.not.bk_idnt_trms .and. mod((l1+l2+l),2).ne.0) CYCLE 	  
!      N_TERMS = N_TERMS + 1
!       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
!     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"	  
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO	  
      ALLOCATE(A_TOP(5,N_TERMS))
      ALLOCATE(expansion_terms(n_r_coll,N_TERMS))
      ALLOCATE(expansion_terms_im(n_r_coll,N_TERMS))	  
	  
      ALLOCATE(weig_terms(N_TERMS))	  
      IF(MYID.EQ.0) PRINT*, "NTERMS",N_TERMS  
	  n_tmp = N_TERMS	  
	  if(N_TERMS.gt.0) then
	  if(allocated(l1s)) deallocate(l1s)
	  if(allocated(n1s)) deallocate(n1s)
	  if(allocated(l2s)) deallocate(l2s)
	  if(allocated(n2s)) deallocate(n2s)
	  if(allocated(ls)) deallocate(ls)
	  allocate(l1s(N_TERMS), n1s(N_TERMS), l2s(N_TERMS), n2s(N_TERMS),
     & ls(N_TERMS))
	  open(1,file = "tau_functions.tmp")
	  do i = 1, N_TERMS
	  read(1,*) l1s(i), n1s(i), l2s(i), n2s(i), ls(i)
	  end do
	  close(1)
	  end if
	  
      N_TERMS = 0
      A_TOP = 0  	  
	  if(.not.bk_idnt_pes) then
	  DO L1 = L1_MIN,L1_MAX
      DO L2 = L2_MIN,L2_MAX! min(L2_MAX,L1)
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      DO NJU1 = 0, min(L1,NJU1_MAX)
      DO NJU2 = -min(L2,NJU2_MAX), min(L2,NJU2_MAX)
      IF(NJU1.eq.0 .and. NJU2.lt.0) CYCLE 	  
      IF(.not.bk_idnt_trms .and. mod((l1+l2+l),2).ne.0) CYCLE 	  
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"
      N_TERMS = N_TERMS + 1
      A_TOP(1,N_TERMS) = L1
      A_TOP(2,N_TERMS) = NJU1
      A_TOP(3,N_TERMS) = L2
      A_TOP(4,N_TERMS) = NJU2
      A_TOP(5,N_TERMS) = L
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO	 
	  
	  else if(bk_idnt_pes) then
	  DO L1 = L1_MIN,L1_MAX
      DO L2 = L2_MIN, min(L2_MAX,L1)
      DO L = abs(L1-L2),min(L1+L2,L_MAX)
      DO NJU1 = 0, min(L1,NJU1_MAX)
      DO NJU2 = -min(L2,NJU2_MAX), min(L2,NJU2_MAX)
      IF(NJU1.eq.0 .and. NJU2.lt.0) CYCLE 	  
	  if(L1.eq.L2 .and. abs(NJU1).lt.abs(NJU2)) cycle
      IF(.not.bk_idnt_trms .and. mod((l1+l2+l),2).ne.0) CYCLE 	  
       IF(max(l,l1,l2)*2 .gt.l1+l2+l )
     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO	 
	  
	  same_term = .false.
	  do i = 1, n_tmp
	  if(l1s(i).eq.L2 .and. n1s(i).eq.NJU1 .and. l2s(i).eq.L1 .and. 
     & n2s(i).eq.NJU1 .and. ls(i).eq.L) then
	  same_term = .true.
	  exit
	  end if
	  end do

	  if(.not.same_term) then
	  N_TERMS = N_TERMS + 1
      A_TOP(1,N_TERMS) = L1
      A_TOP(2,N_TERMS) = NJU1
      A_TOP(3,N_TERMS) = L2
      A_TOP(4,N_TERMS) = NJU2
      A_TOP(5,N_TERMS) = L
	  end if
	  
	  else
	  print*, "Error in terms generation."
	  stop
	  end if

      inquire(file = "tau_functions.tmp", exist = exst)
	  if(exst) call system ("rm tau_functions.tmp")
	  
	  
!      DO L1 = L1_MIN,L1_MAX
!      DO L2 = L2_MIN,L2_MAX
!      DO L = abs(L1-L2),min(L1+L2,L_MAX)
!      DO NJU1 = 0, L1
!      DO NJU2 = -L2, L2
!      IF(NJU1.eq.0 .and. NJU2.lt.0) CYCLE   	  
!      N_TERMS = N_TERMS + 1
!      A_TOP(1,N_TERMS) = L1
!      A_TOP(2,N_TERMS) = NJU1
!      A_TOP(3,N_TERMS) = L2
!      A_TOP(4,N_TERMS) = NJU2
!      A_TOP(5,N_TERMS) = L
!      IF(max(l,l1,l2)*2 .gt.l1+l2+l )
!     &	  STOP"ERROR:TRINAGULAR RULE IN EXPANSION"	  
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO
	  
      N_TERMS_EXP = N_TERMS 
      END SELECT
      IF(MYID.EQ.0) PRINT*,"NUMBER_OF_TERMS",N_TERMS	  
      CALL CALC_EXPANSION_TERMS
      RETURN	  
      END SUBROUTINE CALC_EXPANSION

      SUBROUTINE CALC_EXPANSION_TERMS
      USE VARIABLES
      USE POT_STORE	  
      USE GRID_INTEGR
      USE CONSTANTS
      USE MPI_DATA
      USE MPI
      IMPLICIT NONE
	  integer bk_i,bk_io
	  logical file_exst
      CHARACTER(LEN=19) :: EXP_FILE_NAME = "PES_EXPAN_TERMS.DAT"	  
      INTEGER l1_t,l2_t,l_t,nju1_t,nju2_t,nju_t
      INTEGER i2,i3,i4,i5,i6,i,task_size,i_r_point,i_vib,i_vib1,i_vib2
      INTEGER buffer_size,istat
      REAL*8 R,V,intgrlr,wfr,wfi
      REAL*8 integrant,intgrl,intgrli,weig_b
      REAL*8 buffer(n_r_coll),buffer_w,buffer_i(n_r_coll)
	  real*8 bk_tym,bk_tym_1,bk_tym_2			!Bikram Jan,'19
      INTEGER coll_type_checker, 
     & nb1_checker,na1_checker,ng1_checker,
     & nb2_checker,na2_checker,ng2_checker,nr_checker
      INTEGER i_nr_ini, i_nr_fin
      REAL*8, ALLOCATABLE :: V_TRANSFER(:,:,:,:,:,:)
      REAL*8, ALLOCATABLE :: V_TRANS_BUFF(:,:,:,:,:)
      REAL*8 delta, norm_cost_bk
	  integer kroneker, nmb_terms, l1bk, n1bk, l2bk, n2bk, lbk
	  integer planr_coeff, symmetry_coeff, grd_tot, grd_cntr
	  real*8, allocatable :: V_store(:), V_recalc(:)
	  real*8, allocatable :: wfr_store(:,:),wfi_store(:,:)
	  real*8, allocatable :: wfr_tmp(:), wfi_tmp(:)
	  real*8 v_sum, rms
	  logical same_term
      EXTERNAL delta, kroneker
      CHARACTER(LEN=22) :: header = "V                     "
      CHARACTER(LEN=100) :: fname, bk_header
	  
	  
      intgrlr = 0d0
      intgrli = 0d0
      buffer = 0d0
      i_nr_ini = max(ir_bgn_exp_pnt,1)
      i_nr_fin = min(ir_fin_exp_pnt,n_r_coll)
! Bikram Oct'18 Start:
	  IF(ir_fin_exp_pnt.lt.0 .or. ir_bgn_exp_pnt.lt.0) THEN
      i_nr_ini = 1
      i_nr_fin = n_r_coll	  
      ENDIF
! Bikram End.
      IF(i_nr_ini.gt.n_r_coll .or. i_nr_fin.lt.1) THEN
      IF(myid.eq.0) THEN
      PRINT*,
     & "ERROR: DEFINE PROPERLY GRID RANGE"
      PRINT*,"i_nr_ini",i_nr_ini
      PRINT*,"n_r_coll",n_r_coll  
      ENDIF	
      STOP	  
      ENDIF	  
      conv_unit_r = 1d0
      N_TERMS_EXP = nterms	  
      task_size = n_r_coll*N_TERMS_EXP
      IF(N_TERMS_EXP.ne.nproc) THEN
      IF(MYID.eq.0) WRITE(*,'(a23,a44)')
     & "TERMINATION: N_TERMS IN",
     & " EXPANSION MUST BE EQUAL TO NUMBER OF PROCS."
      STOP	  
      ENDIF	  
      CALL INI_POT_STORAGE
      IF(make_grid_file .and. grid_file_found
     & .and. (coll_type.le.1 .or. coll_type.gt.5) ) THEN   
      IF(myid.eq.0) PRINT*,"EXPANSION TERMS WILL BE COMPUTED FROM GRID"
      IF(coll_type.eq.1) THEN
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) nb1_checker,nr_checker
      READ(1)  V_2_1
      CLOSE(1)
      ENDIF	 	  
      SELECT CASE(coll_type)
!!!! SUBSEQUENT GRID READING FOR R_GRID ONLY FOR ASYM + ASYM
      CASE(6)
      ALLOCATE(	  
     & V_TRANSFER(1,n_gamma1,n_beta2,n_beta1,
     & n_alpha2,i_nr_fin-i_nr_ini+1),
     & V_TRANS_BUFF(1,n_gamma1,n_beta2,n_beta1,
     & n_alpha2))	  
      buffer_size = 1*n_gamma1*n_beta2*
     & n_beta1*n_alpha2*(i_nr_fin-i_nr_ini+1)  

      IF(MYID.EQ.0) THEN	  
      OPEN(1,FILE=potential_file_name,FORM="UNFORMATTED",
     & ACTION="READ",STATUS="OLD")	  
      READ(1) coll_type_checker
      READ(1)	na2_checker,nb1_checker,nb2_checker,ng1_checker,
     & nr_checker	  
      IF(coll_type_checker.ne.coll_type) STOP "WRONG COLL TYPE"
      IF(na2_checker.ne.n_alpha2) STOP "WRONG GRID"
      IF(nb1_checker.ne.n_beta1) STOP "WRONG GRID"
      IF(nb2_checker.ne.n_beta2) STOP "WRONG GRID"
      IF(ng1_checker.ne.n_gamma1) STOP "WRONG GRID"
      IF(nr_checker.ne.n_r_coll) STOP "WRONG GRID"
      PRINT*,"DATA READING STARTED"	 	  
      DO i_r_point=1,n_r_coll	  
      READ(1) V_TRANS_BUFF
      IF(i_r_point.ge.i_nr_ini .and. i_r_point.le.i_nr_fin) THEN
      V_TRANSFER(:,:,:,:,:,i_r_point-i_nr_ini+1) = V_TRANS_BUFF
      IF(i_r_point.eq.1) PRINT*,"DATA READING STARTED"	 
      ENDIF	 
      ENDDO	
      PRINT*, "DATA BEEN READ"
      ENDIF		  
	  

      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      CALL MPI_BCAST(V_TRANSFER, buffer_size, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)		  
      CASE(7)
      ALLOCATE(
     & V_TRANSFER(1,n_gamma1,n_beta2,n_beta1,
     & n_alpha2,i_nr_fin-i_nr_ini+1),
     & V_TRANS_BUFF(1,n_gamma1,n_beta2,n_beta1,
     & n_alpha2))	  
      buffer_size = 1*n_gamma1*n_beta2*
     & n_beta1*n_alpha2*(i_nr_fin-i_nr_ini+1)  

      IF(MYID.EQ.0) THEN	  
      OPEN(1,FILE=potential_file_name,FORM="UNFORMATTED",
     & ACTION="READ",STATUS="OLD")	  
      READ(1) coll_type_checker
      READ(1)	na2_checker,nb1_checker,nb2_checker,ng1_checker,
     & nr_checker	  
      IF(coll_type_checker.ne.coll_type) STOP "WRONG COLL TYPE"
      IF(na2_checker.ne.n_alpha2) STOP "WRONG GRID"
      IF(nb1_checker.ne.n_beta1) STOP "WRONG GRID"
      IF(nb2_checker.ne.n_beta2) STOP "WRONG GRID"
      IF(ng1_checker.ne.n_gamma1) STOP "WRONG GRID"
      IF(nr_checker.ne.n_r_coll) STOP "WRONG GRID"
      PRINT*,"DATA READING STARTED"	 	  
      DO i_r_point=1,n_r_coll	  
      READ(1) V_TRANS_BUFF
      IF(i_r_point.ge.i_nr_ini .and. i_r_point.le.i_nr_fin) THEN
      V_TRANSFER(:,:,:,:,:,i_r_point-i_nr_ini+1) = V_TRANS_BUFF
      IF(i_r_point.eq.1) PRINT*,"DATA READING STARTED"	 
      ENDIF	 
      ENDDO	
      PRINT*, "DATA BEEN READ"
      ENDIF		  
	  

      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      CALL MPI_BCAST(V_TRANSFER, buffer_size, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)		  
      CASE(8)
      ALLOCATE(
     & V_TRANSFER(1,n_gamma1,n_beta2,n_beta1,
     & n_alpha2,i_nr_fin-i_nr_ini+1),
     & V_TRANS_BUFF(1,n_gamma1,n_beta2,n_beta1,
     & n_alpha2))	  
      buffer_size = 1*n_gamma1*n_beta2*
     & n_beta1*n_alpha2*(i_nr_fin-i_nr_ini+1)  

      IF(MYID.EQ.0) THEN	  
      OPEN(1,FILE=potential_file_name,FORM="UNFORMATTED",
     & ACTION="READ",STATUS="OLD")	  
      READ(1) coll_type_checker
      READ(1)	na2_checker,nb1_checker,nb2_checker,ng1_checker,
     & nr_checker	  
      IF(coll_type_checker.ne.coll_type) RETURN !!! ERROR MESSAGE
      IF(na2_checker.ne.n_alpha2) RETURN
      IF(nb1_checker.ne.n_beta1) RETURN
      IF(nb2_checker.ne.n_beta2) RETURN
      IF(ng1_checker.ne.n_gamma1) RETURN
      IF(nr_checker.ne.n_r_coll) RETURN
      PRINT*,"DATA READING STARTED"	 	  
      DO i_r_point=1,n_r_coll	  
      READ(1) V_TRANS_BUFF(1,:,:,:,:)
      IF(i_r_point.ge.i_nr_ini .and. i_r_point.le.i_nr_fin) THEN
      V_TRANSFER(:,:,:,:,:,i_r_point-i_nr_ini+1) = V_TRANS_BUFF
      IF(i_r_point.eq.1) PRINT*,"DATA READING STARTED"	 
      ENDIF	 
      ENDDO	
      PRINT*, "DATA BEEN READ"
      ENDIF		  
	  

      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      CALL MPI_BCAST(V_TRANSFER, buffer_size, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)		  
      CASE(9)
      ALLOCATE(
     & V_TRANSFER(n_gamma2,n_gamma1,n_beta2,n_beta1,
     & n_alpha2,i_nr_fin-i_nr_ini+1),
     & V_TRANS_BUFF(n_gamma2,n_gamma1,n_beta2,n_beta1,
     & n_alpha2))	  
      buffer_size = n_gamma2*n_gamma1*n_beta2*
     & n_beta1*n_alpha2*(i_nr_fin-i_nr_ini+1)  

      IF(MYID.EQ.0) THEN	  
      OPEN(1,FILE=potential_file_name,FORM="UNFORMATTED",
     & ACTION="READ",STATUS="OLD")	  
      READ(1) coll_type_checker
      READ(1)	na2_checker,nb1_checker,nb2_checker,ng1_checker,
     & ng2_checker,nr_checker	  
      IF(coll_type_checker.ne.coll_type) STOP "WRONG COLL TYPE"
      IF(na2_checker.ne.n_alpha2) STOP "WRONG GRID"
      IF(nb1_checker.ne.n_beta1) STOP "WRONG GRID"
      IF(nb2_checker.ne.n_beta2) STOP "WRONG GRID"
      IF(ng1_checker.ne.n_gamma1) STOP "WRONG GRID"
      IF(ng2_checker.ne.n_gamma2) STOP "WRONG GRID"
      IF(nr_checker.ne.n_r_coll) STOP "WRONG GRID"
      PRINT*,"DATA READING STARTED"	 	  
      DO i_r_point=1,n_r_coll	  
      READ(1) V_TRANS_BUFF
      IF(i_r_point.ge.i_nr_ini .and. i_r_point.le.i_nr_fin) THEN
      V_TRANSFER(:,:,:,:,:,i_r_point-i_nr_ini+1) = V_TRANS_BUFF
      IF(i_r_point.eq.1) PRINT*,"DATA READING STARTED"	 
      ENDIF	 
      ENDDO	
      PRINT*, "DATA BEEN READ"
      ENDIF		  
	  

      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      CALL MPI_BCAST(V_TRANSFER, buffer_size, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)		  
      CASE(0)	  
      ALLOCATE(
     & V_TRANSFER(n_gamma2,n_gamma1,n_beta2,n_beta1,
     & n_alpha2,i_nr_fin-i_nr_ini+1),
     & V_TRANS_BUFF(n_gamma2,n_gamma1,n_beta2,n_beta1,
     & n_alpha2))	  
      buffer_size = n_gamma2*n_gamma1*n_beta2*
     & n_beta1*n_alpha2*(i_nr_fin-i_nr_ini+1)  

      IF(MYID.EQ.0) THEN	  
      OPEN(1,FILE=potential_file_name,FORM="UNFORMATTED",
     & ACTION="READ",STATUS="OLD")	  
      READ(1) coll_type_checker
      READ(1)	na2_checker,nb1_checker,nb2_checker,ng1_checker,
     & ng2_checker,nr_checker	  
      IF(coll_type_checker.ne.coll_type) STOP "WRONG COLL TYPE"
      IF(na2_checker.ne.n_alpha2) STOP "WRONG GRID"
      IF(nb1_checker.ne.n_beta1) STOP "WRONG GRID"
      IF(nb2_checker.ne.n_beta2) STOP "WRONG GRID"
      IF(ng1_checker.ne.n_gamma1) STOP "WRONG GRID"
      IF(ng2_checker.ne.n_gamma2) STOP "WRONG GRID"
      IF(nr_checker.ne.n_r_coll) STOP "WRONG GRID"
      PRINT*,"DATA READING STARTED"	 	  
      DO i_r_point=1,n_r_coll	  
      READ(1) V_TRANS_BUFF
      IF(i_r_point.ge.i_nr_ini .and. i_r_point.le.i_nr_fin) THEN
      V_TRANSFER(:,:,:,:,:,i_r_point-i_nr_ini+1) = V_TRANS_BUFF
      IF(i_r_point.eq.1) PRINT*,"DATA READING STARTED"	 
      ENDIF	 
      ENDDO	
      PRINT*, "DATA BEEN READ"
      ENDIF		  
	  

      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      CALL MPI_BCAST(V_TRANSFER, buffer_size, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)		  
      END SELECT   	  
      IF(MYID.EQ.0) PRINT*,"BUFFER HAS BEEN BROADCASTED"
      IF(MYID.eq.0) CLOSE(1)	  
      ENDIF        
      IF(angs_unit) conv_unit_r = a_bohr
      IF(MYID.EQ.0) PRINT*,"EXPANSION TERMS COMPUTING STARTED"
      IF(MYID.EQ.0) PRINT*,"PROGRESS WILL BE PRINTED INTO 
     &PROGRESS_EXPAN.tmp FILE"								!Bikram Jan,'19
      DO i=myid+1,myid+1
	  bk_tym_1=MPI_Wtime()
      IF(i .le. N_TERMS_EXP) THEN	  
      R_GRID: DO i_r_point=i_nr_ini,i_nr_fin
      intgrlr = 0d0
      intgrli = 0d0
      weig_b = 0d0	  
      SELECT CASE(coll_type)
      CASE(1)
	  
	  if(.not.bikram_identical_pes) then

      DO i2=1,n_beta
      beta = xbm + xbl*xb(i2)
      R =  R_COM(i_r_point)*conv_unit_r
      l_t = A_TOP(1,i)
	  
      CALL EXPANSION_TERM
     & (l_t,0,0,0,0,wfr,wfi,
     & beta,0,0,0,0,coll_type,bikram_w3j_fact)
      V = V_2_1(i2,i_r_point)
	  
      integrant	 = (wfr*V)*wb(i2)*dsin(beta)
      intgrlr  = intgrlr + integrant*xbl*2d0*pi	  
      ENDDO
	  
	  else
	  print*, "Something is wrong in PES expansion."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if

! Bikram Start: This piece is to compute the RMSE wrt the original PES
! Bikram End.

      CASE(3)
	  
	  l_t = A_TOP(1,i)
	  nju1_t = A_TOP(2,i)
	  
	  if(.not.bikram_identical_pes) then
	  
	  same_term = .false.
	  if(nju1_t.lt.0) then
	  write(*,'(a, a, 5(i0,a))')"The values of M1 should be ",
     & "always positive. Please check the expansion term: ",
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_term = .true.
	  end if
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  if(same_term) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if

	  else
	  print*, "Something is wrong in PES expansion."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if

	  
      DO i2=1,n_beta
      DO i3=1,n_gamma
      beta = xbm + xbl*xb(i2)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i3)
	  else
	  gamma = xg(i3)
	  end if
! Bikram End.

      R =  R_COM(i_r_point)*conv_unit_r
	  
      CALL EXPANSION_TERM
     & (l_t,0,0,nju1_t,0,wfr,wfi,
     & beta,gamma,0,0,0,coll_type,bikram_w3j_fact)							
      V = V_3_1(i3,i2,i_r_point)
	  
       integrant = wfr*V*wb(i2)*wg(i3)*dsin(beta)
      intgrlr  = intgrlr + integrant*xbl*xgl
      ENDDO
      ENDDO
	  
! Bikram Start: This piece is to compute the RMSE wrt the original PES
! Bikram End.
	  
      CASE(4)
	  
	  l_t = A_TOP(1,i)
	  nju1_t = A_TOP(2,i)	  
	  
	  if(.not.bikram_identical_pes) then
	  
	  same_term = .false.
	  if(nju1_t.lt.0) then
	  write(*,'(a, a, 5(i0,a))')"The values of M1 should be ",
     & "always positive. Please check the expansion term: ",
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_term = .true.
	  end if
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  if(same_term) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if

	  else
	  print*, "Something is wrong in PES expansion."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
      DO i2=1,n_beta
      DO i3=1,n_gamma
      beta = xbm + xbl*xb(i2)
!      gamma = xgm + xgl*xg(i3)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i3)
	  else
	  gamma = xg(i3)
	  end if
! Bikram End.
      R =  R_COM(i_r_point)*conv_unit_r
	  
      CALL EXPANSION_TERM
     & (l_t,0,0,nju1_t,0,wfr,wfi,
     & beta,gamma,0,0,0,coll_type,bikram_w3j_fact)
      V = V_3_1(i3,i2,i_r_point)
       integrant = wfr*V*wb(i2)*wg(i3)*dsin(beta)
      intgrlr  = intgrlr + integrant*xbl*xgl
      ENDDO
      ENDDO

! Bikram Start: This piece is to compute the RMSE wrt the original PES
! Bikram End.
	  
      CASE(5)	  
	  
	  if(.not.bikram_identical_pes) then
! nothing to do here since there is L1, L2, and L only
	  
	  else if(bikram_identical_pes) then
      l1_t = A_TOP(1,i)
      l2_t = A_TOP(2,i)
      l_t = A_TOP(3,i)
	  
	  same_term = .false.
	  if(l1_t.lt.l2_t) then
	  write(*,'(a, a, a, 5(i0,a))')"The values of L2 should be ",
     & "less or equal to L1. ", 
     & "Please check the expansion term: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_term = .true.
	  end if
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  if(same_term) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if	  

	  else
	  print*, "Something is wrong in PES expansion."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if

      DO i2=1,n_beta1
      DO i3=1,n_beta2
      DO i4=1,n_alpha2
      beta = xbm + xbl*xb(i2)
      bbeta = xbbm + xbbl*xbb(i3)
!      gamma = xgm + xgl*xg(i4)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
!      gamma = xgm + xgl*xg(i4)
      aalpha = xaam + xaal*xaa(i4)
	  else
!	  gamma = xg(i4)
	  aalpha = xaa(i4)
	  end if
! Bikram End.

      R =  R_COM(i_r_point)*conv_unit_r
      l1_t = A_TOP(1,i)
      l2_t = A_TOP(2,i)
      l_t = A_TOP(3,i)
  
       CALL EXPANSION_TERM
     & (l_t,l1_t,l2_t,0,0,wfr,wfi,
     & beta,0.0d0,aalpha,bbeta,0.0d0,coll_type,bikram_w3j_fact)
      V = V_2_2(i4,i3,i2,i_r_point)
	  
      integrant	 = (wfr*V)*wb(i2)*wbb(i3)*dsin(beta)*
     & waa(i4)*dsin(bbeta)
      intgrlr  = intgrlr + integrant*xbl*xaal
     & *xbbl*2d0*pi	 
      ENDDO
      ENDDO
      ENDDO
	  
! Bikram Start: This piece is to compute the RMSE wrt the original PES
! Bikram End.
	  
      CASE(7)
	  
      l1_t = A_TOP(1,i)
      nju1_t = A_TOP(2,i)
      l2_t = A_TOP(3,i)
      l_t = A_TOP(4,i) 
      IF(max(l1_t,l2_t,l_t)*2 .gt. l1_t+l2_t+l_t) THEN
      PRINT*,l1_t,l2_t,l_t
      STOP "TRINAGULAR RULE FAILED"
      ENDIF
	  
	  if(.not.bikram_identical_pes) then
	  
	  same_term = .false.
	  if(nju1_t.lt.0) then
	  write(*,'(a, a, 5(i0,a))')"The values of M1 should be ",
     & "always positive. Please check the expansion term: ",
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_term = .true.
	  end if
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  if(same_term) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
	  do nmb_terms = 1, nterms
      l1bk = A_TOP(1,nmb_terms)
      n1bk = A_TOP(2,nmb_terms)
      l2bk = A_TOP(3,nmb_terms)
      lbk = A_TOP(4,nmb_terms)
	  
      do planr_coeff = 1, 2 - KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	 
      if(planr_coeff.eq.2) then
      nju1_t = -nju1_t
      end if
	  
	  if(nmb_terms.ne.i) then
	  if(l1bk.eq.l1_t .and. n1bk.eq.nju1_t .and. l2bk.eq.l2_t .and. 
     & lbk.eq.l_t) then
	  write(*,'(a, a, a, 4(i0,a),5x,4(i0,a))')"Found two same terms.",
     & "Please check the expansion terms: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', l_t,
     & l1bk, ',', n1bk, ',', l2bk, ',', lbk
	  same_term = .true.
	  end if
	  end if
	  end do
	  end do
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  if(same_term) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
	  else
	  print*, "Something is wrong in PES expansion."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  	  	  
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      beta = xbm + xbl*xb(i3)
      bbeta = xbbm + xbbl*xbb(i4)
!      gamma = xgm + xgl*xg(i5)
!      aalpha =  xaam + xaal*xaa(i2)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
	  else
	  gamma = xg(i5)
	  aalpha = xaa(i2)
	  end if
! Bikram End.
      R =  R_COM(i_r_point)*conv_unit_r
	  
      l1_t = A_TOP(1,i)
      nju1_t = A_TOP(2,i)
      l2_t = A_TOP(3,i)
      l_t = A_TOP(4,i)
      IF(max(l1_t,l2_t,l_t)*2 .gt. l1_t+l2_t+l_t) THEN
      PRINT*,l1_t,l2_t,l_t
      STOP "TRINAGULAR RULE FAILED"
      ENDIF
      
      CALL EXPANSION_TERM
     & (l_t,l1_t,l2_t,nju1_t,nju2_t,wfr,wfi,
     & beta,gamma,aalpha,bbeta,ggamma,coll_type,bikram_w3j_fact)
      IF(make_grid_file) THEN
      IF(grid_file_found) THEN	  
      V = V_TRANSFER(1,i5,i4,i3,
     & i2,i_r_point-i_nr_ini+1) 
      ELSE
      IF(MYID.eq.0) PRINT*, "ERROR:GRID_FILE_NOT_FOUND"	  
	  STOP
      ENDIF	  
      ELSE
      V = V_3_2(i5,i4,i3,
     & i2,i_r_point-i_nr_ini+1)	  
      ENDIF	  
      integrant	 = V*wfr*wb(i3)*wg(i5)*dsin(beta)*
     & waa(i2)*wbb(i4)*dsin(bbeta)
      intgrlr  = intgrlr + integrant*xbl*xgl
     & *xaal*xbbl
      integrant	 = V*wfi*wb(i3)*wg(i5)*dsin(beta)*
     & waa(i2)*wbb(i4)*dsin(bbeta)
      intgrli  = intgrli + integrant*xbl*xgl
     & *xaal*xbbl
	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO	 
	  
! Bikram Start: This piece is to compute the RMSE wrt the original PES
! Bikram End.
	  
      CASE(8)
	  
      l1_t = A_TOP(1,i)
      nju1_t = A_TOP(2,i)
      l2_t = A_TOP(3,i)
      l_t = A_TOP(4,i) 
      IF(max(l1_t,l2_t,l_t)*2 .gt. l1_t+l2_t+l_t) THEN
      PRINT*,l1_t,l2_t,l_t
      STOP "TRINAGULAR RULE FAILED"
      ENDIF
	  
	  if(.not.bikram_identical_pes) then
	  
	  same_term = .false.
	  if(nju1_t.lt.0) then
	  write(*,'(a, a, 5(i0,a))')"The values of M1 should be ",
     & "always positive. Please check the expansion term: ",
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_term = .true.
	  end if
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  if(same_term) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
	  do nmb_terms = 1, nterms
      l1bk = A_TOP(1,nmb_terms)
      n1bk = A_TOP(2,nmb_terms)
      l2bk = A_TOP(3,nmb_terms)
      lbk = A_TOP(4,nmb_terms)
	  
      do planr_coeff = 1, 2 - KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	 
      if(planr_coeff.eq.2) then
      nju1_t = -nju1_t
      end if
	  
	  if(nmb_terms.ne.i) then
	  if(l1bk.eq.l1_t .and. n1bk.eq.nju1_t .and. l2bk.eq.l2_t .and. 
     & lbk.eq.l_t) then
	  write(*,'(a, a, a, 4(i0,a),5x,4(i0,a))')"Found two same terms.",
     & "Please check the expansion terms: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', l_t,
     & l1bk, ',', n1bk, ',', l2bk, ',', lbk
	  same_term = .true.
	  end if
	  end if
	  end do
	  end do
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  if(same_term) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
	  else
	  print*, "Something is wrong in PES expansion."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  	  	  
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      beta = xbm + xbl*xb(i3)
      bbeta = xbbm + xbbl*xbb(i4)
!      gamma = xgm + xgl*xg(i5)
!      aalpha =  xaam + xaal*xaa(i2)
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
	  else
	  gamma = xg(i5)
	  aalpha = xaa(i2)
	  end if
! Bikram End.
      R =  R_COM(i_r_point)*conv_unit_r
	  
      l1_t = A_TOP(1,i)
      nju1_t = A_TOP(2,i)
      l2_t = A_TOP(3,i)
      l_t = A_TOP(4,i) 
      IF(max(l1_t,l2_t,l_t)*2 .gt. l1_t+l2_t+l_t) THEN
      PRINT*,l1_t,l2_t,l_t
      STOP "TRINAGULAR RULE FAILED"
      ENDIF
 
      CALL EXPANSION_TERM
     & (l_t,l1_t,l2_t,nju1_t,nju2_t,wfr,wfi,
     & beta,gamma,aalpha,bbeta,0d0,coll_type,bikram_w3j_fact)
      IF(make_grid_file) THEN
      IF(grid_file_found) THEN	  
      V = V_TRANSFER(1,i5,i4,i3,
     & i2,i_r_point-i_nr_ini+1) 
      ELSE
      IF(MYID.eq.0) PRINT*, "ERROR:GRID_FILE_NOT_FOUND"	  
	  STOP
      ENDIF	  
      ELSE
      V = V_3_2(i5,i4,i3,
     & i2,i_r_point-i_nr_ini+1)	  
      ENDIF	
      integrant	 = (V*wfr)*wb(i3)*wg(i5)*dsin(beta)*
     & waa(i2)*wbb(i4)*dsin(bbeta)
      intgrlr  = intgrlr + integrant*xbl*xgl	  
     & *xaal*xbbl
      integrant	 = V*wfi*wb(i3)*wg(i5)*dsin(beta)*
     & waa(i2)*wbb(i4)*dsin(bbeta)
      intgrli  = intgrli + integrant*xbl*xgl
     & *xaal*xbbl
	 
      ENDDO
      ENDDO
      ENDDO
      ENDDO
	  
! Bikram Start: This piece is to compute the RMSE wrt the original PES
! Bikram End.

      CASE(9)
	  
      l1_t = A_TOP(1,i)
      nju1_t = A_TOP(2,i)
      l2_t = A_TOP(3,i)
      nju2_t = A_TOP(4,i)
      l_t = A_TOP(5,i)
      IF(max(l1_t,l2_t,l_t)*2 .gt. l1_t+l2_t+l_t) THEN
      PRINT*,l1_t,l2_t,l_t
      STOP "TRINAGULAR RULE FAILED"
      ENDIF	  
	  
	  if(.not.bikram_identical_pes) then
	  
	  same_term = .false.
	  if(nju1_t.eq.0 .and. nju2_t.lt.0) then
	  write(*,'(a, a, a, 5(i0,a))')"The values of M2 should be ",
     & "positive only for M1 = 0. ", 
     & "Please check the expansion term: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_term = .true.
	  end if
	  if(nju1_t.lt.0) then
	  write(*,'(a, a, 5(i0,a))')"The values of M1 should be ",
     & "always positive. Please check the expansion term: ",
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_term = .true.
	  end if
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  if(same_term) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
	  do nmb_terms = 1, nterms
      l1bk = A_TOP(1,nmb_terms)
      n1bk = A_TOP(2,nmb_terms)
      l2bk = A_TOP(3,nmb_terms)
      n2bk = A_TOP(4,nmb_terms)
      lbk = A_TOP(5,nmb_terms)
	  
      do planr_coeff = 1, 2 - KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	 
      if(planr_coeff.eq.2) then
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      end if
	  
	  if(nmb_terms.ne.i) then
	  if(l1bk.eq.l1_t .and. n1bk.eq.nju1_t .and. l2bk.eq.l2_t .and. 
     & n2bk.eq.nju2_t .and. lbk.eq.l_t) then
	  write(*,'(a, a, a, 5(i0,a),5x,5(i0,a))')"Found two same terms.",
     & "Please check the expansion terms: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t,
     & l1bk, ',', n1bk, ',', l2bk, ',', n2bk, ',', lbk
	  same_term = .true.
	  end if
	  end if
	  end do
	  end do
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  if(same_term) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
	  else
	  print*, "Something is wrong in PES expansion."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
 
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      DO i6=1,n_gamma2
      beta = xbm + xbl*xb(i3)
      bbeta = xbbm + xbbl*xbb(i4)
!      gamma = xgm + xgl*xg(i5)
!      aalpha =  xaam + xaal*xaa(i2)
!      ggamma = xggm + xggl*xgg(i6)	  
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
      ggamma = xggm + xggl*xgg(i6)	  
	  else
	  gamma = xg(i5)
	  aalpha = xaa(i2)
	  ggamma = xgg(i6)
	  end if
! Bikram End.
      R =  R_COM(i_r_point)*conv_unit_r
      l1_t = A_TOP(1,i)
      nju1_t = A_TOP(2,i)
      l2_t = A_TOP(3,i)
      nju2_t = A_TOP(4,i)
      l_t = A_TOP(5,i) 
      IF(max(l1_t,l2_t,l_t)*2 .gt. l1_t+l2_t+l_t) THEN
      PRINT*,l1_t,l2_t,l_t
      STOP "TRINAGULAR RULE FAILED"
      ENDIF
 
      CALL EXPANSION_TERM
     & (l_t,l1_t,l2_t,nju1_t,nju2_t,wfr,wfi,
     & beta,gamma,aalpha,bbeta,ggamma,coll_type,bikram_w3j_fact)
	 
      IF(make_grid_file) THEN
      IF(grid_file_found) THEN	  
      V = V_TRANSFER(i6,i5,i4,i3,
     & i2,i_r_point-i_nr_ini+1)
      ELSE
      IF(MYID.eq.0) PRINT*, "ERROR:GRID_FILE_NOT_FOUND"	  
	  STOP
      ENDIF	  
      ELSE
      CALL V_POT_ASYM_SYM(V,R,0d0,beta,gamma,
     & aalpha,bbeta,ggamma) 
      ENDIF		 
	 	
      integrant	 = V*wfr*wb(i3)*wg(i5)*dsin(beta)*
     & waa(i2)*wbb(i4)*wgg(i6)*dsin(bbeta)
      intgrlr  = intgrlr + integrant*xbl*xgl
     & *xaal*xbbl*xggl*2*pi
      integrant	 = V*wfi*wb(i3)*wg(i5)*dsin(beta)*
     & waa(i2)*wbb(i4)*wgg(i6)*dsin(bbeta)
      intgrli  = intgrli + integrant*xbl*xgl
     & *xaal*xbbl*xggl*2*pi
	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO	  

! Bikram Start: This piece is to compute the RMSE wrt the original PES
! Bikram End.	  

      CASE(0)

      l1_t = A_TOP(1,i)
      nju1_t = A_TOP(2,i)
      l2_t = A_TOP(3,i)
      nju2_t = A_TOP(4,i)
      l_t = A_TOP(5,i) 
      IF(max(l1_t,l2_t,l_t)*2 .gt. l1_t+l2_t+l_t) THEN
      PRINT*,l1_t,l2_t,l_t
      STOP "TRINAGULAR RULE FAILED"
      ENDIF
	  
	  if(.not.bikram_identical_pes) then
	  
	  same_term = .false.
	  if(nju1_t.eq.0 .and. nju2_t.lt.0) then
	  write(*,'(a, a, a, 5(i0,a))')"The values of M2 should be ",
     & "positive only for M1 = 0. ", 
     & "Please check the expansion term: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_term = .true.
	  end if
	  if(nju1_t.lt.0) then
	  write(*,'(a, a, 5(i0,a))')"The values of M1 should be ",
     & "always positive. Please check the expansion term: ",
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_term = .true.
	  end if
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  if(same_term) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
	  do nmb_terms = 1, nterms
      l1bk = A_TOP(1,nmb_terms)
      n1bk = A_TOP(2,nmb_terms)
      l2bk = A_TOP(3,nmb_terms)
      n2bk = A_TOP(4,nmb_terms)
      lbk = A_TOP(5,nmb_terms)
	  
      do planr_coeff = 1, 2 - KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	 
      if(planr_coeff.eq.2) then
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      end if
	  
	  if(nmb_terms.ne.i) then
	  if(l1bk.eq.l1_t .and. n1bk.eq.nju1_t .and. l2bk.eq.l2_t .and. 
     & n2bk.eq.nju2_t .and. lbk.eq.l_t) then
	  write(*,'(a, a, a, 5(i0,a),5x,5(i0,a))')"Found two same terms.",
     & "Please check the expansion terms: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t,
     & l1bk, ',', n1bk, ',', l2bk, ',', n2bk, ',', lbk
	  same_term = .true.
	  end if
	  end if
	  end do
	  end do
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  if(same_term) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
	  else if(bikram_identical_pes) then
	  
	  same_term = .false.
	  if(nju1_t.eq.0 .and. nju2_t.lt.0) then
	  write(*,'(a, a, a, 5(i0,a))')"The values of M2 should be ",
     & "positive only for M1 = 0. ", 
     & "Please check the expansion term: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_term = .true.
	  end if
	  if(nju1_t.lt.0) then
	  write(*,'(a, a, 5(i0,a))')"The values of M1 should be ",
     & "always positive. Please check the expansion term: ",
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_term = .true.
	  end if
	  if(l1_t.lt.l2_t) then
	  write(*,'(a, a, a, 5(i0,a))')"The values of L2 should be ",
     & "less or equal to L1. ", 
     & "Please check the expansion term: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_term = .true.
	  end if
	  if(l1_t.eq.l2_t .and. abs(nju1_t).lt.abs(nju2_t)) then
	  write(*,'(a, a, a, 5(i0,a))')"The values of M1 should be ",
     & "greater or equal to abs(M2) for L1 = L2. ", 
     & "Please check the expansion term: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_term = .true.
	  end if
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  if(same_term) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if	  
	  
	  do nmb_terms = 1, nterms
      l1bk = A_TOP(1,nmb_terms)
      n1bk = A_TOP(2,nmb_terms)
      l2bk = A_TOP(3,nmb_terms)
      n2bk = A_TOP(4,nmb_terms)
      lbk = A_TOP(5,nmb_terms)      
	  
	  do symmetry_coeff = 1, 
     & 2 - KRONEKER(l1_t,l2_t)*KRONEKER(abs(nju1_t),abs(nju2_t))
      if(symmetry_coeff.eq.2) then
      l2_t = A_TOP(1,i)
      nju2_t = A_TOP(2,i)
      l1_t = A_TOP(3,i)
      nju1_t = A_TOP(4,i)
      l_t = A_TOP(5,i)	  
      ENDIF	  
	  
      do planr_coeff = 1, 2 - KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	 
      if(planr_coeff.eq.2) then
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      end if
	  
	  if(nmb_terms.ne.i) then
	  if(l1bk.eq.l1_t .and. n1bk.eq.nju1_t .and. l2bk.eq.l2_t .and. 
     & n2bk.eq.nju2_t .and. lbk.eq.l_t) then
	  write(*,'(a, a, a, 5(i0,a),5x,5(i0,a))')"Found two same terms.",
     & "Please check the expansion terms: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t,
     & l1bk, ',', n1bk, ',', l2bk, ',', n2bk, ',', lbk
	  same_term = .true.
	  end if
	  end if
	  end do
	  end do
	  end do
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  if(same_term) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
	  else
	  print*, "Something is wrong in PES expansion."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  	  
      DO i2=1,n_alpha2
      DO i3=1,n_beta1
      DO i4=1,n_beta2
      DO i5=1,n_gamma1
      DO i6=1,n_gamma2
      beta = xbm + xbl*xb(i3)
      bbeta = xbbm + xbbl*xbb(i4)
!      gamma = xgm + xgl*xg(i5)
!      aalpha =  xaam + xaal*xaa(i2)
!      ggamma = xggm + xggl*xgg(i6)	  
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
      ggamma = xggm + xggl*xgg(i6)	  
	  else
	  gamma = xg(i5)
	  aalpha = xaa(i2)
	  ggamma = xgg(i6)
	  end if
! Bikram End.
      R =  R_COM(i_r_point)*conv_unit_r
	  l1_t = A_TOP(1,i)
      nju1_t = A_TOP(2,i)
      l2_t = A_TOP(3,i)
      nju2_t = A_TOP(4,i)
      l_t = A_TOP(5,i)
 
      CALL EXPANSION_TERM
     & (l_t,l1_t,l2_t,nju1_t,nju2_t,wfr,wfi,
     & beta,gamma,aalpha,bbeta,ggamma,coll_type,bikram_w3j_fact)


      IF(make_grid_file) THEN
      IF(grid_file_found) THEN	  
      V = V_TRANSFER(i6,i5,i4,i3,i2,i_r_point-i_nr_ini+1)
      ELSE
      IF(MYID.eq.0) PRINT*, "ERROR:GRID_FILE_NOT_FOUND"	  
	  STOP
      ENDIF	  
      ELSE
      CALL V_POT_ASYM_ASYM(V,R,0d0,beta,gamma,aalpha,bbeta,ggamma)
      ENDIF		 
	  
      integrant	 = V*wfr*wb(i3)*wg(i5)*dsin(beta)*
     & waa(i2)*wbb(i4)*wgg(i6)*dsin(bbeta)
      intgrlr  = intgrlr + integrant*xbl*xgl
     & *xaal*xbbl*xggl*2*pi
      integrant	 = V*wfi*wb(i3)*wg(i5)*dsin(beta)*
     & waa(i2)*wbb(i4)*wgg(i6)*dsin(bbeta)
      intgrli  = intgrli + integrant*xbl*xgl
     & *xaal*xbbl*xggl*2*pi
	 
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
	  
! Bikram Start: This piece is to compute the RMSE wrt the original PES
	  if(rms_defined .and. i_r_point.eq.rms_r) then
	  if(bikram_rms_ang1.gt.n_alpha2) bikram_rms_ang1 = n_alpha2
	  if(bikram_rms_ang2.gt.n_beta1)  bikram_rms_ang2 = n_beta1
	  if(bikram_rms_ang3.gt.n_gamma1)  bikram_rms_ang3 = n_gamma1
	  
	  grd_tot = bikram_rms_ang1*bikram_rms_ang2*bikram_rms_ang2*
     & bikram_rms_ang3*bikram_rms_ang3
	  allocate(V_store(grd_tot), V_recalc(grd_tot))
	  if(.not.bikram_identical_pes) then
	  allocate(wfr_store(grd_tot,nproc), wfr_tmp(grd_tot), 
     & wfi_store(grd_tot,nproc), wfi_tmp(grd_tot))
	  else if(bikram_identical_pes) then
	  allocate(wfr_store(grd_tot,nproc), wfr_tmp(grd_tot), 
     & wfi_store(grd_tot,nproc), wfi_tmp(grd_tot))
	  else
	  print*, "Something is wrong in PES expansion."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  grd_cntr = 0
	  wfr_tmp = 0
	  wfi_tmp = 0
	  
      DO i2=1,bikram_rms_ang1
      DO i3=1,bikram_rms_ang2
      DO i4=1,bikram_rms_ang2
      DO i5=1,bikram_rms_ang3
      DO i6=1,bikram_rms_ang3
      beta = xbm + xbl*xb(i3)
      bbeta = xbbm + xbbl*xbb(i4)
!      gamma = xgm + xgl*xg(i5)
!      aalpha =  xaam + xaal*xaa(i2)
!      ggamma = xggm + xggl*xgg(i6)	  
	  
! Bikram Start Nov 2021: Making the equidistant grid for angle gamma for this system type	  	  
	  if(.not. bikram_eq_grd_gamma) then
      gamma = xgm + xgl*xg(i5)
      aalpha =  xaam + xaal*xaa(i2)
      ggamma = xggm + xggl*xgg(i6)	  
	  else
	  gamma = xg(i5)
	  aalpha = xaa(i2)
	  ggamma = xgg(i6)
	  end if
! Bikram End.
      R =  R_COM(i_r_point)*conv_unit_r

      IF(make_grid_file) THEN
      IF(grid_file_found) THEN	  
      V = V_TRANSFER(i6,i5,i4,i3,i2,i_r_point-i_nr_ini+1)
      ELSE
      IF(MYID.eq.0) PRINT*, "ERROR:GRID_FILE_NOT_FOUND"	  
	  STOP
      ENDIF	  
      ELSE
      CALL V_POT_ASYM_ASYM(V,R,0d0,beta,gamma,aalpha,bbeta,ggamma)
      ENDIF	
	  
	  grd_cntr = grd_cntr + 1
	  V_store(grd_cntr) = V
 
 	  if(.not.bikram_identical_pes) then
      CALL EXPANSION_TERM
     & (l_t,l1_t,l2_t,nju1_t,nju2_t,wfr,wfi,
     & beta,gamma,aalpha,bbeta,ggamma,coll_type,bikram_w3j_fact)
	  wfr_tmp(grd_cntr) = wfr_tmp(grd_cntr) + wfr*intgrlr
	  wfi_tmp(grd_cntr) = wfi_tmp(grd_cntr) + wfi*intgrlr
      CALL EXPANSION_TERM
     & (l_t,l1_t,l2_t,-nju1_t,-nju2_t,wfr,wfi,
     & beta,gamma,aalpha,bbeta,ggamma,coll_type,bikram_w3j_fact)
	  wfr_tmp(grd_cntr) = wfr_tmp(grd_cntr) + 
     & wfr*(-1)**(l_t+l1_t+l2_t+nju1_t+nju2_t)*intgrlr
	  wfi_tmp(grd_cntr) = wfi_tmp(grd_cntr) + 
     & wfi*(-1)**(l_t+l1_t+l2_t+nju1_t+nju2_t)*intgrlr
	  else if(bikram_identical_pes) then
      CALL EXPANSION_TERM
     & (l_t,l1_t,l2_t,nju1_t,nju2_t,wfr,wfi,
     & beta,gamma,aalpha,bbeta,ggamma,coll_type,bikram_w3j_fact)
	  wfr_tmp(grd_cntr) = wfr_tmp(grd_cntr) + wfr*intgrlr
	  wfi_tmp(grd_cntr) = wfi_tmp(grd_cntr) + wfi*intgrlr
      CALL EXPANSION_TERM
     & (l_t,l1_t,l2_t,-nju1_t,-nju2_t,wfr,wfi,
     & beta,gamma,aalpha,bbeta,ggamma,coll_type,bikram_w3j_fact)
	  wfr_tmp(grd_cntr) = wfr_tmp(grd_cntr) + 
     & wfr*(-1)**(l_t+l1_t+l2_t+nju1_t+nju2_t)*intgrlr
	  wfi_tmp(grd_cntr) = wfi_tmp(grd_cntr) + 
     & wfi*(-1)**(l_t+l1_t+l2_t+nju1_t+nju2_t)*intgrlr
      CALL EXPANSION_TERM
     & (l_t,l2_t,l1_t,nju2_t,nju1_t,wfr,wfi,
     & beta,gamma,aalpha,bbeta,ggamma,coll_type,bikram_w3j_fact)
	  wfr_tmp(grd_cntr) = wfr_tmp(grd_cntr) + 
     & wfr*(-1)**(l1_t+l2_t)*intgrlr
	  wfi_tmp(grd_cntr) = wfi_tmp(grd_cntr) + 
     & wfi*(-1)**(l1_t+l2_t)*intgrlr
      CALL EXPANSION_TERM
     & (l_t,l2_t,l1_t,-nju2_t,-nju1_t,wfr,wfi,
     & beta,gamma,aalpha,bbeta,ggamma,coll_type,bikram_w3j_fact)
	  wfr_tmp(grd_cntr) = wfr_tmp(grd_cntr) + 
     & wfr*(-1)**(l_t+nju1_t+nju2_t)*intgrlr
	  wfi_tmp(grd_cntr) = wfi_tmp(grd_cntr) + 
     & wfi*(-1)**(l_t+nju1_t+nju2_t)*intgrlr
	  end if
	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
	  
	  if(grd_cntr.ne.grd_tot) then
	  write(*,*)"Error in storing V, and Expansion functions"
	  stop
	  end if
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  
      CALL MPI_GATHER(wfr_tmp,grd_tot,MPI_DOUBLE_PRECISION,
     & wfr_store,grd_tot,MPI_DOUBLE_PRECISION,
     & 0,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_GATHER(wfi_tmp,grd_tot,MPI_DOUBLE_PRECISION,
     & wfi_store,grd_tot,MPI_DOUBLE_PRECISION,
     & 0,MPI_COMM_WORLD,ierr_mpi)
	  
	  if(myid.eq.0) then
	  V_recalc = 0.0d0
	  v_sum = 0.0d0
	  do grd_cntr = 1, grd_tot
	  do i2 = 1, nproc
	  V_recalc(grd_cntr) = V_recalc(grd_cntr) + wfr_store(grd_cntr,i2)
	  end do
	  v_sum = v_sum + (V_store(grd_cntr) - V_recalc(grd_cntr))**2d0
	  end do
	  rms = dsqrt(v_sum/grd_tot)
	  print*, "RMSE = ", rms
	  end if
	  end if
! Bikram End.	  

      END SELECT      	  
      buffer(i_r_point)	 =  intgrlr
      buffer_i(i_r_point)	 =  intgrli	  
      IF(i_r_point.eq.1) buffer_w = weig_b      	  
	  bk_tym_2=MPI_Wtime()
	  bk_tym=bk_tym_2-bk_tym_1
	  open(143,file='PROGRESS_EXPAN.tmp',access='append')
	  write(143,'(a19,i5,a1,5x,a15,i4,a1,5x,a17,f15.3,a3)')
     &"Done: R-grid point#",i_r_point,",","expansion term#",(myid+1),
     &",","time since start:",bk_tym,"sec"		!Bikram Jan,'19
	  close(143)		!Bikram Jan,'19
      ENDDO R_GRID
      ENDIF	  
      ENDDO	
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      CALL MPI_GATHER(buffer,n_r_coll,
     & MPI_DOUBLE_PRECISION,expansion_terms,
     &	  n_r_coll,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_GATHER(buffer_i,n_r_coll,
     & MPI_DOUBLE_PRECISION,expansion_terms_im,
     &	  n_r_coll,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)	 
      CALL MPI_GATHER(buffer_w,1,
     & MPI_DOUBLE_PRECISION,weig_terms,
     &	  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)	 
      IF(MYID.eq.0) THEN
! Bikram Oct'18 Start:
      bk_i=0
      inquire(file='PES_EXPAN_TERMS.DAT',exist=file_exst)
	  if(file_exst) then
	  open(1991,file='PES_EXPAN_TERMS.DAT',status='old')
	  do
      read(1991,*,iostat=bk_io)
	  if(bk_io/=0) then
	  exit
	  else
      bk_i=bk_i+1
	  endif
      enddo
	  endif
      close(1991)
! Bikram End.
      OPEN(1,FILE=EXP_FILE_NAME,POSITION = "APPEND",ACTION="WRITE")
      SELECT CASE(coll_type)
      CASE(1)
      DO i=1,N_TERMS_EXP
      l_t= A_TOP(1,i)
	  write(bk_header, '(1(a,i0))')"v_", l_t
! Bikram Oct'18 Start:
      if(bk_i.le.1) then
      IF(i.eq.1) THEN
      WRITE(1,'(a6,1x,a19,3x)',ADVANCE="NO")
     & "R_Bohr",trim(bk_header)
      ELSE
      WRITE(1,'(a19,3x)',ADVANCE="NO")
     & trim(bk_header)
      ENDIF	
	  endif  
      ENDDO
      if(bk_i.le.1) then
	  WRITE(1,*)
	  endif
! Bikram End.	
      DO i_r_point=i_nr_ini,i_nr_fin	  
      DO i=1,N_TERMS_EXP	  
      IF(i.eq.1) THEN	  
      WRITE(1,'(f6.2,1x,e19.12,3x)',ADVANCE="NO")
     & R_COM(i_r_point),expansion_terms(i_r_point,i)
      ELSE
      WRITE(1,'(e19.12,3x)',ADVANCE="NO")
     & expansion_terms(i_r_point,i)
      ENDIF	  
      ENDDO
      WRITE(1,*)	  
      ENDDO 
      CLOSE(1)
      CASE(3)
      DO i=1,N_TERMS_EXP
      l_t= A_TOP(1,i)
      nju_t = A_TOP(2,i)
	  write(bk_header, '(2(a,i0))')"v_", l_t, "_", nju_t
! Bikram Oct'18 Start:
      if(bk_i.le.1) then
      IF(i.eq.1) THEN
      WRITE(1,'(a6,1x,a19,3x)',ADVANCE="NO")
     & "R_Bohr",trim(bk_header)
      ELSE	  
      WRITE(1,'(a19,3x)',ADVANCE="NO")
     & trim(bk_header)
      ENDIF	
	  endif	  
      ENDDO
      if(bk_i.le.1) then
      WRITE(1,*)
      endif
! Bikram End.	  
      DO i_r_point=i_nr_ini,i_nr_fin	  
      DO i=1,N_TERMS_EXP	  
      IF(i.eq.1) THEN	  
      WRITE(1,'(f6.2,1x,e19.12,3x)',ADVANCE="NO")
     & R_COM(i_r_point),expansion_terms(i_r_point,i)
      ELSE
      WRITE(1,'(e19.12,3x)',ADVANCE="NO")
     & expansion_terms(i_r_point,i)
      ENDIF	  
      ENDDO
      WRITE(1,*)	  
      ENDDO 
      CLOSE(1)	  
      CASE(4)
      DO i=1,N_TERMS_EXP
      l_t= A_TOP(1,i)
      nju_t = A_TOP(2,i)
	  write(bk_header, '(2(a,i0))')"v_", l_t, "_", nju_t
! Bikram Oct'18 Start:
      if(bk_i.le.1) then
      IF(i.eq.1) THEN
      WRITE(1,'(a6,1x,a19,3x)',ADVANCE="NO")
     & "R_Bohr",trim(bk_header)
      ELSE	  
      WRITE(1,'(a19,3x)',ADVANCE="NO")
     & trim(bk_header)
      ENDIF
      endif	  
      ENDDO
      if(bk_i.le.1) then
      WRITE(1,*)
      endif
! Bikram End.	  
      DO i_r_point=i_nr_ini,i_nr_fin	  
      DO i=1,N_TERMS_EXP	  
      IF(i.eq.1) THEN	  
      WRITE(1,'(f6.2,1x,e19.12,3x)',ADVANCE="NO")
     & R_COM(i_r_point),expansion_terms(i_r_point,i)
      ELSE
      WRITE(1,'(e19.12,3x)',ADVANCE="NO")
     & expansion_terms(i_r_point,i)
      ENDIF	  
      ENDDO
      WRITE(1,*)	  
      ENDDO 
      CLOSE(1)	  	  
      CASE(5)
!! HEADER PRINTING	  
      DO i=1,N_TERMS_EXP!!!!!!!!!! i.eq.1 or i .eq 29
      l1_t= A_TOP(1,i)
      l2_t = A_TOP(2,i)
      l_t = A_TOP(3,i)	
	  write(bk_header, '(3(a,i0))')"v_", l1_t, "_", l2_t, "_", l_t
! Bikram Oct'18 Start:
      if(bk_i.le.1) then
      IF(i.eq.1) THEN
      WRITE(1,'(a6,1x,a19,3x)',ADVANCE="NO")
     & "R_Bohr",trim(bk_header)
      ELSE	  
      WRITE(1,'(a19,3x)',ADVANCE="NO")
     & trim(bk_header)
      ENDIF
      endif	  
      ENDDO
!! HEADER PRINTING	 
      if(bk_i.le.1) then
      WRITE(1,*)
      endif
! Bikram End.		  
      DO i_r_point=i_nr_ini,i_nr_fin	  
      DO i=1,N_TERMS_EXP	  
      IF(i.eq.1) THEN	  
      WRITE(1,'(f6.2,1x,e19.12,3x)',ADVANCE="NO")
     & R_COM(i_r_point),expansion_terms(i_r_point,i)
      ELSE
      WRITE(1,'(e19.12,3x)',ADVANCE="NO")
     & expansion_terms(i_r_point,i)
      ENDIF	  
      ENDDO
      WRITE(1,*)	  
      ENDDO
      CLOSE(1)
      CASE(7)
!!! HEADER PRINTING	  
      DO i=1,N_TERMS_EXP
      l1_t= A_TOP(1,i)
      nju1_t = A_TOP(2,i)
      l2_t = A_TOP(3,i)
      l_t = A_TOP(4,i)		
	  write(bk_header, '(4(a,i0))')"v_", l1_t, "_", nju1_t, 
     & "_", l2_t, "_", l_t
! Bikram Oct'18 Start:
      if(bk_i.le.1) then
      IF(i.eq.1) THEN
      WRITE(1,'(a6,1x,a19,3x)',ADVANCE="NO")
     & "R_Bohr",trim(bk_header)
      ELSE	  
      WRITE(1,'(a19,3x)',ADVANCE="NO")
     & trim(bk_header)
      ENDIF
      endif	  
      ENDDO
!! HEADER PRINTING	 
      if(bk_i.le.1) then
      WRITE(1,*)
      endif
! Bikram End.	  
      DO i_r_point=i_nr_ini,i_nr_fin  
      DO i=1,N_TERMS_EXP	  
      IF(i.eq.1) THEN	  
      WRITE(1,'(f6.2,1x,e19.12,3x)',ADVANCE="NO")
     & R_COM(i_r_point),expansion_terms(i_r_point,i)
      WRITE(*,'(f6.2,1x,e19.12,3x)',ADVANCE="NO")
     & R_COM(i_r_point),expansion_terms(i_r_point,i)
      ELSE
      WRITE(1,'(e19.12,3x)',ADVANCE="NO")
     & expansion_terms(i_r_point,i)
      WRITE(*,'(e19.12,3x)',ADVANCE="NO")
     & expansion_terms(i_r_point,i)
      ENDIF	  
      ENDDO
      WRITE(1,*)	  
      WRITE(*,*)	  
      ENDDO
      CLOSE(1)	   
      DO i_r_point=i_nr_ini,i_nr_fin  
      DO i=1,N_TERMS_EXP	  
      IF(i.eq.1) THEN	  
      WRITE(*,'(f6.2,1x,e19.12,3x)',ADVANCE="NO")
     & R_COM(i_r_point),expansion_terms_im(i_r_point,i)
      ELSE
      WRITE(*,'(e19.12,3x)',ADVANCE="NO")
     & expansion_terms_im(i_r_point,i)
      ENDIF	  
      ENDDO
      WRITE(*,*)	  
      ENDDO
      CASE(8)
!!! HEADER PRINTING	  
      DO i=1,N_TERMS_EXP
      l1_t= A_TOP(1,i)
      nju1_t = A_TOP(2,i)
      l2_t = A_TOP(3,i)
      l_t = A_TOP(4,i)	
	  write(bk_header, '(4(a,i0))')"v_", l1_t, "_", nju1_t, 
     & "_", l2_t, "_", l_t
! Bikram Oct'18 Start:
      if(bk_i.le.1) then
      IF(i.eq.1) THEN
      WRITE(1,'(a6,1x,a19,3x)',ADVANCE="NO")
     & "R_Bohr",trim(bk_header)
      ELSE	  
      WRITE(1,'(a19,3x)',ADVANCE="NO")
     & trim(bk_header)
      ENDIF
      endif  
      ENDDO
!! HEADER PRINTING	 
      if(bk_i.le.1) then
      WRITE(1,*)
      endif
! Bikram End.		  
      DO i_r_point=i_nr_ini,i_nr_fin  
      DO i=1,N_TERMS_EXP	  
      IF(i.eq.1) THEN	  
      WRITE(1,'(f6.2,1x,e19.12,3x)',ADVANCE="NO")
     & R_COM(i_r_point),expansion_terms(i_r_point,i)
      WRITE(*,'(f6.2,1x,e19.12,3x)',ADVANCE="NO")
     & R_COM(i_r_point),expansion_terms(i_r_point,i)
      ELSE
      WRITE(1,'(e19.12,3x)',ADVANCE="NO")
     & expansion_terms(i_r_point,i)
      WRITE(*,'(e19.12,3x)',ADVANCE="NO")
     & expansion_terms(i_r_point,i)
      ENDIF	  
      ENDDO
      WRITE(1,*)	  
      WRITE(*,*)	  
      ENDDO
      CLOSE(1)	  		  
      DO i_r_point=i_nr_ini,i_nr_fin  
      DO i=1,N_TERMS_EXP	  
      IF(i.eq.1) THEN	  
      WRITE(*,'(f6.2,1x,e19.12,3x)',ADVANCE="NO")
     & R_COM(i_r_point),expansion_terms_im(i_r_point,i)
      ELSE
      WRITE(*,'(e19.12,3x)',ADVANCE="NO")
     & expansion_terms_im(i_r_point,i)
      ENDIF	  
      ENDDO
      WRITE(*,*)	  
      ENDDO
      CASE(9)
!!! HEADER PRINTING	  
      DO i=1,N_TERMS_EXP
      l1_t= A_TOP(1,i)
      nju1_t = A_TOP(2,i)
      l2_t = A_TOP(3,i)
      nju2_t = A_TOP(4,i)
      l_t = A_TOP(5,i)	
	  write(bk_header, '(5(a,i0))')"v_", l1_t, "_", nju1_t, 
     & "_", l2_t, "_", nju2_t, "_", l_t
! Bikram Oct'18 Start:
      if(bk_i.le.1) then
      IF(i.eq.1) THEN
      WRITE(1,'(a6,1x,a19,3x)',ADVANCE="NO")
     & "R_Bohr",trim(bk_header)
      ELSE	  
      WRITE(1,'(a19,3x)',ADVANCE="NO")
     & trim(bk_header)
      ENDIF
      endif  
      ENDDO
!! HEADER PRINTING	 
      if(bk_i.le.1) then
      WRITE(1,*)
      endif
! Bikram End.	  
      DO i_r_point=i_nr_ini,i_nr_fin  
      DO i=1,N_TERMS_EXP	  
      IF(i.eq.1) THEN	  
      WRITE(1,'(f6.2,1x,e19.12,3x)',ADVANCE="NO")
     & R_COM(i_r_point),expansion_terms(i_r_point,i)!,
      ELSE
      WRITE(1,'(e19.12,3x)',ADVANCE="NO")
     & expansion_terms(i_r_point,i)
      ENDIF	  
      ENDDO
      WRITE(1,*)	  
      ENDDO
      CLOSE(1)	  
      CASE(0)
!!! HEADER PRINTING	  
! Bikram Start: Aug 14 2022, Changing header printing formatting and PES_EXPAN_TERMS.DAT file
      if(bk_i.le.1) then
      DO i=1,N_TERMS_EXP
      l1_t= A_TOP(1,i)
      nju1_t = A_TOP(2,i)
      l2_t = A_TOP(3,i)
      nju2_t = A_TOP(4,i)
      l_t = A_TOP(5,i)	
	  write(bk_header, '(5(a,i0))')"v_", l1_t, "_", nju1_t, 
     & "_", l2_t, "_", nju2_t, "_", l_t
! Bikram Oct'18 Start:
      IF(i.eq.1) THEN
      WRITE(1,'(a6,1x,a19,3x)',ADVANCE="NO")
     & "R_Bohr",trim(bk_header)  
      ELSE	  
      WRITE(1,'(a19,3x)',ADVANCE="NO")
     & trim(bk_header)  
      ENDIF
      ENDDO
      WRITE(1,*)
      endif
! Bikram End.		  
      DO i_r_point=i_nr_ini,i_nr_fin  
      DO i=1,N_TERMS_EXP	  
      IF(i.eq.1) THEN	  
      WRITE(1,'(f6.2,1x,e19.12,3x)',ADVANCE="NO")
     & R_COM(i_r_point),expansion_terms(i_r_point,i)!,
      ELSE
      WRITE(1,'(e19.12,3x)',ADVANCE="NO")
     & expansion_terms(i_r_point,i)
      ENDIF	  
      ENDDO
      WRITE(1,*)	  
      ENDDO
      CLOSE(1)
      END SELECT	  
      ENDIF
      IF(MYID.EQ.0) PRINT*,"TERMS HAVE BEEN COMPUTED" 	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
!      CALL MPI_FINALIZE (ierr_mpi)	  
      RETURN
	  
      END SUBROUTINE CALC_EXPANSION_TERMS

      SUBROUTINE READ_EXPANSION_TERMS
      USE VARIABLES
      USE GRID_INTEGR
      USE CONSTANTS
      USE MPI_DATA
      USE MPI	  
      IMPLICIT NONE
      LOGICAL file_exst
      CHARACTER(LEN=19)  :: TERMS_DAT = "PES_EXPAN_TERMS.DAT"	  
      INTEGER i,i_r_point,istat,ich1,ich2
      IF(.not.ALLOCATED(expansion_terms))
     & ALLOCATE(expansion_terms(n_r_coll,nterms))
      IF(.not.ALLOCATED(expansion_terms_der))
     & ALLOCATE(expansion_terms_der(n_r_coll,nterms))	 
	 
       IF(.not.ALLOCATED(A_TOP))
     & ALLOCATE(A_TOP(5,nterms))
       IF(.not.ALLOCATED(expansion_terms_im))
     & ALLOCATE(expansion_terms_im(n_r_coll,nterms))
      IF(MYID.EQ.0) PRINT*, "TERMS ARE BEING READ"	  
      DO i=1,nterms
      SELECT CASE(coll_type)
      CASE(1)
      A_TOP(1,i) = L_TERM(i)
      IF(fine_structure_defined .and. LORB_FINE.gt.0)
     & A_TOP(2,i) = M_TERM(i)  
      CASE(2)
      A_TOP(1,i) = L_TERM(i)
      CASE(3)
      A_TOP(1,i) = L_TERM(i)
      A_TOP(2,i) = M_TERM(i)	 
      CASE(4)
      A_TOP(1,i) = L_TERM(i)
      A_TOP(2,i) = M_TERM(i)	  
      CASE(5)
      A_TOP(1,i) = L1_TERM(i) 
      A_TOP(2,i) = L2_TERM(i)
      A_TOP(3,i) = L_TERM(i)
      CASE(6)
      A_TOP(1,i) = L1_TERM(i) 
      A_TOP(2,i) = L2_TERM(i)
      A_TOP(3,i) = L_TERM(i)	  
      CASE(7)
      A_TOP(1,i) = L1_TERM(i) 
      A_TOP(2,i)  = NJU1_TERM(i) 
      A_TOP(3,i) = L2_TERM(i)
      A_TOP(4,i) = L_TERM(i)	  
      CASE(8)
      A_TOP(1,i) = L1_TERM(i) 
      A_TOP(2,i)  = NJU1_TERM(i) 
      A_TOP(3,i) = L2_TERM(i)
      A_TOP(4,i) = L_TERM(i)	  
      CASE(9)
      A_TOP(1,i) = L1_TERM(i) 
      A_TOP(2,i)  = NJU1_TERM(i) 
      A_TOP(3,i) = L2_TERM(i)
      A_TOP(4,i) = NJU2_TERM(i)
      A_TOP(5,i) = L_TERM(i)	  
      CASE(0)	  
      A_TOP(1,i) = L1_TERM(i) 
      A_TOP(2,i)  = NJU1_TERM(i) 
      A_TOP(3,i) = L2_TERM(i)
      A_TOP(4,i) = NJU2_TERM(i)
      A_TOP(5,i) = L_TERM(i)
      END SELECT	  
      ENDDO

      IF(pot_expansion_defined) THEN
      DO i_r_point=1,n_r_coll	  
      DO i=1,nterms
      CALL USER_DEFINED_TERMS(expansion_terms(i_r_point,i),
     & i,R_COM(i_r_point))
      ENDDO
      ENDDO
      TERM_READ_STAT = .TRUE.
      RETURN 	  
      ENDIF
      INQUIRE( FILE=TERMS_DAT, EXIST=file_exst)
      IF(.not.file_exst) THEN
      IF(MYID.eq.0) PRINT*,TERMS_DAT," NOT FOUND"
      CALL MPI_FINALIZE(ierr_mpi)
      STOP	  
      ENDIF	  
      IF(MYID.eq.0) THEN
	  
      OPEN(1,FILE=TERMS_DAT,ACTION="READ",STATUS="OLD")	  
      READ(1,*,IOSTAT=istat) ! HEADER
      DO i_r_point=1,n_r_coll	  
      DO i=1,nterms
      IF(i.eq.1) THEN	  
      READ(1,'(f6.2,1x,e19.12,3x)',ADVANCE="NO",IOSTAT=istat)
     & R_COM(i_r_point),expansion_terms(i_r_point,i)
      ELSE
      READ(1,'(e19.12,3x)',ADVANCE="NO")
     & expansion_terms(i_r_point,i)!,
      ENDIF	  
      ENDDO
      READ(1,*,IOSTAT=istat) 
      ENDDO
      CLOSE(1)
      TERM_READ_STAT = .TRUE.	  
      IF(istat.gt.0) THEN
      PRINT*,"ERROR IN READING TERMS.DAT"
      TERM_READ_STAT = .FALSE.	  
      ENDIF	 	  
      ENDIF
      CALL MPI_BCAST(TERM_READ_STAT, 1, MPI_LOGICAL,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CALL MPI_BCAST(R_COM, n_r_coll, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(expansion_terms, n_r_coll*nterms
     & , MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	 
      IF(.not.TERM_READ_STAT) THEN
      STOP
      ENDIF
      R_COM = R_COM/conv_unit_r	  
      IF(MYID.EQ.0) PRINT*, "TERMS HAVE BEEN READ"
      RETURN	  
      END SUBROUTINE READ_EXPANSION_TERMS
      !!! ASSYMETRIC + ASSYMETRIX TOP ROUTINE TO EXTRACT EXPAJSION   
	  
      SUBROUTINE READ_USER_TERMS
      USE VARIABLES
      USE GRID_INTEGR
      USE CONSTANTS
      USE MPI_DATA
      IMPLICIT NONE
      INTEGER i,i_r_point
      ALLOCATE(weig_terms(nterms))		  
      IF(.not.ALLOCATED(expansion_terms))
     & ALLOCATE(expansion_terms(n_r_coll,nterms))
      IF(.not.ALLOCATED(A_TOP))  
     &  ALLOCATE(A_TOP(5,nterms))
      IF(.not.ALLOCATED(expansion_terms_im))
     & ALLOCATE(expansion_terms_im(n_r_coll,nterms))
      IF(MYID.EQ.0) PRINT*,
     & " USER DEFINED EXPANSION TERMS WILL BE COMPUTED"	  
      DO i=1,nterms
      SELECT CASE(coll_type)
      CASE(1)
      A_TOP(1,i) = L_TERM(i)
      IF(fine_structure_defined .and. LORB_FINE.gt.0)
     & A_TOP(2,i) = M_TERM(i)  	  
      CASE(2)
      A_TOP(1,i) = L_TERM(i)
      CASE(3)
      A_TOP(1,i) = L_TERM(i)
      A_TOP(2,i) = M_TERM(i)	 
      CASE(4)
      A_TOP(1,i) = L_TERM(i)
      A_TOP(2,i) = M_TERM(i)	  
      CASE(5)
      A_TOP(1,i) = L1_TERM(i) 
      A_TOP(2,i) = L2_TERM(i)
      A_TOP(3,i) = L_TERM(i)
      CASE(6)
      A_TOP(1,i) = L1_TERM(i) 
      A_TOP(2,i) = L2_TERM(i)
      A_TOP(3,i) = L_TERM(i)	  
      CASE(7)
      A_TOP(1,i) = L1_TERM(i) 
      A_TOP(2,i)  = NJU1_TERM(i) 
      A_TOP(3,i) = L2_TERM(i)
      A_TOP(4,i) = L_TERM(i)	  
      CASE(8)
      A_TOP(1,i) = L1_TERM(i) 
      A_TOP(2,i)  = NJU1_TERM(i) 
      A_TOP(3,i) = L2_TERM(i)
      A_TOP(4,i) = L_TERM(i)	  
      CASE(9)
      A_TOP(1,i) = L1_TERM(i) 
      A_TOP(2,i)  = NJU1_TERM(i) 
      A_TOP(3,i) = L2_TERM(i)
      A_TOP(4,i) = NJU2_TERM(i)
      A_TOP(5,i) = L_TERM(i)	  
      CASE(0)	  
      A_TOP(1,i) = L1_TERM(i) 
      A_TOP(2,i)  = NJU1_TERM(i) 
      A_TOP(3,i) = L2_TERM(i)
      A_TOP(4,i) = NJU2_TERM(i)
      A_TOP(5,i) = L_TERM(i)
      END SELECT	  
      ENDDO
  
      RETURN	  
	  
      END SUBROUTINE READ_USER_TERMS	  

      SUBROUTINE EXPANSION_MATRIX_ELEMENT(M_coulp,k,i_r_point,tmp1)
      USE VARIABLES
      USE FINE_STRUCT_LEVELS	  
      USE EIGEN_VECTORS	  
      USE GRID_INTEGR
      USE EXPANSION_MAT_STORAGE
      USE MPI_DATA
      USE MPI
      USE COMPUTE_EXPANSION_VARIABLES	  
      IMPLICIT NONE
      LOGICAL TRIANG_RULE
      INTEGER KRONEKER,i,k,i_r_point,ident_coeff,round
      INTEGER :: i_nr_ini   
      REAL*8 CG,delta,W3JS,	M_coulp,CG_HF,sign_v  
	  real*8 bk_xn,bk_z,bk_dq1
	  real*8 bk_cg1, bk_cg2, CG_bikram
	  real*8 coef1, coef2, coef3, coef4, CG_out, CG_in_all
	  real*8 CG_in1a, CG_in1b, CG_in2a, CG_in2b
	  real*8 CG_in3a, CG_in3b, CG_in4a, CG_in4b
	  integer bk_str_k, tmp1
!	  real*8 tmp_bk_mat(4), tmp_bk_mat1(4)
	  integer bk_par1, bk_par2, sign_correction
!	  character (len = 1) psign1, psign2, psign3, psign4
      EXTERNAL CG,delta,W3JS,TRIANG_RULE,KRONEKER,CG_HF,round,sign_v
	  external CG_bikram
	  logical same_terms
      term_limit_user = nterms	  
      M_coulp = 0d0
      buff= 0d0 	
      i_nr_ini = max(ir_bgn_exp_pnt,1)		  
      stp = ind_mat(1,k)
      stpp = ind_mat(2,k)	  
      SIMPLIFICATION_EXP_MAT = .FALSE.
	  if(bikram_rebalance) tmp1 = 0	  

      IF(.not.TERM_READ_STAT) THEN
      IF(MYID.eq.0 .and. fine_structure_defined) PRINT*, b_fine_coeff !!! DELETE       	  
      IF(MYID.EQ.0) PRINT *, "USER TERM FILE IS BEING READ "	  
      CALL READ_EXPANSION_TERMS!! DELETE AFTER ALLL	  
      ENDIF 	  
      IF(.not.EXPANSION_GRID_DEFINED) THEN

!!      CALL READ_COCO_TERMS	  
      ALLOCATE(TERM_MATRIX_ELEMENT(2,2,4,nterms))
      TERM_MATRIX_ELEMENT = 0d0  
      EXPANSION_GRID_DEFINED = .TRUE.
      IF(MYID.EQ.0) PRINT *, "TERM MATRIX CREATED "  
      ENDIF
      SELECT CASE(coll_type)
      CASE(1)
      IF(.not.fine_structure_defined) THEN	  
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)	  
      m_t = m12(stp)
      IF(m_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"
 	  
      j_p_t = j_ch(channp)
      j_pp_t = j_ch(channpp)
      IF(j_pp_t.ne.j12(stpp)) STOP "ERROR: INI IN BASIS"
      IF(j_p_t.ne.j12(stp)) STOP "ERROR: INI IN BASIS"
	  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for any molecule-molecule collision using expansion
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022
	  if(.not.bikram_identical_pes) then
	  
	  coef1 = dsqrt((2d0*j_p_t + 1d0)/(2d0*j_pp_t + 1d0))
	  
      DO i=1,nterms
      ind_t =i
	  exp_coeff_int = expansion_terms(i_r_point,ind_t)
      l_t= A_TOP(1,ind_t)
      IF(.NOT.TRIANG_RULE(j_pp_t,l_t,j_p_t)) CYCLE
	  coef2 = dsqrt((2d0*l_t + 1d0) / (4d0*pi))
	  CG_j1_l1 = CG_bikram(j_p_t,l_t,j_pp_t,m_t,0,m_t,bikram_w3j_fact)
	  CG_j2_l2 = CG_bikram(j_p_t,l_t,j_pp_t,0,0,0,bikram_w3j_fact)
	  
      M_coulp = M_coulp + coef1 * exp_coeff_int* coef2 *
     & CG_j1_l1 * CG_j2_l2
      ENDDO
	  else
	  print*,"Something is wrong in PES expansion matrix computation."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
! Bikram Start: Feb 2022
! This part is to use recursive method of computing Wigner 3j symbol
	  if(M_coulp.ne.M_coulp) then
	  if(myid.eq.0) then
	  print*, "Computed matrix element is wrong for state:"
      write(*,'(6(a,i0),2x,e19.12)')'ini: ',j_p_t,',',k_ch(channp),
     & ',',eps_p_t, '  fin: ', j_pp_t,',',k_ch(channpp),',',
     & eps_pp_t, M_coulp
	  print*, "Program will stop now."
	  if(bikram_w3j_fact) print*, "Recommended to use recursive
     & algorithm for computing Wigner 3j coefficients."
	  stop
      CALL MPI_FINALIZE (ierr_mpi)
	  end if
	  stop
      CALL MPI_FINALIZE (ierr_mpi)
	  end if
! Bikram End.
	  
      ELSE
      SELECT CASE(LORB_FINE)
      CASE(0)	  
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)
	  m_t = m12(stp)
      IF(m_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"
      j_p_t = j_ch(channp)
      j_pp_t = j_ch(channpp)
      IF(j_pp_t.ne.j12(stpp)) STOP "ERROR: INI IN BASIS"
      IF(j_p_t.ne.j12(stp)) STOP "ERROR: INI IN BASIS"	
      DO i=1,nterms 
      ind_t =i! index_term_in_file(i)	  
      l_t= A_TOP(1,ind_t)
      IF(.NOT.TRIANG_RULE(j_pp_t,l_t,j_p_t)) CYCLE
      
      IF(i_r_point.eq.i_nr_ini) THEN
      TERM_MATRIX_ELEMENT(1,1,1,ind_t) = 0d0	  
      DO f_p_t = -1,1 ! 1==S_FINE	!!! CHANGE IT
      b_p_t = b_fine_coeff(f_p_t+2,channp)
      n_p_t = j_p_t + f_p_t 
  
      DO f_pp_t = -1,1 ! 1==S_FINE	
      b_pp_t = b_fine_coeff(f_pp_t+2,channpp)
      n_pp_t = j_pp_t + f_pp_t 	  
      DO s_p_t = -1,1 !!! S_FINE
      m_s_t = m_t-s_p_t		  
      IF(abs(m_s_t).gt.n_p_t) CYCLE
      IF(abs(m_s_t).gt.n_pp_t) CYCLE
      CG_p_s = CG(n_p_t,1,j_p_t,m_s_t,s_p_t,m_t) ! 1==S_FINE
      CG_pp_s = CG(n_pp_t,1,j_pp_t,m_s_t,s_p_t,m_t) ! 1==S_FINE		  
      TERM_MATRIX_ELEMENT(1,1,1,ind_t) = 
     & TERM_MATRIX_ELEMENT(1,1,1,ind_t) + 
     & 	b_p_t*b_pp_t*CG_p_s*CG_pp_s*
     &  CG(n_p_t,l_t,n_pp_t,m_s_t,0,m_s_t)*
     & CG(n_p_t,l_t,n_pp_t,0,0,0)*	 
     & dsqrt((2d0*n_p_t+1)/(2d0*n_pp_t+1d0))
      IF(k.eq.-14) THEN !!!! CHECK LATER
      WRITE(*,'(i2,1x,i2,1x,i2,1x,i2,1x,f7.4)')
     & f_p_t ,f_pp_t,s_p_t,m_s_t,TERM_MATRIX_ELEMENT(1,1,1,ind_t)
	  WRITE(*,'(f6.3,1x,f6.3,1x,f6.3,1x,f6.3,1x,f6.3)')
     &	 b_p_t,b_pp_t,CG_p_s,CG_pp_s,
     & CG(n_p_t,l_t,n_pp_t,m_s_t,0,m_s_t)*
     & CG(n_p_t,l_t,n_pp_t,0,0,0)*	 
     & dsqrt((2d0*n_p_t+1)/(2d0*n_pp_t+1d0))
      ENDIF	  
      ENDDO
      ENDDO
      ENDDO	  
      	  
      ENDIF		  
      M_coulp = M_coulp + 
     & expansion_terms(i_r_point,ind_t)*
     & TERM_MATRIX_ELEMENT(1,1,1,ind_t)
      ENDDO
      CASE(1)
       channp = indx_chann(stp)
      channpp = indx_chann(stpp)
      eps_p_t = (-1)**par_lorb_ch(channp) !! eps = 0, +, eps = 1 -
	  eps_pp_t = (-1)**par_lorb_ch(channpp)	  
	  m_h_t = m12_h(stp) !!! change m12_h !!par_lorb_ch(:)
      IF(int(m_h_t).ne.int(m12_h(stpp))
     & .or. m_h_t*m12_h(stpp).lt.0 ) STOP "ERROR: INI IN BASIS"
      j_h_p_t = j_h_ch(channp)
      j_h_pp_t = j_h_ch(channpp)
      IF(int(j_h_pp_t).ne.int(j12_h(stpp))) STOP "ERROR: INI IN BASIS"
      IF(int(j_h_p_t).ne.int(j12_h(stp))) STOP "ERROR: INI IN BASIS"
      DO i=1,nterms 
	  ind_t =i! index_term_in_file(i)	  
      l_t= A_TOP(1,ind_t)
!!! must be two or zero
!      IF(.NOT.TRIANG_RULE(j_pp_t,l_t,j_p_t)) CYCLE
!      IF(2*max(j_h_pp_t,l_t,j_h_p_t).gt.
!     & j_h_pp_t+l_t+j_h_p_t) CYCLE
      IF(i_r_point.eq.i_nr_ini) THEN
      TERM_MATRIX_ELEMENT(1,1,1,ind_t) = 0d0
      DO f_p_t = 1,2
      b_p_t = b_fine_coeff(f_p_t,channp) !! f_p_t =1, PI 1/2, f_pp_t=2 = PI=3/2
!      w_h_p =  (f_p_t-0.5d0)!*(-1)**eps_p_t	!!! minus  
      DO f_pp_t = 1,2
      b_pp_t = b_fine_coeff(f_pp_t,channpp)
!      w_h_pp =  (f_pp_t-0.5d0)!*(-1)**eps_pp_t !!!  	plus 
      DO sp_p=1,2
      DO sp_pp=1,2
      IF(sp_pp.ne.sp_p) CYCLE
      sp_p_t = sp_p-1.5d0
      sp_pp_t = sp_pp-1.5d0
      w_h_p = (2*(f_p_t-1) - 0.5d0)*(-1)**sp_p
      w_h_pp = (2*(f_pp_t-1) - 0.5d0)*(-1)**sp_pp
      lambda_p_h = round(w_h_p-sp_p_t)
      lambda_pp_h = round(w_h_pp-sp_pp_t)
      IF(abs(lambda_p_h).ne.1) STOP "ERROR IN 1/2PI FINE"
      IF(abs(lambda_pp_h).ne.1) STOP "ERROR IN 1/2PI FINE" 
	  coeff_1_p = 1
      coeff_1_pp = 1	  
      IF(lambda_p_h.eq.-1) coeff_1_p =  eps_p_t
      IF(lambda_pp_h.eq.-1) coeff_1_pp =  eps_pp_t
      nju_t = A_TOP(2,ind_t)*sign_v(lambda_pp_h-lambda_p_h)	  
      TERM_MATRIX_ELEMENT(1,1,1,ind_t) =
     & TERM_MATRIX_ELEMENT(1,1,1,ind_t) + 
     & sqrt((2*j_h_p_t+1)/(2*j_h_pp_t+1))
     & *CG_HF(j_h_p_t,l_t,j_h_pp_t,m_h_t,0,m_h_t)
     & *CG_HF(j_h_p_t,l_t,j_h_pp_t,w_h_p,nju_t,w_h_pp)/2
     & *coeff_1_p*coeff_1_pp*b_p_t*b_pp_t!! COULD BE SOMETHING HERE!!! to be continued
	!!! TO BE CONTINUED 
      IF(k.eq.17 .and. ind_t.eq.3) THEN
      PRINT*,f_p_t,f_pp_t,sp_p_t
      PRINT*,CG_HF(j_h_p_t,l_t,j_h_pp_t,m_h_t,0,m_h_t),
     & CG_HF(j_h_p_t,l_t,j_h_pp_t,w_h_p,nju_t,w_h_pp)
      PRINT*,coeff_1_p*coeff_1_pp,
     & b_p_t,b_pp_t
      PRINT*,w_h_p,w_h_pp,j_h_p_t+l_t-j_h_pp_t
      PRINT*,
     & sqrt((2*j_h_p_t+1)/(2*j_h_pp_t+1))*
     & CG_HF(j_h_p_t,l_t,j_h_pp_t,m_h_t,0,m_h_t)
     & *CG_HF(j_h_p_t,l_t,j_h_pp_t,w_h_p,nju_t,w_h_pp)/2
     & *coeff_1_p*coeff_1_pp*b_p_t*b_pp_t	  
      PRINT*,"trt"	  
      ENDIF	  
      ENDDO
      ENDDO	  
      ENDDO	 	  
      ENDDO	  
      ENDIF	 
      M_coulp = M_coulp + 
     & expansion_terms(i_r_point,ind_t)*
     & TERM_MATRIX_ELEMENT(1,1,1,ind_t)

	  
      ENDDO
      END SELECT	  
      ENDIF      
      CASE(2)
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)	  
      m_t = m12(stp)
      IF(m_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"
 	  
c      PRINT*, j12_p_t,j12_pp_t,m12_t	  
      j_p_t = j_ch(channp)
      j_pp_t = j_ch(channpp)
      IF(j_pp_t.ne.j12(stpp)) STOP "ERROR: INI IN BASIS"
      IF(j_p_t.ne.j12(stp)) STOP "ERROR: INI IN BASIS"	  
      DO i=1,nterms
      ind_t =i! index_term_in_file(i)	  
      l_t= A_TOP(1,ind_t)
      IF(.NOT.TRIANG_RULE(j_pp_t,l_t,j_p_t)) CYCLE	  
      IF(i_r_point.eq.i_nr_ini) THEN
      TERM_MATRIX_ELEMENT(1,1,1,i) =
     & CG(j_p_t,l_t,j_pp_t,m_t,0,m_t)*
     & CG(j_p_t,l_t,j_pp_t,0,0,0)*
     & dsqrt((2d0*j_p_t+1)/(2d0*j_pp_t+1d0))  	  
      ENDIF
      M_coulp = M_coulp + 
     & expansion_terms(i_r_point,ind_t)*
     & TERM_MATRIX_ELEMENT(1,1,1,ind_t)
      ENDDO
      M_coulp = M_coulp * vib_overlap(channp,channpp)	  
      CASE(3)
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)	  
      m_t = m12(stp)
      IF(m_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"
 	  
! specifying local variables to store quantum numbers of basis
      j_p_t = j_ch(channp)
      j_pp_t = j_ch(channpp)
      k_p_t = k_ch(channp)
      k_pp_t = k_ch(channpp)	  
      eps_p_t = eps_ch(channp)	  
      eps_pp_t = eps_ch(channpp)	  
      IF(j_pp_t.ne.j12(stpp)) STOP "ERROR: INI IN BASIS"
      IF(j_p_t.ne.j12(stp)) STOP "ERROR: INI IN BASIS"
	  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for any molecule-molecule collision using expansion
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022
	  if(.not.bikram_identical_pes) then

	  same_terms = .false.
	  if(i_r_point.eq.1 .and. k.eq.1) then
      DO i = 1,nterms
      ind_t = i
      l_t= A_TOP(1,ind_t)
      nju_t = A_TOP(2,ind_t)
	  
	  if(nju1_t.lt.0) then
	  write(*,'(a, a, 5(i0,a))')"The values of M should be ",
     & "always positive. Please check the expansion term: ",
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_terms = .true.
	  end if
	  end do
	  end if
	  if(same_terms) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if

	  coef1 = dsqrt((2d0*j_p_t + 1d0)/(2d0*j_pp_t + 1d0))

      DO i=1,nterms
      ind_t =i
      exp_coeff_int = expansion_terms(i_r_point,ind_t)
      l_t= A_TOP(1,ind_t)
      nju_t = A_TOP(2,ind_t)
	  
	  coef2 = 1.d0 / dsqrt(2d0*(1d0 + KRONEKER(k_p_t,0)) * 
     &                     2d0*(1d0 + KRONEKER(k_pp_t,0)))
	  coef3 = dsqrt((2d0*l_t + 1d0) / (4d0*pi))
	  CG_j1_l1 = CG_bikram(j_p_t,l_t,j_pp_t,m_t,0,m_t,bikram_w3j_fact)

      DO planr_coeff =1,2-KRONEKER(nju_t,0)
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int = exp_coeff_int*(-1)**nju_t 
      nju_t = -nju_t
      ENDIF
	  
      DO par_p  = 1, 2
      DO par_pp = 1, 2
      k_p_t = (-1)**(par_p-1)*k_ch(channp)
      k_pp_t = (-1)**(par_pp-1)*k_ch(channpp)
      IF(.NOT.TRIANG_RULE(j_pp_t,l_t,j_p_t)) CYCLE
      IF(k_pp_t.ne.k_p_t+nju_t) CYCLE
	  CG_j1_k1 = CG_bikram(j_p_t,l_t,j_pp_t,k_p_t,+nju_t,k_pp_t,
     & bikram_w3j_fact)
	  
      M_coulp = M_coulp + coef1 * coef2 * exp_coeff_int * coef3 * 
     & CG_j1_l1 * (-1)**nju_t * CG_j1_k1 * 
     & ((-1)**eps_p_t)**(par_p-1)*((-1)**eps_pp_t)**(par_pp-1)
      ENDDO
      ENDDO
      ENDDO	  
      ENDDO  	 	  
	  
	  else
	  print*,"Something is wrong in PES expansion matrix computation."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
! Bikram Start: Feb 2022
! This part is to use recursive method of computing Wigner 3j symbol
	  if(M_coulp.ne.M_coulp) then
	  if(myid.eq.0) then
	  print*, "Computed matrix element is wrong for state:"
      write(*,'(6(a,i0),2x,e19.12)')'ini: ',j_p_t,',',k_ch(channp),
     & ',',eps_p_t, '  fin: ', j_pp_t,',',k_ch(channpp),',',
     & eps_pp_t, M_coulp
	  print*, "Program will stop now."
	  if(bikram_w3j_fact) print*, "Recommended to use recursive
     & algorithm for computing Wigner 3j coefficients."
	  stop
      CALL MPI_FINALIZE (ierr_mpi)
	  end if
	  stop
      CALL MPI_FINALIZE (ierr_mpi)
	  end if
! Bikram End.

      CASE(4)
      CALL ASYM_TOP_VECTORS	  
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)	  
      m_t = m12(stp)
      IF(m_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"
 	  
c      PRINT*, j12_p_t,j12_pp_t,m12_t	  
      j_p_t = j_ch(channp)
      j_pp_t = j_ch(channpp)
      IF(j_pp_t.ne.j12(stpp)) STOP "ERROR: INI IN BASIS"
      IF(j_p_t.ne.j12(stp)) STOP "ERROR: INI IN BASIS"
	  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for any molecule-molecule collision using expansion
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022
	  
	  if(.not.bikram_identical_pes) then

	  same_terms = .false.
	  if(i_r_point.eq.1 .and. k.eq.1) then
      DO i = 1,nterms
      ind_t = i
      l_t= A_TOP(1,ind_t)
      nju_t = A_TOP(2,ind_t)
	  
	  if(nju1_t.lt.0) then
	  write(*,'(a, a, 5(i0,a))')"The values of M should be ",
     & "always positive. Please check the expansion term: ",
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_terms = .true.
	  end if
	  end do
	  end if
	  if(same_terms) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
	  coef1 = dsqrt((2d0*j_p_t + 1d0)/(2d0*j_pp_t + 1d0))

      DO i=1,nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,ind_t)
      l_t= A_TOP(1,ind_t)
      nju_t = A_TOP(2,ind_t)
	  
	  coef2 = dsqrt((2d0*l_t + 1d0) / (4d0*pi))
	  CG_j1_l1 = CG_bikram(j_p_t,l_t,j_pp_t,m_t,0,m_t,bikram_w3j_fact)
	  
      DO planr_coeff =1, 2-KRONEKER(nju_t,0)
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int = exp_coeff_int*(-1)**nju_t 
      nju_t = -nju_t
      ENDIF
	  
      DO  k_p_t = -j_p_t,j_p_t
      k_pp_t = k_p_t + nju_t
      IF(abs(k_pp_t).gt.j_pp_t) CYCLE		  
      IF(.NOT.TRIANG_RULE(j_pp_t,l_t,j_p_t)) CYCLE
      coeff_p = M_VECTORS(channp,j_p_t+1+k_p_t)
      coeff_pp = M_VECTORS(channpp,j_pp_t+1+k_pp_t)
	  CG_j1_k1 = CG_bikram(j_p_t,l_t,j_pp_t,k_p_t,nju_t,k_pp_t,
     & bikram_w3j_fact)
 
      M_coulp = M_coulp + coef1 * exp_coeff_int * coef2 *(-1)**nju_t
     & * CG_j1_l1 * CG_j1_k1 * coeff_p * coeff_pp	 							!!!!!! ADD EPSILON  	  
      ENDDO
      ENDDO	  
      ENDDO

	  else
	  print*,"Something is wrong in PES expansion matrix computation."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
! Bikram Start: Feb 2022
! This part is to use recursive method of computing Wigner 3j symbol
	  if(M_coulp.ne.M_coulp) then
	  if(myid.eq.0) then
	  print*, "Computed matrix element is wrong for state:"
      write(*,'(6(a,i0),2x,e19.12)')'ini: ',j_p_t,',',k_ch(channp),
     & ',',eps_p_t, '  fin: ', j_pp_t,',',k_ch(channpp),',',
     & eps_pp_t, M_coulp
	  print*, "Program will stop now."
	  if(bikram_w3j_fact) print*, "Recommended to use recursive
     & algorithm for computing Wigner 3j coefficients."
	  stop
      CALL MPI_FINALIZE (ierr_mpi)
	  end if
	  stop
      CALL MPI_FINALIZE (ierr_mpi)
	  end if
! Bikram End.
	  
      CASE(5)
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)	  

      j12_p_t = j12(stp)	  
      m12_t = m12(stp)
      j12_pp_t = j12(stpp)
      IF(m12_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"   	  

      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for any molecule-molecule collision using expansion
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022
	  

	  if(.not.bikram_identical_pes) then
	  
	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
      DO i=1,nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,i)
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(4d0*pi)
	  
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
	  CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,0,0,0,bikram_w3j_fact)
	  CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,0,0,0,bikram_w3j_fact)
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
	  CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)
	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
	  CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
	  CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
	  CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
	  CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)

	  M_coulp = M_coulp + coef1*coef2*exp_coeff_int*coef3* 
     & CG_j1_k1 * CG_j2_k2 * CG_j1_j2_p * CG_j1_j2_pp *
     & CG_l1_l2 * CG_j1_l1 * CG_j2_l2
      ENDDO	  
      ENDDO  
      ENDDO

	  else if(bikram_identical_pes) then
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for identical molecule-molecule collision using expansion
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022

	  same_terms = .false.
	  if(i_r_point.eq.1 .and. k.eq.1) then
      DO i = 1, nterms
      ind_t = i
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
	  
	  if(l1_t.lt.l2_t) then
	  write(*,'(a, a, a, 5(i0,a))')"The values of L2 should be ",
     & "less or equal to L1. ", 
     & "Please check the expansion term: ", 
     & l1_t, ',', l2_t, ',', l_t
	  same_terms = .true.
	  end if
	  end do
	  end if
	  if(same_terms) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if	  

	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
      DO i=1,nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,i)
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(4d0*pi)
	  sign_correction = 1 
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
	  CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)
	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
	  CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
	  
	  DO symmetry_coeff = 1, 2 - KRONEKER(l1_t,l2_t)
      if(symmetry_coeff.eq.1) then
	  sign_correction = 1
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      else if(symmetry_coeff.eq.2) then
      sign_correction = (-1)**(l1_t+l2_t)
      l1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(1,ind_t)
      l_t = A_TOP(3,ind_t)
	  else
	  if(myid.eq.0 .and. i_r_point.eq.1) 
     & print*, "Error in identical swapping variable."
	  stop
      end if
	  IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
	  
	  CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,0,0,0,bikram_w3j_fact)
	  CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,0,0,0,bikram_w3j_fact)
	  CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
	  CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
	  CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)

	  M_coulp = M_coulp + coef1*coef2*exp_coeff_int*coef3* 
     & CG_j1_k1 * CG_j2_k2 * CG_j1_j2_p * CG_j1_j2_pp *
     & CG_l1_l2 * CG_j1_l1 * CG_j2_l2 * sign_correction
      ENDDO	  
      ENDDO  
      ENDDO
      ENDDO
	  
	  else
	  print*,"Something is wrong in PES expansion matrix computation."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if

	  
! Bikram Start: Feb 2022
! This part is to use recursive method of computing Wigner 3j symbol
	  if(M_coulp.ne.M_coulp) then
	  if(myid.eq.0) then
	  print*, "Computed matrix element is wrong for state:"
      write(*,'(6(a,i0),2x,e19.12)')'ini: ',j_p_t,',',k_ch(channp),
     & ',',eps_p_t, '  fin: ', j_pp_t,',',k_ch(channpp),',',
     & eps_pp_t, M_coulp
	  print*, "Program will stop now."
	  if(bikram_w3j_fact) print*, "Recommended to use recursive
     & algorithm for computing Wigner 3j coefficients."
	  stop
      CALL MPI_FINALIZE (ierr_mpi)
	  end if
	  stop
      CALL MPI_FINALIZE (ierr_mpi)
	  end if
! Bikram End.

! This piece is for identical molecules case with both inversion and exchange parity
	  M_coulp_non_ident = M_coulp
      IF(identical_particles_defined) THEN
      par_p = parity_state(stp)
      par_pp = parity_state(stpp)
      j12_p_t = j12(stp)	  
      j12_pp_t = j12(stpp)
      M_coulp_ident = 0d0	  
      j1_p_t = j2_ch(channp)
      j2_p_t = j1_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)	  

	  if(.not.bikram_identical_pes) then
	  
	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
      DO i=1,nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,i)
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(4d0*pi)
	  
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
	  CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,0,0,0,bikram_w3j_fact)
	  CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,0,0,0,bikram_w3j_fact)
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
	  CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)
	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
	  CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
	  CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
	  CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
	  CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)

	  M_coulp_ident = M_coulp_ident + coef1*coef2*exp_coeff_int*coef3* 
     & CG_j1_k1 * CG_j2_k2 * CG_j1_j2_p * CG_j1_j2_pp *
     & CG_l1_l2 * CG_j1_l1 * CG_j2_l2
      ENDDO	  
      ENDDO  
      ENDDO

	  else if(bikram_identical_pes) then
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for identical molecule-molecule collision using expansion
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022

	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
      DO i=1,nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,i)
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(4d0*pi)
	  sign_correction = 1 
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
	  CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)
	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
	  CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
	  
	  DO symmetry_coeff = 1, 2 - KRONEKER(l1_t,l2_t)
      if(symmetry_coeff.eq.1) then
	  sign_correction = 1
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      else if(symmetry_coeff.eq.2) then
      sign_correction = (-1)**(l1_t+l2_t)
      l1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(1,ind_t)
      l_t = A_TOP(3,ind_t)
	  else
	  if(myid.eq.0 .and. i_r_point.eq.1) 
     & print*, "Error in identical swapping variable."
	  stop
      end if
	  IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
	  
	  CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,0,0,0,bikram_w3j_fact)
	  CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,0,0,0,bikram_w3j_fact)
	  CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
	  CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
	  CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)

	  M_coulp_ident = M_coulp_ident + coef1*coef2*exp_coeff_int*coef3* 
     & CG_j1_k1 * CG_j2_k2 * CG_j1_j2_p * CG_j1_j2_pp *
     & CG_l1_l2 * CG_j1_l1 * CG_j2_l2 * sign_correction
      ENDDO	  
      ENDDO  
      ENDDO
      ENDDO
	  
	  else
	  print*,"Something is wrong in PES expansion matrix computation."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
      
	  M_coulp_ident_1 = M_coulp_ident
      M_coulp = M_coulp + M_coulp_ident*
     & (-1)**j12_p_t*par_p
	 
      M_coulp_ident = 0d0
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j2_ch(channpp)
      j2_pp_t = j1_ch(channpp)	 

	  if(.not.bikram_identical_pes) then
	  
	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
      DO i=1,nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,i)
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(4d0*pi)
	  
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
	  CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,0,0,0,bikram_w3j_fact)
	  CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,0,0,0,bikram_w3j_fact)
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
	  CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)
	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
	  CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
	  CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
	  CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
	  CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)

	  M_coulp_ident = M_coulp_ident + coef1*coef2*exp_coeff_int*coef3* 
     & CG_j1_k1 * CG_j2_k2 * CG_j1_j2_p * CG_j1_j2_pp *
     & CG_l1_l2 * CG_j1_l1 * CG_j2_l2
      ENDDO	  
      ENDDO  
      ENDDO

	  else if(bikram_identical_pes) then
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for identical molecule-molecule collision using expansion
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022

	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
      DO i=1,nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,i)
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(4d0*pi)
	  sign_correction = 1 
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
	  CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)
	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
	  CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
	  
	  DO symmetry_coeff = 1, 2 - KRONEKER(l1_t,l2_t)
      if(symmetry_coeff.eq.1) then
	  sign_correction = 1
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      else if(symmetry_coeff.eq.2) then
      sign_correction = (-1)**(l1_t+l2_t)
      l1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(1,ind_t)
      l_t = A_TOP(3,ind_t)
	  else
	  if(myid.eq.0 .and. i_r_point.eq.1) 
     & print*, "Error in identical swapping variable."
	  stop
      end if
	  IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
	  
	  CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,0,0,0,bikram_w3j_fact)
	  CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,0,0,0,bikram_w3j_fact)
	  CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
	  CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
	  CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)

	  M_coulp_ident = M_coulp_ident + coef1*coef2*exp_coeff_int*coef3* 
     & CG_j1_k1 * CG_j2_k2 * CG_j1_j2_p * CG_j1_j2_pp *
     & CG_l1_l2 * CG_j1_l1 * CG_j2_l2 * sign_correction
      ENDDO	  
      ENDDO  
      ENDDO
      ENDDO
	  
	  else
	  print*,"Something is wrong in PES expansion matrix computation."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if

      M_coulp_ident_2 = M_coulp_ident
	  sign_correction = +1
	  if((-1)**j12_p_t*(-1)**j12_pp_t.lt.0.d0) then
	  if(M_coulp_ident_1*M_coulp_ident_2.gt.0.d0) sign_correction = -1
	  else if((-1)**j12_p_t*(-1)**j12_pp_t.gt.0.d0) then
	  if(M_coulp_ident_1*M_coulp_ident_2.lt.0.d0) sign_correction = -1
	  end if

      M_coulp = M_coulp + M_coulp_ident*
     & (-1)**j12_pp_t*par_pp*sign_correction

      M_coulp_ident = 0d0
      j1_p_t = j2_ch(channp)
      j2_p_t = j1_ch(channp)
      j1_pp_t = j2_ch(channpp)
      j2_pp_t = j1_ch(channpp)

	  if(.not.bikram_identical_pes) then
	  
	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
      DO i=1,nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,i)
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(4d0*pi)
	  
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
	  CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,0,0,0,bikram_w3j_fact)
	  CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,0,0,0,bikram_w3j_fact)
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
	  CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)
	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
	  CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
	  CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
	  CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
	  CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)

	  M_coulp_ident = M_coulp_ident + coef1*coef2*exp_coeff_int*coef3* 
     & CG_j1_k1 * CG_j2_k2 * CG_j1_j2_p * CG_j1_j2_pp *
     & CG_l1_l2 * CG_j1_l1 * CG_j2_l2
      ENDDO	  
      ENDDO  
      ENDDO

	  else if(bikram_identical_pes) then
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for identical molecule-molecule collision using expansion
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022

	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
      DO i=1,nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,i)
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(4d0*pi)
	  sign_correction = 1 
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
	  CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)
	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
	  CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
	  
	  DO symmetry_coeff = 1, 2 - KRONEKER(l1_t,l2_t)
      if(symmetry_coeff.eq.1) then
	  sign_correction = 1
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      else if(symmetry_coeff.eq.2) then
      sign_correction = (-1)**(l1_t+l2_t)
      l1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(1,ind_t)
      l_t = A_TOP(3,ind_t)
	  else
	  if(myid.eq.0 .and. i_r_point.eq.1) 
     & print*, "Error in identical swapping variable."
	  stop
      end if
	  IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
	  
	  CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,0,0,0,bikram_w3j_fact)
	  CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,0,0,0,bikram_w3j_fact)
	  CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
	  CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
	  CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)

	  M_coulp_ident = M_coulp_ident + coef1*coef2*exp_coeff_int*coef3* 
     & CG_j1_k1 * CG_j2_k2 * CG_j1_j2_p * CG_j1_j2_pp *
     & CG_l1_l2 * CG_j1_l1 * CG_j2_l2 * sign_correction
      ENDDO	  
      ENDDO  
      ENDDO
      ENDDO
	  
	  else
	  print*,"Something is wrong in PES expansion matrix computation."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if

	  M_coulp_ident_3 = M_coulp_ident
	  sign_correction = +1
	  if((-1)**j12_p_t*(-1)**j12_pp_t.lt.0.d0) then
	  if(M_coulp_ident_3*M_coulp_non_ident.gt.0.d0) 
     & sign_correction = -1
	  else if((-1)**j12_p_t*(-1)**j12_pp_t.gt.0.d0) then
	  if(M_coulp_ident_3*M_coulp_non_ident.lt.0.d0) 
     & sign_correction = -1
	  end if

      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**j12_p_t*par_p*
     & (-1)**j12_pp_t*par_pp*sign_correction

      M_coulp = M_coulp/
     & dsqrt(2d0*(1d0+delta(j1_p_t,j2_p_t)))
     & /dsqrt(2d0*(1d0+delta(j1_pp_t,j2_pp_t)))
      ENDIF
	  if(k == 1 .and. i_r_point == 1) write(*,'(4(e19.12,1x))')
     & M_coulp_non_ident, M_coulp_ident_3, 
     & M_coulp_ident_2, M_coulp_ident_1
	  return

! This peice below is from the old code for the identical collision partners
! which is not used anymore

      IF(identical_particles_defined) THEN
      par_p = parity_state(stp)
      par_pp = parity_state(stpp)
      j12_p_t = j12(stp)	  
      j12_pp_t = j12(stpp)
      M_coulp_ident = 0d0	  
      j1_p_t = j2_ch(channp)
      j2_p_t = j1_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)	  

      DO i=1,nterms
      IF(i_r_point.eq.i_nr_ini) THEN		  
      ind_t = i
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE		  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      DO mp = -j1_p_t,j1_p_t
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
	  
!      IF(myid.eq.0) PRINT*, "COULPING", M_coulp	  
      TERM_MATRIX_ELEMENT(1,1,2,i) = 
     & TERM_MATRIX_ELEMENT(1,1,2,i) + 
     & dsqrt((2d0*l_t+1d0)/4d0/dacos(-1d0))*
     & CG(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t)*
     & CG(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,m12_t)*
     & dsqrt((2d0*j1_pp_t+1)/(2d0*j1_p_t+1))*
     & dsqrt((2d0*j2_pp_t+1)/(2d0*j2_p_t+1))*
     & dsqrt((2d0*l1_t+1d0)*(2d0*l2_t+1d0))/4d0/pi*
     & CG(l1_t,l2_t,l_t,m_exp,-m_exp,0)*
     & CG(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp)*
     & CG(j1_pp_t,l1_t,j1_p_t,0,0,0)*
     & CG(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,m12_t-mp)*
     & CG(j2_pp_t,l2_t,j2_p_t,0,0,0)
      ENDDO	  
      ENDDO	
      ENDIF 
      M_coulp_ident = M_coulp_ident + 
     & expansion_terms(i_r_point,ind_t)*
     & TERM_MATRIX_ELEMENT(1,1,2,ind_t)	 
      ENDDO	
    
      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**j12_p_t*par_p
	 
      M_coulp_ident = 0d0
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j2_ch(channpp)
      j2_pp_t = j1_ch(channpp)	 
	 
      DO i=1,nterms
      IF(i_r_point.eq.i_nr_ini) THEN		  
      ind_t = i
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE		  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      DO mp = -j1_p_t,j1_p_t
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE	  
c      IF(myid.eq.0) PRINT*, "COULPING", M_coulp	  
      TERM_MATRIX_ELEMENT(1,1,3,ind_t) =
     & TERM_MATRIX_ELEMENT(1,1,3,ind_t) + 
     & dsqrt((2d0*l_t+1d0)/4d0/dacos(-1d0))*
     & CG(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t)*
     & CG(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,m12_t)*
     & dsqrt((2d0*j1_pp_t+1)/(2d0*j1_p_t+1))*
     & dsqrt((2d0*j2_pp_t+1)/(2d0*j2_p_t+1))*
     & dsqrt((2d0*l1_t+1d0)*(2d0*l2_t+1d0))/4d0/pi*
     & CG(l1_t,l2_t,l_t,m_exp,-m_exp,0)*
     & CG(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp)*
     & CG(j1_pp_t,l1_t,j1_p_t,0,0,0)*
     & CG(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,m12_t-mp)*
     & CG(j2_pp_t,l2_t,j2_p_t,0,0,0)

      ENDDO	  
      ENDDO
      ENDIF
      M_coulp_ident = M_coulp_ident + 
     & expansion_terms(i_r_point,ind_t)*
     & TERM_MATRIX_ELEMENT(1,1,3,ind_t)		 
      ENDDO	

      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**j12_pp_t*par_pp

      M_coulp_ident = 0d0
      j1_p_t = j2_ch(channp)
      j2_p_t = j1_ch(channp)
      j1_pp_t = j2_ch(channpp)
      j2_pp_t = j1_ch(channpp)

      DO i=1,nterms
      IF(i_r_point.eq.i_nr_ini) THEN	  
      ind_t = i
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      DO mp = -j1_p_t,j1_p_t
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE	  
c      IF(myid.eq.0) PRINT*, "COULPING", M_coulp	  
      TERM_MATRIX_ELEMENT(1,1,4,ind_t) =
     & TERM_MATRIX_ELEMENT(1,1,4,ind_t) +
     & dsqrt((2d0*l_t+1d0)/4d0/dacos(-1d0))*
     & CG(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t)*
     & CG(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,m12_t)*
     & dsqrt((2d0*j1_pp_t+1)/(2d0*j1_p_t+1))*
     & dsqrt((2d0*j2_pp_t+1)/(2d0*j2_p_t+1))*
     & dsqrt((2d0*l1_t+1d0)*(2d0*l2_t+1d0))/4d0/pi*
     & CG(l1_t,l2_t,l_t,m_exp,-m_exp,0)*
     & CG(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp)*
     & CG(j1_pp_t,l1_t,j1_p_t,0,0,0)*
     & CG(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,m12_t-mp)*
     & CG(j2_pp_t,l2_t,j2_p_t,0,0,0)

      ENDDO	  
      ENDDO
      ENDIF
      M_coulp_ident = M_coulp_ident + 
     & expansion_terms(i_r_point,ind_t)*
     & TERM_MATRIX_ELEMENT(1,1,4,ind_t)	
	  
      ENDDO		  
	  
	 
	 
      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**j12_pp_t*par_pp*par_pp*
     & (-1)**j12_p_t

      M_coulp = M_coulp/
     & dsqrt(2d0*(1d0+delta(j1_p_t,j2_p_t)))
     & /dsqrt(2d0*(1d0+delta(j1_pp_t,j2_pp_t)))

	 
     	 
      ENDIF	
      CASE(6)
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)
c      PRINT*, "CHANNP=",	channp, stp
c      PRINT*,"CHANNPP=",	channpp, stpp	  
      j12_p_t = j12(stp)	  
      m12_t = m12(stp)
      j12_pp_t = j12(stpp)
      IF(m12_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"   	  
c      PRINT*, j12_p_t,j12_pp_t,m12_t	  
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)
      v1_p = v1_ch(channp)
      v2_p = v2_ch(channp)
      v1_pp = v1_ch(channpp)
      v2_pp = v2_ch(channpp)	  
      DO i=1,nterms
      IF(i_r_point.eq.i_nr_ini) THEN	  
      ind_t = i!index_term_in_file(i)
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE		  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      DO mp = -j1_p_t,j1_p_t
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      TERM_MATRIX_ELEMENT(1,1,1,i) = 
     & TERM_MATRIX_ELEMENT(1,1,1,i) + 
     & dsqrt((2d0*l_t+1d0)/4d0/dacos(-1d0))*
     & CG(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t)*
     & CG(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,m12_t)*
     & dsqrt((2d0*j1_pp_t+1)/(2d0*j1_p_t+1))*
     & dsqrt((2d0*j2_pp_t+1)/(2d0*j2_p_t+1))*
     & dsqrt((2d0*l1_t+1d0)*(2d0*l2_t+1d0))/4d0/pi*
     & CG(l1_t,l2_t,l_t,m_exp,-m_exp,0)*
     & CG(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp)*
     & CG(j1_pp_t,l1_t,j1_p_t,0,0,0)*
     & CG(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,m12_t-mp)*
     & CG(j2_pp_t,l2_t,j2_p_t,0,0,0)	 
!      IF(myid.eq.0) PRINT*, "COULPING", 
!     & i_r_point,expansion_terms(i_r_point,1)
      ENDDO	  
      ENDDO
      ENDIF		  
      M_coulp = M_coulp + 
     & expansion_terms(i_r_point,i)*
     & TERM_MATRIX_ELEMENT(1,1,1,i)	  
      ENDDO
	  
c      IF(k.eq.2) PRINT*,i_r_point,  M_coulp  
      IF(identical_particles_defined) THEN
      par_p = parity_state(stp)
      par_pp = parity_state(stpp)
      j12_p_t = j12(stp)	  
      j12_pp_t = j12(stpp)
      M_coulp_ident = 0d0	  
      j1_p_t = j2_ch(channp)
      j2_p_t = j1_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)	  

      DO i=1,nterms
      IF(i_r_point.eq.i_nr_ini) THEN		  
      ind_t = i
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE		  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      DO mp = -j1_p_t,j1_p_t
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
	  
!      IF(myid.eq.0) PRINT*, "COULPING", M_coulp	  
      TERM_MATRIX_ELEMENT(1,1,2,i) = 
     & TERM_MATRIX_ELEMENT(1,1,2,i) + 
     & dsqrt((2d0*l_t+1d0)/4d0/dacos(-1d0))*
     & CG(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t)*
     & CG(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,m12_t)*
     & dsqrt((2d0*j1_pp_t+1)/(2d0*j1_p_t+1))*
     & dsqrt((2d0*j2_pp_t+1)/(2d0*j2_p_t+1))*
     & dsqrt((2d0*l1_t+1d0)*(2d0*l2_t+1d0))/4d0/pi*
     & CG(l1_t,l2_t,l_t,m_exp,-m_exp,0)*
     & CG(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp)*
     & CG(j1_pp_t,l1_t,j1_p_t,0,0,0)*
     & CG(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,m12_t-mp)*
     & CG(j2_pp_t,l2_t,j2_p_t,0,0,0)
      ENDDO	  
      ENDDO	
      ENDIF 
      M_coulp_ident = M_coulp_ident + 
     & expansion_terms(i_r_point,ind_t)*
     & TERM_MATRIX_ELEMENT(1,1,2,ind_t)	 
      ENDDO	
      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**j12_p_t*par_p
	 
      M_coulp_ident = 0d0
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j2_ch(channpp)
      j2_pp_t = j1_ch(channpp)	 
	 
      DO i=1,nterms
      IF(i_r_point.eq.i_nr_ini) THEN		  
      ind_t = i
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE		  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      DO mp = -j1_p_t,j1_p_t
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE	  
c      IF(myid.eq.0) PRINT*, "COULPING", M_coulp	  
      TERM_MATRIX_ELEMENT(1,1,3,ind_t) =
     & TERM_MATRIX_ELEMENT(1,1,3,ind_t) + 
     & dsqrt((2d0*l_t+1d0)/4d0/dacos(-1d0))*
     & CG(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t)*
     & CG(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,m12_t)*
     & dsqrt((2d0*j1_pp_t+1)/(2d0*j1_p_t+1))*
     & dsqrt((2d0*j2_pp_t+1)/(2d0*j2_p_t+1))*
     & dsqrt((2d0*l1_t+1d0)*(2d0*l2_t+1d0))/4d0/pi*
     & CG(l1_t,l2_t,l_t,m_exp,-m_exp,0)*
     & CG(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp)*
     & CG(j1_pp_t,l1_t,j1_p_t,0,0,0)*
     & CG(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,m12_t-mp)*
     & CG(j2_pp_t,l2_t,j2_p_t,0,0,0)

      ENDDO	  
      ENDDO
      ENDIF
      M_coulp_ident = M_coulp_ident + 
     & expansion_terms(i_r_point,ind_t)*
     & TERM_MATRIX_ELEMENT(1,1,3,ind_t)		 
      ENDDO	

      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**j12_pp_t*par_pp

      M_coulp_ident = 0d0
      j1_p_t = j2_ch(channp)
      j2_p_t = j1_ch(channp)
      j1_pp_t = j2_ch(channpp)
      j2_pp_t = j1_ch(channpp)

      DO i=1,nterms
      IF(i_r_point.eq.i_nr_ini) THEN	  
      ind_t = i
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      DO mp = -j1_p_t,j1_p_t
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE	  
c      IF(myid.eq.0) PRINT*, "COULPING", M_coulp	  
      TERM_MATRIX_ELEMENT(1,1,4,ind_t) =
     & TERM_MATRIX_ELEMENT(1,1,4,ind_t) +
     & dsqrt((2d0*l_t+1d0)/4d0/dacos(-1d0))*
     & CG(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t)*
     & CG(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,m12_t)*
     & dsqrt((2d0*j1_pp_t+1)/(2d0*j1_p_t+1))*
     & dsqrt((2d0*j2_pp_t+1)/(2d0*j2_p_t+1))*
     & dsqrt((2d0*l1_t+1d0)*(2d0*l2_t+1d0))/4d0/pi*
     & CG(l1_t,l2_t,l_t,m_exp,-m_exp,0)*
     & CG(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp)*
     & CG(j1_pp_t,l1_t,j1_p_t,0,0,0)*
     & CG(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,m12_t-mp)*
     & CG(j2_pp_t,l2_t,j2_p_t,0,0,0)

      ENDDO	  
      ENDDO
      ENDIF
      M_coulp_ident = M_coulp_ident + 
     & expansion_terms(i_r_point,ind_t)*
     & TERM_MATRIX_ELEMENT(1,1,4,ind_t)	
	  
      ENDDO		  
 
	 
	 
      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**j12_pp_t*par_pp*par_pp*
     & (-1)**j12_p_t

      M_coulp = M_coulp/
     & dsqrt(2d0*(1d0+delta(j1_p_t,j2_p_t)))
     & /dsqrt(2d0*(1d0+delta(j1_pp_t,j2_pp_t)))

	 
     	 
      ENDIF		  
       M_coulp = M_coulp * vib_overlap(channp,channpp)	  
      CASE(7)
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)
      j12_p_t = j12(stp)	  
      m12_t = m12(stp)
      j12_pp_t = j12(stpp)
      IF(m12_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)
      k_p_t = k1_ch(channp)
      k_pp_t = k1_ch(channpp)	  
      eps_p_t = eps1_ch(channp)		  
      eps_pp_t = eps1_ch(channpp)	  

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for any molecule-molecule collision using expansion
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022
	  
	  if(.not.bikram_identical_pes) then

	  same_terms = .false.
	  if(i_r_point.eq.1 .and. k.eq.1) then
      DO i = 1,nterms
      ind_t = i
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      l_t = A_TOP(4,ind_t)
	  
	  if(nju1_t.lt.0) then
	  write(*,'(a, a, 5(i0,a))')"The values of M1 should be ",
     & "always positive. Please check the expansion term: ",
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_terms = .true.
	  end if
	  end do
	  end if
	  if(same_terms) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  	  
	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  coef3 = 1d0 / dsqrt(2d0*(1d0 + KRONEKER(k_p_t,0))*
     &                    2d0*(1d0 + KRONEKER(k_pp_t,0)))
	  
      DO i = 1,nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,ind_t)
      l1_t = A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      l_t = A_TOP(4,ind_t)
	  
	  coef4 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(4d0*pi)
      CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,0,0,0,bikram_w3j_fact)
	  
      DO mp = -j1_p_t, j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)
	  
      DO m_exp = -min(l1_t,l2_t), min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE	  
      CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
      CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
      CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
      CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)	  

      DO planr_coeff = 1, 2-KRONEKER(nju1_t,0)	
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int=exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t)	  
      nju1_t = -nju1_t
      ENDIF	
! This piece below is to compute four terms withn 
! the square bracket of table from the MQCT 2022 manual
!      DO par_p =  1, 2
!      DO par_pp = 1, 2
!      k_p_t = (-1)**(par_p-1)*k1_ch(channp)
!      k_pp_t = (-1)**(par_pp-1)*k1_ch(channpp)
!      IF(k_pp_t.ne.k_p_t-nju1_t) CYCLE
	  
!      CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,k_pp_t,nju1_t,k_p_t,
!     & bikram_w3j_fact)
!     & *((-1)**eps_p_t)**(par_p-1)*((-1)**eps_pp_t)**(par_pp-1)

	  CG_in1a = CG_bikram(j1_pp_t, l1_t, j1_p_t, k_pp_t, nju1_t, k_p_t,
     & bikram_w3j_fact)
	  CG_in2a = CG_bikram(j1_pp_t, l1_t, j1_p_t, k_pp_t, nju1_t,-k_p_t,
     & bikram_w3j_fact)
	  CG_in3a = CG_bikram(j1_pp_t, l1_t, j1_p_t,-k_pp_t, nju1_t, k_p_t,
     & bikram_w3j_fact)
	  CG_in4a = CG_bikram(j1_pp_t, l1_t, j1_p_t,-k_pp_t, nju1_t,-k_p_t,
     & bikram_w3j_fact)
	 
	  CG_in_all = CG_in1a + 
     &            CG_in2a * ((-1)**eps_p_t) +
     &            CG_in3a * ((-1)**eps_pp_t) +
     &            CG_in4a * ((-1)**eps_p_t)*((-1)**eps_pp_t)
	 
      M_coulp = M_coulp + coef1*coef2*coef3*exp_coeff_int*coef4* 
     & CG_j1_j2_p*
     & CG_j1_j2_pp*
     & CG_l1_l2*
     & CG_j1_l1*
     & CG_in_all*!CG_j1_k1*
     & CG_j2_l2*
     & CG_j2_k2
	  
	  if(bikram_rebalance) tmp1 = tmp1 + 1
!      ENDDO	  
!      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
	  
	  else
	  print*,"Something is wrong in PES expansion matrix computation."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
! Bikram Start: Feb 2022
! This part is to use recursive method of computing Wigner 3j symbol
	  if(M_coulp.ne.M_coulp) then
	  if(myid.eq.0) then
	  print*, "Computed matrix element is wrong for state:"
      write(*,'(6(a,i0),2x,e19.12)')'ini: ',j_p_t,',',k_ch(channp),
     & ',',eps_p_t, '  fin: ', j_pp_t,',',k_ch(channpp),',',
     & eps_pp_t, M_coulp
	  print*, "Program will stop now."
	  if(bikram_w3j_fact) print*, "Recommended to use recursive
     & algorithm for computing Wigner 3j coefficients."
	  stop
      CALL MPI_FINALIZE (ierr_mpi)
	  end if
	  stop
      CALL MPI_FINALIZE (ierr_mpi)
	  end if
! Bikram End.

      CASE(8)
      CALL ASYM_TOP_VECTORS 	  
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)	  
      j12_p_t = j12(stp)	  
      m12_t = m12(stp)
      j12_pp_t = j12(stpp)
      IF(m12_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"
	  
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for any molecule-molecule collision using expansion
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022
	  
	  if(.not.bikram_identical_pes) then

	  same_terms = .false.
	  if(i_r_point.eq.1 .and. k.eq.1) then
      DO i = 1,nterms
      ind_t = i
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      l_t = A_TOP(4,ind_t)
	  
	  if(nju1_t.lt.0) then
	  write(*,'(a, a, 5(i0,a))')"The values of M1 should be ",
     & "always positive. Please check the expansion term: ",
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_terms = .true.
	  end if
	  end do
	  end if
	  if(same_terms) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
      DO i = 1, nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,ind_t)
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      l_t = A_TOP(4,ind_t)
	  
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(4d0*pi)
      CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,0,0,0,bikram_w3j_fact)
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
      CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
      CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
      CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
      CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)
	  
      DO planr_coeff = 1, 2-KRONEKER(nju1_t,0)
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int=exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t)	  
      nju1_t = -nju1_t
      ENDIF	  
	  
      DO k1p = -j1_p_t, j1_p_t
      k1pp = k1p - nju1_t
      IF(abs(k1pp).gt.j1_pp_t) CYCLE	  
      CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,k1pp,nju1_t,k1p,
     & bikram_w3j_fact)
      coeff_1_p = M_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF

      M_coulp = M_coulp + coef1 * coef2 * exp_coeff_int * coef3 * 
     & CG_j1_j2_p * CG_j1_j2_pp * CG_l1_l2 *
     & CG_j1_l1 * CG_j1_k1 * CG_j2_l2 * CG_j2_k2 *
     & coeff_1_p * coeff_1_pp
	  
	  if(bikram_rebalance) tmp1 = tmp1 + 1
      ENDDO	  
      ENDDO	  
      ENDDO
      ENDDO
      ENDDO	  	  
	  
	  else
	  print*,"Something is wrong in PES expansion matrix computation."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
! Bikram Start: Feb 2022
! This part is to use recursive method of computing Wigner 3j symbol
	  if(M_coulp.ne.M_coulp) then
	  if(myid.eq.0) then
	  print*, "Computed matrix element is wrong for state:"
      write(*,'(6(a,i0),2x,e19.12)')'ini: ',j_p_t,',',k_ch(channp),
     & ',',eps_p_t, '  fin: ', j_pp_t,',',k_ch(channpp),',',
     & eps_pp_t, M_coulp
	  print*, "Program will stop now."
	  if(bikram_w3j_fact) print*, "Recommended to use recursive
     & algorithm for computing Wigner 3j coefficients."
	  stop
      CALL MPI_FINALIZE (ierr_mpi)
	  end if
	  stop
      CALL MPI_FINALIZE (ierr_mpi)
	  end if
! Bikram End.

      CASE(9)

      CALL ASYM_TOP_VECTORS 	  
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)
      j12_p_t = j12(stp)	  
      m12_t = m12(stp)
      j12_pp_t = j12(stpp)
      IF(m12_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"   	  

      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)
      k_p_t = k2_ch(channp)
      k_pp_t = k2_ch(channpp)	  
      eps_p_t = eps2_ch(channp)		  
      eps_pp_t = eps2_ch(channpp)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for any molecule-molecule collision using expansion
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022
	  
	  if(.not.bikram_identical_pes) then
	  
	  same_terms = .false.
	  if(i_r_point.eq.1 .and. k.eq.1) then
      DO i = 1,nterms
      ind_t = i
      l1_t = A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  
	  if(nju1_t.eq.0 .and. nju2_t.lt.0) then
	  write(*,'(a, a, a, 5(i0,a))')"The values of M2 should be ",
     & "positive only for M1 = 0. ", 
     & "Please check the expansion term: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_terms = .true.
	  end if
	  if(nju1_t.lt.0) then
	  write(*,'(a, a, 5(i0,a))')"The values of M1 should be ",
     & "always positive. Please check the expansion term: ",
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_terms = .true.
	  end if
	  end do
	  end if
	  if(same_terms) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  coef3 = 1d0 / dsqrt(2d0*(1d0 + KRONEKER(k_p_t,0))*
     &                    2d0*(1d0 + KRONEKER(k_pp_t,0)))
	  
	  DO i = 1, nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,ind_t)	 
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  
	  coef4 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(8d0*pi**2)
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)
	  DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
      CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
      CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
      CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
      CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)
	  
      DO planr_coeff = 1, 2-KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int = exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t)	  
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      ENDIF

      DO k1p = -j1_p_t,j1_p_t
      k1pp = k1p - nju1_t
      IF(abs(k1pp).gt.j1_pp_t) CYCLE	  
      CG_j1_k1 = CG(j1_pp_t,l1_t,j1_p_t,k1pp,nju1_t,k1p)	  
      coeff_1_p = M_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF
	  
	  CG_in1a = CG_bikram(j2_pp_t, l2_t, j2_p_t, k_pp_t, nju2_t, k_p_t,
     & bikram_w3j_fact)
	  CG_in2a = CG_bikram(j2_pp_t, l2_t, j2_p_t, k_pp_t, nju2_t,-k_p_t,
     & bikram_w3j_fact)
	  CG_in3a = CG_bikram(j2_pp_t, l2_t, j2_p_t,-k_pp_t, nju2_t, k_p_t,
     & bikram_w3j_fact)
	  CG_in4a = CG_bikram(j2_pp_t, l2_t, j2_p_t,-k_pp_t, nju2_t,-k_p_t,
     & bikram_w3j_fact)
	 
	  CG_in_all = CG_in1a + 
     &            CG_in2a * ((-1)**eps_p_t) +
     &            CG_in3a * ((-1)**eps_pp_t) +
     &            CG_in4a * ((-1)**eps_p_t)*((-1)**eps_pp_t)

	  M_coulp = M_coulp + coef1*coef2*coef3*exp_coeff_int*coef4*
     & CG_j1_j2_p * CG_j1_j2_pp * CG_l1_l2 *
     & CG_j1_l1 * CG_j2_l2 * CG_j1_k1 * CG_in_all *
     & coeff_1_p*coeff_1_pp
	  
	  if(bikram_rebalance) tmp1 = tmp1 + 1
      ENDDO
      ENDDO
      ENDDO	  
      ENDDO
      ENDDO
	  return
	  
	  else
	  print*,"Something is wrong in PES expansion matrix computation."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if


      DO i=1,nterms
      ind_t =i
      exp_coeff_int=expansion_terms(i_r_point,ind_t)	  
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
!      IF(l1_t.ne.3 .and. l1_t.ne.0 ) CYCLE
!      IF(l2_t.ne.0) CYCLE
! HERE IS THE MISTAKE
      DO planr_coeff =1,2-KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	  
      DO par_p = 1,2
      DO par_pp =1,2	  
      symmetry_coeff = 2*par_p+par_pp-2
      k_p_t = (-1)**(par_p-1)*k2_ch(channp)
      k_pp_t = (-1)**(par_pp-1)*k2_ch(channpp)
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int=exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t)	  
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      ENDIF	
      IF(k_pp_t.ne.k_p_t-nju2_t) CYCLE		  
      IF(i_r_point.eq.i_nr_ini) THEN
      TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t) = 0d0	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t)	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE	  
      CG_j1_j2_pp = CG(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,m12_t)
      CG_l1_l2 = CG(l1_t,l2_t,l_t,m_exp,-m_exp,0)
      CG_j1_l1 = CG(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp)
      CG_j2_l2 = CG(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,m12_t-mp)	  
      DO k1p = -j1_p_t,j1_p_t
      k1pp = k1p - nju1_t
      IF(abs(k1pp).gt.j1_pp_t) CYCLE	  
      CG_j1_k1 = CG(j1_pp_t,l1_t,j1_p_t,k1pp,nju1_t,k1p)	  
      coeff_1_p = M_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF	  

      CG_j2_k2 = CG(j2_pp_t,l2_t,j2_p_t,k_pp_t,nju2_t,k_p_t)*	 
!     & eps_p_t**(par_p-1)*eps_pp_t**(par_pp-1)
     & *((-1)**eps_p_t)**(par_p-1)*((-1)**eps_pp_t)**(par_pp-1)
     % /dsqrt(2d0*(1d0 + KRONEKER(k_p_t,0)))
     & /dsqrt(2d0*(1d0 + KRONEKER(k_pp_t,0)))	
	  

      TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t) =
     & TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t)+ 	  
     & CG_j1_j2_p*
     & CG_j1_j2_pp*
     & dsqrt((2d0*j1_pp_t+1)/(2d0*j1_p_t+1))*
     & dsqrt((2d0*j2_pp_t+1)/(2d0*j2_p_t+1))*
     & dsqrt((2d0*l1_t+1d0)*(2d0*l2_t+1d0))/8d0/pi**2*
     & CG_l1_l2*
     & CG_j1_l1*
     & CG_j1_k1*
     & CG_j2_l2*
     & CG_j2_k2*coeff_1_p*coeff_1_pp
  
      ENDDO	  
      ENDDO	  
      ENDDO
      ENDIF
      matrix_exp_coefficent = 
     & TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t)
      M_coulp = M_coulp +
     & exp_coeff_int*matrix_exp_coefficent	 
      ENDDO
      ENDDO
      ENDDO	  
      ENDDO	  
      CASE(0)
      CALL ASYM_TOP_VECTORS
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)
      j12_p_t = j12(stp)
      m12_t = m12(stp)
      j12_pp_t = j12(stpp)
      IF(m12_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)
	  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for any molecule-molecule collision using expansion
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022
	  
	  if(.not.bikram_identical_pes) then
	  
	  same_terms = .false.
	  if(i_r_point.eq.1 .and. k.eq.1) then
      DO i = 1,nterms
      ind_t = i
      l1_t = A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  
	  if(nju1_t.eq.0 .and. nju2_t.lt.0) then
	  write(*,'(a, a, a, 5(i0,a))')"The values of M2 should be ",
     & "positive only for M1 = 0. ", 
     & "Please check the expansion term: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_terms = .true.
	  end if
	  if(nju1_t.lt.0) then
	  write(*,'(a, a, 5(i0,a))')"The values of M1 should be ",
     & "always positive. Please check the expansion term: ",
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_terms = .true.
	  end if
	  end do
	  end if
	  if(same_terms) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
	  DO i=1,nterms
      ind_t =i
      exp_coeff_int=expansion_terms(i_r_point,ind_t)	 
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(8d0*pi**2)
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)
	  DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
      CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
      CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
      CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
      CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)
	  
      DO planr_coeff =1,2-KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	 
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int=exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t)	  
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      ENDIF
	  
	  DO k1p = -j1_p_t,j1_p_t
      k1pp = k1p - nju1_t	  
      IF(abs(k1pp).gt.j1_pp_t) CYCLE	 	  
      CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,k1pp,nju1_t,k1p,
     & bikram_w3j_fact)
      coeff_1_p = M1_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M1_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF	  
      DO k2p =-j2_p_t,j2_p_t
      k2pp = k2p - nju2_t
      IF(abs(k2pp).gt.j2_pp_t) CYCLE	  
      coeff_2_p = M2_VECTORS(channp,j2_p_t+1+k2p)
      CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,k2pp,nju2_t,k2p,
     & bikram_w3j_fact)
      IF(abs(k2pp).le.j2_pp_t) THEN
      coeff_2_pp = M2_VECTORS(channpp,j2_pp_t+1+k2pp)
      ELSE
      coeff_2_pp = 0d0       	  
      ENDIF

	  M_coulp = M_coulp + coef1*coef2*exp_coeff_int*coef3*
     & CG_j1_j2_p * CG_j1_j2_pp * CG_l1_l2 *
     & CG_j1_l1 * CG_j2_l2 * CG_j1_k1 * CG_j2_k2 *
     & coeff_1_p*coeff_2_p*coeff_1_pp*coeff_2_pp

      ENDDO	  
      ENDDO
      ENDDO
      ENDDO	  
      ENDDO
      ENDDO	  
	  
	  else if(bikram_identical_pes) then
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for identical molecule-molecule collision using expansion
! an example would be H2O-H2O, honestly, this derivation was done using water collision
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022

	  same_terms = .false.
	  if(i_r_point.eq.1 .and. k.eq.1) then
      DO i = 1, nterms
      ind_t = i
      l1_t = A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  
	  if(nju1_t.eq.0 .and. nju2_t.lt.0) then
	  write(*,'(a, a, a, 5(i0,a))')"The values of M2 should be ",
     & "positive only for M1 = 0. ", 
     & "Please check the expansion term: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_terms = .true.
	  end if
	  if(nju1_t.lt.0) then
	  write(*,'(a, a, 5(i0,a))')"The values of M1 should be ",
     & "always positive. Please check the expansion term: ",
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_terms = .true.
	  end if
	  if(l1_t.lt.l2_t) then
	  write(*,'(a, a, a, 5(i0,a))')"The values of L2 should be ",
     & "less or equal to L1. ", 
     & "Please check the expansion term: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_terms = .true.
	  end if
	  if(l1_t.eq.l2_t .and. abs(nju1_t).lt.abs(nju2_t)) then
	  write(*,'(a, a, a, 5(i0,a))')"The values of M1 should be ",
     & "less or equal to abs(M2) for L1 = L2. ", 
     & "Please check the expansion term: ", 
     & l1_t, ',', nju1_t, ',', l2_t, ',', nju2_t, ',', l_t
	  same_terms = .true.
	  end if
	  end do
	  end if
	  if(same_terms) then
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
	  
	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
      DO i = 1, nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,ind_t)	 
      l1_t = A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  	  
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(8d0*pi**2)
	  sign_correction = 1 
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)

	  DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
	  
      DO symmetry_coeff = 1, 
     & 2 - KRONEKER(l1_t,l2_t)*KRONEKER(abs(nju1_t),abs(nju2_t))
      if(symmetry_coeff.eq.1) then
	  sign_correction = 1 
      l1_t = A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
      else if(symmetry_coeff.eq.2) then
!      exp_coeff_int = expansion_terms(i_r_point,ind_t)*(-1)**(l1_t+l2_t)	  
      sign_correction = (-1)**(l1_t+l2_t)	  
      l2_t = A_TOP(1,ind_t)
      nju2_t = A_TOP(2,ind_t)
      l1_t = A_TOP(3,ind_t)
      nju1_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  else
	  if(myid.eq.0 .and. i_r_point.eq.1) 
     & print*, "Error in identical swapping variable."
	  stop
      end if
	  
      DO planr_coeff = 1, 2 - KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	 
      IF(planr_coeff.eq.2) THEN
!      exp_coeff_int = exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t) 
      sign_correction = sign_correction
     & *(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t) 
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      ENDIF
	  	  
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
      CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
      CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
      CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)
	  
	  DO k1p = -j1_p_t,j1_p_t
      k1pp = k1p - nju1_t	  
      IF(abs(k1pp).gt.j1_pp_t) CYCLE	 	  
      CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,k1pp,nju1_t,k1p,
     & bikram_w3j_fact)
      coeff_1_p = M1_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M1_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF	  
      DO k2p =-j2_p_t,j2_p_t
      k2pp = k2p - nju2_t
      IF(abs(k2pp).gt.j2_pp_t) CYCLE	  
      coeff_2_p = M2_VECTORS(channp,j2_p_t+1+k2p)
      CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,k2pp,nju2_t,k2p,
     & bikram_w3j_fact)
      IF(abs(k2pp).le.j2_pp_t) THEN
      coeff_2_pp = M2_VECTORS(channpp,j2_pp_t+1+k2pp)
      ELSE
      coeff_2_pp = 0d0       	  
      ENDIF

	  M_coulp = M_coulp + coef1*coef2*exp_coeff_int*coef3*
     & CG_j1_j2_p * CG_j1_j2_pp * CG_l1_l2 *
     & CG_j1_l1 * CG_j2_l2 * CG_j1_k1 * CG_j2_k2 *
     & coeff_1_p*coeff_2_p*coeff_1_pp*coeff_2_pp*sign_correction
	  
	  if(bikram_rebalance) tmp1 = tmp1 + 1
      ENDDO	  
      ENDDO	  
      ENDDO
      ENDDO
      ENDDO	  
      ENDDO
      ENDDO
	  
	  else
	  print*,"Something is wrong in PES expansion matrix computation."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if
	  
! Bikram Start: Feb 2022
! This part is to use recursive method of computing Wigner 3j symbol
	  if(M_coulp.ne.M_coulp) then
	  if(myid.eq.0) then
	  print*, "Computed matrix element is wrong for state:"
      write(*,'(6(a,i0),2x,e19.12)')'ini: ',j_p_t,',',k_ch(channp),
     & ',',eps_p_t, '  fin: ', j_pp_t,',',k_ch(channpp),',',
     & eps_pp_t, M_coulp
	  print*, "Program will stop now."
	  if(bikram_w3j_fact) print*, "Recommended to use recursive
     & algorithm for computing Wigner 3j coefficients."
	  stop
      CALL MPI_FINALIZE (ierr_mpi)
	  end if
	  stop
      CALL MPI_FINALIZE (ierr_mpi)
	  end if
! Bikram End.

! This piece for identical collision partners
      IF(identical_particles_defined) THEN
      M_coulp_non_ident = M_coulp	  

! Bikram Start:
	  bk_par1 = parity_state_bk(stp)
	  bk_par2 = parity_state_bk(stpp)
	  if(stp.eq.stpp .and. j1_ch(channp).eq.j2_ch(channp) .and.
     & ka1_ch(channp).eq.ka2_ch(channp) .and. 
     & kc1_ch(channp).eq.kc2_ch(channp)) then
      bk_par1 = (-1)**j12(stp)
      bk_par2 = (-1)**j12(stpp)
	  end if
! Bikram End.
      par_p = parity_state(stp)
      par_pp = parity_state(stpp)
      j12_p_t = j12(stp)	  
      j12_pp_t = j12(stpp)
      M_coulp_ident = 0d0	  
      j1_p_t = j2_ch(channp)
      j2_p_t = j1_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)	 

	  if(.not.bikram_identical_pes) then
	  
	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
	  DO i=1,nterms
      ind_t =i
      exp_coeff_int=expansion_terms(i_r_point,ind_t)	 
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(8d0*pi**2)
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)
	  DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
      CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
      CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
      CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
      CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)
	  
      DO planr_coeff =1,2-KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	 
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int=exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t)	  
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      ENDIF
	  
	  DO k1p = -j1_p_t,j1_p_t
      k1pp = k1p - nju1_t	  
      IF(abs(k1pp).gt.j1_pp_t) CYCLE	 	  
      CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,k1pp,nju1_t,k1p,
     & bikram_w3j_fact)
      coeff_1_p = M2_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M1_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF	  
      DO k2p =-j2_p_t,j2_p_t
      k2pp = k2p - nju2_t
      IF(abs(k2pp).gt.j2_pp_t) CYCLE	  
      coeff_2_p = M1_VECTORS(channp,j2_p_t+1+k2p)
      CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,k2pp,nju2_t,k2p,
     & bikram_w3j_fact)
      IF(abs(k2pp).le.j2_pp_t) THEN
      coeff_2_pp = M2_VECTORS(channpp,j2_pp_t+1+k2pp)
      ELSE
      coeff_2_pp = 0d0       	  
      ENDIF

	  M_coulp_ident = M_coulp_ident + coef1*coef2*exp_coeff_int*coef3*
     & CG_j1_j2_p * CG_j1_j2_pp * CG_l1_l2 *
     & CG_j1_l1 * CG_j2_l2 * CG_j1_k1 * CG_j2_k2 *
     & coeff_1_p*coeff_2_p*coeff_1_pp*coeff_2_pp

      ENDDO	  
      ENDDO
      ENDDO
      ENDDO	  
      ENDDO
      ENDDO	  
	  
	  else if(bikram_identical_pes) then
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for identical molecule-molecule collision using expansion
! an example would be H2O-H2O, honestly, this derivation was done using water collision
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022

	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
      DO i = 1, nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,ind_t)	 
      l1_t = A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  	  
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(8d0*pi**2)
	  sign_correction = 1 
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)

	  DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
	  
      DO symmetry_coeff = 1, 
     & 2 - KRONEKER(l1_t,l2_t)*KRONEKER(abs(nju1_t),abs(nju2_t))
      if(symmetry_coeff.eq.1) then
	  sign_correction = 1 
      l1_t = A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
      else if(symmetry_coeff.eq.2) then
      sign_correction = (-1)**(l1_t+l2_t)	  
      l2_t = A_TOP(1,ind_t)
      nju2_t = A_TOP(2,ind_t)
      l1_t = A_TOP(3,ind_t)
      nju1_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  else
	  if(myid.eq.0 .and. i_r_point.eq.1) 
     & print*, "Error in identical swapping variable."
	  stop
      end if
	  
      DO planr_coeff = 1, 2 - KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	 
      IF(planr_coeff.eq.2) THEN
      sign_correction = sign_correction
     & *(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t) 
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      ENDIF
	  	  
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
      CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
      CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
      CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)
	  
	  DO k1p = -j1_p_t,j1_p_t
      k1pp = k1p - nju1_t	  
      IF(abs(k1pp).gt.j1_pp_t) CYCLE	 	  
      CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,k1pp,nju1_t,k1p,
     & bikram_w3j_fact)
      coeff_1_p = M2_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M1_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF	  
      DO k2p =-j2_p_t,j2_p_t
      k2pp = k2p - nju2_t
      IF(abs(k2pp).gt.j2_pp_t) CYCLE	  
      coeff_2_p = M1_VECTORS(channp,j2_p_t+1+k2p)
      CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,k2pp,nju2_t,k2p,
     & bikram_w3j_fact)
      IF(abs(k2pp).le.j2_pp_t) THEN
      coeff_2_pp = M2_VECTORS(channpp,j2_pp_t+1+k2pp)
      ELSE
      coeff_2_pp = 0d0       	  
      ENDIF

	  M_coulp_ident = M_coulp_ident + coef1*coef2*exp_coeff_int*coef3*
     & CG_j1_j2_p * CG_j1_j2_pp * CG_l1_l2 *
     & CG_j1_l1 * CG_j2_l2 * CG_j1_k1 * CG_j2_k2 *
     & coeff_1_p*coeff_2_p*coeff_1_pp*coeff_2_pp*sign_correction
	  
	  if(bikram_rebalance) tmp1 = tmp1 + 1
      ENDDO	  
      ENDDO	  
      ENDDO
      ENDDO
      ENDDO	  
      ENDDO
      ENDDO
	  
	  else
	  print*,"Something is wrong in PES expansion matrix computation."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if

      M_coulp_ident_1 = M_coulp_ident
      M_coulp = M_coulp + M_coulp_ident*par_p
     & *bk_par1
	  	  
      M_coulp_ident = 0d0
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j2_ch(channpp)
      j2_pp_t = j1_ch(channpp)
      M_coulp_simplified =M_coulp/
     & dsqrt(1d0+delta(j1_p_t,j2_p_t))
     & /dsqrt(1d0+delta(j1_pp_t,j2_pp_t))	  
      IF(SIMPLIFICATION_EXP_MAT) THEN	  
      M_coulp = M_coulp_simplified
      RETURN	  
      ENDIF

	  if(.not.bikram_identical_pes) then
	  
	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
	  DO i=1,nterms
      ind_t =i
      exp_coeff_int=expansion_terms(i_r_point,ind_t)	 
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(8d0*pi**2)
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)
	  DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
      CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
      CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
      CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
      CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)
	  
      DO planr_coeff =1,2-KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	 
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int=exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t)	  
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      ENDIF
	  
	  DO k1p = -j1_p_t,j1_p_t
      k1pp = k1p - nju1_t	  
      IF(abs(k1pp).gt.j1_pp_t) CYCLE	 	  
      CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,k1pp,nju1_t,k1p,
     & bikram_w3j_fact)
      coeff_1_p = M1_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M2_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF	  
      DO k2p =-j2_p_t,j2_p_t
      k2pp = k2p - nju2_t
      IF(abs(k2pp).gt.j2_pp_t) CYCLE	  
      coeff_2_p = M2_VECTORS(channp,j2_p_t+1+k2p)
      CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,k2pp,nju2_t,k2p,
     & bikram_w3j_fact)
      IF(abs(k2pp).le.j2_pp_t) THEN
      coeff_2_pp = M1_VECTORS(channpp,j2_pp_t+1+k2pp)
      ELSE
      coeff_2_pp = 0d0       	  
      ENDIF

	  M_coulp_ident = M_coulp_ident + coef1*coef2*exp_coeff_int*coef3*
     & CG_j1_j2_p * CG_j1_j2_pp * CG_l1_l2 *
     & CG_j1_l1 * CG_j2_l2 * CG_j1_k1 * CG_j2_k2 *
     & coeff_1_p*coeff_2_p*coeff_1_pp*coeff_2_pp

      ENDDO	  
      ENDDO
      ENDDO
      ENDDO	  
      ENDDO
      ENDDO	  
	  
	  else if(bikram_identical_pes) then
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for identical molecule-molecule collision using expansion
! an example would be H2O-H2O, honestly, this derivation was done using water collision
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022

	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
      DO i = 1, nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,ind_t)	 
      l1_t = A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  	  
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(8d0*pi**2)
	  sign_correction = 1 
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)

	  DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
	  
      DO symmetry_coeff = 1, 
     & 2 - KRONEKER(l1_t,l2_t)*KRONEKER(abs(nju1_t),abs(nju2_t))
      if(symmetry_coeff.eq.1) then
	  sign_correction = 1 
      l1_t = A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
      else if(symmetry_coeff.eq.2) then
      sign_correction = (-1)**(l1_t+l2_t)	  
      l2_t = A_TOP(1,ind_t)
      nju2_t = A_TOP(2,ind_t)
      l1_t = A_TOP(3,ind_t)
      nju1_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  else
	  if(myid.eq.0 .and. i_r_point.eq.1) 
     & print*, "Error in identical swapping variable."
	  stop
      end if
	  
      DO planr_coeff = 1, 2 - KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	 
      IF(planr_coeff.eq.2) THEN
      sign_correction = sign_correction
     & *(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t) 
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      ENDIF
	  	  
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
      CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
      CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
      CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)
	  
	  DO k1p = -j1_p_t,j1_p_t
      k1pp = k1p - nju1_t	  
      IF(abs(k1pp).gt.j1_pp_t) CYCLE	 	  
      CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,k1pp,nju1_t,k1p,
     & bikram_w3j_fact)
      coeff_1_p = M1_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M2_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF	  
      DO k2p =-j2_p_t,j2_p_t
      k2pp = k2p - nju2_t
      IF(abs(k2pp).gt.j2_pp_t) CYCLE	  
      coeff_2_p = M2_VECTORS(channp,j2_p_t+1+k2p)
      CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,k2pp,nju2_t,k2p,
     & bikram_w3j_fact)
      IF(abs(k2pp).le.j2_pp_t) THEN
      coeff_2_pp = M1_VECTORS(channpp,j2_pp_t+1+k2pp)
      ELSE
      coeff_2_pp = 0d0       	  
      ENDIF

	  M_coulp_ident = M_coulp_ident + coef1*coef2*exp_coeff_int*coef3*
     & CG_j1_j2_p * CG_j1_j2_pp * CG_l1_l2 *
     & CG_j1_l1 * CG_j2_l2 * CG_j1_k1 * CG_j2_k2 *
     & coeff_1_p*coeff_2_p*coeff_1_pp*coeff_2_pp*sign_correction
	  
	  if(bikram_rebalance) tmp1 = tmp1 + 1
      ENDDO	  
      ENDDO	  
      ENDDO
      ENDDO
      ENDDO	  
      ENDDO
      ENDDO
	  
	  else
	  print*,"Something is wrong in PES expansion matrix computation."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if

      M_coulp_ident_2 = M_coulp_ident 
	  
      sign_correction = +1
	  if(bk_par1*bk_par2.lt.0.d0) then
	  if(M_coulp_ident_1*M_coulp_ident_2.gt.0.d0) sign_correction = -1
	  else if(bk_par1*bk_par2.gt.0.d0) then
	  if(M_coulp_ident_1*M_coulp_ident_2.lt.0.d0) sign_correction = -1
	  end if
	  M_coulp = M_coulp + M_coulp_ident*par_pp
     & *bk_par2*sign_correction
! Bikram End.
	  

      M_coulp_ident = 0d0
      buff = 0d0	  
      j1_p_t = j2_ch(channp)
      j2_p_t = j1_ch(channp)
      j1_pp_t = j2_ch(channpp)
      j2_pp_t = j1_ch(channpp)	  

	  if(.not.bikram_identical_pes) then
	  
	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
	  DO i=1,nterms
      ind_t =i
      exp_coeff_int=expansion_terms(i_r_point,ind_t)	 
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(8d0*pi**2)
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)
	  DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
      CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
      CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
      CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
      CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)
	  
      DO planr_coeff =1,2-KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	 
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int=exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t)	  
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      ENDIF
	  
	  DO k1p = -j1_p_t,j1_p_t
      k1pp = k1p - nju1_t	  
      IF(abs(k1pp).gt.j1_pp_t) CYCLE	 	  
      CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,k1pp,nju1_t,k1p,
     & bikram_w3j_fact)
      coeff_1_p = M2_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M2_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF	  
      DO k2p =-j2_p_t,j2_p_t
      k2pp = k2p - nju2_t
      IF(abs(k2pp).gt.j2_pp_t) CYCLE	  
      coeff_2_p = M1_VECTORS(channp,j2_p_t+1+k2p)
      CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,k2pp,nju2_t,k2p,
     & bikram_w3j_fact)
      IF(abs(k2pp).le.j2_pp_t) THEN
      coeff_2_pp = M1_VECTORS(channpp,j2_pp_t+1+k2pp)
      ELSE
      coeff_2_pp = 0d0       	  
      ENDIF

	  M_coulp_ident = M_coulp_ident + coef1*coef2*exp_coeff_int*coef3*
     & CG_j1_j2_p * CG_j1_j2_pp * CG_l1_l2 *
     & CG_j1_l1 * CG_j2_l2 * CG_j1_k1 * CG_j2_k2 *
     & coeff_1_p*coeff_2_p*coeff_1_pp*coeff_2_pp

      ENDDO	  
      ENDDO
      ENDDO
      ENDDO	  
      ENDDO
      ENDDO	  
	  
	  else if(bikram_identical_pes) then
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to construct matrix for identical molecule-molecule collision using expansion
! an example would be H2O-H2O, honestly, this derivation was done using water collision
! for more details follow MQCT 2022 manual
! Bikram Start: Aug 2022

	  coef1 = dsqrt((2d0*j1_pp_t + 1d0)/(2d0*j1_p_t + 1d0))
	  coef2 = dsqrt((2d0*j2_pp_t + 1d0)/(2d0*j2_p_t + 1d0))
	  
      DO i = 1, nterms
      ind_t = i
      exp_coeff_int = expansion_terms(i_r_point,ind_t)	 
      l1_t = A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  	  
	  coef3 = dsqrt((2d0*l1_t + 1d0)*(2d0*l2_t + 1d0))/(8d0*pi**2)
	  sign_correction = 1 
	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG_bikram(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t,
     & bikram_w3j_fact)

	  DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      CG_j1_j2_pp = CG_bikram(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,
     & m12_t,bikram_w3j_fact)
	  
      DO symmetry_coeff = 1, 
     & 2 - KRONEKER(l1_t,l2_t)*KRONEKER(abs(nju1_t),abs(nju2_t))
      if(symmetry_coeff.eq.1) then
	  sign_correction = 1 
      l1_t = A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
      else if(symmetry_coeff.eq.2) then
      sign_correction = (-1)**(l1_t+l2_t)	  
      l2_t = A_TOP(1,ind_t)
      nju2_t = A_TOP(2,ind_t)
      l1_t = A_TOP(3,ind_t)
      nju1_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
	  else
	  if(myid.eq.0 .and. i_r_point.eq.1) 
     & print*, "Error in identical swapping variable."
	  stop
      end if
	  
      DO planr_coeff = 1, 2 - KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	 
      IF(planr_coeff.eq.2) THEN
      sign_correction = sign_correction
     & *(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t) 
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      ENDIF
	  	  
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) CYCLE
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) CYCLE
      CG_l1_l2 = CG_bikram(l1_t,l2_t,l_t,m_exp,-m_exp,0,
     & bikram_w3j_fact)
      CG_j1_l1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp,
     & bikram_w3j_fact)
      CG_j2_l2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,
     & m12_t-mp,bikram_w3j_fact)
	  
	  DO k1p = -j1_p_t,j1_p_t
      k1pp = k1p - nju1_t	  
      IF(abs(k1pp).gt.j1_pp_t) CYCLE	 	  
      CG_j1_k1 = CG_bikram(j1_pp_t,l1_t,j1_p_t,k1pp,nju1_t,k1p,
     & bikram_w3j_fact)
      coeff_1_p = M2_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M2_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF	  
      DO k2p =-j2_p_t,j2_p_t
      k2pp = k2p - nju2_t
      IF(abs(k2pp).gt.j2_pp_t) CYCLE	  
      coeff_2_p = M1_VECTORS(channp,j2_p_t+1+k2p)
      CG_j2_k2 = CG_bikram(j2_pp_t,l2_t,j2_p_t,k2pp,nju2_t,k2p,
     & bikram_w3j_fact)
      IF(abs(k2pp).le.j2_pp_t) THEN
      coeff_2_pp = M1_VECTORS(channpp,j2_pp_t+1+k2pp)
      ELSE
      coeff_2_pp = 0d0       	  
      ENDIF

	  M_coulp_ident = M_coulp_ident + coef1*coef2*exp_coeff_int*coef3*
     & CG_j1_j2_p * CG_j1_j2_pp * CG_l1_l2 *
     & CG_j1_l1 * CG_j2_l2 * CG_j1_k1 * CG_j2_k2 *
     & coeff_1_p*coeff_2_p*coeff_1_pp*coeff_2_pp*sign_correction
	  
	  if(bikram_rebalance) tmp1 = tmp1 + 1
      ENDDO	  
      ENDDO	  
      ENDDO
      ENDDO
      ENDDO	  
      ENDDO
      ENDDO
	  
	  else
	  print*,"Something is wrong in PES expansion matrix computation."
      CALL MPI_Abort(MPI_COMM_WORLD,ierr_mpi,ierr_mpi)
	  stop
	  end if

      M_coulp_ident_3 = M_coulp_ident

	  sign_correction = +1
	  if(bk_par1*bk_par2.lt.0.d0) then
	  if(M_coulp_ident_3*M_coulp_non_ident.gt.0.d0) 
     & sign_correction = -1
	  else if(bk_par1*bk_par2.gt.0.d0) then
	  if(M_coulp_ident_3*M_coulp_non_ident.lt.0.d0) 
     & sign_correction = -1
	  end if
	  M_coulp = M_coulp + M_coulp_ident*par_pp*par_p
     & *bk_par1*bk_par2*sign_correction

      M_coulp = M_coulp/
     & dsqrt(2d0*(1d0+delta(j1_p_t,j2_p_t)))
     & /dsqrt(2d0*(1d0+delta(j1_pp_t,j2_pp_t)))
      ENDIF	 	  
	  
      END SELECT 	  
      END SUBROUTINE EXPANSION_MATRIX_ELEMENT

      LOGICAL FUNCTION TRIANG_RULE(a,b,c)
      IMPLICIT NONE
      INTEGER a,b,c	  
      TRIANG_RULE = .TRUE.	  
      IF(2*max(a,b,c).gt.(a+b+c)) TRIANG_RULE = .FALSE.	  
      END FUNCTION TRIANG_RULE	  
     
      INCLUDE "basis_wave_functions.f"