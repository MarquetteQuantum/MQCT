	  subroutine resize
	  use variables
	  use mpi_data
	  use MPI_TASK_TRAJECT
	  integer i,j,k,ii
      integer j_12,m_12,m1,m0	  
	  real*8 del_E_tmp,cntr_tmp
	  real*8,allocatable :: bk_del_e(:)
	  logical same
	  
	  mat_sz_bk = portion_of_MIJ_per_task(2,myid+1) - 
     & portion_of_MIJ_per_task(1,myid+1) +1
	  
	  if(allocated(ind_mat_bk)) deallocate(ind_mat_bk)
	  allocate(ind_mat_bk(2,mat_sz_bk))
	  
	  if(mpi_task_per_traject.eq.1) then
	  ind_mat_bk(:,:) = ind_mat(:,:)
	  else
	  if(bikram_mij_multiprint) then
	  ind_mat_bk(:,:) = ind_mat(:,:)
	  else
	  do i = 1, mat_sz_bk
	  ind_mat_bk(:,i) = ind_mat(:,
     & (i+portion_of_MIJ_per_task(1,myid+1)-1))
	  end do
	  end if
	  end if
	  
      if(allocated(misc_bk)) deallocate(misc_bk)
	  allocate(misc_bk(4,mat_sz_bk),bk_del_e(mat_sz_bk))
	  
	  do k = 1, mat_sz_bk
	  i = ind_mat_bk(1,k)
      j = ind_mat_bk(2,k)
	  
	  j_12 = j12(j)  
      m_12 = m12(i) 
      m1 = m_12+1d0
      m0 = m_12-1d0	  
	  misc_bk(1,k) = sqrt(j_12*(j_12+1d0)-m1*(m1-1d0))
	  misc_bk(2,k) = sqrt(j_12*(j_12+1d0)-m0*(m0+1d0))
	  
	  j_12 = j12(i)  
      m_12 = m12(j)  
      m1 = m_12+1
      m0 = m_12-1
	  misc_bk(3,k) = sqrt(j_12*(j_12+1d0)-m1*(m1-1d0))
	  misc_bk(4,k) = sqrt(j_12*(j_12+1d0)-m0*(m0+1d0))
	  
	  enddo
	  
	  if(allocated(bk_indx)) deallocate(bk_indx)
	  allocate(bk_indx(mat_sz_bk))
	  
	  ph_cntr_bk = 0
	  bk_del_e(1) = Ej(ind_mat_bk(2,1))-Ej(ind_mat_bk(1,1))
	  ph_cntr_bk = ph_cntr_bk + 1
	  bk_indx(1) = ph_cntr_bk
	  cntr_tmp = ph_cntr_bk
	  
      do k = 2, mat_sz_bk
	  i = ind_mat_bk(1,k)
      j = ind_mat_bk(2,k)
	  del_E_tmp = Ej(j)-Ej(i)
	  
	  same = .false.
	  do ii=1,cntr_tmp
	  if(abs(del_E_tmp - bk_del_e(ii)).lt.1.E-12) then
	  same = .true.
	  bk_indx(k) = ii
	  
	  endif
	  enddo
	  
	  if(same) then
	  cycle
	  else
	  ph_cntr_bk = ph_cntr_bk + 1
	  bk_indx(k) = ph_cntr_bk
	  bk_del_e(ph_cntr_bk) = del_E_tmp
	  cntr_tmp = ph_cntr_bk
	  endif
	  enddo
	  
	  if(allocated(bk_delta_E)) deallocate(bk_delta_E)
	  allocate(bk_delta_E(ph_cntr_bk))
	  bk_delta_E(1:ph_cntr_bk) = bk_del_e(1:ph_cntr_bk)
	  
      if(.not. allocated(bk_splint_bfr_mat)) 
     & allocate(bk_splint_bfr_mat(mat_sz_bk))
      if(.not. allocated(bk_splint_afr_mat)) 
     & allocate(bk_splint_afr_mat(mat_sz_bk))
      if(.not. allocated(bk_splint_bfr_dmat)) 
     & allocate(bk_splint_bfr_dmat(mat_sz_bk))
      if(.not. allocated(bk_splint_afr_dmat)) 
     & allocate(bk_splint_afr_dmat(mat_sz_bk))
	  
	  if(allocated(bk_mat_splnt)) deallocate(bk_mat_splnt)
	  if(allocated(bk_dmat_splnt)) deallocate(bk_dmat_splnt)
      allocate(bk_mat_splnt(mat_sz_bk,n_r_coll))
      allocate(bk_dmat_splnt(mat_sz_bk,n_r_coll))
	  
	  bk_mat_splnt = transpose(Mat_el)
	  bk_dmat_splnt = transpose(Mat_el_der)
	  
	  splnt_bfr = -1
	  splnt_afr = -1
	  
	  end subroutine 

      SUBROUTINE DERIVS_BK(t,q,dqdt) 										!Bikram: my derivs routine for all in one
      USE VARIABLES
      USE MPI_DATA
      USE MPI
      USE INTEGRATOR	  
      USE MPI_TASK_TRAJECT	  
      USE CONSTANTS
      USE OLD_MIJ	    
      IMPLICIT NONE
      LOGICAL BELONGS	  
      INTEGER status(MPI_STATUS_SIZE)	  
      REAL*8 dqdt(states_size*2+8), q(states_size*2+8)
      REAL*8 rphase,rfact_phase,iphase,ifact_phase
      REAL*8 d_f,tet,phi,sin_tet,cos_tet
      REAL*8 t,delta,Rrr,dPhidT,dThetadT
      REAL*8 dq_dt_buffer(states_size*2+8)	   
      INTEGER i,j,k,s1,s0
      INTEGER s11,s00
      REAL*8 j_12,m_12,m1,m0	  
      INTEGER i_root,dest,origin,round	  
      REAL*8 aver_poten_temp_0, der_aver_poten_temp_0
! Bikram Start Dec 2019
	  integer tmp_indx
	  REAL*8,allocatable :: bk_mat(:),bk_der_mat(:)
	  real*8 sqrt1,sqrt2						
	  real*8,allocatable :: bk_amp_phi(:),bk_amp_theta(:)
	  real*8,allocatable :: bk_sin_cos(:,:)			
! Bikram End.
      EXTERNAL delta,BELONGS,round


!--------------------------------------------------------------------
! Store the values of R and Theta and their derivatives as local variables
!--------------------------------------------------------------------
      dqdt=0d0
      Rrr = q(2*states_size+1)
      tet = q(2*states_size+3)
      sin_tet = dsin(tet)
      cos_tet = dcos(tet)
! Bikram Start Dec 2019:
	  dThetadT = 0.d0
	  if(bikram_theta) then
      dThetadT = q(2*states_size+4)
     % /q(2*states_size+1)/q(2*states_size+1)/massAB_M
	  endif
! Bikram End.
      dPhidT = q(2*states_size+6)
     % /q(2*states_size+1)/q(2*states_size+1)/massAB_M/
     & sin_tet**2	  
      rphase = 1d0
      iphase = 0d0
!--------------------------------------------------------------------
! Allocate a working array, bk_mat, for the RHS calculations
!--------------------------------------------------------------------
! Bikram Start Dec 2019:
	  allocate(bk_mat(mat_sz_bk),bk_der_mat(mat_sz_bk))
	  allocate(bk_amp_phi(mat_sz_bk),bk_amp_theta(mat_sz_bk))
	  bk_amp_phi = 0.d0
	  if(bikram_theta) bk_amp_theta = 0.d0
!--------------------------------------------------------------------
! Computes all Matrix elements, and their derivative at this value of R
! and stores in the working array
!--------------------------------------------------------------------
	  call potential_data_bk(Rrr, bk_mat, bk_der_mat)
! Bikram End.
!--------------------------------------------------------------------
! Computes phase factors
!--------------------------------------------------------------------
	  allocate(bk_sin_cos(2,ph_cntr_bk))
	  do k= 1,ph_cntr_bk
	  bk_sin_cos(1,k) = dsin(t*bk_delta_E(k))
	  bk_sin_cos(2,k) = dcos(t*bk_delta_E(k))
	  enddo
!--------------------------------------------------------------------
! Loop over the matrix elements for each processor to prepare another 
! working array, bk_amp_phi and bk_amp_theta, for the RHS calculations
!--------------------------------------------------------------------
      DO k = 1, mat_sz_bk
      i = ind_mat_bk(1,k)
      j = ind_mat_bk(2,k)  
!--------------------------------------------------------------------
! Skipping if parity of the state is different from the initial state parity
!--------------------------------------------------------------------	  
	  if(identical_particles_defined .and. ini_st_ident) then
	  if(((-1)**(bk_parity-1)).ne.parity_state(i)) cycle
	  end if
!--------------------------------------------------------------------
! Skipping if coupled state approx defined
!--------------------------------------------------------------------
	  if(coupled_states_defined) then
	  if(m12(s_st).ne.m12(i)) cycle
	  endif
!--------------------------------------------------------------------
! Skipping states indicated in the input file
!--------------------------------------------------------------------  
      IF(states_to_exclude_defined) THEN
      IF(stts_to_excl(i) .or. stts_to_excl(j)) CYCLE
      ENDIF
!--------------------------------------------------------------------
! Store matrix element and its derivative as local variables
!--------------------------------------------------------------------	  
      aver_poten_temp_0 = bk_mat(k)							
	  der_aver_poten_temp_0 = bk_der_mat(k)
!--------------------------------------------------------------------
! Computing phase factor for this transition
!--------------------------------------------------------------------	  
	  d_f = 2d0
      if(i.eq.j) d_f = 1d0
	  tmp_indx = bk_indx(k)
      iphase = bk_sin_cos(1,tmp_indx)
      rphase = bk_sin_cos(2,tmp_indx)
      rfact_phase = -d_f*iphase
      ifact_phase = d_f*rphase
!--------------------------------------------------------------------
! Computing R.H.S. of quantum equations of motion by accumulation
!--------------------------------------------------------------------
       dqdt(i) = dqdt(i)
     & + (q(j)*rfact_phase + q(j+states_size)*ifact_phase)	 
     &   *aver_poten_temp_0/2d0
	 
      dqdt(i+states_size) = dqdt(i+states_size)
     & + (q(j+states_size)*rfact_phase - q(j)*ifact_phase)	 
     &   *aver_poten_temp_0/2d0
	 
      dqdt(j) = dqdt(j)
     & - (q(i)*rfact_phase - q(i+states_size)*ifact_phase)	 
     &   *aver_poten_temp_0/2d0
	 
      dqdt(j+states_size) = dqdt(j+states_size)
     & - (q(i+states_size)*rfact_phase + q(i)*ifact_phase)	 
     &   *aver_poten_temp_0/2d0
!--------------------------------------------------------------------
! Computing R.H.S. of classical equation for P_R by accumulation
!--------------------------------------------------------------------
      dqdt(states_size*2+2) = dqdt(states_size*2+2)
     ^ - ((q(i)*q(j) + q(i+states_size)*q(j+states_size))*rphase
     ^ - (q(i+states_size)*q(j) - q(i)*q(j+states_size))*iphase)
     ^   *d_f*der_aver_poten_temp_0
!--------------------------------------------------------------------
! Coriolis Coupling
!--------------------------------------------------------------------
	  if(.not. coupled_states_defined) then										!Skipping if coupled state approx defined
	  
      j_12 = j12(j)  
      m_12 = m12(i) 
      m1 = m_12+1d0
      m0 = m_12-1d0	
	  
      IF(int(abs(m1)).le.int(j_12))  THEN 
      s1 =  j+1
	  s11 = s1 +states_size 
      ELSE 
	  s1 = 2*states_size+7
	  s11 = 2*states_size+8	   
	  ENDIF
	  IF(int(abs(m0)).le.int(j_12)) then
	  s0 = j-1
	  s00 = s0+states_size	 
      ELSE
      s0 = 2*states_size+7
      s00 = 2*states_size+8	
      ENDIF
!--------------------------------------------------------------------
! To compute Coriolis coupling term for P_Phi and P_Theta and store 
! in the working array bk_amp_phi and bk_amp_theta
!--------------------------------------------------------------------	 
      bk_amp_phi(k) = bk_amp_phi(k)
     & - sin_tet*rfact_phase/2d0*
     &     (misc_bk(1,k)*(q(i)*q(s1) + q(i+states_size)*q(s11)) + 
     & 	      misc_bk(2,k)*(q(i)*q(s0) + q(i+states_size)*q(s00)))
     & - sin_tet*(-ifact_phase)/2d0*
     &     (misc_bk(1,k)*(q(i+states_size)*q(s1) - q(i)*q(s11)) + 	 
     & 	      misc_bk(2,k)*(q(i+states_size)*q(s0) - q(i)*q(s00)))
! Bikram Start Dec 2019:
	  if(bikram_theta) then
	  bk_amp_theta(k) = bk_amp_theta(k)
     & + d_f*rphase/2d0*
     &     (misc_bk(1,k)*(q(i)*q(s1) + q(i+states_size)*q(s11)) - 
     & 	      misc_bk(2,k)*(q(i)*q(s0) + q(i+states_size)*q(s00)))
     & + d_f*(-iphase)/2d0*
     &     (misc_bk(1,k)*(q(i+states_size)*q(s1) - q(i)*q(s11)) - 	 
     & 	      misc_bk(2,k)*(q(i+states_size)*q(s0) - q(i)*q(s00)))
	  endif
! Bikram End.
	 
      j_12 = j12(i)  
      m_12 = m12(j)  
      m1 = m_12+1
      m0 = m_12-1
	  
      IF(int(abs(m1)).le.int(j_12)) THEN 
      s1 =  i+1
      s11 = s1 +states_size
      ELSE 
      s1 = 2*states_size+7
      s11 = 2*states_size+8
      ENDIF
      IF(int(abs(m0)).le.int(j_12)) THEN
      s0 = i-1 
      s00 = s0+states_size
      ELSE
      s0 = states_size*2+7
      s00 = states_size*2+8
      ENDIF
!--------------------------------------------------------------------
! To compute rest of the RHS of equation for P_Phi and P_Theta and store 
! in the working array bk_amp_phi and bk_amp_theta
!--------------------------------------------------------------------	 	
      bk_amp_phi(k) = bk_amp_phi(k)
     & - sin_tet*(-rfact_phase)/2d0*	 
     &     (misc_bk(3,k)*(q(s1)*q(j) + q(s11)*q(j+states_size)) +	 
     & 	      misc_bk(4,k)*(q(s0)*q(j) + q(s00)*q(j+states_size)))
	 
      bk_amp_phi(k) = bk_amp_phi(k)
     & - sin_tet*(-ifact_phase)/2d0*	 
     &     (misc_bk(3,k)*(q(s1)*q(j+states_size) - q(s11)*q(j)) + 	 
     & 	      misc_bk(4,k)*(q(s0)*q(j+states_size) - q(s00)*q(j)))
! Bikram Start Dec 2019:
	  if(bikram_theta) then
	  bk_amp_theta(k) = bk_amp_theta(k)
     & + d_f*(rphase)/2d0*	 
     &     (misc_bk(3,k)*(q(s1)*q(j) + q(s11)*q(j+states_size)) -	 
     & 	      misc_bk(4,k)*(q(s0)*q(j) + q(s00)*q(j+states_size)))
	  
	  bk_amp_theta(k) = bk_amp_theta(k)
     & + d_f*(iphase)/2d0*	 
     &     (misc_bk(3,k)*(q(s1)*q(j+states_size) - q(s11)*q(j)) - 	 
     & 	      misc_bk(4,k)*(q(s0)*q(j+states_size) - q(s00)*q(j)))
	  endif
! Bikram End.
	  endif																		!finishing if loop for coupled state approx.
      ENDDO
	  
	  if(.not. coupled_states_defined) then										!Skipping if coupled state approx defined
!--------------------------------------------------------------------
! RHS as a scalar product of two working arrays that 
! takes care of triple sum over n, n', and m
!--------------------------------------------------------------------
	  dqdt(states_size*2+6) = dot_product(bk_mat, bk_amp_phi)				
	  if(bikram_theta) then
	  dqdt(states_size*2+4) = dot_product(bk_mat, bk_amp_theta)			
	  endif
	  endif																		!finishing if loop for coupled state approx.
   
      DO i=1,traject_roots	  
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	  
      origin = 0
      CALL MPI_Reduce(dqdt, dq_dt_buffer, 2*states_size+8, MPI_REAL8,
     & MPI_SUM, origin,	 
     & comms(i),ierr_mpi)
      dqdt  = dq_dt_buffer
      CALL MPI_Bcast(dqdt,2*states_size+8,
     & MPI_REAL8, origin,comms(i),ierr_mpi)
      ENDIF
      ENDDO	  
	   
	  if(.not. coupled_states_defined) then										!Skipping if coupled state approx defined
	  
	  if(bikram_theta) then														!for spin-rot calculation
	  
	  do i=1,states_size
	  
	  j_12 = j12(i)  
      m_12 = m12(i) 
      m1 = m_12+1d0
      m0 = m_12-1d0	  
	  
	  sqrt1 = sqrt(j_12*(j_12+1d0)-m0*(m0+1d0))					
	  sqrt2 = sqrt(j_12*(j_12+1d0)-m1*(m1-1d0))			
	  
      dqdt(i) = dqdt(i)
     & + dPhidT*m_12*cos_tet*q(i)

      dqdt(i+states_size) = dqdt(i+states_size)
     & - dPhidT*m_12*cos_tet*q(i+states_size)		
	  
      if(i.gt.1) then 
      dqdt(i) = dqdt(i)
     & + (dPhidT*sin_tet*q(i-1+states_size)
     & - dThetadT*q(i-1))*sqrt1/2d0

      dqdt(i+states_size) = dqdt(i+states_size)
     & - (dPhidT*sin_tet*q(i-1)
     & - dThetadT*q(i-1+states_size))*sqrt1/2d0
      endif 		  
	  
	  if(i.lt.states_size) then 
      dqdt(i) = dqdt(i)
     & + (dPhidT*sin_tet*q(i+1+states_size)
     & - dThetadT*q(i+1))*sqrt2/2d0

      dqdt(i+states_size) = dqdt(i+states_size)
     & - (dPhidT*sin_tet*q(i+1)
     & - dThetadT*q(i+1+states_size))*sqrt2/2d0
      endif
      ENDDO
	  
	  else
	  
	  do i=1,states_size
	  
	  j_12 = j12(i)  
      m_12 = m12(i) 
      m1 = m_12+1d0
      m0 = m_12-1d0	  
	  
	  sqrt1 = sqrt(j_12*(j_12+1d0)-m0*(m0+1d0))					
	  sqrt2 = sqrt(j_12*(j_12+1d0)-m1*(m1-1d0))					
	  
      if(i.gt.1) then 
      dqdt(i) = dqdt(i)
     & + dPhidT/2d0*sqrt1*q(i-1+states_size)

      dqdt(i+states_size) = dqdt(i+states_size)
     & - dPhidT/2d0*sqrt1*q(i-1)
      endif 		  
	  
	  if(i.lt.states_size) then 
      dqdt(i) = dqdt(i)
     & + dPhidT/2d0*sqrt2*q(i+1+states_size)

      dqdt(i+states_size) = dqdt(i+states_size)
     & - dPhidT/2d0*sqrt2*q(i+1)
      endif
      ENDDO
	  
	  endif
	  endif																		!finishing if loop for coupled state approx.
!--------------------------------------------------------------------
! The RHS of equations for R, Phi, and Theta
! and the remaining pieces of equations P_Theta and P_Phi
!--------------------------------------------------------------------	  
      dqdt(2*states_size+1) = q(2*states_size+2)/massAB_M						! For R
      dqdt(2*states_size+5) = q(2*states_size+6)								! For Phi
     % /q(2*states_size+1)/q(2*states_size+1)/massAB_M/
     & sin_tet**2	 

      dqdt(states_size*2+2)= dqdt(states_size*2+2)	   							! For P_R
     ^ +q(states_size*2+4)/q(states_size*2+1)
     & /q(states_size*2+1)/massAB_M*
     & q(states_size*2+4)/q(states_size*2+1)+
     ^ q(states_size*2+6)/q(states_size*2+1)
     & /q(states_size*2+1)/massAB_M*
     & q(states_size*2+6)/q(states_size*2+1)/sin_tet**2
! Bikram Start Dec 2019: 
	  if(bikram_theta) then														! For Theta
      dqdt(2*states_size+3) = q(2*states_size+4)
     % /q(2*states_size+1)/q(2*states_size+1)/massAB_M
	 
       dqdt(states_size*2+4)=dqdt(states_size*2+4)								! For P_Theta
     ^ +q(states_size*2+6)/q(states_size*2+1)
     & /q(states_size*2+1)/massAB_M*
     & q(states_size*2+6)
     ^ /sin_tet**3*cos_tet
	  endif
! Bikram End.
	  
      RETURN	 
      END SUBROUTINE DERIVS_BK

      SUBROUTINE potential_data_bk(R,bk_mat,bk_der_mat)
      USE VARIABLES
      USE GRID_INTEGR
      USE CONSTANTS
      USE MPI_TASK_TRAJECT 
      USE MPI_DATA
      USE INTEGRATOR	  
      USE OLD_Mij	  
      IMPLICIT NONE	 
      REAL*8 x,v_m,der_v_m,R,V_COULPING_TERMS
      INTEGER k,k_real,i_exp_term
	  real*8 bk_mat(mat_sz_bk), bk_der_mat(mat_sz_bk)

      x = R
      v_m = 0d0	  
      der_v_m = 0d0
	  
      IF(term_pot_defined) THEN
	  do k =1,mat_sz_bk
      CALL COMPUTE_V_TERMS(R)
      DO i_exp_term=1,nterms
      CALL TERM_PES_MATRIX_ELEMENT(V_COULPING_TERMS,k,i_exp_term)	  
      v_m = v_m	+
     &	V_COULPING_TERMS*V_TERMS(i_exp_term)
      der_v_m = der_v_m	+
     &	V_COULPING_TERMS*V_TERMS_DER(i_exp_term) 	 
      ENDDO		  
	  bk_mat(k)=v_m
	  bk_der_mat(k)=der_v_m
	  enddo
      RETURN	  
      ENDIF	  
	  
      IF(x.lt.R_COM(1)) THEN
	  do k =1,mat_sz_bk
      CALL 	OUT_OF_RANGE(v_m,der_v_m,R,k)  	  
	  bk_mat(k)=v_m
	  bk_der_mat(k)=der_v_m
	  enddo
      RETURN
      ENDIF	  
	  
      IF(x.ge.R_COM(n_r_coll)) then
	  bk_mat(:) = 0.d0
	  bk_der_mat(:) = 0.d0
	  RETURN
	  endif
	  
	  CALL splint_bk(x,bk_mat,bk_der_mat)  
	  
      END SUBROUTINE potential_data_bk   

      SUBROUTINE splint_bk(x,bk_mat,bk_der_mat)
      USE VARIABLES
	  implicit none
      REAL*8 x,y,xaaa(n_r_coll),y2a(n_r_coll),ya(n_r_coll),y_prime
	  integer bk,k_real,bkk
      INTEGER k,khi,klo
	  real*8 chk
      REAL*8 aa,bb,h
	  real*8 cc1,cc2,cc3,cc4,cc5,cc6
	  real*8 m1,m2,dm1,dm2
	  real*8 bk_mat(mat_sz_bk), bk_der_mat(mat_sz_bk)
	  
      klo=1
      khi=n_r_coll
1       if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(R_COM(k).gt.x)then
      khi=k 
      else
      klo=k
      endif
      goto 1
      endif
      h=R_COM(khi)-R_COM(klo)
      if (h.eq.0.) then
      WRITE(*,*) R_COM(n_r_coll),x,n_r_coll
      WRITE(*,*) "DISTANCE",R_COM	  
      STOP 'bad xa input in matrix splint' 
      endif
      aa=(R_COM(khi)-x)/h
      bb=(x-R_COM(klo))/h
	  
	  cc1=aa**3-aa
	  cc2=bb**3-bb
	  cc3=3*aa**2-1d0
	  cc4=3*bb**2-1d0
	  cc5=(h**2)/6.
	  cc6=cc5/h
	  if(splnt_bfr.eq.khi .and. splnt_afr.eq.klo) then
	  else
	  bk_splint_bfr_mat(:) = bk_mat_splnt(:,klo)
	  bk_splint_afr_mat(:) = bk_mat_splnt(:,khi)
	  bk_splint_bfr_dmat(:) = bk_dmat_splnt(:,klo)
	  bk_splint_afr_dmat(:) = bk_dmat_splnt(:,khi)
	  endif
	  
	  bk_mat = aa*bk_splint_bfr_mat + bb*bk_splint_afr_mat 
     &  + (cc1*bk_splint_bfr_dmat + cc2*bk_splint_afr_dmat)*cc5
	  
	  bk_der_mat = (bk_splint_afr_mat - bk_splint_bfr_mat)/h
     & + (-cc3*bk_splint_bfr_dmat + cc4*bk_splint_afr_dmat)*cc6
	  
	  splnt_bfr = khi
	  splnt_afr = klo
	  
      return
      END

	  subroutine resize_adia												!Bikram: my derivs routine for adiabatic trejectories
	  use variables
	  use mpi_data
	  use MPI_TASK_TRAJECT
	  integer i,j,k,ii
	  real*8 del_E_tmp,cntr_tmp
	  real*8,allocatable :: bk_del_e(:)
	  logical same
	  
	  mat_sz_bk = portion_of_MIJ_per_task(2,myid+1) - 
     & portion_of_MIJ_per_task(1,myid+1) +1
	  
!	  if(allocated(ind_mat_bk)) deallocate(ind_mat_bk)
	  allocate(ind_mat_bk(2,mat_sz_bk))
	  
	  if(mpi_task_per_traject.eq.1) then
	  ind_mat_bk(:,:) = ind_mat(:,:)
	  else
	  if(bikram_mij_multiprint) then
	  ind_mat_bk(:,:) = ind_mat(:,:)
	  else
	  do i = 1, mat_sz_bk
	  ind_mat_bk(:,i) = ind_mat(:,
     & (i+portion_of_MIJ_per_task(1,myid+1)-1))
	  end do
	  end if
	  end if
	  
!	  if(allocated(bk_indx)) deallocate(bk_indx)
	  allocate(bk_indx(mat_sz_bk), bk_del_e(mat_sz_bk))
	  
	  ph_cntr_bk = 0
	  bk_del_e(1) = Ej(ind_mat_bk(2,1))-Ej(ind_mat_bk(1,1))
	  ph_cntr_bk = ph_cntr_bk + 1
	  bk_indx(1) = ph_cntr_bk
	  cntr_tmp = ph_cntr_bk
	  
      do k = 2, mat_sz_bk
	  i = ind_mat_bk(1,k)
      j = ind_mat_bk(2,k)
	  del_E_tmp = Ej(j)-Ej(i)
	  
	  same = .false.
	  do ii=1,cntr_tmp
	  if(abs(del_E_tmp - bk_del_e(ii)).lt.1.E-12) then
	  same = .true.
	  bk_indx(k) = ii
	  
	  endif
	  enddo
	  
	  if(same) then
	  cycle
	  else
	  ph_cntr_bk = ph_cntr_bk + 1
	  bk_indx(k) = ph_cntr_bk
	  bk_del_e(ph_cntr_bk) = del_E_tmp
	  cntr_tmp = ph_cntr_bk
	  endif
	  enddo
	  
!	  if(allocated(bk_delta_E)) deallocate(bk_delta_E)
	  allocate(bk_delta_E(ph_cntr_bk))
	  bk_delta_E(1:ph_cntr_bk) = bk_del_e(1:ph_cntr_bk)
	  
!	  if(allocated(bk_splint_bfr_mat)) deallocate(bk_splint_bfr_mat)
!	  if(allocated(bk_splint_afr_mat)) deallocate(bk_splint_afr_mat)
!	  if(allocated(bk_splint_bfr_dmat)) deallocate(bk_splint_bfr_dmat)
!	  if(allocated(bk_splint_afr_dmat)) deallocate(bk_splint_afr_dmat)
      allocate(bk_splint_bfr_mat(mat_sz_bk))
      allocate(bk_splint_afr_mat(mat_sz_bk))
      allocate(bk_splint_bfr_dmat(mat_sz_bk))
      allocate(bk_splint_afr_dmat(mat_sz_bk))
	  
!	  if(allocated(bk_mat_splnt)) deallocate(bk_mat_splnt)
!	  if(allocated(bk_dmat_splnt)) deallocate(bk_dmat_splnt)
      allocate(bk_mat_splnt(mat_sz_bk,n_r_coll))
      allocate(bk_dmat_splnt(mat_sz_bk,n_r_coll))
	  bk_mat_splnt = transpose(Mat_el)
	  bk_dmat_splnt = transpose(Mat_el_der)
	  
	  splnt_bfr = -1
	  splnt_afr = -1
	  
	  end subroutine 
	  
	  subroutine resize_adia_dealloc										!Bikram: my derivs routine for adiabatic trejectories
	  use variables
	  use mpi_data
	  use MPI_TASK_TRAJECT

	  if(allocated(ind_mat_bk)) deallocate(ind_mat_bk)

	  if(allocated(bk_indx)) deallocate(bk_indx)
	  
	  if(allocated(bk_delta_E)) deallocate(bk_delta_E)
	  
	  if(allocated(bk_splint_bfr_mat)) deallocate(bk_splint_bfr_mat)
	  if(allocated(bk_splint_afr_mat)) deallocate(bk_splint_afr_mat)
	  if(allocated(bk_splint_bfr_dmat)) deallocate(bk_splint_bfr_dmat)
	  if(allocated(bk_splint_afr_dmat)) deallocate(bk_splint_afr_dmat)
	  
	  if(allocated(bk_mat_splnt)) deallocate(bk_mat_splnt)
	  if(allocated(bk_dmat_splnt)) deallocate(bk_dmat_splnt)
	  
	  end subroutine 

      SUBROUTINE DERIVS_BK_adia(t,q,dqdt) 									!Bikram: my derivs routine for adiabatic trejectories
      USE VARIABLES
      USE MPI_DATA
      USE MPI
      USE INTEGRATOR	  
      USE MPI_TASK_TRAJECT	  
      USE CONSTANTS
      USE OLD_MIJ	    
      IMPLICIT NONE
      LOGICAL BELONGS	  
      INTEGER status(MPI_STATUS_SIZE)	  
      REAL*8 dqdt(states_size*2+8), q(states_size*2+8)
      REAL*8 dq_dt_buffer(states_size*2+8)
	  INTEGER origin
      REAL*8 rphase,rfact_phase,iphase,ifact_phase
      REAL*8 d_f,t,Rrr,tet,ptet,pphi,dPhidT,dThetadT,yprm(4),sin_tet
      INTEGER i,j,k
      REAL*8 j_12,m_12,m1,m0,aver_poten_temp_0
	  integer tmp_indx
	  REAL*8,allocatable :: bk_mat(:)
	  real*8 sqrt1,sqrt2						
	  real*8,allocatable :: bk_sin_cos(:,:)	


!--------------------------------------------------------------------
! Information for classical DOF are obtained by interpolation of 
! Adiabatic Trajectory data computed previously
!--------------------------------------------------------------------
	  call splint(bk_adia_t, bk_sys_var(:,1), bk_sys_var_der(:,1), 
     & bk_adia_n, t, Rrr, yprm(1))
	  call splint(bk_adia_t, bk_sys_var(:,3), bk_sys_var_der(:,3), 
     & bk_adia_n, t, tet, yprm(2))
	  call splint(bk_adia_t, bk_sys_var(:,4), bk_sys_var_der(:,4), 
     & bk_adia_n, t, ptet, yprm(3))
	  call splint(bk_adia_t, bk_sys_var(:,6), bk_sys_var_der(:,6), 
     & bk_adia_n, t, pphi, yprm(4))
      sin_tet = dsin(tet)
	  dThetadT = 0.d0
	  if(bikram_theta) dThetadT = ptet/ Rrr/ Rrr/ massAB_M
      dPhidT = pphi/ Rrr/ Rrr/ massAB_M/ sin_tet**2	
!--------------------------------------------------------------------
! This part of the code is simillar to the CC-MQCT DERIVS_BK subroutine
! See comments there.
! Except that CS-MQCT approximation is not available within AT-MQCT
!--------------------------------------------------------------------	  
      dqdt=0d0
	  rphase = 1d0
      iphase = 0d0
	  allocate(bk_mat(mat_sz_bk))
	  call potential_data_bk_adia(Rrr,bk_mat)
	  allocate(bk_sin_cos(2,ph_cntr_bk))
	  do k= 1,ph_cntr_bk
	  bk_sin_cos(1,k) = dsin(t*bk_delta_E(k))
	  bk_sin_cos(2,k) = dcos(t*bk_delta_E(k))
	  enddo
	  
      DO k =1,mat_sz_bk
      i = ind_mat_bk(1,k)
      j = ind_mat_bk(2,k)
	  
	  if(identical_particles_defined .and. ini_st_ident) then
	  if(((-1)**(bk_parity-1)).ne.parity_state(i)) cycle
	  end if
	  
      IF(states_to_exclude_defined) THEN
      IF(stts_to_excl(i) .or. stts_to_excl(j)) CYCLE
      ENDIF
	  
      aver_poten_temp_0 = bk_mat(k)
	  d_f = 2d0
      if(i.eq.j) d_f = 1d0
	  tmp_indx = bk_indx(k)
      iphase = bk_sin_cos(1,tmp_indx)
      rphase = bk_sin_cos(2,tmp_indx)
      rfact_phase = -d_f*iphase
      ifact_phase = d_f*rphase

       dqdt(i) = dqdt(i)
     & + (q(j)*rfact_phase + q(j+states_size)*ifact_phase)	 
     &   *aver_poten_temp_0/2d0
	 
      dqdt(i+states_size) = dqdt(i+states_size)
     & + (q(j+states_size)*rfact_phase - q(j)*ifact_phase)	 
     &   *aver_poten_temp_0/2d0
	 
      dqdt(j) = dqdt(j)
     & - (q(i)*rfact_phase - q(i+states_size)*ifact_phase)	 
     &   *aver_poten_temp_0/2d0
	 
      dqdt(j+states_size) = dqdt(j+states_size)
     & - (q(i+states_size)*rfact_phase + q(i)*ifact_phase)	 
     &   *aver_poten_temp_0/2d0  
	  enddo 
	  
	  DO i=1,traject_roots	  
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	  
      origin = 0
      CALL MPI_Reduce(dqdt, dq_dt_buffer, 2*states_size+8, MPI_REAL8,
     & MPI_SUM, origin,	 
     & comms(i),ierr_mpi)
      dqdt  = dq_dt_buffer
      CALL MPI_Bcast(dqdt,2*states_size+8,
     & MPI_REAL8, origin,comms(i),ierr_mpi)
      ENDIF
      ENDDO	  
	  
	  do i=1,states_size
	  j_12 = j12(i)  
      m_12 = m12(i) 
      m1 = m_12+1d0
      m0 = m_12-1d0	  
	  
	  sqrt1 = sqrt(j_12*(j_12+1d0)-m0*(m0+1d0))					
	  sqrt2 = sqrt(j_12*(j_12+1d0)-m1*(m1-1d0))					
	  
      if(i.gt.1) then 
      dqdt(i) = dqdt(i)
     & + dPhidT/2d0*sqrt1*q(i-1+states_size)

      dqdt(i+states_size) = dqdt(i+states_size)
     & - dPhidT/2d0*sqrt1*q(i-1)
      endif 		  
	  
	  if(i.lt.states_size) then 
      dqdt(i) = dqdt(i)
     & + dPhidT/2d0*sqrt2*q(i+1+states_size)

      dqdt(i+states_size) = dqdt(i+states_size)
     & - dPhidT/2d0*sqrt2*q(i+1)
      endif
      ENDDO
      RETURN	 
      END SUBROUTINE DERIVS_BK_adia

      SUBROUTINE potential_data_bk_adia(R,bk_mat)							!Bikram: my derivs routine for adiabatic trejectories
      USE VARIABLES
      USE GRID_INTEGR
      USE CONSTANTS
      USE MPI_TASK_TRAJECT 
      USE MPI_DATA
      USE INTEGRATOR	  
      USE OLD_Mij	  
      IMPLICIT NONE	 
      REAL*8 x,v_m,der_v_m,R,V_COULPING_TERMS
      INTEGER k,k_real,i_exp_term
	  real*8 bk_mat(mat_sz_bk)

      x = R
      v_m = 0d0	  
      der_v_m = 0d0
	  
      IF(term_pot_defined) THEN
	  do k =1,mat_sz_bk
      CALL COMPUTE_V_TERMS(R)
      DO i_exp_term=1,nterms
      CALL TERM_PES_MATRIX_ELEMENT(V_COULPING_TERMS,k,i_exp_term)	  
      v_m = v_m	+
     &	V_COULPING_TERMS*V_TERMS(i_exp_term)
      der_v_m = der_v_m	+
     &	V_COULPING_TERMS*V_TERMS_DER(i_exp_term) 	 
      ENDDO		  
	  bk_mat(k)=v_m
	  enddo
      RETURN	  
      ENDIF	  
	  
      IF(x.lt.R_COM(1)) THEN
	  do k =1,mat_sz_bk
      CALL 	OUT_OF_RANGE(v_m,der_v_m,R,k)  	  
	  bk_mat(k)=v_m
	  enddo
      RETURN
      ENDIF	  
	  
      IF(x.ge.R_COM(n_r_coll)) then
	  bk_mat(:) = 0.d0
	  RETURN
	  endif
	  
	  CALL splint_bk_adia(x,bk_mat)  
	  
      END SUBROUTINE potential_data_bk_adia   

      SUBROUTINE splint_bk_adia(x,bk_mat)									!Bikram: my derivs routine for adiabatic trejectories
      USE VARIABLES
	  implicit none
      REAL*8 x,y,xaaa(n_r_coll),y2a(n_r_coll),ya(n_r_coll),y_prime
	  integer bk,k_real,bkk
      INTEGER k,khi,klo
	  real*8 chk
      REAL*8 aa,bb,h
	  real*8 cc1,cc2,cc5
	  real*8 m1,m2,dm1,dm2
	  real*8 bk_mat(mat_sz_bk)
	  
      klo=1
      khi=n_r_coll
1       if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(R_COM(k).gt.x)then
      khi=k 
      else
      klo=k
      endif
      goto 1
      endif
      h=R_COM(khi)-R_COM(klo)
      if (h.eq.0.) then
      WRITE(*,*) R_COM(n_r_coll),x,n_r_coll
      WRITE(*,*) "DISTANCE",R_COM	  
      STOP 'bad xa input in adiabatic splint' 
      endif
      aa=(R_COM(khi)-x)/h
      bb=(x-R_COM(klo))/h
	  
	  cc1=aa**3-aa
	  cc2=bb**3-bb
	  cc5=(h**2)/6.
	  if(splnt_bfr.eq.khi .and. splnt_afr.eq.klo) then
	  else
	  bk_splint_bfr_mat(:) = bk_mat_splnt(:,klo)
	  bk_splint_afr_mat(:) = bk_mat_splnt(:,khi)
	  bk_splint_bfr_dmat(:) = bk_dmat_splnt(:,klo)
	  bk_splint_afr_dmat(:) = bk_dmat_splnt(:,khi)
	  endif
	  
	  bk_mat = aa*bk_splint_bfr_mat + bb*bk_splint_afr_mat 
     &  + (cc1*bk_splint_bfr_dmat + cc2*bk_splint_afr_dmat)*cc5
	  
	  splnt_bfr = khi
	  splnt_afr = klo
	  
      return
      END

      SUBROUTINE TERM_PES_MATRIX_ELEMENT(M_coulp,k,i_term)
      USE VARIABLES
      USE EIGEN_VECTORS	  
      USE GRID_INTEGR
      USE EXPANSION_MAT_STORAGE
      USE MPI_DATA
      USE COMPUTE_EXPANSION_VARIABLES	  
      IMPLICIT NONE
      LOGICAL TRIANG_RULE
      INTEGER KRONEKER,i,k,i_r_point,ident_coeff,i_term
      REAL*8 CG,delta,W3JS,	M_coulp  
      EXTERNAL CG,delta,W3JS,TRIANG_RULE,KRONEKER
      term_limit_user = nterms	  
      M_coulp = 0d0
      buff= 0d0 
      ind_t = i_term
      i = i_term	  
      stp = ind_mat(1,k)
      stpp = ind_mat(2,k)	  
      SIMPLIFICATION_EXP_MAT = .FALSE.

      IF(.not.EXPANSION_GRID_DEFINED) THEN
      CALL READ_USER_TERMS
!!      CALL READ_COCO_TERMS	  
      ALLOCATE(TERM_MATRIX_ELEMENT(2,2,4,nterms))
      TERM_MATRIX_ELEMENT = 0d0  
      EXPANSION_GRID_DEFINED = .TRUE.
      IF(MYID.EQ.0) PRINT *, "TERM MATRIX CREATED "  
      ENDIF
      SELECT CASE(coll_type)
      CASE(1)
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)	  
      m_t = m12(stp)
      IF(m_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"
 	  
!      PRINT*, j12_p_t,j12_pp_t,m12_t	  
      j_p_t = j_ch(channp)
      j_pp_t = j_ch(channpp)
      IF(j_pp_t.ne.j12(stpp)) STOP "ERROR: INI IN BASIS"
      IF(j_p_t.ne.j12(stp)) STOP "ERROR: INI IN BASIS"	  

      l_t= A_TOP(1,ind_t)
      IF(.NOT.TRIANG_RULE(j_pp_t,l_t,j_p_t)) RETURN	  
      TERM_MATRIX_ELEMENT(1,1,1,i) =
     & CG(j_p_t,l_t,j_pp_t,m_t,0,m_t)*
     & CG(j_p_t,l_t,j_pp_t,0,0,0)*	 
     & dsqrt((2d0*j_p_t+1)/(2d0*j_pp_t+1d0))  	  
      M_coulp = M_coulp + 
     & TERM_MATRIX_ELEMENT(1,1,1,ind_t)
     
      CASE(2)
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)	  
      m_t = m12(stp)
      IF(m_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"
      j_p_t = j_ch(channp)
      j_pp_t = j_ch(channpp)
      IF(j_pp_t.ne.j12(stpp)) STOP "ERROR: INI IN BASIS"
      IF(j_p_t.ne.j12(stp)) STOP "ERROR: INI IN BASIS"	  
      l_t= A_TOP(1,ind_t)
      IF(.NOT.TRIANG_RULE(j_pp_t,l_t,j_p_t)) RETURN	  
      TERM_MATRIX_ELEMENT(1,1,1,i) =
     & CG(j_p_t,l_t,j_pp_t,m_t,0,m_t)*
     & CG(j_p_t,l_t,j_pp_t,0,0,0)*	 
     & dsqrt((2d0*j_p_t+1)/(2d0*j_pp_t+1d0))  	  
      M_coulp = M_coulp + 
     & TERM_MATRIX_ELEMENT(1,1,1,ind_t)
      M_coulp = M_coulp * vib_overlap(channp,channpp)	  
      CASE(3)
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)	  
      m_t = m12(stp)
      IF(m_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"
 	  
      j_p_t = j_ch(channp)
      j_pp_t = j_ch(channpp)
      k_p_t = k_ch(channp)
      k_pp_t = k_ch(channpp)	  
      eps_p_t = eps_ch(channp)		  
      eps_pp_t = eps_ch(channpp)	  
      IF(j_pp_t.ne.j12(stpp)) STOP "ERROR: INI IN BASIS"
      IF(j_p_t.ne.j12(stp)) STOP "ERROR: INI IN BASIS"	  
      exp_coeff_int=1d0	  
      l_t= A_TOP(1,ind_t)
      nju_t = A_TOP(2,ind_t)
      DO planr_coeff =1,2-KRONEKER(nju_t,0)
      DO par_p = 1,2
      DO par_pp =1,2	  
      symmetry_coeff = 2*par_p+par_pp-2
      k_p_t = (-1)**(par_p-1)*k_ch(channp)
      k_pp_t = (-1)**(par_pp-1)*k_ch(channpp)	  
      exp_coeff_int=(-1)**(nju_t*((1+(-1)**planr_coeff)/2))	  
      nju_t = (-1)**(planr_coeff-1)*A_TOP(2,ind_t)
      IF(.NOT.TRIANG_RULE(j_pp_t,l_t,j_p_t)) RETURN
      IF(k_pp_t.ne.k_p_t+nju_t) CYCLE
!      IF(eps_p_t.ne.eps_pp_t) CYCLE	  
      TERM_MATRIX_ELEMENT(planr_coeff,1,symmetry_coeff,i) =
     & CG(j_p_t,l_t,j_pp_t,m_t,0,m_t)*
     & CG(j_p_t,l_t,j_pp_t,k_p_t,nju_t,k_pp_t)*	 
     & dsqrt((2d0*j_p_t+1)*(2d0*
     & l_t+1)/(2d0*j_pp_t+1d0)/4d0/pi)
     & *eps_p_t**(par_p-1)*eps_pp_t**(par_pp-1)/
     % dsqrt(2d0*(1d0 + KRONEKER(k_p_t,0)))/
     & dsqrt(2d0*(1d0 + KRONEKER(k_pp_t,0)))  
      M_coulp = M_coulp + 
     & exp_coeff_int*
     & TERM_MATRIX_ELEMENT(planr_coeff,1,symmetry_coeff,ind_t)
      ENDDO
      ENDDO
      ENDDO	  
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
      exp_coeff_int=1d0	  
      l_t= A_TOP(1,ind_t)
      nju_t = A_TOP(2,ind_t)
      DO planr_coeff =1,2-KRONEKER(nju_t,0)
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int = (-1)**nju_t 
      nju_t = - A_TOP(2,ind_t)
      ENDIF

      TERM_MATRIX_ELEMENT(planr_coeff,1,1,i) = 0d0	  
      DO  k_p_t = -j_p_t,j_p_t
      k_pp_t = k_p_t +nju_t
      IF(abs(k_pp_t).gt.j_pp_t) CYCLE		  
      IF(.NOT.TRIANG_RULE(j_pp_t,l_t,j_p_t)) CYCLE
      coeff_p = M_VECTORS(channp,j_p_t+1+k_p_t)
      coeff_pp = M_VECTORS(channpp,j_pp_t+1+k_pp_t)	  
 
      TERM_MATRIX_ELEMENT(planr_coeff,1,1,i) =
     & TERM_MATRIX_ELEMENT(planr_coeff,1,1,i) + 
     & CG(j_p_t,l_t,j_pp_t,m_t,0,m_t)*
     & CG(j_p_t,l_t,j_pp_t,k_p_t,nju_t,k_pp_t)*	 
     & dsqrt((2d0*j_p_t+1)*(2d0*
     & l_t+1)/(2d0*j_pp_t+1d0)/4d0/pi)
     & *coeff_p*coeff_pp	 !!!!!! ADD EPSILON  	  
      ENDDO
	  
      M_coulp = M_coulp + 
     & exp_coeff_int*
     & TERM_MATRIX_ELEMENT(planr_coeff,1,1,ind_t)
	 
      ENDDO	  
   	  
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
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) RETURN
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) RETURN		  
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
	  
      M_coulp = M_coulp + 
     & TERM_MATRIX_ELEMENT(1,1,1,ind_t)	  

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

      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) RETURN
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) RETURN		  
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
      M_coulp_ident = M_coulp_ident + 
     & TERM_MATRIX_ELEMENT(1,1,2,ind_t)	 
    
      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**j12_p_t*par_p
	 
      M_coulp_ident = 0d0
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j2_ch(channpp)
      j2_pp_t = j1_ch(channpp)	 
	 
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) RETURN
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) RETURN		  
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
      M_coulp_ident = M_coulp_ident + 
     & TERM_MATRIX_ELEMENT(1,1,3,ind_t)		 

      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**j12_pp_t*par_pp

      M_coulp_ident = 0d0
      j1_p_t = j2_ch(channp)
      j2_p_t = j1_ch(channp)
      j1_pp_t = j2_ch(channpp)
      j2_pp_t = j1_ch(channpp)

      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) RETURN
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) RETURN	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      DO mp = -j1_p_t,j1_p_t
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE	  
!      IF(myid.eq.0) PRINT*, "COULPING", M_coulp	  
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
      M_coulp_ident = M_coulp_ident + 
     & TERM_MATRIX_ELEMENT(1,1,4,ind_t)	
	  
	  
	 
	 
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
      j12_p_t = j12(stp)	  
      m12_t = m12(stp)
      j12_pp_t = j12(stpp)
      IF(m12_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"   	  
!      PRINT*, j12_p_t,j12_pp_t,m12_t	  
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)
      v1_p = v1_ch(channp)
      v2_p = v2_ch(channp)
      v1_pp = v1_ch(channpp)
      v2_pp = v2_ch(channpp)	  
      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) RETURN
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) RETURN		  
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
      M_coulp = M_coulp + 
     & TERM_MATRIX_ELEMENT(1,1,1,ind_t)	  
	  
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

      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) RETURN
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) RETURN		  
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
      M_coulp_ident = M_coulp_ident + 
     & TERM_MATRIX_ELEMENT(1,1,2,ind_t)	 

      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**j12_p_t*par_p
	 
      M_coulp_ident = 0d0
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j2_ch(channpp)
      j2_pp_t = j1_ch(channpp)	 
	 

      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) RETURN
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) RETURN		  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      DO mp = -j1_p_t,j1_p_t
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE	  
!      IF(myid.eq.0) PRINT*, "COULPING", M_coulp	  
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

      M_coulp_ident = M_coulp_ident + 
     & TERM_MATRIX_ELEMENT(1,1,3,ind_t)		 
	

      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**j12_pp_t*par_pp

      M_coulp_ident = 0d0
      j1_p_t = j2_ch(channp)
      j2_p_t = j1_ch(channp)
      j1_pp_t = j2_ch(channpp)
      j2_pp_t = j1_ch(channpp)

	  

      l1_t = A_TOP(1,ind_t)
      l2_t = A_TOP(2,ind_t)
      l_t = A_TOP(3,ind_t)
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) RETURN
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) RETURN	  
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

      M_coulp_ident = M_coulp_ident + 
     & TERM_MATRIX_ELEMENT(1,1,4,ind_t)	
	  

 
	 
	 
      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**j12_pp_t*par_pp*par_pp*
     & (-1)**j12_p_t

      M_coulp = M_coulp/
     & dsqrt(2d0*(1d0+delta(j1_p_t,j2_p_t)))
     & /dsqrt(2d0*(1d0+delta(j1_pp_t,j2_pp_t)))

	 
     	 
      ENDIF		  
       M_coulp = M_coulp * vib_overlap(channp,channpp)	  
      CASE(7)
!!!!!!! STILL NOT CORRECT	  
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)	  
!      PRINT*, "CHANNP=",	channp, stp
!      PRINT*,"CHANNPP=",	channpp, stpp	  
      j12_p_t = j12(stp)	  
      m12_t = m12(stp)
      j12_pp_t = j12(stpp)
      IF(m12_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"   	  
!      PRINT*, j12_p_t,j12_pp_t,m12_t	  
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)
      k_p_t = k1_ch(channp)
      k_pp_t = k1_ch(channpp)	  
      eps_p_t = eps1_ch(channp)		  
      eps_pp_t = eps1_ch(channpp)	  
!      PRINT*, j1_p_t,j2_p_t,j1_pp_t,j2_pp_t	  

      exp_coeff_int=1d0
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      l_t = A_TOP(4,ind_t)
!      IF(l1_t.ne.3 .and. l1_t.ne.0 ) CYCLE
!      IF(l2_t.ne.0) CYCLE
      DO planr_coeff =1,2-KRONEKER(nju1_t,0)
      DO par_p = 1,2
      DO par_pp =1,2	  
      symmetry_coeff = 2*par_p+par_pp-2
      k_p_t = (-1)**(par_p-1)*k1_ch(channp)
      k_pp_t = (-1)**(par_pp-1)*k1_ch(channpp)	
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int=exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t)	  
      nju1_t = -nju1_t
      ENDIF	
      IF(k_pp_t.ne.k_p_t+nju1_t) CYCLE		  
      TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t) = 0d0	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t)	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) RETURN
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) RETURN	  
      CG_j1_j2_pp = CG(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,m12_t)
      CG_l1_l2 = CG(l1_t,l2_t,l_t,m_exp,-m_exp,0)
      CG_j1_l1 = CG(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp)
      CG_j2_l2 = CG(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,m12_t-mp)	  
	  
      CG_j1_k1 = CG(j1_pp_t,l1_t,j1_p_t,k_pp_t,-nju1_t,k_p_t)*	 
     & dsqrt((2d0*j1_pp_t+1)/(2d0*j1_p_t+1d0)*(2d0*l1_t+1d0))
     & * eps_p_t**(par_p-1)*eps_pp_t**(par_pp-1)/
     % dsqrt(2d0*(1d0 + KRONEKER(k_p_t,0)))/
     & dsqrt(2d0*(1d0 + KRONEKER(k_pp_t,0)))	  

      CG_j2_k2 = CG(j2_pp_t,l2_t,j2_p_t,0,0,0)*
     & dsqrt((2*j2_pp_t+1)*(2d0*l2_t+1d0)/(2*j2_p_t+1))/4d0/pi	  

      TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t) =
     & TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t)+ 	  
     & CG_j1_j2_p*
     & CG_j1_j2_pp*
     & CG_l1_l2*
     & CG_j1_l1*
     & CG_j1_k1*
     & CG_j2_l2*
     & CG_j2_k2*
     & coeff_1_p*coeff_1_pp

  
      ENDDO	  
      ENDDO

      matrix_exp_coefficent = 
     & TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t)
      M_coulp = M_coulp +
     & exp_coeff_int*matrix_exp_coefficent	 
      ENDDO
      ENDDO
      ENDDO
	  
      CASE(8)
!!!!!!! STILL NOT CORRECT	  
      CALL ASYM_TOP_VECTORS 	  
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)	  
!      PRINT*, "CHANNP=",	channp, stp
!      PRINT*,"CHANNPP=",	channpp, stpp	  
      j12_p_t = j12(stp)	  
      m12_t = m12(stp)
      j12_pp_t = j12(stpp)
      IF(m12_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"   	  
!      PRINT*, j12_p_t,j12_pp_t,m12_t	  
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)
!      PRINT*, j1_p_t,j2_p_t,j1_pp_t,j2_pp_t	  
      exp_coeff_int=1d0
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      l_t = A_TOP(4,ind_t)
!      IF(l1_t.ne.3 .and. l1_t.ne.0 ) CYCLE
!      IF(l2_t.ne.0) CYCLE
      symmetry_coeff = 1
      DO planr_coeff =1,2-KRONEKER(nju1_t,0)
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int=exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t)	  
      nju1_t = -nju1_t
      ENDIF	  
      TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t) = 0d0	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t)	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) RETURN
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) RETURN	  
      CG_j1_j2_pp = CG(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,m12_t)
      CG_l1_l2 = CG(l1_t,l2_t,l_t,m_exp,-m_exp,0)
      CG_j1_l1 = CG(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp)
      CG_j2_l2 = CG(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,m12_t-mp)	  
      DO k1p = -j1_p_t,j1_p_t
      k1pp = k1p +nju1_t
      IF(abs(k1pp).gt.j1_pp_t) CYCLE	  
      CG_j1_k1 = CG(j1_pp_t,l1_t,j1_p_t,k1pp,-nju1_t,k1p)*
     & DSQRT((2*j1_pp_t+1d0)/(2*j1_p_t+1d0)*(2*l1_t+1d0))	  
      coeff_1_p = M_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF	  

      CG_j2_k2 = CG(j2_pp_t,l2_t,j2_p_t,0,0,0)*
     & dsqrt((2*j2_pp_t+1)*(2d0*l2_t+1d0)/(2*j2_p_t+1))/4d0/pi	  

      TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t) =
     & TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t)+ 	  
     & CG_j1_j2_p*
     & CG_j1_j2_pp*
     & CG_l1_l2*
     & CG_j1_l1*
     & CG_j1_k1*
     & CG_j2_l2*
     & CG_j2_k2*
     & coeff_1_p*coeff_1_pp
      ENDDO	  
      ENDDO	  
      ENDDO

      matrix_exp_coefficent = 
     & TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t)
      M_coulp = M_coulp +
     & exp_coeff_int*matrix_exp_coefficent	 
      ENDDO
 
      CASE(9)
!!!!!!! STILL NOT CORRECT	  
      CALL ASYM_TOP_VECTORS 	  
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)	  
!      PRINT*, "CHANNP=",	channp, stp
!      PRINT*,"CHANNPP=",	channpp, stpp	  
      j12_p_t = j12(stp)	  
      m12_t = m12(stp)
      j12_pp_t = j12(stpp)
      IF(m12_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"   	  
!      PRINT*, j12_p_t,j12_pp_t,m12_t	  
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)
      k_p_t = k2_ch(channp)
      k_pp_t = k2_ch(channpp)	  
      eps_p_t = eps2_ch(channp)		  
      eps_pp_t = eps2_ch(channpp)	  
!      PRINT*, j1_p_t,j2_p_t,j1_pp_t,j2_pp_t	  

      exp_coeff_int=1d0  
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

      TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t) = 0d0	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t)	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) RETURN
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) RETURN	  
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
     & eps_p_t**(par_p-1)*eps_pp_t**(par_pp-1)/
     % dsqrt(2d0*(1d0 + KRONEKER(k_p_t,0)))/
     & dsqrt(2d0*(1d0 + KRONEKER(k_pp_t,0)))	
	  

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
     & CG_j2_k2*
     & coeff_1_p*coeff_2_p*coeff_1_pp*coeff_2_pp
  
      ENDDO	  
      ENDDO	  
      ENDDO

      matrix_exp_coefficent = 
     & TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t)
      M_coulp = M_coulp +
     & exp_coeff_int*matrix_exp_coefficent	 
      ENDDO
      ENDDO
      ENDDO	  
 
      CASE(0)
!!!!!!! STILL NOT CORRECT	  
      CALL ASYM_TOP_VECTORS 	  
      channp = indx_chann(stp)
      channpp = indx_chann(stpp)	  
!      PRINT*, "CHANNP=",	channp, stp
!      PRINT*,"CHANNPP=",	channpp, stpp	  
      j12_p_t = j12(stp)	  
      m12_t = m12(stp)
      j12_pp_t = j12(stpp)
      IF(m12_t.ne.m12(stpp)) STOP "ERROR: INI IN BASIS"   	  
!      PRINT*, j12_p_t,j12_pp_t,m12_t	  
      j1_p_t = j1_ch(channp)
      j2_p_t = j2_ch(channp)
      j1_pp_t = j1_ch(channpp)
      j2_pp_t = j2_ch(channpp)
!      PRINT*, j1_p_t,j2_p_t,j1_pp_t,j2_pp_t	  

      exp_coeff_int=1d0	  
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
!      IF(l1_t.ne.3 .and. l1_t.ne.0 ) CYCLE
!      IF(l2_t.ne.0) CYCLE
! HERE IS THE MISTAKE
      DO symmetry_coeff=1,2-KRONEKER(l1_t,l2_t)*KRONEKER(nju1_t,nju2_t)
      IF(.not.identical_particles_defined) THEN
      IF(symmetry_coeff.gt.1) CYCLE	  
      ENDIF 	  
      DO planr_coeff =1,2-KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	  
      IF(symmetry_coeff.eq.2) THEN
      exp_coeff_int=(-1)**(l1_t+l2_t)	  
      l2_t= A_TOP(1,ind_t)
      nju2_t = A_TOP(2,ind_t)
      l1_t = A_TOP(3,ind_t)
      nju1_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)	  
      ENDIF
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int=exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t)	  
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      ENDIF	  

      TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t) = 0d0	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	  
      CG_j1_j2_p = CG(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t)	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) RETURN
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) RETURN	  
      CG_j1_j2_pp = CG(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,m12_t)
      CG_l1_l2 = CG(l1_t,l2_t,l_t,m_exp,-m_exp,0)
      CG_j1_l1 = CG(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp)
      CG_j2_l2 = CG(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,m12_t-mp)	  
      DO k1p = -j1_p_t,j1_p_t
      k1pp = k1p - nju1_t
      IF(abs(k1pp).gt.j1_pp_t) CYCLE	  
      CG_j1_k1 = CG(j1_pp_t,l1_t,j1_p_t,k1pp,nju1_t,k1p)	  
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
      CG_j2_k2 = CG(j2_pp_t,l2_t,j2_p_t,k2pp,nju2_t,k2p)	  
      IF(abs(k2pp).le.j2_pp_t) THEN
      coeff_2_pp = M2_VECTORS(channpp,j2_pp_t+1+k2pp)
      ELSE
      coeff_2_pp = 0d0       	  
      ENDIF
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
     & CG_j2_k2*
     & coeff_1_p*coeff_2_p*coeff_1_pp*coeff_2_pp
      ENDDO	  
      ENDDO	  
      ENDDO	  
      ENDDO

      matrix_exp_coefficent = 
     & TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,1,ind_t)
      M_coulp = M_coulp +
     & exp_coeff_int*matrix_exp_coefficent	 
      ENDDO
      ENDDO	  
  
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
	  

  
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
      DO symmetry_coeff=1,2-KRONEKER(l1_t,l2_t)*KRONEKER(nju1_t,nju2_t)
      DO planr_coeff =1,2-KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)
      exp_coeff_int=1d0	  
      IF(symmetry_coeff.eq.2) THEN
      exp_coeff_int=(-1)**(l1_t+l2_t)	  
      l2_t= A_TOP(1,ind_t)
      nju2_t = A_TOP(2,ind_t)
      l1_t = A_TOP(3,ind_t)
      nju1_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)	  
      ENDIF
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int=exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t)	  
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      ENDIF		  

      TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,2,ind_t) = 0d0	  
      DO mp = -j1_p_t,j1_p_t
      IF(abs(m12_t-mp).gt.j2_p_t) CYCLE	 	  
      CG_j1_j2_p = CG(j1_p_t,j2_p_t,j12_p_t,mp,m12_t-mp,m12_t)	  
      DO m_exp = -min(l1_t,l2_t),min(l1_t,l2_t)
      mpp = mp - m_exp
      IF(abs(mpp).gt.j1_pp_t) CYCLE
      IF(abs(m12_t-mpp).gt.j2_pp_t) CYCLE
      IF(.NOT.TRIANG_RULE(j1_pp_t,l1_t,j1_p_t)) RETURN
      IF(.NOT.TRIANG_RULE(j2_pp_t,l2_t,j2_p_t)) RETURN		  
      CG_j1_j2_pp = CG(j1_pp_t,j2_pp_t,j12_pp_t,mpp,m12_t-mpp,m12_t)
      CG_l1_l2 = CG(l1_t,l2_t,l_t,m_exp,-m_exp,0)
      CG_j1_l1 = CG(j1_pp_t,l1_t,j1_p_t,mpp,m_exp,mp)
      CG_j2_l2 = CG(j2_pp_t,l2_t,j2_p_t,m12_t-mpp,-m_exp,m12_t-mp)	  
      DO k1p = -j1_p_t,j1_p_t
      k1pp = k1p - nju1_t
      IF(abs(k1pp).gt.j1_pp_t) CYCLE	  
      CG_j1_k1 = CG(j1_pp_t,l1_t,j1_p_t,k1pp,nju1_t,k1p)  
      coeff_1_p = M2_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M1_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF		  
      DO k2p =-j2_p_t,j2_p_t
      k2pp = k2p - nju2_t
      IF(abs(k2pp).gt.j2_pp_t) CYCLE	  
      CG_j2_k2 = CG(j2_pp_t,l2_t,j2_p_t,k2pp,nju2_t,k2p)	  
      coeff_2_p = M1_VECTORS(channp,j2_p_t+1+k2p)
      IF(abs(k2pp).le.j2_pp_t) THEN
      coeff_2_pp = M2_VECTORS(channpp,j2_pp_t+1+k2pp)
      ELSE
      coeff_2_pp = 0d0       	  
      ENDIF
      TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,2,ind_t) =
     & TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,2,ind_t) + 	  
     & CG_j1_j2_p*
     & CG_j1_j2_pp*
     & dsqrt((2d0*j1_pp_t+1)/(2d0*j1_p_t+1))*
     & dsqrt((2d0*j2_pp_t+1)/(2d0*j2_p_t+1))*
     & dsqrt((2d0*l1_t+1d0)*(2d0*l2_t+1d0))/8d0/pi**2*
     & CG_l1_l2*
     & CG_j1_l1*
     & CG_j1_k1*
     & CG_j2_l2*
     & CG_j2_k2*
     & coeff_1_p*coeff_2_p*coeff_1_pp*coeff_2_pp
      ENDDO	  
      ENDDO	  
      ENDDO	  
      ENDDO

      matrix_exp_coefficent = 
     & TERM_MATRIX_ELEMENT(symmetry_coeff,1,2,ind_t)
      M_coulp_ident = M_coulp_ident +
     & exp_coeff_int*matrix_exp_coefficent		  
      ENDDO
      ENDDO
      M_coulp_ident_1 = M_coulp_ident
      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**(j12_p_t)*par_p*parity_inversion(channp)
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

 	  
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
      DO symmetry_coeff=1,2-KRONEKER(l1_t,l2_t)*KRONEKER(nju1_t,nju2_t)
      DO planr_coeff =1,2-KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	  
      exp_coeff_int=1d0	  
      IF(symmetry_coeff.eq.2) THEN
      exp_coeff_int=(-1)**(l1_t+l2_t)	  
      l2_t= A_TOP(1,ind_t)
      nju2_t = A_TOP(2,ind_t)
      l1_t = A_TOP(3,ind_t)
      nju1_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)	  
      ENDIF
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int=exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t)	  
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      ENDIF	  

      TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,3,ind_t) = 0d0	  
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
      coeff_1_p = M1_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M2_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF		  
      DO k2p =-j2_p_t,j2_p_t
      k2pp = k2p - nju2_t
      IF(abs(k2pp).gt.j2_pp_t) CYCLE	  
      CG_j2_k2 = CG(j2_pp_t,l2_t,j2_p_t,k2pp,nju2_t,k2p)	  
      coeff_2_p = M2_VECTORS(channp,j2_p_t+1+k2p)
      IF(abs(k2pp).le.j2_pp_t) THEN
      coeff_2_pp = M1_VECTORS(channpp,j2_pp_t+1+k2pp)
      ELSE
      coeff_2_pp = 0d0       	  
      ENDIF

       TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,3,ind_t) =
     & TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,3,ind_t)+ 	  
     & CG_j1_j2_p*
     & CG_j1_j2_pp*
     & dsqrt((2d0*j1_pp_t+1)/(2d0*j1_p_t+1))*
     & dsqrt((2d0*j2_pp_t+1)/(2d0*j2_p_t+1))*
     & dsqrt((2d0*l1_t+1d0)*(2d0*l2_t+1d0))/8d0/pi**2*
     & CG_l1_l2*
     & CG_j1_l1*
     & CG_j1_k1*
     & CG_j2_l2*
     & CG_j2_k2*
     & coeff_1_p*coeff_2_p*coeff_1_pp*coeff_2_pp
      ENDDO	  
      ENDDO	  
      ENDDO	  
      ENDDO

      matrix_exp_coefficent = 
     & TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,3,ind_t)
      M_coulp_ident = M_coulp_ident +
     & exp_coeff_int*matrix_exp_coefficent	  
   
      ENDDO
      ENDDO	  
 
      M_coulp_ident_2 = M_coulp_ident
      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**j12_pp_t*par_pp*parity_inversion(channpp)

      M_coulp_ident = 0d0
      buff = 0d0	  
      j1_p_t = j2_ch(channp)
      j2_p_t = j1_ch(channp)
      j1_pp_t = j2_ch(channpp)
      j2_pp_t = j1_ch(channpp)	  
	  


	  
      l1_t= A_TOP(1,ind_t)
      nju1_t = A_TOP(2,ind_t)
      l2_t = A_TOP(3,ind_t)
      nju2_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)
      DO symmetry_coeff=1,2-KRONEKER(l1_t,l2_t)*KRONEKER(nju1_t,nju2_t)
      DO planr_coeff =1,2-KRONEKER(nju1_t,0)*KRONEKER(nju2_t,0)	  
      exp_coeff_int=1d0	  
      IF(symmetry_coeff.eq.2) THEN
      exp_coeff_int=(-1)**(l1_t+l2_t)	  
      l2_t= A_TOP(1,ind_t)
      nju2_t = A_TOP(2,ind_t)
      l1_t = A_TOP(3,ind_t)
      nju1_t = A_TOP(4,ind_t)
      l_t = A_TOP(5,ind_t)	  
      ENDIF
      IF(planr_coeff.eq.2) THEN
      exp_coeff_int=exp_coeff_int*(-1)**(l1_t+l2_t+l_t+nju1_t+nju2_t)	  
      nju1_t = -nju1_t
      nju2_t = -nju2_t	  
      ENDIF		  

      TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,4,ind_t) = 0d0	  
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
      coeff_1_p = M2_VECTORS(channp,j1_p_t+1+k1p)
      IF(abs(k1pp).le.j1_pp_t) THEN
      coeff_1_pp = M2_VECTORS(channpp,j1_pp_t+1+k1pp)
      ELSE
      coeff_1_pp = 0d0       	  
      ENDIF		  
      DO k2p =-j2_p_t,j2_p_t
      k2pp = k2p - nju2_t
      IF(abs(k2pp).gt.j2_pp_t) CYCLE	  
      CG_j2_k2 = CG(j2_pp_t,l2_t,j2_p_t,k2pp,nju2_t,k2p)	  
      coeff_2_p = M1_VECTORS(channp,j2_p_t+1+k2p)
      IF(abs(k2pp).le.j2_pp_t) THEN
      coeff_2_pp = M1_VECTORS(channpp,j2_pp_t+1+k2pp)
      ELSE
      coeff_2_pp = 0d0       	  
      ENDIF
      TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,4,ind_t) =
     & TERM_MATRIX_ELEMENT(symmetry_coeff,planr_coeff,4,ind_t)+ 	  
     & CG_j1_j2_p*
     & CG_j1_j2_pp*
     & dsqrt((2d0*j1_pp_t+1)/(2d0*j1_p_t+1))*
     & dsqrt((2d0*j2_pp_t+1)/(2d0*j2_p_t+1))*
     & dsqrt((2d0*l1_t+1d0)*(2d0*l2_t+1d0))/8d0/pi**2*
     & CG_l1_l2*
     & CG_j1_l1*
     & CG_j1_k1*
     & CG_j2_l2*
     & CG_j2_k2*
     & coeff_1_p*coeff_2_p*coeff_1_pp*coeff_2_pp
      ENDDO	  
      ENDDO	  
      ENDDO	  
      ENDDO

      matrix_exp_coefficent = 
     & TERM_MATRIX_ELEMENT(symmetry_coeff,1,4,ind_t)
      M_coulp_ident = M_coulp_ident +
     & exp_coeff_int*matrix_exp_coefficent	  
      ENDDO	  
      ENDDO
  
      M_coulp_ident_3 = M_coulp_ident

      M_coulp = M_coulp + M_coulp_ident*
     ^ (-1)**j12_pp_t*par_pp*par_p*
     & (-1)**j12_p_t*parity_inversion(channp)*
     & parity_inversion(channpp)

      M_coulp = M_coulp/
     & dsqrt(2d0*(1d0+delta(j1_p_t,j2_p_t)))
     & /dsqrt(2d0*(1d0+delta(j1_pp_t,j2_pp_t)))
   	  
      ENDIF	  
	  
	  
	  
	  
      END SELECT 	  
      END SUBROUTINE TERM_PES_MATRIX_ELEMENT	  
