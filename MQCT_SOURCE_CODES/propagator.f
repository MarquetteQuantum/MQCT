      MODULE INTEGRATOR !! CONTAINS VARIABLES FOR PROPAGATION
! This module is updated by Bikramaditya Mandal
      IMPLICIT NONE
      LOGICAL :: INTERUPT_PROPAGATION = .FALSE.
      LOGICAL :: term_pot_defined = .FALSE.
      LOGICAL ident_skip 
      LOGICAL :: head_on = .FALSE.	  
      LOGICAL, ALLOCATABLE :: bill_exst(:,:)
      LOGICAL mpi_task_per_proc,traj_run
      INTEGER s_st_size,s_st_check,p_cur,j_cur,m_cur,st_traject,
     & parity_st
      INTEGER :: MAXSTP = 1000
      INTEGER ang_res_elast	  
      INTEGER loop_count,vibration_cnt
      INTEGER*4 idum
      INTEGER mpi_traject	  
      INTEGER tot_number_of_traject,ini_channel
      INTEGER L_MAX_TRAJECT,L_MIN_TRAJECT,
     & J_UP_INT,J_DOWN_INT, j_int_ini,dl_step_integ	  
      INTEGER :: odeint_fail_traject = 0
      INTEGER, ALLOCATABLE :: itraj_myid_l(:,:)	  
      INTEGER itraject,ntraject,n_traject_alloc,id_proc,id_proc_root	  
      INTEGER s_st,i_ener,task_size,i_curr
      INTEGER num_traj_ener,num_traj_prb
      INTEGER sys_var_size
      INTEGER traject_sum_lim
      INTEGER timer_random(8)
      INTEGER j12_s_st,m12_s_st	  
      INTEGER, ALLOCATABLE :: ERROR_INDEX_PROBAB(:)
      INTEGER, ALLOCATABLE :: ERROR_INDEX_ENERGY(:)
      INTEGER, ALLOCATABLE :: err_proc_max_ener(:),err_proc_max_prb(:)
      INTEGER, ALLOCATABLE :: err_traj_max_ener(:),err_traj_max_prb(:)
      INTEGER, ALLOCATABLE :: 
     & err_state_max_ener(:),err_state_max_prb(:),l_scatter(:,:)	  
      REAL*8, ALLOCATABLE :: sys_var(:),deriv_sys(:),probab(:)
     &,probab_J(:,:), total_probab(:),sigma_f(:,:),
     & sigma_m(:,:),deflectangle(:,:),probab_J_all(:,:,:),sigma(:),
     & sigma_m_j(:,:,:), opacity(:,:,:)!, probab_bikram(:)
      REAL*8, ALLOCATABLE :: probab_traj(:)	 
      REAL*8, ALLOCATABLE:: E_coll(:), E_bill(:,:)
      REAL*8, ALLOCATABLE :: err_ener_larg(:),
     &	  err_prb_larg(:)
      REAL*8, ALLOCATABLE ::  error_ener_largest(:),
     & error_prb_largest(:)
      REAL*8 err_prb_tmp, err_ener_tmp,err_id_ener,err_id_prb
      REAL*8, PARAMETER :: ENER_REF = 300d0	  
      REAL*8 E_sct,k_vec,current_vel_ini(3)
      REAL*8 ::  current_error = 0d0
      REAL*8 dE_bllng,exchange_ident
      REAL*8 TIME_BUFFER	  
      REAL*8, PARAMETER :: pi = dacos(-1d0)
      REAL*8 dJ_int_range
      REAL*8 rand1,rand2,time_st,rand_m,rand_j,rand_par
      REAL*8, ALLOCATABLE :: ampl_wf_real(:,:),ampl_wf_imag(:,:),
     & phase_wf(:,:),angle_phase(:,:),angle_scatter(:,:),
     & buffer_elastic(:)
	 
! Bikram Start:
	  REAL*8, ALLOCATABLE :: bk_tym(:,:)
	  REAL*8, ALLOCATABLE :: bk_prd(:,:)
	  REAL*8, ALLOCATABLE :: bk_vib(:,:)
	  REAL*8, ALLOCATABLE :: bk_erre(:,:)
	  REAL*8, ALLOCATABLE :: bk_errp(:,:)
	  REAL*8, ALLOCATABLE :: bk_time(:)
      INTEGER, ALLOCATABLE :: bk_period(:)
      INTEGER, ALLOCATABLE :: bk_vibration(:)
	  REAL*8, ALLOCATABLE :: bk_err_E(:)
	  REAL*8, ALLOCATABLE :: bk_err_P(:)
	  integer bk_d,left_traj,right_traj
! Bikram End.
	 
      REAL*8, ALLOCATABLE :: phase_scatter(:)
      REAL*8, ALLOCATABLE :: def_fnc(:)	 
      REAL*8, ALLOCATABLE :: sigma_elast(:)
      REAL*8, ALLOCATABLE :: mean_values(:,:)
      REAL*8, ALLOCATABLE :: variance_values(:,:)
      REAL*8, ALLOCATABLE :: V_TERMS(:),V_TERMS_DER(:)
      REAL*8, ALLOCATABLE :: ORBIT_TRAJ_DATA(:,:)
      REAL*8, ALLOCATABLE :: TRAJECT_END_DATA(:,:,:,:)
      REAL*8, ALLOCATABLE :: TRAJECT_END_DATA_ALL(:,:,:,:,:)	  
      LOGICAL, ALLOCATABLE :: ORBIT_TRAJ_FLAG(:)
      REAL*8, ALLOCATABLE :: ORBIT_TRAJ_DATA_ALL(:,:,:)
      REAL*8 j_curr_f,m_curr_f	  
      LOGICAL, ALLOCATABLE :: ORBIT_TRAJ_FLAG_ALL(:,:)		  
      REAL*8 dt_step,dt_corrected 	  
      END MODULE INTEGRATOR
      MODULE CHECK_FILE !!! CONTAINS VASRIABLES FOR CHECK_POINT FILE
! This module is updated by Bikramaditya Mandal
      IMPLICIT NONE
      CHARACTER(LEN=14) :: CHECK_FILE_NAME = "CHECK_FILE.DAT"
      CHARACTER(LEN=16) :: TRAJECT_NAME = "TRAJECTORIES.DAT"	  
      INTEGER num_ini_states,num_chann_check_point,ini_st_check_file
      REAL*8, ALLOCATABLE :: all_prob_l(:,:,:),sv_ch_prob_l(:,:,:,:)
      REAL*8, ALLOCATABLE :: sv_ch_prob_l_buffer(:,:,:)	  
      INTEGER, ALLOCATABLE :: all_def_fnc(:,:),sv_ch_def_fnc(:,:,:),
     & sv_ch_def_fnc_buffer(:,:),sv_ch_vib_fnc_buffer(:,:),
     & sv_ch_vib_fnc(:,:,:),all_vib_fnc(:,:)
      REAL*8, ALLOCATABLE :: TRAJECT_DATA_CHECK_FILE(:)
      INTEGER, ALLOCATABLE :: what_computed(:,:),
     & sv_ch_what_computed(:,:,:),sv_ch_what_computed_buffer(:,:) 
!      REAL*8, ALLOCATABLE :: final_traject_info(:,:,:,:)	!!! FINALIZE TOMORROW 
      INTEGER id_trajec_mnt_crl	 
	  integer bk_num_ini_states
      END MODULE CHECK_FILE
      MODULE MONTE_CARLO_STORAGE !!! CONTAINS VARAIBLES FOR MONTE CARLO INTEGRATION
! This module is updated by Bikramaditya Mandal
      IMPLICIT NONE
      INTEGER p_monte_c,tot_number_of_traject_mc,tot_number_to_cut	  
      REAL*8, ALLOCATABLE :: x_mean_storage(:),x2_mean_storage(:)
      REAL*8, ALLOCATABLE :: monte_carlo_err(:,:)
      REAL*8, ALLOCATABLE :: total_probab_check_monte_carlo(:)
      REAL*8, ALLOCATABLE :: x_mean_storage_ch(:),x2_mean_storage_ch(:)
      REAL*8 buffer_mean,buffer_mean2	  
      END MODULE MONTE_CARLO_STORAGE
      MODULE CS_MATRIX
! This module is updated by Bikramaditya Mandal
      IMPLICIT NONE
      REAL*8,ALLOCATABLE :: sys_var_cs(:),deriv_sys_cs(:),Ej_cs_im(:)
      REAL*8, ALLOCATABLE :: Mat_el_cs(:,:,:), Mat_el_cs_der(:,:,:)
      REAL*8, ALLOCATABLE :: mat_buffer(:,:),mat_buffer_der(:,:)	  
      INTEGER states_size_cs,total_size_cs
      INTEGER size_csst_chunk_mpi
      INTEGER residue_csst_mpi
      INTEGER size_csts_chunk_mpi
      INTEGER mat_ts_mpi	  
      INTEGER residue_csts_mpi	  
      INTEGER, ALLOCATABLE :: ind_mat_cs(:) 
      INTEGER, ALLOCATABLE :: ind_state_cs(:)
      INTEGER, ALLOCATABLE :: portion_of_csst_per_task(:,:)	
      INTEGER, ALLOCATABLE :: portion_of_csts_per_task(:,:)
      INTEGER, ALLOCATABLE :: ind_cs_rule(:,:)	  
      LOGICAL, ALLOCATABLE :: integrator_flag(:)
      END MODULE CS_MATRIX		  
      
	  module bk_l_values
! This module is written by Bikramaditya Mandal, Feb 2021
	  implicit none
	  integer chk_par, postv_par, negtv_par, delta_l_lr, l_switch_bk
	  integer l_range1, l_range2
	  logical transfer_prob_spln
	  end module
	  
	  module monte_carlo_sampling
! This module is written by Bikramaditya Mandal, Sept 2021
	  implicit none
	  real*8 J_tot_bk, l_real_bk
	  end module
	  
	  SUBROUTINE PROPAGATE
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE CONSTANTS
      USE CHECK_FILE	  
      USE INTEGRATOR
      USE MPI_DATA
      USE MPI_TASK_TRAJECT
      USE MPI
      USE ERRORS	  
      USE MONTE_CARLO_STORAGE	  
	  use bk_l_values															!Bikram Feb 2021
	  use monte_carlo_sampling													!Bikram Sept 2021
	  use iso_fortran_env,only:output_unit
      IMPLICIT NONE
      LOGICAL sampl_succ	  
      INTEGER s_ini,m_t,j_t,j1_t,j2_t,j12_min,j12_max,j_count,j_summ,st
      INTEGER i,j,k,l_parity,p_parity,ident_max,KRONEKER
      INTEGER l_counter,round 
      REAL*8 delta,rand0,MIDDLE_NUM
	  integer temp_l_counter										    		!Bikram
	  real*8 temp_t_bfr_dotrjct,temp_t_afr_dotrjct,t_prop			    		!Bikram
	  real*8 prop_time_bgn, prop_time_end, prop_time
	  real*8, allocatable :: bk_probab_J(:,:), bk_probab_J_all(:,:,:)
	  logical if_l																!Bikram
	  
! Bikram Start August 2021: Monte-Carlo Adiabatic Calculations
	  CHARACTER(LEN=100) bk_adia_dir1,bk_adia_filepath,bk_adia_dir2
	  integer j_ini_tr, origin, ii, dmm1, dmm2, dmm3, mc_traj
	  integer total_traject_bikram, old_nmb_traj, io_traj_file,tmp_avg
	  integer, allocatable :: bk_s_st(:)
	  real*8, allocatable :: bk_l_real(:)
	  real*8 J_tot, l_real
	  logical belongs, mc_same_traj, mc_traj_file_exst, chk_files_exst
	  logical reduce_nmb_traj
	  character (len = 100) :: chk_files_name, bk_label
	  integer :: chk_file_unit = 123123
	  external belongs
! BIkram End. 
      EXTERNAL delta,rand0,KRONEKER,round
	  
	  
!--------------------------------------------------------------------
! MPI initilization
!--------------------------------------------------------------------
      mpi_task_per_proc = mpi_task_defined
      sys_var_size = 8+states_size*2 	!!!! INTIALIZNG  
      IF(mpi_task_per_proc) mpi_traject = mpi_task_per_traject
      IF(.not.mpi_task_per_proc) mpi_traject = 1
      dl_step_integ =  delta_l_step
      IF(.not.expansion_defined .and. pot_expansion_defined
     & 	.and. .not. mpi_task_defined) THEN !!! IN CASE OF MOLSCAT LIKE POTENTIAL
      term_pot_defined = .TRUE.           !!! ALCCATING WORKING ARRAYS
      ALLOCATE (V_TERMS(nterms),V_TERMS_DER(nterms))	  
      ENDIF	 	  
!!! CHECKING INPUT SIZE OF THE BASIS
      IF(states_size.le.0) STOP"TOTAL NUMBER OF STATES IS < 0"
      IF(myid.eq.0) PRINT*, " "
      IF(myid.eq.0) PRINT*, "SCATTERING CALCULATIONS STARTED"
	  prop_time_bgn = mpi_wtime()
!!! ALLOCATING ARRAYS
!!! ALLOCATION OF POPULATION ARRAY FOR PRINTING TRAJECTORIES
      IF(prn_l_defined) THEN
	  ALLOCATE(probab_traj(number_of_channels))
      probab_traj = 0d0	  
      ENDIF
      ALLOCATE(sigma_elast(states_size))
      sigma_elast = 0d0	  
      ALLOCATE(sys_var(2*states_size+8),deriv_sys(2*states_size+8))
!      ALLOCATE(probab_bikram(number_of_channels))
      ALLOCATE(probab(number_of_channels),
     & total_probab(number_of_channels),sigma(number_of_channels))
      ALLOCATE(E_coll(nmbr_of_enrgs),
     & sigma_f(number_of_channels,nmbr_of_enrgs),
     & E_bill(number_of_channels,nmbr_of_enrgs),
     & bill_exst(number_of_channels,nmbr_of_enrgs))
      ALLOCATE(error_ener_largest(nmbr_of_enrgs),
     & error_prb_largest(nmbr_of_enrgs))
      ALLOCATE(ERROR_INDEX_PROBAB(nproc),
     & ERROR_INDEX_ENERGY(nproc))
      ALLOCATE(err_proc_max_ener(nmbr_of_enrgs))
      ALLOCATE(err_proc_max_prb(nmbr_of_enrgs))
      ALLOCATE(err_state_max_ener(nmbr_of_enrgs))
      ALLOCATE(err_state_max_prb(nmbr_of_enrgs))	  
      ALLOCATE(err_traj_max_ener(nmbr_of_enrgs))
      ALLOCATE(err_traj_max_prb(nmbr_of_enrgs))
      ALLOCATE( monte_carlo_err(number_of_channels,nmbr_of_enrgs))
	  reduce_nmb_traj = .false.												!Bikram
	  if_l = .false.														!Bikram
	  transfer_prob_spln = .false.
	  if(bk_prob_interpolation) transfer_prob_spln = .true.

      ERROR_INDEX_ENERGY = 0d0
!--------------------------------------------------------------------
! Billing energy calculations for all transitions and all energies
!--------------------------------------------------------------------
      DO k=1,nmbr_of_enrgs
      DO j=1,number_of_channels
      dE_bllng = E_ch(j) - E_ch(chann_ini)
      IF(abs(dE_bllng/4d0) .le. U(k) ) THEN
      bill_exst(j,k) = .TRUE.
      ELSE
      bill_exst(j,k) = .FALSE.	  
      ENDIF	  
      E_bill(j,k) = U(k) + dE_bllng/2d0
     & + dE_bllng**2/16d0/U(k)
      ENDDO
      ENDDO
!--------------------------------------------------------------------
! Create file to print resonance trajectories informations
!--------------------------------------------------------------------
! Bikram Oct'18 Start:
      IF(orbit_traj_defined) then
	  if(myid.eq.0) then
      OPEN(2345,FILE="RESONANCE_TRAJECT.out",POSITION = "APPEND")
      WRITE(2345,'(a36)')"RESONANCE TRAJECTORIES IF THEY EXIST"	  
      close(2345)
      endif
      endif
! Bikram End
!--------------------------------------------------------------------
! Decides the value of j12 and m12 of the initial channel to print trajectory
!--------------------------------------------------------------------
      IF(prn_l_defined) THEN
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) STOP "NOT READY"	  
      IF(j12m12_print(1).lt.0) THEN
      IF(coll_type.gt.4 .or. coll_type.eq.0) THEN	  
      j12m12_print(1) = max(j1_ch(chann_ini),j2_ch(chann_ini))
      j12m12_print(2) = 0
      ELSE
      j12m12_print(1) = j_ch(chann_ini)	  
      j12m12_print(2) = 0
      ENDIF	  
      ENDIF	  
      ENDIF
!--------------------------------------------------------------------
! Setting up initial conditions
!--------------------------------------------------------------------
      error_ener_largest = 0d0
      error_prb_largest = 0d0
      monte_carlo_err = 0d0	  
      err_ener_tmp = 0d0
      err_prb_tmp = 0d0	  
      tot_number_of_traject = 0  
      E_coll = U/eVtown/autoeV
      s_ini = chann_indx(chann_ini)
      ini_channel = chann_ini
      ident_max = 1
	   
      IF(identical_particles_defined) THEN
      ident_skip = identical_particles_defined
      exchange_ident = exch_par_w_pl	  
      SELECT CASE(coll_type)		  
      CASE(5)	  

      p_lim_min = 1
      p_lim_max = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
      CASE(6)
      p_lim_min = 1
      p_lim_max = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))*
     & KRONEKER(v1_ch(chann_ini),v2_ch(chann_ini)) 	  
      CASE(0)
      p_lim_min = 1
      p_lim_max = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
     & *KRONEKER(ka1_ch(chann_ini),ka2_ch(chann_ini))
     & *KRONEKER(kc1_ch(chann_ini),kc2_ch(chann_ini))	  
      END SELECT	  

      ident_max = p_lim_max
	  
! Bikram Start Jan 2021:
	  if(ident_max.eq.2d0) then
	  ini_st_ident = .true.
	  else
	  ini_st_ident = .false.
	  end if
! Bikram End.
	  
      ENDIF	
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) THEN
      IF(round(m12_h(s_ini)*2) .ne. 1)  THEN
      PRINT*, myid,m12_h(s_ini)	  
	  STOP"ERROR:INI WENT WRONG 1/2"
      ENDIF	  
      ELSE	  
      IF(m12(s_ini).ne.0) STOP"ERROR:INI WENT WRONG"
      ENDIF	  
      sigma_f = 0d0
      p_cur = p_lim_min
      j_cur = j_min_ind(chann_ini)
      m_cur = 0
      st_traject = 1
      time_st = 0d0
      i_curr = 1
      tot_number_of_traject = nmbr_of_traj
!--------------------------------------------------------------------
! Old Monte-Carlo checkpointing. Disabled by Bikram
!--------------------------------------------------------------------
!!! FOR MONTE CARLO NOT CHECK POINT OR IDENTICAL SYMMETRY LOOP	  
      ! IF(monte_carlo_defined .and. check_point_defined) THEN
      ! CALL READ_CHECK_POINT_MC
      ! ENDIF	  
!!! SETTING UP INTITAL CONDIONS	  
!      PRINT*,myid,check_point_defined	  
      ! IF(check_point_defined .and. .not. monte_carlo_defined) THEN
      ! IF(MYID.EQ.0) PRINT*,"READING FROM CHECKPOINT STARTED"	  
      ! CALL READ_CHECK_POINT	
      ! IF(MYID.EQ.0) PRINT*,"READING FROM CHECKPOINT FINISHED"
      ! CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
! !      STOP 
      ! ENDIF
!!! CHECKPOINT READING	

!--------------------------------------------------------------------
! LOOP OVER ENERGIES
!--------------------------------------------------------------------
      DO i_ener = i_curr, nmbr_of_enrgs
!--------------------------------------------------------------------
! Creating directory to store checkpoint files
!--------------------------------------------------------------------
	  if(write_check_file_defined) then
	  inquire(file = trim(bk_dir44), exist = chk_files_exst)
	  if(chk_files_exst) call system ( "rm -r " // trim(bk_dir44) )
	  call system ( "mkdir -p " // trim(bk_dir44) )
	  write(chk_files_name, '(a,a,i0,a)') trim(bk_dir44), 
     & "/Proc_ID_", myid, ".dat"
	  chk_files_name = trim(chk_files_name)
	  open(newunit = chk_file_unit, file = trim(chk_files_name))
	  end if
	  if(check_point_defined) then
	  write(chk_files_name, '(a,a,i0,a)') trim(bk_dir44), 
     & "/Proc_ID_", myid, ".dat"
	  chk_files_name = trim(chk_files_name)
	  open(newunit = chk_file_unit, file = trim(chk_files_name))
	  end if
      E_sct = E_coll(i_ener)
!--------------------------------------------------------------------
! Computes the number of trajectories per procs group
!--------------------------------------------------------------------  
      IF(monte_carlo_defined) THEN
!--------------------------------------------------------------------
! This is the case of Monte-Carlo
!--------------------------------------------------------------------  
      ident_max  = 1		  
      IF(check_point_defined .and. i_curr.eq.i_ener) THEN
      tot_number_of_traject = nmbr_of_traj - tot_number_to_cut
      ENDIF
      IF(check_point_defined .and. i_curr.lt.i_ener) THEN
      tot_number_of_traject = nmbr_of_traj
      ENDIF	 	  
!!!!! IF MONTE CARLO DEFINED    	  
      IF(mpi_task_defined) THEN	  
      n_traject_alloc = int(tot_number_of_traject/(nproc/mpi_traject))+1
      ELSE
      n_traject_alloc = int(tot_number_of_traject/nproc)+1	  
      ENDIF
!!!! NUMBER OF TRAJECTORIES	
      ntraject = n_traject_alloc-1  
      IF(mpi_task_defined) THEN		  
      IF(myid/mpi_traject
     & .lt. tot_number_of_traject - ntraject*(nproc/mpi_traject) )
     & ntraject = ntraject+1
      ELSE
      IF(myid .lt. tot_number_of_traject - ntraject*nproc )
     & ntraject = ntraject+1	  
      ENDIF
      IF(MYID.EQ.0) WRITE(*,*)
     & "NUMBER OF TRAJECTORIES PER PROC_GROUP",ntraject
      IF(MYID.EQ.0) WRITE(*,*)
     & "TOTAL NUMBER OF TRAJECTORIES",tot_number_of_traject
	 
      ELSE
!--------------------------------------------------------------------
! For regular calculations, it computes the total number of trajectories
!--------------------------------------------------------------------  
      IF(b_impact_defined) THEN
      J_tot_max = int(b_impact_parameter*sqrt(massAB_M*2d0*E_sct))
     & -j_max_ind(chann_ini)	  ! DB change 8/23
      J_tot_min = 0
      ENDIF
       J_DOWN_INT = J_tot_min
       J_UP_INT  = 	J_tot_max
      IF(myid.eq.0)      
     & WRITE(*,
     & '(a19,1x,i4,1x,a2,1x,i4)')
     & "TOTAL J RANGE: FROM",J_DOWN_INT,"TO",J_UP_INT	   
!!!  DEFINING L_MIN AND L_MAX  
      L_MAX_TRAJECT = J_tot_max+j_max_ind(chann_ini)
      L_MIN_TRAJECT = L_MAX_TRAJECT	  
      DO i= J_tot_min, J_tot_max
      DO j= j_min_ind(chann_ini), j_max_ind(chann_ini)
      L_MIN_TRAJECT = min(L_MIN_TRAJECT,
     & abs(i-j) )      	  
      ENDDO		  
      ENDDO		  
      IF(myid.eq.0)      
     & WRITE(*,
     & '(a19,1x,i4,1x,a2,1x,i4)')
     & "TOTAL L RANGE: FROM",L_MIN_TRAJECT,"TO",L_MAX_TRAJECT  
!!!!! TOTAL NUMBER OF TRAJECTORIES      
      tot_number_of_traject = 
     & (L_MAX_TRAJECT - L_MIN_TRAJECT)/delta_l_step + 1

! Bikram Start Feb 2021: This is related to switching the delta_L values in two regions
	  l_switch_bk = -100d0
	  if(bk_b_switch.gt.0.d0) then
	  l_switch_bk = int(bk_b_switch*sqrt(massAB_M*2d0*E_sct))
	  if(myid.eq.0) write(*,'(a,f0.3,a,i0)') 
     & 'B_SWITCH INDICATED, VALUE = ', 
     & bk_b_switch, ' CALCULATED VALUE OF L TO SWITCH = ', l_switch_bk
	  delta_l_lr = bk_dl_lr
	  
	  l_range1 = (l_switch_bk - L_MIN_TRAJECT)/delta_l_step + 1
	  l_range2 = (L_MAX_TRAJECT - l_switch_bk)/delta_l_lr + 1
	  if(myid.eq.0) write(*,'(a,i0)') 
     & '#TRAJECTORIES IN SHORT RANGE = ', l_range1
	  if(myid.eq.0) write(*,'(a,i0)') 
     & '#TRAJECTORIES IN LONG RANGE = ', l_range2
	  
	  tot_number_of_traject = l_range1 + l_range2
	  end if
! Bikram End.

!!!   NUMBER OF TRAJECTORIES PER MPI TASK	  
      IF(mpi_task_defined) THEN	  
      n_traject_alloc = int(tot_number_of_traject/(nproc/mpi_traject))+1
      ELSE
      n_traject_alloc = int(tot_number_of_traject/nproc)+1	 !!!! CHANGE AFTER ALL 
      ENDIF
!!!! NUMBER OF TRAJECTORIES	
      ntraject = n_traject_alloc-1 !!! FIX AFTER ALL 
      IF(mpi_task_defined) THEN		  
      IF(myid/mpi_traject
     & .lt. tot_number_of_traject - ntraject*(nproc/mpi_traject) )
     & ntraject = ntraject+1
      ELSE
      IF(myid .lt. tot_number_of_traject - ntraject*nproc )
     & ntraject = ntraject+1	  
      ENDIF
      IF(MYID.EQ.0) WRITE(*,*)
     & "NUMBER OF TRAJECTORIES PER PROC_GROUP MAX",ntraject
	  IF(MYID.EQ.nproc-1) WRITE(*,*)
     & "NUMBER OF TRAJECTORIES PER PROC_GROUP MIN",ntraject
      IF(MYID.EQ.0) WRITE(*,*)
     & "TOTAL NUMBER OF TRAJECTORIES",tot_number_of_traject	  
      ENDIF
!--------------------------------------------------------------------
! This is to check if #trajectories < (nproc/mpi_pertraj)
! In this case, the code would stop propagation with an error message
!--------------------------------------------------------------------    
! Bikram Jan 2020 Start:
	  if(tot_number_of_traject.lt.(nproc/mpi_traject)) then
	  if(myid.eq.0) then
      write(*,*)' '
      write(*,*)'You Requested:'
      write(*,*)'#trajectories =              ',tot_number_of_traject
      write(*,*)'#processors =                ',nproc
      write(*,*)'#processors per trajectory = ',mpi_traject
      write(*,*)' '
      write(*,*)'It is impossible to distribute the load equally.'
      write(*,*)'Please try reducing the #processors'
      write(*,*)'or increase #trajectories'
      write(*,*)'or increase #processors per trajectory.'
      write(*,*)' '
	  endif
	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  stop
	  endif
! Bikram End.
	  
!!!! IF CHECKPOINT_SAVE DEFINED
      IF(write_check_file_defined .and. .not.monte_carlo_defined) THEN
!--------------------------------------------------------------------
! Old Monte-Carlo checkpointing. Disabled by Bikram
!--------------------------------------------------------------------
      ! IF(fine_structure_defined .and. SPIN_FINE.eq.2) STOP "NOT READY"	  
      ! m_t = 0!!! DEFINING HOW MANY DIFFERENT J12 and M12 WE have
      ! j_t =  j_min_ind(chann_ini)	  
      ! IF(.not.identical_particles_defined) THEN
      ! i = indx_corr(m_t + j_t+1,
      ! & j_t+1,chann_ini) 
      ! ELSE
      ! i = indx_corr_id(p_cur,m_t + j_t+1,
      ! & j_t+1,chann_ini)
      ! ENDIF
      ! j_t =  j_max_ind(chann_ini)
      ! m_t = j_t		  
      ! IF(.not.identical_particles_defined) THEN
      ! j = indx_corr(m_t + j_t+1,
      ! & j_t+1,chann_ini) 
      ! ELSE
      ! j = indx_corr_id(ident_max,m_t + j_t+1,
      ! & j_t+1,chann_ini) 	  
      ! ENDIF	
      ! num_ini_states = j - i	+ 1
      ! ini_st_check_file = i
      ! IF(s_ini.ne.ini_st_check_file .or. num_ini_states.le.0) THEN
      ! IF(MYID.eq.0) PRINT*,"ERROR IN INITIALIZATION OF CHECK FILES"
      ! IF(MYID.eq.0) PRINT*,"ERROR IN CHECK_FILE"      	  
      ! STOP
      ! ENDIF	  
      ! IF (MYID.eq.0 ) THEN
      ! num_chann_check_point = 	number_of_channels
! !!!  THAT IS WHAT WOULD BE WRITTERN IN CHECKPOINT
      ! IF(.not.check_point_defined .or. i_ener.gt.i_curr) THEN	  
      ! ALLOCATE(what_computed(num_ini_states,tot_number_of_traject))	  
      ! ALLOCATE(      
      ! & all_prob_l
      ! & (sys_var_size+6,num_ini_states,tot_number_of_traject))
      ! ALLOCATE(all_def_fnc(num_ini_states,tot_number_of_traject))
      ! ALLOCATE(all_vib_fnc(num_ini_states,tot_number_of_traject))	  
      ! ENDIF	  
 ! !!! CHECKPOINT ARRAYS	 
      ! ALLOCATE(      
      ! & sv_ch_prob_l
      ! & (sys_var_size+6,num_ini_states,n_traject_alloc,nproc))
      ! sv_ch_prob_l = 0d0
      ! ALLOCATE(      
      ! & sv_ch_what_computed
      ! & (num_ini_states,n_traject_alloc,nproc))	  
      ! sv_ch_what_computed  = 0
      ! ALLOCATE(sv_ch_def_fnc(num_ini_states,n_traject_alloc,nproc))
      ! sv_ch_def_fnc = 0
      ! if(allocated(sv_ch_vib_fnc)) deallocate(sv_ch_vib_fnc)		!Bikram Oct'18.	  
      ! ALLOCATE(sv_ch_vib_fnc(num_ini_states,n_traject_alloc,nproc))
      ! sv_ch_vib_fnc = 0	 
      ! ENDIF	  
! !!! ALL PROCESSORS ALLOCATE THESE BUFFER ARRAYS
 
      ! ALLOCATE(sv_ch_prob_l_buffer(
      ! & sys_var_size+6,num_ini_states,n_traject_alloc))
      ! sv_ch_prob_l_buffer = 0d0
      ! ALLOCATE(      
      ! & sv_ch_what_computed_buffer
      ! & (num_ini_states,n_traject_alloc))
      ! sv_ch_what_computed_buffer  = 0
      ! ALLOCATE(sv_ch_def_fnc_buffer(num_ini_states,n_traject_alloc))
      ! ALLOCATE(sv_ch_vib_fnc_buffer(num_ini_states,n_traject_alloc))	  
      ! sv_ch_def_fnc_buffer = 0
      ! sv_ch_vib_fnc_buffer = 0	  
      ! ALLOCATE(TRAJECT_DATA_CHECK_FILE(sys_var_size+6))
      ELSE
      CALL DEFINE_WIDTH(s_ini)	  
      ENDIF		  
!!!!! NUMBER OF TRAJECTORIES HAVE BEEN SETUP	  
      CALL ALLOCATE_ARRAYS
!!! ERRORS FOR EACH ENERGIES	  
      err_id_ener = 0d0
      err_id_prb = 0d0 	  
!!!!! DISTRIBUTE VALUE OF L	  
      IF(.not.monte_carlo_defined) THEN
      ALLOCATE(itraj_myid_l(2,tot_number_of_traject))	  
      l_counter = 0	  
      DO id_proc_root = 0,(nproc-1)/mpi_traject
!!!!!!!! ASSIGNING VALUE OF L
      DO itraject = 1,n_traject_alloc-1  !!! FIX AFTER ALL
!      l_counter = l_counter + 1
	  l_counter = id_proc_root+(itraject-1)*nproc/mpi_traject+1 		!Bikram Oct'19
!! IF YOU WANT TO SEE HOW ITS ASSIGNED put myid.eq.0		  
      IF(myid.eq.-4) THEN
      WRITE(*,*) l_counter,itraject,id_proc_root	  
      ENDIF	
!! ASSIGNING  proccesors and trajectories to values of L
      itraj_myid_l(1,l_counter) = itraject
      itraj_myid_l(2,l_counter) = id_proc_root
      DO i=1,mpi_traject
      id_proc = id_proc_root*mpi_traject+i
!      GOTO 3247	  
      CALL TRAJ_ORB(l_scatter(itraject,id_proc),
     & l_counter,
     & -1,
     & delta_l_step,
     & identical_particles_defined,
     & L_MIN_TRAJECT)	  
!3247      l_scatter(itraject,id_proc) = (l_counter-1)*delta_l_step+
!     & L_MIN_TRAJECT
!! TESTING	 
      IF(myid.eq.-4) THEN
      WRITE(*,*) l_scatter(itraject,id_proc),id_proc	  
      ENDIF
!! TESTING	ASSIGNEMET  
      ENDDO	 
      ENDDO
	  
      IF(id_proc_root
     & .lt. tot_number_of_traject -
     & (n_traject_alloc-1)*(nproc/mpi_traject)) THEN !!! FIX AFTER ALL
!      l_counter = l_counter + 1
!      l_counter = id_proc_root+(itraject-1)/mpi_traject*nproc+1 		!Bikram Oct 2019
      l_counter = l_counter+nproc/mpi_traject					 		!Bikram Dec 2019
	  temp_l_counter=l_counter								!Bikram Oct'19
	  if_l=.TRUE.											!Bikram
      itraj_myid_l(1,l_counter) = n_traject_alloc
      itraj_myid_l(2,l_counter) = id_proc_root
      DO i=1,mpi_traject
      id_proc = id_proc_root*mpi_traject+i
!      GOTO 3248	  
      CALL TRAJ_ORB(l_scatter(n_traject_alloc,id_proc),
     & l_counter,
     & -1,
     & delta_l_step,
     & identical_particles_defined,
     & L_MIN_TRAJECT)	  
!3248      l_scatter(n_traject_alloc,id_proc)=(l_counter-1)*delta_l_step+
!     & L_MIN_TRAJECT
      ENDDO		  
      ENDIF	 
!	  write(*,'(a,1x,i0,1x,i0,1x,i0,1x,i0)')'Lp',myid, l_counter, 
!     & itraj_myid_l(1,l_counter),itraj_myid_l(2,l_counter)
      ENDDO
	  if(if_l) l_counter=temp_l_counter			!Bikram
	  
      IF(l_counter.ne.tot_number_of_traject)
     & STOP "ERROR IN L ASSIGNMENT"
!!! TESTING ENDED
      ENDIF
      TIME_BUFFER=MPI_Wtime()
      IF(myid.eq.0) THEN
      WRITE(*,'(a31,1x,i5)')
     & "TIME IT TOOK TO SETUP ARRAYS, s",
     & int(TIME_BUFFER-TIME_PROPAGATE)	 
      ENDIF	  
!--------------------------------------------------------------------
! Loop over the two values of exchange parity
!--------------------------------------------------------------------
      DO p_parity = p_cur,ident_max
      IF(.not.monte_carlo_defined) THEN 								! NON_MONTECARLO LOOP
!--------------------------------------------------------------------
! This is loop over j12
!--------------------------------------------------------------------
      DO j_t = j_cur, j_max_ind(chann_ini)
      j_int_ini = j_t
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) THEN
      j_curr_f = j_t+0.5d0	  
      ENDIF
!--------------------------------------------------------------------
! This is loop over m12
!--------------------------------------------------------------------
      DO m_t = m_cur, j_t
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) THEN
      m_curr_f = m_cur+0.5d0
      IF(MYID.eq.0) WRITE(*,'(a27,1x,f4.1)')
     & "CALCULATIONS START FOR M = ",m_curr_f	  
      ELSE
      IF(MYID.eq.0) WRITE(*,'(a27,1x,i3)')
     & "CALCULATIONS START FOR M = ",m_t
      ENDIF	 
      IF(.not.identical_particles_defined) THEN
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) THEN
      s_st = indx_corr(round(m_curr_f +j_curr_f+1),
     & round(j_curr_f+0.5d0),chann_ini)  
      ELSE	  
      s_st = indx_corr(m_t + j_t+1,
     & j_t+1,chann_ini)
      ENDIF	 
      ELSE
      s_st = indx_corr_id(p_parity,m_t + j_t+1,
     & j_t+1,chann_ini) 	  
      ENDIF
!--------------------------------------------------------------------
! Checks whether the initial state is defined correctly or not
!--------------------------------------------------------------------
      IF(s_st.ne.0) THEN
!--------------------------------------------------------------------
! The piece below regarding CS-MQCT is not currently working
!--------------------------------------------------------------------
      IF(coupled_states_defined 
!!     & .and.(.not.mpi_task_per_proc)       !!! 531
     & .and. (.not.term_pot_defined) .and. .FALSE.) THEN 
      IF(MYID.eq.0 .and. s_st.eq.s_ini)	
     & PRINT*,"IMPROVED CS WILL BE USED" 
      CALL CS_MATRIX_IDENT	  
      ENDIF
!--------------------------------------------------------------------
! Setting up trajectories
!--------------------------------------------------------------------
      total_probab = 0d0
      sigma = 0d0
      probab_J = 0d0
      probab_J_all = 0d0
      idum = -myid-1
      ampl_wf_real = 0d0	  
      ampl_wf_imag = 0d0
      phase_wf = 0d0
      angle_scatter = 0d0
! Bikram Start:
      bk_tym = 0d0
      bk_prd = 0d0
      bk_vib = 0d0
      bk_erre = 0d0
      bk_errp = 0d0
! Bikram End.
	  
      IF(SPIN_FINE.ne.2) THEN	  
      j12_s_st = j12(s_st) !!! TO BE CONTINUED
      m12_s_st = m12(s_st)
      ENDIF
      IF(identical_particles_defined)THEN
      parity_st = parity_state(s_st)
      ENDIF	  
	  bk_parity = p_parity														!Bikram Nov 2020
!--------------------------------------------------------------------
! This is loop over trajectories
!--------------------------------------------------------------------
      DO itraject=st_traject,ntraject
	  temp_t_bfr_dotrjct=MPI_Wtime()
!--------------------------------------------------------------------
! Dotraject is the actual trajectory propagation subroutine
!--------------------------------------------------------------------
      CALL DOTRAJECT
      IF((myid/mpi_traject)*mpi_traject.eq.myid .and. .not.traj_run)
     & then
	  temp_t_afr_dotrjct=MPI_Wtime()
	  t_prop=temp_t_afr_dotrjct-temp_t_bfr_dotrjct
	  WRITE(*,'(a22,1x,i5,1x,a11,1x,f12.3,1x,a4)')
     & "TRAJECTORY DONE FOR L=",l_scatter(itraject,myid+1),
     & "TIME_SPENT=",t_prop,"SEC."
      endif	 
!--------------------------------------------------------------------
! Old Monte-Carlo checkpointing. Disabled by Bikram
!--------------------------------------------------------------------
      ! IF(write_check_file_defined) THEN!!!! SAVING THE DATA FROM TRAEJCTORY TO SAVE IT IN CHECK FILE
      ! sv_ch_prob_l_buffer
      ! & (:,s_st -ini_st_check_file + 1 ,itraject) = 
      ! & TRAJECT_DATA_CHECK_FILE
      ! sv_ch_def_fnc_buffer
      ! & (s_st -ini_st_check_file + 1 ,itraject) =  loop_count	
      ! sv_ch_vib_fnc_buffer
      ! & (s_st -ini_st_check_file + 1 ,itraject) =  vibration_cnt 
      ! IF(traj_run) THEN 
      ! sv_ch_what_computed_buffer
      ! & (s_st -ini_st_check_file + 1 ,itraject) = -1
      ! ELSE
      ! sv_ch_what_computed_buffer
      ! & (s_st -ini_st_check_file + 1 ,itraject) = 1
      ! !      IF(MYID.eq.12) PRINT*,       sv_ch_what_computed_buffer
      ! !     & (s_st -ini_st_check_file + 1 ,itraject),	itraject 
      ! ENDIF	  
      ! ENDIF
      CALL TIME_CHECKER(INTERUPT_PROPAGATION)
      IF(INTERUPT_PROPAGATION) THEN
      EXIT	  
      ENDIF	  

!!!   GATHERING PROBABILTIES FROM ALL TRAJECTORIES
      DO i=1,number_of_channels		  
      probab_J(i,itraject) = probab(i)
	  if(probab(i).ne.probab(i)) print *, 'probab_error', probab(i)
      ENDDO	  
!      IF(MYID.eq.1 .and. itraject.eq.1) PRINT*,	 probab !!!! DEELET AFER ALL	  
!!!!! FINDING THE LARGEST ERROR FOR MPI TASK	  
      IF(err_id_ener.lt.err_ener_tmp) THEN
      err_id_ener = err_ener_tmp
      ERROR_INDEX_ENERGY(myid+1) = itraject
      num_traj_ener = itraject	  
      ENDIF	  
      IF(err_id_prb.lt.err_prb_tmp) THEN
      err_id_prb = err_prb_tmp
      ERROR_INDEX_PROBAB(myid+1) = itraject
      num_traj_prb = itraject	  
      ENDIF 	  
!      WRITE(*,*) itraject,"err_id_prb,",err_id_prb	  
      ENDDO
!!!!! STOP TRAJECTORY LOOP
!!!! WAITING FOR ALL OTHER PROCESSORS	 
      TIME_DOTRAJECT = MPI_Wtime()
!!!! SAVING DATA IN CHECKPOINT
      CALL TIME_CHECKER(INTERUPT_PROPAGATION)
      IF(INTERUPT_PROPAGATION) THEN
      IF(MYID.EQ.0) PRINT*, 
     & "USER TIME LIMIT REACHED.", "CHECKPOINT IS BEING CREATED"
      task_size	= (sys_var_size+6)*n_traject_alloc*
     &	num_ini_states 	 
      CALL MPI_GATHER( sv_ch_prob_l_buffer,
     & task_size,MPI_DOUBLE_PRECISION,
     & sv_ch_prob_l,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
      IF(MYID.eq.0 .and. ierr_mpi.eq.0) 
     & PRINT*, "QUANTUM/CLASSICAL VARIABLES GATHERED"	 	 
      task_size	= n_traject_alloc*
     &	num_ini_states 	 
      CALL MPI_GATHER( sv_ch_what_computed_buffer,
     & task_size,MPI_INTEGER,
     & sv_ch_what_computed,
     &	  task_size,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
      IF(MYID.eq.0) PRINT*, "CHECK POINT FLAGS GATHERED"		 
      CALL MPI_GATHER( sv_ch_def_fnc_buffer,
     & task_size,MPI_INTEGER,
     & sv_ch_def_fnc,
     &	  task_size,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_GATHER( sv_ch_vib_fnc_buffer,
     & task_size,MPI_INTEGER,
     & sv_ch_vib_fnc,
     &	  task_size,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)	 
      IF(MYID.eq.0) PRINT*, " SCATTERING DATA GATHERED"		 
      IF(MYID.eq.0) PRINT*, "ALL DATA GATHERED"	 
      IF(MYID.EQ.0) THEN
      CALL PRINT_CHECKPOINT  
      ENDIF	  
      IF(MYID.eq.0) PRINT*, "CHEKPOINT DONE. RESTART THE PROGRAM"	 
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	 
      RETURN	  
      ENDIF
!!!!! GATHERING DATA. 
      task_size = number_of_channels*n_traject_alloc	  
	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF(MYID.eq.0) THEN
      IF(.not.identical_particles_defined) THEN	  
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) THEN
      WRITE(*,'(a36,1x,i3,1x,a9,1x,f4.1,1x,a6,1x,f4.1)')
     & "ALL TRAJECTORIES DONE FOR ENER_NUM =",i_ener,"AND M12 =",
     & m12_h(s_st)," J12= ", j12_h(s_st)	  
      ELSE	  
      WRITE(*,'(a36,1x,i3,1x,a9,1x,i4,1x,a6,1x,i4)')
     & "ALL TRAJECTORIES DONE FOR ENER_NUM =",i_ener,"AND M12 =",
     & m12(s_st)," J12= ", j12(s_st)
      ENDIF
      ELSE
      WRITE(*,'(a36,1x,i3,1x,a9,1x,i4,1x,a6,1x,i4,1x,a3,1x,i2)')
     & "ALL TRAJECTORIES DONE FOR ENER_NUM =",i_ener,"AND M12 =",
     & m12(s_st)," J12= ", j12(s_st), "P= ",p_parity	  
	  call flush(output_unit)
      ENDIF	  
      WRITE(*,'(a27,1x,i10)') "TIME IT TOOK TO PROPOGATE,s",
     & int(TIME_DOTRAJECT-TIME_BUFFER)  
	  call flush(output_unit)
      TIME_BUFFER = TIME_DOTRAJECT	 
      ENDIF
	 	  
!!!!   GATHERING DATA to determine the MAXIMUM ENERGY Error
      CALL MPI_GATHER(num_traj_ener,1,MPI_INTEGER,
     & ERROR_INDEX_ENERGY,
     &	  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)

      CALL MPI_GATHER(num_traj_prb,1,MPI_INTEGER,
     & ERROR_INDEX_PROBAB,
     &	  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
	  
      CALL MPI_GATHER(probab_J,task_size,MPI_DOUBLE_PRECISION,
     & probab_J_all,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
!! GATHERING DATA	 
      buffer_elastic = 	ampl_wf_real(:,myid+1) 
      CALL MPI_GATHER(buffer_elastic,
     & n_traject_alloc,MPI_DOUBLE_PRECISION,
     & ampl_wf_real,
     &	  n_traject_alloc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
	 
      buffer_elastic = 	ampl_wf_imag(:,myid+1) 
      CALL MPI_GATHER(buffer_elastic,
     & n_traject_alloc,MPI_DOUBLE_PRECISION,
     & ampl_wf_imag,
     &	  n_traject_alloc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
	 
      buffer_elastic = 	phase_wf(:,myid+1) 
      CALL MPI_GATHER(buffer_elastic,
     & n_traject_alloc,MPI_DOUBLE_PRECISION,
     & phase_wf,
     &	  n_traject_alloc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
	 
      buffer_elastic = 	angle_scatter(:,myid+1) 
      CALL MPI_GATHER(buffer_elastic,
     & n_traject_alloc,MPI_DOUBLE_PRECISION,
     & angle_scatter,
     &	  n_traject_alloc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)

! Bikram Start:
	  buffer_elastic = 	bk_tym(:,myid+1) 
      CALL MPI_GATHER(buffer_elastic,
     & n_traject_alloc,MPI_DOUBLE_PRECISION,
     & bk_tym,
     &	  n_traject_alloc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
	 
	  buffer_elastic = 	bk_prd(:,myid+1) 
      CALL MPI_GATHER(buffer_elastic,
     & n_traject_alloc,MPI_DOUBLE_PRECISION,
     & bk_prd,
     &	  n_traject_alloc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
	 
	  buffer_elastic = 	bk_vib(:,myid+1) 
      CALL MPI_GATHER(buffer_elastic,
     & n_traject_alloc,MPI_DOUBLE_PRECISION,
     & bk_vib,
     &	  n_traject_alloc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
	 
	  buffer_elastic = 	bk_erre(:,myid+1) 
      CALL MPI_GATHER(buffer_elastic,
     & n_traject_alloc,MPI_DOUBLE_PRECISION,
     & bk_erre,
     &	  n_traject_alloc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
	 
	  buffer_elastic = 	bk_errp(:,myid+1) 
      CALL MPI_GATHER(buffer_elastic,
     & n_traject_alloc,MPI_DOUBLE_PRECISION,
     & bk_errp,
     &	  n_traject_alloc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
! Bikram End.
	 
      buffer_elastic = 	angle_phase(:,myid+1) 
      CALL MPI_GATHER(buffer_elastic,
     & n_traject_alloc,MPI_DOUBLE_PRECISION,
     & angle_phase,
     &	  n_traject_alloc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)     
      CALL MPI_GATHER(err_id_ener,1,MPI_DOUBLE_PRECISION,
     & err_ener_larg,
     &	  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_GATHER(err_id_prb,1,MPI_DOUBLE_PRECISION,
     & err_prb_larg,
     &	  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)	 
!      PRINT*, "GATHER DONE ENDED"	 
!--------------------------------------------------------------------
! Broadcasting the trajectory information after gathering
!--------------------------------------------------------------------
      CALL MPI_BCAST(probab_J_all, task_size*nproc, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(angle_scatter, n_traject_alloc*nproc, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
! Bikram Start:
      CALL MPI_BCAST(bk_tym, n_traject_alloc*nproc, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
	  
	  CALL MPI_BCAST(bk_prd, n_traject_alloc*nproc, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
	 
	  CALL MPI_BCAST(bk_vib, n_traject_alloc*nproc, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
	 
	  CALL MPI_BCAST(bk_erre, n_traject_alloc*nproc, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
	  
	  CALL MPI_BCAST(bk_errp, n_traject_alloc*nproc, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
! Bikram End
!--------------------------------------------------------------------
! Computing the maximum energy error and norm conservation error
!--------------------------------------------------------------------	 
      IF(myid.eq.0) THEN
      DO i=1,nproc
      IF(error_ener_largest(i_ener).lt.err_ener_larg(i)) THEN
      error_ener_largest(i_ener)=err_ener_larg(i)
      err_proc_max_ener(i_ener) = i
      err_traj_max_ener(i_ener) = ERROR_INDEX_ENERGY(i)
      err_state_max_ener(i_ener) = s_st	  
      ENDIF	  
      ENDDO	
      DO i=1,nproc
      IF(error_prb_largest(i_ener).lt.err_prb_larg(i)) THEN
      error_prb_largest(i_ener)=err_prb_larg(i)
      err_proc_max_prb(i_ener) = i
      err_traj_max_prb(i_ener) = ERROR_INDEX_PROBAB(i)
      err_state_max_prb(i_ener) = s_st		  
      ENDIF
      ENDDO		  
      ENDIF	
!!!!  CHECK IF EVERYTHING IS ASSIGNED AND GATHERED CORRECTLY	  
      DO j=1,ntraject
      DO i=1,number_of_channels	  
      IF(probab_J_all(i,j,myid+1).ne.probab_J(i,j)) THEN
      PRINT*, "CRITICAL ERROR.MPI FINISH FOR MYID = ", myid
      PRINT*, "ERROR HERE= ", j,i
      PRINT*, "ERROR VALUES= ",probab_J_all(i,j,myid+1),probab_J(i,j) !!!! CHECKING		  
      MPIS_ERROR = 0
      ENDIF  	  
      ENDDO	  
      ENDDO
      CALL MPI_Reduce(MPIS_ERROR,MPIS_ERROR_TOT, 1,
     & MPI_INTEGER,MPI_PROD,0,
     & MPI_COMM_WORLD,ierr_mpi )
      CALL MPI_BCAST(MPIS_ERROR_TOT, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	 
       IF(MPIS_ERROR_TOT.ne.1)THEN
       IF(MYID.eq.0) PRINT*, "ERROR IN MPI PROBABILITY DISTRIBUTION"	   
       CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
       CALL MPI_FINALIZE (ierr_mpi)
       STOP	   
       ENDIF	   
!      PRINT*, "BROADCASTING ENDED"

      IF(MYID.eq.0) PRINT*, "DATA FROM ALL TRAJECTORIES GATHERED"
	  
! Bikram Start:
      IF(orbit_traj_defined) then 
      CALL PRINT_ORBITING_TRAJECTORY
      endif
! Bikram End
	  
      IF(calc_elast_defined) ang_res_elast = ang_res_num	  
      IF(calc_elast_defined .and. SPIN_FINE.ne.2) 						! RIGHT NOW FOR SPIN 1/2 NOT COMPUTE ELASTIC
!--------------------------------------------------------------------
! Computes Elastic integral and differential (if indicated) cross-section, 
! prints deflection function, and phase
!--------------------------------------------------------------------
     & CALL ELAST_CALC(sigma_elast(s_st),ang_res_elast,
     & deflect_fun_defined,diff_cross_defined,number_of_channels)
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) THEN
      IF(MYID.eq.0) CALL INELAST_CALC_FINE(sigma, !!! CALC INELASTIC CROSS SECTION
     & number_of_channels)	 
      IF(MYID.eq.0) THEN	 
      PRINT*,	"INELASTIC CROSS SECTION CALCULATED FOR M", m_curr_f
      PRINT*, " "
      ENDIF	  
      ELSE
!--------------------------------------------------------------------
! Computes the inelastic opacity functions
!--------------------------------------------------------------------
      IF(MYID.eq.0) CALL INELAST_CALC(sigma, number_of_channels)	 
      IF(MYID.eq.0) THEN	 
      PRINT*, "INELASTIC CROSS SECTION CALCULATED FOR M", m_t
      PRINT*, " "
      ENDIF	  
      ENDIF
!--------------------------------------------------------------------
! Implements Billing's correction to cross-sections and changing units to Ang^2
!--------------------------------------------------------------------
      DO j=1,number_of_channels
      sigma(j) = sigma(j)*a_bohr**2*pi*U(i_ener)/E_bill(j,i_ener)
!--------------------------------------------------------------------
! This is for idential colliding partners for inelastic cross-sections
!--------------------------------------------------------------------
      IF(identical_particles_defined) THEN
      SELECT CASE(coll_type)
      CASE(5)       	  
        sigma(j) = sigma(j)*
     & (delta(j1_ch(chann_ini),j2_ch(chann_ini))+1d0)*
     & (delta(j1_ch(j),j2_ch(j))+1d0)
      CASE(6)
        sigma(j) = sigma(j)*
     & (delta(j1_ch(chann_ini),j2_ch(chann_ini))
     & *delta(v1_ch(chann_ini),v2_ch(chann_ini))
     & +1d0)*
     & (delta(j1_ch(j),j2_ch(j))*delta(v1_ch(j),v2_ch(j))
     & +1d0)	  
      CASE(0)
        sigma(j) = sigma(j)*
     & (delta(j1_ch(chann_ini),j2_ch(chann_ini))*
     & delta(ka1_ch(chann_ini),ka2_ch(chann_ini))*
     & delta(kc1_ch(chann_ini),kc2_ch(chann_ini))+1d0)*
     & (delta(j1_ch(j),j2_ch(j))*
     & delta(ka1_ch(j),ka2_ch(j))*
     & delta(kc1_ch(j),kc2_ch(j))+1d0)	  
      END SELECT     	 
      ENDIF
!--------------------------------------------------------------------
! Converting unit of cross-sections A.U. to Ang^2
!--------------------------------------------------------------------
      IF(j.eq.chann_ini) sigma(j) = sigma_elast(s_st)*a_bohr**2
!--------------------------------------------------------------------
! For identical colliding partners Elastic cross-sections
!--------------------------------------------------------------------
      IF(identical_particles_defined.and.j.eq.chann_ini) THEN
      SELECT CASE(coll_type)	  
      CASE(5)       	  
        sigma(j) = sigma(j)*
     & (delta(j1_ch(j),j2_ch(j))+1d0)/2d0
      CASE(6)
        sigma(j) = sigma(j)*
     & (delta(j1_ch(j),j2_ch(j))*delta(v1_ch(j),v2_ch(j))
     & +1d0)/2d0	  
      CASE(0)
        sigma(j) = sigma(j)*
     & (delta(j1_ch(j),j2_ch(j))*
     & delta(ka1_ch(j),ka2_ch(j))*
     & delta(kc1_ch(j),kc2_ch(j))+1d0)/2d0
      END SELECT
      ENDIF	 
      ENDDO
!--------------------------------------------------------------------
! Averaging over the initial m12 states
!--------------------------------------------------------------------
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) THEN
      DO i=1,number_of_channels	  
      sigma_f(i,i_ener)	 =  sigma_f(i,i_ener) +
     & sigma(i)*2d0
     & /((j_max_ind(chann_ini)+j_min_ind(chann_ini)+2d0)  
     & *(j_max_ind(chann_ini)-j_min_ind(chann_ini)+1d0))  
      ENDDO	  
      ELSE
      DO i=1,number_of_channels	  
      sigma_f(i,i_ener)	 =  sigma_f(i,i_ener) +
     & sigma(i)*2d0/(1d0+delta(m_t,0))
     & /((j_max_ind(chann_ini)+j_min_ind(chann_ini)+1d0)  
     & *(j_max_ind(chann_ini)-j_min_ind(chann_ini)+1d0))
      ENDDO
      ENDIF	  
      ENDIF	  
  
      ENDDO
      ENDDO
      ELSE
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) STOP "NOT READY"
!--------------------------------------------------------------------
! Setting up trajectories for Monte-Carlo case
!--------------------------------------------------------------------  
      total_probab = 0d0
      sigma = 0d0
      probab_J = 0d0
      probab_J_all = 0d0
      CALL date_and_time(VALUES=timer_random)	  
      idum = (-myid-1)*timer_random(7)
      ampl_wf_real = 0d0	  
      ampl_wf_imag = 0d0
      phase_wf = 0d0
      angle_scatter = 0d0
! Bikram Start:
	  bk_tym = 0d0
	  bk_prd = 0d0
	  bk_vib = 0d0
	  bk_erre = 0d0
	  bk_errp = 0d0
! Bikram End.
      IF(myid.eq.0) PRINT*, "MONTE_CARLO INTEGRATION IS USED"	  

!--------------------------------------------------------------------
! Checking if Monte-Carlo sampling already exists in the file
!--------------------------------------------------------------------
	  old_nmb_traj = 0
	  total_traject_bikram = tot_number_of_traject
	  mc_traj_file_exst = .false.
	  inquire(file = 'MC_Initial_Cond.dat', 
     & exist = mc_traj_file_exst)
	  if(myid.eq.0 .and. .not. bikram_adiabatic) then	  
!--------------------------------------------------------------------
! Decides how many trajectories still need to be computed
!--------------------------------------------------------------------
      IF(b_impact_defined) THEN
      J_tot_max = int(b_impact_parameter*sqrt(massAB_M*2d0*E_sct))
     & -j_max_ind(chann_ini)	  														! DB change 8/23
      J_tot_min = 0
      ENDIF
      J_DOWN_INT = J_tot_min
      J_UP_INT  = 	J_tot_max
      
      L_MAX_TRAJECT = J_tot_max+j_max_ind(chann_ini)
      L_MIN_TRAJECT = L_MAX_TRAJECT	  
      DO i= J_tot_min, J_tot_max
      DO j= j_min_ind(chann_ini), j_max_ind(chann_ini)
      L_MIN_TRAJECT = min(L_MIN_TRAJECT,
     & abs(i-j) )      	  
      ENDDO		  
      ENDDO		  

!--------------------------------------------------------------------
! this piece is to compute total maximum number of trajectories by 
! taking into account the j12 and m12 sampling.
! For eaxh j12 we have a sampling of m12 from only zero to j12.
! In this way, the range of j12 cancells out. One can double check
! by enabling the commented out lines below.
!--------------------------------------------------------------------
      open(77, file = "tmp_avg.out")
	  tmp_avg = 0
	  write(77,*) 'j_max_ind(chann_ini)', j_max_ind(chann_ini)
	  do i = j_min_ind(chann_ini)+1, j_max_ind(chann_ini)+1
	    tmp_avg = tmp_avg + i
		write(77,*) 'j_max_ind(chann_ini)', j_max_ind(chann_ini)
	  end do
!	  tmp_avg = tmp_avg/(j_max_ind(chann_ini) - j_min_ind(chann_ini) + 1)
	  
	  total_traject_bikram = min(tot_number_of_traject,
     & ((L_MAX_TRAJECT - L_MIN_TRAJECT)/delta_l_step + 1)*tmp_avg)
!     &  (j_max_ind(chann_ini) - j_min_ind(chann_ini) + 1)*tmp_avg)

!--------------------------------------------------------------------
! Bikram Start Dec 2021: This is to read the trajectories 
! that were calculated in the previous run.
!--------------------------------------------------------------------
	  mc_traj_file_exst = .false.
	  inquire(file = 'MC_Initial_Cond.dat', 
     & exist = mc_traj_file_exst)
	  if(mc_traj_file_exst) then
	  open(13,file='MC_Initial_Cond.dat',action='read')
	  read(13, *)
	  old_nmb_traj = 0
	  do
	  read(13,*,iostat = io_traj_file)
	  if(io_traj_file/=0) exit
	  old_nmb_traj = old_nmb_traj + 1
	  end do
	  close(13)
	  
	  write(*,*)'#trajectories found from the previous run = ', 
     & old_nmb_traj
	  if(tot_number_of_traject.ge.old_nmb_traj) then 
      write(*,*) '#trajectories to be calculated in the current 
     & run = ', tot_number_of_traject - old_nmb_traj
	  else
	  write(*,'(a,a,a)')'You want to compute less #trajectories ',
     & 'than already computed in the previous run. ',
     & 'Please increase #trajectories and compute again.'
	  stop
	  end if
	  else
	  old_nmb_traj = 0
	  end if
	  
	  allocate(bk_s_st(tot_number_of_traject))
	  allocate(bk_l_real(tot_number_of_traject))
	  
	  if(mc_traj_file_exst) then
	  open(13,file='MC_Initial_Cond.dat',action='read')
	  read(13, *)
	  do itraject = 1, old_nmb_traj
	  read(13,*)bk_s_st(itraject), bk_l_real(itraject)
	  end do
	  close(13)
	  open(12,file='MC_Sample.dat',action='write',
     & access='append')
	  open(13,file='MC_Initial_Cond.dat',
     & action='write',access='append')
	  else	  
	  open(12,file='MC_Sample.dat',action='write')
	  write(12,'(6(a9),2(a19,2x))') "j_indx", "m_indx", "Parity", 
     & "j12(ini)", "m12(ini)", "l_real", "J_tot", "J_range"
	  open(13,file='MC_Initial_Cond.dat',
     & action='write')
	  write(13,'(2(a8))') "ini_st", "l_real"
	  end if
! Bikram End.
!--------------------------------------------------------------------
! Sampling of the initial conditions either new if this is the first run
! or the remaining piece after the previous run if trajectories are added
!-------------------------------------------------------------------- 
	  do itraject = old_nmb_traj+1, tot_number_of_traject
 192  rand1 = rand0(idum)
      rand2 = rand0(idum)
      rand_par = rand0(idum)
      p_monte_c  = p_lim_min	  
      if(rand_par.gt.0.5) p_monte_c = p_lim_max	  
      rand_j = rand0(idum)	  
      rand_m = rand0(idum)
      rand_j = rand_j*(j_max_ind(chann_ini)+1d0-j_min_ind(chann_ini))
      j_t = int(rand_j)+j_min_ind(chann_ini)
      if(rand_m.le.1d0/(2d0*j_t+1d0)) then
      m_t = 0
      else
      rand_m = rand_m-1d0/(2d0*j_t+1d0)
      rand_m = rand_m*(2d0*j_t+1d0)/2d0
      m_t = int(rand_m) + 1 	  
      end if
	  
	  IF(.not.identical_particles_defined) THEN
      s_st = indx_corr(m_t + j_t+1,
     & j_t+1,chann_ini) 
      ELSE
      s_st = indx_corr_id(p_monte_c,m_t + j_t+1,
     & j_t+1,chann_ini)
      ENDIF
      IF(s_st.le.0) STOP "s_st NOT DEFINED PROPERLY"
	  
	  j_ini_tr = j12(s_st)
	  if(b_impact_defined) then
      k_vec = sqrt(massAB_M*2d0*E_sct)
      J_tot = b_impact_parameter*rand1*k_vec
      dJ_int_range = b_impact_parameter*k_vec
      ELSE
      J_tot = dble(J_tot_max)*rand1 + dble(J_tot_min)*(1d0-rand1)
      dJ_int_range = dble(J_tot_max-J_tot_min)
      ENDIF
! "l_real" - ORBITAL MOMENTA SAMPLING  
      l_real = rand2*abs(J_tot-dble(j_ini_tr)) + 
     & (1d0-rand2)*(J_tot+dble(j_ini_tr))
      l_real = dble(round(l_real))	 
	  
	  mc_same_traj = .false.
	  bk_s_st(itraject) = s_st
	  bk_l_real(itraject) = l_real
	  if(itraject.gt.1) then
	  do mc_traj = 1, itraject-1
	  if(bk_s_st(mc_traj).eq.s_st 
     & .and. int(bk_l_real(mc_traj)).eq.int(l_real)) then
	  mc_same_traj = .true.
	  exit
	  end if
	  end do
	  end if
	  if(itraject.gt.total_traject_bikram) then
	  write(*,'(a,a)') 'The Monte-Carlo random sampling has reached ',
     &'maximum possible combination with current initial conditions.'
	  write(*,'(a,i0)')'Reducing total #trajectories = ', 
     & total_traject_bikram
	  exit
	  end if
	  if(mc_same_traj) go to 192
!--------------------------------------------------------------------
! Creates new file with current sampling for possible future use
!--------------------------------------------------------------------	  
	  write(12,'(6(i7,2x),2(e19.12,2x))') j_t, m_t, p_monte_c, 
     & int(j12(s_st)), int(m12(s_st)), int(l_real),J_tot,dJ_int_range
	  write(13,'(2(i6,2x))') s_st, int(l_real)
	  end do
	  close(12)
	  close(13)
	  end if
!--------------------------------------------------------------------
! Sampling information is broadcasted
!--------------------------------------------------------------------
	  call mpi_bcast(total_traject_bikram, 1, mpi_integer, 0,
     & MPI_COMM_WORLD, ierr_mpi)
	  call mpi_bcast(old_nmb_traj,1,mpi_integer,
     &  0, MPI_COMM_WORLD, ierr_mpi)
	  if(tot_number_of_traject.ge.total_traject_bikram) 
     & reduce_nmb_traj = .true.
	  tot_number_of_traject = total_traject_bikram
	  call MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )

! Bikram Start Dec 2021: Redistributing the trajectories among procs.
	  if(reduce_nmb_traj) then
      IF(mpi_task_defined) THEN	  
      n_traject_alloc = int(tot_number_of_traject/(nproc/mpi_traject))+1
      ELSE
      n_traject_alloc = int(tot_number_of_traject/nproc)+1	  
      ENDIF

      ntraject = n_traject_alloc-1  
      IF(mpi_task_defined) THEN		  
      IF(myid/mpi_traject
     & .lt. tot_number_of_traject - ntraject*(nproc/mpi_traject) )
     & ntraject = ntraject+1
      ELSE
      IF(myid .lt. tot_number_of_traject - ntraject*nproc )
     & ntraject = ntraject+1	  
      ENDIF
      IF(MYID.EQ.0) WRITE(*,*)
     & "New #trajectories per PROC_GROUP = ",ntraject
	  
	  if(allocated(probab_J)) deallocate(probab_J)
	  if(allocated(probab_J_all)) deallocate(probab_J_all)
      allocate(probab_J(number_of_channels,n_traject_alloc))		  
      allocate(probab_J_all(number_of_channels,n_traject_alloc,nproc))
	  end if
! Bikram End.
!--------------------------------------------------------------------
! This is a loop over Monte-Carlo trajectories
!--------------------------------------------------------------------
      DO itraject = st_traject, ntraject
	  
      IF(mpi_task_defined) THEN
      DO i=1,traject_roots		
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	 	    
!	  write(*,'(6(i5,2x))') myid, i, traject_roots, itraject, ntraject
!     & ,(itraject-1)*traject_roots+i
!	  print*, myid, (itraject-1)*traject_roots+i
	  open(12,file='MC_Sample.dat',action='read')
	  read(12, *)
	  
	  do ii = 1, (itraject-1)*traject_roots+i
	  read(12,'(6(i6,2x),2(e19.12,2x))') j_t, m_t, p_monte_c,
     & dmm1, dmm2, dmm3, J_tot_bk,dJ_int_range	 
	  	  
	  IF(.not.identical_particles_defined) THEN
      s_st = indx_corr(m_t + j_t+1,
     & j_t+1,chann_ini) 
      ELSE
      s_st = indx_corr_id(p_monte_c,m_t + j_t+1,
     & j_t+1,chann_ini)
      ENDIF
      IF(s_st.le.0) STOP "s_st NOT DEFINED PROPERLY"
	  	 
	  if(int(j12(s_st)).ne.dmm1) STOP "j12 isn't matching"
	  if(int(m12(s_st)).ne.dmm2) STOP "m12 isn't matching"
	  l_real = dmm3*1.0d0
	  l_real_bk = dmm3*1.0d0
	  bk_parity = p_monte_c
	  end do
	  close(12)
!	  write(*,'(9(i6,2x))') myid, (itraject-1)*traject_roots+i, j_t, 
!     & m_t, p_monte_c, s_st, int(j12(s_st)), int(m12(s_st)), 
!     & int(l_real)
      ENDIF	
      ENDDO		  
      ENDIF	 
! Bikram End.	  
	  
      IF(s_st.le.0) STOP "s_st NOT DEFINED PROPERLY"
	  temp_t_bfr_dotrjct=MPI_Wtime()
	    
! Bikram Start: Reading probabilities for the trajectories in case of addition of trajectories
	  if(mc_traj_file_exst) then
	  if((ii-1).gt.old_nmb_traj) then
      CALL DOTRAJECT
!      write(*,'(4(i0,1x))') myid, ii, old_nmb_traj, dmm3
	  else
!      write(*,'(4(i0,1x))') myid, ii-1, old_nmb_traj, dmm3
	  open(121,file='MC_Partial_XSections.out',action='read')
	  read(121,*) 
	  if(ii.gt.2) then
	  do i = 1, ii-2
	  read(121,*)
	  end do
	  end if
	  do j = 1, number_of_channels
	  read(121,'(e19.12, 2x)',advance='no') probab(j)
	  end do
	  close(121)  
	  end if
	  
	  else
      CALL DOTRAJECT
	  end if
! Bikram End.	  
	  
      IF((myid/mpi_traject)*mpi_traject.eq.myid .and. .not.traj_run)
     & then
	  temp_t_afr_dotrjct=MPI_Wtime()
	  t_prop=temp_t_afr_dotrjct-temp_t_bfr_dotrjct
	  WRITE(*,'(2(a,i5),a,i7,a,f10.3,a)')
     & "TRAJECTORY DONE FOR J12 = ",int(j12(s_st)), ", M12 = ", 
     & int(m12(s_st)), ", L = ", int(l_real),
     & ", TIME_SPENT = ",t_prop," SEC."
      endif	 
      CALL TIME_CHECKER(INTERUPT_PROPAGATION) 										! CHEKING IF WE SHOULD A STOP TRAJECTORY
      IF(INTERUPT_PROPAGATION) THEN
      id_trajec_mnt_crl = itraject	  
      EXIT
      ENDIF
!--------------------------------------------------------------------
! Collects probabilities and error information from Monte-Carlo trajecotries
!--------------------------------------------------------------------
      DO i=1,number_of_channels		  
      probab_J(i,itraject) = probab(i)
      ENDDO
      IF(err_id_ener.lt.err_ener_tmp) THEN
      err_id_ener = err_ener_tmp
      ERROR_INDEX_ENERGY(myid+1) = itraject
      num_traj_ener = itraject	  
      ENDIF	  
      IF(err_id_prb.lt.err_prb_tmp) THEN
      err_id_prb = err_prb_tmp
      ERROR_INDEX_PROBAB(myid+1) = itraject
      num_traj_prb = itraject	  
      ENDIF
      ENDDO
	  
! Bikram Start August 2021: Monte-Carlo AT-MQCT Calculations
	  if(bikram_save_traj) then
      write(bk_adia_dir1, '(a,a,f0.5)') 
     & trim(bk_directory), '/', U(i_ener)
	  bk_adia_dir2 = trim(bk_adia_dir1)
	  write(bk_adia_filepath,'(a,a,a)') './',  trim(bk_adia_dir2), '/'
	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  
	  if(myid.eq.0) then
	  call system("find " // bk_adia_dir2 // 
     & " -name 'adiabatic_traj*' > filename.out")
	  call system("mv filename.out " // bk_adia_dir2 )
	  end if 	
	  end if
! Bikram End.
	  
      task_size = number_of_channels*n_traject_alloc
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )

      IF(MYID.eq.0) WRITE(*,'(a43,1x,i3)')
     & "MONTE CARLO INTEGRATION DONE FOR ENER_NUM =",i_ener

      CALL MPI_GATHER(num_traj_ener,1,MPI_INTEGER,
     & ERROR_INDEX_ENERGY,
     &	  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
	 
      CALL MPI_GATHER(num_traj_prb,1,MPI_INTEGER,
     & ERROR_INDEX_PROBAB,
     &	  1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
	 
      CALL MPI_GATHER(probab_J,task_size,MPI_DOUBLE_PRECISION,
     & probab_J_all,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
	 
!!! GATHERING DATA   	 
      CALL MPI_GATHER(err_id_ener,1,MPI_DOUBLE_PRECISION,
     & err_ener_larg,
     &	  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_GATHER(err_id_prb,1,MPI_DOUBLE_PRECISION,
     & err_prb_larg,
     &	  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)	 

      CALL MPI_BCAST(probab_J_all, task_size*nproc, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      
      IF(myid.eq.0) THEN
      DO i=1,nproc
      IF(error_ener_largest(i_ener).lt.err_ener_larg(i)) THEN
      error_ener_largest(i_ener)=err_ener_larg(i)
      err_proc_max_ener(i_ener) = i
      err_traj_max_ener(i_ener) = ERROR_INDEX_ENERGY(i)
      err_state_max_ener(i_ener) = s_st	  
      ENDIF	  
      ENDDO	
      DO i=1,nproc
      IF(error_prb_largest(i_ener).lt.err_prb_larg(i)) THEN
      error_prb_largest(i_ener)=err_prb_larg(i)
!      PRINT*,"i",i,err_prb_larg(i)	  
      err_proc_max_prb(i_ener) = i
      err_traj_max_prb(i_ener) = ERROR_INDEX_PROBAB(i)
      err_state_max_prb(i_ener) = s_st		  
      ENDIF
      ENDDO		  
      ENDIF	
 
      DO j=1,ntraject
      DO i=1,number_of_channels	  
      IF(probab_J_all(i,j,myid+1).ne.probab_J(i,j)) THEN
      write(*,'(a,1x,3(i0,1x),2(f12.5,1x))') 
     & "CRITICAL ERROR.MPI FINISH FOR MYID = ", 
     & myid,i,j,probab_J_all(i,j,myid+1),probab_J(i,j)
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      CALL MPI_FINALIZE (ierr_mpi)
      ENDIF  	  
      ENDDO	  
      ENDDO
!      PRINT*, "BROADCASTING ENDED", myid
!--------------------------------------------------------------------
! Computes mean value, and variance for Monte-Carlo sampling
!--------------------------------------------------------------------
      i = 0
      IF(ALLOCATED(x2_mean_storage)) THEN
      DEALLOCATE(x2_mean_storage,x_mean_storage,
     & variance_values,mean_values)  
      ENDIF	  
      ALLOCATE(
     & x2_mean_storage(number_of_channels))
      ALLOCATE(
     & x_mean_storage(number_of_channels))
      x2_mean_storage = 0d0
      x_mean_storage = 0d0
      IF(check_point_defined .and. i_ener.eq.i_curr) THEN
      x2_mean_storage = x2_mean_storage_ch
      x_mean_storage = x_mean_storage_ch
      total_probab = total_probab_check_monte_carlo	  
      ENDIF
      IF(.not.check_point_defined)  tot_number_to_cut= 0	  
      DO i=1,number_of_channels
      x_mean_storage(i)  =x_mean_storage(i)+ probab_J_all(i,1,1)
      x2_mean_storage(i)  =x2_mean_storage(i)+ probab_J_all(i,1,1)**2      	  
      ENDDO	
	  
      ALLOCATE(mean_values(number_of_channels,tot_number_of_traject+
     & tot_number_to_cut ))
      ALLOCATE(variance_values
     & (number_of_channels,tot_number_of_traject-1+
     & tot_number_to_cut))
      mean_values = 0d0
      variance_values = 0d0
      i	 = 0
      buffer_mean = 0d0
      buffer_mean2 = 0d0	  
	  
! Bikram Start: Printing Monte-Carlo Informations. Nov 2021
	  if(.not. bikram_save_traj) then
	  if(myid.eq.0) open(1937, file = "MC_Traj_Info.out")
	  if(myid.eq.0) open(1938, file = "MC_Partial_XSections.out")
      DO j=1,number_of_channels
	  write(bk_label, '(a, i0)') "Part_XSections_", j
	  write(1937,'(a5,2x,a19)',advance='no') "Traj#", trim(bk_label)
	  write(1938,'(a19, 2x)',advance='no') trim(bk_label)
	  end do
	  write(1937,*)
	  write(1938,*)
	  end if
! Bikram End.
	  
      DO k=1,nproc
      id_proc = k-1	  
      IF(id_proc.ne.(id_proc/mpi_traject)*mpi_traject) CYCLE	  
      traject_sum_lim = n_traject_alloc - 1  
      IF(mpi_task_defined) THEN	
      IF(id_proc/mpi_traject.lt.tot_number_of_traject-traject_sum_lim
     & *(nproc/mpi_traject)) traject_sum_lim = n_traject_alloc
      ELSE
      IF(id_proc.lt.tot_number_of_traject-traject_sum_lim*nproc)
     & traject_sum_lim = n_traject_alloc  
      ENDIF
      ! IF(INTERUPT_PROPAGATION .and. write_check_file_defined) THEN
      ! traject_sum_lim = min(id_trajec_mnt_crl-1,traject_sum_lim)      
      ! ENDIF		  
	  
      DO itraject=1,traject_sum_lim
      i = i + 1	  
	  !print*, 'i', k, i
      DO j=1,number_of_channels
      total_probab(j) = total_probab(j) + probab_J_all(j,itraject,k)
      mean_values(j,i) = total_probab(j)/
     & dble(i+tot_number_to_cut)
      IF(i.gt.1) THEN	 
      variance_values(j,i-1) = variance_values(j,i-1) + 
     &  (probab_J_all(j,itraject,k)-mean_values(j,i))**2
      ENDIF	  

! Bikram Start: Printing Monte-Carlo Informations. Nov 2021	  
	  if(.not. bikram_save_traj) then
	  if(myid.eq.0) then
	  if(j.eq.1) then 
	  write(1937,'(i5, 2x, e19.12)',advance='no') i, total_probab(j)*
     & dJ_int_range*a_bohr**2*pi*U(i_ener)/E_bill(j,i_ener)/i 
	  else
	  write(1937,'(2x, e19.12)',advance='no') total_probab(j)*
     & dJ_int_range*a_bohr**2*pi*U(i_ener)/E_bill(j,i_ener)/i
	  end if
	  write(1938,'(e19.12, 2x)',advance='no')
     & probab_J_all(j,itraject,k)
	  end if
	  end if
	  
      ENDDO
	  if(.not. bikram_save_traj .and. myid.eq.0) then
	  write(1937,*)
	  write(1938,*)
	  end if
      ENDDO	  
      ENDDO
	  if(.not. bikram_save_traj .and. myid.eq.0) then
	  close(1937)
	  close(1938)
	  end if
! Bikram End.
	  
      IF(INTERUPT_PROPAGATION) THEN
      tot_number_of_traject_mc = i
      tot_number_of_traject = i	  
      CALL PRINT_CHECKPOINT_MC
      ELSE	
      IF(tot_number_of_traject.ne.i) STOP "ERROR IN INITIALIZATION"
      ENDIF	  
      IF(MYID.eq.0) PRINT*, "DATA FROM ALL TRAJECTORIES GATHERED"
      IF(MYID.eq.0) THEN
!      GOTO 124	  
      OPEN(UNIT=3,FILE="MONTE_CARLO_ERROR.out",POSiTION="APPEND")
      WRITE(3,*) "#ENERGY",i_ener      	  
      WRITE(3,'(a10,1x)',ADVANCE="NO")
     &	"N#_traject"
      DO i=1,number_of_channels
      WRITE(3,'(a10,1x,i6,1x,a3,1x,a5,1x)',ADVANCE = "NO")
     &	"#CHNL=",i,"ERR","VALUE"	  
      ENDDO
      WRITE(3,*)	  
      DO itraject = 1,tot_number_of_traject-1
      WRITE(3,'(i9,2x)',ADVANCE="NO") itraject
      DO i=1,number_of_channels
      WRITE(3,'(1x,f9.3,1x,e17.10)',ADVANCE="NO") 
     & dsqrt(variance_values(i,itraject))/
     & mean_values(i,itraject+1)*1d2/
     & dsqrt(dble(itraject+1+tot_number_to_cut)),
     & mean_values(i,itraject+1) 
      IF(itraject.eq.tot_number_of_traject-1 .and. 
     & tot_number_of_traject.gt.4) 
     & monte_carlo_err(i,i_ener) = 
     & MIDDLE_NUM(
     & dsqrt(variance_values(i,itraject-1))/
     & mean_values(i,itraject)*1d2/
     & dsqrt(dble(itraject+tot_number_to_cut)),
	 
     & dsqrt(variance_values(i,itraject))/
     & mean_values(i,itraject+1)*1d2/
     & dsqrt(dble(itraject+1+tot_number_to_cut)),
	 
     & dsqrt(variance_values(i,itraject-2))/
     & mean_values(i,itraject-1)*1d2/
     & dsqrt(dble(itraject-1+tot_number_to_cut)))	 
      ENDDO
      WRITE(3,*)	  
      ENDDO	  
      CLOSE(3)
      	  
c      PRINT*,	"dJ_int_range", dJ_int_range 
124      DO j=1,number_of_channels
!! MONTE CARLO INTEGRATION	 
      sigma(j) = total_probab(j)*dJ_int_range*
     & a_bohr**2*pi*U(i_ener)/E_bill(j,i_ener)
     & /tot_number_of_traject !!! COMMENT IF NOT INTEGRATED
      
      IF(identical_particles_defined) THEN
      SELECT CASE(coll_type)
      CASE(5)       	  
        sigma(j) = sigma(j)*
     & (delta(j1_ch(j),j2_ch(j))+1d0)
      CASE(6)
        sigma(j) = sigma(j)*
     & (delta(j1_ch(j),j2_ch(j))*delta(v1_ch(j),v2_ch(j))
     & +1d0)	 
      CASE(0)
        sigma(j) = sigma(j)*
     & (delta(j1_ch(j),j2_ch(j))*
     & delta(ka1_ch(j),ka2_ch(j))*
     & delta(kc1_ch(j),kc2_ch(j))+1d0)	  
      END SELECT     	 
      ENDIF	 
      IF(j.eq.chann_ini) sigma(j) = sigma_elast(s_st)*a_bohr**2
      IF(identical_particles_defined.and.j.eq.chann_ini) THEN
      SELECT CASE(coll_type)	  
      CASE(5)       	  
        sigma(j) = sigma(j)*
     & (delta(j1_ch(chann_ini),j2_ch(chann_ini))+1d0)*
     & (delta(j1_ch(j),j2_ch(j))+1d0)/4d0
      IF(p_parity.eq.1) sigma(j)=sigma(j)*exch_par_w_pl
      IF(p_parity.eq.2) sigma(j)=sigma(j)*exch_par_w_mn		 
      CASE(0)
        sigma(j) = sigma(j)*
     & (delta(j1_ch(chann_ini),j2_ch(chann_ini))*
     & delta(ka1_ch(chann_ini),ka2_ch(chann_ini))*
     & delta(kc1_ch(chann_ini),kc2_ch(chann_ini))+1d0)*
     & (delta(j1_ch(j),j2_ch(j))*
     & delta(ka1_ch(j),ka2_ch(j))*
     & delta(kc1_ch(j),kc2_ch(j))+1d0)/4d0
      IF(p_parity.eq.1) sigma(j)=sigma(j)*exch_par_w_pl
      IF(p_parity.eq.2) sigma(j)=sigma(j)*exch_par_w_mn	  
      END SELECT  
      ENDIF	 
      ENDDO	  
!       IF(MYID.eq.0) PRINT*, "SIGMA_J_DEFINED",sigma(2) 
      DO i=1,number_of_channels        
      sigma_f(i,i_ener)	 = sigma_f(i,i_ener) +  sigma(i)
      ENDDO

! Bikram July 09, 2023: This is to call the subroutine for new Monte-Carlo Error.
	  call MCerror(tot_number_of_traject, number_of_channels)
! Bikram end.
      ENDIF
      ENDIF
      ENDDO
	  
	  if(check_point_defined .or. write_check_file_defined) 
     & close(chk_file_unit)
!--------------------------------------------------------------------
! Finishing MPI task and printing trajecotry computing time
!--------------------------------------------------------------------
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      IF(MYID.eq.0)PRINT*,
     & "CROSS SECTION CALCULATED FOR ENERGY=",i_ener
      IF(MYID.eq.0)PRINT*," "
	  prop_time_end = mpi_wtime()
	  prop_time = prop_time_end - prop_time_bgn
	  if(myid.eq.0) print*, "Trajectory Propagation Time:", prop_time
      IF(MYID.eq.0)PRINT*," "

! Bikram Start:
	  if(myid.eq.0) then
	  OPEN(1111,FILE="CROSS_SECTIONS.tmp",POSiTION="APPEND")
      WRITE(1111,'(a,f14.6,a)') "U= ",U(i_ener)," cm^-1"
      WRITE(1111,'(2(a6,2x),a19,2x,a14,2x,a19)') 
     & 'ilv', 'flv', 'sigma(U),ANG^2',
     & 'E_coll,cm^-1', 'sigma(E_coll),ANG^2'
      DO i = 1,number_of_channels    	  
      IF(bill_exst(i,i_ener)) 
     & WRITE(1111,'(i6,2x,i6,2x,e19.10,2x,f14.6,2x,e19.10)') 
     & chann_ini, i, sigma_f(i,i_ener)/U(i_ener)*E_bill(i,i_ener), 
     & E_bill(i,i_ener), sigma_f(i,i_ener)
      ENDDO
! Bikram Start May 2020:	  
	  if(bk_nrg_err) then	  
      WRITE(1111,'(a28,1x,e12.5)') "TOTAL P CONSERVATION ERROR,%",
     & error_prb_largest(i_ener)
      WRITE(1111,'(a28,1x,e12.5)') "TOTAL E CONSERVATION ERROR,%",
     &	error_ener_largest(i_ener)
      WRITE(1111,*)
	  endif
! Bikram End.
	  close(1111)
	  endif
! Bikram End.
	  
      CALL DEALLOCATE_ARRAYS 
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      ENDDO
	  
	  return
      END SUBROUTINE PROPAGATE

      SUBROUTINE DOTRAJECT
! This subroutine is updated by Bikramaditya Mandal 
      USE VARIABLES
      USE CONSTANTS	  
      USE INTEGRATOR
      USE CS_MATRIX	  
      USE CHECK_FILE	  
      USE MPI_TASK_TRAJECT
      USE MPI_DATA
      USE MPI		  
	  use monte_carlo_sampling																!Bikram Sept 2021
	  use iso_fortran_env,only:output_unit
      IMPLICIT NONE
      LOGICAL BELONGS,fail_odeint
      LOGICAL :: orbiting_defined = .FALSE.	  
      INTEGER j_ini_tr,period_cnt!(Bikram)
      INTEGER i,j,k	,chann_num,k_st,ind_to_buff 
      INTEGER nok,nbad,i_root,dest,origin,unit_f
      INTEGER summ_range_min,summ_range_max,state_ch,traj_ch
      INTEGER st_sum_rn_max,st_sum_rn_min,ds_state	  
      INTEGER status(MPI_STATUS_SIZE)
	  INTEGER round
      REAL*8 J_tot,l_orb,l_real,teta,dir,dt,def_angle1,def_angle2
      REAL*8 randnum,eror,mom_r,tcur,pot,Eqi,Eqf,Ef,P_total,Eki	  
      REAL*8 rand0,delta,Mjmr,R_st,R_fin,dR_step,R_closest_app
      REAL*8 R,pr,ql,l,phi,qphi
      REAL*8 ort_r(3),ort_teta(3),ort_phi(3)
      REAL*8 velocity_ini(3),velocity_fin(3),scalpr,deflect_angle
      REAL*8 r_offset,vel_cross(3)
      REAL*8 buffer_pot,der_buffer_pot,cos_def,buffer_eq
      REAL*8 Mjmr_cs	
	  integer nmbr_r,nmbr_phi	  															!Bikram
      REAL*8 tmp_R2,tmp_R1,bk_tmp_r,bk_tmp_phi,bk_vel				   						!(Bikram) 
	  real*8 updated_R,updated_phi,f_time,int_1st_stp,bk_int_r		   						!(Bikram) 
	  
      REAL*8 q_cartes(6),buffer_sp(6)	  
!!! DELETE AFTER ALL	  
      COMPLEX*16, ALLOCATABLE:: ampl_inp(:),ampl_out(:),ampl_out2(:)	  
      EXTERNAL delta,Mjmr,scalpr,BELONGS,Mjmr_cs
      EXTERNAL DERIVS,RK4,RKQS,ODEINT,DERIVS_CS,DERIVS_CS_IM,
     & DERIVS_MPI,DERIVS_MPI_CS,BSSTEP,derivs_BK
      CHARACTER(LEN=10) :: istate_word = "imagi_    "  
      CHARACTER(LEN=10) :: rstate_word = "realp_    "	  									! DELETE AFTER ALL
      CHARACTER(LEN=18) :: probab_out = "  _PROBAB_OUT_.out"
      CHARACTER(LEN=16) :: traject_out_file = "TRAJECTORY  .out"
!Bikram Start August 2021:
	  integer dir_length, dm1, dm3
	  CHARACTER(LEN=100) fmt1, dm2, dm5
	  CHARACTER(LEN=:), allocatable :: fmt2
	  logical dm4
! Bikram End.
!Bikram Start April 2020:
      CHARACTER(LEN=100) bk_adia_print,bk_adia_filepath,bk_adia_print1
      CHARACTER(LEN=100) bk_adia_print2
      CHARACTER(LEN=100) bk_temp,bk_temp1, bkpr
      CHARACTER(LEN=100) bk_adia_dir1,bk_adia_dir2
	  logical bk_adia_file_exst, bk_b_chk, bk_adia_file_exst1
	  integer adia_io, ii, nnn
	  real*8 bikram_t
	  real*8, allocatable :: yp1(:), ypn(:), yprm(:)
	  real*8 bk_del, tmp_indx, iphase, rphase
	  REAL*8,allocatable :: bk_matt(:), bk_sin_coss(:,:), tt(:), tr(:)
	  real*8 magic, maxpot, mintm, maxtm, potbox, tmstp, fa, fb, tmchk
	  real*8 tol, bknorm, tau, v_temp, damp_coef, v_max, factor, tmp1
	  real*8 bk_probab, i_max, kin_E, tmp2, tmp3, tmp4
	  real*8, parameter :: pii = 4.0d0*atan(1.0d0)
	  integer :: chk_file_unit = 123123
	  integer stps, counter, stps_cntr, ph_cntr_bkk
	  logical advnc_rk4, at_dir_exst
	  external derivs_BK_adia
! Bikram End.

! Bikram Start May 2020:
	  counter = 0
	  if(bikram_adiabatic) then
	  call resize_adia
	  else
	  call resize
	  endif
	  
	  if(.not. bk_nrg_err .and. .not. bikram_print) then	  
	  else
	  allocate(bk_matt(mat_sz_bk))
	  ph_cntr_bkk = ph_cntr_bk
	  allocate(bk_sin_coss(2,ph_cntr_bk))
	  endif
	  bk_b_chk = .false.	
! Bikram End.

	  tmp_R1=0.d0
	  tmp_R2=0.d0
! "traj_run" - LOGICAL VARIABLE: DEFINES A CONDTION WHEN A TRAJECTORY IS PROPAGATED 
      traj_run = .TRUE.
! "period_cnt" - NUMBER OF PERIODS OF AN ORBITING TRAJECTORY   
      period_cnt = 0
      vibration_cnt = 0		!(Bikram)
!	  write(*,'(i0,1x,i0,1x,i0,1x,i0)') states_size, j12(s_st), 
!     & m12(s_st), parity_state(s_st)
!      dR_step = 0.5 ! NOT USED FOR NOW
! "sys_var"	- A REAL ARRAY WHICH CONTAINS POPULATION AMPLITUDES ai (SEE THE PAPER)
! AND CLASSICAL COORDINATES
! IT IS ORGANIZED AS FOLLOWS: FOR 1=<i=<states_size -  sys_var(i) = RE(ai)
! FOR states_size+1=<i+states_size<states_size*2 - sys_var(i+states_size) = IM(ai)
!(WHERE "states_size" - NUMBER OF QUANTUM STATES) 
! sys_var(1+2*states_size) - IT IS R - INTERMOLECULAR DISTANCE
! sys_var(2+2*states_size) - IT IS Pr - MOMENTA ALONG R**2
! sys_var(3+2*states_size) - IT IS THETA ANGLE
! sys_var(4+2*states_size) - IT IS Ptheta - MOMENTA CONJUGATE TO THETA
! sys_var(5+2*states_size) - IT IS PHI
! sys_var(6+2*states_size) - IT IS Pphi - MOMENTA CONJUGATE TO PHI 
! "probab" - IS AN ARRAY WHICH CONTAINS TRANSITION PROBABLITY, THE SIZE IS NUMBER OF CHANNELS	  
      probab = 0d0
! 	"deriv_sys" - FIRST TIME DERIVATIVE OF  "sys_var" 
      deriv_sys = 0d0
! "j_ini_tr" - THE TOTAL ROTATIONAL MOMENTA OF THE COLLISIONAL SYSTEm (NOT TOTAL ANGULAR)	  
      IF(SPIN_FINE.ne.2) j_ini_tr = j12(s_st)
! 	"k_vec" - THE INTIAL MOMENTA  , "E_sct" - SCATTERING ENRGY IN A.U.
! "massAB_M" - THE REDUCED MASS IN A.U. (ON INPUT - IN A.M.U.)
      k_vec = sqrt(massAB_M*2d0*E_sct)
! "dt" - TIME STEP IN INTEGRATION
! "U(i_ener)" - AN ARRAY WHICH CONTAINS KINETIC ENERGY,
! "ENER_REF" - A REFERENCE ENRGY (300 CM-1)
! "time_step" - TIME STEP ON INPUT 
      dt  = time_step/dsqrt(U(i_ener)/ENER_REF)
! "R_st" - INITISL POSITION DEFINED BY USER ON INPUT	  
      R_st = R_max_dist
! "b_impact_defined" - SAMPLING THROUHG IMPACT PARAMETER
! "rand1" - RANDOM NUMBER BETWEEN 0 AND 1
! SAMPLING FOR MONTE CARLO 
      IF(monte_carlo_defined) THEN 
	  
	  J_tot = J_tot_bk
	  l_real = l_real_bk
	  
!      IF(fine_structure_defined .and. SPIN_FINE.eq.2)STOP "NOT READY"	  
!      IF(b_impact_defined) THEN 
! "J_tot" - TOTAL ANGULAR MOMENTUM OF A TRAJECTORY 	  
!      J_tot = b_impact_parameter*rand1*k_vec
! "dJ_int_range" - MONTE_CARLO INTEGRATION RANGE
!      dJ_int_range = b_impact_parameter*k_vec
!      ELSE
! SAMPLING THROUG MAXIMUM AND MINIMAL TOTAL J	  
!      J_tot = dble(J_tot_max)*rand1 + dble(J_tot_min)*(1d0-rand1)
!      dJ_int_range = dble(J_tot_max-J_tot_min)
!      ENDIF
! "l_real" - ORBITAL MOMENTA SAMPLING  
!      l_real = rand2*abs(J_tot-dble(j_ini_tr)) + 
!     & (1d0-rand2)*(J_tot+dble(j_ini_tr))
!      l_real = dble(round(l_real))	 
      ELSE
      l_real = dble(l_scatter(itraject,myid+1))
      J_tot	 = l_real
      ENDIF
	  
! Bikram Start August 2021: Monte-Carlo testing
      if(monte_carlo_defined .and. myid.eq.-11) then
	  if(myid.eq.0) write(*,'(a)') 'Sampling for Monte-Carlo' 
	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  write(*,'(7(i6,2x),e19.12)') myid, itraject, bk_parity, s_st, 
     & int(j12(s_st)), int(m12(s_st)), int(l_real), J_tot
!	  write(*,'(6(i5,2x))') myid, itraject, int(l_real), int(J_tot), 
!     & int(j12(s_st)), int(m12(s_st))
	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  end if
! Bikram End.
	  
! "l_orb" - CLASSICAL MOMENTUM	  
      l_orb = dsqrt(l_real*(l_real+1d0))
      IF(mpi_task_defined) THEN
! IF MPI_TASK_PER_TRJAECTORY_DEFINED, BCASTING THE INTITIAL CONDITION
      DO i=1,traject_roots
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	  
      origin = 0		  
      CALL MPI_Bcast(J_tot, 1, MPI_REAL8, origin, comms(i),ierr_mpi)
      CALL MPI_Bcast(l_real, 1, MPI_REAL8, origin, comms(i),ierr_mpi)
      CALL MPI_Bcast(l_orb, 1, MPI_REAL8, origin, comms(i),ierr_mpi)
      CALL MPI_Bcast(s_st, 1, MPI_INTEGER, origin, comms(i),ierr_mpi)	  
      ENDIF		  
      ENDDO		  
      ENDIF	  
	  
! Bikram Start August 2021: Monte-Carlo testing
      if(monte_carlo_defined .and. myid.eq.-11) then
	  if(myid.eq.0) write(*,'(a)') 'Broadcasted for Monte-Carlo' 
	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  write(*,'(7(i6,2x),e19.12)') myid, itraject, bk_parity, s_st, 
     & int(j12(s_st)), int(m12(s_st)), int(l_real), J_tot
!	  write(*,'(6(i5,2x))') myid, itraject, int(l_real), int(J_tot), 
!     & int(j12(s_st)), int(m12(s_st))
	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  end if
! Bikram End.
	  
!Bikram Start: May 2020
! IF THE INITIAL SAMPLED POSITION HAPPENS TO BE OUTSIDE THE RANGE, MOVE THE PARTICLE 	  
      IF(l_orb/k_vec.gt.R_st) then 
	  write(*,'(a,i0,a)') 'WARNING: COLLISION IMPACT PARAMETER, b, 
     &EXCEEDS RMAX FOR TRAJECTORY L = ', int(l_real), 
     & '. THIS TRAJECTORY IS NOT PROPAGATED.'
	  bk_b_chk = .true.
!	  R_st = l_orb/k_vec + dt/2d0*k_vec
	  endif
!Bikram End.
! "R_fin" - WHERE TRAJECTORIES STOP	  
      R_fin =  R_st
! "R" - CURRENT R	  
      R = R_st
!	"R_closest_app" - THE DISTANCE OF THE CLOSEST APPROACH  
      R_closest_app = R_st
! "mom_r" - MOMENTUM ALONG R 
      mom_r = -sqrt(k_vec**2-l_orb**2/R_st**2)
! "teta" - 	ANGLE THETA
      teta =  dacos(0d0)
! "phi" - ANGLE PHI	  
      phi = 0d0
! "dir" - SCATTERING DIRECTION ANGLE. dir=0 SIMPLY MEANS THAT SCATTERING STARTS IN (PHI,R) PLANE	  
      dir = 0d0/2d0  !!!!!!!! WE CHANHING NOW!!!!!!!!!!!!!
! "tcur" - CURRENT TIME	  
      tcur = 0d0
! "pot" - THE VALUE OF THE POTENTIAL	  
      pot = 0d0
! "Eqi" - FULL QUANTUM ENERGY	  
      Eqi = 0d0	
! NOW "sys_var" IS INITIALIZNG
! INITIALLY ONLY THE INTIAL STATE "s_st" IS POPULATED
      sys_var = 0d0	  
      sys_var(s_st) = 1d0
	  buffer_eq = 0d0
	  
      IF(vib_mix_state_defined) THEN
      DO i=1,states_size	  
      chann_num = indx_chann(i)
      sys_var(i) = vib_mix_real(chann_num)
     & *delta(m12(i),m12(s_st))
      sys_var(i+states_size) = vib_mix_imag(chann_num)
     & *delta(m12(i),m12(s_st))
      buffer_eq = buffer_eq  
     & + sys_var(i+states_size)**2 
     & + sys_var(i)**2	 
      ENDDO	
      	sys_var= sys_var/sqrt(buffer_eq)  
      ENDIF
      buffer_eq = 0d0	  
      sys_var(1+states_size*2)=R_st
      sys_var(2+states_size*2)=mom_r
      sys_var(3+states_size*2)=teta
      sys_var(4+states_size*2)=l_orb*sin(dir)
      sys_var(5+states_size*2)=phi
      sys_var(6+states_size*2)=l_orb*cos(dir)
! "eror" - ENERGY CONSERVATION ERROR
      eror = 0d0
      current_error = 0d0
! ORTS	  
      ort_r(1) = dcos(sys_var(5+states_size*2))
     ^ *dsin(sys_var(3+states_size*2))
      ort_r(2) = dsin(sys_var(5+states_size*2))
     ^ *sin(sys_var(3+states_size*2)) 
      ort_r(3) = dcos(sys_var(3+states_size*2))
      ort_teta(1) = dcos(sys_var(5+states_size*2))
     ^ *dcos(sys_var(3+states_size*2))
      ort_teta(2) = dsin(sys_var(5+states_size*2))*
     ^ dcos(sys_var(3+states_size*2)) 
      ort_teta(3) = -dsin(sys_var(3+states_size*2))
      ort_phi(1) = -dsin(sys_var(5+states_size*2))
     % *dsin(sys_var(3+states_size*2))
      ort_phi(2) = dcos(sys_var(5+states_size*2))
     $ *dsin(sys_var(3+states_size*2)) 
      ort_phi(3) =0d0
! THE PROGRAM COPIES CLASSICAL VARIABLES FROM "sys_var" FOR FURTHER USE 	  
      R =  sys_var(1+states_size*2)
	  nmbr_r=1+states_size*2								!Bikram
      pr =  sys_var(2+states_size*2)
      ql = sys_var(3+states_size*2)
      l = sys_var(4+states_size*2)
      phi = sys_var(5+states_size*2)
	  nmbr_phi=5+states_size*2								!Bikram
      qphi = sys_var(6+states_size*2)
!!! SDF	  
      deflect_angle = 0d0
      def_angle1 = 0d0
      def_angle2 = 0d0
!!!! k_st identification
      time_st = 0d0
	  if(.not.bikram_mij_multiprint) then
      DO k=1,total_size
      IF(ind_mat(2,k).eq.s_st .and. s_st.eq.ind_mat(1,k)) THEN
	  k_st = k
      IF(coupled_states_defined .and. 
!     & (.not. mpi_task_per_proc) .and. !!! 1422
     & (.not. term_pot_defined)) THEN	  
      DO i=1,states_size_cs
	  IF(ind_state_cs(i).eq.s_st)ind_to_buff = i
      ENDDO
      ENDIF	  
	  EXIT
      ENDIF	  
      ENDDO	  
	  end if
! INTIAL VELOCITY	  
      DO i=1,3
      velocity_ini(i) = ort_r(i)*pr/massAB_M+(
     & ort_teta(i)*l/R/massAB_M
     & +ort_phi(i)*qphi/R/dsin(ql)**2/massAB_M)
      current_vel_ini(i) = velocity_ini(i)
      ENDDO	  
!   loops in case of orbiting
      loop_count = 0 
      ! IF(check_point_defined .and. .not. monte_carlo_defined .and. 
      ! & i_ener.eq. i_curr ) THEN
      ! state_ch = s_st-ini_st_check_file+1
      ! CALL TRAJ_ORB(INT(l_real),
      ! & traj_ch,
      ! & 1,
      ! & dl_step_integ,
      ! & ident_skip,
      ! & L_MIN_TRAJECT)	  
      ! SELECT CASE(what_computed(state_ch,traj_ch ))	  
      ! CASE(1)
      ! sys_var = all_prob_l(1:sys_var_size,state_ch,traj_ch)
      ! loop_count = all_def_fnc(state_ch,traj_ch)
	  ! vibration_cnt =  all_vib_fnc(state_ch,traj_ch)                   !!! DO NOT FORGET 
      ! tcur = all_prob_l(sys_var_size+1,state_ch,traj_ch)
      ! time_st = tcur	  
      ! eror = all_prob_l(sys_var_size+2,state_ch,traj_ch)
      ! current_error = all_prob_l(sys_var_size+3,state_ch,traj_ch)
      ! deflect_angle = all_prob_l(sys_var_size+4,state_ch,traj_ch)
      ! dt_corrected  = all_prob_l(sys_var_size+5,state_ch,traj_ch)
      ! dt_step  = all_prob_l(sys_var_size+6,state_ch,traj_ch)
      ! traj_run = .FALSE.
      ! CASE(-1)
      ! sys_var = all_prob_l(1:sys_var_size,state_ch,traj_ch)
      ! loop_count = all_def_fnc(state_ch,traj_ch)
	  ! vibration_cnt =  all_vib_fnc(state_ch,traj_ch)    !!! DO NOT FORGET 
      ! tcur = all_prob_l(sys_var_size+1,state_ch,traj_ch)
      ! time_st = tcur	  
      ! eror = all_prob_l(sys_var_size+2,state_ch,traj_ch)
      ! current_error = all_prob_l(sys_var_size+3,state_ch,traj_ch)
      ! deflect_angle = all_prob_l(sys_var_size+4,state_ch,traj_ch)
      ! dt_corrected  = all_prob_l(sys_var_size+5,state_ch,traj_ch)
      ! dt_step  = all_prob_l(sys_var_size+6,state_ch,traj_ch)
      ! def_angle2 = deflect_angle + 2d0*acos(-1d0)*loop_count	  
      ! !! REDEFINING ALL VALUES AGAIN	  
      ! ort_r(1) = dcos(sys_var(5+states_size*2))
      ! ^ *dsin(sys_var(3+states_size*2))
      ! ort_r(2) = dsin(sys_var(5+states_size*2))
      ! ^ *sin(sys_var(3+states_size*2)) 
      ! ort_r(3) = dcos(sys_var(3+states_size*2))
      ! ort_teta(1) = dcos(sys_var(5+states_size*2))
      ! ^ *dcos(sys_var(3+states_size*2))
      ! ort_teta(2) = dsin(sys_var(5+states_size*2))*
      ! ^ dcos(sys_var(3+states_size*2)) 
      ! ort_teta(3) = -dsin(sys_var(3+states_size*2))
      ! ort_phi(1) = -dsin(sys_var(5+states_size*2))
      ! % *dsin(sys_var(3+states_size*2))
      ! ort_phi(2) = dcos(sys_var(5+states_size*2))
      ! $ *dsin(sys_var(3+states_size*2)) 
      ! ort_phi(3) =0d0
      ! ! THE PROGRAM COPIES CLASSICAL VARIABLES FROM "sys_var" FOR FURTHER USE 	  
      ! R =  sys_var(1+states_size*2)
      ! pr =  sys_var(2+states_size*2)
      ! ql = sys_var(3+states_size*2)
      ! l = sys_var(4+states_size*2)
      ! phi = sys_var(5+states_size*2)
      ! qphi = sys_var(6+states_size*2)
      ! CASE(0)
      ! END SELECT
      ! ENDIF
!      IF MPI TASKS PER TRAJECORY IS USED, THE DRAGGING MASTER PROC SPREADS OUT THE INTIAL CONSITION TO SLAVES 	  
      IF(mpi_task_defined) THEN
      DO i=1,traject_roots	  
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	  
      origin = 0		  
      CALL MPI_Bcast(sys_var, 2*states_size+8,
     & MPI_REAL8, origin,comms(i),ierr_mpi)
      ENDIF		  
      ENDDO	 
      ENDIF
      DO i = 1,states_size !QUATUM ENERGY COMPUTING
      Eqi = Eqi + (sys_var(i)**2+sys_var(i+states_size)**2)*Ej(i)
      ENDDO	  

      IF(mpi_task_defined) THEN
      summ_range_min = 1!portion_of_MIJ_per_task(1,myid+1)
      summ_range_max = mat_sz_bk!portion_of_MIJ_per_task(2,myid+1)
      st_sum_rn_max = portion_of_state_per_task(2,myid+1) 
      st_sum_rn_min = 	portion_of_state_per_task(1,myid+1)  
      ELSE
      st_sum_rn_max = states_size 
      st_sum_rn_min = 	1 
      summ_range_min =1 
      summ_range_max = mat_sz_bk!total_size
      ENDIF	  
! Bikram Start May 2020:	 
	  if(.not. bk_nrg_err .and. .not. bikram_print) then	  
	  else
	  call potential_data_bk_adia(R,bk_matt)
	  do k= 1,ph_cntr_bkk
	  bk_sin_coss(1,k) = dsin(tcur*bk_delta_E(k))
	  bk_sin_coss(2,k) = dcos(tcur*bk_delta_E(k))
	  enddo
      DO k = summ_range_min,summ_range_max! POT ENERGY COMPUTING
      i = ind_mat_bk(1,k)!ind_mat(1,k)
      IF(coupled_states_defined .and. (m12(s_st) .ne. m12(i))) CYCLE	  
      j = ind_mat_bk(2,k)!ind_mat(2,k)
	  bk_del = 2.d0
	  if(i.eq.j) bk_del = 1.d0
	  tmp_indx = bk_indx(k)
      iphase = bk_sin_coss(1,tmp_indx)
      rphase = bk_sin_coss(2,tmp_indx)
       pot = pot + bk_del*bk_matt(k)*
     & ((sys_var(i)*sys_var(j) + 
     &   sys_var(i+states_size)*sys_var(j+states_size))*rphase -
     &  (sys_var(i+states_size)*sys_var(j) - 
     &   sys_var(i)*sys_var(j+states_size))*iphase)
      ENDDO
	  endif
! Bikram End.
!	  print *, 'chk1', myid  
!	  call flush(output_unit)
!	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
! BROADCASTING THE VALUE OF POTENTIAL	  
      IF(mpi_task_defined) THEN
      DO i=1,traject_roots	  
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	  
      origin = 0		  
      CALL MPI_Reduce(pot,buffer_pot, 1, MPI_REAL8,MPI_SUM,origin,
     & comms(i),ierr_mpi )
      pot = buffer_pot	 
      CALL MPI_Bcast(pot,1,
     & MPI_REAL8, origin, comms(i),ierr_mpi)	 
      ENDIF		  
      ENDDO		  
      ENDIF	  
      Eqi = Eqi + pot !  QUANTUM ENERGY IN THE BEGINING OF COLLISION
!! KINETIC ENERGY BEFORE COLLISION	  
      Eki = pr**2/2d0/massAB_M+l**2/2d0/massAB_M/R**2+	
     & qphi**2/2d0/massAB_M/R**2/sin(ql)**2
	  
! Bikram Start: April 2020, Adiabatic Trajectory Approximation
	  kin_E = pr**2/2d0/massAB_M
	  if(bikram_print) then
      write(bk_temp, '(a,i0,a,i0,a,i0,a,i0,a)') 
     & 'bk_aout_', int(i_ener), '_', int(j12(s_st)), 
     & '_',  int(m12(s_st)), '_', int(l_real), '.out'
	  bk_temp1 = "Bikram_Trajs"
	  call system ( "mkdir -p " // trim(bk_temp1) )
	  bk_temp1 = "./Bikram_Trajs/"
	  open(100002, file = trim(bk_temp1)//
     & trim(bk_temp), action = 'write')
	  endif
	  
! Bikram Start August 2021: Monte-Carlo Adiabatic Calculations
	  if(monte_carlo_defined) then
	  if(bikram_save_traj .and. .not. bk_b_chk) then
      write(bk_adia_dir1, '(a,a,f0.5)') trim(bk_directory), '/', 
     & U(i_ener)
	  bk_adia_dir2 = trim(bk_adia_dir1)
	  inquire(file = trim(bk_adia_dir2)//'/.', exist = at_dir_exst)
	  if(.not. at_dir_exst) then 
	  call system ( "mkdir -p " // bk_adia_dir2 )
	  else
	  endif
!	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  
	  write(bk_adia_filepath,'(a,a,a)') './',  trim(bk_adia_dir2), '/'
      write(bk_adia_print, '(a,i0,a,i0,a,i0,a,i0)') 
     & 'AT_', int(j12(s_st)), 
     & '_',  int(m12(s_st)), '_', int(l_real),
     & '_', bk_parity
      write(bk_adia_print1, '(a,i0,a,i0,a,i0,a,i0)') 
     & 'n_AT_', int(j12(s_st)), 
     & '_', int(m12(s_st)), '_', int(l_real),
     & '_', bk_parity 
!      write(bk_adia_print2,'(a)')'Monte_Carlo_info_adiabatic_traj.out'
	  
	  open(100003, file = trim(bk_adia_filepath)//
     & trim(bk_adia_print), action = 'write', form = 'unformatted')
	  open(100004, file = trim(bk_adia_filepath)//
     & trim(bk_adia_print1), action = 'write', form = 'unformatted')
	  
!	  open(100005, file = 'Monte_Carlo_info_adiabatic_traj.out', 
!     & access = 'APPEND')
!      write(100005,'(a)') bk_adia_print
!	  close(100005)	
!	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  
	  bk_adia_n = 0
	  endif
! Finishing up the Monte-Carlo part.	  
	  else
! Bikram End.
	  
	  if(bikram_save_traj .and. .not. bk_b_chk) then
      write(bk_adia_dir1, '(a,a,f0.5,a,i0,a,i0)') 
     & trim(bk_directory), '/', U(i_ener), '/', int(j12(s_st)), 
     & '/',  int(m12(s_st))
	  bk_adia_dir2 = trim(bk_adia_dir1)
	  inquire(file = trim(bk_adia_dir2)//'/AT_*', exist = at_dir_exst)
	  if(.not. at_dir_exst) then
	  call system ( "mkdir -p " // trim(bk_adia_dir2) )
	  else
	  endif
!	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  
	  write(bk_adia_filepath,'(a,a,a)') './',  trim(bk_adia_dir2), '/'
      write(bk_adia_print, '(a,i0,a,i0,a,i0,a,i0)') 
     & 'AT_', int(j12(s_st)), 
     & '_',  int(m12(s_st)), '_', int(l_real),
     & '_', bk_parity
      write(bk_adia_print1, '(a,i0,a,i0,a,i0,a,i0)') 
     & 'n_AT_', int(j12(s_st)), 
     & '_', int(m12(s_st)), '_', int(l_real),
     & '_', bk_parity 
	  
	  open(100003, file = trim(bk_adia_filepath)//
     & trim(bk_adia_print), action = 'write', form = 'unformatted')
	  open(100004, file = trim(bk_adia_filepath)//
     & trim(bk_adia_print1), action = 'write', form = 'unformatted')
	  bk_adia_n = 0

	  endif
	  endif

! Bikram Start August 2021: Monte-Carlo Adiabatic Calculations
	  if(monte_carlo_defined) then
	  if(bikram_adiabatic .and. .not. bk_b_chk) then
      write(bk_adia_dir1, '(a,a,f0.5)') 
     & trim(bk_directory), '/', U(i_ener)
	  bk_adia_dir2 = trim(bk_adia_dir1)
	  write(bk_adia_filepath,'(a,a,a)') './',  trim(bk_adia_dir2), '/'
	  
!	  if(myid.eq.0)call system("ls -p " // bk_adia_dir2 // 
!     & "> ls.out")
	  if(myid.eq.-10) then
!	  write(dm5,'(a,a)') trim(bk_adia_filepath),"filename.out"
!	  inquire(file=trim(dm5), exist = dm4)
!	  print *, dm4
!	  if(dm4) then
!	  call system("rm " // trim(dm5) )
!	  end if
	  call system("find " // bk_adia_dir2 // 
     & " -name 'adiabatic_traj*' > filename.out")
	  call system("mv filename.out " // bk_adia_dir2 )
	  end if
	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  
      write(bk_adia_print, '(a,i0,a,i0,a,i0,a,i0)') 
     & 'AT_', int(j12(s_st)), 
     & '_',  int(m12(s_st)), '_', int(l_real),
     & '_', bk_parity
      write(bk_adia_print1, '(a,i0,a,i0,a,i0,a,i0)') 
     & 'n_AT_', int(j12(s_st)), 
     & '_',  int(m12(s_st)), '_', int(l_real),
     & '_', bk_parity
  
	  if(myid.eq.-10) then
	  dir_length = len(trim(bk_adia_dir2)) + 1
	  write(fmt1,'(a,i0,a)') "'(a",dir_length,",a)'"
!	  fmt2 = trim(fmt1)
!	  print*, fmt2, fmt2
	  open(12345, file=trim(bk_adia_filepath)// "filename.out")
	  do dm1 = 1, (myid + 1)
	  read(12345, '(a)') dm2
	  end do
	  close(12345)
!	  print *, 'finished_reading'
	  
!	  write(*,'(i0,1x,a)') myid, dm2
	  fmt2 = dm2(1:dir_length)
!	  print*, 'fmt2', myid, fmt2
	  dm3 = len(trim(dm2)) + 1
	  bk_adia_print = dm2(dir_length+1:dm3)
      write(bk_adia_print1, '(a,a)') 'n_', trim(bk_adia_print)
!	  write(*,'(i0,1x,2(a,2x))') myid, bk_adia_print, bk_adia_print1
	  
! Putting the piece together when more than 1 proc is dedicated to a trajectory
	  IF(mpi_task_defined) THEN
      DO i=1,traject_roots	  
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	  
      origin = 0		  
	  print *, comms(i)
!      CALL MPI_Reduce(pot,buffer_pot, 1, MPI_REAL8,MPI_SUM,origin,
!     & comms(i),ierr_mpi )
!      pot = buffer_pot	 
!      CALL MPI_Bcast(pot,1,
!     & MPI_REAL8, origin, comms(i),ierr_mpi)	 
      ENDIF		  
      ENDDO		  
      ENDIF	  
	  print *, myid, bk_adia_print, bk_adia_print1
	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  stop
	  end if
! End.
	  
	  inquire(file=trim(bk_adia_filepath)//
     & trim(bk_adia_print), exist=bk_adia_file_exst)
	  inquire(file=trim(bk_adia_filepath)//
     & trim(bk_adia_print1), exist=bk_adia_file_exst1)
	  if(.not. bk_adia_file_exst) then
	  write(*,'(a,i0,a,i0,a,i0,a,i0)')
     & "The Trajectory File Does Not Exist; 
     & Please Check if it is Printed Properly. j12 = ",
     & int(j12(s_st)), " m12 = ",
     & int(m12(s_st)), " l = ", int(l_real),
     & " P = ", bk_parity
	  stop
	  call MPI_FINALIZE (ierr_mpi)
	  endif
	  if(.not. bk_adia_file_exst1) then
	  write(*,'(a,i0,a,i0,a,i0,a,i0)')
     & "The Trajectory File Size Does Not Exist; 
     & Please Check if it is Printed Properly. j12 = ",
     & int(j12(s_st)), " m12 = ",
     & int(m12(s_st)), " l = ", int(l_real),
     & " P = ", bk_parity
	  stop
	  call MPI_FINALIZE (ierr_mpi)
	  endif
	  
	  open(100004, file = trim(bk_adia_filepath)//
     & trim(bk_adia_print1), status = 'old', action = 'read', 
     & form = 'unformatted')
	  read(100004) bk_adia_n
	  close(100004)
	  
	  nnn = 7
	  allocate(bk_adia_t(bk_adia_n))
	  allocate(bk_sys_var(bk_adia_n,nnn))
	  allocate(bk_sys_var_der(bk_adia_n,nnn))
	  allocate(yp1(nnn), ypn(nnn), yprm(nnn))

	  open(100003, file = trim(bk_adia_filepath)//
     & trim(bk_adia_print), status = 'old', action = 'read', 
     & form = 'unformatted')
	  do i = 1, bk_adia_n
	  read(100003) bk_adia_t(i),bk_sys_var(i,1),bk_sys_var(i,2),
     & bk_sys_var(i,3),bk_sys_var(i,4),bk_sys_var(i,5),bk_sys_var(i,6),
     & bk_sys_var(i,7)
	  enddo
	  close(100003)
	  
	  do i = 1,nnn-1
	  yp1(i) = (bk_sys_var(2,i)-bk_sys_var(1,i))/
     & (bk_adia_t(2)-bk_adia_t(1))
	  ypn(i) = (bk_sys_var(bk_adia_n,i)-bk_sys_var(bk_adia_n-1,i))/
     & (bk_adia_t(bk_adia_n)-bk_adia_t(bk_adia_n-1))
	  call spline(bk_adia_t, bk_sys_var(:,i), bk_adia_n,
     & yp1(i), ypn(i), bk_sys_var_der(:,i))
	  call splint(bk_adia_t, bk_sys_var(:,i), bk_sys_var_der(:,i), 
     & bk_adia_n, tcur, sys_var(i+states_size*2), yprm(i))
	  enddo
	  
	  yp1(7) = (bk_sys_var(2,7)-bk_sys_var(1,7))/
     & (bk_adia_t(2)-bk_adia_t(1))
	  ypn(7) = (bk_sys_var(bk_adia_n,7)-bk_sys_var(bk_adia_n-1,7))/
     & (bk_adia_t(bk_adia_n)-bk_adia_t(bk_adia_n-1))
	  call spline(bk_adia_t, bk_sys_var(:,7), bk_adia_n,
     & yp1(7), ypn(7), bk_sys_var_der(:,7))
	  
	  if(rk4_defined .and. bk_step_size) then
	  maxpot = U(i_ener)/eVtown/autoeV
	  stps = numb_rk4_stps_adia
	  damp_coef = 20.0d0
	  tmp1 = -0.80d0
	  tmp2 = 4.250d0
	  tmp3 = (-1.0d0)* exp(2d0*pii*(U(i_ener)-202.d0)/damp_coef)
	  if(abs(tmp3).gt.1E16) tmp3 = -1E16
	  tmp4 = 0.30d0-tmp1*tmp3/(1.0d0-tmp3)-tmp2*tmp3/(1.0d0-tmp3)**2d0
      tol = bk_rk4_tol_adia*(U(i_ener)/ENER_REF)**tmp4
	  tau = 10.0d0
	  mintm = time_step
	  maxtm = bk_adia_t(bk_adia_n)
	  v_max = maxval(abs(bk_sys_var(:,7)))
	  magic = maxpot*((bk_adia_t(bk_adia_n)*tol)**(4d0))
	  endif
	  endif
! Finishing up the Monte-Carlo part.	  
	  else
! Bikram End.	  
	  if(bikram_adiabatic .and. .not. bk_b_chk) then
      write(bk_adia_dir1, '(a,a,f0.5,a,i0,a,i0)') 
     & trim(bk_directory), '/', U(i_ener), '/', int(j12(s_st)), 
     & '/',  int(m12(s_st))
	  bk_adia_dir2 = trim(bk_adia_dir1)
	  write(bk_adia_filepath,'(a,a,a)') './',  trim(bk_adia_dir2), '/'
      write(bk_adia_print, '(a,i0,a,i0,a,i0,a,i0)') 
     & 'AT_', int(j12(s_st)), 
     & '_',  int(m12(s_st)), '_', int(l_real),
     & '_', bk_parity
      write(bk_adia_print1, '(a,i0,a,i0,a,i0,a,i0)') 
     & 'n_AT_', int(j12(s_st)), 
     & '_',  int(m12(s_st)), '_', int(l_real),
     & '_', bk_parity
	  
	  inquire(file=trim(bk_adia_filepath)//
     & trim(bk_adia_print), exist=bk_adia_file_exst)
	  inquire(file=trim(bk_adia_filepath)//
     & trim(bk_adia_print1), exist=bk_adia_file_exst1)
	  if(.not. bk_adia_file_exst) then
	  write(*,'(a,i0,a,i0,a,i0,a,i0)')
     & "The Trajectory File Does Not Exist; 
     & Please Check if it is Printed Properly. j12 = ",
     & int(j12(s_st)), " m12 = ",
     & int(m12(s_st)), " l = ", int(l_real),
     & " P = ", bk_parity
	  stop
	  call MPI_FINALIZE (ierr_mpi)
	  endif
	  if(.not. bk_adia_file_exst1) then
	  write(*,'(a,i0,a,i0,a,i0,a,i0)')
     & "The Trajectory File Size Does Not Exist; 
     & Please Check if it is Printed Properly. j12 = ",
     & int(j12(s_st)), " m12 = ",
     & int(m12(s_st)), " l = ", int(l_real),
     & " P = ", bk_parity
	  stop
	  call MPI_FINALIZE (ierr_mpi)
	  endif
	  
	  open(100004, file = trim(bk_adia_filepath)//
     & trim(bk_adia_print1), status = 'old', action = 'read', 
     & form = 'unformatted')
	  read(100004) bk_adia_n
	  close(100004)
	  
	  nnn = 7
	  allocate(bk_adia_t(bk_adia_n))
	  allocate(bk_sys_var(bk_adia_n,nnn))
	  allocate(bk_sys_var_der(bk_adia_n,nnn))
	  allocate(yp1(nnn), ypn(nnn), yprm(nnn))

	  open(100003, file = trim(bk_adia_filepath)//
     & trim(bk_adia_print), status = 'old', action = 'read', 
     & form = 'unformatted')
	  do i = 1, bk_adia_n
	  read(100003) bk_adia_t(i),bk_sys_var(i,1),bk_sys_var(i,2),
     & bk_sys_var(i,3),bk_sys_var(i,4),bk_sys_var(i,5),bk_sys_var(i,6),
     & bk_sys_var(i,7)
	  enddo
	  close(100003)
	  
	  do i = 1,nnn-1
	  yp1(i) = (bk_sys_var(2,i)-bk_sys_var(1,i))/
     & (bk_adia_t(2)-bk_adia_t(1))
	  ypn(i) = (bk_sys_var(bk_adia_n,i)-bk_sys_var(bk_adia_n-1,i))/
     & (bk_adia_t(bk_adia_n)-bk_adia_t(bk_adia_n-1))
	  call spline(bk_adia_t, bk_sys_var(:,i), bk_adia_n,
     & yp1(i), ypn(i), bk_sys_var_der(:,i))
	  call splint(bk_adia_t, bk_sys_var(:,i), bk_sys_var_der(:,i), 
     & bk_adia_n, tcur, sys_var(i+states_size*2), yprm(i))
	  enddo
	  
	  yp1(7) = (bk_sys_var(2,7)-bk_sys_var(1,7))/
     & (bk_adia_t(2)-bk_adia_t(1))
	  ypn(7) = (bk_sys_var(bk_adia_n,7)-bk_sys_var(bk_adia_n-1,7))/
     & (bk_adia_t(bk_adia_n)-bk_adia_t(bk_adia_n-1))
	  call spline(bk_adia_t, bk_sys_var(:,7), bk_adia_n,
     & yp1(7), ypn(7), bk_sys_var_der(:,7))
	  
	  if(rk4_defined .and. bk_step_size) then
	  maxpot = U(i_ener)/eVtown/autoeV
	  stps = numb_rk4_stps_adia
	  damp_coef = 20.0d0
	  tmp1 = -0.80d0
	  tmp2 = 4.250d0
	  tmp3 = (-1.0d0)* exp(2d0*pii*(U(i_ener)-202.d0)/damp_coef)
	  if(abs(tmp3).gt.1E16) tmp3 = -1E16
	  tmp4 = 0.30d0-tmp1*tmp3/(1.0d0-tmp3)-tmp2*tmp3/(1.0d0-tmp3)**2d0
      tol = bk_rk4_tol_adia*(U(i_ener)/ENER_REF)**tmp4
	  tau = 10.0d0
	  mintm = time_step
	  maxtm = bk_adia_t(bk_adia_n)
	  v_max = maxval(abs(bk_sys_var(:,7)))
	  magic = maxpot*((bk_adia_t(bk_adia_n)*tol)**(4d0))
	  endif
	  endif
	  endif
	  
	  if(bk_b_chk) then 
	  traj_run = .false.
	  endif
! Bikram End.

! Bikram Start: June 2020
	  if(bikram_print) then
	  bknorm = 0.d0
	  do i = 1, states_size
	  bknorm = bknorm +(sys_var(i)**2+sys_var(i+states_size)**2)
	  enddo
	  bk_probab = 0.d0
	  do i = 1, states_size
      chann_num = indx_chann(i)
      if(chann_num.eq.chann_ini) cycle
	  bk_probab = bk_probab +(sys_var(i)**2+sys_var(i+states_size)**2)
	  enddo
	  write(100002,'(6(e23.16,1x))')
     & sys_var(1+states_size*2)*dcos(sys_var(5+states_size*2)),
     & sys_var(1+states_size*2)*dsin(sys_var(5+states_size*2)),
     & sys_var(1+states_size*2),sys_var(5+states_size*2),
     & sys_var(2+states_size*2),sys_var(6+states_size*2)	  
!	  write(100002,'(f12.5,1x,e23.16)') tcur,bknorm
!	  do i = 1, states_size
!	  if(i.eq.1) then
!	  write(100002,'(f12.5,1x,e23.16,1x,
!!     & e23.16,1x,e23.16,1x,
!     & e23.16)', advance = 'no')
!     & tcur,! bknorm - 1.0d0, bk_probab, dt,! pot*autoeV*eVtown, 
!     & sys_var(i), sys_var(i+states_size)!, abs(1d0 - P_total)*1d2, sys_var(1+states_size*2), tmp_R1, tmp_R2
!!    & sys_var(1+states_size*2),sys_var(5+states_size*2),
!	  else
!	  write(100002,'(1x,e23.16,1x,e23.16)', advance = 'no') 
!     & sys_var(i), sys_var(i+states_size)
!	  endif
!	  enddo
!	  write(100002,*)
!	  do i = 1, states_size
!	  if(i.eq.1) then
!	  write(100002,'(f12.5,1x,e23.16)', advance = 'no')
!     & tcur, sys_var(i)**2 + sys_var(i+states_size)**2
!	  else
!	  write(100002,'(1x,e23.16)', advance = 'no') 
!     & sys_var(i)**2 + sys_var(i+states_size)**2
!	  endif
!	  enddo
!	  write(100002,*)
	  
	  endif
! Bikram End.
	  
! TRAJECTORY STARTED
      DO WHILE(traj_run)
	  counter = counter + 1
! BROADCASTING SYS_VAR	  
      IF(mpi_task_defined) THEN
      DO i=1,traject_roots	  
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	   
      origin = 0		  
      CALL MPI_Bcast(sys_var, 2*states_size+8,
     & MPI_REAL8, origin, comms(i),ierr_mpi)
      ENDIF
      ENDDO
      ENDIF
! "ort_r","ort_teta","ort_phi" - DIRECTION OF VELOCITY IN CARTESIAN COORDIANTES	  
      ort_r(1) = dcos(sys_var(5+states_size*2))
     ^ *dsin(sys_var(3+states_size*2))
      ort_r(2) = dsin(sys_var(5+states_size*2))
     ^ *sin(sys_var(3+states_size*2)) 
      ort_r(3) = dcos(sys_var(3+states_size*2))
      ort_teta(1) = dcos(sys_var(5+states_size*2))
     ^ *dcos(sys_var(3+states_size*2))
      ort_teta(2) = dsin(sys_var(5+states_size*2))*
     ^ dcos(sys_var(3+states_size*2)) 
      ort_teta(3) = -dsin(sys_var(3+states_size*2))
      ort_phi(1) = -dsin(sys_var(5+states_size*2))
     % *dsin(sys_var(3+states_size*2))
      ort_phi(2) = dcos(sys_var(5+states_size*2))
     $ *dsin(sys_var(3+states_size*2)) 
      ort_phi(3) =0d0
! THE PROGRAM COPIES CLASSICAL VARIABLES FROM "sys_var" FOR FURTHER USE 	  
      R =  sys_var(1+states_size*2)
      pr =  sys_var(2+states_size*2)
      ql = sys_var(3+states_size*2)
      l = sys_var(4+states_size*2)
      phi = sys_var(5+states_size*2)
      qphi = sys_var(6+states_size*2)
      IF(R.ne.R) THEN
      PRINT*,"R NOT DEFINED FOR ID,traject",myid,itraject	  
      RETURN
      ENDIF	  
      IF(pr.ne.pr) THEN
      PRINT*,"PR NOT DEFINED FOR ID,traject",myid,itraject	  
      RETURN
      ENDIF	  
      IF(ql.ne.ql) THEN
      PRINT*,"QL NOT DEFINED FOR ID,traject",myid,itraject	  
      RETURN
      ENDIF	  
      IF(l.ne.l) THEN
      PRINT*,"L NOT DEFINED FOR ID,traject",myid,itraject	  
      RETURN
      ENDIF	  
      IF(phi.ne.phi) THEN
      PRINT*,"PHI NOT DEFINED FOR ID,traject",myid,itraject	  
      RETURN
      ENDIF	  
! "velocity_ini" - THE INITIAL VELOCITY	
!! "current_vel_ini" - THE CURRENT VELOCITY 
      DO i=1,3
      current_vel_ini(i) = ort_r(i)*pr/massAB_M+(
     & ort_teta(i)*l/R/massAB_M
     & +ort_phi(i)*qphi/R/dsin(ql)**2/massAB_M)
      ENDDO	  
!!!  COMPUTING DEFLECTION ANGLE ALONG THE TRAJECTORY
      IF(l_real.ne.0d0) THEN
!!!  DETREMINTAION OF ACTUAL ANGLE NOT TAKING ACCOUNT THE PHASE	  
      def_angle1	 =  def_angle2 
      CALL vect_product(vel_cross,current_vel_ini,velocity_ini)
       cos_def = scalpr(current_vel_ini,velocity_ini)/  
     & dsqrt(scalpr(current_vel_ini,current_vel_ini))/
     & dsqrt(scalpr(velocity_ini,velocity_ini))
      IF(abs(cos_def).ge.1) THEN
      def_angle2 = (acos(-1d0)-sign(acos(-1d0),cos_def))/2d0	  
      ELSE	  
      def_angle2 = dacos(cos_def)
      ENDIF	  
      IF(vel_cross(3).lt.0d0) THEN
      def_angle2 = -def_angle2
      ENDIF
!!!! ADDING LOOP IN ATTRATCION REGION	  
      IF(def_angle2*def_angle1.lt.0d0) THEN
      IF(abs(def_angle2-def_angle1).gt.acos(-1d0)) THEN
      loop_count = loop_count + 1	  
      ENDIF	  
      ENDIF
!!! IF "loop_count">0 THE TRAJECTORY IN ATTRATCION REGION	  
       deflect_angle = def_angle2 - 2d0*acos(-1d0)*loop_count	  
      ELSE
!! HEAD ON COLLISION	  
      vel_cross = 0d0
      deflect_angle = dacos(-1d0)	  
      ENDIF	  


     	 
!!!
! "Eqf" - THE QUANTUM ENERGY DURING THE COLLISION	  
      Eqf = 0d0	  
      DO i = st_sum_rn_min,st_sum_rn_max
      IF(coupled_states_defined .and. (m12(s_st) .ne. m12(i))) CYCLE	  
      Eqf = Eqf + (sys_var(i)**2+sys_var(i+states_size)**2)*Ej(i)
      ENDDO
      IF(mpi_task_defined) THEN
      DO i=1,traject_roots	  
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	  
      origin = 0		  
      CALL MPI_Reduce(Eqf,buffer_eq, 1, MPI_REAL8,MPI_SUM,origin,
     & comms(i),ierr_mpi )
      Eqf = buffer_eq	 
      CALL MPI_Bcast(Eqf,1,
     & MPI_REAL8, origin, comms(i),ierr_mpi)	 
      ENDIF		  
      ENDDO		  
      ENDIF		  
!	"Ef" - KINETIC ENERGY DURING THE COLLISION	  
      Ef = pr**2/2d0/massAB_M+l**2/2d0/massAB_M/R**2+	
     & qphi**2/2d0/massAB_M/R**2/sin(ql)**2	
! RECOMPUTING "pot" - THE POTENTIAL ENERGY
      pot = 0d0	
      IF(coupled_states_defined .and. 
!     & (.not. mpi_task_per_proc) .and. !!! 1670
     & (.not. term_pot_defined)) THEN
      DO i=1,states_size_cs
	  
	  DO j=1,i
      sys_var_cs(i) = sys_var(ind_state_cs(i))
      IF(ind_state_cs(i).gt.states_size) STOP "CS FAILED"	  
      sys_var_cs(i+states_size_cs) = 
     & sys_var(ind_state_cs(i)+states_size)	
      sys_var_cs(j) = sys_var(ind_state_cs(j))
      IF(ind_state_cs(j).gt.states_size) STOP "CS FAILED"	  
      sys_var_cs(j+states_size_cs) = 
     & sys_var(ind_state_cs(j)+states_size)	
	 
      pot = pot + 2d0*Mjmr_cs(R,i,j)/(1+delta(i,j))*
     ^   ((sys_var_cs(i)*(sys_var_cs(j)) + 
     ^ sys_var_cs(i+states_size_cs)*(sys_var_cs(j+states_size_cs)))*
     & dcos(tcur*(Ej_cs_im(j)-Ej_cs_im(i)))
     & - (sys_var_cs(i+states_size_cs)*(sys_var_cs(j)) - 
     ^ sys_var_cs(i)*(sys_var_cs(j+states_size_cs)))*
     & dsin(tcur*(Ej_cs_im(j)-Ej_cs_im(i))))	  

      ENDDO
      ENDDO	  
	  
	  ELSE

! Bikram Start May 2020:	  
	  if(.not. bk_nrg_err .and. .not. bikram_print) then	  
	  else
	  call potential_data_bk_adia(R,bk_matt)
	  do k= 1,ph_cntr_bkk
	  bk_sin_coss(1,k) = dsin(tcur*bk_delta_E(k))
	  bk_sin_coss(2,k) = dcos(tcur*bk_delta_E(k))
	  enddo
      DO k = summ_range_min,summ_range_max! POT ENERGY COMPUTING
      i = ind_mat_bk(1,k)!ind_mat(1,k)
      IF(coupled_states_defined .and. (m12(s_st) .ne. m12(i))) CYCLE	  
      j = ind_mat_bk(2,k)!ind_mat(2,k)
	  bk_del = 2.d0
	  if(i.eq.j) bk_del = 1.d0
	  tmp_indx = bk_indx(k)
      iphase = bk_sin_coss(1,tmp_indx)
      rphase = bk_sin_coss(2,tmp_indx)
       pot = pot + bk_del*bk_matt(k)*
     & ((sys_var(i)*sys_var(j) + 
     &   sys_var(i+states_size)*sys_var(j+states_size))*rphase -
     &  (sys_var(i+states_size)*sys_var(j) - 
     &   sys_var(i)*sys_var(j+states_size))*iphase)
      ENDDO
	  endif
! Bikram End.
	  
      IF(mpi_task_defined) THEN
      DO i=1,traject_roots	  
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	 
      origin = 0		  
      CALL MPI_Reduce(pot,buffer_pot, 1, MPI_REAL8,MPI_SUM,origin,
     & comms(i),ierr_mpi )
      pot = buffer_pot
      CALL MPI_Bcast(pot,1,
     & MPI_REAL8, origin, comms(i),ierr_mpi)	 
      ENDIF	  
      ENDDO
      ENDIF
      ENDIF	  
	  
! Bikram Start May 2020:

	  if(bikram_save_traj .and. .not. bk_b_chk) then	  
	  write(100003) tcur,
     & sys_var(1+states_size*2), sys_var(2+states_size*2), 
     & sys_var(3+states_size*2), sys_var(4+states_size*2), 
     & sys_var(5+states_size*2), sys_var(6+states_size*2), pot
	  bk_adia_n = bk_adia_n + 1
	  endif
! Bikram End.
! PROBABILITY CONSERVATION
      P_total = 0d0	  
      DO i = 1,states_size
! "P_total" MUST BE 1. HOWEVER, THERE IS A PROPAGATION ERROR	 
      P_total=P_total+(sys_var(i)**2+sys_var(i+states_size)**2)
      ENDDO
      IF(current_error.lt. abs(1d0 - P_total)) current_error
     &	  = abs(1d0 - P_total)	  
! COMPUTING THE ERROR IN ENRGY CONSERVATION "eror" 

      IF( abs(Eqf+Ef+pot-Eki-Eqi)>eror) eror = 
     & abs(Eqf+Ef+pot-Eki-Eqi)
      IF(R<R_closest_app) R_closest_app = R	! RECOMPUTING THE DISTANCE OF THE CLOSEST APPRAOCH 
! HERE THE MOST TIME CONSUMING PART - HERE A PROPAGATOR IS CALLED
! "rk4" - RUNGE - KUTTA 4TH ORDER
! "odeint" - ADAPTIVE STEP - SIZE CONTROL RUNGE-KUTTA
! "derivs" - CC MQCT, COMPUTES DERIVATIVES OF "sys_var" AND STORES THEM IN "deriv_sys"
! "derivs_cs" - CS MQCT, DOES THE SAME
! CHECKING 
      IF(tcur.eq.0d0) THEN
      dt_step = dt
      dt_corrected = dt	  
      ENDIF	
!!!!! PRINTING OUT THE TRAJECTORY DURING PROPAGATION
!!!!!!! OUTPUT	  
      IF(int(l_real).eq.0 .and. myid.eq.-2) THEN
!      IF(myid.eq.1 .and. itraject.eq. 1) THEN
      unit_f = 2
      IF(myid.eq.1) unit_f=3
      IF(tcur.eq.0)
     & WRITE(unit_f,*)!'(a12,2x,a16)')!2x,a12,2x,a12,2x,a12,2x,a12)')
     & "R","      dt_corrected","      eror"!,"pot","dt_step","tcur"	 
      WRITE
     & (unit_f,'(e17.10,1x,e17.10,1x,e17.10,1x,e17.10)') 	  
!     &  (unit_f,'(f12.5,2x,e12.5,2x,e12.5,2x,e12.5,2x,e12.5,2x,e12.5)')
     & R, pr**2/2d0/massAB_M,l**2/2d0/massAB_M/R**2,	
     & qphi**2/2d0/massAB_M/R**2/sin(ql)**2
      ENDIF			  
! PRINTING ENDED

!!!!

!!!! SWITCH TO CARTESIAN COORDIANTES
      IF(cartes_defined .and. l_real.eq.0d0) THEN
	  head_on = .TRUE.
      buffer_sp =  sys_var(1+states_size*2:6+states_size*2)
	  q_cartes = 0d0
	  q_cartes(1) = R*cos(phi)!*sin(ql) ! x,y,z, px,py,pz
	  q_cartes(2) =0d0 ! R*sin(phi)*sin(ql) ! y
	  q_cartes(3) = 0d0!R*cos(ql)	 !! z
      q_cartes(4) = pr/massAB_M*cos(phi) !current_vel_ini
      sys_var(1+states_size*2:6+states_size*2) = q_cartes
!!!! MAKE MAGIC HERE  !!!! DO 
      ENDIF
!!!!!

! Bikram start Oct 2019
      IF(prn_l_defined) THEN
      IF(fine_structure_defined .and. SPIN_FINE.eq.2)STOP "NOT READY"	  
      IF(int(l_real).eq.prn_l_trj .and.
     &  j12m12_print(1).eq.j12_s_st
     & .and. m12_s_st.eq.j12m12_print(2) ) THEN
      IF(tcur.eq.0) THEN
      OPEN(3452,FILE="TRAJECTORY.out")
      OPEN(3453,FILE="TRAJECTORY_ERRORS.out")	  
      WRITE(3452,'(a17,1x)',ADVANCE="NO") "TIME(a.u.)"
      WRITE(3452,'(a17,1x)',ADVANCE="NO") "R(a.u.)"
      WRITE(3452,'(a17,1x)',ADVANCE="NO") "P_R(a.u.)"
      WRITE(3452,'(a17,1x)',ADVANCE="NO") "Theta(rad)"
      WRITE(3452,'(a17,1x)',ADVANCE="NO") "P_theta(a.u.)"
      WRITE(3452,'(a17,1x)',ADVANCE="NO") "Phi(rad)"
      WRITE(3452,'(a17,1x)',ADVANCE="NO") "P_phi(a.u.)"	 	  
      DO i=1,number_of_channels	  
      WRITE(3452,'(a11,i3,4x)',ADVANCE="NO") "POP_CH#",i
      ENDDO
      WRITE(3452,*)
      WRITE(3453,'(a17,1x)',ADVANCE="NO") "TIME(a.u.)"
      WRITE(3453,'(a17,1x)',ADVANCE="NO") "X"	  
      WRITE(3453,'(a17,1x)',ADVANCE="NO") "Y"	  
      WRITE(3453,'(a17,1x)',ADVANCE="NO") "error_energy(%)"
      WRITE(3453,'(a17,1x)',ADVANCE="NO") "error_probab(%)"
!      WRITE(3453,'(a17,1x)',ADVANCE="NO") "mean_field(cm-1)"
      WRITE(3453,'(a17,1x)',ADVANCE="NO") "poten_term(cm-1)"	  
      WRITE(3453,*)	 
	  endif
      probab_traj = 0d0 
      DO i = 1,states_size
      chann_num = indx_chann(i)
      probab_traj(chann_num) =	probab_traj(chann_num) 
     &		+ (sys_var(i)**2+sys_var(i+states_size)**2)
      ENDDO	  
	  
	  bk_tmp_r=sys_var(1+states_size*2)
	  bk_tmp_phi=sys_var(5+states_size*2)
      WRITE(3453,'(e17.10,1x)',ADVANCE="NO") tcur
      WRITE(3453,'(e17.10,1x)',ADVANCE="NO") bk_tmp_r*cos(bk_tmp_phi)	  
      WRITE(3453,'(e17.10,1x)',ADVANCE="NO") bk_tmp_r*sin(bk_tmp_phi)	  
      WRITE(3453,'(e17.10,1x)',ADVANCE="NO") eror/E_sct*1d2
      WRITE(3453,'(e17.10,1x)',ADVANCE="NO") current_error*1d2
!      WRITE(3453,'(e17.10,1x)',ADVANCE="NO")
!     & Mjmr_cs(R,ind_to_buff,ind_to_buff)*autoeV*eVtown!pot	ind_to_buff	
      WRITE(3453,'(e17.10,1x)',ADVANCE="NO")pot*autoeV*eVtown!pot*Mjmr(R,k_st)autoeV*eVtown	
      WRITE(3452,'(e17.10,1x)',ADVANCE="NO") tcur
      DO i=states_size*2+1,states_size*2+6
      WRITE(3452,'(e17.10,1x)',ADVANCE="NO") sys_var(i)              	  
      ENDDO	  
      DO i=1,number_of_channels
      WRITE(3452,'(e17.10,1x)',ADVANCE="NO") probab_traj(i)              	  
      ENDDO
      WRITE(3452,*)
      WRITE(3453,*)		  
!      ENDIF	 
      ENDIF
      ENDIF	  
	  
! Finding 1st guess for time step of Bikram's propagator
	  bk_vel=dsqrt(current_vel_ini(1)**2+current_vel_ini(2)**2
     & +current_vel_ini(3)**2)
	  bk_int_r=2d0*sqrt(R**2-(sqrt(l_real*(l_real+1d0))/k_vec)**2)
	  int_1st_stp=bk_int_r/bk_vel/10d0
! Bikram End Oct 2019

      fail_odeint = .FALSE.
      IF(coupled_states_defined .and. 
!     & (.not. mpi_task_per_proc) .and. !!! 1775
     & (.not. term_pot_defined) .and. .FALSE.) THEN !! IMPROVED CS IS TURNED OFF
      DO i=1,states_size_cs
      sys_var_cs(i) = sys_var(ind_state_cs(i))
      IF(ind_state_cs(i).gt.states_size) STOP "CS FAILED"	  
      sys_var_cs(i+states_size_cs) = 
     & sys_var(ind_state_cs(i)+states_size)	 
      ENDDO	
      IF(myid.eq.0 .and. tcur.le.0d0) PRINT*, "INI_CS_SEEMS WORKING"  	  
      sys_var_cs(2*states_size_cs+1:2*states_size_cs+8) = 
     &  sys_var(2*states_size+1:2*states_size+8)
  	 
      IF(odeint_defined) THEN 
!      call odeint(sys_var_cs,states_size_cs*2+8,tcur,tcur+dt,
!     & eps_odeint,dt_step,min_t_stp,
!     & nok,nbad,derivs_cs_im,  
!     & rkqs,dt_corrected,fail_odeint)
      dt_step = dt_corrected
      ELSE
!      CALL derivs_cs_im(tcur,sys_var_cs,deriv_sys_cs)
!      CALL rk4(sys_var_cs,deriv_sys_cs, 
!     & states_size_cs*2+8,tcur,dt,sys_var_cs,
!     & derivs_cs_im)	  
      ENDIF	
      DO i=1,states_size_cs
      sys_var(ind_state_cs(i)) = sys_var_cs(i) 
      sys_var(ind_state_cs(i)+states_size) = 		  
     &  sys_var_cs(i+states_size_cs) 
	  
      ENDDO	
      sys_var(2*states_size+1:2*states_size+8) = 	  
     & sys_var_cs(2*states_size_cs+1:2*states_size_cs+8) 
	
	  
      GOTO 2222	  
      ENDIF	 	  
	  
! Bikram Start: April 2020, Adiabatic Trajectory Approximation
!	  traj_run = .FALSE.
!	  if(myid.eq.-1) then
	  if(bikram_adiabatic .and. .not. bk_b_chk) then
	  
	  if(rk4_defined) then
	  
	  if(bk_step_size) then
!	  if(counter.le.2) then
!	  dt = 2000.0d0
!	  else
	  call splint(bk_adia_t, bk_sys_var(:,7), bk_sys_var_der(:,7), 
     & bk_adia_n, tmstp, v_temp, yprm(7))
	  advnc_rk4 = .true.
	  tmstp = tcur
	  tmchk = tmstp
	  potbox = 0.d0
	  stps_cntr = 0
	  
	  do while(advnc_rk4)
	  tmstp = tmchk + tau*(stps_cntr + 0.50d0)
	  stps_cntr = stps_cntr + 1
	  call splint(bk_adia_t, bk_sys_var(:,7), bk_sys_var_der(:,7), 
     & bk_adia_n, tmstp, fb, yprm(7))
	  potbox = potbox + abs(fb)
	  i_max = ((tau*stps_cntr)**4d0)*potbox/stps_cntr
	  if(i_max.ge.magic) then
	  dt = tau*(stps_cntr - 1)
	  if(dt.lt.mintm) dt = mintm
	  advnc_rk4 = .false.
	  else if(tmstp.gt.bk_adia_t(bk_adia_n)) then
	  dt = bk_adia_t(bk_adia_n)-tmchk
	  advnc_rk4 = .false.
	  traj_run = .false.
	  else if(tcur.gt.0.d0) then
	  if(fa/fb.lt.0.d0 .and. fa/fb.eq.fa/fb) then
	  dt = tau*(stps_cntr - 1)
	  if(dt.lt.mintm) dt = mintm
	  fa = fb
	  advnc_rk4 = .false.
	  endif
	  endif
	  fa = fb
	  enddo	  
	  endif
	  
	  if(.not. bk_step_size) then
	  if((tcur+dt).ge.bk_adia_t(bk_adia_n)) then
	  dt = bk_adia_t(bk_adia_n)-tcur
	  traj_run = .false.
	  endif
	  endif
	  
	  if(int(j12(s_st)).eq.6 .and. int(m12(s_st)).eq.6 .and. 
     & int(l_real).eq.382) then
!	  print*, tcur, sys_var(33+states_size)
	  end if
	  call DERIVS_BK_adia(tcur, sys_var(1:states_size*2), 
     & deriv_sys(1:states_size*2))
	  if(int(j12(s_st)).eq.6 .and. int(m12(s_st)).eq.6 .and. 
     & int(l_real).eq.382) then
!	  print*, tcur, sys_var(33+states_size), deriv_sys(33+states_size)
	  end if
      call rk4(sys_var(1:states_size*2), deriv_sys(1:states_size*2), 
     & states_size*2,tcur,dt,sys_var(1:states_size*2),DERIVS_BK_adia)
	  bikram_t = tcur + dt
	  if(int(j12(s_st)).eq.6 .and. int(m12(s_st)).eq.6 .and. 
     & int(l_real).eq.382) then
	  print*, tcur, sys_var(33+states_size), deriv_sys(33+states_size)
	  end if
	  endif
	  
	  IF(odeint_defined) THEN
	  if(bikram_int) then
	  call bk_int_adia(sys_var(1:states_size*2), states_size*2, 
     & tcur, bk_adia_t(bk_adia_n), eps_odeint, int_1st_stp, min_t_stp,
     & nok, nbad, DERIVS_BK_adia,  
     & rkqs, dt_corrected, fail_odeint, f_time, R_fin,
     & vibration_cnt, numb_oscl_prds,
     & period_cnt, numb_orb_prds,
     & bk_adia_t, bk_adia_n, 
     & bk_sys_var(:,1), bk_sys_var_der(:,1),
     & bk_sys_var(:,5), bk_sys_var_der(:,5))
	  else
      call odeint(sys_var(1:states_size*2), states_size*2, 
     & tcur, bk_adia_t(bk_adia_n), eps_odeint, dt_step, min_t_stp,
     & nok, nbad, DERIVS_BK_adia,  
     & rkqs, dt_corrected, fail_odeint)
	  endif
      dt_step = dt_corrected 	
	  traj_run = .false.
	  bikram_t = bk_adia_t(bk_adia_n)
      endif

	  do i = 1, nnn-1
	  call splint(bk_adia_t, bk_sys_var(:,i), bk_sys_var_der(:,i), 
     & bk_adia_n, bikram_t, sys_var(i+states_size*2), yprm(i))
	  enddo

!Switching to normal propagation
	  else
      IF(rk4_defined) THEN
      CALL derivs_BK(tcur,sys_var,deriv_sys)
      CALL rk4(sys_var,deriv_sys,states_size*2+8,tcur,dt,sys_var,
     & derivs_BK)
      ENDIF
!!! CALLING ADAPTIVE STEP_SIZE CONTROL METHOD	  
      IF(odeint_defined) THEN
! Bikram Start Nov 2019: 
	  if(bikram_int) then
	  call bk_int(sys_var,states_size*2+8,tcur,time_lim,
     & eps_odeint,int_1st_stp,min_t_stp,
     & nok,nbad,derivs_BK,  
     & rkqs,dt_corrected,fail_odeint,f_time,R_fin,nmbr_r,
     & vibration_cnt,numb_oscl_prds,
     & nmbr_phi,period_cnt,numb_orb_prds)
	  else
      call odeint(sys_var,states_size*2+8,tcur,tcur+dt,
     & eps_odeint,dt_step,min_t_stp,
     & nok,nbad,derivs_BK,  
     & rkqs,dt_corrected,fail_odeint)
	  endif
! Bikram End.

!!!!!!!!! IN THE END REPLACE GUESS STEP WITH NEW CORRECTED VALUE
      dt_step = dt_corrected 	
      endif
      endif
!      endif 
!Bikram End.

! Bikram Oct '19 start:	  
	  updated_R=sys_var(1+states_size*2)
	  updated_phi=sys_var(5+states_size*2) 
!	  if(bikram_int) traj_run = .FALSE.
! Bikram End	  

! PRINTING TRAJECTORY. FOR NOW IS DISABLED
!!!! SWITCH BACK TO SPHER COORDIANTES
2222  IF(cartes_defined .and. l_real.eq.0d0) THEN
      q_cartes =  sys_var(1+states_size*2:6+states_size*2)
!!!!! MAKE MAGIC HERE  !!!! DO 
      R=abs(q_cartes(1))
      IF(q_cartes(1).ge.0d0)  THEN 
	  phi = 0d0
      ELSE
      phi = acos(-1d0)	  
      ENDIF	  
      ql = acos(-1d0)/2d0
      pr = massAB_M*q_cartes(4)*cos(phi)
      qphi =0d0
      l = 0d0
      sys_var(1+states_size*2) = R
      sys_var(2+states_size*2) = pr
      sys_var(3+states_size*2) = ql
      sys_var(4+states_size*2) = l
      sys_var(5+states_size*2) = phi
      sys_var(6+states_size*2) = qphi	  
      ENDIF
!!!!!  IF ODEINT FAILED	  
      IF(INTERUPT_PROPAGATION .and. fail_odeint) THEN
      time_st = tcur
! "err_prb_tmp" - PROBABILITY CONSERVATION	  
      err_prb_tmp = abs(current_error)*1d2
!      PRINT*,	"err_prb_tmp",err_prb_tmp,myid  
! "err_prb_tmp" - ENERGY CONSERVATION	  
      err_ener_tmp = eror/E_sct*1d2
!     CHECKING THAT PROPAGTION IS INTERUPTED
      PRINT*, myid, "fail_odeint,PROPAGATION ENDED"	  
      RETURN
      ENDIF

!!!!! END	  

! MYID - THE MPI TASK RANK 
! PRINTING TRAJECTORY. FOR NOW IS DISABLED
! MYID - THE MPI TASK RANK
      IF(fail_odeint) THEN
      IF(int(myid/mpi_task_per_traject)*mpi_task_per_traject.eq.myid)
     & THEN
      WRITE(*,'(a26,1x,i5,1x,a11,1x,i4,1x,a7,i5,1x,a10,1x,e17.10,1x,a11,
     & 1x,f9.3)')
     & "ODEINT FAILED FOR MYID =  ", 
     & myid, "#trajectory=",
     & itraject, "#state=",
     & s_st, "at time = ",
     & tcur, "for energy=",
     & U(i_ener)       	  
      ENDIF
      EXIT	  
      ENDIF
	  
! Bikram Start:

! Check if the molecule is trapped and vibrating
	  if(bikram_int)  then
	  else
      if(tcur.gt.(time_st+dt) .and. counter.gt.2) then
	  if((tmp_R2-tmp_R1).lt.0.00d0 .and. (tmp_R1-updated_R).gt.0.00d0) 
     & vibration_cnt=vibration_cnt + 1
	  endif
	  tmp_R2=tmp_R1
      tmp_R1=updated_R
	  endif
! Bikram End.
	  
! RECOMPUTING CURRENT TIME	  
	  if(bikram_int) then
	  else
      tcur = tcur + dt
	  endif
! COMPTUNING NUMBER OF PERIODS FOR AN ORBITING TRAJECTORY	  
	  if(bikram_int)  then
	  else
      period_cnt = int(abs(updated_phi)/2d0/pi)
	  endif
! CHECK IF WE SHOULD EXIT THE TRAJECTORY	  
      IF(numb_orb_prds.le.period_cnt) THEN
      traj_run = .FALSE.
	  orbiting_defined = .TRUE.
      ENDIF

! Bikram Start:
      IF(vibration_cnt.ge.numb_oscl_prds) THEN
      traj_run = .FALSE.
	  orbiting_defined = .TRUE.	  
      ENDIF
! Bikram End.
	  
      IF(tm_lim_defined.and.tcur.ge.time_lim) traj_run=.FALSE.
      IF(tm_lim_defined.and.tcur.ge.time_lim) then
	  print*, "Warning: Time limit has reached. You may want to
     & increase time limit", tcur, int(l_real)
	  end if
      IF(updated_R.gt.R_fin) traj_run = .FALSE.

!!!! SAVING DATA IN CHEKPOINT
      ! IF(write_check_file_defined .and. .not.monte_carlo_defined) THEN
      ! TRAJECT_DATA_CHECK_FILE(1:sys_var_size) = sys_var
      ! TRAJECT_DATA_CHECK_FILE(sys_var_size+1) = tcur
      ! TRAJECT_DATA_CHECK_FILE(sys_var_size+2) = eror
      ! TRAJECT_DATA_CHECK_FILE(sys_var_size+3) = current_error
      ! TRAJECT_DATA_CHECK_FILE(sys_var_size+4) = deflect_angle
      ! TRAJECT_DATA_CHECK_FILE(sys_var_size+5) = dt_corrected
      ! TRAJECT_DATA_CHECK_FILE(sys_var_size+6) = dt_step	  
      ! ENDIF	  
!!!!! CHECK CURRENT TIME FOR CHECKPOINT	  
      IF(mpi_task_defined) THEN
      IF(int(myid/mpi_task_per_traject)*mpi_task_per_traject.eq.myid)
     & THEN
      CALL TIME_CHECKER(INTERUPT_PROPAGATION)
      ENDIF
      DO i=1,traject_roots	  
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	   
      origin = 0		  
      CALL MPI_Bcast(INTERUPT_PROPAGATION, 1,
     & MPI_LOGICAL, origin, comms(i),ierr_mpi)
      ENDIF
      ENDDO
      ELSE
      CALL TIME_CHECKER(INTERUPT_PROPAGATION)
      ENDIF	  
      IF(INTERUPT_PROPAGATION) THEN
      time_st = tcur
! "err_prb_tmp" - PROBABILITY CONSERVATION	  
      err_prb_tmp = abs(current_error)*1d2
!      PRINT*,	"err_prb_tmp",err_prb_tmp,myid  
! "err_prb_tmp" - ENERGY CONSERVATION	  
      err_ener_tmp = eror/E_sct*1d2
!     CHECKING THAT PROPAGTION IS INTERUPTED
      PRINT*, myid, "PROPAGATION ENDED"	  
      RETURN
      ENDIF
	  
! Bikram Start June 2020:
	  if(bikram_print .and. mod(counter,10).eq.0) then
	  bknorm = 0.d0
	  do i = 1, states_size
	  bknorm = bknorm +(sys_var(i)**2+sys_var(i+states_size)**2)
	  enddo
	  bk_probab = 0.d0
	  do i = 1, states_size
      chann_num = indx_chann(i)
      if(chann_num.eq.chann_ini) cycle
	  bk_probab = bk_probab +(sys_var(i)**2+sys_var(i+states_size)**2)
	  enddo	  
	  write(100002,'(6(e23.16,1x))')
     & sys_var(1+states_size*2)*dcos(sys_var(5+states_size*2)),
     & sys_var(1+states_size*2)*dsin(sys_var(5+states_size*2)),
     & sys_var(1+states_size*2),sys_var(5+states_size*2),
     & sys_var(2+states_size*2),sys_var(6+states_size*2)	  
!	  write(100002,'(f12.5,1x,e23.16)') tcur,bknorm
!	  do i = 1, states_size
!	  if(i.eq.1) then
!	  write(100002,'(f12.5,1x,e23.16,1x,
!!     & e23.16,1x,e23.16,1x,
!     & e23.16)', advance = 'no')
!     & tcur,! bknorm - 1.0d0, bk_probab, dt,! pot*autoeV*eVtown, 
!     & sys_var(i), sys_var(i+states_size)!, abs(1d0 - P_total)*1d2, sys_var(1+states_size*2), tmp_R1, tmp_R2
!!    & sys_var(1+states_size*2),sys_var(5+states_size*2),
!	  else
!	  write(100002,'(1x,e23.16,1x,e23.16)', advance = 'no') 
!     & sys_var(i), sys_var(i+states_size)
!	  endif
!	  enddo
!	  write(100002,*)
!	  endif
!	  do i = 1, states_size
!	  if(i.eq.1) then
!	  write(100002,'(f12.5,1x,e23.16)', advance = 'no')
!     & tcur, sys_var(i)**2 + sys_var(i+states_size)**2
!	  else
!	  write(100002,'(1x,e23.16)', advance = 'no') 
!     & sys_var(i)**2 + sys_var(i+states_size)**2
!	  endif
!	  enddo
!	  write(100002,*)
	  
	  endif
! Bikram End. 
	  
      ENDDO
!   END OF THE TRAJECTORY LOOP	  
!	  print *, 'prop_done', myid, l_real 
!	  call flush(output_unit)
!	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
!      IF(INTERUPT_PROPAGATION) RETURN
!	  open(12345,file='#STEPS.out', access = 'append')
!	  write(12345,*)int(l_real), counter
!	  close(12345)
	  
! Bikram Start: April 2020, Adiabatic Trajectory Approximation
	  if(bikram_adiabatic .and. .not. bk_b_chk) then
	  deallocate(bk_adia_t, bk_sys_var, bk_sys_var_der)	
	  deallocate(yp1, ypn, yprm)
	  endif	
	  if(bikram_print) close(100002)
! Bikram End.
	 
! Bikram start Oct 2019:
	  
	  if(bikram_int) tcur=f_time
	  if(.not. check_point_defined) then
	  R =  sys_var(1+states_size*2)
      pr =  sys_var(2+states_size*2)
      ql = sys_var(3+states_size*2)
      l = sys_var(4+states_size*2)
      phi = sys_var(5+states_size*2)
      qphi = sys_var(6+states_size*2)

! "Eqf" - THE QUANTUM ENERGY DURING THE COLLISION	  
      Eqf = 0d0	  
      DO i = st_sum_rn_min,st_sum_rn_max
      IF(coupled_states_defined .and. (m12(s_st) .ne. m12(i))) CYCLE	  
      Eqf = Eqf + (sys_var(i)**2+sys_var(i+states_size)**2)*Ej(i)
      ENDDO
      IF(mpi_task_defined) THEN
      DO i=1,traject_roots	  
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	  
      origin = 0		  
      CALL MPI_Reduce(Eqf,buffer_eq, 1, MPI_REAL8,MPI_SUM,origin,
     & comms(i),ierr_mpi )
      Eqf = buffer_eq	 
      CALL MPI_Bcast(Eqf,1,
     & MPI_REAL8, origin, comms(i),ierr_mpi)	 
      ENDIF		  
      ENDDO		  
      ENDIF		
!	  print *, 'chk3', myid, l_real  
!	  call flush(output_unit)
!	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
!	"Ef" - KINETIC ENERGY DURING THE COLLISION	  
      Ef = pr**2/2d0/massAB_M+l**2/2d0/massAB_M/R**2+	
     & qphi**2/2d0/massAB_M/R**2/sin(ql)**2	
! RECOMPUTING "pot" - THE POTENTIAL ENERGY
      pot = 0d0	
      IF(coupled_states_defined .and. 
!     & (.not. mpi_task_per_proc) .and. !!! 1670
     & (.not. term_pot_defined)) THEN
      DO i=1,states_size_cs
	  
	  DO j=1,i
      sys_var_cs(i) = sys_var(ind_state_cs(i))
      IF(ind_state_cs(i).gt.states_size) STOP "CS FAILED"	  
      sys_var_cs(i+states_size_cs) = 
     & sys_var(ind_state_cs(i)+states_size)	
      sys_var_cs(j) = sys_var(ind_state_cs(j))
      IF(ind_state_cs(j).gt.states_size) STOP "CS FAILED"	  
      sys_var_cs(j+states_size_cs) = 
     & sys_var(ind_state_cs(j)+states_size)	
	 
      pot = pot + 2d0*Mjmr_cs(R,i,j)/(1+delta(i,j))*
     ^   ((sys_var_cs(i)*(sys_var_cs(j)) + 
     ^ sys_var_cs(i+states_size_cs)*(sys_var_cs(j+states_size_cs)))*
     & dcos(tcur*(Ej_cs_im(j)-Ej_cs_im(i)))
     & - (sys_var_cs(i+states_size_cs)*(sys_var_cs(j)) - 
     ^ sys_var_cs(i)*(sys_var_cs(j+states_size_cs)))*
     & dsin(tcur*(Ej_cs_im(j)-Ej_cs_im(i))))	  

      ENDDO
      ENDDO	  
	  ELSE

! Bikram Start May 2020:	  
	  if(.not. bk_nrg_err .and. .not. bikram_print) then	  
	  else
	  call potential_data_bk_adia(R,bk_matt)
	  do k= 1,ph_cntr_bkk
	  bk_sin_coss(1,k) = dsin(tcur*bk_delta_E(k))
	  bk_sin_coss(2,k) = dcos(tcur*bk_delta_E(k))
	  enddo
      DO k = summ_range_min,summ_range_max! POT ENERGY COMPUTING
      i = ind_mat_bk(1,k)!ind_mat(1,k)
      IF(coupled_states_defined .and. (m12(s_st) .ne. m12(i))) CYCLE	  
      j = ind_mat_bk(2,k)!ind_mat(2,k)
	  bk_del = 2.d0
	  if(i.eq.j) bk_del = 1.d0
	  tmp_indx = bk_indx(k)
      iphase = bk_sin_coss(1,tmp_indx)
      rphase = bk_sin_coss(2,tmp_indx)
       pot = pot + bk_del*bk_matt(k)*
     & ((sys_var(i)*sys_var(j) + 
     &   sys_var(i+states_size)*sys_var(j+states_size))*rphase -
     &  (sys_var(i+states_size)*sys_var(j) - 
     &   sys_var(i)*sys_var(j+states_size))*iphase)
      ENDDO
	  deallocate(bk_matt, bk_sin_coss)
	  endif
	  call resize_adia_dealloc
! Bikram End.
	  
      IF(mpi_task_defined) THEN
      DO i=1,traject_roots	  
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	 
      origin = 0		  
      CALL MPI_Reduce(pot,buffer_pot, 1, MPI_REAL8,MPI_SUM,origin,
     & comms(i),ierr_mpi )
      pot = buffer_pot
      CALL MPI_Bcast(pot,1,
     & MPI_REAL8, origin, comms(i),ierr_mpi)	 
      ENDIF	  
      ENDDO
      ENDIF
      ENDIF	    
! Bikram Start May 2020:	
	  if(bikram_save_traj .and. .not. bk_b_chk) then 
	  write(100003) tcur+dt,
     & sys_var(1+states_size*2), sys_var(2+states_size*2), 
     & sys_var(3+states_size*2), sys_var(4+states_size*2), 
     & sys_var(5+states_size*2), sys_var(6+states_size*2), pot
	  bk_adia_n = bk_adia_n + 1
	  write(100004) bk_adia_n
	  close(100003)	
	  close(100004)	
	  endif
! Bikram End.
! COMPUTING THE ERROR IN ENRGY CONSERVATION "eror" 

      IF( abs(Eqf+Ef+pot-Eki-Eqi)>eror) eror = 
     & abs(Eqf+Ef+pot-Eki-Eqi)
! PROBABILITY CONSERVATION
      P_total = 0d0	  
      DO i = 1,states_size
! "P_total" MUST BE 1. HOWEVER, THERE IS A PROPAGATION ERROR	 
      P_total=P_total+(sys_var(i)**2+sys_var(i+states_size)**2)
      ENDDO
      IF(current_error.lt. abs(1d0 - P_total)) current_error
     &	  = abs(1d0 - P_total)	
! COMPUTING FINAL VELOCITY	  
      DO i=1,3
      velocity_fin(i) = ort_r(i)*pr/massAB_M+(
     & ort_teta(i)*l/R/massAB_M
     & + ort_phi(i)*qphi/R/sin(ql)**2/massAB_M)
      ENDDO
      probab = 0d0
      P_total = 0d0
! COMPUTING OPACITIES ( ARRAY "probab" WHICH SIZE IS THE NUMBER OF CHANNELS)
      DO i = 1,states_size
      chann_num = indx_chann(i)
      IF(chann_num.gt.number_of_channels .and. myid.eq.0)
     & PRINT*, "INTIALIZATION WENT WRONG"
      IF(chann_num.lt.1 .and. myid.eq.0)
     & PRINT*, "INTIALIZATION WENT WRONG"	  	 
      IF(monte_carlo_defined) THEN	  
      probab(chann_num) =	probab(chann_num) 
     &		+ (sys_var(i)**2+sys_var(i+states_size)**2)*
     & (2d0*J_tot + 1d0)/k_vec**2
      ELSE
      probab(chann_num) =	probab(chann_num) 
     &		+ (sys_var(i)**2+sys_var(i+states_size)**2)
      ENDIF	  
! "P_total" MUST BE 1. HOWEVER, THERE IS A PROPAGATION ERROR	 
      P_total=P_total+(sys_var(i)**2+sys_var(i+states_size)**2)
      ENDDO	
	  
! Bikram Start:
!      print *, 'chk6', myid, l_real  
!	  call flush(output_unit)
!	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF (no_orbits_defined.and. R.lt.R_max_dist) THEN
      probab = 0.00d0
	  probab(chann_ini) = 1.00d0
	  deflect_angle = 0.00d0
! Bikram End.
	  
      ENDIF	  	  
! "err_prb_tmp" - PROBABILITY CONSERVATION	  
      err_prb_tmp = current_error*1d2
!      PRINT*,	"err_prb_tmp",err_prb_tmp,myid  
! "err_prb_tmp" - ENERGY CONSERVATION	  
      err_ener_tmp = eror/E_sct*1d2
	  
! Bikram start Oct 2019:
      IF(prn_l_defined) THEN
      IF(fine_structure_defined .and. SPIN_FINE.eq.2)STOP "NOT READY"	  
      IF(int(l_real).eq.prn_l_trj .and.
     &  j12m12_print(1).eq.j12_s_st
     & .and. m12_s_st.eq.j12m12_print(2) ) THEN

      probab_traj = 0d0
      DO i = 1,states_size
      chann_num = indx_chann(i)
      probab_traj(chann_num) =	probab_traj(chann_num) 
     &		+ (sys_var(i)**2+sys_var(i+states_size)**2)
      ENDDO	  
	   
	  bk_tmp_r=sys_var(1+states_size*2)
	  bk_tmp_phi=sys_var(5+states_size*2)
      WRITE(3453,'(e17.10,1x)',ADVANCE="NO") tcur
      WRITE(3453,'(e17.10,1x)',ADVANCE="NO") bk_tmp_r*cos(bk_tmp_phi)	  
      WRITE(3453,'(e17.10,1x)',ADVANCE="NO") bk_tmp_r*sin(bk_tmp_phi)	  
      WRITE(3453,'(e17.10,1x)',ADVANCE="NO") eror/E_sct*1d2
      WRITE(3453,'(e17.10,1x)',ADVANCE="NO") current_error*1d2
!      WRITE(3453,'(e17.10,1x)',ADVANCE="NO")
!     & Mjmr_cs(R,ind_to_buff,ind_to_buff)*autoeV*eVtown!pot	ind_to_buff	
      WRITE(3453,'(e17.10,1x)',ADVANCE="NO")pot*autoeV*eVtown!pot*Mjmr(R,k_st)autoeV*eVtown	
      WRITE(3452,'(e17.10,1x)',ADVANCE="NO") tcur
      DO i=states_size*2+1,states_size*2+6
      WRITE(3452,'(e17.10,1x)',ADVANCE="NO") sys_var(i)              	  
      ENDDO	  
      DO i=1,number_of_channels
      WRITE(3452,'(e17.10,1x)',ADVANCE="NO") probab_traj(i)              	  
      ENDDO
      WRITE(3452,*)
      WRITE(3453,*)		  
      CLOSE(3452)
      CLOSE(3453)	  
      ENDIF
      ENDIF	  
! Bikram End
      IF(monte_carlo_defined) RETURN	  
! "ampl_wf_real" AND "ampl_wf_imag" - REAL and IMAGINARY PART OF AMPLITUDE	 
      ampl_wf_real(itraject,myid+1) = sys_var(s_st)!!!!!!! HERE WE SHOULD MOVE
      ampl_wf_imag(itraject,myid+1) = sys_var(s_st+states_size)
! "phase_wf" - SCATTERING PHASE	  
      phase_wf(itraject,myid+1) = 
     & datan2(sys_var(s_st+states_size),sys_var(s_st))
! "angle_scatter" - AN DEFLECTION ANGLE	  
      angle_scatter(itraject,myid+1) = deflect_angle
      angle_phase(itraject,myid+1) = vel_cross(3)/abs(vel_cross(3))	

! Bikram Start:
	  bk_tym(itraject,myid+1) = tcur
	  bk_prd(itraject,myid+1) = period_cnt
	  bk_vib(itraject,myid+1) = vibration_cnt
	  bk_erre(itraject,myid+1) = err_ener_tmp
	  bk_errp(itraject,myid+1) = err_prb_tmp
! Bikram End.
	  
      IF( prn_l_defined.and. prn_l_trj.eq.int(l_real) .and. myid.eq.-1)
     & THEN
      WRITE(probab_out(1:2),'(i2.2)') myid	  
      OPEN(23,FILE=probab_out,POSiTION="APPEND")
      WRITE(23,'(a17,1x,e14.6)') "IMPACT_PARAMETER=", l_real/k_vec
      WRITE(23,'(a17,1x,i4)') "ORBITAL MOMENTUM=", int(l_real)
      WRITE(23,'(a13,1x,i4)') "PROJECTION M=",m12_s_st	  
      IF(identical_particles_defined) THEN
      WRITE(23,'(a8,1x,a6,1x,a17,1x,a6)')
     &  "#CHANNEL","#STATE","TRANS_PROBABILITY","PARITY"	  
      DO i = 1,states_size
      chann_num = indx_chann(i)	  
      WRITE(23,'(i8,1x,i6,1x,e17.9,1x,i6)')chann_num, i,
     & sys_var(i)**2+sys_var(i+states_size)**2,parity_state(i)
      ENDDO	  
      ELSE	  
      WRITE(23,'(a8,1x,a6,1x,a17)')
     &  "#CHANNEL","#STATE","TRANS_PROBABILITY"	  
      DO i = 1,states_size
      chann_num = indx_chann(i)	  
      WRITE(23,'(i8,1x,i6,1x,e17.9)')chann_num, i,
     & sys_var(i)**2+sys_var(i+states_size)**2!,parity_state(i)
      ENDDO
      ENDIF	  
      CLOSE(23)	  
      ENDIF
	  
!!!   SAVING TRAJECTORY ENDS
      ds_state = s_st-ini_st_check_file+1
      TRAJECT_END_DATA(1,itraject,ds_state,i_ener) = l_orb/k_vec
      TRAJECT_END_DATA(2,itraject,ds_state,i_ener) = l_real
      TRAJECT_END_DATA(3,itraject,ds_state,i_ener) = R
      TRAJECT_END_DATA(4,itraject,ds_state,i_ener) = phi
      TRAJECT_END_DATA(5,itraject,ds_state,i_ener) =
     & sys_var(s_st)**2+sys_var(s_st+states_size)**2
      TRAJECT_END_DATA(6,itraject,ds_state,i_ener) = err_ener_tmp
      TRAJECT_END_DATA(7,itraject,ds_state,i_ener) = err_prb_tmp	

      IF(orbit_traj_defined) THEN
!      OPEN(2345,FILE="ORBITING_TRAJECTORY.out",POSiTION="APPEND")

! Bikram Start:
      IF(R.lt.R_max_dist .and. orbiting_defined) THEN
! Bikram End.
	  
      IF(mpi_task_defined) THEN	  
      DO i=1,traject_roots	  
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN	  
      CALL MPI_Comm_rank(comms(i),origin ,ierr_mpi)
      IF(origin.eq.0) THEN
!      WRITE(2345,'(i7,1x,i8,1x,e18.11,1x,e18.11)')myid,itraject,
!     & J_tot/k_vec,l_real
      ORBIT_TRAJ_DATA(1,itraject) = l_orb/k_vec
      ORBIT_TRAJ_DATA(2,itraject) = l_real
      ORBIT_TRAJ_DATA(3,itraject) = R
      ORBIT_TRAJ_DATA(4,itraject) = phi
      ORBIT_TRAJ_DATA(5,itraject) =
     & sys_var(s_st)**2+sys_var(s_st+states_size)**2
      ORBIT_TRAJ_DATA(6,itraject) = err_ener_tmp
      ORBIT_TRAJ_DATA(7,itraject) = err_prb_tmp  	  
      ORBIT_TRAJ_DATA(8,itraject) = period_cnt		  	  
      ORBIT_TRAJ_DATA(9,itraject) = vibration_cnt	  	  
      ORBIT_TRAJ_DATA(10,itraject) = 1.00d0	  	  
      ORBIT_TRAJ_FLAG(itraject)  = .TRUE.	  
      ENDIF
      ENDIF
      ENDDO
!      PRINT*, "orbiting",myid, itraject,
!     & J_tot,l_real,ORBIT_TRAJ_FLAG(itraject)	!!!!! TESTING  
      ELSE
      ORBIT_TRAJ_DATA(1,itraject) = l_orb/k_vec
      ORBIT_TRAJ_DATA(2,itraject) = l_real
      ORBIT_TRAJ_DATA(3,itraject) = R
      ORBIT_TRAJ_DATA(4,itraject) = phi
      ORBIT_TRAJ_DATA(5,itraject) =
     & sys_var(s_st)**2+sys_var(s_st+states_size)**2
      ORBIT_TRAJ_DATA(6,itraject) = err_ener_tmp
      ORBIT_TRAJ_DATA(7,itraject) = err_prb_tmp  	  
      ORBIT_TRAJ_DATA(8,itraject) = period_cnt		  	  
      ORBIT_TRAJ_DATA(9,itraject) = vibration_cnt		  	  
      ORBIT_TRAJ_DATA(10,itraject) = 1.00d0		  	  
      ORBIT_TRAJ_FLAG(itraject)  = .TRUE.		  
      ENDIF	  
      ENDIF
!      CLOSE(2345)
!!!!! TESTING MERELY	

!!!! TESTING MERELY  
      ENDIF

!--------------------------------------------------------------------
! printing trajectory end information to checkpoint files
!--------------------------------------------------------------------
	  if(write_check_file_defined) then
	  write(chk_file_unit) tcur
	  do i = 1, 6
	  write(chk_file_unit) sys_var(i+states_size*2)
	  end do
	  
	  do i = 1, number_of_channels
	  write(chk_file_unit) probab(i)
	  end do
	  
      if(monte_carlo_defined) return 
      
	  write(chk_file_unit) ampl_wf_real(itraject,myid+1), 
     & ampl_wf_imag(itraject,myid+1), phase_wf(itraject,myid+1),
     & angle_scatter(itraject,myid+1), angle_phase(itraject,myid+1),
     & bk_tym(itraject,myid+1), bk_prd(itraject,myid+1),
     & bk_vib(itraject,myid+1), bk_erre(itraject,myid+1),
     & bk_errp(itraject,myid+1),
     & sys_var(s_st), sys_var(s_st+states_size)
	  end if
	  end if

!--------------------------------------------------------------------
! reading trajectory end information from checkpoint 
! files as RESTART = YES indicated
!--------------------------------------------------------------------
	  if(check_point_defined) then
	  read(chk_file_unit) tcur
	  do i = 1, 6
	  read(chk_file_unit) sys_var(i+states_size*2)
	  end do  
	  R =  sys_var(1+states_size*2)
      pr =  sys_var(2+states_size*2)
      ql = sys_var(3+states_size*2)
      l = sys_var(4+states_size*2)
      phi = sys_var(5+states_size*2)
      qphi = sys_var(6+states_size*2)
	  
	  do i = 1, number_of_channels
	  read(chk_file_unit) probab(i)
	  end do
	  
      if(monte_carlo_defined) return 
      
	  read(chk_file_unit) ampl_wf_real(itraject,myid+1), 
     & ampl_wf_imag(itraject,myid+1), phase_wf(itraject,myid+1),
     & angle_scatter(itraject,myid+1), angle_phase(itraject,myid+1),
     & bk_tym(itraject,myid+1), bk_prd(itraject,myid+1),
     & bk_vib(itraject,myid+1), bk_erre(itraject,myid+1),
     & bk_errp(itraject,myid+1),
     & sys_var(s_st), sys_var(s_st+states_size)
	  end if

      RETURN
      	  
      END SUBROUTINE DOTRAJECT
	  
      REAL*8 FUNCTION delta(m,n)
! This function is updated by Bikramaditya Mandal
      INTEGER m,n
      delta=0d0
      IF(m.eq.n) delta = 1d0
      END

      SUBROUTINE OUT_OF_RANGE(V,DV,R,k)
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE GRID_INTEGR
      USE CONSTANTS
      USE MPI_TASK_TRAJECT
      USE MPI_DATA
      USE OLD_MIJ	  
      IMPLICIT NONE
      INTEGER k,k_real	  
      REAL*8 x,y,y_prime,V,DV,R
      REAL*8, PARAMETER :: dx_inf = 1d-4	  
      V = 0d0
      dV = 0d0

      x  = R_COM(1)+dx_inf
      IF(mpi_task_defined) THEN
      k_real = k -portion_of_MIJ_per_task(1,myid+1)+1	  
      ELSE
      k_real = k	  
      ENDIF	  
      CALL splint(R_COM,Mat_el(:,k_real),Mat_el_der(:,k_real)
     & ,n_r_coll,x,y,y_prime)
      DV = y_prime      
      V = DV*(R-x) + y
      END SUBROUTINE OUT_OF_RANGE
	  !! PEREDELAT'
      FUNCTION Mjmr(R,k)
! This function is updated by Bikramaditya Mandal
      USE VARIABLES
      USE GRID_INTEGR
      USE INTEGRATOR	  
      USE CONSTANTS
      USE MPI_TASK_TRAJECT
      USE MPI_DATA
      USE OLD_Mij		  
      IMPLICIT NONE	   
      REAL*8 x,Mjmr,y,y_prime,R,V_COULPING_TERMS
      INTEGER k,k_real,i_exp_term
      x = R
      Mjmr = 0d0
      IF(term_pot_defined) THEN
      CALL COMPUTE_V_TERMS(R)
      DO i_exp_term=1,nterms
      CALL TERM_PES_MATRIX_ELEMENT(V_COULPING_TERMS,k,i_exp_term)	  
      Mjmr = Mjmr	+
     &	V_COULPING_TERMS*V_TERMS(i_exp_term)
      ENDDO		  
      RETURN	  
      ENDIF	  
      IF(x.lt.R_COM(1)) THEN
      CALL 	OUT_OF_RANGE(Mjmr,y_prime,R,k)  
      RETURN
      ENDIF	  
      IF(x.ge.R_COM(n_r_coll)) RETURN
      IF(mpi_task_defined) THEN
      k_real = k - portion_of_MIJ_per_task(1,myid+1)+1	  
      ELSE
      k_real = k	  
      ENDIF
	  
      CALL splint(R_COM,Mat_el(:,k_real),Mat_el_der(:,k_real)
     & ,n_r_coll,x,y,y_prime)

!--------------------------------------------------------------------
! Shifting matrix elements using the last value of R
!--------------------------------------------------------------------
      Mjmr = y - Mat_el(n_r_coll,k_real)
      END FUNCTION Mjmr
	  
      SUBROUTINE INELAST_CALC(sigma_el,nchann)
! This subroutine is updated by Bikramaditya Mandal
      USE INTEGRATOR
      USE MPI_DATA
      USE MPI_TASK_TRAJECT	  
	  use bk_l_values												!Bikram Feb 2021
      IMPLICIT NONE
      LOGICAL EVEN_NUM	  
      LOGICAL, ALLOCATABLE :: nproc_actual(:)	  
      INTEGER nchann,tot_num_of_traj_actual  
      INTEGER i,i_ip,i_traj,k,j_chann,traj_works,j_int_traj,
     & l_int_traj,id_traject 
      REAL*8 sigma_el(nchann)
      REAL*8, PARAMETER :: a_bohr = 5.2917721067d-1
      REAL*8, PARAMETER :: conv_E_cm_Ev = 27.21113845d0*8065.5448d0	  
      REAL*8, ALLOCATABLE :: j_def(:),opac_chann_all(:,:)
      REAL*8 posit_imp,negat_imp,part_cross	  
	  
	  real*8, allocatable :: bk_l_tmp(:)
	  integer tmp_l
	  real*8, allocatable :: bk_prob1(:,:), spln_der(:,:)
	  real*8 spln_yp1, spln_yp2, spln_yp, prob_reslt
	  
      tot_num_of_traj_actual = J_UP_INT-J_DOWN_INT+1	  
      ALLOCATE(j_def(tot_num_of_traj_actual),
     & opac_chann_all(nchann,tot_num_of_traj_actual))
	  	  
      opac_chann_all = 0d0	 
      sigma_el = 0d0
      DO k=1,nchann	  
      DO j_int_traj = J_DOWN_INT,J_UP_INT
      traj_works = j_int_traj - J_DOWN_INT + 1	  
      DO l_int_traj = abs(j_int_traj-j_int_ini),j_int_traj+j_int_ini
      IF(myid.eq.0) THEN
      ENDIF
	  
!      GOTO 3456
      CALL TRAJ_ORB(l_int_traj,
     & id_traject,
     & 1,
     & dl_step_integ,
     & ident_skip,
     & L_MIN_TRAJECT)
	 
!3456      id_traject = (l_int_traj-L_MIN_TRAJECT)/dl_step_integ+1
      i_traj = itraj_myid_l(1,id_traject)	  
      i_ip = itraj_myid_l(2,id_traject)
!      GOTO	5782  

      IF(ident_skip) THEN
      IF(parity_st.eq.1) THEN 
      posit_imp = exchange_ident
      negat_imp = 1d0 - exchange_ident	  
      ELSE
      posit_imp = 1d0 - exchange_ident	
      negat_imp = exchange_ident	  
      ENDIF	  
      ENDIF	  
	  
      id_proc = i_ip*mpi_traject
      part_cross =
     & probab_J_all(k,i_traj,id_proc+1)*(2d0*j_int_traj+1d0)/k_vec**2/
     & (2d0*j_int_ini+1)*pi
	  
!	  write(*,'(i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,i0,1x,f12.5,1x,
!     & f12.5)')k,j_int_traj, J_DOWN_INT,J_UP_INT,l_int_traj,j_int_ini, 
!     & id_traject,probab_J_all(k,i_traj,id_proc+1), part_cross
!      GOTO 2356

      IF(ident_skip) THEN
      IF(EVEN_NUM(l_int_traj)) THEN 
      part_cross = 	part_cross*posit_imp
      ELSE
      part_cross = 	part_cross*negat_imp	  
      ENDIF	  
      ENDIF
	 
2356      sigma_el(k) = sigma_el(k) +  part_cross     

      opac_chann_all(k,traj_works) = opac_chann_all(k,traj_works)
     &	+ part_cross

      j_def(traj_works) = j_int_traj
      ENDDO	  
      ENDDO
      ENDDO
	  

      IF(myid.eq.0) THEN
      OPEN(2,FILE="PARTIAL_CROSS_SECTIONS_J.out",POSiTION="APPEND") !!! CHECPOINT RENEWAL
      WRITE(2,'(a,i0)') "STATE= ", s_st
      WRITE(2,'(a,i0)') "j12= ",j12_s_st	  
      WRITE(2,'(a,i0)') "m12= ",m12_s_st	  
      WRITE(2,'(a,i0,a,f0.3,a)')
     & "ENERGY#",i_ener,", U(i)= ",
     & E_coll(i_ener)*conv_E_cm_Ev," cm^-1"
      DO j_chann = 1,nchann	  
      IF(j_chann.eq.1) WRITE(2,'(a7,5x)',ADVANCE="NO")
     & "J_TOTAL"
      WRITE(2,'(a8,i3,2x)',ADVANCE="NO") "CHANNEL=",j_chann
      ENDDO
      WRITE(2,*)	  
      DO i=1,tot_num_of_traj_actual
      DO j_chann = 1,nchann
      IF(j_chann.eq.1)
     & WRITE(2,'(i5,3x)',ADVANCE="NO")
     &	 INT(j_def(i))
      WRITE(2,'(e12.5,1x)',ADVANCE="NO") opac_chann_all(j_chann,i)
     & *a_bohr**2
      ENDDO
      WRITE(2,*)  
      ENDDO  
      CLOSE(2)
      OPEN(2,FILE="M_PROJ_CROSSECTIONS.out",POSiTION="APPEND")
      WRITE(2,*) "M RESOLVED CROSSECTIONS WITHOUT BILLING CORRECTION"	  
      WRITE(2,"(a10,1x,i4)") "INI_STATE#", s_st
      WRITE(2,'(a10,1x,i4)') "ENERGY#",i_ener
      WRITE(2,*) "j12=",j12_s_st	  
      WRITE(2,*) "m12=",m12_s_st	  
      WRITE(2,'(a8,1x,a19)') "CHANNEL#", "CROSSSECTIONS,ANG^2"
      DO k=1,nchann
      WRITE(2,'(i4,6x,e18.11)') k,sigma_el(k)*a_bohr**2
      ENDDO
      FLUSH(2)	  
      CLOSE(2)	  
      ENDIF	  
      opacity
     & (:,1:tot_num_of_traj_actual,i_ener) = 
     & opac_chann_all(:,:) 
	  
	  DEALLOCATE(opac_chann_all,j_def)	 
	  
! Bikram Start, Feb 2021:
	  if(transfer_prob_spln) then
	  if(.not. allocated(bk_prob1))
     & allocate(bk_prob1(nchann,tot_number_of_traject))
	  if(.not. allocated(bk_l_tmp))
     & allocate(bk_l_tmp(tot_number_of_traject))
	  if(.not. allocated(spln_der))
     & allocate(spln_der(nchann,tot_number_of_traject))
	  
      DO l_int_traj = L_MIN_TRAJECT,L_MAX_TRAJECT
	  
      CALL TRAJ_ORB(l_int_traj,
     & id_traject,
     & 1,
     & dl_step_integ,
     & ident_skip,
     & L_MIN_TRAJECT)
	 
      i_traj = itraj_myid_l(1,id_traject)	  
      i_ip = itraj_myid_l(2,id_traject)
      id_proc = i_ip*mpi_traject 
	  
	  CALL TRAJ_ORB(tmp_l,
     & id_traject,
     & -1,
     & dl_step_integ,
     & ident_skip,
     & L_MIN_TRAJECT)
	  
	  bk_prob1(:,id_traject) = probab_J_all(:,i_traj,id_proc+1)
	  bk_l_tmp(id_traject) = tmp_l*1.d0
	 
      ENDDO
	  
	  if(tot_number_of_traject.ne.size(bk_l_tmp)) then
      print*, 'Size Error in the Inelastic Probability Interpolation'
	  stop
	  end if
	  
	  do k = 1, nchann
	  spln_yp1=(bk_prob1(k,2)-bk_prob1(k,1))/(bk_l_tmp(2)-bk_l_tmp(1))
	  spln_yp2=(bk_prob1(k,tot_number_of_traject)-
     & bk_prob1(k,tot_number_of_traject-1))/
     & (bk_l_tmp(tot_number_of_traject)-
     & bk_l_tmp(tot_number_of_traject-1))
	  call spline(bk_l_tmp,bk_prob1(k,:),tot_number_of_traject,
     & spln_yp1,spln_yp2,spln_der(k,:))
	  end do
	  
      sigma_el = 0d0
      DO k=1,nchann	  
      DO j_int_traj = J_DOWN_INT,J_UP_INT
      DO l_int_traj = abs(j_int_traj-j_int_ini),j_int_traj+j_int_ini

      IF(ident_skip) THEN
      IF(parity_st.eq.1) THEN 
      posit_imp = exchange_ident
      negat_imp = 1d0 - exchange_ident	  
      ELSE
      posit_imp = 1d0 - exchange_ident	
      negat_imp = exchange_ident	  
      ENDIF	  
      ENDIF	  
	  
	  call splint(bk_l_tmp,bk_prob1(k,:),spln_der(k,:),
     & tot_number_of_traject,l_int_traj*1.d0,prob_reslt,spln_yp)
	  
      part_cross =
     & prob_reslt*(2d0*j_int_traj+1d0)/k_vec**2/(2d0*j_int_ini+1)

      IF(ident_skip) THEN
      IF(EVEN_NUM(l_int_traj)) THEN 
      part_cross = 	part_cross*posit_imp
      ELSE
      part_cross = 	part_cross*negat_imp	  
      ENDIF	  
      ENDIF
      sigma_el(k) = sigma_el(k) +  part_cross     

      ENDDO	  
      ENDDO
      ENDDO	  
	  if(allocated(bk_prob1)) deallocate(bk_prob1)
	  if(allocated(bk_l_tmp)) deallocate(bk_l_tmp)
	  if(allocated(spln_der)) deallocate(spln_der)
	  end if
! Bikram End.
	  
      END SUBROUTINE INELAST_CALC
	  
      INTEGER FUNCTION round(a)
! This function is updated by Bikramaditya Mandal
      IMPLICIT NONE
      REAL*8 a,b
      b = abs(a)
      if(b>int(b)+0.5d0) round = int(b) + 1
      if( b<int(b) + 0.5d0) round=int(b)
      if(a<0) round=-round   	  
      END FUNCTION round	  

      LOGICAL FUNCTION BELONGS(A,Arr,N)
! This function is updated by Bikramaditya Mandal
      INTEGER A,N
      INTEGER Arr(N)
      BELONGS = .FALSE.	  
      IF(A.ge.Arr(1) .and. A.le.Arr(N)) BELONGS = .TRUE.
      END 
      SUBROUTINE PRINT_CHECKPOINT
! This subroutine is updated by Bikramaditya Mandal
      USE INTEGRATOR
      USE CHECK_FILE
      USE MPI_DATA	  
      IMPLICIT NONE
      INTEGER nchann,tot_num_of_traj_actual  
      INTEGER i,i_ip,i_traj,k,j_chann,traj_works,j_int_traj,
     & l_int_traj,id_traject
      DO id_traject=1,tot_number_of_traject
      i_traj = itraj_myid_l(1,id_traject)	  
      i_ip = itraj_myid_l(2,id_traject)
      id_proc = i_ip*mpi_traject
      all_prob_l (:,:,id_traject) = sv_ch_prob_l
     & (:,:,i_traj,id_proc+1)
      what_computed(:,id_traject) = sv_ch_what_computed
     & (:,i_traj,id_proc+1)
      all_def_fnc(:,id_traject) = sv_ch_def_fnc(:,i_traj,id_proc+1)
      all_vib_fnc(:,id_traject) = sv_ch_vib_fnc(:,i_traj,id_proc+1)		  
      ENDDO		 
!      WRITE(*,*) id_traject,l_scatter(i_traj,id_proc+1)
      OPEN(33,FILE=CHECK_FILE_NAME,ACTION="WRITE",FORM="UNFORMATTED")
!!!! WRITING INFROMATION ABOUT SYSTEM	  
      WRITE(33) num_ini_states,ini_st_check_file	  
      WRITE(33) L_MIN_TRAJECT,L_MAX_TRAJECT,tot_number_of_traject
      WRITE(33) what_computed
      WRITE(33) all_def_fnc  
      WRITE(33) all_prob_l
      WRITE(33) i_ener
      WRITE(33) all_vib_fnc	  
      CLOSE(33)	  
  
!!!! END WRITING	  	 
  
      END SUBROUTINE PRINT_CHECKPOINT
	  
      SUBROUTINE PRINT_CHECKPOINT_MC
! This subroutine is updated by Bikramaditya Mandal
      USE INTEGRATOR
      USE MPI_DATA
      USE MPI	  
      USE CHECK_FILE
      USE VARIABLES	  
      USE MONTE_CARLO_STORAGE	  
      IMPLICIT NONE
      OPEN(33,FILE=CHECK_FILE_NAME,ACTION="WRITE",FORM="UNFORMATTED")
      WRITE(33) tot_number_of_traject_mc
      WRITE(33) i_ener
      ALLOCATE(x_mean_storage_ch(number_of_channels))
      ALLOCATE(x2_mean_storage_ch(number_of_channels))
      ALLOCATE(total_probab_check_monte_carlo(number_of_channels))
      x_mean_storage_ch=x_mean_storage
      x2_mean_storage_ch=x2_mean_storage
      total_probab_check_monte_carlo = total_probab	
      WRITE(33) x_mean_storage_ch
      WRITE(33) x2_mean_storage_ch
      WRITE(33) total_probab_check_monte_carlo	  
      CLOSE(33) 
!!!! END WRITING	  	 
  
      END SUBROUTINE PRINT_CHECKPOINT_MC	  
      FUNCTION MIDDLE_NUM(a,b,c)
! This function is updated by Bikramaditya Mandal
      REAL*8 a,b,c,	MIDDLE_NUM
      IF(a.le. max(b,c) .and. a.ge.min(b,c)) MIDDLE_NUM = a
      IF(c.le. max(b,a) .and. c.ge.min(b,a)) MIDDLE_NUM = c	
      IF(b.le. max(a,c) .and. b.ge.min(c,a)) MIDDLE_NUM = b		  
      END FUNCTION MIDDLE_NUM	  

      SUBROUTINE READ_CHECK_POINT_MC
! This subroutine is updated by Bikramaditya Mandal
      USE INTEGRATOR
      USE MPI_DATA
      USE MPI
      USE CHECK_FILE
      USE VARIABLES	  
      USE MONTE_CARLO_STORAGE	  
      IMPLICIT NONE
      LOGICAL file_exst	  
      INQUIRE( FILE=CHECK_FILE_NAME, EXIST=file_exst ) !! CHECKING IF FILE EXISTS
      IF(.not.file_exst) THEN
      IF(MYID.eq.0) PRINT*,"CHECK_FILE_NOT_FOUND"	  
      STOP	  
      ENDIF		  
      IF(myid.eq.0) THEN	  
      OPEN(33,FILE=CHECK_FILE_NAME,ACTION="READ",FORM="UNFORMATTED")
      READ(33) tot_number_to_cut
      READ(33) i_curr
      ENDIF
      	  
      ALLOCATE(x_mean_storage_ch(number_of_channels))
      ALLOCATE(x2_mean_storage_ch(number_of_channels))
      ALLOCATE(total_probab_check_monte_carlo(number_of_channels))
      IF(myid.eq.0)	 THEN 
      READ(33) x_mean_storage_ch
      READ(33) x2_mean_storage_ch
      READ(33) total_probab_check_monte_carlo	  
      CLOSE(33) 
!!!! END WRITING	  	 
      ENDIF 
      CALL MPI_BCAST(tot_number_to_cut, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(i_curr, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	 
      CALL MPI_BCAST(x_mean_storage_ch, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(x2_mean_storage_ch, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(total_probab_check_monte_carlo, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi) 
	 	  
      END SUBROUTINE READ_CHECK_POINT_MC

      SUBROUTINE TRAJ_ORB(l_momentum,id_of_traject,flag,dl,ident,
     & lmin)
	  use bk_l_values
! This subroutine is updated by Bikramaditya Mandal
      IMPLICIT NONE 
      INTEGER l_momentum,id_of_traject,flag,dl,lmin
      LOGICAL ident,EVEN_NUM
      IF(abs(flag).ne.1) STOP "ERROR: WRONG FLAG IN TRAJ_ORB"
      IF(.not.ident) THEN	  
      SELECT CASE(flag)
      CASE(1)
!      id_of_traject = (l_momentum-lmin)/dl + 1	  
      id_of_traject = (l_momentum-lmin)/dl + 1
	  if (l_switch_bk.gt.0d0) then
	  if(l_momentum.le.l_switch_bk)then
	  id_of_traject = (l_momentum-lmin)/dl + 1
	  else
	  id_of_traject = (l_momentum-l_switch_bk)/delta_l_lr + 1+l_range1
	  end if
	  end if
      CASE(-1)
!      l_momentum = (id_of_traject-1)*dl + lmin    
	  l_momentum = (id_of_traject-1d0/2d0)*dl + lmin
	  if (l_switch_bk.gt.0d0) then
	  if(id_of_traject.le.l_range1) then
      l_momentum = (id_of_traject-1d0/2d0)*dl + lmin
	  else
	  l_momentum = ((id_of_traject-l_range1)-1d0/2d0)*delta_l_lr 
     & + l_switch_bk
	  end if
	  end if
      END SELECT
      ELSE
      SELECT CASE(flag)
      CASE(1)
      id_of_traject = (l_momentum-lmin)/dl + 1
	  if (l_switch_bk.gt.0d0) then
	  if(l_momentum.le.l_switch_bk)then
	  id_of_traject = (l_momentum-lmin)/dl + 1
	  else
	  id_of_traject = (l_momentum-l_switch_bk)/delta_l_lr + 1+l_range1
	  end if
	  end if
      CASE(-1)
!      l_momentum = (id_of_traject-1)*dl + lmin
      l_momentum = (id_of_traject-1d0/2d0)*dl + lmin	
	  if (l_switch_bk.gt.0d0) then
	  if(id_of_traject.le.l_range1) then
      l_momentum = (id_of_traject-1d0/2d0)*dl + lmin
	  else
	  l_momentum = ((id_of_traject-l_range1)-1d0/2d0)*delta_l_lr 
     & + l_switch_bk
	  end if
	  end if
      END SELECT	  
      ENDIF	  
      END SUBROUTINE TRAJ_ORB	  
	  
      SUBROUTINE READ_CHECK_POINT
! This subroutine is updated by Bikramaditya Mandal
      USE INTEGRATOR
      USE CHECK_FILE
      USE MPI_DATA
      USE MPI
      IMPLICIT NONE
      LOGICAL file_exst	  
      INQUIRE( FILE=CHECK_FILE_NAME, EXIST=file_exst ) !! CHECKING IF FILE EXISTS
      IF(.not.file_exst) THEN
      IF(MYID.eq.0) PRINT*,"CHECK_FILE_NOT_FOUND"	  
      STOP	  
      ENDIF	  
      IF(MYID.EQ.0) THEN	  
      OPEN(33,FILE=CHECK_FILE_NAME,ACTION="READ",FORM="UNFORMATTED")
!!!! READING INFROMATION ABOUT SYSTEM	  
      READ(33) num_ini_states,ini_st_check_file	  
      READ(33) L_MIN_TRAJECT,L_MAX_TRAJECT,tot_number_of_traject
!! CREATING ALLOCATION	  
      ALLOCATE(what_computed(num_ini_states,tot_number_of_traject))	  
      ALLOCATE(      
     & all_prob_l
     & (sys_var_size+6,num_ini_states,tot_number_of_traject))
      ALLOCATE(all_def_fnc(num_ini_states,tot_number_of_traject))
      ALLOCATE(all_vib_fnc(num_ini_states,tot_number_of_traject))	  
      READ(33) what_computed
      READ(33) all_def_fnc  
      READ(33) all_prob_l
      READ(33) i_curr
      READ(33) all_vib_fnc  	  
      CLOSE(33)
      ENDIF
      CALL MPI_BCAST(num_ini_states, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(i_curr, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	 
      CALL MPI_BCAST(L_MIN_TRAJECT, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(L_MAX_TRAJECT, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(tot_number_of_traject, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(ini_st_check_file, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	 
	 
      IF(MYID.NE.0) THEN
      ALLOCATE(what_computed(num_ini_states,tot_number_of_traject))	  
      ALLOCATE(      
     & all_prob_l
     & (sys_var_size+6,num_ini_states,tot_number_of_traject))
      ALLOCATE(all_def_fnc(num_ini_states,tot_number_of_traject))
      ALLOCATE(all_vib_fnc(num_ini_states,tot_number_of_traject))	  
      ENDIF
      task_size  = (sys_var_size+6)* 
     & num_ini_states*tot_number_of_traject
      CALL MPI_BCAST(all_prob_l, task_size, MPI_DOUBLE_PRECISION,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      task_size  = 
     & num_ini_states*tot_number_of_traject
      CALL MPI_BCAST(what_computed, task_size, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(all_def_fnc, task_size, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(all_vib_fnc, task_size, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	 
      END SUBROUTINE READ_CHECK_POINT	  
      SUBROUTINE ALLOCATE_ARRAYS
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE CONSTANTS
	  USE CHECK_FILE
      USE INTEGRATOR
      USE MPI_DATA
      USE MPI
      IMPLICIT NONE
      INTEGER 	 max_numb_trajec 
      IF(i_ener.eq.i_curr) THEN
	  CALL DEFINE_MAX_NUM_TRAJ(max_numb_trajec)
      ALLOCATE(TRAJECT_END_DATA
     & (7,max_numb_trajec,num_ini_states,nmbr_of_enrgs))
      TRAJECT_END_DATA  = 0d0	 
      IF(myid.eq.0) THEN	  
      ALLOCATE(TRAJECT_END_DATA_ALL
     & (7,max_numb_trajec,num_ini_states,nmbr_of_enrgs,nproc))
      ENDIF	 	  
      ENDIF	  
      ALLOCATE(probab_J(number_of_channels,n_traject_alloc))		  
      ALLOCATE(probab_J_all(number_of_channels,n_traject_alloc,nproc))
      ALLOCATE(err_ener_larg(nproc),err_prb_larg(nproc))
!	  write(*,'(a,1x,i0,1x,i0,1x,i0)')'angle_scatter', myid, 
!     & n_traject_alloc,nproc
      ALLOCATE(ampl_wf_real(n_traject_alloc,nproc),
     & ampl_wf_imag(n_traject_alloc,nproc),
     & angle_phase(n_traject_alloc,nproc),
     & l_scatter(n_traject_alloc,nproc),
     & angle_scatter(n_traject_alloc,nproc)

! Bikram Start:
     & ,bk_tym(n_traject_alloc,nproc)
     & ,bk_prd(n_traject_alloc,nproc)
     & ,bk_vib(n_traject_alloc,nproc)
     & ,bk_erre(n_traject_alloc,nproc)
     & ,bk_errp(n_traject_alloc,nproc)
! Bikram End.

     & ,phase_wf(n_traject_alloc,nproc),
     & buffer_elastic(n_traject_alloc))
      ALLOCATE(opacity
     & (number_of_channels,J_UP_INT - J_DOWN_INT + 1,nmbr_of_enrgs))	 
      opacity = 0d0
      angle_scatter = 0d0

! Bikram Start:
	  bk_tym = 0d0
	  bk_prd = 0d0
	  bk_vib = 0d0
	  bk_erre = 0d0
	  bk_errp = 0d0 
! Bikram End.
	  
      buffer_elastic = 0d0
      probab_J = 0d0
      probab_J_all = 0d0
      ampl_wf_real = 0d0
      ampl_wf_imag = 0d0
      l_scatter = 0d0
      phase_wf = 0d0

      	  
      IF(orbit_traj_defined) THEN
      ALLOCATE(ORBIT_TRAJ_DATA(10,n_traject_alloc))
      ALLOCATE(ORBIT_TRAJ_FLAG(n_traject_alloc))
      IF(myid.eq.0) THEN
      ALLOCATE(ORBIT_TRAJ_DATA_ALL(10,n_traject_alloc,nproc))
      ALLOCATE(ORBIT_TRAJ_FLAG_ALL(n_traject_alloc,nproc))
      ENDIF	  
      ORBIT_TRAJ_DATA = 0d0
      DO itraject = 1,ntraject      	  
      ORBIT_TRAJ_FLAG(itraject) = .FALSE.
      ENDDO	  
      ENDIF	  	  
      END SUBROUTINE ALLOCATE_ARRAYS
      SUBROUTINE DEALLOCATE_ARRAYS
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE CONSTANTS
      USE INTEGRATOR
      USE CHECK_FILE
      USE MONTE_CARLO_STORAGE	  
      USE MPI_DATA
      USE MPI
      IF(ALLOCATED(probab_J))DEALLOCATE(probab_J)		  
      IF(ALLOCATED(probab_J_all))DEALLOCATE(probab_J_all)
      IF(ALLOCATED(err_ener_larg))DEALLOCATE(err_ener_larg)
      IF(ALLOCATED(err_prb_larg))DEALLOCATE(err_prb_larg)
!      IF(ALLOCATED(TRAJECT_END_DATA))DEALLOCATE(TRAJECT_END_DATA)
!      IF(ALLOCATED(TRAJECT_END_DATA_ALL).and.myid.eq.0) THEN
!	  DEALLOCATE(TRAJECT_END_DATA_ALL)
!      ENDIF	  
      IF(ALLOCATED(ampl_wf_real)) THEN	  
      DEALLOCATE(ampl_wf_real,angle_phase,
     & ampl_wf_imag,
     & l_scatter,
     & angle_scatter
	 
! Bikram Start:
     & ,bk_tym
     & ,bk_prd
     & ,bk_vib
     & ,bk_erre
     & ,bk_errp
! Bikram End.
	 
     & ,phase_wf,
     & buffer_elastic)
      ENDIF	
      DEALLOCATE(opacity
     & )
      IF(MYID.eq.0) THEN
      IF(ALLOCATED(what_computed)) THEN	  
      DEALLOCATE(what_computed,all_prob_l,
     & all_def_fnc,all_vib_fnc)
      ENDIF
      IF(ALLOCATED( sv_ch_def_fnc)) THEN
      DEALLOCATE(sv_ch_prob_l,sv_ch_what_computed,sv_ch_def_fnc)	  
      ENDIF	  
      ENDIF
      IF(ALLOCATED(sv_ch_prob_l_buffer)) THEN	  
      DEALLOCATE(sv_ch_prob_l_buffer,sv_ch_what_computed_buffer,
     & sv_ch_def_fnc_buffer,sv_ch_vib_fnc_buffer)
      ENDIF
      IF(ALLOCATED(TRAJECT_DATA_CHECK_FILE))
     & DEALLOCATE(TRAJECT_DATA_CHECK_FILE)	  
      IF(ALLOCATED(itraj_myid_l)) DEALLOCATE(itraj_myid_l)
      IF(ALLOCATED(x2_mean_storage)) THEN
      DEALLOCATE(x_mean_storage)
      DEALLOCATE(x2_mean_storage)
      DEALLOCATE(mean_values)
      DEALLOCATE(variance_values)	  
      ENDIF	 
      IF(orbit_traj_defined .and. ALLOCATED(ORBIT_TRAJ_DATA) ) THEN
      DEALLOCATE(ORBIT_TRAJ_DATA)
      DEALLOCATE(ORBIT_TRAJ_FLAG)
      IF(myid.eq.0) DEALLOCATE(ORBIT_TRAJ_DATA_ALL)
      IF(myid.eq.0)DEALLOCATE(ORBIT_TRAJ_FLAG_ALL)	  
      ENDIF	 	  
      END SUBROUTINE DEALLOCATE_ARRAYS
      SUBROUTINE TIME_CHECKER(RETURN_VAL)
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE CHECK_FILE	  
      USE INTEGRATOR	  
      USE MPI
      USE MPI_DATA
      IMPLICIT NONE
      LOGICAL	RETURN_VAL  
      RETURN_VAL = .FALSE.	  
      IF(write_check_file_defined) THEN
      TIME_DOTRAJECT = MPI_Wtime()
      IF((TIME_DOTRAJECT-TIME_BEGIN).gt.60d0*TIME_MIN_CHECK) THEN
      RETURN_VAL = .TRUE.	  
      RETURN
      ENDIF      
      ENDIF
      RETURN	  
      END SUBROUTINE TIME_CHECKER
      
      SUBROUTINE DEFINE_MAX_NUM_TRAJ(max_numb_trajec)
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE CONSTANTS
      USE CHECK_FILE	  
      USE INTEGRATOR
      USE MPI_DATA
      USE MPI
      USE ERRORS	  
      USE MONTE_CARLO_STORAGE	  
      IMPLICIT NONE
      INTEGER max_numb_trajec,j_count
	  integer l_switch_bk_tmp, l1_tmp, l2_tmp										! Bikram Feb 2021
      REAL*8 j_bound_up,j_bound_d	
      REAL*8 max_energ	  
	  IF(monte_carlo_defined) THEN
      max_numb_trajec = int(tot_number_of_traject/(nproc/mpi_traject))+1
      ELSE
	  max_energ = E_coll(1)
      DO 	j_count=1,nmbr_of_enrgs
      IF(	max_energ<  E_coll(j_count)) max_energ = E_coll(j_count)
      ENDDO
!!!! IF NOT MONTE CARLO DEFINED	  
      IF(b_impact_defined) THEN
      j_bound_up = int(b_impact_parameter*sqrt(massAB_M*2d0*max_energ))
      j_bound_d = 0
	  ELSE
      
       j_bound_d = J_tot_min
       j_bound_up  = J_tot_max
      ENDIF
      max_numb_trajec=
     & int(j_bound_up - j_bound_d+jmax_included)/delta_l_step + 1
	  
! Bikram Start: Feb 2021, regarding the L_switch
	  if(bk_b_switch.gt.0.d0) then
	  l_switch_bk_tmp = int(bk_b_switch*sqrt(massAB_M*2d0*max_energ))
	  l1_tmp = (l_switch_bk_tmp-j_bound_d+jmax_included)/delta_l_step + 1
	  l2_tmp = (j_bound_up - l_switch_bk_tmp)/bk_dl_lr + 1
	  max_numb_trajec = l1_tmp + l2_tmp
	  end if
      ENDIF
      END  SUBROUTINE DEFINE_MAX_NUM_TRAJ     
      FUNCTION EVEN_NUM(n)
! This function is updated by Bikramaditya Mandal
      INTEGER n
      LOGICAL EVEN_NUM
      EVEN_NUM = .FALSE.	  
      IF((n/2)*2 .eq. n) EVEN_NUM = .TRUE.	  
      END FUNCTION EVEN_NUM	  
      SUBROUTINE ELAST_CALC(sigma_el,n_angl_size,print_deflect,
     & diff_def,bk_nchanl)
! This subroutine is updated by Bikramaditya Mandal
      USE INTEGRATOR
      USE MPI	  
      USE MPI_DATA
      USE MPI_TASK_TRAJECT	  
	  use bk_l_values												!Bikram Feb 2021
      IMPLICIT NONE
      LOGICAL print_deflect,diff_def	  
      INTEGER k,i_traj,id_traject,i_ip,l_int_traj,traj_works,j_int_traj
      INTEGER n_angl_size
      INTEGER k_angl_1,k_angl_2,d_angl_step,d_angle_div
      INTEGER l1_temp,l2_temp,dl_temp	  
      REAL*8 plgndr,sigma_el
      EXTERNAL 	 plgndr
      COMPLEX*16 f_scatter_ampl(n_angl_size)
      REAL*8 diff_cross(n_angl_size),diff_cross_buffer(n_angl_size)	  
      REAL*8 angle_array(n_angl_size),th_arr
      REAL*8 elast_probab
      REAL*8, PARAMETER :: a_bohr = 5.2917721067d-1
      REAL*8, PARAMETER :: autown = 27.21113845d0*8065.5448d0
	  integer bkn,bk_nchanl			!Bikram Feb 2020
	  real*8 bk_probab_sum_tmp			!Bikram Feb 2020
	  real*8, allocatable :: bk_probab_tmp(:,:)			!Bikram Feb 2020
      logical monte_carlo_defined	  !Bikram Oct'18
	  	  
	  real*8, allocatable :: bk_l_tmp(:), bk_phase(:)
	  integer tmp_l, tmp_l1
	  real*8, allocatable :: bk_prob1(:,:), spln_der(:,:)
	  real*8, allocatable :: bk_def1(:), spln_der_def(:)
	  real*8 spln_yp1, spln_yp2, spln_yp, prob_reslt
	  real*8 def_reslt1, def_reslt2
	  real*8, parameter :: bk_cm_Ev = 27.21113845d0*8065.5448d0
	  	  
      ALLOCATE(
     & phase_scatter(tot_number_of_traject),
     & def_fnc(tot_number_of_traject))
	 
! Bikram Start:
	  ALLOCATE(bk_time(tot_number_of_traject))
	  ALLOCATE(bk_period(tot_number_of_traject))
	  ALLOCATE(bk_vibration(tot_number_of_traject))
	  ALLOCATE(bk_err_E(tot_number_of_traject))
	  ALLOCATE(bk_err_P(tot_number_of_traject))
	  ALLOCATE(bk_probab_tmp(bk_nchanl,tot_number_of_traject))
! Bikram End.
	 
      def_fnc = 0d0	  
      diff_cross_buffer =  0d0
      diff_cross = 0d0	  
      phase_scatter =  0d0
      sigma_el = 0d0

! Bikram Start:
      bk_time = 0d0
      bk_period = 0d0
      bk_vibration = 0d0
      bk_err_E = 0d0
      bk_err_P = 0d0
      bk_probab_tmp = 0d0
! Bikram End.     

!!! CRATE DEF FUNCTION ARRAY	  
      DO l_int_traj = L_MIN_TRAJECT,L_MAX_TRAJECT
!      id_traject = (l_int_traj-L_MIN_TRAJECT)/dl_step_integ+1
      CALL TRAJ_ORB(l_int_traj,
     & id_traject,
     & 1,
     & dl_step_integ,
     & ident_skip,
     & L_MIN_TRAJECT)
	 
      i_traj = itraj_myid_l(1,id_traject)	  
      i_ip = itraj_myid_l(2,id_traject)
      id_proc = i_ip*mpi_traject 
	  
      def_fnc(id_traject)  = angle_scatter(i_traj,id_proc+1)
! Bikram Start:
	  bk_time(id_traject)  = bk_tym(i_traj,id_proc+1)
	  bk_period(id_traject)  = bk_prd(i_traj,id_proc+1)
	  bk_vibration(id_traject)  = bk_vib(i_traj,id_proc+1)
	  bk_err_E(id_traject)  = bk_erre(i_traj,id_proc+1)
	  bk_err_P(id_traject)  = bk_errp(i_traj,id_proc+1)
	  bk_probab_tmp(:,id_traject) = probab_J_all(:,i_traj,id_proc+1)
! Bikram End
	  
      ENDDO	  
	  
	  
!!!! COMPUTING SCATTERING PHASE

      DO id_traject= tot_number_of_traject-1,1,-1 
      CALL TRAJ_ORB(l1_temp,
     & id_traject,
     & -1,
     & dl_step_integ,
     & ident_skip,
     & L_MIN_TRAJECT)	  
      CALL TRAJ_ORB(l2_temp,
     & id_traject+1,
     & -1,
     & dl_step_integ,
     & ident_skip,
     & L_MIN_TRAJECT)
      dl_temp = abs(l2_temp-l1_temp)	 
      phase_scatter(id_traject) = phase_scatter(id_traject+1)
     & - (def_fnc(id_traject)+def_fnc(id_traject+1))/2d0*dl_temp
      i_traj = itraj_myid_l(1,id_traject)	  
      i_ip = itraj_myid_l(2,id_traject)
      id_proc = i_ip*mpi_traject 	  
      elast_probab     =  probab_J_all(ini_channel,i_traj,id_proc+1)
      sigma_el = sigma_el + abs(1d0-abs(dsqrt(elast_probab))
     & *dcmplx(dcos(phase_scatter(id_traject)),
     & dsin(phase_scatter(id_traject))))**2/k_vec**2*pi*(2d0*l1_temp+1)
     & * dl_temp	 
	 
!	  write(*,'(i0,1x,i0,1x,f12.5,1x,f12.5,1x,f12.5,1x,f12.5)')
!     & id_traject,dl_temp,def_fnc(id_traject),def_fnc(id_traject+1),
!     & phase_scatter(id_traject),elast_probab
	  
      ENDDO
!!! SCATTERING AMPLITUTE	  
      f_scatter_ampl = (0d0,0d0)
      IF(diff_def) THEN 
      sigma_el = 0d0	  
      d_angl_step = int(n_angl_size/nproc)
      d_angle_div = n_angl_size - d_angl_step*nproc
      IF(myid.lt.d_angle_div) THEN	  
      k_angl_1 = myid*(d_angl_step+1) + 1
      k_angl_2 = (myid+1)*(d_angl_step+1)
      ELSE
      k_angl_1 = myid*d_angl_step + 1+d_angle_div
      k_angl_2 = (myid+1)*d_angl_step+d_angle_div	  
      ENDIF
      DO k =1,n_angl_size
      th_arr= (dble(k)-1d0)*dacos(-1d0)/(dble(n_angl_size)-1d0)
      angle_array(k) = th_arr
      ENDDO	  
      DO k =k_angl_1,k_angl_2
      th_arr= angle_array(k)
      DO j_int_traj = J_DOWN_INT,J_UP_INT
      DO l_int_traj = abs(j_int_traj-j_int_ini),j_int_traj+j_int_ini
      CALL TRAJ_ORB(l_int_traj,
     & id_traject,
     & 1,
     & dl_step_integ,
     & ident_skip,
     & L_MIN_TRAJECT)	  
!	  if(abs(itraj_myid_l(1,id_traject)).gt.1) 
!     & itraj_myid_l(1,id_traject) = 1
!	  if(abs(itraj_myid_l(2,id_traject)).gt.0) 
!     & itraj_myid_l(2,id_traject) = 0
!      id_traject =(l_int_traj-L_MIN_TRAJECT)/dl_step_integ+1
      i_traj = itraj_myid_l(1,id_traject)	  
      i_ip = itraj_myid_l(2,id_traject)
      id_proc = i_ip*mpi_traject 	  
      elast_probab     =  probab_J_all(ini_channel,i_traj,id_proc+1)	  
      f_scatter_ampl(k) =f_scatter_ampl(k) +
     & plgndr(l_int_traj,0,cos(th_arr))*
     & (1d0-abs(dsqrt(elast_probab))
     & *dcmplx(dcos(phase_scatter(id_traject)),
     & dsin(phase_scatter(id_traject))))*
     & (2d0*j_int_traj+1d0)*dcmplx(0d0,1d0)/k_vec/2d0/
     & (2d0*j_int_ini+1)
      diff_cross(k) = abs(f_scatter_ampl(k))**2	 
      ENDDO
      ENDDO		  
      ENDDO
      CALL MPI_Reduce(diff_cross, diff_cross_buffer,
     & n_angl_size, MPI_REAL8,
     & MPI_SUM, 0,	 
     & MPI_COMM_WORLD,ierr_mpi)
      diff_cross  = diff_cross_buffer 
      CALL MPI_Bcast(diff_cross,n_angl_size,
     & MPI_REAL8, 0,MPI_COMM_WORLD,ierr_mpi) 
      DO k=1,n_angl_size  
      sigma_el = sigma_el +  
     & (diff_cross(k-1)*sin(angle_array(k-1))
     & +diff_cross(k)*sin(angle_array(k))
     & )* dacos(-1d0)*(angle_array(k)-angle_array(k-1))
      ENDDO	
      ENDIF	  
      IF(print_deflect .and. myid.eq.0) THEN
	  
! Bikram Start:
	  open(22,file = 'OPACITY_FUNCTIONS_L.out',position = 'append')
	  open(23,file = 'TOTAL_PROBABILITIES_L.out',position = 'append')
	  
      WRITE(23,'(a,i0)') "STATE= ", s_st
      WRITE(23,'(a,i0)') "j12= ",j12_s_st	  
      WRITE(23,'(a,i0)') "m12= ",m12_s_st	  
      WRITE(23,'(a,i0,a,f0.3,a)')
     & "ENERGY#",i_ener,", U(i)= ",
     & E_coll(i_ener)*bk_cm_Ev," cm^-1"
      WRITE(23,'(a8,7x,a7,1x,a15,1x,a16)',
     & ADVANCE="NO") "B_IMPACT", "     L ", "TOT_ELAST_PROB", 
     & "TOT_INELAST_PROB"
      WRITE(23,*)	  
	  
      WRITE(22,'(a,i0)') "STATE= ", s_st
      WRITE(22,'(a,i0)') "j12= ",j12_s_st	  
      WRITE(22,'(a,i0)') "m12= ",m12_s_st	  
      WRITE(22,'(a,i0,a,f0.3,a)')
     & "ENERGY#",i_ener,", U(i)= ",
     & E_coll(i_ener)*bk_cm_Ev," cm^-1"
      DO bkn = 1, bk_nchanl
      IF(bkn.eq.1) WRITE(22,'(a8,7x,a7,5x)',ADVANCE="NO")
     & "B_IMPACT", "     L "    	 
      WRITE(22,'(a8,i6,1x)',ADVANCE="NO") "CHANNEL=",bkn
      ENDDO
      WRITE(22,*)	  
	  
	  do id_traject = 1, tot_number_of_traject
	  bk_probab_sum_tmp = 0.d0
	  do bkn = 1, bk_nchanl
	  bk_probab_sum_tmp = bk_probab_sum_tmp + 
     & bk_probab_tmp(bkn,id_traject)
	  if(dl_step_integ > 1) then
	  tmp_l1 = int((2d0*id_traject-1)*dl_step_integ/2d0+L_MIN_TRAJECT)
	  else
	  tmp_l1 = int((id_traject-1)*dl_step_integ+L_MIN_TRAJECT)
	  end if
	  if(bkn.eq.1) then
	  write(22,'(e15.8,1x,i6,1x,e15.8,1x)',ADVANCE="NO") 
     & tmp_l1/k_vec, tmp_l1,
     & bk_probab_tmp(bkn,id_traject)
	  else
	  write(22,'(e15.8,1x)',ADVANCE="NO") 
     & bk_probab_tmp(bkn,id_traject)
	  endif
	  enddo
	  write(22,*)
	  
	  write(23,'(e15.8,1x,i6,1x,e15.8,1x,e15.8)')
     & tmp_l1/k_vec, tmp_l1,
     & bk_probab_sum_tmp, 
     & bk_probab_sum_tmp - bk_probab_tmp(ini_channel,id_traject)
	  enddo
	  write(22,*)
	  write(23,*)
	  close(22)
	  close(23)
      OPEN(2,FILE="DEFLECTION_FUNCTION.out",POSiTION = "APPEND")
      WRITE(2,'(a,i0,a,f0.3,a)') 
     & "energy#",i_ener,", (U= ",E_sct*autown,' cm^-1)'	   
      WRITE(2,'(a,i0)') "state= ", s_st
      WRITE(2,'(a,i0)') "j12= ",j12_s_st	  
      WRITE(2,'(a,i0)') "m12= ",m12_s_st
      WRITE(2,'(a5,2x,a14,5x,a17,4x,a15,7x,a17,4x,a15,6x,a13
     & ,4x,a5,3x,a5)')
     & "L_ORB","B_IMPACT(a.u.)","DEFLCT_ANGLE(rad)","SCAT_PHASE(rad)" 	
     & ,"FINISH_TIME(a.u.)","ENERGY_ERROR(%)","NORM_ERROR(%)"
     & ,"LOOPS","OSCIL"	 
      DO id_traject = 1,tot_number_of_traject
	  if(dl_step_integ > 1) then
	  tmp_l1 = int((2d0*id_traject-1)*dl_step_integ/2d0+L_MIN_TRAJECT)
	  else
	  tmp_l1 = int((id_traject-1)*dl_step_integ+L_MIN_TRAJECT)
	  end if
      WRITE(2,'(i4,2x,e15.8,1x,e19.8,2x,e19.8,3x,e19.8,1x,e19.8,1x
     & ,e19.8,3x,i5,3x,i5)')
     & tmp_l1, tmp_l1/k_vec,
     & def_fnc(id_traject),phase_scatter(id_traject)
     & ,bk_time(id_traject),bk_err_E(id_traject)
     & ,bk_err_P(id_traject),bk_period(id_traject)
     & ,bk_vibration(id_traject)
      ENDDO
      WRITE(2,*)	  
      CLOSE(2)
! Bikram End.
	  
      IF(diff_def) THEN  
      OPEN(2,FILE="DIFF_CROSS.out", POSiTION = "APPEND")
      WRITE(2,*) "energy",i_ener	  
      WRITE(2,*) "state = ", s_st	  
      WRITE(2,*) "j12=",j12_s_st	  
      WRITE(2,*) "m12=",m12_s_st
      WRITE(2,'(a19,1x,a19)') "ANGLE,deg","DIFF. CROSS.,angs^2" 	  
      DO k = 1,n_angl_size
      WRITE(2,'(e19.7,1x,e19.7)')
     & angle_array(k)*180d0/dacos(-1d0), diff_cross(k)*a_bohr**2
      ENDDO
      WRITE(2,*)	  
      CLOSE(2)	
      ENDIF	  
      ENDIF
      DEALLOCATE(phase_scatter,def_fnc)

! Bikram Start:
	  DEALLOCATE(bk_time)
	  DEALLOCATE(bk_period)
	  DEALLOCATE(bk_vibration)
	  DEALLOCATE(bk_err_E)
	  DEALLOCATE(bk_err_P)
	  DEALLOCATE(bk_probab_tmp)
! Bikram End

! Bikram Start, Feb 2021:
	  if(transfer_prob_spln) then
	  if(.not. allocated(bk_prob1))
     & allocate(bk_prob1(bk_nchanl,tot_number_of_traject))
	  if(.not. allocated(bk_def1))
     & allocate(bk_def1(tot_number_of_traject))
	  if(.not. allocated(bk_l_tmp))
     & allocate(bk_l_tmp(tot_number_of_traject))
	  if(.not. allocated(spln_der))
     & allocate(spln_der(bk_nchanl,tot_number_of_traject))
	  if(.not. allocated(spln_der_def))
     & allocate(spln_der_def(tot_number_of_traject))
	  if(.not. allocated(bk_phase))
     & allocate(bk_phase(L_MAX_TRAJECT-L_MIN_TRAJECT+1))

      DO l_int_traj = L_MIN_TRAJECT,L_MAX_TRAJECT
	  
      CALL TRAJ_ORB(l_int_traj,
     & id_traject,
     & 1,
     & dl_step_integ,
     & ident_skip,
     & L_MIN_TRAJECT)
	 
      i_traj = itraj_myid_l(1,id_traject)	  
      i_ip = itraj_myid_l(2,id_traject)
      id_proc = i_ip*mpi_traject 
	  
	  CALL TRAJ_ORB(tmp_l,
     & id_traject,
     & -1,
     & dl_step_integ,
     & ident_skip,
     & L_MIN_TRAJECT)
	  
	  bk_prob1(:,id_traject) = probab_J_all(:,i_traj,id_proc+1)
	  bk_def1(id_traject) = angle_scatter(i_traj,id_proc+1)	  
	  bk_l_tmp(id_traject) = tmp_l*1.d0
	 
      ENDDO
	  
	  if(tot_number_of_traject.ne.size(bk_l_tmp)) then
      print*, 'Size Error in the Inelastic Probability Interpolation'
	  stop
	  end if
	  
	  do k = 1, bk_nchanl
	  spln_yp1=(bk_prob1(k,2)-bk_prob1(k,1))/(bk_l_tmp(2)-bk_l_tmp(1))
	  spln_yp2=(bk_prob1(k,tot_number_of_traject)-
     & bk_prob1(k,tot_number_of_traject-1))/
     & (bk_l_tmp(tot_number_of_traject)-
     & bk_l_tmp(tot_number_of_traject-1))
	  call spline(bk_l_tmp,bk_prob1(k,:),tot_number_of_traject,
     & spln_yp1,spln_yp2,spln_der(k,:))
	  
	  spln_yp1=(bk_def1(2)-bk_def1(1))/(bk_l_tmp(2)-bk_l_tmp(1))
	  spln_yp2=(bk_def1(tot_number_of_traject)-
     & bk_def1(tot_number_of_traject-1))/
     & (bk_l_tmp(tot_number_of_traject)-
     & bk_l_tmp(tot_number_of_traject-1))
	  call spline(bk_l_tmp,bk_def1,tot_number_of_traject,
     & spln_yp1,spln_yp2,spln_der_def)
	  end do
	  
	  bk_phase = 0.d0
	  sigma_el = 0.d0
      DO id_traject= L_MAX_TRAJECT-1,L_MIN_TRAJECT,-1
	  call splint(bk_l_tmp,bk_def1,spln_der_def,
     & tot_number_of_traject,id_traject*1.d0,def_reslt1,spln_yp)
	  call splint(bk_l_tmp,bk_def1,spln_der_def,
     & tot_number_of_traject,(id_traject+1)*1.d0,def_reslt2,spln_yp)
	  call splint(bk_l_tmp,bk_prob1(ini_channel,:),
     & spln_der(ini_channel,:),
     & tot_number_of_traject,id_traject*1.d0,prob_reslt,spln_yp)
	  
      bk_phase(id_traject+1) = bk_phase(id_traject+2)
     & - (def_reslt1+def_reslt2)/2d0	
      elast_probab = prob_reslt

      sigma_el = sigma_el + abs(1d0-abs(dsqrt(elast_probab))
     & *dcmplx(dcos(bk_phase(id_traject+1)),
     & dsin(bk_phase(id_traject+1))))**2
     & /k_vec**2*pi*(2d0*id_traject+1)

      ENDDO
	  
	  open(25,file = 'INTERPOLATED_PROBS_DEFLECTION_L.out')
      DO bkn = 1, bk_nchanl
      IF(bkn.eq.1) WRITE(25,'(a8,7x,a7,5x,a12,2x)',ADVANCE="NO")
     & "B_IMPACT", "     L ", "DEF_FUNCTION"    	 
      WRITE(25,'(a8,i6,2x)',ADVANCE="NO") "CHANNEL=",bkn
      ENDDO
      WRITE(25,*)
	  
	  do id_traject = L_MIN_TRAJECT, L_MAX_TRAJECT
	  do bkn = 1, bk_nchanl
	  call splint(bk_l_tmp,bk_def1,spln_der_def,
     & tot_number_of_traject,id_traject*1.d0,def_reslt1,spln_yp)
	  call splint(bk_l_tmp,bk_prob1(bkn,:),
     & spln_der(bkn,:),
     & tot_number_of_traject,id_traject*1.d0,prob_reslt,spln_yp)
	  if(bkn.eq.1) then
	  write(25,'(e15.8,1x,i6,1x,f12.4,5x,e15.8,1x)',ADVANCE="NO") 
     & dble(id_traject + L_MIN_TRAJECT )/k_vec,
     &  int(id_traject + L_MIN_TRAJECT ), def_reslt1,
     & prob_reslt
	  else
	  write(25,'(e15.8,1x)',ADVANCE="NO") prob_reslt
	  endif
	  enddo
	  write(25,*)
	  enddo
	  write(25,*)
	  close(25)
	  
	  if(allocated(bk_prob1)) deallocate(bk_prob1)
	  if(allocated(bk_def1)) deallocate(bk_def1)
	  if(allocated(bk_l_tmp)) deallocate(bk_l_tmp)
	  if(allocated(spln_der)) deallocate(spln_der)
	  if(allocated(spln_der_def)) deallocate(spln_der_def)
	  if(allocated(bk_phase)) deallocate(bk_phase)
	  end if
! Bikram End.
	  
      END SUBROUTINE ELAST_CALC
      SUBROUTINE VARIANCE_COMPUTE(x_values,s_value,n,j)
! This subroutine is updated by Bikramaditya Mandal
      USE MONTE_CARLO_STORAGE
      IMPLICIT NONE
      INTEGER n,i,j
      REAL*8 x_mean,x_values,s_value,x2_mean
      x_mean = 0d0
      x2_mean = 0d0
      x_mean_storage(j) = x_mean_storage(j) + x_values
      x_mean =  x_mean_storage(j) /dble(n+tot_number_to_cut)
      x2_mean_storage(j) = x2_mean_storage(j)+ x_values**2 	  
      x2_mean = x2_mean_storage(j)/dble(n+tot_number_to_cut)
      s_value = dsqrt(x2_mean-x_mean**2)
	  
      END	SUBROUTINE VARIANCE_COMPUTE 

      SUBROUTINE COMPUTE_V_TERMS(R)
! This subroutine is updated by Bikramaditya Mandal
      USE INTEGRATOR
      USE VARIABLES
      USE ERRORS
      USE MPI_DATA	  
      IMPLICIT NONE
      INTEGER i  
      REAL*8 R	    
      IF(.not.term_pot_defined) THEN
      show_error = 
     & "ERROR:USER TERMS POTENTIAL WAS SUPPOSED TO BE DEFINED"
      error_compute_exe = .TRUE.
      IF(MYID.eq.0) PRINT*,show_error	  
      STOP	  
      ENDIF
  
      DO I=1,nterms
      CALL USER_DEFINED_COEFFS(V_TERMS(i),V_TERMS_DER(i),i,R)	  
      ENDDO	  
      END SUBROUTINE COMPUTE_V_TERMS

      SUBROUTINE DEFINE_WIDTH(s_ini)
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE CONSTANTS
      USE CHECK_FILE	  
      USE INTEGRATOR
      USE MPI_DATA
      USE MPI
      USE ERRORS	  
      USE MONTE_CARLO_STORAGE	  
      IMPLICIT NONE
      LOGICAL sampl_succ	  
      INTEGER s_ini,m_t,j_t,j1_t,j2_t,j12_min,j12_max,j_count,j_summ,st
      INTEGER i,j,k,l_parity,p_parity,ident_max,KRONEKER  
	  m_t = 0!!! DEFINING HOW MANY DIFFERENT J12 and M12 WE have
      j_t =  j_min_ind(chann_ini)	  
      IF(.not.identical_particles_defined) THEN
      i = indx_corr(m_t + j_t+1,
     & j_t+1,chann_ini) 
      ELSE
      i = indx_corr_id(p_cur,m_t + j_t+1,
     & j_t+1,chann_ini)
      ENDIF
      j_t =  j_max_ind(chann_ini)
      m_t = j_t		  
      IF(.not.identical_particles_defined) THEN
      j = indx_corr(m_t + j_t+1,
     & j_t+1,chann_ini) 
      ELSE
      j = indx_corr_id(p_cur,m_t + j_t+1,
     & j_t+1,chann_ini) 	  
      ENDIF	
	  if(identical_particles_defined) then
	  if(identical_particles_defined .and. p_lim_max.gt.p_cur) then
	  bk_num_ini_states = indx_corr_id(p_lim_max,m_t + j_t+1,
     & j_t+1,chann_ini)
	  else
	  bk_num_ini_states = indx_corr_id(p_cur,m_t + j_t+1,
     & j_t+1,chann_ini)
	  end if
	  else
	  bk_num_ini_states = indx_corr(m_t + j_t+1,
     & j_t+1,chann_ini)
	 end if
      num_ini_states = bk_num_ini_states - i	+ 1
!	  write(*,'(a,4(1x,i0))')'num_ini_states',myid,i,bk_num_ini_states
!     & ,num_ini_states
      ini_st_check_file = i
      IF(s_ini.ne.ini_st_check_file .or. num_ini_states.le.0) THEN
      IF(MYID.eq.0) PRINT*,
     & "ERROR IN INITIALIZATION OF INI_ARRAY SETUP"	
      IF(myid.eq.0) THEN
      PRINT*,"s_ini=",s_ini
      PRINT*,"ini_st_check_file=",ini_st_check_file
      PRINT*,"num_ini_states=",num_ini_states	 	  
      ENDIF	  
      STOP
      ENDIF	
      END SUBROUTINE DEFINE_WIDTH
      SUBROUTINE PRINT_ORBITING_TRAJECTORY
! This subroutine is updated by Bikramaditya Mandal
      USE INTEGRATOR
      USE MPI
      USE MPI_DATA
      IMPLICIT NONE	
      INTEGER index_traj,index_proc,numb_proc,i
      LOGICAL :: CREATE_FILE = .FALSE.	  
      REAL*8	R_fin, l_real,b_imparameter,phi_fin,surv_prob,
     & err_energ_rec,err_probab_rec,bk_or,bk_os,bk_sw	  
      task_size = 10*n_traject_alloc
!      DO i=1,n_traject_alloc
!      IF(ORBIT_TRAJ_FLAG(i)) THEN
!      PRINT*, myid, i
!      PRINT*,ORBIT_TRAJ_DATA(:,i)	   
!!      ENDIF	  
!      ENDDO	  
      CALL MPI_GATHER(ORBIT_TRAJ_DATA,task_size,MPI_DOUBLE_PRECISION,
     & ORBIT_TRAJ_DATA_ALL,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)	
      task_size = n_traject_alloc
      CALL MPI_GATHER(ORBIT_TRAJ_FLAG,task_size,MPI_LOGICAL,
     & ORBIT_TRAJ_FLAG_ALL,
     &	  task_size,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr_mpi)
!      ORBIT_TRAJ_DATA = 0d0
      IF(myid.eq.0) THEN     	  
      DO index_proc = 0,nproc/mpi_traject-1
      DO index_traj = 1,n_traject_alloc
!      ORBIT_TRAJ_FLAG(index_traj) = .FALSE.	  
      IF(ORBIT_TRAJ_FLAG_ALL(index_traj,
     & index_proc*mpi_traject+1))
     &	CREATE_FILE = .TRUE.
      ENDDO
      ENDDO	
      ENDIF	  

! Bikram Start:	  
	  bk_sw = 0.00d0
	  IF(myid.eq.0) THEN     	  
      DO index_proc = 0,nproc/mpi_traject-1
      DO index_traj = 1,n_traject_alloc
	  numb_proc = 	index_proc*mpi_traject+1
      bk_sw=bk_sw+ORBIT_TRAJ_DATA_ALL(10,index_traj,numb_proc)
      ENDDO
      ENDDO	
      ENDIF
	  
	  if(bk_sw.ge.1.00d0) then
      IF(CREATE_FILE.and. myid.eq.0) THEN!! IN CASE OF ORBITING
      OPEN(2345,FILE="RESONANCE_TRAJECT.out",POSiTION = "APPEND")

!      WRITE(2345,'(a35)')"ORBITING TRAJECTORIES IF THEY EXIST FOR"
      WRITE(2345,"(a9,1x,i3,2x,a10,1x,i6)")
     & "i_energy=",i_ener,"state_ini=",s_st	  
      WRITE(2345,
     & '(a9,1x,a12,2x,a8,2x,a9,2x,a11,2x,a10,2x,a11,2x,a5,2x,a5)')
     & "L_ORBITAL","B_IMPACT","R_FINAL",
     &  "PHI_FINAL", "SURV_PROBAB", "%ERR_ENERG",
     &  "%ERR_PROBAB","LOOPS","OSCIL"	 
      DO index_proc = 0,nproc/mpi_traject-1
      DO index_traj = 1,n_traject_alloc
      ORBIT_TRAJ_FLAG(index_traj) = .FALSE.	  
      IF(ORBIT_TRAJ_FLAG_ALL(index_traj,
     & index_proc*mpi_traject+1)) THEN
      numb_proc = 	index_proc*mpi_traject+1 
      b_imparameter  =	ORBIT_TRAJ_DATA_ALL(1,index_traj,numb_proc) 
      l_real  =	ORBIT_TRAJ_DATA_ALL(2,index_traj,numb_proc)
      R_fin  =  ORBIT_TRAJ_DATA_ALL(3,index_traj,numb_proc)
      phi_fin  =  ORBIT_TRAJ_DATA_ALL(4,index_traj,numb_proc)
      surv_prob  =  ORBIT_TRAJ_DATA_ALL(5,index_traj,numb_proc)	
      err_energ_rec  =  ORBIT_TRAJ_DATA_ALL(6,index_traj,numb_proc)	
      err_probab_rec  =  ORBIT_TRAJ_DATA_ALL(7,index_traj,numb_proc)
      bk_or  =  ORBIT_TRAJ_DATA_ALL(8,index_traj,numb_proc)
      bk_os  =  ORBIT_TRAJ_DATA_ALL(9,index_traj,numb_proc)
!      WRITE(2345,*)	 index_traj,numb_proc
      if(ORBIT_TRAJ_DATA_ALL(10,index_traj,numb_proc).eq.1.00d0) then 
      WRITE(2345,
     & '(i9,1x,f12.6,2x,f8.3,2x,f9.3,2x,f11.8,2x,f10.7,2x,f11.8
     & ,2x,i5,2x,i5)')
     &  int(l_real),b_imparameter,R_fin,
     & phi_fin,surv_prob,err_energ_rec,err_probab_rec
     & ,int(bk_or),int(bk_os)
! Bikram Oct'18 Start:
	  open(1991,file='bikram_mqct.tmp',POSiTION='append')
	  write(1991,*)100
	  close(1991)
! Bikram End
      endif	 
	 
      ENDIF

      ENDDO
      ENDDO	  
	 
      write(2345,*)
      CLOSE(2345)

      ENDIF
      else
      endif
! Bikram End.
	  
	       ORBIT_TRAJ_DATA = 0d0	
      DO i=1,n_traject_alloc
      IF(ORBIT_TRAJ_FLAG(i)) THEN
      ORBIT_TRAJ_FLAG(i) = .FALSE.
      ENDIF	  
      ENDDO	
      END SUBROUTINE PRINT_ORBITING_TRAJECTORY

      SUBROUTINE PRINT_ERRORS
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE CONSTANTS
      USE CHECK_FILE	  
      USE INTEGRATOR
      USE MPI_DATA
      USE MPI
      USE ERRORS	  
      USE MONTE_CARLO_STORAGE	  
      IMPLICIT NONE      
      INTEGER i,j,k,i_e
      REAL*8 energy_error_maximum, probability_error_maximum	  
	  
	  CALL DEFINE_MAX_NUM_TRAJ(n_traject_alloc)
	  IF(MYID.eq.0) THEN
      err_traj_max_prb = 0
      err_proc_max_prb = 0
      err_state_max_prb = 0	
      err_traj_max_ener = 0 
      err_proc_max_ener = 0!-1	 
      err_state_max_ener = 0!	  
      ENDIF 	  
	  task_size = 7*n_traject_alloc*nmbr_of_enrgs*num_ini_states
      CALL MPI_GATHER(TRAJECT_END_DATA,task_size,MPI_DOUBLE_PRECISION,
     & TRAJECT_END_DATA_ALL,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)	
	  IF(MYID.EQ.0) THEN
      OPEN(4,FILE="INTEGRATOR_ERRORS.out",POSiTION="APPEND")	  
      DO i_e=i_curr,i_ener-1
      probability_error_maximum = 0d0
      energy_error_maximum = 0d0	
      DO k= 1, nproc	  
      DO i = 1,num_ini_states
      DO j = 1,n_traject_alloc
	  IF(probability_error_maximum.lt.
     & TRAJECT_END_DATA_ALL(7,j,i,i_e,k)) THEN
      probability_error_maximum =
     & TRAJECT_END_DATA_ALL(7,j,i,i_e,k)
      err_traj_max_prb(i_e) = j 
      err_proc_max_prb(i_e) = k!-1	 
      err_state_max_prb(i_e) = i 
      ENDIF	 
	 
	  IF(energy_error_maximum.lt.
     & TRAJECT_END_DATA_ALL(6,j,i,i_e,k)) THEN
      energy_error_maximum =
     & TRAJECT_END_DATA_ALL(6,j,i,i_e,k)
      err_traj_max_ener(i_e) = j 
      err_proc_max_ener(i_e) = k!-1	 
      err_state_max_ener(i_e) = i!	 
      ENDIF	 
	 
      ENDDO 	  
      ENDDO
      ENDDO
      WRITE(4,*) "FOR ENERGY=",i_e
      WRITE(4,*) "MAXIMUM PROBABILITY ERROR"		  
      WRITE(4,*)"B_IMP",
     & TRAJECT_END_DATA_ALL(1,err_traj_max_prb(i_e),
     & err_state_max_prb(i_e),i_e,err_proc_max_prb(i_e))
      WRITE(4,*)"L_ORB",
     & TRAJECT_END_DATA_ALL(2,err_traj_max_prb(i_e),
     & err_state_max_prb(i_e),i_e,err_proc_max_prb(i_e))	 
      WRITE(4,*)"R_FIN",
     & TRAJECT_END_DATA_ALL(3,err_traj_max_prb(i_e),
     & err_state_max_prb(i_e),i_e,err_proc_max_prb(i_e))	
      WRITE(4,*)"PHI_FIN",
     & TRAJECT_END_DATA_ALL(4,err_traj_max_prb(i_e),
     & err_state_max_prb(i_e),i_e,err_proc_max_prb(i_e))
      WRITE(4,*)"PROBAB_SURV",
     & TRAJECT_END_DATA_ALL(5,err_traj_max_prb(i_e),
     & err_state_max_prb(i_e),i_e,err_proc_max_prb(i_e))
      WRITE(4,*)"ERROR_ENERG",
     & TRAJECT_END_DATA_ALL(6,err_traj_max_prb(i_e),
     & err_state_max_prb(i_e),i_e,err_proc_max_prb(i_e))
      WRITE(4,*)"ERROR_PROBAB",
     & TRAJECT_END_DATA_ALL(7,err_traj_max_prb(i_e),
     & err_state_max_prb(i_e),i_e,err_proc_max_prb(i_e))
      WRITE(4,*)"STATE",err_state_max_prb(i_e)-1+ini_st_check_file
      WRITE(4,*)"#NTRAJECT",err_traj_max_prb(i_e)
      WRITE(4,*)"#NPROC",err_proc_max_prb(i_e)-1
      WRITE(4,*)		 
      WRITE(4,*) "MAXIMUM ENERGY ERROR"		  
      WRITE(4,*)"B_IMP",
     & TRAJECT_END_DATA_ALL(1,err_traj_max_ener(i_e),
     & err_state_max_ener(i_e),i_e,err_proc_max_ener(i_e))
      WRITE(4,*)"L_ORB",
     & TRAJECT_END_DATA_ALL(2,err_traj_max_ener(i_e),
     & err_state_max_ener(i_e),i_e,err_proc_max_ener(i_e))	 
      WRITE(4,*)"R_FIN",
     & TRAJECT_END_DATA_ALL(3,err_traj_max_ener(i_e),
     & err_state_max_ener(i_e),i_e,err_proc_max_ener(i_e))	
      WRITE(4,*)"PHI_FIN",
     & TRAJECT_END_DATA_ALL(4,err_traj_max_ener(i_e),
     & err_state_max_ener(i_e),i_e,err_proc_max_ener(i_e))
      WRITE(4,*)"PROBAB_SURV",
     & TRAJECT_END_DATA_ALL(5,err_traj_max_ener(i_e),
     & err_state_max_ener(i_e),i_e,err_proc_max_ener(i_e))
      WRITE(4,*)"ERROR_ENERG",
     & TRAJECT_END_DATA_ALL(6,err_traj_max_ener(i_e),
     & err_state_max_ener(i_e),i_e,err_proc_max_ener(i_e))
      WRITE(4,*)"ERROR_PROBAB",
     & TRAJECT_END_DATA_ALL(7,err_traj_max_ener(i_e),
     & err_state_max_ener(i_e),i_e,err_proc_max_ener(i_e))
      WRITE(4,*)"STATE",err_state_max_ener(i_e)-1+ini_st_check_file
      WRITE(4,*)"#NTRAJECT",err_traj_max_ener(i_e)
      WRITE(4,*)"#NPROC",err_proc_max_ener(i_e)-1	  
      WRITE(4,*) "-------------------------"	  
      ENDDO	  
      CLOSE(4)	  
      ENDIF	
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  	  
      END SUBROUTINE PRINT_ERRORS	  

      SUBROUTINE INELAST_CALC_FINE(sigma_el,nchann)
! This subroutine is updated by Bikramaditya Mandal
      USE INTEGRATOR
      USE MPI_DATA
      USE MPI_TASK_TRAJECT	  
      IMPLICIT NONE
      LOGICAL EVEN_NUM	  
      LOGICAL, ALLOCATABLE :: nproc_actual(:)	  
      INTEGER nchann,tot_num_of_traj_actual  
      INTEGER i,i_ip,i_traj,k,j_chann,traj_works,j_int_traj,
     & l_int_traj,id_traject 
      REAL*8 sigma_el(nchann)
      REAL*8, PARAMETER :: a_bohr = 5.2917721067d-1
      REAL*8, PARAMETER :: conv_E_cm_Ev = 27.21113845d0*8065.5448d0	  
      REAL*8, ALLOCATABLE :: j_def(:),opac_chann_all(:,:)
      REAL*8 posit_imp,negat_imp,part_cross	  
      tot_num_of_traj_actual = J_UP_INT-J_DOWN_INT+1	  
      ALLOCATE(j_def(tot_num_of_traj_actual),
     & opac_chann_all(nchann,tot_num_of_traj_actual))
      opac_chann_all = 0d0	 
      sigma_el = 0d0
      DO k=1,nchann	  
      DO l_int_traj = L_MIN_TRAJECT,L_MAX_TRAJECT
      traj_works = l_int_traj - L_MIN_TRAJECT + 1	
      IF(myid.eq.0) THEN	  
!      WRITE(*,*) 'j_int_traj',j_int_traj
!      WRITE(*,*) 'l_int_traj',l_int_traj
      ENDIF
!      GOTO 3456
      CALL TRAJ_ORB(l_int_traj,
     & id_traject,
     & 1,
     & dl_step_integ,
     & ident_skip,
     & L_MIN_TRAJECT)
!3456      id_traject = (l_int_traj-L_MIN_TRAJECT)/dl_step_integ+1
      i_traj = itraj_myid_l(1,id_traject)	  
      i_ip = itraj_myid_l(2,id_traject)
!      GOTO	5782  
      IF(ident_skip) THEN
      IF(parity_st.eq.1) THEN 
      posit_imp = exchange_ident
      negat_imp = 1d0 - exchange_ident	  
      ELSE
      posit_imp = 1d0 - exchange_ident	
      negat_imp = exchange_ident	  
      ENDIF	  
      ENDIF	  
!      WRITE(*,*) 'i_ip',i_ip
!      WRITE(*,*) 'i_traj',i_traj	  
      id_proc = i_ip*mpi_traject
!      WRITE(*,*) 'id_proc',id_proc
      part_cross = probab_J_all(k,i_traj,id_proc+1)
     & *(2d0*l_int_traj+1d0)/k_vec**2*pi
!      GOTO 2356
      IF(ident_skip) THEN
      IF(EVEN_NUM(l_int_traj)) THEN 
      part_cross = 	part_cross*posit_imp
      ELSE
      part_cross = 	part_cross*negat_imp	  
      ENDIF	  
      ENDIF
	 
2356      sigma_el(k) = sigma_el(k) +  part_cross     

      opac_chann_all(k,traj_works) = opac_chann_all(k,traj_works)
     &	+ part_cross

      j_def(traj_works) = l_int_traj
      ENDDO	  
      ENDDO

	  

      IF(myid.eq.0) THEN
      OPEN(2,FILE="OPACITY_FUNCTIONS_J.out",POSiTION="APPEND")
      WRITE(2,*) "STATE= ", s_st
      WRITE(2,'(a4,1x,f4.1)') "j12= ",j_curr_f	  
      WRITE(2,'(a4,1x,f4.1)') "m12= ",m_curr_f	  
      WRITE(2,'(a8,1x,i3,a,1x,f8.2,1x,a5)')
     &"#ENERGY=",i_ener,", U(i)= ",E_coll(i_ener)*conv_E_cm_Ev," cm^-1"
      DO j_chann = 1,nchann	  
      IF(j_chann.eq.1) WRITE(2,'(a18,1x,a18,1x)',ADVANCE="NO")
     & "B_impact_parameter","J_TOTAL_momentum"      	 
      WRITE(2,'(a8,1x,i3,2x)',ADVANCE="NO") "CHANNEL=",j_chann
      ENDDO
      WRITE(2,*)	  
      DO i=1,tot_num_of_traj_actual
      DO j_chann = 1,nchann
      IF(j_chann.eq.1)
     & WRITE(2,'(e18.11,3x,i9,6x)',ADVANCE="NO")
     &	 dsqrt((j_def(i)+1)*j_def(i))/k_vec,INT(j_def(i))
      WRITE(2,'(e12.5,1x)',ADVANCE="NO") opac_chann_all(j_chann,i)
     & *a_bohr**2
      ENDDO
      WRITE(2,*)  
      ENDDO  
      CLOSE(2)
      OPEN(2,FILE="M_PROJ_CROSSECTIONS.out",POSiTION="APPEND")
      WRITE(2,*) "M RESOLVED CROSSECTIONS WITHOUT BILLING CORRECTION"	  
      WRITE(2,"(a10,1x,i4)") "INI_STATE#", s_st
      WRITE(2,'(a10,1x,i4)') "ENERGY#",i_ener
      WRITE(2,'(a4,1x,f4.1)') "j12=",j_curr_f	  
      WRITE(2,'(a4,1x,f4.1)') "m12=",m_curr_f		  
      WRITE(2,'(a8,1x,a19)') "CHANNEL#", "CROSSSECTIONS,ANG^2"
      DO k=1,nchann
      WRITE(2,'(i4,6x,e18.11)') k,sigma_el(k)*a_bohr**2
      ENDDO
      FLUSH(2)	  
      CLOSE(2)	  
      ENDIF	  
      opacity
     & (:,1:tot_num_of_traj_actual,i_ener) = 
     & opac_chann_all(:,:) 	 
      END SUBROUTINE INELAST_CALC_FINE
      SUBROUTINE CS_MATRIX_IDENT
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE CONSTANTS
      USE INTEGRATOR
      USE MPI_DATA
      USE MPI
      USE MPI_TASK_TRAJECT		  
      USE CS_MATRIX
      IMPLICIT NONE 
      LOGICAL BELONGS	  
      INTEGER i,j,k,ind_cs,jnd_cs,origin
      EXTERNAL BELONGS	  
      states_size_cs = 0
      DO i=1,states_size
      IF(m12(i).eq.m12(s_st))	states_size_cs = states_size_cs + 1  
      ENDDO	  
      IF(number_of_channels.lt.states_size_cs)
     & STOP "ERROR: CS FAILS"
	  mat_ts_mpi = states_size_cs*(states_size_cs+1)/2	 
      IF(myid.eq.0) PRINT*, "states_size_cs",states_size_cs	 
      IF(ALLOCATED(sys_var_cs)) THEN
      DEALLOCATE(sys_var_cs,deriv_sys_cs,ind_mat_cs,
     & ind_state_cs,Mat_el_cs,Mat_el_cs_der,Ej_cs_im,
     & portion_of_csst_per_task, 
     & portion_of_csts_per_task,ind_cs_rule	 )
      ENDIF	  
      ALLOCATE(sys_var_cs(2*states_size_cs+8))
	  ALLOCATE(deriv_sys_cs(2*states_size_cs+8))
      ALLOCATE(ind_state_cs(states_size_cs))
      ALLOCATE(ind_mat_cs(states_size))
      ALLOCATE(ind_cs_rule(2,mat_ts_mpi))	  
      ALLOCATE(Mat_el_cs(n_r_coll,states_size_cs,states_size_cs))
      ALLOCATE(Mat_el_cs_der(n_r_coll,states_size_cs,states_size_cs))
      ALLOCATE(Ej_cs_im(states_size_cs))
      ALLOCATE(portion_of_csst_per_task(2,nproc))
      ALLOCATE(portion_of_csts_per_task(2,nproc))
      j = 0	  
      DO i=1,states_size
      IF(m12(i).eq.m12(s_st)) THEN
      j = j + 1	  
      ind_state_cs(j) = i
      ind_mat_cs(i) = j	  
      Ej_cs_im(j) = Ej(i)	  
      ENDIF	  
      ENDDO	 
      ind_cs_rule = 0
      ind_cs = 0	  
	  DO i=1,states_size_cs
      DO j=1,i
      ind_cs = ind_cs + 1 
      ind_cs_rule(1,ind_cs) = i
      ind_cs_rule(2,ind_cs) = j      	  
      ENDDO
      ENDDO	 
      IF(ind_cs.ne.mat_ts_mpi) STOP "ERROR:CS TOTAL SIZE WRONG"	  

	  


	  
      IF(mpi_task_defined) THEN
	  
	  
      DO k=1,total_size
      i  = ind_mat(1,k)
      j  = ind_mat(2,k)
      ind_cs = ind_mat_cs(i)
      jnd_cs = ind_mat_cs(j)	  
      IF(m12(i).eq.m12(s_st) .and. myid.eq.0) THEN
      Mat_el_cs_der(:,ind_cs,jnd_cs) = mat_buffer_der(:,k) !!!!
      Mat_el_cs_der(:,jnd_cs,ind_cs) = mat_buffer_der(:,k)	!!!!  PROBLEM HERE  	  
      Mat_el_cs(:,ind_cs,jnd_cs) = mat_buffer(:,k)  !!!!
      Mat_el_cs(:,jnd_cs,ind_cs) = mat_buffer(:,k)  !!!
      ENDIF	  
      ENDDO		  
	  

!!!! BROADCASTING MATRIX
      IF(MYID.eq.0) PRINT*,"DATA CS BROADCASTING"
      task_size	= states_size_cs**2*n_r_coll  
      CALL MPI_Bcast(Mat_el_cs,task_size,
     & MPI_REAL8, 0,MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_Bcast(Mat_el_cs_der,task_size,
     & MPI_REAL8, 0,MPI_COMM_WORLD,ierr_mpi)	 
      IF(MYID.eq.0) PRINT*,"DATA CS BROADCASTED"	  

	    
	  ELSE
      DO k=1,total_size
      i  = ind_mat(1,k)
      j  = ind_mat(2,k)
      ind_cs = ind_mat_cs(i)
      jnd_cs = ind_mat_cs(j)	  
      IF(m12(i).eq.m12(s_st)) THEN
      Mat_el_cs(:,ind_cs,jnd_cs) = Mat_el(:,k)  !!!!
      Mat_el_cs(:,jnd_cs,ind_cs) = Mat_el(:,k)  !!!
      Mat_el_cs_der(:,ind_cs,jnd_cs) = Mat_el_der(:,k) !!!!
      Mat_el_cs_der(:,jnd_cs,ind_cs) = Mat_el_der(:,k)	!!!!  PROBLEM HERE  
      ENDIF	  
      ENDDO	
      ENDIF	

      IF(myid.eq.0) PRINT*, "ARRAYS_CS ALLOCATED FOR",s_st
!      STOP	
	  state_size_check = 0  
      size_csst_chunk_mpi = states_size_cs/mpi_task_per_traject
      residue_csst_mpi = states_size_cs-
     & size_csst_chunk_mpi*mpi_task_per_traject
	  
      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_csst_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)
      IF(k_mpi_proc.gt.nproc .or. k_mpi_proc.le.0) THEN
      IF(MYID.eq.0) PRINT*, "ERROR_CS",k,k_p,mpi_traject_roots(k)
      STOP	  
      ENDIF 	  
      portion_of_csst_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_csst_chunk_mpi+1)
      portion_of_csst_per_task(2,k_mpi_proc) =
     & (k_p)*(size_csst_chunk_mpi+1)
      ENDDO
      state_size_check = state_size_check + size_csst_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 
      IF(k_mpi_proc.gt.nproc .or. k_mpi_proc.le.0) THEN
      IF(MYID.eq.0) PRINT*, "ERROR_CS",k,k_p,mpi_traject_roots(k)
      STOP	  
      ENDIF  	  
      portion_of_csst_per_task(1,k_mpi_proc)
     & = residue_csst_mpi*(size_csst_chunk_mpi+1)+1  
     & + (k_p-1-residue_csst_mpi)*size_csst_chunk_mpi
      portion_of_csst_per_task(2,k_mpi_proc)
     & = residue_csst_mpi*(size_csst_chunk_mpi+1)+	  
     & (k_p-residue_csst_mpi)*size_csst_chunk_mpi
      ENDDO	 
      state_size_check = state_size_check + size_csst_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_csst_per_task(2,k_mpi_proc).ne.states_size_cs) 
     & STOP "WRONG CS ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO	
	  
	  state_size_check = 0
      size_csts_chunk_mpi = mat_ts_mpi
     & /mpi_task_per_traject
      residue_csts_mpi = mat_ts_mpi-
     & size_csts_chunk_mpi*mpi_task_per_traject
	  
      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_csts_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)
      IF(k_mpi_proc.gt.nproc .or. k_mpi_proc.le.0) THEN
      IF(MYID.eq.0) PRINT*, "ERROR_CS",k,k_p,mpi_traject_roots(k)
      STOP	  
      ENDIF 	  
      portion_of_csts_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_csts_chunk_mpi+1)
      portion_of_csts_per_task(2,k_mpi_proc) =
     & (k_p)*(size_csts_chunk_mpi+1)
      ENDDO
      state_size_check = state_size_check + size_csts_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 
      IF(k_mpi_proc.gt.nproc .or. k_mpi_proc.le.0) THEN
      IF(MYID.eq.0) PRINT*, "ERROR_CS",k,k_p,mpi_traject_roots(k)
      STOP	  
      ENDIF  	  
      portion_of_csts_per_task(1,k_mpi_proc)
     & = residue_csts_mpi*(size_csts_chunk_mpi+1)+1  
     & + (k_p-1-residue_csts_mpi)*size_csts_chunk_mpi
      portion_of_csts_per_task(2,k_mpi_proc)
     & = residue_csts_mpi*(size_csts_chunk_mpi+1)+	  
     & (k_p-residue_csts_mpi)*size_csts_chunk_mpi
      ENDDO	 
      state_size_check = state_size_check + size_csts_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_csts_per_task(2,k_mpi_proc).ne.mat_ts_mpi) 
     & STOP "WRONG CS ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO	  

	  IF(MYID.eq.0) PRINT*, "CS ARRAYES MADE"
 
      END SUBROUTINE CS_MATRIX_IDENT
	  	  
      SUBROUTINE potential_data_cs(R,i,j,v_m,der_v_m)
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE GRID_INTEGR
      USE CONSTANTS
      USE MPI_TASK_TRAJECT
      USE MPI_DATA
      USE INTEGRATOR
      USE CS_MATRIX	  
      USE OLD_Mij	  
      IMPLICIT NONE	 
      REAL*8 x,v_m,der_v_m,R,V_COULPING_TERMS
      INTEGER k,k_real,i_exp_term,i,j
      x = R
      v_m = 0d0	  
      der_v_m = 0d0
      IF(x.gt.R_COM(n_r_coll)) RETURN
      CALL splint(R_COM,Mat_el_cs(:,i,j),Mat_el_cs_der(:,i,j)
     & ,n_r_coll,x,v_m,der_v_m) 
      v_m = v_m - Mat_el_cs(n_r_coll,i,j)
!	  v_m = 0d0
!	  der_v_m = 0d0
      END SUBROUTINE potential_data_cs

      FUNCTION Mjmr_cs(R,i,j)
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE GRID_INTEGR
      USE CONSTANTS
      USE MPI_TASK_TRAJECT
      USE MPI_DATA
      USE INTEGRATOR
      USE CS_MATRIX	  
      USE OLD_Mij	  
      IMPLICIT NONE	 
      REAL*8 x,v_m,der_v_m,R,V_COULPING_TERMS,Mjmr_cs
      INTEGER k,k_real,i_exp_term,i,j
      x = R
      Mjmr_cs = 0d0	  
      der_v_m = 0d0
      IF(x.gt.R_COM(n_r_coll)) RETURN
      CALL splint(R_COM,Mat_el_cs(:,i,j),Mat_el_cs_der(:,i,j)
     & ,n_r_coll,x,Mjmr_cs,der_v_m) 
      Mjmr_cs = Mjmr_cs - Mat_el_cs(n_r_coll,i,j)
!	  Mjmr_cs = 0d0
      END FUNCTION Mjmr_cs	  

!--------------------------------------------------------------------
! Bikram July 09, 2023: Making this change
! this subroutine is to find the new Monte-Carlo error values
! this code is written by Carolin Joy and now imcorporated within MQCT by me
!--------------------------------------------------------------------
      subroutine MCerror(num_rows, num_columns)
      implicit none
      integer ::  ii, k, irow, icol, jj, num_rows, num_columns
      real, dimension(:,:), allocatable :: data_1, sum, CumSum, 
     & mean_val, sig_minus, A, err
      
      open(1, file='MC_XSections.out', status='old', action='read')
	  open(20, file = "MC_Error.out", status = 'new', action = 'write')
      allocate(data_1(num_rows,num_columns))
      do irow = 1, num_rows
        read(1,*) (data_1(irow,icol), icol = 1, num_columns) ! array where opacities are stored
      end do
      close(1)
	  
      allocate(sum(num_rows,num_columns))
      allocate(CumSum(num_rows,num_columns))
      allocate(mean_val(num_rows,num_columns))
      allocate(sig_minus(num_rows,num_columns))
      allocate(A(num_rows,num_columns))
      allocate(err(num_rows,num_columns))
	  
      do jj = 1, num_columns
        sum = 0.d0
        do ii = 1, num_rows 
          sum =  sum + data_1(ii,jj) 
          CumSum(ii,jj) =sum(ii,jj)
          mean_val(ii,jj) = CumSum(ii,jj)/ii
          ! write(12,*) CumSum(ii,jj), mean_val(ii,jj)    ! array where cumulative sum & mean values are stored.
        end do
      end do     
        
      do jj = 1, num_columns
        A = 0.d0
        ! No of rows loop
        do ii = 2, num_rows
          do k = 1, ii
            sig_minus(k,jj) = (CumSum(ii,jj) - data_1(k,jj))/(ii-1)  ! calculating sig_minus
            A(ii,jj) = A(ii,jj) + abs(mean_val(ii,jj) - sig_minus(k,jj))**2
          end do
          A(ii,jj) = A(ii,jj) / ii
          A(ii,jj) = (sqrt(A(ii,jj))/ mean_val(ii,jj)) *100
          err(ii,jj) = A(ii,jj)
        end do
        ! write(20,*)
      end do
         
!             write(20,'(f12.1)', Advance ='no') A(ii,jj)
      do ii = 2, num_rows
        do jj = 1, num_columns
            write(20,'(f12.6)', Advance ='no') err(ii,jj)
        end do
        write(20,*)
      end do
      close(20)
      return
      end subroutine

! Bikram End.
