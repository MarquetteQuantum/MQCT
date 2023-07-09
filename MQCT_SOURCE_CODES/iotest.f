      MODULE ERRORS
!!!! TAGS FOR DIFFERENRT TYPES OF ERROR.   
      LOGICAL :: error_file_reading = .FALSE.
      LOGICAL :: error_compute_exe = .FALSE.
      LOGICAL :: error_int_num_read = .FALSE.
      LOGICAL :: error_real_num_read = .FALSE.
      LOGICAL :: error_re_e_num_read = .FALSE.	  
      INTEGER :: MPIS_ERROR = 1
      INTEGER :: MPIS_ERROR_TOT	  
      CHARACTER(LEN=150)show_error  
      END MODULE ERRORS
      MODULE CONSTANTS
      CHARACTER(LEN=22)::USR_INP_LEVEL = "USER_DEFINED_BASIS.DAT" !CONTAINS CONSTANTS AND NAMES
      CHARACTER(LEN=11)	:: MATRIX_NAME_MIJ_UF = "MTRX_UF.DAT"
      CHARACTER(LEN=8)	:: MATRIX_NAME_MIJ = "MTRX.DAT"
      CHARACTER(LEN=17)  :: VIB_MIX_STATE = "VIB_MIX_STATE.DAT"	  
      REAL*8 autoeV, amutoau, eVtown,hartree_time_ps,picosecond,autokelv
      REAL*8 kcal
      REAL*8,PARAMETER ::  a_bohr = 0.52917721067d0 	  
      PARAMETER (autoeV=27.21113845d0,kcal = 627.509469d0)
      PARAMETER(autokelv= 1.1604519d4)
      PARAMETER (amutoau=1822.8885d0, eVtown=8065.5448d0
     &  ,hartree_time_ps=2.418884326505d-5,picosecond=1/hartree_time_ps)
      REAL*8 r_eq_diatom
      REAL*8 mass_vib_diatom,mass_vib_diatom1,mass_vib_diatom2	  
      REAL*8 massAB_M
      REAL*8 conv_unit_e,conv_unit_r
!      PARAMETER(MIJ_ZERO_CUT=1d-12)	  
      CHARACTER(LEN=12) :: potential_file_name = "VGRID_UF.DAT"	  
      END MODULE CONSTANTS
      MODULE POT_STORE
!!!! MODULE CONTAINS PES and W-FUNCTIONS ARRAYS FOR GAUSS_LEGENDRE INTEGRATION.	  
      REAL*8, ALLOCATABLE :: V_3_3(:,:,:,:,:,:)
      REAL*8, ALLOCATABLE :: V_3_3_int_buffer(:,:,:,:,:)
      REAL*8, ALLOCATABLE :: V_2_2(:,:,:,:)
      REAL*8, ALLOCATABLE :: V_3_1(:,:,:)
      REAL*8, ALLOCATABLE :: V_2_1(:,:)
      REAL*8, ALLOCATABLE :: V_3_2(:,:,:,:,:)
      REAL*8, ALLOCATABLE :: V_3_2_int_buffer(:,:,:,:)	  
      REAL*8, ALLOCATABLE :: V_VIB_2(:,:,:)
      REAL*8, ALLOCATABLE :: V_VIB_2_2(:,:,:,:,:,:)
      REAL*8, ALLOCATABLE :: V_VIB_2_2_int_buffer(:,:,:,:,:)		  
      REAL*8, ALLOCATABLE :: wf_vib_part_diatom(:,:)
      REAL*8, ALLOCATABLE :: wf_vib_part_diatom1(:,:)
      REAL*8, ALLOCATABLE :: wf_vib_part_diatom2(:,:)
      END MODULE POT_STORE
      MODULE FACTORIAL
      REAL*8, ALLOCATABLE :: LOGFACT(:)	  
      END MODULE FACTORIAL
      MODULE EIGEN_VECTORS
!!! STORE BASIS FOR ASYMETRIC TOP	  
      REAL*8, ALLOCATABLE :: M1_VECTORS(:,:)
      REAL*8, ALLOCATABLE :: M2_VECTORS(:,:)
      REAL*8, ALLOCATABLE :: M_VECTORS(:,:)	
      LOGICAL :: eig_vec_def = .FALSE.
      REAL*8 A_I, B_I, C_I
      REAL*8 A2_I, B2_I, C2_I
      REAL*8 A1_I, B1_I, C1_I 	  
      END MODULE EIGEN_VECTORS
      MODULE GRID_INTEGR
c! CONTAINS GRIDs FOR GAUSS LEGENDRE INTEGRATION	  
      REAL*8 xam,xbm,xgm 
      REAL*8, ALLOCATABLE :: wa(:),wb(:),wg(:)
      REAL*8, ALLOCATABLE :: xa(:),xb(:),xg(:)	  
      REAL*8 xal,xbl,xgl,alpha,beta,gamma
      REAL*8 xaam,xbbm,xggm
      REAL*8, ALLOCATABLE :: xaa(:),xbb(:),xgg(:)
      REAL*8, ALLOCATABLE :: waa(:),wbb(:),wgg(:)
      REAL*8 xaal,xbbl,xggl,aalpha,bbeta,ggamma
      REAL*8, PARAMETER :: pi = dacos(-1d0)
      REAL*8, ALLOCATABLE :: expansion_terms(:,:)
      REAL*8, ALLOCATABLE :: expansion_terms_der(:,:)	  
      REAL*8, ALLOCATABLE :: expansion_terms_im(:,:)	  
      REAL*8, ALLOCATABLE :: weig_terms(:)
      REAL*8, ALLOCATABLE :: vib_overlap(:,:)	  
      REAL*8, ALLOCATABLE :: r_vb_dt_inegrator(:)
      REAL*8, ALLOCATABLE :: r_vb_dt_integrator1(:)
      REAL*8, ALLOCATABLE :: r_vb_dt_integrator2(:)
      REAL*8, ALLOCATABLE :: EXP_POT_ELEMENT(:),
     & EXP_FORCE_ELEMENT(:)
      INTEGER, ALLOCATABLE :: A_TOP(:,:)
      INTEGER, ALLOCATABLE :: index_term_in_file(:,:)	  
      INTEGER N_TERMS_EXP
      LOGICAL :: EXPANSION_GRID_DEFINED = .FALSE.
      LOGICAL :: GAUSS_LEGENDRE_GRID_DEFINED = .FALSE.
      LOGICAL :: EXPANSION_MATRIX_DEFINED = .FALSE.
      LOGICAL :: TERM_READ_STAT = .FALSE.	  
      END MODULE GRID_INTEGR
      MODULE FINE_STRUCT_LEVELS
      IMPLICIT NONE
      REAL*8, ALLOCATABLE :: b_fine_coeff(:,:)
      INTEGER f_p_t,f_pp_t,s_p_t,s_pp_t
      INTEGER n_p_t,n_pp_t,m_s_t
      REAL*8 b_p_t,b_pp_t
      REAL*8 CG_p_s,CG_pp_s
      END MODULE FINE_STRUCT_LEVELS	  
      MODULE COMPUTE_EXPANSION_VARIABLES
      IMPLICIT	NONE  
      LOGICAL SIMPLIFICATION_EXP_MAT	  !CONTAINS LOCAL VARIABLES FOR COMPUTATION OF POTENTIAL EXPANSION
      INTEGER stp,stpp,ind_t
      INTEGER j12_t,m_t,j_p_t,j_pp_t,v1_p,v2_p,v1_pp,v2_pp	  
      INTEGER j12_p_t,j12_pp_t,j1_p_t,j1_pp_t,j2_p_t,j2_pp_t,
     & p_1,p_2,channp,channpp,m12_t 
      INTEGER ka1_p_t,ka2_p_t,kc1_p_t,kc2_p_t
      INTEGER ka1_pp_t,ka2_pp_t,kc1_pp_t,kc2_pp_t
      INTEGER ka_p_t,ka_pp_t,kc_p_t,kc_pp_t	  
      INTEGER mp,mpp,k1p,k1pp,k2p,k2pp,m_exp
      INTEGER l1_t,l2_t,l_t,nju1_t,nju2_t
      INTEGER nju_t,k_p_t,k_pp_t,eps_p_t,eps_pp_t
      INTEGER j_h_i_p,j_h_i_pp,m_h_i
      INTEGER sp_p,sp_pp
      INTEGER lambda_p_h,lambda_pp_h	  
      REAL*8  sp_p_t,sp_pp_t
      REAL*8 j_h_pp_t,j_h_p_t,m_h_t
      REAL*8 w_h_p,w_h_pp	  
      REAL*8 coeff_1_p,coeff_2_p,coeff_1_pp,coeff_2_pp
      REAL*8 coeff_p,coeff_pp	  
      REAL*8 M_coulp_ident,M_coulp_non_ident
      REAL*8 M_coulp_ident_1,M_coulp_ident_2,M_coulp_ident_3
      REAL*8 M_coulp_simplified, buff  
      REAL*8 switcher
      REAL*8 pot_coff_corr,exp_coeff_int,matrix_exp_coefficent	  
      REAL*8 CG_j1_j2_p,CG_j1_j2_pp,CG_j1_l1,CG_j2_l2,CG_l1_l2,
     & CG_j1_k1, CG_j2_k2
      INTEGER par_p, par_pp,symmetry_coeff,planr_coeff	  
      INTEGER term_limit_user
      END MODULE COMPUTE_EXPANSION_VARIABLES
	  
      MODULE OLD_MIJ
      INTEGER n_r_coll_old,states_size_old,total_size_old !! CONTAINS LOCAL VARIABLES FOR READING PREVIOUSLY SAVED Mij.dat
      INTEGER coll_type_old,number_of_channels_old,st_old
      INTEGER j_old_b,v_old_b,k_old_b,nr_cold,k_old,k1_old_b,
     & k2_old_b,f_old_b,n_old_b	        	  
      INTEGER j1_old_b,ka1_old_b,kc1_old_b,par_old_b,v1_old_b,v2_old_b
     & ,j2_old_b,ka2_old_b,kc2_old_b,i_old,ind_mat_old_1,ind_mat_old_2
      INTEGER eps1_old_b,eps2_old_b,eps_old_b,ka_old_b,kc_old_b
      INTEGER, ALLOCATABLE :: ind_mat_old(:,:), indx_chann_old(:),
     & j12_old(:),m12_old(:),j_ch_old(:),v_ch_old(:),eps_ch_old(:),
     & k_ch_old(:),eps1_ch_old(:),eps2_ch_old(:),ka_ch_old(:),
     & kc_ch_old(:),v1_ch_old(:),v2_ch_old(:),k1_ch_old(:),
     & k2_ch_old(:),
     & j1_ch_old(:),j2_ch_old(:),ka1_ch_old(:),kc1_ch_old(:)
     & ,ka2_ch_old(:),kc2_ch_old(:),parity_states_old(:)
      REAL*8, ALLOCATABLE :: j12_h_old(:),m12_h_old(:),j_h_ch_old(:)
      REAL*8 j_h_old_b
      INTEGER ir2,ir3,ir4,ir5,ir6,buffer_size_V	 
      REAL*8, ALLOCATABLE :: Mat_rest(:,:),Mat_rest_der(:,:)
      REAL*8, ALLOCATABLE :: Mat_el_non_zero(:,:),
     & Mat_el_non_zero_der(:,:)
      INTEGER, ALLOCATABLE :: ind_mat_non_zero(:,:)	 
      INTEGER k_non_zero
      LOGICAL, ALLOCATABLE	:: stts_to_excl(:)	  
      CHARACTER(LEN=12) :: buffer_word_1
      CHARACTER(LEN=15) :: buffer_word_2
      CHARACTER(LEN=17)	:: buffer_word_3
      LOGICAL :: CRITICAL_ERROR = .FALSE.
      REAL*8 TIME_MAT_START,TIME_MAT_FINISH
      REAL*8 deriv_bgn,deriv_end	  
      END MODULE OLD_MIJ
      MODULE MPI_TASK_TRAJECT
      INTEGER, ALLOCATABLE :: mpi_traject_roots(:)!!! MPI PARALLEZATION OF EACH TRAJECTORY
      INTEGER, ALLOCATABLE :: portion_of_MIJ_per_task(:,:)
      INTEGER, ALLOCATABLE :: portion_of_state_per_task(:,:)
      INTEGER, ALLOCATABLE :: portion_of_work_per_task(:,:)		  
      INTEGER traject_roots,size_mij_chunk_mpi,residue_mij_mpi
      INTEGER size_state_chunk_mpi,residue_state_mpi
      INTEGER size_work_chunk_mpi,residue_work_mpi	  
      INTEGER :: total_size_check = 0
      INTEGER :: state_size_check = 0
      INTEGER :: work_size_check = 0		  
      INTEGER total_size_mpi,task_portion_size,tag1,tag2,tag3,k_p
      INTEGER k_mpi_proc
      INTEGER GROUPS_IND, ID_IN_GROUP_INDEX	  
      INTEGER, ALLOCATABLE :: mpi_root_belongs(:)	  
      REAL*8, ALLOCATABLE :: buffer_mpi_portion(:,:)
      REAL*8, ALLOCATABLE  :: dq_dt_mpi(:,:)
      INTEGER, ALLOCATABLE :: groups(:),comms(:)
      INTEGER, ALLOCATABLE :: process_rank(:,:) 
      INTEGER, ALLOCATABLE :: groups_distr(:),comms_distr(:)
      INTEGER, ALLOCATABLE :: process_rank_distr(:,:)
      INTEGER wrld_group,id_proc_in_group	  
      END MODULE MPI_TASK_TRAJECT			  
	  
      MODULE EXPANSION_MAT_STORAGE
      REAL*8 , ALLOCATABLE :: TERM_MATRIX_ELEMENT(:,:,:,:)	  !!! SAVES EXPANSION MATRIX ELEMENTS 
      END MODULE EXPANSION_MAT_STORAGE	  
      MODULE VARIABLES 																			!!! GLOBAL MAIN VARIABLES
c! VARIABLES	  
      LOGICAL rk4_defined
      LOGICAL odeint_defined	  
      LOGICAL scatter_param_defined
      LOGICAL ener_man_defined
      LOGICAL ener_auto_defined
      LOGICAL prn_l_defined	  
      LOGICAL deflect_fun_defined
      LOGICAL orbit_traj_defined
      LOGICAL tm_lim_defined
      LOGICAL b_impact_defined
      LOGICAL check_point_defined
      LOGICAL unformat_defined	  
      LOGICAL atomic_masses_defined
      LOGICAL mpi_task_defined	  
      LOGICAL user_defined
      LOGICAL all_level_included
      LOGICAL energy_defined
      LOGICAL para_ortho_defined
      LOGICAL emax_defined
      LOGICAL rot_const_defined
      LOGICAL sys_type_defined
      LOGICAL channels_defined
      LOGICAL fine_structure_defined
      LOGICAL hyper_fine_stricture_defined	  
      LOGICAL number_of_channels_defined
      LOGICAL identical_particles_defined	  
      LOGICAL ini_chann_defined
      LOGICAL vib_const_defined
      LOGICAL mlc_mlc_chn_num_defined
      LOGICAL mlc_mlc_emax_defined	  
      LOGICAL morse_pot_defined
      LOGICAL calc_elast_defined
      LOGICAL diff_cross_defined	  
      LOGICAL monte_carlo_defined	  
      LOGICAL bikram_int																		!Bikram
      LOGICAL bikram_theta																		!Bikram
      LOGICAL bikram_print																		!Bikram
      LOGICAL bk_prob_interpolation																!Bikram
      LOGICAL bk_nrg_err																		!Bikram
      LOGICAL bk_step_size																		!Bikram
      LOGICAL bikram_mij_shift																	!Bikram
      LOGICAL bikram_mij_multiprint																!Bikram
      LOGICAL bikram_eq_grd_gamma																!Bikram
      LOGICAL bikram_identical_pes																!Bikram
      LOGICAL rms_defined																		!Bikram
      LOGICAL bikram_ident_terms																!Bikram
      LOGICAL bikram_rebalance, bikram_rebalance_comp											!Bikram
      LOGICAL bikram_on																			!Bikram
      LOGICAL bikram_mtrx_path																	!Bikram
      LOGICAL bikram_w3j_fact, bikram_truncate_MTRX												!Bikram
! Bikram Start April 2020:
!	  real*8, allocatable :: bk_adia_t(:), bk_adia_r(:), bk_adia_pr(:)			
!	  real*8, allocatable :: bk_adia_theta(:), bk_adia_ptheta(:)					
!	  real*8, allocatable :: bk_adia_phi(:), bk_adia_pphi(:)						  
!	  real*8, allocatable :: bk_adia_r_der(:), bk_adia_pr_der(:)				
!	  real*8, allocatable :: bk_adia_theta_der(:), bk_adia_ptheta_der(:)					
!	  real*8, allocatable :: bk_adia_phi_der(:), bk_adia_pphi_der(:)				
	  real*8, allocatable :: bk_adia_t(:)					
	  real*8, allocatable :: bk_sys_var(:,:), bk_sys_var_der(:,:)					
	  integer bk_adia_n
      LOGICAL bikram_adiabatic																	!Bikram
      LOGICAL bikram_save_traj																	!Bikram
! Bikram End. 


! Bikram Start Nov 2020:	
	  integer bk_parity
	  logical ini_st_ident
! Bikram End. 

      LOGICAL coupled_states_defined
      LOGICAL user_vib_mat_el_defined
      LOGICAL print_states_defined	  
      LOGICAL expansion_defined
      LOGICAL no_orbits_defined	
      LOGICAL pot_expansion_defined	  
      LOGICAL units_defined
      LOGICAL grid_defined
      LOGICAL cartes_defined
      LOGICAL terms_defined
      LOGICAL angs_unit,au_unit_r,au_unit_e,cm_unit,klv_unit,kal_unit
      LOGICAL matrix_reading_defined
      LOGICAL write_check_file_defined	  
      LOGICAL run_prog_defined
      LOGICAL atom_coord_dist
      LOGICAL make_grid_file
      LOGICAL grid_file_found
      LOGICAL calc_expansion_defined
      LOGICAL read_r_grid_from_file_defined
      LOGICAL print_matrix_defined
      LOGICAL test_expansion_defined
      LOGICAL calc_matrix_defined	  
      LOGICAL, ALLOCATABLE :: Ks_have_to_be_skipped(:)
      LOGICAL states_to_exclude_defined
      LOGICAL terms_onfly_defined
      LOGICAL terms_file_defined
      LOGICAL eq_grid_defined
      LOGICAL vib_mix_state_defined	  
      LOGICAL non_format_out_defined	  
      INTEGER numb_orb_prds
      INTEGER numb_oscl_prds																	!Bikram
      INTEGER numb_rk4_stps_adia																!Bikram
      INTEGER ang_res_num	  
      INTEGER nmbr_of_enrgs	
      INTEGER nmbr_of_traj	  
      INTEGER system_type
      INTEGER prn_l_trj	  
      INTEGER J_tot_max
      INTEGER J_tot_min
      INTEGER delta_l_step	  
	  integer, allocatable :: bk_non_zero_mij_gather(:)											!Bikram May 2022
      INTEGER bk_dl_lr																			!Bikram Feb 2021	  
      INTEGER bk_adiabatic_input																!Bikram   
      INTEGER mtrx_cutoff_r1, mtrx_cutoff_r2, rms_r												!Bikram   
      INTEGER mpi_task_per_traject	  
      INTEGER number_of_channels
      INTEGER TIME_MIN_CHECK
      INTEGER coll_type
      INTEGER chann_ini	  
      INTEGER vmax,vmin,jmax,jmin,jmax_included,vmax_included
      INTEGER vmax1,vmin1,jmax1,jmin1
      INTEGER vmax2,vmin2,jmax2,jmin2
      INTEGER nchann_1,nchann_2
      INTEGER LORB_FINE,SPIN_FINE,SPIN_HYPERFINE	  
      INTEGER, ALLOCATABLE :: j_ch(:),ka_ch(:),kc_ch(:),
     & k_ch(:), eps_ch(:), j1_ch(:),ka1_ch(:),kc1_ch(:),
     & k1_ch(:), eps1_ch(:), j2_ch(:),ka2_ch(:),kc2_ch(:),
     & k2_ch(:), eps2_ch(:), v_ch(:), v1_ch(:), v2_ch(:),
     & j12(:),m12(:),chann_indx(:),indx_chann(:), indx_corr(:,:,:),
     & ind_mat(:,:),parity_state(:),indx_corr_id(:,:,:,:)
     & ,j_max_ind(:),j_min_ind(:),f_ch(:),par_lorb_ch(:)
	  integer, allocatable :: parity_state_bk(:)												!Bikram April 2021
	  integer, allocatable :: parity_state_sign_bk(:)											!Bikram April 2021
	  integer, allocatable ::  p1p2_bk(:,:)
! Bikram Start Dec 2019:
	  integer,allocatable :: ind_mat_bk(:,:),bk_indx(:)						
	  integer mat_sz_bk,ph_cntr_bk,splnt_bfr,splnt_afr												
	  real*8,allocatable :: misc_bk(:,:),bk_delta_E(:)
	  real*8,allocatable :: bk_splint_bfr_mat(:),bk_splint_afr_mat(:)
	  real*8,allocatable :: bk_splint_bfr_dmat(:)
	  real*8,allocatable :: bk_splint_afr_dmat(:)
	  real*8,allocatable :: bk_mat_splnt(:,:),bk_dmat_splnt(:,:)
! Bikram End. 
      INTEGER p_lim_max,p_lim_min	 
      INTEGER vib_diat_atom(2)
      INTEGER symm_top_atom(3)
      INTEGER asymm_top_atom(3)
      INTEGER diat_diat(2)
      INTEGER vib_diat_diat(2)
      INTEGER symm_diat(4)
      INTEGER asymm_diat(4)
      INTEGER asymm_symm(6)
      INTEGER asymm_asymm(6)
      INTEGER j_ini,ka_ini,kc_ini,k_ini, eps_ini,
     & j1_ini,ka1_ini,kc1_ini,k1_ini,
     & j2_ini,ka2_ini,kc2_ini,k2_ini,v1_ini,v2_ini,v_ini,
     & eps1_ini, eps2_ini,f_ini,par_orb_ini
      INTEGER nterms
      INTEGER n_r_coll
      INTEGER n_r_vib,n_r_vib1,	n_r_vib2
      INTEGER n_alpha,n_beta, n_gamma
      INTEGER n_alpha1,n_beta1, n_gamma1 
      INTEGER n_alpha2,n_beta2, n_gamma2
      INTEGER, ALLOCATABLE :: L_TERM(:), M_TERM(:),
     & L1_TERM(:),L2_TERM(:), NJU1_TERM(:),NJU2_TERM(:)
      INTEGER LM(2)
      INTEGER L1L2L(3)
      INTEGER L1NJU1L2L(4)	  
      INTEGER L1L2N1N2L(5)
      INTEGER L_MAX_EXPAN,NJU_MAX_EXPAN,L1_MAX_EXPAN,
     &	L2_MAX_EXPAN,NJU1_MAX_EXPAN,NJU2_MAX_EXPAN,
     & L1_MIN_EXPAN,L2_MIN_EXPAN
      INTEGER ir_bgn_exp_pnt,ir_fin_exp_pnt	 
      INTEGER n_2_pnts(2)
      INTEGER m_elastic_proj_print	  
      INTEGER states_size,total_size
      INTEGER, ALLOCATABLE :: parity_inversion(:)
      INTEGER j12m12_print(2)
	  integer bikram_rms_ang1, bikram_rms_ang2, bikram_rms_ang3									!Bikram 
	  integer bikram_axl_sym1, bikram_axl_sym2													!Bikram 
	  integer bikram_equ_sym1, bikram_equ_sym2													!Bikram 
      REAL*8, ALLOCATABLE :: j_h_ch(:),m12_h(:),j12_h(:)	  
      REAL*8 b_impact_parameter	  
      REAL*8 R_max_dist, R_min_dist
      REAL*8 mass_red
      REAL*8, ALLOCATABLE :: U(:)
      REAL*8 time_step
      REAL*8 min_t_stp
      REAL*8 bk_rk4_tol_adia																	!Bikram
      REAL*8 bk_b_switch																		!Bikram
      REAL*8 time_lim
      REAL*8 eps_odeint
      REAL*8 mnt_crl_intgrt_err
      REAL*8 U_max,U_min, dU
      REAL*8 EMAX1,EMAX2	  
      REAL*8 BE,DE,A,B,C,EMAX,HE_ROT
      REAL*8 BE1,DE1,A1,B1,C1
      REAL*8 BE2,DE2,A2,B2,C2
      REAL*8 WE,XE,WE1,WE2,XE1,XE2
      REAL*8 atomic_masses(4)
      REAL*8 atomic_masses1(2),atomic_masses2(2)	  
      REAL*8 exch_par_w_pl,exch_par_w_mn	  
      REAL*8 morse_a, morse_de, morse_re
      REAL*8 morse12_a(2), morse12_de(2), morse12_re(2)
      REAL*8 r_vib_diatom_min,r_vib_diatom_max
      REAL*8 r_vib_diatom_min12(2),r_vib_diatom_max12(2)	  
      REAL*8 r_vib_diatom_min1,r_vib_diatom_max1
      REAL*8 r_vib_diatom_min2,r_vib_diatom_max2
      REAL*8 gamma_fine,lambda_fine,a_spin_orbit_fine
      REAL*8 gamma_fine_d,lambda_fine_d,lambda_fine_dd
      REAL*8 p_double_fine, q_double_fine		  
      REAL*8, ALLOCATABLE :: E_ch(:),Ej(:),Mat_el(:,:),Mat_el_der(:,:)
     & , R_COM(:),Mij_exp_coulp(:,:)
      REAL*8, ALLOCATABLE :: vib_mix_real(:), vib_mix_imag(:)	 
      REAL*8, ALLOCATABLE :: vibrational_wavefunctions(:,:)
      REAL*8, ALLOCATABLE :: vibrational_wavefunctions_1(:,:)
      REAL*8, ALLOCATABLE :: vibrational_wavefunctions_2(:,:)
      REAL*8, ALLOCATABLE :: jacobian_vib(:),jacobian_vib1(:),
     & jacobian_vib2(:)
      REAL*8, ALLOCATABLE :: r_grid_vib(:),r_grid_vib1(:),
     & r_grid_vib2(:)	 
      REAL*8, ALLOCATABLE :: M_EIGEN_VECTORS_USER(:,:)
      REAL*8 r_unit,e_unit
      REAL*8 atomic_red_mass,atomic_red_mass1,atomic_red_mass2
	  real*8 MIJ_ZERO_CUT																		!Bikram
	  real*8 bikram_cutoff_r1, bikram_cutoff_r2, bikram_rms_r									!Bikram
	  character(len = *),parameter :: bk_directory = "AT_APPROX_TRAJS"							!Bikram
	  character(len = *),parameter :: bk_dir11 = "MATRIX_FILES"									!Bikram
	  character(len = *),parameter :: bk_dir22 = "MATRIX_TRUNCATED"								!Bikram
	  character(len = *),parameter :: bk_dir33 = "REBALANCE_FILES/"								!Bikram
	  character (len=500) :: bk_matrix_path11, bk_dir1, bk_dir2									!Bikram
      CHARACTER(LEN=:), ALLOCATABLE :: label
	  character(len = *),parameter :: bk_dir44 = "CHK_FILES"									!Bikram
      END MODULE VARIABLES	  
      SUBROUTINE INITIALIZATION
! READING INPUT FILE.
      USE MPI_DATA
      USE MPI
      USE ERRORS	  
      USE VARIABLES	  
      IMPLICIT NONE
      LOGICAL file_exst	  
      INTEGER*4 iostatus,mypos,mypos1,file_end
      INTEGER*4 sizein,string_len,row_num,i_row
      CHARACTER(LEN=50) :: filename! = "basisH2H2.dat"
      CHARACTER(LEN=:), ALLOCATABLE :: str, s_temp
      CHARACTER(LEN=1) buffer
      INTEGER*4, ALLOCATABLE :: tab_pos(:)
      INTEGER*4 tab_len
      IF(MYID.EQ.0) THEN	  
      sizein = 0
      iostatus = 0
      mypos = 1
      mypos1 = mypos
      file_end = 0
      tab_len = 0
      string_len = 0
      row_num = 0 
      INQUIRE( FILE="INPUT_NAME.inp", EXIST=file_exst ) !! CHECKING IF FILE EXISTS
      IF(.NOT.file_exst) THEN
      PRINT*, "ERROR: INPUT_NAME.inp NOT FOUND"
      error_file_reading = .TRUE.
      GOTO 3523	  
      ENDIF	  
      OPEN(1,FILE="INPUT_NAME.inp") !!! NAME OF INPUT FILE STORED HERE
      READ(1,'(a)') filename
      CLOSE(1)	  
c      STOP "OKAY"
      INQUIRE( FILE=filename, EXIST=file_exst )
      IF(.NOT.file_exst) THEN
      WRITE(*,'(a7,1x,a1,a20,a1,1x,a10)')
     & "ERROR: ","""",filename,""""," NOT FOUND"
      error_file_reading = .TRUE.
      GOTO 3523		 
      ENDIF
      OPEN(UNIT = 1,FILE=filename,ACTION="READ")
      row_num = 0
      DO
      READ (1, *, END=20)
      row_num = row_num + 1
      END DO
20    CONTINUE
      IF(MYID.eq.0 .and. nproc.eq.1) PRINT*,"rows in Input", !!! FOR USER IFNO
     & row_num
      REWIND(1)
      DO i_row=1,row_num
      DO
      READ(1,'(A)',ADVANCE = "NO", EOR=10) buffer
      IF((buffer.ne.' ') .and. (buffer.ne.'	')) THEN 
      string_len = string_len + 1	  
      ENDIF	  
      ENDDO	  
10    CONTINUE
      ENDDO
      ALLOCATE(CHARACTER(LEN=string_len) :: str)	  
      string_len = 0	  
      REWIND(1)	  
      DO i_row=1,row_num
      DO
      READ(1,'(A)',ADVANCE = "NO", EOR=15) buffer
      IF((buffer.ne.' ') .and. (buffer.ne.'	')) THEN
      string_len = string_len + 1	  
      str(string_len:string_len) = buffer	  !!! INPUT SAVED INTO STRING
      ENDIF	  
      ENDDO	  
15    CONTINUE
      ENDDO	  
	  
      CLOSE(1)
      ENDIF
3523  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi) !! READING ENDS HERE
      CALL MPI_BCAST(error_file_reading, 1, MPI_LOGICAL,0,
     &  MPI_COMM_WORLD,ierr_mpi)  !!! BROADCASTING INPUT ONTO OTHER PROCCESSORS
      IF(error_file_reading) STOP	 
      CALL MPI_BCAST(string_len, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(MYID.NE.0) THEN
      ALLOCATE(CHARACTER(LEN=string_len) :: str)  
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(str, string_len, MPI_CHARACTER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
	 
      CALL INPUT_PARSING(str,string_len) !!! PARSING
  
      CALL OUTPUT !!! PRINT YOUR SYSTEM SETUP

      TIME_START = MPI_Wtime()
      CALL INI_ARRAYS
	  !!! ALLOCATE WORK ARRAYS
      TIME_FINISH = MPI_Wtime()
      TIME_1_ARRAY = TIME_FINISH-TIME_START 	  
      CALL INI_FACT    !!! PREPARE WORK ARRAY WICH CONTAINS FACTORIAL VALUES
      CALL MATRIX_MAKE  !!! MAKING COULPING MATRIX
	  return
      END SUBROUTINE INITIALIZATION
      SUBROUTINE MAKE_INP_STRING(filename,s,n,len_file,t_p,t_l,
     &  s_temp,temp_len) !!! SUBROUTINE CLEANS TABS, SPACES AND LINES FROM INPUT
      IMPLICIT NONE
      INTEGER n,len_file,t_l,temp_len
      CHARACTER(LEN=n) s
      CHARACTER(LEN=len_file) filename
      CHARACTER(LEN=temp_len) s_temp
      INTEGER t_p(t_l)
      INTEGER iostatus,i,p1,p2,j
      CHARACTER(LEN=1) buffer
      OPEN(UNIT=1,FILE = filename)
      p1 = 1
      DO i = 1,t_l
      p2 = t_p(i)
      READ (UNIT=1,FMT = '(a)', IOSTAT=iostatus) s_temp(p1:p2)
c      PRINT*,s_temp(p1:p2)
      IF(p2+1 .le.temp_len ) THEN  
      WRITE(s_temp(p2:p2+1),FMT = '(a)') ''
      p1 = p2 + 1
      ENDIF	  
      ENDDO
      CLOSE(1)
c      STOP "OKAY"	  
      j = 0
      DO i=1,temp_len
      buffer = s_temp(i:i)
      IF(buffer.ne.' ' .and. buffer.ne.'	') THEN
      j = j + 1	  
      s(j:j) = buffer
      ENDIF
       	  
      ENDDO
	  
      END SUBROUTINE MAKE_INP_STRING
      SUBROUTINE INPUT_PARSING(input,length)
!!! DIVIDE YOUR INPUT FILE BY BLOCKS AND READ EACH BLOCK!!! INPUT - user input file, length - its length	  
      IMPLICIT NONE
      INTEGER length	  
      CHARACTER(LEN=length) input
      CHARACTER(LEN=5) :: basis_k_word="BASIS"
      CHARACTER(LEN=6) :: system_k_word="SYSTEM"
      CHARACTER(LEN=9) :: pot_k_word = "POTENTIAL"
      CHARACTER(LEN=3) :: end_key = "END"
      CHARACTER(LEN=1) buffer
      CHARACTER(LEN=:), ALLOCATABLE :: basis_inp
      CHARACTER(LEN=:), ALLOCATABLE :: system_inp
      CHARACTER(LEN=:), ALLOCATABLE :: potential_inp	  
      INTEGER posit,key_sig
      INTEGER sys_pos_beg, sys_pos_end
      INTEGER basis_pos_beg, basis_pos_end
      INTEGER pot_pos_beg, pot_pos_end
      LOGICAL ins_block, key_word
      key_word = .FALSE.
      ins_block = .FALSE.
      posit = 0
      key_sig = 0
      DO WHILE(posit.lt.length)
      posit = posit + 1
      buffer = input(posit:posit)
c      PRINT *,posit,ins_block
      IF(buffer.ne."$" .and. (.not.ins_block)) STOP
     & "ERROR IN INPUT. BLOCK STARTS WITHOUT A KEY WORD" 
      IF(buffer.eq.'$' .and. ins_block) THEN
      buffer = input(posit+1:posit+1)
      IF(input(posit+1:posit+3).ne.end_key) THEN
      STOP "WRONG KEY WORD IN THE END OF INPUT BLOCK"
      ELSE
      IF(key_sig.eq.0) 
     & STOP "CRITICAL ERROR.INPUT BLOCK NOT IDENTIFIED"
      SELECT CASE(key_sig)
      CASE(1)
      sys_pos_end = posit-1  
      CASE(2)
      basis_pos_end = posit-1 	  
      CASE(3)
      pot_pos_end = posit-1 	  
      END SELECT
      key_sig = 0
      posit = posit + 3
      ins_block = .FALSE. 	  
      ENDIF
      ENDIF	  
      IF(buffer.eq.'$' .and. (.not.ins_block)) THEN !!! FINDIGN THE PROPER KEW WORD
      buffer = input(posit+1:posit+1)
      ins_block = .TRUE.
      key_word = .TRUE.	 
      IF(input(posit+1:posit+5).eq.basis_k_word) THEN
      basis_pos_beg = posit+6  
      key_sig = 2
      key_word = .FALSE.
      ENDIF	
      IF(input(posit+1:posit+6).eq.system_k_word) THEN
      sys_pos_beg = posit+7  
      key_sig = 1
      key_word = .FALSE.
      ENDIF
      IF(input(posit+1:posit+9).eq.pot_k_word) THEN
      pot_pos_beg = posit+10
      key_sig = 3
      key_word = .FALSE.
      ENDIF	 
      IF(key_sig.eq.0 .and. key_word) THEN
      STOP
     & "THE BLOCK STARTS WITH WRONG KEY WORD"
      ENDIF	 
      ENDIF	  	  
      ENDDO	  
c      PRINT *, "NO_ERRORS_IN_READING, INI_DONE"
      ALLOCATE(CHARACTER(LEN=basis_pos_end-basis_pos_beg+1) ::
     &	  basis_inp)
      ALLOCATE(CHARACTER(LEN=sys_pos_end-sys_pos_beg+1) ::
     &	  system_inp)
      ALLOCATE(CHARACTER(LEN=pot_pos_end-pot_pos_beg+1) ::
     &	  potential_inp)		 
      basis_inp = input(basis_pos_beg:basis_pos_end)
      system_inp = input(sys_pos_beg:sys_pos_end)
      potential_inp = input(pot_pos_beg:pot_pos_end)	  
c      OPEN(UNIT=1,FILE="JUST_BASIS.dat")	  
c      WRITE(1,*) basis_inp
c      CLOSE(1)
c      OPEN(UNIT=1,FILE="JUST_SYSTEM.dat")	  
c      WRITE(1,*) system_inp
c      CLOSE(1)
c      OPEN(UNIT=1,FILE="JUST_POTENTIAL.dat")	  
c      WRITE(1,*) potential_inp
c      CLOSE(1)	  
      CALL SYSTEM_PARSING(system_inp,LEN(system_inp)) !!! PARSING OF SYSTEM BLOCK
c      STOP "SYSTEM PARSING DONE"	  
      CALL BASIS_PARSING(basis_inp,LEN(basis_inp)) !!! PARSING OF BASIS BLOCK
      CALL POTENTIAL_PARSING(potential_inp,LEN(potential_inp))	  !!! PARISNG OF POTENTIAL BLOCK
      END SUBROUTINE INPUT_PARSING
      SUBROUTINE REAL_NUMBER_READING(inp,len_inp,posit,real_numb)
      IMPLICIT NONE !!! READ A REAL NUMBER IN YOUR INPUT. INP - input string, len_inp - length of the string
      REAL*8 real_numb!! posit - position, where reading starts, real_numb - output
      INTEGER len_inp
      CHARACTER(LEN = len_inp) inp
      INTEGER posit,place, end_pos,bgn_pos,decrement
      CHARACTER(LEN =1) buffer
      LOGICAL ERROR_NUM
      IF(len_inp.le.0) RETURN	  
      place = 1	  
      buffer = inp(place:place)
      ERROR_NUM=.FALSE.
      bgn_pos = 1
      decrement = posit - bgn_pos
      DO WHILE(buffer.ne.',' .and. place.le.len_inp)
      IF((ICHAR(buffer).gt.57 .or. ICHAR(buffer).lt.48)
     & .and. buffer.ne.'.') STOP "ERROR: WRONG FORMAT FOR REAL NUMBER"
      IF(buffer.eq.'.') THEN
      IF(.not.ERROR_NUM) THEN
      ERROR_NUM = .TRUE.
      ELSE	  
      STOP "TOO MANY DOTS"	  !!! SEEKING FOR A DOT
      ENDIF
      ENDIF
      place = place + 1
      buffer = inp(place:place)	  
      ENDDO
      end_pos= place - 1
      IF(.not.ERROR_NUM) PRINT *, " ERROR: REAL NUMBER READING FAILED"
      IF(.not.ERROR_NUM) STOP
      IF((end_pos-bgn_pos).lt.0) THEN
!      PRINT*,inp(bgn_pos:end_pos)	
      STOP "NO NUMBER ON INPUT"
      ENDIF	  
      READ(inp(bgn_pos:end_pos),*) real_numb
!      PRINT*, "rmin", real_numb	  
      posit = end_pos + decrement+2	  !!! SHIFTING THE POSTION OF READING	  
      END SUBROUTINE REAL_NUMBER_READING
      SUBROUTINE INT_NUMBERS_READING(inp,len_inp,posit,int_numb,
     & coun) !! READ INTEGER(S) IN YOUR INPUT
      IMPLICIT NONE
      INTEGER coun
      INTEGER int_numb(coun)
      INTEGER len_inp,i
      LOGICAL minus_def	  
      CHARACTER(LEN = len_inp) inp
      INTEGER posit,place, end_pos,bgn_pos,decrement
      CHARACTER(LEN =1) buffer
! NEW	  
      minus_def = .FALSE.!! MINUS INTEGER OR NOT
! NEW	  
      IF(len_inp.le.0) RETURN
      bgn_pos = 1
      decrement = posit - bgn_pos
      DO i=1,coun
      place = bgn_pos
      buffer = inp(place:place)
      DO WHILE(buffer.ne.',' .and. place.le.len_inp) !SEEKING FOR COMA
! NEW	  
      IF(place.eq.bgn_pos) THEN
      IF(buffer.eq.'-') THEN
      minus_def = .TRUE.
      place = place + 1
      buffer = inp(place:place) 	  
      CYCLE	  
      ENDIF	  
      ENDIF
! NEW	  
      IF(ICHAR(buffer).gt.57 .or. ICHAR(buffer).lt.48
     & ) THEN
      PRINT*,"PROBLEM WITH USER INPUT (SEE BELOW THE PIECE)"
      PRINT*,inp	 
!      PRINT*,"buffer = ",buffer	 
      STOP "ERROR: NON-INTEGER NUMBER. CHECK INPUT" !!! CHECKING THAT ONLY ALLOWED CHARACTERS
      ENDIF	 
      place = place + 1
      buffer = inp(place:place)
      ENDDO
      end_pos= place - 1
      IF((end_pos-bgn_pos).lt.0) STOP "ERROR: NO NUMBER ON INPUT"
      READ(inp(bgn_pos:end_pos),*) int_numb(i)
      bgn_pos = end_pos + 2
      ENDDO
      posit = bgn_pos + decrement
  	  
      END SUBROUTINE INT_NUMBERS_READING
      SUBROUTINE KEY_WORD_BASIS(inp,length,key_word_used,key,place)
      IMPLICIT NONE !!! KEY WORDS IN THE INPUT BASIS. SEE MANUAL
      INTEGER, PARAMETER :: num_key_word = 76
      INTEGER length,posit,key_word_used(num_key_word),
     & key,place,decrement,i
      CHARACTER(LEN=length) inp
      CHARACTER(LEN=9) :: type_k_word="SYS_TYPE="                 									!CASE(1)
      CHARACTER(LEN=12) :: user_k_word="LEVELS_FILE="   		  									!CASE(2)
      CHARACTER(LEN=12) :: e_lvl_k_word="CHNL_ENERGS="            									!CASE(3)
      CHARACTER(LEN=10) :: num_ch_k_word="NMB_CHNLS="             									!CASE(4)
      CHARACTER(LEN=11) :: chann_k_word="CHNLS_LIST="             									!CASE(5)
      CHARACTER(LEN=10) :: ini_k_word="INIT_CHNL="                									!CASE(6)
      CHARACTER(LEN=3) :: be_word="BE="                           									!CASE(7)
      CHARACTER(LEN=3) :: de_word="DE="                           									!CASE(8)
      CHARACTER(LEN=2) :: a_word="A="                             									!CASE(9) 
      CHARACTER(LEN=2) :: b_word="B="                             									!CASE(10)
      CHARACTER(LEN=2) :: c_word="C="                             									!CASE(11)
      CHARACTER(LEN=5) :: em_word="EMAX="                         									!CASE(12)
      CHARACTER(LEN=4) :: be1_word="BE1="                         									!CASE(13)
      CHARACTER(LEN=4) :: de1_word="DE1="                         									!CASE(14)
      CHARACTER(LEN=3) :: a1_word="A1="                           									!CASE(15)
      CHARACTER(LEN=3) :: b1_word="B1="                           									!CASE(16)
      CHARACTER(LEN=3) :: c1_word="C1="                           									!CASE(17)
      CHARACTER(LEN=4) :: be2_word="BE2="						  									!CASE(18)
      CHARACTER(LEN=4) :: de2_word="DE2="						  									!CASE(19)
      CHARACTER(LEN=3) :: a2_word="A2="							  									!CASE(20)
      CHARACTER(LEN=3) :: b2_word="B2="							  									!CASE(21)	
      CHARACTER(LEN=3) :: c2_word="C2="                           									!CASE(22)
      CHARACTER(LEN=3) :: we_word="WE="                           									!CASE(23)
      CHARACTER(LEN=3) :: xe_word="XE="                           									!CASE(24)
      CHARACTER(LEN=4) :: we1_word="WE1="                         									!CASE(25)
      CHARACTER(LEN=4) :: we2_word="WE2="                         									!CASE(26)
      CHARACTER(LEN=4) :: xe1_word="XE1="                         									!CASE(27)
      CHARACTER(LEN=4) :: xe2_word="XE2="                         									!CASE(28)
      CHARACTER(LEN=5) :: vmax_word="VMAX="                       									!CASE(29)
      CHARACTER(LEN=5) :: vmin_word="VMIN="                       									!CASE(30)
      CHARACTER(LEN=5) :: jmax_word="JMAX="                       									!CASE(31)
      CHARACTER(LEN=5) :: jmin_word="JMIN="                       									!CASE(32)
      CHARACTER(LEN=6) :: vmax1_word="VMAX1="                     									!CASE(33)
      CHARACTER(LEN=6) :: vmin1_word="VMIN1="                     									!CASE(34)
      CHARACTER(LEN=6) :: jmax1_word="JMAX1="                     									!CASE(35)
      CHARACTER(LEN=6) :: jmin1_word="JMIN1="                     									!CASE(36)
      CHARACTER(LEN=6) :: vmax2_word="VMAX2="                     									!CASE(37)
      CHARACTER(LEN=6) :: vmin2_word="VMIN2="                     									!CASE(38)
      CHARACTER(LEN=6) :: jmax2_word="JMAX2="                     									!CASE(39)
      CHARACTER(LEN=6) :: jmin2_word="JMIN2="                     									!CASE(40)
      CHARACTER(LEN=6) :: emax1_word="EMAX1="                     									!CASE(41)
      CHARACTER(LEN=6) :: emax2_word="EMAX2="                     									!CASE(42)
      CHARACTER(LEN=6) :: nchnl1_word="NCHL1="                    									!CASE(43)
      CHARACTER(LEN=6) :: nchnl2_word="NCHL2="                    									!CASE(44)
      CHARACTER(LEN=10) :: idnt_word="IDENTICAL="                 									!CASE(45)
      CHARACTER(LEN=14) :: masses_of_atoms="ATOMIC_MASSES="       									!CASE(46)
      CHARACTER(LEN=12)	:: morse_width = "MORSE_WIDTH="           									!CASE(47)
      CHARACTER(LEN=12) :: morse_depth = "MORSE_DEPTH="           									!CASE(48)
      CHARACTER(LEN=13) :: morse_posit = "MORSE_POSITN="          									!CASE(49)
      CHARACTER(LEN=13) :: r_vib_word_min = "RMIN_VIBGRID="   	  									!CASE(50)
      CHARACTER(LEN=13) :: r_vib_word_max = "RMAX_VIBGRID="   	  									!CASE(51)
      CHARACTER(LEN=10)	:: cs_word = "CS_APPROX="                 									!CASE(52)
      CHARACTER(LEN=13) :: print_states_word="PRINT_STATES="      									!CASE(53)
      CHARACTER(LEN=15)	 :: excl_states_word = "EXCLUDE_STATES="  									!CASE(54)
      CHARACTER(LEN=12) :: exch_par_word = "WGHT_POSPAR="         									!CASE(55)
      CHARACTER(LEN=9)   :: para_ortho_word = "SYMMETRY="         									!CASE(56)
	  CHARACTER(LEN=11) :: fine_struct_word = "FINE_STRCT="       									!CASE(57)
	  CHARACTER(LEN=13) :: hyp_fine_struct_word = "HPFINE_STRCT=" 									!CASE(58)
	  CHARACTER(LEN=12) :: spin_rot_cop_word ="SPROT_CPLNG="      									!CASE(59)
	  CHARACTER(LEN=12) :: spin_orb_cop_word ="SPORB_CPLNG="      									!CASE(60)
	  CHARACTER(LEN=8) :: el_spin_word ="EL_SPIN="                									!CASE(61)
	  CHARACTER(LEN=8) :: el_orb_word ="EL_ORB="                  									!CASE(62)
	  CHARACTER(LEN=11) :: spin_spin_cop_word ="SPSP_CPLNG="      									!CASE(63)
	  CHARACTER(LEN=12) :: spin_spind_cop_word ="SPSPD_CPLNG="      								!CASE(64)
	  CHARACTER(LEN=13) :: spin_spindd_cop_word ="SPSPDD_CPLNG="    								!CASE(65)
	  CHARACTER(LEN=13) :: spin_rotd_cop_word ="SPROTD_CPLNG="      								!CASE(66)
      CHARACTER(LEN=10) :: fine_level_word="FINE_LIST="             								!CASE(67)
      CHARACTER(LEN=14) :: ini_fine_level_word="INI_FINE_CHNL="     								!CASE(68)
      CHARACTER(LEN=7) :: he_rot_word="HE_ROT="              										!CASE(69)
      CHARACTER(LEN=6) :: p_double_word="P_DBL="              										!CASE(70)
      CHARACTER(LEN=6) :: q_double_word="Q_DBL="              										!CASE(71)
      CHARACTER(LEN=14) :: par_fine_word="PAR_FINE_LIST="           								!CASE(72)
      CHARACTER(LEN=13) :: par_fine_ini_word="PAR_FINE_INI="        								!CASE(73)
      CHARACTER(LEN=12)	:: vib_mix_state_def_word = "MIXED_STATE="  								!CASE(74)  
      CHARACTER(LEN=10)	:: bk_adia = "AT_APPROX="                									!CASE(75)
      CHARACTER(LEN=13)	:: bk_save_traj = "SAVE_TRAJECT="                								!CASE(76)
      CHARACTER(LEN=1) buffer
      LOGICAL key_used
      IF(length.le.0) RETURN
      key_used = .FALSE.
      posit = 1
      decrement = place - posit
      buffer = 	inp(posit:posit)  
      DO WHILE(buffer.ne."=" .and. posit.lt.length)
      posit = posit + 1
      buffer = 	inp(posit:posit)
      ENDDO
      place = posit + decrement + 1	  
!      PRINT*,"WANNA CHECK",key,inp(1:posit)	  
      IF(inp(1:posit).eq.type_k_word) THEN
      key = 1
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	  
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.user_k_word) THEN
      key = 2
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.e_lvl_k_word) THEN
      key = 3
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.num_ch_k_word) THEN
      key = 4
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.chann_k_word) THEN
      key = 5
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.ini_k_word) THEN
      key = 6
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.be_word) THEN
      key = 7
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.de_word) THEN
      key = 8
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.a_word) THEN
      key = 9
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.b_word) THEN
      key = 10
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.c_word) THEN
      key = 11
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.em_word) THEN
      key = 12
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	
      IF(inp(1:posit).eq.be1_word) THEN
      key = 13
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.de1_word) THEN
      key = 14
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.a1_word) THEN
      key = 15
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.b1_word) THEN
      key = 16
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.c1_word) THEN
      key = 17
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.be2_word) THEN
      key = 18
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.de2_word) THEN
      key = 19
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.a2_word) THEN
      key = 20
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.b2_word) THEN
      key = 21
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.c2_word) THEN
      key = 22
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.we_word) THEN
      key = 23
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.xe_word) THEN
      key = 24
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.we1_word) THEN
      key = 25
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.we2_word) THEN
      key = 26
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.xe1_word) THEN
      key = 27
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.xe2_word) THEN
      key = 28
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.vmax_word) THEN
      key = 29
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.vmin_word) THEN
      key = 30
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	 
      IF(inp(1:posit).eq.jmax_word) THEN
      key = 31
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	 
      IF(inp(1:posit).eq.jmin_word) THEN
      key = 32
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.vmax1_word) THEN
      key = 33
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.vmin1_word) THEN
      key = 34
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	 
      IF(inp(1:posit).eq.jmax1_word) THEN
      key = 35
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	 
      IF(inp(1:posit).eq.jmin1_word) THEN
      key = 36
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.vmax2_word) THEN
      key = 37
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.vmin2_word) THEN
      key = 38
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	 
      IF(inp(1:posit).eq.jmax2_word) THEN
      key = 39
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	 
      IF(inp(1:posit).eq.jmin2_word) THEN
      key = 40
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	
      IF(inp(1:posit).eq.emax1_word) THEN
      key = 41
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.emax2_word) THEN
      key = 42
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.nchnl1_word) THEN
      key = 43
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.nchnl2_word) THEN
      key = 44
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.idnt_word) THEN
      key = 45
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.masses_of_atoms) THEN
      key = 46
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.morse_width) THEN
      key = 47
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.morse_depth) THEN
      key = 48
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.morse_posit) THEN
      key = 49
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.r_vib_word_min) THEN
      key = 50
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	  
      IF(inp(1:posit).eq.r_vib_word_max) THEN
      key = 51
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.cs_word) THEN
      key = 52
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.print_states_word) THEN
      key = 53
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.excl_states_word) THEN
      key = 54
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.exch_par_word) THEN
      key = 55
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.para_ortho_word) THEN
      key = 56
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	
      IF(inp(1:posit).eq.fine_struct_word) THEN
      key = 57
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(inp(1:posit).eq.hyp_fine_struct_word) THEN
      key = 58
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
	  IF(inp(1:posit).eq.spin_rot_cop_word) THEN
      key = 59
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
	  IF(inp(1:posit).eq.spin_orb_cop_word) THEN
      key = 60
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	  
	  IF(inp(1:posit).eq.el_spin_word) THEN
      key = 61
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	  
	  IF(inp(1:posit).eq.el_orb_word) THEN
      key = 62
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	
	  IF(inp(1:posit).eq.spin_spin_cop_word) THEN
      key = 63
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF		  
	  IF(inp(1:posit).eq.spin_spind_cop_word) THEN
      key = 64
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
	  IF(inp(1:posit).eq.spin_spindd_cop_word) THEN
      key = 65
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
	  IF(inp(1:posit).eq.fine_level_word) THEN
      key = 67
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF		  
	  IF(inp(1:posit).eq.spin_rotd_cop_word) THEN
      key = 66
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
	  IF(inp(1:posit).eq.ini_fine_level_word) THEN
      key = 68
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
	  IF(inp(1:posit).eq.he_rot_word) THEN
      key = 69
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	
	  IF(inp(1:posit).eq.p_double_word) THEN
      key = 70
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
	  IF(inp(1:posit).eq.q_double_word) THEN
      key = 71
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
	  IF(inp(1:posit).eq.par_fine_word) THEN
      key = 72
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
	  IF(inp(1:posit).eq.par_fine_ini_word) THEN
      key = 73
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	
	  IF(inp(1:posit).eq.vib_mix_state_def_word) THEN
      key = 74
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	 	        
	  IF(inp(1:posit).eq.bk_adia) THEN
      key = 75
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF	        
	  IF(inp(1:posit).eq.bk_save_traj) THEN
      key = 76
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF	
      key_word_used(key) = 1
      RETURN	  
      ENDIF
      IF(.not.key_used) THEN
      PRINT*, inp(1:posit)	  
      STOP 
     & "ERROR IN BASIS: WORD NOT FOUND OR INPUT ABOVE IS INCORRECT"
      ENDIF

      END SUBROUTINE KEY_WORD_BASIS	  
      SUBROUTINE BASIS_PARSING(basis_inp,len_inp)
      USE VARIABLES !!! HERE WE READ THOSE WORDS. FOR EACH WORD WE ASSIGN A KEY!
      IMPLICIT NONE
      INTEGER len_inp
      CHARACTER(LEN = len_inp) basis_inp 
      INTEGER i
      INTEGER posit
      INTEGER,PARAMETER :: KEY_WORD_NUM = 76
      INTEGER key_words(KEY_WORD_NUM),key
      posit = 1
!!!! INI SETUP	  
      BE = 0d0
      DE = 0d0
	  a_spin_orbit_fine = 0d0 
      lambda_fine =0d0
      lambda_fine_d = 0d0
      lambda_fine_dd = 0d0
      gamma_fine = 0d0
      gamma_fine_d = 0d0
      A = 0d0
      B = 0d0
      C = 0d0
      BE1 = 0d0
      DE1 = 0d0
      A1 = 0d0
      B1 = 0d0
      C1 = 0d0
      BE2 = 0d0
      DE2 = 0d0
      A2 = 0d0
      B2 = 0d0
      C2 = 0d0
      EMAX = 0d0
      EMAX1 = 0d0
      EMAX2 = 0d0	  
      WE = 0d0
      XE = 0d0
      WE1 = 0d0
      XE1 = 0d0
      WE2 = 0d0
      XE2 = 0d0
      jmax = -1
      jmin = -1
      vmax = -1
      vmin = -1
      jmax1 = -1
      jmin1 = -1
      jmax2 = -1
      jmin2 = -1
      vmax1 = -1
      vmax2 = -1
      vmin1 = -1
      vmin2 = -1	  
      number_of_channels  = 0
      bk_adiabatic_input  = 0
      nchann_1 = 0
      nchann_2 = 0 	  
      key_words = 0d0
      jmax_included = 0
      vmax_included = 0	  
      coll_type = -1 !! NOT DEFINED YET
      user_defined = .FALSE.
      energy_defined= .FALSE.
      sys_type_defined= .FALSE.
      channels_defined= .FALSE.
      ini_chann_defined= .FALSE.
      all_level_included=.FALSE.
      rot_const_defined=.FALSE.
      emax_defined =.FALSE.
      number_of_channels_defined = .FALSE.
      vib_const_defined = .FALSE.
      morse_pot_defined = .FALSE.
      mlc_mlc_chn_num_defined = .FALSE.
      mlc_mlc_emax_defined	= .FALSE.
      identical_particles_defined = .FALSE.
      atomic_masses_defined = .FALSE.
      coupled_states_defined = .FALSE.
      bikram_adiabatic = .FALSE.										!Bikram April 2020
      bikram_save_traj = .FALSE.										!Bikram April 2020
      user_vib_mat_el_defined = .FALSE.
      states_to_exclude_defined = .FALSE.
      print_states_defined  = .FALSE.	  
      para_ortho_defined = .FALSE.
      fine_structure_defined  = .FALSE.
      hyper_fine_stricture_defined = .FALSE.
      vib_mix_state_defined = .FALSE.
      atomic_masses = 0d0
      morse_a = 0d0
      morse_de = 0d0
      morse_re = 0d0
      morse12_a = 0d0
      morse12_de = 0d0
      morse12_re = 0d0
      HE_ROT = 0d0 	  
      r_vib_diatom_min = -1d0
      r_vib_diatom_max = -1d0
      r_vib_diatom_min1 = -1d0
      r_vib_diatom_max1 = -1d0
      r_vib_diatom_min2 = -1d0
      r_vib_diatom_max2 = -1d0	  
      exch_par_w_pl = 1d0
      exch_par_w_mn = 0d0
      p_double_fine = 0d0
      q_double_fine = 0d0
      SPIN_FINE = 1
      LORB_FINE = 0	  
!!! INI SETUP	  
      DO WHILE(posit.lt.len_inp)
      CALL KEY_WORD_BASIS(basis_inp(posit:len_inp),
     & len_inp - posit + 1,key_words,key,posit) !!! EVERY TIME READING POSTION IS CHANGING
c      PRINT*, "position changing", posit,"key=",key,
c     & "buffer=",basis_inp(posit:posit)
      IF(key.gt.1) THEN 
      IF(coll_type.lt.0) STOP "ERROR IN SYS_TYPE_READING"
      ENDIF	  
      SELECT CASE(key)
      CASE(1)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,coll_type,
     & 1)
      IF(coll_type.ge.0) sys_type_defined = .TRUE.
      IF(coll_type.ge.20) THEN
      coupled_states_defined = .TRUE.
      coll_type = coll_type - 20      	  
      ENDIF
      CASE(2)	 
      IF(basis_inp(posit:posit+2).eq."YES") user_defined = .TRUE.
      IF(basis_inp(posit:posit+1).eq."NO") user_defined = .FALSE.	  
      IF((basis_inp(posit:posit+2).ne."YES") .and.
     & (basis_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(2,1) !!! IF INCORRECT INPUT - STOP AND SIGNAL ERROR
      IF(user_defined) posit = posit + 4	 
      IF(.not.user_defined) posit = posit + 3
      CASE(3)
      IF(number_of_channels.le.0) CALL ERROR_SIGNALING(3,1)
      ALLOCATE(E_ch(number_of_channels))
      DO i=1,number_of_channels
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,E_ch(i))	
      ENDDO	  
      energy_defined = .TRUE.
      CASE(4)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,number_of_channels,
     & 1)
      IF(number_of_channels.le.0) CALL ERROR_SIGNALING(3,1)
      number_of_channels_defined = .TRUE.
      CASE(5)
      IF(.not.number_of_channels_defined)
     & 	  STOP "ERROR: NUMBER OF CHANNELS NOT SPECIFIED"	  
      SELECT CASE(coll_type)
      CASE(1)
      ALLOCATE(j_ch(number_of_channels)) 
      DO i=1,number_of_channels
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,j_ch(i),1)
      ENDDO 
      CASE(2)
      ALLOCATE(j_ch(number_of_channels))
      ALLOCATE(v_ch(number_of_channels)) 	  
      DO i=1,number_of_channels
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,vib_diat_atom,2)	 
      j_ch(i) = vib_diat_atom(2)
      v_ch(i) = vib_diat_atom(1)
      ENDDO 	  
      CASE(3)
      ALLOCATE(j_ch(number_of_channels))
      ALLOCATE(k_ch(number_of_channels))
      ALLOCATE(eps_ch(number_of_channels))  	  
      DO i=1,number_of_channels
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,symm_top_atom,3)	 
      j_ch(i) = symm_top_atom(1)
      k_ch(i) = symm_top_atom(2)
      eps_ch(i) = symm_top_atom(3)	  
      ENDDO 	 	  
      CASE(4)
      ALLOCATE(j_ch(number_of_channels))
      ALLOCATE(ka_ch(number_of_channels))
      ALLOCATE(kc_ch(number_of_channels))  	  
      DO i=1,number_of_channels
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,asymm_top_atom,3)	 
      j_ch(i) = asymm_top_atom(1)
      ka_ch(i) = asymm_top_atom(2)
      kc_ch(i) = asymm_top_atom(3)	  
      ENDDO	  
      CASE(5)
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))
      DO i=1,number_of_channels
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,diat_diat,2)	 
      j1_ch(i) = diat_diat(1)
      j2_ch(i) = diat_diat(2)
c      STOP "HERE WE ARE DONE"
      ENDDO		  
      CASE(6)
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))
      ALLOCATE(v1_ch(number_of_channels))
      ALLOCATE(v2_ch(number_of_channels))  
      DO i=1,number_of_channels
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,vib_diat_diat,4)	 
      v1_ch(i) = vib_diat_diat(1)
      j1_ch(i) =  vib_diat_diat(2)
      v2_ch(i) = vib_diat_diat(3)
      j2_ch(i) =  vib_diat_diat(4)
      ENDDO			  	  
      CASE(7)
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))
      ALLOCATE(k1_ch(number_of_channels))
      ALLOCATE(eps1_ch(number_of_channels))
      DO i=1,number_of_channels
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,symm_diat,4)	 
      j1_ch(i) = symm_diat(1)
      k1_ch(i) =  symm_diat(2)
      eps1_ch(i) = symm_diat(3)
      j2_ch(i) =  symm_diat(4)	  
      ENDDO			  
      CASE(8)
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))
      ALLOCATE(ka1_ch(number_of_channels))
      ALLOCATE(kc1_ch(number_of_channels))
      DO i=1,number_of_channels
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,asymm_diat,4)	 
      j1_ch(i) = asymm_diat(1)
      ka1_ch(i) =  asymm_diat(2)
      kc1_ch(i) = asymm_diat(3)
      j2_ch(i) =  asymm_diat(4)	  
      ENDDO			  
      CASE(9)
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))
      ALLOCATE(ka1_ch(number_of_channels))
      ALLOCATE(kc1_ch(number_of_channels))
      ALLOCATE(k2_ch(number_of_channels))
      ALLOCATE(eps2_ch(number_of_channels))	  
      DO i=1,number_of_channels
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,asymm_symm,6)	 
      j1_ch(i) = asymm_symm(1)
      ka1_ch(i) =  asymm_symm(2)
      kc1_ch(i) = asymm_symm(3)
      j2_ch(i) = asymm_symm(4)
      k2_ch(i) =  asymm_symm(5)
      eps2_ch(i) = asymm_symm(6)
      ENDDO		  
      CASE(0)
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))
      ALLOCATE(ka1_ch(number_of_channels))
      ALLOCATE(kc1_ch(number_of_channels))
      ALLOCATE(ka2_ch(number_of_channels))
      ALLOCATE(kc2_ch(number_of_channels))	  
      DO i=1,number_of_channels
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,asymm_asymm,6)	 
      j1_ch(i) = asymm_asymm(1)
      ka1_ch(i) =  asymm_asymm(2)
      kc1_ch(i) = asymm_asymm(3)
      j2_ch(i) = asymm_asymm(4)
      ka2_ch(i) =  asymm_asymm(5)
      kc2_ch(i) = asymm_asymm(6)
      ENDDO		  
      END SELECT
      channels_defined = .TRUE.  	  
      CASE(6)
      IF(basis_inp(posit:posit+2).eq."ALL") THEN
      all_level_included = .TRUE.
      posit = posit + 3
      ELSE
      SELECT CASE(coll_type)
      CASE(1)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,j_ini,1)
      CASE(2)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,vib_diat_atom,2)	 
      j_ini = vib_diat_atom(2)
      v_ini = vib_diat_atom(1)	  
      CASE(3)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,symm_top_atom,3)	 
      j_ini = symm_top_atom(1)
      k_ini = symm_top_atom(2)
      eps_ini = symm_top_atom(3)	  
      CASE(4)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,asymm_top_atom,3)	 
      j_ini = asymm_top_atom(1)
      ka_ini = asymm_top_atom(2)
      kc_ini = asymm_top_atom(3)	  
      CASE(5)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,diat_diat,2)	 
      j1_ini = diat_diat(1)
      j2_ini = diat_diat(2)
      CASE(6)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,vib_diat_diat,4)	 
      v1_ini = vib_diat_diat(1)
      j1_ini =  vib_diat_diat(2)
      v2_ini = vib_diat_diat(3)
      j2_ini =  vib_diat_diat(4)	  
      CASE(7)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,symm_diat,4)	 
      j1_ini = symm_diat(1)
      k1_ini =  symm_diat(2)
      eps1_ini = symm_diat(3)
      j2_ini =  symm_diat(4)	  
      CASE(8)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,asymm_diat,4)	 
      j1_ini = asymm_diat(1)
      ka1_ini =  asymm_diat(2)
      kc1_ini = asymm_diat(3)
      j2_ini =  asymm_diat(4)	  
      CASE(9)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,asymm_symm,6)	 
      j1_ini = asymm_symm(1)
      ka1_ini =  asymm_symm(2)
      kc1_ini = asymm_symm(3)
      j2_ini = asymm_symm(4)
      k2_ini =  asymm_symm(5)
      eps2_ini = asymm_symm(6)
      CASE(0)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,asymm_asymm,6)	 
      j1_ini = asymm_asymm(1)
      ka1_ini =  asymm_asymm(2)
      kc1_ini = asymm_asymm(3)
      j2_ini = asymm_asymm(4)
      ka2_ini =  asymm_asymm(5)
      kc2_ini = asymm_asymm(6)
      END SELECT
      ENDIF  
      ini_chann_defined = .TRUE.
      CASE(7)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,BE)	
      rot_const_defined	= .TRUE. 
      CASE(8)
      IF(BE.le.0d0) CALL ERROR_SIGNALING(8,1)
      CALL REAL_E_FMT_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,DE)	  
      CASE(9)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,A)
      rot_const_defined	= .TRUE.	 
      CASE(10)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,B)
      rot_const_defined	= .TRUE.	 	  
      CASE(11)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,C)
      rot_const_defined	= .TRUE.	 	  
      CASE(12)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,EMAX)
      emax_defined = .TRUE.	 
      CASE(13)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,BE1)
      rot_const_defined	= .TRUE.	  
      CASE(14)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,DE1)
      rot_const_defined	= .TRUE.	 
      CASE(15)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,A1)
      rot_const_defined	= .TRUE.	  
      CASE(16)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,B1)
      rot_const_defined	= .TRUE.	 
      CASE(17)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,C1)
      rot_const_defined	= .TRUE.	 
      CASE(18)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,BE2)
      rot_const_defined	= .TRUE.	 
      CASE(19)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,DE2)
      rot_const_defined	= .TRUE.	 
      CASE(20)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,A2)
      rot_const_defined	= .TRUE.	 
      CASE(21)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,B2)
      rot_const_defined	= .TRUE.	 
      CASE(22)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,C2)
c      PRINT*, "C2=","defined",C2	 
      rot_const_defined	= .TRUE.
      CASE(23)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,WE)
      vib_const_defined = .TRUE.
      CASE(24)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,XE)
      vib_const_defined = .TRUE.
      morse_pot_defined = .TRUE.	  
      CASE(25)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,WE1)
      vib_const_defined = .TRUE.
      CASE(26)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,WE2)
      vib_const_defined = .TRUE.
      CASE(27)	  
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,XE1)
      vib_const_defined = .TRUE.
      morse_pot_defined = .TRUE.	  
      CASE(28)	  
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,XE2)
      vib_const_defined = .TRUE.
      morse_pot_defined = .TRUE.	  
      CASE(29)	 	  
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,vmax,
     & 1)
      channels_defined = .TRUE.	 
      CASE(30)	 	  
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,vmin,
     & 1)
      channels_defined = .TRUE.	 
      CASE(31)	 	  
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,jmax,
     & 1)
      channels_defined = .TRUE.	 
      CASE(32)	 	  
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,jmin,
     & 1)
      channels_defined = .TRUE.
      CASE(33)	 	  
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,vmax1,
     & 1)
      channels_defined = .TRUE.	 
      CASE(34)	 	  
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,vmin1,
     & 1)
      channels_defined = .TRUE.	 
      CASE(35)	 	  
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,jmax1,
     & 1)
      channels_defined = .TRUE.	 
      CASE(36)	 	  
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,jmin1,
     & 1)
      channels_defined = .TRUE.
      CASE(37)	 	  
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,vmax2,
     & 1)
      channels_defined = .TRUE.	 
      CASE(38)	 	  
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,vmin2,
     & 1)
      channels_defined = .TRUE.	 
      CASE(39)	 	  
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,jmax2,
     & 1)
      channels_defined = .TRUE.	 
      CASE(40)	 	  
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,jmin2,
     & 1)
      channels_defined = .TRUE.
      CASE(41)	 	  
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,EMAX1)
      IF(EMAX2.GT.0d0) mlc_mlc_emax_defined = .TRUE.
      CASE(42)	 	  
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,EMAX2)
      IF(EMAX1.GT.0d0) mlc_mlc_emax_defined = .TRUE.
      CASE(43)	 	  
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,nchann_1,
     & 1)
      IF(nchann_2.gt.0)
     & mlc_mlc_chn_num_defined = .TRUE.
      CASE(44)	 	  
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,nchann_2,
     & 1)
      IF(nchann_1.gt.0)
     & mlc_mlc_chn_num_defined = .TRUE.
      CASE(45)
      IF(basis_inp(posit:posit+2).eq."YES")
     & identical_particles_defined = .TRUE.
      IF(basis_inp(posit:posit+1).eq."NO")
     & identical_particles_defined = .FALSE.	  
      IF((basis_inp(posit:posit+2).ne."YES") .and.
     & (basis_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(45,1)
      IF(identical_particles_defined) posit = posit + 4	 
      IF(.not.identical_particles_defined) posit = posit + 3
      CASE(46)
      IF(coll_type.eq.2) THEN	  
      DO i = 1,2 
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,atomic_masses(i))
      ENDDO
      ENDIF
      IF(coll_type.eq.6) THEN
      DO i = 1,4 
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,atomic_masses(i))
      ENDDO
      atomic_masses1 = atomic_masses(1:2)
      atomic_masses2 = atomic_masses(3:4) 	  
      ENDIF	  
      atomic_masses_defined = .TRUE.	  
      CASE(47)
      IF(coll_type.eq.2) THEN	  
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,morse_a)
      ENDIF
      IF(coll_type.eq.6) THEN
      DO i=1,2	  
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,morse12_a(i))
      ENDDO	 
      ENDIF	  
      morse_pot_defined = .TRUE.	 
      CASE(48)
      IF(coll_type.eq.2) THEN	  
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,morse_de)
      ENDIF
      IF(coll_type.eq.6) THEN
      DO i=1,2	  
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,morse12_de(i))
      ENDDO	 
      ENDIF	  
      morse_pot_defined = .TRUE.	 
      CASE(49)
      IF(coll_type.eq.2) THEN	  
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,morse_re)
      ENDIF
      IF(coll_type.eq.6) THEN
      DO i=1,2	  
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,morse12_re(i))
      ENDDO	 
      ENDIF	  
      morse_pot_defined = .TRUE.
      CASE(50)
      IF(coll_type.eq.2) THEN	  
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,r_vib_diatom_min)
      ENDIF	  
      IF(coll_type.eq.6) THEN
      DO i = 1,2 
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,r_vib_diatom_min12(i))
      ENDDO
      r_vib_diatom_min1 = r_vib_diatom_min12(1)
      r_vib_diatom_min2 = r_vib_diatom_min12(2)	  
      ENDIF	  
      CASE(51)
      IF(coll_type.eq.2) THEN	  
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,r_vib_diatom_max)
      ENDIF
      IF(coll_type.eq.6) THEN
      DO i = 1,2 
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,r_vib_diatom_max12(i))
      ENDDO
      r_vib_diatom_max1 = r_vib_diatom_max12(1)
      r_vib_diatom_max2 = r_vib_diatom_max12(2)	  
      ENDIF	  	  
      CASE(52)
      IF(basis_inp(posit:posit+2).eq."YES") 
     & coupled_states_defined = .TRUE.
      IF(basis_inp(posit:posit+1).eq."NO") 
     & coupled_states_defined = .FALSE. 
      IF((basis_inp(posit:posit+2).ne."YES") .and.
     & (basis_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(52,1)
      IF(coupled_states_defined) posit = posit + 4	 
      IF(.not.coupled_states_defined) posit = posit + 3
      CASE(53)
      IF(basis_inp(posit:posit+2).eq."YES") 
     & print_states_defined = .TRUE.
      IF(basis_inp(posit:posit+1).eq."NO") 
     & print_states_defined = .FALSE. 
      IF((basis_inp(posit:posit+2).ne."YES") .and.
     & (basis_inp(posit:posit+1).ne."NO"))CALL ERROR_SIGNALING(53,1)
      IF(print_states_defined) posit = posit + 4	 
      IF(.not.print_states_defined) posit = posit + 3
      CASE(54)
      IF(basis_inp(posit:posit+2).eq."YES") 
     & states_to_exclude_defined = .TRUE.
      IF(basis_inp(posit:posit+1).eq."NO") 
     & states_to_exclude_defined = .FALSE. 
! Bikram Oct'18 Start:	 
      IF((basis_inp(posit:posit+2).ne."YES") .and.
     & (basis_inp(posit:posit+1).ne."NO"))CALL ERROR_SIGNALING(54,1)
! Bikram End.
      IF(states_to_exclude_defined) posit = posit + 4	 
      IF(.not.states_to_exclude_defined) posit = posit + 3
      CASE(55)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,exch_par_w_pl)
      IF(exch_par_w_pl.gt.1d0) 	exch_par_w_pl = 1d0
      IF(exch_par_w_pl.lt.0d0) exch_par_w_pl = 0d0
      exch_par_w_mn = 1d0 - exch_par_w_pl   
      IF(.not.identical_particles_defined) THEN
      WRITE(*,'(a53,1x,a19)')
     & "ERROR: SPIN WEIGHT STATITICS MUST BE DEFINED ONLY FOR",
     & "IDENTICAL PARTICLES"
      STOP	 
      ENDIF
      CASE(56)
      IF(basis_inp(posit:posit+2).eq."YES") 
     & para_ortho_defined = .TRUE.
      IF(basis_inp(posit:posit+1).eq."NO") 
     & para_ortho_defined = .FALSE. 
      IF((basis_inp(posit:posit+2).ne."YES") .and.
     & (basis_inp(posit:posit+1).ne."NO"))CALL ERROR_SIGNALING(53,1)
      IF(para_ortho_defined) posit = posit + 4	 
      IF(.not.para_ortho_defined) posit = posit + 3
      CASE(57)
      IF(basis_inp(posit:posit+2).eq."YES") 
     & fine_structure_defined = .TRUE.
      IF(basis_inp(posit:posit+1).eq."NO") 
     & fine_structure_defined = .FALSE. 
      IF((basis_inp(posit:posit+2).ne."YES") .and.
     & (basis_inp(posit:posit+1).ne."NO"))CALL ERROR_SIGNALING(53,1)
      IF(fine_structure_defined) posit = posit + 4	 
      IF(.not.fine_structure_defined) posit = posit + 3
      CASE(58)
      IF(basis_inp(posit:posit+2).eq."YES") 
     & hyper_fine_stricture_defined = .TRUE.
      IF(basis_inp(posit:posit+1).eq."NO") 
     & hyper_fine_stricture_defined = .FALSE. 
      IF((basis_inp(posit:posit+2).ne."YES") .and.
     & (basis_inp(posit:posit+1).ne."NO"))CALL ERROR_SIGNALING(53,1)
      IF(hyper_fine_stricture_defined) posit = posit + 4	 
      IF(.not.hyper_fine_stricture_defined) posit = posit + 3
      CASE(59)
      CALL REAL_E_FMT_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,gamma_fine)	  
      CASE(60)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,a_spin_orbit_fine)	  
      CASE(61)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,SPIN_FINE,
     & 1)	  
      CASE(62)
	  CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,posit,LORB_FINE,
     & 1)	 
      CASE(63)
      CALL REAL_NUMBER_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,lambda_fine)
      CASE(64)
      CALL REAL_E_FMT_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,lambda_fine_d)
      CASE(65)
      CALL REAL_E_FMT_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,lambda_fine_dd)
      CASE(66)
      CALL REAL_E_FMT_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,gamma_fine_d)	
      CASE(67)
      SELECT CASE(coll_type)
      CASE(1)	  
      ALLOCATE(f_ch(number_of_channels)) 
      DO i=1,number_of_channels
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,f_ch(i),1)
	  ENDDO
      END SELECT 
      CASE(68)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,f_ini,1)
      CASE(69)
      CALL REAL_E_FMT_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,HE_ROT)
      CASE(70)
      CALL REAL_E_FMT_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,p_double_fine)	  
      CASE(71)
      CALL REAL_E_FMT_READING(basis_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,q_double_fine)
      CASE(72)	 
      SELECT CASE(coll_type)
      CASE(1)	  
      ALLOCATE(par_lorb_ch(number_of_channels)) 
      DO i=1,number_of_channels
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,par_lorb_ch(i),1)
	  ENDDO
      END SELECT
      CASE(73)
      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
     & len_inp+1-posit
     & ,posit,par_orb_ini,1)
      CASE(74)
      IF(basis_inp(posit:posit+2).eq."YES") 
     & vib_mix_state_defined= .TRUE.
      IF(basis_inp(posit:posit+1).eq."NO") 
     & vib_mix_state_defined = .FALSE. 
      IF((basis_inp(posit:posit+2).ne."YES") .and.
     & (basis_inp(posit:posit+1).ne."NO"))CALL ERROR_SIGNALING(53,1)
      IF(vib_mix_state_defined) posit = posit + 4	 
      IF(.not.vib_mix_state_defined) posit = posit + 3	  
      CASE(75)
      IF(basis_inp(posit:posit+2).eq."YES") 
     & bikram_adiabatic = .TRUE.
      IF(basis_inp(posit:posit+1).eq."NO") 
     & bikram_adiabatic = .FALSE. 
      IF((basis_inp(posit:posit+2).ne."YES") .and.
     & (basis_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(52,1)
      IF(bikram_adiabatic) posit = posit + 4	 
      IF(.not.bikram_adiabatic) posit = posit + 3

!      CALL INT_NUMBERS_READING(basis_inp(posit:len_inp),
!     & len_inp-posit+1,posit,bk_adiabatic_input,
!     & 1)
!      IF(bk_adiabatic_input.lt.0 .or. bk_adiabatic_input.gt.2) 
!     & stop "Calculation should be either CC-MQCT/AT-MQCT"
!	  if(bk_adiabatic_input.eq.1) bikram_save_traj = .true.
!	  if(bk_adiabatic_input.eq.2) bikram_adiabatic = .true.
      CASE(76)
      IF(basis_inp(posit:posit+2).eq."YES") 
     & bikram_save_traj = .TRUE.
      IF(basis_inp(posit:posit+1).eq."NO") 
     & bikram_save_traj = .FALSE. 
      IF((basis_inp(posit:posit+2).ne."YES") .and.
     & (basis_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(52,1)
      IF(bikram_save_traj) posit = posit + 4	 
      IF(.not.bikram_save_traj) posit = posit + 3
      END SELECT 
      ENDDO
      system_type = coll_type	  
      END SUBROUTINE BASIS_PARSING
      SUBROUTINE SYSTEM_PARSING(system_inp,len_inp)
      USE VARIABLES	  !!! READING BLOCK SYSTEM
      IMPLICIT NONE	  
      INTEGER len_inp
      CHARACTER(LEN = len_inp) system_inp 
      INTEGER, PARAMETER :: key_num = 46 
      INTEGER posit,key,key_words(key_num)
      INTEGER lbl_bgn,lbl_end
      INTEGER i	  
      CHARACTER(LEN=1) buffer
      tm_lim_defined = .FALSE.
      rk4_defined = .TRUE.
      odeint_defined = .FALSE.
      scatter_param_defined = .FALSE.
      ener_man_defined = .FALSE.
      ener_auto_defined = .FALSE.
      deflect_fun_defined = .TRUE.
      orbit_traj_defined = .TRUE.
      prn_l_defined = .FALSE.
      b_impact_defined = .FALSE.
      mpi_task_defined = .TRUE.				!Bikram
      check_point_defined = .FALSE.
      write_check_file_defined = .FALSE.
!      read_from_matrix_defined = .FALSE.
      monte_carlo_defined = .FALSE.
	  bikram_int=.FALSE. 			!Bikram
	  bikram_theta=.FALSE. 			!Bikram
	  bikram_print=.false.			!Bikram
	  bk_prob_interpolation=.false.			!Bikram
	  bk_nrg_err=.true.				!Bikram
	  bk_step_size=.false.				!Bikram
      calc_elast_defined = .TRUE.
      no_orbits_defined = .FALSE.
      diff_cross_defined = .FALSE.
      non_format_out_defined = .FALSE.
      cartes_defined = .FALSE.	  
      delta_l_step = 1  
      bk_dl_lr = 1  											!Bikram Feb 2021
      TIME_MIN_CHECK = -1	  
      posit = 1
      nmbr_of_traj = 100
      mpi_task_per_traject = 1	  
      numb_orb_prds	= 1
      numb_oscl_prds = 1				!(Bikram)	  
      numb_rk4_stps_adia = 5			!(Bikram)	  
	  bk_rk4_tol_adia = 0.50d0 			!Bikram
	  bk_b_switch = -10.0d0 			!Bikram
      key_words = 0
      mass_red = 0d0
      R_min_dist = 0d0
      R_max_dist = 0d0
      J_tot_max = 0
      J_tot_min = 0
      U_min = 0d0
      U_max = 0d0
      dU = 0d0
      j12m12_print(1) = -1
      j12m12_print(2) = -1
      	  
      time_lim = -100d0	  
      nmbr_of_enrgs = -1
      prn_l_trj = -1
      b_impact_parameter = -1d0
      ang_res_num = 1000	  
      DO WHILE(posit.lt.len_inp)
      CALL KEY_WORD_SYSTEM(system_inp(posit:len_inp),
     & len_inp - posit + 1,key_words,key,posit)
c      PRINT*, "position changing", posit,"key=",key,
c     & "buffer=",system_inp(posit:posit)	 
      SELECT CASE(key)
      CASE(1)
      buffer = system_inp(posit:posit)	  
      IF(buffer.ne."""") STOP "ERROR:LABEL NOT FOUND"
      posit = posit + 1
      lbl_bgn = posit	  
      buffer = system_inp(posit:posit)	  
      DO WHILE(buffer.ne."""" .and. posit.lt. len_inp)
      posit = posit + 1	
      buffer = system_inp(posit:posit)	  
      ENDDO
      IF(posit.ge.len_inp) STOP
     & "THE SECOND BRACKET LABEL MISSING OR PUT IN THE END" 
      lbl_end = posit-1
      ALLOCATE(CHARACTER(LEN=lbl_end-lbl_bgn+1) :: label)
      label = system_inp(lbl_bgn:lbl_end)
      posit = posit + 1
      buffer = system_inp(posit:posit)
      IF(buffer.ne.",") STOP "COMA IS MISSING"
      posit = posit	+ 1  
      CASE(2)
      CALL REAL_NUMBER_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,mass_red)
      IF(mass_red.le.0d0) STOP "NON-POSITVE REDUCED MASS"
      IF(R_max_dist*mass_red*R_min_dist.ge.0d0)
     &	scatter_param_defined = .TRUE.	 	  
      CASE(3)
      CALL REAL_NUMBER_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,R_min_dist)
c      PRINT*, R_min_dist	 
      IF(R_min_dist.lt.0d0) STOP "NON-POSITIVE RMIN"
      IF(R_max_dist*mass_red*R_min_dist.ge.0d0)
     &	scatter_param_defined = .TRUE.	  
      CASE(4)
      CALL REAL_NUMBER_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,R_max_dist)
      IF(R_max_dist.le.0d0) STOP "NON-POSITIVE RMIN"	 
      IF(R_max_dist*mass_red*R_min_dist.ge.0d0)
     &	scatter_param_defined = .TRUE.	    
      CASE(5)
      IF(system_inp(posit:posit+2).ne."RK4" .and.
     & system_inp(posit:posit+5).ne."ODEINT") CALL ERROR_SIGNALING(5,2)	  
      IF(system_inp(posit:posit+2).eq."RK4") THEN
      rk4_defined = .TRUE.
      posit = posit + 4
      ENDIF
      IF(system_inp(posit:posit+5).eq."ODEINT") THEN
      odeint_defined = .TRUE.
	  rk4_defined=.FALSE.
      posit = posit + 7
      ENDIF	  
	  
      CASE(6)
      CALL INT_NUMBERS_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,J_tot_min,1)
      IF(J_tot_min.lt.0) STOP "NEGATIVE J total max"	 
      CASE(7)
      CALL INT_NUMBERS_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,J_tot_max,1)
      IF(J_tot_max.lt.0) STOP "NEGATIVE J total min"	  
      CASE(8)
      CALL INT_NUMBERS_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,nmbr_of_enrgs,1)
      IF(nmbr_of_enrgs.le.0) 
     & STOP "NUMBER OF COLL ENERGIES IS NON-POSITIVE"
      ALLOCATE(U(nmbr_of_enrgs))	 
      CASE(9)
      IF(nmbr_of_enrgs.lt.1) 
     & STOP "SPECIFY FIRST NUMBER OF COLL ENERGIES"
      DO i = 1,	nmbr_of_enrgs 
      CALL REAL_NUMBER_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,U(i))
      ENDDO
      ener_man_defined = .TRUE.	  
      CASE(10)
      CALL REAL_NUMBER_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,min_t_stp)
      IF(min_t_stp.lt. 0d0) STOP "NON-POSITIVE MINIMUM TIME_STEP"	 
      CASE(11)
      CALL REAL_NUMBER_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,time_step)
      IF(time_step.lt. 0d0) STOP "NON-POSITIVE TIME_STEP"		  
      CASE(12)
      CALL REAL_E_FMT_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,eps_odeint)
      IF(eps_odeint.le.0d0) STOP "NON_POSITIVE ERROR"	 
      odeint_defined = .TRUE.
      CASE(13)
      CALL REAL_NUMBER_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,mnt_crl_intgrt_err)	  
      CASE(14)
      IF(system_inp(posit:posit+2).eq."YES")
     & deflect_fun_defined = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & deflect_fun_defined = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(14,2)	
      IF(deflect_fun_defined) posit = posit + 4	 
      IF(.not.deflect_fun_defined) posit = posit + 3	  
      CASE(15)
      CALL INT_NUMBERS_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,numb_orb_prds,1)	 
      CASE(16)
      CALL REAL_NUMBER_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,U_min)
      ener_auto_defined = .TRUE.
      IF(U_max*U_min*dU.ne. 0d0) 
     & STOP "TOO MANY PARAMETERS DEFINED FOR U"	
      CASE(17)
      CALL REAL_NUMBER_READING(system_inp(posit:len_inp),	  
     & len_inp-posit+1,
     & posit,U_max)
      ener_auto_defined = .TRUE.
      IF(U_max*U_min*dU.ne. 0d0)
     & STOP "TOO MANY PARAMETERS DEFINED FOR U"	  
      CASE(18)
      CALL REAL_NUMBER_READING(system_inp(posit:len_inp),	  
     & len_inp-posit+1,
     & posit,dU)
c      PRINT*,"dU=",dU	 
      CASE(19)
      CALL REAL_E_FMT_READING(system_inp(posit:len_inp),	  
     & len_inp-posit+1,
     & posit,time_lim)
c      PRINT*,time_lim	 
      tm_lim_defined = .TRUE.	 
      CASE(20)
      IF(system_inp(posit:posit+2).eq."YES")
     & orbit_traj_defined = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & orbit_traj_defined = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(20,2)
      IF(orbit_traj_defined) posit = posit + 4	 
      IF(.not.orbit_traj_defined) posit = posit + 3
      CASE(21)
      CALL INT_NUMBERS_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,prn_l_trj,1)
      IF(prn_l_trj.ge.0) prn_l_defined = .TRUE.
      CASE(22)
      CALL REAL_NUMBER_READING(system_inp(posit:len_inp),	  
     & len_inp-posit+1,
     & posit,b_impact_parameter)
      IF(b_impact_parameter.gt.0d0) b_impact_defined = .TRUE.
      CASE(23)
      CALL INT_NUMBERS_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,nmbr_of_traj,1)
      IF(nmbr_of_traj.le.0) STOP
     & "ERROR: NUMBER OF TRAJECTOREIS MUST BE > 0"
      CASE(24)
      CALL INT_NUMBERS_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,mpi_task_per_traject,1)
      IF(mpi_task_per_traject.lt.0) STOP
     & "MPI_TASK IS NON_POSITIVE"
      mpi_task_defined = .TRUE.
      CASE(25)
      IF(system_inp(posit:posit+2).eq."YES")
     & check_point_defined = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & check_point_defined = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(25,2)
      IF(check_point_defined) posit = posit + 4	 
      IF(.not.check_point_defined) posit = posit + 3
      CASE(26)
      CALL INT_NUMBERS_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,TIME_MIN_CHECK,1)
      IF(TIME_MIN_CHECK.le.0) STOP "SUPPLY POSITIVE MAX. TIME"
      write_check_file_defined = .TRUE.
      CASE(27)
      IF(system_inp(posit:posit+2).eq."YES")
     & monte_carlo_defined = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & monte_carlo_defined = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(27,2)
      IF(monte_carlo_defined) posit = posit + 4	 
      IF(.not.monte_carlo_defined) posit = posit + 3
      CASE(28)
      IF(system_inp(posit:posit+2).eq."YES")
     & calc_elast_defined = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & calc_elast_defined = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(28,2)
      IF(calc_elast_defined) posit = posit + 4	 
      IF(.not.calc_elast_defined) posit = posit + 3
      CASE(29)
      CALL INT_NUMBERS_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,delta_l_step,1)
      IF(delta_l_step.lt.1) STOP "ERROR:SUPPLY POSITIVE DELTA L"
      CASE(30)
      CALL INT_NUMBERS_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,ang_res_num,1)
      IF(ang_res_num.lt.10)
     & STOP"ERROR:SUPPLY ANGLE RESOULTION>10"	  
      CASE(31)
      IF(system_inp(posit:posit+2).eq."YES")
     & no_orbits_defined = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & no_orbits_defined = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(28,2)
      IF(no_orbits_defined) posit = posit + 4	 
      IF(.not.no_orbits_defined) posit = posit + 3	  
      CASE(32)
      CALL INT_NUMBERS_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,j12m12_print,2) 
      CASE(33)	 
      IF(system_inp(posit:posit+2).eq."YES")
     & diff_cross_defined = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & diff_cross_defined = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(25,2)
      IF(diff_cross_defined) posit = posit + 4	 
      IF(.not.diff_cross_defined) posit = posit + 3

! Bikram Start:
      CASE(34)
      CALL INT_NUMBERS_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,numb_oscl_prds,1)	 
! Bikram End.
	  CASE(35)
      IF(system_inp(posit:posit+2).eq."YES")
     & non_format_out_defined = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & non_format_out_defined = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(25,2)
      IF(non_format_out_defined) posit = posit + 4	 
      IF(.not.non_format_out_defined) posit = posit + 3
	  CASE(36)
      IF(system_inp(posit:posit+2).eq."YES")
     & cartes_defined = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & cartes_defined = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(25,2)
      IF(cartes_defined) posit = posit + 4	 
      IF(.not.cartes_defined) posit = posit + 3	  	  
! Bikram Start Nov 2019:
	  CASE(37)
      IF(system_inp(posit:posit+2).eq."YES")
     & bikram_int = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & bikram_int = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(37,2)
      IF(bikram_int) posit = posit + 4	 
      IF(.not.bikram_int) posit = posit + 3
	  CASE(38)
      IF(system_inp(posit:posit+2).eq."YES")
     & bikram_theta = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & bikram_theta = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(37,2)
      IF(bikram_theta) posit = posit + 4	 
      IF(.not.bikram_theta) posit = posit + 3
	  CASE(39)
      IF(system_inp(posit:posit+2).eq."YES")
     & bikram_print = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & bikram_print = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(37,2)
      IF(bikram_print) posit = posit + 4	 
      IF(.not.bikram_print) posit = posit + 3
	  CASE(40)
      IF(system_inp(posit:posit+2).eq."YES")
     & bk_nrg_err = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & bk_nrg_err = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(37,2)
      IF(bk_nrg_err) posit = posit + 4	 
      IF(.not.bk_nrg_err) posit = posit + 3
	  CASE(41)
      IF(system_inp(posit:posit+2).eq."YES")
     & bk_step_size = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & bk_step_size = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(37,2)
      IF(bk_step_size) posit = posit + 4	 
      IF(.not.bk_step_size) posit = posit + 3  
      CASE(42)
      CALL INT_NUMBERS_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,numb_rk4_stps_adia,1)	  
      CASE(43)
      CALL REAL_NUMBER_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,bk_rk4_tol_adia)
      IF(bk_rk4_tol_adia.lt. 0d0) then
	  write(*,'(a)') "NEGATIVE TOLERANCE FOR RK4, ADIABATIC"	 
	  STOP 
	  endif
	  CASE(44)
      IF(system_inp(posit:posit+2).eq."YES")
     & bk_prob_interpolation = .TRUE.
      IF(system_inp(posit:posit+1).eq."NO")
     & bk_prob_interpolation = .FALSE.	  
      IF((system_inp(posit:posit+2).ne."YES") .and.
     & (system_inp(posit:posit+1).ne."NO")) CALL ERROR_SIGNALING(37,2)
      IF(bk_prob_interpolation) posit = posit + 4	 
      IF(.not.bk_prob_interpolation) posit = posit + 3
      CASE(45)
      CALL INT_NUMBERS_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,bk_dl_lr,1)
      IF(bk_dl_lr.lt.1) STOP "ERROR:SUPPLY POSITIVE DELTA L LONG-RANGE"
      CASE(46)
      CALL REAL_NUMBER_READING(system_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,bk_b_switch)
      IF(bk_b_switch.le. 0d0) then
	  write(*,'(a)') "NEGATIVE OR ZERO B_SWITCH, PLEASE,
     & PROVIDE CORRECTLY"	 
	  STOP 
	  endif
! Bikram End.
      END SELECT
! Bikram Oct'18 Start:
      if(b_impact_parameter.gt.R_max_dist) then
	  b_impact_parameter=R_max_dist
      PRINT*,
     & "WARNING: MAXIMUM IMPACT PARAMETER SHOULD BE SMALLER THAN RMAX"
	  PRINT*,
     & "MAXIMUM IMPACT PARAMETER WAS SET TO RMAX NOW" 
! Bikram End.	  
	  endif
! Bikram End
      ENDDO 	 

! Bikram Start: Feb 2021, related to L_switch
      if(bk_b_switch.gt.b_impact_parameter) then
	  bk_b_switch=b_impact_parameter
      PRINT*,
     & "WARNING: B_SWICTH SHOULD BE SMALLER THAN MAXIMUM IMPACT
     & PARAMETER"
	  PRINT*,
     & "B_SWICTH WAS SET TO MAXIMUM IMPACT PARAMETER" 
	  endif	 
	  if(bk_b_switch.gt.0.d0) then
	  if(bk_dl_lr.lt.delta_l_step) then
	  write(*,'(a,i0,1x,i0)')
     & "WARNING: LONG-RANGE DELTA_L SHOULD BE LARGER THAN
     & SHORT-RANGE DELTA_L: ", delta_l_step, bk_dl_lr
	  bk_dl_lr = delta_l_step
	  PRINT*,
     & "LONG-RANGE DELTA_L WAS SET TO SHORT-RANGE DELTA_L" 
	  endif	 
	  end if
	  
! Bikram End.
      END SUBROUTINE SYSTEM_PARSING
      SUBROUTINE KEY_WORD_SYSTEM(inp,length,key_word_used,key,place)
      IMPLICIT NONE !!! KEY WORDS FOR SYSTEM. SEE MANUAL
      INTEGER, PARAMETER :: num_key_word = 46 	  
      INTEGER length,posit,key_word_used(num_key_word),
     & key,place,decrement,i
      CHARACTER(LEN=length) inp
      CHARACTER(LEN=6) :: label_word="LABEL="                   			!CASE(1)
      CHARACTER(LEN=9) :: mass_word="MASS_RED="                 			!CASE(2)
      CHARACTER(LEN=5) :: rmin_word="RMIN="                     			!CASE(3)
      CHARACTER(LEN=5) :: rmax_word="RMAX="                     			!CASE(4)
      CHARACTER(LEN=11) :: prop_word="PROPAGATOR="              			!CASE(5)
      CHARACTER(LEN=6) :: jtotl_word="JTOTL="                   			!CASE(6)
      CHARACTER(LEN=6) :: jtotu_word="JTOTU="                   			!CASE(7)
      CHARACTER(LEN=11) :: nb_ener="NMB_ENERGS="                			!CASE(8)
      CHARACTER(LEN=9) :: u_ener_word="U_ENERGY="               			!CASE(9) 
      CHARACTER(LEN=10) :: min_tm_stp="MIN_TMSTP="              			!CASE(10)
      CHARACTER(LEN=10) :: tm_stp="TIME_STEP="                  			!CASE(11)
      CHARACTER(LEN=11) :: eps_word="EPS_ODEINT="               			!CASE(12)
      CHARACTER(LEN=11) :: pres_cross_sec_word="EPS_MONCAR="    			!CASE(13)
      CHARACTER(LEN=10) :: dflc_word="DFLC_FNCT="               			!CASE(14)
      CHARACTER(LEN=10) :: num_of_period_word="NMB_LOOPS="      			!CASE(15)
      CHARACTER(LEN=5) :: u_min_word="UMIN="                    			!CASE(16)
      CHARACTER(LEN=5) :: u_max_word="UMAX="                    			!CASE(17)
      CHARACTER(LEN=3) :: du_word="DU="						    			!CASE(18)
      CHARACTER(LEN=9) :: time_lim_word="TIME_LIM="			    			!CASE(19)
      CHARACTER(LEN=11) :: orbit_word="SHOW_ORBIT="							!CASE(20)
      CHARACTER(LEN=10) :: max_nmb_tr_word="PRN_TRJCT="		    			!CASE(21)
      CHARACTER(LEN=8) :: b_impact_word="B_IMPCT="		        			!CASE(22)
      CHARACTER(LEN=9) :: n_taject_word="NMB_TRAJ="		        			!CASE(23)
      CHARACTER(LEN=12):: mpi_proc_word="MPI_PERTRAJ="          			!CASE(24)
      CHARACTER(LEN=8):: continue_word="RESTART="               			!CASE(25)
      CHARACTER(LEN=12):: time_lim_check_file_word="CHECK_POINT="           !CASE(26)
      CHARACTER(LEN=12):: monte_carlo_word="MONTE_CARLO="       			!CASE(27)
      CHARACTER(LEN=13):: calc_elast_word="CALC_ELASTIC="       			!CASE(28)
      CHARACTER(LEN=3)  :: delta_l_step_word  = "DL="	        			!CASE(29)
      CHARACTER(LEN=8) :: ang_res_word = "ANG_RES="             			!CASE(30)
      CHARACTER(LEN=13)	:: no_orbits_word = "NO_RESONANCE="     			!CASE(31)
      CHARACTER(LEN=7)	:: prn_j12m12_word = "PRN_JM="          			!CASE(32)	 
      CHARACTER(LEN=11)	:: diff_cross_word = "DIFF_CROSS="      			!CASE(33)
	  CHARACTER(LEN=10) :: num_of_oscl_word="NMB_OSCIL="        			!CASE(34)		!Bikram
	  CHARACTER(LEN=15) :: non_format_out_word="NONFORM_OUTPUT="			!CASE(35)
	  CHARACTER(LEN=7) :: cartes_word="CARTES="								!CASE(36)		  
      CHARACTER(LEN=12):: bikram_int_word="SINGLE_STEP="        			!CASE(37)		!Bikram
      CHARACTER(LEN=12):: bikram_theta_word="POLAR_ANGLE="       			!CASE(38)		!Bikram
      CHARACTER(LEN=9):: bikram_traj_print="BK_PRINT="       				!CASE(39)		!Bikram
      CHARACTER(LEN=13):: bikram_nrg_err="ENERGY_ERROR="        			!CASE(40)		!Bikram
      CHARACTER(LEN=13):: bikram_step_size="AT_ADAPT_STP="      			!CASE(41)		!Bikram
      CHARACTER(LEN=10) :: num_of_stps_rk4_adia="NMB_STEPS="    			!CASE(42)		!Bikram
      CHARACTER(LEN=11) :: rk4_tol_adia="AT_ADAPTOL="               		!CASE(43)		!Bikram
      CHARACTER(LEN=12):: bikram_prob_interpolation="PROB_SPLINE="       	!CASE(44)		!Bikram
      CHARACTER(LEN=6) :: bikram_delta_l_LR="DL_LR="                 		!CASE(45)		!Bikram
      CHARACTER(LEN=9) :: bikram_b_switch="B_SWITCH="                 		!CASE(46)		!Bikram
      CHARACTER(LEN=1) buffer
      LOGICAL key_used
      IF(length.le.0) RETURN
      key_used = .FALSE.	  
      posit = 1
      decrement = place - posit
      buffer = 	inp(posit:posit)  
      DO WHILE(buffer.ne."=" .and. posit.lt.length)
      posit = posit + 1
      buffer = 	inp(posit:posit)
c      PRINT*, "WANNA BELIEVE",buffer	  
      ENDDO
c      PRINT*,"WANNA CHECK",posit,inp(1:posit)	  
      IF(inp(1:posit).eq.label_word) THEN
      key = 1
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.mass_word) THEN
      key = 2
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.rmin_word) THEN
      key = 3
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.rmax_word) THEN
      key = 4
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.prop_word) THEN
      key = 5
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.jtotl_word) THEN
      key = 6
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.jtotu_word) THEN
      key = 7
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.nb_ener) THEN
      key = 8
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.u_ener_word) THEN
      key = 9
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.min_tm_stp) THEN
      key = 10
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.tm_stp) THEN
      key = 11
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.eps_word) THEN
      key = 12
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	
      IF(inp(1:posit).eq.pres_cross_sec_word) THEN
      key = 13
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.dflc_word) THEN
      key = 14
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.num_of_period_word) THEN
      key = 15
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.u_min_word) THEN
      key = 16
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.u_max_word) THEN
      key = 17
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.du_word) THEN
      key = 18
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.time_lim_word) THEN
      key = 19
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.orbit_word) THEN
      key = 20
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.max_nmb_tr_word) THEN
      key = 21
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.b_impact_word) THEN
      key = 22
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.n_taject_word) THEN
      key = 23
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.mpi_proc_word ) THEN
      key = 24
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.continue_word ) THEN
      key = 25
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.time_lim_check_file_word ) THEN
      key = 26
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF		  
      IF(inp(1:posit).eq.monte_carlo_word) THEN
      key = 27
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.calc_elast_word) THEN
      key = 28
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.delta_l_step_word) THEN
      key = 29
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.ang_res_word) THEN
      key = 30
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.no_orbits_word) THEN
      key = 31
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.prn_j12m12_word) THEN
      key = 32
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.diff_cross_word) THEN
      key = 33
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF

! Bikram Start:
      IF(inp(1:posit).eq.num_of_oscl_word) THEN
      key = 34
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
! Bikram End.
      IF(inp(1:posit).eq.non_format_out_word) THEN
      key = 35
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	
      IF(inp(1:posit).eq.cartes_word) THEN
      key = 36
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
! Bikram Start Nov 2019:
	  IF(inp(1:posit).eq.bikram_int_word) THEN
      key = 37
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
	  IF(inp(1:posit).eq.bikram_theta_word) THEN
      key = 38
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
	  IF(inp(1:posit).eq.bikram_traj_print) THEN
      key = 39
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
	  IF(inp(1:posit).eq.bikram_nrg_err) THEN
      key = 40
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
	  IF(inp(1:posit).eq.bikram_step_size) THEN
      key = 41
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.num_of_stps_rk4_adia) THEN
      key = 42
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.rk4_tol_adia) THEN
      key = 43
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
	  IF(inp(1:posit).eq.bikram_prob_interpolation) THEN
      key = 44
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.bikram_delta_l_LR) THEN
      key = 45
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.bikram_b_switch) THEN
      key = 46
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
! Bikram End.
	  
      IF(.not.key_used) THEN
      PRINT*,inp(1:posit)	  
      STOP 
     & "ERROR IN SYSTEM: WORD NOT FOUND OR INPUT ABOVE IS INCORRECT"
      ENDIF	 
      place = posit + decrement + 1	  
      END SUBROUTINE KEY_WORD_SYSTEM
      SUBROUTINE REAL_E_FMT_READING(inp,len_inp,posit,e_numb)
      USE ERRORS!!! READING FORMATED REAL NUMBER
      USE MPI_DATA	  
      IMPLICIT NONE
      REAL*8 e_numb
      INTEGER len_inp
      CHARACTER(LEN = len_inp) inp
      INTEGER posit,place, end_pos,bgn_pos,decrement
      INTEGER pos_e,pos_sign,power,val_pow
      CHARACTER(LEN =1) buffer
      LOGICAL dot_num
      LOGICAL e_num
      LOGICAL sign_num,minus_def
      LOGICAL ERROR_NUM	  
      IF(len_inp.le.0) THEN 
	  show_error = "ERROR: NON DATA ON INPUT WHILE READING E. FORMAT"
      IF(MYID.eq.0) THEN
      PRINT*,show_error
      PRINT*, "CHECK YOUR INPUT"
      ENDIF	  
      STOP	  
      ENDIF	  
      place = 1
      power = 1	  
      buffer = inp(place:place)
      dot_num=.FALSE.
      e_num = .FALSE.
      sign_num = .FALSE.
      bgn_pos = 1
      decrement = posit - bgn_pos
      DO WHILE(buffer.ne.',' .and. place.le.len_inp)
      IF(place.eq.bgn_pos) THEN
      IF(buffer.eq.'-') THEN
      minus_def = .TRUE.
      place = place + 1
      buffer = inp(place:place) 	  
      CYCLE	  
      ENDIF	  
      ENDIF	  
      IF((ICHAR(buffer).gt.57 .or. ICHAR(buffer).lt.48)
     & .and. buffer.ne.'.'  .and. buffer.ne.'E' .and.
     &	 buffer.ne.'+' .and. buffer.ne.'-')
     &	 STOP "NONREAL_NUMBER_FOR_E_REAL_FMT"
      IF(buffer.eq.'.') THEN
      IF(.not.dot_num) THEN	  
      dot_num = .TRUE.
      ELSE
      STOP "TOO MANY DOTS"	  
      ENDIF
      ENDIF
      IF(buffer.eq.'E') THEN
      pos_e = place
      IF(.not.e_num) THEN	  
      e_num = .TRUE.
      ELSE
      STOP "TOO MANY E'S"	  
      ENDIF
      ENDIF	  	  
      IF(buffer.eq.'+' .or. buffer.eq.'-' ) THEN
      pos_sign = place
      IF(buffer.eq.'-')	power = -1  
      IF((pos_sign-pos_e).ne.1) STOP "WRONG_E_SIGN_POS"	  
      IF(.not.sign_num) THEN	  
      sign_num = .TRUE.
      ELSE
      STOP "TOO MANY SIGNS"	  
      ENDIF
      ENDIF	  
      place = place + 1
      buffer = inp(place:place)	  
      ENDDO
      end_pos= place - 1
      ERROR_NUM = dot_num.and.e_num.and.sign_num	  
      IF(.not.ERROR_NUM) PRINT *, " NOT REAL NUM IN ",place+decrement
      IF(.not.ERROR_NUM) STOP
      IF(bgn_pos.gt.pos_e-1) STOP "WRONG E FMT INP"
      IF(end_pos.lt.pos_sign+1) STOP "WRONG SING FMT INP"	  
      READ(inp(bgn_pos:pos_e-1),*) e_numb
      READ(inp(pos_sign+1:end_pos),*) val_pow
      e_numb = e_numb*10d0**(power*val_pow)	  
      posit = end_pos + decrement+2	 	  
      END SUBROUTINE REAL_E_FMT_READING
      SUBROUTINE POTENTIAL_PARSING(potential_inp,len_inp)
      USE VARIABLES	  !!!! READING POTENTIAL BLOCK
      IMPLICIT NONE
      INTEGER len_inp
      CHARACTER(LEN = len_inp) potential_inp
      INTEGER, PARAMETER :: key_num = 57
      INTEGER posit,key_words(key_num)
      INTEGER key
      INTEGER i
      CHARACTER(LEN=1) buffer
	  logical bk_exst
      grid_defined = .FALSE.	  
      expansion_defined = .FALSE.
      matrix_reading_defined = .FALSE.
      run_prog_defined = .TRUE.  
      posit = 1
      key_words = 0
      r_unit = 1.0
      au_unit_r = .TRUE.
      angs_unit	= .FALSE.
      au_unit_e = .TRUE.
      cm_unit = .FALSE.
      klv_unit = .FALSE.
      kal_unit = .FALSE.
      atom_coord_dist = .FALSE.
      make_grid_file = .FALSE.
      grid_file_found = .FALSE.
      calc_expansion_defined = .FALSE.	
      read_r_grid_from_file_defined = .FALSE.
      print_matrix_defined = .FALSE.
      pot_expansion_defined = .FALSE.
      unformat_defined = .TRUE.
      test_expansion_defined = .FALSE.
      terms_onfly_defined = .FALSE.
      terms_file_defined = .FALSE.
      eq_grid_defined = .TRUE.	
      calc_matrix_defined = .TRUE.	  
	  bikram_mij_shift = .TRUE.																!Bikram
	  bikram_mij_multiprint = .FALSE.														!Bikram
	  bikram_eq_grd_gamma = .TRUE.															!Bikram
	  bikram_identical_pes = .FALSE.														!Bikram
	  rms_defined = .FALSE.																	!Bikram
	  bikram_ident_terms = .FALSE.															!Bikram
	  bikram_rebalance = .FALSE.															!Bikram
	  bikram_rebalance_comp = .FALSE.														!Bikram
	  bikram_w3j_fact = .FALSE.																!Bikram
	  bikram_truncate_MTRX = .FALSE.														!Bikram
	  bikram_on = .FALSE.																	!Bikram
	  bikram_mtrx_path = .FALSE.															!Bikram
      e_unit = 1.0
      n_r_coll = 0
      nterms = 0
      L_MAX_EXPAN = -1
      NJU_MAX_EXPAN = -1
      L1_MAX_EXPAN = -1
      L2_MAX_EXPAN = -1
      L1_MIN_EXPAN = 0
      L2_MIN_EXPAN = 0	  
      NJU1_MAX_EXPAN = -1
      NJU2_MAX_EXPAN = -1
      ir_bgn_exp_pnt = -1
      ir_fin_exp_pnt = -1
      bikram_rms_ang1 = 20																	!Bikram
      bikram_rms_ang2 = 10																	!Bikram
      bikram_rms_ang3 = 20																	!Bikram
      bikram_axl_sym1 = 1																	!Bikram
      bikram_axl_sym2 = 1																	!Bikram
      bikram_equ_sym1 = 1																	!Bikram
      bikram_equ_sym2 = 1																	!Bikram
      m_elastic_proj_print = 0	  
	  MIJ_ZERO_CUT = 1.d-12																	!Bikram
	  bikram_cutoff_r1 = 6.0d0																!Bikram
	  bikram_cutoff_r2 = 8.0d0																!Bikram
	  bikram_rms_r = 8.0d0																	!Bikram
      DO WHILE(posit.lt.len_inp)
      CALL KEY_WORD_POTENTIAL(potential_inp(posit:len_inp),
     & len_inp - posit + 1,key_words,key,posit)
c      PRINT*, "position changing", posit,"key=",key,
c     & "buffer=",potential_inp(posit:posit)	 
      SELECT CASE(key)
      CASE(1)
      IF(potential_inp(posit:posit+3).eq."YES,")
     & expansion_defined = .TRUE.
      IF(potential_inp(posit:posit+2).eq."NO,")
     & expansion_defined = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO")) STOP "USER_DEF_ERROR"
      IF(expansion_defined) posit = posit + 4	 
      IF(.not.expansion_defined) posit = posit + 3
      CASE(2)
      IF(potential_inp(posit:posit+3).eq."A.U.")
     & au_unit_e = .TRUE.
      IF(potential_inp(posit:posit+3).eq."CM-1") THEN
      cm_unit = .TRUE.
      au_unit_e = .FALSE. 	  
      e_unit  =  27.21113845d0*8065.5448d0       	  
      ENDIF
      IF(potential_inp(posit:posit+3).eq."KLVN") THEN
      klv_unit = .TRUE.
      au_unit_e = .FALSE. 	  
      ENDIF
      IF(potential_inp(posit:posit+3).eq."KCAL") THEN
      kal_unit = .TRUE.
      au_unit_e = .FALSE. 	  
      ENDIF		  
      IF((potential_inp(posit:posit+3).ne."A.U.") .and.
     & (potential_inp(posit:posit+3).ne."CM-1").and.
     & (potential_inp(posit:posit+3).ne."KLVN")
     & .and. (potential_inp(posit:posit+3).ne."KCAL"))
     &	STOP "USER_DEF_ERROR"
      posit = posit + 5	 
      CASE(3)
      IF(potential_inp(posit:posit+3).eq."A.U.")
     & au_unit_r = .TRUE.
      IF(potential_inp(posit:posit+3).eq."ANGS") THEN
      angs_unit = .TRUE.
      au_unit_r = .FALSE. 	  
      r_unit  =  1d0/0.52917721067d0     	  
      ENDIF	  
      IF((potential_inp(posit:posit+3).ne."A.U.") .and.
     & (potential_inp(posit:posit+3).ne."ANGS"))
     & CALL ERROR_SIGNALING(3,3)
      posit = posit + 5
      CASE(4)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,nterms,1)
!      IF(.not.expansion_defined) STOP "DONT SPECIFY TERMS NUMBER"
      CASE(5)
      terms_defined = .TRUE.
      IF(nterms.le.0) STOP "NON-POSITIVE NUMBER OF TERMS"	  
      SELECT CASE(system_type)
      CASE(1)
      IF(.not.fine_structure_defined .or. LORB_FINE.eq.0) THEN	  
      ALLOCATE(L_TERM(nterms))
      L_TERM = 0
      DO i=1,nterms	  
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,L_TERM(i),1)
      ENDDO	 
      ELSE
      ALLOCATE(L_TERM(nterms))
      ALLOCATE(M_TERM(nterms))	  
      L_TERM = 0
      M_TERM = 0	  
      DO i=1,nterms	  
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,LM,2)
      L_TERM(i) = LM(1)
      M_TERM(i) = LM(2)	  
      ENDDO 	  
      ENDIF	  
      CASE(2)
      ALLOCATE(L_TERM(nterms))
      L_TERM = 0
      DO i=1,nterms	  
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,L_TERM(i),1)
      ENDDO	 
      CASE(3)
      ALLOCATE(L_TERM(nterms))
      ALLOCATE(M_TERM(nterms))	  
      L_TERM = 0
      M_TERM = 0	  
      DO i=1,nterms	  
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,LM,2)
      L_TERM(i) = LM(1)
      M_TERM(i) = LM(2)	  
      ENDDO
      CASE(4)
      ALLOCATE(L_TERM(nterms))
      ALLOCATE(M_TERM(nterms))	  
      L_TERM = 0
      M_TERM = 0	  
      DO i=1,nterms	  
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,LM,2)
      L_TERM(i) = LM(1)
      M_TERM(i) = LM(2)	  
      ENDDO
      CASE(5)
      ALLOCATE(L1_TERM(nterms))
      ALLOCATE(L2_TERM(nterms))
      ALLOCATE(L_TERM(nterms))	  
      L1_TERM = 0
      L2_TERM = 0
      L_TERM = 0	  
      DO i=1,nterms	  
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,L1L2L,3)
      L1_TERM(i) = L1L2L(1)
      L2_TERM(i) = L1L2L(2)
      L_TERM(i) = L1L2L(3)	  
      ENDDO
      CASE(6)
      ALLOCATE(L1_TERM(nterms))
      ALLOCATE(L2_TERM(nterms))
      ALLOCATE(L_TERM(nterms))	  
      L1_TERM = 0
      L2_TERM = 0
      L_TERM = 0	  
      DO i=1,nterms	  
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,L1L2L,3)
      L1_TERM(i) = L1L2L(1)
      L2_TERM(i) = L1L2L(2)
      L_TERM(i) = L1L2L(3)	  
      ENDDO
      CASE(7)
      ALLOCATE(L1_TERM(nterms))
      ALLOCATE(NJU1_TERM(nterms))	  
      ALLOCATE(L2_TERM(nterms))
      ALLOCATE(L_TERM(nterms))	  
      L1_TERM = 0
      NJU1_TERM = 0	  
      L2_TERM = 0
      L_TERM = 0	  
      DO i=1,nterms	  
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,L1NJU1L2L,4)
      L1_TERM(i) = L1NJU1L2L(1)
      NJU1_TERM(i) = L1NJU1L2L(2)     	  
      L2_TERM(i) = L1NJU1L2L(3)
      L_TERM(i) = L1NJU1L2L(4)
      ENDDO	  
      CASE(8)
      ALLOCATE(L1_TERM(nterms))
      ALLOCATE(NJU1_TERM(nterms))	  
      ALLOCATE(L2_TERM(nterms))
      ALLOCATE(L_TERM(nterms))	  
      L1_TERM = 0
      NJU1_TERM = 0	  
      L2_TERM = 0
      L_TERM = 0	  
      DO i=1,nterms	  
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,L1NJU1L2L,4)
      L1_TERM(i) = L1NJU1L2L(1)
      NJU1_TERM(i) = L1NJU1L2L(2)     	  
      L2_TERM(i) = L1NJU1L2L(3)
      L_TERM(i) = L1NJU1L2L(4)
      ENDDO
      CASE(9)
      ALLOCATE(L1_TERM(nterms))
      ALLOCATE(L2_TERM(nterms))
      ALLOCATE(L_TERM(nterms))
      ALLOCATE(NJU1_TERM(nterms))
      ALLOCATE(NJU2_TERM(nterms))	  
      L1_TERM = 0
      L2_TERM = 0
      L_TERM = 0
      NJU1_TERM = 0
      NJU2_TERM = 0	  
      DO i=1,nterms	  
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,L1L2N1N2L,5)
      L1_TERM(i) = L1L2N1N2L(1)
      L2_TERM(i) = L1L2N1N2L(3)
      NJU1_TERM(i) = L1L2N1N2L(2) 
      NJU2_TERM(i) = L1L2N1N2L(4) 	  
      L_TERM(i) = L1L2N1N2L(5)	  
      ENDDO	  	  
      CASE(0)
      ALLOCATE(L1_TERM(nterms))
      ALLOCATE(L2_TERM(nterms))
      ALLOCATE(L_TERM(nterms))
      ALLOCATE(NJU1_TERM(nterms))
      ALLOCATE(NJU2_TERM(nterms))	  
      L1_TERM = 0
      L2_TERM = 0
      L_TERM = 0
      NJU1_TERM = 0
      NJU2_TERM = 0	  
      DO i=1,nterms	  
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,L1L2N1N2L,5)
      L1_TERM(i) = L1L2N1N2L(1)
      L2_TERM(i) = L1L2N1N2L(3)
      NJU1_TERM(i) = L1L2N1N2L(2) 
      NJU2_TERM(i) = L1L2N1N2L(4) 	  
      L_TERM(i) = L1L2N1N2L(5)	  
      ENDDO	  
      END SELECT	  
      CASE(6)	  
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,n_r_coll,1)
      IF(n_r_coll.le.0) STOP "NON_POSITIVE NUMBER OF R-GRID PNTS"
      CASE(7)
      IF(system_type.eq.2) THEN	  
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,n_r_vib,1)
      grid_defined = .TRUE.
c      PRINT*,n_r_vib,grid_defined	!!!!!!!!!! DELETE  
      ENDIF
      IF(system_type.eq.5 .or. system_type.eq.6 ) THEN	  
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,n_2_pnts,2)
      n_r_vib1 = n_2_pnts(1)
      n_r_vib2 = n_2_pnts(2)
      grid_defined = .TRUE.	  
      ENDIF
      IF(.not.grid_defined) STOP "CHECK COLLISION TYPE"  	  
      CASE(8)
      IF(system_type.le.4 .and. system_type.ge.1) THEN
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,n_alpha,1) 
      grid_defined = .TRUE. 	 
      ENDIF
      IF(system_type.ge.5 .or. system_type.eq.0) THEN
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,n_2_pnts,2)
      n_alpha1 = n_2_pnts(1)
      n_alpha2 = n_2_pnts(2)
      grid_defined = .TRUE. 	  
      ENDIF
      IF(.not.grid_defined) STOP "CHECK COLLISION TYPE" 	  
      CASE(9)
      IF(system_type.le.4 .and. system_type.ge.1) THEN
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,n_beta,1)
      grid_defined = .TRUE.	 
      ENDIF
      IF(system_type.ge.5 .or. system_type.eq.0) THEN
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,n_2_pnts,2)
      n_beta1 = n_2_pnts(1)
      n_beta2 = n_2_pnts(2)
      grid_defined = .TRUE.	  
      ENDIF
      IF(.not.grid_defined) STOP "CHECK COLLISION TYPE" 	  
      CASE(10)
      IF(system_type.le.4 .and. system_type.ge.3) THEN
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,n_gamma,1)
      grid_defined = .TRUE.
      ENDIF	  
      IF(system_type.ge.5 .or. system_type.eq.0) THEN
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,n_2_pnts,2)
      n_gamma1 = n_2_pnts(1)
      n_gamma2 = n_2_pnts(2)
      grid_defined=.TRUE.	  
      ENDIF
      IF(.not.grid_defined) STOP "CHECK COLLISION TYPE"
      CASE(11)
      IF(potential_inp(posit:posit+2).eq."YES")
     & matrix_reading_defined = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & matrix_reading_defined = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &	 CALL ERROR_SIGNALING(11,3)
      IF(matrix_reading_defined) posit = posit + 4	 
      IF(.not.matrix_reading_defined) posit = posit + 3
      CASE(12)
      IF(potential_inp(posit:posit+2).eq."YES")
     & run_prog_defined = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & run_prog_defined = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO")) 
     & CALL ERROR_SIGNALING(12,3)
      IF(run_prog_defined) posit = posit + 4	 
      IF(.not.run_prog_defined) posit = posit + 3	  
      CASE(13)
      IF(potential_inp(posit:posit+2).eq."YES")
     & atom_coord_dist = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & atom_coord_dist = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     & CALL ERROR_SIGNALING(13,3)
      IF(atom_coord_dist) posit = posit + 4	 
      IF(.not.atom_coord_dist) posit = posit + 3
      CASE(14)
      IF(potential_inp(posit:posit+2).eq."YES")
     & make_grid_file = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & make_grid_file = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     & CALL ERROR_SIGNALING(14,3)
      IF(make_grid_file) posit = posit + 4	 
      IF(.not.make_grid_file) posit = posit + 3
      CASE(15)
      IF(potential_inp(posit:posit+2).eq."YES")
     & calc_expansion_defined = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & calc_expansion_defined = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     & CALL ERROR_SIGNALING(15,3)
      IF(calc_expansion_defined) posit = posit + 4	 
      IF(.not.calc_expansion_defined) posit = posit + 3
      CASE(16)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,L_MAX_EXPAN,1)
      CASE(17)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,NJU_MAX_EXPAN,1)
      CASE(18)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,L1_MAX_EXPAN,1)
      CASE(19)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,L2_MAX_EXPAN,1)
      CASE(20)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,NJU1_MAX_EXPAN,1)
      CASE(21)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,NJU2_MAX_EXPAN,1)
      CASE(22)
      IF(potential_inp(posit:posit+2).eq."YES")
     & read_r_grid_from_file_defined = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & read_r_grid_from_file_defined = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     & CALL ERROR_SIGNALING(22,3)
      IF(read_r_grid_from_file_defined) posit = posit + 4	 
      IF(.not.read_r_grid_from_file_defined) posit = posit + 3	  
      CASE(23)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,L1_MIN_EXPAN,1)	  
      CASE(24)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,L2_MIN_EXPAN,1)
      CASE(25)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,ir_bgn_exp_pnt,1)
      CASE(26)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,ir_fin_exp_pnt,1)
      CASE(27)
      IF(potential_inp(posit:posit+2).eq."YES")
     & print_matrix_defined = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & print_matrix_defined = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     & CALL ERROR_SIGNALING(27,3)	 
      IF(print_matrix_defined) posit = posit + 4	 
      IF(.not.print_matrix_defined) posit = posit + 3
      CASE(28)
      IF(potential_inp(posit:posit+2).eq."YES")
     & terms_onfly_defined = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & terms_onfly_defined = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(27,3)	 
      IF(terms_onfly_defined) posit = posit + 4	 
      IF(.not.terms_onfly_defined) posit = posit + 3
      CASE(29)
      IF(potential_inp(posit:posit+2).eq."YES")
     & unformat_defined = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & unformat_defined = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(29,3)	
      IF(unformat_defined) posit = posit + 4	 
      IF(.not.unformat_defined) posit = posit + 3
      CASE(30)
      IF(potential_inp(posit:posit+2).eq."YES")
     & test_expansion_defined = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & test_expansion_defined = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     & CALL ERROR_SIGNALING(30,3)	
      IF(test_expansion_defined) posit = posit + 4	 
      IF(.not.test_expansion_defined) posit = posit + 3
      CASE(31)
      IF(potential_inp(posit:posit+2).eq."YES")
     & terms_file_defined = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & terms_file_defined = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     & CALL ERROR_SIGNALING(30,3)	
      IF(terms_file_defined) posit = posit + 4	 
      IF(.not.terms_file_defined) posit = posit + 3
      CASE(32)	  
      IF(potential_inp(posit:posit+2).eq."YES")
     & eq_grid_defined = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & eq_grid_defined = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     & CALL ERROR_SIGNALING(30,3)	
      IF(eq_grid_defined) posit = posit + 4	 
      IF(.not.eq_grid_defined) posit = posit + 3  
      CASE(33)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,m_elastic_proj_print,1)
      CASE(34)
      IF(potential_inp(posit:posit+2).eq."YES")
     & calc_matrix_defined = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & calc_matrix_defined = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(29,3)	
      IF(calc_matrix_defined) posit = posit + 4	 
      IF(.not.calc_matrix_defined) posit = posit + 3	 
      CASE(35)
      IF(potential_inp(posit:posit+2).eq."YES")
     & bikram_mij_shift = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & bikram_mij_shift = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(29,3)	
      IF(bikram_mij_shift) posit = posit + 4	 
      IF(.not.bikram_mij_shift) posit = posit + 3	  
      CASE(36)	  
	  call REAL_E_FMT_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,MIJ_ZERO_CUT)	 
      CASE(37)
      IF(potential_inp(posit:posit+2).eq."YES")
     & bikram_mij_multiprint = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & bikram_mij_multiprint = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(29,3)	
      IF(bikram_mij_multiprint) posit = posit + 4	 
      IF(.not.bikram_mij_multiprint) posit = posit + 3	  
      CASE(38)
      IF(potential_inp(posit:posit+2).eq."YES")
     & bikram_eq_grd_gamma = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & bikram_eq_grd_gamma = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(29,3)	
      IF(bikram_eq_grd_gamma) posit = posit + 4	 
      IF(.not.bikram_eq_grd_gamma) posit = posit + 3	   
      CASE(39)
      IF(potential_inp(posit:posit+2).eq."YES")
     & bikram_on = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & bikram_on = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(29,3)	
      IF(bikram_on) posit = posit + 4	 
      IF(.not.bikram_on) posit = posit + 3	  	   
      CASE(40)
      IF(potential_inp(posit:posit+2).eq."YES")
     & bikram_mtrx_path = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & bikram_mtrx_path = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(29,3)	
      IF(bikram_mtrx_path) posit = posit + 4	 
      IF(.not.bikram_mtrx_path) posit = posit + 3
	  
!Changing matrix file path 
!so that user don't need to copy the large matrix file in each directory.
	  if(bikram_mtrx_path) then
	  inquire(file = "Matrix_Path.DAT", exist = bk_exst)
	  if(.not.bk_exst) then
	  write(*,'(a)')"File containing matrix path is not provided."
	  stop
	  
	  else
	  open(111,file="Matrix_Path.DAT",status="old",action = "read")
	  read(111,'(a)') bk_matrix_path11
	  close(111)
	  end if
	  bk_matrix_path11 = trim(bk_matrix_path11)
	  end if
      CASE(41)
      IF(potential_inp(posit:posit+2).eq."YES")
     & bikram_w3j_fact = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & bikram_w3j_fact = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(29,3)	
      IF(bikram_w3j_fact) posit = posit + 4	 
      IF(.not.bikram_w3j_fact) posit = posit + 3	  
      CASE(42)
      IF(potential_inp(posit:posit+2).eq."YES")
     & bikram_truncate_MTRX = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & bikram_truncate_MTRX = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(29,3)	
      IF(bikram_truncate_MTRX) posit = posit + 4	 
      IF(.not.bikram_truncate_MTRX) posit = posit + 3	    
      CASE(43)	  
	  call REAL_NUMBER_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,bikram_cutoff_r1)
	  if(bikram_cutoff_r1.lt.R_min_dist .or. 
     & bikram_cutoff_r1.gt.R_max_dist) then
	  print*, "CUTOFF_R1 is <Rmin or >Rmax", bikram_cutoff_r1, 
     & R_min_dist, R_max_dist
	  stop
	  end if
      CASE(44)	  
	  call REAL_NUMBER_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,bikram_cutoff_r2)	
	  if(bikram_cutoff_r2.lt.R_min_dist .or. 
     & bikram_cutoff_r2.gt.R_max_dist) then
	  print*, "CUTOFF_R2 is <Rmin or >Rmax", bikram_cutoff_r2, 
     & R_min_dist, R_max_dist
	  stop
	  end if
	  CASE(45)
      IF(potential_inp(posit:posit+2).eq."YES")
     & bikram_identical_pes = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & bikram_identical_pes = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(29,3)	
      IF(bikram_identical_pes) posit = posit + 4	 
      IF(.not.bikram_identical_pes) posit = posit + 3	   
	  CASE(46)
      IF(potential_inp(posit:posit+2).eq."YES")
     & rms_defined = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & rms_defined = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(29,3)	
      IF(rms_defined) posit = posit + 4	 
      IF(.not.rms_defined) posit = posit + 3	   
	  CASE(47)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,bikram_rms_ang1,1)
      CASE(48)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,bikram_rms_ang2,1)
      CASE(49)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,bikram_rms_ang3,1)
      CASE(50)	  
	  call REAL_NUMBER_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,bikram_rms_r)
	  if(bikram_rms_r.lt.R_min_dist .or. 
     & bikram_rms_r.gt.R_max_dist) then
	  print*, "EXP_RMS_R is <Rmin or >Rmax", bikram_rms_r, 
     & R_min_dist, R_max_dist
	  stop
	  end if
	  CASE(51)
      IF(potential_inp(posit:posit+2).eq."YES")
     & bikram_ident_terms = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & bikram_ident_terms = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(29,3)	
      IF(bikram_ident_terms) posit = posit + 4	 
      IF(.not.bikram_ident_terms) posit = posit + 3	   
      CASE(52)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,bikram_axl_sym1,1)
      CASE(53)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,bikram_axl_sym2,1)
      CASE(54)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,bikram_equ_sym1,1)
      CASE(55)
      CALL INT_NUMBERS_READING(potential_inp(posit:len_inp),
     & len_inp-posit+1,
     & posit,bikram_equ_sym2,1)
	  CASE(56)
      IF(potential_inp(posit:posit+2).eq."YES")
     & bikram_rebalance = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & bikram_rebalance = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(29,3)	
      IF(bikram_rebalance) posit = posit + 4	 
      IF(.not.bikram_rebalance) posit = posit + 3	   
	  CASE(57)
      IF(potential_inp(posit:posit+2).eq."YES")
     & bikram_rebalance_comp = .TRUE.
      IF(potential_inp(posit:posit+1).eq."NO")
     & bikram_rebalance_comp = .FALSE.	  
      IF((potential_inp(posit:posit+2).ne."YES") .and.
     & (potential_inp(posit:posit+1).ne."NO"))
     &  CALL ERROR_SIGNALING(29,3)	
      IF(bikram_rebalance_comp) posit = posit + 4	 
      IF(.not.bikram_rebalance_comp) posit = posit + 3	   

      END SELECT
      ENDDO	
	  
	  if(expansion_defined) terms_file_defined = .true.						!Bikram 
      END SUBROUTINE POTENTIAL_PARSING
      SUBROUTINE KEY_WORD_POTENTIAL(inp,length,key_word_used,key,place)
      IMPLICIT NONE 															!!! KEY WORDS FOR POTENTIAL
      INTEGER, PARAMETER :: num_key_word = 57
      INTEGER length,posit,key_word_used(num_key_word),
     & key,place,decrement,i
      CHARACTER(LEN=length) inp
      CHARACTER(LEN=15) :: expan_word="READ_EXPANSION="         				!CASE(1)
      CHARACTER(LEN=8) :: energy_units_word="E_UNITS="     						!CASE(2) CM-1 ,A.U. KCAL, KLVN
      CHARACTER(LEN=8) :: r_units_word="R_UNITS="          						!CASE(3) A.U., ANGS
      CHARACTER(LEN=10) :: nmb_terms_word="NMB_TERMS="     						!CASE(4)
      CHARACTER(LEN=6) :: terms_word="TERMS="              						!CASE(5)
      CHARACTER(LEN=6) :: size_R_word="GRD_R="             						!CASE(6)
      CHARACTER(LEN=8) :: size_r_vib_word="GRD_VIB="       						!CASE(7)
      CHARACTER(LEN=9) :: size_thet_word="GRD_ANG1="       						!CASE(8)
      CHARACTER(LEN=9) :: size_bet_word="GRD_ANG2="        						!CASE(9) 
      CHARACTER(LEN=18) :: size_gam_word="GRD_ANG3="       						!CASE(10)
      CHARACTER(LEN=10)  :: matrix_read_word = "READ_MTRX=" 					!CASE(11)
      CHARACTER(LEN=9)  :: run_prog_word = "PROG_RUN=" 							!CASE(12)
      CHARACTER(LEN=13) :: dist_coord_word = "PES_ADDITIVE="					!CASE(13)
      CHARACTER(LEN=12)	:: make_grid_file_word ="VGRID_FILE="   				!CASE(14)
      CHARACTER(LEN=15) :: calculate_expansion_word = "CALC_EXPANSION=" 		!CASE(15)
      CHARACTER(LEN=6)	:: L_word="L_MAX="  									!CASE(16)
      CHARACTER(LEN=6)	:: NJU_word="M_MAX="  									!CASE(17)
      CHARACTER(LEN=7)	:: L1_word="L1_MAX="  									!CASE(18)
      CHARACTER(LEN=7)	:: L2_word="L2_MAX="  									!CASE(19)
      CHARACTER(LEN=7)	:: NJU1_word="M1_MAX="  								!CASE(20)
      CHARACTER(LEN=7)	:: NJU2_word="M2_MAX="  								!CASE(21)
      CHARACTER(LEN=12 ):: r_grid_file_def = "RGRID_FILE=" 						!CASE(22)
      CHARACTER(LEN=7)	:: L1_mword="L1_MIN="  									!CASE(23)
      CHARACTER(LEN=7)	:: L2_mword="L2_MIN="  									!CASE(24)
      CHARACTER(LEN=7)	:: ir_bgn_word = "IR_BGN="								!CASE(25)
      CHARACTER(LEN=7)	:: ir_fin_word = "IR_FIN=" 								!CASE(26)
      CHARACTER(LEN=10)	:: print_matrix_defined_word = "SAVE_MTRX=" 			!CASE(27)
      CHARACTER(LEN=14) :: ini_pot_expansion_word="TERMS_ONFLY=" 				!CASE(28)
      CHARACTER(LEN=9) :: unformatted_reading_word ="UNFORMAT="					!CASE(29)
      CHARACTER(LEN=15)	:: test_exp_word = "PRINT_DIAGONAL=" 					!CASE(30)			!Bikram Oct'18
      CHARACTER(LEN=11) :: terms_file_defined_word  =	"TERMS_FILE=" 			!CASE(31)
      CHARACTER(LEN=11)	:: eq_grid_word = "RGRID_EQDS="               			!CASE(32)
      CHARACTER(LEN=12)	:: elast_m_proj_word = "ELAST_MPROJ="               	!CASE(33) 	
      CHARACTER(LEN=9)	:: calc_matrix_word = "MIJ_CALC="               		!CASE(34) 	      	  
      CHARACTER(LEN=10)	:: bk_mij_shift = "MIJ_SHIFT="               			!CASE(35) 	   		!Bikram Dec 2019   	  
      CHARACTER(LEN=11)	:: bk_mij_cutoff = "MIJ_CUTOFF="               			!CASE(36) 	   		!Bikram Dec 2019   	  
      CHARACTER(LEN=12)	:: bk_mij_multiprint = "PARALLEL_IO="               	!CASE(37) 	   		!Bikram Nov 2020   	  
      CHARACTER(LEN=12)	:: bk_eq_grd_gamma = "EQ_ANG_GRID="               		!CASE(38) 	   		!Bikram Nov 2021   	  
      CHARACTER(LEN=15)	:: bk_on = "ORTHONORMALITY="               				!CASE(39) 	   		!Bikram Nov 2021   	  
      CHARACTER(LEN=10)	:: bk_mtrx_path = "MTRX_PATH="              			!CASE(40) 	   		!Bikram Nov 2021   	  
      CHARACTER(LEN=14)	:: bk_w3j_fact = "W3J_RECURSIVE="               		!CASE(41) 	   		!Bikram Feb 2022   	  
      CHARACTER(LEN=11)	:: bk_truncate_mtrx = "RETRUNCATE="               		!CASE(42) 	   		!Bikram May 2022   	  
      CHARACTER(LEN=10)	:: bk_cutoff_r1 = "CUTOFF_R1="               			!CASE(43) 	   		!Bikram June 2022   	  
      CHARACTER(LEN=10)	:: bk_cutoff_r2 = "CUTOFF_R2="               			!CASE(44) 	   		!Bikram June 2022   	  
      CHARACTER(LEN=14)	:: bk_identical_pes = "IDENTICAL_PES="               	!CASE(45) 	   		!Bikram Aug 2022   	  
      CHARACTER(LEN=9)	:: bk_calc_rms = "CALC_RMS="               				!CASE(46) 	   		!Bikram Aug 2022
      CHARACTER(LEN=9)	:: bk_rms_ang1 = "RMS_GRD1="							!CASE(47) 	   		!Bikram Aug 2022
      CHARACTER(LEN=9)	:: bk_rms_ang2 = "RMS_GRD2="							!CASE(48) 	   		!Bikram Aug 2022
      CHARACTER(LEN=9)	:: bk_rms_ang3 = "RMS_GRD3="							!CASE(49) 	   		!Bikram Aug 2022
      CHARACTER(LEN=6)	:: bk_rms_r = "RMS_R="									!CASE(50) 	   		!Bikram Aug 2022
      CHARACTER(LEN=10)	:: bk_ident_terms = "ODD_L1L2L="						!CASE(51) 	   		!Bikram Aug 2022
      CHARACTER(LEN=9)	:: bk_axl_sym1 = "AXL_SYM1="							!CASE(52) 	   		!Bikram Sep 2022
      CHARACTER(LEN=9)	:: bk_axl_sym2 = "AXL_SYM2="							!CASE(53) 	   		!Bikram Sep 2022
      CHARACTER(LEN=9)	:: bk_equ_sym1 = "EQU_SYM1="							!CASE(54) 	   		!Bikram Sep 2022
      CHARACTER(LEN=9)	:: bk_equ_sym2 = "EQU_SYM2="							!CASE(55) 	   		!Bikram Sep 2022
      CHARACTER(LEN=10)	:: bk_rebalance = "REBALANCE ="							!CASE(56) 	   		!Bikram Oct 2022
      CHARACTER(LEN=12)	:: bk_reblnc_comp = "BALANCED_MIJ="						!CASE(57) 	   		!Bikram Oct 2022
      CHARACTER(LEN=1) buffer
      LOGICAL key_used
      IF(length.le.0) RETURN
      key_used = .FALSE.	  
      posit = 1
      decrement = place - posit
      buffer = 	inp(posit:posit)  
      DO WHILE(buffer.ne."=" .and. posit.lt.length)
      posit = posit + 1
      buffer = 	inp(posit:posit)
!      PRINT*, "WANNA BELIEVE",buffer	  
      ENDDO
!      PRINT*,"WANNA CHECK",posit,inp(1:posit)	  
      IF(inp(1:posit).eq.expan_word) THEN
      key = 1
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.energy_units_word) THEN
      key = 2
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.r_units_word) THEN
      key = 3
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.nmb_terms_word) THEN
      key = 4
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.terms_word) THEN
      key = 5
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.size_R_word) THEN
      key = 6
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.size_r_vib_word) THEN
      key = 7
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.size_thet_word) THEN
      key = 8
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.size_bet_word) THEN
      key = 9
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.size_gam_word) THEN
      key = 10
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.matrix_read_word) THEN
      key = 11
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.run_prog_word) THEN
      key = 12
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.dist_coord_word) THEN
      key = 13
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	
      IF(inp(1:posit).eq.make_grid_file_word) THEN
      key = 14
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.calculate_expansion_word) THEN
      key = 15
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.L_word) THEN
      key = 16
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.NJU_word) THEN
      key = 17
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.L1_word) THEN
      key = 18
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.L2_word) THEN
      key = 19
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.NJU1_word) THEN
      key = 20
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.NJU2_word) THEN
      key = 21
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.r_grid_file_def) THEN
      key = 22
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.L1_mword) THEN
      key = 23
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.L2_mword) THEN
      key = 24
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.ir_bgn_word) THEN
      key = 25
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.ir_fin_word) THEN
      key = 26
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.print_matrix_defined_word) THEN
      key = 27
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.ini_pot_expansion_word) THEN
      key = 28
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.unformatted_reading_word) THEN
      key = 29
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.test_exp_word) THEN
      key = 30
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.terms_file_defined_word) THEN
      key = 31
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.eq_grid_word) THEN
      key = 32
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF 
      IF(inp(1:posit).eq.elast_m_proj_word) THEN
      key = 33
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF 
      IF(inp(1:posit).eq.calc_matrix_word) THEN
      key = 34
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.bk_mij_shift) THEN
      key = 35
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.bk_mij_cutoff) THEN
      key = 36
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.bk_mij_multiprint) THEN
      key = 37
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.bk_eq_grd_gamma) THEN
      key = 38
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.bk_on) THEN
      key = 39
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.bk_mtrx_path) THEN
      key = 40
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.bk_w3j_fact) THEN
      key = 41
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.bk_truncate_mtrx) THEN
      key = 42
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.bk_cutoff_r1) THEN
      key = 43
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.bk_cutoff_r2) THEN
      key = 44
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.bk_identical_pes) THEN
      key = 45
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.bk_calc_rms) THEN
      key = 46
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.bk_rms_ang1) THEN
      key = 47
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.bk_rms_ang2) THEN
      key = 48
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.bk_rms_ang3) THEN
      key = 49
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
      IF(inp(1:posit).eq.bk_rms_r) THEN
      key = 50
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF	  
      IF(inp(1:posit).eq.bk_ident_terms) THEN
      key = 51
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
	  IF(inp(1:posit).eq.bk_axl_sym1) THEN
      key = 52
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
	  IF(inp(1:posit).eq.bk_axl_sym2) THEN
      key = 53
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
	  IF(inp(1:posit).eq.bk_equ_sym1) THEN
      key = 54
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
	  IF(inp(1:posit).eq.bk_equ_sym2) THEN
      key = 55
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
	  IF(inp(1:posit).eq.bk_rebalance) THEN
      key = 56
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF
	  IF(inp(1:posit).eq.bk_reblnc_comp) THEN
      key = 57
      key_used = .TRUE.	  
      IF(key_word_used(key).eq.1) THEN
      PRINT*,inp(1:posit)	  
      STOP "ERROR:THIS WORD IS ALREADY USED"
      ENDIF
      key_word_used(key) = 1
      ENDIF

      IF(.not.key_used) THEN
      PRINT*,inp(1:posit)	  
      STOP 
     & "ERROR IN POTENTIAL: WORD NOT FOUND OR INPUT ABOVE IS INCORRECT"
      ENDIF	  
      place = posit + decrement + 1	    
      END SUBROUTINE KEY_WORD_POTENTIAL
      SUBROUTINE OUTPUT
      USE CONSTANTS	  !!!! ANALYZING OUTPUT - CREATING ARRAYS-BASIS SET
      USE VARIABLES
      USE FINE_STRUCT_LEVELS	  
      USE MPI_DATA
      USE MPI	  
      IMPLICIT NONE
      LOGICAL exst,EXST_STATE	  
      INTEGER i,k, nlvl,decr,j_t,k_t,eps_t,nchann_est,v_t,e_t,n_prev
      INTEGER nchann1_est,nchann2_est,i_ground
      INTEGER nlvl1,nlvl2,n_prev1,n_prev2	  
      INTEGER, PARAMETER :: max_nmb_chann = 100000,j_max_allowed = 10000
     & ,v_max_allowed = 1000
      REAL*8, ALLOCATABLE :: e_temp(:),e1_temp(:),e2_temp(:)
      INTEGER, ALLOCATABLE :: j_ch_tmp(:),k_ch_tmp(:),v_ch_temp(:),
     & eps_ch_tmp(:), ja(:),ka_ch_tmp(:),kc_ch_tmp(:),indch_ABC(:,:,:)
      INTEGER, ALLOCATABLE :: k1_ch_tmp(:),eps1_ch_tmp(:),ka1_ch_tmp(:),
     & kc1_ch_tmp(:)
      INTEGER, ALLOCATABLE :: k2_ch_tmp(:),eps2_ch_tmp(:),ka2_ch_tmp(:),
     & kc2_ch_tmp(:)	 
      INTEGER KRONEKER,ka_t,kc_t,dk_t,nchnl_ini_guess,
     & buffer_ch1,buffer_ch2
      INTEGER, ALLOCATABLE :: j1_ch_tmp(:),j2_ch_tmp(:),v1_ch_temp(:),
     & 	v2_ch_temp(:),ja1(:),ja2(:)  
      INTEGER j1_t,j2_t,v1_t,v2_t,stat_of_file
      INTEGER, ALLOCATABLE :: SORT_STORAGE(:,:)
      INTEGER TIME_DATA_EXE(8)
      CHARACTER(LEN=4) :: am_pm
      CHARACTER(LEN=1) :: par_sym	  
      REAL*8 Ech1,Ech2,Eground
      REAL*8 rot_const_inp(11)
      REAL*8 j_f_p_b	   
      INTEGER rot_state_inp(6)
	  integer kappa1, kappa2
	  character (len=100) :: bk_temp_dir					!Bikram
      EXTERNAL KRONEKER,EXST_STATE	  
      IF(.not.scatter_param_defined) 
     & STOP "ERROR: MASS AND RANGE NO DEFINED"
      IF(MYID.eq.0) THEN	 
      OPEN(UNIT=1, FILE = "USER_INPUT_CHECK.out")
      WRITE(1,'(2(a))')
     &"|--------------------------------------------------",
     & "--------------------------------------------------|"
	  write(1,'(3(a))')
     & "|                                              ", 
     & "MQCT_2022", 
     & "                                             |"
	  write(1,'(3(a))')
     & "|                               ", 
     & "THE PROGRAM IS TO PERFORM ROTATIONALLY", 
     & "                               |"
	  write(1,'(3(a))')
     & "|                        ",
     & "AND VIBRATIONALLY INELASTIC SCATTERING CALCULATIONS.",
     & "                        |"
	  write(1,'(3(a))')
     & "|                                 ",
     & "TO CITE MQCT_2022 PLEASE REFER TO:",
     & "                                 |"
      write(1,'(3(a))')
     & "|                        ",
     & "B. MANDAL, C. JOY, D. BOSTAN, A. ENG and D. BABIKOV,",
     & "                        |"
      write(1,'(3(a))')
     & "|                          ",
     & "JOURNAL OF CHEMICAL THEORY AND COMPUTATION, 2022",
     & "                          |"
      write(1,'(2(a))')
     &"|--------------------------------------------------",
     & "--------------------------------------------------|"
      WRITE(1,*)
      CALL date_and_time(VALUES=TIME_DATA_EXE)
      IF(TIME_DATA_EXE(5).lt.12 .and.TIME_DATA_EXE(5).gt.0) THEN
      am_pm = 'a.m.'  	  
      ENDIF
      IF(TIME_DATA_EXE(5).eq.12) THEN
      am_pm = 'p.m.'  	  
      ENDIF
      IF(TIME_DATA_EXE(5).eq.0) THEN
      am_pm = 'a.m.'
      TIME_DATA_EXE(5) = 12  	  
      ENDIF	
      IF(TIME_DATA_EXE(5).gt.12) THEN
      am_pm = 'p.m.'
      TIME_DATA_EXE(5) =TIME_DATA_EXE(5)- 12  	  
      ENDIF		  
      WRITE(1,'(a22,i2,a1,i2,a1,i4,a4,i3,a1,a4,a1,i2,a28)')
     & "THIS RUN PERFORMED ON ",
     & TIME_DATA_EXE(2),'/',TIME_DATA_EXE(3),'/'
     & ,TIME_DATA_EXE(1),' at ',TIME_DATA_EXE(5)," ",am_pm," ",
     & TIME_DATA_EXE(6)," minutes, LOCAL MACHINE TIME" 	 
      ENDIF	  
      IF(ener_auto_defined) THEN
      IF(nmbr_of_enrgs.le.1) STOP
     & "ERROR:YOU NEED MORE THAN ONE COLLISIONAL ENERGY."
      IF(dU.eq.0d0) THEN
      U(1) = U_min
      IF(U_min.eq.0) STOP "ERROR: Umin IS NON_POSITIVE"
      U(nmbr_of_enrgs) = U_max
      IF(U_max.eq.0) STOP "ERROR: Umax IS NON_POSITIVE"
      dU = (U_max-U_min)/(nmbr_of_enrgs-1)
      IF(dU.lt.0d0) STOP "ERROR:Umax < Umin"	  
      DO i=2,nmbr_of_enrgs-1
      U(i) = U(1) + dU*(i-1d0)
      ENDDO	  
      ENDIF
      IF(U_max.eq.0d0) THEN
      U(1) = U_min
      IF(U_min.eq.0) STOP "ERROR: Umin IS NON-POSITIVE"
      IF(dU.le.0d0) STOP "ERROR:dU IS NON-POSITIVE"	  
      DO i=2,nmbr_of_enrgs
      U(i) = U(i-1) + dU	  
      U_max = U(nmbr_of_enrgs)
      ENDDO	  
      ENDIF
      IF(U_min.eq.0d0) THEN
      U(nmbr_of_enrgs) = U_max
      IF(U_max.eq.0) STOP "ERROR: Umax IS NON-POSITIVE"
      IF(dU.le.0d0) STOP "ERROR:dU IS NON-POSITIVE"	  
      DO i=nmbr_of_enrgs,2,-1
      U(i-1) = U(i) - dU	  
      U_min = U(1)
      ENDDO  
      IF(U_min.eq.0) STOP "ERROR: Umin IS NON-POSITIVE"	  
      ENDIF			  
      ENDIF
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a14,1x,a)') "YOUR JOB NAME:",label
      WRITE(1,*)	  
      WRITE(1,'(a28,1x,f7.3)') "REDUCED MASS,A.M.U.: MASS = ",mass_red
      WRITE(1,*)
      WRITE(1,'(a20,1x,a11,1x,f6.2,1x,a9,f6.2,1x,a4)')
     & "THE POTENTIAL RANGE:",
     & "FROM Rmin =",R_min_dist,'TO Rmax =', R_max_dist,"BOHR"
      WRITE(1,*)
      WRITE(1,'(a27)')
     & "COLLISIONAL ENERGIES, CM-1:"
      DO i=1,nmbr_of_enrgs
      WRITE(1,'(a2,i4,a2,f9.3,a1)') "E(",i,")=",U(i),";" 	  
      ENDDO	  
 	 
      WRITE(1,*)
      IF(J_tot_min.gt.J_tot_max) THEN
      WRITE(1,'(a45)') "ERROR: JTOT_MAX MUST BE NO LESS THAN JTOT_MIN"	  
      STOP
      ENDIF
      ENDIF
      IF(MYID.eq.0) THEN	  
      IF(b_impact_defined) THEN
      WRITE(1,'(a30,f10.3)') "MAX IMPACT PARAMETER IS,a.u = ",
     & b_impact_parameter	  
      ELSE	  
      WRITE(1,'(a36,i3,1x,a2,1x,i3)')
     & "TOTAL ANGULAR MOMENTA INCLUDED: FROM"
     &,J_tot_min,"TO",J_tot_max
      WRITE(1,*)
      ENDIF	
      IF(delta_l_step.ge.1) THEN
      WRITE(1,'(a36,1x,i3)')
     & "ORBITAL ANGULAR MOMENTUM STEP_SIZE:=", delta_l_step
      ELSE
      PRINT*,"DL IS NOT POSITIVE"	  
      ENDIF	  
      IF(.not.deflect_fun_defined) THEN
!      WRITE(1,'(a30)') "PRINT DEFLECTION FUNCTION: YES"	  
!      ELSE
      WRITE(1,'(a29)') "PRINT DEFLECTION FUNCTION: NO"	  
      ENDIF
      WRITE(1,*)
      IF(monte_carlo_defined)	  
     & WRITE(1,'(a33,f7.2)') "MONTE-CARLO INTEGRATION ERROR, %:",
     & mnt_crl_intgrt_err
      WRITE(1,*)
	  
!      IF(orbit_traj_defined) THEN
!      WRITE(1,'(a60)') 
!     & "IF ORBITING TRAJECTORIES EXIST, PRINT IMPACT PARAMETERS: YES"
!      ELSE
!      WRITE(1,'(a59)')
!     & "IF ORBITING TRAJECTORIES EXIST, PRINT IMPACT PARAMETERS: NO"	  
!      ENDIF
!      WRITE(1,*)
      IF(tm_lim_defined) WRITE(1,'(a,1x,e11.4)')
     & "ALL TRAJECTORIES WILL BE STOPPED AFTER TIME_LIMT=",time_lim
      WRITE(1,*)

! Bikram Start:
      WRITE(1,'(a34,i4,1x,a10,i5,1x,a12)') 
     & "TRAPPED TRAJECTIRY WILL STOP AFTER",
     & numb_orb_prds,
     & "ROUNDS AND",numb_oscl_prds,"OSCILLATIONS"
! Bikram End.

      WRITE(1,*)
      IF(rk4_defined) odeint_defined = .FALSE.	  
      IF(odeint_defined) WRITE(1,'(a16,1x,a38)') "PROPAGATOR TYPE:",
     &	                              "ADAPTIVE STEP SIZE CONTROL RUNGE-K
     &UTTA"
      IF(rk4_defined) WRITE(1,'(a16,1x,a22)') "PROPAGATOR TYPE:",
     & "RUNGE-KUTTA 4-TH ORDER"
      IF(odeint_defined) WRITE(1,'(a66,1x,f12.2)')"MAXIMUM TIME STEP IN 
     &ADAPTIVE STEP SIZE CONTROL SCHEME TIME_STEP = ",time_step
      IF(odeint_defined) WRITE(1,'(a66,1x,f12.2)')"MINIMUM TIME STEP IN 
     &ADAPTIVE STEP SIZE CONTROL SCHEME MIN_TSTEP = ",min_t_stp
      IF(rk4_defined) WRITE(1,'(a34,1x,f12.2)')"TIME STEP IN RUNGE-KUTT
     &A 4TH ORDER",time_step	
      IF(odeint_defined) WRITE(1,'(a27,e11.4)')
     & "ADAPTIVE STEP SIZE ACCURACY",eps_odeint	
      WRITE(1,*)
      IF(number_of_channels_defined.and.emax_defined) THEN
      WRITE(1,'(a81)')                         "ERROR: NUMBER OF CHANNEL
     &S AND CLOSED CHANNEL ENERGY ARE NOT INDEPENDENT VARIABLES"	 
      STOP
      ENDIF
      IF(coupled_states_defined) 
     & WRITE(1,'(a41)') "COUPLED STATES APPROXIMATION WILL BE USED"	  
      IF(bikram_adiabatic) WRITE(1,'(a47)') 
     & "ADIABATIC TRAJECTORY APPROXIMATION WILL BE USED"	 
      if(bikram_save_traj) then
	  do i = 1, nmbr_of_enrgs
	  write(bk_temp_dir, '(a,a,f0.5)') bk_directory, '/', U(i)
	  call system ( "mkdir -p " // trim(bk_directory) )
	  write(bk_temp_dir, '(a,a,f0.5,a)') bk_directory, '/', U(i), '/*'
	  if(.not. monte_carlo_defined) 
     & call system ( "rm -r " // trim(bk_temp_dir) )
	  enddo
      WRITE(1,'(a,a)') 
     & "APPROXIMATED TRAJECTORIES WILL BE SAVED IN ", bk_directory
	  endif
      ENDIF  
      IF(delta_l_step.le.0) STOP
! Bikram Start Jan,19:
      IF(.not. number_of_channels_defined .and. .not. emax_defined)THEN
      number_of_channels=6
	  number_of_channels_defined = .TRUE.
	  jmin=0
	  jmin1=0
	  jmin2=0 
      ENDIF	  
! Bikram End.
      SELECT CASE(system_type)
      CASE(1)
      IF(MYID.eq.0) 	  
     & WRITE(1,'(a23,1x,a25)')"THE COLLISIONAL SYSTEM:",
     & "RIGID DIATOMIC TOP + ATOM"
      IF(user_defined) THEN
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a7,2x,a60)') "WARNING:",
     & "THE WAVEFUNCTIONS AND ENERGY LEVELS MUST BE SUPPLIED BY USER"
      IF(rot_const_defined)WRITE(1,'(a27)')"ROT. CONST WILL BE INGNORED"
      ENDIF	  
      INQUIRE( FILE=USR_INP_LEVEL, EXIST=exst )
      IF(.not.exst)	THEN
      IF(MYID.eq.0)PRINT*,USR_INP_LEVEL, " NOT FOUND"
      STOP	  
      ENDIF	  
      OPEN(325,FILE=USR_INP_LEVEL,
     & STATUS="OLD",ACTION="READ")
      IF(number_of_channels.le.0) THEN
      IF(MYID.eq.0) PRINT*,"NUMBER OF CHANNELS NOT SPECIFIED"
      STOP	  
      ENDIF
      IF(n_r_vib.le.0) THEN
      IF(MYID.eq.0) PRINT*,"INTEGRATION GRID NOT DEFINED"
      STOP	  
      ENDIF	  
      IF(.not.channels_defined) THEN	  
      ALLOCATE(j_ch(number_of_channels))
      ENDIF
      IF(.not.energy_defined) THEN	  
      ALLOCATE(E_ch(number_of_channels))
      ENDIF
      READ(325,*) !!! HEADER
      jmax_included = 0
      DO i=1,number_of_channels
      READ(325,'(i4,1x,i3,1x,e17.10)',
     & IOSTAT = stat_of_file)
     & decr,j_ch(i), E_ch(i)
      IF(jmax_included.lt.j_ch(i)) jmax_included= j_ch(i)
      IF(stat_of_file.gt.0 .or. i.ne.decr) THEN
      PRINT*,"WRONG CHANNELS IN USER INPUT FILE"	  
      STOP
      ENDIF 	  
      ENDDO
      CLOSE(325)
      channels_defined = .TRUE.
      energy_defined = .TRUE.	
	   	  	 
      ELSE
      IF(BE.eq.0) STOP "ERROR:SUPPLY POSITVE BE"
      IF(DE.GT.BE) STOP"ERROR: CHECK INPUT-DE SHOULD MUCH LESS THAN BE"	  
      IF(MYID.eq.0) THEN
      WRITE(1,'(a32,f10.4)') "ROTATIONAL CONSTANT, CM-1, BE = ",Be	  
      WRITE(1,'(a23,f12.7)') "DISTORTION, CM-1, DE = ",DE
      IF(fine_structure_defined) THEN
      WRITE(1,"(a21,i2)") "MULTIPLIICITY 2S+1 = ",SPIN_FINE !!! WRITE A CONDTION FOR HALF SPIN
      WRITE(1,"(a22,i2)") "ELECTRONIC MOMENTUM = ",LORB_FINE		  
      WRITE(1,'(a31,f10.4)') "SPIN_ORBIT_COUPLING, CM-1, A = ",
     & a_spin_orbit_fine
      WRITE(1,'(a35,f10.4)') "SPIN_SPIN_COUPLING, CM-1, Lambda = ",
     & lambda_fine      	 
      WRITE(1,'(a38,e11.4)') "SPIN_SPIN_COUPLING_D, CM-1, LambdaD = ",
     & lambda_fine_d
      WRITE(1,'(a40,e11.4)')"SPIN_SPIN_COUPLING_DD, CM-1, LambdaDD = ",
     & lambda_fine_dd
      WRITE(1,'(a33,e11.4)') "SPIN_ROT_COUPLING, CM-1, Gamma = ",
     & gamma_fine
      WRITE(1,'(a36,e11.4)') "SPIN_ROT_COUPLING_D, CM-1, GammaD = ",
     & gamma_fine_d
      WRITE(1,'(a34,e11.4)') "P_LAMBDA_DOUBLING, CM-1, GammaD = ",
     & p_double_fine
      WRITE(1,'(a34,e11.4)') "Q_LAMBDA_DOUBLING, CM-1, GammaD = ",
     & q_double_fine    	 
      ENDIF 
      ENDIF  
      IF(fine_structure_defined) THEN
      ALLOCATE(b_fine_coeff(SPIN_FINE,number_of_channels))
      ENDIF	  
      IF(jmax.ge.0 .and. jmin.ge.0) THEN
      IF(jmax.lt.jmin) STOP "ERROR: JMAX MUST BE > OR = THAN JMIN"	  
      number_of_channels_defined = .TRUE.
      energy_defined = .TRUE.
      channels_defined = .TRUE.
      number_of_channels = jmax - jmin + 1
      ALLOCATE(j_ch(number_of_channels))
      ALLOCATE(E_ch(number_of_channels)) 	  
      DO i=1,number_of_channels
      j_ch(i) = i - 1+jmin
      E_ch(i) = Be*j_ch(i)*(j_ch(i)+1d0) -
     &	  De*(j_ch(i)*(j_ch(i)+1d0))**2	  
      ENDDO	  
      ENDIF	  
      IF(number_of_channels_defined) THEN
      IF(energy_defined.and. .not.channels_defined) STOP 
     & "ERROR:FOR MANUALLY DEFINED ENERGY YOU SHOULD SPECIFY CHANNELS"	  
      IF(.not.energy_defined) ALLOCATE(E_ch(number_of_channels))

      IF(.not.channels_defined) THEN
      ALLOCATE(j_ch(number_of_channels)) 	  
      DO i=1,number_of_channels
      j_ch(i) = i - 1+jmin
      E_ch(i) = Be*j_ch(i)*(j_ch(i)+1d0) -
     &	  De*(j_ch(i)*(j_ch(i)+1d0))**2	  
      ENDDO
      ELSE
      DO i=1,number_of_channels
      IF(.not.energy_defined) 
     &      E_ch(i) = Be*j_ch(i)*(j_ch(i)+1d0) -
     &	  De*(j_ch(i)*(j_ch(i)+1d0))**2
      IF(fine_structure_defined) THEN
      rot_const_inp(1) = a_spin_orbit_fine
      rot_const_inp(2) = Be
      rot_const_inp(3) = De
      rot_const_inp(4) = lambda_fine
      rot_const_inp(5) = gamma_fine
      rot_const_inp(6) = lambda_fine_d
      rot_const_inp(7) = gamma_fine_d
      rot_const_inp(8) = lambda_fine_dd
      rot_const_inp(9) = HE_ROT	
      rot_const_inp(10) = p_double_fine
      rot_const_inp(11) = q_double_fine	  
      rot_state_inp(1) = LORB_FINE
      rot_state_inp(2) = SPIN_FINE
      rot_state_inp(3) = 0
      rot_state_inp(4) = f_ch(i)
      IF (LORB_FINE.ne.0) rot_state_inp(5) = par_lorb_ch(i)
      rot_state_inp(6) = j_ch(i)	  
      CALL FINE_ENERGY_DIATOMIC(E_ch(i),
     & rot_const_inp, rot_state_inp,
     & b_fine_coeff(:,i),SPIN_FINE)	 
      ENDIF	  
      ENDDO	  
      ENDIF

      channels_defined = .TRUE.
      energy_defined = .TRUE.	  
      ENDIF	
      IF(emax_defined.and. .not.channels_defined) THEN
      nchann_est = int(sqrt(EMAX/Be))+1
      i = nchann_est - 1	  
      IF(Be*i*(i+1d0) - De*(i*(i+1d0))**2 .lt. EMAX) THEN
      DO WHILE(Be*i*(i+1d0) - De*(i*(i+1d0))**2 .lt. EMAX)
      i = i +1
      ENDDO
      nchann_est = i - 1	  
      ENDIF
      IF(Be*i*(i+1d0) - De*(i*(i+1d0))**2 .GT. EMAX) THEN
      DO WHILE(Be*i*(i+1d0) - De*(i*(i+1d0))**2 .GT. EMAX)
      i = i -1
      ENDDO
      nchann_est = i  
      ENDIF
      number_of_channels = nchann_est+1
      IF(number_of_channels_defined) STOP"ERROR: DEFINE N_CHNLS OR EMAX"
      ALLOCATE(j_ch(number_of_channels))
      ALLOCATE(E_ch(number_of_channels)) 	  
      DO i=1,number_of_channels
      j_ch(i) = i - 1
      E_ch(i) = Be*j_ch(i)*(j_ch(i)+1d0) -
     &	  De*(j_ch(i)*(j_ch(i)+1d0))**2	  
      ENDDO
      channels_defined = .TRUE.
      energy_defined = .TRUE.	  
      ENDIF
      IF(MYID.eq.0) THEN	  
      WRITE(1,*)
      ENDIF	  
      jmax_included	 = 0
      IF(para_ortho_defined) CALL SYMMETRY_SORT_PARA_ORTHO	  
      IF(MYID.eq.0) THEN
      IF(.not.fine_structure_defined) THEN	  
      WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"	  
      WRITE(1,'(a5,9x,a2,9x,a2)')"#N = "
     & ,"J","E"
      ELSE
      IF(LORB_FINE.eq.0) THEN	  
      WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"	  
      WRITE(1,'(a5,9x,a2,9x,a2,9x,a2)')"#N = "
     & ,"J","F","E"	  
      ELSE
      WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"	  
      WRITE(1,'(a5,9x,a2,9x,a2,9x,a2,9x,a2,9x,a2)')"#N = "
     & ,"J","N","F","P","E"	  
      ENDIF	  
	  
      ENDIF

	  
      ENDIF
      IF(.not.fine_structure_defined) THEN		  
      DO i=1,number_of_channels
      IF(jmax_included.lt.j_ch(i)) jmax_included = j_ch(i)	  
      IF(MYID.eq.0) THEN      
      WRITE(1,'(i5,i11,f11.3)')i,
     & j_ch(i),E_ch(i)
      ENDIF	 
      ENDDO
      ELSE
      IF (LORB_FINE.eq.0) THEN	  
	  
      DO i=1,number_of_channels
      IF(jmax_included.lt.j_ch(i)) jmax_included = j_ch(i)	  
      IF(MYID.eq.0) THEN      
      WRITE(1,'(i5,i11,i11,f11.3)')i,
     & j_ch(i),f_ch(i),E_ch(i)
      ENDIF	 
      ENDDO	
      ELSE
      ALLOCATE(j_h_ch(number_of_channels))
      Eground	 =   E_ch(1)
      DO i=2,number_of_channels
      IF(E_ch(i).lt.Eground) Eground = E_ch(i)	  
      ENDDO
      DO i=1,number_of_channels
      E_ch(i)  =  E_ch(i)	- 	Eground 
      ENDDO		  
      DO i=1,number_of_channels
      IF(jmax_included.lt.j_ch(i)) jmax_included = j_ch(i)	  
 	  j_f_p_b  =  j_ch(i) - (-1)**f_ch(i)/2d0
      j_h_ch(i) = j_f_p_b 	  
	  par_sym  = 'e'
      IF(par_lorb_ch(i).eq.1)  par_sym  = 'f'
      IF(MYID.eq.0) THEN 	  
      WRITE(1,'(i5,f11.1,i11,i11,a11,f11.3)')i,
     & j_h_ch(i),j_ch(i),f_ch(i),par_sym,E_ch(i)
      ENDIF	 
      ENDDO		  
      ENDIF



	  
      ENDIF	  
      ENDIF	
      IF(MYID.eq.0) THEN	  
      WRITE(1,*)
      ENDIF	  
      IF(ini_chann_defined .and. .not.all_level_included) THEN
      chann_ini = 0	  
      DO i=1,number_of_channels
      IF(fine_structure_defined) THEN	
      IF(LORB_FINE.eq.0) THEN	  
      IF(j_ini.eq.j_ch(i) .and. 
     &  f_ini.eq.f_ch(i) ) chann_ini = i
      ELSE
      IF(j_ini.eq.j_ch(i) .and. 
     &  f_ini.eq.f_ch(i) .and.
     & par_lorb_ch(i).eq.par_orb_ini ) chann_ini = i	  
      ENDIF	  
      ELSE
      IF(j_ini.eq.j_ch(i)) chann_ini = i	  
      ENDIF
      ENDDO	
      IF(chann_ini.eq.0) THEN
      IF(MYID.eq.0)CALL FLUSH(1)	  
      STOP "ERROR: INITIAL CHANNEL WRONG"
      ENDIF
      IF(.not.fine_structure_defined) THEN 	  
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a17,a4,i4)') "INITIAL STATE IS:",
     & "J = ",j_ini
      ENDIF	
      ELSE
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a17,a4,i4,a4,i4)') "INITIAL STATE IS:",
     & "J = ",j_ini,"F = ",f_ini
!      WRITE(*,*) b_fine_coeff	 
      ENDIF		  
      ENDIF	  
      ENDIF	 
	  
      IF(.not.ini_chann_defined .and. channels_defined .and.
     & energy_defined) THEN
	  i_ground = 1
      DO i=1,number_of_channels
      IF(E_ch(i_ground).gt.E_ch(i)) i_ground = i   
      ENDDO	 
      chann_ini = i_ground
      ini_chann_defined = .TRUE.	  
      ENDIF
 	  
      IF(MYID.eq.0) WRITE(1,*)
      IF(MYID.eq.0) THEN	  
      IF(angs_unit)WRITE(1,'(a35)')"THE POTENTIAL R-UNITS ARE ANGSTROMS"
      IF(au_unit_r)WRITE(1,'(a31)') "THE POTENTIAL R-UNITS ARE BOHRS"
      IF(cm_unit)WRITE(1,'(a31)')"THE POTENTIAL ENERGY IS IN CM-1"
      IF(au_unit_e)WRITE(1,'(a34)')"THE POTENTIAL ENERGY IS IN HARTREE"
      IF(klv_unit)WRITE(1,'(a33)')"THE POTENTIAL ENERGY IS IN KELVIN"
      IF(kal_unit)WRITE(1,'(a35)')"THE POTENTIAL ENERGY IS IN KCAL/MOL"	  
      ENDIF      
      IF(.not.angs_unit .and..not.au_unit_r )   
     & STOP "ERROR:R_UNITS NOT DEFINED"
      IF(.not.cm_unit .and..not.au_unit_e.and..not.klv_unit
     & .and. .not.kal_unit)   
     & STOP "ERROR:ENERGY NOT DEFINED"
      IF(MYID.eq.0)WRITE(1,*) 	 
      IF(expansion_defined) THEN
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a42)') "USER SUPPLIED POTENTIAL TERMS WILL BE USED"
      WRITE(1,*)
      WRITE(1,'(a45,i4)')
     & "NUMBER OF POTENTIAL TERMS INCLUDED N_TERMS = ",nterms
      WRITE(1,'(a38)')
     & "THE POTENTIAL TERMS ARE, (LAMBDA):"  
      ENDIF
      IF(MYID.eq.0) THEN	  
      DO i=1,nterms
      WRITE(1,'(a1,i3,a3)',ADVANCE ="NO")
     & '(',L_TERM(i),'); '	  
      ENDDO
      IF(fine_structure_defined.and.SPIN_FINE.eq.2) THEN
      WRITE(1,'(a40)')
     & "THE POTENTIAL TERMS ARE, (NJU = 0 or 2):"
      DO i=1,nterms
      WRITE(1,'(a1,i3,a3)',ADVANCE ="NO")
     & '(',M_TERM(i),'); '	  
      ENDDO	 
      ENDIF	 	  
      WRITE(1,*)
      WRITE(1,*)	  
      WRITE(1,'(a33,i3)') "NUMBER OF R-GRID POINTS IS N_R = ",
     & n_r_coll
      ENDIF 	 
      ELSE
      IF(MYID.eq.0) THEN
      IF(.not.matrix_reading_defined)	  
     & WRITE(1,'(a97)') "COULPING MATRIX WILL BE COMPUTED USING NUMERICAL
     &  INTEGRATION OF POTENTIAL OVER THE WAVEFUNCTIONS"  
      ENDIF
      ENDIF	  
      CASE(2)
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a23,1x,a29)')"THE COLLISIONAL SYSTEM:",
     & "VIBRATING DIATOMIC TOP + ATOM"
      ENDIF	 
      IF(user_defined) THEN
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a7,2x,a60)') "WARNING:",
     & "THE WAVEFUNCTIONS AND ENERGY LEVELS MUST BE SUPPLIED BY USER"
      IF(rot_const_defined)WRITE(1,'(a27)')"ROT. CONST WILL BE INGNORED"
      IF(vib_const_defined)WRITE(1,'(a27)')"VIB. CONST WILL BE INGNORED"
      ENDIF	  
      INQUIRE( FILE=USR_INP_LEVEL, EXIST=exst )
      IF(.not.exst)	THEN
      IF(MYID.eq.0)PRINT*,USR_INP_LEVEL, " NOT FOUND"
      STOP	  
      ENDIF
	  
      IF(vib_mix_state_defined) THEN
      IF(MYID.eq.0) THEN
      OPEN(325,FILE=USR_INP_LEVEL,FORM="UNFORMATTED",
     & STATUS="OLD",ACTION="READ")
      IF(number_of_channels.le.0) THEN
      IF(MYID.eq.0) PRINT*,"NUMBER OF CHANNELS NOT SPECIFIED"
      STOP	  
      ENDIF
      IF(n_r_vib.le.0) THEN
      IF(MYID.eq.0) PRINT*,"INTEGRATION GRID NOT DEFINED"
      STOP	  
      ENDIF	  
	  
      IF(.not.channels_defined) THEN	  
      ALLOCATE(j_ch(number_of_channels))
      ALLOCATE(v_ch(number_of_channels))
      ENDIF
	  
      IF(.not.energy_defined) THEN	  
      ALLOCATE(E_ch(number_of_channels))
      ENDIF
	  
      READ(325)j_ch
      READ(325)v_ch
      READ(325)E_ch	  
      jmax_included = 0
      vmax_included = 0	  
      DO i=1,number_of_channels
      IF(jmax_included.lt.j_ch(i)) jmax_included= j_ch(i)
      IF(vmax_included.lt.v_ch(i)) vmax_included= v_ch(i)	  
      ENDDO
      ALLOCATE(
     & vibrational_wavefunctions(n_r_vib,
     & number_of_channels))

      ALLOCATE(
     & jacobian_vib(n_r_vib))

      ALLOCATE(r_grid_vib(n_r_vib))
	  
      READ(325)r_grid_vib
      READ(325)jacobian_vib	  
      READ(325) vibrational_wavefunctions
 
      CLOSE(325)	  
      ENDIF
	  
      IF(MYID.gt.0) THEN
      IF(.not.channels_defined) THEN	  
      ALLOCATE(j_ch(number_of_channels))
      ALLOCATE(v_ch(number_of_channels))
      ENDIF
	  
      IF(.not.energy_defined) THEN	  
      ALLOCATE(E_ch(number_of_channels))
      ENDIF	  
      ALLOCATE(
     & vibrational_wavefunctions(n_r_vib,
     & number_of_channels))

      ALLOCATE(
     & jacobian_vib(n_r_vib))

      ALLOCATE(r_grid_vib(n_r_vib))	  
	  !!!! TO BE CONTINUED
      ENDIF	
      channels_defined = .TRUE.
      energy_defined = .TRUE.
      CALL MPI_BCAST(j_ch,
     & number_of_channels, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(v_ch,
     & number_of_channels, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(j_ch,
     & number_of_channels, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(E_ch,
     & number_of_channels, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(r_grid_vib,
     & n_r_vib, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(jacobian_vib,
     & n_r_vib, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(vibrational_wavefunctions,
     & n_r_vib*number_of_channels, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	 
      ELSE 
      OPEN(325,FILE=USR_INP_LEVEL,
     & STATUS="OLD",ACTION="READ")
      IF(number_of_channels.le.0) THEN
      IF(MYID.eq.0) PRINT*,"NUMBER OF CHANNELS NOT SPECIFIED"
      STOP	  
      ENDIF
      IF(n_r_vib.le.0) THEN
      IF(MYID.eq.0) PRINT*,"INTEGRATION GRID NOT DEFINED"
      STOP	  
      ENDIF	  
      IF(.not.channels_defined) THEN	  
      ALLOCATE(j_ch(number_of_channels))
      ALLOCATE(v_ch(number_of_channels))

      ENDIF
      IF(.not.energy_defined) THEN	  
      ALLOCATE(E_ch(number_of_channels))
      ENDIF
      READ(325,*) !!! HEADER
      jmax_included = 0
      vmax_included = 0	  
      DO i=1,number_of_channels
      READ(325,'(i4,1x,i3,1x,i3,1x,e17.10)',
     & IOSTAT = stat_of_file)
     & decr,j_ch(i), v_ch(i), E_ch(i)
      IF(jmax_included.lt.j_ch(i)) jmax_included= j_ch(i)
      IF(vmax_included.lt.v_ch(i)) vmax_included= v_ch(i)	  
      IF(stat_of_file.gt.0 .or. i.ne.decr) THEN
      IF(MYID.eq.0) PRINT*,"WRONG CHANNELS IN USER INPUT FILE"	  
      STOP
      ENDIF 	  
      ENDDO
      ALLOCATE(
     & vibrational_wavefunctions(n_r_vib,
     & number_of_channels))

      ALLOCATE(
     & jacobian_vib(n_r_vib))

      ALLOCATE(r_grid_vib(n_r_vib))
	  
c      ! WRITTING A GRID	 
      READ(325,*) ! A HEADER
      DO i=1,n_r_vib! LOOP OVER VIBRATIONAL GRID
      READ(325,*)r_grid_vib(i),jacobian_vib(i)	  
!      READ(325,*)r_vb_dt_integrator1(i),wgg(i)
	  ! r_vb_dt_inegrator(i)-
! r-value  for i-point of the grid, wg(i) - weight(jacobian)
      ENDDO
      
      READ(325,*)! "WAVEFUNCTIONS ARE LISTED BELOW"	  
      DO k=1,number_of_channels !! k - channel number
c      !! WRITING THE WAVEFUNCTIONS	  
      READ(325,*)!	  "MOLECULE#1,CHANNEL=",k !!! A HEADER
      DO i=1,n_r_vib
      READ(325,*) vibrational_wavefunctions(i,k)	  
!      READ(325,*)wf_vib_part_diatom1(i,k) ! VALUE OF WAVEFUNCTION FOR k-channel at i-point of the grid
      ENDDO
      READ(325,*) !! LEAVE IT BLANK	  
      ENDDO
 
      CLOSE(325)
      ENDIF	  
      channels_defined = .TRUE.
      energy_defined = .TRUE.	
	   	  
      ELSE
      IF(BE.le.0) STOP "ERROR:SUPPLY POSITVE BE"
      IF(DE.GT.BE) STOP"ERROR: CHECK INPUT-DE SHOULD MUCH LESS THAN BE"	
      IF(morse_pot_defined) THEN
!      IF(morse_re.le.0)														!Dulat: This piece is moved down because Morse R needs atomic masses
!     & STOP "ERROR: YOU MUST SUPPLY R_EQ FOR MORSE POTENTIAL"
      IF(atomic_masses(1)*atomic_masses(2).le.0) STOP
     & "ERROR: SUPPLY MASS OF ATOMS"
      atomic_red_mass = atomic_masses(1)*atomic_masses(2)/
     & (atomic_masses(1)+atomic_masses(2))
      IF(MYID.eq.0) THEN	 
      WRITE(1,'(a31,f12.7)') "MOLECULAR REDUCED MASS,A.M.U = ",
     & atomic_red_mass
      ENDIF	 
      mass_vib_diatom = atomic_masses(1)+atomic_masses(2)	 

! Dulat Bostan start:
	  if(morse_re.le.0) then
	  if(myid.eq.0) print*, 
     & "Equilibrium bond distance is not specified in the input"
	  morse_re = dsqrt(1.d0/2.d0/(atomic_red_mass*amutoau)
     & /(Be/eVtown/autoeV))
	  if(myid.eq.0) print*, "Equilibrium bond distance: morse_re = ", 
     & morse_re, "a.u."
	  else
	  if(myid.eq.0) print*, "Equilibrium bond distance: morse_re = ", 
     & morse_re, "a.u."
	  endif      
! Dulat Bostan end	  
	  
      IF(morse_a.gt.0 .and. morse_de.gt.0 .and. .not. vib_const_defined
     & ) THEN
!      WE = 1d0/morse_a															!Dulat: disabled this line and corrected equation in the next line
      WE = morse_a
     & /dsqrt(atomic_red_mass*amutoau/2d0/morse_de/eVtown/autoeV)
      XE = WE**2/4d0/morse_de
      ENDIF
      IF(r_vib_diatom_max.le.0d0 .or. r_vib_diatom_min.lt.0d0) STOP
     & "ERROR: SUPPLY MORSE POTENTIAL RANGE."
      IF(MYID.eq.0) THEN	 
      WRITE(1,'(a33,f10.4)')
     & "RMIN FOR MORSE POTENTIAL, a.u. = ",r_vib_diatom_min	 
      WRITE(1,'(a33,f10.4)')
     & "RMAX FOR MORSE POTENTIAL, a.u. = ",r_vib_diatom_max
      ENDIF	 
      ENDIF	  
      IF(WE.le.0 ) STOP "ERROR:SUPPLY POSITVE WE"
      IF(XE.le.0 .or.XE.ge.WE ) STOP "ERROR:SUPPLY POSITVE XE<WE"
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a32,f10.4)') "ROTATIONAL CONSTANT, CM-1, BE = ",Be	  
      WRITE(1,'(a23,f12.7)') "DISTORTION, CM-1, DE = ",DE	 
      WRITE(1,'(a33,f10.4)') "VIBRATIONAL CONSTANT, CM-1, WE = ",WE	  
      WRITE(1,'(a24,f12.7)') "ANHARMONICITY, CM-1, XE = ",XE
      ENDIF	  
      IF(emax_defined .and. .not.channels_defined) THEN
      nchann_est = (int(EMAX/WE)+2)*(int(DSQRT(EMAX/BE))+2)
      ALLOCATE(j_ch_tmp(nchann_est),v_ch_temp(nchann_est),
     & e_temp(nchann_est),ja(nchann_est))
      n_prev = 0	 
      DO v_t = 0,int(EMAX/WE)+1
      DO j_t = 0,int(DSQRT(EMAX/BE))+1
      n_prev = n_prev + 1	  
      j_ch_tmp(n_prev) = j_t
      v_ch_temp(n_prev) = v_t
      CALL energy_diatom(e_temp(n_prev),v_t,j_t,we,xe,be,de)   	  
      ENDDO
      ENDDO
      CALL	hpsort(nchann_est,e_temp,ja)
      nlvl = 0	  
      DO WHILE(.not.channels_defined)
      nlvl = nlvl + 1
      IF(e_temp(nlvl).gt.EMAX) THEN
      number_of_channels = 
     & nlvl - 1
      channels_defined = .TRUE.	 
      ENDIF	 
      ENDDO
      IF(number_of_channels.le.0)STOP"ERROR:ZERO NUMBER OF CHANNELS"
      IF(ALLOCATED(j_ch)) DEALLOCATE (j_ch)
      IF(ALLOCATED(v_ch)) DEALLOCATE (v_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE (E_ch)	  
      ALLOCATE(j_ch(number_of_channels))
      ALLOCATE(v_ch(number_of_channels))	  
      ALLOCATE(E_ch(number_of_channels)) 	  
      ALLOCATE(j_ch(number_of_channels),
     & v_ch(number_of_channels),E_ch(number_of_channels))
      DO nlvl = 1,number_of_channels
      n_prev = ja(nlvl)	  
      j_ch(nlvl) = j_ch_tmp(n_prev)
      v_ch(nlvl) = v_ch_temp(n_prev)
      E_ch(nlvl) = e_temp(nlvl)
      ENDDO
      DEALLOCATE(j_ch_tmp,v_ch_temp,e_temp,ja)
      energy_defined = .TRUE.
      channels_defined = .TRUE.	  
      ENDIF
  
      IF(jmax.ge.0 .and. jmin.ge.0 .and. vmax.ge.0 .and. vmin.ge.0)
     & THEN
      IF(jmax.lt.jmin) STOP "ERROR: JMAX MUST BE > OR = THAN JMIN"
      IF(vmax.lt.vmin) STOP "ERROR: VMAX MUST BE > OR = THAN VMIN"	  
      IF(number_of_channels_defined) 
     & STOP "ERROR: NUMBER_OF_CHANNELS ALREADY DEFINED"
      IF(energy_defined) 
     & STOP "ERROR: ENRGIES ALREADY DEFINED"	 
      number_of_channels_defined = .TRUE.
      energy_defined = .TRUE.
      number_of_channels = (jmax - jmin + 1)*(vmax - vmin + 1)
      IF(ALLOCATED(j_ch)) DEALLOCATE (j_ch)
      IF(ALLOCATED(v_ch)) DEALLOCATE (v_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE (E_ch)	  
      ALLOCATE(j_ch(number_of_channels))
      ALLOCATE(v_ch(number_of_channels))	  
      ALLOCATE(E_ch(number_of_channels)) 	  
      i =  0 	  
      DO v_t=vmin,vmax
      DO j_t=jmin,jmax
      i = i + 1	  
      j_ch(i) = j_t
      v_ch(i) = v_t	  
      CALL energy_diatom(E_ch(i),v_t,j_t,we,xe,be,de)
      ENDDO
      ENDDO
      IF(i.ne.number_of_channels) STOP"ERROR:PROGRAM ERROR"	
      	  

      ELSE
      IF(channels_defined .and. .not. energy_defined) THEN 	  
      IF(number_of_channels.le.0) 
     & STOP"ERROR:NUMBER_OF_CHANNLES MUST BE >0"
      IF(.not.energy_defined) THEN	 
      ALLOCATE(E_ch(number_of_channels))	 
      DO i=1,number_of_channels	  
      j_t = j_ch(i)
      v_t = v_ch(i)
      IF(jmax_included.lt.j_t)	jmax_included = j_t
      IF(vmax.lt.v_t) vmax = v_t  
      CALL energy_diatom(E_ch(i),v_t,j_t,we,xe,be,de)
      ENDDO
      energy_defined = .TRUE.	  
      ENDIF	  
      ENDIF
      IF(number_of_channels_defined .and. .not.channels_defined) THEN
      nchann_est = (number_of_channels+1)**2
      n_prev = 0	 
      DO v_t = 0,number_of_channels
      IF(n_prev.gt.nchann_est) EXIT	  
      DO j_t = 0,number_of_channels
      IF(n_prev.gt.nchann_est) EXIT		  
      n_prev = n_prev + 1	  
      j_ch_tmp(n_prev) = j_t
      v_ch_temp(n_prev) = v_t
      CALL energy_diatom(e_temp(n_prev),v_t,j_t,we,xe,be,de)   	  
      ENDDO
      ENDDO
      CALL	hpsort(nchann_est,e_temp,ja)
      IF(ALLOCATED(j_ch)) DEALLOCATE (j_ch)
      IF(ALLOCATED(v_ch)) DEALLOCATE (v_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE (E_ch)	  
      ALLOCATE(j_ch(number_of_channels))
      ALLOCATE(v_ch(number_of_channels))	  
      ALLOCATE(E_ch(number_of_channels)) 	  
      ALLOCATE(j_ch(number_of_channels),
     & v_ch(number_of_channels),E_ch(number_of_channels))
      DO nlvl = 1,number_of_channels
      n_prev = ja(nlvl)	  
      j_ch(nlvl) = j_ch_tmp(n_prev)
      v_ch(nlvl) = v_ch_temp(n_prev)
      E_ch(nlvl) = e_temp(nlvl)
      ENDDO
      DEALLOCATE(j_ch_tmp,v_ch_temp,e_temp,ja)
      energy_defined = .TRUE.
      channels_defined = .TRUE.		  
      ENDIF	  
	  
	  
      ENDIF	  
      ENDIF

      IF(atom_coord_dist) THEN
      IF(MYID.eq.0) THEN	  
      WRITE(1,*)	  
      WRITE(1,'(a28)') "PAIR WISE POTENTIAL DEFINED."
      WRITE(1,'(a47,f10.4,1x,a10,f10.4)')
     & "MASS OF ATOMS IN THE MOLECULE, A.M.U.: MASS1 = ",
     & atomic_masses(1),", MASS2 = ",atomic_masses(2)
      WRITE(1,*)
      ENDIF	  
      ENDIF	
      IF(para_ortho_defined) CALL SYMMETRY_SORT_PARA_ORTHO	  
      IF(MYID.eq.0)WRITE(1,*)
      IF(channels_defined.and.energy_defined)
     & THEN
      IF(MYID.eq.0) THEN	 
      WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"	  
      WRITE(1,'(a5,9x,a2,9x,a2,9x,a2)')"#N = "
     & ,"J","V","E"
      ENDIF	 
      DO i=1,	number_of_channels
      IF(jmax_included.lt.j_ch(i)) jmax_included = j_ch(i)
      IF(vmax_included.lt.v_ch(i)) vmax_included = v_ch(i)	  
      IF(MYID.eq.0)WRITE(1,'(i5,i11,i11,f11.3)')i,
     & j_ch(i),v_ch(i),E_ch(i)
      ENDDO	  
      IF(MYID.eq.0)WRITE(1,*)
      IF(ini_chann_defined .and. .not.all_level_included) THEN
      chann_ini	 = 0  
      DO i=1,number_of_channels
      IF(j_ini.eq.j_ch(i) .and. v_ini.eq.v_ch(i) ) chann_ini = i
      ENDDO	
      IF(chann_ini.eq.0) STOP "ERROR: INITIAL CHANNEL WRONG"
      IF(MYID.eq.0)WRITE(1,'(a17,a4,i4,a6,i4)') "INITIAL STATE IS:",
     & "J = ",j_ini,"; V = ", v_ini
      ENDIF
      ENDIF
      IF(.not.channels_defined.or. .not.energy_defined) 
     & STOP "ERROR IN INPUT: ENERGIES AND LEVELS ARE NOT DEFINED"
      IF(.not.ini_chann_defined .and. channels_defined .and.
     & energy_defined) THEN
	  i_ground = 1
      DO i=1,number_of_channels
      IF(E_ch(i_ground).gt.E_ch(i)) i_ground = i   
      ENDDO	 
      chann_ini = i_ground
      ini_chann_defined = .TRUE.	  
      ENDIF	 
	  IF(vib_mix_state_defined) THEN
      INQUIRE( FILE=VIB_MIX_STATE, EXIST=exst )
      IF(.not.exst)	THEN
      IF(MYID.eq.0)PRINT*,VIB_MIX_STATE, " NOT FOUND"
      STOP	  
      ENDIF	
      ALLOCATE(vib_mix_real(number_of_channels),
     & vib_mix_imag(number_of_channels))
      OPEN(22,FILE=VIB_MIX_STATE,ACTION="READ",STATUS="OLD")
      READ(22,*)	  
	  DO i=1,number_of_channels
      READ(22,*) k,	vib_mix_real(i),vib_mix_imag(i)  
      ENDDO
      CLOSE(22)
      IF(MYID.eq.0) PRINT*,"MIXED_INI_STATE INITILIALZID"	  
      ENDIF	  
      IF(MYID.eq.0) THEN	 
      WRITE(1,*)  	  
      IF(angs_unit)WRITE(1,'(a35)')"THE POTENTIAL R-UNITS ARE ANGSTROMS"
      IF(au_unit_r)WRITE(1,'(a31)') "THE POTENTIAL R-UNITS ARE BOHRS"
      IF(cm_unit)WRITE(1,'(a31)')"THE POTENTIAL ENERGY IS IN CM-1"
      IF(au_unit_e)WRITE(1,'(a34)')"THE POTENTIAL ENERGY IS IN HARTREE"
      IF(klv_unit)WRITE(1,'(a33)')"THE POTENTIAL ENERGY IS IN KELVIN"
      IF(kal_unit)WRITE(1,'(a35)')"THE POTENTIAL ENERGY IS IN KCAL/MOL"	  
      ENDIF
      IF(.not.angs_unit .and..not.au_unit_r )   
     & STOP "ERROR:R_UNITS NOT DEFINED"
      IF(.not.cm_unit .and..not.au_unit_e.and..not.klv_unit
     & .and. .not.kal_unit)   
     & STOP "ERROR:ENERGY NOT DEFINED"
      IF(MYID.eq.0) WRITE(1,*) 	 
      IF(expansion_defined) THEN
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a42)') "USER SUPPLIED POTENTIAL TERMS WILL BE USED"
      WRITE(1,*)
      WRITE(1,'(a45,i4)')
     & "NUMBER OF POTENTIAL TERMS INCLUDED N_TERMS = ",nterms
      WRITE(1,'(a38)')
     & "THE POTENTIAL TERMS ARE, (LAMBDA):"
      DO i=1,nterms
      WRITE(1,'(a1,i3,a3)',ADVANCE ="NO")
     & '(',L_TERM(i),'); '	  
      ENDDO	
      WRITE(1,*)
      WRITE(1,*)	  
      WRITE(1,'(a33,i3)') "NUMBER OF R-GRID POINTS IS N_R = ",
     & n_r_coll
      WRITE(1,'(a39,i3)') "NUMBER OF R_VIB-GRID POINTS IS N_VIB = ",
     & n_r_vib
      ENDIF	 
      ELSE
      IF(MYID.eq.0) THEN	  
      IF(.not.matrix_reading_defined)	  
     & WRITE(1,'(a97)')"COULPING MATRIX WILL BE COMPUTED USING NUMERICAL
     &  INTEGRATION OF POTENTIAL OVER THE WAVEFUNCTIONS"
      ENDIF
      ENDIF	  
      CASE(3)
      IF(MYID.eq.0)	  
     & WRITE(1,'(a23,1x,a26)')"THE COLLISIONAL SYSTEM:",
     & "RIGID SYMMETRIC TOP + ATOM"
      IF(user_defined) THEN
      IF(MYID.eq.0)	  
     & WRITE(1,'(a7,2x,a60)') "WARNING:",
     & "THE WAVEFUNCTIONS AND ENERGY LEVELS MUST BE SUPPLIED BY USER"
      IF(rot_const_defined)WRITE(1,'(a27)')"ROT. CONST WILL BE INGNORED"
      INQUIRE( FILE=USR_INP_LEVEL, EXIST=exst )
      IF(.not.exst)	THEN
      IF(MYID.eq.0)PRINT*,USR_INP_LEVEL, " NOT FOUND"
      STOP	  
      ENDIF	  
      OPEN(325,FILE=USR_INP_LEVEL,
     & STATUS="OLD",ACTION="READ")
      IF(number_of_channels.le.0) THEN
      PRINT*,"NUMBER OF CHANNELS NOT SPECIFIED"
      STOP	  
      ENDIF
      IF(.not.channels_defined) THEN	  
      ALLOCATE(j_ch(number_of_channels))
      ALLOCATE(k_ch(number_of_channels))
      ALLOCATE(eps_ch(number_of_channels))	  
      ENDIF
      IF(.not.energy_defined) THEN	  
      ALLOCATE(E_ch(number_of_channels))
      ENDIF
      READ(325,*) !HEADER
      jmax_included = 0	  
      DO i=1,number_of_channels
      READ(325,'(i5,1x,i4,1x,i4,1x,i2,1x,e17.10)',
     & IOSTAT = stat_of_file)
     & decr,j_ch(i), k_ch(i),eps_ch(i),E_ch(i)
      IF(jmax_included.lt.j_ch(i))jmax_included=j_ch(i)	 
      IF(stat_of_file.gt.0 .or. decr.ne.i) THEN
      PRINT*,"WRONG CHANNELS IN USER INPUT FILE"	  
      STOP
      ENDIF 	  
      ENDDO
      CLOSE(325)
      channels_defined = .TRUE.
      energy_defined = .TRUE.	  	  
      ELSE
      IF(A+B.eq.0) STOP "ERROR: SUPPLY A OR B"
      IF(A.eq.0) A=B
      B=A
      IF(C.eq.0) STOP "ERROR : SYPPLY C"
       IF(MYID.eq.0)	  
     & WRITE(1,'(a40,f7.2,a6,f7.2,a2)')
     & "ROTATIONAL CONSTANTS ARE, CM-1: A = B = ",A
     &, "; C = ", C,";"
      IF(number_of_channels_defined) THEN
      IF(energy_defined.and. .not.channels_defined) STOP 
     & "ERROR:FOR MANUALLY DEFINED ENERGY YOU SHOULD SPECIFY CHANNELS"	  
      IF(.not.energy_defined) ALLOCATE(E_ch(number_of_channels))	  
      IF(.not.channels_defined) THEN
c      PRINT*, "CHANNELS NOT DEFIEND = ",number_of_channels
      ALLOCATE(j_ch(number_of_channels),k_ch(number_of_channels),
     & eps_ch(number_of_channels))
      j_ch = 0
      k_ch = 0
      eps_ch = 0
      decr = 1
      IF(max_nmb_chann.lt.number_of_channels) STOP "TOO MANY CHANNELS"	  
      nlvl = number_of_channels*decr	  
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
c      PRINT*,"nlvl==",nlvl	  
      ALLOCATE(j_ch_tmp(nlvl),
     & k_ch_tmp(nlvl),ja(nlvl),
     & eps_ch_tmp(nlvl),e_temp(nlvl))
      i = 0	 
      DO j_t = 0,j_max_allowed
      DO k_t = 0,j_t
      DO eps_t = 0,1-KRONEKER(k_t,0)
      i = i+1    
      IF(i.gt.nlvl) EXIT
      j_ch_tmp(i) = j_t
      k_ch_tmp(i) = k_t
      eps_ch_tmp(i) = eps_t
      e_temp(i) = A*(j_t+1d0)*j_t - (A-C)*k_t**2 
      ENDDO
      IF(i.gt.nlvl) EXIT	  
      ENDDO	  
      IF(i.gt.nlvl) EXIT	  
      ENDDO
!      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
!      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,number_of_channels	  
      IF(E_ch(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,number_of_channels
      j_ch(i) = j_ch_tmp(ja(i))
      k_ch(i) = k_ch_tmp(ja(i))
      eps_ch(i) = eps_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = number_of_channels*decr      	  
      DEALLOCATE(k_ch_tmp,j_ch_tmp,ja,e_temp,eps_ch_tmp)	  
      ENDDO 
      energy_defined = .TRUE.	  
      ENDIF	  
      jmax_included = 0
      IF(MYID.eq.0)	  
     & WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"	  
      IF(MYID.eq.0)	  
     & WRITE(1,'(a5,9x,a2,9x,a2,9x,a2,9x,a2)')"#N = "
     & ,"J","K","P", "E"	  
      DO i=1,number_of_channels
      IF(j_ch(i).lt.k_ch(i))
     & STOP
     & "ERROR: J MUST BE NO LESS THAN K. CHECK THE INPUT CHANNELS."
      IF(k_ch(i).eq.0 .and. eps_ch(i).ne.0)
     &	  STOP"ERROR:K=0, SO PARITY MUST BE 0"
      IF(abs(eps_ch(i)).gt.1) STOP "ERROR: PARITY MUST BE -1 OR 1"
      IF(jmax_included.lt.j_ch(i)) jmax_included = j_ch(i)
      ENDDO	  
      ALLOCATE(indch_ABC(2,2*jmax_included+1,jmax_included+1))
      indch_ABC = 0
      DO i=1,number_of_channels
      IF(indch_ABC(eps_ch(i)+1,j_ch(i)+k_ch(i)+1,j_ch(i)+1).eq.0)
     & THEN
      indch_ABC(eps_ch(i)+1,j_ch(i)+k_ch(i)+1,j_ch(i)+1) = i 	 
      ELSE
      IF(MYID.eq.0)	  
     & WRITE(1,'(a42,1x,a4,i3,a4,1x,i3,a9,1x,i3)')
     & "ERROR: YOU ENTERED TWICE THE SAME CHANNEL:"," j =",
     & j_ch(i), " k =",k_ch(i)," PARITY =", eps_ch(i)
      STOP "ERROR: CHECK INPUT. ENTERED THE SAME CHANNEL TWICE"	 
      ENDIF
      IF(.not.energy_defined) 
     & E_ch(i) = A*(j_ch(i)+1d0)*j_ch(i) - (A-C)*k_ch(i)**2 
      ENDDO
      energy_defined = .TRUE.	  
      ENDIF
	  
	  
	  
      IF(emax_defined.and. .not.channels_defined) THEN
c      PRINT*,"EMAX=",EMAX	  
      j_t = int(abs(min(sqrt(EMAX/C),(EMAX/A-1)/2d0)))+1
      nchnl_ini_guess = (j_t+1)**2	  
c      PRINT*,"j_ini",j_t	  
      nchann_est = nchnl_ini_guess
      n_prev = -1
      nlvl = 0
      decr = 1
      DO WHILE(nlvl.ne.n_prev)
      n_prev =  nlvl	  
      nchann_est = nchnl_ini_guess*decr
c      PRINT*,"nchann_est = ",nchann_est	  
      ALLOCATE(j_ch_tmp(nchann_est),k_ch_tmp(nchann_est),
     & eps_ch_tmp(nchann_est),e_temp(nchann_est),ja(nchann_est))
      nlvl = 0
      DO i=0,j_max_allowed
      DO k_t=0,i
      DO eps_t = 0,1-KRONEKER(k_t,0)
      nlvl = nlvl + 1
      IF(nlvl.gt.nchann_est) EXIT	  
      j_ch_tmp(nlvl) = i
      k_ch_tmp(nlvl) = k_t
      eps_ch_tmp(nlvl) = eps_t
      e_temp(nlvl) = A*(i+1d0)*i - (A-C)*k_t**2 	  
      ENDDO
      IF(nlvl.gt.nchann_est) EXIT		  
      ENDDO
      IF(nlvl.gt.nchann_est) EXIT		  
      ENDDO
c      PRINT*,"nlvl",nlvl	  
      CALL	hpsort(nchann_est,e_temp,ja)
      i = 1	  
      DO WHILE(e_temp(i).lt.EMAX)
      i = i + 1
c      PRINT*,"E",i,e_temp(i)	  
      ENDDO
      nlvl = i - 1	  
c      PRINT*,'nlvl',nlvl      	  
      decr = decr + 1 
      DEALLOCATE(j_ch_tmp,k_ch_tmp,
     & eps_ch_tmp,e_temp,ja)	  
      ENDDO
       number_of_channels = nlvl
      ALLOCATE(j_ch(number_of_channels),k_ch(number_of_channels),
     & eps_ch(number_of_channels),E_ch(number_of_channels)) 
       ALLOCATE(j_ch_tmp(nchann_est),k_ch_tmp(nchann_est),
     & eps_ch_tmp(nchann_est),e_temp(nchann_est),ja(nchann_est))
	 
      nlvl = 0
      DO i=0,j_max_allowed
      DO k_t=0,i
      DO eps_t = 0,1-KRONEKER(k_t,0)
      nlvl = nlvl + 1
      IF(nlvl.gt.nchann_est) EXIT	  
      j_ch_tmp(nlvl) = i
      k_ch_tmp(nlvl) = k_t
      eps_ch_tmp(nlvl) = eps_t
      e_temp(nlvl) = A*(i+1d0)*i - (A-C)*k_t**2 	  
      ENDDO
      IF(nlvl.gt.nchann_est) EXIT		  
      ENDDO
      IF(nlvl.gt.nchann_est) EXIT		  
      ENDDO
c      PRINT*,"number_of_channels = ",number_of_channels     
      CALL	hpsort(nchann_est,e_temp,ja)	  
      DO i=1,number_of_channels
      j_ch(i) = j_ch_tmp(ja(i))
      k_ch(i) = k_ch_tmp(ja(i))
      eps_ch(i) = eps_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      DEALLOCATE(j_ch_tmp,k_ch_tmp,
     & eps_ch_tmp,e_temp,ja)
      energy_defined = .TRUE.
      IF(MYID.eq.0)	  
     & WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"	  
       IF(MYID.eq.0)	  
     & WRITE(1,'(a5,9x,a2,9x,a2,9x,a2,9x,a2)')"#N = "
     & ,"J","K","P", "E"
      DO i=1,number_of_channels	 
      IF(.not.energy_defined) 
     & E_ch(i) = A*(j_ch(i)+1d0)*j_ch(i) - (A-C)*k_ch(i)**2 
      ENDDO	 
      ENDIF		  
      ENDIF
      IF(para_ortho_defined) CALL SYMMETRY_SORT_PARA_ORTHO
      DO i=1,number_of_channels	 
      IF(MYID.eq.0)	  
     & WRITE(1,'(i5,i11,i11,i11,f11.3)')i,
     & j_ch(i),k_ch(i),eps_ch(i),E_ch(i)
      ENDDO	 	 
      IF(ini_chann_defined .and. .not.all_level_included) THEN
      chann_ini = 0
      jmax_included = 0	  
      DO i=1,number_of_channels
      IF(j_ini.eq.j_ch(i) .and. k_ini.eq.k_ch(i) .and. 
     & eps_ini.eq.eps_ch(i)) chann_ini = i
      IF(jmax_included.lt.j_ch(i)) jmax_included = j_ch(i)	 
      ENDDO	
      IF(chann_ini.eq.0) STOP "ERROR: INITIAL CHANNEL WRONG"
       IF(MYID.eq.0)	  
     & WRITE(1,'(a17,a4,i4,a6,i4,a6,i4)') "INITIAL STATE IS:",
     & "J = ",j_ini,"; K = ", k_ini,"; P =	",eps_ini  
      ENDIF
      IF(.not.ini_chann_defined .and. channels_defined .and.
     & energy_defined) THEN
	  i_ground = 1
      DO i=1,number_of_channels
      IF(E_ch(i_ground).gt.E_ch(i)) i_ground = i   
      ENDDO	 
      chann_ini = i_ground
      ini_chann_defined = .TRUE.	  
      ENDIF	  
      IF(MYID.EQ.0) THEN	  
      WRITE(1,*)
      IF(angs_unit)WRITE(1,'(a35)')"THE POTENTIAL R-UNITS ARE ANGSTROMS"
      IF(au_unit_r)WRITE(1,'(a31)') "THE POTENTIAL R-UNITS ARE BOHRS"
      IF(cm_unit)WRITE(1,'(a31)')"THE POTENTIAL ENERGY IS IN CM-1"
      IF(au_unit_e)WRITE(1,'(a34)')"THE POTENTIAL ENERGY IS IN HARTREE"
      IF(klv_unit)WRITE(1,'(a33)')"THE POTENTIAL ENERGY IS IN KELVIN"
      IF(kal_unit)WRITE(1,'(a35)')"THE POTENTIAL ENERGY IS IN KCAL/MOL"
      ENDIF	  
      IF(.not.angs_unit .and..not.au_unit_r )   
     & STOP "ERROR:R_UNITS NOT DEFINED"
      IF(.not.cm_unit .and..not.au_unit_e.and..not.klv_unit
     & .and. .not.kal_unit)   
     & STOP "ERROR:ENERGY NOT DEFINED"
      IF(MYID.eq.0)	  
     & WRITE(1,*) 	 
      IF(expansion_defined) THEN
      IF(myid.eq.0) THEN	  
      WRITE(1,'(a42)') "USER SUPPLIED POTENTIAL TERMS WILL BE USED"
      WRITE(1,*)
      WRITE(1,'(a45,i4)')
     & "NUMBER OF POTENTIAL TERMS INCLUDED N_TERMS = ",nterms
      WRITE(1,'(a38)')
     & "THE POTENTIAL TERMS ARE, (LAMBDA,NJU):"
      DO i=1,nterms
      WRITE(1,'(a1,i3,a1,i3,a3)',ADVANCE ="NO")
     & '(',L_TERM(i),',',M_TERM(i),'); '	  
      ENDDO	
      WRITE(1,*)
      WRITE(1,*)	  
      WRITE(1,'(a33,i3)') "NUMBER OF R-GRID POINTS IS N_R = ",
     & n_r_coll
      ENDIF	 
      ELSE
      IF(.not.matrix_reading_defined .and. myid.eq.0)	  
     & WRITE(1,'(a97)')"COULPING MATRIX WILL BE COMPUTED USING NUMERICAL
     &  INTEGRATION OF POTENTIAL OVER THE WAVEFUNCTIONS"	  
      ENDIF
       IF(MYID.eq.0)	  
     & WRITE(1,*)
      	  
      CASE(4)
      IF(MYID.eq.0) 	  
     & WRITE(1,'(a23,1x,a27)')"THE COLLISIONAL SYSTEM:",
     & "RIGID ASYMMETRIC TOP + ATOM"
      IF(user_defined) THEN
      IF(MYID.eq.0)	  
     & WRITE(1,'(a7,2x,a60)') "WARNING:",
     & "THE WAVEFUNCTIONS AND ENERGY LEVELS MUST BE SUPPLIED BY USER"
      IF(rot_const_defined)WRITE(1,'(a27)')"ROT. CONST WILL BE INGNORED"
      INQUIRE( FILE=USR_INP_LEVEL, EXIST=exst )
      IF(.not.exst)	THEN
      IF(MYID.eq.0)PRINT*,USR_INP_LEVEL, " NOT FOUND"
      STOP	  
      ENDIF	 
	  OPEN(325,FILE=USR_INP_LEVEL,
     & STATUS="OLD",ACTION="READ")
      IF(number_of_channels.le.0) THEN
      PRINT*,"NUMBER OF CHANNELS NOT SPECIFIED"
      STOP	  
      ENDIF
      IF(.not.channels_defined) THEN	  
      ALLOCATE(j_ch(number_of_channels))
      ALLOCATE(k_ch(number_of_channels))
      ALLOCATE(eps1_ch(number_of_channels))	  
      ENDIF
      IF(.not.energy_defined) THEN	  
      ALLOCATE(E_ch(number_of_channels))
      ENDIF
	  READ(325,*) !HEADER
      jmax_included = 0	  
      DO i=1,number_of_channels
      READ(325,*,IOSTAT = stat_of_file)
     & decr,j_ch(i), ka_ch(i),kc_ch(i),E_ch(i)
      IF(jmax_included.lt.j_ch(i)) jmax_included= j_ch(i)
      IF(stat_of_file.gt.0 .or. decr.ne.i) THEN
      PRINT*,"WRONG CHANNELS IN USER INPUT FILE"	  
      STOP
      ENDIF 	  
      ENDDO
      READ(325,*) !! HEADER FOR BASIS
      ALLOCATE(M_EIGEN_VECTORS_USER
     & (number_of_channels,jmax_included*2+1))
      M_EIGEN_VECTORS_USER = 0d0	 
      DO i=1,number_of_channels
      READ(325,*) decr,M_EIGEN_VECTORS_USER
     & (i,1:1+2*j_ch(i))
      IF(i.ne.decr .or. stat_of_file.gt.0) STOP "ERROR IN USER FILE"	 
      ENDDO	  
      CLOSE(325)
      channels_defined = .TRUE.
      energy_defined = .TRUE.	  	  	  
	  
      ELSE
      IF(A*B*C.le.0d0)
     & STOP"ERROR: ALL A,B AND C MUST BE POSITIVE. CHECK INPUT"
       IF(MYID.eq.0) 	  
     & WRITE(1,'(a35,f7.2,a7,f7.2,a6,f7.2,a2)')
     & "ROTATIONAL CONSTANTS ARE, CM-1: A =",A," ; B = ",
     &  B,"; C = ", C,";"	 
      IF(number_of_channels_defined) THEN
      IF(energy_defined.and. .not.channels_defined) STOP 
     & "ERROR:FOR MANUALLY DEFINED ENERGY YOU SHOULD SPECIFY CHANNELS"
      IF(.not.energy_defined)  ALLOCATE(E_ch(number_of_channels))
      IF(.not.channels_defined) THEN
      ALLOCATE(j_ch(number_of_channels),ka_ch(number_of_channels),
     & kc_ch(number_of_channels))
      j_ch = 0
      ka_ch = 0
      kc_ch = 0
      decr = 1
      IF(max_nmb_chann.lt.number_of_channels) STOP "TOO MANY CHANNELS"	  
      nlvl = number_of_channels*decr	  
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
c      PRINT*,"nlvl==",nlvl	  
      ALLOCATE(j_ch_tmp(nlvl),
     & ka_ch_tmp(nlvl),ja(nlvl),
     & kc_ch_tmp(nlvl),e_temp(nlvl))
      i = 0	 
      DO j_t = 0,j_max_allowed
      DO dk_t = -j_t,j_t
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"	  
      i = i+1    
      IF(i.gt.nlvl) EXIT
      j_ch_tmp(i) = j_t
      ka_ch_tmp(i) = ka_t
      kc_ch_tmp(i) = kc_t
c      PRINT*,"j_t,ka_t,kc_t="	, j_t, ka_t,kc_t	  
      CALL EIGEN_VAL(A,B,C,j_t,ka_t,kc_t,e_temp(i))
c      PRINT*,"e_temp,="	, e_temp(i) 
      ENDDO
      IF(i.gt.nlvl) EXIT	  
      ENDDO
c      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
c      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,number_of_channels	  
      IF(E_ch(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,number_of_channels
      j_ch(i) = j_ch_tmp(ja(i))
      ka_ch(i) = ka_ch_tmp(ja(i))
      kc_ch(i) = kc_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = number_of_channels*decr      	  
      DEALLOCATE(ka_ch_tmp,j_ch_tmp,ja,e_temp,kc_ch_tmp)	  
      ENDDO 
      energy_defined = .TRUE.	  
      ENDIF
      jmax_included = 0
       IF(MYID.eq.0) 	  
     & WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"	  
       IF(MYID.eq.0) 	  
     & WRITE(1,'(a5,9x,a2,9x,a2,9x,a2,9x,a2)')"#N = "
     & ,"J","KA","KC", "E"	  
      DO i=1,number_of_channels
      IF(j_ch(i).lt.max(ka_ch(i),kc_ch(i)))
     & STOP
     & "ERROR: J MUST BE NO LESS THAN K. CHECK THE INPUT CHANNELS."
      IF(.not. EXST_STATE(j_ch(i),ka_ch(i),kc_ch(i))) THEN 
       IF(MYID.eq.0) 	  
     & WRITE(1,'(a48,1x,i4,1x,i4,1x,i4)')
     & "ERROR: CHECK CHANNELS. THIS STATE DOES NOT EXIST",
     & j_ch(i),ka_ch(i),kc_ch(i)
      STOP
      ENDIF
      IF(jmax_included.lt.j_ch(i)) jmax_included = j_ch(i)
      ENDDO	  
      ALLOCATE(indch_ABC
     & (2*jmax_included+1,2*jmax_included+1,jmax_included+1))
      indch_ABC = 0
      DO i=1,number_of_channels
      IF(indch_ABC(j_ch(i)+kc_ch(i)+1,j_ch(i)+ka_ch(i)+1,j_ch(i)+1) 
     & .eq.0)
     & THEN
      indch_ABC(j_ch(i)+kc_ch(i)+1,j_ch(i)+ka_ch(i)+1,j_ch(i)+1) = i 	 
      ELSE
       IF(MYID.eq.0) 	  
     & WRITE(1,'(a42,1x,a5,i3,a4,1x,i3,a5,1x,i3)')
     & "ERROR: YOU ENTERED TWICE THE SAME CHANNEL:"," j =",
     & j_ch(i), " ka =",ka_ch(i)," kc =", kc_ch(i)
      STOP "ERROR: CHECK INPUT. ENTERED THE SAME CHANNEL TWICE"	 
      ENDIF
      IF(.not.energy_defined) THEN
      CALL EIGEN_VAL(A,B,C,j_ch(i),ka_ch(i),kc_ch(i),E_ch(i))	  
      ENDIF	 
       IF(MYID.eq.0) 	  
     & WRITE(1,'(i5,i11,i11,i11,f11.3)')i,
     & j_ch(i),ka_ch(i),kc_ch(i),E_ch(i) 	  
      ENDDO
      energy_defined = .TRUE.
      channels_defined = .TRUE.	 
      ENDIF
 
      IF(emax_defined.and. .not.channels_defined) THEN
c      PRINT*,"EMAX=",EMAX	  
      j_t = int(abs(min(sqrt(EMAX/C),(EMAX/A-1)/2d0)))+1
c      PRINT*,"j_ini",j_t
      nchnl_ini_guess =(j_t+1)**2  	  
      nchann_est = nchnl_ini_guess
      n_prev = -1
      nlvl = 0
      decr = 1
      DO WHILE(nlvl.ne.n_prev)
      n_prev =  nlvl	  
      nchann_est = nchnl_ini_guess*decr
c      PRINT*,"nchann_est = ",nchann_est	  
      ALLOCATE(j_ch_tmp(nchann_est),ka_ch_tmp(nchann_est),
     & kc_ch_tmp(nchann_est),e_temp(nchann_est),ja(nchann_est))
      nlvl = 0
      DO j_t=0,j_max_allowed
      DO dk_t=-j_t,j_t
      nlvl = nlvl + 1
      IF(nlvl.gt.nchann_est) EXIT	  
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"
      j_ch_tmp(nlvl) = j_t
      ka_ch_tmp(nlvl) = ka_t
      kc_ch_tmp(nlvl) = kc_t
      CALL EIGEN_VAL(A,B,C,j_t,ka_t,kc_t,e_temp(nlvl))

	  
      ENDDO
      IF(nlvl.gt.nchann_est) EXIT		  
      ENDDO
c      PRINT*,"nlvl",nlvl	  
      CALL	hpsort(nchann_est,e_temp,ja)
      i = 1	  
      DO WHILE(e_temp(i).lt.EMAX)
      i = i + 1
c      PRINT*,"E",i,e_temp(i)	  
      ENDDO
      nlvl = i - 1	  
c      PRINT*,'nlvl',nlvl      	  
      decr = decr + 1 
      DEALLOCATE(j_ch_tmp,ka_ch_tmp,
     & kc_ch_tmp,e_temp,ja)	  
      ENDDO
       number_of_channels = nlvl
c      PRINT*,"number_of_channels = ",number_of_channels 	   
      ALLOCATE(j_ch(number_of_channels),ka_ch(number_of_channels),
     & kc_ch(number_of_channels),E_ch(number_of_channels)) 
       ALLOCATE(j_ch_tmp(nchann_est),ka_ch_tmp(nchann_est),
     & kc_ch_tmp(nchann_est),e_temp(nchann_est),ja(nchann_est))
	 
      nlvl = 0
      DO j_t=0,j_max_allowed
      DO dk_t=-j_t,j_t
      nlvl = nlvl + 1
      IF(nlvl.gt.nchann_est) EXIT	  
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst)
     &	  STOP"ERROR:ASSYM. TOP INITALIZATION FAILED FOR EMAX_DEFINED"
      j_ch_tmp(nlvl) = j_t
      ka_ch_tmp(nlvl) = ka_t
      kc_ch_tmp(nlvl) = kc_t
      CALL EIGEN_VAL(A,B,C,j_t,ka_t,kc_t,e_temp(nlvl))	  
      ENDDO
      IF(nlvl.gt.nchann_est) EXIT		  
      ENDDO
c      PRINT*,"number_of_channels_again = ",number_of_channels     
      CALL	hpsort(nchann_est,e_temp,ja)	  
      DO i=1,number_of_channels
      j_ch(i) = j_ch_tmp(ja(i))
      ka_ch(i) = ka_ch_tmp(ja(i))
      kc_ch(i) = kc_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      DEALLOCATE(j_ch_tmp,ka_ch_tmp,
     & kc_ch_tmp,e_temp,ja)
      energy_defined = .TRUE.	 
       IF(MYID.eq.0) 	  
     & WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"	  
       IF(MYID.eq.0) 	  
     & WRITE(1,'(a5,9x,a2,9x,a2,9x,a2,9x,a2)')"#N = "
     & ,"J","KA","KC", "E"
      jmax_included = 0	 
      DO i=1,number_of_channels
      IF(jmax_included.lt.j_ch(i)) jmax_included = j_ch(i)	  
      IF(.not.energy_defined) THEN
      CALL EIGEN_VAL(A,B,C,j_ch(i),ka_ch(i),kc_ch(i),E_ch(i))	  
      ENDIF	
       IF(MYID.eq.0) 	  
     & WRITE(1,'(i5,i11,i11,i11,f11.3)')i,
     & j_ch(i),ka_ch(i),kc_ch(i),E_ch(i)
      ENDDO
      energy_defined = .TRUE.
      channels_defined = .TRUE.	  
      ENDIF		  
      ENDIF
       IF(user_defined) THEN	  
       IF(MYID.eq.0) 	  
     & WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"	  
       IF(MYID.eq.0) 	  
     & WRITE(1,'(a5,9x,a2,9x,a2,9x,a2,9x,a2)')"#N = "
     & ,"J","KA","KC", "E"
      jmax_included = 0
      IF(para_ortho_defined) CALL SYMMETRY_SORT_PARA_ORTHO	  
      DO i=1,number_of_channels
      IF(jmax_included.lt.j_ch(i)) jmax_included = j_ch(i)	  
      IF(MYID.eq.0) 	  
     & WRITE(1,'(i5,i11,i11,i11,f11.3)')i,
     & j_ch(i),ka_ch(i),kc_ch(i),E_ch(i)
      ENDDO
      ENDIF	  
      IF(ini_chann_defined .and. .not.all_level_included) THEN
      chann_ini = 0
      DO i=1,number_of_channels
      IF(j_ini.eq.j_ch(i) .and. ka_ini.eq.ka_ch(i) .and. 
     & kc_ini.eq.kc_ch(i)) chann_ini = i
      ENDDO	
      IF(chann_ini.eq.0) STOP "ERROR: INITIAL CHANNEL WRONG"
       IF(MYID.eq.0) 	  
     & WRITE(1,'(a17,a4,i4,a7,i4,a7,i4)') "INITIAL STATE IS:",
     & "J = ",j_ini,"; KA = ", ka_ini,"; KC = ",kc_ini  
      ENDIF
      IF(.not.ini_chann_defined .and. channels_defined .and.
     & energy_defined) THEN
	  i_ground = 1
      DO i=1,number_of_channels
      IF(E_ch(i_ground).gt.E_ch(i)) i_ground = i   
      ENDDO	 
      chann_ini = i_ground
      ini_chann_defined = .TRUE.	  
      ENDIF	  
      IF(myid.eq.0) THEN	  
      WRITE(1,*)
      IF(angs_unit)WRITE(1,'(a35)')"THE POTENTIAL R-UNITS ARE ANGSTROMS"
      IF(au_unit_r)WRITE(1,'(a31)') "THE POTENTIAL R-UNITS ARE BOHRS"
      IF(cm_unit)WRITE(1,'(a31)')"THE POTENTIAL ENERGY IS IN CM-1"
      IF(au_unit_e)WRITE(1,'(a34)')"THE POTENTIAL ENERGY IS IN HARTREE"
      IF(klv_unit)WRITE(1,'(a33)')"THE POTENTIAL ENERGY IS IN KELVIN"
      IF(kal_unit)WRITE(1,'(a35)')"THE POTENTIAL ENERGY IS IN KCAL/MOL"	  
      IF(.not.angs_unit .and..not.au_unit_r )   
     & STOP "ERROR:R_UNITS NOT DEFINED"
      IF(.not.cm_unit .and..not.au_unit_e.and..not.klv_unit
     & .and. .not.kal_unit)   
     & STOP "ERROR:ENERGY NOT DEFINED"
      WRITE(1,*)
      ENDIF	  
      IF(expansion_defined) THEN
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a42)') "USER SUPPLIED POTENTIAL TERMS WILL BE USED"
      WRITE(1,*)
      WRITE(1,'(a45,i4)')
     & "NUMBER OF POTENTIAL TERMS INCLUDED N_TERMS = ",nterms
      WRITE(1,'(a38)')
     & "THE POTENTIAL TERMS ARE, (LAMBDA,NJU):"
      DO i=1,nterms
      WRITE(1,'(a1,i3,a1,i3,a3)',ADVANCE ="NO")
     & '(',L_TERM(i),',',M_TERM(i),'); '	  
      ENDDO	
      WRITE(1,*)
      WRITE(1,*)	  
      WRITE(1,'(a33,i3)') "NUMBER OF R-GRID POINTS IS N_R = ",
     & n_r_coll
      ENDIF	 
      ELSE
      IF(.not.matrix_reading_defined .and. myid.eq.0)	  
     & WRITE(1,'(a97)')"COULPING MATRIX WILL BE COMPUTED USING NUMERICAL
     &  INTEGRATION OF POTENTIAL OVER THE WAVEFUNCTIONS"	  
      ENDIF	  
      CASE(5)
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a23,1x,a32)')"THE COLLISIONAL SYSTEM:",
     & "RIGID DIATOMIC TOP + DIATOMC TOP"
      IF(identical_particles_defined) THEN
      WRITE(1,*)	  
   	  WRITE(1,'(a57)')
     & "THE COLLISIONAL PARTNERS ARE DEFINED AS INDINSTIGUISHABLE"
      WRITE(1,*)	 
      ENDIF
      ENDIF	  
      IF(user_defined) THEN
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a7,2x,a60)') "WARNING:",
     & "THE WAVEFUNCTIONS AND ENERGY LEVELS MUST BE SUPPLIED BY USER"
      IF(rot_const_defined)WRITE(1,'(a27)')"ROT. CONST WILL BE INGNORED"
      ENDIF	  
      INQUIRE( FILE=USR_INP_LEVEL, EXIST=exst )
      IF(.not.exst)	THEN
      IF(MYID.eq.0)PRINT*,USR_INP_LEVEL, " NOT FOUND"
      STOP	  
      ENDIF	  
      OPEN(325,FILE=USR_INP_LEVEL,
     & STATUS="OLD",ACTION="READ")
      IF(number_of_channels.le.0) THEN
      IF(MYID.eq.0) PRINT*,"NUMBER OF CHANNELS NOT SPECIFIED"
      STOP	  
      ENDIF
      IF(n_r_vib1*n_r_vib2.le.0) THEN
      IF(MYID.eq.0) PRINT*,"INTEGRATION GRID NOT DEFINED"
      STOP	  
      ENDIF	  
      IF(.not.channels_defined) THEN	  
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))
      ENDIF
      IF(.not.energy_defined) THEN	  
      ALLOCATE(E_ch(number_of_channels))
      ENDIF
      READ(325,*) !!! HEADER	  
      DO i=1,number_of_channels
      READ(325,'(i4,1x,i3,1x,i3,1x,e17.10)',
     & IOSTAT = stat_of_file)
     & decr,  j1_ch(i), 
     & j2_ch(i), E_ch(i)
      IF(stat_of_file.gt.0 .or. i.ne. decr) THEN
      IF(jmax_included.lt.j1_ch(i))	jmax_included = j1_ch(i)
      IF(jmax_included.lt.j2_ch(i))	jmax_included = j2_ch(i)  
	  
      PRINT*,"WRONG CHANNELS IN USER INPUT FILE"	  
      STOP
      ENDIF 	  
      ENDDO

      CLOSE(325)
      channels_defined = .TRUE.
      energy_defined = .TRUE.	
	  	 
      ELSE
      IF(BE.ne.0 .or. DE.ne.0) WRITE(1,'(a49)')
     & "WARNING: BE AND DE ARE FOR SYS_TYPE #1, NOT FOR #5"
      IF(BE1.eq.0) STOP "ERROR:SUPPLY POSITVE BE1"
      IF(DE1.GT.BE1)
     &	  STOP"ERROR: CHECK INPUT-DE1 SHOULD MUCH LESS THAN BE1"
      IF(identical_particles_defined) THEN
      BE2 = BE1
      DE2 = DE1	  
      ENDIF		 
      IF(BE2.eq.0) STOP "ERROR:SUPPLY POSITVE BE2"
      IF(DE2.GT.BE2)
     &	 STOP"ERROR: CHECK INPUT-DE2 SHOULD MUCH LESS THAN BE2"
      IF(MYID.eq.0) THEN		 
      WRITE(1,'(a49,f10.4)')
     & "ROTATIONAL CONSTANT FOR MOLECULE #1, CM-1, BE1 = ",Be1 	  
      WRITE(1,'(a48,f12.7)')
     & "DISTORTION FOR MOLECULE #1, CM-1, DE1 = ",DE1
      WRITE(1,'(a49,f10.4)')
     & "ROTATIONAL CONSTANT FOR MOLECULE #2, CM-1, BE2 = ",Be2 	  
      WRITE(1,'(a48,f12.7)')
     & "DISTORTION FOR MOLECULE #2, CM-1, DE2 = ",DE2
      ENDIF	 
      IF(jmax1.ge.0 .and. jmin1.ge.0 .and. 
     & jmax2.ge.0 .and. jmin2.ge.0) THEN
      IF(jmax1.lt.jmin1) STOP "ERROR: JMAX1 MUST BE > OR = THAN JMIN1"
      IF(jmax2.lt.jmin2) STOP "ERROR: JMAX1 MUST BE > OR = THAN JMIN1"	  
      number_of_channels_defined = .TRUE.
      energy_defined = .TRUE.
      channels_defined = .TRUE.
      number_of_channels = (jmax1 - jmin1 + 1)*(jmax2 - jmin2 + 1)
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))	  
      ALLOCATE(E_ch(number_of_channels))
      i = 0	  
      DO j1_t = jmin1,jmax1
      DO j2_t = jmin2,jmax2
      i = i + 1
      j1_ch(i) = j1_t
      j2_ch(i) = j2_t	  
      E_ch(i) = Be1*j1_t*(j1_t+1d0) -
     &	  De1*(j1_t*(j1_t+1d0))**2 + Be2*j2_t*(j2_t+1d0) -
     &	  De2*(j2_t*(j2_t+1d0))**2 	  
      ENDDO	  
      ENDDO
	  
      ENDIF	 
      IF(number_of_channels_defined
     &  .and. .not.channels_defined) THEN
      nchann_est = (number_of_channels+1)**2	   
      ALLOCATE(j1_ch_tmp(nchann_est),j2_ch_tmp(nchann_est),
     & e_temp(nchann_est),ja(nchann_est))
      i = 0
      DO j1_t =0,number_of_channels
      DO j2_t =0,number_of_channels
      i = i + 1	  
      j1_ch_tmp(i) = j1_t
      j2_ch_tmp(i) = j2_t 	  
      e_temp(i) = Be1*j1_t*(j1_t+1d0) -
     &	  De1*(j1_t*(j1_t+1d0))**2 + Be2*j2_t*(j2_t+1d0) -
     &	  De2*(j2_t*(j2_t+1d0))**2 	  	  
      ENDDO
      ENDDO
      CALL hpsort(nchann_est,e_temp,ja)
      IF(ALLOCATED(j1_ch)) DEALLOCATE (j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE (j2_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE (E_ch)		  
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))	  
      ALLOCATE(E_ch(number_of_channels))	  
      ALLOCATE(j1_ch(number_of_channels),j2_ch(number_of_channels),
     & E_ch(number_of_channels))
      DO i=1,number_of_channels
      E_ch(i) = e_temp(i)
      j1_ch(i) = j1_ch_tmp(ja(i))
      j2_ch(i) = j2_ch_tmp(ja(i))	  
      ENDDO	  
      DEALLOCATE(j1_ch_tmp,j2_ch_tmp,e_temp,ja)	 
      energy_defined = .TRUE.
      channels_defined = .TRUE.	  
	  
      ELSE
      IF(.not.energy_defined .and.channels_defined) THEN
      IF(number_of_channels.le.0)
     &  STOP"ERROR:N_CHANNELS MUST >0"	  
      ALLOCATE(E_ch(number_of_channels))
      DO i =1,number_of_channels
      j1_t = j1_ch(i)
      j2_t = j2_ch(i)  	  
      E_ch(i) = Be1*j1_t*(j1_t+1d0) -
     &	  De1*(j1_t*(j1_t+1d0))**2 + Be2*j2_t*(j2_t+1d0) -
     &	  De2*(j2_t*(j2_t+1d0))**2 	  	  
      ENDDO	  
      ENDIF
      energy_defined = .TRUE.	  
      ENDIF	 
      IF(mlc_mlc_chn_num_defined .and. .not.channels_defined) THEN
      IF(number_of_channels_defined) WRITE(1,'(a58)')
     & "WARNING: TOTAL N_CHANNELS AND THEIR LEVELS WILL BE IGNORED"
      IF(mlc_mlc_emax_defined)
     & STOP"ERROR:DEFINE EMAX1,2 OR EITHER NCHNL1,2"
      number_of_channels_defined = .TRUE.
      energy_defined = .TRUE.
      channels_defined = .TRUE.
      number_of_channels = nchann_1*nchann_2
      IF(ALLOCATED(j1_ch)) DEALLOCATE (j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE (j2_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE (E_ch)	  
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))	  
      ALLOCATE(E_ch(number_of_channels))
      i = 0	  
      DO j1_t = 0,nchann_1
      DO j2_t = 0,nchann_2
      i = i + 1
      j1_ch(i) = j1_t
      j2_ch(i) = j2_t	  
      E_ch(i) = Be1*j1_t*(j1_t+1d0) -
     &	  De1*(j1_t*(j1_t+1d0))**2 + Be2*j2_t*(j2_t+1d0) -
     &	  De2*(j2_t*(j2_t+1d0))**2 	  
      ENDDO	  
      ENDDO	 
      ENDIF
      IF(mlc_mlc_emax_defined.and. .not.channels_defined) THEN
c      PRINT *,	"mlc_mlc_emax_defined"  
      IF(emax_defined) WRITE(1,'(a38)')
     & "WARNING: VALUE OF EMAX WILL BE IGNORED"
      nchann_est = int(sqrt(EMAX1/Be1))
      i = nchann_est - 1	  
      IF(Be1*i*(i+1d0) - De1*(i*(i+1d0))**2 .lt. EMAX1) THEN
      DO WHILE(Be1*i*(i+1d0) - De1*(i*(i+1d0))**2 .lt. EMAX1)
      i = i +1
      ENDDO
      nchann_est = i - 1	  
      ENDIF
      IF(Be1*i*(i+1d0) - De1*(i*(i+1d0))**2 .GT. EMAX1) THEN
      DO WHILE(Be1*i*(i+1d0) - De1*(i*(i+1d0))**2 .GT. EMAX1)
      i = i -1
      ENDDO
      nchann_est = i  
      ENDIF	 
      nchann_1 = 	nchann_est+1

      nchann_est = int(sqrt(EMAX2/Be2))
      i = nchann_est - 1	  
      IF(Be2*i*(i+1d0) - De2*(i*(i+1d0))**2 .lt. EMAX2) THEN
      DO WHILE(Be2*i*(i+1d0) - De2*(i*(i+1d0))**2 .lt. EMAX2)
      i = i +1
      ENDDO
      nchann_est = i - 1	  
      ENDIF
      IF(Be2*i*(i+1d0) - De2*(i*(i+1d0))**2 .GT. EMAX2) THEN
      DO WHILE(Be2*i*(i+1d0) - De2*(i*(i+1d0))**2 .GT. EMAX2)
      i = i -1
      ENDDO
      nchann_est = i  
      ENDIF	 
      nchann_2 = 	nchann_est+1
      number_of_channels_defined = .TRUE.
      energy_defined = .TRUE.
      channels_defined = .TRUE.
      number_of_channels = nchann_1*nchann_2
      IF(ALLOCATED(j1_ch)) DEALLOCATE (j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE (j2_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE (E_ch)		  
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))	  
      ALLOCATE(E_ch(number_of_channels))
      i = 0	  
      DO j1_t = 0,nchann_1-1
      DO j2_t = 0,nchann_2-1
      i = i + 1
      j1_ch(i) = j1_t
      j2_ch(i) = j2_t	  
      E_ch(i) = Be1*j1_t*(j1_t+1d0) -
     &	  De1*(j1_t*(j1_t+1d0))**2 + Be2*j2_t*(j2_t+1d0) -
     &	  De2*(j2_t*(j2_t+1d0))**2 	  
      ENDDO	  
      ENDDO	 
      WRITE(1,*)	  
      ENDIF 
      IF (emax_defined.and. .not.channels_defined) THEN
      nchann_est = MAX(int(sqrt(EMAX/Be1))+1,int(sqrt(EMAX/Be2))+1)**2
      ALLOCATE(j1_ch_tmp(nchann_est),j2_ch_tmp(nchann_est),
     & e_temp(nchann_est),ja(nchann_est))
      i = 0	 
      DO j1_t =0,nchann_est
      IF(i.gt.nchann_est) EXIT		  
      DO j2_t =0,nchann_est
      IF(i.gt.nchann_est) EXIT		  
      i = i +1	  
      j1_ch_tmp(i) = j1_t
      j2_ch_tmp(i) = j2_t 	  
      e_temp(i) = Be1*j1_t*(j1_t+1d0) -
     &	  De1*(j1_t*(j1_t+1d0))**2 + Be2*j2_t*(j2_t+1d0) -
     &	  De2*(j2_t*(j2_t+1d0))**2 	  	  
      ENDDO
      ENDDO	
      CALL	hpsort(nchann_est,e_temp,ja)
      nlvl = 0	  
      DO WHILE(.not.channels_defined)
      nlvl = nlvl + 1
      IF(e_temp(nlvl).gt.EMAX) THEN
      number_of_channels = 
     & nlvl - 1
      channels_defined = .TRUE.	 
      ENDIF	 
      ENDDO
      IF(number_of_channels.le.0)STOP"ERROR:ZERO NUMBER OF CHANNELS"
      IF(ALLOCATED(j1_ch)) DEALLOCATE (j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE (j2_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE (E_ch)	  
      ALLOCATE(j_ch(number_of_channels))
      ALLOCATE(v_ch(number_of_channels))	  
      ALLOCATE(E_ch(number_of_channels)) 	  
      ALLOCATE(j_ch(number_of_channels),
     & v_ch(number_of_channels),E_ch(number_of_channels))
      DO nlvl = 1,number_of_channels
      n_prev = ja(nlvl)	  
      j_ch(nlvl) = j_ch_tmp(n_prev)
      v_ch(nlvl) = v_ch_temp(n_prev)
      E_ch(nlvl) = e_temp(nlvl)
      ENDDO
      DEALLOCATE(j_ch_tmp,v_ch_temp,e_temp,ja)
      energy_defined = .TRUE.
      channels_defined = .TRUE.	
	  
      ENDIF	 
      ENDIF
      jmax_included = 0
      IF(para_ortho_defined) CALL SYMMETRY_SORT_PARA_ORTHO      	  
      IF(myid.eq.0) WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"	  
      IF(myid.eq.0)     WRITE(1,'(a5,9x,a2,9x,a2,9x,a2)')"#N = "
     & ,"J1","J2","E"
      DO i=1,	number_of_channels
      IF(jmax_included.lt.j1_ch(i)) jmax_included = j1_ch(i)
      IF(jmax_included.lt.j2_ch(i)) jmax_included = j2_ch(i)
      IF(identical_particles_defined) THEN
      IF(j1_ch(i).lt.j2_ch(i)) THEN
      IF(myid.eq.0)WRITE(1,'(a42)')
     & "ERROR: FOR IDENT. PARTICLES MUST BE J1>=J2"
      STOP	  
      ENDIF	  
      ENDIF	  
      IF(myid.eq.0) WRITE(1,'(i5,i11,i11,f11.3)')i,
     & j1_ch(i),j2_ch(i),E_ch(i)
      ENDDO	  
      IF(myid.eq.0) WRITE(1,*)	  
      IF(ini_chann_defined .and. .not.all_level_included) THEN
      chann_ini = 0	  
      DO i=1,number_of_channels
      IF(j1_ini.eq.j1_ch(i) .and.j2_ini.eq.j2_ch(i)) chann_ini = i
      ENDDO	
      IF(chann_ini.eq.0) STOP "ERROR: INITIAL CHANNEL WRONG"
      IF(myid.eq.0) WRITE(1,'(a17,a8,i4,1x,i4)') "INITIAL STATE IS:",
     & "J1,J2 = ",j1_ini,j2_ini 	  
      ENDIF
      IF(.not.ini_chann_defined .and. channels_defined .and.
     & energy_defined) THEN
	  i_ground = 1
      DO i=1,number_of_channels
      IF(E_ch(i_ground).gt.E_ch(i)) i_ground = i   
      ENDDO	 
      chann_ini = i_ground
      ini_chann_defined = .TRUE.	  
      ENDIF	  
      IF(myid.eq.0) THEN	  
      WRITE(1,*)
      IF(angs_unit)WRITE(1,'(a35)')"THE POTENTIAL R-UNITS ARE ANGSTROMS"
      IF(au_unit_r)WRITE(1,'(a31)') "THE POTENTIAL R-UNITS ARE BOHRS"
      IF(cm_unit)WRITE(1,'(a31)')"THE POTENTIAL ENERGY IS IN CM-1"
      IF(au_unit_e)WRITE(1,'(a34)')"THE POTENTIAL ENERGY IS IN HARTREE"
      IF(klv_unit)WRITE(1,'(a33)')"THE POTENTIAL ENERGY IS IN KELVIN"
      IF(kal_unit)WRITE(1,'(a35)')"THE POTENTIAL ENERGY IS IN KCAL/MOL"
      ENDIF	  
      IF(.not.angs_unit .and..not.au_unit_r )   
     & STOP "ERROR:R_UNITS NOT DEFINED"
      IF(.not.cm_unit .and..not.au_unit_e.and..not.klv_unit
     & .and. .not.kal_unit)   
     & STOP "ERROR:ENERGY NOT DEFINED"
      IF(myid.eq.0)WRITE(1,*) 	 
      IF(expansion_defined .and. myid.eq.0) THEN
      WRITE(1,'(a42)') "USER SUPPLIED POTENTIAL TERMS WILL BE USED"
      WRITE(1,*)
      WRITE(1,'(a45,i4)')
     & "NUMBER OF POTENTIAL TERMS INCLUDED N_TERMS = ",nterms
      WRITE(1,'(a38)')
     & "THE POTENTIAL TERMS ARE, (L1,L2,L):"
      DO i=1,nterms
      IF(max(L1_TERM(i),L1_TERM(i),L_TERM(i))*2 .gt.
     & L2_TERM(i)+L1_TERM(i)+L_TERM(i)) THEN 
      WRITE(1,'(a47)')
     &"WARNING: L1, L2 and L must follow triangular rule"
      WRITE(1,'(a3,1x,i3,1x,a3,1x,i3,1x,a3,1x,i3)') 
     & "L1=",L1_TERM(i),"L2=",L1_TERM(i),"L=",L_TERM(i)      	  
      ENDIF	  
      WRITE(1,'(a1,i3,1x,i3,1x,i3,a3)',ADVANCE ="NO")
     & '(',L1_TERM(i),L2_TERM(i),L_TERM(i),'); '	  
      ENDDO	
      WRITE(1,*)
      WRITE(1,*)	  
      WRITE(1,'(a33,i3)') "NUMBER OF R-GRID POINTS IS N_R = ",
     & n_r_coll	  
      ELSE
      IF(.not.matrix_reading_defined .and. myid.eq.0)	  
     & WRITE(1,'(a97)')"COULPING MATRIX WILL BE COMPUTED USING NUMERICAL
     &  INTEGRATION OF POTENTIAL OVER THE WAVEFUNCTIONS"	  
      ENDIF	  
	  
      CASE(6)
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a23,1x,a36)')"THE COLLISIONAL SYSTEM:",
     & "VIBRATING DIATOMIC TOP + DIATOMC TOP"
      IF(identical_particles_defined) THEN
      WRITE(1,*)	  
   	  WRITE(1,'(a57)')
     & "THE COLLISIONAL PARTNERS ARE DEFINED AS INDINSTIGUISHABLE"
      WRITE(1,*)	 
      ENDIF
      ENDIF	  
      IF(user_defined) THEN
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a7,2x,a60)') "WARNING:",
     & "THE WAVEFUNCTIONS AND ENERGY LEVELS MUST BE SUPPLIED BY USER"
      IF(rot_const_defined)WRITE(1,'(a27)')"ROT. CONST WILL BE INGNORED"
      IF(vib_const_defined)WRITE(1,'(a27)')"VIB. CONST WILL BE INGNORED"
      ENDIF	  
      INQUIRE( FILE=USR_INP_LEVEL, EXIST=exst )
      IF(.not.exst)	THEN
      IF(MYID.eq.0)PRINT*,USR_INP_LEVEL, " NOT FOUND"
      STOP	  
      ENDIF	  
      OPEN(325,FILE=USR_INP_LEVEL,
     & STATUS="OLD",ACTION="READ")
      IF(number_of_channels.le.0) THEN
      IF(MYID.eq.0) PRINT*,"NUMBER OF CHANNELS NOT SPECIFIED"
      STOP	  
      ENDIF
      IF(n_r_vib1*n_r_vib2.le.0) THEN
      IF(MYID.eq.0) PRINT*,"INTEGRATION GRID NOT DEFINED"
      STOP	  
      ENDIF	  
      IF(.not.channels_defined) THEN	  
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(v1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))
      ALLOCATE(v2_ch(number_of_channels))	  
      ENDIF
      IF(.not.energy_defined) THEN	  
      ALLOCATE(E_ch(number_of_channels))
      ENDIF
      READ(325,*) !!! HEADER	  
      DO i=1,number_of_channels
      READ(325,'(i4,1x,i3,1x,i3,1x,i3,1x,i3,1x,e17.10)',
     & IOSTAT = stat_of_file)
     & decr,  j1_ch(i), v1_ch(i),
     & j2_ch(i), v2_ch(i),E_ch(i)
      IF(jmax_included.lt.j1_ch(i))	jmax_included = j1_ch(i)
      IF(jmax_included.lt.j2_ch(i))	jmax_included = j2_ch(i) 
      IF(vmax_included.lt.v1_ch(i))	vmax_included = v1_ch(i)
      IF(vmax_included.lt.v2_ch(i))	vmax_included = v2_ch(i) 	  
      IF(stat_of_file.gt.0 .or. i.ne. decr) THEN
      PRINT*,"WRONG CHANNELS IN USER INPUT FILE"	  
      STOP
      ENDIF 	  
      ENDDO
      ALLOCATE(
     & vibrational_wavefunctions_1(n_r_vib1,
     & number_of_channels))
      ALLOCATE(
     & vibrational_wavefunctions_2(n_r_vib2,
     & number_of_channels))
      ALLOCATE(
     & jacobian_vib1(n_r_vib1))
      ALLOCATE(
     & jacobian_vib2(n_r_vib2))
      ALLOCATE(r_grid_vib1(n_r_vib1))
      ALLOCATE(r_grid_vib2(n_r_vib2))	  
c      ! WRITTING A GRID	 
      READ(325,*) ! A HEADER
      DO i=1,n_r_vib1! LOOP OVER VIBRATIONAL GRID
      READ(325,*)r_grid_vib1(i),jacobian_vib1(i)	  
!      READ(325,*)r_vb_dt_integrator1(i),wgg(i)
	  ! r_vb_dt_inegrator(i)-
! r-value  for i-point of the grid, wg(i) - weight(jacobian)
      ENDDO
      READ(325,*)! "R2(i),J2(i)" ! A HEADER
      DO i=1,n_r_vib2! LOOP OVER VIBRATIONAL GRID
      READ(325,*)r_grid_vib2(i),jacobian_vib2(i)	  
!      READ(325,*)r_vb_dt_integrator2(i),wbb(i)
	  ! r_vb_dt_inegrator(i)-
! r-value  for i-point of the grid, wg(i) - weight(jacobian)
      ENDDO	  
      READ(325,*)! "WAVEFUNCTIONS ARE LISTED BELOW"	  
      DO k=1,number_of_channels !! k - channel number
c      !! WRITING THE WAVEFUNCTIONS	  
      READ(325,*)!	  "MOLECULE#1,CHANNEL=",k !!! A HEADER
      DO i=1,n_r_vib1
      READ(325,*) vibrational_wavefunctions_1(i,k)	  
!      READ(325,*)wf_vib_part_diatom1(i,k) ! VALUE OF WAVEFUNCTION FOR k-channel at i-point of the grid
      ENDDO
      READ(325,*) !! LEAVE IT BLANK	  
      ENDDO

      DO k=1,number_of_channels !! k - channel number
c      !! WRITING THE WAVEFUNCTIONS	  
      READ(325,*)! "MOLECULE#2,CHANNEL=",k !!! A HEADER
      DO i=1,n_r_vib2
      READ(325,*) vibrational_wavefunctions_2(i,k)	  
!      READ(456,*)wf_vib_part_diatom2(i,k) ! VALUE OF WAVEFUNCTION FOR k-channel at i-point of the grid
      ENDDO
      READ(325,*) !! LEAVE IT BLANK	  
      ENDDO	  
 
      CLOSE(325)
      channels_defined = .TRUE.
      energy_defined = .TRUE.	
	  
      ELSE
      IF(morse_pot_defined) THEN
!      IF(morse12_re(1).le.0 .or. morse12_re(2).le.0 )								!Dulat: moved these lines further down because Morse R needs atomic masses
!     & STOP "ERROR: YOU MUST SUPPLY R_EQ FOR MORSE POTENTIAL"
      IF(atomic_masses(1)*atomic_masses(2).le.0) STOP
     & "ERROR: SUPPLY MASS OF ATOMS"
      atomic_red_mass1 = atomic_masses(1)*atomic_masses(2)/
     & (atomic_masses(1)+atomic_masses(2))
      atomic_red_mass2 = atomic_masses(3)*atomic_masses(4)/
     & (atomic_masses(3)+atomic_masses(4))	 
      IF(MYID.eq.0) THEN	 
      WRITE(1,'(a34,f12.7)') "#1 MOLECULAR REDUCED MASS,A.M.U = ",
     & atomic_red_mass1
      WRITE(1,'(a34,f12.7)') "#2 MOLECULAR REDUCED MASS,A.M.U = ",
     & atomic_red_mass2	 
      ENDIF	 
      mass_vib_diatom1 = atomic_masses(1)+atomic_masses(2)
      mass_vib_diatom2 = atomic_masses(3)+atomic_masses(4)		  
	  
! Dulat Bostan 08/08/2022
	  if(morse12_re(1).le.0) then
	  if(MYID.eq.0) print*, "Molecule #1"
	  if(MYID.eq.0) print*, 
     & "Equilibrium bond distance is not specified in the input"
	  morse12_re(1) = dsqrt(1.d0/2.d0/(atomic_red_mass1*amutoau)
     & /(Be1/eVtown/autoeV))
	  if(MYID.eq.0) print*, "Equilibrium bond distance:", 
     & morse12_re(1), "a.u."
	  else
	  if(MYID.eq.0) print*, "Equilibrium bond distance:", 
     & morse12_re(1), "a.u."
	  endif
	  
	  if(morse12_re(2).le.0) then
	  if(MYID.eq.0) print*, "Molecule #2"
	  if(MYID.eq.0) print*, 
     & "Equilibrium bond distance is not specified in the input"
	  morse12_re(2) = dsqrt(1.d0/2.d0/(atomic_red_mass2*amutoau)
     & /(Be2/eVtown/autoeV))
	  if(MYID.eq.0) print*, "Equilibrium bond distance:", 
     & morse12_re(2), "a.u."
	  else
	  if(MYID.eq.0) print*, "Equilibrium bond distance:", 
     & morse12_re(2), "a.u."
	  endif
! Dulat Bostan - end	  
	  
      IF(morse12_a(1).gt.0 .and. morse12_de(1).gt.0 .and.
     & .not. vib_const_defined) THEN
!      WE1 = 1d0/morse12_a(1)													!Dulat: Disabled this line and corrected equation in the next line
      WE1 = morse12_a(1)
     & /dsqrt(atomic_red_mass1*amutoau/2d0/morse12_de(1)/eVtown/autoeV)
      XE1 = WE1**2/4d0/morse12_de(1)
      ENDIF
      IF(morse12_a(2).gt.0 .and. morse12_de(2).gt.0 .and.
     & .not. vib_const_defined) THEN
!      WE2 = 1d0/morse12_a(2)													!Dulat: Disabled this line and corrected equation in the next line
      WE2 = morse12_a(2)
     & /dsqrt(atomic_red_mass1*amutoau/2d0/morse12_de(2)/eVtown/autoeV)
      XE2 = WE2**2/4d0/morse12_de(2)
      ENDIF	  
      IF(r_vib_diatom_max12(1).le.0d0 .or. r_vib_diatom_min12(1).lt.0d0)
     & STOP
     & "ERROR: SUPPLY MORSE POTENTIAL RANGE."
      IF(r_vib_diatom_max12(2).le.0d0 .or. r_vib_diatom_min12(2).lt.0d0)
     & STOP
     & "ERROR: SUPPLY MORSE POTENTIAL RANGE."	 
      IF(MYID.eq.0) THEN	 
      WRITE(1,'(a33,f10.4,1x,f10.4)')
     & "RMIN FOR MORSE POTENTIAL, a.u. = ",r_vib_diatom_min12(1),
     & r_vib_diatom_min12(2)	 
      WRITE(1,'(a33,f10.4,1x,f10.4)')
     & "RMAX FOR MORSE POTENTIAL, a.u. = ",r_vib_diatom_max12(1),
     & r_vib_diatom_max12(2)
      ENDIF	 
      ENDIF
      IF(identical_particles_defined) THEN
      BE2 = BE1
      DE2 = DE1
      WE2 = WE1
      XE2 = XE1	  
      ENDIF	  
      IF(BE1.le.0) STOP "ERROR:SUPPLY POSITVE BE1"
      IF(DE1.GT.BE1)
     & STOP"ERROR: CHECK INPUT-DE1 SHOULD MUCH LESS THAN BE1"	
      IF(WE1.le.0 ) STOP "ERROR:SUPPLY POSITVE WE1"
      IF(XE1.le.0 .or.XE1.ge.WE1 ) STOP "ERROR:SUPPLY POSITIVE XE1<WE1"
!      IF(XE1.le.0 .or.XE1.ge.1 ) STOP "ERROR:SUPPLY POSITVE XE1<1"				!Dulat: Disabled this line
	  
      IF(BE2.le.0) STOP "ERROR:SUPPLY POSITVE BE2"
      IF(DE2.GT.BE2)
     & STOP"ERROR: CHECK INPUT-DE2 SHOULD MUCH LESS THAN BE2"	
      IF(WE2.le.0 ) STOP "ERROR:SUPPLY POSITVE WE2"
      IF(XE2.le.0 .or.XE2.ge.WE2 ) STOP "ERROR:SUPPLY POSITIVE XE2<WE2"
!      IF(XE2.le.0 .or.XE1.ge.1 ) STOP "ERROR:SUPPLY POSITVE XE2<1"				!Dulat: Disabled this line
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a32,f10.4)') "ROTATIONAL CONSTANT, CM-1, BE1 = ",Be1	  
      WRITE(1,'(a23,f12.7)') "DISTORTION, CM-1, DE1 = ",DE1	 
      WRITE(1,'(a33,f10.4)') "VIBRATIONAL CONSTANT, CM-1, WE1 = ",We1	  
      WRITE(1,'(a24,f12.7)') "ANHARMONICITY, <1, XE1 = ",XE1
      WRITE(1,'(a32,f10.4)') "ROTATIONAL CONSTANT, CM-1, BE2 = ",Be2	  
      WRITE(1,'(a23,f12.7)') "DISTORTION, CM-1, DE2 = ",DE2	 
      WRITE(1,'(a33,f10.4)') "VIBRATIONAL CONSTANT, CM-1, WE2 = ",We2	  
      WRITE(1,'(a24,f12.7)') "ANHARMONICITY, <1, XE2 = ",XE2
      ENDIF
      IF(mlc_mlc_emax_defined .or. emax_defined) THEN
      IF(emax_defined) THEN
      EMAX1=EMAX
      EMAX2=EMAX
      ENDIF
      IF(EMAX1*EMAX2.lt.0) STOP "ERROR IN EMAX DEFINITION"	  
      nchann_est = (int(EMAX1/WE1)+2)*(int(DSQRT(EMAX1/BE1))+2)*
     & (int(EMAX2/WE2)+2)*(int(DSQRT(EMAX2/BE2))+2)
!      IF(MYID.eq.0) PRINT*, "nchann_est",nchann_est	 
      ALLOCATE(j1_ch_tmp(nchann_est),v1_ch_temp(nchann_est),
     & e1_temp(nchann_est),ja(nchann_est),ja1(nchann_est))
      ALLOCATE(j2_ch_tmp(nchann_est),v2_ch_temp(nchann_est),
     & e2_temp(nchann_est),e_temp(nchann_est),ja2(nchann_est))	 
      n_prev = 0	 
      DO v1_t = 0,int(EMAX1/WE1)+1
      DO j1_t = 0,int(DSQRT(EMAX1/BE1))+1
      DO v2_t = 0,int(EMAX2/WE2)+1
      DO j2_t = 0,int(DSQRT(EMAX2/BE2))+1  
      n_prev = n_prev + 1	  
      j1_ch_tmp(n_prev) = j1_t
      v1_ch_temp(n_prev) = v1_t
      j2_ch_tmp(n_prev) = j2_t
      v2_ch_temp(n_prev) = v2_t
      CALL energy_diatom(e1_temp(n_prev),v1_t,j1_t,we1,xe1,be1,de1)
      CALL energy_diatom(e2_temp(n_prev),v2_t,j2_t,we2,xe2,be2,de2) 
      e_temp(n_prev) = e1_temp(n_prev) + e2_temp(n_prev)   	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      IF(emax_defined) THEN	  
      CALL	hpsort(nchann_est,e_temp,ja)
      nlvl = 0	  
      DO WHILE(.not.channels_defined)
      nlvl = nlvl + 1
      IF(nlvl.gt.nchann_est) THEN
      PRINT*, "SOMETHING WENT WORNG"
      STOP	  
      ENDIF	  
      IF(e_temp(nlvl).gt.EMAX) THEN
      number_of_channels = 
     & nlvl - 1
      channels_defined = .TRUE.	 
      ENDIF	 
      ENDDO
      ELSE
      CALL	hpsort(nchann_est,e1_temp,ja1)
      nlvl = 0	  
      DO WHILE(.not.channels_defined)
      nlvl = nlvl + 1
      IF(e1_temp(nlvl).gt.EMAX1) THEN
      nchann_1 = 
     & nlvl - 1
      channels_defined = .TRUE.	 
      ENDIF	 
      ENDDO
      channels_defined = .FALSE.
      CALL	hpsort(nchann_est,e2_temp,ja2)
      nlvl = 0	  
      DO WHILE(.not.channels_defined)
      nlvl = nlvl + 1
      IF(e2_temp(nlvl).gt.EMAX2) THEN
      nchann_2 = 
     & nlvl - 1
      channels_defined = .TRUE.	 
      ENDIF	 
      ENDDO
      number_of_channels = nchann_1*nchann_2	  
	  
      ENDIF	  
      IF(number_of_channels.le.0)STOP"ERROR:ZERO NUMBER OF CHANNELS"	  
!      PRINT*,  number_of_channels 
      IF(ALLOCATED(j1_ch)) DEALLOCATE (j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE (j2_ch)
      IF(ALLOCATED(v1_ch)) DEALLOCATE (j1_ch)
      IF(ALLOCATED(v2_ch)) DEALLOCATE (j2_ch)	  
      IF(ALLOCATED(E_ch)) DEALLOCATE (E_ch)
      ALLOCATE(j1_ch(number_of_channels),j2_ch(number_of_channels),
     & v1_ch(number_of_channels),v2_ch(number_of_channels),
     & E_ch(number_of_channels))
      IF(emax_defined) THEN	 
      DO nlvl = 1,number_of_channels
      n_prev = ja(nlvl)	  
      j1_ch(nlvl) = j1_ch_tmp(n_prev)
      v1_ch(nlvl) = v1_ch_temp(n_prev)
      j2_ch(nlvl) = j2_ch_tmp(n_prev)
      v2_ch(nlvl) = v2_ch_temp(n_prev)	  
      E_ch(nlvl) = e_temp(nlvl)
      ENDDO
      ELSE
      i = 0	  
      DO nlvl1 = 1,nchann_1
      DO nlvl2 = 1,nchann_2	
      i = i + 1		  
      n_prev1 = ja1(nlvl1)
      n_prev2 = ja2(nlvl2)	  
      j1_ch(i) = j1_ch_tmp(n_prev1)
      v1_ch(i) = v1_ch_temp(n_prev1)
      j2_ch(i) = j2_ch_tmp(n_prev2)
      v2_ch(i) = v2_ch_temp(n_prev2)	  
      E_ch(i) = e1_temp(nlvl1) + e2_temp(nlvl2)
      ENDDO
      ENDDO
      IF(i.ne.number_of_channels) STOP "ERROR IN SORTING"	  
      ENDIF
      DEALLOCATE(j1_ch_tmp,v1_ch_temp,e1_temp,ja1)
      DEALLOCATE(j2_ch_tmp,v2_ch_temp,e2_temp,ja2)
      DEALLOCATE(e_temp,ja)
      IF(MYID.eq.0) PRINT*, "ARRAYS DEALLOCATED"	  
      energy_defined = .TRUE.
      channels_defined = .TRUE.	  
      ENDIF

	  
      IF(jmax1.ge.0 .and. jmin1.ge.0 .and. vmax1.ge.0 .and. vmin1.ge.0 
     & .and.
     & jmax2.ge.0 .and. jmin2.ge.0 .and. vmax2.ge.0 .and. vmin2.ge.0 )
     & THEN
      number_of_channels_defined = .TRUE.
      energy_defined = .TRUE.
      channels_defined = .TRUE.
      number_of_channels = (jmax1 - jmin1 + 1)*(vmax1 - vmin1 + 1)*
     & (jmax2 - jmin2 + 1)*(vmax2 - vmin2 + 1)
      IF(ALLOCATED(j1_ch)) DEALLOCATE (j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE (j2_ch)
      IF(ALLOCATED(v1_ch)) DEALLOCATE (j1_ch)
      IF(ALLOCATED(v2_ch)) DEALLOCATE (j2_ch)	  
      IF(ALLOCATED(E_ch)) DEALLOCATE (E_ch)	 
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(v1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))
      ALLOCATE(v2_ch(number_of_channels))	  
      ALLOCATE(E_ch(number_of_channels))
      i =  0 	  
      DO v1_t=vmin1,vmax1
      DO v2_t=vmin2,vmax2      	  
      DO j1_t=jmin1,jmax1
      DO j2_t=jmin2,jmax2	  
      i = i + 1	  
      j1_ch(i) = j1_t
      v1_ch(i) = v1_t
      j2_ch(i) = j2_t
      v2_ch(i) = v2_t		  
      CALL energy_diatom(Ech1,v1_t,j1_t,we1,xe1,be1,de1)
      CALL energy_diatom(Ech2,v2_t,j2_t,we2,xe2,be2,de2)
      E_ch(i) = Ech1 + Ech2	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      IF(i.ne.number_of_channels) STOP"ERROR:PROGRAM ERROR"
      ELSE
      IF(.not.energy_defined .and. channels_defined) THEN	  
      IF(number_of_channels.le.0) 
     & STOP"ERROR:NUMBER_OF_CHANNLES MUST BE >0"
      number_of_channels_defined = .TRUE.
 	 
      ALLOCATE(E_ch(number_of_channels))
  
      DO i=1,number_of_channels
      j1_t = j1_ch(i) 
      v1_t = v1_ch(i)  
      j2_t = j2_ch(i)  
      v2_t = v2_ch(i) 		  
      CALL energy_diatom(Ech1,v1_t,j1_t,we1,xe1,be1,de1)
      CALL energy_diatom(Ech2,v2_t,j2_t,we2,xe2,be2,de2)
      E_ch(i) = Ech1 + Ech2	  
      ENDDO
      energy_defined = .TRUE.
      ENDIF
      ENDIF
      IF(number_of_channels_defined .and. .not. channels_defined) THEN
      nchann_est = (number_of_channels+1)**4
!      IF(MYID.eq.0) PRINT*, "nchann_est",nchann_est	 
      ALLOCATE(j1_ch_tmp(nchann_est),v1_ch_temp(nchann_est),
     & e1_temp(nchann_est),ja(nchann_est),ja1(nchann_est))
      ALLOCATE(j2_ch_tmp(nchann_est),v2_ch_temp(nchann_est),
     & e2_temp(nchann_est),e_temp(nchann_est),ja2(nchann_est))	 
      n_prev = 0	 
      DO v1_t = 0,number_of_channels
      DO j1_t = 0,number_of_channels
      DO v2_t = 0,number_of_channels
      DO j2_t = 0,number_of_channels 
      n_prev = n_prev + 1	  
      j1_ch_tmp(n_prev) = j1_t
      v1_ch_temp(n_prev) = v1_t
      j2_ch_tmp(n_prev) = j2_t
      v2_ch_temp(n_prev) = v2_t
      CALL energy_diatom(e1_temp(n_prev),v1_t,j1_t,we1,xe1,be1,de1)
      CALL energy_diatom(e2_temp(n_prev),v2_t,j2_t,we2,xe2,be2,de2) 
      e_temp(n_prev) = e1_temp(n_prev) + e2_temp(n_prev)   	  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      CALL hpsort(nchann_est,e_temp,ja)
      IF(ALLOCATED(j1_ch)) DEALLOCATE (j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE (j2_ch)
      IF(ALLOCATED(v1_ch)) DEALLOCATE (v1_ch)
      IF(ALLOCATED(v2_ch)) DEALLOCATE (v2_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE (E_ch)
      ALLOCATE(j1_ch(number_of_channels),j2_ch(number_of_channels)
     &  ,v1_ch(number_of_channels),v2_ch(number_of_channels), 
     & E_ch(number_of_channels))	  
      DO i=1,number_of_channels
      E_ch(i) = e_temp(i)
      j1_ch(i) = j1_ch_tmp(ja(i))
      j2_ch(i) = j2_ch_tmp(ja(i))
      v1_ch(i) = v1_ch_temp(ja(i))
      v2_ch(i) = v2_ch_temp(ja(i))		  
      ENDDO	  
      energy_defined = .TRUE.
      channels_defined = .TRUE.
      DEALLOCATE(j1_ch_tmp,v1_ch_temp,
     & e1_temp,ja,ja1)
      DEALLOCATE(j2_ch_tmp,v2_ch_temp,
     & e2_temp,e_temp,ja2)	 
      ENDIF	  
      IF(mlc_mlc_chn_num_defined .and. .not. channels_defined) THEN
      nchann1_est = (nchann_1+1)**2
!      IF(MYID.eq.0) PRINT*, "nchann_est",nchann_est	 
      ALLOCATE(j1_ch_tmp(nchann1_est),v1_ch_temp(nchann1_est),
     & e1_temp(nchann1_est),ja1(nchann1_est))
      n_prev = 0	 
      DO v1_t = 0,nchann_1
      DO j1_t = 0,nchann_1
      n_prev = n_prev + 1	  
      j1_ch_tmp(n_prev) = j1_t
      v1_ch_temp(n_prev) = v1_t
      CALL energy_diatom(e1_temp(n_prev),v1_t,j1_t,we1,xe1,be1,de1)
      ENDDO
      ENDDO	
      CALL hpsort(nchann1_est,e1_temp,ja1)
      nchann2_est = (nchann_2+1)**2	  
      ALLOCATE(j2_ch_tmp(nchann2_est),v2_ch_temp(nchann2_est),
     & e2_temp(nchann2_est),ja2(nchann2_est))
      n_prev = 0	 
      DO v2_t = 0,nchann_2
      DO j2_t = 0,nchann_2
      n_prev = n_prev + 1	  
      j2_ch_tmp(n_prev) = j2_t
      v2_ch_temp(n_prev) = v2_t
      CALL energy_diatom(e2_temp(n_prev),v2_t,j2_t,we2,xe2,be2,de2)
      ENDDO
      ENDDO	
      CALL hpsort(nchann2_est,e2_temp,ja2)
      IF(ALLOCATED(j1_ch)) DEALLOCATE (j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE (j2_ch)
      IF(ALLOCATED(v1_ch)) DEALLOCATE (j1_ch)
      IF(ALLOCATED(v2_ch)) DEALLOCATE (j2_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE (E_ch)
      ALLOCATE(j1_ch(number_of_channels),j2_ch(number_of_channels)
     &  ,v1_ch(number_of_channels),v2_ch(number_of_channels), 
     & E_ch(number_of_channels))
      n_prev =0
      DO i=1,nchann_1
      DO nlvl = 1,nchann_2
      n_prev = n_prev + 1	  
      E_ch(n_prev) = e_temp(i)
      j1_ch(n_prev) = j1_ch_tmp(ja1(i))
      j2_ch(n_prev) = j2_ch_tmp(ja2(nlvl))
      v1_ch(n_prev) = v1_ch_temp(ja1(i))
      v2_ch(n_prev) = v2_ch_temp(ja2(nlvl))
      E_ch(n_prev) = e1_temp(i) + e2_temp(nlvl)	  
      ENDDO	 
      ENDDO
      energy_defined = .true.
      channels_defined = .true.
      DEALLOCATE(j1_ch_tmp,j2_ch_tmp,v2_ch_temp,v1_ch_temp)
      DEALLOCATE(e1_temp,e2_temp,ja1,ja2)		  
      ENDIF	  
      ENDIF 
      IF(.not.energy_defined .or. .not. channels_defined) STOP
     & "ERROR IN SYSTEM INI"
      IF(atom_coord_dist) THEN
      IF(MYID.eq.0) THEN	  
      WRITE(1,*)	  
      WRITE(1,'(a28)') "PAIR WISE POTENTIAL DEFINED."
      WRITE(1,'(a49,f10.4,1x,a10,f10.4)')
     & "MASS OF ATOMS IN THE MOLECULE#1, A.M.U.: MASS1 = ",
     & atomic_masses(1),", MASS2 = ",atomic_masses(2)
      WRITE(1,*)
      WRITE(1,'(a49,f10.4,1x,a10,f10.4)')
     & "MASS OF ATOMS IN THE MOLECULE#2, A.M.U.: MASS1 = ",
     & atomic_masses(3),", MASS2 = ",atomic_masses(4)
      WRITE(1,*)	  
      ENDIF	  
      ENDIF
      IF(myid.eq.0) THEN	  
      WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"	  
      WRITE(1,'(a5,9x,a2,9x,a2,9x,a2,9x,a2,9x,a2)')"#N = "
     & ,"J1","V1","J2","V2","E"
      ENDIF	 
      jmax_included = 0d0
      IF(para_ortho_defined) CALL SYMMETRY_SORT_PARA_ORTHO	  
      DO i=1,number_of_channels
      IF(jmax_included.lt.j1_ch(i)) jmax_included = j1_ch(i)
      IF(jmax_included.lt.j2_ch(i)) jmax_included = j2_ch(i)	  
      IF(MYID.eq.0) WRITE(1,'(i5,i11,i11,i11,i11,f11.3)')i,
     & j1_ch(i),v1_ch(i), j2_ch(i),v2_ch(i),E_ch(i)
      ENDDO	  
      IF(MYID.eq.0) WRITE(1,*)
      IF(ini_chann_defined .and. .not.all_level_included) THEN
      chann_ini	 = 0  
      DO i=1,number_of_channels
      IF(j1_ini.eq.j1_ch(i) .and. v1_ini.eq.v1_ch(i)
     & .and. j2_ini.eq.j2_ch(i) .and. v2_ini.eq.v2_ch(i)) 
     & chann_ini = i
      ENDDO	
      IF(chann_ini.eq.0) STOP "ERROR: INITIAL CHANNEL WRONG"
      IF(MYID.eq.0)	  
     & WRITE(1,'(a17,a4,i4,a6,i4,a4,i4,a6,i4)') "INITIAL STATE IS:",
     & "J1 = ",j1_ini,"; V1 = ",v1_ini,"J2 = ",j2_ini,"; V2 = ",v2_ini
      ENDIF
      IF(.not.ini_chann_defined .and. channels_defined .and.
     & energy_defined) THEN
	  i_ground = 1
      DO i=1,number_of_channels
      IF(E_ch(i_ground).gt.E_ch(i)) i_ground = i   
      ENDDO	 
      chann_ini = i_ground
      ini_chann_defined = .TRUE.	  
      ENDIF  
      IF(MYID.eq.0) WRITE(1,*)
      IF(atomic_masses_defined) THEN
      atomic_red_mass1 = atomic_masses1(1)*atomic_masses1(2)/
     & (atomic_masses1(1)+atomic_masses1(2))
      IF(atomic_red_mass1.ne.atomic_red_mass1) THEN
      STOP "ERROR IN MASSES OF VIBRATING DIATOM #1"	  
      ENDIF
      atomic_red_mass2 = atomic_masses2(1)*atomic_masses2(2)/
     & (atomic_masses2(1)+atomic_masses2(2))
      IF(atomic_red_mass2.ne.atomic_red_mass2) THEN
      STOP "ERROR IN MASSES OF VIBRATING DIATOM #2"		  
      ENDIF
      IF(atomic_red_mass2*atomic_red_mass1.le.0d0) THEN
      STOP "ERROR:SUPPLY POSSITIVE MASSES FOR DIATOMS"	  
      ENDIF	  
      ENDIF
      IF(myid.eq.0) THEN	  
      IF(angs_unit)WRITE(1,'(a35)')"THE POTENTIAL R-UNITS ARE ANGSTROMS"
      IF(au_unit_r)WRITE(1,'(a31)') "THE POTENTIAL R-UNITS ARE BOHRS"
      IF(cm_unit)WRITE(1,'(a31)')"THE POTENTIAL ENERGY IS IN CM-1"
      IF(au_unit_e)WRITE(1,'(a34)')"THE POTENTIAL ENERGY IS IN HARTREE"
      IF(klv_unit)WRITE(1,'(a33)')"THE POTENTIAL ENERGY IS IN KELVIN"
      IF(kal_unit)WRITE(1,'(a35)')"THE POTENTIAL ENERGY IS IN KCAL/MOL"
      ENDIF	  
      IF(.not.angs_unit .and..not.au_unit_r )   
     & STOP "ERROR:R_UNITS NOT DEFINED"
      IF(.not.cm_unit .and..not.au_unit_e.and..not.klv_unit
     & .and. .not.kal_unit)   
     & STOP "ERROR:ENERGY NOT DEFINED"
      IF(myid.eq.0) WRITE(1,*)
      IF(expansion_defined) THEN
      IF(myid.eq.0) THEN	  
      WRITE(1,'(a42)') "USER SUPPLIED POTENTIAL TERMS WILL BE USED"
      WRITE(1,*)
      WRITE(1,'(a45,i4)')
     & "NUMBER OF POTENTIAL TERMS INCLUDED N_TERMS = ",nterms
      WRITE(1,'(a38)')
     & "THE POTENTIAL TERMS ARE, (L1,L2,L):"
      ENDIF	 
      DO i=1,nterms
      IF(max(L1_TERM(i),L1_TERM(i),L_TERM(i))*2 .gt.
     & L2_TERM(i)+L1_TERM(i)+L_TERM(i)) THEN
      IF(myid.eq.0) THEN	 
      WRITE(1,'(a47)')"ERROR: L1, L2 and L must follow triangular rule"
      WRITE(1,'(a3,1x,i3,1x,a3,1x,i3,1x,a3,1x,i3)') 
     & "L1=",L1_TERM(i),"L2=",L1_TERM(i),"L=",L_TERM(i)
      ENDIF	 
      ENDIF	  
      IF(myid.eq.0) WRITE(1,'(a1,i3,1x,i3,1x,i3,a3)',ADVANCE ="NO")
     & '(',L1_TERM(i),L2_TERM(i),L_TERM(i),'); '	  
      ENDDO	
      IF(myid.eq.0) THEN	  
      WRITE(1,*)
      WRITE(1,*)	  
      WRITE(1,'(a33,i3)') "NUMBER OF R-GRID POINTS IS N_R = ",
     & n_r_coll
      ENDIF	 
      ELSE
      IF(.not.matrix_reading_defined) THEN
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a98)')"COULPING MATRIX WILL BE COMPUTED USING NUMERICAL 
     &  INTEGRATION OF POTENTIAL OVER THE WAVEFUNCTIONS"
      WRITE(1,'(a58,i3)')
     & "NUMBER OF THETA-GRID OF #1 MOLECULE POINTS IS N_THETHA1 = ",
     & n_beta1
      WRITE(1,'(a58,i3)')
     & "NUMBER OF THETA-GRID OF #2 MOLECULE POINTS IS N_THETHA2 = ",
     & n_beta2
      WRITE(1,'(a51,i3)')
     & "NUMBER OF PHI-GRID OF  MOLECULES POINTS IS N_PHI = ",
     & n_gamma1 	 
      WRITE(1,'(a54,i3)')
     & "NUMBER OF R_VIB-GRID OF #1 MOLECULE POINTS IS N_VIB = ",
     & n_r_vib1
      WRITE(1,'(a54,i3)')
     & "NUMBER OF R_VIB-GRID OF #2 MOLECULE POINTS IS N_VIB = ",
     & n_r_vib2
      ENDIF 
      ELSE
      IF(MYID.eq.0) 
     & WRITE(1,'(a97)')"COULPING MATRIX WILL BE READ MIJ.dat FILE"	  
      ENDIF	  
      ENDIF 
      CASE(7)
      IF(MYID.eq.0)WRITE(1,*) "SYMMETRIC TOP + DIATOMC TOP"
      IF(user_defined) THEN	  
      IF(MYID.eq.0)	  
     & WRITE(1,'(a7,2x,a60)') "WARNING:",
     & "THE WAVEFUNCTIONS AND ENERGY LEVELS MUST BE SUPPLIED BY USER"
      	 
      IF(rot_const_defined)WRITE(1,'(a27)')"ROT. CONST WILL BE INGNORED"
      INQUIRE( FILE=USR_INP_LEVEL, EXIST=exst )
      IF(.not.exst)	THEN
      IF(MYID.eq.0)PRINT*,USR_INP_LEVEL, " NOT FOUND"
      STOP	  
      ENDIF	  
      OPEN(325,FILE=USR_INP_LEVEL,
     & STATUS="OLD",ACTION="READ")
      IF(number_of_channels.le.0) THEN
      PRINT*,"NUMBER OF CHANNELS NOT SPECIFIED"
      STOP	  
      ENDIF
      IF(.not.channels_defined) THEN	  
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(k1_ch(number_of_channels))
      ALLOCATE(eps1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))	  
      ENDIF
      IF(.not.energy_defined) THEN	  
      ALLOCATE(E_ch(number_of_channels))
      ENDIF
      READ(325,*) !HEADER
      jmax_included = 0	  
      DO i=1,number_of_channels
      READ(325,'(i5,1x,i4,1x,i4,1x,i2,1x,i3,e17.10)',
     & IOSTAT = stat_of_file)
     & decr,j1_ch(i), k1_ch(i),eps1_ch(i),j2_ch(i),E_ch(i)
      IF(jmax_included.lt.j1_ch(i))jmax_included=j1_ch(i)
      IF(jmax_included.lt.j2_ch(i))jmax_included=j2_ch(i)	 
	  
      IF(stat_of_file.gt.0 .or. decr.ne.i) THEN
      PRINT*,"WRONG CHANNELS IN USER INPUT FILE"	  
      STOP
      ENDIF 	  
      ENDDO
      CLOSE(325)
      channels_defined = .TRUE.
      energy_defined = .TRUE.
      ELSE
      IF(A+B.eq.0) STOP "ERROR: SUPPLY A OR B"
      IF(A.eq.0) A=B
      B=A
      IF(C.eq.0) STOP "ERROR : SYPPLY C"
      IF(Be.le.0) STOP "ERROR: SUPPLY BE"
      IF(DE.gt.BE) STOP "ERROR: DE MUST BE MUCH LESS THAN BE"
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a57,f7.2,a6,f7.2,a2)')
     & "ROTATIONAL CONSTANTS OF SYMMETRIC TOP ARE, CM-1: A = B = ",A
     &, "; C = ", C,";"
      WRITE(1,'(a46,f10.4)')
     &	  "ROT.CONSTANT OF DIATOMIC MOLECULE, CM-1, BE = ",Be	  
      WRITE(1,'(a23,f12.7)') "DISTORTION, CM-1, DE = ",DE
      ENDIF	  
      IF(number_of_channels_defined) THEN
      IF(energy_defined.and. .not.channels_defined) STOP 
     & "ERROR:FOR MANUALLY DEFINED ENERGY YOU SHOULD SPECIFY CHANNELS"
      IF(.not.energy_defined) ALLOCATE(E_ch(number_of_channels))
      IF(.not.channels_defined) THEN
      decr = 1	 
      nlvl = number_of_channels*decr
      ALLOCATE(j1_ch_tmp(nlvl),
     & k1_ch_tmp(nlvl),ja1(nlvl**2),
     & eps1_ch_tmp(nlvl),e1_temp(nlvl))
      ALLOCATE(j2_ch_tmp(nlvl),e2_temp(nlvl))	  	 
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
c      PRINT*,"nlvl==",nlvl	  
      ALLOCATE(j_ch_tmp(nlvl),
     & k_ch_tmp(nlvl),ja(nlvl),
     & eps_ch_tmp(nlvl),e_temp(nlvl))
	 
      i = 0	 
      DO j_t = 0,j_max_allowed
      DO k_t = 0,j_t
      DO eps_t = 0,1-KRONEKER(k_t,0)
      i = i+1    
      IF(i.gt.nlvl) EXIT
      j_ch_tmp(i) = j_t
      k_ch_tmp(i) = k_t
      eps_ch_tmp(i) = eps_t
      e_temp(i) = A*(j_t+1d0)*j_t - (A-C)*k_t**2 
      ENDDO
      IF(i.gt.nlvl) EXIT	  
      ENDDO	  
      IF(i.gt.nlvl) EXIT	  
      ENDDO
c      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
c      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,number_of_channels	  
      IF(e1_temp(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,number_of_channels
      j1_ch_tmp(i) = j_ch_tmp(ja(i))
      k1_ch_tmp(i) = k_ch_tmp(ja(i))
      eps1_ch_tmp(i) = eps_ch_tmp(ja(i))
      e1_temp(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = number_of_channels*decr      	  
      DEALLOCATE(k_ch_tmp,j_ch_tmp,ja,e_temp,eps_ch_tmp)	  
      ENDDO

      decr = 1	 
      nlvl = number_of_channels*decr
      channels_defined = .FALSE.	  
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
c      PRINT*,"nlvl==",nlvl	  
      ALLOCATE(j_ch_tmp(nlvl),e_temp(nlvl),ja(nlvl))
	 
      i = 0	 
      DO j_t = 0,j_max_allowed
      i = i+1    
      IF(i.gt.nlvl) EXIT
      j_ch_tmp(i) = j_t
      e_temp(i) =   Be*j_t*(j_t+1d0) -
     &	  De*(j_t*(j_t+1d0))**2 
      ENDDO

c      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
c      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,number_of_channels	  
      IF(e2_temp(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,number_of_channels
      j2_ch_tmp(i) = j_ch_tmp(ja(i))
      e2_temp(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = number_of_channels*decr      	  
      DEALLOCATE(j_ch_tmp,ja,e_temp)	  
      ENDDO	  
	  
      ALLOCATE(SORT_STORAGE(2,number_of_channels**2))
      ALLOCATE(e_temp(number_of_channels**2))
      n_prev = 0 	  
      DO i=1,number_of_channels
      DO nlvl=1,number_of_channels
      n_prev = n_prev +1 
      SORT_STORAGE(1,n_prev)  = i 
      SORT_STORAGE(2,n_prev)  = nlvl
      e_temp(n_prev) = e1_temp(i) + e2_temp(nlvl)
      ENDDO
      ENDDO	  
      CALL hpsort(number_of_channels**2,e_temp,ja1)
      IF(ALLOCATED(j1_ch)) DEALLOCATE(j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE(j2_ch)	  
      IF(ALLOCATED(k1_ch)) DEALLOCATE(k1_ch)
      IF(ALLOCATED(eps1_ch)) DEALLOCATE(eps1_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE(E_ch)
      ALLOCATE(j1_ch(number_of_channels),k1_ch(number_of_channels),
     & eps1_ch(number_of_channels),j2_ch(number_of_channels),
     & E_ch(number_of_channels) )	  
      DO i=1,number_of_channels
      nlvl = SORT_STORAGE(1,ja1(i))
      n_prev = SORT_STORAGE(2,ja1(i))	  
      j1_ch(i) = j1_ch_tmp(nlvl)
      k1_ch(i) = k1_ch_tmp(nlvl)
      eps1_ch(i) = eps1_ch_tmp(nlvl)
      j2_ch(i) = j2_ch_tmp(n_prev)
      E_ch(i) = e_temp(i)	  
      ENDDO
      DEALLOCATE(SORT_STORAGE,ja1,j1_ch_tmp,k1_ch_tmp,
     & eps1_ch_tmp,j2_ch_tmp,e_temp)	  
      channels_defined = .TRUE.
      energy_defined = .TRUE.	  
      ENDIF
	  
	  
	  
      IF(channels_defined .and. .not.energy_defined ) THEN
      DO i=1,number_of_channels	  
      E_ch(i) = A*(j1_ch(i)+1d0)*j1_ch(i) - (A-C)*k1_ch(i)**2 +
     & Be*j2_ch(i)*(j2_ch(i)+1d0) -
     &	  De*(j2_ch(i)*(j2_ch(i)+1d0))**2
      ENDDO
	  
      energy_defined = .TRUE.     	  
      ENDIF
!!!!! PRINTING
      ENDIF
      IF(mlc_mlc_chn_num_defined.and.mlc_mlc_emax_defined) 
     & STOP "ERROR:EITHER EMAX1,2 OR NCHLS1,2 MUST BE SPECIFIED"
      IF(mlc_mlc_chn_num_defined .and. .not.channels_defined) THEN
      number_of_channels = nchann_1*nchann_2
      IF(ALLOCATED(j1_ch)) DEALLOCATE(j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE(j2_ch)	  
      IF(ALLOCATED(k1_ch)) DEALLOCATE(k1_ch)
      IF(ALLOCATED(eps1_ch)) DEALLOCATE(eps1_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE(E_ch)  	  
      ALLOCATE(j1_ch(number_of_channels),k1_ch(number_of_channels),
     & eps1_ch(number_of_channels),j2_ch(number_of_channels),
     & E_ch(number_of_channels) )
      j1_ch = 0
      k1_ch = 0
      eps1_ch = 0
      decr = 1
      IF(max_nmb_chann.lt.number_of_channels) STOP "TOO MANY CHANNELS"	  
      nlvl = nchann_1*decr	  
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
c      PRINT*,"nlvl==",nlvl	  
      ALLOCATE(j_ch_tmp(nlvl),
     & k_ch_tmp(nlvl),ja(nlvl),
     & eps_ch_tmp(nlvl),e_temp(nlvl))
      i = 0	 
      DO j_t = 0,j_max_allowed
      DO k_t = 0,j_t
      DO eps_t = 0,1-KRONEKER(k_t,0)
      i = i+1    
      IF(i.gt.nlvl) EXIT
      j_ch_tmp(i) = j_t
      k_ch_tmp(i) = k_t
      eps_ch_tmp(i) = eps_t
      e_temp(i) = A*(j_t+1d0)*j_t - (A-C)*k_t**2 
      ENDDO
      IF(i.gt.nlvl) EXIT	  
      ENDDO	  
      IF(i.gt.nlvl) EXIT	  
      ENDDO
c      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
c      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,nchann_1	  
      IF(E_ch(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,nchann_1
      j1_ch(i) = j_ch_tmp(ja(i))
      k1_ch(i) = k_ch_tmp(ja(i))
      eps1_ch(i) = eps_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = nchann_1*decr      	  
      DEALLOCATE(k_ch_tmp,j_ch_tmp,ja,e_temp,eps_ch_tmp)	  
      ENDDO 


      DO nlvl = 1,nchann_2
      DO i=1,nchann_1
      j1_ch(i+nchann_1*(nlvl-1)) = j1_ch(i)
      k1_ch(i+nchann_1*(nlvl-1)) = k1_ch(i)
      eps1_ch(i+nchann_1*(nlvl-1)) = eps1_ch(i)
      j2_ch(i+nchann_1*(nlvl-1)) = nlvl-1	  
      E_ch(i+nchann_1*(nlvl-1))=E_ch(i)+Be*(nlvl-1)*((nlvl-1)+1d0)-
     &	  De*((nlvl-1)*((nlvl-1)+1d0))**2
      ENDDO	
      ENDDO

      energy_defined = .TRUE.	  
      ENDIF	  	  
      IF(mlc_mlc_emax_defined.and. .not.channels_defined) THEN
      IF(EMAX1.lt.0d0 .or. EMAX2.lt.0)
     & STOP"ERROR:EMAX1,2 MUST BE POSITIVE"	  
c      PRINT*,"EMAX=",EMAX	  
      j_t = int(abs(min(sqrt(EMAX1/C),(EMAX1/A-1)/2d0)))
      nchnl_ini_guess = (j_t+1)**2	  
c      PRINT*,"j_ini",j_t	  
      nchann_est = nchnl_ini_guess
      n_prev = -1
      nlvl = 0
      decr = 1
      DO WHILE(nlvl.ne.n_prev)
      n_prev =  nlvl	  
      nchann_est = nchnl_ini_guess*decr
c      PRINT*,"nchann_est = ",nchann_est	  
      ALLOCATE(j_ch_tmp(nchann_est),k_ch_tmp(nchann_est),
     & eps_ch_tmp(nchann_est),e_temp(nchann_est),ja(nchann_est))
      nlvl = 0
      DO i=0,j_max_allowed
      DO k_t=0,i
      DO eps_t = 0,1-KRONEKER(k_t,0)
      nlvl = nlvl + 1
      IF(nlvl.gt.nchann_est) EXIT	  
      j_ch_tmp(nlvl) = i
      k_ch_tmp(nlvl) = k_t
      eps_ch_tmp(nlvl) = eps_t
      e_temp(nlvl) = A*(i+1d0)*i - (A-C)*k_t**2 	  
      ENDDO
      IF(nlvl.gt.nchann_est) EXIT		  
      ENDDO
      IF(nlvl.gt.nchann_est) EXIT		  
      ENDDO
c      PRINT*,"nlvl",nlvl	  
      CALL	hpsort(nchann_est,e_temp,ja)
      i = 1	  
      DO WHILE(e_temp(i).lt.EMAX1)
      i = i + 1
c      PRINT*,"E",i,e_temp(i)	  
      ENDDO
      nlvl = i - 1	  
c      PRINT*,'nlvl',nlvl      	  
      decr = decr + 1 
      DEALLOCATE(j_ch_tmp,k_ch_tmp,
     & eps_ch_tmp,e_temp,ja)	  
      ENDDO
      nchann_1 = nlvl
      buffer_ch1 = nchann_est
	  
      nchann_est = int(sqrt(EMAX2/Be))
      i = nchann_est  
      IF(Be*i*(i+1d0) - De*(i*(i+1d0))**2 .lt. EMAX2) THEN
      DO WHILE(Be*i*(i+1d0) - De*(i*(i+1d0))**2 .lt. EMAX2)
      i = i +1
      ENDDO
      nchann_est = i - 1	  
      ENDIF
      IF(Be*i*(i+1d0) - De*(i*(i+1d0))**2 .GT. EMAX2) THEN
      DO WHILE(Be*i*(i+1d0) - De*(i*(i+1d0))**2 .GT. EMAX2)
      i = i -1
      ENDDO
      nchann_est = i  
      ENDIF
      nchann_2 = nchann_est	+ 1 
      number_of_channels = nchann_1*nchann_2
      ALLOCATE(j1_ch(number_of_channels),k1_ch(number_of_channels),
     & eps1_ch(number_of_channels),E_ch(number_of_channels),
     & j2_ch(number_of_channels)) 
       ALLOCATE(j_ch_tmp(buffer_ch1),k_ch_tmp(buffer_ch1),
     & eps_ch_tmp(buffer_ch1),e_temp(buffer_ch1),ja(buffer_ch1))
	 
      nlvl = 0
      DO i=0,j_max_allowed
      DO k_t=0,i
      DO eps_t = 0,1-KRONEKER(k_t,0)
      nlvl = nlvl + 1
      IF(nlvl.gt.buffer_ch1) EXIT	  
      j_ch_tmp(nlvl) = i
      k_ch_tmp(nlvl) = k_t
      eps_ch_tmp(nlvl) = eps_t
      e_temp(nlvl) = A*(i+1d0)*i - (A-C)*k_t**2 	  
      ENDDO
      IF(nlvl.gt.buffer_ch1) EXIT		  
      ENDDO
      IF(nlvl.gt.buffer_ch1) EXIT		  
      ENDDO
c      PRINT*,"number_of_channels = ",number_of_channels     
      CALL	hpsort(nchann_est,e_temp,ja)	  
      DO i=1,nchann_1
      j1_ch(i) = j_ch_tmp(ja(i))
      k1_ch(i) = k_ch_tmp(ja(i))
      eps1_ch(i) = eps_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      DEALLOCATE(j_ch_tmp,k_ch_tmp,
     & eps_ch_tmp,e_temp,ja)
      energy_defined = .TRUE.
      DO nlvl = 1,nchann_2
      DO i=1,nchann_1
      j1_ch(i+nchann_1*(nlvl-1)) = j1_ch(i)
      k1_ch(i+nchann_1*(nlvl-1)) = k1_ch(i)
      eps1_ch(i+nchann_1*(nlvl-1)) = eps1_ch(i)
      j2_ch(i+nchann_1*(nlvl-1)) = nlvl-1	  
      E_ch(i+nchann_1*(nlvl-1))=E_ch(i)+Be*(nlvl-1)*((nlvl-1)+1d0)-
     &	  De*((nlvl-1)*((nlvl-1)+1d0))**2
      ENDDO	
      ENDDO
	  
      channels_defined = .TRUE.
      energy_defined = .TRUE.
      number_of_channels_defined = .TRUE.
      ENDIF	  
	  
	  

      IF(emax_defined .and. .not. channels_defined) THEN
      nchann1_est = DSQRT(EMAX/min(A,C))+2
      nchann2_est = DSQRT(EMAX/BE)+2
      nchann_est = nchann1_est*nchann2_est
      ALLOCATE(j1_ch_tmp(nchann1_est),
     & k1_ch_tmp(nchann1_est),eps1_ch_tmp(nchann1_est),
     & e1_temp(nchann1_est))
      ALLOCATE(j2_ch_tmp(nchann2_est),e2_temp(nchann2_est),
     & e_temp(nchann_est),ja1(nchann_est),SORT_STORAGE(2,nchann_est))
      n_prev = 0 	 
      DO j_t=0,j_max_allowed
      DO k_t=0,j_t
      DO eps_t = 0,1-KRONEKER(k_t,0)
      IF(n_prev.gt.nchann1_est)	EXIT  
      n_prev = n_prev  + 1
      j1_ch_tmp(n_prev) =  j_t
      k1_ch_tmp(n_prev) =  k_t
      eps1_ch_tmp(n_prev) =  eps_t
      e1_temp(n_prev) = A*(j_t+1d0)*j_t - (A-C)*k_t**2
      ENDDO
      IF(n_prev.gt.nchann1_est)	EXIT 	  
      ENDDO
      IF(n_prev.gt.nchann1_est)	EXIT 	  
      ENDDO	  
	  nlvl = 0 
	  DO j_t = 0,	j_max_allowed
      nlvl = nlvl + 1
      IF(nlvl.gt.nchann2_est)	  EXIT
      j2_ch_tmp(nlvl) = j_t
      e2_temp(nlvl)  = Be*j_t*(j_t+1d0) -
     &	  De*(j_t*(j_t+1d0))**2     	  
      ENDDO 
      n_prev = 0	  
      DO i=1,nchann1_est
      DO nlvl=1,nchann2_est
      n_prev = n_prev + 1
      e_temp(n_prev) = e1_temp(i) + e2_temp(nlvl)
      SORT_STORAGE(1,n_prev) = i
      SORT_STORAGE(2,n_prev) = nlvl	  
      ENDDO
      ENDDO	

      CALL hpsort(nchann_est,e_temp,ja1)
      i = 1
      DO WHILE(e_temp(i).lt.EMAX .and. i.lt.nchann_est)
       i = i+ 1	       	  
      ENDDO
      number_of_channels = i-1
      IF(number_of_channels.lt.0) THEN
	  PRINT*,"ERROR:NUMBER OF CHANNELS",number_of_channels
      STOP
      ENDIF	
      IF(ALLOCATED(j1_ch)) DEALLOCATE(j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE(j2_ch)	  
      IF(ALLOCATED(k1_ch)) DEALLOCATE(k1_ch)
      IF(ALLOCATED(eps1_ch)) DEALLOCATE(eps1_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE(E_ch)	  
      ALLOCATE(j1_ch(number_of_channels),
     & j2_ch(number_of_channels),k1_ch(number_of_channels),
     & eps1_ch(number_of_channels),E_ch(number_of_channels))	 
      DO i=1,number_of_channels
      nlvl = SORT_STORAGE(1,ja1(i))
      n_prev = SORT_STORAGE(2,ja1(i))	  
      j1_ch(i) = j1_ch_tmp(nlvl)
      k1_ch(i) = k1_ch_tmp(nlvl)
      eps1_ch(i) = eps1_ch_tmp(nlvl)
      j2_ch(i) = j2_ch_tmp(n_prev)
      E_ch(i) = e_temp(i)	  
      ENDDO 
      channels_defined = .TRUE.
      energy_defined = .TRUE.	  
      ENDIF	  
	  
      ENDIF

	  
 !! PRINTINMG

      jmax_included = 0
      IF(para_ortho_defined) CALL SYMMETRY_SORT_PARA_ORTHO	  
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"	  
      WRITE(1,'(a5,9x,a2,9x,a2,9x,a2,9x,a2,9x,a2)')"#N = "
     & ,"J1","K","P","J2","E"
      ENDIF	 
      DO i=1,number_of_channels
      IF(j1_ch(i).lt.k1_ch(i))
     & STOP
     & "ERROR: J MUST BE NO LESS THAN K. CHECK THE INPUT CHANNELS."
      IF(k1_ch(i).eq.0 .and. eps1_ch(i).ne.0)
     &	  STOP"ERROR:K=0, SO PARITY MUST BE POSITIVE"
      IF(abs(eps1_ch(i)).gt.1) STOP "ERROR: PARITY MUST BE -1 OR 1"
      IF(jmax_included.lt.j1_ch(i)) jmax_included = j1_ch(i)
      IF(jmax_included.lt.j2_ch(i)) jmax_included = j2_ch(i)	  
      ENDDO	  

      DO i=1,number_of_channels
      IF(MYID.eq.0) WRITE(1,'(i5,i11,i11,i11,i11,e23.16)')i,
     & j1_ch(i),k1_ch(i),eps1_ch(i),j2_ch(i),E_ch(i) 	  
      ENDDO
  
	  
	  

      IF(ini_chann_defined .and. .not.all_level_included) THEN
      chann_ini = 0	  
      DO i=1,number_of_channels
      IF(j1_ini.eq.j1_ch(i) .and. k1_ini.eq.k1_ch(i) .and. 
     & eps1_ini.eq.eps1_ch(i).and.j2_ini.eq.j2_ch(i) )chann_ini = i
      ENDDO	
      IF(chann_ini.eq.0) STOP "ERROR: INITIAL CHANNEL WRONG"
      IF(MYID.eq.0)
     & WRITE(1,'(a17,a5,i4,a6,i4,a6,i4,a7,i4)') "INITIAL STATE IS:",
     & "J1 = ",j1_ini,"; K = ", k1_ini,"; P =	",eps1_ini,
     & "; J2 =	", j2_ini
      ENDIF
      IF(.not.ini_chann_defined .and. channels_defined .and.
     & energy_defined) THEN
	  i_ground = 1
      DO i=1,number_of_channels
      IF(E_ch(i_ground).gt.E_ch(i)) i_ground = i   
      ENDDO	 
      chann_ini = i_ground
      ini_chann_defined = .TRUE.	  
      ENDIF	  
      IF(myid.eq.0) THEN	  
      WRITE(1,*)
      IF(angs_unit)WRITE(1,'(a35)')"THE POTENTIAL R-UNITS ARE ANGSTROMS"
      IF(au_unit_r)WRITE(1,'(a31)') "THE POTENTIAL R-UNITS ARE BOHRS"
      IF(cm_unit)WRITE(1,'(a31)')"THE POTENTIAL ENERGY IS IN CM-1"
      IF(au_unit_e)WRITE(1,'(a34)')"THE POTENTIAL ENERGY IS IN HARTREE"
      IF(klv_unit)WRITE(1,'(a33)')"THE POTENTIAL ENERGY IS IN KELVIN"
      IF(kal_unit)WRITE(1,'(a35)')"THE POTENTIAL ENERGY IS IN KCAL/MOL"
      ENDIF	  
      IF(.not.angs_unit .and..not.au_unit_r )   
     & STOP "ERROR:R_UNITS NOT DEFINED"
      IF(.not.cm_unit .and..not.au_unit_e.and..not.klv_unit
     & .and. .not.kal_unit)   
     & STOP "ERROR:ENERGY NOT DEFINED"
      IF(myid.eq.0) WRITE(1,*) 	 
      IF(expansion_defined .and. myid.eq.0) THEN
      WRITE(1,'(a42)') "USER SUPPLIED POTENTIAL TERMS WILL BE USED"
      WRITE(1,*)
      WRITE(1,'(a45,i4)')
     & "NUMBER OF POTENTIAL TERMS INCLUDED N_TERMS = ",nterms
      WRITE(1,'(a40)')
     & "THE POTENTIAL TERMS ARE, (L1,NJU1,L2,L):"
      DO i=1,nterms
      WRITE(1,'(a1,i3,a1,i3,a1,i3,a1,i3,a3)',ADVANCE ="NO")
     & '(',L1_TERM(i),',',NJU1_TERM(i),',',L2_TERM(i),
     & ',',L_TERM(i),'); '	  
      ENDDO		
      WRITE(1,*)
      WRITE(1,*)	  
      WRITE(1,'(a33,i3)') "NUMBER OF R-GRID POINTS IS N_R = ",
     & n_r_coll	  
      ELSE
      IF(.not.matrix_reading_defined .and. myid.eq.0)	  
     & WRITE(1,'(a97)')"COULPING MATRIX WILL BE COMPUTED USING NUMERICAL
     &  INTEGRATION OF POTENTIAL OVER THE WAVEFUNCTIONS"  
      ENDIF	  
      IF(myid.eq.0)  WRITE(1,*)   
      CASE(8)
      IF(MYID.eq.0)WRITE(1,'(a23,1x,a29)')"THE COLLISIONAL SYSTEM:",
     & "ASYMMETRIC TOP + DIATOMIC TOP"
      IF(user_defined) THEN
      IF(MYID.eq.0)	  
     & WRITE(1,'(a7,2x,a60)') "WARNING:",
     & "THE WAVEFUNCTIONS AND ENERGY LEVELS MUST BE SUPPLIED BY USER"
      	 
      IF(rot_const_defined)WRITE(1,'(a27)')"ROT. CONST WILL BE INGNORED"
      INQUIRE( FILE=USR_INP_LEVEL, EXIST=exst )
      IF(.not.exst)	THEN
      IF(MYID.eq.0)PRINT*,USR_INP_LEVEL, " NOT FOUND"
      STOP	  
      ENDIF	  
      OPEN(325,FILE=USR_INP_LEVEL,
     & STATUS="OLD",ACTION="READ")
      IF(number_of_channels.le.0) THEN
      PRINT*,"NUMBER OF CHANNELS NOT SPECIFIED"
      STOP	  
      ENDIF
      IF(.not.channels_defined) THEN	  
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(ka1_ch(number_of_channels))
      ALLOCATE(kc1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))	  
      ENDIF
      IF(.not.energy_defined) THEN	  
      ALLOCATE(E_ch(number_of_channels))
      ENDIF
      READ(325,*) !HEADER
      jmax_included = 0	  
      DO i=1,number_of_channels
      READ(325,*,IOSTAT = stat_of_file)
     & decr,j1_ch(i), ka1_ch(i),kc1_ch(i),j2_ch(i),E_ch(i)
      IF(jmax_included.lt.j1_ch(i))jmax_included=j1_ch(i)
      IF(jmax_included.lt.j2_ch(i))jmax_included=j2_ch(i)	 
	  
      IF(stat_of_file.gt.0 .or. decr.ne.i) THEN
      PRINT*,"WRONG CHANNELS IN USER INPUT FILE"	  
      STOP
      ENDIF 	  
      ENDDO
      READ(325,*) !! HEADER FOR BASIS
      ALLOCATE(M_EIGEN_VECTORS_USER
     & (number_of_channels,jmax_included*2+1))
      M_EIGEN_VECTORS_USER = 0d0	 
      DO i=1,number_of_channels
      READ(325,*) decr,M_EIGEN_VECTORS_USER
     & (i,1:1+2*j1_ch(i))
      IF(i.ne.decr .or. stat_of_file.gt.0) STOP "ERROR IN USER FILE"	 
      ENDDO	  	  
      CLOSE(325)
      channels_defined = .TRUE.
      energy_defined = .TRUE.
      ELSE
      IF(A*B*C.le.0d0)
     & STOP"ERROR: ALL A,B AND C MUST BE POSITIVE. CHECK INPUT"
      IF(myid.eq.0)WRITE(1,'(a50,f7.2,a7,f7.2,a6,f7.2,a2)')
     & "ROTATIONAL CONSTANTS OF #1 MOLECULE ARE, CM-1: A ="
     & ,A," ; B = ",
     &  B,"; C = ", C,";"	 
      IF(Be.le.0) STOP "ERROR: SUPPLY BE"
      IF(DE.gt.BE) STOP "ERROR: DE MUST BE MUCH LESS THAN BE"
      IF(myid.eq.0) THEN	  
      WRITE(1,'(a46,f10.4)')
     &	  "ROT.CONSTANT OF DIATOMIC MOLECULE, CM-1, BE = ",Be	  
      WRITE(1,'(a23,f12.7)') "DISTORTION, CM-1, DE = ",DE
      ENDIF
      IF(mlc_mlc_emax_defined .and. .not.channels_defined) THEN
      IF(EMAX1.lt.0d0 .or. EMAX2.lt.0)
     & STOP"ERROR:EMAX1,2 MUST BE POSITIVE"	  
c      PRINT*,"EMAX=",EMAX1,EMAX2	  
      j_t = int(abs(min(sqrt(EMAX1/C),(EMAX1/A-1)/2d0)))
      nchnl_ini_guess = (j_t+1)**2	  
c      PRINT*,"j_ini",j_t	  
      nchann_est = nchnl_ini_guess
      n_prev = -1
      nlvl = 0
      decr = 1
      DO WHILE(nlvl.ne.n_prev)
      n_prev =  nlvl	  
      nchann_est = nchnl_ini_guess*decr
c      PRINT*,"nchann_est = ",nchann_est	  
      ALLOCATE(j_ch_tmp(nchann_est),ka_ch_tmp(nchann_est),
     & kc_ch_tmp(nchann_est),e_temp(nchann_est),ja(nchann_est))
      nlvl = 0
      DO j_t=0,j_max_allowed
      DO dk_t=-j_t,j_t
      nlvl = nlvl + 1
      IF(nlvl.gt.nchann_est) EXIT	  
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"
      j_ch_tmp(nlvl) = j_t
      ka_ch_tmp(nlvl) = ka_t
      kc_ch_tmp(nlvl) = kc_t
      CALL EIGEN_VAL(A,B,C,j_t,ka_t,kc_t,e_temp(nlvl))
      ENDDO
      IF(nlvl.gt.nchann_est) EXIT		  
      ENDDO
      CALL	hpsort(nchann_est,e_temp,ja)
      i = 1	  
      DO WHILE(e_temp(i).lt.EMAX1)
      i = i + 1
c      PRINT*,"E",i,e_temp(i)	  
      ENDDO
      nlvl = i - 1	  
c      PRINT*,'nlvl',nlvl      	  
      decr = decr + 1 
      DEALLOCATE(j_ch_tmp,ka_ch_tmp,
     & kc_ch_tmp,e_temp,ja)	  
      ENDDO
      nchann_1 = nlvl
      buffer_ch1 = nchann_est
c      PRINT*,nchann_1,buffer_ch1	  
      nchann_est = int(sqrt(EMAX2/Be))
      i = nchann_est  
      IF(Be*i*(i+1d0) - De*(i*(i+1d0))**2 .le. EMAX2) THEN
      DO WHILE(Be*i*(i+1d0) - De*(i*(i+1d0))**2 .le. EMAX2)
      i = i +1
      ENDDO
      nchann_est = i - 1	  
      ENDIF
      IF(Be*i*(i+1d0) - De*(i*(i+1d0))**2 .GT. EMAX2) THEN
      DO WHILE(Be*i*(i+1d0) - De*(i*(i+1d0))**2 .GT. EMAX2)
      i = i -1
      ENDDO
      nchann_est = i  
      ENDIF
      nchann_2 = nchann_est+1	  
      number_of_channels = nchann_1*nchann_2
c      PRINT*,number_of_channels,nchann_2,buffer_ch1
      ALLOCATE(j1_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),E_ch(number_of_channels),
     & j2_ch(number_of_channels)) 
       ALLOCATE(j_ch_tmp(buffer_ch1),ka_ch_tmp(buffer_ch1),
     & kc_ch_tmp(buffer_ch1),e_temp(buffer_ch1),ja(buffer_ch1))
	 
      nlvl = 0
      DO j_t=0,j_max_allowed
      DO dk_t=-j_t,j_t
      nlvl = nlvl + 1
      IF(nlvl.gt.buffer_ch1) EXIT	  
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst)
     &	  STOP"ERROR:ASSYM. TOP INITALIZATION FAILED FOR EMAX_DEFINED"
      j_ch_tmp(nlvl) = j_t
      ka_ch_tmp(nlvl) = ka_t
      kc_ch_tmp(nlvl) = kc_t
      CALL EIGEN_VAL(A,B,C,j_t,ka_t,kc_t,e_temp(nlvl))	  
      ENDDO
      IF(nlvl.gt.buffer_ch1) EXIT		  
      ENDDO
c      PRINT*,"number_of_channels_again = ",number_of_channels     
      CALL	hpsort(buffer_ch1,e_temp,ja)	  
      DO i=1,nchann_1
      j1_ch(i) = j_ch_tmp(ja(i))
      ka1_ch(i) = ka_ch_tmp(ja(i))
      kc1_ch(i) = kc_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      DEALLOCATE(j_ch_tmp,ka_ch_tmp,
     & kc_ch_tmp,e_temp,ja)
      energy_defined = .TRUE.	 
      DO nlvl = 1,nchann_2
      DO i=1,nchann_1
      j1_ch(i+nchann_1*(nlvl-1)) = j1_ch(i)
      ka1_ch(i+nchann_1*(nlvl-1)) = ka1_ch(i)
      kc1_ch(i+nchann_1*(nlvl-1)) = kc1_ch(i)
      j2_ch(i+nchann_1*(nlvl-1)) = nlvl-1	  
      E_ch(i+nchann_1*(nlvl-1))=E_ch(i)+Be*(nlvl-1)*((nlvl-1)+1d0)-
     &	  De*((nlvl-1)*((nlvl-1)+1d0))**2
      ENDDO	
      ENDDO
	  
      channels_defined = .TRUE.
      energy_defined = .TRUE.
      number_of_channels_defined = .TRUE.
      ENDIF	  	  
      IF(number_of_channels_defined) THEN
      IF(energy_defined.and. .not.channels_defined) STOP 
     & "ERROR:FOR MANUALLY DEFINED ENERGY YOU SHOULD SPECIFY CHANNELS"
      IF(.not.channels_defined) THEN 
      decr = 1	 
      nlvl = number_of_channels*decr
      ALLOCATE(j1_ch_tmp(nlvl),
     & ka1_ch_tmp(nlvl),ja1(nlvl**2),
     & kc1_ch_tmp(nlvl),e1_temp(nlvl))
      ALLOCATE(j2_ch_tmp(nlvl),e2_temp(nlvl))	  	 
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
      ALLOCATE(j_ch_tmp(nlvl),
     & ka_ch_tmp(nlvl),ja(nlvl),
     & kc_ch_tmp(nlvl),e_temp(nlvl))
      i = 0	 
      DO j_t=0,j_max_allowed
      DO dk_t=-j_t,j_t
      i = i + 1
      IF(i.gt.nlvl) EXIT	  
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"
      j_ch_tmp(i) = j_t
      ka_ch_tmp(i) = ka_t
      kc_ch_tmp(i) = kc_t
      CALL EIGEN_VAL(A,B,C,j_t,ka_t,kc_t,e_temp(i))
      ENDDO
      IF(i.gt.nlvl) EXIT		  
      ENDDO
!      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
!      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,number_of_channels	  
      IF(e1_temp(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,number_of_channels
      j1_ch_tmp(i) = j_ch_tmp(ja(i))
      ka1_ch_tmp(i) = ka_ch_tmp(ja(i))
      kc1_ch_tmp(i) = kc_ch_tmp(ja(i))
      e1_temp(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = number_of_channels*decr      	  
      DEALLOCATE(ka_ch_tmp,j_ch_tmp,ja,e_temp,kc_ch_tmp)	  
      ENDDO
!      IF(myid.eq.5) PRINT*,"CHANNELS DEFINED",e1_temp
      decr = 1	 
      nlvl = number_of_channels*decr
      channels_defined	= .FALSE.  
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
!      IF(Myid.eq.5)PRINT*,"nlvl==",nlvl	  
      ALLOCATE(j_ch_tmp(nlvl),e_temp(nlvl),ja(nlvl))
	 
      i = 0	 
      DO j_t = 0,j_max_allowed
      i = i+1    
      IF(i.gt.nlvl) EXIT
      j_ch_tmp(i) = j_t
      e_temp(i) =   Be*j_t*(j_t+1d0) -
     &	  De*(j_t*(j_t+1d0))**2 
      ENDDO

!      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
!      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,number_of_channels	  
      IF(e2_temp(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,number_of_channels
      j2_ch_tmp(i) = j_ch_tmp(ja(i))
      e2_temp(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = number_of_channels*decr      	  
      DEALLOCATE(j_ch_tmp,ja,e_temp)	  
      ENDDO	  
!      IF(myid.eq.5) PRINT*,myid,"CHANNELS DEFINED 2",e2_temp	  
      ALLOCATE(SORT_STORAGE(2,number_of_channels**2))
      ALLOCATE(e_temp(number_of_channels**2))
      SORT_STORAGE = 0
!      IF(myid.eq.5)PRINT*,myid," SORT_STORAGE",SORT_STORAGE		  
      n_prev = 0 	  
      DO i=1,number_of_channels
      DO nlvl=1,number_of_channels
      n_prev = n_prev +1
!      IF(myid.eq.4) PRINT*,myid,i,nlvl,n_prev 	  
      SORT_STORAGE(1,n_prev)  = i 
      SORT_STORAGE(2,n_prev)  = nlvl    	  
      e_temp(n_prev) = e1_temp(i) + e2_temp(nlvl)
!      IF(myid.eq.9) PRINT*, e_temp(n_prev),e1_temp(i),e2_temp(nlvl)	  
      ENDDO
	  
      ENDDO
!      IF(myid.eq.1)PRINT*,SORT_STORAGE	  
!      STOP	  
      CALL hpsort(number_of_channels**2,e_temp,ja1)
!      IF(myid.eq.0) PRINT*,"ENERYG SORTED DEFINED",e_temp		  
      IF(ALLOCATED(j1_ch)) DEALLOCATE(j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE(j2_ch)	  
      IF(ALLOCATED(ka1_ch)) DEALLOCATE(ka1_ch)
      IF(ALLOCATED(kc1_ch)) DEALLOCATE(kc1_ch)
      IF(ALLOCATED(E_ch))  DEALLOCATE(E_ch)


!      IF(myid.eq.0) PRINT*,"JS SORTED DEFINED",ja1		  
      ALLOCATE(j1_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),j2_ch(number_of_channels),
     & E_ch(number_of_channels) )
!      IF(myid.eq.0) PRINT*,"ARRAYS CREATED"	
!      STOP	  
      DO i=1,number_of_channels
      nlvl = SORT_STORAGE(1,ja1(i))
      n_prev = SORT_STORAGE(2,ja1(i))	  
      j1_ch(i) = j1_ch_tmp(nlvl)
      ka1_ch(i) = ka1_ch_tmp(nlvl)
      kc1_ch(i) = kc1_ch_tmp(nlvl)
      j2_ch(i) = j2_ch_tmp(n_prev)
      E_ch(i) = e_temp(i)	  
      ENDDO
      DEALLOCATE(SORT_STORAGE,ja1,j1_ch_tmp,ka1_ch_tmp,
     & kc1_ch_tmp,j2_ch_tmp,e_temp)	  
      channels_defined = .TRUE.
      energy_defined = .TRUE.
!      IF(MYID.eq.0) PRINT*,"ALL DEFINED"	  
      ENDIF
      IF(channels_defined .and. .not.energy_defined ) THEN
      IF(.not.ALLOCATED(E_ch)) ALLOCATE(E_ch(number_of_channels))	  
      DO i=1,number_of_channels
      j_t = j1_ch(i)
      ka_t = ka1_ch(i)
      kc_t = kc1_ch(i)	  
      CALL EIGEN_VAL(A,B,C,j_t,ka_t,kc_t,E_ch(i))	  
      E_ch(i) = E_ch(i) +
     & Be*j2_ch(i)*(j2_ch(i)+1d0) -
     &	  De*(j2_ch(i)*(j2_ch(i)+1d0))**2
      ENDDO
      energy_defined = .TRUE.     	  
      ENDIF
!!!!! PRINTING
      ENDIF
      IF(mlc_mlc_chn_num_defined.and.mlc_mlc_emax_defined) 
     & STOP "ERROR:EITHER EMAX1,2 OR NCHLS1,2 MUST BE SPECIFIED"
      IF(mlc_mlc_chn_num_defined.and. .not.channels_defined) THEN
      number_of_channels = nchann_1*nchann_2
      IF(MYID.eq.0) WRITE(1,'(a26,1x,i4)')
     & "#CHANNELS FOR 1ST MOLECULE",nchann_1
      IF(MYID.eq.0) WRITE(1,'(a26,1x,i4)')
     & "#CHANNELS FOR 2ST MOLECULE",nchann_2	  
      ALLOCATE(j1_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),j2_ch(number_of_channels),
     & E_ch(number_of_channels) )
      j1_ch = 0
      ka1_ch = 0
      kc1_ch = 0
      decr = 1
      IF(max_nmb_chann.lt.number_of_channels) STOP "TOO MANY CHANNELS"	  
      nlvl = nchann_1*decr	  
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
      ALLOCATE(j_ch_tmp(nlvl),ka_ch_tmp(nlvl),
     & kc_ch_tmp(nlvl),e_temp(nlvl),ja(nlvl))		  
      i = 0	 
      DO j_t = 0,j_max_allowed
      DO dk_t = -j_t,j_t
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"	  
      i = i+1    
      IF(i.gt.nlvl) EXIT
      j_ch_tmp(i) = j_t
      ka_ch_tmp(i) = ka_t
      kc_ch_tmp(i) = kc_t
c      PRINT*,"j_t,ka_t,kc_t="	, j_t, ka_t,kc_t	  
      CALL EIGEN_VAL(A,B,C,j_t,ka_t,kc_t,e_temp(i))
c      PRINT*,"e_temp,="	, e_temp(i) 
      ENDDO
      IF(i.gt.nlvl) EXIT	  
      ENDDO
c      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
c      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,nchann_1	  
      IF(E_ch(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,nchann_1
      j1_ch(i) = j_ch_tmp(ja(i))
      ka1_ch(i) = ka_ch_tmp(ja(i))
      kc1_ch(i) = kc_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = nchann_1*decr      	  
      DEALLOCATE(ka_ch_tmp,j_ch_tmp,ja,e_temp,kc_ch_tmp)	  
      ENDDO 

      ALLOCATE(j_ch_tmp(number_of_channels),
     & ka_ch_tmp(number_of_channels),
     & kc_ch_tmp(number_of_channels),
     & e_temp(number_of_channels))
      DO nlvl = 1,nchann_2
      DO i=1,nchann_1
      j1_ch(i+nchann_1*(nlvl-1)) = j1_ch(i)
      ka1_ch(i+nchann_1*(nlvl-1)) = ka1_ch(i)
      kc1_ch(i+nchann_1*(nlvl-1)) = kc1_ch(i)
      j2_ch(i+nchann_1*(nlvl-1)) = nlvl-1	  
      E_ch(i+nchann_1*(nlvl-1))=E_ch(i)+Be*(nlvl-1)*((nlvl-1)+1d0)-
     &	  De*((nlvl-1)*((nlvl-1)+1d0))**2	  
!      j_ch_tmp(i+nchann_1*(nlvl-1)) = j1_ch(i)
!      ka_ch_tmp(i+nchann_1*(nlvl-1)) = ka1_ch(i)
!      kc_ch_tmp(i+nchann_1*(nlvl-1)) = kc1_ch(i)
 !     j2_ch(i+nchann_1*(nlvl-1)) = nlvl-1	  
!      e_temp(i+nchann_1*(nlvl-1))=E_ch(i)+Be*(nlvl-1)*((nlvl-1)+1d0)-
!     &	  De*((nlvl-1)*((nlvl-1)+1d0))**2
      ENDDO	
      ENDDO
!      DEALLOCATE(j1_ch,ka1_ch,kc1_ch,E_ch)	
!      ALLOCATE(j1_ch(number_of_channels),ka1_ch(number_of_channels),
!     & kc1_ch(number_of_channels),j2_ch(number_of_channels),
!     & E_ch(number_of_channels) )
!      j1_ch =  j_ch_tmp	 
 !     ka1_ch  =	  ka_ch_tmp
!      kc1_ch = 	   kc_ch_tmp
!	  E_ch = e_temp
!      DEALLOCATE(j_ch_tmp,ka_ch_tmp,kc_ch_tmp,e_temp)		  
      number_of_channels_defined = .TRUE.
      energy_defined = .TRUE.	  
      ENDIF	

       IF(emax_defined .and. .not. channels_defined) THEN
      nchann1_est = INT(2*(DSQRT(EMAX/min(min(A,C),B))+1)**2)
      nchann2_est = DSQRT(EMAX/BE)+2
      nchann_est = nchann1_est*nchann2_est
!      IF(myid.eq.0) PRINT*,nchann_est,nchann1_est,nchann2_est	  
      ALLOCATE(j1_ch_tmp(nchann1_est),
     & ka1_ch_tmp(nchann1_est),kc1_ch_tmp(nchann1_est),
     & e1_temp(nchann1_est))
      ALLOCATE(j2_ch_tmp(nchann2_est),e2_temp(nchann2_est),
     & e_temp(nchann_est),ja1(nchann_est),SORT_STORAGE(2,nchann_est))
      i = 0 	 
      DO j_t = 0,j_max_allowed
      DO dk_t = -j_t,j_t
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"	  
      i = i+1    
      IF(i.gt.nchann1_est) EXIT
      j1_ch_tmp(i) = j_t
      ka1_ch_tmp(i) = ka_t
      kc1_ch_tmp(i) = kc_t
!      IF(myid.eq.0)PRINT*,"j_t,ka_t,kc_t="	, j_t, ka_t,kc_t	  
      CALL EIGEN_VAL(A,B,C,j_t,ka_t,kc_t,e1_temp(i))
!      IF(myid.eq.0)PRINT*,"e1_temp,="	, e1_temp(i) 
      ENDDO
      IF(i.gt.nchann1_est) EXIT	  
      ENDDO 
	  nlvl = 0 
	  DO j_t = 0,	j_max_allowed
      nlvl = nlvl + 1
      IF(nlvl.gt.nchann2_est)	  EXIT
      j2_ch_tmp(nlvl) = j_t
      e2_temp(nlvl)  = Be*j_t*(j_t+1d0) -
     &	  De*(j_t*(j_t+1d0))**2
!      IF(myid.eq.0)PRINT*,"e2_temp,="	, e2_temp(nlvl)	 
      ENDDO
 	  
      n_prev = 0	  
      DO i=1,nchann1_est
      DO nlvl=1,nchann2_est
      IF(n_prev.gt.nchann_est) EXIT	  
      n_prev = n_prev + 1
      e_temp(n_prev) = e1_temp(i) + e2_temp(nlvl)
      SORT_STORAGE(1,n_prev) = i
      SORT_STORAGE(2,n_prev) = nlvl	  
      ENDDO
      IF(n_prev.gt.nchann_est) EXIT	 	  
      ENDDO	
!      IF(myid.eq.0) PRINT*,SORT_STORAGE
      CALL hpsort(nchann_est,e_temp,ja1)
!      IF(myid.eq.0) PRINT*,e_temp	  
      i = 1
      DO WHILE(e_temp(i).lt.EMAX)
       i = i+ 1	       	  
      ENDDO
      number_of_channels = i-1
!      IF(myid.eq.0) PRINT*,"number_of_channels",number_of_channels	  
      IF(number_of_channels.lt.0) THEN
!	  PRINT*,"ERROR:NUMBER OF CHANNELS",number_of_channels
      STOP
      ENDIF	
      IF(ALLOCATED(j1_ch)) DEALLOCATE(j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE(j2_ch)	  
      IF(ALLOCATED(ka1_ch)) DEALLOCATE(k1_ch)
      IF(ALLOCATED(kc1_ch)) DEALLOCATE(kc1_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE(E_ch)	  
      ALLOCATE(j1_ch(number_of_channels),
     & j2_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),E_ch(number_of_channels))	 
      DO i=1,number_of_channels
      nlvl = SORT_STORAGE(1,ja1(i))
      n_prev = SORT_STORAGE(2,ja1(i))	  
      j1_ch(i) = j1_ch_tmp(nlvl)
      ka1_ch(i) = ka1_ch_tmp(nlvl)
      kc1_ch(i) = kc1_ch_tmp(nlvl)
      j2_ch(i) = j2_ch_tmp(n_prev)
      E_ch(i) = e_temp(i)	  
      ENDDO 
      channels_defined = .TRUE.
      energy_defined = .TRUE.
      DEALLOCATE(j1_ch_tmp,ka1_ch_tmp,j2_ch_tmp,kc1_ch_tmp,e_temp)	  
      ENDIF 	  
      ENDIF


 !! PRINTINMG
      jmax_included = 0
      IF(myid.eq.0) THEN	  
      WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"	  
      WRITE(1,'(a5,9x,a2,9x,a2,9x,a2,9x,a2,9x,a2)')"#N = "
     & ,"J1","KA","KC", "J2","E"
      ENDIF	
      IF(para_ortho_defined) CALL SYMMETRY_SORT_PARA_ORTHO	  
      DO i=1,number_of_channels
      IF(j1_ch(i).lt.max(ka1_ch(i),kc1_ch(i)))
     & STOP
     & "ERROR: J MUST BE NO LESS THAN K. CHECK THE INPUT CHANNELS."
      IF(.not. EXST_STATE(j1_ch(i),ka1_ch(i),kc1_ch(i))) THEN 
      WRITE(1,'(a48,1x,i4,1x,i4,1x,i4)')
     & "ERROR: CHECK CHANNELS. THIS STATE DOES NOT EXIST",
     & j1_ch(i),ka1_ch(i),kc1_ch(i)
      STOP
      ENDIF
      IF(jmax_included.lt.j1_ch(i)) jmax_included = j1_ch(i)
      IF(jmax_included.lt.j2_ch(i)) jmax_included = j2_ch(i)	  
      ENDDO	  
      DO i=1,number_of_channels
      IF(myid.eq.0)WRITE(1,'(i5,i11,i11,i11,i11,f11.3)')i,
     & j1_ch(i),ka1_ch(i),kc1_ch(i),j2_ch(i),E_ch(i) 	  
      ENDDO
      energy_defined = .TRUE.	  

      IF(ini_chann_defined .and. .not.all_level_included) THEN
      chann_ini = 0	  
      DO i=1,number_of_channels
      IF(j1_ini.eq.j1_ch(i) .and. ka1_ini.eq.ka1_ch(i) .and. 
     & kc1_ini.eq.kc1_ch(i).and.j2_ini.eq.j2_ch(i) )chann_ini = i
      ENDDO	
      IF(chann_ini.eq.0) THEN
      IF(MYID.eq.0) PRINT*,"ERROR: INITIAL CHANNEL WRONG",chann_ini
      STOP	   
      ENDIF	  
      IF(myid.eq.0)
     & WRITE(1,'(a17,a5,i4,a7,i4,a7,i4,a7,i4)') "INITIAL STATE IS:",
     & "J1 = ",j1_ini,"; KA = ", ka1_ini,"; KC =	",kc1_ini,
     & "; J2 =	", j2_ini
      ENDIF
      IF(.not.ini_chann_defined .and. channels_defined .and.
     & energy_defined) THEN
	  i_ground = 1
      DO i=1,number_of_channels
      IF(E_ch(i_ground).gt.E_ch(i)) i_ground = i   
      ENDDO	 
      chann_ini = i_ground
      ini_chann_defined = .TRUE.	  
      ENDIF	  
      IF(myid.eq.0) THEN	  
      WRITE(1,*)
      IF(angs_unit)WRITE(1,'(a35)')"THE POTENTIAL R-UNITS ARE ANGSTROMS"
      IF(au_unit_r)WRITE(1,'(a31)') "THE POTENTIAL R-UNITS ARE BOHRS"
      IF(cm_unit)WRITE(1,'(a31)')"THE POTENTIAL ENERGY IS IN CM-1"
      IF(au_unit_e)WRITE(1,'(a34)')"THE POTENTIAL ENERGY IS IN HARTREE"
      IF(klv_unit)WRITE(1,'(a33)')"THE POTENTIAL ENERGY IS IN KELVIN"
      IF(kal_unit)WRITE(1,'(a35)')"THE POTENTIAL ENERGY IS IN KCAL/MOL"
      ENDIF	  
      IF(.not.angs_unit .and..not.au_unit_r )   
     & STOP "ERROR:R_UNITS NOT DEFINED"
      IF(.not.cm_unit .and..not.au_unit_e.and..not.klv_unit
     & .and. .not.kal_unit)   
     & STOP "ERROR:ENERGY NOT DEFINED"
      IF(myid.eq.0)      WRITE(1,*) 	 
      IF(expansion_defined .and. myid.eq.0) THEN
      WRITE(1,'(a42)') "USER SUPPLIED POTENTIAL TERMS WILL BE USED"
      WRITE(1,*)
      WRITE(1,'(a45,i4)')
     & "NUMBER OF POTENTIAL TERMS INCLUDED N_TERMS = ",nterms
      WRITE(1,'(a40)')
     & "THE POTENTIAL TERMS ARE, (L1,NJU1,L2,L):"
      DO i=1,nterms
      WRITE(1,'(a1,i3,a1,i3,a1,i3,a1,i3,a3)',ADVANCE ="NO")
     & '(',L1_TERM(i),',',NJU1_TERM(i),',',L2_TERM(i),
     & ',',L_TERM(i),'); '	  
      ENDDO	
      WRITE(1,*)
      WRITE(1,*)	  
      WRITE(1,'(a33,i3)') "NUMBER OF R-GRID POINTS IS N_R = ",
     & n_r_coll	  
      ELSE
      IF(.not.matrix_reading_defined .and. myid.eq.0)	  
     & WRITE(1,'(a97)')"COULPING MATRIX WILL BE COMPUTED USING NUMERICAL
     &  INTEGRATION OF POTENTIAL OVER THE WAVEFUNCTIONS"	  
      ENDIF	  
      IF(myid.eq.0)      WRITE(1,*)   	 
      CASE(9)
      IF(MYID.eq.0)WRITE(1,'(a23,1x,a30)')"THE COLLISIONAL SYSTEM:",
     & "ASYMMETRIC TOP + SYMMETRIC TOP"
	  IF(user_defined) THEN
      IF(MYID.eq.0)	  
     & WRITE(1,'(a7,2x,a60)') "WARNING:",
     & "THE WAVEFUNCTIONS AND ENERGY LEVELS MUST BE SUPPLIED BY USER"
      	 
      IF(rot_const_defined)WRITE(1,'(a27)')"ROT. CONST WILL BE INGNORED"
      INQUIRE( FILE=USR_INP_LEVEL, EXIST=exst )
      IF(.not.exst)	THEN
      IF(MYID.eq.0)PRINT*,USR_INP_LEVEL, " NOT FOUND"
      STOP	  
      ENDIF	  
      OPEN(325,FILE=USR_INP_LEVEL,
     & STATUS="OLD",ACTION="READ")
      IF(number_of_channels.le.0) THEN
      PRINT*,"NUMBER OF CHANNELS NOT SPECIFIED"
      STOP	  
      ENDIF
      IF(.not.channels_defined) THEN	  
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(ka1_ch(number_of_channels))
      ALLOCATE(kc1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))
      ALLOCATE(k2_ch(number_of_channels))
      ALLOCATE(eps2_ch(number_of_channels))
	  
      ENDIF
      IF(.not.energy_defined) THEN	  
      ALLOCATE(E_ch(number_of_channels))
      ENDIF
      READ(325,*) !HEADER
      jmax_included = 0	  
      DO i=1,number_of_channels
      READ(325,*,IOSTAT = stat_of_file)
     & decr,j1_ch(i), ka1_ch(i),kc1_ch(i),
     & j2_ch(i),k2_ch(i),eps2_ch(i),E_ch(i)
      IF(jmax_included.lt.j1_ch(i))jmax_included=j1_ch(i)
      IF(jmax_included.lt.j2_ch(i))jmax_included=j2_ch(i)	 
	  
      IF(stat_of_file.gt.0 .or. decr.ne.i) THEN
      PRINT*,"WRONG CHANNELS IN USER INPUT FILE"	  
      STOP
      ENDIF 	  
      ENDDO
      READ(325,*) !! HEADER FOR BASIS
      ALLOCATE(M_EIGEN_VECTORS_USER
     & (number_of_channels,jmax_included*2+1))
      M_EIGEN_VECTORS_USER = 0d0	 
      DO i=1,number_of_channels
      READ(325,*) decr,M_EIGEN_VECTORS_USER
     & (i,1:1+2*j1_ch(i))
      IF(i.ne.decr .or. stat_of_file.gt.0) STOP "ERROR IN USER FILE"	 
      ENDDO	  	  
      CLOSE(325)
      channels_defined = .TRUE.
      energy_defined = .TRUE.
      ELSE
      IF(A1*B1*C1.le.0d0)
     & STOP"ERROR: ALL A1,B1 AND C1 MUST BE POSITIVE. CHECK INPUT"
      IF(myid.eq.0)WRITE(1,'(a51,f7.2,a8,f7.2,a6,f8.2,a2)')
     & "ROTATIONAL CONSTANTS OF #1 MOLECULE ARE, CM-1: A1 ="
     & ,A1," ; B1 = ",
     &  B1,"; C1 = ", C1,";"
      IF(A2+B2.eq.0) STOP "ERROR: SUPPLY A2 OR B2"
      IF(A2.eq.0) A2=B2
      B2=A2
      IF(C2.eq.0) STOP "ERROR : SYPPLY C2"
      IF(myid.eq.0)WRITE(1,'(a57,f7.2,a7,f7.2,a2)')
     & "ROTATIONAL CONSTANTS OF #2 MOLECULE ARE, CM-1: A2 = B2 = ",A2
     &, "; C2 = ", C2,";"
	 
      IF(number_of_channels_defined) THEN
      IF(energy_defined.and. .not.channels_defined) STOP 
     & "ERROR:FOR MANUALLY DEFINED ENERGY YOU SHOULD SPECIFY CHANNELS"
      IF(.not.energy_defined) ALLOCATE(E_ch(number_of_channels))
      IF(.not.channels_defined) THEN 
      decr = 1	 
      nlvl = number_of_channels*decr
      ALLOCATE(j1_ch_tmp(nlvl),
     & ka1_ch_tmp(nlvl),ja1(nlvl**2),
     & kc1_ch_tmp(nlvl),e1_temp(nlvl))
      ALLOCATE(j2_ch_tmp(nlvl),
     & k2_ch_tmp(nlvl),
     & eps2_ch_tmp(nlvl),e2_temp(nlvl))	  	 
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
      ALLOCATE(j_ch_tmp(nlvl),
     & ka_ch_tmp(nlvl),ja(nlvl),
     & kc_ch_tmp(nlvl),e_temp(nlvl))
      i = 0	 
      DO j_t=0,j_max_allowed
      DO dk_t=-j_t,j_t
      i = i + 1
      IF(i.gt.nlvl) EXIT	  
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"
      j_ch_tmp(i) = j_t
      ka_ch_tmp(i) = ka_t
      kc_ch_tmp(i) = kc_t
      CALL EIGEN_VAL(A1,B1,C1,j_t,ka_t,kc_t,e_temp(i))
      ENDDO
      IF(i.gt.nlvl) EXIT		  
      ENDDO
!      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
!      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,number_of_channels	  
      IF(e1_temp(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,number_of_channels
      j1_ch_tmp(i) = j_ch_tmp(ja(i))
      ka1_ch_tmp(i) = ka_ch_tmp(ja(i))
      kc1_ch_tmp(i) = kc_ch_tmp(ja(i))
      e1_temp(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = number_of_channels*decr      	  
      DEALLOCATE(ka_ch_tmp,j_ch_tmp,ja,e_temp,kc_ch_tmp)	  
      ENDDO
!      IF(myid.eq.5) PRINT*,"CHANNELS DEFINED",e1_temp
      decr = 1	 
      nlvl = number_of_channels*decr
      channels_defined	= .FALSE.  
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
!      PRINT*,"nlvl==",nlvl	  
      ALLOCATE(j_ch_tmp(nlvl),
     & k_ch_tmp(nlvl),ja(nlvl),
     & eps_ch_tmp(nlvl),e_temp(nlvl))
	 
      i = 0	 
      DO j_t = 0,j_max_allowed
      DO k_t = 0,j_t
      DO eps_t = 0,1-KRONEKER(k_t,0)
      i = i+1    
      IF(i.gt.nlvl) EXIT
      j_ch_tmp(i) = j_t
      k_ch_tmp(i) = k_t
      eps_ch_tmp(i) = eps_t
      e_temp(i) = A2*(j_t+1d0)*j_t - (A2-C2)*k_t**2 
      ENDDO
      IF(i.gt.nlvl) EXIT	  
      ENDDO	  
      IF(i.gt.nlvl) EXIT	  
      ENDDO
!      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
!      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,number_of_channels	  
      IF(e1_temp(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,number_of_channels
      j2_ch_tmp(i) = j_ch_tmp(ja(i))
      k2_ch_tmp(i) = k_ch_tmp(ja(i))
      eps2_ch_tmp(i) = eps_ch_tmp(ja(i))
      e2_temp(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = number_of_channels*decr      	  
      DEALLOCATE(k_ch_tmp,j_ch_tmp,ja,e_temp,eps_ch_tmp)	  
      ENDDO  
!      IF(myid.eq.5) PRINT*,myid,"CHANNELS DEFINED 2",e2_temp	  
      ALLOCATE(SORT_STORAGE(2,number_of_channels**2))
      ALLOCATE(e_temp(number_of_channels**2))
      SORT_STORAGE = 0
!      IF(myid.eq.5)PRINT*,myid," SORT_STORAGE",SORT_STORAGE		  
      n_prev = 0 	  
      DO i=1,number_of_channels
      DO nlvl=1,number_of_channels
      n_prev = n_prev +1
!      IF(myid.eq.4) PRINT*,myid,i,nlvl,n_prev 	  
      SORT_STORAGE(1,n_prev)  = i 
      SORT_STORAGE(2,n_prev)  = nlvl    	  
      e_temp(n_prev) = e1_temp(i) + e2_temp(nlvl)
!      IF(myid.eq.9) PRINT*, e_temp(n_prev),e1_temp(i),e2_temp(nlvl)	  
      ENDDO
	  
      ENDDO
!      IF(myid.eq.1)PRINT*,SORT_STORAGE	  
!      STOP	  
      CALL hpsort(number_of_channels**2,e_temp,ja1)
!      IF(myid.eq.0) PRINT*,"ENERYG SORTED DEFINED",e_temp		  
      IF(ALLOCATED(j1_ch)) DEALLOCATE(j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE(j2_ch)	  
      IF(ALLOCATED(ka1_ch)) DEALLOCATE(ka1_ch)
      IF(ALLOCATED(kc1_ch)) DEALLOCATE(kc1_ch)
      IF(ALLOCATED(k2_ch)) DEALLOCATE(k2_ch)
      IF(ALLOCATED(eps2_ch)) DEALLOCATE(eps2_ch)	  
      IF(ALLOCATED(E_ch))  DEALLOCATE(E_ch)


!      IF(myid.eq.0) PRINT*,"JS SORTED DEFINED",ja1		  
      ALLOCATE(j1_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),j2_ch(number_of_channels),
     & k2_ch(number_of_channels),eps2_ch(number_of_channels),
     & E_ch(number_of_channels) )
!      IF(myid.eq.0) PRINT*,"ARRAYS CREATED"	
!      STOP	  
      DO i=1,number_of_channels
      nlvl = SORT_STORAGE(1,ja1(i))
      n_prev = SORT_STORAGE(2,ja1(i))	  
      j1_ch(i) = j1_ch_tmp(nlvl)
      ka1_ch(i) = ka1_ch_tmp(nlvl)
      kc1_ch(i) = kc1_ch_tmp(nlvl)
      j2_ch(i) = j2_ch_tmp(n_prev)
      k2_ch(i) = k2_ch_tmp(n_prev)
      eps2_ch(i) = eps2_ch_tmp(n_prev)	  
      E_ch(i) = e_temp(i)	  
      ENDDO
      DEALLOCATE(SORT_STORAGE,ja1,j1_ch_tmp,ka1_ch_tmp,
     & kc1_ch_tmp,j2_ch_tmp,e_temp,k2_ch_tmp,eps2_ch_tmp)	  
      channels_defined = .TRUE.
      energy_defined = .TRUE.
!      IF(MYID.eq.0) PRINT*,"ALL DEFINED"	  
      ENDIF
      IF(channels_defined .and. .not.energy_defined ) THEN
      DO i=1,number_of_channels
      j_t = j1_ch(i)
      ka_t = ka1_ch(i)
      kc_t = kc1_ch(i)	  
      CALL EIGEN_VAL(A1,B1,C1,j_t,ka_t,kc_t,E_ch(i))	  
      E_ch(i) = E_ch(i) + A2*(j2_ch(i)+1d0)*j2_ch(i)-
     & (A2-C2)*k2_ch(i)**2
      ENDDO
	  
	  
      energy_defined = .TRUE.     	  
      ENDIF
!!!!! PRINTING
      ENDIF
      IF(mlc_mlc_chn_num_defined.and.mlc_mlc_emax_defined) 
     & STOP "ERROR:EITHER EMAX1,2 OR NCHLS1,2 MUST BE SPECIFIED"
      IF(mlc_mlc_chn_num_defined .and. .not. channels_defined) THEN
      number_of_channels = nchann_1*nchann_2
      ALLOCATE(j1_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),j2_ch(number_of_channels),
     & E_ch(number_of_channels),k2_ch(number_of_channels),
     & eps_ch(number_of_channels))
      j1_ch = 0
      ka1_ch = 0
      kc1_ch = 0
      decr = 1
      IF(max_nmb_chann.lt.number_of_channels) STOP "TOO MANY CHANNELS"	  
      nlvl = nchann_1*decr	  
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
      i = 0	 
      DO j_t = 0,j_max_allowed
      DO dk_t = -j_t,j_t
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"	  
      i = i+1    
      IF(i.gt.nlvl) EXIT
      j_ch_tmp(i) = j_t
      ka_ch_tmp(i) = ka_t
      kc_ch_tmp(i) = kc_t
c      PRINT*,"j_t,ka_t,kc_t="	, j_t, ka_t,kc_t	  
      CALL EIGEN_VAL(A1,B1,C1,j_t,ka_t,kc_t,e_temp(i))
c      PRINT*,"e_temp,="	, e_temp(i) 
      ENDDO
      IF(i.gt.nlvl) EXIT	  
      ENDDO
c      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
c      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,nchann_1	  
      IF(E_ch(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,nchann_1
      j1_ch(i) = j_ch_tmp(ja(i))
      ka1_ch(i) = ka_ch_tmp(ja(i))
      kc1_ch(i) = kc_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = nchann_1*decr      	  
      DEALLOCATE(ka_ch_tmp,j_ch_tmp,ja,e_temp,kc_ch_tmp)	  
      ENDDO 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      decr = 1
      nlvl = nchann_2*decr
      channels_defined = .FALSE.	  
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
c      PRINT*,"nlvl==",nlvl	  
      ALLOCATE(j_ch_tmp(nlvl),
     & k_ch_tmp(nlvl),ja(nlvl),
     & eps_ch_tmp(nlvl),e_temp(nlvl))
      i = 0	 
      DO j_t = 0,j_max_allowed
      DO k_t = 0,j_t
      DO eps_t = 0,1-KRONEKER(k_t,0)
      i = i+1    
      IF(i.gt.nlvl) EXIT
      j_ch_tmp(i) = j_t
      k_ch_tmp(i) = k_t
      eps_ch_tmp(i) = eps_t
      e_temp(i) = A2*(j_t+1d0)*j_t - (A2-C2)*k_t**2 
      ENDDO
      IF(i.gt.nlvl) EXIT	  
      ENDDO	  
      IF(i.gt.nlvl) EXIT	  
      ENDDO
c      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
c      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,nchann_2	  
      IF(E_ch(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,nchann_2!!!!!!!!!!!!!!!!!!!!!!!!!!HERE
      j2_ch((i-1)*nchann_1+1) = j_ch_tmp(ja(i))
      k2_ch((i-1)*nchann_1+1) = k_ch_tmp(ja(i))
      eps2_ch((i-1)*nchann_1+1) = eps_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = nchann_1*decr      	  
      DEALLOCATE(k_ch_tmp,j_ch_tmp,ja,e_temp,eps_ch_tmp)	  
      ENDDO 	  
	  
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO nlvl = 1,nchann_2
      DO i=1,nchann_1
      j1_ch(i+nchann_1*(nlvl-1)) = j1_ch(i)
c      PRINT*,"j1",i,j1_ch(i)	  
      ka1_ch(i+nchann_1*(nlvl-1)) = ka1_ch(i)
c      PRINT*,"ka1",i,ka1_ch(i)	  
      kc1_ch(i+nchann_1*(nlvl-1)) = kc1_ch(i)
c      PRINT*,"kc1",i,kc1_ch(i)		  
      j2_ch(i+nchann_1*(nlvl-1)) =  j2_ch(1+nchann_1*(nlvl-1))
c      PRINT*,"j2",i,j2_ch(i)			  
      k2_ch(i+nchann_1*(nlvl-1)) =  k2_ch(1+nchann_1*(nlvl-1))
c      PRINT*,"k2",i,k2_ch(i)		  
      eps2_ch(i+nchann_1*(nlvl-1)) =  eps2_ch(1+nchann_1*(nlvl-1))
c      PRINT*,"p2",i,eps2_ch(i)		  
      j_t = j1_ch(i)
      ka_t = ka1_ch(i)
      kc_t = kc1_ch(i)	  
      CALL EIGEN_VAL(A1,B1,C1,j_t,ka_t,kc_t,Ech1)
      j_t = j2_ch(1+nchann_1*(nlvl-1))
      k_t =  k2_ch(1+nchann_1*(nlvl-1))		  
      Ech2 =  A2*(j_t+1d0)*j_t-
     & (A2-C2)*k_t**2
      E_ch(i+nchann_1*(nlvl-1)) = Ech2 + Ech1
c      PRINT*,i, E_ch(i)	 
      ENDDO	
      ENDDO
      number_of_channels_defined = .TRUE.
      energy_defined = .TRUE.	  
      ENDIF
      IF(mlc_mlc_emax_defined.and. .not.channels_defined) THEN
      IF(EMAX1.lt.0d0 .or. EMAX2.lt.0)
     & STOP"ERROR:EMAX1,2 MUST BE POSITIVE"	  
c      PRINT*,"EMAX=",EMAX1,EMAX2	  
      j_t = int(abs(min(sqrt(EMAX1/C1),(EMAX1/A1-1)/2d0)))
      nchnl_ini_guess = (j_t+1)**2	  
c      PRINT*,"j_ini",j_t	  
      nchann_est = nchnl_ini_guess
      n_prev = -1
      nlvl = 0
      decr = 1
      DO WHILE(nlvl.ne.n_prev)
      n_prev =  nlvl	  
      nchann_est = nchnl_ini_guess*decr
c      PRINT*,"nchann_est = ",nchann_est	  
      ALLOCATE(j_ch_tmp(nchann_est),ka_ch_tmp(nchann_est),
     & kc_ch_tmp(nchann_est),e_temp(nchann_est),ja(nchann_est))
      nlvl = 0
      DO j_t=0,j_max_allowed
      DO dk_t=-j_t,j_t
      nlvl = nlvl + 1
      IF(nlvl.gt.nchann_est) EXIT	  
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"
      j_ch_tmp(nlvl) = j_t
      ka_ch_tmp(nlvl) = ka_t
      kc_ch_tmp(nlvl) = kc_t
      CALL EIGEN_VAL(A1,B1,C1,j_t,ka_t,kc_t,e_temp(nlvl))
      ENDDO
      IF(nlvl.gt.nchann_est) EXIT		  
      ENDDO
      CALL	hpsort(nchann_est,e_temp,ja)
      i = 1	  
      DO WHILE(e_temp(i).lt.EMAX1)
      i = i + 1
c      PRINT*,"E",i,e_temp(i)	  
      ENDDO
      nlvl = i - 1	  
c      PRINT*,'nlvl',nlvl      	  
      decr = decr + 1 
      DEALLOCATE(j_ch_tmp,ka_ch_tmp,
     & kc_ch_tmp,e_temp,ja)	  
      ENDDO
      nchann_1 = nlvl
      buffer_ch1 = nchann_est
c      PRINT *,"nchann_1",nchann_1,"buffer_ch1",buffer_ch1 	  
!!!!! FOR SYMMETRIC TOP	  
      j_t = int(abs(min(sqrt(EMAX2/C2),(EMAX2/A2-1)/2d0)))
      nchnl_ini_guess = (j_t+1)**2	  
c      PRINT*,"j_ini",j_t	  
      nchann_est = nchnl_ini_guess
      n_prev = -1
      nlvl = 0
      decr = 1
      DO WHILE(nlvl.ne.n_prev)
      n_prev =  nlvl	  
      nchann_est = nchnl_ini_guess*decr
c      PRINT*,"nchann_est = ",nchann_est	  
      ALLOCATE(j_ch_tmp(nchann_est),k_ch_tmp(nchann_est),
     & eps_ch_tmp(nchann_est),e_temp(nchann_est),ja(nchann_est))
      nlvl = 0
      DO i=0,j_max_allowed
      DO k_t=0,i
      DO eps_t = 0,1-KRONEKER(k_t,0)
      nlvl = nlvl + 1
      IF(nlvl.gt.nchann_est) EXIT	  
      j_ch_tmp(nlvl) = i
      k_ch_tmp(nlvl) = k_t
      eps_ch_tmp(nlvl) = eps_t
      e_temp(nlvl) = A2*(i+1d0)*i - (A2-C2)*k_t**2 	  
      ENDDO
      IF(nlvl.gt.nchann_est) EXIT		  
      ENDDO
      IF(nlvl.gt.nchann_est) EXIT		  
      ENDDO
c      PRINT*,"nlvl",nlvl	  
      CALL	hpsort(nchann_est,e_temp,ja)
      i = 1	  
      DO WHILE(e_temp(i).lt.EMAX2)
      i = i + 1
c      PRINT*,"E",i,e_temp(i)	  
      ENDDO
      nlvl = i - 1	  
c      PRINT*,'nlvl',nlvl      	  
      decr = decr + 1 
      DEALLOCATE(j_ch_tmp,k_ch_tmp,
     & eps_ch_tmp,e_temp,ja)	  
      ENDDO
      nchann_2 = nlvl
      buffer_ch2 = nchann_est
c      PRINT *,"nchann_2",nchann_2,"buffer_ch2",buffer_ch2  
      number_of_channels = nchann_1*nchann_2
	  !end
c      PRINT*,number_of_channels,nchann_2,buffer_ch1
      ALLOCATE(j1_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),E_ch(number_of_channels),
     & j2_ch(number_of_channels),eps2_ch(number_of_channels),
     & k2_ch(number_of_channels))	 
       ALLOCATE(j_ch_tmp(buffer_ch1),ka_ch_tmp(buffer_ch1),
     & kc_ch_tmp(buffer_ch1),e_temp(buffer_ch1),ja(buffer_ch1))
	 
      nlvl = 0
      DO j_t=0,j_max_allowed
      DO dk_t=-j_t,j_t
      nlvl = nlvl + 1
      IF(nlvl.gt.buffer_ch1) EXIT	  
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst)
     &	  STOP"ERROR:ASSYM. TOP INITALIZATION FAILED FOR EMAX_DEFINED"
      j_ch_tmp(nlvl) = j_t
      ka_ch_tmp(nlvl) = ka_t
      kc_ch_tmp(nlvl) = kc_t
      CALL EIGEN_VAL(A1,B1,C1,j_t,ka_t,kc_t,e_temp(nlvl))	  
      ENDDO
      IF(nlvl.gt.buffer_ch1) EXIT		  
      ENDDO
c      PRINT*,"number_of_channels_again = ",number_of_channels     
      CALL	hpsort(buffer_ch1,e_temp,ja)	  
      DO i=1,nchann_1
      j1_ch(i) = j_ch_tmp(ja(i))
      ka1_ch(i) = ka_ch_tmp(ja(i))
      kc1_ch(i) = kc_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      DEALLOCATE(j_ch_tmp,ka_ch_tmp,
     & kc_ch_tmp,e_temp,ja)
      energy_defined = .TRUE.
c      PRINT*,"j1_ch",j1_ch	  
!!!!!!!!!!!!! for assymetric top
      ALLOCATE(j_ch_tmp(buffer_ch2),k_ch_tmp(buffer_ch2),
     & eps_ch_tmp(buffer_ch2),e_temp(buffer_ch2),ja(buffer_ch2))
      nlvl = 0
      DO i=0,j_max_allowed
      DO k_t=0,i
      DO eps_t = 0,1-KRONEKER(k_t,0)
      nlvl = nlvl + 1
      IF(nlvl.gt.buffer_ch2) EXIT	  
      j_ch_tmp(nlvl) = i
      k_ch_tmp(nlvl) = k_t
      eps_ch_tmp(nlvl) = eps_t
      e_temp(nlvl) = A2*(i+1d0)*i - (A2-C2)*k_t**2 	  
      ENDDO
      IF(nlvl.gt.buffer_ch2) EXIT		  
      ENDDO
      IF(nlvl.gt.buffer_ch2) EXIT		  
      ENDDO
c      PRINT*,"nlvl",nlvl	  
      CALL	hpsort(buffer_ch2,e_temp,ja)
      i = 1	  
      DO WHILE(e_temp(i).lt.EMAX2)
      i = i + 1
c      PRINT*,"E",i,e_temp(i)	  
      ENDDO
      nlvl = i - 1	  
c      PRINT*,'nlvl',nlvl      	  
  

      DO i=1,nchann_2!!!!!!!!!!!!!!!!!!!!!!!!!!HERE
      j2_ch((i-1)*nchann_1+1) = j_ch_tmp(ja(i))
      k2_ch((i-1)*nchann_1+1) = k_ch_tmp(ja(i))
      eps2_ch((i-1)*nchann_1+1) = eps_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      DEALLOCATE(j_ch_tmp,k_ch_tmp,
     & eps_ch_tmp,e_temp,ja)	
c      PRINT *,"j2=",j2_ch  
      DO nlvl = 1,nchann_2
      DO i=1,nchann_1
      j1_ch(i+nchann_1*(nlvl-1)) = j1_ch(i)
c      PRINT*,"j1",i,j1_ch(i)	  
      ka1_ch(i+nchann_1*(nlvl-1)) = ka1_ch(i)
c      PRINT*,"ka1",i,ka1_ch(i)	  
      kc1_ch(i+nchann_1*(nlvl-1)) = kc1_ch(i)
c      PRINT*,"kc1",i,kc1_ch(i)		  
      j2_ch(i+nchann_1*(nlvl-1)) =  j2_ch(1+nchann_1*(nlvl-1))
c      PRINT*,"j2",i,j2_ch(i)			  
      k2_ch(i+nchann_1*(nlvl-1)) =  k2_ch(1+nchann_1*(nlvl-1))
c      PRINT*,"k2",i,k2_ch(i)		  
      eps2_ch(i+nchann_1*(nlvl-1)) =  eps2_ch(1+nchann_1*(nlvl-1))
c      PRINT*,"p2",i,eps2_ch(i)		  
      j_t = j1_ch(i)
      ka_t = ka1_ch(i)
      kc_t = kc1_ch(i)	  
      CALL EIGEN_VAL(A1,B1,C1,j_t,ka_t,kc_t,Ech1)
      j_t = j2_ch(1+nchann_1*(nlvl-1))
      k_t =  k2_ch(1+nchann_1*(nlvl-1))		  
      Ech2 =  A2*(j_t+1d0)*j_t-
     & (A2-C2)*k_t**2
      E_ch(i+nchann_1*(nlvl-1)) = Ech2 + Ech1
c      PRINT*,i, E_ch(i)	 
      ENDDO	
      ENDDO
c      PRINT	  
      channels_defined = .TRUE.
      energy_defined = .TRUE.
      number_of_channels_defined = .TRUE.
      ENDIF	  
      ENDIF
	  
 !! PRINTINMG
       IF(emax_defined .and. .not. channels_defined) THEN
      nchann1_est = INT(2*(DSQRT(EMAX/min(min(A1,C1),B1))+1)**2)
      nchann2_est =DSQRT(EMAX/min(min(A2,C2),B2))+2
      nchann_est = nchann1_est*nchann2_est
      ALLOCATE(j1_ch_tmp(nchann1_est),
     & ka1_ch_tmp(nchann1_est),kc1_ch_tmp(nchann1_est),
     & e1_temp(nchann1_est),k2_ch_tmp(nchann2_est),
     & eps_ch_tmp(nchann2_est))
      ALLOCATE(j2_ch_tmp(nchann2_est),e2_temp(nchann2_est),
     & e_temp(nchann_est),ja1(nchann_est),SORT_STORAGE(2,nchann_est))
      i = 0 	 
      DO j_t = 0,j_max_allowed
      DO dk_t = -j_t,j_t
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"	  
      i = i+1    
      IF(i.gt.nchann1_est) EXIT
      j1_ch_tmp(i) = j_t
      ka1_ch_tmp(i) = ka_t
      kc1_ch_tmp(i) = kc_t
c      PRINT*,"j_t,ka_t,kc_t="	, j_t, ka_t,kc_t	  
      CALL EIGEN_VAL(A,B,C,j_t,ka_t,kc_t,e1_temp(i))
c      PRINT*,"e_temp,="	, e_temp(i) 
      ENDDO
      IF(i.gt.nchann1_est) EXIT	  
      ENDDO 

      n_prev = 0	  
      DO j_t=0,j_max_allowed
      DO k_t=0,j_t
      DO eps_t = 0,1-KRONEKER(k_t,0)
      IF(n_prev.gt.nchann2_est)	EXIT  
      n_prev = n_prev  + 1
      j2_ch_tmp(n_prev) =  j_t
      k2_ch_tmp(n_prev) =  k_t
      eps2_ch_tmp(n_prev) =  eps_t
      e2_temp(n_prev) = A*(j_t+1d0)*j_t - (A-C)*k_t**2
      ENDDO
      IF(n_prev.gt.nchann2_est)	EXIT 	  
      ENDDO
      IF(n_prev.gt.nchann2_est)	EXIT 	  
      ENDDO	  
	  
	  
      n_prev = 0	  
      DO i=1,nchann1_est
      DO nlvl=1,nchann2_est
      n_prev = n_prev + 1
      e_temp(n_prev) = e1_temp(i) + e2_temp(nlvl)
      SORT_STORAGE(1,n_prev) = i
      SORT_STORAGE(2,n_prev) = nlvl	  
      ENDDO
      ENDDO	

      CALL hpsort(nchann_est,e_temp,ja1)
      i = 1
      DO WHILE(e_temp(i).lt.EMAX .and. i.lt.max_nmb_chann)
       i = i+ 1	       	  
      ENDDO
      number_of_channels = i
      IF(number_of_channels.lt.0) THEN
	  PRINT*,"ERROR:NUMBER OF CHANNELS",number_of_channels
      STOP
      ENDIF	
      IF(ALLOCATED(j1_ch)) DEALLOCATE(j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE(j2_ch)	  
      IF(ALLOCATED(k1_ch)) DEALLOCATE(k1_ch)
      IF(ALLOCATED(eps1_ch)) DEALLOCATE(eps1_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE(E_ch)	  
      ALLOCATE(j1_ch(number_of_channels),
     & j2_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),E_ch(number_of_channels))	 
      DO i=1,number_of_channels
      nlvl = SORT_STORAGE(1,ja1(i))
      n_prev = SORT_STORAGE(2,ja1(i))	  
      j1_ch(i) = j1_ch_tmp(nlvl)
      ka1_ch(i) = ka1_ch_tmp(nlvl)
      kc1_ch(i) = kc1_ch_tmp(nlvl)
      j2_ch(i) = j2_ch_tmp(n_prev)
      E_ch(i) = e_temp(i)	  
      ENDDO 
      channels_defined = .TRUE.
      energy_defined = .TRUE.
      DEALLOCATE(j1_ch_tmp,ka1_ch_tmp,j2_ch_tmp,kc1_ch_tmp,e_temp)	  
      ENDIF
      jmax_included = 0
      IF(para_ortho_defined) CALL SYMMETRY_SORT_PARA_ORTHO	  
      IF(myid.eq.0) THEN
      WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"	  
      WRITE(1,'(a5,9x,a2,9x,a2,9x,a2,9x,a2,9x,a2,9x,a2,9x,a2)')"#N = "
     & ,"J1","KA1","KC1", "J2","K2","P2","E"
      ENDIF	 
      DO i=1,number_of_channels
      IF(j2_ch(i).lt.k2_ch(i))
     & STOP
     & "ERROR: J2 MUST BE NO LESS THAN K2. CHECK THE INPUT CHANNELS."
      IF(k2_ch(i).eq.0 .and. eps2_ch(i).ne.0)
     &	  STOP"ERROR:K2=0, SO PARITY MUST BE 0"
      IF(abs(eps2_ch(i)).gt.1) STOP "ERROR: PARITY MUST BE -1 OR 1"	  
      IF(j1_ch(i).lt.max(ka1_ch(i),kc1_ch(i)))
     & STOP
     & "ERROR: J1 MUST BE NO LESS THAN K1. CHECK THE INPUT CHANNELS."
      IF(.not. EXST_STATE(j1_ch(i),ka1_ch(i),kc1_ch(i))) THEN 
      WRITE(1,'(a48,1x,i4,1x,i4,1x,i4)')
     & "ERROR: CHECK CHANNELS. THIS STATE DOES NOT EXIST",
     & j1_ch(i),ka1_ch(i),kc1_ch(i)
      STOP
      ENDIF
      IF(jmax_included.lt.j1_ch(i)) jmax_included = j1_ch(i)
      IF(jmax_included.lt.j2_ch(i)) jmax_included = j2_ch(i)	  
      ENDDO	  
      DO i=1,number_of_channels
      IF(myid.eq.0) THEN	  
      WRITE(1,'(i5,i11,i11,i11,i11,i11,i11,f11.3)')i,
     & j1_ch(i),ka1_ch(i),kc1_ch(i),j2_ch(i),
     & k2_ch(i),eps2_ch(i),E_ch(i)
      ENDIF	 
      ENDDO
      energy_defined = .TRUE.	  

      IF(ini_chann_defined .and. .not.all_level_included) THEN
      chann_ini = 0	  
      DO i=1,number_of_channels
      IF(j1_ini.eq.j1_ch(i) .and. ka1_ini.eq.ka1_ch(i) .and. 
     & kc1_ini.eq.kc1_ch(i).and.j2_ini.eq.j2_ch(i).and.
     & k2_ini.eq.k2_ch(i).and.eps2_ini.eq.eps2_ch(i))
     & chann_ini = i
      ENDDO	
      IF(chann_ini.eq.0) STOP "ERROR: INITIAL CHANNEL WRONG"
      IF(myid.eq.0) 	  
     & WRITE(1,'(a17,a5,i4,a8,i4,a8,i4,a7,i4,a7,i4,a7,i4)')
     & "INITIAL STATE IS:",
     & "J1 = ",j1_ini,"; KA1 = ", ka1_ini,"; KC1 =	",kc1_ini,
     & "; J2 =	", j2_ini, "; K2 = ", k2_ini, "; P2 = ", eps2_ini
      ENDIF
      IF(.not.ini_chann_defined .and. channels_defined .and.
     & energy_defined) THEN
	  i_ground = 1
      DO i=1,number_of_channels
      IF(E_ch(i_ground).gt.E_ch(i)) i_ground = i   
      ENDDO	 
      chann_ini = i_ground
      ini_chann_defined = .TRUE.	  
      ENDIF	  
      IF(myid.eq.0) THEN  
      WRITE(1,*)
      IF(angs_unit)WRITE(1,'(a35)')"THE POTENTIAL R-UNITS ARE ANGSTROMS"
      IF(au_unit_r)WRITE(1,'(a31)') "THE POTENTIAL R-UNITS ARE BOHRS"
      IF(cm_unit)WRITE(1,'(a31)')"THE POTENTIAL ENERGY IS IN CM-1"
      IF(au_unit_e)WRITE(1,'(a34)')"THE POTENTIAL ENERGY IS IN HARTREE"
      IF(klv_unit)WRITE(1,'(a33)')"THE POTENTIAL ENERGY IS IN KELVIN"
      IF(kal_unit)WRITE(1,'(a35)')"THE POTENTIAL ENERGY IS IN KCAL/MOL"
      ENDIF	  
      IF(.not.angs_unit .and..not.au_unit_r )   
     & STOP "ERROR:R_UNITS NOT DEFINED"
      IF(.not.cm_unit .and..not.au_unit_e.and..not.klv_unit
     & .and. .not.kal_unit)   
     & STOP "ERROR:ENERGY NOT DEFINED"
      IF(myid.eq.0) WRITE(1,*) 	 
      IF(expansion_defined .and. myid.eq.0) THEN
      WRITE(1,'(a42)') "USER SUPPLIED POTENTIAL TERMS WILL BE USED"
      WRITE(1,*)
      WRITE(1,'(a45,i4)')
     & "NUMBER OF POTENTIAL TERMS INCLUDED N_TERMS = ",nterms
      WRITE(1,'(a45)')
     & "THE POTENTIAL TERMS ARE, (L1,NJU1,L2,NJU2,L):"
      DO i=1,nterms
      WRITE(1,'(a1,i3,a1,i3,a1,i3,a1,i3,a1,i3,a3)',ADVANCE ="NO")
     & '(',L1_TERM(i),',',NJU1_TERM(i),',',L2_TERM(i),
     & ',',NJU2_TERM(i),',',L_TERM(i),'); '	  
      ENDDO		
      WRITE(1,*)
      WRITE(1,*)	  
      WRITE(1,'(a33,i3)') "NUMBER OF R-GRID POINTS IS N_R = ",
     & n_r_coll	  
      ELSE
      IF(.not.matrix_reading_defined .and. myid.eq.0)	  
     & WRITE(1,'(a97)')"COULPING MATRIX WILL BE COMPUTED USING NUMERICAL
     &  INTEGRATION OF POTENTIAL OVER THE WAVEFUNCTIONS"	  
      ENDIF	  
      IF(myid.eq.0) WRITE(1,*)  
      CASE(0)
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a23,1x,a31)')"THE COLLISIONAL SYSTEM:",
     & "ASYMMETRIC TOP + ASYMMETRIC TOP"
      ENDIF	 
      IF(identical_particles_defined) THEN
      IF(MYID.eq.0) THEN	  
      WRITE(1,*)	  
   	  WRITE(1,'(a57)')
     & "THE COLLISIONAL PARTNERS ARE DEFINED AS INDINSTIGUISHABLE"
      WRITE(1,*)
      ENDIF	  
      ENDIF	 
      IF(user_defined) THEN
      IF(MYID.eq.0)	  
     & WRITE(1,'(a7,2x,a60)') "WARNING:",
     & "THE WAVEFUNCTIONS AND ENERGY LEVELS MUST BE SUPPLIED BY USER"
      	 
      IF(rot_const_defined)WRITE(1,'(a27)')"ROT. CONST WILL BE INGNORED"
      INQUIRE( FILE=USR_INP_LEVEL, EXIST=exst )
      IF(.not.exst)	THEN
      IF(MYID.eq.0)PRINT*,USR_INP_LEVEL, " NOT FOUND"
      STOP	  
      ENDIF	  
      OPEN(325,FILE=USR_INP_LEVEL,
     & STATUS="OLD",ACTION="READ")
      IF(number_of_channels.le.0) THEN
      PRINT*,"NUMBER OF CHANNELS NOT SPECIFIED"
      STOP	  
      ENDIF
      IF(.not.channels_defined) THEN	  
      ALLOCATE(j1_ch(number_of_channels))
      ALLOCATE(ka1_ch(number_of_channels))
      ALLOCATE(kc1_ch(number_of_channels))
      ALLOCATE(j2_ch(number_of_channels))
      ALLOCATE(ka2_ch(number_of_channels))
      ALLOCATE(kc2_ch(number_of_channels))
	  
      ENDIF
      IF(.not.energy_defined) THEN	  
      ALLOCATE(E_ch(number_of_channels))
      ENDIF
      READ(325,*) !HEADER
      jmax_included = 0	  
      DO i=1,number_of_channels
      READ(325,*,IOSTAT = stat_of_file)
     & decr,j1_ch(i), ka1_ch(i),kc1_ch(i),
     & j2_ch(i),ka2_ch(i),kc2_ch(i),E_ch(i)
      IF(jmax_included.lt.j1_ch(i))jmax_included=j1_ch(i)
      IF(jmax_included.lt.j2_ch(i))jmax_included=j2_ch(i)	 
	  
      IF(stat_of_file.gt.0 .or. decr.ne.i) THEN
      PRINT*,"WRONG CHANNELS IN USER INPUT FILE"	  
      STOP
      ENDIF 	  
      ENDDO
      READ(325,*) !! HEADER FOR BASIS
      ALLOCATE(M_EIGEN_VECTORS_USER
     & (number_of_channels,jmax_included*4+2))
      M_EIGEN_VECTORS_USER = 0d0	 
      DO i=1,number_of_channels
      READ(325,*) decr,M_EIGEN_VECTORS_USER
     & (i,1:2+2*(j1_ch(i)+j2_ch(i)))
      IF(i.ne.decr .or. stat_of_file.gt.0) STOP "ERROR IN USER FILE"	 
      ENDDO	  	  
      CLOSE(325)
      channels_defined = .TRUE.
      energy_defined = .TRUE.      
      ELSE
      IF(identical_particles_defined) THEN
      A2 = A1
      B2 = B1
      C2 = C1
  
      ENDIF	 	  
      IF(A1*B1*C1.le.0d0)
     & STOP"ERROR: ALL A1,B1 AND C1 MUST BE POSITIVE. CHECK INPUT"
      IF(MYID.eq.0) THEN      
      WRITE(1,'(a51,f7.2,a8,f7.2,a6,f8.2,a2)')
     & "ROTATIONAL CONSTANTS OF #1 MOLECULE ARE, CM-1: A1 ="
     & ,A1," ; B1 = ",
     &  B1,"; C1 = ", C1,";"
      ENDIF	 
      IF(A2*B2*C2.le.0d0)
     & STOP"ERROR: ALL A2,B2 AND C2 MUST BE POSITIVE. CHECK INPUT"
      IF(MYID.eq.0) THEN
      WRITE(1,'(a51,f7.2,a8,f7.2,a6,f8.2,a2)')
     & "ROTATIONAL CONSTANTS OF #2 MOLECULE ARE, CM-1: A2 ="
     & ,A2," ; B2 = ",
     &  B2,"; C2 = ", C2,";"
      ENDIF
       IF(emax_defined .and. .not. channels_defined) THEN
      nchann1_est = INT(2*(DSQRT(EMAX/min(min(A1,C1),B1))+1)**2)
      nchann2_est = INT(2*(DSQRT(EMAX/min(min(A2,C2),B2))+1)**2)
      nchann_est = nchann1_est*nchann2_est
      ALLOCATE(j1_ch_tmp(nchann1_est),
     & ka1_ch_tmp(nchann1_est),kc1_ch_tmp(nchann1_est),
     & e1_temp(nchann1_est),k2_ch_tmp(nchann2_est),
     & eps_ch_tmp(nchann2_est))
      ALLOCATE(j2_ch_tmp(nchann2_est),e2_temp(nchann2_est),
     & e_temp(nchann_est),ja1(nchann_est),SORT_STORAGE(2,nchann_est))
      i = 0 	 
      DO j_t = 0,j_max_allowed
      DO dk_t = -j_t,j_t
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"	  
      i = i+1    
      IF(i.gt.nchann1_est) EXIT
      j1_ch_tmp(i) = j_t
      ka1_ch_tmp(i) = ka_t
      kc1_ch_tmp(i) = kc_t
c      PRINT*,"j_t,ka_t,kc_t="	, j_t, ka_t,kc_t	  
      CALL EIGEN_VAL(A,B,C,j_t,ka_t,kc_t,e1_temp(i))
c      PRINT*,"e_temp,="	, e_temp(i) 
      ENDDO
      IF(i.gt.nchann1_est) EXIT	  
      ENDDO 

      n_prev = 0	  
      DO j_t=0,j_max_allowed
      DO k_t=0,j_t
      DO eps_t = 0,1-KRONEKER(k_t,0)
      IF(n_prev.gt.nchann2_est)	EXIT  
      n_prev = n_prev  + 1
      j2_ch_tmp(n_prev) =  j_t
      k2_ch_tmp(n_prev) =  k_t
      eps2_ch_tmp(n_prev) =  eps_t
      e2_temp(n_prev) = A*(j_t+1d0)*j_t - (A-C)*k_t**2
      ENDDO
      IF(n_prev.gt.nchann2_est)	EXIT 	  
      ENDDO
      IF(n_prev.gt.nchann2_est)	EXIT 	  
      ENDDO	  
	  
	  
      n_prev = 0	  
      DO i=1,nchann1_est
      DO nlvl=1,nchann2_est
      n_prev = n_prev + 1
      e_temp(n_prev) = e1_temp(i) + e2_temp(nlvl)
      SORT_STORAGE(1,n_prev) = i
      SORT_STORAGE(2,n_prev) = nlvl	  
      ENDDO
      ENDDO	

      CALL hpsort(nchann_est,e_temp,ja1)
      i = 1
      DO WHILE(e_temp(i).lt.EMAX .and. i.lt.max_nmb_chann)
       i = i+ 1	       	  
      ENDDO
      number_of_channels = i
      IF(number_of_channels.lt.0) THEN
	  PRINT*,"ERROR:NUMBER OF CHANNELS",number_of_channels
      STOP
      ENDIF	
      IF(ALLOCATED(j1_ch)) DEALLOCATE(j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE(j2_ch)	  
      IF(ALLOCATED(k1_ch)) DEALLOCATE(k1_ch)
      IF(ALLOCATED(eps1_ch)) DEALLOCATE(eps1_ch)
      IF(ALLOCATED(E_ch)) DEALLOCATE(E_ch)	  
      ALLOCATE(j1_ch(number_of_channels),
     & j2_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),E_ch(number_of_channels))	 
      DO i=1,number_of_channels
      nlvl = SORT_STORAGE(1,ja1(i))
      n_prev = SORT_STORAGE(2,ja1(i))	  
      j1_ch(i) = j1_ch_tmp(nlvl)
      ka1_ch(i) = ka1_ch_tmp(nlvl)
      kc1_ch(i) = kc1_ch_tmp(nlvl)
      j2_ch(i) = j2_ch_tmp(n_prev)
      E_ch(i) = e_temp(i)	  
      ENDDO 
      channels_defined = .TRUE.
      energy_defined = .TRUE.
      DEALLOCATE(j1_ch_tmp,ka1_ch_tmp,j2_ch_tmp,kc1_ch_tmp,e_temp)	  
      ENDIF

      IF(mlc_mlc_emax_defined.and. .not.channels_defined) THEN
      IF(EMAX1.lt.0d0 .or. EMAX2.lt.0)
     & STOP"ERROR:EMAX1,2 MUST BE POSITIVE"	  
c      PRINT*,"EMAX=",EMAX1,EMAX2	  
      j_t = int(abs(min(sqrt(EMAX1/C1),(EMAX1/A1-1)/2d0)))
      nchnl_ini_guess = (j_t+1)**2	  
c      PRINT*,"j_ini",j_t	  
      nchann_est = nchnl_ini_guess
      n_prev = -1
      nlvl = 0
      decr = 1
      DO WHILE(nlvl.ne.n_prev)
      n_prev =  nlvl	  
      nchann_est = nchnl_ini_guess*decr
c      PRINT*,"nchann_est = ",nchann_est	  
      ALLOCATE(j_ch_tmp(nchann_est),ka_ch_tmp(nchann_est),
     & kc_ch_tmp(nchann_est),e_temp(nchann_est),ja(nchann_est))
      nlvl = 0
      DO j_t=0,j_max_allowed
      DO dk_t=-j_t,j_t
      nlvl = nlvl + 1
      IF(nlvl.gt.nchann_est) EXIT	  
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"
      j_ch_tmp(nlvl) = j_t
      ka_ch_tmp(nlvl) = ka_t
      kc_ch_tmp(nlvl) = kc_t
      CALL EIGEN_VAL(A1,B1,C1,j_t,ka_t,kc_t,e_temp(nlvl))
      ENDDO
      IF(nlvl.gt.nchann_est) EXIT		  
      ENDDO
      CALL	hpsort(nchann_est,e_temp,ja)
      i = 1	  
      DO WHILE(e_temp(i).lt.EMAX1)
      i = i + 1
c      PRINT*,"E",i,e_temp(i)	  
      ENDDO
      nlvl = i - 1	  
c      PRINT*,'nlvl',nlvl      	  
      decr = decr + 1 
      DEALLOCATE(j_ch_tmp,ka_ch_tmp,
     & kc_ch_tmp,e_temp,ja)	  
      ENDDO
      nchann_1 = nlvl
      buffer_ch1 = nchann_est
c      PRINT *,"nchann_1",nchann_1,"buffer_ch1",buffer_ch1 	  
!!!!! FOR SYMMETRIC TOP	  
      j_t = int(abs(min(sqrt(EMAX2/C2),(EMAX2/A2-1)/2d0)))
      nchnl_ini_guess = (j_t+1)**2	  
c      PRINT*,"j_ini",j_t	  
      nchann_est = nchnl_ini_guess
      n_prev = -1
      nlvl = 0
      decr = 1
      DO WHILE(nlvl.ne.n_prev)
      n_prev =  nlvl	  
      nchann_est = nchnl_ini_guess*decr
c      PRINT*,"nchann_est = ",nchann_est	  
      ALLOCATE(j_ch_tmp(nchann_est),ka_ch_tmp(nchann_est),
     & kc_ch_tmp(nchann_est),e_temp(nchann_est),ja(nchann_est))
      nlvl = 0
      DO j_t=0,j_max_allowed
      DO dk_t=-j_t,j_t
      nlvl = nlvl + 1
      IF(nlvl.gt.nchann_est) EXIT	  
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"
      j_ch_tmp(nlvl) = j_t
      ka_ch_tmp(nlvl) = ka_t
      kc_ch_tmp(nlvl) = kc_t
      CALL EIGEN_VAL(A2,B2,C2,j_t,ka_t,kc_t,e_temp(nlvl))
      ENDDO
      IF(nlvl.gt.nchann_est) EXIT		  
      ENDDO
      CALL	hpsort(nchann_est,e_temp,ja)
      i = 1	  
      DO WHILE(e_temp(i).lt.EMAX2)
      i = i + 1
c      PRINT*,"E",i,e_temp(i)	  
      ENDDO
      nlvl = i - 1	  
c      PRINT*,'nlvl',nlvl      	  
      decr = decr + 1 
      DEALLOCATE(j_ch_tmp,ka_ch_tmp,
     & kc_ch_tmp,e_temp,ja)	  
      ENDDO
      nchann_2 = nlvl
      buffer_ch2 = nchann_est
c      PRINT *,"nchann_1",nchann_1,"buffer_ch1",buffer_ch1 	  
!!!!! FOR SYMMETRIC TOP	  
c      PRINT *,"nchann_2",nchann_2,"buffer_ch2",buffer_ch2  
      number_of_channels = nchann_1*nchann_2
	  !end
c      PRINT*,number_of_channels,nchann_2,buffer_ch1
      ALLOCATE(j1_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),E_ch(number_of_channels),
     & j2_ch(number_of_channels),ka2_ch(number_of_channels),
     & kc2_ch(number_of_channels))	 
       ALLOCATE(j_ch_tmp(buffer_ch1),ka_ch_tmp(buffer_ch1),
     & kc_ch_tmp(buffer_ch1),e_temp(buffer_ch1),ja(buffer_ch1))
	 
      nlvl = 0
      DO j_t=0,j_max_allowed
      DO dk_t=-j_t,j_t
      nlvl = nlvl + 1
      IF(nlvl.gt.buffer_ch1) EXIT	  
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst)
     &	  STOP"ERROR:ASSYM. TOP INITALIZATION FAILED FOR EMAX_DEFINED"
      j_ch_tmp(nlvl) = j_t
      ka_ch_tmp(nlvl) = ka_t
      kc_ch_tmp(nlvl) = kc_t
      CALL EIGEN_VAL(A1,B1,C1,j_t,ka_t,kc_t,e_temp(nlvl))	  
      ENDDO
      IF(nlvl.gt.buffer_ch1) EXIT		  
      ENDDO
c      PRINT*,"number_of_channels_again = ",number_of_channels     
      CALL	hpsort(buffer_ch1,e_temp,ja)	  
      DO i=1,nchann_1
      j1_ch(i) = j_ch_tmp(ja(i))
      ka1_ch(i) = ka_ch_tmp(ja(i))
      kc1_ch(i) = kc_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      DEALLOCATE(j_ch_tmp,ka_ch_tmp,
     & kc_ch_tmp,e_temp,ja)
      energy_defined = .TRUE.
c      PRINT*,"j1_ch",j1_ch	  
!!!!!!!!!!!!! for assymetric top
       ALLOCATE(j_ch_tmp(buffer_ch2),ka_ch_tmp(buffer_ch2),
     & kc_ch_tmp(buffer_ch2),e_temp(buffer_ch2),ja(buffer_ch2))
      nlvl = 0
      DO j_t=0,j_max_allowed
      DO dk_t=-j_t,j_t
      nlvl = nlvl + 1
      IF(nlvl.gt.buffer_ch2) EXIT	  
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst)
     &	  STOP"ERROR:ASSYM. TOP INITALIZATION FAILED FOR EMAX_DEFINED"
      j_ch_tmp(nlvl) = j_t
      ka_ch_tmp(nlvl) = ka_t
      kc_ch_tmp(nlvl) = kc_t
      CALL EIGEN_VAL(A2,B2,C2,j_t,ka_t,kc_t,e_temp(nlvl))	  
      ENDDO
      IF(nlvl.gt.buffer_ch2) EXIT		  
      ENDDO
c      PRINT*,"number_of_channels_again = ",number_of_channels     
      CALL	hpsort(buffer_ch2,e_temp,ja)	  
      DO i=1,nchann_2
      j2_ch(1+nchann_1*(i-1)) = j_ch_tmp(ja(i))
      ka2_ch(1+nchann_1*(i-1)) = ka_ch_tmp(ja(i))
      kc2_ch(1+nchann_1*(i-1)) = kc_ch_tmp(ja(i))
      E_ch(1+nchann_1*(i-1)) = e_temp(i)	  
      ENDDO	 
      DEALLOCATE(j_ch_tmp,ka_ch_tmp,
     & kc_ch_tmp,e_temp,ja)
c      PRINT *,"j2=",j2_ch  
      DO nlvl = 1,nchann_2
      DO i=1,nchann_1
      j1_ch(i+nchann_1*(nlvl-1)) = j1_ch(i)
c      PRINT*,"j1",i,j1_ch(i)	  
      ka1_ch(i+nchann_1*(nlvl-1)) = ka1_ch(i)
c      PRINT*,"ka1",i,ka1_ch(i)	  
      kc1_ch(i+nchann_1*(nlvl-1)) = kc1_ch(i)
c      PRINT*,"kc1",i,kc1_ch(i)		  
      j2_ch(i+nchann_1*(nlvl-1)) =  j2_ch(1+nchann_1*(nlvl-1))
c      PRINT*,"j2",i,j2_ch(i)			  
      ka2_ch(i+nchann_1*(nlvl-1)) =  ka2_ch(1+nchann_1*(nlvl-1))
c      PRINT*,"k2",i,k2_ch(i)		  
      kc2_ch(i+nchann_1*(nlvl-1)) =  kc2_ch(1+nchann_1*(nlvl-1))
c      PRINT*,"p2",i,eps2_ch(i)		  
      j_t = j1_ch(i)
      ka_t = ka1_ch(i)
      kc_t = kc1_ch(i)	  
      CALL EIGEN_VAL(A1,B1,C1,j_t,ka_t,kc_t,Ech1)
      j_t = j2_ch(i+nchann_1*(nlvl-1))
      ka_t = ka2_ch(i+nchann_1*(nlvl-1))
      kc_t = kc2_ch(i+nchann_1*(nlvl-1))	  
      CALL EIGEN_VAL(A2,B2,C2,j_t,ka_t,kc_t,Ech2)
      E_ch(i+nchann_1*(nlvl-1)) = Ech2 + Ech1
c      PRINT*,i, E_ch(i)	 
      ENDDO	
      ENDDO
c      PRINT	  
      channels_defined = .TRUE.
      energy_defined = .TRUE.
      number_of_channels_defined = .TRUE.
      ENDIF		  
      IF(number_of_channels_defined) THEN
      IF(energy_defined.and. .not.channels_defined) STOP 
     & "ERROR:FOR MANUALLY DEFINED ENERGY YOU SHOULD SPECIFY CHANNELS"
      IF(.not.energy_defined) ALLOCATE(E_ch(number_of_channels))
      IF(.not.channels_defined) THEN 
      decr = 1	 
      nlvl = number_of_channels*decr
      ALLOCATE(j1_ch_tmp(nlvl),
     & ka1_ch_tmp(nlvl),ja1(nlvl**2),
     & kc1_ch_tmp(nlvl),e1_temp(nlvl))
      ALLOCATE(j2_ch_tmp(nlvl),
     & ka2_ch_tmp(nlvl),
     & kc2_ch_tmp(nlvl),e2_temp(nlvl))	  	 
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
      ALLOCATE(j_ch_tmp(nlvl),
     & ka_ch_tmp(nlvl),ja(nlvl),
     & kc_ch_tmp(nlvl),e_temp(nlvl))
      i = 0	 
      DO j_t=0,j_max_allowed
      DO dk_t=-j_t,j_t
      i = i + 1
      IF(i.gt.nlvl) EXIT	  
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"
      j_ch_tmp(i) = j_t
      ka_ch_tmp(i) = ka_t
      kc_ch_tmp(i) = kc_t
      CALL EIGEN_VAL(A1,B1,C1,j_t,ka_t,kc_t,e_temp(i))
      ENDDO
      IF(i.gt.nlvl) EXIT		  
      ENDDO
!      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
!      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,number_of_channels	  
      IF(e1_temp(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,number_of_channels
      j1_ch_tmp(i) = j_ch_tmp(ja(i))
      ka1_ch_tmp(i) = ka_ch_tmp(ja(i))
      kc1_ch_tmp(i) = kc_ch_tmp(ja(i))
      e1_temp(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = number_of_channels*decr      	  
      DEALLOCATE(ka_ch_tmp,j_ch_tmp,ja,e_temp,kc_ch_tmp)	  
      ENDDO
!      IF(myid.eq.5) PRINT*,"CHANNELS DEFINED",e1_temp
      decr = 1	 
      nlvl = number_of_channels*decr
      channels_defined	= .FALSE.  
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
      ALLOCATE(j_ch_tmp(nlvl),
     & ka_ch_tmp(nlvl),ja(nlvl),
     & kc_ch_tmp(nlvl),e_temp(nlvl))
      i = 0	 
      DO j_t=0,j_max_allowed
      DO dk_t=-j_t,j_t
      i = i + 1
      IF(i.gt.nlvl) EXIT	  
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"
      j_ch_tmp(i) = j_t
      ka_ch_tmp(i) = ka_t
      kc_ch_tmp(i) = kc_t
      CALL EIGEN_VAL(A2,B2,C2,j_t,ka_t,kc_t,e_temp(i))
      ENDDO
      IF(i.gt.nlvl) EXIT		  
      ENDDO
!      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
!      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,number_of_channels	  
      IF(e2_temp(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,number_of_channels
      j2_ch_tmp(i) = j_ch_tmp(ja(i))
      ka2_ch_tmp(i) = ka_ch_tmp(ja(i))
      kc2_ch_tmp(i) = kc_ch_tmp(ja(i))
      e2_temp(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = number_of_channels*decr      	  
      DEALLOCATE(ka_ch_tmp,j_ch_tmp,ja,e_temp,kc_ch_tmp)	  
      ENDDO 
!      IF(myid.eq.5) PRINT*,myid,"CHANNELS DEFINED 2",e2_temp	  
      ALLOCATE(SORT_STORAGE(2,number_of_channels**2))
      ALLOCATE(e_temp(number_of_channels**2))
      SORT_STORAGE = 0
!      IF(myid.eq.5)PRINT*,myid," SORT_STORAGE",SORT_STORAGE		  
      n_prev = 0 	  
      DO i=1,number_of_channels
      DO nlvl=1,number_of_channels
      n_prev = n_prev +1
!      IF(myid.eq.4) PRINT*,myid,i,nlvl,n_prev 	  
      SORT_STORAGE(1,n_prev)  = i 
      SORT_STORAGE(2,n_prev)  = nlvl    	  
      e_temp(n_prev) = e1_temp(i) + e2_temp(nlvl)
!      IF(myid.eq.9) PRINT*, e_temp(n_prev),e1_temp(i),e2_temp(nlvl)	  
      ENDDO
	  
      ENDDO
!      IF(myid.eq.1)PRINT*,SORT_STORAGE	  
!      STOP	  
      CALL hpsort(number_of_channels**2,e_temp,ja1)
!      IF(myid.eq.0) PRINT*,"ENERYG SORTED DEFINED",e_temp		  
      IF(ALLOCATED(j1_ch)) DEALLOCATE(j1_ch)
      IF(ALLOCATED(j2_ch)) DEALLOCATE(j2_ch)	  
      IF(ALLOCATED(ka1_ch)) DEALLOCATE(ka1_ch)
      IF(ALLOCATED(kc1_ch)) DEALLOCATE(kc1_ch)
      IF(ALLOCATED(ka2_ch)) DEALLOCATE(ka2_ch)
      IF(ALLOCATED(kc2_ch)) DEALLOCATE(kc2_ch)	  
      IF(ALLOCATED(E_ch))  DEALLOCATE(E_ch)


!      IF(myid.eq.0) PRINT*,"JS SORTED DEFINED",ja1		  
      ALLOCATE(j1_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),j2_ch(number_of_channels),
     & ka2_ch(number_of_channels),kc2_ch(number_of_channels),
     & E_ch(number_of_channels) )
!      IF(myid.eq.0) PRINT*,"ARRAYS CREATED"	
!      STOP	  
      DO i=1,number_of_channels
      nlvl = SORT_STORAGE(1,ja1(i))
      n_prev = SORT_STORAGE(2,ja1(i))	  
      j1_ch(i) = j1_ch_tmp(nlvl)
      ka1_ch(i) = ka1_ch_tmp(nlvl)
      kc1_ch(i) = kc1_ch_tmp(nlvl)
      j2_ch(i) = j2_ch_tmp(n_prev)
      ka2_ch(i) = ka2_ch_tmp(n_prev)
      kc2_ch(i) = kc2_ch_tmp(n_prev)	  
      E_ch(i) = e_temp(i)	  
      ENDDO
      DEALLOCATE(SORT_STORAGE,ja1,j1_ch_tmp,ka1_ch_tmp,
     & kc1_ch_tmp,j2_ch_tmp,e_temp,ka2_ch_tmp,kc2_ch_tmp)	  
      channels_defined = .TRUE.
      energy_defined = .TRUE.
!      IF(MYID.eq.0) PRINT*,"ALL DEFINED"	  
      ENDIF
      IF(channels_defined .and. .not.energy_defined ) THEN
      DO i=1,number_of_channels
      j_t = j1_ch(i)
      ka_t = ka1_ch(i)
      kc_t = kc1_ch(i)	  
      CALL EIGEN_VAL(A1,B1,C1,j_t,ka_t,kc_t,Ech1)
      j_t = j2_ch(i)
      ka_t = ka2_ch(i)
      kc_t = kc2_ch(i)		  
      CALL EIGEN_VAL(A2,B2,C2,j_t,ka_t,kc_t,Ech2)
      E_ch(i) = Ech1 + Ech2  
      ENDDO
	  
	  
      energy_defined = .TRUE.     	  
      ENDIF
!!!!! PRINTING
      ENDIF
      IF(mlc_mlc_chn_num_defined.and.mlc_mlc_emax_defined) 
     & STOP "ERROR:EITHER EMAX1,2 OR NCHLS1,2 MUST BE SPECIFIED"
      IF(mlc_mlc_chn_num_defined) THEN
      number_of_channels = nchann_1*nchann_2
      ALLOCATE(j1_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),j2_ch(number_of_channels),
     & ka2_ch(number_of_channels),
     & kc2_ch(number_of_channels),E_ch(number_of_channels))
      j1_ch = 0
      ka1_ch = 0
      kc1_ch = 0
      j2_ch = 0
      ka2_ch = 0
      kc2_ch = 0	
      decr = 1
      IF(max_nmb_chann.lt.number_of_channels) STOP "TOO MANY CHANNELS"	  
      nlvl = nchann_1*decr	  
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
      i = 0	 
      DO j_t = 0,j_max_allowed
      DO dk_t = -j_t,j_t
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"	  
      i = i+1    
      IF(i.gt.nlvl) EXIT
      j_ch_tmp(i) = j_t
      ka_ch_tmp(i) = ka_t
      kc_ch_tmp(i) = kc_t
c      PRINT*,"j_t,ka_t,kc_t="	, j_t, ka_t,kc_t	  
      CALL EIGEN_VAL(A1,B1,C1,j_t,ka_t,kc_t,e_temp(i))
c      PRINT*,"e_temp,="	, e_temp(i) 
      ENDDO
      IF(i.gt.nlvl) EXIT	  
      ENDDO
c      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
c      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,nchann_1	  
      IF(E_ch(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,nchann_1
      j1_ch(i) = j_ch_tmp(ja(i))
      ka1_ch(i) = ka_ch_tmp(ja(i))
      kc1_ch(i) = kc_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = nchann_1*decr      	  
      DEALLOCATE(ka_ch_tmp,j_ch_tmp,ja,e_temp,kc_ch_tmp)	  
      ENDDO 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      decr = 1
      IF(max_nmb_chann.lt.number_of_channels) STOP "TOO MANY CHANNELS"	  
      nlvl = nchann_2*decr	  
      DO WHILE(nlvl.lt.max_nmb_chann .and. .not.channels_defined)
      i = 0	 
      DO j_t = 0,j_max_allowed
      DO dk_t = -j_t,j_t
      CALL STATE(j_t,dk_t,ka_t,kc_t,exst)
      IF(.not.exst) STOP"ERROR:ASSYM. TOP INITALIZATION FAILED"	  
      i = i+1    
      IF(i.gt.nlvl) EXIT
      j_ch_tmp(i) = j_t
      ka_ch_tmp(i) = ka_t
      kc_ch_tmp(i) = kc_t
c      PRINT*,"j_t,ka_t,kc_t="	, j_t, ka_t,kc_t	  
      CALL EIGEN_VAL(A2,B2,C2,j_t,ka_t,kc_t,e_temp(i))
c      PRINT*,"e_temp,="	, e_temp(i) 
      ENDDO
      IF(i.gt.nlvl) EXIT	  
      ENDDO
c      PRINT *, "ENERGY LEVELS",e_temp
      CALL	hpsort(nlvl,e_temp,ja)
c      PRINT *, "ja hpsort",ja	  
      IF(decr.gt.1) THEN
      channels_defined = .TRUE.
      DO i=1,nchann_2	  
      IF(E_ch(i).ne.e_temp(i)) THEN
      channels_defined = .FALSE.  	 
      ENDIF
      ENDDO
      ENDIF	  
      DO i=1,nchann_2
      j2_ch(i) = j_ch_tmp(ja(i))
      ka2_ch(i) = ka_ch_tmp(ja(i))
      kc2_ch(i) = kc_ch_tmp(ja(i))
      E_ch(i) = e_temp(i)	  
      ENDDO	 
      decr = decr + 1
      nlvl = nchann_2*decr      	  
      DEALLOCATE(ka_ch_tmp,j_ch_tmp,ja,e_temp,kc_ch_tmp)	  
      ENDDO 
	  
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO nlvl = 1,nchann_2
      DO i=1,nchann_1
      j1_ch(i+nchann_1*(nlvl-1)) = j1_ch(i)
c      PRINT*,"j1",i,j1_ch(i)	  
      ka1_ch(i+nchann_1*(nlvl-1)) = ka1_ch(i)
c      PRINT*,"ka1",i,ka1_ch(i)	  
      kc1_ch(i+nchann_1*(nlvl-1)) = kc1_ch(i)
c      PRINT*,"kc1",i,kc1_ch(i)		  
      j2_ch(i+nchann_1*(nlvl-1)) =  j2_ch(1+nchann_1*(nlvl-1))
c      PRINT*,"j2",i,j2_ch(i)			  
      ka2_ch(i+nchann_1*(nlvl-1)) =  ka2_ch(1+nchann_1*(nlvl-1))
c      PRINT*,"k2",i,k2_ch(i)		  
      kc2_ch(i+nchann_1*(nlvl-1)) =  kc2_ch(1+nchann_1*(nlvl-1))
c      PRINT*,"p2",i,eps2_ch(i)		  
      j_t = j1_ch(i)
      ka_t = ka1_ch(i)
      kc_t = kc1_ch(i)	  
      CALL EIGEN_VAL(A1,B1,C1,j_t,ka_t,kc_t,Ech1)
      j_t = j2_ch(i)
      ka_t = ka2_ch(i)
      kc_t = kc2_ch(i)	  
      CALL EIGEN_VAL(A2,B2,C2,j_t,ka_t,kc_t,Ech2)
      E_ch(i+nchann_1*(nlvl-1)) = Ech2 + Ech1
c      PRINT*,i, E_ch(i)	 
      ENDDO	
      ENDDO
      number_of_channels_defined = .TRUE.
      energy_defined = .TRUE.	  
      ENDIF	  	  
      ENDIF
  
 !! PRINTINMG
      IF(.not.user_defined) THEN
      jmax_included = 0
      IF(para_ortho_defined) CALL SYMMETRY_SORT_PARA_ORTHO	  
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a24)') "CHANNELS ENERGIES, CM-1:"
      if(.not.identical_particles_defined) then
	  WRITE(1,'(a5,9x,a2,9x,a3,9x,a3,9x,a2,9x,a3,9x,a3,9x,a2)')"#N = "
     & ,"J1","KA1","KC1", "J2","KA2","KC2","E"
	  else
	  WRITE(1,'(a5,9x,a2,9x,a3,9x,a3,9x,a2,9x,a3,9x,a3,
     & a9,a9,a9,a9,9x,a2)')
     & "#N = ","J1","KA1","KC1", "J2","KA2","KC2",
     & "p1","p2","kappa1","kappa2","E"
	  end if
      ENDIF      
      DO i=1,number_of_channels
      IF(j2_ch(i).lt.max(ka2_ch(i),kc2_ch(i))) THEN
      PRINT *,j2_ch(i),ka2_ch(i),kc2_ch(i)	  
      STOP
     & "ERROR: J2 MUST BE NO LESS THAN K2. CHECK THE INPUT CHANNELS."	  
      ENDIF	
      IF(.not. EXST_STATE(j2_ch(i),ka2_ch(i),kc2_ch(i))) THEN 
      IF(MYID.eq.0) WRITE(1,'(a48,1x,i4,1x,i4,1x,i4)')
     & "ERROR: CHECK CHANNELS. THIS STATE DOES NOT EXIST",
     & j2_ch(i),ka2_ch(i),kc2_ch(i)
      STOP
      ENDIF       
      IF(j1_ch(i).lt.max(ka1_ch(i),kc1_ch(i)))
     & STOP
     & "ERROR: J1 MUST BE NO LESS THAN K1. CHECK THE INPUT CHANNELS."
      IF(.not. EXST_STATE(j1_ch(i),ka1_ch(i),kc1_ch(i))) THEN 
      IF(MYID.eq.0) WRITE(1,'(a48,1x,i4,1x,i4,1x,i4)')
     & "ERROR: CHECK CHANNELS. THIS STATE DOES NOT EXIST",
     & j1_ch(i),ka1_ch(i),kc1_ch(i)
      STOP
      ENDIF
      IF(jmax_included.lt.j1_ch(i)) jmax_included = j1_ch(i)
      IF(jmax_included.lt.j2_ch(i)) jmax_included = j2_ch(i)	  
      ENDDO	  
      DO i=1,number_of_channels
      IF(MYID.eq.0) THEN
      if(.not.identical_particles_defined) then
      WRITE(1,'(i5,i11,i12,i12,i11,i12,i12,f11.3)')i,
     & j1_ch(i),ka1_ch(i),kc1_ch(i),j2_ch(i),
     & ka2_ch(i),kc2_ch(i),E_ch(i)
	  else
	  call ASYM_TOP_VECTORS
	  call bikram_kappa(ka1_ch(i),kc1_ch(i),ka2_ch(i),kc2_ch(i),
     & kappa1, kappa2)
      WRITE(1,'(i5,i11,i12,i12,i11,i12,i12,i9,i9,i9,i9,f11.3)')i,
     & j1_ch(i),ka1_ch(i),kc1_ch(i),
     & j2_ch(i),ka2_ch(i),kc2_ch(i),
     & p1p2_bk(1,i),p1p2_bk(2,i),(-1)**kappa1,(-1)**kappa2,E_ch(i)
	  end if
      ENDIF	 
      ENDDO
      energy_defined = .TRUE.	  
      ENDIF
      IF(ini_chann_defined .and. .not.all_level_included) THEN
      chann_ini = 0	  
      DO i=1,number_of_channels
      IF(j1_ini.eq.j1_ch(i) .and. ka1_ini.eq.ka1_ch(i) .and. 
     & kc1_ini.eq.kc1_ch(i).and.j2_ini.eq.j2_ch(i).and.
     & ka2_ini.eq.ka2_ch(i).and.kc2_ini.eq.kc2_ch(i))
     & chann_ini = i
      ENDDO	
      IF(chann_ini.eq.0) STOP "ERROR: INITIAL CHANNEL WRONG"
      IF(MYID.eq.0) THEN      
      WRITE(1,'(a17,a5,i4,a8,i4,a8,i4,a7,i4,a8,i4,a8,i4)')
     & "INITIAL STATE IS:",
     & "J1 = ",j1_ini,"; KA1 = ", ka1_ini,"; KC1 =	",kc1_ini,
     & "; J2 =	", j2_ini, "; KA2 = ", ka2_ini, "; KC2 = ", eps2_ini
      ENDIF      
      ENDIF
      IF(.not.ini_chann_defined .and. channels_defined .and.
     & energy_defined) THEN
	  i_ground = 1
      DO i=1,number_of_channels
      IF(E_ch(i_ground).gt.E_ch(i)) i_ground = i   
      ENDDO	 
      chann_ini = i_ground
      ini_chann_defined = .TRUE.	  
      ENDIF	  
      IF(MYID.eq.0) WRITE(1,*)
      IF(MYID.eq.0) THEN	  
      IF(angs_unit)WRITE(1,'(a35)')"THE POTENTIAL R-UNITS ARE ANGSTROMS"
      IF(au_unit_r)WRITE(1,'(a31)') "THE POTENTIAL R-UNITS ARE BOHRS"
      IF(cm_unit)WRITE(1,'(a31)')"THE POTENTIAL ENERGY IS IN CM-1"
      IF(au_unit_e)WRITE(1,'(a34)')"THE POTENTIAL ENERGY IS IN HARTREE"
      IF(klv_unit)WRITE(1,'(a33)')"THE POTENTIAL ENERGY IS IN KELVIN"
      IF(kal_unit)WRITE(1,'(a35)')"THE POTENTIAL ENERGY IS IN KCAL/MOL"	  
      ENDIF     
      IF(.not.angs_unit .and..not.au_unit_r )   
     & STOP "ERROR:R_UNITS NOT DEFINED"
      IF(.not.cm_unit .and..not.au_unit_e.and..not.klv_unit
     & .and. .not.kal_unit)   
     & STOP "ERROR:ENERGY NOT DEFINED"
      IF(MYID.eq.0)WRITE(1,*) 	 
      IF(expansion_defined) THEN
      IF(MYID.eq.0) THEN	  
      WRITE(1,'(a42)') "USER SUPPLIED POTENTIAL TERMS WILL BE USED"
      WRITE(1,*)
      WRITE(1,'(a45,i4)')
     & "NUMBER OF POTENTIAL TERMS INCLUDED N_TERMS = ",nterms
      WRITE(1,'(a45)')
     & "THE POTENTIAL TERMS ARE, (L1,NJU1,L2,NJU2,L):"
      DO i=1,nterms
      WRITE(1,'(a1,i3,a1,i3,a1,i3,a1,i3,a1,i3,a3)',ADVANCE ="NO")
     & '(',L1_TERM(i),',',NJU1_TERM(i),',',L2_TERM(i),
     & ',',NJU2_TERM(i),',',L_TERM(i),'); '	  
      ENDDO	
      WRITE(1,*)
      WRITE(1,*)	  
      WRITE(1,'(a33,i3)') "NUMBER OF R-GRID POINTS IS N_R = ",
     & n_r_coll
      ENDIF	 
      ELSE
      IF(.not.matrix_reading_defined .and. MYID.eq.0)	  
     & WRITE(1,'(a97)')"COULPING MATRIX WILL BE COMPUTED USING NUMERICAL
     &  INTEGRATION OF POTENTIAL OVER THE WAVEFUNCTIONS"	  
      ENDIF	  
      IF(MYID.eq.0)WRITE(1,*)  	  
      END SELECT
      IF(.not.(energy_defined.and.channels_defined)) THEN
      PRINT*,"ERROR:BASIS SET NOT DEFINED. PROGRAM WILL STOP"
      STOP      	  
      ENDIF      

!Changing matrix file path 
!so that user don't need to copy the large matrix file in each directory.
	  if(bikram_mtrx_path) then
	  if(myid.eq.0) write(*,'(a,a,a)') "Matrix will be read/write ",
     & "from following directory: ", trim(bk_matrix_path11)
	  write(bk_temp_dir,'(a,a,a)')trim(bk_matrix_path11),'/',
     & trim(bk_dir11)
	  bk_dir1 = trim(bk_temp_dir)
	  write(bk_temp_dir,'(a,a,a)')trim(bk_matrix_path11),'/',
     & trim(bk_dir22)
	  bk_dir2 = trim(bk_temp_dir)
	  else
	  bk_dir1 = trim(bk_dir11)
	  bk_dir2 = trim(bk_dir22)
	  end if

      IF(MYID.eq.0) THEN	  
      IF(matrix_reading_defined) WRITE(1,'(a33)')
     & "MATRIX Mij WILL BE READ FROM FILE"
      IF(run_prog_defined) WRITE(1,'(a41)')
     & "SCATTERING CALCULATIONS WILL BE PERFORMED"	  
      WRITE(1,*)
      WRITE(1,'(a33)') "THE SYSTEM SETUP DONE SUCCESFULLY."
      CALL FLUSH(1)	  
      CLOSE(1)
      ENDIF	  
      END SUBROUTINE OUTPUT	        	  
      SUBROUTINE SYMMETRY_SORT_PARA_ORTHO
      USE CONSTANTS	  !!!! ANALYZING OUTPUT - CREATING ARRAYS-BASIS SET
      USE VARIABLES
      USE MPI_DATA	  
      IMPLICIT NONE
      INTEGER numb_chann_symm, i
      INTEGER symm_coeff,symm_coeff1,symm_coeff2	  
      INTEGER, ALLOCATABLE :: chan_incl_sym(:)
      INTEGER, ALLOCATABLE :: buffer_sym(:,:)
      REAL*8, ALLOCATABLE :: energy_sym(:)
      ALLOCATE(chan_incl_sym(number_of_channels))
      ALLOCATE(energy_sym(number_of_channels))	  
      SELECT CASE(coll_type)
      CASE(1)
      ALLOCATE(buffer_sym(1,number_of_channels)) 
      buffer_sym(1,:) = j_ch
      energy_sym = E_ch	  
      chan_incl_sym = 0	  
      numb_chann_symm = 0
      DO i = 1,number_of_channels
      symm_coeff = j_ch(i)-j_ini	  
      IF((symm_coeff/2)*2 .eq. symm_coeff) THEN
      numb_chann_symm = numb_chann_symm + 1	 	  
      chan_incl_sym(numb_chann_symm) = i
      ENDIF 	  
      ENDDO
      DEALLOCATE(j_ch,E_ch)
      number_of_channels = numb_chann_symm
      ALLOCATE(j_ch(number_of_channels),E_ch(number_of_channels))	  
      DO i=1,number_of_channels
      IF(chan_incl_sym(i).eq.0) STOP "ERROR IN PARA_ORTHO_SYMMETRY"	  
      j_ch(i) = buffer_sym(1,chan_incl_sym(i))
      E_ch(i) = energy_sym(chan_incl_sym(i))	  
      ENDDO
      CASE(2)	  
      ALLOCATE(buffer_sym(2,number_of_channels)) 
      buffer_sym(1,:) = j_ch
      buffer_sym(2,:) = v_ch	  
      energy_sym = E_ch	  
      chan_incl_sym = 0	  
      numb_chann_symm = 0
      DO i = 1,number_of_channels
      symm_coeff = j_ch(i)-j_ini	  
      IF((symm_coeff/2)*2 .eq. symm_coeff) THEN
      numb_chann_symm = numb_chann_symm + 1	 	  
      chan_incl_sym(numb_chann_symm) = i
      ENDIF 	  
      ENDDO
      DEALLOCATE(j_ch,v_ch,E_ch)
      number_of_channels = numb_chann_symm
      ALLOCATE(j_ch(number_of_channels),
     & v_ch(number_of_channels), E_ch(number_of_channels))	  
      DO i=1,number_of_channels
      IF(chan_incl_sym(i).eq.0) STOP "ERROR IN PARA_ORTHO_SYMMETRY"	  
      j_ch(i) = buffer_sym(1,chan_incl_sym(i))
      v_ch(i) = buffer_sym(2,chan_incl_sym(i))	  
      E_ch(i) = energy_sym(chan_incl_sym(i))	  
      ENDDO
      CASE(3)
      ALLOCATE(buffer_sym(3,number_of_channels)) 
      buffer_sym(1,:) = j_ch
      buffer_sym(2,:) = k_ch 
      buffer_sym(3,:) = eps_ch
      energy_sym = E_ch	  
      chan_incl_sym = 0	  
      numb_chann_symm = 0
      DO i = 1,number_of_channels
      symm_coeff = eps_ch(i) - eps_ini  
      IF((symm_coeff/2)*2 .eq. symm_coeff) THEN
      numb_chann_symm = numb_chann_symm + 1	 	  
      chan_incl_sym(numb_chann_symm) = i
      ENDIF 	  
      ENDDO
      DEALLOCATE(j_ch,k_ch,eps_ch,E_ch)
      number_of_channels = numb_chann_symm
      ALLOCATE(j_ch(number_of_channels),k_ch(number_of_channels),
     & eps_ch(number_of_channels),E_ch(number_of_channels))	  
      DO i=1,number_of_channels
      IF(chan_incl_sym(i).eq.0) STOP "ERROR IN PARA_ORTHO_SYMMETRY"	  
      j_ch(i) = buffer_sym(1,chan_incl_sym(i))
      k_ch(i) = buffer_sym(2,chan_incl_sym(i))
      eps_ch(i) = buffer_sym(3,chan_incl_sym(i))
      E_ch(i) = energy_sym(chan_incl_sym(i))	  
      ENDDO	  
      CASE(4)
      ALLOCATE(buffer_sym(3,number_of_channels)) 
      buffer_sym(1,:) = j_ch
      buffer_sym(2,:) = ka_ch 
      buffer_sym(3,:) = kc_ch
      energy_sym = E_ch	  
      chan_incl_sym = 0	  
      numb_chann_symm = 0
      DO i = 1,number_of_channels
      symm_coeff = (ka_ch(i)-kc_ch(i))-
     & (ka_ini-kc_ini)	  
      IF((symm_coeff/2)*2 .eq. symm_coeff) THEN
      numb_chann_symm = numb_chann_symm + 1	 	  
      chan_incl_sym(numb_chann_symm) = i
      ENDIF 	  
      ENDDO
      DEALLOCATE(j_ch,ka_ch,kc_ch,E_ch)
      number_of_channels = numb_chann_symm
      ALLOCATE(j_ch(number_of_channels),ka_ch(number_of_channels),
     & kc_ch(number_of_channels),E_ch(number_of_channels))	  
      DO i=1,number_of_channels
      IF(chan_incl_sym(i).eq.0) STOP "ERROR IN PARA_ORTHO_SYMMETRY"	  
      j_ch(i) = buffer_sym(1,chan_incl_sym(i))
      ka_ch(i) = buffer_sym(2,chan_incl_sym(i))
      kc_ch(i) = buffer_sym(3,chan_incl_sym(i))
      E_ch(i) = energy_sym(chan_incl_sym(i))	  
      ENDDO
      CASE(5)
      ALLOCATE(buffer_sym(2,number_of_channels))
      buffer_sym(1,:) = j1_ch
      buffer_sym(2,:) = j2_ch         
      energy_sym = E_ch	  
      chan_incl_sym = 0	  
      numb_chann_symm = 0
      DO i = 1,number_of_channels
      symm_coeff1 = j1_ch(i) - j1_ini
      symm_coeff2 = j2_ch(i) - j2_ini
      IF((symm_coeff1/2)*2 .eq. symm_coeff1) THEN
      IF((symm_coeff2/2)*2 .eq. symm_coeff2) THEN	  
      numb_chann_symm = numb_chann_symm + 1	 	  
      chan_incl_sym(numb_chann_symm) = i
      ENDIF
      ENDIF 	  
      ENDDO
      DEALLOCATE(j1_ch,j2_ch,E_ch)
      number_of_channels = numb_chann_symm
      ALLOCATE(j1_ch(number_of_channels),E_ch(number_of_channels),
     &  j2_ch(number_of_channels))	  
      DO i=1,number_of_channels
      IF(chan_incl_sym(i).eq.0) STOP "ERROR IN PARA_ORTHO_SYMMETRY"	  
      j1_ch(i) = buffer_sym(1,chan_incl_sym(i))
      j2_ch(i) = buffer_sym(2,chan_incl_sym(i))	  
      E_ch(i) = energy_sym(chan_incl_sym(i))	  
      ENDDO	  	  
      CASE(6)
      ALLOCATE(buffer_sym(4,number_of_channels))
      buffer_sym(1,:) = j1_ch
      buffer_sym(2,:) = j2_ch
      buffer_sym(3,:) = v1_ch
      buffer_sym(4,:) = v2_ch  	  
      energy_sym = E_ch	  
      chan_incl_sym = 0	  
      numb_chann_symm = 0
      DO i = 1,number_of_channels
      symm_coeff1 = j1_ch(i) - j1_ini
      symm_coeff2 = j2_ch(i) - j2_ini
      IF((symm_coeff1/2)*2 .eq. symm_coeff1) THEN
      IF((symm_coeff2/2)*2 .eq. symm_coeff2) THEN	  
      numb_chann_symm = numb_chann_symm + 1	 	  
      chan_incl_sym(numb_chann_symm) = i
      ENDIF
      ENDIF 	  
      ENDDO
      DEALLOCATE(j1_ch,v1_ch,v2_ch,j2_ch,E_ch)
      number_of_channels = numb_chann_symm
      ALLOCATE(j1_ch(number_of_channels),E_ch(number_of_channels),
     &  j2_ch(number_of_channels),v1_ch(number_of_channels),
     &	v2_ch(number_of_channels))	  
      DO i=1,number_of_channels
      IF(chan_incl_sym(i).eq.0) STOP "ERROR IN PARA_ORTHO_SYMMETRY"	  
      j1_ch(i) = buffer_sym(1,chan_incl_sym(i))
      j2_ch(i) = buffer_sym(2,chan_incl_sym(i))
      v1_ch(i) = buffer_sym(3,chan_incl_sym(i))
      v2_ch(i) = buffer_sym(4,chan_incl_sym(i))	  
      E_ch(i) = energy_sym(chan_incl_sym(i))	  
      ENDDO		  
      CASE(7)
      ALLOCATE(buffer_sym(4,number_of_channels))
      buffer_sym(1,:) = j1_ch
      buffer_sym(2,:) = k1_ch 
      buffer_sym(3,:) = eps1_ch
      buffer_sym(4,:) = j2_ch         
      energy_sym = E_ch	  
      chan_incl_sym = 0	  
      numb_chann_symm = 0
      DO i = 1,number_of_channels
      symm_coeff1 = eps1_ch(i)-eps1_ini
      symm_coeff2 = j2_ch(i) - j2_ini
      IF((symm_coeff1/2)*2 .eq. symm_coeff1) THEN
      IF((symm_coeff2/2)*2 .eq. symm_coeff2) THEN	  
      numb_chann_symm = numb_chann_symm + 1	 	  
      chan_incl_sym(numb_chann_symm) = i
      ENDIF
      ENDIF 	  
      ENDDO
      DEALLOCATE(j1_ch,k1_ch,eps1_ch,j2_ch,E_ch)
      number_of_channels = numb_chann_symm
      ALLOCATE(j1_ch(number_of_channels),k1_ch(number_of_channels),
     & eps1_ch(number_of_channels),E_ch(number_of_channels),
     &  j2_ch(number_of_channels))	  
      DO i=1,number_of_channels
      IF(chan_incl_sym(i).eq.0) STOP "ERROR IN PARA_ORTHO_SYMMETRY"	  
      j1_ch(i) = buffer_sym(1,chan_incl_sym(i))
      k1_ch(i) = buffer_sym(2,chan_incl_sym(i))
      eps1_ch(i) = buffer_sym(3,chan_incl_sym(i))
      j2_ch(i) = buffer_sym(4,chan_incl_sym(i))	  
      E_ch(i) = energy_sym(chan_incl_sym(i))	  
      ENDDO	  
      CASE(8)
      ALLOCATE(buffer_sym(4,number_of_channels))
      buffer_sym(1,:) = j1_ch
      buffer_sym(2,:) = ka1_ch 
      buffer_sym(3,:) = kc1_ch
      buffer_sym(4,:) = j2_ch         
      energy_sym = E_ch	  
      chan_incl_sym = 0	  
      numb_chann_symm = 0
      DO i = 1,number_of_channels
      symm_coeff1 = (ka1_ch(i)-kc1_ch(i))-
     & (ka1_ini-kc1_ini)
      symm_coeff2 = j2_ch(i) - j2_ini
      IF((symm_coeff1/2)*2 .eq. symm_coeff1) THEN
      IF((symm_coeff2/2)*2 .eq. symm_coeff2) THEN	  
      numb_chann_symm = numb_chann_symm + 1	 	  
      chan_incl_sym(numb_chann_symm) = i
      ENDIF
      ENDIF 	  
      ENDDO
      DEALLOCATE(j1_ch,ka1_ch,kc1_ch,j2_ch,E_ch)
      number_of_channels = numb_chann_symm
      ALLOCATE(j1_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),E_ch(number_of_channels),
     &  j2_ch(number_of_channels))	  
      DO i=1,number_of_channels
      IF(chan_incl_sym(i).eq.0) STOP "ERROR IN PARA_ORTHO_SYMMETRY"	  
      j1_ch(i) = buffer_sym(1,chan_incl_sym(i))
      ka1_ch(i) = buffer_sym(2,chan_incl_sym(i))
      kc1_ch(i) = buffer_sym(3,chan_incl_sym(i))
      j2_ch(i) = buffer_sym(4,chan_incl_sym(i))	  
      E_ch(i) = energy_sym(chan_incl_sym(i))	  
      ENDDO	  
      CASE(9)
      ALLOCATE(buffer_sym(6,number_of_channels))
      buffer_sym(1,:) = j1_ch
      buffer_sym(2,:) = ka1_ch 
      buffer_sym(3,:) = kc1_ch
      buffer_sym(4,:) = j2_ch
      buffer_sym(5,:) = k2_ch 
      buffer_sym(6,:) = eps2_ch	  
      energy_sym = E_ch	  
      chan_incl_sym = 0	  
      numb_chann_symm = 0
      DO i = 1,number_of_channels
      symm_coeff = (ka1_ch(i)-kc1_ch(i))-
     & (ka1_ini-kc1_ini)	  
      IF((symm_coeff/2)*2 .eq. symm_coeff) THEN
      numb_chann_symm = numb_chann_symm + 1	 	  
      chan_incl_sym(numb_chann_symm) = i
      ENDIF 	  
      ENDDO
      DEALLOCATE(j1_ch,ka1_ch,kc1_ch,
     & j2_ch,k2_ch,eps2_ch,E_ch)
      number_of_channels = numb_chann_symm
      ALLOCATE(j1_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),E_ch(number_of_channels),
     &  j2_ch(number_of_channels),k2_ch(number_of_channels),
     & eps2_ch(number_of_channels) )	  
      DO i=1,number_of_channels
      IF(chan_incl_sym(i).eq.0) STOP "ERROR IN PARA_ORTHO_SYMMETRY"	  
      j1_ch(i) = buffer_sym(1,chan_incl_sym(i))
      ka1_ch(i) = buffer_sym(2,chan_incl_sym(i))
      kc1_ch(i) = buffer_sym(3,chan_incl_sym(i))
      j2_ch(i) = buffer_sym(4,chan_incl_sym(i))
      k2_ch(i) = buffer_sym(5,chan_incl_sym(i))
      eps2_ch(i) = buffer_sym(6,chan_incl_sym(i))	  
      E_ch(i) = energy_sym(chan_incl_sym(i))
      ENDDO	  
      CASE(0)
      buffer_sym(1,:) = j1_ch
      buffer_sym(2,:) = ka1_ch 
      buffer_sym(3,:) = kc1_ch
      buffer_sym(4,:) = j2_ch
      buffer_sym(5,:) = ka2_ch 
      buffer_sym(6,:) = kc2_ch	  
      energy_sym = E_ch	  
      chan_incl_sym = 0	  
      numb_chann_symm = 0
      DO i = 1,number_of_channels
      symm_coeff1 = (ka1_ch(i)-kc1_ch(i))-
     & (ka1_ini-kc1_ini)
      symm_coeff2 = (ka2_ch(i)-kc2_ch(i))-
     & (ka2_ini-kc2_ini)	 
      IF(((symm_coeff1/2)*2 .eq. symm_coeff1) .and.
     & ((symm_coeff2/2)*2 .eq. symm_coeff2)) THEN
      numb_chann_symm = numb_chann_symm + 1	 	  
      chan_incl_sym(numb_chann_symm) = i
      ENDIF 	  
      ENDDO
      DEALLOCATE(j1_ch,ka1_ch,kc1_ch,
     & j2_ch,ka2_ch,kc2_ch,E_ch)
      number_of_channels = numb_chann_symm
      ALLOCATE(j1_ch(number_of_channels),ka1_ch(number_of_channels),
     & kc1_ch(number_of_channels),E_ch(number_of_channels),
     &  j2_ch(number_of_channels),ka2_ch(number_of_channels),
     & kc2_ch(number_of_channels) )	  
      DO i=1,number_of_channels
      IF(chan_incl_sym(i).eq.0) STOP "ERROR IN PARA_ORTHO_SYMMETRY"	  
      j1_ch(i) = buffer_sym(1,chan_incl_sym(i))
      ka1_ch(i) = buffer_sym(2,chan_incl_sym(i))
      kc1_ch(i) = buffer_sym(3,chan_incl_sym(i))
      j2_ch(i) = buffer_sym(4,chan_incl_sym(i))
      ka2_ch(i) = buffer_sym(5,chan_incl_sym(i))
      kc2_ch(i) = buffer_sym(6,chan_incl_sym(i))	  
      E_ch(i) = energy_sym(chan_incl_sym(i))	  
      ENDDO	  
      END SELECT	  
      END SUBROUTINE SYMMETRY_SORT_PARA_ORTHO
