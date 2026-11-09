	  module serial_parallel
! This module is written by Bikramaditya Mandal
! The purpose for this module is to make the code running as sequential
	  implicit none
	  integer definer, bk_funit
	  logical parallel_defined, definer_file_exst
	  end module

      MODULE MPI_DATA
! This module is updated by Bikramaditya Mandal
      IMPLICIT NONE	  
      INTEGER myid,nproc,ierr_mpi,id_cnt
      REAL*8 TIME_START,TIME_FINISH,TIME_SPEND
      REAL*8,ALLOCATABLE ::	  TOTAL_TIME(:), TOTAL_TIME_INI(:)
      REAL*8 TIME_BEGIN,TIME_CURRENT,TIME_1_MAT,TIME_2_MAT,TIME_3_MAT
      REAL*8 TIME_1_ARRAY
      REAL*8 TIME_INI,TIME_PROPAGATE,TIME_DOTRAJECT	  
! R outside the Mij grid [RMIN,RMAX], set in potential_data_bk/_adia
! (mqct_eqtns.f) and reported once at the end of the run (below).
! R_below_grid = lowest R below RMIN (stays 1d30 if R never went below)
! R_above_grid = highest R above RMAX (stays -1d30 if R never went above)
! R_grid_min, R_grid_max = RMIN, RMAX of the grid (0 until R leaves it)
      REAL*8 :: R_below_grid = 1d30, R_above_grid = -1d30
      REAL*8 :: R_grid_min = 0d0, R_grid_max = 0d0
      END MODULE MPI_DATA
	  
      PROGRAM MQCT_HEAD
! This module is updated by Bikramaditya Mandal
	  use serial_parallel
      USE MPI_DATA
      USE MPI
      REAL*8 R_lowest, R_highest, RMIN_grid, RMAX_grid
	  
! Bikram July 2020:
! This part is to create conditions so that code can run in sequential
	  definer = 0
	  parallel_defined = .true.
      inquire(file='serial_parallel_input',exist=definer_file_exst)
	  if(definer_file_exst) then
	  open(100001, file='serial_parallel_input')
	  read(100001,*) definer
	  close(100001)
	  endif
	  if(definer.eq.1) parallel_defined = .false.
	  if(.not. parallel_defined) then
	  ierr_mpi = 0
	  nproc = 1
	  myid = 0
	  endif
	  
	  if(parallel_defined) then
      CALL MPI_INIT( ierr_mpi ) !/MPI INITIALIZATION
      IF(ierr_mpi.ne.0) THEN 
      WRITE(*,'(a22,1x,i5)')'ERROR: MPI_INI_FAIL = ',ierr_mpi 
      STOP
      ENDIF
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr_mpi ) !/MPI NUMBER OF PROCESSORS INI
      IF(ierr_mpi.ne.0) THEN 
      WRITE(*,'(a22,1x,i5)')'ERROR: MPI_INI_FAIL = ',ierr_mpi 
      STOP
      ENDIF
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr_mpi ) !/ MPI PROCESSORS ID INI
      IF(ierr_mpi.ne.0) THEN 
      WRITE(*,'(a22,1x,i5)')'ERROR: MPI_INI_FAIL = ',ierr_mpi 
      STOP
      ENDIF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )     !/ MPI CHECKING IF ALL PROCESSORS READY
      IF(ierr_mpi.ne.0) THEN 
      WRITE(*,'(a22,1x,i5)')'ERROR: MPI_INI_FAIL = ',ierr_mpi 
      STOP
      ENDIF
	  endif
      IF(MYID.EQ.0) PRINT*, "MPI INITIALIZED"
c/ READING INPUT AND SETUP MATRIX Mij(R)
      IF(MYID.EQ.0) PRINT*, 
     & "PARSERING OF USER INPUT"
      if(parallel_defined) TIME_BEGIN = MPI_Wtime()
      CALL INITIALIZATION
!      CALL CPU_TIME(TIME_START)	  
      if(parallel_defined) TIME_START = 	MPI_Wtime()  ! TIME START
      TIME_INI = TIME_START-TIME_BEGIN ! TIME SPENT ON INITIALIZATION
      ALLOCATE(TOTAL_TIME_INI(nproc))
      if(parallel_defined) then 
	  CALL MPI_GATHER(TIME_INI,1,MPI_DOUBLE_PRECISION,
     & TOTAL_TIME_INI,
     &	  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)	  
	  else
	  TOTAL_TIME_INI(:) = TIME_INI
	  endif
c/ SAMPLING AND PROPAGATION OF MQCT EQUATIONS
      if(parallel_defined) TIME_PROPAGATE = MPI_Wtime() 	  
      CALL PROPAGATE
!      CALL CPU_TIME(TIME_SPEND)	  
      if(parallel_defined) TIME_FINISH = MPI_Wtime()  ! TIME FINISH
      TIME_SPEND = TIME_FINISH-TIME_START	  
      ALLOCATE(TOTAL_TIME(nproc)) ! GATHERING TIMES FROM ALL TASKS
      if(parallel_defined) CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      if(parallel_defined) then
	  CALL MPI_GATHER(TIME_SPEND,1,MPI_DOUBLE_PRECISION,
     & TOTAL_TIME,
     &	  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
	  else
	  !TOTAL_TIME(:) = TIME_SPEND
	  endif
! R outside the Mij grid: proc 0 collects from all processes
! the lowest R (MPI_MIN), the highest R and the grid edges (MPI_MAX),
! and prints each warning only once.
      R_lowest  = R_below_grid
      R_highest = R_above_grid
      RMIN_grid = R_grid_min
      RMAX_grid = R_grid_max
      if(parallel_defined) then
      CALL MPI_REDUCE(R_below_grid, R_lowest, 1, MPI_DOUBLE_PRECISION,
     & MPI_MIN, 0, MPI_COMM_WORLD, ierr_mpi)
      CALL MPI_REDUCE(R_above_grid, R_highest, 1, MPI_DOUBLE_PRECISION,
     & MPI_MAX, 0, MPI_COMM_WORLD, ierr_mpi)
      CALL MPI_REDUCE(R_grid_min, RMIN_grid, 1, MPI_DOUBLE_PRECISION,
     & MPI_MAX, 0, MPI_COMM_WORLD, ierr_mpi)
      CALL MPI_REDUCE(R_grid_max, RMAX_grid, 1, MPI_DOUBLE_PRECISION,
     & MPI_MAX, 0, MPI_COMM_WORLD, ierr_mpi)
      endif
      IF(MYID.EQ.0 .and. R_lowest.lt.RMIN_grid)
     & WRITE(*,'(a,f9.4,a,f9.4,a)')
     & ' WARNING: R went below RMIN =',RMIN_grid,
     & ' bohr of the Mij grid (lowest R =',R_lowest,
     & '). Mij was linearly extrapolated there.'
      IF(MYID.EQ.0 .and. R_highest.gt.RMAX_grid)
     & WRITE(*,'(a,f9.4,a,f9.4,a)')
     & ' WARNING: R went above RMAX =',RMAX_grid,
     & ' bohr of the Mij grid (highest R =',R_highest,
     & '). Mij used the C6/R^6 tail there.'
      CALL PRINT_OUTPUT ! PRINTING OUTPUT
      DEALLOCATE(TOTAL_TIME)	 !! MEMORY FREE
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi )	  
      CALL MPI_FINALIZE (ierr_mpi)		  
      END
	  

