      SUBROUTINE PRINT_OUTPUT
      USE INTEGRATOR
      USE VARIABLES
      USE MONTE_CARLO_STORAGE	  
      USE MPI_DATA
      USE MPI
      IMPLICIT NONE
      LOGICAL file_exst	  
      INTEGER i,stat_tmp
      INTEGER HOURS, MINUTES
      REAL*8  SECONDS
      REAL*8 MAX_CLOCK_TIME
      INTEGER HOURS_MAX, MINUTES_MAX
      REAL*8  SECONDS_MAX
      INTEGER ID_MAX,N_TRJ_MONT_CARLO_NEED
      REAL*8 scnds
      REAL*8, ALLOCATABLE :: mont_carlo_max(:)	  
      INTEGER mnts,hrs,i_u
      INTEGER act_nmb_of_traject  
      CALL PRINT_ERRORS 	  
!      INQUIRE( FILE="CROSS_SECTIONS.out", EXIST=file_exst )
      IF(non_format_out_defined) THEN
	  CALL NONFORMAT_PRINTING
      IF(myid.eq.0) PRINT*, "DATA_PRINTED IN NON_FORMAT"	  
      RETURN
      ENDIF	
 	  IF(myid.eq.0) THEN
      IF(monte_carlo_defined) THEN
      ALLOCATE(mont_carlo_max(nmbr_of_enrgs))
      mont_carlo_max = 0d0	  
      ENDIF	 
	
! Bikram Oct'18 Start:
      inquire(file='bikram_mqct.tmp',exist=file_exst)
	  if(file_exst) then
	  open(1991,file='bikram_mqct.tmp')
	  close(1991,status='delete')
	  else
	  open(1992,file='RESONANCE_TRAJECT.out')
	  close(1992,status='delete')
	  endif	
      INQUIRE(FILE="CROSS_SECTIONS.tmp", EXIST=file_exst )	
      IF(file_exst) open(1111,FILE="CROSS_SECTIONS.tmp")
      close(1111,status='delete')	
! Bikram End.
      OPEN(1,FILE="CROSS_SECTIONS.out",POSITION="APPEND")
      DO i_u=i_curr, i_ener-1
      WRITE(1,'(a,f14.6,a)') "U= ",U(i_ener)," cm^-1"
      IF(.not.monte_carlo_defined) THEN	  
      WRITE(1,'(2(a6,2x),a19,2x,a14,2x,a19)') 
     & 'ilv', 'flv', 'sigma(U),ANG^2',
     & 'E_coll,cm^-1', 'sigma(E_coll),ANG^2'
      DO i = 1,number_of_channels    	  
      IF(bill_exst(i,i_u))
     & WRITE(1,'(i6,2x,i6,2x,e19.10,2x,f14.6,2x,e19.10)') 
     & chann_ini, i, sigma_f(i,i_u)/U(i_u)*E_bill(i,i_u), 
     & E_bill(i,i_u), sigma_f(i,i_u)
      ENDDO
	  
      ELSE 
      WRITE(1,'(2(a6,2x),a19,2x,a14,2x,a19)') 
     & 'ilv', 'flv', 'sigma(U),ANG^2',
     & 'E_coll,cm^-1', 'sigma(E_coll),ANG^2'
      DO i = 1,number_of_channels    	  
      IF(bill_exst(i,i_u)) 
     & WRITE(1,'(i6,2x,i6,2x,e19.10,2x,f14.6,2x,e19.10)')
     & chann_ini,i, sigma_f(i,i_u)/U(i_u)*E_bill(i,i_u), 
     & E_bill(i,i_u),sigma_f(i,i_u)
      IF(mont_carlo_max(i_u).lt.monte_carlo_err(i,i_u)
     & .and.bill_exst(i,i_u))
     & 	mont_carlo_max(i_u) = monte_carlo_err(i,i_u)  
      ENDDO
      IF(mnt_crl_intgrt_err.gt.0d0)	THEN  
      N_TRJ_MONT_CARLO_NEED  =	int(dble(tot_number_of_traject)*
     & (mont_carlo_max(i_u)/mnt_crl_intgrt_err)**2)
      WRITE(1,'(a33,1x,a25,1x,i9)')
     & "NUMBER OF TRAJECTORIES NEEDED FOR",
     & "THE DESIRABLE ACCURACY IS",
     & N_TRJ_MONT_CARLO_NEED	
      ENDIF	 
      ENDIF	  

! Bikram Start May 2020:
      WRITE(1,'(a28,1x,e12.5)') "TOTAL P CONSERVATION ERROR,%",
     & error_prb_largest(i_u)	  
	  if(bk_nrg_err) then
      WRITE(1,'(a28,1x,e12.5)') "TOTAL E CONSERVATION ERROR,%",
     &	error_ener_largest(i_u)
      WRITE(1,*)
	  endif
! Bikram End.
	  
      ENDDO
      ENDIF      
      IF(myid.eq.0) THEN
      TIME_SPEND = 0d0
      MAX_CLOCK_TIME = 0d0	  
      DO id_cnt=1,nproc
      TIME_SPEND = TIME_SPEND + TOTAL_TIME(id_cnt)
      IF(MAX_CLOCK_TIME.LE.TOTAL_TIME(id_cnt)) THEN
      MAX_CLOCK_TIME=TOTAL_TIME(id_cnt)
      ID_MAX = id_cnt - 1	  
      ENDIF	  
      ENDDO	
      act_nmb_of_traject	=
     & nproc  
      MINUTES = INT(INT(TIME_SPEND/dble(act_nmb_of_traject))/60)
      SECONDS = TIME_SPEND/act_nmb_of_traject - 60d0*MINUTES	  
      HOURS = INT(MINUTES/60)
      MINUTES = MINUTES-60*HOURS

      MINUTES_MAX = INT(INT(MAX_CLOCK_TIME)/60)
      SECONDS_MAX = MAX_CLOCK_TIME - 60d0*MINUTES_MAX	  
      HOURS_MAX = INT(MINUTES_MAX/60)
      MINUTES_MAX = MINUTES_MAX-60*HOURS_MAX
      IF(TIME_SPEND>10d0) THEN	  
      WRITE(*,'(a38,i8)')"TOTAL TIME SPENT ON PROPAGATION,sec = ",
     & INT(TIME_SPEND)
      ELSE
      WRITE(*,'(a38,f6.3)')"TOTAL TIME SPENT ON PROPAGATION,sec = ",
     & TIME_SPEND	  
      ENDIF	  
      WRITE(1,'(a31,i3)')"TOTAL NUMBER OF TRAJECTORIES = ",
     & tot_number_of_traject!/mpi_task_per_traject
      WRITE(*,'(a53,1x,i3,1x,a2,1x,i3,2x,a2,1x,f4.1,1x,a2)')
     & "THE AVERAGE TIME PER MPI TASK SPENT ON PROPAGATION = ",
     & HOURS," h",MINUTES," m", SECONDS, " s"
!      WRITE(1,'(a22,i4,a36,1x,i3,1x,a2,1x,i3,2x,a2,1x,f4.1,1x,a2)')
!     & "THE LONGEST MPI TASK #",ID_MAX ,
!     & " TOOK PROPAGATION WALL CLOCK TIME = ",
!     & HOURS_MAX," h",MINUTES_MAX," m", SECONDS_MAX, " s"
      CALL TIME_CONVERT(TIME_INI,hrs,mnts,scnds)
      WRITE(*,'(a28,a24,1x,i3,1x,a2,1x,i3,2x,a2,1x,f4.1,1x,a2)')
     & "TIME SPENT ON INITIALIZATION",
     & " TOOK WALL CLOCK TIME = ",
     & hrs," h",mnts," m", scnds, " s"	  
      IF(matrix_reading_defined) THEN
      CALL TIME_CONVERT(TIME_1_MAT,hrs,mnts,scnds)	  
      WRITE(*,'(a29,a24,1x,i3,1x,a2,1x,i3,2x,a2,1x,f4.1,1x,a2)')
     & "TIME SPENT ON MIJ.DAT READING",
     & " TOOK WALL CLOCK TIME = ",
     & hrs," h",mnts," m", scnds, " s"	
      ENDIF	 
      IF(mpi_task_defined) THEN
      CALL TIME_CONVERT(TIME_2_MAT,hrs,mnts,scnds)	  
      WRITE(*,'(a30,a24,1x,i3,1x,a2,1x,i3,2x,a2,1x,f4.1,1x,a2)')
     & "TIME SPENT ON MIJ DISTRIBUTION",
     & " TOOK WALL CLOCK TIME = ",
     & hrs," h",mnts," m", scnds, " s"
      ENDIF
      CALL TIME_CONVERT(TIME_3_MAT,hrs,mnts,scnds)	  
      WRITE(*,'(a30,a24,1x,i3,1x,a2,1x,i3,2x,a2,1x,f4.1,1x,a2)')
     & "TIME SPENT ON MIJ REALLOCATION",
     & " TOOK WALL CLOCK TIME = ",
     & hrs," h",mnts," m", scnds, " s"
      CLOSE(1)	 
      ENDIF

      IF(INTERUPT_PROPAGATION) THEN
      IF(MYID.eq.0) THEN
      OPEN(1,FILE="CROSS_SECTIONS.out",POSITION="APPEND")	  
      WRITE(1,'(a30,1x,a33,1x,a28,1x,a20)')
     ^ "ATTENTION: NOT ALL CALULATIONS"
     & ,"WERE COMPLETED DUE TO TIME LIMIT.", 
     & "CHECKPOINT FILE WAS CREATED."
     & ,"RESTART THE PROGRAM."
      CLOSE(1)
      WRITE(*,'(a30,1x,a33,1x,a28,1x,a20)')
     ^ "ATTENTION: NOT ALL CALULATIONS"
     & ,"WERE COMPLETED DUE TO TIME LIMIT.", 
     & "CHECKPOINT FILE WAS CREATED."
     & ,"RESTART THE PROGRAM."	  
      ENDIF	  
      RETURN      
      ENDIF
      IF(MYID.eq.0) PRINT*, "ALL WORK DONE" !!! WORK IS DONE
      END SUBROUTINE PRINT_OUTPUT

      SUBROUTINE TIME_CONVERT(time,h,m,s)
      IMPLICIT NONE
      REAL*8 time,s
      INTEGER h,m
      m = int(time/60d0)
      s = time - m*60d0
      h = int(m/60)
      m = m-60*h 	  
      END SUBROUTINE TIME_CONVERT	  
	 
      SUBROUTINE NONFORMAT_PRINTING
      USE INTEGRATOR
      USE VARIABLES
      USE MONTE_CARLO_STORAGE	  
      USE MPI_DATA
      USE MPI
      IMPLICIT NONE
      LOGICAL file_exst	  
      INTEGER i
      INTEGER HOURS, MINUTES
      REAL*8  SECONDS
      REAL*8 MAX_CLOCK_TIME
      INTEGER HOURS_MAX, MINUTES_MAX
      REAL*8  SECONDS_MAX
      INTEGER ID_MAX,N_TRJ_MONT_CARLO_NEED
      REAL*8 scnds
      REAL*8, ALLOCATABLE :: mont_carlo_max(:)	  
      INTEGER mnts,hrs,i_u
      INTEGER act_nmb_of_traject
      IF(myid.eq.0) THEN
      IF(monte_carlo_defined) THEN
      ALLOCATE(mont_carlo_max(nmbr_of_enrgs))
      mont_carlo_max = 0d0	  
      ENDIF	  
      OPEN(1,FILE="CROSS_SECTIONS.out",POSITION="APPEND")
      DO i_u=i_curr, i_ener-1
!      WRITE(1,'(a3,1x,f10.3)') "U=",U(i_u)
      IF(.not.monte_carlo_defined) THEN	  
!      WRITE(1,'(a3,2x,a3,2x,a13,2x,a17)') 
!     & 'ilv','flv','E_coll,cm^-1','cross sect.,ANG^2'  
      DO i = 1,number_of_channels    	  
      IF(bill_exst(i,i_u)) 
     & WRITE(1,'(i3,2x,i3,2x,f13.3,2x,e17.10)') chann_ini, i,
     &	  E_bill(i,i_u),sigma_f(i,i_u)	  
      ENDDO
      ELSE
!      WRITE(1,'(a3,2x,a3,2x,a13,2x,a17,1x,a19)') 
!     & 'ilv','flv','E_coll,cm^-1','cross sect.,ANG^2',
!     & 'monte_carlo_error,%'  
      DO i = 1,number_of_channels    	  
      IF(bill_exst(i,i_u)) 
     & WRITE(1,'(i3,2x,i3,2x,f13.3,2x,e17.10,10x,f10.4)')chann_ini,i,
     & E_bill(i,i_u),sigma_f(i,i_u),monte_carlo_err(i,i_u)
      IF(mont_carlo_max(i_u).lt.monte_carlo_err(i,i_u)
     & .and.bill_exst(i,i_u))
     & 	mont_carlo_max(i_u) = monte_carlo_err(i,i_u)  
      ENDDO
      IF(mnt_crl_intgrt_err.gt.0d0)	THEN  
      N_TRJ_MONT_CARLO_NEED  =	int(dble(tot_number_of_traject)*
     & (mont_carlo_max(i_u)/mnt_crl_intgrt_err)**2)
!      WRITE(1,'(a33,1x,a25,1x,i9)')
!     & "NUMBER OF TRAJECTORIES NEEDED FOR",
!     & "THE DESIRABLE ACCURACY IS",
!     & N_TRJ_MONT_CARLO_NEED	
      ENDIF	 
      ENDIF	  
!      WRITE(1,'(a28,1x,e12.5)') "TOTAL P CONSERVATION ERROR,%",
!     & error_prb_largest(i_u)
!      WRITE(1,'(a28,1x,e12.5)') "TOTAL E CONSERVATION ERROR,%",
!     &	error_ener_largest(i_u)
      WRITE(1,*)
      ENDDO
      ENDIF  

	  
      END SUBROUTINE NONFORMAT_PRINTING	  	 