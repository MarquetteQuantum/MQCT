28    SUBROUTINE INI_ARRAYS !!! remeber to change int( to round(
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA	  
      IMPLICIT NONE
      INTEGER i,st,j_count,j_summ,m_count,j1_count,j2_count
      INTEGER p_count, percent_counter
      INTEGER st1,st2,KRONEKER,p_lim_max_ini,round	 
	  INTEGER bk_global_parity(number_of_channels)								! Bikram April 2021
	  INTEGER bk_kappa(number_of_channels), kappa1, kappa2						! Bikram April 2021
      LOGICAL EVEN_NUM					  							 
      EXTERNAL KRONEKER,round	  
      st = 0
      IF(MYID.eq.0) PRINT *, "ARRAY_INI_STARTED"
      ALLOCATE(chann_indx(number_of_channels))
      IF(.not.identical_particles_defined) THEN	  
      p_lim_max = 1
      p_lim_min = 1	
      ELSE
      SELECT CASE(coll_type)		  
      CASE(5)	  

      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
      CASE(6)
      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))*
     & KRONEKER(v1_ch(chann_ini),v2_ch(chann_ini)) 	  
      CASE(0)
      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
     & *KRONEKER(ka1_ch(chann_ini),ka2_ch(chann_ini))
     & *KRONEKER(kc1_ch(chann_ini),kc2_ch(chann_ini))
      END SELECT
      ENDIF	
      IF(identical_particles_defined) THEN
      IF(coll_type.eq.5 .or. coll_type.eq.6 .or. coll_type.eq.0) THEN
      IF(myid.eq.0) THEN
      PRINT*,"EXCHANGE SYMMETRIES FOR IDENTICAL PARTICLES:"
      IF(exch_par_w_pl.eq.1d0) THEN
      PRINT*, "POSITIVE EXCHANGE PARITY ONLY"
      ENDIF	  
      IF(exch_par_w_pl.lt.1d0 .and. exch_par_w_pl.gt.0d0) THEN
      PRINT*, "BOTH EXCHANGE PARITIES"
      ENDIF	  
      IF(exch_par_w_pl.eq.0d0) THEN
      PRINT*, "NEGATIVE EXCHANGE PARITY ONLY"
      ENDIF	  	  
      ENDIF	 
      ELSE
      IF(myid.eq.0) 
     & PRINT*,
     & "ERROR: EXCHANGE SYMMETRY DEFINED FOR NON-IDENTICAL PARTICLES"
      STOP	
      ENDIF	  
      ENDIF	  
      chann_indx = 0!!! CREATING CHANNEL -> STATE	  
	  
      DO i=1,number_of_channels
      SELECT CASE(coll_type)
      CASE(1)
      IF(fine_structure_defined .and. LORB_FINE.gt.0 .and.
     & SPIN_FINE.eq.2 ) THEN
      DO j_count = 1, round(j_h_ch(i)*2) + 1	  
      st = st + 1
      IF(j_count*2-2 .eq. round(j_h_ch(i)*2) + 1) chann_indx(i) = st !!!! SLOWLY MODIFIYNG
      ENDDO	
      ELSE
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      IF(j_count.eq.0) chann_indx(i) = st
      ENDDO	  
      ENDIF	  
      CASE(2)
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      IF(j_count.eq.0) chann_indx(i) = st
      ENDDO
      CASE(3)
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      IF(j_count.eq.0) chann_indx(i) = st
      ENDDO
      CASE(4)
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      IF(j_count.eq.0) chann_indx(i) = st
      ENDDO
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i)))
     &	  chann_indx(i) = st
      ENDDO
      ENDDO
      ELSE
      p_lim_min = 1
      p_lim_max = 2-KRONEKER(j1_ch(i),j2_ch(i)) !min(2-KRONEKER(j1_ch(i),j2_ch(i)),p_lim_max_ini)
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch(i),j2_ch(i))	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and. p_count .eq.p_lim_min  )
     &	  chann_indx(i) = st
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i)))
     &	  chann_indx(i) = st
      ENDDO
      ENDDO
      ELSE
      p_lim_min = 1
      p_lim_max = 2-KRONEKER(j1_ch(i),j2_ch(i))
     & *KRONEKER(v1_ch(i),v2_ch(i))
!      p_lim_max = min(2-KRONEKER(j1_ch(i),j2_ch(i))*
!     & KRONEKER(v1_ch(i),v2_ch(i)),p_lim_max_ini) 
      DO p_count=p_lim_min,p_lim_max	
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and.  p_count .eq.p_lim_min  )
     &	  chann_indx(i) = st
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	  
      CASE(7)
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i)))
     &	  chann_indx(i) = st
      ENDDO
      ENDDO
      CASE(8)
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i)))
     &	  chann_indx(i) = st
      ENDDO
      ENDDO
      CASE(9)
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i)))
     &	  chann_indx(i) = st
      ENDDO
      ENDDO
      CASE(0)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i)))
     &	  chann_indx(i) = st
      ENDDO
      ENDDO
      ELSE
      p_lim_min = 1
!      p_lim_max = min(2-KRONEKER(j1_ch(i),j2_ch(i))
!     & *KRONEKER(ka1_ch(i),ka2_ch(i))*KRONEKER(kc1_ch(i),kc2_ch(i)),
!     & p_lim_max_ini)  
      p_lim_max = 2-KRONEKER(j1_ch(i),j2_ch(i))
     & *KRONEKER(ka1_ch(i),ka2_ch(i))*KRONEKER(kc1_ch(i),kc2_ch(i))
	  
      DO p_count=p_lim_min,p_lim_max	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and.  p_count .eq.1  )
     &	  chann_indx(i) = st
!      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3)')
!     & "st"," ch",st,i,
!     & "j1",j1_ch(i),"j2",j2_ch(i),"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	   	  
      END SELECT	  
      ENDDO
      states_size = st
      IF(MYID.eq.0) PRINT *, "STATES_SIZE_FOUND",states_size
      IF(MYID.eq.0) PRINT *, "number_of_channels",number_of_channels
      IF(MYID.eq.0) PRINT *, "jmax_included",jmax_included
	  
      ALLOCATE(indx_chann(states_size))
      ALLOCATE(j_min_ind(number_of_channels))
      ALLOCATE(j_max_ind(number_of_channels))	  
      ALLOCATE(m12(states_size))
      ALLOCATE(j12(states_size))
      ALLOCATE(Ej(states_size))
      IF(fine_structure_defined .and. LORB_FINE.gt.0 .and.
     & SPIN_FINE.eq.2 ) THEN
      ALLOCATE(j12_h(states_size))
      ALLOCATE(m12_h(states_size))
      ENDIF	 
      ALLOCATE(indx_corr(jmax_included*4+1,
     & jmax_included*2+1,number_of_channels))
      indx_corr = 0	 
      IF(identical_particles_defined) THEN	  
      ALLOCATE(parity_state_bk(states_size))											!Bikram April 2021
      ALLOCATE(parity_state_sign_bk(states_size))										!Bikram April 2021
      ALLOCATE(parity_state(states_size))
      ALLOCATE(indx_corr_id(2,4*jmax_included+1,
     & jmax_included*2+1,number_of_channels))
      indx_corr_id = 0	 
      ENDIF
      IF(MYID.eq.0) PRINT *, "BASIS ARRAYS CREATED"
	  
! Bikram Start April 2021:
	  IF(identical_particles_defined) then
      SELECT CASE(coll_type)
	  CASE(0)
	  do i = 1, number_of_channels
	  call ASYM_TOP_VECTORS
	  bk_global_parity(i) = parity_inversion(i)
	  call bikram_kappa(ka1_ch(i),kc1_ch(i),ka2_ch(i),kc2_ch(i),
     & kappa1, kappa2)
	  bk_kappa(i) = (-1d0)**(kappa1 + kappa2)
!	  print*, bk_global_parity(i), parity_inversion(i)
	  END do
      END SELECT
	  END IF
! Bikram End.

      indx_chann = 0
      m12 = 0
      j12 = 0	  
c      PRINT*,"states_size",states_size
c      PRINT*,"states",chann_indx
      st = 0	  
      DO i=1,number_of_channels
	  
      SELECT CASE(coll_type)
      CASE(1)
      IF(fine_structure_defined .and. LORB_FINE.gt.0 .and.
     & SPIN_FINE.eq.2 ) THEN
      DO j_count = 1, round(j_h_ch(i)*2) + 1	  
      st = st + 1
      indx_chann(st) = i
      j12_h(st) = j_h_ch(i)
      m12_h(st) = j_count-1 - j_h_ch(i)	
      indx_corr(round(m12_h(st) +j_h_ch(i)+1),round(j_h_ch(i)+0.5d0),i)
     &	  = st	!!!!! STOPED HERE 
      j_min_ind(i) = int(j_h_ch(i))
      j_max_ind(i) = int(j_h_ch(i))	
      Ej(st) = E_ch(i)	  
      ENDDO		 
      ELSE	 
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      indx_chann(st) = i
      j12(st) = j_ch(i)
      m12(st) = j_count
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_ch(i)+1,j_ch(i)+1,i) = st
      ENDDO
      j_min_ind(i) = j_ch(i)
      j_max_ind(i) = j_ch(i)
      ENDIF	  
      CASE(2)
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      indx_chann(st) = i
      j12(st) = j_ch(i)
      m12(st) = j_count
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_ch(i)+1,j_ch(i)+1,i) = st
      ENDDO
      j_min_ind(i) = j_ch(i)
      j_max_ind(i) = j_ch(i)	  
      CASE(3)
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      indx_chann(st) = i
      j12(st) = j_ch(i)
      m12(st) = j_count
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_ch(i)+1,j_ch(i)+1,i) = st
      ENDDO
      j_min_ind(i) = j_ch(i)
      j_max_ind(i) = j_ch(i)		  
      CASE(4)
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      indx_chann(st) = i
      j12(st) = j_ch(i)
      m12(st) = j_count
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_ch(i)+1,j_ch(i)+1,i) = st
      ENDDO
      j_min_ind(i) = j_ch(i)
      j_max_ind(i) = j_ch(i)	  
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(*,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_summ+1,j_summ+1,i) = st	  
      ENDDO
      ENDDO
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)		  
      ELSE
      p_lim_min = 1
      p_lim_max = 2-KRONEKER(j1_ch(i),j2_ch(i)) !min(2-KRONEKER(j1_ch(i),j2_ch(i)),p_lim_max_ini)		  
      DO p_count=p_lim_min,p_lim_max	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
      j12(st) = j_summ
      m12(st) = j_count
	  
      Ej(st) = E_ch(i)
      parity_state(st) = (-1)**(p_count-1)
! Dulat start: correcting the parity for the identical states	  
	  parity_db = p_count																	
	  IF(j1_ch(i).eq.j2_ch(i).and. .not.EVEN_NUM(j_summ)) then
	  parity_state(st) = -1
	  parity_db = 2																			
	  ENDIF
! Dulat END	 
c      IF((j1_ch(i)-j2_ch(i)).eq.0) THEN
c      parity_state(st) = (-1)**j_summ	  
c      ENDIF
      indx_corr_id(parity_db,																! Dulat
     & j_count+j_summ+1,j_summ+1,i) = st
c      PRINT*,st,j12(st),m12(st),      indx_corr_id(p_count,
c     & j_count+j_summ+1,j_summ+1,i),parity_state(st),indx_chann(st),
c     & ej(st)

	 
      ENDDO
      ENDDO
      ENDDO	 
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(1,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_summ+1,j_summ+1,i) = st	  
      ENDDO
      ENDDO
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)		  
      ELSE
      p_lim_max = 2-KRONEKER(j1_ch(i),j2_ch(i))*
     & KRONEKER(v1_ch(i),v2_ch(i))	        
!	 p_lim_max = min(2-KRONEKER(j1_ch(i),j2_ch(i))*
!     & KRONEKER(v1_ch(i),v2_ch(i)),p_lim_max_ini) 	  
      DO p_count=p_lim_min,p_lim_max	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(1,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      parity_state(st) = (-1)**(p_count-1)
! Dulat start: correcting the parity for the identical states	  
	  parity_db = p_count																	
	  IF(j1_ch(i).eq.j2_ch(i).and. .not.EVEN_NUM(j_summ)) then
	  parity_state(st) = -1
	  parity_db = 2																			
	  ENDIF
! Dulat END				
c      IF((j1_ch(i)-j2_ch(i)).eq.0) THEN
c      parity_state(st) = (-1)**j_summ	  
c      ENDIF
      indx_corr_id(parity_db,
     & j_count+j_summ+1,j_summ+1,i) = st	  
      ENDDO
      ENDDO
      ENDDO	 
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)	  
      ENDIF	 
      CASE(7)
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(1,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_summ+1,j_summ+1,i) = st  
      ENDDO
      ENDDO
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)	  
      CASE(8)
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(1,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_summ+1,j_summ+1,i) = st  
      ENDDO
      ENDDO
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)	  
      CASE(9)
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(1,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_summ+1,j_summ+1,i) = st  
      ENDDO
      ENDDO
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)	  
      CASE(0)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ	
      st = st + 1
      indx_chann(st) = i
c      WRITE(1,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_summ+1,j_summ+1,i) = st	  
      ENDDO
      ENDDO
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)		  
	  
      ELSE
!      p_lim_max = min(2-KRONEKER(j1_ch(i),j2_ch(i))
!     & *KRONEKER(ka1_ch(i),ka2_ch(i))*KRONEKER(kc1_ch(i),kc2_ch(i)),
!     & p_lim_max_ini)
	  p_lim_max = 2-KRONEKER(j1_ch(i),j2_ch(i))
     & *KRONEKER(ka1_ch(i),ka2_ch(i))*KRONEKER(kc1_ch(i),kc2_ch(i))
	  
      DO p_count=p_lim_min,p_lim_max
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(1,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      parity_state_bk(st) = (-1)**j_summ*bk_global_parity(i)*bk_kappa(i)
!	  IF(m12(st).eq.0d0) write(*,'(i3,1x,6(i0),3(1x,i0),1x,f12.5)')st, 
!     & j1_ch(i), ka1_ch(i), kc1_ch(i), j2_ch(i), ka2_ch(i), kc2_ch(i), 
!     & j12(st), m12(st), parity_state(st), Ej(st)
      parity_state(st) = (-1)**(p_count-1)
	  
!	  IF(j1_ch(i).eq.j2_ch(i) .and. ka1_ch(i).eq.ka2_ch(i) .and.
!     & kc1_ch(i).eq.kc2_ch(i)) 
!     & parity_state_bk(st) = parity_state_bk(st)*(-1)**j_summ
	  parity_state_sign_bk(st) = parity_state(st)*parity_state_bk(st)
	  
!	  write(*,'(i0,i0,i0,i0,i0,i0,1x,i0)')j1_ch(i), ka1_ch(i),
!     & kc1_ch(i), j2_ch(i), ka2_ch(i), kc2_ch(i), parity_state(st)
	  
c      IF((j1_ch(i)-j2_ch(i)).eq.0) THEN
c      parity_state(st) = (-1)**j_summ	  
c      ENDIF

! Dulat start: correcting the parity for the identical states	  
	  parity_db = p_count																	
	  IF(j1_ch(i).eq.j2_ch(i)) then
	  IF(ka1_ch(i).eq.ka2_ch(i) .and. kc1_ch(i).eq.kc2_ch(i)) then
	  IF(.not.EVEN_NUM(j_summ)) then
	  parity_state(st) = -1
	  parity_db = 2																			
	  ENDIF
	  ENDIF
	  ENDIF
! Dulat END		

      indx_corr_id(parity_db,
     & j_count+j_summ+1,j_summ+1,i) = st	  
      ENDDO
      ENDDO
      ENDDO	 
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)	  
      ENDIF
      END SELECT      
	  
      ENDDO
!      PRINT*, "INITIATION FINISHED"	
	  
c       PRINT*, j_min_ind,j_max_ind   
c       PRINT*,m12,j12
c c      PRINT*,Ej
c       PRINT*,parity_state
c       PRINT*,  indx_corr_id    	  
c      STOP "HERE" 	  
      total_size = 0	  
      DO st1=1,states_size
      DO st2=1,st1
      IF(identical_particles_defined) THEN
      IF(parity_state(st1).ne.parity_state(st2)) CYCLE  
      ENDIF	  
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined) THEN
      IF(int(m12_h(st1)).eq.int(m12_h(st2)).and.
     & m12_h(st2)*m12_h(st1).gt.0 ) total_size = total_size + 1	  
      ELSE	  
      IF(m12(st1).eq.m12(st2)) total_size = total_size + 1
      ENDIF	  
      ENDDO
      ENDDO
!      PRINT*, "TOTAL_SIZE",total_size        !IT WAS CHANGED
!      STOP	  
	  IF(bikram_mij_multiprint .and. matrix_reading_defined .and. 
     & .not. print_matrix_defined) then
! not to allocate and store large arrays for parallel reading
	  st = 0
      DO st1=1,states_size
      DO st2=1,st1
      IF(identical_particles_defined) THEN
      IF(parity_state(st1).ne.parity_state(st2)) CYCLE  
      ENDIF	  	  
      
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined) THEN
      IF(int(m12_h(st1)).eq.int(m12_h(st2)) .and.
     & m12_h(st2)*m12_h(st1).gt.0 )THEN
      st  = st+1	  
      ind_mat(1,st) = st1
      ind_mat(2,st) = st2	  
      ENDIF
	  
      ELSE	  
	  IF(m12(st1).eq.m12(st2)) THEN
      st  = st+1
      ENDIF
      ENDIF	  
      ENDDO
      ENDDO
	  
	  ELSE
      ALLOCATE(ind_mat(2,total_size))
	  
!      IF(MYID.eq.0) PRINT *, "TOTAL_SIZE_NOW",total_size
!      STOP	  
	  st = 0
      DO st1=1,states_size
      DO st2=1,st1
      IF(identical_particles_defined) THEN
      IF(parity_state(st1).ne.parity_state(st2)) CYCLE  
      ENDIF	  	  
      
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined) THEN
      IF(int(m12_h(st1)).eq.int(m12_h(st2)) .and.
     & m12_h(st2)*m12_h(st1).gt.0 )THEN
      st  = st+1	  
      ind_mat(1,st) = st1
      ind_mat(2,st) = st2	  
      ENDIF	  
	  
      ELSE	  
	  IF(m12(st1).eq.m12(st2)) THEN
      st  = st+1	  
      ind_mat(1,st) = st1
      ind_mat(2,st) = st2	  
      ENDIF
      ENDIF	  
      ENDDO
!	  IF(myid.eq.0) print*,"progress", st1
      ENDDO
!	  PRINT *, myid, total_size, st
!	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
	  END IF
	  
!      IF(st.ne.total_size) STOP"ERROR: INI FAILED"
      IF(st.ne.total_size) st=total_size
	  IF(.not.bikram_mij_multiprint) then
      IF(mpi_task_defined) THEN
      IF(MYID.EQ.0) THEN	  
      ALLOCATE(Mat_el(n_r_coll,total_size))
      ALLOCATE(Mat_el_der(n_r_coll,total_size))
      ENDIF	  
      ELSE
      ALLOCATE(Mat_el(n_r_coll,total_size))
      ALLOCATE(Mat_el_der(n_r_coll,total_size))	  
      ENDIF	  
	  END IF
      ALLOCATE(R_COM(n_r_coll))
c      IF(identical_particles_defined) THEN	  
c      DO st=1,states_size
c      i = indx_chann(st)   
c      WRITE(*,'(i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
c    & i,st,j1_ch(i),j2_ch(i),j12(st),m12(st),
c     &l_state(st),parity_state(st)	  
c      ENDDO
c      ENDIF	  
c     STOP	
      IF(MYID.eq.0 .and. print_states_defined) THEN  !!! COMMENT FOR A WHILE
      PRINT*, "ARRAYS_INITIALIZED"
      OPEN(4,FILE="STATES.out")

      SELECT CASE(coll_type)
      CASE(1)
      IF(SPIN_FINE.eq.2 .and. LORB_FINE.gt.0) THEN
      WRITE(4,'(a8,1x,a8,1x,a4,1x,a5,1x,1x,a4)')
     & "#STATE","#CHANNEL","J12","M12","J"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,f4.1,1x,f5.1,1x,1x,f4.1)')
     & st,i,j12_h(st),m12_h(st),j_h_ch(i)	  
      ENDDO	  
      ELSE	  
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,1x,i3)')
     & st,i,j12(st),m12(st),j_ch(i)	  
      ENDDO
      ENDIF	  
      CASE(2)
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","V","J"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),v_ch(i),j_ch(i)	  
      ENDDO
      CASE(3)
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J","K","EPS"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j_ch(i),k_ch(i),eps_ch(i)	  
      ENDDO
      CASE(4)
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J","KA","KC"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j_ch(i),ka_ch(i),kc_ch(i)	  
      ENDDO	  
      CASE(5)
      IF(.NOT.identical_particles_defined) THEN 	  
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","J2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i),parity_state(st)
      ENDDO	 
      ENDIF
      CASE(6)
      IF(.NOT.identical_particles_defined) THEN 	  
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","V1","J1","V2","J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),v1_ch(i),j1_ch(i),v2_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","V1","J1","V2","J2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),v1_ch(i),j1_ch(i),v2_ch(i),j2_ch(i),
     & parity_state(st)	 
      ENDDO 
      ENDIF
      CASE(7)
      WRITE(4,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","K1","EPS", "J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),k1_ch(i),eps1_ch(i),j2_ch(i)	  
      ENDDO	  
      CASE(8)
      WRITE(4,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","KA1","KC1","J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i),j2_ch(i)	  
      ENDDO	  
      CASE(9)
      WRITE(4,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","KA1","KC1",
     & "J2","K2","EPS"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),k2_ch(i),eps2_ch(i)	 
      ENDDO	  
      CASE(0)
      IF(.NOT.identical_particles_defined) THEN  
      WRITE(4,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","KA1","KC1",
     & "J2","KA2","KC2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i)	 
      ENDDO
      ELSE
      WRITE(4,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a2)
     & ')
     & "#STATE","#CHANNEL","J12","M12","J1","KA1","KC1",
     & "J1","KA2","KC2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i2)
     & ')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i),parity_state(st)
      ENDDO	 
      ENDIF	 
      END SELECT
      CLOSE(4)
      OPEN(6,FILE="STATE_Mij_INDX.out")
      WRITE(6,"(a9,1x,a4,1x,a4)") "Mij_INDEX", "ST_1", "ST_2"
      DO i=1,total_size
      WRITE(6,"(i9,1x,i4,1x,i4)") i,ind_mat(1,i),ind_mat(2,i)	  
      ENDDO	  
      CLOSE(6)
      ENDIF
!      PRINT *, "barrier", myid
!	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
!	  IF(MYID.eq.0) PRINT *, "done"
      END SUBROUTINE INI_ARRAYS 
      SUBROUTINE MATRIX_MAKE
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      USE MPI_TASK_TRAJECT	  
      USE POT_STORE
      USE CS_MATRIX	 	  
      USE OLD_MIJ	  
      IMPLICIT NONE
      LOGICAL file_exst
      CHARACTER(LEN=22) :: R_GRID_DAT = "USER_DEFINED_RGRID.DAT"	  
      INTEGER ISTAT	  
      INTEGER i,st_1,st_2,k,k_st_mpi,k_fn_mpi,chunk_mpi_size,task_size
      INTEGER L_NJU(8)
      INTEGER mean_size_1,mean_size_2,i_nr_fin,i_nr_ini
      INTEGER status(MPI_STATUS_SIZE)	  
      INTEGER k_check,i_check,n_point_dis_1,n_point_dis_2	  
      REAL*8 intgeral,R_dist,dR_step
      REAL*8 TIME_WRITING_MATRIX
	  REAL*8 bk_time_bgn, bk_time_fin,bk_time															!Bikram Jan,'19	  
      REAL*8, ALLOCATABLE :: Mat_el_temp(:,:),Mat_el_der_temp(:,:),
     & Mat_el_res(:),Mat_el_der_res(:),Mat_el_resf(:,:)
     &,Mat_el_der_resf(:,:)
      REAL*8,ALLOCATABLE :: Mat_el_r_temp(:,:), Mat_el_der_r_temp(:,:)
	  REAL*8, allocatable :: bk_mat_temp1(:,:)															!Bikram Nov 2020
	  REAL*8, allocatable :: deviation_r(:)																!Bikram Aug 2022
	  REAL*8 bk_mat_temp(n_r_coll)																		!Bikram Nov 2020
	  INTEGER mij_remander, bk_mij_addition, bk_mat_counter, mt_chk										!Bikram Dec 2020
	  logical mij_k_skip																				!Bikram Dec 2020
	  REAL*8 bgn_tym, END_tym, calc_tym																	!Bikram Dec 2020
	  REAL*8 tot_tym_calc, tot_tym_prn, tot_tym_rd														!Bikram Dec 2020
	  INTEGER mij_counter,bk_non_zero_counter,bk_non_zero_counter_old
	  CHARACTER (LEN=28) :: dm5
	  CHARACTER (LEN=21) :: dm6
	  CHARACTER (LEN=9) :: dm7
	  CHARACTER (LEN=500) :: bk_matrix_path1, bk_matrix_path2
	  CHARACTER (LEN=500) :: bk_matrix_path3, bk_matrix_path4
	  CHARACTER (LEN=500) :: bk_matrix_path5
	  REAL*8, allocatable :: bk_mat_array(:)
	  INTEGER file_counter,fc_old1,k_st,k_fn,ii,dmm1,dmm2,dmm3(2),dmm4
	  INTEGER dmm5(8), proc_cntr, MPI_Req(nproc), prg_cntr, prct_cntr_rcv
	  INTEGER percent_counter, prcnt, prct_cntr(nproc), cyc
	  INTEGER, allocatable :: bk_ind_tmp(:,:)
	  INTEGER, allocatable :: file_old(:,:), file_old1(:,:), cyc_cntr(:)
	  REAL*8 MIJ_ZERO_CUT_old, mtrx_cutoff_chk1, mtrx_cutoff_chk2
	  REAL*8 bk_tym1, bk_tym2, bk_tym
	  INTEGER load_cntr, dm11, dm_bgn, dm_END, k_rstrt_chk, k_start
	  CHARACTER (LEN=500) :: load_file
	  REAL*8, allocatable :: load_mij(:)
	  logical cut_r, mpi_test_flag
      LOGICAL, ALLOCATABLE :: K_SKIPPED_BY_ROUTINE(:)
      INTEGER, allocatable :: bk_k_gather(:)
		REAL*8, allocatable :: Mat_el_R_array(:),exp_result_array(:)
		INTEGER :: local_k, local_cyc
		INTEGER :: local_n_r_coll
		
		REAL*8 :: stack_result_array(10)  ! Fixed size on stack
		INTEGER :: copy_size
		
      LOGICAL :: K_SKIPPED_RES = .FALSE.
      LOGICAL BELONGS, bk_mtrx_re
      EXTERNAL BELONGS
		
		INTERFACE
        SUBROUTINE EXPANSION_MATRIX_ELEMENT(M_coulp_array,k,tmp1)
          USE VARIABLES
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: k
          REAL*8, INTENT(OUT) :: M_coulp_array(*)
			 INTEGER:: tmp1
        END SUBROUTINE
      END INTERFACE
		
	  cut_r = .false.
	  
! Prints timing information for array allocation												   
      IF(myid.eq.0) WRITE(*,'(a32,1x,f17.3)')
     & "ARRAYS ALLOCATION TOOK TIME, s =",
     & TIME_1_ARRAY	  
! Validates matrix size (total_size > 0)										   
      IF(myid.eq.0) PRINT*,"TOTAL_SIZE_OF_M = ",total_size
      IF(total_size.le.0) STOP "ERROR:ZERO DIM Mij MATRIX"
! Checks IF matrix calculation is properly defined													
      IF(.not. calc_matrix_defined) THEN
      IF(MYID.eq.0) PRINT*,"PROGRAM STOPED"	  
      CALL MPI_FINALIZE(ierr_mpi)
      STOP      	  
      ENDIF	  

      conv_unit_r = 1d0
      conv_unit_e = 1d0
      mean_size_1 = total_size
      mean_size_2 = 0
      i_nr_ini = max(ir_bgn_exp_pnt,1)
      i_nr_fin = max(min(ir_fin_exp_pnt,n_r_coll),i_nr_ini)
! Determines the R-grid range for expansion points													
      IF(ir_fin_exp_pnt.lt.0 .or. ir_bgn_exp_pnt.lt.0) THEN
      i_nr_ini = 1
      i_nr_fin = n_r_coll	  
      ENDIF
      IF(i_nr_ini.gt.n_r_coll .or. i_nr_fin.lt.1) THEN
      IF(myid.eq.0) THEN
      PRINT*,
     & "ERROR: DEFINE PROPERLY GRID RANGE"
      PRINT*,"i_nr_ini",i_nr_ini
      PRINT*,"n_r_coll",n_r_coll	  
      ENDIF	
      STOP	  
      ENDIF	 	  
!Sets up unit conversions for energy and distance
      IF(angs_unit) conv_unit_r = a_bohr
      IF(cm_unit)  conv_unit_e =1d0/eVtown/autoeV 	  
      IF(klv_unit)	conv_unit_e =1d0/autokelv/autoeV
      IF(kal_unit) conv_unit_e = 1d0/kcal	  
      massAB_M = mass_red*amutoau
      Ej = Ej/eVtown/autoeV
     	  
      IF(expansion_defined) THEN	  
      IF(terms_file_defined ) THEN
      pot_expansion_defined = .FALSE.
      ELSE
      pot_expansion_defined = .TRUE.	  
      ENDIF
      IF(terms_onfly_defined) THEN
      expansion_defined = .FALSE.
      pot_expansion_defined = .TRUE.	  
      ENDIF	  
      ENDIF
	  
      IF(.not.expansion_defined .and. pot_expansion_defined .and.
     & .not. mpi_task_defined ) THEN
      DEALLOCATE(Mat_el,Mat_el_der)
      CALL EXCLUDE_STATES	  
      RETURN
      ELSE
      IF(pot_expansion_defined)expansion_defined = .TRUE.	  
      ENDIF	  
      IF(expansion_defined) THEN
      make_grid_file = .FALSE.	  
      ENDIF	  
      ALLOCATE(Ks_have_to_be_skipped(total_size))	  
      IF(n_r_coll.le.1) STOP "ERROR: TOO FEW POINTS FOR COLL GRID"
! EXPANSION NOT DEFINED,	  
      IF(.not. expansion_defined .or. pot_expansion_defined .or.
     & calc_expansion_defined) THEN
! Handles three scenarios for creating the radial distance grid:																
      IF(.not.read_r_grid_from_file_defined) THEN	 
      IF(.not.eq_grid_defined) THEN	  
! 1. Creates a non-unIForm "optimized" grid with three sections -  Dense points near min distance, regular spacing in the middle and adjusted spacing towards maximum																																									 
      IF(MYID.eq.0) PRINT*," MAKING AN ""OPTIMIZED"" GRID	"  
      n_point_dis_1 = int(n_r_coll/3)
      n_point_dis_2 = 	int(n_r_coll/2)+1  
      R_COM(1) = R_min_dist
      R_COM(n_r_coll) = R_max_dist	  
      dR_step = (R_max_dist-R_min_dist)/(n_r_coll-1d0)	  
      DO i=1,n_point_dis_1
      R_COM(i) = R_min_dist+ dR_step*(i-1d0)/3d0    
      ENDDO
      DO i=n_point_dis_1+1,n_point_dis_2
      R_COM(i) = R_COM(n_point_dis_1)+ dR_step*(i-n_point_dis_1)    
      ENDDO
      dR_step = (R_max_dist-R_COM(n_point_dis_2))/
     & (n_r_coll-n_point_dis_2)
      DO i=n_point_dis_2,n_r_coll
      R_COM(i) = R_COM(n_point_dis_2)+ dR_step*(i-n_point_dis_2)    
      ENDDO
      ELSE
! 2. Making an equidistant grid between R_min and R_max														
      IF(MYID.eq.0) PRINT*," MAKING AN EQUIDIST GRID	" 	  
      dR_step = (R_max_dist-R_min_dist)/(n_r_coll-1d0)	  
      DO i=1,n_r_coll
      R_COM(i) = R_min_dist+ dR_step*(i-1d0)   
      ENDDO	  
      ENDIF	  
      ENDIF

! 3. Reads R-grid from "USER_DEFINED_RGRID.DAT" and broadcasts grid to all MPI processess																						 
	  IF(read_r_grid_from_file_defined) THEN
      INQUIRE( FILE=R_GRID_DAT, EXIST=file_exst)
      IF(file_exst) THEN	  
      IF(MYID.eq.0) THEN
      OPEN(23,FILE=R_GRID_DAT,STATUS="OLD",ACTION="READ")
      READ(23,*,IOSTAT=ISTAT)	  
      DO i=1,n_r_coll
      READ(23,*,IOSTAT=ISTAT) i_old,R_COM(i)
      ENDDO
      R_COM = R_COM!/conv_unit_r	  				!This conversion has been disable because the file is always in Bohr unit
      CLOSE(23)
      IF(ISTAT.gt.0) STOP "ERROR:USER DEFINED R_GRID FILE NOT FOUND"
      ENDIF	
      CALL MPI_BCAST(R_COM, n_r_coll, MPI_DOUBLE_PRECISION,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      ELSE
      IF(MYID.eq.0) PRINT*," R_GRID FILE NOT FOUND" 	  
      IF(MYID.eq.0) PRINT*," MAKING AN ""OPTIMIZED"" GRID	"  
      n_point_dis_1 = int(n_r_coll/3)
      n_point_dis_2 = 	int(n_r_coll/2)+1  
      R_COM(1) = R_min_dist
      R_COM(n_r_coll) = R_max_dist	  
      dR_step = (R_max_dist-R_min_dist)/(n_r_coll-1d0)	  
      DO i=1,n_point_dis_1
      R_COM(i) = R_min_dist+ dR_step*(i-1d0)/3d0    
      ENDDO
      DO i=n_point_dis_1+1,n_point_dis_2
      R_COM(i) = R_COM(n_point_dis_1)+ dR_step*(i-n_point_dis_1)    
      ENDDO
      dR_step = (R_max_dist-R_COM(n_point_dis_2))/
     & (n_r_coll-n_point_dis_2)
      DO i=n_point_dis_2,n_r_coll
      R_COM(i) = R_COM(n_point_dis_2)+ dR_step*(i-n_point_dis_2)    
      ENDDO	  
      ENDIF	  
	  
      ELSE
      INQUIRE( FILE=R_GRID_DAT, EXIST=file_exst )	  
      IF (MYID.EQ.0 .and. .not.file_exst) THEN
!     IF GRID DOES NOT EXIST CREATE ITS OWN R-GRID FILE	  
      OPEN(23,FILE=R_GRID_DAT,STATUS="NEW",ACTION="WRITE")	  
      WRITE(23,'(a3,1x,a12)') '#iR','R_COM(iR)'	  
      DO i=1,n_r_coll
      WRITE(23,'(i3,1x,e19.12)') i,R_COM(i)
!      IF(i.ne.i_old) STOP "ERROR IN R_GRID_FILE"	  
      ENDDO
      CLOSE(23)	  
      ENDIF	  
      ENDIF	 
	  
! finding R for computing RMS
	  allocate(deviation_r(n_r_coll))
	  do i = 1, n_r_coll
	  deviation_r(i) = abs(R_COM(i) - bikram_rms_r)
	  END do
	  rms_r = minloc(deviation_r,1)
	  deallocate(deviation_r)
      ENDIF
	  
!     CREATING/READING V-POTENTIAL FILE ON THE TOTAL ANGLE(S)/R-GRID 	  
      IF(make_grid_file .and. (coll_type.gt.5 .or. coll_type.le.1))
     & THEN	   
! 	Checks for existing potential grid files											 
      INQUIRE(FILE=potential_file_name,EXIST=grid_file_found)       
      IF(grid_file_found) THEN
      IF(MYID.eq.0)PRINT*,"GRID WILL BE READ FROM FILE" 
      ELSE
      IF(MYID.eq.0)PRINT*,
     & "GRID WILL BE STORED IN FILE FOR FURTHER USE"
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi ) 	   	   
! If not found, calls INI_POT_STORAGE to create them														
      CALL INI_POT_STORAGE
!!    CREATING A GRID	  
      IF(MYID.eq.0) PRINT*,"GRID HAD BEEN STORED IN FILE"
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi ) 	   
      ENDIF	   
      ENDIF
!      STOP

!     COMPUTING POTENTIAL EXPANSION	  
! If expansion calculation is requested									   
      IF(calc_expansion_defined) THEN
	  bk_time_bgn=MPI_Wtime()		!Bikram Jan,'19
      IF(MYID.EQ.0) WRITE(*,*)"EXPANSION WILL BE COMPUTED"
! Sets up angular momentum quantum numbers										  
      L_NJU  = 0 	  
      L_NJU(1) = L_MAX_EXPAN
      L_NJU(2) = L1_MAX_EXPAN
      L_NJU(3) = L2_MAX_EXPAN
      L_NJU(4) = NJU_MAX_EXPAN
      L_NJU(5) = NJU1_MAX_EXPAN
      L_NJU(6) = NJU2_MAX_EXPAN
      L_NJU(7) = L1_MIN_EXPAN
      L_NJU(8) = L2_MIN_EXPAN	  
! Calls CALC_EXPANSION to compute potential expansion coefficients																	
      CALL CALC_EXPANSION(coll_type,n_r_coll,L_NJU,
     & bikram_ident_terms,bikram_identical_pes, 
     & bikram_axl_sym1, bikram_axl_sym2, 
     & bikram_equ_sym1, bikram_equ_sym2)
! Times the calculation and exits
      bk_time_fin=MPI_Wtime()		!Bikram Jan,'19
	  bk_time=bk_time_fin-bk_time_bgn		!Bikram Jan,'19
	  IF(myid.eq.0)print*,"TIME TOOK TO COMPUTE EXPANSION TERMS = ",
     & bk_time,"SEC"	 !Bikram Jan,'19
!	  IF (myid.eq.0) then
!	  INQUIRE(FILE="PROGRESS_EXPAN.tmp", EXIST=file_exst )	
!      IF(file_exst) open(1111,FILE="PROGRESS_EXPAN.tmp")
!      close(1111,status='delete')
!	  ENDIF
      IF(MYID.eq.0)PRINT*,"EXPANSION HAS BEEN COMPUTED,PROGRAM HAS ",
     & "BEEN STOPED"
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi ) 	   
      STOP	  
      ENDIF
! ASSIGNING THE SIZE OF A CHUNK
! Divides the total matrix work among MPI processes
! Each process gets a chunk: chunk_mpi_size = total_size/nproc
      chunk_mpi_size = total_size/nproc
! Calculates start (k_st_mpi) and END (k_fn_mpi) indices for each process																		 
      k_st_mpi = myid*chunk_mpi_size + 1
      k_fn_mpi = (myid+1)*chunk_mpi_size 	  
!      IF(myid.eq.nproc-1) k_fn_mpi = total_size
      task_size = (k_fn_mpi-k_st_mpi +1)*n_r_coll
      IF(MYID.eq.0 .and. .not.matrix_reading_defined
     & .and. chunk_mpi_size.gt.0 )
     & WRITE(*,'(a47,1x,i0)')
     & "THE SIZE OF THE Mij COMPUTED BY EACH PROCESSOR=",
     & chunk_mpi_size
      IF(chunk_mpi_size.eq.0 .and. myid.eq.0 ) WRITE(*,*)
     & "NUMBER OF CPUS > SIZE OF Mij"	  
       ! MATRIX READING OR CALCULATINGS 
      TIME_MAT_START = MPI_Wtime()	   
		
! Matrix Computing (matrix_reading_defined = .FALSE.)
      IF(.not.matrix_reading_defined) THEN
	  
!	  tot_tym_calc = 0.d0											!Bikram
!	  tot_tym_prn = 0.d0											!Bikram
!	  tot_tym_rd = 0.d0												!Bikram

!!-------------------------------------------------------------------------------------------------	  
!! Bikram Start Dec 2020:
	  	  
! Handles parallel matrix computation with file output																
	  IF (print_matrix_defined .and. bikram_mij_multiprint) then
	  IF(allocated(Mat_el)) deallocate(Mat_el)
	  IF(allocated(Mat_el_der)) deallocate(Mat_el_der)
	  IF(.not.check_point_defined) then
	  IF(myid.eq.0) call system ( "rm -r " // trim(bk_dir1) )
	  IF(myid.eq.0) call system ( "rm -r " // trim(bk_dir2) )
	  IF(myid.eq.0) call system ( "mkdir -p " // trim(bk_dir2) )
	  END IF
	  IF(bikram_rebalance) then
	  IF(myid.eq.0) call system ( "rm -r " // trim(bk_dir33) )
	  IF(myid.eq.0) call system ( "mkdir -p " // trim(bk_dir33) )
	  END IF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  
	  IF(total_size.ge.nproc) then
	  
      chunk_mpi_size = total_size/nproc
      k_st_mpi = myid*chunk_mpi_size + 1
      k_fn_mpi = (myid+1)*chunk_mpi_size 	  	 
	  mij_remander = 0
	  IF(total_size.gt.(chunk_mpi_size*nproc)) 
     & mij_remander = total_size - chunk_mpi_size*nproc
	  IF(myid.lt.mij_remander) then	  
      k_st_mpi = k_st_mpi + myid
      k_fn_mpi = k_fn_mpi + (myid + 1)
	  ELSE	  
      k_st_mpi = k_st_mpi + mij_remander
      k_fn_mpi = k_fn_mpi + mij_remander
	  END IF
	  
!!!!  COMPUTING MATRIX	 
      IF(MYID.EQ.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "COMPUTING MATRIX ELEMENTS STARTED"	  
! Computes matrix elements in chunks of 1000	 
	  bk_mat_counter = 0
	  bk_non_zero_counter = 0
	  allocate(bk_mat_temp1(n_r_coll,1000))
	  allocate(bk_non_zero_mij_gather(nproc))
	  IF(write_check_file_defined .or. check_point_defined) 
     & allocate(bk_k_gather(nproc))
	  IF(bikram_rebalance .or. bikram_rebalance_comp) 
     & allocate(cyc_cntr(1000))
	  k_rstrt_chk = 0
	  
! finding #r for matrix truncation
	  IF(.not.cut_r) then
      CALL	INTEGRATOR_TEMP(intgeral,k_st_mpi,1,cyc) 
! It finds two truncation points (mtrx_cutoff_r1 and mtrx_cutoff_r2)
! Finding the closest grid points to bikram_cutoff_r1 and bikram_cutoff_r2
! Using minloc to find indices where R_COM values are closest to these cutoff values
	  allocate(deviation_r(n_r_coll))
	  do ii = 1, n_r_coll
! bikram_cutoff_r1 = 6.0d0 &  bikram_cutoff_r2 = 8.0d0															  
	  deviation_r(ii) = abs(R_COM(ii) - bikram_cutoff_r1)
	  END do
	  mtrx_cutoff_r1 = minloc(deviation_r,1)
	  deallocate(deviation_r)
	  allocate(deviation_r(n_r_coll))
	  do ii = 1, n_r_coll
	  deviation_r(ii) = abs(R_COM(ii) - bikram_cutoff_r2)
	  END do
	  mtrx_cutoff_r2 = minloc(deviation_r,1)
	  deallocate(deviation_r)
	  cut_r = .true.
	  IF(myid.eq.0) write(*,'(2(a,i0,a,f0.3))') 
     & "Truncation #1 at #R = ", mtrx_cutoff_r1, ", R = ", 
     & R_COM(mtrx_cutoff_r1), ", and #2 at #R = ",mtrx_cutoff_r2,  
     & ", R = ",R_COM(mtrx_cutoff_r2)
	  MIJ_ZERO_CUT = MIJ_ZERO_CUT/eVtown/autoeV
	  END IF

!----------------------------------------------------------------------------------------
! Bikram: this piece is to read those matrix elements to be computed	
!----------------------------------------------------------------------------------------
	  IF(bikram_rebalance_comp) then
	  open(1, file = trim(bk_dir33)//"Load.dat", status = 'old', 
     & action = 'read')
	  dm11 = 0
	  do k = 1, myid+1
	  read(1,*) load_cntr
	  dm11 = dm11 + 1
	  END do
	  close(1)
	  IF(dm11-1 /= myid) stop "correct load not read"
	  allocate(load_mij(load_cntr))
	  write(load_file,"(a,a,i0,a)")trim(bk_dir33),"Proc_",myid,".DAT"
	  load_file = trim(load_file)
	  open(1, file = trim(load_file), status = 'old', action = 'read')
	  do dm11 = 1, load_cntr
	  read(1,*) load_mij(dm11)
	  IF(load_mij(dm11) == 0) stop "Error laoding distributed k"
	  END do
	  close(1)
	  
	  TIME_1_MAT = MPI_Wtime()
	  IF(myid.eq.0) call bk_print_matrix_info
	  TIME_2_MAT = MPI_Wtime()
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  
	  bk_tym1 = MPI_Wtime()
	  percent_counter = 1
	  dm_bgn = load_mij(1)
	  dm_END = load_mij(load_cntr)
	  
	  IF(check_point_defined) then
	  IF(myid.eq.0) print*, "Matrix caculations restarting."
	  open(11,file="Restart_info.DAT", action = 'read')
	  do i = 1, myid+1
	  read(11,*) k_rstrt_chk
	  END do
	  close(11)
	  k_st_mpi = k_rstrt_chk
	  IF(k_st_mpi == 0) k_st_mpi = 1
	  ELSE
	  k_st_mpi = 1
	  END IF

! =====================================================================
! PARALLEL MATRIX COMPUTATION LOOP - LOAD BALANCED VERSION
! =====================================================================
! This loop processes matrix elements that have been distributed across
! MPI processes for load balancing. Each process computes a subset of
! the total matrix elements.
	  
	  DO  dm11 = k_st_mpi, load_cntr
! ------------------------------------------------------------------
! LOAD BALANCING INDEX MAPPING
! ------------------------------------------------------------------
! dm11: Local loop counter (ranges from k_st_mpi to load_cntr)
! load_mij(dm11): Contains the actual global matrix element index
!                 that this process should compute
	  k = load_mij(dm11)
      bk_mat_counter = bk_mat_counter + 1
! ! Store global matrix index for file writing	  
	   cyc_cntr(bk_mat_counter) = k

! Print progress of matrix computation	 
	  IF(mod((dm11 - k_st_mpi), int((load_cntr - k_st_mpi)/10)) == 0) then
	  bk_tym2 = MPI_Wtime()
	  bk_tym = bk_tym2 - bk_tym1
	  write(*,'(2(a, i5),a,f12.3)') "proc_id = ",myid,", Progress = ",
     & int(dble(dm11 - k_st_mpi)/dble(load_cntr - k_st_mpi)*100.d0), 
     & "%, Time(sec.) = ", bk_tym
	  percent_counter = percent_counter + 1
	  END IF
	  IF(dm11 == load_cntr) then
	  bk_tym2 = MPI_Wtime()
	  bk_tym = bk_tym2 - bk_tym1
	  write(*,'(2(a, i5),a,f12.3)') "proc_id = ",myid,", Progress = ",
     & int(100.d0), "%, Time(sec.) = ", bk_tym
	  END IF
! Bikram End.

! Matrix computation 	  	  
		IF (.NOT. expansion_defined) THEN
! Direct integration case: loop over R points	  	  
      DO i=1,n_r_coll
      CALL	INTEGRATOR_TEMP(intgeral,k,i,cyc)
      bk_mat_temp1(i,bk_mat_counter) = intgeral*conv_unit_e
      ENDDO
		ELSE
! Expansion-based computation
!		write(*,*) "#0.1"
      ALLOCATE(Mat_el_R_array(n_r_coll))
      CALL EXPANSION_MATRIX_ELEMENT(Mat_el_R_array,k,cyc)
      bk_mat_temp1(:,bk_mat_counter) = Mat_el_R_array(:)*conv_unit_e
      DEALLOCATE(Mat_el_R_array)
      ENDIF
		
!	   cyc_cntr(bk_mat_counter) = k
! File writing when buffer is full or at END
	   bgn_tym = MPI_Wtime()
	   IF(bk_mat_counter.eq.1000 .or. dm11.eq.load_cntr) then
	   IF(write_check_file_defined) then
      IF(abs(bgn_tym-bk_tym1) < TIME_MIN_CHECK*60.d0) then
	   call bk_print_matrix (dm_bgn, dm_END, k,
     & bk_mat_counter, bk_mat_temp1, mt_chk, cyc_cntr)
	   bk_mat_counter = 0
	   k_rstrt_chk = dm11 + 1
	   ELSE
	   exit
	   ENDIF
	  
	   ELSE
	   call bk_print_matrix (dm_bgn, dm_END, k,
     & bk_mat_counter, bk_mat_temp1, mt_chk, cyc_cntr) 
	   bk_mat_counter = 0
	   END IF
	   END IF

	  IF(write_check_file_defined .and. 
     & abs(bgn_tym-bk_tym1) > TIME_MIN_CHECK*60.d0) exit
      ENDDO
	  
	  IF(write_check_file_defined .and. k_rstrt_chk < load_cntr) then
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  CALL MPI_GATHER(k_rstrt_chk,1,MPI_INTEGER,bk_k_gather,1,
     & MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
	  IF(myid.eq.0) then
	  open(11,file="Restart_info.DAT")
	  do i = 1, nproc
	  write(11,*) bk_k_gather(i)
	  END do
	  close(11)
	  write(*,'(a,a)')
     & "Matrix calculations did not finish due to time limit. ",
     & "Please restart. Program will stop now."
	  END IF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      STOP
	  END IF
	  
	  bk_tym2 = MPI_Wtime()
	  bk_tym = bk_tym2 - bk_tym1
	  open(111, file = "Time_per_proc.out", access = 'appEND')
	  write(111,'(a,i5,a,f12.5)')"#proc = ",myid,", time = ",bk_tym
	  close(111)
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )

	  IF(myid.eq.0) then
	  open(1, file = trim(bk_dir33)//"Load.dat", status = 'old', 
     & action = 'read')
	  dm11 = 0
	  do i = 1, nproc
	  read(1,*) load_cntr
	  dm11 = dm11 + load_cntr
	  END do
	  close(1)
	  write(bk_matrix_path5, '(a,a)') trim(bk_dir2), 
     & "/MATRIX_NONZERO_INDEX.DAT"
	  open(11,file=trim(bk_matrix_path5))
	  write(11,'(a9,i16)') "#Files = ", nproc
	  write(11,'(a28,i16)') "Non-Zero #Matrix Elements = ", dm11
	  write(11,'(a21,e19.12,a6)') "The Matrix Cut-Off = ", 
     & MIJ_ZERO_CUT*eVtown*autoeV, " cm^-1"
	  write(11,'(a16,2x,a16)') '          file_#', '   #non_zero_Mij'
	  
	  open(1, file = trim(bk_dir33)//"Load.dat", status = 'old', 
     & action = 'read')
	  do i = 1, nproc
	  read(1,*) load_cntr
	  write(11,'(i16,2x,i16)') i, load_cntr
	  END do
	  close(11)
	  write(bk_matrix_path5, '(4(a))') trim(bk_dir33), 
     & "/MATRIX_FILE_INDEX_NEW.DAT ", trim(bk_dir2), 
     & "/MATRIX_FILE_INDEX.DAT"
	  call system ( "cp " // trim(bk_matrix_path5) )
	  END IF
	  
      TIME_WRITING_MATRIX = TIME_2_MAT - TIME_1_MAT
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      TIME_MAT_FINISH= MPI_Wtime()
      IF(MYID.EQ.0) WRITE(*,'(a57,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij READING/COMPUTING/SAVING ON DISK"
     & ,",s",(TIME_MAT_FINISH-TIME_MAT_START)/nproc
	  IF(print_matrix_defined.and. myid.eq.0) 
     & WRITE(*,'(a47,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij SAVING ON DISK"
     & ,",s",TIME_WRITING_MATRIX	 
      IF(MYID.EQ.0) PRINT*, "ALL WORK ON MATRIX IS DONE"
      IF(run_prog_defined) THEN
	  IF(myid.eq.0) print*, "You indicated Parallel Matrix printing. "
     & , "Start trajectory calculation separately. ",
     & "Program will stop now."
      STOP	  
	  ELSE
	  stop
      ENDIF

	  return
	  END IF
	  
	  TIME_1_MAT = MPI_Wtime()
	  IF(myid.eq.0) call bk_print_matrix_info
	  TIME_2_MAT = MPI_Wtime()
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  
	  bk_tym1 = MPI_Wtime()
! Bikram Start: this is to print progress of matrix computation
	  percent_counter = 1
	  k_start = k_st_mpi
	  IF(check_point_defined) then
	  IF(myid.eq.0) print*, "Matrix caculations restarting."
	  open(11,file="Restart_info.DAT", action = 'read')
	  do i = 1, myid+1
	  read(11,*) k_rstrt_chk
	  END do
	  close(11)
	  k_start = k_rstrt_chk
	  IF(k_start == 0) k_start = k_st_mpi
	  END IF

      DO  k = k_start, k_fn_mpi
	  
	  IF((k_fn_mpi - k_start) .ge. 10) then
	  IF(mod((k - k_start), int((k_fn_mpi - k_start)/10)) == 0) then
	  bk_tym2 = MPI_Wtime()
	  bk_tym = bk_tym2 - bk_tym1
	  write(*,'(2(a, i5),a,f12.3)') "proc_id = ",myid,", Progress = ",
     & int(dble(k - k_start)/dble(k_fn_mpi - k_start)*100.d0), 
     & "%, Time(sec.) = ", bk_tym
	  percent_counter = percent_counter + 1
	  END IF
	  END IF
	  IF(k == k_fn_mpi) then
	  bk_tym2 = MPI_Wtime()
	  bk_tym = bk_tym2 - bk_tym1
	  write(*,'(2(a, i5),a,f12.3)') "proc_id = ",myid,", Progress = ",
     & int(100.d0), "%, Time(sec.) = ", bk_tym	  
	  END IF
! Bikram End.
	  
	  bk_mat_counter = bk_mat_counter + 1
	  IF(bk_mat_counter.eq.1) mt_chk = k
	  mij_k_skip = .false.
	  
! calling matrix calculation for 2 values of R to check 
! whether to compute the entrire array along R
! Skip this IF its expansion 

      IF (.NOT. expansion_defined) THEN
      CALL	INTEGRATOR_TEMP(intgeral,k,mtrx_cutoff_r1,cyc)
	  mtrx_cutoff_chk1 = intgeral*conv_unit_e
	  IF(bikram_rebalance) cyc_cntr(bk_mat_counter) = cyc
      CALL	INTEGRATOR_TEMP(intgeral,k,mtrx_cutoff_r2,cyc)
	  mtrx_cutoff_chk2 = intgeral*conv_unit_e
	  IF(max(abs(mtrx_cutoff_chk1), abs(mtrx_cutoff_chk2)).lt.
     & MIJ_ZERO_CUT) then
      bk_mat_temp1(:,bk_mat_counter) = 0d0
      mij_k_skip = .TRUE.
	   bk_non_zero_counter = bk_non_zero_counter + 1
      ENDIF	  
      ENDIF
	  
      IF(.not.mij_k_skip) then 	  
	   IF (.NOT. expansion_defined) THEN
! Direct integration case: loop over R points  
      DO i=1,n_r_coll
	   IF(mij_k_skip) CYCLE
	   IF(i.eq.mtrx_cutoff_r1) then
			bk_mat_temp1(i,bk_mat_counter) = mtrx_cutoff_chk1
	   ELSE IF(i.eq.mtrx_cutoff_r2) then
			bk_mat_temp1(i,bk_mat_counter) = mtrx_cutoff_chk2
	   ELSE
!		Write(*,*) '#1.0'
      CALL	INTEGRATOR_TEMP(intgeral,k,i,cyc)	  	
      bk_mat_temp1(i,bk_mat_counter) = intgeral*conv_unit_e
	   END IF
		ENDDO
		
		ELSE
! Expansion-based computation (more efficient)
	   ALLOCATE(Mat_el_R_array(n_r_coll))
!	   Write(*,*) '#1.1'
	   CALL EXPANSION_MATRIX_ELEMENT(Mat_el_R_array,k,cyc)
	   bk_mat_temp1(:,bk_mat_counter) = Mat_el_R_array(:)*conv_unit_e
	   DEALLOCATE(Mat_el_R_array)
      ENDIF
	  
! Check cutoff after computing the matrix elements
      IF(max(abs(bk_mat_temp1(mtrx_cutoff_r1,bk_mat_counter)),
     &     abs(bk_mat_temp1(mtrx_cutoff_r2,bk_mat_counter))).lt.
     &     MIJ_ZERO_CUT) then
	   bk_mat_temp1(:,bk_mat_counter) = 0d0
	   mij_k_skip = .TRUE.
	   bk_non_zero_counter = bk_non_zero_counter + 1
	   ENDIF
	   ENDIF
		
	  bgn_tym = MPI_Wtime()
	  IF(bk_mat_counter.eq.1000 .or. k.eq.k_fn_mpi) then
	  IF(write_check_file_defined) then
      IF(abs(bgn_tym-bk_tym1) < TIME_MIN_CHECK*60.d0) then
	  call bk_print_matrix (k_st_mpi, k_fn_mpi, k,
     & bk_mat_counter, bk_mat_temp1, mt_chk, cyc_cntr)
	  bk_mat_counter = 0
	  k_rstrt_chk = k + 1
	  ELSE
	  exit
!	  END_tym = MPI_Wtime()
!	  calc_tym = END_tym - bgn_tym
!	  tot_tym_prn = tot_tym_prn + calc_tym
	  END IF
	  
	  ELSE
	  call bk_print_matrix (k_st_mpi, k_fn_mpi, k,
     & bk_mat_counter, bk_mat_temp1, mt_chk, cyc_cntr)
	  bk_mat_counter = 0
	  END IF
	  END IF
	  
	  IF(write_check_file_defined .and. 
     & abs(bgn_tym-bk_tym1) > TIME_MIN_CHECK*60.d0) exit
      ENDDO
	  	  
	  IF(write_check_file_defined .and. k_rstrt_chk < k_fn_mpi) then
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  CALL MPI_GATHER(k_rstrt_chk,1,MPI_INTEGER,bk_k_gather,1,
     & MPI_INTEGER,0,MPI_COMM_WORLD,ierr_mpi)
	  IF(myid.eq.0) then
	  open(11,file="Restart_info.DAT")
	  do i = 1, nproc
	  write(11,*) bk_k_gather(i)
	  END do
	  close(11)
	  write(*,'(a,a)')
     & "Matrix calculations did not finish due to time limit. ",
     & "Please restart. Program will stop now."
	  END IF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      STOP
	  END IF
	  
	  IF(bikram_rebalance) then
	  bk_tym2 = MPI_Wtime()
	  bk_tym = bk_tym2 - bk_tym1
	  open(111, file = "Time_per_proc.out", access = 'appEND')
	  write(111,'(a,i5,a,f12.5)')"#proc = ",myid,", time = ",bk_tym
	  close(111)
	  END IF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  
      CALL MPI_GATHER(k_fn_mpi-k_st_mpi+1-bk_non_zero_counter,1,
     & MPI_INTEGER,bk_non_zero_mij_gather,1,MPI_INTEGER,0,
     & MPI_COMM_WORLD,ierr_mpi)
	  IF(myid.eq.0) then
	  bk_non_zero_counter = 0
	  do i = 1, nproc
	  bk_non_zero_counter = bk_non_zero_counter + 
     & bk_non_zero_mij_gather(i)
	  END do
	  write(bk_matrix_path5, '(a,a)') trim(bk_dir2), 
     & "/MATRIX_NONZERO_INDEX.DAT"
	  open(11,file=trim(bk_matrix_path5))
	  write(11,'(a9,i16)') "#Files = ", nproc
	  write(11,'(a28,i16)') "Non-Zero #Matrix Elements = ", 
     & bk_non_zero_counter
	  write(11,'(a21,e19.12,a6)') "The Matrix Cut-Off = ", 
     & MIJ_ZERO_CUT*eVtown*autoeV, " cm^-1"
	  write(11,'(a16,2x,a16)') '          file_#', '   #non_zero_Mij'
	  do i = 1, nproc
	  write(11,'(i16,2x,i16)') i, bk_non_zero_mij_gather(i)
	  END do
	  close(11)
	  END IF
	  
	  
	  ELSE
	  IF(myid.eq.0) then
	  write(*,'(a)')"The #procs is larger than #elements in matrix."
	  write(*,'(a)')"It is not optimal to have many files with 1 Mij."
	  write(*,'(a)')"Program will stop now."
	  END IF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	 
	  stop
	  return
	  END IF
	  
	  IF(myid.eq.0) then
	  call bk_matrix_combination(total_size,nproc)
	  IF(bikram_rebalance) then
	  write(bk_matrix_path5, '(4(a))') trim(bk_dir2), 
     & "/MATRIX_NONZERO_INDEX.DAT ", trim(bk_dir33), "/"
	  call system ( "cp " // trim(bk_matrix_path5) )
	  write(bk_matrix_path5, '(4(a))') trim(bk_dir2), 
     & "/MATRIX_FILE_INDEX.DAT ", trim(bk_dir33), "/"
	  call system ( "cp " // trim(bk_matrix_path5) )
	  END IF
	  END IF
	  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  
      TIME_WRITING_MATRIX = TIME_2_MAT - TIME_1_MAT
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      TIME_MAT_FINISH= MPI_Wtime()
      IF(MYID.EQ.0) WRITE(*,'(a57,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij READING/COMPUTING/SAVING ON DISK"
     & ,",s",(TIME_MAT_FINISH-TIME_MAT_START)/nproc
	  IF(print_matrix_defined.and. myid.eq.0) 
     & WRITE(*,'(a47,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij SAVING ON DISK"
     & ,",s",TIME_WRITING_MATRIX	 
      IF(MYID.EQ.0) PRINT*, "ALL WORK ON MATRIX IS DONE"
      IF(run_prog_defined) THEN
	  IF(myid.eq.0) print*, "You indicated Parallel Matrix printing. "
     & , "Start trajectory calculation separately. ",
     & "Program will stop now."
      STOP	  
	  ELSE
	  stop
      ENDIF
	  
	  return
	  
	  END IF

!! Bikram End.
!!-------------------------------------------------------------------------------------------------	  
	  
	  
	  
	  
!!!!  COMPUTING MATRIX	  
      IF(myid.eq.0) PRINT*,"MATRIX_INI_STARTED"
!      PRINT*, "HELLO"!DELETE	  
      IF(chunk_mpi_size.ge.0) THEN
!! ALLOCATING CHUNK_SIZE ARRAYS
      IF(chunk_mpi_size.gt.0) THEN	  
      ALLOCATE(Mat_el_temp(n_r_coll,k_fn_mpi-k_st_mpi +1)
     &, Mat_el_der_temp(n_r_coll,k_fn_mpi-k_st_mpi +1))
      ALLOCATE(K_SKIPPED_BY_ROUTINE(k_fn_mpi-k_st_mpi +1))
      DO k=k_st_mpi,k_fn_mpi
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.FALSE.	  
      ENDDO
      ENDIF	  
	  
!!!!! ALLOCATION OF SKIPPING ARRAYS
!!!! THIS IS THE CASE WHEN FOR LARGE SYSTEMS IT IS NOT POSSIBLE TO
!! OR PRACTICAL TO STORE IN THE MEMORY THE ENTIRE GRID
	  
      IF(make_grid_file .and. (coll_type.gt.5 .or. coll_type.eq.0))
     & THEN
      IF(coll_type.eq.1) THEN
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) n_beta,nr_cold
      ALLOCATE(V_2_1(n_beta,nr_cold))
      READ(1)  V_2_1
      CLOSE(1)
      ENDIF	  
      IF(MYID.EQ.0) THEN
      PRINT*, "MATRIX WILL BE COMPUTED USING V-GRID FILE"	  
      SELECT CASE(coll_type)
      CASE(6)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) n_beta1,n_beta2,n_gamma1,n_r_vib1,
     & n_r_vib2, nr_cold
      ALLOCATE(V_VIB_2_2(n_gamma1,n_beta2,n_beta1,n_r_vib2,
     & n_r_vib1,nr_cold)) 
      CASE(7)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) n_alpha2,n_beta1,n_beta2,n_gamma1,nr_cold 
      ALLOCATE(
     & V_3_2(n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))		  
      CASE(8)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) n_alpha2,n_beta1,n_beta2,n_gamma1,nr_cold 
      ALLOCATE(
     & V_3_2(n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  
      CASE(9)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat)
     &	n_alpha2,n_beta1,n_beta2,n_gamma1,n_gamma2,nr_cold 
      ALLOCATE(
     & V_3_3(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  	  
      CASE(0)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat)
     & n_alpha2,n_beta1,n_beta2,n_gamma1,n_gamma2,nr_cold 
      ALLOCATE(
     & V_3_3(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  
      END SELECT
      PRINT*,"VIJ_BUFFER_ALLOCATED"	
!!!!! READS AND PASSE V_ij BUFFER FOR FURTHER BROADCASTING,
!!!! ONLY FOR CASE OF ASYMETRIC + ASYMETRIC	
      IF(istat.gt.0) THEN
      CRITICAL_ERROR=.TRUE.
      PRINT*,"ERROR in V_GRID"	  
      ENDIF	  
      ENDIF
      CALL MPI_BCAST(CRITICAL_ERROR, 1, MPI_LOGICAL,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(CRITICAL_ERROR) STOP	 
      CALL MPI_BCAST(n_alpha2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_beta1, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_beta2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_gamma1, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_gamma2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_r_vib1, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_r_vib2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	 
      CALL MPI_BCAST(nr_cold, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(nr_cold.lt.n_r_coll) THEN 
      IF(MYID.EQ.0) THEN
      CRITICAL_ERROR=.TRUE.
      PRINT*,"ERROR in V_GRID"
      ENDIF
      STOP	  
      ENDIF	  
	  
      SELECT CASE(coll_type)
      CASE(6)
      ALLOCATE(V_VIB_2_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_r_vib2,
     & n_r_vib1)) 
      CASE(7)
      ALLOCATE(V_3_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_alpha2))	 	  
      CASE(8)
      ALLOCATE(V_3_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      CASE(9)
      ALLOCATE(
     & V_3_3_int_buffer(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      CASE(0)
!!!  MINI BUFFERS FOR EACH R-VALUE TO AVOID MEMORY LIMITATIONS	  
      ALLOCATE(
     & V_3_3_int_buffer(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2))
!      PRINT*,"BUFFER_ALLOCATED FOR=",MYID
      END SELECT	
		
! First, read all the grid data
      DO i=1,n_r_coll
      IF(MYID.EQ.0) THEN	  
      SELECT CASE(coll_type)
      CASE(6)

      READ(1) 
     & V_VIB_2_2_int_buffer
      V_VIB_2_2(:,:,:,:,:,i)=
     & V_VIB_2_2_int_buffer	  

      CASE(7)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_2_int_buffer
      V_3_2(:,:,:,:,i)=
     & V_3_2_int_buffer 	  
      CASE(8)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_2_int_buffer
      V_3_2(:,:,:,:,i)=
     & V_3_2_int_buffer 

      CASE(9)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_3_int_buffer
      V_3_3(:,:,:,:,:,i)=
     & V_3_3_int_buffer 
 
      CASE(0)	  

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_3_int_buffer
      V_3_3(:,:,:,:,:,i)=
     & V_3_3_int_buffer	 
      END SELECT
      IF(i.ge.i_nr_ini .and. i.le. i_nr_fin) THEN	  
	   PRINT*,"BUFFER HAS BEEN READ FOR Ri = ",i
      ENDIF	  
      ENDIF
		
      SELECT CASE(coll_type)
      CASE(6)
      buffer_size_V = n_r_vib1*n_r_vib2*n_beta1*n_beta2*n_gamma1
      CASE(7)
      buffer_size_V = n_gamma1*n_beta2*n_beta1*n_alpha2		  
      CASE(8)
      buffer_size_V = n_gamma1*n_beta2*n_beta1*n_alpha2	  
      CASE(9)
      buffer_size_V = n_gamma2*n_gamma1*n_beta2*n_beta1*n_alpha2	  
      CASE(0)	  
      buffer_size_V = n_gamma2*n_gamma1*n_beta2*n_beta1*n_alpha2
      END SELECT	
		
      IF(MYID.EQ.0) PRINT*,"BUFFER OF VIJ SIZE= ",
     & buffer_size_V!,size(V_3_3_int_buffer)	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
		
!!!!! BROADCASTING MINI _BUFFER	  
      IF(MYID.EQ.0) PRINT*,"BROADCASTING OF VIJ_BUFFER STARTED"
!      IF(MYID.EQ.0) WRITE(3,*)	V_3_3_int_buffer
      SELECT CASE(coll_type)
      CASE(6)
      CALL MPI_BCAST(V_VIB_2_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CASE(7)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(8)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(9)
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(0)  
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      END SELECT	 
		
      IF(MYID.EQ.0) PRINT*,"THE BUFFER HAS BEEN BROADCASTED"
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      ENDDO ! End of first loop over i for reading grid data
		
! SECOND LOOP: Matrix computation with proper separation
      DO  k=k_st_mpi,k_fn_mpi
		
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE	 
      
      IF (.NOT. expansion_defined) THEN
! Direct integration case: loop over R points
         DO i=1,n_r_coll
            IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE
!				write(*,*) "#2.1"
            CALL INTEGRATOR_TEMP(intgeral,k,i,cyc)
            Mat_el_temp(i,k-k_st_mpi+1) = intgeral*conv_unit_e
         ENDDO
      ELSE
! Expansion-based computation (more efficient)
         ALLOCATE(Mat_el_R_array(n_r_coll))
!			write(*,*) "#2.2"
         CALL EXPANSION_MATRIX_ELEMENT(Mat_el_R_array,k,cyc)
         Mat_el_temp(:,k-k_st_mpi+1) = Mat_el_R_array(:)*conv_unit_e
         DEALLOCATE(Mat_el_R_array)
      END IF
      ENDDO ! End of loop over k
		
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
		
! finding #r for matrix truncation
	  IF(.not.cut_r) then
	  allocate(deviation_r(n_r_coll))
	  do ii = 1, n_r_coll
	  deviation_r(ii) = abs(R_COM(ii) - bikram_cutoff_r1)
	  END do
	  mtrx_cutoff_r1 = minloc(deviation_r,1)
	  deallocate(deviation_r)
	  allocate(deviation_r(n_r_coll))
	  do ii = 1, n_r_coll
	  deviation_r(ii) = abs(R_COM(ii) - bikram_cutoff_r2)
	  END do
	  mtrx_cutoff_r2 = minloc(deviation_r,1)
	  deallocate(deviation_r)
	  cut_r = .true.
	  IF(myid.eq.0) write(*,'(2(a,i0,a,f0.3))') 
     & "Truncation #1 at #R = ", mtrx_cutoff_r1, ", R = ", 
     & R_COM(mtrx_cutoff_r1), ", and #2 at #R = ",mtrx_cutoff_r2, 
     & ", R = ",R_COM(mtrx_cutoff_r2)
	  MIJ_ZERO_CUT = MIJ_ZERO_CUT/eVtown/autoeV
	  END IF

! ============================================================================
! MATRIX ELEMENT PROCESSING: Two main branches 
! ============================================================================	  
! Branch 1: Process matrix elements loaded from grid file
! Check IF elements fall below threshold and mark them for exclusion
      DO  k=k_st_mpi,k_fn_mpi
! Skip elements that have already been marked for exclusion
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE	  
! Check IF matrix elements at specIFic cutoff rows are below threshold
! Compare the maximum absolute value of elements at two reference rows
! against the zero cutoff criterion (MIJ_ZERO_CUT)
      IF(max(abs(Mat_el_temp(mtrx_cutoff_r1,k-k_st_mpi+1)),
     & abs(Mat_el_temp(mtrx_cutoff_r2,k-k_st_mpi+1))).LT.
     & MIJ_ZERO_CUT) THEN    !!! EXCLUDE SOME MATRIX ELEMENTS IF THEY ARE BELOW CERTAIN VALUE
! Zero out matrix elements and derivatives that are below threshold
      Mat_el_temp(:,k-k_st_mpi+1) = 0d0
      Mat_el_der_temp(:,k-k_st_mpi+1) = 0d0	  
! Mark this k-index as skipped for future iterations
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.TRUE.	  
      ENDIF	  
	   END do
	  
      ELSE
		
!!!! IF NOT THE GRID FORM FILE DEFINED COMPUTE MATRIX ELEMENTS
      IF(MYID.EQ.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "COMPUTING MATRIX ELEMENTS STARTED"	  
	   bk_tym1 = MPI_Wtime()	!Bikram
!! Main computation loop over k-indices assigned to this MPI process		
      DO  k=k_st_mpi,k_fn_mpi	  
! Bikram July 9, 2023: Display computation progress at 10% intervals
	   IF((k_fn_mpi - k_st_mpi) .ge. 10) then
	   IF(mod((k - k_st_mpi), int((k_fn_mpi - k_st_mpi)/10)) == 0) then
	   bk_tym2 = MPI_Wtime()
	   bk_tym = bk_tym2 - bk_tym1
	   write(*,'(2(a, i5),a,f12.3)') "proc_id = ",myid,", Progress = ",
     & int(dble(k - k_st_mpi)/dble(k_fn_mpi - k_st_mpi)*100.d0), 
     & "%, Time(sec.) = ", bk_tym
	   percent_counter = percent_counter + 1
	   END IF
	   END IF
! Print final 100% completion message		
	   IF(k == k_fn_mpi) then
	   bk_tym2 = MPI_Wtime()
	   bk_tym = bk_tym2 - bk_tym1
	   write(*,'(2(a, i5),a,f12.3)') "proc_id = ",myid,", Progress = ",
     & int(100.d0), "%, Time(sec.) = ", bk_tym	  
	   END IF
! End of progress reporting section
! Bikram END.	  
	  
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE !!! SKIP SOME ELEMENTS IF NESSEACRY	
! ================================================================
! MATRIX ELEMENT COMPUTATION: Two methods - DI & Expansion		
! ================================================================
		IF (.NOT. expansion_defined) THEN
! Direct integration case: loop over R points
      DO i=1,n_r_coll
! Skip points outside the defined range [i_nr_ini, i_nr_fin]		
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
! Compute matrix element using numerical integration
      CALL	INTEGRATOR_TEMP(intgeral,k,i,cyc)	
! Store result with unit conversion (conv_unit_e)		
      Mat_el_temp(i,k-k_st_mpi+1) = intgeral*conv_unit_e
      ENDDO
		
! Method 2: Use expansion-based computation (more efficient)
		ELSE			
! Allocate temporary array for expansion results	
		Allocate(exp_result_array(n_r_coll))
! Compute matrix elements using expansion method
		CALL EXPANSION_MATRIX_ELEMENT(exp_result_array,k,cyc)	
!		write(1544,*) exp_result_array/ 219474.6d0
! Store result with unit conversion (conv_unit_e)	
		Mat_el_temp(:,k-k_st_mpi+1) = exp_result_array(:)*conv_unit_e
		DEALLOCATE(exp_result_array)
		ENDIF
		
! finding #r for matrix truncation
	  IF(.not.cut_r) then
	  allocate(deviation_r(n_r_coll))
	  do ii = 1, n_r_coll
	  deviation_r(ii) = abs(R_COM(ii) - bikram_cutoff_r1)
	  END do
	  mtrx_cutoff_r1 = minloc(deviation_r,1)
	  deallocate(deviation_r)
	  allocate(deviation_r(n_r_coll))
	  do ii = 1, n_r_coll
	  deviation_r(ii) = abs(R_COM(ii) - bikram_cutoff_r2)
	  END do
	  mtrx_cutoff_r2 = minloc(deviation_r,1)
	  deallocate(deviation_r)
	  cut_r = .true.
	  IF(myid.eq.0) write(*,'(2(a,i0,a,f0.3))') 
     & "Truncation #1 at #R = ", mtrx_cutoff_r1, ", R = ", 
     & R_COM(mtrx_cutoff_r1), ", and #2 at #R = ",mtrx_cutoff_r2, 
     & ", R = ",R_COM(mtrx_cutoff_r2)
	  MIJ_ZERO_CUT = MIJ_ZERO_CUT/eVtown/autoeV
	  END IF
	  
      IF(max(abs(Mat_el_temp(mtrx_cutoff_r1,k-k_st_mpi+1)),
     & abs(Mat_el_temp(mtrx_cutoff_r2,k-k_st_mpi+1))).LT.
     & MIJ_ZERO_CUT) THEN !!!! THIS IS EXCLUSION CONDITION
      Mat_el_temp(:,k-k_st_mpi+1) = 0d0
      Mat_el_der_temp(:,k-k_st_mpi+1) = 0d0	  
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.TRUE.	  
      ENDIF	  
	  
      ENDDO
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      ENDIF	
!!!      SPLININING STARTED
      IF(MYID.EQ.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "SPLINING OF Mij STARTED"	  
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE
      deriv_bgn = (Mat_el_temp(2,k-k_st_mpi+1)
     & - Mat_el_temp(1,k-k_st_mpi+1))/
     & (R_COM(2)-R_COM(1)) 
      deriv_END = (Mat_el_temp(n_r_coll,k-k_st_mpi+1)
     & - Mat_el_temp(n_r_coll-1,k-k_st_mpi+1))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1)) 	  
      CALL  spline(R_COM,Mat_el_temp(:,k-k_st_mpi+1),
     & n_r_coll,deriv_bgn,deriv_END,
     & Mat_el_der_temp(:,k-k_st_mpi+1))	 
      ENDDO	 
      IF(MYID.EQ.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "SPLINING OF Mij FINISHED"	 
!!! GATHERING MATRIX IN Mij.dat IN ONE PROCESSOR	 

      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF(chunk_mpi_size.gt.0) THEN	  
      CALL MPI_GATHER(Mat_el_temp,task_size,MPI_DOUBLE_PRECISION,Mat_el,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 	  
      CALL MPI_GATHER(Mat_el_der_temp,task_size,MPI_DOUBLE_PRECISION,
     & Mat_el_der,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
      ENDIF	 
      ENDIF		 
!!!! COMPUTING THE REST OF THE MATRIX MIJ
      IF(expansion_defined) THEN
      CALL READ_EXPANSION_TERMS	  
      ENDIF
		
      IF(total_size.gt.nproc*chunk_mpi_size) THEN
      IF(MYID.EQ.0) PRINT*,"COMPUTING THE RESIDUE OF Mij"
      ALLOCATE(Mat_el_res(n_r_coll),Mat_el_der_res(n_r_coll))
      ALLOCATE(Mat_el_resf(n_r_coll,nproc),
     & Mat_el_der_resf(n_r_coll,nproc))
      Mat_el_res = 0d0
      Mat_el_der_res = 0d0
		
		IF (.NOT. expansion_defined) THEN
! Non-expansion case: loop over R points
      DO i=1,n_r_coll
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
      IF(K_SKIPPED_RES) CYCLE
!!! AGAIN IF READING FROM GRID DEFINED!	  
      IF(make_grid_file .and. (coll_type.gt.5 .or. coll_type.eq.0))
     & THEN
      IF(MYID.EQ.0) THEN
		
      SELECT CASE(coll_type)
      CASE(6)
      V_VIB_2_2_int_buffer(:,:,:,:,:) = 
     & V_VIB_2_2(:,:,:,:,:,i)
      CASE(7)
      V_3_2_int_buffer(:,:,:,:) = V_3_2(:,:,:,:,i)	  
      CASE(8)
      V_3_2_int_buffer(:,:,:,:) = V_3_2(:,:,:,:,i)	  
      CASE(0)	  
      V_3_3_int_buffer(:,:,:,:,:) = V_3_3(:,:,:,:,:,i)
      END SELECT	  
      ENDIF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
		
      SELECT CASE(coll_type)
      CASE(6)
      CALL MPI_BCAST(V_VIB_2_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CASE(7)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(8)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(9)
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(0)	  
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      PRINT*, "CHUNK BUFFER HAS BEEN BROADCASTED FOR IR=",i	 
      END SELECT	 
      ENDIF	  
		
      IF(myid.lt.total_size-nproc*chunk_mpi_size) THEN
      k=nproc*chunk_mpi_size+1+myid	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i,cyc)	
!		Write(*,*) '#3'
      Mat_el_res(i) = intgeral*conv_unit_e
      ENDIF	  
      ENDDO
		
	   ELSE
! Expansion case: process all R values at once for each k
		IF(myid.lt.total_size-nproc*chunk_mpi_size) THEN
		k = nproc*chunk_mpi_size + 1 + myid
		IF(.NOT. K_SKIPPED_RES) THEN
		ALLOCATE(Mat_el_R_array(n_r_coll))
!		Write(*,*) '#3.1'
		CALL EXPANSION_MATRIX_ELEMENT(Mat_el_R_array,k,cyc)
		Mat_el_res(:) = Mat_el_R_array(:)*conv_unit_e	
		DEALLOCATE(Mat_el_R_array)
		ENDIF
		ENDIF
		ENDIF
		
! finding #r for matrix truncation
	  IF(.not.cut_r) then
	  allocate(deviation_r(n_r_coll))
	  do ii = 1, n_r_coll
	  deviation_r(ii) = abs(R_COM(ii) - bikram_cutoff_r1)
	  END do
	  mtrx_cutoff_r1 = minloc(deviation_r,1)
	  deallocate(deviation_r)
	  allocate(deviation_r(n_r_coll))
	  do ii = 1, n_r_coll
	  deviation_r(ii) = abs(R_COM(ii) - bikram_cutoff_r2)
	  END do
	  mtrx_cutoff_r2 = minloc(deviation_r,1)
	  deallocate(deviation_r)
	  cut_r = .true.
	  IF(myid.eq.0) write(*,'(2(a,i0,a,f0.3))') 
     & "Truncation #1 at #R = ", mtrx_cutoff_r1, ", R = ", 
     & R_COM(mtrx_cutoff_r1), ", and #2 at #R = ",mtrx_cutoff_r2, 
     & ", R = ",R_COM(mtrx_cutoff_r2)
	  MIJ_ZERO_CUT = MIJ_ZERO_CUT/eVtown/autoeV
	  END IF
	  
      IF(max(abs(Mat_el_res(mtrx_cutoff_r1)),
     & abs(Mat_el_res(mtrx_cutoff_r2))).LT.MIJ_ZERO_CUT) THEN
      Mat_el_res(:) = 0d0
      Mat_el_der_res(:) = 0d0	  
      K_SKIPPED_RES=.TRUE.	  
      ENDIF
	  
!!! SPLINING THE REST OF MIJ
      IF(MYID.eq.0) PRINT*,"SPLINING RESIDUE OF Mij"
      IF(.NOT. K_SKIPPED_RES) THEN
      deriv_bgn = (Mat_el_res(2)
     & - Mat_el_res(1))/
     & (R_COM(2)-R_COM(1)) 
      deriv_END = (Mat_el_res(n_r_coll)
     & - Mat_el_res(n_r_coll-1))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1)) 	  	  
      CALL  spline(R_COM,Mat_el_res,n_r_coll,
     & deriv_bgn,deriv_END,
     & Mat_el_der_res)	  
      ENDIF

!!!!!  GATHERING THE REST	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
!	  print *, 'calc_alex', myid, tot_tym_calc						!Bikram
      CALL MPI_GATHER(Mat_el_res,n_r_coll,MPI_DOUBLE_PRECISION,
     & Mat_el_resf,
     &	  n_r_coll,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 	  
      CALL MPI_GATHER(Mat_el_der_res,n_r_coll,MPI_DOUBLE_PRECISION,
     & Mat_el_der_resf,
     &	  n_r_coll,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
!!!!!  GATHERING THE REST      
      IF(MYID.EQ.0) THEN
      DO i=1,total_size-nproc*chunk_mpi_size
      Mat_el(:,i+nproc*chunk_mpi_size) = Mat_el_resf(:,i)
      Mat_el_der(:,i+nproc*chunk_mpi_size) = Mat_el_der_resf(:,i)	  
      ENDDO	  
      ENDIF	  
      ENDIF
	  TIME_1_MAT = MPI_Wtime()
      IF(MYID.EQ.0 .and. print_matrix_defined) THEN
	  IF(.not. bikram_mij_multiprint) then
      IF(.not.unformat_defined) THEN	  
!	  bgn_tym = MPI_Wtime()											!Bikram
      CALL PRINT_MIJ
!	  END_tym = MPI_Wtime()											!Bikram
!	  calc_tym = END_tym - bgn_tym									!Bikram
!	  tot_tym_prn = tot_tym_prn + calc_tym							!Bikram
      WRITE(*,'(a35,1x,a2,a11,a2)')
     & "MATRIX HAS BEEN SAVED INTO THE FILE"
     & ,"""",MATRIX_NAME_MIJ,""""
      ELSE	 
      CALL PRINT_MIJ_USER
      WRITE(*,'(a35,1x,a2,a11,a2)')
     & "MATRIX HAS BEEN SAVED INTO THE FILE"
     & ,"""",MATRIX_NAME_MIJ_UF,""""	  
      ENDIF	  
	  
	  ELSE
	  print*, "Something wrong in Parallel_IO"
	  stop
	  ENDIF
!	  print *, 'prn_alex', myid, tot_tym_prn						!Bikram
      ENDIF
	  
!!!!! MATRIX BROADCASTING	  
!      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
!	  print *, myid, size(Mat_el)
!      CALL MPI_BCAST(Mat_el, n_r_coll*total_size, MPI_REAL8,0,
!     &  MPI_COMM_WORLD,ierr_mpi)
!	  IF (print_matrix_defined .and. bikram_mij_multiprint) then
!	  do i = (myid+1), total_size, nproc
!	  print *, i, myid, nproc
!	  call bk_print_matrix (i)
!	  print *, 'finished', i
!	  END do
!	  ENDIF
	  
	  TIME_2_MAT = MPI_Wtime()
      TIME_WRITING_MATRIX = TIME_2_MAT - TIME_1_MAT		  
!!! WAITING FOR OTHER PROCESSORS
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF(.not.run_prog_defined) THEN
      TIME_MAT_FINISH= MPI_Wtime()
      IF(MYID.EQ.0) WRITE(*,'(a57,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij READING/COMPUTING/SAVING ON DISK"
     & ,",s",(TIME_MAT_FINISH-TIME_MAT_START)/nproc
	  IF(print_matrix_defined.and. myid.eq.0) 
     & WRITE(*,'(a47,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij SAVING ON DISK"
     & ,",s",TIME_WRITING_MATRIX	 
      IF(MYID.EQ.0) PRINT*, "ALL WORK ON MATRIX IS DONE"
      STOP	  
      ENDIF

! Bikram Start Dec 2019:	  
	  IF(bikram_mij_shIFt) then
	  DO i=1,n_r_coll
	  Mat_el(i,:) = Mat_el(i,:) - Mat_el(n_r_coll,:)
      ENDDO 
	  ENDIF
! Bikram End.
	  
	  IF(coupled_states_defined .and. myid.eq.0) THEN
      ALLOCATE(mat_buffer(n_r_coll,total_size))
      ALLOCATE(mat_buffer_der(n_r_coll,total_size))
	  mat_buffer = Mat_el
	  mat_buffer_der = Mat_el_der
      ENDIF		  
      IF(.not.mpi_task_defined) THEN
!!!!! MATRIX BROADCASTING	  
      CALL MPI_BCAST(Mat_el, n_r_coll*total_size, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CALL MPI_BCAST(Mat_el_der, n_r_coll*total_size, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(R_COM, n_r_coll, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
!!!!!  CHECKING IF Mij.dat HAS BEEN BROADCASTED CORRECTLY	 
      DO  k=k_st_mpi,k_fn_mpi
      DO i=1,n_r_coll	  
      IF(Mat_el_temp(i,k-k_st_mpi+1).ne.Mat_el(i,k)) THEN
      CALL MPI_FINALIZE (ierr_mpi)
      PRINT*,R_COM
!!!!!! IF NOT STOP THE PROGRAM	  
      STOP "ERROR:MATRIX_INI_WRONG"	  
      ENDIF
      IF(Mat_el_der_temp(i,k-k_st_mpi+1).ne.Mat_el_der(i,k)) THEN
      CALL MPI_FINALIZE (ierr_mpi)
!!!!!! IF NOT STOP THE PROGRAM	  
      STOP "ERROR:MATRIX_INI_WRONG_IN_DER"	  
      ENDIF	  
      ENDDO
      ENDDO 
      ENDIF
      IF(mpi_task_defined) THEN
!!!!! IF MPI_TASKS DEFINED THEN IT WILL REALLOCATE Mij.dat	  
      IF(MYID.eq.0) THEN	  
      ALLOCATE(Mat_el_non_zero(n_r_coll,total_size))
      ALLOCATE(Mat_el_non_zero_der(n_r_coll,total_size))
      Mat_el_non_zero = Mat_el
      Mat_el_non_zero_der = Mat_el_der
      DEALLOCATE(Mat_el,Mat_el_der)	  
      ENDIF	  
!      TIME_MAT_START = MPI_Wtime()	  
      IF(myid.eq.0) PRINT*,"MPI TASK PER TRAJECTORY WILL BE USED"
      IF(myid.eq.0)	
     & WRITE(*,'(a53,1x,i4)')
     & "MPI TASKS WHICH ARE ASSOCIATED WITH ONE TRAJECTORY = ",
     & mpi_task_per_traject	 
      traject_roots = nproc/mpi_task_per_traject
      IF(traject_roots*mpi_task_per_traject.ne.nproc) THEN
      STOP "mpi_task_number must be a delimeter of nproc"	  
      ENDIF	  
      ALLOCATE(mpi_traject_roots(traject_roots))
      ALLOCATE(portion_of_MIJ_per_task(2,nproc))
      ALLOCATE(portion_of_state_per_task(2,nproc))
      ALLOCATE(portion_of_work_per_task(2,nproc))	  
      ALLOCATE(mpi_root_belongs(nproc))	  
!!!!!!!!!!!!!   REMAKE	  
      size_mij_chunk_mpi = total_size/mpi_task_per_traject
      residue_mij_mpi = total_size-
     & size_mij_chunk_mpi*mpi_task_per_traject
	 
      size_state_chunk_mpi = states_size/mpi_task_per_traject
      residue_state_mpi = states_size-
     & size_state_chunk_mpi*mpi_task_per_traject
	 
      size_work_chunk_mpi = (2*states_size+8)/mpi_task_per_traject !!!! DO IN FUTURE
      residue_work_mpi = 2*states_size+8-!!!! DO IN FUTURE
     & size_work_chunk_mpi*mpi_task_per_traject	!!!! DO IN FUTURE
	 
      DO k=1,traject_roots
      mpi_traject_roots(k) = (k-1)*mpi_task_per_traject
      DO k_p=1,mpi_task_per_traject	  
      mpi_root_belongs(k_p+mpi_traject_roots(k))=mpi_traject_roots(k)
      ENDDO	  
      ENDDO
      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_mij_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_MIJ_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_mij_chunk_mpi+1)
      portion_of_MIJ_per_task(2,k_mpi_proc) =
     & (k_p)*(size_mij_chunk_mpi+1)
      ENDDO
      total_size_check = total_size_check + size_mij_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_MIJ_per_task(1,k_mpi_proc)
     & = residue_mij_mpi*(size_mij_chunk_mpi+1)+1  
     & + (k_p-1-residue_mij_mpi)*size_mij_chunk_mpi
      portion_of_MIJ_per_task(2,k_mpi_proc)
     & = residue_mij_mpi*(size_mij_chunk_mpi+1)+	  
     & (k_p-residue_mij_mpi)*size_mij_chunk_mpi
      ENDDO	 
      total_size_check = total_size_check + size_mij_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_MIJ_per_task(2,k_mpi_proc).ne.total_size) 
     & STOP "WRONG MATRIX ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO


      	  
      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_state_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_state_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_state_chunk_mpi+1)
      portion_of_state_per_task(2,k_mpi_proc) =
     & (k_p)*(size_state_chunk_mpi+1)
      ENDDO
      state_size_check = state_size_check + size_state_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_state_per_task(1,k_mpi_proc)
     & = residue_state_mpi*(size_state_chunk_mpi+1)+1  
     & + (k_p-1-residue_state_mpi)*size_state_chunk_mpi
      portion_of_state_per_task(2,k_mpi_proc)
     & = residue_state_mpi*(size_state_chunk_mpi+1)+	  
     & (k_p-residue_state_mpi)*size_state_chunk_mpi
      ENDDO	 
      state_size_check = state_size_check + size_state_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_state_per_task(2,k_mpi_proc).ne.states_size) 
     & STOP "WRONG MATRIX ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO	


      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_work_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_work_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_work_chunk_mpi+1)
      portion_of_work_per_task(2,k_mpi_proc) =
     & (k_p)*(size_work_chunk_mpi+1)
      ENDDO
      work_size_check = work_size_check + size_work_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_work_per_task(1,k_mpi_proc)
     & = residue_work_mpi*(size_work_chunk_mpi+1)+1  
     & + (k_p-1-residue_work_mpi)*size_work_chunk_mpi
      portion_of_work_per_task(2,k_mpi_proc)
     & = residue_work_mpi*(size_work_chunk_mpi+1)+	  
     & (k_p-residue_work_mpi)*size_work_chunk_mpi
      ENDDO	 
      work_size_check = work_size_check + size_work_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_work_per_task(2,k_mpi_proc).ne.states_size*2+8) 
     & STOP "WRONG MATRIX ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO

	  
	  
!      PRINT*,"STATES",portion_of_state_per_task	  
      IF(state_size_check.ne.states_size) THEN
      IF(myid.eq.0) WRITE(*,*)state_size_check,states_size	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF	
	  
      IF(total_size_check.ne.total_size) THEN
      IF(myid.eq.0) WRITE(*,*)total_size_check,total_size	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF
		  
      IF(work_size_check.ne.states_size*2+8) THEN
      IF(myid.eq.0) WRITE(*,*)work_size_check,states_size*2+8	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF	  
	  
      total_size_mpi = portion_of_MIJ_per_task(2,myid+1) - 
     & portion_of_MIJ_per_task(1,myid+1) + 1	  
      ALLOCATE(Mat_el(n_r_coll,total_size_mpi))
      ALLOCATE(Mat_el_der(n_r_coll,total_size_mpi))
      CALL MPI_Comm_group(MPI_COMM_WORLD, wrld_group,ierr_mpi)
      ALLOCATE(process_rank_distr(mpi_task_per_traject,traject_roots))
      ALLOCATE(comms_distr(mpi_task_per_traject),
     & groups_distr(mpi_task_per_traject))   	  
      DO i=1,traject_roots	  
      DO k=1,mpi_task_per_traject
      process_rank_distr(k,i) = k - 1 + (i-1)*mpi_task_per_traject     	  
      ENDDO
      ENDDO	

      DO i=1,mpi_task_per_traject      	  
      CALL MPI_Group_incl(wrld_group, traject_roots,
     & process_rank_distr(i,:), groups_distr(i),ierr_mpi)
      CALL MPI_Comm_create(MPI_COMM_WORLD,groups_distr(i),
     & comms_distr(i),ierr_mpi)
      ENDDO	  
!      PRINT*,myid,total_size_mpi	  
      tag1 = 1
      tag2 = 2
      IF(MYID.eq.0) PRINT*,"Mij ROOT PROC DISTRIBUTION STARTED"	  
      IF(MYID.eq.0) THEN
      DO k=2,mpi_task_per_traject
      total_size_mpi= portion_of_MIJ_per_task(2,k) - 
     & portion_of_MIJ_per_task(1,k) + 1	
      task_portion_size = total_size_mpi*n_r_coll	 
      ALLOCATE(buffer_mpi_portion(n_r_coll,total_size_mpi))
      buffer_mpi_portion = Mat_el_non_zero(:,
     & portion_of_MIJ_per_task(1,k):portion_of_MIJ_per_task(2,k)) 	  
      CALL MPI_SEND(buffer_mpi_portion,
     & task_portion_size, MPI_REAL8, k-1, 
     &  tag1, MPI_COMM_WORLD, ierr_mpi)
      buffer_mpi_portion = 	 Mat_el_non_zero_der (:,
     & portion_of_MIJ_per_task(1,k):portion_of_MIJ_per_task(2,k))
      CALL MPI_SEND(buffer_mpi_portion,
     & task_portion_size, MPI_REAL8, k-1, 
     &  tag2, MPI_COMM_WORLD, ierr_mpi)
      DEALLOCATE(buffer_mpi_portion) 
      ENDDO
      total_size_mpi= portion_of_MIJ_per_task(2,1) - 
     & portion_of_MIJ_per_task(1,1) + 1		  
      Mat_el=Mat_el_non_zero(:,
     & portion_of_MIJ_per_task(1,1):portion_of_MIJ_per_task(2,1))
      Mat_el_der=Mat_el_non_zero_der(:,
     & portion_of_MIJ_per_task(1,1):portion_of_MIJ_per_task(2,1))
      DEALLOCATE(Mat_el_non_zero,Mat_el_non_zero_der)	 
      ELSE
      IF(myid.le.mpi_task_per_traject-1) THEN	  
      task_portion_size = total_size_mpi*n_r_coll	  
      CALL MPI_RECV(Mat_el, task_portion_size, MPI_REAL8, 
     & 0, tag1, MPI_COMM_WORLD, status, ierr_mpi)	  
      CALL MPI_RECV(Mat_el_der, task_portion_size, MPI_REAL8, 
     & 0, tag2, MPI_COMM_WORLD, status, ierr_mpi)
      ENDIF	 
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
      IF(myid.eq.0)	 PRINT*, "Mij ROOT PROC DISTRIBUTION DONE"  
      IF(myid.eq.0) PRINT*,"Mij ALL PROC DISTRIBUTION STARTED" 	  
      DO i=1,mpi_task_per_traject
      IF(i-1 .eq. myid - int(myid/mpi_task_per_traject)
     & *mpi_task_per_traject) THEN
      CALL MPI_Comm_rank(comms_distr(i),id_proc_in_group ,ierr_mpi)
      IF(id_proc_in_group.ne.myid/mpi_task_per_traject) PRINT*,
     & "ERROR IN COMUNICATIONS ASSIGNEMNET_1_distr",id_proc_in_group,
     & i-1
      IF(myid.lt.mpi_task_per_traject .and. id_proc_in_group.ne.0 )
     & PRINT*,"ERROR IN COMUNICATIONS ASSIGNEMNET_2_distr"
      task_portion_size = total_size_mpi*n_r_coll	  
      CALL MPI_BCAST(Mat_el,task_portion_size, MPI_REAL8,0,
     &  comms_distr(i),ierr_mpi)
      CALL MPI_BCAST(Mat_el_der,task_portion_size, MPI_REAL8,0,
     &  comms_distr(i),ierr_mpi)
      ENDIF	 
      ENDDO
      IF(myid.eq.0) PRINT*,"Mij ALL PROC DISTRIBUTION DONE"	  
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)	  
      IF(MYID.EQ.0) PRINT *,"MPI_COMMUNICATORS_CREATION STARTED"
      ALLOCATE(comms(traject_roots),groups(traject_roots))     	  
      ALLOCATE(process_rank(traject_roots,mpi_task_per_traject))		  
      DO i=1,traject_roots	  
      DO k=1,mpi_task_per_traject
      process_rank(i,k) = k - 1 + (i-1)*mpi_task_per_traject     	  
      ENDDO
      ENDDO	

      DO i=1,traject_roots      	  
      CALL MPI_Group_incl(wrld_group, mpi_task_per_traject,
     & process_rank(i,:), groups(i),ierr_mpi)
      CALL MPI_Comm_create(MPI_COMM_WORLD,groups(i),comms(i),ierr_mpi)
      ENDDO

      DO i=1,traject_roots
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN
      CALL MPI_Comm_rank(comms(i),id_proc_in_group ,ierr_mpi)	  
      IF(mpi_traject_roots(i)+id_proc_in_group.ne.myid) STOP
     & "WRONG GROUP MPI ASSIGNEMENT"	  
	  
      ENDIF		  
      ENDDO	 
      IF(MYID.EQ.0) PRINT *,"MPI_COMMUNICATORS_CREATION DONE"
      TIME_MAT_FINISH = MPI_Wtime()
      TIME_2_MAT = TIME_MAT_FINISH - TIME_MAT_START
      CALL MPI_BCAST(R_COM, n_r_coll, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(ind_mat, 2*total_size, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      ENDIF		
!	  call resize
      ELSE
	  

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! Bikram Start May 2022:
	  IF (bikram_mij_multiprint .and. bikram_truncate_MTRX) then
	  IF(myid.eq.0) print*, "Retruncating starts"
	  write(bk_matrix_path2,'(a,a)')trim(bk_dir1),"/MATRIX_INFO.DAT"
	  bk_matrix_path2 = trim(bk_matrix_path2)
	  bk_mtrx_re = .false.
	  inquire(file = trim(bk_matrix_path2), exist = bk_mtrx_re)
	  IF(.not.bk_mtrx_re) then
	  IF(myid.eq.0) write(*,'(a,a)') 
     & "Please make sure that MATRIX_FILES directory is provided. ",
     & "Program will stop now."
	  stop
	  END IF
	  IF(myid.eq.0) then
	  call system ( "rm -r " // trim(bk_dir2) )
	  call system ( "mkdir -p " // trim(bk_dir2) )
	  END IF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  
	  IF(.not.unformat_defined) then
	  open(1,file = trim(bk_matrix_path2))      
	  read(1,*)
	  read(1,*)
      read(1,'(a17,1x,i6)') dm7, dmm1
      read(1,'(a12,1x,i9)') dm7, dmm2	  
      read(1,*)
      read(1,*)
	  do ii = 1, dmm1
	  read(1,*)
	  END do
      read(1,*)
	  do ii = 1, dmm2
	  read(1,*)
	  END do
      read(1,*)
	  do ii = 1, n_r_coll
	  read(1,*) dmm1, R_COM(ii)
	  IF(dmm1.ne.ii) then
	  print*, "Error in truncation while reading R-grid"
	  END IF
	  END do
	  close(1)
	  
	  ELSE
	  open(1,file = bk_matrix_path2, form = "unformatted")
	  read(1,iostat=istat) dmm4
	  read(1,iostat=istat) dmm4
	  read(1,iostat=istat) dmm1 
	  read(1,iostat=istat) dmm2
	  read(1,iostat=istat) dmm4
	  IF(dmm4.ne.n_r_coll) then
	  print*,"Error in unformatted matrix truncation", dmm4, n_r_coll
	  stop
	  END IF
	  
      select case(coll_type)
      case(1)
      do ii = 1, dmm1
      read(1,iostat=istat)dmm2,dmm2,dmm2,dmm2,dmm2
      END do	  
      case(2)
      do ii = 1, dmm1
      read(1,iostat=istat)dmm2,dmm2,dmm2,dmm2,dmm2,dmm2
      END do
      case(3)
      do ii = 1, dmm1
      read(1,iostat=istat)dmm2,dmm2,dmm2,dmm2,dmm2,dmm2,dmm2
      END do
      case(4)
      do ii = 1, dmm1
      read(1,iostat=istat)dmm2,dmm2,dmm2,dmm2,dmm2,dmm2,dmm2
      END do
      case(5)
      do ii = 1, dmm1
      IF(.not.identical_particles_defined) then	 
      read(1,iostat=istat)dmm2,dmm2,dmm2,dmm2,dmm2,dmm2
      ELSE
      read(1,iostat=istat)dmm2,dmm2,dmm2,dmm2,dmm2,dmm2,dmm2
      END IF
      END do
      case(6)
      do ii = 1, dmm1
      IF(.not.identical_particles_defined) then	 
      read(1,iostat=istat)dmm2,dmm2,dmm2,dmm2,dmm2,dmm2,dmm2,dmm2
      ELSE
      read(1,iostat=istat)dmm2,dmm2,dmm2,dmm2,dmm2,dmm2,dmm2,dmm2,dmm2
      END IF
      END do
      case(7)
      do ii = 1, dmm1
      read(1,iostat=istat)dmm2,dmm2,dmm2,dmm2,dmm2,dmm2,dmm2,dmm2
      END do
      case(8)
	  do ii = 1, dmm1
      read(1,iostat=istat)dmm2,dmm2,dmm2,dmm2,dmm2,dmm2,dmm2,dmm2
      END do
      case(9)
      do ii = 1, dmm1
      read(1,iostat=istat)dmm2,dmm2,dmm2,dmm2,
     & dmm2,dmm2,dmm2,dmm2,dmm2,dmm2
      END do
      case(0)
      do ii = 1, dmm1
      IF(.NOT.identical_particles_defined) then
      read(1,iostat=istat)dmm5(1),dmm5(2),dmm5(3),dmm5(4),dmm5(5),
     & dmm5(6),dmm5(7),dmm5(8)
      ELSE
      read(1,iostat=istat)dmm2,dmm2,dmm2,dmm2,
     & dmm2,dmm2,dmm2,dmm2,dmm2,dmm2,dmm2	  
      END IF
      END do
	  
      END select
	  do ii = 1, n_r_coll
	  read(1,iostat=istat) dmm1, R_COM(ii)
	  IF(dmm1.ne.ii) then
	  print*, "Error in truncation while reading R-grid"
	  END IF
	  END do
	  close(1)
	  END IF
	  
! finding #r for matrix truncation
	  IF(.not.cut_r) then
	  allocate(deviation_r(n_r_coll))
	  do ii = 1, n_r_coll
	  deviation_r(ii) = abs(R_COM(ii) - bikram_cutoff_r1)
	  END do
	  mtrx_cutoff_r1 = minloc(deviation_r,1)
	  deallocate(deviation_r)
	  allocate(deviation_r(n_r_coll))
	  do ii = 1, n_r_coll
	  deviation_r(ii) = abs(R_COM(ii) - bikram_cutoff_r2)
	  END do
	  mtrx_cutoff_r2 = minloc(deviation_r,1)
	  deallocate(deviation_r)
	  cut_r = .true.
	  IF(myid.eq.0) write(*,'(2(a,i0,a,f0.3))') 
     & "Truncation #1 at #R = ", mtrx_cutoff_r1, ", R = ", 
     & R_COM(mtrx_cutoff_r1), ", and #2 at #R = ",mtrx_cutoff_r2, 
     & ", R = ",R_COM(mtrx_cutoff_r2)
	  MIJ_ZERO_CUT = MIJ_ZERO_CUT/eVtown/autoeV
	  END IF
	  
	  allocate(bk_non_zero_mij_gather(nproc))
	  write(bk_matrix_path2,'(a,a)')trim(bk_dir1),
     & "/MATRIX_FILE_INDEX.DAT"
	  bk_matrix_path2 = trim(bk_matrix_path2)
	  
	  open(1,file = trim(bk_matrix_path2))
	  read(1,'(a15,i10)') dm7, fc_old1
	  
	  IF(fc_old1.ne.nproc) then
	  IF(myid.eq.0) write(*,'(a,i0,a,i0,a)') "The #procs, ", nproc, 
     &", is not equal to the #files, ", fc_old1, ", to be truncated."
	  IF(myid.eq.0) write(*,'(a,a)') "Please make sure to use equal ",
     & "#procs. Program will stop now."
	  stop
	  return
	  END IF
	  
	  read(1,*)
	  IF(myid.eq.0) then
	  read(1,'(3(i16,2x))') i, k_st, k_fn
	  IF((i-1).ne.myid)print*,"File index number does not match",myid
	  ELSE
	  do ii = 1, myid
	  read(1,*)
	  END do
	  read(1,'(3(i16,2x))') i, k_st, k_fn
	  IF((i-1).ne.myid)print*,"File index number does not match",myid	  
	  END IF
	  close(1)
	  
	  write(bk_matrix_path2,'(a,a)')trim(bk_dir1),
     & "/MATRIX_NONZERO_INDEX.DAT"
	  bk_matrix_path2 = trim(bk_matrix_path2)
	  
	  open(1,file = bk_matrix_path2)
	  read(1,'(a9,i16)') dm7, file_counter
	  read(1,'(a28,i16)') dm5, bk_non_zero_counter_old
	  read(1,'(a21,e19.12)') dm6, MIJ_ZERO_CUT_old
	  read(1,*)
	  allocate(file_old(2,file_counter))
	  do i = 1, file_counter
	  read(1,'(i16,2x,i16)')file_old(1,i),file_old(2,i)
	  END do
	  close(1)
	  IF(file_old(1,myid+1).ne.myid+1) then
	  print*, "Error reading MATRIX_NONZERO_INDEX.DAT"
	  stop
	  END IF
	  
	  write(bk_matrix_path2,'(a,a,i0,a,i0,a)') trim(bk_dir1), 
     & "/MIJ_", k_st, "_", k_fn, ".DAT"
	  bk_matrix_path2 = trim(bk_matrix_path2)
	  write(bk_matrix_path3,'(a,a,i0,a,i0,a)') trim(bk_dir1), 
     & "/Ind_", k_st, "_", k_fn, ".DAT"
	  bk_matrix_path3 = trim(bk_matrix_path3)
	  write(bk_matrix_path4,'(a,a,i0,a,i0,a)') trim(bk_dir2), 
     & "/MIJ_", k_st, "_", k_fn, ".DAT"
	  bk_matrix_path4 = trim(bk_matrix_path4)
	  write(bk_matrix_path5,'(a,a,i0,a,i0,a)') trim(bk_dir2), 
     & "/Ind_", k_st, "_", k_fn, ".DAT"
	  bk_matrix_path5 = trim(bk_matrix_path5)
	  
	  allocate(bk_mat_array(n_r_coll))
	  bk_non_zero_counter = 0
	  
	  IF(.not.unformat_defined) then
	  open(11,file = bk_matrix_path2)
	  open(13,file = bk_matrix_path3)
	  open(12,file = bk_matrix_path4)
	  open(14,file = bk_matrix_path5)
	  do i = 1, file_old(2,myid+1)!k_fn-k_st+1
      do ii=1,n_r_coll	  
      read(11,'(i16,1x,i8,1x,e19.12)')dmm1,dmm2,
     & bk_mat_array(ii)
      END do
	  read(13,*) dmm3(1), dmm3(2)
	  
	  IF(max(abs(bk_mat_array(mtrx_cutoff_r1)),
     & abs(bk_mat_array(mtrx_cutoff_r2))).gt.MIJ_ZERO_CUT) then
	  bk_non_zero_counter = bk_non_zero_counter + 1
      do ii=1,n_r_coll	  
      write(12,'(i16,1x,i8,1x,e19.12)') dmm1,ii,
     & bk_mat_array(ii)
	  END do
	  write(14,*) dmm3(1), dmm3(2)
      END IF
      END do
      close(11)
      close(13)
      close(12)
      close(14)
	  
	  ELSE
	  open(11,file = bk_matrix_path2, form = "unformatted")
	  open(13,file = bk_matrix_path3, form = "unformatted")
	  open(12,file = bk_matrix_path4, form = "unformatted")
	  open(14,file = bk_matrix_path5, form = "unformatted")
	  do i = 1, file_old(2,myid+1) !k_fn-k_st+1
      read(11) dmm1
	  read(11) bk_mat_array(:)
	  read(13) dmm3(:)
	  
!	  IF(abs(bk_mat_array(1)).lt.MIJ_ZERO_CUT) then
	  IF(max(abs(bk_mat_array(mtrx_cutoff_r1)),
     & abs(bk_mat_array(mtrx_cutoff_r2))).gt.MIJ_ZERO_CUT) then
!	  print*, myid, i, bk_mat_array(mtrx_cutoff_r1),
!     & bk_mat_array(mtrx_cutoff_r2), MIJ_ZERO_CUT
	  bk_non_zero_counter = bk_non_zero_counter + 1
      WRITE(12) dmm1
      WRITE(12) bk_mat_array(:)
	  write(14)dmm3(:)
      END IF
      END do
      close(11)
      close(13)
      close(12)
      close(14)
	  
	  END IF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	 
	  
      CALL MPI_GATHER(bk_non_zero_counter,1,
     & MPI_INTEGER,bk_non_zero_mij_gather,1,MPI_INTEGER,0,
     & MPI_COMM_WORLD,ierr_mpi)
	  IF(myid.eq.0) then
	  bk_non_zero_counter = 0
	  do i = 1, nproc
	  bk_non_zero_counter = bk_non_zero_counter + 
     & bk_non_zero_mij_gather(i)
	  END do
	  
	  write(bk_matrix_path4,'(a,a)') trim(bk_dir2),
     & "/MATRIX_NONZERO_INDEX.DAT"
	  bk_matrix_path4 = trim(bk_matrix_path4)
	  
	  open(11,file=bk_matrix_path4)
	  write(11,'(a9,i16)') "#Files = ", nproc
	  write(11,'(a28,i16)') "Non-Zero #Matrix Elements = ", 
     & bk_non_zero_counter
	  write(11,'(a21,e19.12,a6)') "The Matrix Cut-Off = ", 
     & MIJ_ZERO_CUT*eVtown*autoeV, " cm^-1"
	  write(11,'(a16,2x,a16)') '          file_#', '   #non_zero_Mij'
	  do i = 1, nproc
	  write(11,'(i16,2x,i16)') i, bk_non_zero_mij_gather(i)
	  END do
	  close(11)
	  
	  write(bk_matrix_path4,'(4(a))') trim(bk_dir1),
     & "/MATRIX_FILE_INDEX.DAT ", trim(bk_dir2),
     & "/MATRIX_FILE_INDEX.DAT"
	  bk_matrix_path4 = trim(bk_matrix_path4)
	  call system ( "cp " // trim(bk_matrix_path4) )
	  
	  write(bk_matrix_path4,'(4(a))') trim(bk_dir1),
     & "/MATRIX_INFO.DAT ", trim(bk_dir2), "/MATRIX_INFO.DAT"
	  bk_matrix_path4 = trim(bk_matrix_path4)
	  call system ( "cp " // trim(bk_matrix_path4) )
	  
	  print*, "Matrix has been truncated. Program will stop now."
	  END IF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	 
	  stop
	  return	  
	  END IF
	  
!! Bikram End.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!-------------------------------------------------------------------------------------------------	  
!! Bikram Start Dec 2020:

	  IF (bikram_mij_multiprint .and. .not.bikram_truncate_MTRX) then
	  IF(allocated(Mat_el)) deallocate(Mat_el)
	  IF(allocated(Mat_el_der)) deallocate(Mat_el_der)
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	 

	  IF (myid.eq.0) then
      call bk_read_matrix_info
      END IF
      CALL MPI_BCAST(R_COM, n_r_coll, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
	  
      call MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      TIME_MAT_FINISH = MPI_Wtime()
      TIME_1_MAT = TIME_MAT_FINISH - TIME_MAT_START	  
      CALL MPI_BCAST(CRITICAL_ERROR, 1, MPI_LOGICAL,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(CRITICAL_ERROR) THEN
      IF(MYID.EQ.0) PRINT*,"CRITICAL ERROR IN MATRIX READING"
      STOP 	  
      ENDIF	  
      CALL MPI_BCAST(total_size_old, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	

	  write(bk_matrix_path2,'(a,a)')trim(bk_dir2),
     & "/MATRIX_NONZERO_INDEX.DAT"
	  bk_matrix_path2 = trim(bk_matrix_path2)
	  
	  open(1,file = bk_matrix_path2)
	  read(1,'(a9,i16)') dm7, file_counter
	  read(1,'(a28,i16)') dm5, bk_non_zero_counter_old
	  read(1,'(a21,e19.12)') dm6, MIJ_ZERO_CUT_old
	  read(1,*)
	  allocate(file_old(2,file_counter))
	  do i = 1, file_counter
	  read(1,'(i16,2x,i16)')file_old(1,i),file_old(2,i)
	  END do
	  close(1)
	  IF(myid.eq.0) print*,"Current size of the non-zero Matrix = ",
     & bk_non_zero_counter_old
	  IF(MIJ_ZERO_CUT_old.ne.MIJ_ZERO_CUT) then
	  print*, "The truncated Matrix does not match the Cut-Off Value."
	  print*, "Please check and provide correct truncated Matrix."
	  print*, MIJ_ZERO_CUT_old, MIJ_ZERO_CUT
	  stop
	  return
	  END IF
	  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This piece is to add additional matrix elements, 
! need to work on this for the case of Parallel_I/O
	  IF(total_size.gt.total_size_old) then
	  bk_mij_addition = total_size - total_size_old
	  IF(bk_mij_addition.gt.nproc) then
	  
      chunk_mpi_size = bk_mij_addition/nproc
      k_st_mpi = myid*chunk_mpi_size + 1
      k_fn_mpi = (myid+1)*chunk_mpi_size 	 
	  mij_remander = 0
	  IF(bk_mij_addition.gt.(chunk_mpi_size*nproc)) 
     & mij_remander = bk_mij_addition - chunk_mpi_size*nproc
	  IF(myid.lt.mij_remander) then	 
	  k_st_mpi = k_st_mpi + myid + total_size_old
      k_fn_mpi = k_fn_mpi + (myid + 1) + total_size_old
	  ELSE	  
      k_st_mpi = k_st_mpi + mij_remander + total_size_old
      k_fn_mpi = k_fn_mpi + mij_remander + total_size_old
	  END IF
	  
!!!!  COMPUTING MATRIX	 
      IF(MYID.EQ.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "COMPUTING ADDITIONAL MATRIX ELEMENTS STARTED"	  
	  
	  bk_mat_counter = 0
	  bk_non_zero_counter = 0
	  allocate(bk_mat_temp1(n_r_coll,1000))
	  allocate(bk_non_zero_mij_gather(nproc))
	  
! finding #r for matrix truncation
	  IF(.not.cut_r) then
      CALL	INTEGRATOR_TEMP(intgeral,k_st_mpi,1,cyc)
	  allocate(deviation_r(n_r_coll))
	  do ii = 1, n_r_coll
	  deviation_r(ii) = abs(R_COM(ii) - bikram_cutoff_r1)
	  END do
	  mtrx_cutoff_r1 = minloc(deviation_r,1)
	  deallocate(deviation_r)
	  allocate(deviation_r(n_r_coll))
	  do ii = 1, n_r_coll
	  deviation_r(ii) = abs(R_COM(ii) - bikram_cutoff_r2)
	  END do
	  mtrx_cutoff_r2 = minloc(deviation_r,1)
	  deallocate(deviation_r)
	  cut_r = .true.
	  IF(myid.eq.0) write(*,'(2(a,i0,a,f0.3))') 
     & "Truncation #1 at #R = ", mtrx_cutoff_r1, ", R = ", 
     & R_COM(mtrx_cutoff_r1), ", and #2 at #R = ",mtrx_cutoff_r2, 
     & ", R = ",R_COM(mtrx_cutoff_r2)
	  MIJ_ZERO_CUT = MIJ_ZERO_CUT/eVtown/autoeV
	  END IF

! Bikram Start: this is to print progress of matrix computation
	  percent_counter = 1
	  bk_tym1 = MPI_Wtime()
      DO  k=k_st_mpi,k_fn_mpi
	  
	  IF((k_fn_mpi - k_st_mpi) .gt. 10) then
	  IF(mod((k - k_st_mpi), int((k_fn_mpi - k_st_mpi)/10)) == 0) then
	  bk_tym2 = MPI_Wtime()
	  bk_tym = bk_tym2 - bk_tym1
	  write(*,'(2(a, i5),a,f12.3)') "proc_id = ",myid,", Progress = ",
     & int(dble(k - k_st_mpi)/dble(k_fn_mpi - k_st_mpi)*100.d0), 
     & "%, Time(sec.) = ", bk_tym
	  percent_counter = percent_counter + 1
	  END IF
	  END IF
	  IF(k == k_fn_mpi) then
	  bk_tym2 = MPI_Wtime()
	  bk_tym = bk_tym2 - bk_tym1
	  write(*,'(2(a, i5),a,f12.3)') "proc_id = ",myid,", Progress = ", 
     & int(100.d0), "%, Time(sec.) = ", bk_tym
	  END IF
! Bikram End.
	  
	  bk_mat_counter = bk_mat_counter + 1
	  IF(bk_mat_counter.eq.1) mt_chk = k
	  
	  mij_k_skip = .false.
	  
! calling matrix calculation for 2 values of R to check 
! whether to compute the entrire array along R
      IF (.NOT. expansion_defined) THEN
!       Cutoff check first
      CALL	INTEGRATOR_TEMP(intgeral,k,mtrx_cutoff_r1,cyc)
	   mtrx_cutoff_chk1 = intgeral*conv_unit_e
		
      CALL	INTEGRATOR_TEMP(intgeral,k,mtrx_cutoff_r2,cyc)
	   mtrx_cutoff_chk2 = intgeral*conv_unit_e

	  IF(max(abs(mtrx_cutoff_chk1), abs(mtrx_cutoff_chk2)).lt.
     & MIJ_ZERO_CUT) then
      bk_mat_temp1(:,bk_mat_counter) = 0d0
      mij_k_skip = .TRUE.
	   bk_non_zero_counter = bk_non_zero_counter + 1
      ENDIF	  
	  
      IF(.not.mij_k_skip) then  
      DO i=1,n_r_coll
	   IF(mij_k_skip) CYCLE
	   IF(i.eq.mtrx_cutoff_r1) then
	   bk_mat_temp1(i,bk_mat_counter) = mtrx_cutoff_chk1
	   ELSE IF(i.eq.mtrx_cutoff_r2) then
	   bk_mat_temp1(i,bk_mat_counter) = mtrx_cutoff_chk2
	   ELSE
      CALL INTEGRATOR_TEMP(intgeral,k,i,cyc)	
	   bk_mat_temp1(i,bk_mat_counter) = intgeral*conv_unit_e
	   END IF
      ENDDO
	   END IF
		
		ELSE
! Expansion: compute all R at once, no cutoff check needed
		ALLOCATE(Mat_el_R_array(n_r_coll))
		CALL EXPANSION_MATRIX_ELEMENT(Mat_el_R_array,k,cyc)
		bk_mat_temp1(:,bk_mat_counter) = Mat_el_R_array(:)*conv_unit_e
		DEALLOCATE(Mat_el_R_array)
		
		IF(max(abs(bk_mat_temp1(mtrx_cutoff_r1,bk_mat_counter)),
     &       abs(bk_mat_temp1(mtrx_cutoff_r2,bk_mat_counter))).lt.
     &       MIJ_ZERO_CUT) then
		bk_mat_temp1(:,bk_mat_counter) = 0d0
		mij_k_skip = .TRUE.
		bk_non_zero_counter = bk_non_zero_counter + 1
		ENDIF
		ENDIF
	  
	   IF(bk_mat_counter.eq.1000 .or. k.eq.k_fn_mpi) then
	   call bk_print_matrix (k_st_mpi, k_fn_mpi, k,
     & bk_mat_counter, bk_mat_temp1, mt_chk, cyc_cntr)
	   bk_mat_counter = 0
	   ENDIF
      ENDDO
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	    
      CALL MPI_GATHER(k_fn_mpi-k_st_mpi+1-bk_non_zero_counter,1,
     & MPI_INTEGER,bk_non_zero_mij_gather,1,MPI_INTEGER,0,
     & MPI_COMM_WORLD,ierr_mpi)
	  IF(myid.eq.0) then
	  bk_non_zero_counter = 0
	  do i = 1, nproc
	  bk_non_zero_counter = bk_non_zero_counter + 
     & bk_non_zero_mij_gather(i)
	  END do
	  open(11,file = bk_matrix_path2)
	  write(11,'(a9,i16)') "#Files = ", nproc + file_counter
	  write(11,'(a28,i16)') "Non-Zero #Matrix Elements = ", 
     & bk_non_zero_counter + bk_non_zero_counter_old
	  write(11,'(a21,e19.12,a6)') "The Matrix Cut-Off = ", 
     & MIJ_ZERO_CUT*eVtown*autoeV, " cm^-1"
	  write(11,'(a16,2x,a16)') '          file_#', '   #non_zero_Mij'
	  do i = 1, file_counter
	  write(11,'(i16,2x,i16)') file_old(1,i),file_old(2,i)
	  IF(file_old(1,i).ne.i) stop "Error in additional matrix print"
	  END do
	  do i = 1, nproc
	  write(11,'(i16,2x,i16)')i+file_counter,bk_non_zero_mij_gather(i)
	  END do
	  close(11)
	  
	  write(bk_matrix_path2,'(a,a)')trim(bk_dir2),
     & "/MATRIX_FILE_INDEX.DAT"
	  bk_matrix_path2 = trim(bk_matrix_path2)
!	  write(bk_matrix_path4,'(a,a)')trim(bk_dir2),
!     & "/MATRIX_COMBINE.sh"
	  bk_matrix_path4 = trim(bk_matrix_path4)
	  
	  allocate(file_old1(3,file_counter))
	  open(1,file = trim(bk_matrix_path2))
	  read(1,*)
	  read(1,*)
	  do i = 1, file_counter
	  read(1,'(3(i16,2x))')
     & file_old1(1,i),file_old1(2,i),file_old1(3,i)
	  END do
	  close(1)
	  open(11,file = trim(bk_matrix_path2))
!	  open(111,file = trim(bk_matrix_path4), access = "appEND")
	  write(11,'(a15,i10)')"Total #files = ", nproc + file_counter
	  write(11,'(3(a16,2x))') '            #Mij', '      #Mij_begin'
     & ,'        #Mij_END'
	  do i = 1, file_counter
	  write(11,'(3(i16,2x))')
     & file_old1(1,i),file_old1(2,i),file_old1(3,i)
	  END do
	  do i = 1, nproc
      k_st_mpi = (i-1)*chunk_mpi_size + 1
      k_fn_mpi = i*chunk_mpi_size 	 
	  mij_remander = 0
	  IF(bk_mij_addition.gt.(chunk_mpi_size*nproc)) 
     & mij_remander = bk_mij_addition - chunk_mpi_size*nproc
	  IF((i-1).lt.mij_remander) then	 
	  k_st_mpi = k_st_mpi + i-1 + total_size_old
      k_fn_mpi = k_fn_mpi + i + total_size_old
	  ELSE	  
      k_st_mpi = k_st_mpi + mij_remander + total_size_old
      k_fn_mpi = k_fn_mpi + mij_remander + total_size_old
	  END IF
	  write(11,'(3(i16,2x))') i+file_counter, k_st_mpi, k_fn_mpi
	  
!	  write(111, '(a,i0,a,i0,a)') 
!     & 'cat MIJ_',k_st_mpi,'_',k_fn_mpi,'.DAT >> ../MTRX.DAT'
	  END do
	  close(11)
!	  close(111)
	  END IF  
	  
	  ELSE	  
	  IF(myid.eq.0) then
	  write(*,'(a)')"The #procs is larger than #elements in matrix."
	  write(*,'(a)')"It is not optimal to have many files with 1 Mij."
	  write(*,'(a)')"Program will stop now."
	  END IF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	 
	  stop
	  return
	  END IF
	  
	  TIME_1_MAT = MPI_Wtime()
	  IF(myid.eq.0) call bk_print_matrix_info
	  
	  TIME_2_MAT = MPI_Wtime()
      TIME_WRITING_MATRIX = TIME_2_MAT - TIME_1_MAT		  
!!! WAITING FOR OTHER PROCESSORS
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      TIME_MAT_FINISH= MPI_Wtime()
      IF(MYID.EQ.0) WRITE(*,'(a57,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij READING/COMPUTING/SAVING ON DISK"
     & ,",s",(TIME_MAT_FINISH-TIME_MAT_START)/nproc
	  IF(print_matrix_defined.and. myid.eq.0) 
     & WRITE(*,'(a47,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij SAVING ON DISK"
     & ,",s",TIME_WRITING_MATRIX	 
      IF(MYID.EQ.0) PRINT*, "ALL WORK ON MATRIX IS DONE"
      IF(run_prog_defined) THEN
	  IF(myid.eq.0) then
	  print*,"You indicated Parallel_I/O."
	  print*,"Trajs can't be computed in the same run with Matrix."
	  print*,"Please run again. Program will stop now."
	  END IF
      ENDIF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      STOP	  
	  return
	  END IF
	  	  
	  total_size = bk_non_zero_counter_old
	  
	  IF(myid.eq.0) PRINT*,"MPI TASK PER TRAJECTORY WILL BE USED"
      IF(myid.eq.0)	
     & WRITE(*,'(a53,1x,i4)')
     & "MPI TASKS WHICH ARE ASSOCIATED WITH ONE TRAJECTORY = ",
     & mpi_task_per_traject	 
      traject_roots = nproc/mpi_task_per_traject
      IF(traject_roots*mpi_task_per_traject.ne.nproc) THEN
      STOP "mpi_task_number must be a delimiter of nproc"	  
      ENDIF	  
      ALLOCATE(mpi_traject_roots(traject_roots))
      ALLOCATE(portion_of_MIJ_per_task(2,nproc))
      ALLOCATE(portion_of_state_per_task(2,nproc))
      ALLOCATE(portion_of_work_per_task(2,nproc))	  
      ALLOCATE(mpi_root_belongs(nproc))	  
!!!!!!!!!!!!!   REMAKE	  
      size_mij_chunk_mpi = total_size/mpi_task_per_traject
      residue_mij_mpi = total_size-
     & size_mij_chunk_mpi*mpi_task_per_traject
      size_state_chunk_mpi = states_size/mpi_task_per_traject
      residue_state_mpi = states_size-
     & size_state_chunk_mpi*mpi_task_per_traject
      size_work_chunk_mpi = (2*states_size+8)/mpi_task_per_traject !!!! DO IN FUTURE
      residue_work_mpi = 2*states_size+8-!!!! DO IN FUTURE
     & size_work_chunk_mpi*mpi_task_per_traject	!!!! DO IN FUTURE
	 
      DO k=1,traject_roots
      mpi_traject_roots(k) = (k-1)*mpi_task_per_traject
      DO k_p=1,mpi_task_per_traject	  
      mpi_root_belongs(k_p+mpi_traject_roots(k))=mpi_traject_roots(k)
      ENDDO	  
      ENDDO
      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_mij_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_MIJ_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_mij_chunk_mpi+1)
      portion_of_MIJ_per_task(2,k_mpi_proc) =
     & (k_p)*(size_mij_chunk_mpi+1)
      ENDDO
      total_size_check = total_size_check + size_mij_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_MIJ_per_task(1,k_mpi_proc)
     & = residue_mij_mpi*(size_mij_chunk_mpi+1)+1  
     & + (k_p-1-residue_mij_mpi)*size_mij_chunk_mpi
      portion_of_MIJ_per_task(2,k_mpi_proc)
     & = residue_mij_mpi*(size_mij_chunk_mpi+1)+	  
     & (k_p-residue_mij_mpi)*size_mij_chunk_mpi
      ENDDO	 
      total_size_check = total_size_check + size_mij_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_MIJ_per_task(2,k_mpi_proc).ne.total_size) 
     & STOP "WRONG MATIX ASSIGNEMENT"  
      ENDIF	  
      ENDIF
	  
      ENDDO	
      IF(total_size_check.ne.total_size) THEN
      IF(myid.eq.0) WRITE(*,*)total_size_check,total_size	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF	 

      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_state_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_state_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_state_chunk_mpi+1)
      portion_of_state_per_task(2,k_mpi_proc) =
     & (k_p)*(size_state_chunk_mpi+1)
      ENDDO
      state_size_check = state_size_check + size_state_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_state_per_task(1,k_mpi_proc)
     & = residue_state_mpi*(size_state_chunk_mpi+1)+1  
     & + (k_p-1-residue_state_mpi)*size_state_chunk_mpi
      portion_of_state_per_task(2,k_mpi_proc)
     & = residue_state_mpi*(size_state_chunk_mpi+1)+	  
     & (k_p-residue_state_mpi)*size_state_chunk_mpi
      ENDDO	 
      state_size_check = state_size_check + size_state_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_state_per_task(2,k_mpi_proc).ne.states_size) 
     & STOP "WRONG MATRIX ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO	  
	  
!      PRINT*,"STATES",portion_of_state_per_task	  
      IF(state_size_check.ne.states_size) THEN
      IF(myid.eq.0) WRITE(*,*)state_size_check,states_size	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF	

      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_work_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_work_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_work_chunk_mpi+1)
      portion_of_work_per_task(2,k_mpi_proc) =
     & (k_p)*(size_work_chunk_mpi+1)
      ENDDO
      work_size_check = work_size_check + size_work_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_work_per_task(1,k_mpi_proc)
     & = residue_work_mpi*(size_work_chunk_mpi+1)+1  
     & + (k_p-1-residue_work_mpi)*size_work_chunk_mpi
      portion_of_work_per_task(2,k_mpi_proc)
     & = residue_work_mpi*(size_work_chunk_mpi+1)+	  
     & (k_p-residue_work_mpi)*size_work_chunk_mpi
      ENDDO	 
      work_size_check = work_size_check + size_work_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_work_per_task(2,k_mpi_proc).ne.states_size*2+8) 
     & STOP "WRONG MATRIX ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO	  
	  
      IF(work_size_check.ne.states_size*2+8) THEN
      IF(myid.eq.0) WRITE(*,*)work_size_check,states_size*2+8	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF
	  
      total_size_mpi = portion_of_MIJ_per_task(2,myid+1) - 
     & portion_of_MIJ_per_task(1,myid+1) + 1	  
      ALLOCATE(Mat_el(n_r_coll,total_size_mpi))
      ALLOCATE(Mat_el_der(n_r_coll,total_size_mpi))
      CALL MPI_Comm_group(MPI_COMM_WORLD, wrld_group,ierr_mpi)
      ALLOCATE(process_rank_distr(mpi_task_per_traject,traject_roots))
      ALLOCATE(comms_distr(mpi_task_per_traject),
     & groups_distr(mpi_task_per_traject))   	  
      DO i=1,traject_roots	  
      DO k=1,mpi_task_per_traject
      process_rank_distr(k,i) = k - 1 + (i-1)*mpi_task_per_traject     	  
      ENDDO
      ENDDO	

      DO i=1,mpi_task_per_traject      	  
      CALL MPI_Group_incl(wrld_group, traject_roots,
     & process_rank_distr(i,:), groups_distr(i),ierr_mpi)
      CALL MPI_Comm_create(MPI_COMM_WORLD,groups_distr(i),
     & comms_distr(i),ierr_mpi)
      ENDDO	  
!      PRINT*,myid,total_size_mpi	  
!      tag1 = 1
!      tag2 = 2
      IF(MYID.eq.0) PRINT*,"Mij ROOT PROC DISTRIBUTION STARTED"	  
!      IF(MYID.eq.0) THEN
!      DO k=2,mpi_task_per_traject
!      total_size_mpi= portion_of_MIJ_per_task(2,k) - 
!     & portion_of_MIJ_per_task(1,k) + 1	
!      task_portion_size = total_size_mpi*n_r_coll	 
!      ALLOCATE(buffer_mpi_portion(n_r_coll,total_size_mpi))
!      buffer_mpi_portion = Mat_el_non_zero(:,
!     & portion_of_MIJ_per_task(1,k):portion_of_MIJ_per_task(2,k)) 	  
!      CALL MPI_SEND(buffer_mpi_portion,
!     & task_portion_size, MPI_REAL8, k-1, 
!     &  tag1, MPI_COMM_WORLD, ierr_mpi)
!      buffer_mpi_portion = 	 Mat_el_non_zero_der (:,
!     & portion_of_MIJ_per_task(1,k):portion_of_MIJ_per_task(2,k))
!      CALL MPI_SEND(buffer_mpi_portion,
!     & task_portion_size, MPI_REAL8, k-1, 
!     &  tag2, MPI_COMM_WORLD, ierr_mpi)
!      DEALLOCATE(buffer_mpi_portion) 
!      ENDDO
!      total_size_mpi= portion_of_MIJ_per_task(2,1) - 
!     & portion_of_MIJ_per_task(1,1) + 1		  
!      Mat_el=Mat_el_non_zero(:,
!     & portion_of_MIJ_per_task(1,1):portion_of_MIJ_per_task(2,1))
!      Mat_el_der=Mat_el_non_zero_der(:,
!     & portion_of_MIJ_per_task(1,1):portion_of_MIJ_per_task(2,1))
!      DEALLOCATE(Mat_el_non_zero,Mat_el_non_zero_der)	 
!      ELSE
!      IF(myid.le.mpi_task_per_traject-1) THEN	  
!      task_portion_size = total_size_mpi*n_r_coll	  
!      CALL MPI_RECV(Mat_el, task_portion_size, MPI_REAL8, 
!     & 0, tag1, MPI_COMM_WORLD, status, ierr_mpi)	  
!      CALL MPI_RECV(Mat_el_der, task_portion_size, MPI_REAL8, 
!     & 0, tag2, MPI_COMM_WORLD, status, ierr_mpi)
!      ENDIF	 
!      ENDIF
!	  
!	  IF(myid.eq.0) call bk_matrix_splitting
!	  call MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
	  
!	  allocate(K_SKIPPED_BY_ROUTINE(total_size_mpi))
!	  IF(allocated(ind_mat)) deallocate(ind_mat)
!	  allocate(ind_mat(2,total_size_mpi))
!	  K_SKIPPED_BY_ROUTINE = .false.
	  allocate(bk_ind_tmp(2,total_size_mpi))
	  call bk_read_matrix
     & (total_size_mpi, Mat_el, bk_ind_tmp)
!     & (total_size_mpi, Mat_el, ind_mat)
!	  Mat_el(:,k) = bk_mat_temp(:)
	  IF(allocated(ind_mat)) deallocate(ind_mat)
	  allocate(ind_mat(2,total_size_mpi))
	  ind_mat(:,:) = bk_ind_tmp(:,:)
!	  mij_counter = 0
!	  do k = 1, total_size_mpi
!	  IF(ABS(Mat_el(1,k)).LT.MIJ_ZERO_CUT) THEN
!      Mat_el_der(:,k) = 0d0	  
!      K_SKIPPED_BY_ROUTINE(k) = .TRUE.	  
!	  mij_counter = mij_counter + 1
!      ENDIF
!	  END do
!	  print*,myid, '#Non_zero_Mij = ', total_size_mpi - mij_counter
!	  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
	  IF(myid.eq.0) print*,"Mij_Splining_Started"
	  
	  DO  k = 1, total_size_mpi
!      IF(K_SKIPPED_BY_ROUTINE(k)) CYCLE
      deriv_bgn = (Mat_el(2,k) - Mat_el(1,k))/
     & (R_COM(2)-R_COM(1)) 
      deriv_END = (Mat_el(n_r_coll,k) - Mat_el(n_r_coll-1,k))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1))
      CALL spline(R_COM,Mat_el(:,k),
     & n_r_coll,deriv_bgn,deriv_END,Mat_el_der(:,k))	 
      ENDDO
	  
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
	  IF(myid.eq.0) print*,"Mij_Splining_finished"
	  

      IF(myid.eq.0)	 PRINT*, "Mij ROOT PROC DISTRIBUTION DONE"  
      IF(myid.eq.0) PRINT*,"Mij ALL PROC DISTRIBUTION STARTED" 	  
      DO i=1,mpi_task_per_traject
      IF(i-1 .eq. myid - int(myid/mpi_task_per_traject)
     & *mpi_task_per_traject) THEN
      CALL MPI_Comm_rank(comms_distr(i),id_proc_in_group ,ierr_mpi)
      IF(id_proc_in_group.ne.myid/mpi_task_per_traject) PRINT*,
     & "ERROR IN COMUNICATIONS ASSIGNEMNET_1_distr",id_proc_in_group,
     & i-1
      IF(myid.lt.mpi_task_per_traject .and. id_proc_in_group.ne.0 )
     & PRINT*,"ERROR IN COMUNICATIONS ASSIGNEMNET_2_distr"
      task_portion_size = total_size_mpi*n_r_coll	  
      CALL MPI_BCAST(Mat_el,task_portion_size, MPI_REAL8,0,
     &  comms_distr(i),ierr_mpi)
      CALL MPI_BCAST(Mat_el_der,task_portion_size, MPI_REAL8,0,
     &  comms_distr(i),ierr_mpi)
      ENDIF	 
      ENDDO
      IF(myid.eq.0) PRINT*,"Mij ALL PROC DISTRIBUTION DONE"	  
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)	  
      IF(MYID.EQ.0) PRINT *,"MPI_COMMUNICATORS_CREATION STARTED"
      ALLOCATE(comms(traject_roots),groups(traject_roots))     	  
      ALLOCATE(process_rank(traject_roots,mpi_task_per_traject))		  
      DO i=1,traject_roots	  
      DO k=1,mpi_task_per_traject
      process_rank(i,k) = k - 1 + (i-1)*mpi_task_per_traject     	  
      ENDDO
      ENDDO	
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      DO i=1,traject_roots   
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      CALL MPI_Group_incl(wrld_group, mpi_task_per_traject,
     & process_rank(i,:), groups(i),ierr_mpi)
      CALL MPI_Comm_create(MPI_COMM_WORLD,groups(i),comms(i),ierr_mpi)
      ENDDO
      DO i=1,traject_roots
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN
      CALL MPI_Comm_rank(comms(i),id_proc_in_group ,ierr_mpi)	  
      IF(mpi_traject_roots(i)+id_proc_in_group.ne.myid) STOP
     & "WRONG GROUP MPI ASSIGNEMENT"	  
	  
      ENDIF		  
      ENDDO	 
      IF(MYID.EQ.0) PRINT *,"MPI_COMMUNICATORS_CREATION DONE"
      TIME_MAT_FINISH = MPI_Wtime()
      TIME_2_MAT = TIME_MAT_FINISH - TIME_MAT_START	  
	  
!	  ENDIF
	  return
	  	  
	  END IF
!! Bikram End.
!!-------------------------------------------------------------------------------------------------	  
	  
	  
	  
	  	  
!!!   READING Mij FROM A FILE	  
!      TIME_MAT_START = MPI_Wtime()	  
      IF (myid.eq.0) THEN
      IF(.not.unformat_defined) THEN
!	  bgn_tym = MPI_Wtime()											!Bikram
      CALL READ_MIJ !!!! TESTING
!	  END_tym = MPI_Wtime()											!Bikram
!	  calc_tym = END_tym - bgn_tym									!Bikram
!	  tot_tym_rd = tot_tym_rd + calc_tym							!Bikram
      ELSE 
      CALL READ_MIJ_USER	  
      ENDIF	
!	  print *, 'rd_alex', myid, tot_tym_rd							!Bikram
      ENDIF
!      STOP	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      TIME_MAT_FINISH = MPI_Wtime()
      TIME_1_MAT = TIME_MAT_FINISH - TIME_MAT_START	  
      CALL MPI_BCAST(CRITICAL_ERROR, 1, MPI_LOGICAL,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(CRITICAL_ERROR) THEN
      IF(MYID.EQ.0) PRINT*,"CRITICAL ERROR IN MATRIX READING"
      STOP 	  
      ENDIF	  
      CALL MPI_BCAST(total_size_old, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
!      IF(MYID.EQ.0) PRINT*,"TOTAL_OLD_SIZE",total_size_old,total_size
!!!   COMPUTING OF NEW MATRIX BEGINS
!!!        EXPANDING ON R_GRID
      IF(total_size.eq.total_size_old) THEN
      IF(ir_bgn_exp_pnt.gt.0 .and. ir_fin_exp_pnt.gt.0 .and.
     &  .not. (ir_bgn_exp_pnt.eq.1 .and. ir_fin_exp_pnt.eq.
     & n_r_coll ) ) THEN
      IF(myid.eq.0) THEN	  
      ALLOCATE(Mat_el_r_temp(n_r_coll,total_size))
      ALLOCATE(Mat_el_der_r_temp(n_r_coll,total_size))
      ENDIF	  
!!!!  COMPUTING MATRIX	  
      IF(myid.eq.0) PRINT*,"MATRIX_INI__R_ADD_STARTED"
!      PRINT*, "HELLO"!DELETE	  
      IF(chunk_mpi_size.ge.0) THEN
!! ALLOCATING CHUNK_SIZE ARRAYS
      IF(chunk_mpi_size.gt.0) THEN	  
      ALLOCATE(Mat_el_temp(n_r_coll,k_fn_mpi-k_st_mpi +1)
     &, Mat_el_der_temp(n_r_coll,k_fn_mpi-k_st_mpi +1))
      ALLOCATE(K_SKIPPED_BY_ROUTINE(k_fn_mpi-k_st_mpi +1))
      DO k=k_st_mpi,k_fn_mpi
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.FALSE.	  
      ENDDO
      ENDIF	  
	  
!!!!! ALLOCATION OF SKIPPING ARRAYS
!!!! THIS IS THE CASE WHEN FOR LARGE SYSTEMS IT IS NOT POSSIBLE TO
!! OR PRACTICAL TO STORE IN THE MEMORY THE ENTIRE GRID
	  
      IF(make_grid_file .and. (coll_type.gt.5 .or. coll_type.eq.0))
     & THEN
      IF(MYID.EQ.0) THEN
      PRINT*, "MATRIX WILL BE COMPUTED FROM GRID"	  
      SELECT CASE(coll_type)
      CASE(6)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) n_beta1,n_beta2,n_gamma1,n_r_vib1,
     & n_r_vib2, nr_cold
      ALLOCATE(V_VIB_2_2(n_gamma1,n_beta2,n_beta1,n_r_vib2,
     & n_r_vib1,nr_cold)) 
      CASE(7)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) n_alpha2,n_beta1,n_beta2,n_gamma1,nr_cold 
      ALLOCATE(
     & V_3_2(n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))		  
      CASE(8)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) n_alpha2,n_beta1,n_beta2,n_gamma1,nr_cold 
      ALLOCATE(
     & V_3_2(n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  
      CASE(9)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat)
     &	n_alpha2,n_beta1,n_beta2,n_gamma1,n_gamma2,nr_cold 
      ALLOCATE(
     & V_3_3(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  	  
      CASE(0)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat)
     & n_alpha2,n_beta1,n_beta2,n_gamma1,n_gamma2,nr_cold 
      ALLOCATE(
     & V_3_3(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  
      END SELECT
      PRINT*,"VIJ_BUFFER_ALLOCATED"	
!!!!! READS AND PASSE V_ij BUFFER FOR FURTHER BROADCASTING,
!!!! ONLY FOR CASE OF ASYMETRIC + ASYMETRIC	
      IF(istat.gt.0) THEN
      CRITICAL_ERROR=.TRUE.
      PRINT*,"ERROR in V_GRID"	  
      ENDIF	  
      ENDIF
      CALL MPI_BCAST(CRITICAL_ERROR, 1, MPI_LOGICAL,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(CRITICAL_ERROR) STOP	 
      CALL MPI_BCAST(n_alpha2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_beta1, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_beta2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_gamma1, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_gamma2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_r_vib1, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_r_vib2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	 
      CALL MPI_BCAST(nr_cold, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(nr_cold.lt.n_r_coll) THEN 
      IF(MYID.EQ.0) THEN
      CRITICAL_ERROR=.TRUE.
      PRINT*,"ERROR in V_GRID"
      ENDIF
      STOP	  
      ENDIF	  
	  
      SELECT CASE(coll_type)
      CASE(6)
      ALLOCATE(V_VIB_2_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_r_vib2,
     & n_r_vib1)) 
      CASE(7)
      ALLOCATE(V_3_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_alpha2))	 	  
      CASE(8)
      ALLOCATE(V_3_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      CASE(9)
      ALLOCATE(
     & V_3_3_int_buffer(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      CASE(0)
!!!  MINI BUFFERS FOR EACH R-VALUE TO AVOID MEMORY LIMITATIONS	  
      ALLOCATE(
     & V_3_3_int_buffer(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2))
!      PRINT*,"BUFFER_ALLOCATED FOR=",MYID
      END SELECT	 

      DO i=1,n_r_coll
      IF(MYID.EQ.0) THEN	  
      SELECT CASE(coll_type)
      CASE(6)

      READ(1) 
     & V_VIB_2_2_int_buffer
      V_VIB_2_2(:,:,:,:,:,i)=
     & V_VIB_2_2_int_buffer	  

      CASE(7)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_2_int_buffer
      V_3_2(:,:,:,:,i)=
     & V_3_2_int_buffer 	  
      CASE(8)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_2_int_buffer
      V_3_2(:,:,:,:,i)=
     & V_3_2_int_buffer 

      CASE(9)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_3_int_buffer
      V_3_3(:,:,:,:,:,i)=
     & V_3_3_int_buffer 
 
      CASE(0)	  

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_3_int_buffer
      V_3_3(:,:,:,:,:,i)=
     & V_3_3_int_buffer	 
      END SELECT	  
      PRINT*,"BUFFER HAS BEEN READ FOR Ri = ",i
	  
      ENDIF
      SELECT CASE(coll_type)
      CASE(6)
      buffer_size_V = n_r_vib1*n_r_vib2*n_beta1*n_beta2*n_gamma1
      CASE(7)
      buffer_size_V = n_gamma1*n_beta2*n_beta1*n_alpha2		  
      CASE(8)
      buffer_size_V = n_gamma1*n_beta2*n_beta1*n_alpha2	  
      CASE(9)
      buffer_size_V = n_gamma2*n_gamma1*n_beta2*n_beta1*n_alpha2	  
      CASE(0)	  
      buffer_size_V = n_gamma2*n_gamma1*n_beta2*n_beta1*n_alpha2
      END SELECT	  
      IF(MYID.EQ.0) PRINT*,"BUFFER OF VIJ SIZE= ",
     & buffer_size_V!,size(V_3_3_int_buffer)	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
!!!!! BROADCASTING MINI _BUFFER	  
      IF(MYID.EQ.0) PRINT*,"BROADCASTING OF VIJ_BUFFER STARTED"
!      IF(MYID.EQ.0) WRITE(3,*)	V_3_3_int_buffer
      SELECT CASE(coll_type)
      CASE(6)
      CALL MPI_BCAST(V_VIB_2_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CASE(7)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(8)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(9)
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(0)  
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      END SELECT	 
      IF(MYID.EQ.0) PRINT*,"THE BUFFER HAS BEEN BROADCASTED"
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i,cyc) !!! CALCULATING VALUE OVER INTEGRATION FOR EACH R
			
!		Write(*,*) '#5'
!      PRINT*,st_1,st_2,intgeral	  
      Mat_el_temp(i,k-k_st_mpi+1) = intgeral*conv_unit_e
!      IF(ABS(Mat_el_temp(1,k-k_st_mpi+1)).LT.MIJ_ZERO_CUT) THEN    !!! EXCLUDE SOME MATRIX ELEMENTS IF THEY ARE BELOW CERTAIN VALUE
      IF(max(abs(Mat_el_temp(mtrx_cutoff_r1,k-k_st_mpi+1)),
     & abs(Mat_el_temp(mtrx_cutoff_r2,k-k_st_mpi+1))).LT.
     & MIJ_ZERO_CUT) THEN    !!! EXCLUDE SOME MATRIX ELEMENTS IF THEY ARE BELOW CERTAIN VALUE
      Mat_el_temp(:,k-k_st_mpi+1) = 0d0
      Mat_el_der_temp(:,k-k_st_mpi+1) = 0d0	  
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.TRUE.	  
      ENDIF	  
!      PRINT*,st_1,st_2,	 intgeral 
      ENDDO
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      ENDDO
	  
      ELSE
!!!! IF NOT THE GRID FORM FILE DEFINED COMPUTE MATRIX ELEMENTS
      IF(MYID.EQ.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "COMPUTING MATRIX ELEMENTS STARTED"	  
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE !!! SKIP SOME ELEMENTS IF NESSEACRY		  
      DO i=1,n_r_coll 
      IF(i.gt.i_nr_fin .or. i.lt.i_nr_ini) CYCLE	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i,cyc)
			
		Write(*,*) '#6'
!!! COMPUTING AND STROING IN THE TEMPORARY ARRAY ASSIGNED FOR EACN ID_PROC	  
      Mat_el_temp(i,k-k_st_mpi+1) = intgeral*conv_unit_e
!      IF(ABS(Mat_el_temp(1,k-k_st_mpi+1)).LT.MIJ_ZERO_CUT) THEN !!!! THIS IS EXCLUSION CONDITION
      IF(max(abs(Mat_el_temp(mtrx_cutoff_r1,k-k_st_mpi+1)),
     & abs(Mat_el_temp(mtrx_cutoff_r2,k-k_st_mpi+1))).LT.
     & MIJ_ZERO_CUT) THEN !!!! THIS IS EXCLUSION CONDITION
      Mat_el_temp(:,k-k_st_mpi+1) = 0d0
      Mat_el_der_temp(:,k-k_st_mpi+1) = 0d0	  
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.TRUE.	  
      ENDIF	  
      ENDDO
      ENDDO
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      ENDIF	
!!!      SPLININING STARTED
      IF(MYID.EQ.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "SPLINING OF Mij STARTED"	  
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE
      deriv_bgn = (Mat_el_temp(2,k-k_st_mpi+1)
     & - Mat_el_temp(1,k-k_st_mpi+1))/
     & (R_COM(2)-R_COM(1)) 
      deriv_END = (Mat_el_temp(n_r_coll,k-k_st_mpi+1)
     & - Mat_el_temp(n_r_coll-1,k-k_st_mpi+1))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1)) 	  
      CALL  spline(R_COM,Mat_el_temp(:,k-k_st_mpi+1),
     & n_r_coll,deriv_bgn,deriv_END,
     & Mat_el_der_temp(:,k-k_st_mpi+1))	 
      ENDDO	 
      IF(MYID.EQ.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "SPLINING OF Mij FINISHED"	 
!!! GATHERING MATRIX IN Mij.dat IN ONE PROCESSOR	 

      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF(chunk_mpi_size.gt.0) THEN	  
      CALL MPI_GATHER(Mat_el_temp,task_size,MPI_DOUBLE_PRECISION,
     & Mat_el_r_temp,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 	  
      CALL MPI_GATHER(Mat_el_der_temp,task_size,MPI_DOUBLE_PRECISION,
     & Mat_el_der_r_temp,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
      ENDIF	 
      ENDIF
      IF(MYID.eq.0) THEN	  
      Mat_el(i_nr_ini:i_nr_fin,:) = Mat_el_r_temp(i_nr_ini:i_nr_fin,:)
      Mat_el_der(i_nr_ini:i_nr_fin,:) = 
     & Mat_el_der_r_temp(i_nr_ini:i_nr_fin,:) 
      DEALLOCATE(Mat_el_r_temp,Mat_el_der_r_temp)
      ENDIF	  
!!!! COMPUTING THE REST OF THE MATRIX MIJ
      IF(expansion_defined) THEN
      CALL READ_EXPANSION_TERMS	  
      ENDIF
      IF(total_size.gt.nproc*chunk_mpi_size) THEN
      IF(MYID.EQ.0) PRINT*,"COMPUTING THE RESIDUE OF Mij"
      ALLOCATE(Mat_el_res(n_r_coll),Mat_el_der_res(n_r_coll))
      ALLOCATE(Mat_el_resf(n_r_coll,nproc),
     & Mat_el_der_resf(n_r_coll,nproc))
      Mat_el_res = 0d0
      Mat_el_der_res = 0d0
      DO i=1,n_r_coll
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
      IF(K_SKIPPED_RES) CYCLE
!!! AGAIN IF READING FROM GRID DEFINED!	  
      IF(make_grid_file .and. (coll_type.gt.5 .or. coll_type.eq.0))
     & THEN
      IF(MYID.EQ.0) THEN
      SELECT CASE(coll_type)
      CASE(6)
      V_VIB_2_2_int_buffer(:,:,:,:,:) = 
     & V_VIB_2_2(:,:,:,:,:,i)
      CASE(7)
      V_3_2_int_buffer(:,:,:,:) = V_3_2(:,:,:,:,i)	  
      CASE(8)
      V_3_2_int_buffer(:,:,:,:) = V_3_2(:,:,:,:,i)	  
      CASE(0)	  
      V_3_3_int_buffer(:,:,:,:,:) = V_3_3(:,:,:,:,:,i)
      END SELECT	  
      ENDIF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      SELECT CASE(coll_type)
      CASE(6)
      CALL MPI_BCAST(V_VIB_2_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CASE(7)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(8)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(9)
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(0)	  
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      PRINT*, "CHUNK BUFFER HAS BEEN BROADCASTED FOR IR=",i	 
      END SELECT	 
      ENDIF	  
!      R_dist =  R_COM(i)*conv_unit_r
      IF(myid.lt.total_size-nproc*chunk_mpi_size) THEN
      k=nproc*chunk_mpi_size+1+myid	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i,cyc)
	  	
		Write(*,*) '#7'
!      PRINT*,"chunk test",intgeral,k,i	  
      Mat_el_res(i) = intgeral*conv_unit_e
!      IF(ABS(Mat_el_res(1)).LT.MIJ_ZERO_CUT) THEN
      IF(max(abs(Mat_el_res(mtrx_cutoff_r1)),
     & abs(Mat_el_res(mtrx_cutoff_r2))).LT.MIJ_ZERO_CUT) THEN
      Mat_el_res(:) = 0d0
      Mat_el_der_res(:) = 0d0	  
      K_SKIPPED_RES=.TRUE.	  
      ENDIF
      ENDIF	  
      ENDDO
!!! SPLINING THE REST OF MIJ
      IF(MYID.eq.0) PRINT*,"SPLINING RESIDUE OF Mij"
      IF(.NOT. K_SKIPPED_RES) THEN
      deriv_bgn = (Mat_el_res(2)
     & - Mat_el_res(1))/
     & (R_COM(2)-R_COM(1)) 
      deriv_END = (Mat_el_res(n_r_coll)
     & - Mat_el_res(n_r_coll-1))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1)) 	  	  
      CALL  spline(R_COM,Mat_el_res,n_r_coll,
     & deriv_bgn,deriv_END,
     & Mat_el_der_res)	  
      ENDIF

!!!!!  GATHERING THE REST	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      CALL MPI_GATHER(Mat_el_res,n_r_coll,MPI_DOUBLE_PRECISION,
     & Mat_el_resf,
     &	  n_r_coll,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 	  
      CALL MPI_GATHER(Mat_el_der_res,n_r_coll,MPI_DOUBLE_PRECISION,
     & Mat_el_der_resf,
     &	  n_r_coll,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
!!!!!  GATHERING THE REST      
      IF(MYID.EQ.0) THEN
      DO i=1,total_size-nproc*chunk_mpi_size
      Mat_el(i_nr_ini:i_nr_fin,i+nproc*chunk_mpi_size)
     & = Mat_el_resf(i_nr_ini:i_nr_fin,i)
      Mat_el_der(i_nr_ini:i_nr_fin,i+nproc*chunk_mpi_size) =
     & Mat_el_der_resf(i_nr_ini:i_nr_fin,i)	  
      ENDDO	  
      ENDIF	  
      ENDIF	 
	  TIME_1_MAT = MPI_Wtime() 
      IF(MYID.EQ.0 .and. print_matrix_defined) THEN
      IF(.not.unformat_defined) THEN	  
      CALL PRINT_MIJ
      PRINT*,
     &  "MATRIX HAS BEEN SAVED INTO THE FILE ""MTRX.DAT"" "	  
      ELSE	  
      CALL PRINT_MIJ_USER
      PRINT*,"MATRIX HAS BEEN SAVED INTO THE FILE ""MTRX_UF.DAT"" "	  
      ENDIF		  
      ENDIF		 
	  TIME_2_MAT = MPI_Wtime()
      TIME_WRITING_MATRIX = TIME_2_MAT - TIME_1_MAT	  
	 
! Bikram Start Dec 2019:	  
	  IF(bikram_mij_shIFt) then
	  DO i=1,n_r_coll
	  Mat_el(i,:) = Mat_el(i,:) - Mat_el(n_r_coll,:)
      ENDDO 
	  ENDIF
! Bikram End.

      ENDIF
      ENDIF  
      IF(total_size.gt.total_size_old) THEN
      ALLOCATE(Mat_rest(n_r_coll,total_size-total_size_old))
      ALLOCATE(Mat_rest_der(n_r_coll,total_size-total_size_old))	  
      chunk_mpi_size = (total_size-total_size_old)/nproc
!      IF(MYID.EQ.0) PRINT*,"CHUNK SIZE",chunk_mpi_size	  
      k_st_mpi = myid*chunk_mpi_size + 1+total_size_old
      k_fn_mpi = (myid+1)*chunk_mpi_size+total_size_old
!      PRINT*,"INTEGRATION RANGE",myid,k_st_mpi,k_fn_mpi 	  
!      IF(myid.eq.nproc-1) k_fn_mpi = total_size
      task_size = (k_fn_mpi-k_st_mpi +1)*n_r_coll	  
      IF(myid.eq.0) PRINT*,"MATRIX_ADD_STARTED"
      IF(chunk_mpi_size.ge.0) THEN
!      PRINT*, "HELLO"!DELETE
      IF(chunk_mpi_size.gt.0)	THEN  
      ALLOCATE(Mat_el_temp(n_r_coll,k_fn_mpi-k_st_mpi +1)
     &, Mat_el_der_temp(n_r_coll,k_fn_mpi-k_st_mpi +1))
      ALLOCATE(K_SKIPPED_BY_ROUTINE(k_fn_mpi-k_st_mpi +1))
      DO k=k_st_mpi,k_fn_mpi
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.FALSE.	  
      ENDDO
      ENDIF	  
      IF(make_grid_file .and. (coll_type.gt.5 .or. coll_type.eq.0))
     & THEN
      IF(MYID.EQ.0) THEN
      PRINT*, "MATRIX WILL BE COMPUTED FROM GRID"	  
      SELECT CASE(coll_type)
      CASE(6)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1)
      READ(1) n_beta1,n_beta2,n_gamma1,n_r_vib1,
     & n_r_vib2, nr_cold
      ALLOCATE(V_VIB_2_2(n_gamma1,n_beta2,n_beta1,n_r_vib2,
     & n_r_vib1,nr_cold)) 
      CASE(7)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1)
      READ(1) n_alpha2,n_beta1,n_beta2,n_gamma1,nr_cold 
      ALLOCATE(
     & V_3_2(n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  
      CASE(8)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1)
      READ(1) n_alpha2,n_beta1,n_beta2,n_gamma1,nr_cold 
      ALLOCATE(
     & V_3_2(n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  
      CASE(9)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1)
      READ(1) n_alpha2,n_beta1,n_beta2,n_gamma1,n_gamma2,nr_cold 
      ALLOCATE(
     & V_3_3(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  	  
      CASE(0)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1)
      READ(1) n_alpha2,n_beta1,n_beta2,n_gamma1,n_gamma2,nr_cold 
      ALLOCATE(
     & V_3_3(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  
      END SELECT
      PRINT*,"VIJ_BUFFER_ALLOCATED"	
!!!!! READS AND PASSE V_ij BUFFER FOR FURTHER BROADCASTING,
!!!! ONLY FOR CASE OF ASYMETRIC + ASYMETRIC	  
      ENDIF
      CALL MPI_BCAST(n_alpha2, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_beta1, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_beta2, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_gamma1, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_gamma2, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_r_vib1, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_r_vib2, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	 
      CALL MPI_BCAST(nr_cold, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(nr_cold.lt.n_r_coll) THEN 
      IF(MYID.EQ.0)
     & PRINT*,"ERROR: WRONG FILE V_IJ GRID"
!!!! BROADCASTING ENDED NOW CHECKING	  
      STOP
      ENDIF
      SELECT CASE(coll_type)
      CASE(6)
      ALLOCATE(V_VIB_2_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_r_vib2,
     & n_r_vib1)) 
      CASE(7)
      ALLOCATE(V_3_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      CASE(8)
      ALLOCATE(V_3_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      CASE(9)
      ALLOCATE(
     & V_3_3_int_buffer(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      CASE(0)
!!!  MINI BUFFERS FOR EACH R-VALUE TO AVOID MEMORY LIMITATIONS	  
      ALLOCATE(
     & V_3_3_int_buffer(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2))
!      PRINT*,"BUFFER_ALLOCATED FOR=",MYID
      END SELECT	 

      DO i=1,n_r_coll
      IF(MYID.EQ.0) THEN	  
      SELECT CASE(coll_type)
      CASE(6)

      READ(1) 
     & V_VIB_2_2_int_buffer
      V_VIB_2_2(:,:,:,:,:,i)=
     & V_VIB_2_2_int_buffer	  

      CASE(7)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_2_int_buffer
      V_3_2(:,:,:,:,i)=
     & V_3_2_int_buffer	 	 
  
      CASE(8)
 
!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_2_int_buffer
      V_3_2(:,:,:,:,i)=
     & V_3_2_int_buffer	 

      CASE(9)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_3_int_buffer
      V_3_3(:,:,:,:,:,i)=
     & V_3_3_int_buffer		 
	  
      CASE(0)	  

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_3_int_buffer
      V_3_3(:,:,:,:,:,i)=
     & V_3_3_int_buffer	 
      END SELECT	  
      PRINT*,"BUFFER HAS BEEN READ FOR Ri = ",i
	  
      ENDIF
      SELECT CASE(coll_type)
      CASE(6)
      buffer_size_V = n_r_vib1*n_r_vib2*n_beta1*n_beta2*n_gamma1
      CASE(7)
      buffer_size_V = n_gamma1*n_beta2*n_beta1*n_alpha2	  
      CASE(8)
      buffer_size_V = n_gamma1*n_beta2*n_beta1*n_alpha2	  
      CASE(9)
      buffer_size_V = n_gamma2*n_gamma1*n_beta2*n_beta1*n_alpha2	  
      CASE(0)	  
      buffer_size_V = n_gamma2*n_gamma1*n_beta2*n_beta1*n_alpha2
      END SELECT	  
      IF(MYID.EQ.0) PRINT*,"BUFFER OF VIJ SIZE= ",
     & buffer_size_V!,size(V_3_3_int_buffer)	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
!!!!! BROADCASTING MINI _BUFFER	  
      IF(MYID.EQ.0) PRINT*,"BROADCASTING OF VIJ_BUFFER STARTED"
!      IF(MYID.EQ.0) WRITE(3,*)	V_3_3_int_buffer
      SELECT CASE(coll_type)
      CASE(6)
      CALL MPI_BCAST(V_VIB_2_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CASE(7)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(8)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(9)
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(0)  
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      END SELECT	 
      IF(MYID.EQ.0) PRINT*,"THE BUFFER HAS BEEN BROADCASTED"
      IF(MYID.eq.0 .and. chunk_mpi_size.gt.0 .and. i.eq.1 )
     & PRINT*, "ADD_MATRIX IS BEING COMPUTED ON V_GRID"	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i,cyc) !!! CALCULATING VALUE OVER INTEGRATION FOR EACH R
			
		Write(*,*) '#8'
!      PRINT*,st_1,st_2,intgeral	  
      Mat_el_temp(i,k-k_st_mpi+1) = intgeral*conv_unit_e
!      IF(ABS(Mat_el_temp(1,k-k_st_mpi+1)).LT.MIJ_ZERO_CUT) THEN    !!! EXCLUDE SOME MATRIX ELEMENTS IF THEY ARE BELOW CERTAIN VALUE
      IF(max(abs(Mat_el_temp(mtrx_cutoff_r1,k-k_st_mpi+1)),
     & abs(Mat_el_temp(mtrx_cutoff_r2,k-k_st_mpi+1))).LT.
     & MIJ_ZERO_CUT) THEN    !!! EXCLUDE SOME MATRIX ELEMENTS IF THEY ARE BELOW CERTAIN VALUE
      Mat_el_temp(:,k-k_st_mpi+1) = 0d0
      Mat_el_der_temp(:,k-k_st_mpi+1) = 0d0	  
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.TRUE.	  
      ENDIF	  
!      PRINT*,st_1,st_2,	 intgeral 
      ENDDO
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      ENDDO	  
	  
      ELSE	  
!!!   REGULAR COMPUTING
      IF(MYID.eq.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "ADD_MATRIX IS BEING COMPUTED"	  
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE 	  
      DO i=1,n_r_coll
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
	  
c      R_dist =  R_COM(i)*conv_unit_r
c      PRINT*,st_1,st_2,R_dist
c      STOP	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i,cyc)
			
		Write(*,*) '#9'
c      PRINT*,st_1,st_2,intgeral	  
      Mat_el_temp(i,k-k_st_mpi+1) = intgeral*conv_unit_e
!      IF(ABS(Mat_el_temp(1,k-k_st_mpi+1)).LT.MIJ_ZERO_CUT) THEN
      IF(max(abs(Mat_el_temp(mtrx_cutoff_r1,k-k_st_mpi+1)),
     & abs(Mat_el_temp(mtrx_cutoff_r2,k-k_st_mpi+1))).LT.
     & MIJ_ZERO_CUT) THEN
      Mat_el_temp(:,k-k_st_mpi+1) = 0d0
      Mat_el_der_temp(:,k-k_st_mpi+1) = 0d0	  
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.TRUE.	  
      ENDIF	  
c      PRINT*,st_1,st_2,	 intgeral 
      ENDDO
  
      ENDDO
      ENDIF	  
!!!   END REGULAR COMPUTING
!! SPLINING	 
      IF(MYID.eq.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "ADD_MATRIX IS BEING SPLINED" 
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE
      deriv_bgn = (Mat_el_temp(2,k-k_st_mpi+1)
     & - Mat_el_temp(1,k-k_st_mpi+1))/
     & (R_COM(2)-R_COM(1)) 
      deriv_END = (Mat_el_temp(n_r_coll,k-k_st_mpi+1)
     & - Mat_el_temp(n_r_coll-1,k-k_st_mpi+1))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1))	  
      CALL  spline(R_COM,Mat_el_temp(:,k-k_st_mpi+1),
     & n_r_coll,deriv_bgn,deriv_END,
     & Mat_el_der_temp(:,k-k_st_mpi+1))		  
      ENDDO	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF(chunk_mpi_size.gt.0) THEN	  
      CALL MPI_GATHER(Mat_el_temp,task_size,MPI_DOUBLE_PRECISION,
     & Mat_rest,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 	  
      CALL MPI_GATHER(Mat_el_der_temp,task_size,MPI_DOUBLE_PRECISION,
     & Mat_rest_der,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
      ENDIF	 
      ENDIF		 
!!!   COMPUTING CHUNKS RESIDUE
      IF(expansion_defined) THEN
      CALL READ_EXPANSION_TERMS	  
      ENDIF
      IF(total_size-total_size_old
     & .gt.nproc*chunk_mpi_size) THEN
      IF(myid.eq.0) PRINT*,"REST OF THE CHUNK COMPUTING STARTED"
      IF(myid.eq.0) PRINT*,"THE SIZE OF CHUNK PER EACH PROCESSOR= ",
     & chunk_mpi_size	  
      ALLOCATE(Mat_el_res(n_r_coll),Mat_el_der_res(n_r_coll))
      ALLOCATE(Mat_el_resf(n_r_coll,nproc),
     & Mat_el_der_resf(n_r_coll,nproc))
      Mat_el_res = 0d0
      Mat_el_der_res = 0d0	  
      DO i=1,n_r_coll
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
      IF(make_grid_file .and. (coll_type.gt.5 .or. coll_type.eq.0))
     & THEN
      IF(MYID.EQ.0) THEN
      SELECT CASE(coll_type)
      CASE(6)
      V_VIB_2_2_int_buffer(:,:,:,:,:) = 
     & V_VIB_2_2(:,:,:,:,:,i)
      CASE(7)
      V_3_2_int_buffer(:,:,:,:) = V_3_2(:,:,:,:,i)	  
      CASE(8)
      V_3_2_int_buffer(:,:,:,:) = V_3_2(:,:,:,:,i)	  
      CASE(0)	  
      V_3_3_int_buffer(:,:,:,:,:) = V_3_3(:,:,:,:,:,i)
      END SELECT	  
      ENDIF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      SELECT CASE(coll_type)
      CASE(6)
      CALL MPI_BCAST(V_VIB_2_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CASE(7)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(8)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(9)
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(0)	  
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      PRINT*, "CHUNK BUFFER HAS BEEN BROADCASTED FOR IR=",i	 
      END SELECT	 
      ENDIF	  	  
      IF(K_SKIPPED_RES) CYCLE	  
c      R_dist =  R_COM(i)*conv_unit_r
      IF(myid.lt.total_size-total_size_old
     & -nproc*chunk_mpi_size) THEN
      k=nproc*chunk_mpi_size+1+myid+total_size_old 	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i,cyc)
			
		Write(*,*) '#10'
      Mat_el_res(i) = intgeral*conv_unit_e
!      IF(ABS(Mat_el_res(1)).LT.MIJ_ZERO_CUT) THEN
      IF(max(abs(Mat_el_res(mtrx_cutoff_r1)),
     & abs(Mat_el_res(mtrx_cutoff_r2))).LT.MIJ_ZERO_CUT) THEN
      Mat_el_res(:) = 0d0
      Mat_el_der_res(:) = 0d0	  
      K_SKIPPED_RES=.TRUE.	  
      ENDIF
      ENDIF	  
      ENDDO
c      PRINT*,st_1,st_2,intgeral	  
      IF(.NOT.K_SKIPPED_RES) THEN
      deriv_bgn = (Mat_el_res(2)
     & - Mat_el_res(1))/
     & (R_COM(2)-R_COM(1)) 
      deriv_END = (Mat_el_res(n_r_coll)
     & - Mat_el_res(n_r_coll-1))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1))	 	  
      CALL  spline(R_COM,Mat_el_res,n_r_coll,
     & deriv_bgn,deriv_END,
     & Mat_el_der_res)	  
      ENDIF

      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
!      STOP		  
!!!!!!!!!!!!! SEND RECIEVE
      CALL MPI_GATHER(Mat_el_res,n_r_coll,MPI_DOUBLE_PRECISION,
     & Mat_el_resf,
     &	  n_r_coll,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 	  
      CALL MPI_GATHER(Mat_el_der_res,n_r_coll,MPI_DOUBLE_PRECISION,
     & Mat_el_der_resf,
     &	  n_r_coll,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
      IF(MYID.EQ.0) THEN
      DO i=1,total_size-total_size_old-nproc*chunk_mpi_size
      Mat_rest(:,i+nproc*chunk_mpi_size) = Mat_el_resf(:,i)
      Mat_rest_der(:,i+nproc*chunk_mpi_size) = Mat_el_der_resf(:,i)	  
      ENDDO	  
      ENDIF	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	 
	 
!!!!!!!!! SEND RECIEVE	 
	  
      ENDIF
      IF(MYID.EQ.0) THEN	  
      DO i=total_size_old+1,total_size
      Mat_el(:,i) = Mat_rest(:,i-total_size_old)
      Mat_el_der(:,i) = Mat_rest_der(:,i-total_size_old)	  
      ENDDO
      ENDIF	  
	  
	  
      IF(MYID.EQ.0 .and. print_matrix_defined) THEN
      IF(.not.unformat_defined) THEN	  
      CALL PRINT_MIJ
      PRINT*,
     &  "MATRIX HAS BEEN SAVED INTO THE FILE ""MTRX.DAT"" "	  
      ELSE	  
      CALL PRINT_MIJ_USER
      PRINT*,"MATRIX HAS BEEN SAVED INTO THE FILE ""MTRX_UF.DAT"" "	  
      ENDIF		
      ENDIF	  
!      IF(MYID.EQ.0) CALL PRINT_MREST
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)	  
	  
! Bikram Start Dec 2019:	  	  
	  IF(bikram_mij_shIFt) then
	  DO i=1,n_r_coll
	  Mat_el(i,:) = Mat_el(i,:) - Mat_el(n_r_coll,:)
      ENDDO 
	  ENDIF
! Bikram End.
	  
      ENDIF
	  IF(run_prog_defined .and. 
     & coupled_states_defined .and. myid.eq.0) THEN
      ALLOCATE(mat_buffer(n_r_coll,total_size))
      ALLOCATE(mat_buffer_der(n_r_coll,total_size))
	  mat_buffer = Mat_el
	  mat_buffer_der = Mat_el_der
      ENDIF	
!      TIME_MAT_START = MPI_Wtime()

! finding #r for matrix truncation
	  IF(.not.cut_r) then
	  allocate(deviation_r(n_r_coll))
	  do ii = 1, n_r_coll
	  deviation_r(ii) = abs(R_COM(ii) - bikram_cutoff_r1)
	  END do
	  mtrx_cutoff_r1 = minloc(deviation_r,1)
	  deallocate(deviation_r)
	  allocate(deviation_r(n_r_coll))
	  do ii = 1, n_r_coll
	  deviation_r(ii) = abs(R_COM(ii) - bikram_cutoff_r2)
	  END do
	  mtrx_cutoff_r2 = minloc(deviation_r,1)
	  deallocate(deviation_r)
	  cut_r = .true.
	  IF(myid.eq.0) write(*,'(2(a,i0,a,f0.3))') 
     & "Truncation #1 at #R = ", mtrx_cutoff_r1, ", R = ", 
     & R_COM(mtrx_cutoff_r1), ", and #2 at #R = ",mtrx_cutoff_r2, 
     & ", R = ",R_COM(mtrx_cutoff_r2)
	  MIJ_ZERO_CUT = MIJ_ZERO_CUT/eVtown/autoeV
	  END IF

      IF(MYID.EQ.0) THEN
      DO k=1,total_size 
!      IF(ABS(Mat_el(1,k)).lt.MIJ_ZERO_CUT) THEN
      IF(max(abs(Mat_el(mtrx_cutoff_r1,k)),
     & abs(Mat_el(mtrx_cutoff_r2,k))).lt.MIJ_ZERO_CUT) THEN
      Ks_have_to_be_skipped(k) = .TRUE.
      mean_size_1 = mean_size_1 - 1	  
!   	  IF(ind_mat(1,k).eq.ind_mat(2,k)) THEN
 !     PRINT*, "ERROR GLOBAL IN MIJ",k
!      run_prog_defined = .FALSE.	  
!      CALL MPI_BCAST(run_prog_defined,1, MPI_LOGICAL,0,
!     &  MPI_COMM_WORLD,ierr_mpi)	  
!      ENDIF	  
	  
      ELSE
      Ks_have_to_be_skipped(k) = .FALSE.     
      ENDIF	 
      ENDDO	 
      ENDIF	
      IF(MYID.EQ.-1) THEN
      OPEN(345,FILE="MAT_CHECK.out")	  
      DO i=1,n_r_coll	  
      WRITE(345,'(e19.12,1x,e19.12)')R_COM(i),Mat_el(i,1)
      ENDDO
      CLOSE(345)	  
      ENDIF
      IF(MYID.EQ.0)WRITE(*,'(a35,1x,i9)')
     & "THE SIZE OF NON_ZERO PART OF Mij = ",mean_size_1	  
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
      IF(.not. run_prog_defined) THEN
      TIME_MAT_FINISH= MPI_Wtime()
      IF(MYID.EQ.0) WRITE(*,'(a57,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij READING/COMPUTING/SAVING ON DISK"
     & ,",s",(TIME_MAT_FINISH-TIME_MAT_START)/nproc
	  IF(print_matrix_defined.and. myid.eq.0) 
     & WRITE(*,'(a47,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij SAVING ON DISK"
     & ,",s",TIME_WRITING_MATRIX	 
      IF(MYID.eq.0)PRINT*, "ALL WORK ON MATRIX IS DONE"	  
      CALL MPI_FINALIZE (ierr_mpi)	  
      STOP
	  
      ELSE
      IF(MYID.EQ.0) THEN	  
      ALLOCATE(Mat_el_non_zero(n_r_coll,mean_size_1))
      ALLOCATE(Mat_el_non_zero_der(n_r_coll,mean_size_1))	  
      ALLOCATE(ind_mat_non_zero(2,mean_size_1))
      k_non_zero = 0	  
      DO k=1,total_size
      IF(Ks_have_to_be_skipped(k)) CYCLE
      k_non_zero = k_non_zero + 1
      ind_mat_non_zero(:,k_non_zero) = ind_mat(:,k)
      Mat_el_non_zero(:,k_non_zero) = Mat_el(:,k)
      Mat_el_non_zero_der(:,k_non_zero) = Mat_el_der(:,k)	  
      ENDDO 
      IF(k_non_zero.ne.mean_size_1) CRITICAL_ERROR = .TRUE.
      total_size = mean_size_1	  
      DEALLOCATE(Mat_el)
      DEALLOCATE(Mat_el_der)	  
      DEALLOCATE(ind_mat)
      IF(.NOT.mpi_task_defined) THEN	  
      ALLOCATE(Mat_el(n_r_coll,total_size))	  
      ALLOCATE(Mat_el_der(n_r_coll,total_size))
      ENDIF	  
      ALLOCATE(ind_mat(2,total_size))
      IF(.NOT.mpi_task_defined) THEN  
      Mat_el = Mat_el_non_zero
      Mat_el_der = Mat_el_non_zero_der
      ENDIF	  
      ind_mat = ind_mat_non_zero
      	  
      IF(.NOT.mpi_task_defined) THEN
      DEALLOCATE(Mat_el_non_zero,Mat_el_non_zero_der,ind_mat_non_zero)
      ENDIF 	  
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)	  
	  
      CALL MPI_BCAST(total_size, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CALL MPI_BCAST(CRITICAL_ERROR, 1, MPI_LOGICAL,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(CRITICAL_ERROR) STOP "CRITICAL ERROR"
      IF(MYID.NE.0) THEN
      DEALLOCATE(ind_mat)
      ALLOCATE(ind_mat(2,total_size))
      IF(.NOT.mpi_task_defined) THEN	  
      DEALLOCATE(Mat_el)
      DEALLOCATE(Mat_el_der)
      ALLOCATE(Mat_el(n_r_coll,total_size))	  
      ALLOCATE(Mat_el_der(n_r_coll,total_size))
      ENDIF
      ENDIF	  
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)	  
      TIME_MAT_FINISH = MPI_Wtime()
      TIME_3_MAT = TIME_MAT_FINISH - TIME_MAT_START 	  
      IF(mpi_task_defined) THEN
!      TIME_MAT_START = MPI_Wtime()	  
      IF(myid.eq.0) PRINT*,"MPI TASK PER TRAJECTORY WILL BE USED"
      IF(myid.eq.0)	
     & WRITE(*,'(a53,1x,i4)')
     & "MPI TASKS WHICH ARE ASSOCIATED WITH ONE TRAJECTORY = ",
     & mpi_task_per_traject	 
      traject_roots = nproc/mpi_task_per_traject
      IF(traject_roots*mpi_task_per_traject.ne.nproc) THEN
      STOP "mpi_task_number must be a delimeter of nproc"	  
      ENDIF	  
      ALLOCATE(mpi_traject_roots(traject_roots))
      ALLOCATE(portion_of_MIJ_per_task(2,nproc))
      ALLOCATE(portion_of_state_per_task(2,nproc))
      ALLOCATE(portion_of_work_per_task(2,nproc))	  
      ALLOCATE(mpi_root_belongs(nproc))	  
!!!!!!!!!!!!!   REMAKE	  
      size_mij_chunk_mpi = total_size/mpi_task_per_traject
      residue_mij_mpi = total_size-
     & size_mij_chunk_mpi*mpi_task_per_traject
      size_state_chunk_mpi = states_size/mpi_task_per_traject
      residue_state_mpi = states_size-
     & size_state_chunk_mpi*mpi_task_per_traject
      size_work_chunk_mpi = (2*states_size+8)/mpi_task_per_traject !!!! DO IN FUTURE
      residue_work_mpi = 2*states_size+8-!!!! DO IN FUTURE
     & size_work_chunk_mpi*mpi_task_per_traject	!!!! DO IN FUTURE
      DO k=1,traject_roots
      mpi_traject_roots(k) = (k-1)*mpi_task_per_traject
      DO k_p=1,mpi_task_per_traject	  
      mpi_root_belongs(k_p+mpi_traject_roots(k))=mpi_traject_roots(k)
      ENDDO	  
      ENDDO
      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_mij_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_MIJ_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_mij_chunk_mpi+1)
      portion_of_MIJ_per_task(2,k_mpi_proc) =
     & (k_p)*(size_mij_chunk_mpi+1)
      ENDDO
      total_size_check = total_size_check + size_mij_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_MIJ_per_task(1,k_mpi_proc)
     & = residue_mij_mpi*(size_mij_chunk_mpi+1)+1  
     & + (k_p-1-residue_mij_mpi)*size_mij_chunk_mpi
      portion_of_MIJ_per_task(2,k_mpi_proc)
     & = residue_mij_mpi*(size_mij_chunk_mpi+1)+	  
     & (k_p-residue_mij_mpi)*size_mij_chunk_mpi
      ENDDO	 
      total_size_check = total_size_check + size_mij_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_MIJ_per_task(2,k_mpi_proc).ne.total_size) 
     & STOP "WRONG MATIX ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO	
      IF(total_size_check.ne.total_size) THEN
      IF(myid.eq.0) WRITE(*,*)total_size_check,total_size	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF	  

      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_state_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_state_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_state_chunk_mpi+1)
      portion_of_state_per_task(2,k_mpi_proc) =
     & (k_p)*(size_state_chunk_mpi+1)
      ENDDO
      state_size_check = state_size_check + size_state_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_state_per_task(1,k_mpi_proc)
     & = residue_state_mpi*(size_state_chunk_mpi+1)+1  
     & + (k_p-1-residue_state_mpi)*size_state_chunk_mpi
      portion_of_state_per_task(2,k_mpi_proc)
     & = residue_state_mpi*(size_state_chunk_mpi+1)+	  
     & (k_p-residue_state_mpi)*size_state_chunk_mpi
      ENDDO	 
      state_size_check = state_size_check + size_state_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_state_per_task(2,k_mpi_proc).ne.states_size) 
     & STOP "WRONG MATRIX ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO	  
	  
!      PRINT*,"STATES",portion_of_state_per_task	  
      IF(state_size_check.ne.states_size) THEN
      IF(myid.eq.0) WRITE(*,*)state_size_check,states_size	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF	

      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_work_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_work_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_work_chunk_mpi+1)
      portion_of_work_per_task(2,k_mpi_proc) =
     & (k_p)*(size_work_chunk_mpi+1)
      ENDDO
      work_size_check = work_size_check + size_work_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_work_per_task(1,k_mpi_proc)
     & = residue_work_mpi*(size_work_chunk_mpi+1)+1  
     & + (k_p-1-residue_work_mpi)*size_work_chunk_mpi
      portion_of_work_per_task(2,k_mpi_proc)
     & = residue_work_mpi*(size_work_chunk_mpi+1)+	  
     & (k_p-residue_work_mpi)*size_work_chunk_mpi
      ENDDO	 
      work_size_check = work_size_check + size_work_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_work_per_task(2,k_mpi_proc).ne.states_size*2+8) 
     & STOP "WRONG MATRIX ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO	  
	  
      IF(work_size_check.ne.states_size*2+8) THEN
      IF(myid.eq.0) WRITE(*,*)work_size_check,states_size*2+8	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF
	  
      total_size_mpi = portion_of_MIJ_per_task(2,myid+1) - 
     & portion_of_MIJ_per_task(1,myid+1) + 1	  
      ALLOCATE(Mat_el(n_r_coll,total_size_mpi))
      ALLOCATE(Mat_el_der(n_r_coll,total_size_mpi))
      CALL MPI_Comm_group(MPI_COMM_WORLD, wrld_group,ierr_mpi)
      ALLOCATE(process_rank_distr(mpi_task_per_traject,traject_roots))
      ALLOCATE(comms_distr(mpi_task_per_traject),
     & groups_distr(mpi_task_per_traject))   	  
      DO i=1,traject_roots	  
      DO k=1,mpi_task_per_traject
      process_rank_distr(k,i) = k - 1 + (i-1)*mpi_task_per_traject     	  
      ENDDO
      ENDDO	

      DO i=1,mpi_task_per_traject      	  
      CALL MPI_Group_incl(wrld_group, traject_roots,
     & process_rank_distr(i,:), groups_distr(i),ierr_mpi)
      CALL MPI_Comm_create(MPI_COMM_WORLD,groups_distr(i),
     & comms_distr(i),ierr_mpi)
      ENDDO	  
!      PRINT*,myid,total_size_mpi	  
      tag1 = 1
      tag2 = 2
      IF(MYID.eq.0) PRINT*,"Mij ROOT PROC DISTRIBUTION STARTED"	  
      IF(MYID.eq.0) THEN
      DO k=2,mpi_task_per_traject
      total_size_mpi= portion_of_MIJ_per_task(2,k) - 
     & portion_of_MIJ_per_task(1,k) + 1	
      task_portion_size = total_size_mpi*n_r_coll	 
      ALLOCATE(buffer_mpi_portion(n_r_coll,total_size_mpi))
      buffer_mpi_portion = Mat_el_non_zero(:,
     & portion_of_MIJ_per_task(1,k):portion_of_MIJ_per_task(2,k)) 	  
      CALL MPI_SEND(buffer_mpi_portion,
     & task_portion_size, MPI_REAL8, k-1, 
     &  tag1, MPI_COMM_WORLD, ierr_mpi)
      buffer_mpi_portion = 	 Mat_el_non_zero_der (:,
     & portion_of_MIJ_per_task(1,k):portion_of_MIJ_per_task(2,k))
      CALL MPI_SEND(buffer_mpi_portion,
     & task_portion_size, MPI_REAL8, k-1, 
     &  tag2, MPI_COMM_WORLD, ierr_mpi)
      DEALLOCATE(buffer_mpi_portion) 
      ENDDO
      total_size_mpi= portion_of_MIJ_per_task(2,1) - 
     & portion_of_MIJ_per_task(1,1) + 1		  
      Mat_el=Mat_el_non_zero(:,
     & portion_of_MIJ_per_task(1,1):portion_of_MIJ_per_task(2,1))
      Mat_el_der=Mat_el_non_zero_der(:,
     & portion_of_MIJ_per_task(1,1):portion_of_MIJ_per_task(2,1))
      DEALLOCATE(Mat_el_non_zero,Mat_el_non_zero_der)	 
      ELSE
      IF(myid.le.mpi_task_per_traject-1) THEN	  
      task_portion_size = total_size_mpi*n_r_coll	  
      CALL MPI_RECV(Mat_el, task_portion_size, MPI_REAL8, 
     & 0, tag1, MPI_COMM_WORLD, status, ierr_mpi)	  
      CALL MPI_RECV(Mat_el_der, task_portion_size, MPI_REAL8, 
     & 0, tag2, MPI_COMM_WORLD, status, ierr_mpi)
      ENDIF	 
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
      IF(myid.eq.0)	 PRINT*, "Mij ROOT PROC DISTRIBUTION DONE"  
      IF(myid.eq.0) PRINT*,"Mij ALL PROC DISTRIBUTION STARTED" 	  
      DO i=1,mpi_task_per_traject
      IF(i-1 .eq. myid - int(myid/mpi_task_per_traject)
     & *mpi_task_per_traject) THEN
      CALL MPI_Comm_rank(comms_distr(i),id_proc_in_group ,ierr_mpi)
      IF(id_proc_in_group.ne.myid/mpi_task_per_traject) PRINT*,
     & "ERROR IN COMUNICATIONS ASSIGNEMNET_1_distr",id_proc_in_group,
     & i-1
      IF(myid.lt.mpi_task_per_traject .and. id_proc_in_group.ne.0 )
     & PRINT*,"ERROR IN COMUNICATIONS ASSIGNEMNET_2_distr"
      task_portion_size = total_size_mpi*n_r_coll	  
      CALL MPI_BCAST(Mat_el,task_portion_size, MPI_REAL8,0,
     &  comms_distr(i),ierr_mpi)
      CALL MPI_BCAST(Mat_el_der,task_portion_size, MPI_REAL8,0,
     &  comms_distr(i),ierr_mpi)
      ENDIF	 
      ENDDO
      IF(myid.eq.0) PRINT*,"Mij ALL PROC DISTRIBUTION DONE"	  
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)	  
      IF(MYID.EQ.0) PRINT *,"MPI_COMMUNICATORS_CREATION STARTED"
      ALLOCATE(comms(traject_roots),groups(traject_roots))     	  
      ALLOCATE(process_rank(traject_roots,mpi_task_per_traject))		  
      DO i=1,traject_roots	  
      DO k=1,mpi_task_per_traject
      process_rank(i,k) = k - 1 + (i-1)*mpi_task_per_traject     	  
      ENDDO
      ENDDO	

      DO i=1,traject_roots      	  
      CALL MPI_Group_incl(wrld_group, mpi_task_per_traject,
     & process_rank(i,:), groups(i),ierr_mpi)
      CALL MPI_Comm_create(MPI_COMM_WORLD,groups(i),comms(i),ierr_mpi)
      ENDDO

      DO i=1,traject_roots
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN
      CALL MPI_Comm_rank(comms(i),id_proc_in_group ,ierr_mpi)	  
      IF(mpi_traject_roots(i)+id_proc_in_group.ne.myid) STOP
     & "WRONG GROUP MPI ASSIGNEMENT"	  
	  
      ENDIF		  
      ENDDO	 
      IF(MYID.EQ.0) PRINT *,"MPI_COMMUNICATORS_CREATION DONE"
      TIME_MAT_FINISH = MPI_Wtime()
      TIME_2_MAT = TIME_MAT_FINISH - TIME_MAT_START	  
      ENDIF	  
! mpi_task_per_traject
! mpi_task_defined
      IF(.NOT.mpi_task_defined) THEN 	  
      CALL MPI_BCAST(Mat_el, n_r_coll*total_size, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(Mat_el_der, n_r_coll*total_size, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(Ks_have_to_be_skipped,total_size, MPI_LOGICAL,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      ENDIF
      CALL MPI_BCAST(R_COM, n_r_coll, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(ind_mat, 2*total_size, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
!	  call resize
      ENDIF
      IF(states_to_exclude_defined) THEN
      ALLOCATE(stts_to_excl(states_size))	  
      IF(myid.eq.0) THEN
      CALL EXCLUDE_STATES	  
      ENDIF	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      	  
      CALL MPI_BCAST(stts_to_excl,states_size,MPI_LOGICAL,0,
     & MPI_COMM_WORLD,ierr_mpi)
      ENDIF	 	  
!      ALLOCATE(dq_dt_mpi(size_mij_chunk_mpi,states_size*2+8))
      TIME_MAT_FINISH= MPI_Wtime()
      IF(myid.eq.0) THEN	  
      WRITE(*,'(a57,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij READING/COMPUTING/SAVING ON DISK"
     & ,",s",TIME_MAT_FINISH-TIME_MAT_START
	  IF(print_matrix_defined.and. myid.eq.0) 
     & WRITE(*,'(a47,1x,a2,f10.1)')  
     & "TIME SPENT ON MATRIX Mij SAVING ON DISK"
     & ,",s",TIME_WRITING_MATRIX
      ENDIF	 
      IF(.not.run_prog_defined) STOP	  
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
!      STOP	  
      END SUBROUTINE MATRIX_MAKE
      SUBROUTINE PRINT_MIJ
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      IMPLICIT NONE
      INTEGER st,i,k	  
      PRINT*,"SYSTEM_SETUP_DONE"
      OPEN(1,FILE=MATRIX_NAME_MIJ,ACTION="WRITE")
      WRITE(1,'(a15,1x,i1)') "COLLISION_TYPE=",coll_type
      WRITE(1,'(a17,x,i4)') "NUMBER_OF_CHANLS=",number_of_channels
      WRITE(1,'(a17,1x,i6)') "NUMBER_OF_STATES=",states_size
      WRITE(1,'(a12,1x,i9)') "MATRIX_SIZE=",total_size
      WRITE(1,'(a12,1x,i4)') "R_GRID_SIZE=",n_r_coll	  
      SELECT CASE(coll_type)
      CASE(1)
      IF(.not.fine_structure_defined) THEN	  
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j_ch(i)	  
      ENDDO
      ELSE
      IF(SPIN_FINE.ne.2) THEN       	  
      WRITE(1,'(a8,1x,a8,1x,a3,1x,a3,1x,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J","F","N"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i3,1x,i3,1x,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j_ch(i),2+f_ch(i),j_ch(i)+f_ch(i)	  
      ENDDO
      ELSE
      WRITE(1,'(a8,1x,a8,1x,a4,1x,a5,1x,1x,a4,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J","F","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,f4.1,1x,f5.1,1x,1x,f4.1,1x,i3,1x,i3)')
     & st,i,j12_h(st),m12_h(st),j_h_ch(i),f_ch(i),par_lorb_ch(i)	  !!!! TO MODIFY
      ENDDO	  
      ENDIF	  
	  
      ENDIF	  
      CASE(2)
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","V","J"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),v_ch(i),j_ch(i)	  
      ENDDO
      CASE(3)
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J","K", "EPS"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j_ch(i),k_ch(i),eps_ch(i)	  
      ENDDO
      CASE(4)
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J","KA", "KC"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j_ch(i),ka_ch(i),kc_ch(i)	  
      ENDDO      	  
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J1","J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J1","J2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i),parity_state(st)	  
      ENDDO	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined) THEN	  
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","V1","V2","J1","J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),v1_ch(i),v2_ch(i),j1_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","V1","V2","J1","J2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),v1_ch(i),v2_ch(i),j1_ch(i),j2_ch(i),
     & parity_state(st)	 
      ENDDO	  
      ENDIF	  
      CASE(7)
      WRITE(1,
     & '(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J1","K","EPS",
     & "J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j1_ch(i),k1_ch(i),eps1_ch(i)
     & ,j2_ch(i) 
      ENDDO	  
      CASE(8)
      WRITE(1,
     & '(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J1","K1A","K1C",
     & "J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i) 
      ENDDO		  
      CASE(9)
      WRITE(1,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","K1A","K1C",
     & "J2","K2","EPS"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),k2_ch(i),eps2_ch(i) 
      ENDDO		  
      CASE(0)
      IF(.NOT.identical_particles_defined) THEN  
      WRITE(1,
     & '(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J1","KA1","KC1",
     & "J2","KA2","KC2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i)	 
      ENDDO
      ELSE
      WRITE(1,
     & '(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a2)
     & ')
     & "#STATE","#CHANNEL","J12","M12","J1","KA1","KC1",
     & "J2","KA2","KC2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i2)
     & ')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i),parity_state(st)
      ENDDO	 
      ENDIF	 
      END SELECT
      WRITE(1,"(a16,1x,a8,1x,a8)") "MIJ_INDEX", "ST_1", "ST_2"
      DO i=1,total_size
      WRITE(1,"(i16,1x,i8,1x,i8)") i,ind_mat(1,i),ind_mat(2,i)	  
      ENDDO		  
      WRITE(1,'(a5,1x,a19)') "#iR","R_COM(iR)"	  
      DO i=1,n_r_coll
      WRITE(1,'(i5,1x,e19.12)')i,R_COM(i)	  
      ENDDO
      WRITE(1,'(a16,1x,a8,1x,a19)') "#k_matel","#iR",
     & "Mij(iR,k_matel)"!,
!     & "Mij_der(iR,k_matel)"	  
      DO k=1,total_size
      DO i=1,n_r_coll	  
      WRITE(1,'(i16,1x,i8,1x,e19.12)')k,i,Mat_el(i,k)!,
!     & Mat_el_der(i,k)	  
      ENDDO
      ENDDO
      CLOSE(1)
!      OPEN(333,file="LAST_MIJ.DAT")	   !!!!! DELETE AFTER ALL
!	  DO k=1,total_size !!!!! DELETE AFTER ALL
!	  IF(ind_mat(1,k).eq.ind_mat(2,k)) THEN
!	  WRITE(333,*) k,Mat_el(n_r_coll,k)
!      ENDIF
!      ENDDO	
!      CLOSE(333)	  
      IF(test_expansion_defined) THEN 
      CALL PRINT_ELASTIC_MIJ
      PRINT*, "ELASTIC_TERMS_FILE_CREATED"  
      ENDIF	  
      END SUBROUTINE PRINT_MIJ

      SUBROUTINE READ_MIJ !!! TO DO
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      USE OLD_MIJ
      IMPLICIT NONE
      REAL*8 Mij_buffer	  
      INTEGER i,k,st,tot_siz_read,st_siz_read,j_count,j_summ,p_count
      LOGICAL file_exst, bk_exst
	  LOGICAL EVEN_NUM		! Dulat 9/25/2023
	  CHARACTER(LEN = 500) bk_matrix_path1, bk_matrix_path2
      INTEGER KRONEKER,p_lim_max_ini	  
      EXTERNAL KRONEKER	  
	  
! Bikram Start Nov 2021: Changing matrix file path 
!so that user don't need to copy the large matrix file in each directory.
	  IF(bikram_mtrx_path) then
	  inquire(file = "Matrix_Path.DAT", exist = bk_exst)
	  IF(.not.bk_exst) then
	  write(*,'(a)')"File containing matrix path is not provided."
	  stop
	  ELSE
	  open(111,file="Matrix_Path.DAT",status="old",action = "read")
	  read(111,'(a)') bk_matrix_path1
	  close(111)
	  write(*,'(a,a,a)') 'Matrix file path is defined. ',
     & 'Matrix reading from ', trim(bk_matrix_path1)
	  write(bk_matrix_path2,'(a,a,a)')trim(bk_matrix_path1),'/',
     & trim(MATRIX_NAME_MIJ)
	  END IF
	  END IF
	  
	  IF(.not.bikram_mtrx_path) then	  
	  INQUIRE( FILE=MATRIX_NAME_MIJ, EXIST=file_exst )
	  ELSE
	  INQUIRE( FILE=trim(bk_matrix_path2), EXIST=file_exst )
	  END IF
! Bikram End.
	  
      IF(.not.file_exst) THEN 
      PRINT*, "ERROR: MTRX.DAT NOT FOUND"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF
      IF(.not.identical_particles_defined) THEN	  
      p_lim_max = 1
      p_lim_min = 1	
      ELSE
      SELECT CASE(coll_type)		  
      CASE(5)	  

      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
      CASE(6)
      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))*
     & KRONEKER(v1_ch(chann_ini),v2_ch(chann_ini)) 	  
      CASE(0)
      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
     & *KRONEKER(ka1_ch(chann_ini),ka2_ch(chann_ini))
     & *KRONEKER(kc1_ch(chann_ini),kc2_ch(chann_ini))	  
      END SELECT		  
      ENDIF	  
	  
! Bikram Start Nov 2021: 
	  IF(.not.bikram_mtrx_path) then	  
	  OPEN(1,FILE=MATRIX_NAME_MIJ,STATUS="OLD",ACTION="READ")
	  ELSE
	  OPEN(1,FILE=trim(bk_matrix_path2), STATUS="OLD",ACTION="READ")
	  END IF
! Bikram End.

      PRINT*,"MTRX.DAT READING STARTED"	  
      READ(1,'(a15,1x,i1)') buffer_word_2,coll_type_old
      IF(coll_type_old.ne.coll_type) THEN
      PRINT*, "ERROR: SYSTEM IS DIFFERENT"
      CRITICAL_ERROR = .TRUE.	  
      ENDIF	 
      READ(1,'(a17,x,i4)') buffer_word_3,number_of_channels_old	 
      READ(1,'(a17,1x,i6)') buffer_word_3,states_size_old
      READ(1,'(a12,1x,i9)') buffer_word_1,total_size_old	  
      READ(1,'(a12,1x,i4)') buffer_word_1,n_r_coll_old
      IF(n_r_coll_old.ne.n_r_coll) THEN
      PRINT*, "ERROR:WRONG GRID"
      CRITICAL_ERROR = .TRUE.	  
      RETURN	  
      ENDIF	  
      ALLOCATE(j12_old(states_size_old),m12_old(states_size_old))
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined) THEN
      ALLOCATE(j12_h_old(states_size_old),m12_h_old(states_size_old))
      ENDIF	  
      ALLOCATE(indx_chann_old(states_size_old))
      IF(identical_particles_defined) 
     & ALLOCATE(parity_states_old(states_size_old))  
      SELECT CASE(coll_type)
      CASE(1)
      READ(1,*)
      ALLOCATE(j_ch_old(number_of_channels_old))
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined)
     & ALLOCATE(j_h_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      IF(fine_structure_defined) THEN	!!!! MODIFY
      IF(SPIN_FINE.ne.2) THEN	  
	  READ(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,f_old_b,n_old_b
      ELSE
	  READ(1,'(i8,1x,i8,1x,f4.1,1x,f5.1,1x,f4.1,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),
     & j12_h_old(st),m12_h_old(st),j_h_old_b,f_old_b,par_old_b	  
!!!! TO ADD	  
      ENDIF
      ELSE
	  READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined) THEN
      j_h_ch_old(indx_chann_old(st_old)) = j_h_old_b	
      ELSE	  
      j_ch_old(indx_chann_old(st_old)) = j_old_b  
      ENDIF
      ENDDO		  
      CASE(2)
      READ(1,*)
      ALLOCATE(v_ch_old(number_of_channels_old),
     & j_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),v_old_b,j_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      v_ch_old(indx_chann_old(st_old)) = v_old_b	  
      ENDDO
      CASE(3)
      READ(1,*)
      ALLOCATE(j_ch_old(number_of_channels_old),
     & k_ch_old(number_of_channels_old),
     & eps_ch_old(number_of_channels_old) )	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,k_old_b,eps_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      k_ch_old(indx_chann_old(st_old)) = k_old_b
      eps_ch_old(indx_chann_old(st_old)) = eps_old_b	  
      ENDDO	  
      CASE(4)
      READ(1,*)
      ALLOCATE(j_ch_old(number_of_channels_old),
     & ka_ch_old(number_of_channels_old),
     & kc_ch_old(number_of_channels_old) )	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,ka_old_b,kc_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      ka_ch_old(indx_chann_old(st_old)) = ka_old_b
      kc_ch_old(indx_chann_old(st_old)) = kc_old_b	  
      ENDDO	  	  
      CASE(5)
      READ(1,*)
      IF(.not.identical_particles_defined) THEN	  
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & j1_old_b,j2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b		 
      ENDDO
      ELSE
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & j1_old_b,j2_old_b,par_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Muj FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      parity_states_old(st_old) = par_old_b	  
      ENDDO	  
      ENDIF	  
      CASE(6)
      READ(1,*)
      IF(.not.identical_particles_defined) THEN	  
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))
      ALLOCATE(v1_ch_old(number_of_channels_old),
     & v2_ch_old(number_of_channels_old))	 
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & v1_old_b,v2_old_b,j1_old_b,j2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      v1_ch_old(indx_chann_old(st_old)) = v1_old_b
      v2_ch_old(indx_chann_old(st_old)) = v2_old_b		  
      ENDDO
      ELSE
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & v1_old_b,v2_old_b,j1_old_b,j2_old_b,par_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      v1_ch_old(indx_chann_old(st_old)) = v1_old_b
      v2_ch_old(indx_chann_old(st_old)) = v2_old_b		  
      parity_states_old(st_old) = par_old_b	  
      ENDDO	  
      ENDIF	  
      CASE(7)
      READ(1,*)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & k1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),	 
     & eps1_ch_old(number_of_channels_old) )	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j1_old_b,k1_old_b,eps1_old_b,j2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b	  
      k1_ch_old(indx_chann_old(st_old)) = k1_old_b
      eps1_ch_old(indx_chann_old(st_old)) = eps1_old_b
      ENDDO	  
      CASE(8)
      READ(1,*)	  
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old))
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j1_old_b,ka1_old_b,kc1_old_b,
     & j2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
	  
      ENDDO	 	 
      CASE(9)
      READ(1,*)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old)
     & ,k2_ch_old(number_of_channels_old),
     & eps2_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old
      READ(1,
     & '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,k2_old_b,eps2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      k2_ch_old(indx_chann_old(st_old)) = k2_old_b
      eps2_ch_old(indx_chann_old(st_old)) = eps2_old_b
      ENDDO	  
      CASE(0)
      READ(1,*)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old)
     & ,ka2_ch_old(number_of_channels_old),
     & kc2_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old
      IF(.NOT.identical_particles_defined) THEN	  
      READ(1,
     & '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,ka2_old_b,kc2_old_b
      ELSE
      READ(1,
     & '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i2)
     & ')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,ka2_old_b,kc2_old_b,par_old_b	  
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      ka2_ch_old(indx_chann_old(st_old)) = ka2_old_b
      kc2_ch_old(indx_chann_old(st_old)) = kc2_old_b
      IF(identical_particles_defined)
     & parity_states_old(st_old) = par_old_b	  
      ENDDO	 
      END SELECT
      st = 0    
!!!! TO BE CONTINUED	  
      DO i=1,min(number_of_channels_old,number_of_channels)
      SELECT CASE(coll_type)
      CASE(1)
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) THEN
      DO j_count = 1, int(j_h_ch_old(i)*2) + 1	  
      st = st + 1
      IF(j_count*2 .eq. int(j_h_ch_old(i)*2)+1
     & .and. chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF	 
      ENDDO	  
      ELSE	  
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	  
      ENDDO
      ENDIF	  
	  
      CASE(2)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(3)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(4)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
      ENDDO
      ENDDO
      ELSE
      p_lim_max= 2-KRONEKER(j1_ch(i),j2_ch(i)) !min(p_lim_max_ini,2-KRONEKER(j1_ch(i),j2_ch(i))) 
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch(i),j2_ch(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and. p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  


	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	
      ENDDO
      ENDDO
      ELSE
      p_lim_max= 2-KRONEKER(j1_ch(i),j2_ch(i))*
     & KRONEKER(v1_ch(i),v2_ch(i))  											   										 
!      p_lim_max=min(p_lim_max_ini,2-KRONEKER(j1_ch(i),j2_ch(i))*
!     & KRONEKER(v1_ch(i),v2_ch(i))) 		  
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch(i),j2_ch(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and. p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO 
  

	  
      ENDIF	  
      CASE(7)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(8)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(9)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(0)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      ELSE
!      p_lim_max=min(p_lim_max_ini,2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
!     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
!     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i)))	  
	  p_lim_max = 2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i))
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
!     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
!     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     &	  abs(j1_ch_old(i)-j2_ch_old(i))
     & .and.  p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	   	  
      END SELECT	  
      ENDDO
! HERE WE CHECK THE ORDER !! FOR FINE SPIN=2 WE JUST OVERJUM
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined) GOTO 1994
      DO st=1,min(states_size_old,states_size)
      IF(j12(st).ne.j12_old(st)) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: J12s ARE DIFFERENT"
      RETURN	  
      ENDIF	  
      IF(m12(st).ne.m12_old(st)) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR IN INI: CHECK M12"
      RETURN	  
      ENDIF	  
      SELECT CASE(coll_type)
      CASE(1)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*,"ERROR IN INI : CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      CASE(2)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      IF(v_ch_old(indx_chann_old(st)).ne.v_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI : CHECK V"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	    
      CASE(3)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      IF(k_ch_old(indx_chann_old(st)).ne.k_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK K"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
      IF(eps_ch_old(indx_chann_old(st)).ne.eps_ch(indx_chann(st)))
     & THEN
      PRINT*, "ERROR IN INI: CHECK EPS"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      CASE(4)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka_ch_old(indx_chann_old(st)).ne.ka_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc_ch_old(indx_chann_old(st)).ne.kc_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(5)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(6)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(v1_ch_old(indx_chann_old(st)).ne.v1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK V1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(v2_ch_old(indx_chann_old(st)).ne.v2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK V2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(7)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(k1_ch_old(indx_chann_old(st)).ne.k1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK K1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(eps1_ch_old(indx_chann_old(st)).ne.eps1_ch(indx_chann(st)))
     & THEN
      PRINT*, "ERROR IN INI: CHECK EPS1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(8)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(9)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(k2_ch_old(indx_chann_old(st)).ne.k2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK K2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(eps2_ch_old(indx_chann_old(st)).ne.eps2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK EPS2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(0)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka2_ch_old(indx_chann_old(st)).ne.ka2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc2_ch_old(indx_chann_old(st)).ne.kc2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(identical_particles_defined) THEN	 
      IF(parity_states_old(st)
     & .ne.parity_state(st))
     & STOP "ERROR IN INI: CHECK PARITY"
      ENDIF	 
      END SELECT	  
	  
	  
      ENDDO	 	  

! HERE WE CHECK THE ORDER	  
!      PRINT*, "WE HAVE CKECED"	  
1994      READ(1,*)
!      PRINT*, buffer_word_3
      DO i=1,total_size_old
      READ(1,"(i16,1x,i8,1x,i8)") i_old,ind_mat_old_1,ind_mat_old_2
      IF(i.ne.i_old) CRITICAL_ERROR = .TRUE.
      IF(i.le.total_size) THEN	  
      IF(ind_mat_old_1.ne.ind_mat(1,i)) CRITICAL_ERROR = .TRUE.
      IF(ind_mat_old_2.ne.ind_mat(2,i))	CRITICAL_ERROR = .TRUE.
      ENDIF	  
      ENDDO		  
      READ(1,*)	  
      DO i=1,n_r_coll
      READ(1,'(i5,1x,e19.12)')i_old,R_COM(i)	  
      ENDDO
      READ(1,*)	  
      DO  k=1,min(total_size,total_size_old)
      DO i=1,n_r_coll
      READ(1,'(i16,1x,i8,1x,e19.12)')k_old,i_old,
     & Mat_el(i,k) 	  
      IF(k_old.ne.k) PRINT*,"ERROR IN MIJ: k is wrong"
      IF(i.ne.i_old) PRINT*,"ERROR IN MIJ : ir is wrong"
      ENDDO

! Bikram Start Dec 2019:
	  IF(bikram_mij_shIFt) then
	  do i=1,n_r_coll
	  Mat_el(i,k) = Mat_el(i,k) - Mat_el(n_r_coll,k)		!Bikram
	  ENDdo
	  ENDIF
! Bikram End.

      deriv_bgn = (Mat_el(2,k)
     & - Mat_el(1,k))/
     & (R_COM(2)-R_COM(1)) 
      deriv_END = (Mat_el(n_r_coll,k)
     & - Mat_el(n_r_coll-1,k))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1)) 	 
      CALL  spline(R_COM,Mat_el(:,k),n_r_coll,
     & deriv_bgn,deriv_END,
     & Mat_el_der(:,k))
 
      ENDDO
      CLOSE(1)
	  
      PRINT*,"FORMATTED READING DONE"
      IF(test_expansion_defined)CALL PRINT_ELASTIC_MIJ
      IF(CRITICAL_ERROR) PRINT*,"ERROR IN ASSIGNMENT"	  
      RETURN	  
      END SUBROUTINE READ_MIJ	  

      SUBROUTINE EXCLUDE_STATES
      USE VARIABLES
      USE OLD_MIJ	  
      IMPLICIT NONE
      INTEGER i_state,st_user_exclude,file_stat
      LOGICAL exst 	  
      CHARACTER(LEN=21)::STATES_TO_EXCLUDE="STATES_TO_EXCLUDE.DAT"
      INQUIRE( FILE=STATES_TO_EXCLUDE, EXIST=exst )
      DO i_state = 1,states_size
      stts_to_excl(i_state) = .FALSE.	  
      ENDDO	  
      IF(.not.exst) RETURN	  
      OPEN(UNIT=234,FILE=STATES_TO_EXCLUDE,STATUS="OLD",
     & ACTION = "READ")	  
      READ(234,*,IOSTAT=file_stat) !HEADER
      IF(file_stat.ne.0) THEN
      CLOSE(234)
      RETURN
      ENDIF	  
      DO WHILE(file_stat.eq.0)
      READ(234,*,IOSTAT=file_stat) st_user_exclude
      IF(st_user_exclude.gt.states_size) THEN
      CLOSE(234)	  
      RETURN
      ENDIF	   
      stts_to_excl(st_user_exclude) = .TRUE.	  
      ENDDO	  
      RETURN	  
      END SUBROUTINE EXCLUDE_STATES	  
      SUBROUTINE READ_MIJ_USER
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      USE OLD_MIJ
      IMPLICIT NONE
      REAL*8 Mij_buffer	  
      INTEGER i,k,st,tot_siz_read,st_siz_read,j_count,j_summ,p_count
      INTEGER istat,p_lim_max_ini	  
      LOGICAL file_exst, bk_exst
	  LOGICAL EVEN_NUM 			! Dulat 9/25/2023										
	  CHARACTER(LEN = 500) bk_matrix_path1, bk_matrix_path2
      INTEGER KRONEKER
      CHARACTER (LEN=12) ::
     & MIJ_FILE_NAME_N ="Mij_UF_N.dat"
      CHARACTER (LEN=10) ::
     & MIJ_FILE_NAME ="Mij_UF.dat"       	  
      EXTERNAL KRONEKER	  
	  
! Bikram Start Nov 2021: Changing matrix file path 
!so that user don't need to copy the large matrix file in each directory.
	  IF(bikram_mtrx_path) then
	  inquire(file = "Matrix_Path.DAT", exist = bk_exst)
	  IF(.not.bk_exst) then
	  write(*,'(a)')"File containing matrix path is not provided."
	  stop
	  ELSE
	  open(111,file="Matrix_Path.DAT",status="old",action = "read")
	  read(111,'(a)') bk_matrix_path1
	  close(111)
	  write(*,'(a,a,a)') 'Matrix file path is defined. ',
     & 'Matrix reading from ', trim(bk_matrix_path1)
	  write(bk_matrix_path2,'(a,a,a)')trim(bk_matrix_path1),'/',
     & trim(MATRIX_NAME_MIJ_UF)
	  END IF
	  END IF
	  
	  IF(.not.bikram_mtrx_path) then	  
	  INQUIRE( FILE=MATRIX_NAME_MIJ_UF, EXIST=file_exst )
	  ELSE
	  INQUIRE( FILE=trim(bk_matrix_path2), EXIST=file_exst )
	  END IF
! Bikram End.
	  
!      INQUIRE( FILE=MATRIX_NAME_MIJ_UF, EXIST=file_exst )
      IF(.not.file_exst) THEN
      PRINT*,"ERROR: MTRX_UF.DAT NOT FOUND"
      CRITICAL_ERROR = .TRUE.	  
      RETURN	  
      ENDIF
      IF(.not.identical_particles_defined) THEN	  
      p_lim_max = 1
      p_lim_min = 1	
      ELSE
      SELECT CASE(coll_type)		  
      CASE(5)	  

      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
      CASE(6)
      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))*
     & KRONEKER(v1_ch(chann_ini),v2_ch(chann_ini)) 	  
      CASE(0)
      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
     & *KRONEKER(ka1_ch(chann_ini),ka2_ch(chann_ini))
     & *KRONEKER(kc1_ch(chann_ini),kc2_ch(chann_ini))	  
      END SELECT		  
      ENDIF		   
	  
! Bikram Start Nov 2021: 
	  IF(.not.bikram_mtrx_path) then	  
	  OPEN(1,FILE=MATRIX_NAME_MIJ_UF,STATUS="OLD",ACTION="READ",
     & form = 'unformatted')
	  ELSE
	  OPEN(1,FILE=trim(bk_matrix_path2), STATUS="OLD",ACTION="READ",
     & form = 'unformatted')
	  END IF
! Bikram End.

!      OPEN(1,FORM="UNFORMATTED",
!     % FILE=MATRIX_NAME_MIJ_UF,STATUS="OLD",ACTION="READ")
      PRINT*,"MTRX_UF.DAT READING STARTED"	  
      READ(1,IOSTAT=istat) coll_type_old
      IF(coll_type_old.ne.coll_type) THEN 
      PRINT*, "ERROR: SYSTEM IS DIFFERENT"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
      READ(1,IOSTAT=istat) number_of_channels_old	 
      READ(1,IOSTAT=istat) states_size_old
      READ(1,IOSTAT=istat) total_size_old	  
      READ(1,IOSTAT=istat)n_r_coll_old
      IF(n_r_coll_old.ne.n_r_coll) THEN
      n_r_coll_old = n_r_coll
!      CRITICAL_ERROR = .TRUE.
      PRINT*,"WRONG GRID", n_r_coll_old, n_r_coll
!      RETURN	  
      ENDIF	  
      IF(istat.ne.0) THEN
      CRITICAL_ERROR = .TRUE.
      PRINT*,"ERROR in ",MATRIX_NAME_MIJ_UF
      RETURN	  
      ENDIF	  
!      IF(total_size_old.lt.total_size) THEN
!      CALL SYSTEM('cp '//
!     & MIJ_FILE_NAME//' '//MIJ_FILE_NAME_N)	  
!      ENDIF
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) GOTO 1992	  
      ALLOCATE(j12_old(states_size_old),m12_old(states_size_old))
      ALLOCATE(indx_chann_old(states_size_old))
      IF(identical_particles_defined) 
     & ALLOCATE(parity_states_old(states_size_old))  
      SELECT CASE(coll_type)
      CASE(1)
      ALLOCATE(j_ch_old(number_of_channels_old))     
      DO st=1,states_size_old
      READ(1)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	  
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      ENDDO		  
      CASE(2)
      ALLOCATE(v_ch_old(number_of_channels_old),
     & j_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),v_old_b,j_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      v_ch_old(indx_chann_old(st_old)) = v_old_b	  
      ENDDO
      CASE(3)
      ALLOCATE(k_ch_old(number_of_channels_old),
     & j_ch_old(number_of_channels_old),
     & eps_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,k_old_b,eps_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      k_ch_old(indx_chann_old(st_old)) = k_old_b
      eps_ch_old(indx_chann_old(st_old)) = eps_old_b
      ENDDO	  
      CASE(4)	  
      ALLOCATE(ka_ch_old(number_of_channels_old),
     & j_ch_old(number_of_channels_old),
     & kc_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,ka_old_b,kc_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      ka_ch_old(indx_chann_old(st_old)) = ka_old_b
      kc_ch_old(indx_chann_old(st_old)) = kc_old_b
      ENDDO	  
      CASE(5)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))
      DO st=1,states_size_old
      IF(.not.identical_particles_defined) THEN	 
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & j1_old_b,j2_old_b
      ELSE
      READ(1,IOSTAT=istat)	  
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & j1_old_b,j2_old_b,par_old_b
      parity_states_old(st_old) = par_old_b	 
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b		 
      ENDDO
      CASE(6)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))
      ALLOCATE(v1_ch_old(number_of_channels_old),
     & v2_ch_old(number_of_channels_old))	 
      DO st=1,states_size_old
      IF(.not.identical_particles_defined) THEN	 
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & v1_old_b,j1_old_b,v2_old_b,j2_old_b
      ELSE
      READ(1,IOSTAT=istat)	  
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & v1_old_b,j1_old_b,v2_old_b,j2_old_b,par_old_b
      parity_states_old(st_old) = par_old_b	 
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      v1_ch_old(indx_chann_old(st_old)) = v1_old_b
      v2_ch_old(indx_chann_old(st_old)) = v2_old_b	  
      ENDDO
      CASE(7)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & k1_ch_old(number_of_channels_old),
     & eps1_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old

      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,k1_old_b,eps1_old_b
     & ,j2_old_b
	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      k1_ch_old(indx_chann_old(st_old)) = k1_old_b
      eps1_ch_old(indx_chann_old(st_old)) = eps1_old_b
      ENDDO		  
      CASE(8)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old

      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b
	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      ENDDO	 	  	  
      CASE(9)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old)
     & ,k2_ch_old(number_of_channels_old),
     & eps2_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old

      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,k2_old_b,eps2_old_b
	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      k2_ch_old(indx_chann_old(st_old)) = k2_old_b
      eps2_ch_old(indx_chann_old(st_old)) = eps2_old_b
      ENDDO	 	  
      CASE(0)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old)
     & ,ka2_ch_old(number_of_channels_old),
     & kc2_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old
      IF(.NOT.identical_particles_defined) THEN	  
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,ka2_old_b,kc2_old_b
      ELSE
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,ka2_old_b,kc2_old_b,par_old_b	  
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      ka2_ch_old(indx_chann_old(st_old)) = ka2_old_b
      kc2_ch_old(indx_chann_old(st_old)) = kc2_old_b
      IF(identical_particles_defined)
     & parity_states_old(st_old) = par_old_b	  
      ENDDO	 
      END SELECT
      st = 0       
      DO i=1,min(number_of_channels_old,number_of_channels)
      SELECT CASE(coll_type)
      CASE(1)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	  
      ENDDO
      CASE(2)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(3)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(4)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
      ENDDO
      ENDDO
      ELSE
      p_lim_max= 2-KRONEKER(j1_ch(i),j2_ch(i)) !min(p_lim_max_ini,2-KRONEKER(j1_ch(i),j2_ch(i)))	  	  
	  
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch(i),j2_ch(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and. p_count .eq.1  ) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	
      ENDDO
      ENDDO
      ELSE
      p_lim_max= 2-KRONEKER(j1_ch(i),j2_ch(i))*
     & KRONEKER(v1_ch(i),v2_ch(i))
!	 p_lim_max=min(p_lim_max_ini,2-KRONEKER(j1_ch(i),j2_ch(i))*
!     & KRONEKER(v1_ch(i),v2_ch(i)))
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch(i),j2_ch(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and. p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO 
  

	  
      ENDIF	  
      CASE(7)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(8)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(9)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(0)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      ELSE
      p_lim_max=2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i))      
!	 p_lim_max=min(p_lim_max_ini,2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
!     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
!     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i)))		
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
!     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
!     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     &	  abs(j1_ch_old(i)-j2_ch_old(i))
     & .and.  p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	   	  
      END SELECT	  
      ENDDO
! HERE WE CHECK THE ORDER 
      DO st=1,min(states_size_old,states_size)
      IF(j12(st).ne.j12_old(st)) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: J12s ARE DIFFERENT"
      RETURN	  
      ENDIF	  
      IF(m12(st).ne.m12_old(st)) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR IN INI: CHECK M12"
      RETURN	  
      ENDIF	  
      SELECT CASE(coll_type)
      CASE(1)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*,"ERROR IN INI : CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      CASE(2)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      IF(v_ch_old(indx_chann_old(st)).ne.v_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI : CHECK V"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	    
      CASE(3)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      IF(k_ch_old(indx_chann_old(st)).ne.k_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK K"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
      IF(eps_ch_old(indx_chann_old(st)).ne.eps_ch(indx_chann(st)))
     & THEN
      PRINT*, "ERROR IN INI: CHECK EPS"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      CASE(4)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka_ch_old(indx_chann_old(st)).ne.ka_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc_ch_old(indx_chann_old(st)).ne.kc_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(5)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(6)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(v1_ch_old(indx_chann_old(st)).ne.v1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK V1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(v2_ch_old(indx_chann_old(st)).ne.v2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK V2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(7)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(k1_ch_old(indx_chann_old(st)).ne.k1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK K1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(eps1_ch_old(indx_chann_old(st)).ne.eps1_ch(indx_chann(st)))
     & THEN
      PRINT*, "ERROR IN INI: CHECK EPS1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(8)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(9)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(k2_ch_old(indx_chann_old(st)).ne.k2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK K2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(eps2_ch_old(indx_chann_old(st)).ne.eps2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK EPS2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(0)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka2_ch_old(indx_chann_old(st)).ne.ka2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc2_ch_old(indx_chann_old(st)).ne.kc2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(identical_particles_defined) THEN	 
      IF(parity_states_old(st)
     & .ne.parity_state(st))
     & STOP "ERROR IN INI: CHECK PARITY"
      ENDIF	 
      END SELECT	  
	  
	  
      ENDDO	  

! HERE WE CHECK THE ORDER	  
1992  DO i=1,total_size_old
      READ(1) i_old,ind_mat_old_1,ind_mat_old_2
      IF(i.ne.i_old) STOP "ERROR IN INI: CHECK FILE"
      IF(i.le.total_size) THEN	  
      IF(ind_mat_old_1.ne.ind_mat(1,i)) THEN
      PRINT*, "ERROR IN INI:IND_MAT_1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
      IF(ind_mat_old_2.ne.ind_mat(2,i))	THEN
      PRINT*,"ERROR IN INI:IND_MAT_2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
      ENDIF	  
      ENDDO		  
      DO i=1,n_r_coll
      READ(1)i_old,R_COM(i)	  
      ENDDO
      DO  k=1,min(total_size,total_size_old)
      READ(1) k_old
!      DO i=1,n_r_coll
      READ(1) Mat_el(:,k) 	  
      IF(k_old.ne.k) THEN
      PRINT*,"ERROR IN MIJ: k is wrong"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
!      IF(i.ne.i_old) PRINT*,"ERROR IN MIJ : ir is wrong"
!      ENDDO
! SPLINING	  


! Bikram Start Dec 2019:
	  IF(bikram_mij_shIFt) then
	  do i=1,n_r_coll
	  Mat_el(i,k) = Mat_el(i,k) - Mat_el(n_r_coll,k)		!Bikram
	  ENDdo
	  ENDIF
! Bikram End.

      deriv_bgn = (Mat_el(2,k)
     & - Mat_el(1,k))/
     & (R_COM(2)-R_COM(1)) 
      deriv_END = (Mat_el(n_r_coll,k)
     & - Mat_el(n_r_coll-1,k))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1))
!      deriv_END = 0d0
!      deriv_bgn = 0d0	  
      CALL  spline(R_COM,Mat_el(:,k),n_r_coll,
     & deriv_bgn,deriv_END,
     & Mat_el_der(:,k))
      ENDDO
      CLOSE(1)
      IF(test_expansion_defined)CALL PRINT_ELASTIC_MIJ	  
      PRINT*,"UNFORMATTED READING DONE"
      !!! TESTING
!      CALL SYSTEM('rm '// MIJ_FILE_NAME_N)
!      OPEN(333,file="LAST_MIJ.DAT")	            !!!!! DELETE AFTER ALL
!	  DO k=1,total_size
!	  IF(ind_mat(1,k).eq.ind_mat(2,k)) THEN
!	  WRITE(333,*) k,Mat_el(n_r_coll,k)
!      ENDIF
!      ENDDO	
!      CLOSE(333)		  
      RETURN	  
      !!! TESTING	  
      END SUBROUTINE READ_MIJ_USER 
!!!!!!!!!!!!! PREPARE IT FOR MPI TASK PER TRAJECTORY
     
      SUBROUTINE PRINT_MIJ_USER
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      USE OLD_MIJ
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      IMPLICIT NONE
      INTEGER st,i,k
      CHARACTER (LEN=12) ::
     & MIJ_FILE_NAME_N ="Mij_UF_N.dat"
      CHARACTER (LEN=10) ::
     & MIJ_FILE_NAME ="Mij_UF.dat"    	  
      IF(MYID.ne.0) RETURN	  
      PRINT*,"SYSTEM_SETUP_DONE"
      OPEN(111,FORM="UNFORMATTED",FILE=MATRIX_NAME_MIJ_UF,
     & ACTION="WRITE")
      WRITE(111)coll_type
      WRITE(111)number_of_channels
      WRITE(111)states_size
      WRITE(111)total_size
      WRITE(111)n_r_coll
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) GOTO 1993	  
      SELECT CASE(coll_type)
      CASE(1)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j_ch(i)	  
      ENDDO		  
      CASE(2)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),v_ch(i),j_ch(i)
      ENDDO	
      CASE(3)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j_ch(i),k_ch(i),eps_ch(i)
      ENDDO
      CASE(4)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j_ch(i),ka_ch(i),kc_ch(i)
      ENDDO		  
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i),parity_state(st)	  
      ENDDO	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined)	THEN  
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),v1_ch(i),j1_ch(i),v2_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),v1_ch(i),j1_ch(i),v2_ch(i),j2_ch(i),
     & parity_state(st)	  
      ENDDO	  
      ENDIF
      CASE(7)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j1_ch(i),k1_ch(i),eps1_ch(i)
     & ,j2_ch(i)	 
      ENDDO
      CASE(8)	  
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i)	 
      ENDDO
      CASE(9)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),kc2_ch(i),eps2_ch(i)	 
      ENDDO	  
      CASE(0)
      IF(.NOT.identical_particles_defined) THEN  
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i)	 
      ENDDO
      ELSE
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i),parity_state(st)
      ENDDO	 
      ENDIF	 
      END SELECT
1993  DO i=1,total_size
      WRITE(111) i,ind_mat(1,i),ind_mat(2,i)	  
      ENDDO		  
      DO i=1,n_r_coll
      WRITE(111)i,R_COM(i)	  
      ENDDO
      DO k=1,total_size
      WRITE(111) k	  
!      DO i=1,n_r_coll	  
      WRITE(111) Mat_el(:,k)!,
!      ENDDO
      ENDDO
      CLOSE(111)
!      CALL SYSTEM('rm '// MIJ_FILE_NAME_N)
!      OPEN(333,file="LAST_MIJ.DAT")	            !!!!! DELETE AFTER ALL
!	  DO k=1,total_size
!	  IF(ind_mat(1,k).eq.ind_mat(2,k)) THEN
!	  WRITE(333,*) k,Mat_el(n_r_coll,k)
!      ENDIF
!      ENDDO	
!      CLOSE(333)	
      IF(test_expansion_defined)CALL PRINT_ELASTIC_MIJ  
      END SUBROUTINE PRINT_MIJ_USER
      SUBROUTINE PRINT_ELASTIC_MIJ
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      USE OLD_MIJ
      IMPLICIT NONE
      INTEGER i,j,k
      INTEGER non_zero_size	  
      REAL*8, ALLOCATABLE :: Mij_elast(:,:),
     & Mij_elast_cs(:,:)
      INTEGER, ALLOCATABLE :: index_elastic_corr(:)
      CHARACTER(LEN=22) :: elast_mij_out = "ELASTIC_ELEMENTS  .out"	  
      ALLOCATE(Mij_elast(n_r_coll,states_size))
      non_zero_size = 0	  
      DO k=1,total_size
      IF(ind_mat(1,k).eq.ind_mat(2,k)) THEN
      IF(m12(ind_mat(1,k)).eq.m_elastic_proj_print)
     & non_zero_size = non_zero_size +1	  
      Mij_elast(:,ind_mat(1,k)) = Mat_el(:,k)	  
      ENDIF	  
      ENDDO
      ALLOCATE(Mij_elast_cs(n_r_coll,non_zero_size))
	  ALLOCATE(index_elastic_corr(non_zero_size))
      non_zero_size = 0	  
      DO k=1,total_size
      IF(ind_mat(1,k).eq.ind_mat(2,k)) THEN
      IF(m12(ind_mat(1,k)).eq.m_elastic_proj_print) THEN
      non_zero_size = non_zero_size +1	  
      Mij_elast_cs(:,non_zero_size) = Mat_el(:,k)
      index_elastic_corr(non_zero_size) = k	  
      ENDIF	  
      ENDIF	  
      ENDDO	  
      WRITE(elast_mij_out(17:18),'(i2.2)')m_elastic_proj_print	  
      OPEN(234,FILE=elast_mij_out)
       DO j=1,non_zero_size
      IF(j.eq.1) WRITE(234,"(a6)",ADVANCE="NO") "R_COM" 
      k = 	index_elastic_corr(j)  
      WRITE(234,"(2x,i4,2x,i3,1x,i3,2x)",ADVANCE="NO") k,
     & ind_mat(1,k),ind_mat(2,k)
      ENDDO
      WRITE(234,*)	  
	  
      DO i=1,n_r_coll
      DO j=1,non_zero_size
      IF(j.eq.1) WRITE(234,"(f6.2)",ADVANCE="NO") R_COM(i) 	  
      WRITE(234,"(e17.5)",ADVANCE="NO") Mij_elast_cs(i,j)*autoeV*eVtown
      ENDDO
      WRITE(234,*)	  
      ENDDO	  
      CLOSE(234)	  
	  
      END SUBROUTINE PRINT_ELASTIC_MIJ	  
	  
	   SUBROUTINE bk_print_matrix_info
! This subroutine is created by Bikramaditya Mandal, Nov 2020
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      IMPLICIT NONE
      INTEGER st,i,k
	  CHARACTER (LEN=100) :: bk_dir_temp_1,bk_dir_temp_2
	  CHARACTER (LEN=100) :: bk_dir_temp_5
	  
      PRINT*,"SYSTEM_SETUP_DONE"
	  IF(.not.unformat_defined) then
	  write(bk_dir_temp_1, '(a,a)') trim(bk_dir2), "/MATRIX_INFO.DAT"
	  bk_dir_temp_2 = trim(bk_dir_temp_1)
	  OPEN(1,FILE=bk_dir_temp_2,ACTION="WRITE")
      WRITE(1,'(a15,1x,i1)') "COLLISION_TYPE=",coll_type
      WRITE(1,'(a17,x,i4)') "NUMBER_OF_CHANLS=",number_of_channels
      WRITE(1,'(a17,1x,i6)') "NUMBER_OF_STATES=",states_size
      WRITE(1,'(a12,1x,i9)') "MATRIX_SIZE=",total_size
      WRITE(1,'(a12,1x,i4)') "R_GRID_SIZE=",n_r_coll	  
      SELECT CASE(coll_type)
      CASE(1)
      IF(.not.fine_structure_defined) THEN	  
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j_ch(i)	  
      ENDDO
      ELSE
      IF(SPIN_FINE.ne.2) THEN       	  
      WRITE(1,'(a8,1x,a8,1x,a3,1x,a3,1x,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J","F","N"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i3,1x,i3,1x,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j_ch(i),2+f_ch(i),j_ch(i)+f_ch(i)	  
      ENDDO
      ELSE
      WRITE(1,'(a8,1x,a8,1x,a4,1x,a5,1x,1x,a4,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J","F","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,f4.1,1x,f5.1,1x,1x,f4.1,1x,i3,1x,i3)')
     & st,i,j12_h(st),m12_h(st),j_h_ch(i),f_ch(i),par_lorb_ch(i)	  !!!! TO MODIFY
      ENDDO	  
      ENDIF	  
	  
      ENDIF	  
      CASE(2)
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","V","J"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),v_ch(i),j_ch(i)	  
      ENDDO
      CASE(3)
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J","K", "EPS"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j_ch(i),k_ch(i),eps_ch(i)	  
      ENDDO
      CASE(4)
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J","KA", "KC"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j_ch(i),ka_ch(i),kc_ch(i)	  
      ENDDO      	  
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J1","J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J1","J2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i),parity_state(st)	  
      ENDDO	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined) THEN	  
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","V1","V2","J1","J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),v1_ch(i),v2_ch(i),j1_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      WRITE(1,'(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","V1","V2","J1","J2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),v1_ch(i),v2_ch(i),j1_ch(i),j2_ch(i),
     & parity_state(st)	 
      ENDDO	  
      ENDIF	  
      CASE(7)
      WRITE(1,
     & '(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J1","K","EPS",
     & "J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j1_ch(i),k1_ch(i),eps1_ch(i)
     & ,j2_ch(i) 
      ENDDO	  
      CASE(8)
      WRITE(1,
     & '(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J1","K1A","K1C",
     & "J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i) 
      ENDDO		  
      CASE(9)
      WRITE(1,
     & '(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J1","K1A","K1C",
     & "J2","K2","EPS"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),k2_ch(i),eps2_ch(i) 
      ENDDO		  
      CASE(0)
      IF(.NOT.identical_particles_defined) THEN  
      WRITE(1,
     & '(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8)')
     & "#STATE","#CHANNEL","J12","M12","J1","KA1","KC1",
     & "J2","KA2","KC2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i)	 
      ENDDO
      ELSE
      WRITE(1,
     & '(a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a8,1x,a2)
     & ')
     & "#STATE","#CHANNEL","J12","M12","J1","KA1","KC1",
     & "J2","KA2","KC2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i2)
     & ')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i),parity_state(st)
      ENDDO	 
      ENDIF	 
      END SELECT
      WRITE(1,"(a16,1x,a8,1x,a8)") "MIJ_INDEX", "ST_1", "ST_2"
      DO i=1,total_size
      WRITE(1,"(i16,1x,i8,1x,i8)") i,ind_mat(1,i),ind_mat(2,i)	  
      ENDDO		  
      WRITE(1,'(a5,1x,a19)') "#iR","R_COM(iR)"	  
      DO i=1,n_r_coll
      WRITE(1,'(i5,1x,e19.12)')i,R_COM(i)	  
      ENDDO
	  WRITE(1,*)
	  close(1)
	  
	  call getcwd(bk_dir_temp_5)
	  call chdir(trim(bk_dir2))
!	  call system ( "chmod +x MATRIX_COMBINE.sh" )
!	  call system ( "./MATRIX_COMBINE.sh" )
	  call chdir(trim(bk_dir_temp_5))
!	  print *, 'need work to do in combine'

	  ELSE
	  write(bk_dir_temp_1, '(a,a)') trim(bk_dir2), "/MATRIX_INFO.DAT"
	  bk_dir_temp_2 = trim(bk_dir_temp_1)
	  OPEN(1,FILE=bk_dir_temp_2,ACTION="WRITE",form="unformatted")
      WRITE(1)coll_type
      WRITE(1)number_of_channels
      WRITE(1)states_size
      WRITE(1)total_size
      WRITE(1)n_r_coll
      SELECT CASE(coll_type)
      CASE(1)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1)
     & st,i,j12(st),m12(st),j_ch(i)	  
      ENDDO		  
      CASE(2)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1)
     & st,i,j12(st),m12(st),v_ch(i),j_ch(i)
      ENDDO	
      CASE(3)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1)
     & st,i,j12(st),m12(st),j_ch(i),k_ch(i),eps_ch(i)
      ENDDO
      CASE(4)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1)
     & st,i,j12(st),m12(st),j_ch(i),ka_ch(i),kc_ch(i)
      ENDDO		  
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1)
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1)
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i),parity_state(st)	  
      ENDDO	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined)	THEN  
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1)
     & st,i,j12(st),m12(st),v1_ch(i),j1_ch(i),v2_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1)
     & st,i,j12(st),m12(st),v1_ch(i),j1_ch(i),v2_ch(i),j2_ch(i),
     & parity_state(st)	  
      ENDDO	  
      ENDIF
      CASE(7)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1)
     & st,i,j12(st),m12(st),j1_ch(i),k1_ch(i),eps1_ch(i)
     & ,j2_ch(i)	 
      ENDDO
      CASE(8)	  
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1)
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i)	 
      ENDDO
      CASE(9)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1)
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),kc2_ch(i),eps2_ch(i)	 
      ENDDO	  
      CASE(0)
      IF(.NOT.identical_particles_defined) THEN  
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1)
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i)	 
      ENDDO
      ELSE
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1)
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i),parity_state(st)
      ENDDO	 
      ENDIF	 
      END SELECT
!      DO i=1,total_size
!      WRITE(1) i,ind_mat(1,i),ind_mat(2,i)	  
!      ENDDO		  
      DO i=1,n_r_coll
      WRITE(1)i,R_COM(i)
      ENDDO

	  END IF
	  
      IF(test_expansion_defined) THEN 
      CALL PRINT_ELASTIC_MIJ
      PRINT*, "ELASTIC_TERMS_FILE_CREATED"  
      ENDIF	 
      END SUBROUTINE
	  
	   SUBROUTINE bk_print_matrix(bk_n1, bk_n2, 
     & bk_nn, bk_ncount, bk_mat_array, mat_chk, tmp1)
!--------------------------------------------------------------------
! This subroutine is created by Bikramaditya Mandal, Nov 2020
! This subroutine is called after 1000 matrix elements are computed
! and are ready to be stored in the file.
!--------------------------------------------------------------------
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      IMPLICIT NONE
      INTEGER st,ii,kk
	  INTEGER bk_nn, bk_n1, bk_n2, bk_ncount, mat_chk
	  REAL*8 bk_mat_array(n_r_coll,bk_ncount)
	  INTEGER tmp1(bk_ncount)
	  CHARACTER (LEN=1000) :: bk_dir_tmp, bk_dir_mtrx, bk_dir_ind
	  CHARACTER (LEN=1000) :: bk_rebalance

	  IF(.not.unformat_defined) then
	  write(bk_dir_tmp, '(a,a,i0,a,i0,a)') 
     & trim(bk_dir2), "/MIJ_",bk_n1,"_",bk_n2,".DAT"
	  bk_dir_mtrx = trim(bk_dir_tmp)
	  OPEN(1,FILE=trim(bk_dir_mtrx),ACTION="WRITE",access="appEND")
	  write(bk_dir_tmp, '(a,a,i0,a,i0,a)') 
     & trim(bk_dir2), "/Ind_",bk_n1,"_",bk_n2,".DAT"
	  bk_dir_ind = trim(bk_dir_tmp)
	  OPEN(11,FILE=trim(bk_dir_ind),ACTION="WRITE",access="appEND")
	  IF(bikram_rebalance) then
	  write(bk_dir_tmp, '(a,a,i0,a,i0,a)') 
     & trim(bk_dir33), "/Info_",bk_n1,"_",bk_n2,".DAT"
	  bk_rebalance = trim(bk_dir_tmp)
	  open(12,file=trim(bk_rebalance),ACTION="WRITE",access="appEND")
	  END IF
	  
! Then your original loop with bounds checking											  
	  IF(bikram_rebalance_comp) then
	  do kk = 1, bk_ncount
      DO ii = 1, n_r_coll	  
      WRITE(1,'(i16,1x,i8,1x,e19.12)') tmp1(kk),ii,bk_mat_array(ii,kk)
      ENDDO
	  write(11,*)ind_mat( 1, tmp1(kk) ), ind_mat( 2, tmp1(kk) )
      ENDDO
      CLOSE(1)
      CLOSE(11)
	  return
	  END IF
	  
	  do kk = 1, bk_ncount
	  IF(max(abs(bk_mat_array(mtrx_cutoff_r1,kk)),
     & abs(bk_mat_array(mtrx_cutoff_r2,kk))).gt.MIJ_ZERO_CUT) then
      DO ii=1,n_r_coll	  
      WRITE(1,'(i16,1x,i8,1x,e19.12)')(kk-1)+mat_chk,ii,
     & bk_mat_array(ii,kk)
      ENDDO
	  write(11,*)ind_mat(1,(kk-1)+mat_chk), ind_mat(2,(kk-1)+mat_chk)
	  IF(bikram_rebalance) write(12,*)(kk-1)+mat_chk, tmp1(kk)
	  END IF
      ENDDO
      CLOSE(1)
      CLOSE(11)
      IF(bikram_rebalance) CLOSE(12)
	  
	  ELSE
	  write(bk_dir_tmp, '(a,a,i0,a,i0,a)') 
     & trim(bk_dir2), "/MIJ_",bk_n1,"_",bk_n2,".DAT"
	  bk_dir_mtrx = trim(bk_dir_tmp)
	  OPEN(1,FILE=bk_dir_mtrx,ACTION="WRITE",form="unformatted"
     & ,access="appEND")
	  write(bk_dir_tmp, '(a,a,i0,a,i0,a)') 
     & trim(bk_dir2), "/Ind_",bk_n1,"_",bk_n2,".DAT"
	  bk_dir_ind = trim(bk_dir_tmp)
	  OPEN(11,FILE=bk_dir_ind,ACTION="WRITE",form="unformatted"
     & ,access="appEND")
	  
	  IF(bikram_rebalance_comp) then
	  do kk = 1, bk_ncount
      WRITE(1) tmp1(kk)
      WRITE(1) bk_mat_array(:,kk)
	  write(11)ind_mat( :, tmp1(kk) )
      ENDDO
      CLOSE(1)
      CLOSE(11)
	  return
	  END IF
	  
	  do kk = 1, bk_ncount
	  IF(max(abs(bk_mat_array(mtrx_cutoff_r1,kk)),
     & abs(bk_mat_array(mtrx_cutoff_r2,kk))).gt.MIJ_ZERO_CUT) then
      WRITE(1) kk - 1 + mat_chk
      WRITE(1) bk_mat_array(:,kk)
	  write(11)ind_mat(:,(kk-1)+mat_chk)
	  END IF
      ENDDO
      CLOSE(1)
      CLOSE(11)
	  END IF
      END SUBROUTINE
	  
	   SUBROUTINE bk_read_matrix(bk_nn,
     & bk_mat_array, bk_ind_mat)
! This subroutine is created by Bikramaditya Mandal, Nov 2020 
! to read matrix files for PARALLEL_IO input
!--------------------------------------------------------------------
! Parallel_IO matrix calculations involve 4 major files. All these files 
! are in directory "MATRIX_TRUNCATED". One type starts with "Ind_*", and 
! another type with "MIJ_*". The number of these files are equal to number of procs.
! Two other important file names are MATRIX_FILE_INDEX.DAT and MATRIX_NONZERO_INDEX.DAT
!--------------------------------------------------------------------

      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      IMPLICIT NONE
      INTEGER st, ii, iii, bk_nn, bk_nnn, bk_st, bk_fn, jj, jjj
	  INTEGER file_counter, file_nmb, file_io, file_nmb1, file_nmb2
	  INTEGER dm1, dm2, dm3, st_store, bk_non_zero_counter_old,fc_old
	  REAL*8 dm4, MIJ_ZERO_CUT_old, dm9(n_r_coll)
	  INTEGER, allocatable ::file_nmb_bgn(:),file_nmb_END(:),nmb_m(:)
	  REAL*8 bk_mat_array(n_r_coll, bk_nn)
	  REAL*8 db_time1, db_time2, db_time										! Dulat 12/2/2024 - matrix reading time check
	  INTEGER bk_ind_mat(2, bk_nn), wrk_id, bk_nz, bk_dm, bk_dm1
	  INTEGER proc_start, proc_remain, bk_tmp_ind_mat(2), dm10(2)
	  CHARACTER(LEN=500) bk_matrix_path1, bk_matrix_path2
	  CHARACTER (LEN=500) :: bk_dir_temp_1,bk_dir_temp_2
	  CHARACTER (LEN=500) :: bk_dir_temp_3,bk_dir_temp_4
	  CHARACTER (LEN=28) :: dm5
	  CHARACTER (LEN=21) :: dm6
	  CHARACTER (LEN=9) :: dm7
	  CHARACTER (LEN=15) :: dm8
	  logical bk_exst


      db_time1 = mpi_wtime()													! Dulat 12/2/2024 - matrix reading time check																													
	  IF(myid.eq.0) write(*,*)"Parallel Matrix Reading Started"
!--------------------------------------------------------------------
! Read information about the PARALLEL_IO output performed during the previous run
!--------------------------------------------------------------------
	  write(bk_dir_temp_1, '(a,a)') trim(bk_dir2),
     & "/MATRIX_NONZERO_INDEX.DAT"
	  bk_dir_temp_2 = trim(bk_dir_temp_1)
	  open(1,file = bk_dir_temp_2, action = 'read')
	  read(1,'(a9,i16)') dm7, file_counter										! total number of files of Ind* and MIJ*
	  read(1,'(a28,i16)') dm5, bk_non_zero_counter_old							! total number of non-zero matrix elemnts
	  read(1,'(a21,e19.12)') dm6, MIJ_ZERO_CUT_old								! cut-off value used for truncation of this matrix
	  read(1,*)
!--------------------------------------------------------------------
! Check whether the cut-off value is corect. If not, then stop.
!--------------------------------------------------------------------
	  IF(MIJ_ZERO_CUT_old.ne.MIJ_ZERO_CUT) then
	  print*, "The truncated Matrix does not match the Cut-Off Value."
	  print*, "Please check and provide correct truncated Matrix."
	  print*, MIJ_ZERO_CUT_old, MIJ_ZERO_CUT
	  stop
	  return
	  END IF
!--------------------------------------------------------------------
! Read the number of non-zero matrix elements in each individual MIJ_* file
!--------------------------------------------------------------------
	  allocate(nmb_m(file_counter))
	  do ii = 1, file_counter
	  read(1,'(3(i16,2x))') dm2, nmb_m(ii)
	  IF(dm2.ne.ii) then
	  print*, "Error in getting matrix index #", ii, dm2
	  stop
	  END IF
	  END do
	  close(1)
!--------------------------------------------------------------------
! Computing the number of non-zero matrix elements needed to read 
! by each procs depENDing on MPI_PERTRAJ
!--------------------------------------------------------------------
	  wrk_id = myid - int(myid/mpi_task_per_traject)
     & *mpi_task_per_traject
	  bk_nz = bk_non_zero_counter_old/mpi_task_per_traject
	  IF(wrk_id.lt.mod(bk_non_zero_counter_old,
     & mpi_task_per_traject)) bk_nz = bk_nz + 1
!--------------------------------------------------------------------
! Read information about which matrix elements are stored in each file
! Check that the info is consistent. If not, then stop
!--------------------------------------------------------------------
	  write(bk_dir_temp_1, '(a,a)') trim(bk_dir2),
     & "/MATRIX_FILE_INDEX.DAT"
	  bk_dir_temp_2 = trim(bk_dir_temp_1)
	  open(1,file = bk_dir_temp_2, action = 'read')
	  read(1,'(a15,i10)') dm8, fc_old
	  read(1,'(a15,i10)') 
	  IF(fc_old.ne.file_counter) then
	  print*, "#Files do not match in matrix directories"
	  stop
	  END IF
!--------------------------------------------------------------------
! Read indices for first and last matrix elements stored in this file
!--------------------------------------------------------------------
	  allocate(file_nmb_bgn(file_counter),file_nmb_END(file_counter))
	  do ii = 1, file_counter
	  read(1,'(3(i16,2x))') dm2, file_nmb_bgn(ii), file_nmb_END(ii)
	  IF(dm2.ne.ii) then
	  print*, "Error in getting matrix index #", ii, dm2
	  stop
	  END IF
	  END do
	  close(1)
!--------------------------------------------------------------------
! Now each procs detrermines which files it should read, and which part of each file
!--------------------------------------------------------------------	  
	  proc_start = 1
	  IF(wrk_id.gt.0) then
	  do ii = 0, wrk_id-1
	  bk_dm = bk_non_zero_counter_old/mpi_task_per_traject							! determines the size of the chunk of the total matrix it needs to read
	  IF(ii.lt.mod(bk_non_zero_counter_old,
     & mpi_task_per_traject)) bk_dm = bk_dm + 1
	  proc_start = proc_start + bk_dm
	  END do
	  END IF
	  
	  bk_st = 0
	  do ii = 1, file_counter														! determines how many files contain that information
!	  IF(bk_st.ge.file_nmb_bgn(ii) .and. 
!     & bk_st.le.file_nmb_END(ii)) then
!      file_nmb1 = ii
!	  exit
!	  END IF
	  bk_st = bk_st + nmb_m(ii)														
	  IF(bk_st.ge.proc_start) then
      file_nmb1 = ii																! determine the number where the reading starts by each procs
!	  bk_dm1 = bk_st - (proc_start - 1)
	  exit
	  END IF
	  END do   

	  dm1 = 0
	  IF(wrk_id.gt.0)  then															
	  IF(file_nmb1.gt.1) then	  
	  bk_dm1 = 0
	  do ii = 1, file_nmb1-1
	  bk_dm1 = bk_dm1 + nmb_m(ii)
	  END do
	  dm1 = (proc_start - 1) - bk_dm1
	  ELSE
	  dm1 = proc_start - 1
	  END IF
	  END IF
	  
	  bk_st = 0
	  do ii = 1, file_counter
!	  IF(bk_fn.ge.file_nmb_bgn(ii) .and. 
!     & bk_fn.le.file_nmb_END(ii)) then
!      file_nmb2 = ii
!	  exit
!	  END IF
	  bk_st = bk_st + nmb_m(ii)
	  IF(bk_st.ge.proc_start+bk_nn-1) then
      file_nmb2 = ii																! determine the number where the reading ENDs by each procs
	  exit
	  END IF
	  END do
	  file_nmb = file_nmb2 - file_nmb1 + 1
	  st_store = 0   
	  
	  do jj = 1, file_nmb															! using start and END numbers make the file names to read
	  write(bk_dir_temp_1, '(a,a,i0,a,i0,a)') trim(bk_dir2),
     & '/MIJ_', file_nmb_bgn(jj - 1 + file_nmb1),'_',
     & file_nmb_END(jj - 1 + file_nmb1),'.DAT'
	  bk_dir_temp_2 = trim(bk_dir_temp_1)
	  write(bk_dir_temp_3, '(a,a,i0,a,i0,a)') trim(bk_dir2),
     & '/Ind_', file_nmb_bgn(jj - 1 + file_nmb1),'_',
     & file_nmb_END(jj - 1 + file_nmb1),'.DAT'
	  bk_dir_temp_4 = trim(bk_dir_temp_3)
!--------------------------------------------------------------------
! This is the actual reading of the matrix elements by individual processors
!--------------------------------------------------------------------	  
	  IF(jj.eq.1) then
	  
	  IF(.not.unformat_defined) then
	  OPEN(1,FILE=bk_dir_temp_2,ACTION="read")
	  OPEN(11,FILE=bk_dir_temp_4,ACTION="read")
	  IF(dm1.gt.0) then
	  do jjj = 1, dm1
      do ii=1,n_r_coll	  
	  read(1,'(i16,1x,i8,1x,e19.12)') dm2, dm3, dm4
!	  IF(myid.eq.2) write(*,'(i16,1x,i8,1x,e19.12)') dm2, dm3, dm4
!	  IF(dm2.eq.bk_st .and. dm3.eq.1) then
!	  write(*,'(a,a,2(2x,i0))') 
!     & "Something is wrong in Matrix reading in file ",
!     & bk_dir_temp_2, bk_st, bk_fn
!	  END IF
	  END do
	  read(11,*)
	  END do
	  END IF
	  
	  ELSE
	  OPEN(1,FILE=bk_dir_temp_2,ACTION="read",form="unformatted")
	  OPEN(11,FILE=bk_dir_temp_4,ACTION="read",form="unformatted")
	  IF(dm1.gt.0) then
	  do jjj = 1, dm1
      READ(1) dm2
      READ(1) dm9(:)
      READ(11) dm10(:)
	  END do
	  END IF
	  
	  END IF
	  
	  IF(.not.unformat_defined) then
	  do st = 1, min0(bk_nn, nmb_m(file_nmb1)-dm1)
      do ii=1,n_r_coll	  
      read(1,'(i16,1x,i8,1x,e19.12)')bk_nnn,iii,bk_mat_array(ii, st)
!      write(*,'(i2,1x,i16,1x,i8,1x,e19.12)')myid,bk_nnn,
!     & iii,bk_mat_array(ii, st)
	  END do
	  read(11,*) bk_ind_mat(1,st), bk_ind_mat(2,st)
!	  bk_ind_mat(:,st) = ind_mat(:,bk_nnn)
	  END do

	  ELSE
	  do st = 1, min0(bk_nn, nmb_m(file_nmb1)-dm1)
      read(1) bk_nnn
	  read(1) bk_mat_array(:, st)
	  read(11) bk_ind_mat(:,st)
!	  bk_ind_mat(:,st) = ind_mat(:,bk_nnn)
	  END do
	  
	  END IF
      CLOSE(1)
      CLOSE(11)
	  st_store = st
	  
	  ELSE
	  
	  IF(.not.unformat_defined) then
	  OPEN(1,FILE=bk_dir_temp_2,ACTION="read")
	  OPEN(11,FILE=bk_dir_temp_4,ACTION="read")
	  do st = st_store, min0(bk_nn, 
     & (st_store-1)+nmb_m(file_nmb1+jj-1))
      do ii=1,n_r_coll	  
      read(1,'(i16,1x,i8,1x,e19.12)')bk_nnn,iii,bk_mat_array(ii, st)
!      write(*,'(i2,1x,i16,1x,i8,1x,e19.12)')myid, bk_nnn,iii,
!     & bk_mat_array(ii, st)
	  END do
	  read(11,*) bk_ind_mat(1,st), bk_ind_mat(2,st)
	  END do
	  
	  ELSE
	  OPEN(1,FILE=bk_dir_temp_2,ACTION="read",form="unformatted")
	  OPEN(11,FILE=bk_dir_temp_4,ACTION="read",form="unformatted")
	  do st = st_store, min0(bk_nn, 
     & (st_store-1)+nmb_m(file_nmb1+jj-1))
      read(1) bk_nnn
	  read(1) bk_mat_array(:, st)
	  read(11) bk_ind_mat(:,st)
!	  bk_ind_mat(:,st) = ind_mat(:,bk_nnn)
	  END do
	  
	  END IF
      CLOSE(1)
      CLOSE(11)
	  st_store = st
	  END IF
	  END do
	  
	  do st = 1, bk_nn
	  IF(bikram_mij_shIFt) then
	  do ii=1,n_r_coll
	  bk_mat_array(ii, st) = bk_mat_array(ii, st) - 
     & bk_mat_array(n_r_coll, st)		!Bikram
	  ENDdo
	  ENDIF
	  END do
	  
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
	  db_time2 = mpi_wtime()												!Dulat 12/2/2024 - matrix reading time check
	  db_time = db_time2 - db_time1											!Dulat 12/2/2024 - matrix reading time check
	  IF(myid.eq.0) write(*,*)"Matrix Reading Finished", db_time
	  return
	  	  
      END SUBROUTINE

	   SUBROUTINE bk_read_matrix_info
! This subroutine is created by Bikramaditya Mandal, Nov 2020
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      USE OLD_MIJ
      IMPLICIT NONE
      INTEGER j_count,j_summ,p_count
      INTEGER st,i,k,istat
      INTEGER KRONEKER,p_lim_max_ini	  
      EXTERNAL KRONEKER
      LOGICAL file_exst
	  LOGICAL EVEN_NUM			! Dulat 9/25/2023									   
	  CHARACTER(LEN=500) :: bk_matrix_path1, bk_matrix_path2
	  CHARACTER(LEN=500) :: bk_dir_temp_1,bk_dir_temp_2
	  logical bk_exst
	  
	  write(bk_dir_temp_1, '(a,a)') trim(bk_dir2), "/MATRIX_INFO.DAT"
	  bk_dir_temp_2 = trim(bk_dir_temp_1)
	  
      INQUIRE( FILE=bk_dir_temp_2, EXIST=file_exst )
      IF(.not.file_exst) THEN 
      PRINT*, "ERROR: MATRIX_INFO.DAT NOT FOUND"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF
      IF(.not.identical_particles_defined) THEN	  
      p_lim_max = 1
      p_lim_min = 1	
      ELSE
      SELECT CASE(coll_type)		  
      CASE(5)	  

      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
      CASE(6)
      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))*
     & KRONEKER(v1_ch(chann_ini),v2_ch(chann_ini)) 	  
      CASE(0)
      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
     & *KRONEKER(ka1_ch(chann_ini),ka2_ch(chann_ini))
     & *KRONEKER(kc1_ch(chann_ini),kc2_ch(chann_ini))	  
      END SELECT		  
      ENDIF	  
	  
	  IF(.not.unformat_defined) then
      OPEN(1,FILE=bk_dir_temp_2,STATUS="OLD",ACTION="READ")
      PRINT*,"Matrix Info Reading Started"	  
      READ(1,'(a15,1x,i1)') buffer_word_2,coll_type_old
      IF(coll_type_old.ne.coll_type) THEN
      PRINT*, "ERROR: SYSTEM IS DIFFERENT"
      CRITICAL_ERROR = .TRUE.	
      ENDIF	 
      READ(1,'(a17,x,i4)') buffer_word_3,number_of_channels_old	 
      READ(1,'(a17,1x,i6)') buffer_word_3,states_size_old
      READ(1,'(a12,1x,i9)') buffer_word_1,total_size_old	  
      READ(1,'(a12,1x,i4)') buffer_word_1,n_r_coll_old
      IF(n_r_coll_old.ne.n_r_coll) THEN
      PRINT*, "ERROR:WRONG GRID"
      CRITICAL_ERROR = .TRUE.	  
      RETURN	  
      ENDIF	  
      ALLOCATE(j12_old(states_size_old),m12_old(states_size_old))
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined) THEN
      ALLOCATE(j12_h_old(states_size_old),m12_h_old(states_size_old))
      ENDIF	  
      ALLOCATE(indx_chann_old(states_size_old))
      IF(identical_particles_defined) 
     & ALLOCATE(parity_states_old(states_size_old))  
      SELECT CASE(coll_type)
      CASE(1)
      READ(1,*)
      ALLOCATE(j_ch_old(number_of_channels_old))
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined)
     & ALLOCATE(j_h_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      IF(fine_structure_defined) THEN	!!!! MODIFY
      IF(SPIN_FINE.ne.2) THEN	  
	  READ(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,f_old_b,n_old_b
      ELSE
	  READ(1,'(i8,1x,i8,1x,f4.1,1x,f5.1,1x,f4.1,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),
     & j12_h_old(st),m12_h_old(st),j_h_old_b,f_old_b,par_old_b	  
!!!! TO ADD	  
      ENDIF
      ELSE
	  READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined) THEN
      j_h_ch_old(indx_chann_old(st_old)) = j_h_old_b	
      ELSE	  
      j_ch_old(indx_chann_old(st_old)) = j_old_b  
      ENDIF
      ENDDO		  
      CASE(2)
      READ(1,*)
      ALLOCATE(v_ch_old(number_of_channels_old),
     & j_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),v_old_b,j_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      v_ch_old(indx_chann_old(st_old)) = v_old_b	  
      ENDDO
      CASE(3)
      READ(1,*)
      ALLOCATE(j_ch_old(number_of_channels_old),
     & k_ch_old(number_of_channels_old),
     & eps_ch_old(number_of_channels_old) )	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,k_old_b,eps_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      k_ch_old(indx_chann_old(st_old)) = k_old_b
      eps_ch_old(indx_chann_old(st_old)) = eps_old_b	  
      ENDDO	  
      CASE(4)
      READ(1,*)
      ALLOCATE(j_ch_old(number_of_channels_old),
     & ka_ch_old(number_of_channels_old),
     & kc_ch_old(number_of_channels_old) )	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,ka_old_b,kc_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      ka_ch_old(indx_chann_old(st_old)) = ka_old_b
      kc_ch_old(indx_chann_old(st_old)) = kc_old_b	  
      ENDDO	  	  
      CASE(5)
      READ(1,*)
      IF(.not.identical_particles_defined) THEN	  
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & j1_old_b,j2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b		 
      ENDDO
      ELSE
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & j1_old_b,j2_old_b,par_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Muj FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      parity_states_old(st_old) = par_old_b	  
      ENDDO	  
      ENDIF	  
      CASE(6)
      READ(1,*)
      IF(.not.identical_particles_defined) THEN	  
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))
      ALLOCATE(v1_ch_old(number_of_channels_old),
     & v2_ch_old(number_of_channels_old))	 
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & v1_old_b,v2_old_b,j1_old_b,j2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      v1_ch_old(indx_chann_old(st_old)) = v1_old_b
      v2_ch_old(indx_chann_old(st_old)) = v2_old_b		  
      ENDDO
      ELSE
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & v1_old_b,v2_old_b,j1_old_b,j2_old_b,par_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      v1_ch_old(indx_chann_old(st_old)) = v1_old_b
      v2_ch_old(indx_chann_old(st_old)) = v2_old_b		  
      parity_states_old(st_old) = par_old_b	  
      ENDDO	  
      ENDIF	  
      CASE(7)
      READ(1,*)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & k1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),	 
     & eps1_ch_old(number_of_channels_old) )	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j1_old_b,k1_old_b,eps1_old_b,j2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b	  
      k1_ch_old(indx_chann_old(st_old)) = k1_old_b
      eps1_ch_old(indx_chann_old(st_old)) = eps1_old_b
      ENDDO	  
      CASE(8)
      READ(1,*)	  
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old))
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j1_old_b,ka1_old_b,kc1_old_b,
     & j2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
	  
      ENDDO	 	 
      CASE(9)
      READ(1,*)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old)
     & ,k2_ch_old(number_of_channels_old),
     & eps2_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old
      READ(1,
     & '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,k2_old_b,eps2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      k2_ch_old(indx_chann_old(st_old)) = k2_old_b
      eps2_ch_old(indx_chann_old(st_old)) = eps2_old_b
      ENDDO	  
      CASE(0)
      READ(1,*)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old)
     & ,ka2_ch_old(number_of_channels_old),
     & kc2_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old
      IF(.NOT.identical_particles_defined) THEN	  
      READ(1,
     & '(2(i8,1x),8(i8,1x))')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,ka2_old_b,kc2_old_b
      ELSE
      READ(1,
     & '(i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i8,1x,i2)
     & ')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,ka2_old_b,kc2_old_b,par_old_b	  
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	 
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      ka2_ch_old(indx_chann_old(st_old)) = ka2_old_b
      kc2_ch_old(indx_chann_old(st_old)) = kc2_old_b
      IF(identical_particles_defined)
     & parity_states_old(st_old) = par_old_b	  
      ENDDO	 
      END SELECT
	  
      st = 0    
!!!! TO BE CONTINUED	  
      DO i=1,min(number_of_channels_old,number_of_channels)
      SELECT CASE(coll_type)
      CASE(1)
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) THEN
      DO j_count = 1, int(j_h_ch_old(i)*2) + 1	  
      st = st + 1
      IF(j_count*2 .eq. int(j_h_ch_old(i)*2)+1
     & .and. chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF	 
      ENDDO	  
      ELSE	  
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	  
      ENDDO
      ENDIF	  
	  
      CASE(2)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(3)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(4)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
      ENDDO
      ENDDO
      ELSE
      p_lim_max= 2-KRONEKER(j1_ch(i),j2_ch(i)) !min(p_lim_max_ini,2-KRONEKER(j1_ch(i),j2_ch(i)))	  
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch(i),j2_ch(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and. p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	
      ENDDO
      ENDDO
      ELSE
      p_lim_max= 2-KRONEKER(j1_ch(i),j2_ch(i))*
     & KRONEKER(v1_ch(i),v2_ch(i))      
!	 p_lim_max=min(p_lim_max_ini,2-KRONEKER(j1_ch(i),j2_ch(i))*
!     & KRONEKER(v1_ch(i),v2_ch(i)))    
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch(i),j2_ch(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and. p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO 
  

	  
      ENDIF	  
      CASE(7)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(8)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(9)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(0)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      ELSE
      p_lim_max= 2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i))	      
!	 p_lim_max=min(p_lim_max_ini,2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
!     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
!     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i)))	
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
!     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
!     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     &	  abs(j1_ch_old(i)-j2_ch_old(i))
     & .and.  p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	   	  
      END SELECT	  
      ENDDO
	  
! HERE WE CHECK THE ORDER !! FOR FINE SPIN=2 WE JUST OVERJUM
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined) GOTO 1994
      DO st=1,min(states_size_old,states_size)
      IF(j12(st).ne.j12_old(st)) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: J12s ARE DIFFERENT"
      RETURN	  
      ENDIF	  
      IF(m12(st).ne.m12_old(st)) THEN
      CRITICAL_ERROR = .TRUE.	
      PRINT*, "ERROR IN INI: CHECK M12"
      RETURN	  
      ENDIF	 
      SELECT CASE(coll_type)
      CASE(1)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*,"ERROR IN INI : CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      CASE(2)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      IF(v_ch_old(indx_chann_old(st)).ne.v_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI : CHECK V"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	    
      CASE(3)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      IF(k_ch_old(indx_chann_old(st)).ne.k_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK K"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
      IF(eps_ch_old(indx_chann_old(st)).ne.eps_ch(indx_chann(st)))
     & THEN
      PRINT*, "ERROR IN INI: CHECK EPS"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      CASE(4)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka_ch_old(indx_chann_old(st)).ne.ka_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc_ch_old(indx_chann_old(st)).ne.kc_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(5)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(6)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(v1_ch_old(indx_chann_old(st)).ne.v1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK V1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(v2_ch_old(indx_chann_old(st)).ne.v2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK V2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(7)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(k1_ch_old(indx_chann_old(st)).ne.k1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK K1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(eps1_ch_old(indx_chann_old(st)).ne.eps1_ch(indx_chann(st)))
     & THEN
      PRINT*, "ERROR IN INI: CHECK EPS1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(8)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(9)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(k2_ch_old(indx_chann_old(st)).ne.k2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK K2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(eps2_ch_old(indx_chann_old(st)).ne.eps2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK EPS2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(0)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka2_ch_old(indx_chann_old(st)).ne.ka2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc2_ch_old(indx_chann_old(st)).ne.kc2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(identical_particles_defined) THEN	 
      IF(parity_states_old(st)
     & .ne.parity_state(st))
     & STOP "ERROR IN INI: CHECK PARITY"
      ENDIF	 
      END SELECT	  	  	  
      ENDDO	 
! HERE WE CHECK THE ORDER	  
!      PRINT*, "WE HAVE CKECED"	  
1994      READ(1,*)
!      PRINT*, buffer_word_3
      DO i=1,total_size_old
      READ(1,"(i16,1x,i8,1x,i8)") i_old,ind_mat_old_1,ind_mat_old_2
      IF(i.ne.i_old) CRITICAL_ERROR = .TRUE.
!      IF(i.le.total_size) THEN	  
!      IF(ind_mat_old_1.ne.ind_mat(1,i)) CRITICAL_ERROR = .TRUE.
!      IF(ind_mat_old_2.ne.ind_mat(2,i))	CRITICAL_ERROR = .TRUE.
!      ENDIF	  
      ENDDO		  
      READ(1,*)	  
      DO i=1,n_r_coll
      READ(1,'(i5,1x,e19.12)')i_old,R_COM(i)	  
      ENDDO
	  
	  ELSE
      OPEN(1,FILE=bk_dir_temp_2,STATUS="OLD",ACTION="READ",
     & form="unformatted")
      PRINT*,"Matrix_unformatted Info Reading Started"	  
      READ(1,IOSTAT=istat) coll_type_old
      IF(coll_type_old.ne.coll_type) THEN 
      PRINT*, "ERROR: SYSTEM IS DIFFERENT"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
      READ(1,IOSTAT=istat) number_of_channels_old	 
      READ(1,IOSTAT=istat) states_size_old
      READ(1,IOSTAT=istat) total_size_old	  
      READ(1,IOSTAT=istat)n_r_coll_old
      IF(n_r_coll_old.ne.n_r_coll) THEN
      n_r_coll_old = n_r_coll
!      CRITICAL_ERROR = .TRUE.
      PRINT*,"WRONG GRID", n_r_coll_old, n_r_coll
!      RETURN	  
      ENDIF	  
      IF(istat.ne.0) THEN
      CRITICAL_ERROR = .TRUE.
      PRINT*,"ERROR in ",trim(bk_dir_temp_2)
      RETURN	  
      ENDIF	  
      ALLOCATE(j12_old(states_size_old),m12_old(states_size_old))
      ALLOCATE(indx_chann_old(states_size_old))
      IF(identical_particles_defined) 
     & ALLOCATE(parity_states_old(states_size_old))  
      SELECT CASE(coll_type)
      CASE(1)
      ALLOCATE(j_ch_old(number_of_channels_old))     
      DO st=1,states_size_old
      READ(1)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	  
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      ENDDO		  
      CASE(2)
      ALLOCATE(v_ch_old(number_of_channels_old),
     & j_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),v_old_b,j_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      v_ch_old(indx_chann_old(st_old)) = v_old_b	  
      ENDDO
      CASE(3)
      ALLOCATE(k_ch_old(number_of_channels_old),
     & j_ch_old(number_of_channels_old),
     & eps_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,k_old_b,eps_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      k_ch_old(indx_chann_old(st_old)) = k_old_b
      eps_ch_old(indx_chann_old(st_old)) = eps_old_b
      ENDDO	  
      CASE(4)	  
      ALLOCATE(ka_ch_old(number_of_channels_old),
     & j_ch_old(number_of_channels_old),
     & kc_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,ka_old_b,kc_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      ka_ch_old(indx_chann_old(st_old)) = ka_old_b
      kc_ch_old(indx_chann_old(st_old)) = kc_old_b
      ENDDO	  
      CASE(5)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))
      DO st=1,states_size_old
      IF(.not.identical_particles_defined) THEN	 
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & j1_old_b,j2_old_b
      ELSE
      READ(1,IOSTAT=istat)	  
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & j1_old_b,j2_old_b,par_old_b
      parity_states_old(st_old) = par_old_b	 
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b		 
      ENDDO
      CASE(6)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))
      ALLOCATE(v1_ch_old(number_of_channels_old),
     & v2_ch_old(number_of_channels_old))	 
      DO st=1,states_size_old
      IF(.not.identical_particles_defined) THEN	 
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & v1_old_b,j1_old_b,v2_old_b,j2_old_b
      ELSE
      READ(1,IOSTAT=istat)	  
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & v1_old_b,j1_old_b,v2_old_b,j2_old_b,par_old_b
      parity_states_old(st_old) = par_old_b	 
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      v1_ch_old(indx_chann_old(st_old)) = v1_old_b
      v2_ch_old(indx_chann_old(st_old)) = v2_old_b	  
      ENDDO
      CASE(7)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & k1_ch_old(number_of_channels_old),
     & eps1_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old

      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,k1_old_b,eps1_old_b
     & ,j2_old_b
	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      k1_ch_old(indx_chann_old(st_old)) = k1_old_b
      eps1_ch_old(indx_chann_old(st_old)) = eps1_old_b
      ENDDO		  
      CASE(8)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old

      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b
	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      ENDDO	 	  	  
      CASE(9)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old)
     & ,k2_ch_old(number_of_channels_old),
     & eps2_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old

      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,k2_old_b,eps2_old_b
	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      k2_ch_old(indx_chann_old(st_old)) = k2_old_b
      eps2_ch_old(indx_chann_old(st_old)) = eps2_old_b
      ENDDO	 	  
      CASE(0)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old)
     & ,ka2_ch_old(number_of_channels_old),
     & kc2_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old
      IF(.NOT.identical_particles_defined) THEN	  
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,ka2_old_b,kc2_old_b
      ELSE
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,ka2_old_b,kc2_old_b,par_old_b	  
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      ka2_ch_old(indx_chann_old(st_old)) = ka2_old_b
      kc2_ch_old(indx_chann_old(st_old)) = kc2_old_b
      IF(identical_particles_defined)
     & parity_states_old(st_old) = par_old_b	  
      ENDDO	 
      END SELECT
      st = 0       
      DO i=1,min(number_of_channels_old,number_of_channels)
      SELECT CASE(coll_type)
      CASE(1)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	  
      ENDDO
      CASE(2)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(3)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(4)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
      ENDDO
      ENDDO
      ELSE
      p_lim_max= 2-KRONEKER(j1_ch(i),j2_ch(i)) !min(p_lim_max_ini,2-KRONEKER(j1_ch(i),j2_ch(i)))	  
	  
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch(i),j2_ch(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and. p_count .eq.1  ) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	
      ENDDO
      ENDDO
      ELSE
      p_lim_max= 2-KRONEKER(j1_ch(i),j2_ch(i))*
     & KRONEKER(v1_ch(i),v2_ch(i))      
!	 p_lim_max=min(p_lim_max_ini,2-KRONEKER(j1_ch(i),j2_ch(i))*
!     & KRONEKER(v1_ch(i),v2_ch(i)))	  
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch(i),j2_ch(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and. p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO 
  

	  
      ENDIF	  
      CASE(7)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(8)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(9)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(0)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      ELSE
      p_lim_max= 2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i))     
!	 p_lim_max=min(p_lim_max_ini,2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
!     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
!     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i)))	
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
!     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
!     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     &	  abs(j1_ch_old(i)-j2_ch_old(i))
     & .and.  p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	   	  
      END SELECT	  
      ENDDO
! HERE WE CHECK THE ORDER 
      DO st=1,min(states_size_old,states_size)
      IF(j12(st).ne.j12_old(st)) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: J12s ARE DIFFERENT"
      RETURN	  
      ENDIF	  
      IF(m12(st).ne.m12_old(st)) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR IN INI: CHECK M12"
      RETURN	  
      ENDIF	  
      SELECT CASE(coll_type)
      CASE(1)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*,"ERROR IN INI : CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      CASE(2)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      IF(v_ch_old(indx_chann_old(st)).ne.v_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI : CHECK V"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	    
      CASE(3)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      IF(k_ch_old(indx_chann_old(st)).ne.k_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK K"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
      IF(eps_ch_old(indx_chann_old(st)).ne.eps_ch(indx_chann(st)))
     & THEN
      PRINT*, "ERROR IN INI: CHECK EPS"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      CASE(4)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka_ch_old(indx_chann_old(st)).ne.ka_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc_ch_old(indx_chann_old(st)).ne.kc_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(5)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(6)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(v1_ch_old(indx_chann_old(st)).ne.v1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK V1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(v2_ch_old(indx_chann_old(st)).ne.v2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK V2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(7)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(k1_ch_old(indx_chann_old(st)).ne.k1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK K1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(eps1_ch_old(indx_chann_old(st)).ne.eps1_ch(indx_chann(st)))
     & THEN
      PRINT*, "ERROR IN INI: CHECK EPS1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(8)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(9)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(k2_ch_old(indx_chann_old(st)).ne.k2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK K2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(eps2_ch_old(indx_chann_old(st)).ne.eps2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK EPS2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(0)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka2_ch_old(indx_chann_old(st)).ne.ka2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc2_ch_old(indx_chann_old(st)).ne.kc2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(identical_particles_defined) THEN	 
      IF(parity_states_old(st)
     & .ne.parity_state(st))
     & STOP "ERROR IN INI: CHECK PARITY"
      ENDIF	 
      END SELECT	  
	  
	  
      ENDDO	  

! HERE WE CHECK THE ORDER	  
!	  DO i=1,total_size_old
!      READ(1) i_old,ind_mat_old_1,ind_mat_old_2
!      IF(i.ne.i_old) STOP "ERROR IN INI: CHECK FILE"
!      IF(i.le.total_size) THEN	  
!      IF(ind_mat_old_1.ne.ind_mat(1,i)) THEN
!      PRINT*, "ERROR IN INI:IND_MAT_1"
!      CRITICAL_ERROR = .TRUE.
!      RETURN	  
!      ENDIF	  
!      IF(ind_mat_old_2.ne.ind_mat(2,i))	THEN
!      PRINT*,"ERROR IN INI:IND_MAT_2"
!      CRITICAL_ERROR = .TRUE.
!      RETURN	  
!      ENDIF	  
!      ENDIF	  
!      ENDDO		  
      DO i=1,n_r_coll
      READ(1)i_old,R_COM(i)
      ENDDO

	  END IF

      END SUBROUTINE

	   SUBROUTINE bk_matrix_combination(mat_size, tot_proc)
	  use MPI
	  use VARIABLES
! This subroutine is created by Bikramaditya Mandal, Nov 2020
      IMPLICIT NONE
	  INTEGER mat_size, tot_proc, chunk, k_st, k_fn
	  INTEGER ierr_mpi, mij_remainder, i
	  CHARACTER (LEN=100) :: bk_dir_temp_1,bk_dir_temp_2
	  CHARACTER (LEN=100) :: bk_dir_temp_3,bk_dir_temp_4
	  CHARACTER (LEN=255) :: bk_dir_temp_5
	  logical file_exst
	  
	  write(bk_dir_temp_3, '(a,a)') trim(bk_dir2), 
     & "/MATRIX_FILE_INDEX.DAT"
	  bk_dir_temp_4 = trim(bk_dir_temp_3)
	  open(11,file = trim(bk_dir_temp_4))
	  write(11,'(a15,i10)')"Total #files = ", tot_proc
	  write(11,'(3(a16,2x))') '            #Mij', '      #Mij_begin'
     & ,'        #Mij_END'
	  
	  mij_remainder = 0
	  do i = 1, tot_proc
	  chunk = mat_size/tot_proc
      k_st = (i-1)*chunk + 1
      k_fn = i*chunk 	  
	  IF(mat_size.gt.(chunk*tot_proc)) 
     & mij_remainder = mat_size - chunk*tot_proc
	  
	  IF((i-1).lt.mij_remainder) then	  
      k_st = k_st + i - 1
      k_fn = k_fn + i
	  ELSE	  
      k_st = k_st + mij_remainder
      k_fn = k_fn + mij_remainder
	  END IF
	  
	  write(11,'(3(i16,2x))') i, k_st, k_fn
	  write(bk_dir_temp_1, '(a,a,i0,a,i0,a)') 
     & trim(bk_dir2), '/MIJ_',k_st,'_',k_fn,'.DAT'
	  bk_dir_temp_2 = trim(bk_dir_temp_1)
	  INQUIRE(FILE=bk_dir_temp_2, EXIST=file_exst )
      IF(.not.file_exst) THEN 
      PRINT*, "ERROR: MATRIX FILES NOT FOUND ", trim(bk_dir_temp_2)
	  END IF	  
	  
	  END do
	  
	  close(11)
	  
!	  write(bk_dir_temp_3, '(a,a)') trim(bk_dir2), 
!     & "/MATRIX_COMBINE.sh"
!	  bk_dir_temp_4 = trim(bk_dir_temp_3)
!	  open(11,file = bk_dir_temp_4)
!	  write(11,*)'#!/bin/sh'
!	  write(11,*)'cat MATRIX_INFO.DAT >> ../MTRX.DAT'
	  
	  mij_remainder = 0
	  do i = 1, tot_proc
	  chunk = mat_size/tot_proc
      k_st = (i-1)*chunk + 1
      k_fn = i*chunk 	  
	  IF(mat_size.gt.(chunk*tot_proc)) 
     & mij_remainder = mat_size - chunk*tot_proc
	  
	  IF((i-1).lt.mij_remainder) then	  
      k_st = k_st + i - 1
      k_fn = k_fn + i
	  ELSE	  
      k_st = k_st + mij_remainder
      k_fn = k_fn + mij_remainder
	  END IF
	  
	  write(bk_dir_temp_1, '(a,a,i0,a,i0,a)') 
     & trim(bk_dir2), '/MIJ_',k_st,'_',k_fn,'.DAT'
	  bk_dir_temp_2 = trim(bk_dir_temp_1)
	  INQUIRE(FILE=bk_dir_temp_2, EXIST=file_exst )
      IF(.not.file_exst) THEN 
      PRINT*, "ERROR: MATRIX FILES NOT FOUND ", trim(bk_dir_temp_2)
	  return
	  stop
	  call MPI_FINALIZE(ierr_mpi)
	  END IF	  
	  
	  write(bk_dir_temp_1, '(a,i0,a,i0,a)') 
     & 'cat MIJ_',k_st,'_',k_fn,'.DAT >> ../MTRX.DAT'
	  bk_dir_temp_2 = trim(bk_dir_temp_1)
!	  write(11,*) bk_dir_temp_2
	  
	  END do
	  
!	  close(11)
	  
      END SUBROUTINE	
	  
      subroutine bikram_kappa(bk_ka1, bk_kc1, bk_ka2, bk_kc2, 
     & kappa1, kappa2)
! This subroutine is written by Bikramaditya Mandal
! This determine the sign of kappa1 and kappa2 which is used to compute total parity.
      implicit none
      INTEGER bk_ka1, bk_kc1, bk_ka2, bk_kc2
	  INTEGER symm_coeff1, symm_coeff2
      INTEGER k12, kappa1, kappa2

      symm_coeff1 = bk_ka1 - bk_kc1
	  IF(mod(abs(symm_coeff1),2).eq.0d0 .and. 
     & mod(abs(symm_coeff1-1),2).eq.1d0) then
	  kappa1 = 0d0
	  ELSE IF(mod(abs(symm_coeff1),2).eq.1d0 .and. 
     & mod(abs(symm_coeff1-1),2).eq.0d0) then
	  kappa1 = 1d0
	  ELSE
	  write(*,'(a)')'Something is wrong in Computation of Kappa1'
	  stop
	  END IF
	  
      symm_coeff2 = bk_ka2 - bk_kc2	 
	  IF(mod(abs(symm_coeff2),2).eq.0d0 .and. 
     & mod(abs(symm_coeff2-1),2).eq.1d0) then
	  kappa2 = 0d0
	  ELSE IF(mod(abs(symm_coeff2),2).eq.1d0 .and. 
     & mod(abs(symm_coeff2-1),2).eq.0d0) then
	  kappa2 = 1d0
	  ELSE
	  write(*,'(a)')'Something is wrong in Computation of Kappa2'
	  stop
	  END IF
!	  k12 = kappa1 + kappa2
	  
      END subroutine

