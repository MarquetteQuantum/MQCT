      SUBROUTINE WF_ASYM_TOP_ASYM_TOP(wfr,wfi,i_state,
     & alpha1,beta1,gamma1,alpha2,beta2,gamma2)
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      IMPLICIT NONE
      REAL*8 wfr,wfi,dr1,dr2,di1,di2,clb_grd
      INTEGER i_state,chann,j1_t,ka1_t,kc1_t,m1_t,m12_t,j12_t
     &,j2_t,ka2_t,kc2_t,m2_t
      REAL*8 alpha1,beta1,gamma1,alpha2,beta2,gamma2
      REAL*8 CG	  
      EXTERNAL  CG
!!! WAVEFUNCTIONS OF ASYM_TOP + ASYM_TOP	  
      wfr = 0d0
      wfi = 0d0	  
      chann = indx_chann(i_state) !!! INI OF CHANNEL
!!!  INI OF QUNATUM NUMBERS
      j12_t = j12(i_state)	!! TOTAL j12  
      m12_t = m12(i_state)  !! TOTAL M12		  
      j1_t = j1_ch(chann)   !! TOTAL J1
      j2_t = j2_ch(chann)   !! TOTAL J2
      ka1_t = ka1_ch(chann) !! KA1
      ka2_t = ka2_ch(chann) !! KA2
      kc1_t = kc1_ch(chann) !! KC1
      kc2_t = kc2_ch(chann) !! KC2
       IF(max(j12_t,j1_t,j2_t)*2 .gt.j12_t+j1_t+j2_t )
     &	  STOP"ERROR:TRINAGULAR RULE"
      DO m1_t = -j1_t,j1_t
      m2_t = m12_t - m1_t
      IF(abs(m2_t).le.j2_t) THEN
      CALL WF_ASYM_TOP(A1,B1,C1,j1_t,1,chann,m1_t,dr1,di1,
     & alpha1,beta1,gamma1)!!! WAVEFUNCTION OF THE FIRST MOLECULE
      CALL WF_ASYM_TOP(A2,B2,C2,j2_t,2,chann,m2_t,dr2,di2,
     & alpha2,beta2,gamma2) !!! WAVEFUNCTION OF THE SECOND MOLECULE
      clb_grd	= CG(j1_t,j2_t,j12_t,m1_t,m2_t,m12_t)
       wfr =wfr +  clb_grd*(dr1*dr2 - 
     &	di1*di2)
      wfi =wfi +  clb_grd*(dr1*di2 + 
     &	di1*dr2)	 
      ENDIF	  
      ENDDO
      END SUBROUTINE WF_ASYM_TOP_ASYM_TOP
      SUBROUTINE WF_ASYM_TOP_SYM_TOP(wfr,wfi,i_state,
     & alpha1,beta1,gamma1,alpha2,beta2,gamma2)
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      IMPLICIT NONE
      REAL*8 wfr,wfi,dr1,dr2,di1,di2,clb_grd
      INTEGER i_state,chann,j1_t,ka1_t,kc1_t,m1_t,m12_t,j12_t
     &,j2_t,k2_t,eps2_t,m2_t
      REAL*8 alpha1,beta1,gamma1,alpha2,beta2,gamma2
      REAL*8 CG	  
      EXTERNAL CG
      wfr = 0d0
      wfi = 0d0
! !!! WAVE FUNCTION OF ASYM + SYM TOP	  
      chann = indx_chann(i_state) !!!! INI CHANNEL
      j12_t = j12(i_state)	!!!   
      m12_t = m12(i_state)		  
      j1_t = j1_ch(chann)
      j2_t = j2_ch(chann)
      ka1_t = ka1_ch(chann)
      k2_t = k2_ch(chann)
      kc1_t = kc1_ch(chann)
      eps2_t = eps2_ch(chann)
	  	  
      IF(max(j12_t,j1_t,j2_t)*2 .gt.j12_t+j1_t+j2_t )
     &	  STOP"ERROR:TRINAGULAR RULE"
      DO m1_t = -j1_t,j1_t
      m2_t = m12_t - m1_t
      IF(abs(m2_t).le.j2_t) THEN
      CALL WF_ASYM_TOP(A1,B1,C1,j1_t,1,chann,m1_t,dr1,di1,
     & alpha1,beta1,gamma1)
      CALL WF_SYM_TOP(j2_t,k2_t,eps2_t,m2_t,dr2,di2,
     & alpha2,beta2,gamma2)
      clb_grd	= CG(j1_t,j2_t,j12_t,m1_t,m2_t,m12_t)
       wfr =wfr +  clb_grd*(dr1*dr2 - 
     &	di1*di2)
      wfi =wfi +  clb_grd*(dr1*di2 + 
     &	di1*dr2)	 
      ENDIF	  
      ENDDO	  
      END SUBROUTINE WF_ASYM_TOP_SYM_TOP	  
      SUBROUTINE WF_ASYM_TOP(A,B,C,J,mol_ind,
     & chnn_nmb,M,wfr,wfi,alpha,beta,gamma)
! This subroutine is updated by Bikramaditya Mandal
      USE EIGEN_VECTORS 	 
      IMPLICIT NONE 
      REAL*8 A,B,C,wfr,wfi,alpha,beta,gamma 
      INTEGER J,KA,KC,i,M,DK,chnn_nmb,mol_ind
      REAL*8 coeff(2*J+1)
      REAL*8 dwr,dwi,pi
      CALL ASYM_TOP_VECTORS	  
      pi = dacos(-1d0)
      coeff = 0d0
      IF(mol_ind.eq.1) THEN	  
      IF(A.eq.A1_I .and. B.eq.B1_I .and. C.eq.C1_I)
     & coeff = M1_VECTORS(chnn_nmb,1:2*J+1)
      IF(A.eq.A_I .and. B.eq.B_I .and. C.eq.C_I) then
      coeff = M_VECTORS(chnn_nmb,1:2*J+1)
	  end if
      ELSE	 
      IF(A.eq.A2_I .and. B.eq.B2_I .and. C.eq.C2_I)
     & coeff = M2_VECTORS(chnn_nmb,1:2*J+1)
      ENDIF 
      wfr = 0d0
      wfi = 0d0
      DO i=1,2*J+1
      DK = i-J-1
      CALL Djkm(dwr,dwi,J,DK,M,alpha,beta,gamma)
      dwr = dwr*dsqrt((2d0*dble(J)+1d0)/8d0/pi**2)
      dwi = dwi*dsqrt((2d0*dble(J)+1d0)/8d0/pi**2)	  
      wfr = wfr + dwr*coeff(i)
      wfi = wfi + dwi*coeff(i)	  
      ENDDO	  
!      IF(wfr.eq.0d0 .and. wfi.eq.0d0) STOP "ERROR:ASYM TOP MATRIX"	  			!Disabled because the wavefunction can be Zero at nodes
      END SUBROUTINE WF_ASYM_TOP 
	  
      SUBROUTINE WF_SYM_TOP(J,K,P,M,wfr,wfi,alpha,beta,gamma)
! This subroutine is updated by Bikramaditya Mandal
      IMPLICIT NONE
      REAL*8 wfr,wfi,alpha,beta,gamma 
      INTEGER J,K,P,i,M
      REAL*8 coeff(2*J+1)
      REAL*8 dwr,dwi,pi!, dmmr1, dmmr2, dmmi1, dmmi2
      pi = dacos(-1d0)
      wfr = 0d0
      wfi = 0d0
      CALL Djkm(dwr,dwi,J,K,M,alpha,beta,gamma)
      wfr = wfr + dwr
      wfi = wfi + dwi
!	  dmmr1 = dwr
!	  dmmi1 = dwi
      CALL Djkm(dwr,dwi,J,-K,M,alpha,beta,gamma)
      wfr = wfr + dwr*(1d0-2*P)
      wfi = wfi + dwi*(1d0-2*P)
!	  dmmr2 = dwr
!	  dmmi2 = dwi
      IF(K.eq.0) THEN
      wfr = wfr/2d0*dsqrt((2d0*dble(J)+1d0)/8d0/pi**2)  
      wfi = wfi/2d0*dsqrt((2d0*dble(J)+1d0)/8d0/pi**2) 
      IF(P.ne.0) then
	  print*, "ERROR: WRONG LEVEL DEFINITION FOR SYMMETRIC TOP", K, P
	  stop
	  end if
      ELSE	  
      wfr = wfr/dsqrt(2d0)*dsqrt((2d0*dble(J)+1d0)/8d0/pi**2) 	  
      wfi = wfi/dsqrt(2d0)*dsqrt((2d0*dble(J)+1d0)/8d0/pi**2) 
      ENDIF	
!	  write(*,'(4(i4,2x),9(f12.5,2x))')J, K, P, M, alpha, gamma, beta, 
!     & dmmr1, dmmr2, wfr, dmmi1, dmmi2, wfi
      END SUBROUTINE WF_SYM_TOP 
	  
      SUBROUTINE WF_DIAT_TOP(J,M,theta,phi,wfr,wfi)
! This subroutine is updated by Bikramaditya Mandal
      USE FACTORIAL	  
      IMPLICIT NONE
      INTEGER J,M
      REAL*8 theta,phi,wfr,wfi,plgndr,pi
      EXTERNAL plgndr
      pi = dacos(-1d0)
	  
      wfr = plgndr(J,abs(M),dcos(theta))
     & *dsqrt((2d0*J+1d0)/4d0/pi*
     & EXP(LOGFACT(J-abs(M)) - LOGFACT(J+abs(M))))*
     & dcos(M*phi)!*(-1)**M
      wfi = plgndr(J,abs(M),dcos(theta))
     & *dsqrt((2d0*J+1d0)/4d0/pi*
     & EXP(LOGFACT(J-abs(M)) - LOGFACT(J+abs(M))))*
     & dsin(M*phi)!*(-1)**M
!      RETURN
      IF(M<0) THEN
      wfr = wfr*(-1)**M
      wfi = wfi*(-1)**M	  
      ENDIF	  
      END SUBROUTINE WF_DIAT_TOP
	  
      SUBROUTINE WF_DIAT_DIAT(J12,M12,J1,J2,theta1,theta2,phi,wfr,wfi)
! This subroutine is updated by Bikramaditya Mandal
      USE FACTORIAL	  
      IMPLICIT NONE
      INTEGER J12,M12,J1,J2,M1,M2
      REAL*8 wfr1,wfr2,wfi1,wfi2
      REAL*8 theta1,theta2,phi,wfr,wfi,plgndr,pi,CG,clb_grd,phi1,phi2
      EXTERNAL  plgndr,CG
      pi = dacos(-1d0)
      IF(max(j12,j1,j2)*2 .gt.j12+j1+j2)
     &	  STOP"ERROR:TRINAGULAR RULE"
      phi1 = 0d0
      phi2 = phi
      wfr = 0d0
      wfi = 0d0	  
      DO M1=-J1,J1
      M2 = M12 - M1
      IF(ABS(M2).le.J2) THEN	  
      CALL WF_DIAT_TOP(J1,M1,theta1,phi1,wfr1,wfi1)
      CALL WF_DIAT_TOP(J2,M2,theta2,phi2,wfr2,wfi2)	  
	 
      clb_grd	= CG(j1,j2,j12,m1,m2,m12)	 
      wfr = wfr +  clb_grd*(wfr1*wfr2 - wfi1*wfi2)
      wfi = wfi +  clb_grd*(wfr1*wfi2 + wfi1*wfr2)	 
      ENDIF	 
      ENDDO  	 
      END SUBROUTINE WF_DIAT_DIAT	  
	  
      SUBROUTINE WF_DIAT_DIAT_IDENT(
     & J12,M12,J1,J2,theta1,theta2,phi,wfr,wfi,parity_ident)
! This subroutine is updated by Bikramaditya Mandal
      USE FACTORIAL	  
      IMPLICIT NONE
      INTEGER J12,M12,J1,J2,M1,M2,parity_ident
      REAL*8 wfr1,wfr2,wfi1,wfi2,delta
      REAL*8 theta1,theta2,phi,wfr,wfi,plgndr,pi,CG,clb_grd,phi1,phi2
      EXTERNAL  plgndr,CG,delta
      pi = dacos(-1d0)
      IF(max(j12,j1,j2)*2 .gt.j12+j1+j2)
     &	  STOP"ERROR:TRINAGULAR RULE"
      phi2 = 0d0
      phi1 = phi
      wfr = 0d0
      wfi = 0d0	  
      CALL WF_DIAT_DIAT(J12,M12,J1,J2,theta1,theta2,phi,wfr1,wfi1)
      CALL WF_DIAT_DIAT(J12,M12,J2,J1,theta1,theta2,phi,wfr2,wfi2)
      wfr = wfr1+wfr2*parity_ident*(-1)**J12
      wfi = wfi1+wfi2*parity_ident*(-1)**J12      	  
      wfr = wfr/dsqrt(2d0*(1d0+delta(J1,J2)))
      wfi = wfi/dsqrt(2d0*(1d0+delta(J1,J2)))


	  
      END SUBROUTINE WF_DIAT_DIAT_IDENT
	  
      SUBROUTINE WF_ASYM_TOP_ASYM_TOP_IDENT(wfr,wfi,i_state,
     & alpha1,beta1,gamma1,alpha2,beta2,gamma2)
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      IMPLICIT NONE
      REAL*8 wfr,wfi,dr1,dr2,di1,di2,clb_grd,delta
      INTEGER i_state,chann,j1_t,ka1_t,kc1_t,m1_t,m12_t,j12_t
     &,j2_t,ka2_t,kc2_t,m2_t,parity_ident
      REAL*8 alpha1,beta1,gamma1,alpha2,beta2,gamma2
      REAL*8 CG,wfr1,wfi1,wfr2,wfi2	  
      EXTERNAL  CG,delta
      INTEGER p_inv	  
      parity_ident = parity_state(i_state)
      wfr = 0d0
      wfi = 0d0
      wfr1 = 0d0
      wfi1 = 0d0
      wfr2 = 0d0
      wfi2 = 0d0		  
      CALL WF_ASYM_TOP_ASYM_TOP(wfr1,wfi1,i_state,
     & alpha1,beta1,gamma1,alpha2,beta2,gamma2)
  
      chann = indx_chann(i_state)
      p_inv = parity_inversion(chann)
      j12_t = j12(i_state)	  
      m12_t = m12(i_state)		  
      j1_t = j2_ch(chann)
      j2_t = j1_ch(chann)
       IF(max(j12_t,j1_t,j2_t)*2 .gt.j12_t+j1_t+j2_t )
     &	  STOP"ERROR:TRINAGULAR RULE"
      DO m1_t = -j1_t,j1_t
      m2_t = m12_t - m1_t
      IF(abs(m2_t).le.j2_t) THEN
      CALL WF_ASYM_TOP(A1,B1,C1,j1_t,2,chann,m1_t,dr1,di1,
     & alpha1,beta1,gamma1)
      CALL WF_ASYM_TOP(A2,B2,C2,j2_t,1,chann,m2_t,dr2,di2,
     & alpha2,beta2,gamma2)
      clb_grd	= CG(j1_t,j2_t,j12_t,m1_t,m2_t,m12_t)
      wfr2 =wfr2 +  clb_grd*(dr1*dr2 - 
     &	di1*di2)
      wfi2 =wfi2 +  clb_grd*(dr1*di2 + 
     &	di1*dr2)	 
      ENDIF	  
      ENDDO
      wfr = (wfr1+wfr2*parity_ident*(-1)**j12_t*p_inv)/
     & dsqrt(2d0*(1d0+delta(j1_t,j2_t)))
      wfi = (wfi1+wfi2*parity_ident*(-1)**j12_t*p_inv)/
     & dsqrt(2d0*(1d0+delta(j1_t,j2_t)))
      END SUBROUTINE WF_ASYM_TOP_ASYM_TOP_IDENT
	  
      SUBROUTINE WF_VIB_DIAT_DIAT
     & (state_num,
     & i_vib1,i_vib2,alpha,beta,gamma,wfr,wfi )
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE POT_STORE	  
      IMPLICIT NONE
      REAL*8 dr1,dr2,wfr,wfi
      REAL*8 alpha,beta,gamma
      INTEGER channel_state,i_vib1,i_vib2,v1_chann,v2_chann	  
      INTEGER state_num,j12_state,m12_state,j1_state,j2_state
      j12_state = j12(state_num)
      m12_state = m12(state_num)
      channel_state = indx_chann(state_num)
      j1_state = j1_ch(channel_state)
      j2_state = j2_ch(channel_state)
      v1_chann = v1_ch(channel_state)
      v2_chann = v2_ch(channel_state)	  
      CALL WF_DIAT_DIAT
     & (j12_state,m12_state,
     & j1_state,j2_state,beta,gamma,alpha,wfr,wfi)
      dr1 = wf_vib_part_diatom1(i_vib1,channel_state)
      dr2 = wf_vib_part_diatom2(i_vib2,channel_state)	  
      wfr = dr1*dr2*wfr
      wfi = dr1*dr2*wfi	  
      END SUBROUTINE WF_VIB_DIAT_DIAT

      SUBROUTINE WF_VIB_DIAT_DIAT_IDENT
     & (state_num,
     & i_vib1,i_vib2,alpha,beta,gamma,wfr,wfi,parity_ident )
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE POT_STORE	  
      IMPLICIT NONE
      REAL*8 dr1,dr2,wfr,wfi,wfr1,wfr2,wfi1,wfi2
      REAL*8 alpha,beta,gamma,delta
      INTEGER channel_state,i_vib1,i_vib2,v1_chann,v2_chann,
     & parity_ident		  
      INTEGER state_num,j12_state,m12_state,j1_state,j2_state
      EXTERNAL delta 
      j12_state = j12(state_num)
      m12_state = m12(state_num)
      channel_state = indx_chann(state_num)
      j1_state = j1_ch(channel_state)
      j2_state = j2_ch(channel_state)
      v1_chann = v1_ch(channel_state)
      v2_chann = v2_ch(channel_state)	  
      CALL WF_DIAT_DIAT(j12_state,m12_state,j1_state,j2_state,
     & beta,gamma,alpha,wfr1,wfi1)
      dr1 = wf_vib_part_diatom1(i_vib1,channel_state)
      dr2 = wf_vib_part_diatom2(i_vib2,channel_state)	 
      wfr1 = wfr1*dr1*dr2
      wfi1 = wfi1*dr1*dr2 	 
      CALL WF_DIAT_DIAT(j12_state,m12_state,j2_state,j1_state,
     & beta,gamma,alpha,wfr2,wfi2)
      dr1 = wf_vib_part_diatom2(i_vib1,channel_state)
      dr2 = wf_vib_part_diatom1(i_vib2,channel_state)
      wfr2 = wfr2*dr1*dr2
      wfi2 = wfi2*dr1*dr2	  
      wfr = wfr1+wfr2*parity_ident*(-1)**j12_state
      wfi = wfi1+wfi2*parity_ident*(-1)**j12_state      	  
      wfr = wfr/
     & dsqrt(2d0*(1d0+delta(j1_state,j2_state)
     & *delta(v1_chann,v2_chann)))
      wfi = wfi
     & /dsqrt(2d0*(1d0+delta(j1_state,j2_state)
     & *delta(v1_chann,v2_chann)))

	  
      END SUBROUTINE WF_VIB_DIAT_DIAT_IDENT
      
      SUBROUTINE WF_SYM_TOP_DIAT_TOP(wfr,wfi,i_state,
     & alpha1,beta1,gamma1,alpha2,beta2,gamma2)
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      IMPLICIT NONE
      REAL*8 wfr,wfi,dr1,dr2,di1,di2,clb_grd
      INTEGER i_state,chann,m12_t,j12_t
     &,j2_t,m2_t,j1_t,k1_t,eps1_t,m1_t
      REAL*8 alpha1,beta1,gamma1,alpha2,beta2,gamma2
      REAL*8 CG	  
      EXTERNAL CG
      wfr = 0d0
      wfi = 0d0	  
      chann = indx_chann(i_state)
      j12_t = j12(i_state)	  
      m12_t = m12(i_state)		  
      j1_t = j1_ch(chann)
      j2_t = j2_ch(chann)
      k1_t = k1_ch(chann)
      eps1_t = eps1_ch(chann)
	  	  
      IF(max(j12_t,j1_t,j2_t)*2 .gt.j12_t+j1_t+j2_t )
     &	  STOP"ERROR:TRINAGULAR RULE"
      DO m1_t = -j1_t,j1_t
      m2_t = m12_t - m1_t
      IF(abs(m2_t).le.j2_t) THEN
      CALL WF_SYM_TOP(j1_t,k1_t,eps1_t,m1_t,dr1,di1,
     & alpha1,beta1,gamma1)
      CALL WF_DIAT_TOP(j2_t,m2_t,beta2,alpha2,dr2,di2)
      clb_grd	= CG(j1_t,j2_t,j12_t,m1_t,m2_t,m12_t)
       wfr =wfr +  clb_grd*(dr1*dr2 - 
     &	di1*di2)
      wfi =wfi +  clb_grd*(dr1*di2 + 
     &	di1*dr2)	 
      ENDIF	  
      ENDDO	  	 
      END SUBROUTINE WF_SYM_TOP_DIAT_TOP
      SUBROUTINE WF_ASYM_TOP_DIAT_TOP(wfr,wfi,i_state,
     & alpha1,beta1,gamma1,alpha2,beta2,gamma2)
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      IMPLICIT NONE
      REAL*8 wfr,wfi,dr1,dr2,di1,di2,clb_grd
      INTEGER i_state,chann,j1_t,ka1_t,kc1_t,m1_t,m12_t,j12_t
     &,j2_t,m2_t
      REAL*8 alpha1,beta1,gamma1,alpha2,beta2,gamma2
      REAL*8 CG	  
      EXTERNAL CG
      wfr = 0d0
      wfi = 0d0	  
      chann = indx_chann(i_state)
      j12_t = j12(i_state)	  
      m12_t = m12(i_state)		  
      j1_t = j1_ch(chann)
      j2_t = j2_ch(chann)
      ka1_t = ka1_ch(chann)
      kc1_t = kc1_ch(chann)
      IF(max(j12_t,j1_t,j2_t)*2 .gt.j12_t+j1_t+j2_t )
     &	  STOP"ERROR:TRINAGULAR RULE"
      DO m1_t = -j1_t,j1_t
      m2_t = m12_t - m1_t
      IF(abs(m2_t).le.j2_t) THEN
      CALL WF_ASYM_TOP(A,B,C,j1_t,1,chann,m1_t,dr1,di1,
     & alpha1,beta1,gamma1)
      CALL WF_DIAT_TOP(j2_t,m2_t,beta2,alpha2,dr2,di2)
      clb_grd	= CG(j1_t,j2_t,j12_t,m1_t,m2_t,m12_t)
       wfr =wfr +  clb_grd*(dr1*dr2 - 
     &	di1*di2)
      wfi =wfi +  clb_grd*(dr1*di2 + 
     &	di1*dr2)	 
      ENDIF	  
      ENDDO	  	 
      END SUBROUTINE WF_ASYM_TOP_DIAT_TOP

      SUBROUTINE WF_DIAT_TOP_FINE(wfr,wfi,
     & spin_inp,i_state_inp,theta,phi,m_fine_proj)
! This subroutine is updated by Bikramaditya Mandal
      USE VARIABLES
      USE FINE_STRUCT_LEVELS
      IMPLICIT NONE
      INTEGER i_state_inp	  
      INTEGER i_chann_inp,spin_inp,m_fine_proj,m_n_proj
      REAL*8 theta,phi	 
      REAL*8 wfr(2*spin_inp+1),	wfi(2*spin_inp+1)
	  INTEGER i_spin_proj
      INTEGER j_fine_inp,f_fine_inp,n_fine_inp	 
      REAL*8 CG,clb_grd,buffer_r,buffer_i
 
      EXTERNAL CG
	   i_chann_inp =  indx_chann(i_state_inp)
      j_fine_inp = j_ch(i_chann_inp)	
      wfr = 0d0
      wfi = 0d0
      DO f_fine_inp=-1,1 !!! CHANGE IT
      n_fine_inp = j_fine_inp+f_fine_inp !!! VALUE ADDING HALF INTEGER
      DO i_spin_proj=-spin_inp,spin_inp 
      m_n_proj = m_fine_proj-i_spin_proj
      IF(n_fine_inp.lt.abs(m_n_proj)) CYCLE	 
      clb_grd	= CG(n_fine_inp,spin_inp,j_fine_inp,
     &  m_n_proj,i_spin_proj,m_fine_proj)
      CALL WF_DIAT_TOP(n_fine_inp,m_n_proj,theta,phi
     &  ,buffer_r,buffer_i)
      wfr(i_spin_proj+spin_inp+1) = wfr(i_spin_proj+spin_inp+1)+
     & buffer_r*clb_grd*
     &  b_fine_coeff(f_fine_inp+2,i_chann_inp)
      wfi(i_spin_proj+spin_inp+1) = wfi(i_spin_proj+spin_inp+1)+
     & buffer_i*clb_grd*
     & b_fine_coeff(f_fine_inp+2,i_chann_inp)   
      ENDDO	
      ENDDO	 
	 
     	 
      END SUBROUTINE WF_DIAT_TOP_FINE	 
