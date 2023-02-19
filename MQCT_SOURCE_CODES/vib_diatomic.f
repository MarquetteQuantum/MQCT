      SUBROUTINE energy_diatom(E,v,j,we,xe,be,de)
      IMPLICIT NONE	  
      INTEGER v,j,vmax
      REAL*8 E,we,xe,be,de
    
      IF(xe.gt.we) STOP "ERROR: IT'S NOT A  DIATOMIC MOLECULE:XE>we"
      vmax = int((1d0-xe/we)/2d0/(xe/we))	  
      IF(v.gt.vmax) THEN
       
      WRITE(*,'(a46,i3,1x,a6,i3)')
     & "ERROR: VIBRATION MODE V  IS MORE THAN V_MAX = ",
     & vmax, ", V = ", v
      ENDIF	 
      E = we*(v+0.5d0) - xe*(v+0.5d0)**2 + Be*dble(j*(j+1)) - 
     & de*dble(j*(j+1))**2
      
      END SUBROUTINE energy_diatom 
	  
      SUBROUTINE VIBRATIONAL_WAVEFUNCTION
     & (wpsir,r,n_grid,n_vib,nju_red,De,a_eq,r_eq)
      USE FACTORIAL	 
      IMPLICIT NONE

      INTEGER n_grid,i,n_vib
      REAL*8 wpsir(n_grid),r(n_grid),nju_red,De,a_eq,r_eq
      REAL*8 z(n_grid),x(n_grid),xe,lambda,ev,Nn,alhpa,gammln
      REAL*8 lag(1:n_grid,0:n_vib),alpha
      EXTERNAL gammln
      lambda = dsqrt(2d0*nju_red*De)/a_eq
c      PRINT*,"lambda",lambda
      alpha = 2d0*lambda-2d0*n_vib-1
       
c      PRINT*,"alpha",alpha	  
      Nn = 0.5d0*(LOGFACT(n_vib)+log(2d0*lambda-2*n_vib-1)-
     & gammln(2d0*lambda-n_vib))
c      PRINT*,"Nn",Nn,
c     * log(2d0*lambda-2*n_vib-1),gammln(2d0*lambda-n_vib)
c      STOP
c      PRINT*,"r_1",r(1)
      xe = a_eq*r_eq	  
      DO i=1,n_grid	 
      x(i) = a_eq*r(i)
      z(i) = 2d0*lambda*exp(xe-x(i))
      ENDDO
c cc     PRINT*,"r1",r(1)
!      STOP	  
      CALL lf_function ( n_grid, n_vib, alpha, z, lag )
      DO i=1,n_grid	 
      wpsir(i) = Nn+log(z(i))*(alpha/2d0) -z(i)/2d0!+ log(lag(i,n_vib)
      wpsir(i) = exp(wpsir(i))*lag(i,n_vib)	  
      ENDDO	
      wpsir = wpsir*dsqrt(a_eq)	  
	  
	  
      END SUBROUTINE VIBRATIONAL_WAVEFUNCTION

      SUBROUTINE FINE_ENERGY_DIATOMIC(E,rot_const,
     & rot_state,vector_cff,SPIN )
      IMPLICIT NONE
      REAL*8 rot_const(11)
      INTEGER rot_state(6)	  
      REAL*8 E,A,B,lambda,gamma,p,q
      REAL*8 gammaD,lambdaDD, lambdaD,D	,H  	  
      INTEGER L,SPIN,N,F,SIGMA,J
      REAL*8 vector_cff(SPIN)
      REAL*8 E1,E2,E3,a11,a12,a22,norm,j_h
      REAL*8 m11,m12,m22,m33,m44,m34,e11,e22
      A = rot_const(1)! = a_spin_orbit_fine
      B = rot_const(2)! = Be
      D = rot_const(3)! = De
      lambda = rot_const(4)! = lambda_fine
      gamma = rot_const(5)! = gamma_fine
      lambdaD = rot_const(6)! = lambda_fine_d
      gammaD = rot_const(7)! = gamma_fine_d
      lambdaDD = rot_const(8)! = lambda_fine_dd
      H = rot_const(9)! = lambda_fine_dd
      p = rot_const(10)! = p_double_fine
      q = rot_const(11)! = q_double_fine		  
      L = rot_state(1)! = LORB_FINE
      SPIN = rot_state(2)! = SPIN_FINE
      N = rot_state(3)! = 0
      F = rot_state(4)! = f_ch(i)
      SIGMA = rot_state(5)! = 0
      J = rot_state(6)! = j_ch(i)
      vector_cff  =0d0	  
      SELECT CASE(L)
      CASE(0) !!! Hund's b with departure	  
      !!! J is input,F is input, SPIN is input, S=1 ADD CASE FURTHER
      vector_cff = 0d0	  
      E2 = B*(J+1)*J+2d0/3d0*lambda-gamma - D*(J+1)**2*J**2-
     & gammaD*(J+1)*J + 2d0/3d0*lambdaD*J*(J+1) 
     & + 2d0/3d0*lambdaDD*J**2*(J+1)**2  + H*(J+1)**3*J**3!!! MAIN DIAGONAL

	 
      a11 = B*(J+1)*(J+2)-2d0/3d0*lambda*(J+2d0)/(2*J+1)-
     & gamma*(J+2) - D*(J+1)**2*(J+2)**2 - gammaD*(J+1)*(J+2)**2
     &  - 2d0/3d0*lambdaD*(J+2d0)/(2*J+1)*(J+1)*(J+2)
     & -2d0/3d0*lambdaDD*(J+2d0)/(2*J+1)*(J+1)**2*(J+2)**2
     &   + H*(J+2)**3*(J+1)**3	 !!! N = J+1
	 
      a22 = B*J*(J-1)-2d0/3d0*lambda*(J-1d0)/(2*J+1)+
     & gamma*(J-1)-D*(J-1)**2*J**2+gammaD*(J-1)**2*J
     & - 2d0/3d0*lambdaD*(J-1d0)/(2*J+1)*(J-1)*J
     & - 2d0/3d0*lambdaDD*(J-1d0)/(2*J+1)*(J-1)**2*J**2
     & + H*(J-1)**3*J**3
	 !!! N = J-1
      a12 = 2*lambda*sqrt((J+1d0)*J/(2*J+1d0)**2)+
     & lambdaD*sqrt((J+1d0)*J/(2*J+1d0)**2)*
     & ((J-1)*J+(J+1)*(J+2))+
     & lambdaDD*sqrt((J+1d0)*J/(2*J+1d0)**2)*
     & ((J-1)**2*J**2+(J+1)**2*(J+2)**2)
	  IF(j.eq.0) a22 = 0d0
      E1 = (a11+a22)/2d0 - sqrt((a11-a22)**2+4*a12**2)/2
      E3 = (a11+a22)/2d0 + sqrt((a11-a22)**2+4*a12**2)/2
      n = j + f
!      if(mod(n,2).ne.1 .or. n.lt. 0) THEN
!      PRINT*,"STATE FOR THE GIVEN FINE LEVEL DOES NOT EXIST" !!! THINK ABOUT IT
!      STOP	  
!      ENDIF
      if(f.eq.0) then
	  E = E2
	  vector_cff(2) = 1d0
      endif	  
      if(f.eq.-1) then
      E = E1
	  vector_cff(1) = 1d0
      IF(j.ne.0) then	  
      vector_cff(3) =  (E1-a22)/a12
      endif	  
      norm = 
     & sqrt(vector_cff(1)**2 + vector_cff(2)**2 + vector_cff(3)**2)
      vector_cff = vector_cff/norm	  
      endif
      if(f.eq.1) then
	  E = E3
      IF(j.ne.0) then	  
	  vector_cff(1) = (E3-a11)/a12
      endif	  
      vector_cff(3) = 1d0
      norm =
     & sqrt(vector_cff(1)**2 + vector_cff(2)**2 + vector_cff(3)**2)
      vector_cff = vector_cff/norm 	  
      endif	  
      CASE(1)
      j_h = j - (-1d0)**F/2d0
!!! SPIN = 1/2	
      m11 = A/2 + B*(j_h*(j_h+1)-3.d0/4.0d0)    
      m33 = m11 !!! P3/2
	  m22 = -A/2 + B*(j_h*(j_h+1)+5.d0/4.0d0)
     & - gamma -1/2d0*(p+2*q)*(j_h+0.5d0) !P 1/2 +
	  m44 = -A/2 + B*(j_h*(j_h+1)+5.d0/4.0d0)
     & - gamma +1/2d0*(p+2*q)*(j_h+0.5d0)!P 1/2 -
      m12 = - sqrt((j_h+3.0/2.0)*(j_h-1.0/2.0))
     & *((B-gamma/2)-q/2*(j_h+0.5d0))
      m34 = - sqrt((j_h+3.0/2.0)*(j_h-1.0/2.0))
     & *((B-gamma/2)+q/2*(j_h+0.5d0))
      if(F.eq.2 .and. j.eq.1) then
      vector_cff(1) = 1.0
      vector_cff(2) = 0.0	  
      select case(sigma)
      case(0)
      E = m22	  
      case(1)
      E = m44	  
      end select	  
      return
      endif	
       !!!!vector_cff(1) - P1/2
       !!!!vector_cff(1) - P3/2	   
      select case(sigma)
      case(0)
      E11 = (m11 + m22 - sqrt((m11-m22)**2 + m12**2*4))/2
      E22 = (m11 + m22 + sqrt((m11-m22)**2 + m12**2*4))/2
	  if(F.eq.1) then
      vector_cff(1) = 1d0
      vector_cff(2) = (E11 - m22)/m12
	  norm = sqrt(vector_cff(1)**2 + vector_cff(2)**2)
	  vector_cff = vector_cff/norm
	  E = E11   
      endif
	  if(F.eq.2) then
      vector_cff(1) = 1d0
      vector_cff(2) = (E22 - m22)/m12
	  norm = sqrt(vector_cff(1)**2 + vector_cff(2)**2)
	  vector_cff = vector_cff/norm
	  E = E22   
      endif	  
      case(1)
	  if(f.eq.1) then
      E11 = (m33 + m44 - sqrt((m33-m44)**2 + m34**2*4))/2
      E22 = (m33 + m44 + sqrt((m33-m44)**2 + m34**2*4))/2
      vector_cff(1) = 1d0
      vector_cff(2) = (E11 - m44)/m34
	  norm = sqrt(vector_cff(1)**2 + vector_cff(2)**2)
	  vector_cff = vector_cff/norm	
	  E = E11 	  
      endif
	  if(f.eq.2) then
      E11 = (m33 + m44 - sqrt((m33-m44)**2 + m34**2*4))/2
      E22 = (m33 + m44 + sqrt((m33-m44)**2 + m34**2*4))/2
      vector_cff(1) = 1d0
      vector_cff(2) = (E22 - m44)/m34
	  norm = sqrt(vector_cff(1)**2 + vector_cff(2)**2)
	  vector_cff = vector_cff/norm
	  E = E22 	  
      endif	  	  
      end select   
	  
      return
	  
      END SELECT	  
	 
      END SUBROUTINE FINE_ENERGY_DIATOMIC	 