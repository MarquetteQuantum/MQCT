      program fit5d
C -- NORMALISATION MOLSCAT POUR LES COEFFICIENTS  ----
C---------------------------------------------------------
      implicit none
      integer mxlmb, i, ityp
      parameter(mxlmb=10000)
      integer icntrl, nlam, lam(4,mxlmb), p1, q1, p2, p
      integer mxlam
      character dir*80
      real*8 R, COEF(mxlmb),
     $     thx, phx, thpx, phpx, tnormed, y, t
      INTEGER :: N_POINTS = 200
      INTEGER I_R
      REAL*8 :: R_MIN = 3d0
      REAL*8 :: R_MAX = 30d0
      REAL*8 DR_GRID
      REAL*8, ALLOCATABLE :: TERMS(:,:)	  
      logical corr_9d_flag, r12_flag
c-------------------------------- INITIALIZATION
      icntrl=-1
      dir='./H2Ov000_H2v0/'
!      corr_9d_flag=.true.
!      r12_flag=.true.
c      write(6,'(a,$)')
c     $     'Enter corr_9d_flag (t/f) and r12_flag (t/f) ? '
c      read(5,*) corr_9d_flag, r12_flag
      r=0
      do i=1,mxlmb
         coef(i)=0
      enddo
      mxlam = mxlmb
      call potenl(icntrl, 
     $     mxlam, nlam, lam, r, coef,ityp)
      write(6,*) 'NLAM', nlam,'  MXLAM', MXLAM

c-------------------------------- PES EVALUATION
      icntrl=0
      r=2.5
      DR_GRID = (R_MAX-R_MIN)/(N_POINTS-1d0)
      ALLOCATE(TERMS(N_POINTS,nlam))	  
      do I_R=1,N_POINTS
!         write(6,'(a,$)') 'Distance and angles (degrees) ? '
!         read(5,*,end=999) r, thx, phx, thpx, phpx
      thx= 0d0
      phx = 0d0
      thpx = 0d0
      phpx = 0d0	  
      r = R_MIN + DR_GRID*(I_R-1)
         call potenl(icntrl, 
     $        mxlam, nlam, lam, r, coef,ityp)
         y=0
         do i=1,nlam
            p1=lam(1,i)
            q1=lam(2,i)
            p2=lam(3,i)
            p =lam(4,i)
            y=y+coef(i)*t(p1,q1,p2,p,thx,phx,thpx,phpx)
         TERMS(I_R,i) = coef(i)			 
         enddo
!         write(6,'(5f10.3,f15.3)') r, thx, phx, thpx, phpx, y
      enddo
      OPEN(1,FILE="TERMS_H2O_H2.DAT",ACTION="WRITE")
      DO i=1,nlam
      IF(i.eq.1) THEN
      WRITE(1,'(a6,1x,i3,a1,1x,i3,a1,1x,i3,a1,1x,i3,a1,3x)',
     & ADVANCE = "NO" )
     & "R,a.u.",lam(1,i),',',lam(2,i),',',lam(3,i),',',lam(4,i),','
      ELSE	 
      WRITE(1,'(i3,a1,1x,i3,a1,1x,i3,a1,1x,i3,a1,3x)',
     & ADVANCE = "NO" )
     & lam(1,i),',',lam(2,i),',',lam(3,i),',',lam(4,i),','
      ENDIF	 
      ENDDO
      WRITE(1,*)	  
      DO I_R=1,N_POINTS	  
      DO i=1,nlam
      IF(i.eq.1) THEN	  
      WRITE(1,'(f6.3,1x,e19.12,3x)',ADVANCE="NO")
     &  R_MIN + DR_GRID*(I_R-1),TERMS(I_R,i)
      ELSE
      WRITE(1,'(e19.12,3x)',ADVANCE="NO")
     & TERMS(I_R,i)!,
      ENDIF	  
      ENDDO
      WRITE(1,*) 
      ENDDO
      CLOSE(1)	  
      end

