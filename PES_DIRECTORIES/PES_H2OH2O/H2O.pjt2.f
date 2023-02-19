      SUBROUTINE POTS(V,Q1,Q2,THETA)
C     Potential PJT2 due Polyansky, Jensen and Tennyson, 
C     J. Chem. Phys., 105, 6490-6497 (1996)
C     Update of Polyansky, Jensen and Tennyson, J Chem Phys 101, 7651 (1994))
C     Units: Hartree and Bohr
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     RZ = OH equilibrium value
C     RHO = equilibrium value of pi - bond angle(THETA)
      DATA TOANG/0.5291772/, CMTOAU/219474.624/
      DATA X1/1.0/
      DATA RHO1    /    75.50035308/
      DATA FA1     /     .00000000/
      DATA FA2     /18902.44193433/
      DATA FA3     /  1893.99788146/
      DATA FA4     /  4096.73443772/
      DATA FA5     /-1959.60113289/
      DATA FA6     /  4484.15893388/
      DATA FA7     /  4044.55388819/
      DATA FA8     / -4771.45043545/
      DATA FA9     /     0.00000000/
      DATA FA10    /     0.00000000/
      DATA RZ    /     .95792059/
      DATA A     /    2.22600000/
      DATA F1A1    /  -6152.40141181/
      DATA F2A1    / -2902.13912267/
      DATA F3A1    / -5732.68460689/
      DATA F4A1    /  953.88760833/
      DATA F11     / 42909.88869093/
      DATA F1A11   /  -2767.19197173/
      DATA F2A11   /  -3394.24705517/
      DATA F3A11   /     .00000000/
      DATA F13     /  -1031.93055205/
      DATA F1A13   /  6023.83435258/
      DATA F2A13   /     .00000000/
      DATA F3A13   /     .00000000/
      DATA F111    /     .00000000/
      DATA F1A111  /   124.23529382/
      DATA F2A111  /  -1282.50661226/
      DATA F113    /  -1146.49109522/
      DATA F1A113  /  9884.41685141/
      DATA F2A113  /  3040.34021836/ 
      DATA F1111   /  2040.96745268/
      DATA FA1111  /  .00000000/
      DATA F1113   /  -422.03394198/
      DATA FA1113  /-7238.09979404/
      DATA F1133   /     .00000000/
      DATA FA1133  /     .00000000/
      DATA F11111  / -4969.24544932/
      DATA f111111/  8108.49652354/
      DATA F71   /  90.00000000/
 
      data c1/50.0/,c2/10.0/,beta1/22.0/,beta2/13.5/,gammas/0.05/,
     *     gammaa/0.10/,delta/0.85/,rhh0/1.40/
                 RHO=RHO1*3.141592654/180.000000000
      fa11=0.0
      f1a3=f1a1
      f2a3=f2a1
      f3a3=f3a1
      f4a3=f4a1
      f33=f11
      f1a33=f1a11
      f2a33=f2a11
      f333=f111
      f1a333=f1a111
      f2a333=f2a111
      f133=f113
      f1a133=f1a113
      f2a133=f2a113
      f3333=f1111
      fa3333=fa1111
      f1333=f1113
      fa1333=fa1113
      f33333=f11111
      f333333 =f111111
      f73     =f71
 
C     Find value for DR and DS
      DR = TOANG*Q1 - RZ
      DS = TOANG*Q2 - RZ
 
C     Transform to Morse coordinates
      Y1 = X1 - EXP(-A * DR)
      Y3 = X1 - EXP(-A * DS)
 
C     transform to Jensens angular coordinate
      CORO = DCOS(THETA) + DCOS(RHO)
 
C     Now for the potential
      V0=(FA2+FA3*CORO+FA4*CORO**2+FA6*CORO**4+FA7*CORO**5)*CORO**2
      V0=V0+(FA8*CORO**6+FA5*CORO**3+FA9*CORO**7+FA10*CORO**8 )*CORO**2
      V0=V0+(                                    FA11*CORO**9 )*CORO**2
      FE1= F1A1*CORO+F2A1*CORO**2+F3A1*CORO**3+F4A1*CORO**4
      FE3= F1A3*CORO+F2A3*CORO**2+F3A3*CORO**3+F4A3*CORO**4
      FE11= F11+F1A11*CORO+F2A11*CORO**2
      FE33= F33+F1A33*CORO+F2A33*CORO**2
      FE13= F13+F1A13*CORO
      FE111= F111+F1A111*CORO+F2A111*CORO**2
      FE333= F333+F1A333*CORO+F2A333*CORO**2
      FE113= F113+F1A113*CORO+F2A113*CORO**2
      FE133= F133+F1A133*CORO+F2A133*CORO**2
      FE1111= F1111+FA1111*CORO
      FE3333= F3333+FA3333*CORO
      FE1113= F1113+FA1113*CORO
      FE1333= F1333+FA1333*CORO
      FE1133=       FA1133*CORO
      FE11111=F11111
      FE33333=F33333
      FE111111=F111111
      FE333333=F333333
      FE71    =F71
      FE73    =F73
      V   = V0 +  FE1*Y1+FE3*Y3
     1         +FE11*Y1**2+FE33*Y3**2+FE13*Y1*Y3
     2         +FE111*Y1**3+FE333*Y3**3+FE113*Y1**2*Y3
     3         +FE133*Y1*Y3**2
     4         +FE1111*Y1**4+FE3333*Y3**4+FE1113*Y1**3*Y3
     5         +FE1333*Y1*Y3**3+FE1133*Y1**2*Y3**2
     6         +FE11111*Y1**5+FE33333*Y3**5
     7         +FE111111*Y1**6+FE333333*Y3**6
     8         +FE71    *Y1**7+FE73    *Y3**7
C     modification by Choi & Light, J. Chem. Phys., 97, 7031 (1992).
      sqrt2=sqrt(2.0)
      xmup1=sqrt2/3.0+0.5
      xmum1=xmup1-x1
      term=2.0*xmum1*xmup1*q1*q2*cos(theta)
      r1=toang*sqrt((xmup1*q1)**2+(xmum1*q2)**2-term)
      r2=toang*sqrt((xmum1*q1)**2+(xmup1*q2)**2-term)
      rhh=sqrt(q1**2+q2**2-2.0*q1*q2*cos(theta))
      rbig=(r1+r2)/sqrt2
      rlit=(r1-r2)/sqrt2
 
      alpha=(x1-tanh(gammas*rbig**2))*(x1-tanh(gammaa*rlit**2))
      alpha1=beta1*alpha
      alpha2=beta2*alpha
      drhh=toang*(rhh-delta*rhh0)
      DOLEG=     (1.4500-THETA)
C     IF (THETA.LE.0.64  ) V=0.1E17
C     IF((DR.LE.-0.4).AND.(THETA.LE.1.1)) V=0.1E17
C     IF((DS.LE.-0.4).AND.(THETA.LE.1.1)) V=0.1E17
C     IF (DS.LE. 0.0  ) V=0.1E17
      v = v + c1*exp(-alpha1*drhh) + c2*exp(-alpha2*drhh)
 
C     Convert to Hartree
      V=V/CMTOAU
      RETURN
      END
