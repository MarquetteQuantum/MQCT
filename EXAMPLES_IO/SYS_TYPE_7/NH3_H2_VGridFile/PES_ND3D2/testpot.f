c program testpot
c computes nd3-h2 potential in Jacobi coordinates
c First make dummy call to nd3h2pot.
c Then ransform Jacobi coordinates to coordinates defined in Maret et al.,
c   MNRAS, 399, 425–431 (2009).
c Subroutine nd3h2pot first transforms coordinates of ND3-H2 to those of
c   NH3-H2 and then calls nh3h2pot.
      implicit none
      real*8 r, thx, phx, thxp, phxp, energy, x
      real*8 theta, phi, thetap, phip, rad2deg, au2cm, pi
      real*8 alpha1,beta1,gamma1,alpha2,beta2, b1, b2, g1, a2
      parameter(au2cm=219474.63137098d0)
c
      pi=dacos(-1d0)
      rad2deg=180d0/pi
c call first to initialize
      call nd3h2pot(5d0, 1d0, 1d0, 1d0, 1d0, x)
c Distance in a_0 at minimum of -266.7 cm-1
      r = 6.216763d0
c Euler angles in degrees (at equilibrium geometry)
      b1 = 0d0
      g1 = 0d0
      b2 = 0d0
      a2 = 0d0
c convert angles to radians
      alpha1 = 0d0
      beta1 = b1/rad2deg
      gamma1 = g1/rad2deg
      beta2 = b2/rad2deg
      alpha2 = a2/rad2deg
c  The Euler angles used here differ from those used in the potential
c  Coordinate transformation
      call convert(alpha1,beta1,gamma1,alpha2,beta2,
     &              theta,phi,thetap,phip)
c Angles defined in Maret et al., MNRAS, 399, 425–431 (2009)
c convert to degrees
      thx = theta*rad2deg
      phx = phi*rad2deg
      thxp = thetap*rad2deg
      phxp = phip*rad2deg
c call with angles in degrees
      call nd3h2pot(r, thx, phx, thxp, phxp, energy)
c energy in cm-1, should be -266.7 cm-1
      write(6,'(5f12.6,f20.12)') r, b1, g1, b2, a2, energy
c divide by au2cm to get Hartree, not done here
      end
