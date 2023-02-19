c program to generate the angles theta, phi, theta', phi'
c used by Valiron et al. for H2O-H2
c from Euler angles in the two-angle embedded BF frame
      subroutine convert(alphaA,betaA,gammaA,alphaB,betaB,
     &                   theta,phi,thetap,phip)
      implicit real*8 (a-h,o-z)
      dimension rBG(3), rB(3), rmat(3,3)
      pi = dacos(-1d0)
      rad2deg = 180d0/pi
c      theta = -betaA
c      phi   = -gammaA
      theta = betaA
      phitmp   = pi-gammaA
      if (phitmp .ge. 0d0) then
        phi = phitmp
      else
        phi = phitmp + 2d0*pi
      endif
      rBG(1) = dsin(betaB)*dcos(alphaB)
      rBG(2) = dsin(betaB)*dsin(alphaB)
      rBG(3) = dcos(betaB)
c inverse Euler rotation (check with transposed)
c      call eulerrot(rmat,-gammaA,-betaA,-alphaA)
      call eulerrot(rmat,alphaA,betaA,gammaA)
      do i = 1,3
        rB(i) = 0d0
        do j = 1,3
c          rB(i) = rB(i) + rmat(i,j)*rBG(j)
          rB(i) = rB(i) + rmat(j,i)*rBG(j)
        enddo
      enddo
      cost = rB(3)
      if (rB(1)**2+rB(2)**2 .lt. 1d-20) then
        if (cost .gt. 0d0) then
          thetap = 0d0
        else
          thetap = pi
        endif
        phip = 0d0
      else
        phitmp = datan2(rB(2),rB(1))
        cosp = dcos(phitmp)
        sinp = dsin(phitmp)
        if (dabs(sinp) .ge. 1d-10) then
          sint = rB(2)/sinp
        else
          sint = rB(1)/cosp
        endif
        thetatmp = datan2(sint,cost)
        if (thetatmp .ge. 0d0) then
          thetap = thetatmp
          if (phitmp .ge. 0d0) then
            phip = phitmp
          else
            phip = phitmp + 2d0*pi
          endif
        else
          thetap = -thetatmp
          phip = phitmp + pi
        endif
      endif
      end

      subroutine disteul(rdist,alpha1,beta1,gamma1,alpha2,beta2, dist)
      implicit none
      real*8 sitA(3,3), sitB(3,2), rot(3,3), sum
      real*8 rdist,alpha1,beta1,gamma1,alpha2,beta2,gamma2
      real*8 dist(3,2), coor1(3,3), coor2(3,2)
      integer i, j, k

c     cartesian coordinates of the sites
      data sitA/
     +  0.0d0,              0.0d0,  .1255334885d0, ! O
     + -1.45365196228170d0, 0.0d0, -.9961538357d0, ! H1
     +  1.45365196228170d0, 0.0d0, -.9961538357d0/ ! H2
      data sitB/
     +  0.0d0, 0.0d0,  0.725d0, ! H3
     +  0.0d0, 0.0d0, -0.725d0/ ! H4

c     rotate monomer A

      call eulerrot(rot,alpha1,beta1,gamma1)
      do i=1,3
         do k=1,3
            sum=0d0
            do j=1,3
               sum=sum+rot(i,j)*sitA(j,k)
            enddo
            coor1(i,k)=sum
         enddo
      enddo

      gamma2 = 0d0
      call eulerrot(rot,alpha2,beta2,gamma2)
      do i=1,3
         do k=1,2
            sum=0d0
            do j=1,3
               sum=sum+rot(i,j)*sitB(j,k)
            enddo
            coor2(i,k)=sum
         enddo
      enddo

      do k=1,2
         coor2(3,k)=coor2(3,k)+rdist
      enddo

c write cartesian coordinates as Molden input file
c     open(unit=10,file='cart.dat',form='formatted')
c
c     write(10,'(a)') '[Molden Format]'
c     write(10,'(a)') '[Atoms] AU     '
c     write(10,'(a,3f12.6)') 'O1  1  8 ', (coor1(i,1), i= 1, 3)
c     write(10,'(a,3f12.6)') 'H1  2  1 ', (coor1(i,2), i= 1, 3)
c     write(10,'(a,3f12.6)') 'H2  3  1 ', (coor1(i,3), i= 1, 3)
c     write(10,'(a,3f12.6)') 'H3  4  1 ', (coor2(i,1), i= 1, 3)
c     write(10,'(a,3f12.6)') 'H4  5  1 ', (coor2(i,2), i= 1, 3)
c
c     close(unit=10)

c     compute distance matrix

      do i=1,3
         do j=1,2
            sum=0d0
            do k=1,3
               sum=sum+(coor1(k,i)-coor2(k,j))**2
            end do
            dist(i,j)=dsqrt(sum)
         end do
      end do

      end

      subroutine distgren(rdist,theta,phi,thetap,phip,dist)
      implicit none
      real*8 sitA(3,3), sitB(3,2), rot(3,3)
      real*8 rdist,theta,phi,thetap,phip,sum
      real*8 dist(3,2), coor1(3,3), coor2(3,2)
      integer i, j, k

c     cartesian coordinates of the sites
      data sitA/
     +  0.0d0,              0.0d0,  .1255334885d0, ! O
     + -1.45365196228170d0, 0.0d0, -.9961538357d0, ! H1
     +  1.45365196228170d0, 0.0d0, -.9961538357d0/ ! H2
      data sitB/
     +  0.0d0, 0.0d0,  0.725d0, ! H3
     +  0.0d0, 0.0d0, -0.725d0/ ! H4

      do i=1,3
         do k=1,3
            coor1(i,k)=sitA(i,k)
         end do
      end do

      call eulerrot(rot,phip,thetap,0d0)
      do i=1,3
         do k=1,2
            sum=0d0
            do j=1,3
               sum=sum+rot(i,j)*sitB(j,k)
            end do
            coor2(i,k)=sum
         end do
      end do

      do k=1,2
         coor2(1,k)=coor2(1,k)+rdist*dsin(theta)*dcos(phi)
         coor2(2,k)=coor2(2,k)+rdist*dsin(theta)*dsin(phi)
         coor2(3,k)=coor2(3,k)+rdist*dcos(theta)
      end do

c     compute distance matrix

      do i=1,3
         do j=1,2
            sum=0d0
            do k=1,3
               sum=sum+(coor1(k,i)-coor2(k,j))**2
            end do
            dist(i,j)=dsqrt(sum)
         end do
      end do

      end

      subroutine eulerrot(rot,alpha,beta,gamma)
c
c calculate rotation matrix given the euler angles (in radians)
c
c note
c       alpha  corresponds to polar coordinate phi
c       beta   corresponds to polar coordinate theta
c       gamma  corresponds to the "spin" angle chi
c
c       rot(alpha,beta,gamma)=Rz(alpha)Ry(beta)Rz(gamma)
c
      implicit real*8 (a-h,o-z)
      real*8 rot(3,3)
c
      cosa=dcos(alpha)
      sina=dsin(alpha)
      cosb=dcos(beta)
      sinb=dsin(beta)
      cosc=dcos(gamma)
      sinc=dsin(gamma)
c
      rot(1,1)=cosb*cosa*cosc-sina*sinc
      rot(2,1)=cosb*sina*cosc+cosa*sinc
      rot(3,1)=-sinb*cosc
      rot(1,2)=-cosb*cosa*sinc-sina*cosc
      rot(2,2)=-cosb*sina*sinc+cosa*cosc
      rot(3,2)=sinb*sinc
      rot(1,3)=sinb*cosa
      rot(2,3)=sinb*sina
      rot(3,3)=cosb
      return
      end

