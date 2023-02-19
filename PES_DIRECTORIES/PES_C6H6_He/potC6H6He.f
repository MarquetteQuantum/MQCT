
	implicit real*8(a-h,o-z)
	theta=30.0d0
	phi=90.0d0
	R=2.0d0
	do1 i=1,1000
	call V(R,theta,phi,pes)
	write(*,*) R,pes
	R=R+0.01d0
1	continue
	stop
	end


 	 subroutine V(R,theta,phi,pes)
c******************************************************
c R (angstr.), theta, phi (degree): RareGas coordinates
c pes: potential energy in meV
c****************************************************** 
	implicit real*8(a-h,o-z)
	common/blk0/iflag
	common/blk1/xxx(12),yyy(12),zzz(12)
	common/blk2/xbarbnd(12),ybarbnd(12),zbarbnd(12)
	common/blk3/xbnd2(12),ybnd2(12),zbnd2(12)
        common/blk4/xbnd1(12),ybnd1(12),zbnd1(12)
        common/blk5/rmperpch,rmparch,epsperpch,epsparch
        common/blk6/rmperpcc,rmparcc,epsperpcc,epsparcc
        
	parameter (PI=3.1415926535897932d0)
         
c set iflag=0 the first call to initialize the pes parameters
        if(iflag.eq.0) call init

        thetarad=theta*pi/180.d0
        phirad=phi*pi/180.d0
	xAr=R*dsin(phirad)*dcos(thetarad)
	yAr=R*dsin(phirad)*dsin(thetarad)
	zAr=R*dcos(phirad)
        
	pes=0.d0

	do i=1,12
        xbnd1(i)=xAr-xbarbnd(i)
        ybnd1(i)=yAr-ybarbnd(i)
        zbnd1(i)=zAr-zbarbnd(i)
	rrab=dsqrt(xbnd1(i)**2+ybnd1(i)**2+zbnd1(i)**2)
	if (i.le.6) then
           den=(rrab*1.39d0)
         else
           den=(rrab*1.09d0)
         endif
        anum=xbnd1(i)*xbnd2(i)+ybnd1(i)*ybnd2(i)+zbnd1(i)*zbnd2(i)
        costhb2=(anum/den)**2
        sinthb2=1-(anum/den)**2
        if(i.le.6) then
        eps=epsparcc*costhb2+epsperpcc*sinthb2
        rm=rmparcc*costhb2+rmperpcc*sinthb2
        pes=pes+pmms(eps,rm,rrab)
        else
        eps=(epsparch*costhb2+epsperpch*sinthb2)
        rm=rmparch*costhb2+rmperpch*sinthb2
        pes=pes+pmms(eps,rm,rrab)
        endif
        enddo

        return
        end

	function pmms(eps,rm,r)
	implicit real*8(a-h,o-z)
        an=6.00d0
        am=10.0d0+4.0d0*(r/rm)**2
	pmms=eps*(an/(am-an))*(rm/r)**am-(am/(am-an))*eps*(rm/r)**an
	return
	end

	subroutine init
        implicit real*8(a-h,o-z)
	common/blk0/iflag
	common/blk1/xxx(12),yyy(12),zzz(12)
	common/blk2/xbarbnd(12),ybarbnd(12),zbarbnd(12)
	common/blk3/xbnd2(12),ybnd2(12),zbnd2(12)
        common/blk4/xbnd1(12),ybnd1(12),zbnd1(12)
        common/blk5/rmperpch,rmparch,epsperpch,epsparch
	common/blk6/rmperpcc,rmparcc,epsperpcc,epsparcc

	iflag=1

c C6H6 geometry
	DATA xxx/1.20377531126d0,1.20377531126d0,0.00d0,-1.20377531126d0,
     <-1.20377531126d0,0.00d0,2.14774300139d0,2.14774300139d0,0.00d0,
     <-2.14774300139d0,-2.14774300139d0,0.00d0/
	DATA yyy/-0.695,0.695,1.390d0,0.695,-0.695,-1.390d0,-1.240d0,
     <1.240d0,2.480d0,1.240d0,-1.240d0,-2.480d0/
	DATA zzz/0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,
     <0.0d0,0.0d0,0.0d0/
       
c MMS parameters for C6H6-He
	
	rmperpcc=3.583
	rmparcc=4.005
	epsperpcc=0.860
	epsparcc=0.881
	rmperpch=3.234
	rmparch=3.480
	epsperpch=1.364
	epsparch=1.016


        do i=1,12
        if(i.le.5) then
        xbarbnd(i)=(xxx(i+1)+xxx(i))*0.5d0
        ybarbnd(i)=(yyy(i+1)+yyy(i))*0.5d0
        zbarbnd(i)=(zzz(i+1)+zzz(i))*0.5d0
        xbnd2(i)=xxx(i+1)-xxx(i)
        ybnd2(i)=yyy(i+1)-yyy(i)
        zbnd2(i)=zzz(i+1)-zzz(i)
        else
        if(i.eq.6) then
          xbarbnd(i)=(xxx(1)+xxx(6))*0.5d0
          ybarbnd(i)=(yyy(1)+yyy(6))*0.5d0
          zbarbnd(i)=(zzz(1)+zzz(6))*0.5d0
          xbnd2(i)=xxx(1)-xxx(6)
          ybnd2(i)=yyy(1)-yyy(6)
          zbnd2(i)=zzz(1)-zzz(6)
          else
           xbarbnd(i)=(xxx(i)+xxx(i-6))*0.5d0
           ybarbnd(i)=(yyy(i)+yyy(i-6))*0.5d0
           zbarbnd(i)=(zzz(i)+zzz(i-6))*0.5d0
           xbnd2(i)=xxx(i)-xxx(i-6)
           ybnd2(i)=yyy(i)-yyy(i-6)
           zbnd2(i)=zzz(i)-zzz(i-6)
           endif
         endif        
        enddo
	return
	end
