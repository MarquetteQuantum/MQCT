      SUBROUTINE USER_DEFINED_PES(V,R,rvib,rvib2,alpha,beta,gamma,
     &										   aalpha,bbeta,ggamma)
!     INPUT:  R - distance betweenn COMs, rvib - vibrational coordinate
! 		      alpha,beta,gamma - Euler's angles of the first molecule
!   		  aalpha,bbeta, ggamma - Euler's angles of the second molecule
!     OUTPUT: V - value of the potential
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 V,R,alpha,beta,gamma,aalpha,bbeta,ggamma,t,rvib,rvib2
      REAL*8 R1,COG
      DIMENSION rbohr(3)
	   call SUPATATNEW(R,1.1007d0,beta,V)
	   V = V*8.06573d0
      
      END SUBROUTINE USER_DEFINED_PES


!----------!
! OPTION 2 !
!----------! 
!	  USE KEYWORD "EXPANSION=YES" TO INITIATE THIS OPTION
      SUBROUTINE USER_DEFINED_TERMS(T,I,R)
!     THIS SUBROTUNE COMPUTES RADIAL COEFFICENTS OF THE PES EXPANSION AT A GIVEN DISTANCE
!     INPUT:  R - distance between COMs of particles, I - TERM NUMBER
!     OUTPUT: T - value of coefficent 	  
      IMPLICIT NONE	  
      REAL*8 T,R
      INTEGER I
!     USER MUST INSERT A CALL OF AN EXTERNAL SUBROUTINE HERE.
!     DELETE THE "STOP" COMMAND BELOW IF THE SUBROUTINE SUPPLIED.
!     IN CASE IF USER FORGOT TO SUPPLY THE SUBTOUITNE,
!     BUT THE MAIN PROGRAM REQUIRES IT, THEN STOP:
      STOP "ERROR: USER_DEFINED_TERMS IS NOT SUPPLIED"
	  END SUBROUTINE USER_DEFINED_TERMS
!----------!
! OPTION 3 !
!----------! 
!	  USE KEYWORDS "EXPANSION=YES, TERMS_FILE=YES" TO INITIATE THIS OPTION
! 	  SIMILAR TO OPTION 2, BUT NO SUBROUTINE IS REQUIRED
!     USER SHOULD PROVIDE THE FILE EXPAN_PES_TERMS.DAT 
!     IN THE MAIN PROGRAM DIRECTORY CONTAINING THE COEFFICEINS 
!     OF POTENTIAL EXPANSION PRECOMPUTED EXTERNALLY.
! 	  SEE EXAMPLE FILES SUPPLIED WITH THE CODE.
!----------!
! OPTION 4 !
!----------! 
!	  USE KEYWORDS "EXPANSION=YES, TERMS_ONFLY=YES" TO INITIATE THIS OPTION
      SUBROUTINE USER_DEFINED_COEFFS(T,DTDR,I,R) 
!     THIS SUBROUTINE COMPUTES RADIAL COEFFICENTS OF THE PES EXPANSION 
!     AND THEIR DERIVATIVES AT A GIVEN DISTANCE R
      IMPLICIT NONE
!     INPUT : R - distance between COMs of particles, I - TERM NUMBER
!     OUTPUT: T - value of coefficent, DTDR - its radial derivative 	  
      REAL*8 T,R,DTDR 
      INTEGER I
!     USER MUST INCERT A CALL OF AN EXTERNAL SUBROUTINE HERE.
!     DELETE THE "STOP" COMMAND BELOW IF THE SUBROUTINE IS SUPPLIED.
!     IN CASE IF USER FORGOT TO SUPPLY THE SUBTOUITNE,
!     BUT THE MAIN PROGRAM REQUIRES IT, THEN STOP:	
      STOP "ERROR: USER_DEFINED_COEFFS IS NOT SUPPLIED"
      END SUBROUTINE USER_DEFINED_COEFFS 

c--------------------------------------------------------------------------
c  FROM HERE IF YOU USE ANOTHER PROGRAM  ----------------
c--------------------------------------------------------------------------


      SUBROUTINE SUPATATNEW(R,rr2,THETARADB,VTOT1)

      IMPLICIT REAL*8 (A-H,O-Z)
c   R distanza intermolecolare in Angstr. (distanza tra O e centro di massa N2)
c   rr2 distanza interna (r_N2) molecola N2 in angstr
c   THETARADB angolo di jacobi tra asse diatomo e distanza intermolecolare in radianti
c   VTOT1 potenziale di interazione in meV   per los stato Pi
c  
c  valore distanza equilibrio r_N2=1.1007 angstr.




      COMMON/PARAM2/epsi1,rm1

      COMMON/PARAM5/dist

      dimension xxx(6),yyy(6),zzz(6),dist(42),q(6)



      parameter(pigreco = 3.1415926535897932d0)

      real*8 nlj,mlj






c quadrupolo in unità atomiche

c quadrupolo O tripleto P
c
c valore ab initio mio
c aug-cc-pV5Z spazio occ,4,2,2,0,4,2,2,0; ACPF

         QA=0.4751d0  ! lungo x e y (lungo z -> -QA*2) 


c carica su O

c---------------------------------------------------------------------
        q(3)=0.d0


c-------------------------------------------------
c  nuova parte elettrostatica per lo stretching
c   secondo Pirani
c-------------------------------------------------
              

c  seconda molecola (N2)


c calcolo MAXBART aug-cc-pV5Z spazio occ,4,2,2,0,4,2,2,0; ACPF
c preso da nostro JCC 2011



c       QB=-1.11495d0
 
       QB=-1.11495d0+1.8359d0*(rr2-1.1007d0)
     #+0.1900d0*(rr2-1.1007d0)**2+0.2496d0*(rr2-1.1007d0)**3
     #-3.6399d0*(rr2-1.1007d0)**4
     #+2.2262d0*(rr2-1.1007d0)**5

c cariche parziali

         qbpos=QB/(2.d0*(rr2/0.529177d0/2.d0)**2)

c         write(*,*) QB,qbpos
c         stop





c valore relativo a distribuzione CILINDRICA
c         qbpos=0.3238d0

         qbneg=-qbpos/2.d0



c coordinate del CM dell'O
        xxx(3)=0.d0
        yyy(3)=0.d0
        zzz(3)=-R/2.d0



c fissa le coordinate del diatomo 2
                rrmed2=rr2/2.d0
                xxx(4)=rrmed2*dsin(thetaradb)
                xxx(5)=-rrmed2*dsin(thetaradb)
                yyy(4)=0.d0
                yyy(5)=0.d0
                zzz(4)=r/2.d0+rrmed2*dcos(thetaradb)
                zzz(5)=r/2.d0-rrmed2*dcos(thetaradb)


c coordinate del CM del N2
        xxx(6)=0.d0
        yyy(6)=0.d0
        zzz(6)=R/2.d0



c cariche su N2

c---------------------------------------------------------------------
        q(4)=qbpos 

        q(5)=qbpos  

        q(6)=4.d0*qbneg  ! CM

c----------------------------------------------------------------





c--------------------------------------------------------------------------------------
c*** Distanze interatomiche tra 2 molecole  per calcolare parte vdW e elettrostatica


           index=0

           do ii=3,3
                 do hh=4,5

                   index=index+1
                     
       dist(index)=dsqrt((xxx(ii)-xxx(hh))**2.d0+(yyy(ii)-yyy(hh))**2.d0
     <    +(zzz(ii)-zzz(hh))**2.d0)
               
                 enddo
           enddo


      CALL EPSIRMATAT(rr2)


c coppia O--N

c                epsi1=7.4d0
c                rm1=3.30d0



c Calcola potenziale (vdW)



                pes1=0.d0


                

c---------------------------------------------------------------
            do j=1,2

               epsi=epsi1
               rm=rm1


            nlj=6
            mlj=8.0+4.0*(dist(j)/rm)**2


            pes1=pes1+epsi*(nlj/(mlj-nlj))*(rm/dist(j))**mlj-
     #(mlj/(mlj-nlj))*epsi*(rm/dist(j))**nlj

                enddo

c---------------------------------------------------------------
 




c***** Calcolo della parte elettrostatica ***********


c-----------------------------------------
c NUOVO
c-----------------------------------------


c          vel=0.d0

c             index=0

c           do ii2=3,3
c                 do hh2=4,6

c                 index=index+1

c           vel=vel+q(ii2)*q(hh2)/dist(index)

c                 enddo
c           enddo


                thetarada=0.d0
                phirad=0.d0



c-------------------------------------------------------------------
c espressione armoniche sferiche

          ANG224=15.d0/dsqrt(630.d0)*
     #(3.d0/8.d0*(dsin(thetarada)**2)
     #*(dsin(thetaradb)**2)*2.d0*p1(2*phirad)-6.d0*p1(thetarada)
     #*dsin(thetarada)*p1(thetaradb)*dsin(thetaradb)*2.d0*p1(phirad)
     #+6.d0*p2(thetarada)*p2(thetaradb))  

c--------------------------------------------------------------------



c        VQUAD=vel*0.529177d0*27212.d0

        VQUAD=QA*QB*dsqrt(14.d0/5.d0)*27212.d0*0.529177d0**5/R**5.d0



            vtot1=
     #pes1
     #+vquad*ANG224


          return
          end



      SUBROUTINE EPSIRMATAT(rrb1)

      IMPLICIT REAL*8 (A-H,O-Z)                                


       COMMON/PARAM2/epsi1,rm1
       COMMON/PARAM3/C6medio1,alfamedioA,alfamedioB1
     #,alfamedioAeff,alfamedioB1eff



           c13=1.d0/3.d0


                ELMOLA=5.25d0   ! per O 

                ELMOLB=4.4d0    ! per N
 


                              
                 amediaA=0.8d0  ! polarizzabilità O tripletto


 
              
                 alfamedioA=amediaA  ! per O


c-----------------------------------------------------------------------
                call  POLARN2(rrB1,aparB,aperB,amediaB,deltaB)                


              
                 alfamedioB1=amediaB*0.5d0



c----------------------------------------------------------
                      









c primo monomero

c O (tripletto P)


          alfamedioAeff=alfamedioA


c seconda molecola


c N2

c          alfapareffB1new=alfapareffB1

c          alfaparB1new=alfaparB1

c          alfaperpeffB1new=alfaperpeffB1

c          alfaperpB1new=alfaperpB1




          alfamedioB1eff=alfamedioB1

c------------------------------------------------------------------


c--------------------------------
C coefficienti best fit


c           FIT1=1.006d0
c           FIT2=1.d0
 
c           FITC61=1.1d0


c valori originali

           FIT1=1.046d0

           FITC6=1.d0

           FITEPS1=1.00896d0
  

c---------------------------------




c Rm1  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      rm1=(alfamedioAeff**(1.d0/3.d0)+alfamedioB1eff**(1.d0/3.d0))/
     #((alfamedioAeff*alfamedioB1eff)**0.095d0)

         rm1=rm1*1.767d0






c correzione per best fit

           rm1=rm1/FIT1





c----------------------------------------------------------------

c C6 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        c6medio1=(alfamedioA*alfamedioB1)/
     #(dsqrt(alfamedioA/ELMOLA)+dsqrt(alfamedioB1/ELMOLB))
      
        c6medio1=c6medio1*15700.d0






c correzione per best fit
        c6medio1=c6medio1/FITC6




c epsi ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


c 1 c-----------------------------------------------------------
         epsi1=0.72d0*c6medio1/rm1**6.d0


         epsi1=epsi1/FITEPS1



    
          return
          end
c-----------------------------------------------------
c========espressioni dei pol. di Legendre===================


       REAL*8 FUNCTION P1(x)
       IMPLICIT REAL*8(A-H,O-Z)
          p1 = dcos(x)
       RETURN
       END

       REAL*8 FUNCTION P2(x)
       IMPLICIT REAL*8(A-H,O-Z)
          p2 = 0.5d0*(3.0d0*dcos(x)*dcos(x) - 1.0d0)
       RETURN
       END

       REAL*8 FUNCTION P3(x)
       IMPLICIT REAL*8(A-H,O-Z)
          p3 = 0.5d0*(5.d0*(dcos(x))**3-3.d0*dcos(x))
       RETURN
       END


      SUBROUTINE POLARN2(rr,apar,aper,amedia,delta)
c programma calcolo polarizzabilita' N2
c in funzione della distanza intermolecolare

      implicit real*8(a-l,o-z)
      dimension abond(70,10)

c--------------------------------------
c ordine di legame - bo
c distanza di legame - bl
c alfa atomo 1 - alf1
c alfa atomo 2 - alf2
c elettroni di valenza atomo 1 - evt1
c elettroni di valenza atomo 2 - evt2
c elettroni di non legame atomo 1 - enb1
c elettroni di non legame atomo 2 - enb2

      bo=3.0d0
      bl=1.1007d0
      alf1=1.100d0
      alf2=1.100d0
      evt1=5.00d0
      evt2=5.00d0
      enb1=2.00d0
      enb2=2.00d0

c--------------------------------------
C COSTANTI
      c13=0.333333333d0
      c12=0.500000000d0
      c16=1.0d0/6.0d0
      c14=1.0d0/4.0d0
      c23=2.0d0/3.0d0
      c32=1.d0/c23
C PARAMETRI
      a1=alf1/evt1
      a2=alf2/evt2
      at=a1+a2
      alfa=alf1**c13+alf2**c13
      a13=a1**c13
      a23=a2**c13
c-------------------------------------
c-------------------------------------
c oppure distanza fissa di equilibrio
c  rr=bl
c-------------------------------------
c fattore di contrazione
         fc=rr/alfa
c fattore di spegnimento
         fs=rr/bl
c fattore novita'
         fn=rr/(a13+a23)
c termini esponenziali
        esp1=2.000d0*fc-1.000d0

c
c        fe1=dexp(-1.000d0*fs*esp1)
c 2o cambio
        fe1=dexp(-1.000d0*(fs-1.d0)*2.d0*fc)


        esp2=0.7500d0*(fs-1.00d0)*(1.00d0-fc)*(fs**c12)*alfa**c12
        fe2=dexp(esp2)
        esp3=2.00d0*(1.00d0-fs)/bo
        fe3=dexp(esp3)

cc
c        esp4=dabs(esp1)
c 1o cambio
        esp4=esp1


        fe4=dexp(-0.75d0*esp4*(fs**c12)*alfa**c12)
c---------------------------------------
c alfa media legame
c---------------------------------------

c fattore campana
        cex=(fs-1.0d0)*fe2/3.0d0
        camp=1.0d0+cex
c ordine di legame effettivo

c        b2=fs+1.00d0
c        b2=1.00d0/b2
c        b1=dabs(bo-2.0d0)
c        b1=3.0d0*(bo**b1)/4.0d0
c 3o cambio

        b2=4.d0*(fc**2)+1.00d0
        b2=(fs+2.00d0)/b2
        b1=((bo-2.d0)**2)-1.d0/8.d0
        b1=1.0d0*(bo**b1)/4.0d0


        bof=bo-b1*b2*fe1

c---------------------------------------
c polarizzabilita media (amedia)

        alfnb=a1*enb1+a2*enb2
        alfv=at*bof
        amedia=(alfv+alfnb)*camp
c---------------------------------------
c anisotropia (delta)
c---------------------------------------
        pre1=c32*alfa**c12*fn**2.d0
        pre21=fe3*(enb1+enb2)/2.0d0+(enb1*enb2)**c12
        pre21=(pre21/(evt1+evt2))
        pre24=bo**((bo-4.d0)/(bo+1.0d0))
        pre2=(1.0d0+pre24*pre21)**(3.d0)
        pre4=2.0d0*fc/(bo**c13)-dabs(alf1-alf2)/(alf1+alf2)
        delta=fe4*amedia*pre1*pre4/pre2
c---------------------------------------
c parallela--perpendicolare
c---------------------------------------
        apar=amedia+c23*delta
        aper=amedia-c13*delta
c---------------------------------------
 1    continue
      return
      end
