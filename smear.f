      SUBROUTINE SMEAR(I,F,CTE)                
      IMPLICIT REAL*8(A-H,O-Z)             
      real*8 P_OLD(4000,5)
      REAL*4 GAUSS,MEAN,SD                 
      COMMON/JETT/P_OLD
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      
      MEAN=P(I,4)
      IF(MEAN.LT.0.1)GO TO 2
      SD=dsqrt(F**2*MEAN + (CTE*MEAN)**2)
      FAC=ABS(GAUSS(MEAN,SD))/MEAN
c     fac=dsqrt(fac**2+cte**2)
      P_OLD(I,4)=P(I,4)
      P_OLD(I,5)=P(I,5)
      DO 1 J=1,3
         P_OLD(I,J)=P(I,J)   
 1       P(I,J)=P(I,J)*FAC
      P(I,4)=DSQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)
   2  RETURN
      END   
 
      
      SUBROUTINE SMEAR_A(I,F)                
      implicit none             
      REAL*8 fac,pi,alpha,dphi,cphi,phi,cteta,teta
      real*8 pj(4000,5),f,sphi,steta,dteta
c      real*8 aux,pjf(8,5),ran,dcteta,ptes(8,5)
      REAL*4 GAUSS,MEAN,SD,ran2                 
      integer iseed,i
      DATA iseed/125763/
      COMMON/JETT/PJ
c      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      
      pi=3.14159d0
      MEAN=0.                         
      SD=F         
c      fac=f             
      FAC=ABS(GAUSS(MEAN,SD))
      alpha=ran2(iseed)
      alpha=alpha*2.d0*pi
      dphi=fac*dcos(alpha) 
      dteta=fac*dsin(alpha)          
      cphi=pj(i,2)/pj(i,5)
      sphi=pj(i,3)/pj(i,5)
      if (sphi.ge.0.d0) then
        phi=dphi+dacos(cphi)
      else     
        if (cphi.ge.0.d0) then
           phi=dphi+dasin(sphi)
        else
           phi=dphi + 2.d0*pi - dacos(cphi)
        endif
      endif

      cteta=pj(i,4)/pj(i,1)
      steta=pj(i,5)/pj(i,1)
      teta=dteta+dacos(cteta)
      pj(i,5)=pj(i,1)*(steta*dcos(dteta)+dsin(dteta)*cteta)
      pj(i,2)=pj(i,5)*(cphi*dcos(dphi)-sphi*dsin(dphi)) 
      pj(i,3)=pj(i,5)*(sphi*dcos(dphi)+dsin(dphi)*cphi)
      pj(i,4)=pj(i,1)*(cteta*dcos(dteta)-steta*dsin(dteta))
c      pjf(i,5)=pj(i,1)*dsin(teta)
c      pjf(i,2)=pj(i,5)*dcos(phi) 
c      pjf(i,3)=pj(i,5)*dsin(phi) 
c      pjf(i,4)=pj(i,1)*dcos(teta)
c      do j=2,5
c      aux= dabs((ptes(i,j)-pjf(i,j))/pjf(i,j))
c      if (aux.gt.1d-4) then
c        write(15,*) 'Warning j,ptes,pj',j,ptes(i,j),pjf(i,j),pj(i,j),aux
c      endif
c      enddo
      return
      END   
      
      FUNCTION GAUSS(FMEAN,SD)
c      real*8 ran
      real*4 ran2
      DATA iseed/-1/
    1 V1=1-2*RAN2(iseed)          
      V2=RAN2(iseed)              
      S=V1**2+V2**2           
      IF(S.GE.1) GO TO 1      
      X=V1*SQRT(-2*ALOG(S)/S) 
      GAUSS=X*SD+FMEAN        
      RETURN                  
      END                     


      FUNCTION RAN2(IDUM)
C
C MACHINE-INDEPENDENT RANDOM NUMBER GENERATOR, GENERAL PURPOSE VERSION
C                     OK AS LONG AS >= 32 BITS
C                      INITIAL IDUM MUST BE < 0
      COMMON/RANDM/RDM(31),RM1,RM2,IA1,IC1,M1,IX1,IA2,IC2,M2,IX2,
     1       IA3,IC3,M3,IX3
      DATA IA1,IC1,M1 /1279,351762,1664557/
      DATA IA2,IC2,M2 /2011,221592,1048583/
      DATA IA3,IC3,M3 /15091,6171,29201/
      IF (IDUM.GE.0) GO TO 2
C INITIALIZATION
      IX1=MOD(-IDUM,M1)
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IX1,M2)
      IX1=MOD(IA1*IX1+IC1,M1)
      IX3=MOD(IX1,M3)
      RM1=1./FLOAT(M1)
      RM2=1./FLOAT(M2)
      DO 1 J=1,31
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
1     RDM(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
C GENERATE NEXT NUMBER IN SEQUENCE
2     IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(31*IX3)/M3
      RAN2=RDM(J)
      RDM(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
C OMIT FOLLOWING STATEMENT IF FUNCTION ARGUMENTS PASSED BY VALUE:
      IDUM=IX1
      RETURN
      END



