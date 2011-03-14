* **  *****************************************************************
* **  "diff_sub_23102010": is the subroutine differences with MBMs version with 
* **  svn.r57 version. Please check for the differences. However few sugggestions:
* ** 1. integer kfstdaug(1000),lnstdaug(1000) :  kfstdaug(4000),lnstdaug(4000)
* ** 2. double precision pstdaugm(5,1000),pxprng1(3),pxprng2(3): 4000
* ** 3. character*16 namestdaug(1000) : 4000 
* **
* **
* **  *****************************************************************

      subroutine chisignalvertex(jj,chifin,chiline,kfchifin,nfchi,
     $     show,enocut,munocut,taunocut,
     $     pfullmu,pfulltau,pfullp1,pfullp2)
C      nmuevt in commons
      implicit double precision(a-h,o-z)
      include 'pylocalcom.h' 
      EXTERNAL PYDATA
*******************-- local variables --********************************
      double precision pwdec(5,5,10)
      integer orig(5)
      integer chiline(5),ki1w(5),mworig(5)
      integer kfwdec(5,10)
      character*16 nwamendec(5,10)
      logical signalvertex(5)
      logical show
      integer chifin(5),kfchifin(5,10)
      integer nw2jj(10),nwmu(5),nwtau(5),nfchi(5,10)
      logical enocut(10),munocut(10),taunocut(10)
      double precision  pfullmu(10,8),pfulltau(10,8),pfulljet1(10,8)
      double precision pfulljet2(10,8),pfullp1(10,8),pfullp2(10,8)
***************__ Commons __**************
      double precision lum
      common/sig/sig_sm,lum
      integer nevgen,nevgenwmu,nevgenwtau,nevgenwe ! spdas not consistent with other common
      common/porodskands/nevgen,nevgenwmu,nevgenwtau,nevgenwe ! spdas not consistent with other common
      integer nffin(5,100)
      common/dauglines/nffin
      common/muevt/nmuevt
      common/tauevt/ntauevt ! DR
      common/eventnumber/iev,nchilinejj
C... INPUT: chifin(jj),chiline(jj),kfchifin(jj,ii),show,
      signalvertex(jj)=.False.
      nw2jj(jj)=0
      do 202 ii=1,chifin(jj)    ! loop over neutralino final states
c     write(*,*)'chi_10->',kfchifin(jj,ii),ii,iev
         
C......(II.a) CHECK IF W-> qq'
         wevt: if(abs(kfchifin(jj,ii)).eq.24)then ! W id number is 24
            orig(1)=chiline(jj)
            kfwpm=abs(kfchifin(jj,ii))
            call decaychain(orig,1,kfwpm,
     $           iiw,ki1w,mworig,nwamendec,kfwdec,pwdec)
C: common output: nffin(iiw,jjj) see below for jjj meaning
            
            if(show)call writedecay(orig,1,kfwpm,
     $           iiw,ki1w,mworig,nwamendec,kfwdec,pwdec)
            
c.... Check what are the W daughters
            do jjj=1,mworig(iiw)
               if(abs(kfwdec(iiw,jjj)).le.6)then ! id number <=6 for quarks!
C.... Check that the reconstructed vertex satisfy the cuts on jets:
C     jet cuts on W daughters, nffin is in the commons of decaychain
c     if(pyp(nffin(iiw,jjj),10).gt.15..and.
c     $                          pyp(nffin(iiw,jjj),19).lt.4.9)
c     $                
                  nw2jj(jj)=nw2jj(jj)+1
c                  do jjjj=1,2
c                     ppartonw(jjjj,jjj)=pyp(nffin(iiw,jjj),11+jjjj*4)
c                  end do
c                  if(show)print*,'phi,eta:',ppartonw(1,1),ppartonw(2,1)

               end if
            end do
C...  comment out for all possible W decays
            if(nw2jj(jj).eq.2)then
C...  uncomment out for all possible W decays
c     if(nw2jj(jj).ge.0)
                signalvertex(jj)=.True.
                  do jjjj=1,4
                     pfullp1(jj,jjjj)=pyp(nffin(iiw,1),jjjj)
                     pfullp2(jj,jjjj)=pyp(nffin(iiw,2),jjjj)
                  end do
                  pfullp1(jj,5)=pyp(nffin(iiw,1),10)
                  pfullp1(jj,6)=pyp(nffin(iiw,1),19)
                  pfullp1(jj,7)=1
                  pfullp1(jj,8)=pyp(nffin(iiw,1),15)
                  pfullp2(jj,5)=pyp(nffin(iiw,2),10)
                  pfullp2(jj,6)=pyp(nffin(iiw,2),19)
                  pfullp2(jj,7)=1
                  pfullp2(jj,8)=pyp(nffin(iiw,2),15)
             end if
         end if wevt
 202  CONTINUE
      
C.... (II.b) CHECK IF ~chi_10 -> (W -> qq') tau, or (W -> qq') mu  
C.... Check that the W companion is a proper muon or tau
      if(signalvertex(jj))then 
         do 203 ii=1,chifin(jj) ! loop over neutralino final states
            if(abs(kfchifin(jj,ii)).eq.11)then ! e id # is 
C: DEBUG: change to .True.
               signalvertex(jj)=.false. 
               enocut(jj)=.True.
            elseif(abs(kfchifin(jj,ii)).eq.13)then ! mu id # is 13
               nmuevt=nmuevt+1
               munocut(jj)=.True.
               do jjj=1,4
                  pfullmu(jj,jjj)=pyp(nfchi(jj,ii),jjj)
               end do
               pfullmu(jj,5)=pyp(nfchi(jj,ii),10)
               pfullmu(jj,6)=pyp(nfchi(jj,ii),19)
               pfullmu(jj,8)=pyp(nfchi(jj,ii),15)
               pfullmu(jj,7)=1
C...  Cuts for proper leptons
C...  BEGIN comment out for pure chi_10 -> mu W -> mu q q' events
c               qrand=rand(iev+jj)
               if(pyp(nfchi(jj,ii),10).lt.20.or.
     $              abs(pyp(nfchi(jj,ii),19)).gt.2.5)then
!     $              .or.qrand.lt.0.05)then ! 95% efficiency in muons
C     .DEBUG
                  signalvertex(jj)=.false.
               else
C...  END comment out for pure chi_10 -> mu W -> mu q q' events
                  nwmu(jj)=nwmu(jj)+1
C...  comment the end if for pure chi_10 -> mu W -> mu q q' events
               end if
            elseif(abs(kfchifin(jj,ii)).eq.15)then ! tau id # is 15
               ntauevt=ntauevt+1 ! DR
C...  comment as before for pure chi_10 -> tau W -> tau q q' events
               taunocut(jj)=.True.
               do jjj=1,4
                  pfulltau(jj,jjj)=pyp(nfchi(jj,ii),jjj)
               end do
               pfulltau(jj,5)=pyp(nfchi(jj,ii),10)
               pfulltau(jj,6)=pyp(nfchi(jj,ii),19)
               pfulltau(jj,8)=pyp(nfchi(jj,ii),15)
               pfulltau(jj,7)=1
               if(pyp(nfchi(jj,ii),10).lt.20.or.
     $              abs(pyp(nfchi(jj,ii),19)).gt.2.5)then ! BR(tau -> 3-prong) * 85% efficiency in taus =0.15*0.85
C...........for one-prong BR(tau -> 1-prong)=85% and efficiency is 50% 
                  signalvertex(jj)=.false.
               else
                  nwtau(jj)=nwtau(jj)+1
               end if
c               write(*,*)''
            elseif(abs(kfchifin(jj,ii)).eq.24)then
C...  just the W
            else
               write(*,*)'ERROR wrong W companion',
     $              kfchifin(jj,ii)
               stop
            end if
 203     CONTINUE
      end if
      end
***********************************************************************
      subroutine chisignalvertex2(jj,chifin,chiline,kfchifin,nfchi,
     $     show,enocut,munocut,taunocut,pfulle,
     $     pfullmu,pfulltau,pfullp1,pfullp2)    ! pfullp3 is not passing: chnaged by spdas

C      nmuevt in commons
      implicit double precision(a-h,o-z)
      include 'pylocalcom.h' 
      EXTERNAL PYDATA
*******************-- local variables --********************************
      double precision pwdec(5,5,10)
      integer orig(5)
      integer chiline(5),ki1w(5),mworig(5)
      integer kfwdec(5,10)
      character*16 nwamendec(5,10)
      logical signalvertex(5)
      logical show
      integer chifin(5),kfchifin(5,10)
      integer nw2jj(10),nwe(5),nwmu(5),nwtau(5),nfchi(5,10)
      logical enocut(10),munocut(10),taunocut(10)
      double precision pfulle(10,8),pfullp3(10,8)
      double precision  pfullmu(10,8),pfulltau(10,8),pfulljet1(10,8)
      double precision pfulljet2(10,8),pfullp1(10,8),pfullp2(10,8)
***************__ Commons __**************
      double precision lum
      common/sig/sig_sm,lum
      integer nevgen,nevgenwmu,nevgenwtau,nevgenwe
      common/porodskands/nevgen,nevgenwmu,nevgenwtau,nevgenwe
      integer nffin(5,100)
      common/dauglines/nffin
      common/muevt/nmuevt
      common/eevt/neevt
      common/tauevt/ntauevt ! DR
      common/eventnumber/iev,nchilinejj
C... INPUT: chifin(jj),chiline(jj),kfchifin(jj,ii),show,
      signalvertex(jj)=.False.
      nw2jj(jj)=0
      do 202 ii=1,chifin(jj)    ! loop over neutralino final states
c     write(*,*)'chi_10->',kfchifin(jj,ii),ii,iev
         
C......(II.a) CHECK IF W-> qq'
         wevt: if(abs(kfchifin(jj,ii)).eq.24)then ! W id number is 24
            orig(1)=chiline(jj)
            kfwpm=abs(kfchifin(jj,ii))
            call decaychain(orig,1,kfwpm,
     $           iiw,ki1w,mworig,nwamendec,kfwdec,pwdec)
C: common output: nffin(iiw,jjj) see below for jjj meaning
            
            if(show)call writedecay(orig,1,kfwpm,
     $           iiw,ki1w,mworig,nwamendec,kfwdec,pwdec)
            
c.... Check what are the W daughters
            do jjj=1,mworig(iiw)
               if(abs(kfwdec(iiw,jjj)).le.6)then ! id number <=6 for quarks!
C.... Check that the reconstructed vertex satisfy the cuts on jets:
C     jet cuts on W daughters, nffin is in the commons of decaychain
c     if(pyp(nffin(iiw,jjj),10).gt.15..and.
c     $                          pyp(nffin(iiw,jjj),19).lt.4.9)
c     $                
                  nw2jj(jj)=nw2jj(jj)+1
c                  do jjjj=1,2
c                     ppartonw(jjjj,jjj)=pyp(nffin(iiw,jjj),11+jjjj*4)
c                  end do
c                  if(show)print*,'phi,eta:',ppartonw(1,1),ppartonw(2,1)

               end if
            end do
C...  comment out for all possible W decays
            if(nw2jj(jj).eq.2)then
C...  uncomment out for all possible W decays
c     if(nw2jj(jj).ge.0)
                signalvertex(jj)=.True.
                  do jjjj=1,4
                     pfullp1(jj,jjjj)=pyp(nffin(iiw,1),jjjj)
                     pfullp2(jj,jjjj)=pyp(nffin(iiw,2),jjjj)
                  end do
                  pfullp1(jj,5)=pyp(nffin(iiw,1),10)
                  pfullp1(jj,6)=pyp(nffin(iiw,1),19)
                  pfullp1(jj,7)=1
                  pfullp1(jj,8)=pyp(nffin(iiw,1),15)
                  pfullp2(jj,5)=pyp(nffin(iiw,2),10)
                  pfullp2(jj,6)=pyp(nffin(iiw,2),19)
                  pfullp2(jj,7)=1
                  pfullp2(jj,8)=pyp(nffin(iiw,2),15)
             end if
         end if wevt
 202  CONTINUE
      
C.... CHECK IF ~chi_10 -> (W -> qq') tau, or (W -> qq') mu, or (W -> qq') e
C.... Check that the W companion is a proper e or muon or tau
      if(signalvertex(jj))then 
         do 203 ii=1,chifin(jj) ! loop over neutralino final states
            if(abs(kfchifin(jj,ii)).eq.11)then ! e id # is 
C: DEBUG: change to .True.
               neevt=neevt+1
               enocut(jj)=.True.
                do jjj=1,4
                  pfulle(jj,jjj)=pyp(nfchi(jj,ii),jjj)
               end do
               pfulle(jj,5)=pyp(nfchi(jj,ii),10)
               pfulle(jj,6)=pyp(nfchi(jj,ii),19)
               pfulle(jj,8)=pyp(nfchi(jj,ii),15)
               pfulle(jj,7)=1
C...  Cuts for proper leptons
C...  BEGIN comment out for pure chi_10 -> mu W -> mu q q' events
c               qrand=rand(iev+jj)
               if(pyp(nfchi(jj,ii),10).lt.20.or.
     $              abs(pyp(nfchi(jj,ii),19)).gt.2.5)then
!     $              .or.qrand.lt.0.05)then ! 95% efficiency in muons
C     .DEBUG
                  signalvertex(jj)=.false.
               else
C...  END comment out for pure chi_10 -> mu W -> mu q q' events
                  nwe(jj)=nwe(jj)+1
C...  comment the end if for pure chi_10 -> mu W -> mu q q' events
               end if
             elseif(abs(kfchifin(jj,ii)).eq.13)then ! mu id # is 13
               nmuevt=nmuevt+1
               munocut(jj)=.True.
               do jjj=1,4
                  pfullmu(jj,jjj)=pyp(nfchi(jj,ii),jjj)
               end do
               pfullmu(jj,5)=pyp(nfchi(jj,ii),10)
               pfullmu(jj,6)=pyp(nfchi(jj,ii),19)
               pfullmu(jj,8)=pyp(nfchi(jj,ii),15)
               pfullmu(jj,7)=1
C...  Cuts for proper leptons
C...  BEGIN comment out for pure chi_10 -> mu W -> mu q q' events
c               qrand=rand(iev+jj)
               if(pyp(nfchi(jj,ii),10).lt.20.or.
     $              abs(pyp(nfchi(jj,ii),19)).gt.2.5)then
!     $              .or.qrand.lt.0.05)then ! 95% efficiency in muons
C     .DEBUG
                  signalvertex(jj)=.false.
               else
C...  END comment out for pure chi_10 -> mu W -> mu q q' events
                  nwmu(jj)=nwmu(jj)+1
C...  comment the end if for pure chi_10 -> mu W -> mu q q' events
               end if
            elseif(abs(kfchifin(jj,ii)).eq.15)then ! tau id # is 15
               ntauevt=ntauevt+1 ! DR
C...  comment as before for pure chi_10 -> tau W -> tau q q' events
               taunocut(jj)=.True.
               do jjj=1,4
                  pfulltau(jj,jjj)=pyp(nfchi(jj,ii),jjj)
               end do
               pfulltau(jj,5)=pyp(nfchi(jj,ii),10)
               pfulltau(jj,6)=pyp(nfchi(jj,ii),19)
               pfulltau(jj,8)=pyp(nfchi(jj,ii),15)
               pfulltau(jj,7)=1
               if(pyp(nfchi(jj,ii),10).lt.20.or.
     $              abs(pyp(nfchi(jj,ii),19)).gt.2.5)then ! BR(tau -> 3-prong) * 85% efficiency in taus =0.15*0.85
C...........for one-prong BR(tau -> 1-prong)=85% and efficiency is 50% 
                  signalvertex(jj)=.false.
               else
                  nwtau(jj)=nwtau(jj)+1
               end if
c               write(*,*)''
            elseif(abs(kfchifin(jj,ii)).eq.24)then
C...  just the W
            else
               write(*,*)'ERROR wrong W companion',
     $              kfchifin(jj,ii)
               stop
            end if
 203     CONTINUE
      end if
      end
***********************************************************************
******    ROUTINE TO CHECKING 3-BODY DECAYS   *************************
***********************************************************************      
      subroutine chisignalvertex3(jj,chifin,chiline,kfchifin,nfchi,
     $     show,enocut,munocut,taunocut,pfulle,
     $     pfullmu,pfulltau,pfulltau2,pfullp1,pfullp2)    ! pfullp3 is not passing: chnaged by spdas

C      nmuevt in commons
      implicit double precision(a-h,o-z)
      include 'pylocalcom.h' 
      EXTERNAL PYDATA
*******************-- local variables --********************************
      double precision pwdec(5,5,10)
      integer orig(5)
      integer chiline(5),ki1w(5),mworig(5)
      integer kfwdec(5,10)
      character*16 nwamendec(5,10)
      logical signalvertex(5)
      logical show
      integer chifin(5),kfchifin(5,10)
      integer nw2jj(10),nwe(5),nwmu(5),nwtau(5),nfchi(5,10)
      logical enocut(10),munocut(10),taunocut(10)
      double precision pfulle(10,8),pfullp3(10,8)
      double precision  pfullmu(10,8),pfulltau(10,8),pfulljet1(10,8)
      double precision pfulljet2(10,8),pfullp1(10,8),pfullp2(10,8)
      double precision pfulltau2(10,8)
***************__ Commons __**************
      double precision lum,tautag
      common/sig/sig_sm,lum
      integer nevgen,nevgenwmu,nevgenwtau,nevgenwe
      common/porodskands/nevgen,nevgenwmu,nevgenwtau,nevgenwe
      integer nffin(5,100)
      common/dauglines/nffin
      common/muevt/nmuevt
      common/eevt/neevt
      common/tauevt/ntauevt ! DR
      common/jetevt/njetevt !MBM
      common/eventnumber/iev,nchilinejj
C... INPUT: chifin(jj),chiline(jj),kfchifin(jj,ii),show,
      signalvertex(jj)=.False.
      nw2jj(jj)=0
      do 202 ii=1,chifin(jj)    ! loop over neutralino final states
c     write(*,*)'chi_10->',kfchifin(jj,ii),ii,iev
         
         wevt: if(abs(kfchifin(jj,ii)).eq.15)then ! tau id number is 15
            signalvertex(jj)=.True.
            do jjj=1,4
               pfulltau(jj,jjj)=pyp(nfchi(jj,ii),jjj)
            end do
            pfulltau(jj,5)=pyp(nfchi(jj,ii),10)
            pfulltau(jj,6)=pyp(nfchi(jj,ii),19)
            pfulltau(jj,8)=pyp(nfchi(jj,ii),15)
            pfulltau(jj,7)=1
            tautag = pyp(nfchi(jj,ii),19)
         end if wevt
 202  CONTINUE
      
C.... CHECK IF ~chi_10 -> tau + l + nu, and
C.... Check that the tau companion is a proper e or muon 
      if(signalvertex(jj))then 
         do 203 ii=1,chifin(jj) ! loop over neutralino final states
            if(abs(kfchifin(jj,ii)).eq.11)then ! e id # is 11 
C: DEBUG: change to .True.
               neevt=neevt+1
               enocut(jj)=.True.
                do jjj=1,4
                  pfulle(jj,jjj)=pyp(nfchi(jj,ii),jjj)
               end do
               pfulle(jj,5)=pyp(nfchi(jj,ii),10)
               pfulle(jj,6)=pyp(nfchi(jj,ii),19)
               pfulle(jj,8)=pyp(nfchi(jj,ii),15)
               pfulle(jj,7)=1
C...  Cuts for proper leptons
               if(pyp(nfchi(jj,ii),10).lt.20.or.
     $              abs(pyp(nfchi(jj,ii),19)).gt.2.5)then
                  signalvertex(jj)=.false.
               else
                  nwe(jj)=nwe(jj)+1
               end if
            elseif(abs(kfchifin(jj,ii)).eq.13)then ! mu id # is 13
               nmuevt=nmuevt+1
               munocut(jj)=.True.
               do jjj=1,4
                  pfullmu(jj,jjj)=pyp(nfchi(jj,ii),jjj)
               end do
               pfullmu(jj,5)=pyp(nfchi(jj,ii),10)
               pfullmu(jj,6)=pyp(nfchi(jj,ii),19)
               pfullmu(jj,8)=pyp(nfchi(jj,ii),15)
               pfullmu(jj,7)=1
C...  Cuts for proper leptons
               if(pyp(nfchi(jj,ii),10).lt.20.or.
     $              abs(pyp(nfchi(jj,ii),19)).gt.2.5)then
!     $              .or.qrand.lt.0.05)then ! 95% efficiency in muons
                  signalvertex(jj)=.false.
               else
                  nwmu(jj)=nwmu(jj)+1
               end if
            elseif(abs(kfchifin(jj,ii)).eq.15)then ! tau id # is 15
               if(tautag.eq.pyp(nfchi(jj,ii),19))then
                  ntauevt=ntauevt   ! picking same tau
               else
                  ntauevt=ntauevt+1 ! a different tau
                  taunocut(jj)=.True.
                  do jjj=1,4
                     pfulltau2(jj,jjj)=pyp(nfchi(jj,ii),jjj)
                  end do
                  pfulltau2(jj,5)=pyp(nfchi(jj,ii),10)
                  pfulltau2(jj,6)=pyp(nfchi(jj,ii),19)
                  pfulltau2(jj,8)=pyp(nfchi(jj,ii),15)
                  pfulltau2(jj,7)=1
c                  if(pyp(nfchi(jj,ii),10).lt.20.or.
c     $              abs(pyp(nfchi(jj,ii),19)).gt.2.5)then ! BR(tau -> 3-prong) * 85% efficiency in taus =0.15*0.85
C...........for one-prong BR(tau -> 1-prong)=85% and efficiency is 50% 
c                  signalvertex(jj)=.false.
C               else
C                  nwtau(jj)=nwtau(jj)+1
C               end if
c               write(*,*)''
               endif
            elseif(abs(kfchifin(jj,ii)).eq.12.or.abs(kfchifin(jj,ii)).eq
     $              .14.or.abs(kfchifin(jj,ii)).eq.16)then  ! neutrinos are ok
            else
c...  Chi -> tau + jets  3-body decays
               njetevt=njetevt+1
            end if
 203     CONTINUE
      end if
      end
*******************************************************************
      subroutine writedecay(orig,iiorig,kfdecay,iidecay,kdecay,fin,
     $     finname,kffin,pfin)
********************************************************************
*  Search for the decays of the particle with KF=kfdecay. This
* unstable particle is the daughter of any of iiorig-th particle occupying 
* the position i=orig(ii) (i=1,n and ii=1,iiorig) in the event history 
* part (k(i,1)=21) of the pythia table. If orig(1)=-1, the parent of the
* unstable particle is no taken into account 
*  INPUT
*  orig(1) ... orig(iiorig) parents of instable particle wich have:
*  kfdecay as KF code
*  OUTPUT
*  iidecay: number of instable particles
*  kdecay(1) ... kdecay(iidecay) line number in in the event history
*              part (k(i,1)=21) of the pythia table for each of the 
*              inestable 1...iidecay particles
*  fin(1) ... fin(iidecay) number of final states for each of the
*                          inestable 1...iidecay particles
*  finname(1,fin(1))...finname(iidecay,fin(iidecay)): Names of the
*              final states for each of the inestable 
*              1...iidecay particles
*  kffin(1,fin(1)),...kffin(iidecay,finn(iidecay)) KF codes of the final 
*  s         tates for each of the inestable 1...iidecay particles
*  pfin(i,iidecay,fin(iidecay)) momentum i for each inestable particle
*        iidecay, into  
*  OUTPUT FROM THE COMMONS:
*  nffin(iidecay,fin(iidecay)) table position of final states for each of
*                              the inestable 1...iidecay particles
**********************************************************
      Implicit double precision(a-h, o-z)
      integer orig(5),kdecay(5),fin(5),kffin(5,*)
      character*16 finname(5,*)
      character*8 ckfchi,cchifin,cchiname,ckfchifin,cpchifin
      character*8 cchiline
      character*5 ciichi
      character*7 prtcle
      double precision pfin(5,5,*)
**************--EXTERNAL VARIABLES--************************************
C...  Program variables
      integer nffin(5,100)
      common/dauglines/nffin
      common/eventnumber/iev,nchilinejj
C...  Pythia variables
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
**********************************************************************
      call pyname(kfdecay,prtcle)
      if(kfdecay.eq.1000022)then
         ckfchi='kfchi'
         ciichi='iichi'
         cchiline='chiline'
         cchifin='chifin'
         cchiname='chiname'
         ckfchifin='kfchifin'
         cpchifin='pchifin'
      else
         ckfchi='kfwpm'
         ciichi='iiw'
         cchiline='ki1w'
         cchifin='mworig'
         cchiname='nwamende'
         ckfchifin='kfwdec'
         cpchifin='pwdec'
      end if

      if(kfdecay.eq.1000022)then
         write(*,*)'****EVENT NUMBER',iev,' BEGIN HERE*************'
         write(*,*)'BEGIN NEUTRALINO DECAY'
         write(*,*)'See decay table later'
         write(*,*)'INPUT:'

         write(*,*)' chiorig: not used for cascade initiator:',
     $        orig(1)
      else
         write(*,*)'*BEGIN detailed analysis of W decay from chi_10:',
     $        nchilinejj
         write(*,*)'INPUT:'
         write(*,*)' chiline(jj): the jj neutralino is the parent',
     $        orig(1)

      end if
      write(*,*)' iichiorig: not used for cascade initiator:',
     $     iiorig
      write(*,*)ckfchi,': KF code of ',prtcle,':',kfdecay
      write(*,*)'OUTPUT'
      write(*,*)ciichi,': number of ',prtcle,' found in the table:'
     $     ,iidecay
      write(*,*)cchiline,'(',ciichi,')',': line number in the table of'
      write(*,*)'   each ',prtcle,':',(kdecay(jj),jj=1,iidecay) 
      write(*,*)'   you can crosscheck with the next decay table'
      write(*,*)cchifin,'(',ciichi,'): number of final states of each'
      write(*,*)'   each',prtcle,(fin(jj),jj=1,iidecay)
      write(*,*)cchiname,'(',ciichi,',',cchifin,'(',ciichi,
     $     ')): names of final'
      write(*,*)'   states of each ',prtcle,':'
      write(*,*)ciichi,'=1: ',(finname(1,jj),jj=1,fin(1))
      write(*,*)ciichi,'=2: ',(finname(2,jj),jj=1,fin(2))
      write(*,*)ckfchifin,'(',ciichi,',',cchifin,'(',ciichi,
     $     ')): kf code of final'
c      write(*,*)' kfchifin(iichi,chifin(iichi)): kf code of final'
      write(*,*)'   states of each ',prtcle,':'
      write(*,*)ciichi,'=1: ',(kffin(1,jj),jj=1,fin(1))
      write(*,*)ciichi,'=2: ',(kffin(2,jj),jj=1,fin(2))
      write(*,*)' nffin(iichi,chifin(iichi)): line number in table'
      write(*,*)'    of final states of each ',prtcle,':'
      write(*,*)'   iichi=1: ',(nffin(1,jj),jj=1,fin(1))
      write(*,*)'   iichi=2: ',(nffin(2,jj),jj=1,fin(2))
      write(*,*)cpchifin,'(i,',ciichi,',',cchifin,'(',ciichi,
     $     ')): p_i of final'
      write(*,*)'   states of each ',prtcle,':'
      write(*,*)'   i=1,',ciichi,'=1:',(pfin(1,1,jj),jj=1,fin(1))
      write(*,*)'   i=1,',ciichi,'=2:',(pfin(1,2,jj),jj=1,fin(2))
      write(*,*)'Note that p(nffin(1,1),i): i momentum of the first'
      write(*,*)'daughter of the first ',prtcle,', is equal to '
      write(*,*)cpchifin,'(i,1,1). In fact:'
      write(*,*)'(',(p(nffin(1,1),jj),jj=1,3),")=",
     $     '(',(pfin(jj,1,1),jj=1,3),')'
      
      end
**********************************************************
C...  INPUT:
c...   jj: jj-th neutralino
c...   chifin(jj): number of neutralino final states.
c...   kfchifin(jj,ii): kf number of ii neutralino daughter
c...   nfchi(jj,ii): line number of ii neutralino daughter
c...   vchi(jj,1-4): displaced vertex position of neutralino
c...   RRsili: transverse pixel resolution
c...   BBsili: longitudinal pixel resolution
c...   show: be verbose if TRUE
C...  OUTPUT
      subroutine tauanalysis(jj,chifin,kfchifin,nfchi,vchi,
     $     RRsili,BBsili,show,show_tau_event,efficiency,
     $     sphere_line,tau_event,nprong,nprongh,nprongmu,
     $     pfullprng1,pfullprng2,pfullprng3,vtau)
C      nmuevt in commons

      implicit double precision(a-h,o-z)
      include 'pylocalcom.h' 
      EXTERNAL PYDATA
*******************-- local variables --********************************
      double precision ptdec(5,5,100)
      double precision vchi(10,4)
      double precision distvtx(10),dprng(10)
      double precision p1(3),p2(3)
      double precision v1(3),v2(3)
      double precision pprong(3)
      integer ki1t(5),mtorig(5)
      integer kftdec(5,10)
      character*16 ntamendec(5,10)
      logical show_tau_event,show
      logical tau_event,sphere_line(10)
      logical efficiency(5)
      integer chifin(5),kfchifin(5,10)
      integer nfchi(5,10)
c$$$      integer nwfin(5,10)
      integer kfstdaug(1000),lnstdaug(1000)
      double precision pstdaugm(5,1000),pxprng1(3),pxprng2(3)
      logical allfst,ptcut(5),rcut(5),etcut(5)
      character*16 namestdaug(1000)
      character*8 allchrgd
      integer nprong(10),nprongh(10),nprongmu(10)
      double precision v1m(3),v2m(3),p2m(3),pxprng3(3),pcone(3)
      double precision pfullprng1(10,8),pfullprng2(10,8)
      double precision pfullprng3(10,8),vtau(10,4)
***************__ Commons __**************
      double precision lum
      common/sig/sig_sm,lum
      integer nevgen,nevgenwmu,nevgenwtau,nevgenwe ! spdas not consistent with other common
      common/porodskands/nevgen,nevgenwmu,nevgenwtau,nevgenwe ! spdas not consistent with other common
      integer nffin(5,100)
      common/dauglines/nffin
      common/muevt/nmuevt
      common/eventnumber/iev,nchilinejj
      common/teste/ievv,ptcut,rcut,etcut     !MBM
      distvtx(jj)=-1d0
      efficiency(jj)=.False.
      if(show)write(*,*)'tau decay'

      do 206 ii=1,chifin(jj)    ! loop over neutralino final states
         if(show)write(*,*)'~chi_10->',kfchifin(jj,ii),ii,iev
C.....(VI) IDENTIFYING tau DAUGHTERS
C.....Storing tau momentum

         tauanalysisif: if(abs(kfchifin(jj,ii)).eq.15)then !.and.nwtau(jj).gt.0)then ! tau id number is 15
            tau_event=.True.
            if(show)call pylist(3)
***********************************************************************
c...  Find what are the tau daughters
*=====================================================================
C.... (A) By using directly fulldecaychain with full fragmented table
            if(show)write(*,*)'tau line number:',nfchi(jj,ii),
     $           ', in previous fragmented table.'
C...  INPUT:
C...  nfchi(jj,ii): line number in table of tau coming from 
C...  ~chi_10 -> tau W
C..   *********************************************
C...  allfst: .False. Only charged final stable states (FSS)
C...  .True.  All FSS
            allfst=.True. 
C,,,  *********************************************
C...  OUTPUT:
C...  mtorig(1); number of final stable states (FSS) of tau
C...  kfstdaug(1)...kfstdaug(mtorig(1)): kf of FSS 
C...  lnstdaug(1)...lnstdaug(mtorig(1)): line number in table of FSS
C...  namestdaug(1)...namestdaug(mtorig(1)): names of FSS 
C...  pstdaug(1...5,1)...pstdaug(1...5,mtorig(1)): full momentum of FSS
            newiline=nfchi(jj,ii) !nfchi(jj,ii) ! 100
            if(show)write(*,*)'Finding the daughters of kf',
     $           k(newiline,2),' in line number',newiline
            iit=1               ! number of taus to be analysed 
            ki1t(iit)=newiline  ! line number in the table for tau
            call fulldecaychain(newiline,allfst,mtorig(iit),
     $           kfstdaug,lnstdaug,namestdaug,pstdaugm)
            ki1t(iit)=newiline
C...  Set tau daughters output variables:

            do ikik=1,mtorig(iit)
               kftdec(iit,ikik)=kfstdaug(ikik)
               nffin(iit,ikik)=lnstdaug(ikik)
               ntamendec(iit,ikik)=namestdaug(ikik)
               do ikiki=1,5 
                  ptdec(ikiki,iit,ikik)=pstdaugm(ikiki,ikik)
               enddo
            enddo
*=====================================================================
c...  (B) By using decaychainnew:
c     call decaychainnew(chiline(jj),1,abs(kfchifin(jj,ii)), ! used after pyexec
c     $                    iit,ki1t,mtorig,ntamendec,kftdec,ptdec)
*=====================================================================
            if(allfst)then
               allchrgd='  all'
            else
               allchrgd='charged'
            end if
            if(show)then
               write(*,*)'Total number of ',allchrgd,
     $              'stable daugthers:',mtorig(iit),' of line',
     $              ki1t(iit)
               write(*,*)'kf of daug',
     $              (kftdec(iit,jjjj),jjjj=1,mtorig(iit))
               write(*,*)'line number of daug',
     $              (nffin(iit,jjjj),jjjj=1,mtorig(iit))
               write(*,*)'name of daug: ',
     $              (ntamendec(iit,jjjj),jjjj=1,mtorig(iit))
               write(*,*)'p_1 of daug ',
     $              (ptdec(1,iit,jjjj),jjjj=1,mtorig(iit))
            endif
***********************************************************************
c.... Check what are the tau daughters
            ptcut(jj) = .False.  !MBM
            rcut(jj) = .False.  !MBM
            iprong=0
            nprongh(jj)=0
            nprongmu(jj)=0

* ** spdas: mtorig: number of the final statble states of taus.
* **  The below block fill up the "jj^th Neutralinos daughters LSP->\tauW"
* ** 
* ** 

            do jjj=1,mtorig(iit)
               if(abs(kchg(pycomp(kftdec(1,jjj)),1)).gt.0)then
                  iprong=iprong+1
                  if(iprong.eq.1)then
                     do jjjj=1,3
                        pxprng1(jjjj)=pyp(nffin(iit,jjj),jjjj) !MBM
                     enddo
                     ptprng1 = pyp(nffin(iit,jjj),10) !MBM  
                     etaprng1 = pyp(nffin(iit,jjj),19) !MBM 
                     phiprng1 = pyp(nffin(iit,jjj),15) !MBM
                     do jjjj=1,4
                        pfullprng1(jj,jjjj)=pyp(nffin(iit,jjj),jjjj)
                     end do
                     pfullprng1(jj,5)=pyp(nffin(iit,jjj),10) ! jj^th Neutralino 
                     pfullprng1(jj,6)=pyp(nffin(iit,jjj),19)
                     pfullprng1(jj,8)=pyp(nffin(iit,jjj),15)
                     if(abs(k(nffin(iit,jjj),2)).eq.11)then
                        pfullprng1(jj,7)=1
                     else if(abs(k(nffin(iit,jjj),2)).eq.13)then
                        pfullprng1(jj,7)=2
                     else
                        pfullprng1(jj,7)=3 ! other than e/mu decays i.e., taus hadronic decays.
                     end if
                  else if(iprong.eq.2)then
                     do jjjj=1,3
                        pxprng2(jjjj)=pyp(nffin(iit,jjj),jjjj) !MBM
                     enddo
                     ptprng2 = pyp(nffin(iit,jjj),10) !MBM 
                     etaprng2 = pyp(nffin(iit,jjj),19) !MBM
                     phiprng2 = pyp(nffin(iit,jjj),15) !MBM
                     do jjjj=1,4
                        pfullprng2(jj,jjjj)=pyp(nffin(iit,jjj),jjjj)
                     end do
                     pfullprng2(jj,5)=pyp(nffin(iit,jjj),10)
                     pfullprng2(jj,6)=pyp(nffin(iit,jjj),19)
                     pfullprng2(jj,8)=pyp(nffin(iit,jjj),15)
                     pfullprng2(jj,7)=1
                     pfullprng1(jj,7)=0
                  else if(iprong.eq.3)then
                     do jjjj=1,3
                        pxprng3(jjjj)=pyp(nffin(iit,jjj),jjjj) !MBM
                     enddo
                     ptprng3 = pyp(nffin(iit,jjj),10) !MBM
                     etaprng3 = pyp(nffin(iit,jjj),19) !MBM
                     phiprng3 = pyp(nffin(iit,jjj),15) !MBM
c.....check pt prong cuts MBM
                      if(ptprng1.gt.2.d0.and.ptprng2.gt.2.d0
     $                  .and.ptprng3.gt.2.d0)then
                        if(ptprng1.gt.9.d0.or.ptprng2.gt.9.d0
     $                     .or.ptprng3.gt.9.d0)then
                           ptcut(jj)=.True.
                        endif
                     endif
c.....check tau cone cut MBM - if prongs are within DeltaR < 0.2,
                     do jjjj=1,3
                        pcone(jjjj)=0.d0
                     enddo
                     do jjjj=1,3
                        pcone(jjjj)=pcone(jjjj)+pxprng1(jjjj)+
     $                  pxprng2(jjjj)+pxprng3(jjjj)
                     enddo
                     ptau=dsqrt(pcone(1)**2+pcone(2)**2+pcone(3)**2)
                     pttau=dsqrt(pcone(1)**2+pcone(2)**2)
                     phitau=dacos(pcone(1)/pttau)*
     $                       pcone(2)/dabs(pcone(2))
                     etatau=log((ptau+pcone(3))/(ptau-pcone(3)))/2.d0
                     deltar1 = dsqrt((etaprng1-etatau)**2
     $                         +(phiprng1-phitau)**2)
                     deltar2 = dsqrt((etaprng2-etatau)**2
     $                         +(phiprng2-phitau)**2)
                     deltar3 = dsqrt((etaprng3-etatau)**2
     $                         +(phiprng3-phitau)**2)
                     if(deltar1.lt.0.2d0.and.deltar2.lt.0.2d0
     $                  .and.deltar3.lt.0.2d0)rcut(jj)=.True.
c.....end of check MBM
                     do jjjj=1,4
                        pfullprng3(jj,jjjj)=pyp(nffin(iit,jjj),jjjj)
                     end do
                     pfullprng3(jj,5)=pyp(nffin(iit,jjj),10)
                     pfullprng3(jj,6)=pyp(nffin(iit,jjj),19)
                     pfullprng3(jj,8)=pyp(nffin(iit,jjj),15)
                     pfullprng3(jj,7)=1
                     pfullprng1(jj,7)=0
                  end if
               end if
            end do

            if(iprong.gt.3)pfullprng3(jj,7)=-1
            nprong(jj)=iprong
            pfullprng2(jj,7)=iprong
            if(nprong(jj).ne.1.and.nprong(jj).ne.3)then
               write(*,*)'WARNING: more the 3 prongs in tau decay'
            end if

            if(iprong.eq.1)then
               do jjj=1,mtorig(iit)
                  if(abs(kchg(pycomp(kftdec(1,jjj)),1)).gt.0
     $                 .and.abs(kftdec(1,jjj)).gt.15)then !MBM2
                     nprongh(jj)=nprongh(jj)+1
                  elseif(abs(kftdec(1,jjj)).eq.13)then
                     nprongmu(jj)=nprongmu(jj)+1
                  endif
               end do
            end if
C.... (VII) TAU EFFICIENCIES
            qrand=rand(iev+jj)
            if(iprong.eq.1)then
               if(qrand.lt.0.5d0)efficiency(jj)=.True.
            elseif(iprong.eq.3)then
               if(qrand.lt.0.85d0)efficiency(jj)=.True.
            end if
            if(show)write(*,*)iprong,'-prong efficiency',qrand,
     $           efficiency(jj)
            do jjj=1,3
               pprong(jjj)=0d0
            end do
C.... (VI.a) Tau daughter properties 
            ki1t(iit)=newiline

            do 207 jjj=1,mtorig(iit)
               if(iprong.eq.1)then
                  if(abs(kchg(pycomp(kftdec(iit,jjj)),1)).gt.0)
     $                 then
                     if(show_tau_event)then
                        write(*,*)'p of',jjj,' tau daughter',
     $                       (p(nffin(iit,jjj),jjjj),jjjj=1,3),
     $                       pyp(nffin(iit,jjj),8)
                        write(*,*)'ptau',iit,(p(ki1t(iit),jjjj),
     $                       jjjj=1,3),pyp(ki1t(iit),8)
                     end if
                     c1prng=0d0
                     do iii=1,3
                        c1prng=c1prng+p(nffin(iit,jjj),iii)*
     $                       p(ki1t(iit),iii)
                     end do
                     c1prng=c1prng/(pyp(nffin(iit,jjj),8)*
     $                    pyp(ki1t(iit),8) )
                     if(show)write(*,*)'c1prng',c1prng
                     do jjjj=1,4
c....................define vprong
                        if(jjjj.le.3)pprong(jjjj)=
     $                       p(nffin(iit,jjj),jjjj)
                     end do
                     cprng=c1prng
                  end if
               end if
               if(iprong.eq.3)then
                  if(abs(kchg(pycomp(kftdec(1,jjj)),1)).gt.0)
     $                 then
                     if(show_tau_event)then
                        write(*,*)'p of',jjj,' tau daughter',
     $                       (p(nffin(iit,jjj),jjjj),jjjj=1,3),
     $                       pyp(nffin(iit,jjj),8)
                     end if
                     c3prng=0d0
                     do iii=1,4
                        if(iii.le.3)pprong(iii)=pprong(iii)
     $                       +p(nffin(iit,jjj),iii)
                     end do
                  end if
               end if
 207        CONTINUE

            if(iprong.eq.3)then
               p3pm=dsqrt(pprong(1)**2+pprong(2)**2
     $              +pprong(3)**2)
               c3prng=0d0
               do iii=1,3
                  c3prng=c3prng+pprong(iii)*
     $                 p(ki1t(iit),iii)
               end do
               c3prng=c3prng/( p3pm*pyp(ki1t(iit),8) )
               if(show_tau_event)then
                  write(*,*)'pprong',( pprong(iii),iii=1,3 ),p3pm,
     $                 c3prng,'total momentum of 3 prongs'
                  write(*,*)'ptau',iit,(p(ki1t(iit),jjjj),
     $                 jjjj=1,3),pyp(ki1t(iit),8)
                  write(*,*)'c3prng',c3prng
               end if
               cprng=c3prng
            end if
c     if(iprong.eq.1)write(66,*),c1prng
c     if(iprong.eq.3)write(67,*),c3prng
c.... For one inestable particle v(ii,5) is different from zero
c.... v is the point where particle was produced and vp where it
c.... will decay. 
            do jjjj=1,3
C.... v of tau (obtained from v(~chi_10,5))
               v1(jjjj)=vchi(jj,jjjj)
C.... v of prong (obtained directly from the prong)
C.... Note that is the same that the obtained from v(tau,5)
               v2(jjjj)=v(nffin(iit,1),jjjj)
               vtau(jj,jjjj)=v2(jjjj)
               p1(jjjj)=p(nfchi(jj,ii),jjjj) ! p of prong
               p2(jjjj)=pprong(jjjj) ! p of tau
            end do
            if(show)then
               write(*,*)'v tau',(v1(jjjj),jjjj=1,3)
               write(*,*)'v prong',(v2(jjjj),jjjj=1,3)
               write(*,*)'p tau',(p1(jjjj),jjjj=1,3)
               write(*,*)'p prong',(p2(jjjj),jjjj=1,3)
            end if
            if(iprong.eq.1.or.iprong.eq.3)then
C.... There are two aproachs to stabish if the tau prong intersect the neutralino
C.... displaced vertex sphere (sphere because the pixel resolution RRsili)
C.... 1) Two dimensional approach using angle between tau and prong
               distvtx(jj)=dsqrt((v1(1)-v2(1))**2
     $              +(v1(2)-v2(2))**2+(v1(3)-v2(3))**2)
               dprng(jj)=distvtx(jj)*dtan(dacos(cprng))
               if(show)write(*,*)'distance tau-prong<',RRsili,
     $              '?:',dprng(jj),distvtx(jj)
C.....2)Check intersection between prong line and resolution ellipsoid
C..... The formula for line-sphere intersection can be used, provided that
C.....  a) center of sphere, origin of line, direction of line be reescaled
C.....     by the axis of ellipsoid
C.....  b) r=1
C.....See mathematica/line-sphere*.nb for details 
               AA=RRsili ! RRsili
               BB=BBsili !RRsili
               RR=AA
               do jjjj=1,3
                  if(jjjj.eq.3)RR=BB
                  v1m(jjjj)=v1(jjjj)/RR
                  v2m(jjjj)=v2(jjjj)/RR
                  p2m(jjjj)=p2(jjjj)/RR
               end do
c               spheline=SLintersection(v1,v2,p2,RRsili)
               spheline=SLintersection(v1m,v2m,p2m,1d0)
               if(spheline.gt.0)sphere_line(jj)=.True.
               if(show)write(*,*)'intersection line sphere?',
     $              sphere_line(jj)
            else
               write(*,*)'ERROR: failed tau reconstruction'
               write(*,*)'fix tau -> ',
     $              (ntamendec(1,jjjj),jjjj=1,mtorig(1)),
     $              ' at event ',iev
c     call pylist(3)
            end if
            
         end if  tauanalysisif                ! end if over neutralino final states ii

C     :   OUTPUT: efficiency(jj),sphere_line(jj)
 206  CONTINUE
      end 
**********************************************************
**********************************************************
C...  INPUT:
c...   jj: jj-th neutralino
c...   chifin(jj): number of neutralino final states.
c...   kfchifin(jj,ii): kf number of ii neutralino daughter
c...   nfchi(jj,ii): line number of ii neutralino daughter
c...   vchi(jj,1-4): displaced vertex position of neutralino
c...   RRsili: transverse pixel resolution
c...   BBsili: longitudinal pixel resolution
c...   show: be verbose if TRUE
C...  OUTPUT
      subroutine tauanalysis2(jj,chifin,kfchifin,nfchi,vchi,
     $     RRsili,BBsili,show,show_tau_event,efficiency,
     $     sphere_line,tau_event,nprong,nprongh,npronge,
     $     nprongmu,pfullprng1,pfullprng2,pfullprng3,vtau)
C      nmuevt in commons
      implicit double precision(a-h,o-z)
      include 'pylocalcom.h' 
      EXTERNAL PYDATA
*******************-- local variables --********************************
      double precision ptdec(5,5,100)
      double precision vchi(10,4)
      double precision distvtx(10),dprng(10)
      double precision p1(3),p2(3)
      double precision v1(3),v2(3)
      double precision pprong(3)
      integer ki1t(5),mtorig(5)
      integer kftdec(5,10)
      character*16 ntamendec(5,10)
      logical show_tau_event,show
      logical tau_event,sphere_line(10)
      logical efficiency(5)
      integer chifin(5),kfchifin(5,10)
      integer nfchi(5,10)
c$$$      integer nwfin(5,10)
      integer kfstdaug(1000),lnstdaug(1000)
      double precision pstdaugm(5,1000),pxprng1(3),pxprng2(3)
      logical allfst,ptcut(5),rcut(5),etcut(5)
      character*16 namestdaug(1000)
      character*8 allchrgd
      integer nprong(10),nprongh(10),nprongmu(10),npronge(10)
      double precision v1m(3),v2m(3),p2m(3),pxprng3(3),pcone(3)
      double precision pfullprng1(10,8),pfullprng2(10,8)
      double precision pfullprng3(10,8),vtau(10,4)
***************__ Commons __**************
      double precision lum
      common/sig/sig_sm,lum
      integer nevgen,nevgenwmu,nevgenwtau,nevgenwe
      common/porodskands/nevgen,nevgenwmu,nevgenwtau,nevgenwe
      integer nffin(5,100)
      common/dauglines/nffin
      common/eevt/neevt
      common/muevt/nmuevt
      common/eventnumber/iev,nchilinejj
      common/teste/ievv,ptcut,rcut,etcut     !MBM
      distvtx(jj)=-1d0
      efficiency(jj)=.False.
      if(show)write(*,*)'tau decay'

      do 206 ii=1,chifin(jj)    ! loop over neutralino final states
         if(show)write(*,*)'~chi_10->',kfchifin(jj,ii),ii,iev
C.....(VI) IDENTIFYING tau DAUGHTERS
C.....Storing tau momentum

         tauanalysisif: if(abs(kfchifin(jj,ii)).eq.15)then !.and.nwtau(jj).gt.0)then ! tau id number is 15
            tau_event=.True.
            if(show)call pylist(3)
***********************************************************************
c...  Find what are the tau daughters
*=====================================================================
C.... (A) By using directly fulldecaychain with full fragmented table
            if(show)write(*,*)'tau line number:',nfchi(jj,ii),
     $           ', in previous fragmented table.'
C...  INPUT:
C...  nfchi(jj,ii): line number in table of tau coming from 
C...  ~chi_10 -> tau W
C..   *********************************************
C...  allfst: .False. Only charged final stable states (FSS)
C...  .True.  All FSS
            allfst=.True. 
C,,,  *********************************************
C...  OUTPUT:
C...  mtorig(1); number of final stable states (FSS) of tau
C...  kfstdaug(1)...kfstdaug(mtorig(1)): kf of FSS 
C...  lnstdaug(1)...lnstdaug(mtorig(1)): line number in table of FSS
C...  namestdaug(1)...namestdaug(mtorig(1)): names of FSS 
C...  pstdaug(1...5,1)...pstdaug(1...5,mtorig(1)): full momentum of FSS
            newiline=nfchi(jj,ii) !nfchi(jj,ii) ! 100
            if(show)write(*,*)'Finding the daughters of kf',
     $           k(newiline,2),' in line number',newiline
            iit=1               ! number of taus to be analysed 
            ki1t(iit)=newiline  ! line number in the table for tau
            call fulldecaychain(newiline,allfst,mtorig(iit),
     $           kfstdaug,lnstdaug,namestdaug,pstdaugm)
            ki1t(iit)=newiline
C...  Set tau daughters output variables:

            do ikik=1,mtorig(iit)
               kftdec(iit,ikik)=kfstdaug(ikik)
               nffin(iit,ikik)=lnstdaug(ikik)
               ntamendec(iit,ikik)=namestdaug(ikik)
               do ikiki=1,5 
                  ptdec(ikiki,iit,ikik)=pstdaugm(ikiki,ikik)
               enddo
            enddo
*=====================================================================
c...  (B) By using decaychainnew:
c     call decaychainnew(chiline(jj),1,abs(kfchifin(jj,ii)), ! used after pyexec
c     $                    iit,ki1t,mtorig,ntamendec,kftdec,ptdec)
*=====================================================================
            if(allfst)then
               allchrgd='  all'
            else
               allchrgd='charged'
            end if
            if(show)then
               write(*,*)'Total number of ',allchrgd,
     $              'stable daugthers:',mtorig(iit),' of line',
     $              ki1t(iit)
               write(*,*)'kf of daug',
     $              (kftdec(iit,jjjj),jjjj=1,mtorig(iit))
               write(*,*)'line number of daug',
     $              (nffin(iit,jjjj),jjjj=1,mtorig(iit))
               write(*,*)'name of daug: ',
     $              (ntamendec(iit,jjjj),jjjj=1,mtorig(iit))
               write(*,*)'p_1 of daug ',
     $              (ptdec(1,iit,jjjj),jjjj=1,mtorig(iit))
            endif
***********************************************************************
c.... Check what are the tau daughters
            ptcut(jj) = .False.  !MBM
            rcut(jj) = .False.  !MBM
            iprong=0
            nprongh(jj)=0
            npronge(jj)=0
            nprongmu(jj)=0

            do jjj=1,mtorig(iit)
               if(abs(kchg(pycomp(kftdec(1,jjj)),1)).gt.0)then
                  iprong=iprong+1
                  if(iprong.eq.1)then
                     do jjjj=1,3
                        pxprng1(jjjj)=pyp(nffin(iit,jjj),jjjj) !MBM
                     enddo
                     ptprng1 = pyp(nffin(iit,jjj),10) !MBM
                     etaprng1 = pyp(nffin(iit,jjj),19) !MBM
                     phiprng1 = pyp(nffin(iit,jjj),15) !MBM
                     do jjjj=1,4
                        pfullprng1(jj,jjjj)=pyp(nffin(iit,jjj),jjjj)
                     end do
                     pfullprng1(jj,5)=pyp(nffin(iit,jjj),10)
                     pfullprng1(jj,6)=pyp(nffin(iit,jjj),19)
                     pfullprng1(jj,8)=pyp(nffin(iit,jjj),15)
                     if(abs(k(nffin(iit,jjj),2)).eq.11)then
                        pfullprng1(jj,7)=1
                     else if(abs(k(nffin(iit,jjj),2)).eq.13)then
                        pfullprng1(jj,7)=2
                     else
                        pfullprng1(jj,7)=3
                     end if
                  else if(iprong.eq.2)then
                     do jjjj=1,3
                        pxprng2(jjjj)=pyp(nffin(iit,jjj),jjjj) !MBM
                     enddo
                     ptprng2 = pyp(nffin(iit,jjj),10) !MBM
                     etaprng2 = pyp(nffin(iit,jjj),19) !MBM
                     phiprng2 = pyp(nffin(iit,jjj),15) !MBM
                     do jjjj=1,4
                        pfullprng2(jj,jjjj)=pyp(nffin(iit,jjj),jjjj)
                     end do
                     pfullprng2(jj,5)=pyp(nffin(iit,jjj),10)
                     pfullprng2(jj,6)=pyp(nffin(iit,jjj),19)
                     pfullprng2(jj,8)=pyp(nffin(iit,jjj),15)
                     pfullprng2(jj,7)=1
                     pfullprng1(jj,7)=0
                  else if(iprong.eq.3)then
                     do jjjj=1,3
                        pxprng3(jjjj)=pyp(nffin(iit,jjj),jjjj) !MBM
                     enddo
                     ptprng3 = pyp(nffin(iit,jjj),10) !MBM
                     etaprng3 = pyp(nffin(iit,jjj),19) !MBM
                     phiprng3 = pyp(nffin(iit,jjj),15) !MBM
c.....check pt prong cuts MBM
                      if(ptprng1.gt.2.d0.and.ptprng2.gt.2.d0
     $                  .and.ptprng3.gt.2.d0)then
                        if(ptprng1.gt.9.d0.or.ptprng2.gt.9.d0
     $                     .or.ptprng3.gt.9.d0)then
                           ptcut(jj)=.True.
                        endif
                     endif
c.....check tau cone cut MBM - if prongs are within DeltaR < 0.2,
                     do jjjj=1,3
                        pcone(jjjj)=0.d0
                     enddo
                     do jjjj=1,3
                        pcone(jjjj)=pcone(jjjj)+pxprng1(jjjj)+
     $                  pxprng2(jjjj)+pxprng3(jjjj)
                     enddo
                     ptau=dsqrt(pcone(1)**2+pcone(2)**2+pcone(3)**2)
                     pttau=dsqrt(pcone(1)**2+pcone(2)**2)
                     phitau=dacos(pcone(1)/pttau)*
     $                       pcone(2)/dabs(pcone(2))
                     etatau=log((ptau+pcone(3))/(ptau-pcone(3)))/2.d0
                     deltar1 = dsqrt((etaprng1-etatau)**2
     $                         +(phiprng1-phitau)**2)
                     deltar2 = dsqrt((etaprng2-etatau)**2
     $                         +(phiprng2-phitau)**2)
                     deltar3 = dsqrt((etaprng3-etatau)**2
     $                         +(phiprng3-phitau)**2)
                     if(deltar1.lt.0.2d0.and.deltar2.lt.0.2d0
     $                  .and.deltar3.lt.0.2d0)rcut(jj)=.True.
c.....end of check MBM
                     do jjjj=1,4
                        pfullprng3(jj,jjjj)=pyp(nffin(iit,jjj),jjjj)
                     end do
                     pfullprng3(jj,5)=pyp(nffin(iit,jjj),10)
                     pfullprng3(jj,6)=pyp(nffin(iit,jjj),19)
                     pfullprng3(jj,8)=pyp(nffin(iit,jjj),15)
                     pfullprng3(jj,7)=1
                     pfullprng1(jj,7)=0
                  end if
               end if
            end do

            if(iprong.gt.3)pfullprng3(jj,7)=-1
            nprong(jj)=iprong
            pfullprng2(jj,7)=iprong
            if(nprong(jj).ne.1.and.nprong(jj).ne.3)then
               write(*,*)'WARNING: more the 3 prongs in tau decay'
            end if

            if(iprong.eq.1)then
               do jjj=1,mtorig(iit)
                  if(abs(kchg(pycomp(kftdec(1,jjj)),1)).gt.0
     $                 .and.abs(kftdec(1,jjj)).gt.15)then !MBM2
                     nprongh(jj)=nprongh(jj)+1
                  elseif(abs(kftdec(1,jjj)).eq.11)then
                     npronge(jj)=npronge(jj)+1
                  elseif(abs(kftdec(1,jjj)).eq.13)then
                     nprongmu(jj)=nprongmu(jj)+1
                  endif
               end do
            end if
C.... (VII) TAU EFFICIENCIES
            qrand=rand(iev+jj)
            if(iprong.eq.1)then
               if(qrand.lt.0.5d0)efficiency(jj)=.True.
            elseif(iprong.eq.3)then
               if(qrand.lt.0.85d0)efficiency(jj)=.True.
            end if
            if(show)write(*,*)iprong,'-prong efficiency',qrand,
     $           efficiency(jj)
            do jjj=1,3
               pprong(jjj)=0d0
            end do
C.... (VI.a) Tau daughter properties 
            ki1t(iit)=newiline

            do 207 jjj=1,mtorig(iit)
               if(iprong.eq.1)then
                  if(abs(kchg(pycomp(kftdec(iit,jjj)),1)).gt.0)
     $                 then
                     if(show_tau_event)then
                        write(*,*)'p of',jjj,' tau daughter',
     $                       (p(nffin(iit,jjj),jjjj),jjjj=1,3),
     $                       pyp(nffin(iit,jjj),8)
                        write(*,*)'ptau',iit,(p(ki1t(iit),jjjj),
     $                       jjjj=1,3),pyp(ki1t(iit),8)
                     end if
                     c1prng=0d0
                     do iii=1,3
                        c1prng=c1prng+p(nffin(iit,jjj),iii)*
     $                       p(ki1t(iit),iii)
                     end do
                     c1prng=c1prng/(pyp(nffin(iit,jjj),8)*
     $                    pyp(ki1t(iit),8) )
                     if(show)write(*,*)'c1prng',c1prng
                     do jjjj=1,4
c....................define vprong
                        if(jjjj.le.3)pprong(jjjj)=
     $                       p(nffin(iit,jjj),jjjj)
                     end do
                     cprng=c1prng
                  end if
               end if
               if(iprong.eq.3)then
                  if(abs(kchg(pycomp(kftdec(1,jjj)),1)).gt.0)
     $                 then
                     if(show_tau_event)then
                        write(*,*)'p of',jjj,' tau daughter',
     $                       (p(nffin(iit,jjj),jjjj),jjjj=1,3),
     $                       pyp(nffin(iit,jjj),8)
                     end if
                     c3prng=0d0
                     do iii=1,4
                        if(iii.le.3)pprong(iii)=pprong(iii)
     $                       +p(nffin(iit,jjj),iii)
                     end do
                  end if
               end if
 207        CONTINUE

            if(iprong.eq.3)then
               p3pm=dsqrt(pprong(1)**2+pprong(2)**2
     $              +pprong(3)**2)
               c3prng=0d0
               do iii=1,3
                  c3prng=c3prng+pprong(iii)*
     $                 p(ki1t(iit),iii)
               end do
               c3prng=c3prng/( p3pm*pyp(ki1t(iit),8) )
               if(show_tau_event)then
                  write(*,*)'pprong',( pprong(iii),iii=1,3 ),p3pm,
     $                 c3prng,'total momentum of 3 prongs'
                  write(*,*)'ptau',iit,(p(ki1t(iit),jjjj),
     $                 jjjj=1,3),pyp(ki1t(iit),8)
                  write(*,*)'c3prng',c3prng
               end if
               cprng=c3prng
            end if
c     if(iprong.eq.1)write(66,*),c1prng
c     if(iprong.eq.3)write(67,*),c3prng
c.... For one inestable particle v(ii,5) is different from zero
c.... v is the point where particle was produced and vp where it
c.... will decay. 
            do jjjj=1,3
C.... v of tau (obtained from v(~chi_10,5))
               v1(jjjj)=vchi(jj,jjjj)
C.... v of prong (obtained directly from the prong)
C.... Note that is the same that the obtained from v(tau,5)
               v2(jjjj)=v(nffin(iit,1),jjjj)
               vtau(jj,jjjj)=v2(jjjj)
               p1(jjjj)=p(nfchi(jj,ii),jjjj) ! p of prong
               p2(jjjj)=pprong(jjjj) ! p of tau
            end do
            if(show)then
               write(*,*)'v tau',(v1(jjjj),jjjj=1,3)
               write(*,*)'v prong',(v2(jjjj),jjjj=1,3)
               write(*,*)'p tau',(p1(jjjj),jjjj=1,3)
               write(*,*)'p prong',(p2(jjjj),jjjj=1,3)
            end if
            if(iprong.eq.1.or.iprong.eq.3)then
C.... There are two aproachs to stabish if the tau prong intersect the neutralino
C.... displaced vertex sphere (sphere because the pixel resolution RRsili)
C.... 1) Two dimensional approach using angle between tau and prong
               distvtx(jj)=dsqrt((v1(1)-v2(1))**2
     $              +(v1(2)-v2(2))**2+(v1(3)-v2(3))**2)
               dprng(jj)=distvtx(jj)*dtan(dacos(cprng))
               if(show)write(*,*)'distance tau-prong<',RRsili,
     $              '?:',dprng(jj),distvtx(jj)
C.....2)Check intersection between prong line and resolution ellipsoid
C..... The formula for line-sphere intersection can be used, provided that
C.....  a) center of sphere, origin of line, direction of line be reescaled
C.....     by the axis of ellipsoid
C.....  b) r=1
C.....See mathematica/line-sphere*.nb for details 
               AA=RRsili ! RRsili
               BB=BBsili !RRsili
               RR=AA
               do jjjj=1,3
                  if(jjjj.eq.3)RR=BB
                  v1m(jjjj)=v1(jjjj)/RR
                  v2m(jjjj)=v2(jjjj)/RR
                  p2m(jjjj)=p2(jjjj)/RR
               end do
c               spheline=SLintersection(v1,v2,p2,RRsili)
               spheline=SLintersection(v1m,v2m,p2m,1d0)
               if(spheline.gt.0)sphere_line(jj)=.True.
               if(show)write(*,*)'intersection line sphere?',
     $              sphere_line(jj)
            else
               write(*,*)'ERROR: failed tau reconstruction'
               write(*,*)'fix tau -> ',
     $              (ntamendec(1,jjjj),jjjj=1,mtorig(1)),
     $              ' at event ',iev
c     call pylist(3)
            end if
            
         end if  tauanalysisif                ! end if over neutralino final states ii

C     :   OUTPUT: efficiency(jj),sphere_line(jj)
 206  CONTINUE
      end 
**********************************************************
      subroutine chipartonline(ndaug,kfdaug,nldaug,ndaugnew,nlnew)
C...   Which of the jets comes from 
C...   INPUT:
C...    ndaug: number of neutralino daughters
C...    v1(*): displaced vertex of each neutralino
C...    rpix: ratio of first silicon cilinder
C...    kfdaug(*): KF code of neutralino daughters
C...    nldaug(*): line number in table of neutralino daughters
C...   OUTPUT: 
C...    ndaugnew: number of neutralino daughters at parton level
C...    nlnew(ndaugnew): line in table of each of neutralino daughters at parton level
      implicit double precision(a-h,o-z)
C.... INPUT
      integer kfdaug(*),nldaug(*)
      integer nlnew(20)
C...  Internal variables
C... PYTHIA COMMONS
C.....The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
******************PREAMBLE**************************
C...  Here I need to know what are the quarks which are associated to W
*************************************************
C...  Obtain the line numbers of final (parton level) neutalino daughters
      ndaugnew=0
      do i=1,ndaug
         if(abs(kfdaug(i)).ge.11.and.abs(kfdaug(i)).le.15)then 
            ndaugnew=ndaugnew+1
            nlnew(ndaugnew)=nldaug(i)
         end if
         if(k(nldaug(i),1).eq.21)then
            do j=i+1,n
               if(k(j,3).eq.nldaug(i).and.
     $              k(j,1).eq.21)then
                  ndaugnew=ndaugnew+1
                  nlnew(ndaugnew)=j
               end if
            end do
         end if
      end do
*************************************************
      END
**********************************************************
**********************************************************
      subroutine hitcilinder(iline,v1,rpix,hpixmax)
C...   Which of the jets comes from 
C...   INPUT:
C...    iline: line number in table of track to be analysed
C...    v1(*): displaced vertex of each neutralino
C...    rpix: ratio of first silicon cilinder
C...   OUTPUT: 
C.. .   hpixmax: intersection of p->=Sum_i p(iline,i) with cilinder
      implicit double precision(a-h,o-z)
C...  INPUT
      double precision v1(3)
C...  Internal variables
      double precision pl(3),plu(3)
      logical show
C... PYTHIA COMMONS
C.....The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
******************PREAMBLE**************************
      show=.False.
C...  Here I need to know what are the quarks which are associated to W
*************************************************
C...  Obtain the line numbers of final (parton level) neutalino daughters
      do i=1,3
         pl(i)=p(iline,i)
      end do
      call univec(pl,plu)
      if(show)then
         write(*,*)'v',(v1(i),i=1,3)
         write(*,*)'P_L',(pl(i),i=1,3)
         write(*,*)'P_LU',(plu(i),i=1,3)
      end if
      a=-( plu(1)*v1(1) + plu(2)*v1(2) +  
     $     dsqrt(rpix**2*(plu(1)**2 + plu(2)**2) -
     $     (plu(2)*v1(1)-plu(1)*v1(2))**2) )/((plu(1)**2 +
     $     plu(2)**2))
      if(a.lt.0)
     $    a=( -plu(1)*v1(1) - plu(2)*v1(2) + 
     $        dsqrt(rpix**2*(plu(1)**2 + plu(2)**2) -
     $        (plu(2)*v1(1)-plu(1)*v1(2))**2) )/((plu(1)**2 +
     $        plu(2)**2))
      
      hpixmax=v1(3)+a*plu(3)
      if(show)write(*,*)'a',a

      END
**********************************************************
**********************************************************
      subroutine  invmassfromall(v1,RRsili,minvchi)
      implicit double precision(a-h,o-z)
      double precision v1(3),v2(3),minvchi(10),pii(4),p2(3)
      parameter(nnpart=1000)   
      double precision ppartm(4,nnpart)
C...  PYTHIA COMMONS
C.....The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      
*************
      npart=0
C...  Store vertex and momentum of each stable particle in the table
      do ii=1,n
         if(k(ii,1).lt.11)then
            do jjj=1,4
               if(jjj.lt.4)then
                  v2(jjj)=v(ii,jjj)
                  p2(jjj)=p(ii,jjj)
               end if
               pii(jjj)=p(ii,jjj)
            end do
            spheline=SLintersection(v1,v2,p2,RRsili)
            if(spheline.gt.0)then !.and.
c     $                 abs(kchg(pycomp(k(ii,2)),1)).gt.0)then
               npart=npart+1
               do jjj=1,4
                  ppartm(jjj,npart)=pii(jjj)
               end do
            endif
         end if
      end do
c...  Invariant mass
      minvchi(jj)=DInvariantMass(npart,ppartm)
      end
******************************************************
      subroutine triggers_old(triggerps)
      implicit double precision(a-h,o-z)
      logical triggerps,show
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
******************************************************
C...  (V)  TRIGGERS CUTS FROM POROD & SKANDS
      show=.False.
******************************************************
      njets=0
      nleptons=0
      nbottom=0
      njetscheck=0
      nlcheck=0
         do ii=1,n
C..   jets 
            if(k(ii,1).eq.21)then ! check only parton level
               if(abs(k(ii,2)).le.6.or.abs(k(ii,2)).eq.15.or.
     $              abs(k(ii,2)).eq.21)then
                  njetscheck=njetscheck+1
c                  if(show)write(*,*)'q gluon tau',njetscheck
                  if(pyp(ii,10).gt.100.and.abs(pyp(ii,19)).lt.4.9d0)then
                     njets=njets+1
                     if(show)write(*,*)'jet p_t>100 and eta<4.9',
     $                    pyp(ii,10),pyp(ii,19),njets
                  end if
               end if
               if(abs(k(ii,2)).eq.11.or.abs(k(ii,2)).eq.13)then
                  nlcheck=nlcheck+1
                  if(show)write(*,*)'e or mu',nlcheck+1
                  if(pyp(ii,10).gt.20.and.abs(pyp(ii,19)).lt.2.5d0)then
                     nleptons=nleptons+1
                     if(show)write(*,*)'e or mu p_t>20 and eta<2.5',
     $                    pyp(ii,10),pyp(ii,19),nleptons
                  end if
               end if
               if(abs(k(ii,2)).eq.5.and.abs(pyp(ii,19)).lt.4.9d0)
     $              nbottom=nbottom+1
            end if
         end do
C...evento false because porod skands triggers
         if(njets.ge.4)triggerps=.True.
         if(show)write(*,*)'tigger njets>4',triggerps
         if(njets.ge.2.and.nleptons.ge.1)triggerps=.True.
         if(show)write(*,*)'tigger njets>2,nl>1',triggerps
         if(njets.ge.1.and.nleptons.ge.2)triggerps=.True.
         if(show)write(*,*)'tigger njets>1,nl>2',triggerps
         if(nbottom.gt.0)triggerps=.False.
         if(show)write(*,*)'tigger nbot=0',triggerps

      end
******************************************************
******************************************************
      subroutine triggersdv(trigger,Drimin,Drimax,njet,
     $     ET_isola)
      implicit double precision(a-h,o-z)
      logical trigger,show
      integer trigg(10)
C...  INPUT/OUTPUT
      double precision DRimin(15),DRimax(15)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C.....Three Pythia functions return integers, so need declaring.
      INTEGER PYK,PYCHGE,PYCOMP
C.....Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
******************************************************
C...  (V)  TRIGGERS CUTS FROM POROD & SKANDS
      show=.False.
******************************************************
      trigger=.False. 
      ntrigger=5 ! number of different triggers
      n2e=0
      etmx=0d0; etmy=0d0
      do i=1,10
         trigg(i)=0
      end do
      do 1012 ii=1,n
         iitable: if(k(ii,1).lt.11)then ! check for stable particles
            if(abs(k(ii,2)).eq.11)then
               if(pyp(ii,10).gt.15d0.and.abs(pyp(ii,19)).lt.2.5d0)then
                  if(show)write(*,*)'electron pt>15',ii
                  newline=ii
                  nkfdaug=abs(k(ii,2))
                  call isolai(newline,Drimin,Drimax,nkfdaug,ET,pt)
                  if(et.lt.ET_isola)then
                     if(show)write(*,*)'isolated electron pt=',
     $                    pyp(ii,10),'>15',ii
                     n2e=n2e+1
                     if(n2e.eq.2)then
                        trigg(2)=1 ! two isolated electrons pt>15
                        trigger=.True. 
                        goto 1013
                     end if
                     if(pyp(ii,10).gt.20d0)then
                        trigg(1)=1 ! one isolated electron pt>20
                        trigger=.True.
                        goto 1013
                     end if
                  end if
               end if
            end if
            if(abs(k(ii,2)).eq.13)then
               if(pyp(ii,10).gt.6d0.and.abs(pyp(ii,19)).lt.2.5d0)then
                  if(show)write(*,*)'muon pt>6 and eta<2.5',ii
                  newline=ii
                  nkfdaug=abs(k(ii,2))
                  call isolai(newline,Drimin,Drimax,nkfdaug,ET,pt)
                  if(et.lt.ET_isola)then
                  if(show)write(*,*)'isolated muon pt>6  and eta<2.5',ii
                     trigg(3)=1  ! one isolated muon pt>6
                     trigger=.True.
                     goto 1013
                  end if
               end if
            end if
            if(dabs(pyp(ii,19)).gt.5)then
               etmx =etmx+
     $              pyp(ii,1)/pyp(ii,10)*pyp(ii,4)*dsin(pyp(ii,13))
               etmy =etmy+
     $              pyp(ii,2)/pyp(ii,10)*pyp(ii,4)*dsin(pyp(ii,13))
            else
               if(k(ii,2).ne.22.and.kchg(pycomp(k(ii,2)),1).eq.0)then
                  etmx =etmx+
     $                 pyp(ii,1)/pyp(ii,10)*pyp(ii,4)*dsin(pyp(ii,13))
                  etmy =etmy+
     $                 pyp(ii,2)/pyp(ii,10)*pyp(ii,4)*dsin(pyp(ii,13))
               end if
            end if
         end if iitable
 1012 CONTINUE
 1013 CONTINUE ! exit loop
      if(trigger)then
         if(show)write(*,*)'trigger',(trigg(i),i=1,ntrigger)
         return
      end if
      do ii=1,njet
         if(pyp(n+ii,10).gt.100d0.and.abs(pyp(i,19)).lt.5.d0)then 
            if(show)write(*,*)'jet with pt>100',ii
            trigg(4)=1          ! jet with pt> 100
            trigger=.True.
            goto 1014
         end if
      end do
 1014 CONTINUE
      if(trigger)then
         if(show)write(*,*)'trigger',(trigg(i),i=1,ntrigger)
         return
      end if

      if (show)write(*,*)'etmx',etmx,etmy
      etmis=dsqrt(etmx**2+etmy**2)
      if(etmis.gt.100.d0)then
         if(show)write(*,*)'E_T missing > 1000'
         trigger=.true.
      end if
      if(trigger)then
         if(show)write(*,*)'trigger',(trigg(i),i=1,ntrigger)
         return
      end if
      if(show)write(*,*)'trigger',(trigg(i),i=1,ntrigger)
      end
******************************************************
C.............................................................................
C....... SUBROUTINE THAT DETERMINES WETHER AN EVENT PASSES THE EXPERIMENTAL 
C.......               TRIGGERS OR NOT - SUITED FOR LHC LEVEL 1
C.............................................................................

      subroutine trigger_lhc(njet,trigger)
C.....INPUT: njet
C.....Output: trigger
      implicit double precision(a-h,o-z)

      logical trigger
      logical nodebug

c      common/trgr/trigger,etmx,etmy,njet,nleptgr,njettgr,nmistgr,ptx,pty
c     &               nevel4,nevmu4
C.....Some Pythia functions return integers, so need declaring.
      INTEGER PYCOMP

      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)

c..... trigger over leptons
c..... we require either a muon with pt > 6 GeV,
c..... or an electron or photon with pt > 20 GeV or 
c..... or two electrons or photons with pt > 15 GeV or 
      nodebug=.False. 
      trigger=.False. 

c     nleptgr = 0
      nem = 0
      etmx=0d0
      etmy=0d0
      do 666 i=1,n
         if(abs(k(i,2)).eq.13.and.k(i,1).eq.1)then
            if(pyp(i,10).gt.6.d0.and.dabs(pyp(i,19)).lt.3.d0)then
               trigger=.true.
c               nleptgr=nleptgr+1
            endif
         endif
         if(trigger)goto 667

         if((abs(k(i,2)).eq.11.or.k(i,2).eq.22).and.k(i,1).eq.1)then
            if(pyp(i,10).gt.20.d0.and.dabs(pyp(i,19)).lt.3.d0)then 
               trigger=.true.
c               nleptgr=nleptgr+1
            elseif(pyp(i,10).gt.15.d0.and.dabs(pyp(i,19)).lt.3.d0)then
               nem = nem+1
            endif
         endif

         if(nem.ge.2)trigger=.true.
         if(trigger)goto 667

c.... triger over missing transversal energy > 100 GeV
c      nmistgr = 0
c.... Missing energy from  stable particles close to beam pipe
         if(k(i,1).lt.11.and.dabs(pyp(i,19)).gt.5)then 
            etmx = etmx+pyp(i,1)/pyp(i,10)*pyp(i,4)*dsin(pyp(i,13))
            etmy = etmy+pyp(i,2)/pyp(i,10)*pyp(i,4)*dsin(pyp(i,13))
c... Missing energy from stable neutral particles
         elseif(k(i,1).lt.11.and.kchg(pycomp(k(i,2)),1).eq.0)then
            etmx = etmx+pyp(i,1)/pyp(i,10)*pyp(i,4)*dsin(pyp(i,13))
            etmy = etmy+pyp(i,2)/pyp(i,10)*pyp(i,4)*dsin(pyp(i,13))
         endif
 666  CONTINUE
 667  RETURN
      etmis=dsqrt(etmx**2+etmy**2)
      if(etmis.gt.100.d0)then
         trigger=.true.
      endif
      if(trigger)RETURN

c..... trigger over jets with pt > 100 GeV
      
c      njettgr = 0
      if(njet.gt.0)then
         do i=n+1,n+njet
            if(pyp(i,10).gt.100.d0.and.dabs(pyp(i,19)).lt.5.d0)then 
               trigger=.true.
c               njettgr=njettgr+1
            endif
         enddo
      endif

      return
      end

* ** This is added on 23.10.2010 from sv.r57 version

****************************************************
      subroutine hitcilinder_update(iev,plep_hitcil,pjw_hitcil,v1,
     +rpix,hpix,xVec_r_mod_t_pos,plep_theta,pquark_nWdecay)

C...   Which of the jets comes from 
C...   INPUT:
C...    iline: line number in table of track to be analysed
* **    nlnew(jj),jj=1,3 are the LSP daughters...
C...    v1(*): displaced vertex of each neutralino
C...    rpix: ratio of first silicon cilinder

* ** spd: iline is calculated from "chipartonline" subroutines...
* ** However, I think the "iline" should be modified from LSPs 
* ** 3-body decay including the qq'. 
* ** p,plu is the vector and unit-vector...

C...   OUTPUT: 
C.. .   hpixmax: intersection of p->=Sum_i p(iline,i) with cilinder
      implicit double precision(a-h,o-z)
C...  INPUT
      integer nlnew(3),iev
      double precision v1(3)
      double precision pjw_hitcil(4,100)
      double precision plep_hitcil(10,8)
      double precision pquark_nWdecay(10,12)
      double precision Sol_1(3),Sol_2(3)
C...  Internal variables
      double precision pl(3),plu(3)
      logical show
C... PYTHIA COMMONS
C.....The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)

* ** this is the common block to pass the LSP info from main programme...

      common/kinematicsLSP/decayLtmp,decayLps,x_boostLSP,x_bLSP
     +       ,x_thetaLSP,x_yLSP,x_etaLSP,x_cal_yLSP,x_cal_etaLSP
     +       ,x_totE_LSP,x_magP_LSP,x_mass_LSP

******************PREAMBLE**************************
      show=.True. ! .False.  spd: it was FALSE
C...  Here I need to know what are the quarks which are associated to W
*************************************************
C...  Obtain the line numbers of final (parton level) neutalino daughters
* ** here the l,jj momentum should be stored ...
* ** 1st check: the parton level daughters...(this is been set...)
*           write(*,*)'in the subroutine '
*           write(*,*)'lep-1-mom',(plep_hitcil(1,kk),kk=1,4)
*           write(*,*)'lep-1-pt ',plep_hitcil(1,5)
*           write(*,*)'lep-1-eta',plep_hitcil(1,6)
*           write(*,*)'lep-1-phi',plep_hitcil(1,8)
*           do jjj=1,4
*           write(*,*)'jet-',jjj,'-mom',(pjw_hitcil(kk,jjj),kk=1,4)
*           write(*,*)'quark-',jjj,'-mom',(pquark_nWdecay(jjj,kk),kk=1,4)
*           write(*,*)'quark-',jjj,'-pt',pquark_nWdecay(jjj,5)
*           write(*,*)'quark-',jjj,'-Eta',pquark_nWdecay(jjj,6)
*           write(*,*)'quark-',jjj,'-Phi',pquark_nWdecay(jjj,8)
*           write(*,*)'quark-',jjj,'-theta',pquark_nWdecay(jjj,12)
*           enddo 
*        write(*,*)'inside hitcil subroutine...' 
* ** w/o factor of 5 the rpix and hpix are 0.02d0,0.5d0 in mm respectively.

      aa=0.02d0*5d0  ! ellipsoid production  delta_xy
      ba=0.50d0*5d0  ! ellipsoid production  delta_z

      rpix=aa
      hpix=ba

      do i=1,3
         pl(i)=0.d0
         Sol_1(i)=0.d0
         Sol_2(i)=0.d0
      end do

         xMomentum=0.d0
         pl_totE=0.d0
         x_mLSP=0.d0
         x_y_Reco=0.d0
         x_eta_Reco=0.d0

      xDisplaced=0.d0
      xMomentum1=0.d0

      xt_1=0.d0
      xt_2=0.d0

      xDot_plu_r0=0.d0
      xVec_r_1=0.d0
      xVec_r_2=0.d0
      xVec_r_3=0.d0
      xVec_r_mod_t_pos =0.d0
      xVec_r_mod_t_neg =0.d0

      xA=0.d0
      xB=0.d0
      xC=0.d0
      xDis=0.d0
      xDelR_jj=0.d0

      do i=1,3
         pl(i)=plep_hitcil(1,i)+pjw_hitcil(i,1)+pjw_hitcil(i,2)
      end do

      xLSPmass= 0.d0

      xLSPmass=(pjw_hitcil(4,1) + pjw_hitcil(4,2) + plep_hitcil(1,4))**2
     +     -   (pjw_hitcil(3,1) + pjw_hitcil(3,2) + plep_hitcil(1,3))**2
     +     -   (pjw_hitcil(2,1) + pjw_hitcil(2,2) + plep_hitcil(1,2))**2
     +     -   (pjw_hitcil(1,1) + pjw_hitcil(1,2) + plep_hitcil(1,1))**2
       
      if ( xLSPmass.lt.0.d0) go to 999

      xLSPmass=dsqrt(xLSPmass)
*      write(*,*)'LSP masses',xLSPmass
****************************************************
* ** Uncertainties and Error propagation...
* ** Jets... 
* ** IMPORTANT; Please note that for our first approximations, 
* ** we have to use the parton level 4-momentum and smeared 
* ** according to the parametrizations...
*           write(*,*)'jet-',jjj,'-mom',(pjw_hitcil(kk,jjj),kk=1,4)
*           write(*,*)'quark-',jjj,'-mom',(pquark_nWdecay(jjj,kk),kk=1,4)

* ** If this is on then quark level. If off then jet level.
***           do kk=1,4
***              pjw_hitcil(kk,2)= pquark_nWdecay(2,kk)
***              pjw_hitcil(kk,1)= pquark_nWdecay(1,kk)
***           enddo

***       do kkj=1,2 ! quarks 1 and 2.
***          pnfjw_hitcil(1,kkj)= pquark_nWdecay(kkj,8)  !phi 
***          pnfjw_hitcil(2,kkj)= pquark_nWdecay(kkj,6)  !eta
***       enddo

* ******************************************************************
* ** LSP -mass, meomentum reconstruction...

         xMomentum=dsqrt(pl(1)**2+pl(2)**2+pl(3)**2)
         pl_totE=plep_hitcil(1,4)+pjw_hitcil(4,1)+pjw_hitcil(4,2)

         x_mLSP=dsqrt(pl_totE**2-xMomentum**2)
         x_y_Reco=0.5d0*(log((pl_totE+pl(3))/(pl_totE-pl(3))))
         x_eta_Reco=0.5d0*(log((xMomentum+pl(3))/(xMomentum-pl(3))))

*        write(*,*)'totE,abs(p),m LSP:',x_totE_LSP,x_magP_LSP,x_mass_LSP
*        write(*,*)'totE,abs(p),m from LSP-daughters:',pl_totE,xMomentum
*     +                                               ,x_mLSP  
*        write(*,*)'Y,eta directly from LSP:',x_yLSP,x_etaLSP
*        write(*,*)'Y,eta_Reco from LSP-daughters:',x_y_Reco,x_eta_Reco

      call univec(pl,plu)
     
* ** In this case vl(1,2,3)=Displaced VERTEX ...

      xDisplaced=dsqrt(v1(1)**2+v1(2)**2+v1(3)**2)
      xMomentum1=dsqrt(plu(1)**2+plu(2)**2+plu(3)**2)

*      if(show)then
*        write(*,*)'aa,ba in subroutine:',rpix,hpix 
*        write(*,*)'v',(v1(i),i=1,3),'mod(v):Length:=>',xDisplaced
*        write(*,*)'P_L',(pl(i),i=1,3),'magnitude',xMomentum
*        write(*,*)'P_LU',(plu(i),i=1,3),'magnitude1',xMomentum1
*      end if

      xA=pl(1)**2+pl(2)**2
      xB=2.d0*(v1(1)*pl(1)+ v1(2)*pl(2))
      xC=v1(1)**2+v1(2)**2-rpix**2

      xDis=xB**2-4.d0*xA*xC

      if(xDis.lt.0) then
***      write(*,*)'No solution...evt',iev
      go to 999
      endif

      xt_1=(-1.d0*xB+dsqrt(xDis))/(2.d0*xA)
      xt_2=(-1.d0*xB-dsqrt(xDis))/(2.d0*xA)

      Sol_1(1)=v1(1)+pl(1)*xt_1
      Sol_2(1)=v1(1)+pl(1)*xt_2

      Sol_1(2)=v1(2)+pl(2)*xt_1
      Sol_2(2)=v1(2)+pl(2)*xt_2

      Sol_1(3)=v1(3)+pl(3)*xt_1
      Sol_2(3)=v1(3)+pl(3)*xt_2

*      write(*,*)'1st intersection co-ordinates:',(Sol_1(j),j=1,3)
*      write(*,*)'2nd intersection co-ordinates:',(Sol_2(j),j=1,3)

      if(abs(Sol_1(3))  .gt. hpix    ) then 
        write(*,*)'z1-inters outside cylinder, event',iev
        go to 999
      endif
      if(abs(Sol_2(3))  .gt. hpix    ) then 
        write(*,*)'z2-inters outside cylinder, event',iev
        go to 999
      endif

***      write(*,*)'1st: rpix uncer',(Sol_1(1)**2+Sol_1(2)**2-rpix**2)
***      write(*,*)'2nd: rpix uncer',(Sol_2(1)**2+Sol_2(2)**2-rpix**2)
      if(dabs(Sol_1(1)**2+Sol_1(2)**2-rpix**2).gt.0.0001d0) go to 999
      if(dabs(Sol_2(1)**2+Sol_2(2)**2-rpix**2).gt.0.0001d0) go to 999

      xDelta_1=(Sol_1(1)+Sol_2(1))/2.d0
      xDelta_2=(Sol_1(2)+Sol_2(2))/2.d0
      xDelta_3=(Sol_1(3)+Sol_2(3))/2.d0
*      write(*,*)'middle point: x,y,z',xDelta_1,xDelta_2,xDelta_3

      xDist=dsqrt((Sol_1(1)-Sol_2(1))**2+(Sol_1(2)-Sol_2(2))**2+ 
     +      (Sol_1(3)-Sol_2(3))**2)

      if (xDist.gt.2.d0*(dsqrt(rpix**2+hpix**2))) then 
      write(*,*)'LSP do not hit the inside ellipsoid event',iev ! Rapidity might be larger...
      go to 999
      endif
      
*      write(*,*)'1st way:t_1,t_2==>',xt_1,xt_2
* ** Here assuming 3D-geometry:
* **  plu.dot.R=0 (as the momentum vector is perpendicular to the min-distance)
* ** ...

      xDot_plu_r0=plu(1)*v1(1)+plu(2)*v1(2)+plu(3)*v1(3) 

*      write(*,*)'2nd way:t(only one: + or -ve)==>',xDot_plu_r0

      xVec_r_1=v1(1)-xDot_plu_r0*plu(1)
      xVec_r_2=v1(2)-xDot_plu_r0*plu(2)
      xVec_r_3=v1(3)-xDot_plu_r0*plu(3)

      xVec_r_mod_t_pos=dsqrt(xVec_r_1**2+xVec_r_2**2+xVec_r_3**2)

* *** By convention: t is positive.... so this is the correct values:
*      write(*,*)'UPDATE...'
*      write(*,*)'t:Positive:delta(x,y,z),mod(r)',xVec_r_1,xVec_r_2,
*     +           xVec_r_3,xVec_r_mod_t_pos
****************************************************
*      xVec_r_1=v1(1)+xDot_plu_r0*plu(1)
*      xVec_r_2=v1(2)+xDot_plu_r0*plu(2)
*      xVec_r_3=v1(3)+xDot_plu_r0*plu(3)
*      xVec_r_mod_t_neg=dsqrt(xVec_r_1**2+xVec_r_2**2+xVec_r_3**2)
*      write(*,*)'t:Negative:delta(x,y,z),mod(r)',xVec_r_1,xVec_r_2,
*     +           xVec_r_3,xVec_r_mod_t_neg
****************************************************
*      write(*,*)'decayLtmp,xDelta_1,xDelta_2,xDelta_3,xDist,
*     +xVec_r_mod_t_pos,x_y_Reco,x_eta_Reco'
*      write(*,*)decayLtmp,xDelta_1,xDelta_2,xDelta_3,xDist,
*     +xVec_r_mod_t_pos,x_y_Reco,x_eta_Reco
****************************************************

      write(12,108)iev,decayLtmp,xDelta_1,xDelta_2,xDelta_3,xDist      !1-6 
     +,xVec_r_mod_t_pos,x_boostLSP,x_bLSP,x_thetaLSP,x_yLSP,x_etaLSP   !7-12 
     +,x_cal_yLSP,x_cal_etaLSP,x_totE_LSP,x_magP_LSP,x_mass_LSP        !13-17
     +,pl_totE,xMomentum,x_mLSP,x_y_Reco,x_eta_Reco,xDelR_jj,xLSPmass  !18-24
     +,(xDist/decayLtmp)                                               !25

108   format(i9,1x,24(1x,e14.7))
12    format(28(1x,e14.7))  
24    format(18(1x,e14.7)) 
999   RETURN
      END
