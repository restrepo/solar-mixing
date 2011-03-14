C..... CROSS SECTION AND final N_events in fb 

      subroutine sigandnev 

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      include 'pylocalcom.h'
      double precision lum
      CHARACTER CHAP*16
      CHARACTER CHAPM*16
      common/sig/sig_sm,lum
      integer nevgen,nevgenwmu,nevgenwtau,nevgenwe,mucut(7),taucut(7)
      integer p1cut(0:7),p1hcut(0:7),p3cut(0:7),ecut(7)
      common/porodskands/nevgen,nevgenwmu,nevgenwtau,nevgenwe,nevgen3be,
     $     nevgen3bmu
      common/solar/nev_mutau,nev_muqq,nev_tauqq,nev_etau,nev_eqq
      common/tauisola/nevgenwtau1p,nevgenwtau3p,nevgenwtau1ph,
     $     nevgenwtaue,nevgenwtaumu,nevgenwtau4p
      common/eevt/neevt
      common/muevt/nmuevt
      common/aftercuts/ecut,mucut,taucut,p1cut,p1hcut,p3cut
      common/tauevt/ntauevt ! DR
************************************************************************
      integer subname(500,5)
      common/crossname/subname
**********************************************
      open(21,file='result.out',status='unknown') 
      open(22,file='aftercut.out',status='unknown') 

      xsolar_e = 0.d0
      xsolar_mu = 0.d0
      r_solar  = 0.d0

      write(21,*)'------------------------------------Event Input'
      write(21,*)'sqrt(s)=',rs,'GeV' 
      write(21,*)'Luminosity=',lum,'fb^-1' 
      if (mssmrp.eq.2) write(21,*)'interface with SPhenoRP.spc'
      if (npdflib.eq.0) write(21,*)'PDF defaults Pythia'
      write(21,*)'SM cross-section: in this analysis NO ',sig_sm
      write(*,*)'msel,nsub',msel,nsub
      write(21,*)'msel,nsub',msel,nsub
      write(21,*)'------------------------------------'
      write(21,*)'events: total Simulated, after cuts'
      write(21,*)ievmax,trilep

C***************CURRENTLY USED ******************************
      if(MSEL.eq.39)then
         sigmaall=0.d0
         rmainxsec=0.d0
         imainxsec=0
         do i=201,301 ! please note that although we have summed not all processes are ON
            sigmaall=sigmaall+xsec(i,3)
            if(xsec(i,3).gt.rmainxsec)then
               rmainxsec=xsec(i,3)
               imainxsec=i
            endif
            if(i.eq.230)then !  ! fi fj ->LSP_2+Chargino_1 production
               CALL PYNAME(KFPR(i,1),CHAP)
               CALL PYNAME(KFPR(i,2),CHAPM)
               write(21,*)'ixsection (fb) ',CHAP,'+ ',CHAPM,
     $              xsec(i,3)*1.E12
            endif
         enddo
         CALL PYNAME(KFPR(imainxsec,1),CHAP)
         CALL PYNAME(KFPR(imainxsec,2),CHAPM)
        write(21,*)'Main xsection (fb) ',CHAP,'+ ',CHAPM,
     $        rmainxsec*1.E12 ! mb -> pb(10^9) ->fb(10^12)
        write(21,*)'sigmaall (fb), efficiency:',sigmaall*1.E12,
     $        real(inev)/real(ievmax)
        write(21,*)'Xsection (fb):Total,Nev(as 2*chi_10)=2*Total*Lumi',
     $        sigmaall*1.E12,
     $        2.d0*sigmaall*1.E12*lum ! DR
C...BRAT(IDC) Is branching for IDC. IDC=1733   for chi_10 -> W-mu+
C...To obtain IDC: call pylist(12)
* ** spdas: Mutiplication with 0.6796 = Br(W->hadronic)
        write(21,*)'Nev for chi_10->e W,W->q qp, no cuts:'
        write(21,*)'BRAT(4412):LSP->eW possibly from SLHA',BRAT(4412)
***********************************************************************
* ** spdas: The BRAT(4412) is basically from the SPhenoRP.spc file. 
* ** However at the moment it is takeing care properly and this
* ** is the reason it gives us ZERO. However, the obtained numbers are ok.
        write(21,*)'expected',
     $        2d0*sigmaall*1.E12*lum*(BRAT(4412)*2d0)*0.6796,
     $        ' obtained',
     $        sigmaall*1.E12*lum*real(neevt)/real(ievmax)
        write(21,*)'Nev for chi_10->mu W,W->q qp, no cuts:'
        write(21,*)'expected', 
     $        2d0*sigmaall*1.E12*lum*(BRAT(4412)*2d0)*0.6796,
     $        ' obtained',
     $        sigmaall*1.E12*lum*real(nmuevt)/real(ievmax)
        write(21,*)'Nev for chi_10->tau W,W->q qp, no cuts:', ! DR
     $        sigmaall*1.E12*lum*real(ntauevt)/real(ievmax) ! DR
***********************************************************************
        write(21,*)'Below: GenNeve w FINAL CUTS(1-6) & w no TAU'
        write(21,*)'GenNevtau w FINAL CUTS(1-6) w res sphrd crossing'
        write(21,*)'GenNeve GenNevtau',
     $        sigmaall*1.E12*real(nevgenwe)/real(ievmax)*lum,
     $        sigmaall*1.E12*real(nevgenwtau)/real(ievmax)*lum
        write(21,*)'Below: GenNevmu w FINAL CUTS(1-6) & w no TAU'
        write(21,*)'GenNevmu,GenNevtau',
     $        sigmaall*1.E12*real(nevgenwmu)/real(ievmax)*lum,
     $        sigmaall*1.E12*real(nevgenwtau)/real(ievmax)*lum
        write(21,*)'GenNevtau1p,GenNevtau1ph,GenNevtaue,',
     $              'GenNevtaumu,GenNevtau3p,4p+ in below'
        write(21,129)sigmaall*1.E12*real(nevgenwtau1p)/real(ievmax)*lum,
     $        sigmaall*1.E12*real(nevgenwtau1ph)/real(ievmax)*lum,
     $        sigmaall*1.E12*real(nevgenwtaue)/real(ievmax)*lum,
     $        sigmaall*1.E12*real(nevgenwtaumu)/real(ievmax)*lum,
     $        sigmaall*1.E12*real(nevgenwtau3p)/real(ievmax)*lum,
     $        sigmaall*1.E12*real(nevgenwtau4p)/real(ievmax)*lum
        write(21,*)'GenNev3be,GenNev3bmu'
        write(21,*)sigmaall*1.E12*real(nevgen3be)/real(ievmax)*lum,
     $        sigmaall*1.E12*real(nevgen3bmu)/real(ievmax)*lum
        write(21,*)'------------------------------------'
      else
        write(21,*)'for msel=',msel,'use sigandevold'
C....for other cross specific cross section see sigandevold
      end if

      write(22,*)"Number of Events:e,mu,tau,tau-1p,tau-1p-Had,tau-3p"
      write(22,123)sigmaall*1.E12*lum*real(neevt)/real(ievmax),
     $           sigmaall*1.E12*lum*real(nmuevt)/real(ievmax),
     $           sigmaall*1.E12*lum*real(ntauevt)/real(ievmax),
     $           sigmaall*1.E12*lum*real(p1cut(0))/real(ievmax),
     $           sigmaall*1.E12*lum*real(p1hcut(0))/real(ievmax),
     $           sigmaall*1.E12*lum*real(p3cut(0))/real(ievmax)

      write(22,*)"Cut,NumberEvents:e,mu,tau,tau-1p,tau-1p-Had,tau-3p"
      do ii=1,7
         write(22,124)ii,sigmaall*1.E12*lum*real(ecut(ii))
     $              /real(ievmax),sigmaall*1.E12*lum*real(mucut(ii))
     $              /real(ievmax),sigmaall*1.E12*lum*real(taucut(ii))
     $              /real(ievmax),
     $              sigmaall*1.E12*lum*real(p1cut(ii))/real(ievmax),
     $              sigmaall*1.E12*lum*real(p1hcut(ii))/real(ievmax),
     $              sigmaall*1.E12*lum*real(p3cut(ii))/real(ievmax)
      enddo

*      write(*,*)'solar - number in the subroutine program:'
*      write(*,*)'Evt: mutau,muqq,tauqq',nev_mutau,nev_muqq,nev_tauqq
*      write(*,*)'Evt: etau,eqq,tauqq',nev_etau,nev_eqq,nev_tauqq

      xsolar_e  = nev_etau  - nev_eqq*0.159d0  - nev_tauqq*0.156
      xsolar_mu = nev_mutau - nev_muqq*0.159d0 - nev_tauqq*0.156

      r_solar= xsolar_e/xsolar_mu

      write(*,*)'xsolar_e,xsolar_mu,r_solar',xsolar_e,xsolar_mu,r_solar

123   format("TOTAL:",6(1x,f12.2))
124   format("AFTER CUT:",i3,1x,6(1x,'&',f12.2))
129   format("GenNev:tau1p,tau1ph,taue,taumu,tau3p,4p+",6(1x,f12.2))

      close(21)
      close(22) ! spdas

      return 
      end
********************************** 
