***********************************************************************
*********   CALCULATES THE 3-GENERATION R-PARITY VIOLATING  ***********
*********                SUSY CROSS SECTIONS                ***********
*********              LOOKING FOR E/TAU RATIO              ***********
***********************************************************************
C::WARNING: open files must have less than 8 characters (without extens)
C. DEBUG: smearing currently commented
C         trigger currently commented
C         comment jet isolation
C         use input/output in subroutines intead commons
      implicit double precision(a-h,o-z)
      include 'pylocalcom.h' 
      EXTERNAL PYDATA
*******************-- local variables --********************************
      integer i
      double precision vp(3)
      double precision pchifin(5,5,10)
      double precision vchi(10,4),vchim(10,4)
      double precision chirt(100)
      double precision v1(3)
      integer chiline(5)
      character*16 chiname(5,10)
      character cuts*3
      logical ldec,evento
      logical triggerps
      logical show_e_event,show_mu_event,show_tau_event,show
      logical tau_event,sphere_line(10)
      logical efficiency(5)
      integer chiorig(5),iichiorig,kfchi,iichi,chifin(5),kfchifin(5,10)
      integer nfchi(5,10),nprongmu(10),npronge(10)
      logical allfst
      double precision DRmin(20)
      integer ecut(7)
      integer kfdaug(10),nldaug(10),mucut(7),taucut(7),p1cut(0:7)
      integer nlnew(20),nprong(10),nprongh(10),p3cut(0:7),p1hcut(0:7)
      logical ldecjj(10),ldeceta(10),ldecmin(10),ldecmax(10)
      logical tauevt1p(5),tauevt3p(5)
      logical tauevt1ph(5),tauevt4p(5),tauevtmu(5),tauevte(5)
      logical tauevt1p3b(5),tauevt3p3b(5)
      logical tauevt1ph3b(5),tauevt4p3b(5),tauevtmu3b(5),tauevte3b(5)
      logical trigger,aux(7)
      double precision pwparton(4,10),drmintot(10)
      double precision DRimin(15),DRimax(15)
      logical enocut(10),munocut(10),taunocut(10)
      logical enocut3b(10),munocut3b(10),taunocut3b(10)
      logical usecut(13),ldecfix,auxcut,usecut3b
      logical porodskandsl
      double precision pps(3),ppsu(3)
      double precision massinvw(10)
      integer njetwjj(10)
      double precision az(10,5),etltau(10),pfulle(10,8)
      double precision  pfullmu(10,8),pfulltau(10,8),pfulljet1(10,8)
      double precision pfulljet2(10,8),pfullp1(10,8),pfullp2(10,8)
      double precision pfulle3b(10,8),pfullmu3b(10,8),pfulltau3b(10,8)
      double precision pfullp13b(10,8),pfullp23b(10,8),pfulltau2(10,8)
      double precision plep_hitcil(10,8),pjw_new(4,100) ! spdas
      double precision ppartonw(4,10),ptmpj(8,100),pfullchi(10,8)
      double precision pfullprng1(10,8),pfullprng2(10,8)
      double precision pfullprng3(10,8),vtau(10,4)
      double precision ltrack
      double precision pquark_nWdecay(10,12)
      logical finalcute,finalcute3b,finalcutem
      logical finalcuts(20),finalcuts3b,finalcutmu,finalcutmu3b
      logical finalcuttau,ptcut(5),rcut(5)
      integer nfinalcuts(20)
      logical fullwrite,partialcuts,finalcutmum,finalcuttaum,debug
**************__ Commons __**************
      double precision lum
      common/sig/sig_sm,lum
      integer nevgen,nevgenwmu,nevgenwtau,nevgenwe
      common/porodskands/nevgen,nevgenwmu,nevgenwtau,nevgenwe,nevgen3be,
     $     nevgen3bmu
      common/solar/nev_mutau,nev_muqq,nev_tauqq,nev_etau,nev_eqq 
      common/tauisola/nevgenwtau1p,nevgenwtau3p,nevgenwtau1ph,
     $     nevgenwtaue,nevgenwtaumu,nevgenwtau4p
      integer nffin(5,100)
      common/dauglines/nffin
      common/eevt/neevt
      common/muevt/nmuevt
      common/tauevt/ntauevt ! DR
      common/jetevt/njetevt !MBM
      common/eventnumber/iev,nchilinejj
      common/teste/ievv,ptcut,rcut     !MBM
      common/aftercuts/ecut,mucut,taucut,p1cut,p1hcut,p3cut

      common/kinematicsLSP/decayLtmp,decayLps,x_boostLSP,x_bLSP
     +       ,x_thetaLSP,x_yLSP,x_etaLSP,x_cal_yLSP,x_cal_etaLSP
     +       ,x_totE_LSP,x_magP_LSP,x_mass_LSP

      logical cluster
      common/debugcluster/cluster
c      common/jetprop/ptmpj(8,100)
c      common/wpartonlevel/ppartonw(4,10)
*******************-- general input --**********************************

      open(10,file='input_tev.dat',status='old') 
      read(10,*)ievmax          ! number of simulations
      read(10,*)mssmrp          ! option: 0 sm, 1 mssm, 2 expl. rpv
      read(10,*)msel            ! process: you: 0, or full mssm: 39
      read(10,*)npdflib         ! 0: default pythia, 1: default isajet 
      read(10,*)rs              ! dsqrt(s)=14000.d0 (LHC)
      read(10,*)cuts            ! cut selection: 'hc1', 'sc1', etc
      read(10,*)sig_sm          ! SM total cross section
      read(10,*)lum             ! luminosity
      read(10,*)nhisto          ! histograms: (1) yes; (0) no
c      read(10,*)jetcut          ! cut on the number of jets
      read(10,*)nseed            ! seed for pythia mrpy(1)
      close(10)
      if(nseed.eq.0)nseed=time()
C     1 prong tau: 1189718480
      mrpy(1)=nseed ! time() ! 25 seed for random number generation in pythia 25core 3-p 26 1-p
C...  PAW output: 
      morefilter=4 !-1: no output; 0: all, 1: W mu or tau; 2: W mu; 3 W tau; 4: j j mu or tau, 5: j j tau_prongs; 6: all but not mj mass cuts  ! spdas:  DR version it was 4. The fort.92 is generating because changing 5->4.
********************-- chosen processes --********************************
      write(91,*)mrpy(1) !nseed
c.....open subprocesses subpr(1)...subpr(nsub). 
         
      if(msel.eq.0)then
         call opensub           !read file isub.dat 
      end if

      write(12,*)'iev,decayLtmp,xDelta_1,xDelta_2,xDelta_3,xDist
     +,xVec_r_mod_t_pos,x_boostLSP,x_bLSP,x_thetaLSP,x_yLSP,x_etaLSP
     +,x_cal_yLSP,x_cal_etaLSP,x_totE_LSP,x_magP_LSP,x_mass_LSP
     +,pl_totE,xMomentum,x_mLSP,x_y_Reco,x_eta_Reco,xDelR_jj,xDistByLp'
c.....susy initilization: susy model and susy input parameters

c.....reading decay tables and parameters according to Les Houches

      imss(1)=11                ! switches on the les houches file
      imss(21)=10               ! read the table for couplings and mixings
      imss(22)=10               ! read the table for branchings
      if(mssmrp.eq.1)open(10,file='SPheno.spc',status='old') 
      if(mssmrp.eq.2)open(10,file='SPhenoRP.spc',status='old') 

c.....LHC initialization 

      call setlhc
c.....initialization for the lhc. 

      call pyinit('cms','p+','p+',rs)

c.....LEP constraints

      call lepcuts

c.....writing the final decay table

      call pyupda(1,97)         ! output in file fort.97

C...obtain IDC number
* ** However for the processes that are fed from the SLHA check those speacial numbers.
* **       call pylist(12) 
c      stop

C...write Branching for IDC BRAT(IDC). IDC=1733   for chi_10 -> W-mu+
c      write(*,*)"hola",BRAT(1733)
c.....pymaxi was called in pyinit. call again if any change

c      if(mssmrp.gt.0) call pymaxi 

*******************-- event generation --*******************************
C.....initialization of logical variables
      show_e_event=.False.
      show_mu_event=.False.
      show_tau_event=.False.
      show=.False.
      debug=.False.
      if(show_e_event.or.show_mu_event.or.show_tau_event)show=.True.
* ** For the default intialization (the three lines above) "show" is always False. spdas
      tau_event=.False.
      porodskandsl=.False.

* **       porodskandsl=.True. ! to extract the numbers from hep-ph/0401077
c.....initialization of counters
c.....event generation loop
      ievprint=100
      nevtcut=0
      nevgenwe=0
      nevgenwmu=0
      nevgen3be=0
      nevgen3bmu=0
      nevgenwtau=0
      nevgenwtau1p=0
      nevgenwtau1ph=0
      nevgenwtaue=0
      nevgenwtaumu=0
      nevgenwtau3p=0
      nevgenwtau4p=0
      decayL=0d0
      nenocut=0
      nenocut3b=0
      nmunocut=0
      nmunocut3b=0
      ntaunocut=0
      p1cut(0)=0
      p1hcut(0)=0
      p3cut(0)=0
      ncluster=0
* ** To extract the Solar Number:
      nev_mutau = 0
      nev_muqq  = 0
      nev_tauqq = 0

      nev_etau  = 0
      nev_eqq   = 0

      do ii=1,7
         ecut(ii)=0
         mucut(ii)=0
         taucut(ii)=0
         p1cut(ii)=0
         p1hcut(ii)=0
         p3cut(ii)=0
      enddo

      xlspmass         = 0.d0
      xLp              = 0.d0
      xLp_new          = 0.d0
      xabsP_LSP_reco   = 0.d0
      xP_LSP_reco_px   = 0.d0
      xP_LSP_reco_py   = 0.d0
      xP_LSP_reco_pz   = 0.d0

******************All distances in mm ********************
C.... Detector Geometry
      RRsili=0.012d0              ! transverse pixel resolution
      BBsili=0.077d0              ! longitudinal pixel resolution
c.... Testing the neutralino decay vertex in lhc
c.... spheroid x^2/aa^2+y^2/aa^2+z^2/ba^2=1
c.... Inner detector is 7m long and 1.15m radius see REFERENCES
c.... BUG: we cannot include the track detector system here!
c     aa=20d-3*5d0; ba=0.5d0*5d0 !LHC
c     aa=10d-3*5d0; ba=0.025d0*5d0 !Tevatron
      if(porodskandsl)then
         aa=0.02d0*5d0; ba=0.5d0*5d0 !LHC ! 
      else
        aa=0.02*5d0; ba=0.5d0*5d0 !LHC ! 
c        aa=RRsili*5d0; ba=BBsili*5d0 !LHC ! 
      end if
c     aa=0.012*5d0; ba=0.066d0*5d0 !ATLAS
      rtrack=900d0                ! 1.1m  !track system effective radius
      ltrack=3000d0               ! 7m wide.  3.5m dectector lenght from collision point
      trackmaxeta=2.5D0         ! see inner_detector.svg 1.5D0 ! maximum neutralino pseudorapidity
      
C.... ATLAS internal detector a) Pixel detector, b) Semiconductor tracker (SCT), 
C.... c) transition radiation tracker (TRT)
      
C.... Dimensions and resolutions R: recta; V: vertex; P_i: puntos
      rpix=48.d0                ! ratio of first pixel detector (FPD)
      hpix=380.d0               ! length of FPD. To be conservative (real hpix=400 mm)
C.... Information about leptona and tau isolation
      do i=1,3 ! e and mu
         Drimin(9+i*2)=0d0
         Drimax(9+i*2)=0.3d0
      end do
      Drimax(15)=0.4d0          ! tau annulus, outer radio DR
      Drimin(15)=0.1d0          ! tau annulus, inner radio DR
      ET_isola=5d0              ! maximum ET inside isolation cone (l) or annulus (tau)
      neevt=0
      nmuevt=0
      ntauevt=0
      numberofcuts=7
      do i=1,numberofcuts
         nfinalcuts(i)=0
      end do
C.... END Detector Geometry
      ievdebug=6167
      do 200 iev=1,ievmax
         ievv=iev
         if(debug.and.iev.eq.ievdebug)print*,iev
         if(mod(iev,1000).eq.0) write(*,*) 'now at event number',iev
c     do ii=1,dl
c.....Fim dos testes
         call pyevnt           
         

         if(debug.and.iev.eq.ievdebug)print*,iev,'pyevent'
         call pyedit(21) !Store a copy of the table!
         if(debug.and.iev.eq.ievdebug)print*,iev,'pyedit(21)'
         evento=.true.
         tau_event=.False.
C.....Intializarion of variables
C...  ntuple variables
c         if(iev.eq.ievprint)call pylist(3)
C.....HISTORY OF THE EVENT
         chiorig(1)=-1 ! not used for the cascade initiator
         iichiorig=1
         kfchi=1000022
c....neutralino decay output: iichi: number of kfchi found; 
c    chiline(1)...chiline(iichi) table position for each kfchi
c    chifin(1)...chifin(iichi) number of final state of each kfchi
c    chiname(1)...chiname(iichi)  and kfchifin(1)...kfchifin(iichi) are names and KF of final states
         do ii=1,5
            sphere_line(ii)=.False.                    
            tauevt1p(ii)=.False.; tauevt3p(ii)=.False.
            tauevt1ph(ii)=.False.; tauevt4p(ii)=.False.
            tauevtmu(ii)=.False.;tauevte(ii)=.False.
            tauevt1p3b(ii)=.False.; tauevt3p3b(ii)=.False.
            tauevt1ph3b(ii)=.False.; tauevt4p3b(ii)=.False.
            tauevtmu3b(ii)=.False.;tauevte3b(ii)=.False.
            enocut(ii)=.False.; munocut(ii)=.False.
            taunocut(ii)=.False.
            enocut3b(ii)=.False.; munocut3b(ii)=.False.
            taunocut3b(ii)=.False.
         end do
         triggerps=.False.
C... Obtain parton decay chain from informative lines: k(i,1)=21
C.... decaychainnew in decaychain.f  uses the full table         

         call decaychain(chiorig,iichiorig,kfchi,iichi,chiline,chifin,
     $        chiname,kfchifin,pchifin)
***********************************************************************
         if(debug.and.iev.eq.ievdebug)print*,iev,'decaychain'
         do ii=1,iichi
            do jj=1,chifin(ii)
               nfchi(ii,jj)=nffin(ii,jj)
            end do
         end do

c...  explanation of the subroutine for the neutralino case
         if(show)call writedecay(chiorig,iichiorig,kfchi,iichi,chiline,
     $        chifin,chiname,kfchifin,pchifin)
         if(iev.eq.51541)call writedecay(chiorig,iichiorig,kfchi,iichi,
     $        chiline,chifin,chiname,kfchifin,pchifin)
C...
C.....(II) ANALYSIS OF ~chi_10 DAUGHTERS: check if ~chi_10->W lepton
C....Roughly Porod Skands cuts
C...Determine if at least one neutralino decays into W
C   In such a case determines if the W is a good signal vertex
C      Signal vertex
C      a)Check that the event satisfy the trigger cuts
C      b) chi_10 -> W mu, or W tau, check that mu and tau are
C         charged leptons
C      c) W -> jj but not b bar. Check that j is in fact a jet
C      d)  
         do 201 jj=1,iichi      ! BEGIN loop over each LSP neutralino
            hpixmax=1d30
C...RESET important particle propierties
C...  pfullXXX(jj,1-4): momentum
C...  pfullXXX(jj,5): p_T
C...  pfullXXX(jj,6): eta
C...  pfullXXX(jj,7): good event 1, either 0
C...  pfullXXX(jj,8): theta
            do jjj=1,8
               pfullchi(jj,jjj)=-1d30 ! ridiculous large value.
               pfulle(jj,jjj)=-1d30 ! ridiculous large value.
               pfulle3b(jj,jjj)=-1d30 ! ridiculous large value.
               pfullmu(jj,jjj)=-1d30 ! ridiculous large value.
               pfullmu3b(jj,jjj)=-1d30 ! ridiculous large value.
               pfulltau(jj,jjj)=-1d30 ! ridiculous large value.
               pfulltau2(jj,jjj)=-1d30 ! ridiculous large value.
               pfulltau3b(jj,jjj)=-1d30 ! ridiculous large value.
               pfullp1(jj,jjj)=-1d30 ! ridiculous large value.
               pfullp13b(jj,jjj)=-1d30 ! ridiculous large value.
               pfullp2(jj,jjj)=-1d30 ! ridiculous large value.
               pfullp23b(jj,jjj)=-1d30 ! ridiculous large value.
               pfullprng1(jj,jjj)=-1d30
               pfullprng2(jj,jjj)=-1d30
               pfullprng3(jj,jjj)=-1d30
            end do
            pfullchi(jj,7)=0     ! ridiculous large value.
            pfulle(jj,7)=0     ! ridiculous large value.* spdas jjj=>7
            pfulle3b(jj,7)=0     ! ridiculous large value.* spdas jjj=>7
            pfullmu(jj,7)=0     ! ridiculous large value.
            pfullmu3b(jj,7)=0     ! ridiculous large value.
            pfulltau(jj,7)=0     ! ridiculous large value.
            pfulltau2(jj,7)=0     ! ridiculous large value.
            pfulltau3b(jj,7)=0     ! ridiculous large value.
            pfullp1(jj,7)=0     ! ridiculous large value.
            pfullp13b(jj,7)=0     ! ridiculous large value.
            pfullp2(jj,7)=0    ! ridiculous large value.
            pfullp23b(jj,7)=0    ! ridiculous large value.
            pfullprng1(jj,7)=0
            pfullprng2(jj,7)=0
            pfullprng3(jj,7)=0

           do kk=1,4
             do kkk=1,100
                pjw_new(kk,kkk) = 0.d0
             enddo
           enddo

            do jjj=1,3
               vtau(jj,jjj)=-1d30
            end do
            do jjj=1,4
               pfullchi(jj,jjj)=pyp(chiline(jj),jjj)
            end do
            pfullchi(jj,5)=pyp(chiline(jj),10)
            pfullchi(jj,6)=pyp(chiline(jj),19)
***            write(*,*)'flavor---', iabs(k(chiline(jj),2)) ! it is neutralino 
            pfullchi(jj,8)=pyp(chiline(jj),15)
C... check that at least one neutralino decays into mu W or tau W
C...  and the other neutralino have a proper jet or charged lepton
C... INPUT: chifin(jj),chiline(jj),kfchifin(jj,ii),show,

            if(chifin(jj).eq.2)then ! BEGIN TWO BODY DECAYS

               call chisignalvertex2(jj,chifin,chiline,kfchifin,nfchi,
     $     show,enocut,munocut,taunocut,pfulle,
     $              pfullmu,pfulltau,pfullp1,pfullp2)
               if(enocut(jj))nenocut=nenocut+1
               if(munocut(jj))nmunocut=nmunocut+1
               if(taunocut(jj))ntaunocut=ntaunocut+1

            elseif(chifin(jj).eq.3)then 
C.... Three body decays

               call chisignalvertex3(jj,chifin,chiline,kfchifin,nfchi,
     $     show,enocut3b,munocut3b,taunocut3b,pfulle3b,
     $              pfullmu3b,pfulltau3b,pfulltau2,pfullp13b,pfullp23b)
               if(enocut3b(jj))nenocut3b=nenocut3b+1
               if(munocut3b(jj))nmunocut3b=nmunocut3b+1

            end if
            if(debug.and.iev.eq.ievdebug)print*,iev,'chisig'
            if(show)write(*,*)'~chi_10->W(->jj) e, mu or tau',
     $           enocut(jj),munocut(jj),taunocut(jj)
C.....(III) CHECK if ~chi_10 DECAYS OUTSIDE MINIMAL ELLIPSOID AND 
C...END WARNING
C****************************************************
C...  We are now analyzing jj lightest neutralinos... 
C*****************************************************
C...  Simulating track system
C...  Check the psudorapidity of each neutralino
C.... track system pseudorapidity cut for lhc
C.... eta= -log(tan(theta/2)) -> theta=2*atan(exp(-eta))
C.... theta=90 -> eta ->0 ; theta=0 -> eta -> infinity <-> (eta=3)
C...  jj is the jj neutralino 
C...       chiline(jj) is the line in table of each neutralino
c.... decay width PMAS(PYCOMP(1000022),2), 
C.... below is the decay lenght in mm for each neutralino [(hbar*c/Gamma): PMAS(PYCOMP(1000022),4)]
c.... average at the end
            decayLtmp=PMAS(PYCOMP(1000022),4)*
     $           dsqrt(p(chiline(jj),4)**2/p(chiline(jj),5)**2-1.d0) 
            decayL=decayL+decayLtmp
            decayLps=PMAS(PYCOMP(1000022),4) ! without boost
            
          x_bLSP=dsqrt(pyp(chiline(jj),4)**2/pyp(chiline(jj),5)**2-1.d0)
            
            x_totE_LSP=p(chiline(jj),4)
            x_magP_LSP=pyp(chiline(jj),8)
            x_mass_LSP=p(chiline(jj),5)
            
            x_thetaLSP=pyp(chiline(jj),13)
            x_yLSP=pyp(chiline(jj),17)
            x_etaLSP= pyp(chiline(jj),19)
            
          x_cal_yLSP=0.5d0*(log((pyp(chiline(jj),4)+pyp(chiline(jj),3))/
     +           (pyp(chiline(jj),4)-pyp(chiline(jj),3))))
        x_cal_etaLSP=0.5d0*(log((pyp(chiline(jj),8)+pyp(chiline(jj),3))/
     +           (pyp(chiline(jj),8)-pyp(chiline(jj),3))))
            
*           write(*,*)'*************'
*           write(*,*)'for evt:',iev
*           write(*,*)'decayLtmp,totaL,Lps',decayLtmp,decayL,decayLps
*           write(*,*)'xboostLSP using P and PYP',xboostLSP,x_bLSP
***         x_boostLSP=dsqrt(p(chiline(jj),4)**2/p(chiline(jj),5)**2-1d0)
*           write(*,*)'theta,Y,eta-of LSP:',x_thetaLSP,x_yLSP,x_etaLSP
*           write(*,*)'Calculated Y,eta-of LSP:',x_cal_yLSP,x_cal_etaLSP
*           write(*,*)'*************'

        if (iev.le.3)write(*,*)'check:1st: iev,L0,Llab,Ltotal',
     +       iev,decayLps,decayLtmp,decayL
******************-- smearing all momenta --****************************
*******                  LHC energy resolution                **********
****************** if mstu(53)=0 then smear.f***************************
        
        if(mstu(53).ne.0) go to 345    

        do ijk=1,n
           if(k(ijk,1).lt.11)then
              if((iabs(k(ijk,2)).eq.11.or.iabs(k(ijk,2)).eq.13.or.
     $             iabs(k(ijk,2)).eq.15.or.iabs(k(ijk,2)).eq.22)
     $             .and.dabs(pyp(ijk,19)).lt.5.d0)then
                 call smear(ijk,0.1d0,0.01d0)
              else
                 if(dabs(pyp(ijk,19)).le.3)call smear(ijk,0.50d0,0.03d0)
                 if(dabs(pyp(ijk,19)).gt.3)call smear(ijk,1.d0,0.10d0)
              endif
           endif
        enddo

 345    continue
        
C...  To find the displaced vertex we will obtain the v of one of the 
C.... daughters
            do 2010 ii=1,n
C...kfchifin(jj,1): kf number of the first daughter of jj neutralino
C...nfchi(jj,1): line number in informative part of the table, 
C...             K(ii,1)=21, (v(ii,*)=0) of the first daughter of 
C...             jj neutralino
C...The first  daughter of jj neutralino reappear in the table but with
C...            K(ii,1)< or = 11 and v(ii,*) different from zero 
C... k(ii,3) is ... DEBUG
               if(k(ii,2).eq.kfchifin(jj,1).and.k(ii,3).eq.nfchi(jj,1))
     $              then
                  do i=1,3
c.... for one inestable particle v(ii,5) is different from zero
c.... v is the point where particle was produced and vp where it
c.... will decay.  v(ii,5)=0 if the particle is stable
C.... Therefore, vp should be obtained also from v of some of the 
C.....neutrino daughters 
                     vp(i) = v(ii,i)
                     vchi(jj,i)=vp(i)
                     vchim(jj,i)=vp(i)
                     rndmd=sqrt(vp(1)**2+vp(2)**2+vp(3)**2)/decayLtmp ! obtain random distribution
                  end do
               end if
 2010       CONTINUE

            if(porodskandsl)then ! don't use the boost
               do ii=1,3
                  pps(ii)=p(chiline(jj),ii) ! obtain neutralino momentum
               enddo 
               call univec(pps,ppsu)
               decayLps=rndmd*decayLps
               do ii=1,3
                  vp(ii)=decayLps*ppsu(ii)
                  vchi(jj,ii)=vp(ii)
               end do
            end if
            if(show)write(*,*)'vp of neutralino n',jj,':',vp
c.... (IV) CHECKING IF TRACKS HIT ENOUGH PIXEL DETECTORS 
            do i=1,3
               az(jj,i)=10000d0
            end do
            hitcil: if(enocut(jj).or.munocut(jj).or.taunocut(jj))then 
               if(show)then
                  write(*,*)'check if the tracks Hit enough pixel det.'
               end if
C...   Check if neutralino daughters at parton level hit enough pixel detectors
               do ii=1,chifin(jj)
                  kfdaug(ii)=kfchifin(jj,ii) ! kf of W and lepton in ~chi_10 -> W lepton 
                  nldaug(ii)=nfchi(jj,ii)
               end do
               nchidaug=chifin(jj)   ! =2: number of ~chi_10 daugthers in  W and lepton
               do jjj=1,3
                  v1(jjj)=vchi(jj,jjj) ! displaced vertex position
               enddo
C... Obtain line numbers in table for final states of neutralino at parton level
C...   eg: line numbers for q q' lepton in ~chi_10 -> W lepton -> q q' lepton

               call chipartonline(nchidaug,kfdaug,nldaug,ndaugnew,nlnew)

***               write(*,*)'check:daughters...',nlnew(1),nlnew(2),nlnew(3)

               if(debug.and.iev.eq.ievdebug)print*,iev,
     $              'chipartonline'
C... OUTPUT: ndaugnew: number of parton level neutralino daughters
C...         nlnew(1) ... nlnew(ndaugnew): line in table for q q' and lepton
               do i=1,ndaugnew ! loop over neutralino parton level daugthers
                  iline=nlnew(i)
c                  call hitcilinder(iline,v1,rpix,hpixmax)
                  hpix=0.9d0*hpixmax ! deactivate this cut
                  az(jj,i)=hpixmax
                  if(show)then 
                     write(*,*)'abs(hpixmax)',ndaugnew
                  end if
                  if(abs(hpixmax)/abs(hpixmax).ne.1)then
                     write(*,*)'ERROR, NAN displaced vertex position'
                     stop
                 end if
               end do
C.... Parton level invariant mass of W.
               npwparton=0
               do i=1,ndaugnew  ! loop over neutralino parton level daugthers
                  if(abs(k(nlnew(i),2)).le.6)then
                     npwparton=npwparton+1
                     do jjj=1,4
                        pwparton(jjj,npwparton)=p(nlnew(i),jjj)
                     end do
                  end if
               end do
               rminvtotparton=DInvariantMass(npwparton,pwparton)
***********************************************************************************
            end if hitcil

 201     CONTINUE               ! END loop over each LSP neutralino

         if(show)then
            if(iev.eq.ievprint)call pylist(3)
         end if

****************************************
C.... Switch on fragmentation:
         call pyexec

         nniev=80 !13830
         if (iev.eq.nniev) call pylist(3) !spdas: to get the BRAT or IDC number

         if(debug.and.iev.eq.ievdebug)print*,iev,'pyexec'
****************************************
         do 204 jj=1,iichi
C...INPUT: jj,chifin(jj),kfchifin(jj,ii),nwtau(jj),nfchi(jj,ii),show
            if(taunocut(jj))then
               call tauanalysis2(jj,chifin,kfchifin,nfchi,vchi,
     $              RRsili,BBsili,show,show_tau_event,efficiency,
     $              sphere_line,tau_event,nprong,nprongh,npronge,
     $              nprongmu,pfullprng1,pfullprng2,pfullprng3,vtau)
              if(nprong(jj).eq.1)then
                  tauevt1p(jj)=.True.
                  if(nprongh(jj).gt.0)tauevt1ph(jj)=.True.
                  if(npronge(jj).gt.0)tauevte(jj)=.True.
                  if(nprongmu(jj).gt.0)tauevtmu(jj)=.True.
               else if(nprong(jj).eq.3.and.ptcut(jj).and.rcut(jj))then   !MBM
                  tauevt3p(jj)=.True.
               else
                  tauevt4p(jj)=.True.
               end if
            end if              ! end if over neutralino: jj
            if(enocut3b(jj).or.munocut3b(jj))then
               call tauanalysis2(jj,chifin,kfchifin,nfchi,vchi,
     $              RRsili,BBsili,show,show_tau_event,efficiency,
     $              sphere_line,tau_event,nprong,nprongh,
     $              npronge,nprongmu,pfullprng1,pfullprng2,
     $              pfullprng3,vtau)
              if(nprong(jj).eq.1)then
                  tauevt1p3b(jj)=.True.
                  if(nprongh(jj).gt.0)tauevt1ph3b(jj)=.True.
                  if(npronge(jj).gt.0)tauevte3b(jj)=.True.
                  if(nprongmu(jj).gt.0)tauevtmu3b(jj)=.True.
               else if(nprong(jj).eq.3.and.ptcut(jj).and.rcut(jj))then   !MBM
                  tauevt3p3b(jj)=.True.
               else
                  tauevt4p3b(jj)=.True.
               end if
            end if              ! end if over neutralino: jj
  204     CONTINUE

         if(debug.and.iev.eq.ievdebug)print*,iev,'tauanalysis'
C... Aislamniento de leptones y taus
         do jj=1,iichi
            etltau(jj)=10000d0
            if(enocut(jj).or.munocut(jj).or.taunocut(jj).or.enocut3b(jj)
     &         .or.munocut3b(jj))then ! end if in 6666
               do ii=1,chifin(jj)
                  kfdaug(ii)=kfchifin(jj,ii)
                  nldaug(ii)=nfchi(jj,ii)
               end do
               nchidaug=chifin(jj)
               call isola(nchidaug,kfdaug,nldaug,DRimin,DRimax,ET,pt)

               if(show)write(*,*)'E_T=',ET
c this cut replace efficiencies en leptons and taus!
               etltau(jj)=et
            end if
         end do
         if(debug.and.iev.eq.ievdebug)print*,iev,'isola'
C...OUTPUT
C...  whether isolated or no
         njet=0
         call pycell(njet) ! find the number of jets (see setlhc.f)
         if(debug.and.iev.eq.ievdebug)print*,iev,'pycell'
C.... (IX) Implement of the cuts
C.....Full cuts:
C...  (V)  TRIGGERS CUTS FROM POROD & SKANDS
         if(show)write(*,*)'***BEGIN TRIGGERS*****'
         call triggers_old(triggerps)
         if(show)write(*,*)'***END TRIGGERS*****'
         if(show)write(*,*)'****EVENT NUMBER',iev,' END HERE*******'
C... Drimin and Drimax was calculated before by call isola
c         call  triggersdv(trigger,Drimin,Drimax,njet,
c     $        ET_isola)
c         print*,'trigger',trigger
c         trigger=.False.
         ntrigger=0
         call trigger_lhc(njet,trigger) ! trigger initialized .False. inside

         if(debug.and.iev.eq.ievdebug)print*,iev,'triggers'
         if(trigger)ntrigger=1
c         print*,trigger
         do 208 jj=1,iichi
         if(debug.and.iev.eq.ievdebug)print*,iev,'208',jj,iichi
c... reset paw variables
            massinvw(jj)=-1d0
            njetwjj(jj)=-1
            drmintot(jj)=-1d0
            do jjj=1,6
               pfulljet1(jj,jjj)=-1e-30
               pfulljet2(jj,jjj)=-1e-30
               plep_hitcil(jj,jjj)=0.d0 ! spdas

                do kjj=1,10
                   pquark_nWdecay(kjj,jjj)=0.d0
                   pquark_nWdecay(kjj,9)=0.d0
                   pquark_nWdecay(kjj,10)=0.d0
                   pquark_nWdecay(kjj,11)=0.d0
                   pquark_nWdecay(kjj,12)=0.d0
               end do

            end do

            pfulljet1(jj,7)  = 0 ! good event integer values
            pfulljet2(jj,7)  = 0 ! good event integer values
            plep_hitcil(jj,7)= 0 ! good event integer values

            pfulljet1(jj,8)   = 0.d0
            pfulljet2(jj,8)   = 0.d0
            plep_hitcil(jj,8) = 0.d0  

C... BEGIN DEBUG
C...  displaced vertex of each neutralino:
            do jjj=1,3
               v1(jjj)=vchi(jj,jjj)
            enddo
C... Invariant mass from all stable particles crossing displaced vertex
C... resolution sphere
c            call invmassfromall(v1,RRsili,minvchi) ! core dumped!
******************************************************
            do ii=1,chifin(jj)
               kfdaug(ii)=kfchifin(jj,ii)
               nldaug(ii)=nfchi(jj,ii)
            end do
C,,, Invariant mass of the event chi_10 -> (mu or tau) W
            imevt: if(enocut(jj).or.munocut(jj).or.taunocut(jj))then 
               nchidaugjj=chifin(jj)
C... DEBUG: not used ******************************
C...  invariant mass of each part and the event :
               allfst=.True.
c               call invariantmass(nchidaugjj,kfdaug,nldaug,v1,
c     $                 allfst,rminvpar,rminvtot)
         if(debug.and.iev.eq.ievdebug)print*,iev,'invariantmass',jj
C... DEBUG ******************************
               allfst=.True.
c               write(98,*)rminvpar,rminvtot
               call invariantmassjet(nchidaugjj,kfdaug,nldaug,v1,
     $              allfst,njet,rminvtotj)

         if(debug.and.iev.eq.ievdebug)print*,iev,'invariantmassjet'
               jetwreconstruction: if(njet.gt.0)then
C...********* invariantmassjetW ********
C... Calculates the jets by using only the W daughters
C... INPUT 
C...    nchidaugjj: number of neutralino daughters
C...    kfdaug(*): KF code of neutralino daughters
C...    nldaug(*): line number in table of neutralino daughters
C...    allfst: wheter consider only charged tracks (.false.) or all
C...    njet: number of jets in the event
C...    ppartonw: momentum info about partons

                  ppartonw(1,1)=pfullp1(jj,8) ! phi
                  ppartonw(2,1)=pfullp1(jj,6) ! eta

                  ppartonw(1,2)=pfullp2(jj,8) ! phi
                  ppartonw(2,2)=pfullp2(jj,6) ! eta

                  if(debug.and.iev.eq.ievdebug)print*,iev,jj

         if(iev.eq.51541)call pylist(1)
                  call invariantmassjetW(nchidaugjj,kfdaug,nldaug,v1,
     $                 allfst,njet,ppartonw,njetw,DRmin,rminvtotjnew,
     $                 rminvtotjw,ptmpj)

                  if(debug.and.iev.eq.ievdebug)print*,
     $                 iev,'invariantmassjetW',jj

****************************************
* ** NEW START  with 1-lepton + 2-jet
* ** spdas take the jet information from here...
                  do kkk=1,4
                     pjw_new(kkk,1)= ptmpj(kkk,1)
                     pjw_new(kkk,2)= ptmpj(kkk,2)
                  end do

*      write(*,*) 'check:LSP daughters.lt.1 and no call hitcil'
              if (ndaugnew.lt.1) go to 234

      nWdecay  = 0          ! W parton decay product...
      xLSPmass = 0.d0
      xWmass   = 0.d0

           do i=1,ndaugnew ! loop over neutralino parton level daugthers
                  iline=nlnew(i)

              nlepsel_sec=0
              plep_theta=0.d0
***      write(*,*)'check:LSPs daughters',i,'line,flav',iline,k(iline,2)

* ** taking only muon...
              if(abs(k(iline,2)).eq.13) then ! spdas:added abs
                 nlepsel_sec=nlepsel_sec+1
                 il_lep=iline

                do kk=1,4
                   plep_hitcil(nlepsel_sec,kk)=p(il_lep,kk)
                enddo
                   plep_hitcil(nlepsel_sec,5)=pyp(il_lep,10) !p_T
                   plep_hitcil(nlepsel_sec,6)=pyp(il_lep,19) !eta
                   plep_hitcil(nlepsel_sec,8)=pyp(il_lep,15) !phi
                   plep_theta=pyp(il_lep,13) !theta
*      write(*,*)'check:LSPs lep:',(plep_hitcil(nlepsel_sec,kk),kk=1,4)
              endif

*      write(*,*)'check:LSPs daughters:',i,'il,flav',iline,k(iline,2)
              if(abs(k(iline,2)).le.6) then ! W decay daughters
                     nWdecay=nWdecay+1
*      write(*,*)'check:LSPsinside le6:',i,'il,flav',iline,k(iline,2)
                        do jjj=1,4
                           pquark_nWdecay(nWdecay,jjj)=p(nlnew(i),jjj)
*      write(*,*)'check:q-',nWdecay,'4mom', pquark_nWdecay(nWdecay,jjj)
                        enddo
                           pquark_nWdecay(nWdecay,5)=pyp(nlnew(i),10) !p_T
                           pquark_nWdecay(nWdecay,6)=pyp(nlnew(i),19) !eta
                           pquark_nWdecay(nWdecay,8)=pyp(nlnew(i),15) !phi
                           pquark_nWdecay(nWdecay,12)=pyp(nlnew(i),13) !theta
*      write(*,*)'check:q-',nWdecay,'5,6,8,12', pquark_nWdecay(nWdecay,5),
*     + pquark_nWdecay(nWdecay,6),pquark_nWdecay(nWdecay,8),
*     + pquark_nWdecay(nWdecay,12)

              endif ! W decay daughters

            enddo ! loop over neutralino parton level daugthers

      if( nlepsel_sec .lt. 1) go to 234 ! at least one muon 
      if( nWdecay     .lt. 2) go to 234 !W decay products at least parton level-2 particle 
      if( njetw       .lt. 2) go to 234 ! at least two jet

      xLSPmass = (pjw_new(4,1) + pjw_new(4,2) + plep_hitcil(1,4))**2
     +         - (pjw_new(3,1) + pjw_new(3,2) + plep_hitcil(1,3))**2
     +         - (pjw_new(2,1) + pjw_new(2,2) + plep_hitcil(1,2))**2
     +         - (pjw_new(1,1) + pjw_new(1,2) + plep_hitcil(1,1))**2

      xWmass =   (pjw_new(4,1) + pjw_new(4,2))**2
     +         - (pjw_new(3,1) + pjw_new(3,2))**2
     +         - (pjw_new(2,1) + pjw_new(2,2))**2
     +         - (pjw_new(1,1) + pjw_new(1,2))**2

      if ( xLSPmass .lt. 0.d0) go to 234
      if ( xWmass   .lt. 0.d0) go to 234

      xLSPmass=dsqrt(xLSPmass)
      xWmass=dsqrt(xWmass)

      xP_LSP_reco_px= pjw_new(1,1) + pjw_new(1,2) + plep_hitcil(1,1)
      xP_LSP_reco_py= pjw_new(2,1) + pjw_new(2,2) + plep_hitcil(1,2)
      xP_LSP_reco_pz= pjw_new(3,1) + pjw_new(3,2) + plep_hitcil(1,3)

      xabsP_LSP_reco= dabs(dsqrt(xP_LSP_reco_px**2 +
     +                           xP_LSP_reco_py**2 +
     +                           xP_LSP_reco_pz**2   ))

      if(xabsP_LSP_reco .lt. 0.0000001) go to 234

      xLp=PMAS(PYCOMP(1000022),1)*decayLtmp/xabsP_LSP_reco
      xLp_new = decayLtmp/xabsP_LSP_reco

*      write(*,*)'mass-input,reco,L,absP,xLp,xLp_new',
*     +  PMAS(PYCOMP(1000022),1),xLSPmass,decayLtmp,xabsP_LSP_reco,
*     +  xLp,xLp_new
********************************************************
*            decayLtmp=PMAS(PYCOMP(1000022),4)*
*     $           dsqrt(p(chiline(jj),4)**2/p(chiline(jj),5)**2-1d0)
* PMAS(PYCOMP(1000022),4): is the rest-frame decay length.
* p(chiline(jj),{4,5}): (E,m)
* decayLtmp= is the boosted decay length.
* Thus ==>  L_lab= L_0 \times (p \over m)
********************************************************
********************************************************
*           write(*,*)'check:in the main before calling....'
*           write(*,*)'LSP daughters Nlep',nlepsel_sec
*           write(*,*)'lep1-mom',(plep_hitcil(1,kk),kk=1,4)
*           write(*,*)'LSP daughters Njet',njetw
*           write(*,*)'jet1-mom',(pjw_hitcil(kk,1),kk=1,4)
*           write(*,*)'jet2-mom',(pjw_hitcil(kk,2),kk=1,4)
*           write(*,*)'jet3-mom',(pjw_hitcil(kk,3),kk=1,4)
*           write(*,*)'jet4-mom',(pjw_hitcil(kk,4),kk=1,4)
*           write(*,*)'jetNjetW-mom',(pjw_hitcil(kk,njetw),kk=1,4)
*  pjw_new  is replaced for  pjw_hitcil
********************************************************
           xVec_r_mod_t_pos=0.d0

* ** We required at least 2-jets and 1-lepton before calling this "hitcilinder_update"

           call hitcilinder_update(iev,plep_hitcil,pjw_new,v1,
     +rpix,hpix,xVec_r_mod_t_pos,plep_theta,pquark_nWdecay)

*      write(*,*)'check:after hitcil:L,del(r)',decayLtmp,xVec_r_mod_t_pos

234   continue
      if(cluster)then 
         ncluster=ncluster+1
         evento=.False.
      end if

C... OUTPUT:
C...    njetw: number fo jets from pythia table with only W daughters 
C...    DRmin(1)...DRmin(njetw): Minimum distance between each W jet with
C...                             the others jet
C...    rminvtotjnew: invariant W mass calculated from the full jet system
C...    rminvtotjw: invariant W mass calculated from the jets coming only from W system
C...    ptmpj: momentum information about jets
C...CUTS for jets
                  drmintot(jj)=1000d0
                  massinvw(jj)=rminvtotjnew
                  njetwjj(jj)=njetw
                  do jjj=1,njetw
                     if(DRmin(jjj).lt.drmintot(jj))
     $                    drmintot(jj)=DRmin(jjj)
                  enddo
c               if(rminvtotjnew.gt.60d0.and.rminvtotjnew.lt.100d0
c     $              .and.drmintot(jj).gt.0.8d0)invmasjetw(jj)=.True.
C...Jet propiertes
                  do jjj=1,8
                     pfulljet1(jj,jjj)=ptmpj(jjj,1)
                     pfulljet2(jj,jjj)=ptmpj(jjj,2)
                  end do
               else
                  rminvtotjnew=-1d0
                  massinvw(jj)=rminvtotjnew
                  njetwjj(jj)=0
                  do jjj=1,8
                     pfulljet1(jj,jjj)=0d0
                     pfulljet2(jj,jjj)=0d0
                  end do
               end if jetwreconstruction
            end if imevt
 208     CONTINUE
         if(debug.and.iev.eq.ievdebug)print*,iev,'208 end'

C...   We accept events with two signal vertex!
C... DEBUG:  Add electron counter!
C... preparing ldec cut:
         ldec=.True.
         do jj=1,iichi
            ldecjj(jj)=.False.
            ellips2=vchi(jj,1)**2/aa**2+vchi(jj,2)**2/aa**2
     $           +vchi(jj,3)**2/ba**2
            if((ellips2.gt.1) !     minimum ellipsoid
     $           .and.(dsqrt(vchi(jj,1)**2+vchi(jj,2)**2).lt.rtrack) ! inside inner detector
     $           .and.(dabs(vchi(jj,3)).lt.ltrack))then ! inside inner detector
c     $           .and.dabs(pfullchi(jj,6)).lt.trackmaxeta)then ! small enough pseudorapidity
               ldecjj(jj)=.True.
            end if
            if(ldecjj(jj))pfullchi(jj,7)=1
            if(.not.ldecjj(jj))ldec=.False.
         end do
C...defining the cuts 
* ** spdas: The selection cuts are used for each Neutralino in the events. 
* ** We are seeing that the selection cuts are also loop with the variables "jj".

         do 213 jj=1,iichi
            do jjj=1,numberofcuts 
               finalcuts(jjj)=.False.
            end do
            finalcuts3b=.False.

            nusecut=8 ! 0 no cuts ! 8 full cuts ! 12 porod-skands cuts: set porod_skands=.True.
* ** if 7 then  usecut(7) becomes the usecut(0)=presence of e/mu/tau in the event.
C************ Final cuts ***********
C...   1) neutralino properties: minimum ellipsoid, inner detector and psuedorapidity
            finalcuts(1)=ldecjj(jj)
C...   2) Level 1 triggers
            finalcuts(2)=trigger
            finalcute=int(pfulle(jj,7)).eq.1
     $           .and.pfulle(jj,5).gt.20.
     $           .and.abs(pfulle(jj,6)).lt.trackmaxeta
            finalcute3b=int(pfulle3b(jj,7)).eq.1
     $           .and.pfulle3b(jj,5).gt.20.
     $           .and.abs(pfulle3b(jj,6)).lt.trackmaxeta
            finalcutmu=int(pfullmu(jj,7)).eq.1
     $           .and.pfullmu(jj,5).gt.20.
     $           .and.abs(pfullmu(jj,6)).lt.trackmaxeta
            finalcutmu3b=int(pfullmu3b(jj,7)).eq.1
     $           .and.pfullmu3b(jj,5).gt.20.
     $           .and.abs(pfullmu3b(jj,6)).lt.trackmaxeta
            finalcuttau=int(pfulltau(jj,7)).eq.1
     $           .and.pfulltau(jj,5).gt.20.
     $           .and.abs(pfulltau(jj,6)).lt.trackmaxeta
C...  3) neutralino two body decays to W e or W mu or W tau. mu or tau with p_T>20 GeV and eta<2.5 
            finalcuts(3)=(finalcute).or.(finalcutmu).or.(finalcuttau)
            finalcuts3b=(finalcute3b).or.(finalcutmu3b)
C...  4) check if parton level neutralino daughters:
C...      (~chi1_0 -> (mu or tau) q q'), hit enough pixel detectors
            finalcuts(4)=dabs(pfulljet1(jj,6)).lt.trackmaxeta
     $           .and.dabs(pfulljet2(jj,6)).lt.trackmaxeta !dabs(az(jj,1)).lt.hpix
c     $           .and.dabs(az(jj,2)).lt.hpix
c     $           .and.dabs(az(jj,3)).lt.hpix
C...  5) check mu or tau isolation, E_T<ET_isola GeV inside annulus DR
            finalcuts(5)=etltau(jj).lt.ET_isola
C...  6) invariant mass of W from two jets 
            finalcuts(6)=njetwjj(jj).eq.2.and.
     $           massinvw(jj).gt.60..and.massinvw(jj).lt.100.
C...  7) tau events with p_prong pounting out to neutralino and nprong 1 or 3
            finalcuts(7)=taunocut(jj).and.sphere_line(jj)
     $           .and.nprong(jj).le.3
C...Counters
            do jjj=1,numberofcuts
               if(finalcuts(jjj))nfinalcuts(jjj)=nfinalcuts(jjj)+1
            end do
            aux(1)=finalcuts(1)
            aux(2)=finalcuts(1).and.finalcuts(2)
            aux(3)=finalcuts(1).and.finalcuts(2).and.finalcuts(3)
            aux(4)=finalcuts(1).and.finalcuts(2).and.finalcuts(3).and.
     $             finalcuts(4)
            aux(5)=finalcuts(1).and.finalcuts(2).and.finalcuts(3).and.
     $             finalcuts(4).and.finalcuts(5)
            aux(6)=finalcuts(1).and.finalcuts(2).and.finalcuts(3).and.
     $             finalcuts(4).and.finalcuts(5).and.finalcuts(6)
            aux(7)=finalcuts(1).and.finalcuts(2).and.finalcuts(3).and.
     $             finalcuts(4).and.finalcuts(5).and.finalcuts(6).and.
     $             finalcuts(7)
            
            if(tauevt1p(jj))p1cut(0)=p1cut(0)+1    ! 1-prong tau decay
            if(tauevt1ph(jj))p1hcut(0)=p1hcut(0)+1 ! 1-prong tau decay hadronically 
            if(tauevt3p(jj))p3cut(0)=p3cut(0)+1    ! 3-prong tau decay 

* ** spdas: "ecut(N)" is the cuts for presence of electron (at least one? check ) 
* ** from jj^th Neutralino with the "aux(N)" cuts. Similarly, for the 
* ** other cuts, like muon, tau, tau1p,tau-1pHadronic, tau-3prong.

            do iiii=1,7
               if(enocut(jj).and.aux(iiii))ecut(iiii)=ecut(iiii)+1 
               if(munocut(jj).and.aux(iiii))mucut(iiii)=mucut(iiii)+1
               if(taunocut(jj).and.aux(iiii))taucut(iiii)=taucut(iiii)+1
               if(tauevt1p(jj).and.aux(iiii))p1cut(iiii)=p1cut(iiii)+1
               if(tauevt1ph(jj).and.aux(iiii))p1hcut(iiii)=
     $            p1hcut(iiii)+1
               if(tauevt3p(jj).and.aux(iiii))p3cut(iiii)=p3cut(iiii)+1
            enddo
* ** spdas: As events with e/mu do not have the the tau-prong associated 
* ** with it; thus applying aux(7) is meaningless and this is the reason 
* ** finally removed the finalcuts(7) from those events.
* ** check how ecut(6) and ecut(7) values in the output
            if(enocut(jj).and.(.not.aux(7)))ecut(7)=ecut(6)
            if(munocut(jj).and.(.not.aux(7)))mucut(7)=mucut(6)
C***********************************
            ldecfix=ldecjj(jj) !ldec ! ldecjj(jj)
            auxcut=taunocut(jj) 
            usecut(0)=munocut(jj).or.taunocut(jj).or.enocut(jj)
            usecut(1)=finalcuts(1)  ! small efficiency for ldec
            usecut(2)=finalcuts(2)  ! high efficiency
            usecut(3)=finalcuts(3)  ! replaces usecut(0). DEBUG: include e ! high efficiency
            usecut(4)=finalcuts(4)  ! high efficiency
            usecut(5)=finalcuts(5)  ! small efficiency, and, -> nmu >ntau
            usecut(6)=finalcuts(6)  ! small efficency. -> nmu > ntau 
            if(nusecut.eq.7)then 
               auxcut=finalcuts(7)  !tau instersecton of neutralino sphere
               usecut(nusecut)=usecut(0) ! small efficiency for tau: -> nmu >> ntau 
            end if
            if(nusecut.eq.8)then ! Full cut
               usecut(nusecut)=finalcuts(1).and.finalcuts(2)
     $              .and.finalcuts(3).and.finalcuts(4)
     $              .and.finalcuts(5).and.finalcuts(6)
               auxcut=finalcuts(7) ! tau daughters from neutralino resolution sphere
               usecut3b=finalcuts(1).and.finalcuts(2)
     $              .and.finalcuts3b.and.finalcuts(5)
     $              .and.finalcuts(7) 
            end if
C.... DEBUG CUTS: ************************************+
            if(nusecut.eq.9)then  ! porod_skands cuts without efficiences
               usecut(nusecut)=finalcuts(3).and.finalcuts(1)
     $              .and.finalcuts(4).and.finalcuts(2) 
               auxcut=finalcuts(7) ! tau daughters from neutralino resolution sphere
            end if
            if(nusecut.eq.10)then  ! porod_skands cuts without efficiences
               usecut(nusecut)=finalcuts(3).and.finalcuts(1)
     $              .and.finalcuts(4).and.finalcuts(5).and.finalcuts(2)
               auxcut=finalcuts(7) ! tau daughters from neutralino resolution sphere
            end if
            usecut(11)=finalcuts(3).and.finalcuts(5)
            if(nusecut.eq.12)then  ! porod_skands cuts without efficiences
               usecut(nusecut)=ldec
     $              .and.finalcuts(3)
     $              .and.finalcuts(4)
     $              .and.triggerps
               auxcut=taunocut(jj) ! tau daughters from neutralino resolution sphere
            end if
C.... DEBUG CUTS: ************************************+

            cutscheck: if(usecut(nusecut))then  !finalcuts(1-6)
               if(auxcut)then   ! if tau event
                  nevgenwtau=nevgenwtau+1
                  if(tauevt1p(jj))then   ! tau-1
                     nevgenwtau1p=nevgenwtau1p+1
                     if(tauevt1ph(jj))nevgenwtau1ph=nevgenwtau1ph+1
                     if(tauevte(jj))nevgenwtaue=nevgenwtaue+1
                     if(tauevtmu(jj))nevgenwtaumu=nevgenwtaumu+1
                  end if ! tau-1
                  if(tauevt3p(jj))nevgenwtau3p=nevgenwtau3p+1
                  if(tauevt4p(jj))nevgenwtau4p=nevgenwtau4p+1
C...  DEBUG: try separatly 1 and 3 prongs
C...  DEBUG: do not consider 1 prong leptonic tau????
c     print*,'tau event: prong, 1p,3p',nprong(jj),
c     $                    tauevt1p(jj),tauevt3p(jj)
               else if(enocut(jj))then ! else if tau event: it satisfies 1-6 finalcuts but _not_ auxcut.
                  nevgenwe=nevgenwe+1 ! e event with finalcuts(1-6)
               else if(munocut(jj))then 
                  nevgenwmu=nevgenwmu+1 ! mu event with finalcuts(1-6)
C... DEBUG. include electron counter
               end if  ! if tau event
C..DEBUG
               do ii=1,3
                  vp(ii)=vchim(jj,ii) ! is that necesssary? 
               end do 
            end if cutscheck

            cut3b: if(usecut3b)then
               if(enocut3b(jj))nevgen3be=nevgen3be+1
               if(munocut3b(jj))nevgen3bmu=nevgen3bmu+1
            endif cut3b

* ** ******************************************************************
* ** spdas: to get the information for Solar angle section: solar

* **  think about finalcuttau ? correct or not? also invoke of finalcuts(4)?
* ** finalcuts(6)=Njet>2 and Mjj-MW<20
* ** tauevt1p(jj)/tauevt1ph(jj)/tauevt3p(jj)

      if(finalcuts(1).and.finalcuts(2).and.finalcutmu
     $  .and.finalcuttau.and.finalcuts(5).and.finalcuts(7) ) then 
            nev_mutau=nev_mutau+1 
      endif

      if(finalcuts(1).and.finalcuts(2).and.finalcutmu
     $  .and.finalcuts(6) ) then 
            nev_muqq=nev_muqq+1 
      endif

      if(finalcuts(1).and.finalcuts(2).and.finalcuttau
     $  .and.finalcuts(5).and.finalcuts(6).and.finalcuts(7) ) then 
            nev_tauqq=nev_tauqq+1 
      endif

      if(finalcuts(1).and.finalcuts(2).and.finalcute
     $  .and.finalcuttau.and.finalcuts(5).and.finalcuts(7) ) then 
            nev_etau=nev_etau+1 
      endif

      if(finalcuts(1).and.finalcuts(2).and.finalcute
     $  .and.finalcuts(6) ) then 
            nev_eqq=nev_eqq+1 
      endif

            if(.not.usecut(0))evento=.False.

 213     CONTINUE

* ** ******************************************************************
* ** spdas: Right after exiting this 213 loop (Neutralino event characteresitics) 
* ** all the numbers from decay parameters from all the Number of Neutralino present
* ** in the events.

********************check cuts***********************
* ** spdas: finalcutem and finalcute are same; however 
* ** construcnting other cuts using those might be different.

         if(evento)nevtcut=nevtcut+1
         do jj=1,iichi
            finalcutem=int(pfulle(jj,7)).eq.1  
     $           .and.pfulle(jj,5).gt.20.
     $           .and.abs(pfulle(jj,6)).lt.trackmaxeta
            finalcutmum=int(pfullmu(jj,7)).eq.1
     $           .and.pfullmu(jj,5).gt.20.
     $           .and.abs(pfullmu(jj,6)).lt.trackmaxeta
            finalcuttaum=int(pfulltau(jj,7)).eq.1
     $           .and.pfulltau(jj,5).gt.20.
     $           .and.abs(pfulltau(jj,6)).lt.trackmaxeta
            partialcuts=ldecjj(jj) ! finalcuts(1)
     $           .and.trigger   !finalcuts(2)
     $           .and.((finalcutem).or.(finalcutmum).or.(finalcuttaum)) !finalcuts(3)
     $           .and.dabs(pfulljet1(jj,6)).lt.trackmaxeta
     $           .and.dabs(pfulljet2(jj,6)).lt.trackmaxeta !finalcuts(4)
     $           .and.etltau(jj).lt.ET_isola !finalcuts(5) 

            if(morefilter.eq.0)then
               fullwrite=.True.
            elseif(morefilter.eq.1)then
               fullwrite=enocut(jj).or.munocut(jj).or.taunocut(jj)
            elseif(morefilter.eq.7)then
               fullwrite=enocut(jj)
            elseif(morefilter.eq.2)then
               fullwrite=munocut(jj)
            elseif(morefilter.eq.3)then
               fullwrite=taunocut(jj)
            elseif(morefilter.eq.4)then
               fullwrite=usecut(nusecut)
            elseif(morefilter.eq.5)then
               fullwrite=partialcuts.and.njetwjj(jj).eq.2.and.
     $              massinvw(jj).gt.60..and.massinvw(jj).lt.100.
     $              .and.taunocut(jj).and.sphere_line(jj)
     $              .and.nprong(jj).le.3
            elseif(morefilter.eq.6)then
               fullwrite=partialcuts.and.njetwjj(jj).eq.2
            elseif(morefilter.eq.-1)then
               fullwrite=.False.
            end if
            nsphereline=0
            if(sphere_line(jj))nsphereline=1
            if(fullwrite)write(92,*)massinvw(jj),njet,njetwjj(jj),
     $           drmintot(jj),
     $           (pfulle(jj,jjj),jjj=1,8), ! DR's version at the end introduced. spdas
     $           (pfullmu(jj,jjj),jjj=1,8),(pfulltau(jj,jjj),jjj=1,8),
     $           (pfullchi(jj,jjj),jjj=1,8),
     $           (pfulljet1(jj,jjj),jjj=1,8),
     $           (pfulljet2(jj,jjj),jjj=1,8),
     $           (pfullp1(jj,jjj),jjj=1,8),(pfullp2(jj,jjj),jjj=1,8),
     $           (pfullprng1(jj,jjj),jjj=1,8),
     $           (pfullprng2(jj,jjj),jjj=1,8),
     $           (pfullprng3(jj,jjj),jjj=1,8),
     $           (vchi(jj,jjj),jjj=1,3),(vtau(jj,jjj),jjj=1,3),
     $           ntrigger,(az(jj,jjj),jjj=1,3),
     $           nsphereline,etltau(jj),
     $           xLSPmass,xWmass,xLp,xLp_new,xabsP_LSP_reco,
     $           x_bLSP,decayLtmp,x_thetaLSP,x_yLSP,x_etaLSP
*** spdas: other variables that can be useful: decayLtmp,rminvtotjnew,rminvtotjw
         end do

 200  continue ! END OF LOOP over events

      trilep=nevtcut
      inev=nevtcut
* ** spdas: at the moment we are not working on "porodskandsl" parton level cuts.
        write(*,*)'nevgenwe,nevgenwmu,nevgenwtau3p:  (95,95 and 85%)'
        write(*,*)'b4 eff.inclusion:',nevgenwe,nevgenwmu,nevgenwtau3p
      if(porodskandsl)then
         nevgenwe=int(real(nevgenwe)*0.95)
         nevgenwmu=int(real(nevgenwmu)*0.95)
         nevgenwtau3p=int(real(nevgenwtau3p)*0.85)
      end if
      write(*,*)'after eff.inclusion:',nevgenwe,nevgenwmu,nevgenwtau3p

      write(*,*)'solar - number in the main program:'
      write(*,*)'Evt: mutau,muqq,tauqq',nev_mutau,nev_muqq,nev_tauqq
      write(*,*)'Evt: etau,eqq,tauqq',nev_etau,nev_eqq,nev_tauqq

      call sigandnev

c...There are two neutralinos 
      write(*,*)'Final: decaylength and Mass::',decayL/(2.d0*ievmax),
     $     PMAS(PYCOMP(1000022),1)
c      write(*,*)real(inev),real(ievmax),real(inev)/real(ievmax)
c      call pystat(1)
      write(*,132)nevgenwe,nevgenwmu,nevgenwtau
      write(*,133)nevgenwtau1p,nevgenwtau1ph,nevgenwtau3p,nevgenwtau4p

      write(*,139)float(mucut(6))/float(taucut(6)) 
      write(*,141)dsqrt(1.d0/float(mucut(6))+1.d0/float(taucut(6)))
C.DEBUG
      write(*,*)'------------------------------------------------------'
      write(*,134)nenocut,nmunocut,ntaunocut
      write(*,*)'tan^2(theta_atm) from the above mu/tau',
     +         float(nmunocut)/float(ntaunocut)
      write(*,*)'------------------------------------------------------'
      write(*,*)'cut used:',nusecut
      if(nusecut.eq.8) write(*,*)'Full CUTS implementation'
      write(*,*)'------------------------------------------------------'
      write(*,135)(nfinalcuts(i),i=1,numberofcuts)
      write(*,136)(ecut(i),i=1,6)
      write(*,137)(mucut(i),i=1,6)
      write(*,138)(taucut(i),i=1,7)
      write(*,151)(float(mucut(jj))/float(taucut(jj)),jj=1,6)
      write(*,*)'------------------------------------------------------'
      write(*,*)'0 (smear.f) or 2 (Pythia in-built) mstu(53)=>',mstu(53)
      if(porodskandsl)write(*,*)'Porod-Peter-Parton level'


132   format("nwe,nwmu,wtau:",3(1x,i12))
133   format("nwtau 1p, 1ph,nwtau 3p, 4p+:",4(1x,i12))
134   format("TOTAL:nwe,nwmu,nwtau :",3(1x,i12))
135   format("partialcuts          :",7(1x,i12))
136   format("progressivecuts   (e):",6(1x,i12))
137   format("progressivecuts  (mu):",6(1x,i12))
138   format("progressivecuts (tau):",7(1x,i12))
139   format("tan^2_theta{Atm}:",f9.3)
141   format("Delta(R)\over R :",f9.3)
151   format("mu/tau ratios for 1-6 cuts:",6(1x,f9.3)) 

      close(10) ! this is the closing of the SPhenoRP.spc file

      end
****************
c......FLAG if a line intersect a sphere: if less than zero
C......no intersection
      
      double precision function SLintersection(v1,v2,p2,r)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      double precision v1(3),v2(3),p2(3)
      double precision a,b,c,r
      integer i
C.....See mathematica/line-sphere*.nb for details 
* http://www.csee.umbc.edu/~olano/435f02/ray-sphere.html
* http://www.devmaster.net/wiki/Ray-sphere_intersection
******************
      a=0d0;b=0d0;c=0d0
      do i=1,3
         a=a+p2(i)*p2(i)
         b=b+2d0*p2(i)*(v2(i)-v1(i))
         c=c+(v2(i)-v1(i))*(v2(i)-v1(i))
      enddo
      c=c-r**2
      SLintersection=b**2-4d0*a*c
      end      
c********************************
c......Calculates the invariant mass of npart four-momentums
C...  INPUT: 
C...  npart: number of momentums
C...  ppart(1,4),...,ppart(npart,4): Each one of the four-momentums
C...  WARNING "* >= npart" in "double precision ppart(4,*)..." below
      
      double precision function DInvariantMass(npart,ppart)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      double precision ppart(4,*),gmn(4),pminv(4)
******************
      aux=0d0
      gmn(1)=-1d0;gmn(2)=gmn(1);gmn(3)=gmn(1);gmn(4)=1d0
      do i=1,4
         pminv(i)=0
c.... loop over npart: number of tracks of the displaced vertex
         do j=1,npart
            pminv(i)=pminv(i)+ppart(i,j)
         end do
         aux=aux+gmn(i)*pminv(i)**2
      end do
      if(aux.gt.0)then
         DInvariantMass=dsqrt(aux)
      else
         DInvariantMass=-1d0
      end if
      end
*****************************************      
*****************************************      
      subroutine invariantmass(ndaug,kfdaug,nldaug,v1,allfst,
     $     rminvpar,rminvtot)
C...  Calculate (1) the invariant mass of neutralino from the full set of their
C...                stable daughters for 
C...                  (A)  ~chi_10 -> nu_i b barb
C...                  (B)  ~chi_10 -> (mu or tau) W
C...                The nu_i in (A) and mu in (B) belong to the set of stable 
C...                 daughters.
C...            (2) The invariant mass of the heavier unstable (HU), parton level, 
C...                neutralino daughter, from the full set of HU stable daughters.
C...                HU is "b" for (A) or "W" for (B)
C...                
C...   INPUT:
C...    ndaug: number of neutralino daughters
C...    kfdaug(*): KF code of neutralino daughters
C...    nldaug(*): line number in table of neutralino daughters
C...    allfst: wheter consider only charged tracks (.false.) or all
C...   OUTPUT: rminvpar: invariant mass of main daughter
C...           rminvtot: invariant mass of the event
      implicit double precision(a-h,o-z)
      integer kfdaug(*),nldaug(*)
      double precision v1(*) 
      logical allfst
C...Internal variables
      double precision pminvtot(4,10),pminvpar(4,1)
      integer kfstdaug(4000),lnstdaug(4000)
      double precision pstdaugm(5,1000)
      character*16 namestdaug(1000)
C... PYTHIA COMMONS
C.....The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C.....Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
******************PREAMBLE**************************

      ndaugnew=ndaug
C.....Whether of to analyze:
C======== (A) ~chi_10 -> nu_i b barb =============
      if(ndaug.eq.3)then
C          convert the process to a 2 body one ~chi_10 -> nu_i b
C          The string from b will include barb in pythia:
         ndaugnew=2
         inonstable=5 ! flag main unstable particle at parton level
         do ii=1,ndaug
            if(abs(kfdaug(ii)).gt.6)then
               istable=kfdaug(ii) ! flag neutrino as stable
               if(ii.eq.3)then ! never happens!
c... neutrino information in ndaug=3, pass it to 2 (
                  kfdaug(2)=kfdaug(3)
                  nldaug(2)=nldaug(3)
                  stop
               end if
            end if
         end do
C=================================================
      else
C======== (B) ~chi_10 -> (mu or tau) W ==========
         inonstable=24 ! flag main unstable particle at parton level
         istable=13 ! flag muon as stable
      endif
***************************************************
      do 209 ii=1,ndaugnew    ! loop over neutralino final states
         if(abs(kfdaug(ii)).eq.istable)then ! nu (A) or mu (B) contribution
            do jjjj=1,4
               pminvpar(jjjj,1)=p(nldaug(ii),jjjj)
               pminvtot(jjjj,ii)=p(nldaug(ii),jjjj)
            end do
         else  ! find b and bbar (A) or tau and W (B) daughters
C...              bbar is in the same string as b and no require analysis.
            newline=nldaug(ii)
            call fulldecaychain(newline,allfst,istdaug,
     $           kfstdaug,lnstdaug,namestdaug,pstdaugm)
            nmtorig=istdaug
            do jjjj=1,4
               pminvtot(jjjj,ii)=0d0
            end do
            do jjj=1,nmtorig
               do jjjj=1,4
                  pminvtot(jjjj,ii)=pminvtot(jjjj,ii)+
     $                 pstdaugm(jjjj,jjj)
               end do
            end do
            do  jjjj=1,4
               pminvpar(jjjj,1)=pminvtot(jjjj,ii)
            enddo
            rminvchi=DInvariantMass(1,pminvpar) 
            if(abs(kfdaug(ii)).eq.inonstable)then ! b (A) or W (B)
               rminvpar=rminvchi
            endif
         end if        
         
 209  CONTINUE
      rminvtot=DInvariantMass(ndaugnew,pminvtot)
      END
*****************************************      
      subroutine invariantmassjet(ndaug,kfdaug,nldaug,v1,allfst,njet,
     $     rminvtot)
C...   Which of the jets comes from 
C...   INPUT:
C...    ndaug: number of neutralino daughters
C...    kfdaug(*): KF code of neutralino daughters
C...    nldaug(*): line number in table of neutralino daughters
C...    allfst: wheter consider only charged tracks (.false.) or all
C...    njet: number of jets in the event
C...   OUTPUT: rminvpar: invariant mass of main daughter
C...           rminvtot: invariant mass of the event
      implicit double precision(a-h,o-z)
      integer kfdaug(*),nldaug(*)
      double precision v1(*) 
      logical allfst
C...Internal variables
      integer nlwdaug(10),nminj(2),nminjj(2)
      double precision pw(4,1),pq(4,10),pj(4,100),pjj(4,100,100) 
      double precision ptmp(4,100),thmin(2)
      logical show
C... PYTHIA COMMONS
C.....The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
******************PREAMBLE**************************
      show=.False.
      inonstable=24             ! flag main unstable particle at parton level
C...  Here I need to know what are the quarks which are associated to W
***************************************************
C...   obtain W daughters form informative [k(i,3)=21] lines!

      do kkk=1,100
         do kk=1,4
             pj(kk,kkk)=0.d0 ! just to make sure for initilizations added by spdas
         enddo
      enddo

      do 209 ii=1,ndaug         ! loop over neutralino final states
         if(abs(kfdaug(ii)).eq.inonstable)then ! nu (A) or mu (B) contribution
            nwline=nldaug(ii)
            if(k(nldaug(ii),1).eq.21)then 
               nwdaug=0
               do i=nldaug(ii)+1,n
                  if(k(i,3).eq.nldaug(ii))then 
                     if(k(i,2).ne.k(nldaug(ii),2))then
                        nwdaug=nwdaug+1
                        nlwdaug(nwdaug)=i
                     else
                        goto 210 ! exit loop
                     end if
                  end if
               end do
 210           CONTINUE
            end if
         end if
 209  CONTINUE
C...  Compares jet momentum with parent and partons momentum
C...  calculates both the differences and the angle between the vectors
C.....  define each mometum
      do i=1,4
         pw(i,1)=p(nwline,i)
         do ii=1,nwdaug
            pq(i,ii)=p(nlwdaug(ii),i)
         end do
         do ii=1,njet
            pj(i,ii)=p(n+ii,i)
         end do
C...  Calculate the momentum of each pair of jets
C...  There are (njet^2-njet)/2 possibilities:
         njj=0
         do ii=1,njet
            do jj=1,njet
               if(jj.gt.ii)then
                  njj=njj+1
                  pjj(i,ii,jj)=pj(i,ii)+pj(i,jj)
               end if
            end do
         end do
      end do
      pi=dacos(-1d0)
      do i=1,2
         thmin(i)=pi
         nminj(i)=0
         nminjj(i)=0
      end do
      thminj=pi
      do 211 i=1,njet
c         do j=1,nwdaug
c            write(*,*)'Pjet',i,'q',j,' : ',vecA_vecB(pj,pq,i,j)
c         end do
         do j=1,nwdaug ! check each quark against each jet
            thetajq=theta(pj,pq,i,j)
            if(thetajq.gt.pi)thetajq=2d0*pi-thetajq
            if(thetajq.lt.thmin(1))then 
               thmin(2)=thmin(1)
               nminj(2)=nminj(1)
               thmin(1)=thetajq
               nminj(1)=i
            end if
            if(thetajq.lt.thmin(2).and.thetajq.gt.thmin(1))then
               thmin(2)=thetajq
               nminj(2)=i
            end if 
c            write(*,*)'Tjet',i,'q',j,' : ',thetajq,thmin(1),thmin(2),
c     $           nminj(1),nminj(2)
         end do
         do j=1,njet
            if(j.gt.i)then
               do ii=1,3
                  ptmp(ii,1)=pjj(ii,i,j)
               end do
               thetajq=theta(pw,ptmp,1,1)
c               write(*,*)'WJtjet',i,',',j,' : ',thetajq,thminj
               if(thetajq.lt.thminj)then
                  thminj=thetajq
                  nminjj(1)=i;nminjj(2)=j
               end if
            end if
         end do
 211  CONTINUE
      ntmp=nminj(1)
      if(nminj(1).gt.nminj(2))then
         nminj(1)=nminj(2)
         nminj(2)=ntmp
      end if
      if(show)then
         write(*,*)'min',thmin(1),thmin(2),thminj
         write(*,*)'minjet',nminj(1),nminj(2)
         write(*,*)'minjetjj',nminjj(1),nminjj(2)
      end if
      if(nminj(1).eq.nminjj(1).and.nminj(2).eq.nminjj(2))then
         if(show)then
            write(*,*)'W -> jet(',nminj(1),') jet(',nminj(2),')'
            write(*,*)'invariant mass from quarks:',
     $           DInvariantMass(nwdaug,pq)
            write(*,*)'invariant mass from jets:'
         end if
         do j=1,4
            ptmp(j,1)=pjj(j,nminj(1),nminj(2))
         end do
         rminvtot=DInvariantMass(1,ptmp)
         if(show)write(*,*)rminvtot
      else
         if(show)write(*,*)'failed W mass reconstruction from jets!'
         rminvtot=-1d0
      end if
C...  Assume tha only one jet was created:
      if(rminvtot.lt.0.and.nminj(1).eq.nminj(2))then
         do j=1,4
            ptmp(j,1)=pj(j,nminj(1))
         end do
c         write(*,*)'ptmp',(ptmp(j,nminj(1)),j=1,4)
         rminvtot=DInvariantMass(1,ptmp)
c         write(*,*)'invariant mass from 1 jets:',rminvtot
      end if
      END
**********************************************************
      subroutine invariantmassjetW(ndaug,kfdaug,nldaug,v1,allfst,njet,
     $     ppartonw,njetw,DRmin,rminvtot,rminvtotw,ptmpj)
C...   Which of the jets comes from 
C...   INPUT:
C...    ndaug: number of neutralino daughters
C...    kfdaug(*): KF code of neutralino daughters
C...    nldaug(*): line number in table of neutralino daughters
C...    allfst: wheter consider only charged tracks (.false.) or all
C...    njet: number of jets in the event
C...   OUTPUT: 
C...    DRmin(1)...DRmin(njetw): Minimum distance between each W jet with
C...                             the others jet: no longer neccesary
C...    rminvtotj: invariant W mass calculated from the full jet system
      implicit double precision(a-h,o-z)
      integer kfdaug(*),nldaug(*)
      double precision v1(*) 
      logical allfst
      double precision DRmin(*),DRminw(100),DRminwpart(100)
      integer ndrminw(100),ndrminwpart(100)
C...Internal variables
      double precision pminvtot(4,10),pminvpar(4,1)
      integer kfstdaug(4000),lnstdaug(4000)
      double precision pstdaugm(5,1000)
      character*16 namestdaug(1000)
      integer nlwdaug(10),nminj(10),nminjj(2)
      double precision pw(4,1),pq(4,10),pj(4,100),pjj(4,100,100)
      double precision pjother(3,100)
      double precision ptmp(4,100),thmin(10),pnfj(4,100),pjw(4,100)
      double precision pnfjw(4,100),pnfj2(4,100)
      logical wdaug,show
      integer kold(4000),knp(100),ndrmin(20)
      double precision drj1(20),drj2(20),drj3(20),drwjii(20),drwjij(20)
      logical jetw_eq_jet
      double precision ppartonw(4,10),ptmpj(8,100)
C...  Input Commons
c      common/wpartonlevel/ppartonw(4,10)
C...  Output Commons
c      common/jetprop/ptmpj(8,100)
C... PYTHIA COMMONS
C.....The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C.....Three Pythia functions return integers, so need declaring.
      INTEGER PYK,PYCHGE,PYCOMP
C.....Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
C.....Parameters.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      common/eventnumber/iev,nchilinejj
***   common/teste/ievv  ! spdas: it was not consistent with other subroutines...
      common/teste/ievv,ptcut,rcut     !MBM
******************PREAMBLE**************************
      inonstable=24             ! flag main unstable particle at parton level
C...  Here I need to know what are the quarks which are associated to W
***************************************************
c      if(iev.eq.6167)then
c         show=.True.
c      else
         show=.False.
c      end if
      maxwjet=4
C...   obtain momentum of each jet

      njetfull=njet
      do ii=1,njetfull
         do i=1,4
            pj(i,ii)=p(n+ii,i)
         end do
         pjother(1,ii)=pyp(n+ii,10)
         pjother(2,ii)=pyp(n+ii,19)
         pjother(3,ii)=pyp(n+ii,15)
      end do
      do ii=1,2
         do i=1,8
            ptmpj(i,ii)=-1e-30
         end do
         ptmpj(7,ii)=0
      end do
C..   obtain eta and phi for each jet
      do i=1,2
         do ii=1,njetfull
            pnfj(i,ii)=pyp(n+ii,11+i*4) !phi and eta
         end do
      end do
      do i=1,njetfull
         knp(i)=k(n+i,4) ! number of particles assigned to jet i
      end do
C...   obtain W daughters 
      do 209 ii=1,ndaug         ! loop over neutralino final states
         if(abs(kfdaug(ii)).eq.inonstable)then ! nu (A) or mu (B) contribution
            newline=nldaug(ii)
            call fulldecaychain(newline,allfst,istdaug,
     $           kfstdaug,lnstdaug,namestdaug,pstdaugm)
C...  DEBUG:
         end if
 209  CONTINUE

C...  Preparing calculation of jets by using only the W daughters
      do i=1,n ! flag no W stable daughters as unestable particles
         kold(i)=-1
         if(k(i,1).lt.11)then
            wdaug=.False.
            do j=1,istdaug
               if(i.eq.lnstdaug(j))wdaug=.True.
            end do
            if(.not.wdaug)then
               kold(i)=k(i,1)
               k(i,1)=11

            end if
         end if
      end do
c      if(show)call pylist(2)  
      call pycell(njetw) ! uses only W stable daughters

      show=.true. ! to see the effects until next show
      show=.false. ! not write the onwards info
      if(show)then
         write(*,*)'CELLjet *************'
         call pylist(2)
         write(*,*)'njetsfull,njetw',njetfull,njetw
         write(*,*)'phi-',(pyp(n+kk,15),kk=1,njetw)
         write(*,*)'eta-',(pyp(n+kk,19),kk=1,njetw)
         write(*,*)'lines',(lnstdaug(i),i=1,istdaug)
      endif

      show=.false. ! not write the onwards info

C...  store momentun of jets coming from W:
      do i=1,4
         do ii=1,njetw
            pjw(i,ii)=p(n+ii,i)
         end do
      end do
      do i=1,2
         do ii=1,njetw
            pnfjw(i,ii)=pyp(n+ii,11+i*4) !phi and eta
         end do
      end do
      if(show)write(*,*)'phi,eta::',pnfjw(1,1),pnfjw(2,1)
C...  Pick out the originals jets
C.... (A) by using the angle betwen the jets
      pi=dacos(-1d0)
c      do i=1,njetw
c
c      end do
      do i=1,njetw
         thmin(i)=pi
         if(show)write(*,*)'** angle between jw and j *******'
         do j=1,njetfull
            thetajq=theta(pjw,pj,i,j)
            if(thetajq.gt.pi)thetajq=2d0*pi-thetajq
            if(thetajq.lt.thmin(i))then 
               thmin(i)=thetajq
               nminj(i)=j
            end if
            if(show)write(*,*)'jW:',
     $           i,' j:',j,' : ',thetajq,thmin(i),nminj(i)
         end do
      end do
      thmax=0d0
      do i=1,njetw
         if(thmax.lt.thmin(i))thmax=thmin(i)
      end do
      if(show)then
         do i=1,njetw 
            write(*,*)'CELLjet: jetw', i,' is: jet',nminj(i),thmin(i)
         end do
      end if
C...  (B) By using DeltaR Calculation
C...  (B.1) with W jets

      do i=1,njetw
         DRminw(i)=1000D0
      end do

      do j=1,njetfull 
         do i=1,njetw
            deltarij=DeltaR(pnfjw,pnfj,i,j)
            if(deltarij.lt.DRminw(i))then
               DRminw(i)=deltarij
               ndrminw(i)=j

            end if
         end do        
      end do
C...  (B.2) with W partons

      do i=1,2
         DRminwpart(i)=1000D0
      end do

      deltarparton=dsqrt((ppartonw(1,1)-ppartonw(1,2))**2+
     $     (ppartonw(2,1)-ppartonw(2,2))**2)
      do j=1,njetfull 
         do i=1,2
            deltarpartij=DeltaR(ppartonw,pnfj,i,j)

            if(deltarpartij.lt.DRminwpart(i))then
               DRminwpart(i)=deltarpartij
               ndrminwpart(i)=j
            end if
            if(show)print*,'parton:',deltarpartij,' parton',i,' jet',j,
     $           DRminwpart(i),ndrminwpart(i)
         end do        
      end do
      if(show)then
         do i=1,2
            write(*,*)'parton: parton', i,' is: jet',ndrminwpart(i)
     $           ,DRminwpart(i)
         end do
         write(*,*)'DRparton',deltarparton,njetw,ppartonw(1,1),
     $        ppartonw(2,1)
      end if


C... Distance of no W jets to W jets
      do i=1,njetw
         DRmin(i)=1000D0
      end do
      do i=1,5
         if(i.eq.njetw)then
            if(njetw.lt.5)then
               do j=njetw+1,4
                  nminj(j)=-1
               end do
            else
               write(*,*)'WARNING: invariantmassjetw: too many jets:',
     $              njetw
c               stop
            end if
         end if
      end do
C... Analysis only of there at most maxwjet=4 wjets
      do j=1,njetfull 
         do i=1,njetw
            if(j.ne.nminj(1).and.j.ne.nminj(2).and.j.ne.nminj(3).and.
     $           j.ne.nminj(maxwjet))then
               deltarij=DeltaR(pnfjw,pnfj,i,j)
               if(deltarij.lt.DRmin(i))then
                  DRmin(i)=deltarij
                  ndrmin(i)=j
               end if
            end if
         end do        
      end do
      if(show)then
         write(*,*)'DRmin i=1,njetw',(DRmin(i),i=1,njetw)
         write(*,*)'ndrmin i=1,njetw',(ndrmin(i),i=1,njetw)
      end if


C...Invariant W mass calculation from jets
C...(A) By using W jets:
      rminvtotw=DInvariantMass(njetw,pjw)
      if(show)write(*,*),'inv mass from W jets',rminvtotw
C...(B) By using jets identified with W jets:
c...  W jets identified when the identification with theta is equal to  idendification with DeltaR
      jetw_eq_jet=.True. 
      do i=1,njetw
         if(nminj(i).ne.ndrminw(i))jetw_eq_jet=.False. !  this is the isolation-matching of full-jet AND W-jet
         if(njetw.gt.maxwjet)jetw_eq_jet=.False.  

* ** spdas:  as maxwjet=4 thus we are taking events only containing less than 4-jets; 
* ** however in practice I think this is somehow brute force; it should 
* ** be better to re-ordering  the above clause; as the isolation is much important 
* ** than the number of jets. 

      end do
      if(jetw_eq_jet)then 
         do i=1,njetw
            do j=1,4
               ptmp(j,i)=pj(j,nminj(i))
               ptmpj(j,i)=pj(j,nminj(i))
            end do

* ** spdas: Since we are assigning jets originating from W decays 
* ** kinematics  here, as "pjother" contains the info for the full 
* ** jet-system thus the above "if" is picking the correct jets from W.

            ptmpj(5,i)=pjother(1,nminj(i))
            ptmpj(6,i)=pjother(2,nminj(i))
            ptmpj(8,i)=pjother(3,nminj(i))
            ptmpj(7,i)=1
         end do
         rminvtot=DInvariantMass(njetw,ptmp)
      else
         rminvtot=-1D0 
      end if
      if(show)write(*,*),'inv mass from full jets',rminvtot

C...(C) with jets identfied from partons
      do i=1,2
         do j=1,4
            ptmp(j,i)=pj(j,ndrminwpart(i))
         end do
      end do
      rminvtotpar=DInvariantMass(2,ptmp)
      if(show)write(*,*),'inv mass from jets(partons)',rminvtotpar


      if(show)then
         write(*,*)'DRminw: i=1...njetw=',njetw,':', 
     $        (DRminw(i),i=1,njetw)
         do i=1,njetw
            write(*,*)'i=',i,' # particl j/ # particl jw',
     $           real(k(n+i,4))/real(knp(nminj(i))),
     $           ' ndremin',ndrminw(i),' nminj',nminj(i)
         end do
      end if

      call jetisola(show,nminj,njetw,njetfull,pnfjw,pnfj)
      deltarjj=-1000d0
      if(njetw.eq.2)then
         do j=1,2
            pnfj2(j,nminj(2))=pnfj(j,nminj(2))
         end do
         deltarjj=DeltaR(pnfj,pnfj2,nminj(1),nminj(2))
      end if
C... Contamination

c      do i=njetw+1,5
c         DRmin(i)=1000D0
c      end do
c      if(njetw.gt.0)write(93,*)njetfull,njetw,
c     $     (drj1(i),i=1,15),(drj2(i),i=1,15),
c     $     (drj3(i),i=1,15),(drwjii(i),i=1,3),(drwjij(i),i=1,3),
c     $     (thmin(i),i=1,3),deltarmin,thmax,(DRmin(i),i=1,3),
c     $     real(k(n+1,4))/real(knp(nminj(1))),
c     $     DInvariantMass(njetw,pjw)/rminvtot,rminvtot

C.... recover original table
      do i=1,n
         if(kold(i).ne.-1)k(i,1)=kold(i)
      end do

      call pycell(njet)

      END
**********************************************************
      subroutine isola(ndaug,kfdaug,nldaug,DRimin,DRimax,ET,pt)
C...   Which of the jets comes from 
C...   INPUT:
C...    ndaug: number of neutralino daughters
C...    kfdaug(*): KF code of neutralino daughters
C...    nldaug(*): line number in table of neutralino daughters
C...   OUTPUT: 
C...     ET: Transverse energy inside the cone for the lepton in
C...          ~chi_10 -> W+ l- (l=e,mu,tau)
C...         In the case of tau is the maximum ET between 
C...         minimum and maximum cone
      implicit double precision(a-h,o-z)
      integer kfdaug(*),nldaug(*)
C...  Internal variables
      double precision pminvtot(4,10),pminvpar(4,1)
      integer kfstdaug(4000),lnstdaug(4000)
      double precision pstdaugm(5,1000)
      character*16 namestdaug(1000)
      integer nlwdaug(10),nminj(2),nminjj(2),nlep(3),nlepb(3)
      double precision pw(4,1),pq(4,10),pj(4,100),pjj(4,100,100) 
      double precision ptmp(4,100),thmin(2),pnfl(4,1),pnfp(4,1)
      double precision DRimin(15),DRimax(15)
      integer nlnew(20)
      logical show,hadron,allfst,tau1phadron
      double precision dricheck(100)
      double precision pmuta(3),abspmuta,abspmutat,tmpp,tmpn ! DR
C... PYTHIA COMMONS
C.....The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C.....Three Pythia functions return integers, so need declaring.
      INTEGER PYK,PYCHGE,PYCOMP
C.....Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
******************PREAMBLE**************************
      show=.False.
      tau1phadron=.False.
C...  Here I need to know what are the quarks which are associated to W
*************************************************
C...  Replace nldaug(i) for no informative lines nlnew(i)
      do i=1,ndaug
         if(k(nldaug(i),1).eq.21)then
            do j=i+1,n
               if(k(j,3).eq.nldaug(i).and.
     $              abs(k(j,2)).eq.abs(kfdaug(i)))then
                  nlnew(i)=j
                  goto 210
               end if
            end do
 210        CONTINUE
         end if
      end do
*************************************************
      do i=1,3
         nlep(i)=9+i*2
c         Drimin(9+i*2)=0d0
c         Drimax(9+i*2)=0.3d0
      end do
c      Drimax(15)=0.4d0
c      Drimin(15)=0.1d0
C...   obtain W daughters from informative [k(i,3)=21] lines!
      do 209 ii=1,ndaug         ! loop over neutralino final states
         if(abs(kfdaug(ii)).eq.nlep(1).or.abs(kfdaug(ii)).eq.nlep(2)
     $        .or.abs(kfdaug(ii)).eq.nlep(3))then 
            if(k(nlnew(ii),1).lt.11)then ! mu or electron: stable particle DR
               nwline=nlnew(ii) ! line in which the particle is
               nkfdaug=kfdaug(ii) ! KF of lepton for DRimin(KF(lepton))
               call isolai(nwline,Drimin,Drimax,nkfdaug,ET,pt)
            else                ! tau
C...  In case of a tau, we need to check for each charged track!
C... For 3 prong the total momentum is summed up  ! DR
               newline=nlnew(ii) ! line in which the lepton is
               allfst=.False. ! only charged particles DR
               call fulldecaychain(newline,allfst,istdaug,
     $              kfstdaug,lnstdaug,namestdaug,pstdaugm)
               allfst=.True.
               do i=1,2
                  pnfl(i,1)=pyp(lnstdaug(1),11+i*4) !phi and eta
               end do
C... begin changes ! DR
               do i=1,3
                  pmuta(i)=0d0
               end do
               do i=1,istdaug
                  do j=1,3
                     pmuta(j)=pyp(lnstdaug(i),j) !phi and eta
                  end do
               end do
               if(istdaug.gt.1)then
                  abspmuta=dsqrt(pmuta(1)**2+pmuta(2)**2+pmuta(3)**2)
                  abspmutat=dsqrt(pmuta(1)**2+pmuta(2)**2)
C..   Azimuthal angle phi is defined between -pi and pi. 
C...  See http://en.wikipedia.org/wiki/Spherical_coordinate_system
                  pnfl(1,1)=datan2(pmuta(2),pmuta(1)) ! phi
                  pnfl(2,1)=(1d0/2d0)*
     $                 dlog((abspmuta+pmuta(1))/(abspmuta-pmuta(1))) ! eta
               end if
C... end changes ! DR
               etmax=0d0
               ptmax=0d0
               do i=1,istdaug
                  nwline=lnstdaug(i) ! line in which the particle is
                  nkfdaug=abs(kfdaug(ii)) ! KF of lepton for DRimin(KF(lepton))

                  call isolai(nwline,Drimin,Drimax,nkfdaug,ET,pt)
                  if(ET.gt.etmax)etmax=ET
                  if(pt.gt.ptmax)ptmax=pt
               end do 
               ET=etmax
               pt=ptmax
               if(istdaug.eq.1)then
                  if(abs(k(lnstdaug(1),2)).gt.15)tau1phadron=.True.
               end if

            end if
         end if
 209  CONTINUE
      END
**********************************************************
************************************
      subroutine isolai(nwline,Drimin,Drimax,nkfdaug,ET,pt)
C...  INPUT:
C...      nwline: line in table of particle for isolation check
C...      Drimin(abs(nkfdaug)). Minimum cone
C...      Drimax(abs(nkfdaug)). Maximum cone
C...      nkfdaug: =abs(k(nwline,2)) for nkfdaug=11,13
C...               =15 for tau. k(nwline,2) in this case are for
C...                       tau stable daughters
C...  OUTPUT
C...      ET: Transver Energy of other particles between 
C...          minimum and maximum cone
      implicit double precision(a-h,o-z)
C...  INPUT/OUTPUT
      double precision DRimin(15),DRimax(15)
C.... Internal variables
      double precision pnfl(4,1),pnfp(4,1)
      logical hadron
C... PYTHIA COMMONS
C.....The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C.....Three Pythia functions return integers, so need declaring.
      INTEGER PYK,PYCHGE,PYCOMP
C.....Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
******************PREAMBLE**************************
      if(nkfdaug.eq.11.or.nkfdaug.eq.13)then
         if(abs(k(nwline,2)).ne.nkfdaug)then
            write(*,*)'ERROR: Bad nkfdaug input in isolai',
     $           nkfdaug,abs(k(nwline,2))
c            stop
         end if
      end if
      do i=1,2
         pnfl(i,1)=pyp(nwline,11+i*4) !phi and eta
      end do
      et=0d0
      ptmax=0d0
c     call pylist(2)
      do i=1,n
         if(k(i,1).lt.11.and.i.ne.nwline)then
            hadron=.True.
            if(abs(k(i,2)).ge.11.and.abs(k(i,2)).le.16)
     $           hadron=.False. ! exclude leptons
            if(abs(k(i,2)).eq.22)hadron=.False. ! exclude photons
            do j=1,2
               pnfp(j,1)=pyp(i,11+j*4) !phi and eta
            end do
            eti=pyp(i,4)*dsin(pyp(i,13))
            pti=pyp(i,10)
            DRi=DeltaR(pnfl,pnfp,1,1)
            if(Dri.gt.Drimin(abs(nkfdaug)).and.
     $           DRi.lt.Drimax(abs(nkfdaug)).and.hadron)then
               et=et+eti
               if(pti.gt.ptmax)ptmax=pti
            end if

         end if
      end do
      pt=ptmax
**********************************************
      end
**********************************************
      subroutine jetisola(show,nminj,njetw,njetfull,pnfjw,pnfj)
      implicit double precision(a-h,o-z)
C...Internal variables
      integer nminj(10)
      double precision pnfjw(4,100),pnfj(4,100)
      double precision pminvtot(4,10),pminvpar(4,1)
      double precision pstdaugm(5,1000)
      logical show
      double precision drj1(20),drj2(20),drj3(20),drwjii(20),drwjij(20)
C... PYTHIA COMMONS
******************PREAMBLE**************************
      inonstable=24             ! flag main unstable particle at parton level
C...  Here I need to know what are the quarks which are associated to W
***************************************************
C...  Fine analysis of jets. check that others no W jets are 0.8 away
C...  The analysis will be done only until 4 jets from W
      if(show)write(*,*)'CELLjet: other jet info:'
      if(show)write(*,*)'CELLjet: nminj(i=1,4) main:',
     $     (nminj(i),i=1,njetw),', others',(nminj(i),i=njetw+1,4)
      do i=1,20
         drj1(i)=-1d0;drj2(i)=-1d0;drj3(i)=-1d0
         drwjii(i)=-1d0;drwjij(i)=-1d0
      end do
      deltarmin=10000d0
      do 210 j=1,njetfull
         if(nminj(1).eq.j.or.nminj(2).eq.j.or.nminj(3).eq.j
     $        .or.nminj(4).eq.j)then
            do i=1,njetw
               if(nminj(i).eq.j)then
                  if(show)write(*,*)'CELLjet: *DRwjetwii*', i,
     $                 ' with jet',j,DeltaR(pnfjw,pnfj,i,j)
                  drwjii(i)=DeltaR(pnfjw,pnfj,i,j)
               else
                  if(show)write(*,*)'CELLjet: DRwjetwij', i,' with jet',
     $                 j,DeltaR(pnfjw,pnfj,i,j)
                  drwjij(i)=DeltaR(pnfjw,pnfj,i,j)
               end if
            end do
         else
            do i=1,njetw
               deltarij=DeltaR(pnfjw,pnfj,i,j)
               if(show)write(*,*)'CELLjet: DR jetw', i,' with jet',j,
     $              deltarij
               if(i.eq.1)drj1(j)=deltarij
               if(i.eq.2)drj2(j)=deltarij
               if(i.eq.3)drj3(j)=deltarij
               if(deltarij.lt.deltarmin)deltarmin=deltarij
            end do
         end if   
 210  CONTINUE

      end
**********************************************************
      double precision function vecA_vecB(A,B,n,m)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      double precision A(4,*),B(4,*),dacos
*****************************
      vecA_vecB=dsqrt( (A(1,n)-B(1,m))**2 +(A(2,n)-B(2,m))**2+
     $     (A(3,n)-B(3,m))**2)
      end
      
**********************************************************
**********************************************************
      double precision function theta(A,B,n,m)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      double precision A(4,*),B(4,*)
*****************************
C...  Rare bug when theta=0d0: dacos(1d0)=NAN !!!
      theta=abs( 
     $     acos(real(
     $            ( A(1,n)*B(1,m)+A(2,n)*B(2,m)+A(3,n)*B(3,m) )/
     $            ( dsqrt( A(1,n)**2+A(2,n)**2+A(3,n)**2 )*
     $              dsqrt( B(1,m)**2+B(2,m)**2+B(3,m)**2 )
     $            )
     $          ))
     $        )
      end
**********************************************************
      double precision function DeltaR(A,B,n,m)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      double precision A(4,*),B(4,*)
*****************************
C...  Rare bug when theta=0d0: dacos(1d0)=NAN !!!
      DeltaR=dsqrt((A(1,n)-B(1,m))**2+(A(2,n)-B(2,m))**2)
      end
******************************************
**********************************************************
      subroutine univec(A,Au)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      double precision A(4),Au(4)
*****************************
      Am=dsqrt(A(1)*A(1)+A(2)*A(2)+A(3)*A(3))
      do i=1,3
         Au(i)=A(i)/Am
      end do
      end
****************************************** 
