C.....LHC INITIALIZATION 

      subroutine setlhc

      implicit none
      double precision paru,parj,brat
      integer mstu,mstj,i,mdcy,mdme,kfdp
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      double precision parp,pari
      integer mstp,msti
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      integer npdflib
      common/cpdflib/npdflib
      integer trilep,inev,ievmax,isub
      common/threelepc/trilep,inev,ievmax,isub

      mstu(22)  = 20000         ! max number of errors printed (20000)
      MSTP(111) = 0             ! switches on(1)/off(0) fragmentation
      MSTP(61)  = 1             ! on(1-2)/off(0) inicial state qcd and qed rad.
      MSTP(71)  = 1             ! on(1)/off(0) final state qcd/qed rad.
      MSTJ(41)  = 1             ! only qcd(1), not qed(2)
      MSTP(43)  = 3             ! 1 = gg ; 2 = ZZ; 3 = Z/g
      MSTP(7)   = 6             ! generating tops

C.....Using Rick Field's Tune A for Pythia (see P. Skands comments)
c.....Multiple interactions are strongly recommended for LHC

      MSTP(81)  = 1 !1          ! switches on(1)/off(0) multiple interactions
      MSTP(82)  = 4             ! structure of multiple interactions (D=4)  
      PARP(83)  = 0.5d0         ! fraction of hadronic matter (D=0.5)
      PARP(84)  = 0.4           ! core radius (D=0.2)
      PARP(89)  = 1800.d0       ! reference energy scale (D=1800)
      PARP(90)  = 0.25d0        ! power of energy rescaling (D=0.16)
      PARP(67)  = 4.d0          ! Q2 scale for hard scattering (D=1)
      PARP(82)  = 2.0d0         ! regularization scale P_T cutoff (D=2)
      PARP(85)  = 0.9d0         ! probability of gg mult. interact. (D=0.33)
      PARP(86)  = 0.95d0        ! probability of loop gg mult. inter. (D=0.66)

c..... Z -> tau + tau 

      if(isub.eq.23)then
         do i=174,189
            MDME(i,1) = 0
         enddo
         MDME(186,1) = 1        ! Z -> tau + tau 
      endif

C..... setting pythia to LHC toy calorimeter 

c..... using pyclus
c      mstu(46) = 2        ! (D=1) jet reconstruction choice
c      paru(44) = 24.0d0    ! (D=2.5) djoin
c      paru(45) = 0.05d0   ! (D=0.05d0) see mstu(46)=4
c      mstu(47) = 0        ! (D=1) minimum number of jets
c      mstu(41) = 2        ! see below
c..... using pycell
*************************
      mstu(52) = 64       !(D=24)number of azimuthal bins (deltaphi=2*pi/30)
      paru(54) = 0.4d0  !0.15d0   ! (0.7)Delta R=sqrt(delta eta**2+delta phi**2)
      paru(55) = 0.5d0    ! (0.5)see mstu(53)
*********************
      paru(51) = 5.0d0    !(5.0)maximum absolute pseudorapidity
      mstu(51) = 100      !(200)number of pseudorapidity bins(deltaeta=10/50)
      mstu(54) = 3        ! (D=1) P information in pycell
      mstu(53) = 0        !smearing of the energy paru(55)*sqrt(E) for 2
      paru(52) = 2        ! Minimum energy cluster (D=1.5) ????
      paru(53) = 20.d0    ! (25) (D=7) E_T minimum of the jet. If less hadron is jetless
      mstu(41) = 2        ! (D=2) all partons/particles except neutrinos and..

C.....Set proton distribution

      MSTP(51)=7                ! CTEQ 5L (leading order) (D=7,CTEQ 5L)
      if(npdflib.eq.1)then
         mstp(52) = 2           !user in charge
         mstp(51) = 1000 * 4 + 32 !CTEQ4L, 
c         mstp(51) = 1000 * 4 + 17 !CTEQ2L
      endif

      return
      end














