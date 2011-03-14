C.....All real arithmetic in double precision.
C.....Three Pythia functions return integers, so need declaring.
      INTEGER PYK,PYCHGE,PYCOMP
C.....Parameter statement to help give large particle numbers
C.....(left- and righthanded SUSY, excited fermions).
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KEXCIT=4000000)
C.....Commonblocks.
C.....The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C.....Parameters.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C.....Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
C.....Decay information.
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
C.....Selection of hard scattering subprocesses.
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
C.....Parameters. 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
C.....Supersymmetry parameters.
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
C.....Information on neutralino, chargino and sfermion
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &     SFMIX(16,4)
C..... Cross section information
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
C..... other stuff
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYDATR/MRPY(6),RRPY(100)
************************---COMMONS---***********************************
      integer trilep,inev,ievmax,isub
      common/threelepc/trilep,inev,ievmax,isub
c      common/daugprop/pfullmu(10,8),pfulltau(10,8),pfulljet1(10,8),
c     $     pfulljet2(10,8),pfullp1(10,8),pfullp2(10,8)
************************************************************************
      integer subpr(500),nsub
      common/mayprocess/subpr,nsub
**************--EXTERNAL VARIABLES--************************************
      real*8 pyp
      external pyp
      

