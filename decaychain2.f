
      subroutine decaychain(orig,iiorig,kfdecay,iidecay,kdecay,fin,
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
      Implicitdouble precision(a-h, o-z)
      integer orig(5),kdecay(5),fin(5),kffin(5,*)
      character*16 finname(5,*),orpar
      logical decay
      double precision pfin(5,5,*)
**************--EXTERNAL VARIABLES--************************************
      real*8 pyp
      external pyp
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer nffin(5,100)
      common/dauglines/nffin
***********************************************************************
      do ii=1,5
         kdecay(ii)=-1
         fin(ii)=0
      enddo
      iidecay=0
      decay=.false.
      Do i=1,n

         if(k(i,1).eq.21)then ! informative lines
c....search if the particle comes from a previous one
            do iii=1,iiorig 
               if(orig(iii).eq.k(i,3).or.orig(iii).eq.-1)then
                  if(abs(k(i,2)).eq.kfdecay)decay=.true.
               endif
            enddo
C....search the decay products of this particle
            if(decay)then ! store the line number particle
               iidecay=iidecay+1
               kdecay(iidecay)=i
               decay=.false.
            else  ! search for decay products of it
               do iiidecay=1,iidecay 
                  if(k(i,3).eq.kdecay(iiidecay))then
                     fin(iiidecay)=fin(iiidecay)+1   
                     call pyname(k(i,2),
     $                    finname(iiidecay,fin(iiidecay)))
                     kffin(iiidecay,fin(iiidecay))=k(i,2)
                     nffin(iiidecay,fin(iiidecay))=i
                     do jj=1,5
                        pfin(jj,iiidecay,fin(iiidecay))=p(i,jj)
                     enddo
                  endif
               enddo
            endif
         else
            goto 207
         endif
      enddo
 207  continue
      if(iidecay.gt.0.and.fin(iidecay).gt.0)then
         call pyname(kfdecay,orpar)
         do i=1,iidecay
c            write(*,*)orpar,'->',(finname(i,jj),jj=1,fin(i))
         enddo
      endif
      END
***********************************************************
      subroutine decaychainnew(orig,iiorig,kfdecay,iidecay,kdecay,fin,
     $     finname,kffin,pfin)
********************************************************************
*  Search for the decays of the particle with KF=kfdecay. This
* unstable particle is the daughter of any of iiorig-th particle occupying 
* the position i=orig(ii) (i=1,n and ii=1,iiorig) in the event history 
* part (k(i,1)=21) of the pythia table. If orig(1)=-1, the parent of the
* unstable particle is no taken into account 
*  INPUT
*  orig(1) ... orig(iiorig) parents of instable particle wich have:
*  kfdecay: as KF code
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
      character*16 finname(5,*),orpar
      logical decay
      double precision pfin(5,5,*)
      integer kfstdaug(4000),lnstdaug(4000)
      double precision pstdaug(5,1000)
      logical allfst
      character*16 namestdaug(1000)

**************--EXTERNAL VARIABLES--************************************
      real*8 pyp
      external pyp
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      integer nffin(5,100)
      common/dauglines/nffin
      do ii=1,5
         kdecay(ii)=-1
         fin(ii)=0
      enddo
      iidecay=0
      decay=.false.
      Do i=1,n

         if(k(i,1).eq.21)then ! informative lines
c....search if the particle comes from a previous one
            do iii=1,iiorig 
               if(orig(iii).eq.k(i,3).or.orig(iii).eq.-1)then
                  if(abs(k(i,2)).eq.kfdecay)then
                     iidecay=iidecay+1
                     kdecay(iidecay)=i
C... Subroutine fulldecaychain
C... INPUT:
C... kdecay(iidecay): line number in table of kfdecay coming from
C... orig(iii)
C..  *********************************************
C... allfst: .False. Only charged final stable states (FSS)
C...         .True.  All FSS
                     allfst=.True. 
C,,, *********************************************
C... OUTPUT:
C...  mtorig(1); number of final stable states (FSS) of tau
C...  kfstdaug(1)...kfstdaug(istdaug): kf of FSS 
C...  lnstdaug(1)...lnstdaug(istdaug): line number in table of FSS
C...  namestdaug(1)...namestdaug(istdaug): names of FSS 
C...  pstdaug(1...5,1)...pstdaug(1...5,istdaug): full momentum of FSS
                     newline=kdecay(iidecay)
                     call fulldecaychain(newline,allfst,
     $                    istdaug,kfstdaug,lnstdaug,namestdaug,
     $                    pstdaug)
                     do ikik=1,istdaug
                        kffin(iidecay,ikik)=kfstdaug(ikik)
                        nffin(iidecay,ikik)=lnstdaug(ikik)
                        finname(iidecay,ikik)=namestdaug(ikik)
                        do ikiki=1,5 
                           pfin(ikiki,iidecay,ikik)=
     $                          pstdaug(ikiki,ikik)
                        end do
                     end do
                     fin(iidecay)=istdaug
                  end if
               end if
            end do
         else
            goto 207
         endif
      enddo
 207  continue
      if(iidecay.gt.0.and.fin(iidecay).gt.0)then
         call pyname(kfdecay,orpar)
         do i=1,iidecay
c            write(*,*)orpar,'->',(finname(i,jj),jj=1,fin(i))
         enddo
      endif
c      write(*,*)'kdecay',iiorig,iidecay,kdecay(iidecay)

      END
***************************************************
C...  Full Decay Chain 
C...  Find all the daughters of the initial particle at line iline
C...    of Pythia table.
      subroutine fulldecaychain(iline,allfst,istdaug,kfstdaug,
     $     lnstdaug,namestdaug,pstdaug)
*************************************************************
C...  INPUT
C...  iline: line in the table of the initial unstable particle
C...  allfst: .False. only stable charged final states. 
C...           .True. All final stable states
C...  OUTPUT: 
C...  istdaug: Number of final stable states (FSS)
C...  kfstdaug(1)...kfstdaug(istdaug): kf of FSS 
C...  lnstdaug(1)...lnstdaug(istdaug): line number in table of FSS
C...  namestdaug(1)...namestdaug(istdaug): names of FSS 
C...  pstdaug(1...5,1)...pstdaug(1...5,istdaug): full momentum of FSS
C...  IT IS RECOMMENDED TO DEFINE THE VECTOR WITH A DIMENSION 
C.... SUFFICIENTLY HIGH, FOR EXAMPLE:
C.... kfstdaug(1000) for an initial quark
*****************************************************************
      implicit double precision(a-h,o-z)
      integer kfstdaug(*),lnstdaug(*)
c.... The variable length declaration (*) can only be used for the last 
C....  dimension of an array! 
      double precision pstdaug(5,*)
      character*16 namestdaug(*)
      logical allfst
      logical enddecay
C... Local variables
      parameter (nenddecay=200)
      integer nstart(nenddecay),nstop(nenddecay),nstrstasto(10)
      integer kfstdaugw(4000),lnstdaugw(4000)
      double precision pstdaugw(5,1000)
      character*16 namestdaugw(1000)
C.....Three Pythia functions return integers, so need declaring.
      INTEGER PYCOMP
C.... PYTHIA COMMONS
C.....The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C.....Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      common/eventnumber/iev,nchilinejj
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
***********************************
      logical cluster
      common/debugcluster/cluster


      enddecay=.False.
      cluster=.False.
      newiline=iline
C...  Skipping informative part of the pythia table [k(i,1)=21]
C...  The same particle appears again in the table with the same
C...  flavor (K(iline,2)), and with origin iline 
      if(k(iline,1).eq.21)then 
         do i=iline+1,n
            if(k(i,3).eq.iline.and.k(i,2).eq.k(iline,2))newiline=i
         end do
      end if
C...  Check if initial particle is unstable
      if(k(newiline,1).lt.11)then
         write(*,*)'EROR:subroutine fulldecaychain:'
         write(*,*)'No daughters for stable particle: KF=',k(iline,2),
     $        ' in line: ',iline,' of table!'
         stop
      end if      
****************************************
C... The set of daughters of father, k(iline,2), appears in blocks in
C... the pythia table: from k(*,4) and k(*,5). In each block the initial
C... and the final position of the next block can be extracted.
C... Info of first block of daughters:
      nstart(1)=k(newiline,4)
      nstop(1)=k(newiline,5)
      if(iev.eq.13830)write(*,*)'filhos',iev,newiline,nstart(1),nstop(1)
************************************
C.... Preparing loop over number of blocks generated by newiline in table:
      nenddecaym=nenddecay !usually nenddecaym< 10 is enough. See below
      enddecay=.False.
      istdaug=0 ! initialize number of daughters
C...  Loop over each block of daugthers.
C.... The loop is exit when all the particles in the block are stable.
      do 208 ikik=1,nenddecaym

         if(.not.enddecay)then  ! 2081 if

*****************************
C... check number of strings in block
            nstring=0
      if(iev.eq.13830)write(*,*)'nstart ',ikik,nstart(ikik),nstop(ikik)
c      if(iev.eq.13830)write(*,*)'kstring ',k(nstart(ikik),4)
            kstring=k(nstart(ikik),4)
C... Check if the block contains a string:
            if(abs(k(k(nstart(ikik),4),2)).eq.92)then
      if(iev.eq.13830)write(*,*)'nstring < 10',nstring+1
               nstring=nstring+1
      if(iev.eq.13830)write(*,*)'nstrstato ',kstring
               nstrstasto(nstring)=kstring
            end if
C... Loop over the block ikik: 
            do iki=nstart(ikik),nstop(ikik)
      if(iev.eq.13830)write(*,*)'iki ',iki,nstart(ikik),nstop(ikik),
     $              ikik,iev
      if(iev.eq.13830)write(*,*)'kstringnew ',k(iki,4) !,k(502,4),k(502,1)
               kstringnew=k(iki,4)
C... Check if inside the block there is any string:
               if(abs(k(k(iki,4),2)).eq.92)then
                  if(kstring.ne.kstringnew)then
c      if(iev.eq.13830)write(*,*)'nstring ',nstring+1
                     nstring=nstring+1
c      if(iev.eq.13830)write(*,*)'kstring ',kstring
                     kstring=kstringnew
c      if(iev.eq.13830)write(*,*)'nstrstrato ',kstring
                     nstrstasto(nstring)=kstring
                  end if
               end if
C... Check if there are clusters. To be implemented as an special
C... case similar to various string. For the while just ignore this
C... events
               If(abs(k(k(iki,4),2)).eq.91)then
                  write(*,*)'cluster:',iev,nstart(ikik),nstop(ikik)
                  cluster=.True.
                  goto 2082 ! Exit loop
               end if

            enddo

            if(nstring.le.1)then ! end if at 888

C...Normal case: particles decays to 0 or 1 string
********************************************************
C...fulldecaychainsimple: obtain info from stable particles for block:
C..."ikik", whith "nstart(ikik)" and "nstop(ikik)".
C... Also give "nstart(ikik+1)" and "nstop(ikik+1)" with "nonstable.gt.0"
C... for the next ikik iteration,
C... OR "nonstable.eq.0" to exit loop. 
C... Where "nonstable" is number of unstable particles in each ikik block.
C... In each "ikik" iteration the number, of stable daughters: "istdaug" is
C... increased, in addition to all the other info from stable particles
               call fulldecaychainsimple(ikik,nstart,nstop,allfst,
     $              nonstable,istdaug,kfstdaug,lnstdaug,namestdaug,
     $              pstdaug)
               if(nonstable.eq.0)then
                  enddecay=.True.
                  goto 2082     ! Exit loop
               end if

            else ! nstring > 1
*******************************************
C... EXCEPTIONS: Place for special cases
*==================(I)====================
C...(I) Two or more strings involved
C.... print WARNING, search daughters of each string and exit loop
c               write(*,*)'unstable particle decays in',nstring,
c     $              ' strings at lines',(nstrstasto(iki),iki=1,nstring)
               nnstart=0
               istdaug=0
c               write(*,*)'finding the daughters of',nstring,
c     $              ' strings...'
C...Finding the daughters of each string: (nstring > 1)
               do 2010 jjj=1,nstring
****************************************
C...  The set of daughters of string at line, nstrstasto(jjj), appears in blocks in
C...  the pythia table: from k(*,4) and k(*,5). In each block the initial
C...  and the final position of the next block can be extracted.
C...  Info of first block of daughters:
                  nstart(1)=nstrstasto(jjj)
                  nstop(1)=nstrstasto(jjj)
************************************
                  nenddecaym=nenddecay !usually nenddecaym< 10 is enough. See below
                  enddecay=.False.
                  istdaugw=0    ! initialize number of daughters
C...  Loop over each block of daugthers.
C.... The loop is exit when all the particles in the block are stable.
                  if(jjj.gt.1)nnstart=nnstart+istdaug
                  do ikiki=1,nenddecaym
                     if(.not.enddecay)then 

                        call fulldecaychainsimple(ikiki,nstart,nstop,
     $                       allfst,nonstable,istdaugw,kfstdaugw,
     $                       lnstdaugw,namestdaugw,pstdaugw)
                        istdaug=istdaugw
                        if(nonstable.eq.0)then
                           enddecay=.True.
                           goto 2013 ! Exit loop
                        end if
                     end if           

                  end do
 2013             CONTINUE
C...  Set tau daughters output variables:
                  do ikiki=1,istdaug
                     kfstdaug(ikiki+nnstart)=kfstdaugw(ikiki)
                     lnstdaug(ikiki+nnstart)=lnstdaugw(ikiki)
                     namestdaug(ikiki+nnstart)=
     $                    namestdaugw(ikiki)
                     do ikikii=1,5 
                        pstdaug(ikikii,ikiki+nnstart)=
     $                       pstdaugw(ikikii,ikiki)
                     enddo
                  enddo

 2010          CONTINUE
               istdaug=nnstart+istdaug
               goto 2082        ! Exit loop
*==================(I)====================
********************************************************
            end if
         end if                 ! 2081 end if

 208  CONTINUE
 2082 CONTINUE
c      If(.not.cluster)write(*,*)'No cluster',iev
      end
***************************************************
C...  Full Decay Chain 
C...  Find all the daughters of the initial particle at line iline
C...    of Pythia table.
      subroutine fulldecaychainsimple(ikik,nstart,nstop,allfst,
     $     nonstable,istdaug,kfstdaug,lnstdaug,namestdaug,pstdaug)
*************************************************************
C...  INPUT
C...  ikik: number of each block. See main loop in parent subroutine:
C...        do 208 ikik=1,nenddecaym
C...        This loop continue if the analysed block contains unstable 
C...        particles (nonstable.gt.0) and exit if nonstable.eq.0
C...  nstart(ikik): Initial line in pythia table for ikik block
C...  nstop(ikik): Final line in pythia table for ikik block
C...  allfst: .False. only stable charged final states. 
C...           .True. All final stable states
C...  istdaug: number of stable particles of previous block
C...  OUTPUT: 
C...  nstart(ikik+1): Initial line in pythia table for 
C...  nstop(ikik+1): Final line in pythia table for next ikik iteration
C...  nonstable: Number of unstable particles in the ikik block
C...  istdaug: Cumulative number of final stable states (FSS) 
C...  kfstdaug(1)...kfstdaug(istdaug): Cumulative kf of FSS 
C...  lnstdaug(1)...lnstdaug(istdaug): Cumulative line number in table of FSS
C...  namestdaug(1)...namestdaug(istdaug): Cumulative names of FSS 
C...  pstdaug(1...5,1)...pstdaug(1...5,istdaug): Cumulative full momentum of FSS
C...  IT IS RECOMMENDED TO DEFINE THE VECTOR WITH A DIMENSION 
C.... SUFFICIENTLY HIGH, FOR EXAMPLE:
C.... kfstdaug(4000) for an initial quark
*****************************************************************
      implicit double precision(a-h,o-z)
      integer kfstdaug(*),lnstdaug(*)
c.... The variable length declaration (*) can only be used for the last 
C....  dimension of an array! 
      double precision pstdaug(5,*)
      character*16 namestdaug(*)
      logical allfst
      logical startd
C... Local variables
      integer nstart(*),nstop(*)
C.....Three Pythia functions return integers, so need declaring.
      INTEGER PYCOMP
C.... PYTHIA COMMONS
C.....The event record.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C.....Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      common/teste/ievv
      common/eventnumber/iev,nchilinejj
***********************************
      istaug=0
      startd=.True.
C..   debug:
*************************************
C...  Intialize next block:
      nstart(ikik+1)=0
      nstop(ikik+1)=0
*************************************
      nonstable=0 
c      write(*,*)'ikik in',ikik,nstart(ikik),nstart(ikik)
      do 209 iki=nstart(ikik),nstop(ikik)
C...  Warn if a string is found
         if(abs(k(iki,2)).eq.92)then
c            write(*,*)'******************************************'
c            write(*,*)'subroutine fulldecaychain:'
c            write(*,*)'Inital particle belongs to a string at i=',
c     $           iki,' in which should be other initial particles'
c            write(*,*)'******************************************'
         end if

***************************************
C...  Info of next block of daughters:
         if(k(iki,4).gt.0.and.startd)then
            nstart(ikik+1)=k(iki,4)
            startd=.False.
         end if
         stopc=k(iki,5)
         if(stopc.gt.nstop(ikik+1))nstop(ikik+1)=k(iki,5)
*******************************************************
***********************************************
C...  Extacting info from stable daughters
         if(k(iki,1).lt.11)then ! stable daughter
            if(allfst)then
               istdaug=istdaug+1
               kfstdaug(istdaug)=k(iki,2)
               call pyname(kfstdaug(istdaug),
     $              namestdaug(istdaug))
               lnstdaug(istdaug)=iki
               do jkjk=1,5
                  pstdaug(jkjk,istdaug)=
     $                 p(lnstdaug(istdaug),jkjk)
               end do
           else
               if(abs(kchg(pycomp(k(iki,2)),1)).gt.0)then
                  istdaug=istdaug+1
                  kfstdaug(istdaug)=k(iki,2)
                  lnstdaug(istdaug)=iki
                  call pyname(kfstdaug(istdaug),
     $                 namestdaug(istdaug))
                  do jkjk=1,5
                     pstdaug(jkjk,istdaug)=
     $                    p(lnstdaug(istdaug),jkjk)
                  end do
               end if
            end if
         else
            nonstable=nonstable+1
         endif
********************************************
 209  CONTINUE
      end
