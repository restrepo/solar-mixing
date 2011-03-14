*********************************************************************
******  MSSM constraints from LEP  -  ICHEP2000 P. Igo-Kemenes  *****
*********************************************************************

      subroutine lepcuts

      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      include 'pylocalcom.h'

*******************-- LOCAL VARIABLES --********************************
      CHARACTER CHAP*16
************************************************************************

c.....LSP charged or colored

      write(*,*)'stop = ',PMAS(PYCOMP(1000006),1)
      write(*,*)'g~   = ',PMAS(PYCOMP(1000021),1)
      write(*,*)'ur~  = ',PMAS(PYCOMP(2000002),1)
      write(*,*)'dr~  = ',PMAS(PYCOMP(2000001),1)

      open(21,file='result.out',status='unknown')

      if(PMAS(PYCOMP(1000024),1).lt.PMAS(PYCOMP(1000022),1))then
         write(21,*)'point theoretically excluded.'
         write(21,*)'chargino is the LSP - M = ',PMAS(PYCOMP(1000024),1)
         stop
      endif

      if(PMAS(PYCOMP(1000005),1).lt.PMAS(PYCOMP(1000022),1))then
         write(21,*)'point theoretically excluded.'
         write(21,*)'sbottom is the LSP - M = ',PMAS(PYCOMP(1000005),1)
         stop
      endif

      if(PMAS(PYCOMP(1000006),1).lt.PMAS(PYCOMP(1000022),1))then
         write(21,*)'point theoretically excluded.'
         write(21,*)'stop is the LSP - M = ',PMAS(PYCOMP(1000006),1)
         stop
      endif

      if(PMAS(PYCOMP(1000015),1).lt.PMAS(PYCOMP(1000022),1))then
         write(21,*)'point theoretically excluded.'
         write(21,*)'stau is the LSP - M = ',PMAS(PYCOMP(1000015),1)
         stop
      endif

c..... LEP limits

      CALL PYNAME(1000006,CHAP)
      if(PMAS(PYCOMP(1000006),1).lt.95.d0)then
         write(21,*)'SUSY point excluded by LEP.'         
         write(21,*)CHAP,PMAS(PYCOMP(1000006),1) ! ~t_1 mass
         stop 
      endif

      CALL PYNAME(1000005,CHAP) 
      if(PMAS(PYCOMP(1000005),1).lt.85.d0)then
         write(21,*)'SUSY point excluded by LEP.'         
         write(21,*)CHAP,PMAS(PYCOMP(1000005),1) ! ~b_1 mass
         stop
      endif

      CALL PYNAME(1000024,CHAP) 
      if(PMAS(PYCOMP(1000024),1).lt.95.d0)then
         write(21,*)'SUSY point excluded by LEP.'         
         write(21,*)CHAP,PMAS(PYCOMP(1000024),1) ! ~chi_1+ mass
         stop
      endif

      CALL PYNAME(1000022,CHAP) 
      if(PMAS(PYCOMP(1000022),1).lt.42.d0)then
         write(21,*)'SUSY point excluded by LEP.'         
         write(21,*)CHAPP,PMAS(PYCOMP(1000022),1) ! ~chi_10 mass
c         stop
      endif

      CALL PYNAME(1000015,CHAP) 
      if(PMAS(PYCOMP(1000015),1).lt.79.d0)then
         write(21,*)'SUSY point excluded by LEP.'         
         write(21,*)CHAP,PMAS(PYCOMP(1000015),1) ! ~tau_1- mass   
         stop
      endif

      CALL PYNAME(25,CHAP) 
      if(PMAS(PYCOMP(25),1).lt.95.d0)then  ! Tata's limit
         write(21,*)'SUSY point excluded by LEP.'         
         write(21,*)CHAP,PMAS(PYCOMP(25),1) ! h0 mass
         stop
      endif

      close(21)

************************************************************************
      return
      end
************************************************************************













