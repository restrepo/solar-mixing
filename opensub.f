
      subroutine opensub        !output:manysub(thissub): thissub 1 to n

      implicit none 
      integer msel,mselpd,msub,kfin
      double precision ckin
      common/pysubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      integer subpr(500),nsub,i
      common/mayprocess/subpr,nsub

      open(9,file='isub.dat',status='unknown')
      nsub=1
 333  read(9,*,end=8888)subpr(nsub)
      msub(subpr(nsub))=1
      nsub=nsub+1
      goto 333
 8888 continue
      nsub=nsub-1
      write(*,*)'subprocess switched on: ',(subpr(i),i=1,nsub)
      return
      end


