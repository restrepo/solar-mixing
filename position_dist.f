      subroutine position_dist(v1,v2,p1,p2,d2,vpos1,vpos2)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      double precision v1(3),v2(3),p1(3),p2(3),z1i(2),z2i(4)
      double precision x1i(2),y1i(2),x2i(4),y2i(4)
      double precision vpos1(3),vpos2(3)
c...  I have checked that the four solutions are equal
      p1x=p1(1);p1y=p1(2);p1z=p1(3) ! u   |PQ*v|
      p2x=p2(1);p2y=p2(2);p2z=p2(3) ! v d=------  
      xp1=v1(1);yp1=v1(2);zp1=v1(3) ! P    |v|
      xp2=v2(1);yp2=v2(2);zp2=v2(3) ! Q
      a= 2*p1z**2*p2x*p2z*xp1 - 2*p1z**2*p2x*p2z*xp2 + 
     -  2*p1z**2*p2y*p2z*yp1 - 2*p1z**2*p2y*p2z*yp2 + 
     -  2*p1y**2*(p2x**2 + p2z**2)*zp1 + 
     -  2*p1x**2*(p2y**2 + p2z**2)*zp1 + 
     -  2*p1z**2*p2x**2*zp2 + 2*p1z**2*p2y**2*zp2 - 
     -  2*p1x*p1z*(p2y**2*(xp1 - xp2) + 
     -     p2x*p2y*(-yp1 + yp2) + 
     -     p2z*(p2z*(xp1 - xp2) + p2x*(zp1 + zp2))) - 
     -  2*p1y*(2*p1x*p2x*p2y*zp1 + 
     -     p1z*(p2x*p2y*(-xp1 + xp2) + p2x**2*(yp1 - yp2) + 
     -        p2z*(p2z*(yp1 - yp2) + p2y*(zp1 + zp2))))
      b=-4*(p1z**2*(p2x**2 + p2y**2) - 2*p1x*p1z*p2x*p2z - 
     -     2*p1y*p2y*(p1x*p2x + p1z*p2z) + 
     -     p1y**2*(p2x**2 + p2z**2) + 
     -     p1x**2*(p2y**2 + p2z**2))*
     -   (-(d2*p1z**2*(p2x**2 + p2y**2 + p2z**2)) + 
     -     (-2*p1x*p1y*p2x*p2y + p1y**2*(p2x**2 + p2z**2) + 
     -        p1x**2*(p2y**2 + p2z**2))*zp1**2 + 
     -     p1z**2*(p2z**2*
     -         (xp1**2 - 2*xp1*xp2 + xp2**2 + 
     -           (yp1 - yp2)**2) + 
     -        2*p2x*p2z*(xp1 - xp2)*zp2 - 
     -        2*p2y*(yp1 - yp2)*
     -         (p2x*(xp1 - xp2) - p2z*zp2) + 
     -        p2y**2*(xp1**2 - 2*xp1*xp2 + xp2**2 + 
     -           zp2**2) + 
     -        p2x**2*(yp1**2 - 2*yp1*yp2 + yp2**2 + zp2**2))
     -       - 2*p1z*zp1*
     -      (p1x*(p2y**2*(xp1 - xp2) + 
     -           p2x*p2y*(-yp1 + yp2) + 
     -           p2z*(p2z*xp1 - p2z*xp2 + p2x*zp2)) + 
     -        p1y*(p2x*p2y*(-xp1 + xp2) + 
     -           p2x**2*(yp1 - yp2) + 
     -           p2z*(p2z*yp1 - p2z*yp2 + p2y*zp2)))) + 
     -  4*(p1y**2*(p2x**2 + p2z**2)*zp1 + 
     -      p1x**2*(p2y**2 + p2z**2)*zp1 + 
     -      p1z**2*(p2x*p2z*(xp1 - xp2) + p2x**2*zp2 + 
     -         p2y*(p2z*yp1 - p2z*yp2 + p2y*zp2)) - 
     -      p1x*p1z*(p2y**2*(xp1 - xp2) + 
     -         p2x*p2y*(-yp1 + yp2) + 
     -         p2z*(p2z*(xp1 - xp2) + p2x*(zp1 + zp2))) - 
     -      p1y*(2*p1x*p2x*p2y*zp1 + 
     -         p1z*(p2x*p2y*(-xp1 + xp2) + 
     -            p2x**2*(yp1 - yp2) + 
     -            p2z*(p2z*(yp1 - yp2) + p2y*(zp1 + zp2)))))
     -     **2
      c=2*(p1z**2*(p2x**2 + p2y**2) - 2*p1x*p1z*p2x*p2z - 
     -    2*p1y*p2y*(p1x*p2x + p1z*p2z) + 
     -    p1y**2*(p2x**2 + p2z**2) + 
     -    p1x**2*(p2y**2 + p2z**2))
      if(b.eq.0)then
         aoverb=1d-5
      else
         aoverb=abs(a/b)
      endif
      if(aoverb.gt.1d-6)then
         if(b.lt.0)b=0d0
         z1i(1)=(a+dsqrt(b))/c
         z1i(2)=(a-dsqrt(b))/c
      else
         write(*,*)'ERROR: z1 is complex',b/a
      endif

      ii=0
c.....I will take only one solution
      do i=1,1  !2
         z1=z1i(i)
         ap=-2*p1z**2*p2x*p2z*xp1 + 2*p1z**2*p2x*p2z*xp2 - 
     -        2*p1z**2*p2y*p2z*yp1 + 2*p1z**2*p2y*p2z*yp2 - 
     -        2*p1x*p1z*p2x*p2z*z1 - 2*p1y*p1z*p2y*p2z*z1 - 
     -        2*p1z**2*p2z**2*z1 + 2*p1x*p1z*p2x*p2z*zp1 + 
     -        2*p1y*p1z*p2y*p2z*zp1 - 2*p1z**2*p2x**2*zp2 - 
     -        2*p1z**2*p2y**2*zp2
         bp=(2*p1z**2*p2x*p2z*xp1 - 2*p1z**2*p2x*p2z*xp2 + 
     -        2*p1z**2*p2y*p2z*yp1 - 2*p1z**2*p2y*p2z*yp2 + 
     -        2*p1x*p1z*p2x*p2z*z1 + 2*p1y*p1z*p2y*p2z*z1 + 
     -        2*p1z**2*p2z**2*z1 - 2*p1x*p1z*p2x*p2z*zp1 - 
     -        2*p1y*p1z*p2y*p2z*zp1 + 2*p1z**2*p2x**2*zp2 + 
     -        2*p1z**2*p2y**2*zp2)**2 - 
     -        4*(-(p1z**2*p2x**2) - p1z**2*p2y**2 - 
     -        p1z**2*p2z**2)*
     -        (d2*p1z**2*p2z**2 - p1z**2*p2z**2*xp1**2 + 
     -        2*p1z**2*p2z**2*xp1*xp2 - p1z**2*p2z**2*xp2**2 - 
     -        p1z**2*p2z**2*yp1**2 + 2*p1z**2*p2z**2*yp1*yp2 - 
     -        p1z**2*p2z**2*yp2**2 - 2*p1x*p1z*p2z**2*xp1*z1 + 
     -        2*p1x*p1z*p2z**2*xp2*z1 - 
     -        2*p1y*p1z*p2z**2*yp1*z1 + 
     -        2*p1y*p1z*p2z**2*yp2*z1 - p1x**2*p2z**2*z1**2 - 
     -        p1y**2*p2z**2*z1**2 - p1z**2*p2z**2*z1**2 + 
     -        2*p1x*p1z*p2z**2*xp1*zp1 - 
     -        2*p1x*p1z*p2z**2*xp2*zp1 + 
     -        2*p1y*p1z*p2z**2*yp1*zp1 - 
     -        2*p1y*p1z*p2z**2*yp2*zp1 + 
     -        2*p1x**2*p2z**2*z1*zp1 + 
     -        2*p1y**2*p2z**2*z1*zp1 - p1x**2*p2z**2*zp1**2 - 
     -        p1y**2*p2z**2*zp1**2 - 
     -        2*p1z**2*p2x*p2z*xp1*zp2 + 
     -        2*p1z**2*p2x*p2z*xp2*zp2 - 
     -        2*p1z**2*p2y*p2z*yp1*zp2 + 
     -        2*p1z**2*p2y*p2z*yp2*zp2 - 
     -        2*p1x*p1z*p2x*p2z*z1*zp2 - 
     -        2*p1y*p1z*p2y*p2z*z1*zp2 + 
     -        2*p1x*p1z*p2x*p2z*zp1*zp2 + 
     -        2*p1y*p1z*p2y*p2z*zp1*zp2 - 
     -        p1z**2*p2x**2*zp2**2 - p1z**2*p2y**2*zp2**2)
         cp=2*(-(p1z**2*p2x**2) - p1z**2*p2y**2 - p1z**2*p2z**2)
         x1i(i)=-((-(p1z*xp1) - p1x*z1 + p1x*zp1)/p1z)
         y1i(i)=-((-(p1z*yp1) - p1y*z1 + p1y*zp1)/p1z)
         if(abs(ap/bp).gt.1d-6)then
            if(bp.lt.0)bp=0d0
            ii=ii+1
            z2i(ii)=(ap+dsqrt(bp))/cp
            z2=z2i(ii)
            x2i(ii)=-((-(p2z*xp2) - p2x*z2 + p2x*zp2)/p2z)
            y2i(ii)=-((-(p2z*yp2) - p2y*z2 + p2y*zp2)/p2z)

            ii=ii+1
            z2i(ii)=(ap-dsqrt(bp))/cp
            z2=z2i(ii)
            x2i(ii)=-((-(p2z*xp2) - p2x*z2 + p2x*zp2)/p2z)
            y2i(ii)=-((-(p2z*yp2) - p2y*z2 + p2y*zp2)/p2z)
         else
            ii=ii+1
            write(*,*)'ERROR: z2',ii,ii+1,' is complex',ap/bp
         endif
      enddo
      ii=0
c      do i=1,2
c         x1=x1i(i);y1=y1i(i);z1=z1i(i) 
c         ii=ii+1
c         x2=x2i(ii);y2=y2i(ii);z2=z2i(ii) 
c         dist=dsqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
c         write(*,*)'z1',ii,dist,dsqrt(d2)
c         ii=ii+1
c         x2=x2i(ii);y2=y2i(ii);z2=z2i(ii) 
c         dist=dsqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
c         write(*,*)'z1',ii,dist,dsqrt(d2)
c      enddo
      x1=x1i(1);y1=y1i(1);z1=z1i(1) 
      x2=x2i(1);y2=y2i(1);z2=z2i(1) 
      vpos1(1)=x1;vpos1(2)=y1;vpos1(3)=z1
      vpos2(1)=x2;vpos2(2)=y2;vpos2(3)=z2
      end
