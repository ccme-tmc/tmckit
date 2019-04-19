      subroutine latdist(lattyp,latc,p1,p2,dist)
!! this subroutine calculate the distance between a pair of points in a lattice structure 
      implicit none
      character*4:: lattyp
      real*8::      latc(6),p1(3),p2(3),dist

      integer:: i
      real*8:: br2(3,3),help(3),dif(3)
      real*8:: pi,singam,cosgam 
      logical:: ortho

      pi=4.d0*atan(1.d0)
      cosgam=cos(latc(6)/180.d0*pi)
      singam=sin(latc(6)/180.d0*pi)

      call dlatgen(lattyp,latc(4:6),br2,ortho)
      dif=p2-p1

      do i=1,3
        if(dif(i).gt.0.5) then
          dif(i)=dif(i)-1.0
        elseif(dif(i).lt.-0.5) then
          dif(i)=dif(i)+1.0
        endif 
      enddo

      if(.not.ortho) then 
        help=dif
        if(lattyp(1:1).eq.'R') then
          dif(1)=help(1)*BR2(1,1)+help(2)*BR2(2,1)+help(3)*BR2(3,1)
          dif(2)=help(1)*BR2(1,2)+help(2)*BR2(2,2)+help(3)*BR2(3,2)
          dif(3)=help(1)*BR2(1,3)+help(2)*BR2(2,3)+help(3)*BR2(3,3)
        elseif(lattyp(1:3).eq.'CXZ') then
          dif(1)=help(1)*singam
          dif(2)=(help(1)*cosgam*latc(1)+help(2)*latc(2))/latc(2)
          dif(3)=help(3)
        else
          dif(1)=(help(1)*BR2(1,1)*latc(1)+help(2)*BR2(2,1)*latc(2)+    &
     &           help(3)*BR2(3,1)*latc(3))/latc(1)
          dif(2)=(help(1)*BR2(1,2)*latc(1)+help(2)*BR2(2,2)*latc(2)+    &
     &           help(3)*BR2(3,2)*latc(3))/latc(2)
          dif(3)=(help(1)*BR2(1,3)*latc(1)+help(2)*BR2(2,3)*latc(2)+    &
     &           help(3)*BR2(3,3)*latc(3))/latc(3)
        endif
      endif 
      dist=sqrt( sum(  (dif(1:3)*latc(1:3))**2 ) )
      end subroutine 
