      subroutine linfit(x,y,n,b,ferr) 
      implicit none 
      integer,intent(in):: n
      real(8),intent(in):: x(n),y(n)
      real(8),intent(out):: b,ferr 
      real(8):: mx,my,mxy,mx2,a
      mx=sum(x)/n
      my=sum(y)/n
      mx2=sum(x*x)/n
      mxy=sum(x*y)/n
      b=(mxy-mx*my)/(mx2-mx*mx)
      a=my-b*mx
      ferr=sqrt(sum((y-a-b*x)**2)/n)
      end subroutine 

