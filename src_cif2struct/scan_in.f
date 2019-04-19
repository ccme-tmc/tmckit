 subroutine scan_in(filein,nato)

   use struct, only : tobohr,lattyp,sgnum,HRtransf,hex2rho,rho2hex

   integer       nato
   character*80  filein,line
   real*8        orig(3),x(3,nato)
   real*8        a,b,c,alfa,beta,gamma
   real*8        aa,bb,cc
   character*50  sgname,aname(nato)
   integer       n,nat
   real*8        pi

      pi=acos(-1.0d0)

      rewind 1
      open(unit=1,file=filein,status='old')
      read(1,'(a)') line
      do i=1,80
         if (line(i:i).eq.'b') then
            tobohr=1.0d0
            goto 1
         endif
         if (line(i:i).eq.'a') then
            tobohr=1.0d0/0.5291772083d0
            goto 1
         endif
      enddo
1     continue
      n=0
      read(1,*) (orig(i),i=1,3)
      read(1,*) a,b,c,alfa,beta,gamma
      read(1,*) sgname
2     n=n+1
      read(1,*,end=3) aname(n)
      if(aname(n).eq.' ') go to 3
      read(1,*) (x(i,n),i=1,3) 
      goto 2
3     continue
      close (1)

      nat=n-1

     call test_sgname(sgname,.false.)

     if (sgnum.eq.167.or.sgnum.eq.166.or.sgnum.eq.161.or.sgnum.eq.160.or.&
         sgnum.eq.155.or.sgnum.eq.148.or.sgnum.eq.146) then
        HRtransf=.true.
        small=1.0d-6
        if ((abs(alfa-90.0d0).lt.small).and.&  
            (abs(beta-90.0d0).lt.small).and.&  
            (abs(gamma-120.0d0).lt.small)) then
!           do nothing
        else
           if (lattyp.ne.'R  ') then
              write(*,'(2a)') 'lattyp should be R and is: ',trim(lattyp)
              stop
           endif   
           alfa=alfa*pi/180.0d0
           aa=a*2.0d0*cos((pi-alfa)/2.0d0)
           bb=aa
           cc=3.0d0*sqrt(a**2-(aa**2)/3.0d0)
           a=aa
           b=bb
           c=cc
           alfa=90.0d0
           beta=90.0d0
           gamma=120.0d0
           do i=1,nat
              x(1:3,i)=matmul(x(1:3,i),rho2hex)
              do j=1,3
                 if (x(j,i).lt.0.0d0) x(j,i)=x(j,i)+1.0d0
                 if (x(j,i).gt.1.0d0) x(j,i)=x(j,i)-1.0d0
              enddo
           enddo

        endif
     endif

        open(unit=1,status='scratch')        
!     open(unit=1,file='rolask')        
     write(1,'(3f12.8)') orig(1:3)
     write(1,'(6f14.8)') a,b,c,alfa,beta,gamma
     sgname=''''//trim(sgname)//''''
     write(1,'(a)') sgname
     do i=1,nat
        call insert_string(aname(i),'     ')
        aname(i)=''''//trim(aname(i))//''''         
        write(1,'(a)') aname(i)
        write(1,'(3f18.8)') (x(j,i),j=1,3)  
     enddo
     rewind 1


      end subroutine scan_in
