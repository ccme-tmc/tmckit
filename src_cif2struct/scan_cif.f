 subroutine scan_cif(filein,usespg,nato)

   use struct, only: lattyp,sgnum,HRtransf,hex2rho,rho2hex

   include       'ciftbx.cmn'

   logical usespg
   character*80 filein
   logical       f1,f2,f3,f4
   character*32  sgname,sseting
   real*8          a,b,c,siga,sigb,sigc
   real*8          alfa,beta,gamma
   real*8          x(3,nato),sx,sy,sz
   integer       i,j,nat,nsym,isgnum,n
   character*10  aname(nato)
   character*100 symname(nato)
   real*8      org(3),fract,small,sgnum1
   real*8      aa,bb,cc
   real*8      pi
   logical RHtransf,RHexeption


     pi=acos(-1.0d0)

     f1 = init_(1,2,3,6)  
     if(.not.ocif_(filein)) then
        write(*,'(a///)')  ' >>>>>>>>> CIF cannot be opened'
        stop
     endif

     if(.not.data_(' ')) then
        write(*,'(/a/)')   ' >>>>>>> No data_ statement found'
        stop
     endif 
     
     f1 = numd_('_cell_length_a', a, siga)
     f2 = numd_('_cell_length_b', b, sigb)
     f3 = numd_('_cell_length_c', c, sigc)
     if(.not.(f1.and.f2.and.f3)) then
        write(*,'(a)') ' Cell dimension(s) missing!'
        write(*,*) f1,f2,f3
        stop
     endif  

     f1 = numd_('_cell_angle_alpha', alfa, siga)
     f2 = numd_('_cell_angle_beta', beta, sigb)
     f3 = numd_('_cell_angle_gamma', gamma, sigc)
     if(.not.(f1.and.f2.and.f3)) then
        write(*,'(a)') ' Cell angle(s) missing!'
        stop
     endif  

     i=0
160  continue
     i=i+1
     f1 = char_('_atom_site_type_symbol', aname(i))
     if (.not.f1) f1 = char_('_atom_site_label', aname(i))

     f2 = numd_('_atom_site_fract_x', x(1,i), sx)
     f3 = numd_('_atom_site_fract_y', x(2,i), sy)
     f4 = numd_('_atom_site_fract_z', x(3,i), sz)
     if(.not.(f1.and.f2.and.f3.and.f4)) then
        write(*,'(a)') ' Same of the atomic data missing!'
        stop
     endif  
     if(loop_)  goto 160   
     nat=i

     f1 = char_('_symmetry_cell_setting', sseting)
     if (f1.and.sseting.eq.'hexagonal') then
        small=2.0d-5 
        do i=1,nat
           do j=1,2
              do jfract=1,5
!              do fract=1.0d0/6.0d0,5.0d0/6.0d0,1.0d0/6.0d0
                 fract=jfract/6.d0
                 if (abs(x(j,i)-fract).lt.small) x(j,i)=fract
                 if (abs(x(j,i)+fract).lt.small) x(j,i)=-fract
              enddo
           enddo
        enddo
     endif

     sgnum1=0.0d0 
     sgname='                                            '

     f1 = numd_('_symmetry_Int_Tables_number', sgnum1,siga)     
     f2 = char_('_symmetry_space_group_name_H-M', sgname)

     i=0     
170  continue
     i=i+1
     f3 = char_('_symmetry_equiv_pos_as_xyz', symname(i))
     if(loop_)  goto 170   
     nsym=i

     if (.not.(f3.or.f2)) then
        write(*,'(a)') ' Space group name nor symmetry operations are not given!'
        stop
     endif

     isgnum=int(sgnum1)
!     if (.not.f2) call  getsgname(isgnum,sgname)
     call test_sgname(sgname,f3)
     if (lattyp.eq.'   ') lattyp='P  '
     if (f1.and.(.not.f3).and.(isgnum.ne.sgnum)) then
        write(*,'(a,i4,a,i4,2a)') 'Space group number from cif:',isgnum,&
             ' and from internal tables:',sgnum,' are different!',& 
             ' Check space group name.'
        stop
     endif
     
     if (isgnum.ne.0) sgnum=isgnum

     RHtransf=.false.
     RHexeption=.false.
     if (sgnum.eq.167.or.sgnum.eq.166.or.sgnum.eq.161.or.&
          sgnum.eq.160.or.sgnum.eq.155.or.sgnum.eq.167.or.&
          sgnum.eq.148.or.sgnum.eq.146) then
        small=1.0d-6
        if ((abs(alfa-90.0d0).lt.small).and.&  
            (abs(beta-90.0d0).lt.small).and.&  
            (abs(gamma-120.0d0).lt.small)) then
           lattyp='H  '
           RHexeption=.true.
        else
           lattyp='R  '
           RHtransf=.true.
           aa=a*2.0d0*cos((pi-alfa*pi/180.0d0)/2.0d0)
           bb=aa
           cc=3.0d0*sqrt(a**2-(aa**2)/3.0d0)
           a=aa
           b=bb
           c=cc
           alfa=90.0d0
           beta=90.0d0
           gamma=120.0d0
        endif
     endif

     if (.not.f3) then
        if (RHexeption) HRtransf=.true.
        if (RHtransf) then
           HRtransf=.true.
           do i=1,nat
              x(1:3,i)=matmul(x(1:3,i),rho2hex)
              do j=1,3
                 if (x(j,i).lt.0.0d0) x(j,i)=x(j,i)+1.0d0
                 if (x(j,i).gt.1.0d0) x(j,i)=x(j,i)-1.0d0
              enddo
           enddo
        endif
        org(1:3)=0.0d0
        f1 = numd_('_symmetry_origin_shift_x', org(1),siga)
        f2 = numd_('_symmetry_origin_shift_y', org(2),siga)
        f3 = numd_('_symmetry_origin_shift_z', org(3),siga)
        usespg=.true.
        open(unit=1,status='scratch')        
!        open(unit=1,file='rolask')        
        write(1,'(3f14.8)') org(1:3)
        write(1,'(6f14.8)') a,b,c,alfa,beta,gamma
        sgname=''''//trim(sgname)//''''
        write(1,'(a)') sgname
        do i=1,nat
           call insert_string(aname(i),'  ')
           aname(i)=''''//trim(aname(i))//''''         
           write(1,'(a)') aname(i)
           write(1,'(3f14.8)') (x(j,i),j=1,3)  
        enddo
        rewind 1
     else
        call gen_struct_cif(nat,x,aname,nsym,symname,sgname,&
                            a,b,c,alfa,beta,gamma)
        usespg=.false. 
     endif

 end subroutine scan_cif

 subroutine insert_string(aname,str)
 
   character aname*(*),str*(*)
   integer n,i,m
   character*12 numbers
   
   numbers='0123456789+-'
   m=len_trim(aname)
   i1=scan(aname,numbers,.false.)
   i2=scan(aname,numbers,.true.)
!     write(*,*) 'aname:',aname,i1,i2
   if (i1.ne.0) aname=aname(1:i1-1)//str//aname(i1:i2)

 end subroutine insert_string

 subroutine gen_struct_cif(nat1,x,aname1,nsym,symname,sgname1,&
                           a,b,c,alfa,beta,gamma)

   use struct, only: nat,sgname,pos,aname,aa,bb,cc,alpha,mult,lattyp

   character*32  sgname1
   real*8          a,b,c
   real*8          alfa,beta,gamma
   real*8          x(3,1000)
   integer       i,nat1,nsym
   character*10  aname1(*)
   character*100 symname(*)
   real*8,allocatable :: xn(:,:)
   integer,allocatable :: indequiv(:)
   real*8        x1(3),pi
   integer j,index
 
   pi=acos(-1.0d0)

   allocate (xn(3,nsym),indequiv(nsym))

   aa=a
   bb=b
   cc=c
   alpha(1)=alfa
   alpha(2)=beta
   alpha(3)=gamma

   nat=nat1
   
   index=0
   do i=1,nat
      aname(i)=aname1(i)
      call insert_string(aname(i),'  ')
      x1(1:3)=x(1:3,i) 
      do j=1,nsym
         call apply_cifsymop(symname(j),x1,xn(1:3,j),j)         
      enddo
      call gen_equiv(nsym,xn,lattyp,indequiv,mult(i))
      do j=1,mult(i)
         index=index+1
         pos(1:3,index)=xn(1:3,indequiv(j))       
      enddo
   enddo
      
 end subroutine gen_struct_cif

 subroutine apply_cifsymop(symname,x,xn,irec)

   integer irec
   character symname*(*)
   real*8    x(3),xn(3)
   integer n,i,j,lp
   character*10 part(3)
   real*8 a,b
   character symname1*50

   n=len_trim(symname)
   do jj=1,n
        if(symname(jj:jj).eq.'X')symname(jj:jj)='x'
        if(symname(jj:jj).eq.'Y')symname(jj:jj)='y'
        if(symname(jj:jj).eq.'Z')symname(jj:jj)='z'
   enddo

   i=index(symname,',',.false.)
   j=index(symname,',',.true.)
   part(1)=adjustl(symname(1:i-1))
   part(2)=adjustl(symname(i+1:j-1))
   part(3)=adjustl(symname(j+1:n))

!!$   if (irec.eq.1) then       
!!$      if (trim(part(1)).ne.'x'.and.&
!!$          trim(part(2)).ne.'y'.and.&
!!$          trim(part(3)).ne.'z') then
!!$         write(*,'(/,a)') 'First symmetry operation in _symmetry_equiv_pos_as_xyz should be unity'
!!$         write(*,'(a)')symname
!!$         stop 
!!$      endif
!!$   endif

   do i=1,3

      xn(i)=0.0d0  
      lp=0

      j=index(part(i),'-x')
      if (j.ne.0) then
         xn(i)=xn(i)-x(1)
         lp=lp+2 
      else
         j=index(part(i),'x')
         if (j.ne.0) then
            xn(i)=xn(i)+x(1)
            lp=lp+1 
            if (j.ne.1) lp=lp+1
         endif
      endif
      j=index(part(i),'-y')
      if (j.ne.0) then
         xn(i)=xn(i)-x(2) 
         lp=lp+2 
      else
         j=index(part(i),'y')
         if (j.ne.0) then
            xn(i)=xn(i)+x(2)
            lp=lp+1
            if (j.ne.1) lp=lp+1
         endif
      endif
      j=index(part(i),'-z')
      if (j.ne.0) then
         xn(i)=xn(i)-x(3) 
         lp=lp+2 
      else
         j=index(part(i),'z')
         if (j.ne.0) then
            xn(i)=xn(i)+x(3)
            lp=lp+1
            if (j.ne.1) lp=lp+1
         endif
      endif
      j=index(part(i),'/')
      if (j.ne.0) then
         if (j.eq.1.or.j.eq.len_trim(part(i))) go to 1
         read(part(i)(j-1:j-1),*,err=1) a
         read(part(i)(j+1:j+1),*,err=1) b
         if (j.eq.2) then
            xn(i)=xn(i)+a/b                  
            lp=lp+3
         else 
            if (part(i)(j-2:j-2).eq.'-') then
               xn(i)=xn(i)-a/b                  
               lp=lp+4
            else
               xn(i)=xn(i)+a/b                  
               lp=lp+4
            endif
         endif
         
      endif

      if (lp.ne.len_trim(part(i))) goto 1  
  
   enddo

   do i=1,3
      if (xn(i).lt.0.0d0) xn(i)=xn(i)+1.0d0
      if (xn(i).gt.1.0d0) xn(i)=xn(i)-1.0d0
   enddo
   
   do i=1,3
      if (xn(i).lt.0.0d0) xn(i)=xn(i)+1.0d0
      if (xn(i).gt.1.0d0) xn(i)=xn(i)-1.0d0
   enddo
   
   do i=1,3
      if (xn(i).lt.0.0d0) xn(i)=xn(i)+1.0d0
      if (xn(i).gt.1.0d0) xn(i)=xn(i)-1.0d0
   enddo
   
   return
1  continue
   write(*,'(/,a,i4,a,i2)') 'wrong syntax in _symmetry_equiv_pos_as_xyz:    record '&
        ,irec,'     component',i
   stop   

 end subroutine apply_cifsymop

 subroutine gen_equiv(nsym,xn,lattyp,indequiv,mult)

   integer nsym
   real*8  xn(3,nsym)
   character lattyp*(*)
   integer   indequiv(nsym),mult
   integer i,j,k
   logical thesame    

   mult=0
   do i=1,nsym
      do j=1,mult
         k=indequiv(j)
         if (thesame(lattyp,xn(1:3,k),xn(1:3,i))) goto 1
      enddo
      mult=mult+1
      indequiv(mult)=i
1     continue
   enddo

 end subroutine gen_equiv

 logical function thesame(lattyp,x1,x2)

    character lattyp*(*)
    real*8    x1(3),x2(3)
    real*8 tv(3,1000),dx(3),dxl,small,small2
    integer ind,i,j,k


    small=1.0d-5
    small2=small/2.0d0

    ind=0
    do i=-1,1
       do j=-1,1
          do k=-1,1
             ind=ind+1
             tv(1,ind)=i
             tv(2,ind)=j
             tv(3,ind)=k
          enddo 
        enddo
     enddo

     if (lattyp(1:1).eq.'B') then

        do i=-1,1,2
           do j=-1,1,2
              do k=-1,1,2
                 ind=ind+1
                 tv(1,ind)=0.5d0*i
                 tv(2,ind)=0.5d0*j
                 tv(3,ind)=0.5d0*k
              enddo
           enddo
        enddo

     else if (lattyp(1:1).eq.'F') then

        do i=-1,1,2
           do j=-1,1,2
              ind=ind+1
              tv(1,ind)=0.5d0*i
              tv(2,ind)=0.5d0*j
              tv(3,ind)=0.0d0
           enddo
        enddo
        do i=-1,1,2
           do j=-1,1,2
              ind=ind+1
              tv(1,ind)=0.5d0*i
              tv(2,ind)=0.0d0
              tv(3,ind)=0.5d0*j
           enddo
        enddo
        do i=-1,1,2
           do j=-1,1,2
              ind=ind+1
              tv(1,ind)=0.0d0
              tv(2,ind)=0.5d0*j
              tv(3,ind)=0.5d0*j
           enddo
        enddo

     else if (lattyp(1:3).eq.'CXY') then
          
        do i=-1,1,2
           do j=-1,1,2
              ind=ind+1
              tv(1,ind)=0.5d0*i
              tv(2,ind)=0.5d0*j
              tv(3,ind)=0.0d0
           enddo
        enddo

     else if (lattyp(1:3).eq.'CXZ') then

        do i=-1,1,2
           do j=-1,1,2
              ind=ind+1
              tv(1,ind)=0.5d0*i
              tv(2,ind)=0.0d0
              tv(3,ind)=0.5d0*j
           enddo
        enddo

     else if (lattyp(1:3).eq.'CYZ') then

        do i=-1,1,2
           do j=-1,1,2
              ind=ind+1
              tv(1,ind)=0.0d0
              tv(2,ind)=0.5d0*j
              tv(3,ind)=0.5d0*j
           enddo
        enddo

     endif

     thesame=.false.
     do i=1,ind
        dx(1:3)=x1(1:3)-x2(1:3)-tv(1:3,i)
        dx=dx+10.0d0
!!$        do j=1,3
!!$           if (dx(j).lt.0.0d0) dx(j)=dx(j)+1.0d0
!!$           if (dx(j).gt.1.0d0) dx(j)=dx(j)-1.0d0
!!$        enddo
        dx=mod(dx+small2,1.0d0)-small2
        if (abs(dx(1)-1.0d0).lt.small2) dx(1)=0.0d0
        if (abs(dx(2)-1.0d0).lt.small2) dx(2)=0.0d0
        if (abs(dx(3)-1.0d0).lt.small2) dx(3)=0.0d0

        dxl=sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
        if (dxl.lt.small) then
           thesame=.true.
           return
        endif
     enddo

  end function thesame

  subroutine inversa(a,ainv)
    implicit real*8 (a-h,o-z)
    dimension a(3,3),ainv(3,3)
    det= a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1) &
         +a(1,3)*a(2,1)*a(3,2)-a(3,1)*a(2,2)*a(1,3) &
         -a(1,1)*a(3,2)*a(2,3)-a(2,1)*a(1,2)*a(3,3)
    ainv(1,1) =(   a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
    ainv(2,1) =( - a(2,1) * a(3,3) + a(2,3) * a(3,1) ) / det
    ainv(3,1) =(   a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
    ainv(1,2) =( - a(1,2) * a(3,3) + a(1,3) * a(3,2) ) / det
    ainv(2,2) =(   a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
    ainv(3,2) =( - a(1,1) * a(3,2) + a(1,2) * a(3,1) ) / det
    ainv(1,3) =(   a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det
    ainv(2,3) =( - a(1,1) * a(2,3) + a(1,3) * a(2,1) ) / det
    ainv(3,3) =(   a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det
    return
  end subroutine inversa
