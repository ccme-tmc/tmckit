 module struct

   CHARACTER*80             :: title
   CHARACTER*3              :: lattyp,irel
   REAL*8                   :: aa,bb,cc,pia(3),alpha(3)
   REAL*8                   :: vol
   
   INTEGER                  :: nat,ndif
   INTEGER,ALLOCATABLE      :: mult(:),jrj(:),iatnr(:),isplit(:)
   REAL*8,ALLOCATABLE       :: r0(:),dx(:),rmt(:),zz(:)
   CHARACTER*10,ALLOCATABLE :: aname(:)
   REAL*8,allocatable       :: pos(:,:)
   REAL*8, allocatable      :: rotloc(:,:,:)
   integer iord
   integer sgnum
   character*8 sgname
   real*8      tobohr
   logical     HRtransf
   real*8      hex2ort(3,3),ort2rho(3,3),hex2rho(3,3),rho2hex(3,3)

contains

 subroutine init_struct(nato)

   tobohr=1.0d0/0.5291772083d0

    allocate (mult(nato),jrj(nato),iatnr(nato),isplit(nato))  
    allocate (r0(nato),dx(nato),rmt(nato),zz(nato))
    allocate (aname(nato),pos(3,50*nato),rotloc(3,3,nato))  
    do i=1,nato
       iatnr(i)=-i
    enddo
   isplit(1:nato)=15 
   jrj(1:nato)=781
   r0(1:nato)=5.0d-6
   rmt(1:nato)=2.0d0
   rotloc(1:3,1:3,1:nato)=0.0d0
   rotloc(1,1,1:nato)=1.0d0
   rotloc(2,2,1:nato)=1.0d0
   rotloc(3,3,1:nato)=1.0d0
   iord=0
   title='blebleble'

   hex2ort(1,1)=0.0d0
   hex2ort(1,2)=1.0d0
   hex2ort(1,3)=0.0d0
   hex2ort(2,1)=sqrt(3.0d0)/2.0d0
   hex2ort(2,2)=-0.5d0
   hex2ort(2,3)=0.0d0
   hex2ort(3,1)=0.0d0
   hex2ort(3,2)=0.0d0
   hex2ort(3,3)=1.0d0
   ort2rho(1,1)=1.0d0/sqrt(3.0d0)
   ort2rho(1,2)=1.0d0/sqrt(3.0d0)
   ort2rho(1,3)=-2.0d0/sqrt(3.0d0)
   ort2rho(2,1)=-1.0d0
   ort2rho(2,2)=1.0d0
   ort2rho(2,3)=0.0d0
   ort2rho(3,1)=1.0d0
   ort2rho(3,2)=1.0d0
   ort2rho(3,3)=1.0d0
   hex2rho=matmul(hex2ort,ort2rho)
   call inversa(hex2rho,rho2hex)

   HRtransf=.false.

 end subroutine init_struct

 subroutine write_struct
  
   implicit none
   
   integer                :: index,i,j,j1,j2,m,jatom
   
   aa=aa*tobohr
   bb=bb*tobohr
   cc=cc*tobohr

   if (HRtransf) then
      i=0
      do jatom=1,nat
         do m=1,mult(jatom)                                     
            i=i+1                                            
            pos(1:3,i)=matmul(pos(1:3,i),hex2rho)
            do j=1,3
               if (pos(j,i).lt.0.0d0) pos(j,i)=pos(j,i)+1.0d0
               if (pos(j,i).gt.1.0d0) pos(j,i)=pos(j,i)-1.0d0
            enddo
         enddo
      enddo
   endif

   write(20,1000) title
   write(20,1010) lattyp,nat,sgnum,sgname
   write(20,1020) aa,bb,cc,alpha(1),alpha(2),alpha(3)
   index=0
   do jatom=1,nat
      index=index+1
      write(20,1030) iatnr(jatom),( pos(j,index),j=1,3 ), &
           mult(jatom),isplit(jatom) 
      do m=1,mult(jatom)-1                                     
         index=index+1                                            
         write(20,1031) iatnr(jatom),( pos(j,index),j=1,3)         
      enddo
      write(20,1050) aname(jatom),jrj(jatom),r0(jatom),rmt(jatom), &
           zz(jatom)
      write(20,1051) ((rotloc(i,j,jatom),i=1,3),j=1,3)                
   enddo
   write(20,1151) iord
     
1000 format(a80)                                                       
1010 format(a3,1x,'LATTICE,NONEQUIV.ATOMS ',i3,2x,i3,1x,a8/,&
            'MODE OF CALC=RELA unit=bohr')                                 
1020 format(6f10.6)                                          
1030 FORMAT('ATOM',I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8,/,&
             '          MULT=',i2,'          ISPLIT=',i2)          
1031 FORMAT(4X,I4,': X=',F10.8,' Y=',F10.8,' Z=',F10.8)                   
1050 format(a10,' NPT=',i5,'  R0=',f10.9,' RMT=',f10.5,'   Z:',f10.5)
1051 format('LOCAL ROT MATRIX:   ',3f10.7,/,20x,3f10.7,/,20x,3f10.7)
1151 format(i4,'      NUMBER OF SYMMETRY OPERATIONS')

   end subroutine write_struct
   
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

 end module struct
