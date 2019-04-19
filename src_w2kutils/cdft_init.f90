! This program initializes the constrained DFT calculations 
      program cdft_init 
      implicit none 
      integer:: np_orb,lorb,iat_imp,nat,nat_orb,icmplx,iop_occ,   &
     &          iop_neutral,iop_other,iats_orb(1:100)  
      real*8::  nel_orb,temp
      character(len=80)::casename

! occmax: maximal oupation for s=-1 at given kappa 
      integer   :: iat,norb,n,kappa,iorb
      integer   :: occmax,occ,jflag
      real*8    :: shift,newshift,neup,oup(0:3),odn(0:3),nedn,md

!! variables used for *.in1      
      integer   :: l,ndiff,napw,napwl
      real*8    :: etrial,el,de 
      character(len=4)::switch

!! variables used for *.in2
      real*8    :: emin,neval,esepermin,eseper0

      character :: cmplxflag
      integer   :: uinc,uincup_p,uincdn_p,uincup_m,uincdn_m,            &
     &             uin1,uin2,uin1_c,uin2_p,uin2_m,                      &
     &             uin2_j,uincup_j,uincdn_j,uerr 
      character(len=80)::finc,fincup_p,fincdn_p,fincup_m,fincdn_m,      &
     &             fin1,fin2,fin1_c,fin2_p,fin2_m,fin2_j,fincup_j,fincdn_j 
      character(len=80)::str,fname="cdft_init"

!! readin all necessary input parameters
!!   Line 1  casename 
!!        2  iop_occ (0/1/2), iop_netural(0/1), iop_other (0/1/2), icmplx (0/1)
!!        3  np_orb,lorb,nel_orb 
!!        4  nat,nat_orb,iats_orb(1:nat_orb)
!!        5  newshift, temp  

      uerr=999
      open(unit=uerr,file=trim(fname)//".error",action='write') 
      write(uerr,*) "Error in "//trim(fname)

      write(6,*) "Case name?"
      read *,casename 
      write(6,'(a5,a)') " =>",trim(casename)

!! iop_occ: control occupation configurations
!!    0 -- the original scheme used in Ansimov-Gunnarsson (PRB 43, 7570 (1991))
!!    1 -- low  spin configuration 
!!    2 -- high spin configuration 

!! iop_neutral : Option for charge neutralization
!!    0 -- put additional charge to valence 
!!    1 -- neutralized by background charge 

!! iop_other: option to treat localized electrons on the other target atoms
!!    0 -- kept in valence 
!!    1 -- put into core, half up, half dn
!!    2 -- put into core, all up
!! icmplx: whether doing complex calculations 
!!    0 -- real 
!!    1 -- complex

!! newshift : shift used in inc for open-core calculations  
!! temp: smearing temperature used in in2, neglected if negative 

      write(6,*) "iop_occ, iop_neutral, iop_other, icmplx?"
      read *,iop_occ, iop_neutral, iop_other, icmplx
      write(6,'(a5,4i5)') " =>",iop_occ,iop_neutral, iop_other, icmplx

      if(icmplx.eq.0) then 
        cmplxflag=' '
      else 
        cmplxflag='c'
      endif 

      write(6,*) "np_orb,lorb,nel_orb?"
      read *,np_orb,lorb,nel_orb
      write(6,'(a5,2i5,f6.2)') " =>",np_orb,lorb,nel_orb
   
      write(6,*) "nat,nat_orb, and iats_orb? "
      read *,nat,nat_orb,iats_orb(1:nat_orb)
      write(6,'(a5,100i5)') " =>",nat,nat_orb,iats_orb(1:nat_orb)
      iat_imp=iats_orb(1)

      write(6,*) "newshift (inc), smearing temperature(in2) ?"
      read *,newshift,temp 
      write(6,'(a5,2f10.4)') " =>",newshift,temp

!! maximal occupation number of the shell (excluding the spin degeneracy) 
      occmax=2*lorb+1 
      write(6,*) "Maximal occupation (excluding spin):",occmax
      if(nel_orb.ge.2*occmax) then 
        write(6,*)"WARNING: nel_orb is bigger than the physical max.!"
      endif 

!! Setup the oupation configuration for n+0.5 and n-0.5


!! for the case of nel_orb.eq.1, iop_occ=0 or 1 can give negative or completely zero occupation numbers
!!
      if(nel_orb.le.1.d0 .and.iop_occ.eq.0 ) then 
         write(6,*) "ERROR: nel_orb <=1 is not doable for iop_occ=0 " 
         stop 
       endif 

!! Scheme - 0

      if(iop_occ.eq.0) then 
        write(6,*) "iop_occ=0: Anisimov-Gunnarsson CDFT scheme"
        oup(0)=nel_orb+1
        odn(0)=nel_orb-1
        oup(1)=nel_orb+1 
        odn(1)=nel_orb
        oup(2)=nel_orb+1 
        odn(2)=nel_orb-2
        jflag=0
        write(6,'(a,2f8.3)') "neup,nedn=",nel_orb/2,nel_orb/2
      else
        write(6,*) "iop_occ!=0: New CDFT scheme"
        if(iop_occ.lt.0) then 
          write(6,*) "iop_occ<0: customized configuration"
          neup = abs(iop_occ)
          nedn = nel_orb-neup
        elseif(iop_occ.eq.1) then  
          write(6,*) "iop_occ=1: Low Spin configuration"
          neup = nel_orb*0.5d0
          nedn = nel_orb*0.5d0
        else 
          write(6,*) "iop_occ=2: High Spin configuration"
          if(nel_orb.le.occmax) then
            neup = nel_orb
            nedn = 0.d0
          else
            neup = occmax
            nedn = nel_orb-occmax
          endif
        endif 
        write(6,'(a,2f8.3)') "neup,nedn=",neup,nedn

        oup(1)=neup*2.d0 + 1.d0
        odn(1)=nedn*2.d0

        oup(2)=neup*2.d0 - 1.d0
        odn(2)=nedn*2.d0

        if(min(minval(oup(1:2)),minval(odn(1:2))).lt.0.d0) then 
          write(6,*) "ERROR: negative occupation occurs"
          stop
        endif 
         
        !! determine the occ configuration for the calculation of J
        if(neup.ge.1.0.and.abs(neup-nedn-1.0).gt.1.e-5) then 
          oup(0)=neup*2.d0-2.d0
          odn(0)=nedn*2.d0+1.d0
          jflag=-1
        else
          oup(0)=neup*2.d0+2.d0
          odn(0)=nedn*2.d0-1.d0
          jflag=1
        endif 
        write(6,'(a,2f6.2)') " neup,nedn=",neup,nedn
        if(min(oup(0),odn(0)).lt.0.0) then 
          write(6,*) "***"
          write(6,*) "WARNING: the configuration you defined is not &
     &appropirate for determining J"
          write(6,*) "***"
        endif  
      endif


! set occupation for local states on other atoms

      if(iop_other.eq.0) then 
        oup(3)=0.d0
        odn(3)=0.d0 
      elseif(iop_other.eq.1) then 
        oup(3)=nel_orb
        odn(3)=nel_orb
      else 
        oup(3)=nel_orb*2
        odn(3)=0.d0
      endif 

      write(6,'(a,2f6.2)') "Occ_up/dn for np (+0.5) :",oup(1)/2,odn(1)/2
      write(6,'(a,2f6.2)') "Occ_up/dn for nm (-0.5) :",oup(2)/2,odn(2)/2
      write(6,'(a,2f6.2)') "Occ_up/dn for nj        :",oup(0)/2,odn(0)/2
      write(6,'(a,2f6.2)') "Occ_up/dn for others    :",oup(3)/2,odn(3)/2


!! Set I/O files  
      uinc=10
      uincup_p=11
      uincdn_p=12
      uincup_m=13
      uincdn_m=14

      finc=trim(casename)//'.inc'
      fincup_p=trim(casename)//'.incup_cdft_p' 
      fincdn_p=trim(casename)//'.incdn_cdft_p' 
      fincup_m=trim(casename)//'.incup_cdft_m' 
      fincdn_m=trim(casename)//'.incdn_cdft_m' 
      open(unit=uinc,file=trim(finc),action='read') 
      open(unit=uincup_p,file=trim(fincup_p),action='write') 
      open(unit=uincdn_p,file=trim(fincdn_p),action='write') 
      open(unit=uincup_m,file=trim(fincup_m),action='write') 
      open(unit=uincdn_m,file=trim(fincdn_m),action='write') 

      uin1=15
      uin2=16
      uin1_c=17
      uin2_p=18
      uin2_m=19

      fin1=trim(casename)//'.in1'//trim(cmplxflag)
      fin2=trim(casename)//'.in2'//trim(cmplxflag)
      fin1_c=trim(casename)//'.in1'//trim(cmplxflag)//'_cdft'
      fin2_p=trim(casename)//'.in2'//trim(cmplxflag)//'_cdft_p'
      fin2_m=trim(casename)//'.in2'//trim(cmplxflag)//'_cdft_m'

      open(unit=uin1,file=trim(fin1),action='read')
      open(unit=uin2,file=trim(fin2),action='read')
      open(unit=uin1_c,file=trim(fin1_c),action='write')
      open(unit=uin2_p,file=trim(fin2_p),action='write')
      open(unit=uin2_m,file=trim(fin2_m),action='write')

      uin2_j=23
      uincup_j=21
      uincdn_j=22 
      fin2_j=trim(casename)//'.in2'//trim(cmplxflag)//'_cdft_j'
      fincup_j=trim(casename)//'.incup_cdft_j'
      fincdn_j=trim(casename)//'.incdn_cdft_j'
      open(unit=uin2_j,file=trim(fin2_j),action='write')
      open(unit=uincup_j,file=trim(fincup_j),action='write')
      open(unit=uincdn_j,file=trim(fincdn_j),action='write')
      
!
! generate inc for n+0.5 and n-0.5
!
      do iat=1,nat
        read(uinc,*) norb,shift 
        if((iat.eq.iat_imp) .or. (iop_other.gt.0 .and.                  &
     &      l_belong(iat,iats_orb,nat_orb)) ) then 
          write(uincup_p,100) norb+2,newshift
          write(uincdn_p,100) norb+2,newshift
          write(uincup_m,100) norb+2,newshift
          write(uincdn_m,100) norb+2,newshift
          write(uincup_j,100) norb+2,newshift
          write(uincdn_j,100) norb+2,newshift
        else 
          write(uincup_p,100) norb,shift
          write(uincdn_p,100) norb,shift
          write(uincup_m,100) norb,shift
          write(uincdn_m,100) norb,shift
          write(uincup_j,100) norb,shift
          write(uincdn_j,100) norb,shift
        endif 

        do iorb=1,norb 
          read(uinc,*) n,kappa,occ
          write(uincup_p,101) n,kappa,occ
          write(uincdn_p,101) n,kappa,occ
          write(uincup_m,101) n,kappa,occ
          write(uincdn_m,101) n,kappa,occ
          write(uincup_j,101) n,kappa,occ
          write(uincdn_j,101) n,kappa,occ
        enddo 

        !! add additional core electrons 
        if(iat.eq.iat_imp) then 
          write(uincup_p,102) np_orb,   lorb, occlsj1(oup(1),lorb,1) 
          write(uincup_p,102) np_orb,-lorb-1, occlsj1(oup(1),lorb,2) 

          write(uincdn_p,102) np_orb,   lorb, occlsj1(odn(1),lorb,1) 
          write(uincdn_p,102) np_orb,-lorb-1, occlsj1(odn(1),lorb,2) 

          write(uincup_m,102) np_orb,   lorb, occlsj1(oup(2),lorb,1) 
          write(uincup_m,102) np_orb,-lorb-1, occlsj1(oup(2),lorb,2) 

          write(uincdn_m,102) np_orb,   lorb, occlsj1(odn(2),lorb,1) 
          write(uincdn_m,102) np_orb,-lorb-1, occlsj1(odn(2),lorb,2) 

          write(uincup_j,102) np_orb,   lorb, occlsj1(oup(0),lorb,1) 
          write(uincup_j,102) np_orb,-lorb-1, occlsj1(oup(0),lorb,2) 

          write(uincdn_j,102) np_orb,   lorb, occlsj1(odn(0),lorb,1) 
          write(uincdn_j,102) np_orb,-lorb-1, occlsj1(odn(0),lorb,2) 

        elseif(iop_other.gt.0 .and. l_belong(iat,iats_orb,nat_orb) ) then 
          write(uincup_p,102) np_orb,   lorb, occlsj1(oup(3),lorb,1)
          write(uincup_p,102) np_orb,-lorb-1, occlsj1(oup(3),lorb,2)

          write(uincdn_p,102) np_orb,   lorb, occlsj1(odn(3),lorb,1)
          write(uincdn_p,102) np_orb,-lorb-1, occlsj1(odn(3),lorb,2)

          write(uincup_m,102) np_orb,   lorb, occlsj1(oup(3),lorb,1)
          write(uincup_m,102) np_orb,-lorb-1, occlsj1(oup(3),lorb,2)

          write(uincdn_m,102) np_orb,   lorb, occlsj1(odn(3),lorb,1)
          write(uincdn_m,102) np_orb,-lorb-1, occlsj1(odn(3),lorb,2)

          write(uincup_j,102) np_orb,   lorb, occlsj1(oup(3),lorb,1)
          write(uincup_j,102) np_orb,-lorb-1, occlsj1(oup(3),lorb,2)

          write(uincdn_j,102) np_orb,   lorb, occlsj1(odn(3),lorb,1)
          write(uincdn_j,102) np_orb,-lorb-1, occlsj1(odn(3),lorb,2)
        endif 
      enddo 

!
!! generate new *.in1 with frozen impurity d/f states 
!
      !! the first two lines are simply duplicated
      read(uin1,1000) str 
      write(uin1_c,1001) trim(str)
      read(uin1,1000) str 
      write(uin1_c,1001) trim(str)

      do iat=1,nat 
        read(uin1,*) etrial,ndiff,napw 
        write(uin1_c,201) etrial,ndiff,napw 
        do iorb=1,ndiff 
          read(uin1,202) l,el,de,switch,napwl
          if( l.eq.lorb .and. ( iat.eq.iat_imp .or.(iop_other.gt.0      &
     &       .and. l_belong(iat,iats_orb,nat_orb)) ) )  then 
            write(uin1_c,203) l,el+20.d0,0.d0,'CONT',napwl
          else 
            write(uin1_c,203) l,el,de,switch,napwl
          endif 
        enddo 
      enddo 
      read(uin1,1000) str
      write(uin1_c,1001) trim(str)

!
!! generate new *.in2 file 
!
     
      read(uin2,1000) str
      write(uin2_p,1000) str
      write(uin2_m,1000) str
      write(uin2_j,1000) str

      read(uin2,*) emin,neval,esepermin,eseper0 
      if(iop_other.gt.0) then 
        neval = neval - nel_orb*nat_orb 
      else 
        neval = neval - nel_orb 
      endif 
       
      if(iop_neutral.eq.0) then 
        write(uin2_p,301) emin,neval - 0.5d0,esepermin,eseper0
        write(uin2_m,301) emin,neval + 0.5d0,esepermin,eseper0
        if(iop_occ.eq.0) then 
          write(uin2_j,301) emin,neval,esepermin,eseper0
        else 
          write(uin2_j,301) emin,neval+0.5d0,esepermin,eseper0
        endif 
      else 
        write(uin2_p,301) emin,neval,esepermin,eseper0
        write(uin2_m,301) emin,neval,esepermin,eseper0
        write(uin2_j,301) emin,neval,esepermin,eseper0
      endif 

      read(uin2,1000) str
      if(temp.gt.0.d0) then 
        write(uin2_p,302) 'TEMP ',0.002
        write(uin2_m,302) 'TEMP ',0.002
        write(uin2_j,302) 'TEMP ',0.002
      else 
        write(uin2_p,1000) str
        write(uin2_m,1000) str
        write(uin2_j,1000) str
      endif

      !! duplicate all remaining lines 
       
  800 read(uin2,1000,END=900) str 
      write(uin2_p,1000) str
      write(uin2_m,1000) str
      write(uin2_j,1000) str
      goto 800
  900 CONTINUE

!! generate *.cdft_conf file 
      open(8,file=trim(casename)//".cdft_conf",action='write')
      write(8,401) iop_occ,iat_imp,np_orb,lorb,nel_orb   
      write(8,402) occlsj2(oup(1),lorb),occlsj2(odn(1),lorb)  
      write(8,403) occlsj2(oup(2),lorb),occlsj2(odn(2),lorb)  
      write(8,404) occlsj2(oup(0),lorb),occlsj2(odn(0),lorb),jflag   
      close(8) 

      close(uinc)
      close(uincup_p)
      close(uincdn_p)
      close(uincup_m)
      close(uincdn_m)
      close(uin2_j) 
      close(uincup_j)
      close(uincdn_j)

      close(uin1)
      close(uin2)
      close(uin1_c)
      close(uin2_p)
      close(uin2_m)
      close(uerr,status='delete')  

 100  format(i2,f5.2,5x,'NUMBER OF ORBITALS (EXCLUDING SPIN), SHIFT')
 101  format(i1,',',i2,',',i1,15x,'( N,KAPPA,OCCUP)')
 102  format(i1,',',i2,',',f7.5,15x,'( N,KAPPA,OCCUP)')
 1000 format(a80)
 1001 format(a)

 201  format(f7.5,i5,i3,6x,'(GLOBAL E-PARAMETER WITH n OTHER CHOICES, ' &
           ,'global APW/LAPW)')
 202  format(i2,2f10.5,a4,i2)
 203  format(i2,f8.2,f11.3,1x,a4,i2)

 301  format(2f10.1,2f5.2,16x,'EMIN, NE, ESEPERMIN, ESEPER0')
 302  format(a5,f10.5,10x,'(GAUSS,ROOT,TEMP,TETRA,ALL      eval)')
      

 401  format(4i5,f6.2,6x,"# iop_occ,iat_imp,n,l,nel_orb")
 402  format(4f8.4,8x,"# Occ of added core states for np")
 403  format(4f8.4,8x,"# Occ of added core states for nm")
 404  format(4f8.4,i4,4x,"# Occ of added core states for nj and jflag")

      contains 

      function occlsj1(ne,l,is) result(occ)
      real(8)::ne
      integer::l,is
      real*8:: occ
      if(is.eq.1) then
        occ=max(min(ne,2.d0*l),0.0001d0)
      else
        occ=max(0.d0,ne-2.d0*l)
      endif
      end function 

      function occlsj2(ne,l)  result(occ)
      real(8)::ne
      integer::l
      real*8:: occ(2)
      occ(1)=max(min(ne,2.d0*l),0.0001d0) 
      occ(2)=max(0.d0,ne-2.d0*l)
      end function 

      function l_belong (ia,ias,na) result(lb)
      integer:: na,ia,ias(na)
      logical:: lb
      integer:: i
      
      lb=.false.
      do i=1,na
        if(ias(i).eq.ia) then 
          lb=.true.
          return 
        endif 
      enddo
      end function 
    
      end program           
        
             
