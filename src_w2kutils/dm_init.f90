      program dm_init 
      implicit none 
      integer,parameter:: labc=4 
      integer::iop_sp,isp,l,i,j,il,ia,ierr,natorb,itmp,iop_dm,nel
      integer::fid(2)
      real(8) :: nl
      real(8)   :: dm_input(-3:3)
      complex(8):: rdm(-3:-3,-3:-3) 
      real(8),allocatable:: nls(:,:,:)
      complex(8),allocatable::dmat(:,:,:,:,:) 
      integer,allocatable:: nlorb(:),iatom(:),lorb(:,:)
      character(len=80):: casename,fname
      character(len=2) ::spflag(2)

      write(6,*) "====================================================="
      write(6,*) "*                     dm_init                       *"
      write(6,*) "====================================================="

      write(6,*) 'Case name?'
      read *,casename
      write(6,*) "<< ", trim(casename)
      write(6,*) 'iop_sp  ( 0 - runafm, 1 - runsp_c, 2 - runsp )?'
      read *,iop_sp 
      write(6,*) "<< ", iop_sp
      write(6,*) "iop_dm"
      read *,iop_dm
      write(6,*) "<< ",iop_dm
!
! iop_dm determines how to initialize density matrix  
!    iop_dm == 0  -- random dmat 
!              1  -- random diagonal dmat with trace read from input 
!              2  -- diagonal dmat by reading the diagonal elements from input
!              3  -- diagonal dmat by reading diagonal elements from and input and add some small random purturbation 
!
      open(10,file=trim(casename)//'.inorb',iostat=ierr)
      if(ierr.ne.0) then 
        write(6,*) "ERROR in dm_init: fail to open ",trim(casename)&
     &//'.inorb' 
        stop
      endif 
      read(10,*) itmp,natorb,itmp
      read(10,*) 

      allocate (iatom(natorb),nlorb(natorb),lorb(labc,natorb), & 
     &          dmat(-3:3,-3:3,labc,natorb,2),nls(labc,natorb,2))
      dmat=0.d0
      
      do i=1,natorb
        read(10,*) iatom(i),nlorb(i),(lorb(il,i),il=1,nlorb(i))
        write(6,202) iatom(i),(lorb(il,i),il=1,nlorb(i))
        if(nlorb(i).gt.labc) stop 'nlorb gt labc'
      enddo

      call init_random_seed
!
! generate nl's randomly 
!
      spflag(1)='up'
      spflag(2)='dn'
      fid(1)=11
      fid(2)=12

!! set dmat 
      do isp=1,2
        fname=trim(casename)//'.dmat'//spflag(isp)
        open(fid(isp),file=fname,iostat=ierr)
        if(ierr.ne.0) then
          write(6,*) "ERROR in dm_init: fail to open ",trim(fname)
          stop
        endif

        do ia=1,natorb 
          do il=1,nlorb(ia)

            l=lorb(il,ia) 
            if(isp.eq.1) then 
              call set_dmat 
            endif

            if(isp.eq.2) then 
              if(iop_sp.eq.0) then   !! runafm  
                dmat(:,:,il,ia,2)=dmat(:,:,il,natorb-ia+1,1) 
                nls(il,ia,2)=nls(il,natorb-ia+1,1)
              elseif(iop_sp.eq.1) then !! runsp_c
                dmat(:,:,il,ia,2)=dmat(:,:,il,ia,1)
                nls(il,ia,2)=nls(il,ia,1)
              else        !! runsp
                call set_dmat
              endif
            endif  

          enddo  ! il
        enddo  ! iat
      enddo ! isp

!! write dmat file 
      do isp=1,2
        rewind(fid(isp))
        do ia=1,natorb
          do il=1,nlorb(ia)
            l=lorb(il,ia) 
            write(fid(isp),100) iatom(ia)
            write(fid(isp),102) l, 0.d0,0.d0, 0.d0
            do i=-l,l
              write(fid(isp),101) (dmat(i,j,il,ia,isp),j=-l,l)
            enddo
          enddo
        enddo
        close(fid(isp))
      enddo

!! print dmat to stdout
      do ia=1,natorb
        do il=1,nlorb(ia)
          l=lorb(il,ia)
          write(6,*) 
          write(6,300) ia,l,nls(il,ia,1:2)
          do isp=1,2
            write(6,*) 
            write(6,*) "spin ",spflag(isp)
            write(6,301) (j,j=-l,l)
            write(6,*) 
            do i=-l,l
              write(6,302) i,real(dmat(i,-l:l,il,ia,isp))
            enddo
          enddo 
        enddo
      enddo 
      close(10)
      deallocate(iatom,nlorb,lorb,dmat,nls)

100   format(i5,' atom density matrix')
101   format(2(2e16.8,2x))
102   format(i5,3f10.6, ' L, Lx,Ly,Lz in global orthogonal system')
300   format("# dmat for ia=",i3," l=",i1," Trace up:",f5.2," dn:",f5.2)
301   format(7x,7i5)
302   format(i5,2x,7f5.2)
202   format(' Vorb applied to atom',i4,' orbit. numbers',3i4)

      contains 
        subroutine set_dmat()

        if(iop_dm.eq.0) then !! random Hermertian densty matrix
          call random_number(nl)
          nl=nl*(2.d0*l+1.d0)
          call dmat_rand(dmat(-l:l,-l:l,il,ia,isp),l,nl)

        elseif(iop_dm.eq.1) then  !! random diagonal density matrix
          read *,nl
          call dmat_rand(rdm(-l:l,-l:l),l,nl)
          do i=-l,l
            dmat(i,i,il,ia,isp)=rdm(i,i)
          enddo
        elseif(iop_dm.eq.2.or.iop_dm.eq.3) then   !! read the diagonal dm from input and set off-diagonal as zero
          read *,dm_input(-l:l)
          rdm=0.d0
          if(iop_dm.eq.3) then                  !! read diagonal dm from input and then add a small random matrix 
            call dmat_rand(rdm(-l:l,-l:l),l,0.01)
          endif   
          nl=0.d0
          do i=-l,l
            dmat(i,i,il,ia,isp)=dmat(i,i,il,ia,isp)+dm_input(i)+rdm(i,i)
            nl=nl+dmat(i,i,il,ia,isp)
          enddo
        endif
        nls(il,ia,isp)=nl
        end subroutine 
      end program 

      subroutine dmat_rand(dmat,l,nl)  
!! generate an randon Hermitian matrix with given trace nl
      implicit none 
      integer,intent(in) :: l
      real(8),intent(in) :: nl 
      complex(8),intent(out) :: dmat(-l:l,-l:l)

      integer::i,j 
      real(8)::tr
      real(8),allocatable::re(:,:),im(:,:)

      allocate(re(-l:l,-l:l),im(-l:l,-l:l))
      call random_number(re) 
      call random_number(im) 
      dmat=cmplx(re,im,8)

      dmat=dmat+transpose(conjg(dmat))
      tr=0.d0
      do i=-l,l
        tr=tr+real(dmat(i,i))
      enddo 
      dmat=dmat*nl/tr
      deallocate(re,im)
      endsubroutine  

      subroutine dmat_hundrule(dmat,l,nl,nel,isp) 
      integer,intent(in)::nel,isp
      real(8),intent(out)::nl
      complex(8),intent(out)::dmat(-l:l,-l:l)
      integer::i

      dmat=0.d0 
      nl=0.d0
      if(isp.eq.1) then 
        if(nel.le.2*l+1) then 
          do i=1,nel
            dmat(-l+i-1,-l+i-1)=1.d0
          enddo 
          nl=1.d0*nel
        else 
          do i=-l,l
            dmat(i,i)=1.d0
          enddo
          nl=2.d0*l+1
        endif 
      else
        if(nel.gt.2*l+1) then 
          do i=1,nel-(2*l+1),-1
            dmat(l-i+1,l-i+1)=1.d0
          enddo
          nl=dble(nel-2*l-1)
        endif 
      endif 

      end subroutine 

      SUBROUTINE init_random_seed()
      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
      END SUBROUTINE

