  subroutine HerEigen(A,D,N)
    implicit none
    integer:: N
    complex(8):: A(N,N)
    real(8) :: D(N)
    complex(8) :: Work(2*N)
    real(8) :: RWORK(3*N-2)
    integer::info,i,j

    call ZHEEV( 'V', 'L', N, A, N, D, WORK, 2*N, RWORK, INFO )
    if(info/=0 ) then 
      write(6,*) "ERROR HerEigen,info=",info 
      stop
    endif 
  end subroutine 

  subroutine SymEigen(a,d,n)
    implicit none
    integer::N
    real(8) ::a(N,N),d(N)
    real(8) ::Work(N*6)
    integer::info,lwork
    character::JOBZ='V',UPLO='U'

    lwork=6*N
    call DSYEV( JOBZ, UPLO, N, A, N, d, WORK, LWORK, INFO )
    if(info/=0 ) then
      write(6,*) "ERROR HerEigen,info=",info
      stop
    endif
  end subroutine 

  subroutine LinEq(a,b,n)
    implicit none
    integer::n
    real(8)::a(n,n),b(n)
    integer::IPIV(N),info
    real(8)::work(N*2)

    call DSYSV('u', N, 1, A, N, IPIV, B,N, WORK,N, INFO )
    if(info/=0 ) then
      write(6,*) "ERROR HerEigen,info=",info
      stop
    endif

  end subroutine 

  subroutine InvSym(A,NM,N)
    implicit none
    integer::N,NM
    real(8)::a(NM,N)
    integer::IPIV(N),info,i,j
    real(8)::work(N)

    call DSYTRF( 'U', N,  A, NM, IPIV, WORK, N, INFO )
    if(info/=0) then 
      write(6,*) "ERROR InvSym,info=",info
      stop
    endif 

    call DSYTRI( 'U', N,  A, NM, IPIV, WORK, INFO )
    if(info/=0) then 
      write(6,*) "ERROR InvSym,info=",info
      stop
    endif 
    do i=1,N
      do j=1,i-1
        A(i,j)=A(j,i)
      enddo
    enddo

  end subroutine 

  subroutine InvSVD(A,NM,N,tolsvd)
    implicit none
    integer::N,NM
    real(8)::a(NM,N),tolsvd

    INTEGER LWORK,IWORK(8*N),INFO,I,J,k
    REAL(8),allocatable::U(:,:),VT(:,:),S(:),WORK(:),IKERN(:,:)

    LWORK=8*N*N
    ALLOCATE( U(N,N), VT(N,N), S(N),WORK(LWORK),STAT=INFO)
    IF(INFO.NE.0) THEN
       PRINT *,"InvSymSVD: Fail to allocate arrays"
       STOP
    ENDIF

    CALL DGESVD('A','A',N, N,A,NM,S,U,N,VT,N,WORK,LWORK,INFO)
    IF(INFO.NE.0) THEN
       WRITE(6,*) "VXCOPM: Error in calling DGESDD ---",INFO
       STOP
    ENDIF
    
!    write(6,*) "Diagonal Values:"
    do i=1,N
!      write(6,'(e12.4)') s(i)
      if(s(i).lt.tolsvd) then 
         s(i)=0.d0
      else 
         s(i)=1.d0/s(i)
      endif 
    enddo 
    do j=1,N
      do i=1,N
        A(i,j)=0.d0
        do k=1,N
          A(i,j)=A(i,j)+u(i,k)*s(k)*vt(k,j)    
        enddo
      enddo 
    enddo
    deallocate(u,vt,s,work)
  end subroutine 
