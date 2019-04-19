!! This file contains subroutines used for least square fitting 
!!   LSQPNF -- least square polynomial fitting 
!!   LSQFIT -- general subroutines for least square fitting 

      SUBROUTINE LSQPNF(C,ERR,Y, X, Nx,NF)
!  Least SQuare PolyNomial Fitting
!        Y(x)= \sum_n  C_k x^n
!  Input:
!        X(1:Nx) --- x coordinates
!        Y(1:Nx) --- x coordinates
!        Nx   --- Number of data to be fitted
!        NF   --- order of fitting,
!  Output:
!        C(NF+1)   --- fitted coefficients
!        ERR    --- fitting error, defined as
!            SQRT[ \sum_i ( Y(x_i) - \sum_n  C_n f_n(x_i) )^2 ]
      IMPLICIT NONE
      INTEGER*4 NF,Nx
      REAL*8    C(0:NF),Y(Nx),X(Nx),ERR
      INTEGER*4 n
      REAL*8    F(Nx,0:NF)

      F(:,0)=1.0
      DO N=1,NF
        F(:,N)=F(:,N-1)*X
      ENDDO
      CALL LSQFIT(C,ERR,Y,F,Nx,NF+1)
      END SUBROUTINE

      SUBROUTINE LSQFIT(C,ERR,Y,F,Nx,NF)
!  A general subroutine for Least SQuare fitting
!        Y(x)= \sum_k  C_k f_k(x)
!  Input:
!        Nx   --- Number of data to be fitted
!        NF     --- Number of unknow parameters
!        F(Nx,NF) --- Values of basis functions
!            F(i,n)= f_n( x_i )
!  Output:
!        C(NF)   --- values of fitting parameters
!        ERR    --- fitting error, defined as
!            SQRT[ \sum_i ( Y(x_i) - \sum_n  C_n f_n(x_i) )^2 ]

      IMPLICIT NONE
      INTEGER*4 NF,Nx
      REAL*8    C(NF),Y(Nx),F(Nx,NF),ERR
      REAL*8    A(NF,NF),B(NF),WORK(NF*8)
      INTEGER*4 i,n,m,IPIV(NF),INFO

! calculate coefficient matrix and inhomogenity
      DO n=1,NF
         C(n)=SUM(Y*F(:,n))
         DO m=1,n
            A(m,n)=SUM(F(:,m)*F(:,n))
            A(n,m)=A(m,n)
         ENDDO
      ENDDO

! Call LAPACK routine to solve the linear equation
      CALL DSYSV('L',NF,1,A,NF,IPIV, C, NF,WORK,NF*8,INFO)       
      IF(INFO.NE.0) THEN
         WRITE(6,*) "LSQFIT: Error in calling DSYSV ---",INFO
         STOP
      ENDIF
! Calculate the fitting error
      ERR=0.0
      DO I=1,Nx
         ERR=ERR+ (Y(I)- SUM(C*F(I,:)))**2
      ENDDO

      ERR= SQRT(ERR/Nx)

      RETURN
      END SUBROUTINE

      subroutine linfit(x,y,n,a,b,ferr)
!!    Linear Fitting 
      implicit none
      integer,intent(in):: n
      real(8),intent(in):: x(n),y(n)
      real(8),intent(out):: a,b,ferr
      real(8):: mx,my,mxy,mx2
      mx=sum(x)/n
      my=sum(y)/n
      mx2=sum(x*x)/n
      mxy=sum(x*y)/n
      b=(mxy-mx*my)/(mx2-mx*mx)
      a=my-b*mx
      ferr=sqrt(sum((y-a-b*x)**2)/n)
      end subroutine

