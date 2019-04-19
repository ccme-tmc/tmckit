      subroutine randspin(spins,nat)
      implicit none
      integer,intent(in):: nat
      integer,intent(out):: spins(nat)
      real*8:: rand(nat)
      INTEGER              :: isize,idate(8)
      INTEGER,ALLOCATABLE  :: iseed(:)

      CALL DATE_AND_TIME(VALUES=idate)
      CALL RANDOM_SEED(SIZE=isize)
      ALLOCATE( iseed(isize) )
      CALL RANDOM_SEED(GET=iseed)
      iseed = iseed * (idate(8)-500)      ! idate(8) contains millisecond
      CALL RANDOM_SEED(PUT=iseed)

      CALL RANDOM_NUMBER(rand)

      spins(1:nat)=nint(sign(1.d0,rand(1:nat)-0.5d0))
      if( maxval(abs(spins(1:nat))).ne.1 .or. &
     &    minval(abs(spins(1:nat))).ne.1)  then
        write(6,*) "ERROR in randspin: vales other than 1 or -1 found"
      endif
      DEALLOCATE( iseed )
  
      end subroutine 
     
      

