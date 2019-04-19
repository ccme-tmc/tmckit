      subroutine skiplines(fid,n)
      implicit none 
      integer:: fid,n
      integer:: i

      do i=1,n
        read(fid,*) 
      enddo
      end subroutine

