module DirectMadMod

use structMod
use calcmad_func

PUBLIC nD_GetDirectMadelungPotential


CONTAINS

!> Calculate Madelung constant at site of atom OR non-atom position in 2D or 3D cases
real*8 function nD_GetDirectMadelungPotential(task,atom,pos)
    implicit none
    !Declare
    type(taskT),pointer,intent(in) :: task
    integer,intent(in),optional :: atom !< The index of atom. If set to 0 , then use fSite to calculate . Can be neglected
    real*8,intent(in),optional :: pos(3) !< The position of non-atom site. Only used when nIndex = 0. This argument can be neglected if not used
    !local var
    
    
    !local var
    real*8 dLast !< Record last result in converge iteration
    real*8 fR(3)
    integer nRmax,i,j,k,a2 ! vector*3, atom index(for 2nd atom)
    integer nVmax(3)  ! The maximum value of a,b,c-axis ( may change as dimension )
    real*8 dRabs
    real*8 dPotR !< The total potential of G/R part
    real*8 :: dSite(3) ! Store the position used
    integer :: nIndex ! Store the atom index. 0

    type(structT),pointer :: struct

    struct => task%struct

    !Convert input paramter to actual used parameter

    nIndex = 0 
    if ( present(atom)) then
        nIndex = atom
    end if

    if ( nIndex .eq. 0) then
        if ( present(pos)) then
            dSite(1:3) = pos(1:3)
        else
           stop "Must specify one of atom index or non-atom position coordinate when use nD_GetEwaldMadelungPotential" 
        end if
    else
        dSite(1:3) = struct%fCoord(1:3,nIndex)
    end if
    

    !
    !if ( struct%nDimension .ne. 2) then
    !    stop "nD_GetAtomEwaldMadelungPotential must be used only in 2D cases"
    !end if

    if ( nLogLevel >= 2) then
        write(6,*) "----------------------------------------"
    end if

    !Add Real Space
    dPotR=-10
    do nRmax=2,nSpaceMax,3
        dLast=dPotR
        dPotR=0
        nVmax = 0
        do i=1,struct%nDimension
            nVmax(i) = nRmax
        end do
        do i=-nVmax(1),nVmax(1)
            do j=-nVmax(2),nVmax(2)
                do k=-nVmax(3),nVmax(3)
                    fR=i*struct%fCellSize(1:3,1)+j*struct%fCellSize(1:3,2)+k*struct%fCellSize(1:3,3)
                        do a2= 1,struct%nAtom
                            if ( ( nIndex .eq. a2 ) .AND. ( i .eq. 0) .AND.( j .eq. 0) .AND. (k .eq. 0) )  then ! not for same atom in same lattice
                                cycle
                            end if
                            dRabs=VectorToAbs(dSite-struct%fCoord(1:3,a2)+fR)
                            !write(6,*) fTemp,struct%fCharge(a2)*DERFC(task%fAlpha*fTemp)/fTemp,dPotR,DERFC(task%fAlpha*fTemp)
                            dPotR=dPotR+struct%fCharge(a2)/dRabs
                            !fResultD=fResultD+struct%fCharge(a2)/fTemp
                         end do
                     end do
             end do
        end do
        if ( nRmax > nSpaceMax ) then
            write(6,*) "Rmax overflow! R-space does not converge."
        end if
            

        !write(6,'("Rmax: ",i4,"  ",f14.7," Direct: ",f14.7)') nRmax,dPotR ,fResultD-dSelf
        if ( nLogLevel >= 3) then
        write(6,'(" Rmax: ",i4,"  ",f14.7)') nRmax,dPotR
        end if
        if ( abs((dPotR-dLast)/dPotR) .lt. task%fConvLimit .OR. (abs(dPotR) .lt. task%fConvLimit) ) then
            exit
       end if             
    end do
    if ( nLogLevel >= 2) then
        write(6,*) "R-Space Energy: ",dPotR 
        write(6,*) "----------------------------------------"
    end if
              
    nD_GetDirectMadelungPotential=dPotR

end function


end module
