!> This module contains struct information
MODULE structMod

use calcmad_func

implicit NONE

!> Represent the structure of crystal in calculation
type :: structT
    real*8 :: fCoord(3,nMaxAtomCount) !< The cartesian coordinate of atoms
    real*8 :: fCharge(nMaxAtomCount) !< The charge of each atom
    real*8 :: fCellSize(3,3) !< The lattice vector of cell
    integer :: nDimension !< The dimension of system. If set to 2, then c-axis of cell is treated as non-repeat
    integer :: nAtom !< The number of atoms in the cell

    integer :: nInclude(nMaxAtomCount) 
    !< Preserved paramter 1 ( now used to indicate whether an atom is included in calculation

    integer :: nInclude2(nMaxAtomCount)     !< Preserved paramter 2

    integer :: nMFfactor !< The struct unit in one cell, used to calculate Madelung constant 

    real*8 :: fVPUC !< The volume of the unit cell. Calculated after read automatically
    real*8 :: fRSize(3,3) !< The reciporal lattice vector of cell. Calculated 
    real*8 :: fSPUC !< The area of a-b plane of unit cell. Used in 2D. Note use this value require the c-axis perpendicular to a-b plane. 
    

end type

!> Represent all parameters necessary in calculation
type :: taskT
    type(structT) :: struct !< The struct

    integer :: nCoreList(nMaxAtomCount) !< List of Core atom ( atom positions that will be calculated potential )
    integer :: nCoreCount !< The number of nCoreList

    integer :: nSiteCount !< The number of sites
    real*8 :: fSiteCoord(3,nMaxNonAtomPosCount) !< List of position where potentialshould be calculated but not position of atoms

    real*8 :: fConvLimit !< The convergence threshold in summation 
    !! when the  difference between two iteration is less then limit*(value of
    !! this iteration), treated as converge

    real*8 :: fAlpha !< The alpha in Ewald method

    character(len=3)  :: stCalcMode !< The cases of calculation ( 3D1 or 2D1 )
    character(len=10) :: stCalcMethod !< Method of calculation, Ewald or Direct
end type

CONTAINS

!> Read calculation information
SUBROUTINE ReadTask(fFile,task)
    !
    IMPLICIT NONE
    !
    !Declaration
    integer :: fFile
    type(taskT),pointer :: task
    !
    !Local
    !
    !Crystal Parameter
    real*8 fa,fb,fc,fAngleA,fAngleB,fAngleC
    !Cartesian Coordinate
    real*8 fbx,fby,fcx,fcy,fcz
    real*8 fCellAngle(3)
    real*8 fCellSizeABC(3)
    !Internal Coordinate and Charge(manually set)
    real*8 fCoordABC(3,nMaxAtomCount)
    real*8 fSiteCoordABC(3,nMaxNonAtomPosCount)
    integer :: i,j

    !Format control
    character(255) :: line,name_atom
    integer :: err

    !For easy access
    type(structT),pointer :: struct

    struct => task%struct

    read(fFile,*) task%stCalcMode,task%stCalcMethod
    read(fFile,*) struct%nMFfactor,task%fAlpha,task%fConvLimit,nLogLevel
    
    write(6,*) "Calculation lattice: ",task%stCalcMode
    write(6,*) "Calculation Method: ",task%stCalcMethod

    read(task%stCalcMode,fmt='(i1)') struct%nDimension
       
    read(fFile,*) fa,fb,fc
    read(fFile,*) fAngleA,fAngleB,fAngleC
    write(6,'(" Crystal Parameter: ",f10.4,f10.4,f10.4)') fa,fb,fc
    write(6,'(" Angles: ",f10.4,f10.4,f10.4)')  fAngleA,fAngleB,fAngleC  
    !Read atom site list
    read(2,*) struct%nAtom
    j=1
    do i = 1,struct%nAtom,1
    !The first item can be used as a comment if it is not a float
        read(2,'(a)') line
        read(line,*,iostat=err) name_atom,fCoordABC(1:3,i),struct%fCharge(i),struct%nInclude(i),struct%nInclude2(i)
        if (err .ne. 0) then
          if (i .eq. 1) write(6,'(a)') "Input file format: without atom names"
          read(line,*,iostat=err) fCoordABC(1:3,i),struct%fCharge(i),struct%nInclude(i),struct%nInclude2(i)
        else
          if (i .eq. 1) write(6,'(a)') "Input file format: with atom names"
        end if

        if ( struct%nInclude(i) .ne. 0 ) then
            task%nCoreList(j)=i
            j=j+1
        endif
    end do
    task%nCoreCount = j-1
    
    !Read non-atom site list
    read(2,*) task%nSiteCount
    do i = 1,task%nSiteCount,1
        read(2,*) fSiteCoordABC(1,i),fSiteCoordABC(2,i),fSiteCoordABC(3,i)
    end do

  
    !Convert to Radians
    !fAngleA=fAngleA*pi/180
    fAngleA=AngleToRadian(fAngleA)
    fAngleB=AngleToRadian(fAngleB)
    fAngleC=AngleToRadian(fAngleC)
    !write(6,'("Angle in radians: ",3(f10.4))') fAngleA,fAngleB,fAngleC
    
    fCellSizeABC(1)=fa
    fCellSizeABC(2)=fb
    fCellSizeABC(3)=fc
    
    fCellAngle(1)=fAngleA
    fCellAngle(2)=fAngleB
    fCellAngle(3)=fAngleC

    !Convert a-b-c Coord to 90-90-90 x-y-z Coordinate
    !direction: x=a, xy plane = ab plane
    fbx=fb*cos(fAngleC)
    fby=fb*sin(fAngleC)
    fcx=fc*cos(fAngleB)
    fcy=fc*(cos(fAngleA)-cos(fAngleB)*cos(fAngleC))/sin(fAngleC)
    fcz=sqrt((fc*sin(fAngleB))**2-fcy**2)
    !Genereate Cartesian vector a/b/c
    struct%fCellSize(1:3,1)=(/fa,DBLE(0.0),DBLE(0.0)/)
    struct%fCellSize(1:3,2)=(/fbx,fby,DBLE(0.0)/)
    struct%fCellSize(1:3,3)=(/fcx,fcy,fcz/)
    write(6,'(" a=",3(f14.7))') struct%fCellSize(1:3,1)
    write(6,'(" b=",3(f14.7))') struct%fCellSize(1:3,2)
    write(6,'(" c=",3(f14.7))') struct%fCellSize(1:3,3)
    !Convert internal Coord to x-y-z
    write(6,'("      No.      x             y             z         Charge")')
    do i=1,struct%nAtom
        struct%fCoord(1,i)=fCoordABC(1,i)*fa+fCoordABC(2,i)*fbx+fCoordABC(3,i)*fcx
        struct%fCoord(2,i)=fCoordABC(2,i)*fby+fCoordABC(3,i)*fcy
        struct%fCoord(3,i)=fCoordABC(3,i)*fcz
        write(6,'(" Atom",i3,":",3(f14.7)," ",f6.3)') i,struct%fCoord(1:3,i),struct%fCharge(i)
        !write(6,*) sqrt(fCoord(1,i)**2+fCoord(2,i)**2+fCoord(3,i)**2)
    enddo
    
    do i=1,task%nSiteCount
        task%fSiteCoord(1,i)=fSiteCoordABC(1,i)*fa+fSiteCoordABC(2,i)*fbx+fSiteCoordABC(3,i)*fcx
        task%fSiteCoord(2,i)=fSiteCoordABC(2,i)*fby+fSiteCoordABC(3,i)*fcy
        task%fSiteCoord(3,i)=fSiteCoordABC(3,i)*fcz            
    end do        

    !Pre-calculated some properies for later use
    struct%fRSize=GetReciprocalVector(struct%fCellSize)
    struct%fVPUC=GetCellVolume(struct%fCellSize)
    struct%fSPUC = struct%fVPUC / struct%fCellSize(3,3)

END SUBROUTINE

END MODULE
