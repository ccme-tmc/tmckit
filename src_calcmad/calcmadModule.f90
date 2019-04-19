module calcmad_func
         
        IMPLICIT none

        !Const
        real*8 pi,euler
        parameter (pi=3.14159265358979323846d0)
        parameter (euler=0.577215664901532860606d0)
        
        real*8 fElectronChargeUnit,fElectronCharge,fElectronMass,fElectronChargeMass,fLightSpeed,fEpsilon0,fPlanck,fRydberg2eV,fBohr
        parameter (fElectronChargeUnit=1.602176487D-19)
        parameter ( fElectronCharge=-1.602176487D-19)
        parameter ( fElectronChargeMass =-1.758820150D11)
        parameter (fElectronMass=9.10938215D-31)
        parameter ( fLightSpeed=299792458D0)
        parameter( fEpsilon0=8.85418781762038985053656303171D-12)
        parameter ( fPlanck =  6.6260684D-34)
        parameter ( fRydberg2eV = 13.6056942325D0)
        parameter ( fBohr = 5.291772108D-11)
        
        !Some restrict to floating point number
        real*8,parameter :: dMaxExponent = 100d0 !< The maximum absolute value of x in e^-x which is not treated as zero, to avoid overflow of double precision type
        
        !Common Parameter 
        integer nSpaceMax
        parameter (nSpaceMax=1000)!space vector maximum cutoff( in 3-D);include both real and reciprocal
        
        integer nLogLevel !< The details of log. 4 is most verbose, 0 least

        !Max number of atom
        integer,parameter :: nMaxAtomCount=10000
        integer,parameter :: nMaxNonAtomPosCount=10000



        !Function

        !real*8 AngleToRadian
        !interface
        !    function CrossProduct(fVector1,fVector2)
        !    implicit none
        !    real*8 fVector1(3),fVector2(3)
        !    real*8 CrossProduct(3)
        !    end function
        !end interface       
        
        
        
       
        
        

        !File Manipulate Functions

        Contains

        subroutine GetCaseName(stName)
        IMPLICIT none
        character(len=*),intent(out) :: stName

        !variable
        integer*4 nStatus,nPos
        character(len=255) stPath
        integer*4 getcwd

        nStatus = getcwd(stPath)
        nPos = index(stPath,"/",.True.)
        stName = stPath(nPos+1:len(stPath))
        end subroutine

        !Math Functions

        real*8 function AngleToRadian(x)
        IMPLICIT none
        real*8 x
        !real*8 AngleToRadian
        real*8 pi
        parameter (pi=3.14159265358979323846D0)
        !write(6,'(f14.7)') x
        AngleToRadian=x*pi/180.0
        end function

        real*8 function VectorToAbs(x)
        IMPLICIT none
        !real*8 VectorToAbs
        real*8 x(3)
        real*8 results
        results=sqrt(x(1)**2+x(2)**2+x(3)**2)
        VectorToAbs=results
        end function
        
        real*8 function VectorToAbs2D(x)
        IMPLICIT none
        !real*8 VectorToAbs
        real*8 x(2)
        real*8 results
        results=sqrt(x(1)**2+x(2)**2)
        VectorToAbs2D=results
        end function        
        
        real*8 function GetCellVolume(fCellSize)
        !Get volume of a cell
        IMPLICIT none
        real*8 fCellSize(3,3)
        real*8 results
        results=dot_product(fCellSize(1:3,1),CrossProduct(fCellSize(1:3,2),fCellSize(1:3,3)))
        GetCellVolume = results
        end function
        
        real*8 function Get2DCellVolume(fCellSize)
        IMPLICIT none
        real*8 fCellSize(3,2)
        real*8 results
        results=dot_product(fCellSize(1:3,1),fCellSize(1:3,2))
        Get2DCellVolume = results
        end function        
        
        function CrossProduct(fVector1,fVector2)
        IMPLICIT none
        !Get cross product of 2 3-D vector
        real*8 fVector1(3),fVector2(3)
        real*8 CrossProduct(3)
        CrossProduct(1)=fVector1(2)*fVector2(3)-fVector1(3)*fVector2(2)
        CrossProduct(2)=fVector1(3)*fVector2(1)-fVector1(1)*fVector2(3)
        CrossProduct(3)=fVector1(1)*fVector2(2)-fVector1(2)*fVector2(1)
        end function

        function Get2DReciprocalVector(fCellSize)
        IMPLICIT none
        !Get Reciprocal Vector from Real vector
        real*8 fCellSize(3,2)
        real*8 fRSize(3,2)
        real*8 Get2DReciprocalVector(3,2)
        real*8 fTemp
        !local var
        fTemp=Get2DCellVolume(fCellSize)
        fRSize(1:3,1)=2*pi/fTemp*fCellSize(1:3,2)
        fRSize(1:3,2)=2*pi/fTemp*fCellSize(1:3,1)
        Get2DReciprocalVector=fRSize
        end function

        function GetReciprocalVector(fCellSize)
        IMPLICIT none
        !Get Reciprocal Vector from Real vector
        real*8 fCellSize(3,3)
        real*8 fRSize(3,3)
        real*8 GetReciprocalVector(3,3)
        real*8 fTemp
        !local var
        fTemp=GetCellVolume(fCellSize)
        fRSize(1:3,1)=2*pi/fTemp*CrossProduct(fCellSize(1:3,2),fCellSize(1:3,3))
        fRSize(1:3,2)=2*pi/fTemp*CrossProduct(fCellSize(1:3,3),fCellSize(1:3,1))
        fRSize(1:3,3)=2*pi/fTemp*CrossProduct(fCellSize(1:3,1),fCellSize(1:3,2))        
        GetReciprocalVector=fRSize
        end function

        function GetMinimumImageDistance(fCellSize,fCoord)
        IMPLICIT none
        !Get coord by minium image convetion(<=0.5)
        real*8 fCellSize(3,3),fCoord(3),GetMinimumImageDistance(3)
        real*8 i
        do i=1,3
            if ( fCoord(i) > 0.5) then
                GetMinimumImageDistance(i) = 1-fCoord(i)
            else
                GetMinimumImageDistance(i) = fCoord(i) 
            end if
        end do
        end function
        

        real*8 function GetMinimumAtomDistance(fCellSize,fCoord,nAtom)
        !Get minimum distance of 2 atom in cell
        !If only one Atom, return minimum in fCellSize
        implicit none
        real*8 fCellSize(3,3)
        integer nAtom
        real*8 fCoord(3,100)
        
        real*8 i,j
        real*8 fTemp,fMin
        integer fMinIndex(2) ! Atom Index record

        fMin=100
        if ( nAtom  .gt. 1) then
            do i=1,nAtom
                do j=i+1,nAtom
                    !fTemp=VectorToAbs(GetMinimumImageDistance(fCoord(1:3,i)-fCoord(1:3,j)))
                    fTemp=VectorToAbs(fCoord(1:3,i)-fCoord(1:3,j))
                    if ( fTemp < fMin) then
                        fMin = fTemp
                        fMinIndex(1)=i
                        fMinIndex(2)=j
                    end if
                end do
            end do
        else
            do i = 1,3
                fTemp=VectorToAbs(fCellSize(1:3,i))
                if ( fTemp < fMin) then
                    fMin = fTemp
                end if
            end do
        end if
        write(6,'(" Minimum fund between ",i3, " and ",i3)') fMinIndex(1),fMinIndex(2)
        GetMinimumAtomDistance=fMin
        end function

        !>  Get three number na,nb,nc, to make sure na*|a| ~= nb*|b| ~= nc*|c|
        !!  ~= n*max(|a|,|b|,|c|) or n*min(|a|,|b|,|c|) 
        !!  Used to get homogenous grid in space from arbitary base vectors
        subroutine GetEqualDistanceTimes(fCellSize,n,bMax,nC)
        !
        IMPLICIT NONE
        !
        !Delclaration
        !
        real*8,intent(in) :: fCellSize(3,3) !< The base vectors
        integer,intent(in) :: n !< n*longest vector is the reference length 
        integer,intent(out) :: nC(3) !< na,nb,nc
        logical,intent(in) :: bMax !< If true, use longest base as reference; if false, useshortest
        !
        !Local
        !
        real*8 :: dL(3),dMax
        integer ::i

        do i = 1,3
            dL(i) = VectorToAbs(fCellSize(1:3,i))
        end do
        if ( bMax) then
            dMax = MaxVal(dL)
        else
            dMax = MinVal(dL)
        end if

        do i = 1,3
            nC(i) = NINT(n*dMax/dL(i))
            if ( nC(i) .eq. 0) then
                nC(i) = 1
            end if
        end do

        end subroutine

end module
