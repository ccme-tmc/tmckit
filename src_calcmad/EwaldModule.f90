!> This module contains functions of different kind of Ewald summation
MODULE EwaldMod

use calcmad_func
use structMod

CONTAINS

!Reference : Kolafa. et.al 1992 Molecular Simulation
real*8 function GetEwaldKError(fAlpha,nKmax,fLength,fVolume)
    !
    implicit none
    !definition
    real*8 fAlpha,fLength,fVolume
    integer nKmax
    
    real*8 fResult
    
    fResult = fAlpha*exp(-(pi*nKmax/fAlpha/fLength)**2)/(pi**2*nKmax**1.5)
    GetEwaldKError = fResult
end function
    
real*8 function GetEwaldRError(fAlpha,nRmax,fLength,fVolume)
    !
    implicit none
    !definition
    real*8 fAlpha,fLength,fVolume
    integer nRmax
    
    real*8 fResult
    fResult = sqrt(nRmax/(2*fVolume))*exp(-fAlpha**2*sqrt(nRmax*1.0*fLength))/(fAlpha**2*(nRmax*1.0*fLength)**2)
    GetEwaldRError = fResult
end function        
    
real*8 function GetEwaldErrorDifference(fAlpha,nRmax,nKmax,fLength,fVolume)
    implicit none
    real*8 fAlpha,fLength,fVolume
    integer nRmax,nKmax
    
    GetEwaldErrorDifference = GetEwaldKError(fAlpha,nKmax,fLength,fVolume)-GetEwaldRError(fAlpha,nRmax,fLength,fVolume)
end function
   
    
real*8 function OptimizeEwaldAlpha(fCellSize)
    !Get optimized alpha, try make r-space and k-space error eqaul at R/G=3 cutoff
    !use 2-way split
    !
    IMPLICIT NONE       
    
    !definition
    real*8 fCellSize(3,3)
    
    !local
    real*8 fKError,fRError,fMinError,fMaxError,fError
    real*8 fAlphaMin,fAlphaMax,fAlpha
    real*8 fVolume,fLength,fTemp ! abc^(1/3)
    integer nTestSpaceMax !Rmax and Kmax
    fAlphaMin=0.01
    fAlphaMax=10
    fAlpha = 1
    nTestSpaceMax = 5
    
    !set length to maximum in 3 axis
    fVolume=GetCellVolume(fCellSize)
    fLength = VectorToAbs(fCellSize(1:3,1))
    fTemp = VectorToAbs(fCellSize(1:3,2))
    if ( fTemp .gt. fLength ) then
        fLength = fTemp
    end if
    fTemp = VectorToAbs(fCellSize(1:3,3))
    if ( fTemp .gt. fLength ) then
        fLength = fTemp
    end if        
    !Converge limit:  |Kerr(x)/Rerr(x)-1|  < 0.000001
    
    fMinError = GetEwaldErrorDifference(fAlphaMin,nTestSpaceMax,nTestSpaceMax,fLength,fVolume)
    fMaxError = GetEwaldErrorDifference(fAlphaMax,nTestSpaceMax,nTestSpaceMax,fLength,fVolume)
    fError = 1
    
    if ( fMinError * fMaxError .gt. 0 ) then  !Error : Min and max doesn't contain the minimum value
         OptimizeEwaldAlpha=0
         return
    else
        do while ( abs(fError) .gt. 0.000001 )
            fAlpha = (fAlphaMax+fAlphaMin)/2
            fError = GetEwaldErrorDifference(fAlpha,nTestSpaceMax,nTestSpaceMax,fLength,fVolume)
            !write(6,*) fAlpha,":",fError
            if ( fError*fMinError .lt. 0 ) then
                fAlphaMax = fAlpha
                fMaxError = fError
            else
                fAlphaMin = fAlpha
                fMinError = fError
            end if
        end do
    end if
    OptimizeEwaldAlpha = fAlpha
    
end function


real*8 function GetEwaldMadelung(task)
    !
    IMPLICIT NONE
    !Declare
    type(taskT),pointer :: task
    
    !local var
    real*8 fRSize(3,3) ! reciprocal space vector
    real*8 fTemp,fTemp2,fTemp3,fResult,fResult2,fLastResult
    real*8 fG(3)
    real*8 fR(3)
    integer nGmax,nRmax,i,j,k,a1,a2 ! vector*3, atom index*2
    real*8 fVPUC ! volume of PUC
    complex*8 cOneParticleSum ! sum of  G * R for all atom

    type(structT),pointer :: struct
    struct => task%struct
   
    !calc vector
    !fTemp=dot_product(struct%fCellSize(1:3,1),CrossProduct(struct%fCellSize(1:3,2),struct%fCellSize(1:3,3)))
    fTemp=GetCellVolume(struct%fCellSize)
    fRSize(1:3,1)=2*pi/fTemp*CrossProduct(struct%fCellSize(1:3,2),struct%fCellSize(1:3,3))
    fRSize(1:3,2)=2*pi/fTemp*CrossProduct(struct%fCellSize(1:3,3),struct%fCellSize(1:3,1))
    fRSize(1:3,3)=2*pi/fTemp*CrossProduct(struct%fCellSize(1:3,1),struct%fCellSize(1:3,2))
    
    fVPUC=fTemp
    
    !write(6,*) "PUC Volumn: ", fVPUC
    !write(6,*) "Converge Standard: ", task%fConvLimit
    !write(6,*) "Ewald alpha parameter: ", task%fAlpha
    if ( nLogLevel >= 2) then
        write(6,*) "----------------------------------------"
    end if
    !Converge
    fResult=-10
    do nGmax=1,nSpaceMax,2
        fLastResult=fResult
        fResult=0
        do i=-nGmax,nGmax
            do j=-nGmax,nGmax
                do k=-nGmax,nGmax
                    if ( (i .eq. 0) .AND.( j .eq. 0) .AND. (k .eq. 0) ) then ! G!=0,0,0
                        cycle
                    end if
                    fG=i*fRSize(1:3,1)+j*fRSize(1:3,2)+k*fRSize(1:3,3) 
                    ! for each G
                    fTemp=exp(-VectorToAbs(fG)**2/4/task%fAlpha**2)*2*pi/fVPUC/(VectorToAbs(fG)**2) !exp(-G^2/4a^2) *2pi/V
                    
                    !Summation Method 1
                    
                    !Summation of a1, then sum * a2
                    
                    fResult2=0
                    
                         !if ( (i .eq. -3) .AND.( j .eq. -3) .AND. (k .eq. -3) ) then 
                         !   write(6,'("e-G^2 Energy",es)') fTemp
                         !   write(6,*) "ABS(G)",VectorToAbs(fG)
                        !end if
                    
                    cOneParticleSum = (0.0,0.0)
                    do a1=1,struct%nAtom
                        !sum atom1 first
                        fTemp2 = dot_product(fG,struct%fCoord(1:3,a1))
                        cOneParticleSum=cOneParticleSum+struct%fCharge(a1)*cos(fTemp2)-(0.0d0,1.0d0)*sin(fTemp2)
                    end do
                    
                    do a2=1,struct%nAtom
                        !Multiply atom2 to sum of atom1
                        fTemp2 = dot_product(fG,struct%fCoord(1:3,a2))
                        fTemp3 = struct%fCharge(a2)*REAL(cOneParticleSum*(cos(fTemp2)+(0.0d0,1.0d0)*sin(fTemp2)))
                        fResult2= fResult2 + fTemp3
                        !if ( (i .eq. -3) .AND.( j .eq. -3) .AND. (k .eq. -3) ) then 
                        !    write(6,'("Atom",i1,"Energy",es)') a2,fTemp3
                        !end if
                    end do
                    
                    !Summation Method 2
                    
                    !Summation of a1*a2
                    
                    fResult2=0
                    do a1=1,struct%nAtom
                        do a2= 1,struct%nAtom
                            !if ( a1 .eq. a2 ) then ! different atom only
                            !    cycle
                            !end if
                            fTemp2=cos(dot_product(fG,struct%fCoord(1:3,a1)-struct%fCoord(1:3,a2)))
                            fResult2 = fResult2 + struct%fCharge(a1)*struct%fCharge(a2)*fTemp2
                        end do
                    end do
                    
                    fResult2=fResult2*fTemp
                    !write(6,*) "k-point ",i,j,k,": ",fResult2
                    fResult = fResult + fResult2
                 end do
             end do
        end do
        if ( nLogLevel >= 3) then
            write(6,*) "Gmax: ",nGmax,"    ",fResult 
        end if
        if ( (abs((fResult-fLastResult)/fResult) .lt. task%fConvLimit) .OR. (abs(fResult-fLastResult) .lt. 1d-20 )) then
            exit
       end if             
    end do
    !fResult = fResult*2*pi/fVPUC
    if ( nLogLevel >= 2) then
       write(6,*) "K-Space Energy: ",fResult 
    end if
    !Del complemtary
    fResult2=0
    do a1=1,struct%nAtom
       fResult2=fResult2-struct%fCharge(a1)**2/sqrt(pi)*task%fAlpha
    end do
    if ( nLogLevel >= 3) then
       write(6,*) "Self exclude: ",fResult2
    end if
    fResult = fResult + fResult2
    !Add other interact in G=0
    fResult2=-10
    do nRmax=1,nSpaceMax,3
        fLastResult=fResult2
        fResult2=0
        do i=-nRmax,nRmax
            do j=-nRmax,nRmax
                do k=-nRmax,nRmax
                    fR=i*struct%fCellSize(1:3,1)+j*struct%fCellSize(1:3,2)+k*struct%fCellSize(1:3,3)
                    do a1= 1,struct%nAtom
                        do  a2= 1,struct%nAtom
                             if (( a1 .eq. a2 ) .AND.  (i .eq. 0) .AND.( j .eq. 0) .AND. (k .eq. 0) )  then ! not for same atom in same lattice
                                 cycle
                             end if
                             fTemp=VectorToAbs(struct%fCoord(1:3,a1)-struct%fCoord(1:3,a2)+fR)
                            !fResult2=fResult2+struct%fCharge(a1)*struct%fCharge(a2)*(1-ERF(task%fAlpha*fTemp))/fTemp
                            fResult2=fResult2+struct%fCharge(a1)*struct%fCharge(a2)*DERFC(task%fAlpha*fTemp)/fTemp
                         end do
                     end do
                 end do
             end do
        end do
        fResult2 = fResult2/2.0
        if ( nLogLevel >= 3) then
            write(6,*) "Rmax: ",nRmax,"    ",fResult2 
        end if
        if ( abs((fResult2-fLastResult)/fResult2) .lt. task%fConvLimit .OR. (abs(fResult2) .lt. task%fConvLimit) ) then
            exit
       end if             
     end do
     if ( nLogLevel >= 2) then
        write(6,*) "R-Space Energy: ",fResult2 
        write(6,*) "----------------------------------------"
     end if
     
     fResult= fResult+fResult2

     GetEwaldMadelung=fResult                    
    end function


!> Calculate Madelung constant at site of atom  OR non-atom position in 2D or 3D cases
real*8 function nD_GetEwaldMadelungPotential(task,atom,pos)
    implicit none
    !Declare
    type(taskT),pointer,intent(in) :: task
    integer,intent(in),optional :: atom !< The index of atom. If set to 0 , then use fSite to calculate . Can be neglected
    real*8,intent(in),optional :: pos(3) !< The position of non-atom site. Only used when nIndex = 0. This argument can be neglected if not used
    !local var
    
    
    !local var
    real*8 dLast !< Record last result in converge iteration
    real*8 fG(3)
    real*8 fR(3)
    integer nGmax,nRmax,i,j,k,a2 ! vector*3, atom index(for 2nd atom)
    integer nVmax(3)  ! The maximum value of a,b,c-axis ( may change as dimension )
    real*8 dSelf
    real*8 dRabs
    real*8 dPotG,dPotR !< The total potential of G/R part
    real*8 dPotG2 !< G potential at specific point
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
    !Converge
    dPotG=-10
    do nGmax=1,nSpaceMax,2
        dLast=dPotG
        dPotG=0
        nVmax = 0
        do i=1,struct%nDimension
            nVmax(i) = nGmax
        end do
        do i=-nVmax(1),nVmax(1)
            do j=-nVmax(2),nVmax(2)
                do k=-nVmax(3),nVmax(3)
                    if ( (i .eq. 0) .AND. (j.eq.0)) then
                        if ( struct%nDimension .eq. 3 .AND. k .eq. 0) then
                            !Skip G=(0,0,0) in 3D cases
                            cycle
                        else
                            if ( struct%nDimension .eq. 2 ) then
                                dPotG = dPotG + D2D_EwaldPotentialK_00z(struct,task%fAlpha,dSite)
                                cycle
                            end if
                        end if

                    end if

                    fG=i*struct%fRSize(1:3,1)+j*struct%fRSize(1:3,2)+k*struct%fRSize(1:3,3) 

                    if ( struct%nDimension .eq. 3 ) then
                        dPotG2 = D3D_EwaldPotentialK(struct,task%fAlpha,fG,dSite)
                    else
                        dPotG2 =D2D_EwaldPotentialK(struct,task%fAlpha,fG,dSite)
                    end if
                    

                    if ( nLogLevel .ge. 4 .and. abs(i) .le. 1 .and. abs(j) .le. 1 ) then
                        write(6,fmt='("k-point",3i4,": ",d14.7)') i,j,k,dPotG2
                    end if
                    dPotG = dPotG + dPotG2
                end do
            end do
        end do
        if ( nLogLevel >= 3) then
            write(6,*) "Gmax: ",nGmax,"    ",dPotG 
        end if
        !write(6,*) "converge:",dPotG,dLast,abs((dPotG-dLast)/dPotG)
        if ( isnan(dPotG)) then
            stop "NaN found !" 
        end if 
        !Deal with dPotG == 0 ( This may happen at nVmax=1,1,1 for Oh system )
        !Sometimes a very small value may occur due to floating point precision 
        !Always go to next cycle 
        if ( dPotG == 0d0 ) then
          cycle
        end if
        if ( (abs((dPotG-dLast)/dPotG) .lt. task%fConvLimit) .OR. (abs(dPotG-dLast) .lt. 1d-20 )) then
            exit
        end if             
    end do
     
    if ( nLogLevel >= 2) then
        write(6,*) "K-Space Energy: ",dPotG 
    end if
     
    !Del self energy
    dSelf = 0D0
    if ( nIndex .ne. 0) then
        dSelf=-struct%fCharge(nIndex)/sqrt(pi)*task%fAlpha*2
        if ( nLogLevel >= 2) then
            write(6,*) "Self exclude: ",dSelf
        end if
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
                            if (( nIndex .eq. a2 ) .AND.  (i .eq. 0) .AND.( j .eq. 0) .AND. (k .eq. 0) )  then ! not for same atom in same lattice
                                cycle
                            end if
                            dRabs=VectorToAbs(dSite-struct%fCoord(1:3,a2)+fR)
                            !write(6,*) fTemp,struct%fCharge(a2)*DERFC(task%fAlpha*fTemp)/fTemp,dPotR,DERFC(task%fAlpha*fTemp)
                            dPotR=dPotR+struct%fCharge(a2)*DERFC(task%fAlpha*dRabs)/dRabs
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
              
    nD_GetEwaldMadelungPotential=dPotG+dPotR+dSelf

end function

!> Calculate 3D Ewald G-space part Potential
!! We do not need to care about whether the position has an atom or not here
real*8 function D3D_EwaldPotentialK(struct,fAlpha,dG,dSite)
    !
    IMPLICIT NONE
    !
    !Declaration
    !
    type(structT),pointer,intent(in) :: struct !< The structure of lattice
    real*8 :: fAlpha !< The alpha value in calculation
    real*8 :: dG(3) !< The G-vector . It cannot be (0,0,0)!
    real*8 :: dSite(3) !< The position where potential should be calculated
    !
    !Local
    !
    real*8 dCoef1,d_eikr,dTot
    integer :: i

    dCoef1=exp(-VectorToAbs(dG)**2/4/fAlpha**2)*4*pi/struct%fVPUC/(VectorToAbs(dG)**2) !exp(-G^2/4a^2) *2pi/V
                    
    dTot=0
    do i=1,struct%nAtom
        d_eikr=cos(dot_product(dG,dSite-struct%fCoord(1:3,i)))
        dTot = dTot + struct%fCharge(i)*d_eikr
    end do
    
    dTot=dTot*dCoef1

    D3D_EwaldPotentialK = dTot

end function

!> Calculate 2D Ewald G-space part Potential
!! We do not need to care about whether the position has an atom or not here
real*8 function D2D_EwaldPotentialK(struct,fAlpha,dG,dSite)
    !
    IMPLICIT NONE
    !
    !Declaration
    !
    type(structT),pointer,intent(in) :: struct !< The structure of lattice
    real*8,intent(in) :: fAlpha !< The alpha value in calculation
    real*8,intent(in) :: dG(3) !< The G-vector. It must not be (0,0,x)
    real*8,optional,intent(in) :: dSite(3) !< The position where potential should be calculated
    !
    !Local
    !
    real*8 dCoef1,d_eikr,dTot
    real*8 dRz,dGAbs,dIntegrate
    integer :: i

    dGAbs = VectorToAbs(dG)
    dCoef1 = pi / struct%fSPUC/dGAbs
    
    
    dTot=0
    do i=1,struct%nAtom
!If it is 2D case, multiply the limit of integration of the
!third dimensional 
!Assume k_{xy} = k
!Equantion : \pi/2*k*(
!exp(k*r)erfc(k/2a+a*r)+exp(-k*r)erfc(k/2a-a*r)
        dRz = dSite(3)-struct%fCoord(3,i) !z-component 
        d_eikr=cos(dot_product(dG,dSite-struct%fCoord(1:3,i))) ! e^{iGR},we just directly use GR ( include xyz ) as z is always 0
        !Numeric problem : exp and erfc are very easy to overflow !
        !As erfc(z) ~= exp(-z^2)*(C+R(1/2))/z,it will approach 0 at the
        !limit of dGAbs -> \infty
        !So we safely ignore when dGAbs is large enough

        if ( dGAbs*dRz .gt. dMaxExponent ) then
        !if ( .false. ) then
            dIntegrate = 0D0
        else 
            dIntegrate = exp(dGAbs*dRz)*DERFC(dGAbs/2/fAlpha+fAlpha*dRz)+exp(-dGAbs*dRz)*DERFC(dGAbs/2/fAlpha-fAlpha*dRz)
        end if
        if ( nLogLevel .ge. 5 ) then
            write(6,*) dGAbs,dRz,d_eikr,dIntegrate
        end if
        dTot = dTot + struct%fCharge(i)*d_eikr*dIntegrate
    end do
    
    dTot=dTot*dCoef1
    
    D2D_EwaldPotentialK = dTot

end function

!> Calculate 2D Ewald G-space part Potential, summation of all G of (0,0,z) to
!! the limit of c-> \infty
real*8 function D2D_EwaldPotentialK_00z(struct,fAlpha,dSite)
    !
    IMPLICIT NONE
    !
    !Declaration
    !
    type(structT),pointer,intent(in) :: struct !< The structure of lattice
    real*8,intent(in) :: fAlpha !< The alpha value in calculation
    real*8,optional,intent(in) :: dSite(3) !< The position where potential should be calculated
    !
    !Local
    !
    real*8 dCoef1,dTot
    real*8 dRz,dIntegrate
    integer :: i

    ! G=(0,0,z) has contribution to 2D results
    ! it is not nelegible
    ! but we need do it at limit c->\infty
    ! it is an integration contains singularity for each atom alone !
    ! But the singularity can be removed as \sum_iQ_i=0
    ! And what left is what we need
    ! Refer Parry 1975 ( Erratum )
    ! 
    !
    dCoef1 = 2*pi / struct%fSPUC
    dTot = 0d0
    do i=1,struct%nAtom
        dRz = dSite(3)-struct%fCoord(3,i) !z-component 
        !write(6,*) dRz
        dIntegrate =-dRz*DERF(fAlpha*dRz)-1/sqrt(pi)/fAlpha*exp(-fAlpha**2*dRz**2)
        !write(6,*) dIntegrate
        !dIntegrate = 2*pi/fVPUC*struct%fCellSize(3,3)*dIntegrate
        !write(6,*) dIntegrate
        dTot = dTot + struct%fCharge(i)*dIntegrate
    end do
    dTot = dTot * dCoef1
    if ( nLogLevel .ge. 4 ) then
        write(6,fmt='("k-point",2i4,a4,": ",d14.7)') 0,0,"x",dTot
    end if
    D2D_EwaldPotentialK_00z = dTot

end function


END MODULE
