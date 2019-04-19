program CalcMad
         
    use calcmad_func
    use DirectMadMod
    use EwaldMod
           
    IMPLICIT none

    type(taskT),pointer :: task
    type(structT),pointer :: struct
    
    integer i
    real*8 fMinDistance
    real*8 fNewAlpha
    !Loop 
    real*8 dTotalEnergy
    real*8 dTmp ! Contains value before unit conversion
    real*8 dVr ! Potential at r
    real*8 dEatom ! Energy of atom
    
    !Test Parameter
    real*8 fNaCl
    
    !File Definition
    !2 input
    !3 output
    character(len=255) :: stName

    !Format
    character(len=34) ::  sFormat1='(" Core Atom List: ",????????(i) )'

    !command line parameter
    integer nArg

    !Read input file
    !open(2,file="calcmad.input",status="old")
    !open(3,file="calcmad.output",status="unknown")

    nArg = iargc()
    if ( nArg .eq. 2) then
        call getarg(1,stName)
        write(6,*) "Read from file "//Trim(stName)
        open(2,file=Trim(stName),status="old")
        call getarg(2,stName)
        open(3,file=Trim(stName),status="unknown")
    else
        call GetCaseName(stName)
        write(6,*) "Read from file "//Trim(stName)//".inmad"
        open(2,file=Trim(stName)//".inmad",status="old")
        open(3,file=Trim(stName)//".outputmad",status="unknown")
    end if

    !Init paramters
    allocate(task)
    struct => task%struct

    !Read Lattice Parameter
    call ReadTask(2,task)
    
    !Get min distance between atoms
    fMinDistance = GetMinimumAtomDistance(struct%fCellSize,struct%fCoord,struct%nAtom)
    write(6,*) "Minimum Distance: ", fMinDistance
    
    !output core atom list
    write(sFormat1(22:29),'(i8)'),task%nCoreCount !define count
    write(6,sFormat1) task%nCoreList(1:task%nCoreCount)
    
    
    !show parameter
    write(6,*) "Converge Standard: ", task%fConvLimit

    if ( task%stCalcMethod == "Ewald" ) then
      !Ewald Alpha test
      fNewAlpha = OptimizeEwaldAlpha(struct%fCellSize)
      write(6,*) "Ewald Alpha Optimized Estimate: ",fNewAlpha
      
      if ( task%fAlpha .eq. 0) then
          write(6,*) "Use new Ewald Alpha."
          task%fAlpha = fNewAlpha
      else
      write(6,*) "Ewald alpha parameter: ", task%fAlpha
      end if
    end if
    
        
    !Calculate Madelung constant only in 3D mode 
    if ( struct%nDimension .eq. 3 ) then
        fNaCl=-1.747564594633182190636  
        !3D mode
        !Ewald Method ( Total)
        write(6,*) "---------------------Overall Madelung Constant-------------------"
        dTmp = GetEwaldMadelung(task)
        dTotalEnergy = dTmp/4/pi/fEpsilon0*fElectronChargeUnit/fBohr
        write(6,*) "Electrostatic Energy: ",dTotalEnergy
        write(6,*) "Madelung Constant:", dTmp*fMinDistance/struct%nMFfactor
        write(3,'(g20.12,a1,g20.12,a1,g20.12)') dTmp*fMinDistance/struct%nMFfactor,CHAR(09),dTotalEnergy,CHAR(09),dTotalEnergy/fRydberg2eV
        !write(6,*) "Eletrostatic Energy Error: ",fNaCl/fMinDistance*nAtom/2-dTotalEnergy
    end if
    
    write(6,fmt=*) "------------ Atom-Site Madelung Potential--------"
    write(3,*) task%nCoreCount
    do i =1,task%nCoreCount
        write(6,fmt='(" Atom ",i4)') task%nCoreList(i)
        if ( task%stCalcMethod == "Ewald") then
            dTmp = nD_GetEwaldMadelungPotential(task,atom=task%nCoreList(i))
        else
            dTmp = nD_GetDirectMadelungPotential(task,atom=task%nCoreList(i))
        end if
        dVr = dTmp/4/pi/fEpsilon0*fElectronChargeUnit/fBohr
        dEatom = dVr*struct%fCharge(task%nCoreList(i))
        write(6,fmt='(" Electrostatic Potential : ",g20.12," V")') dVr
        write(6,fmt='(" Electrostatic Energy    : ",g20.12," eV")') dEatom
        write(6,fmt='(" Madelung Constant       : ",g20.12)') dTmp*fMinDistance
        write(3,'(g20.12,a1,g20.12,a1,g20.12,a1,g20.12)') dVr,CHAR(09),dVr/fRydberg2eV,CHAR(09),dEatom,CHAR(09),dEatom/fRydberg2eV
    end do
    write(6,*) "--------------------Non-atom Madelung Potential--------------------"
    write(3,*) task%nSiteCount
    do i =1,task%nSiteCount
        write(6,fmt='(" Site ",i4)') i 
        if ( task%stCalcMethod == "Ewald") then
            dVr=nD_GetEwaldMadelungPotential(task,pos=task%fSiteCoord(1:3,i))
        else
            dVr=nD_GetDirectMadelungPotential(task,pos=task%fSiteCoord(1:3,i))
        end if
        dVr = dVr/4/pi/fEpsilon0*fElectronChargeUnit/fBohr
        write(6,fmt='(" Electrostatic Potential : ",g20.12," V")') dVr
        write(3,'(g20.12,a1,g20.12)') dVr,CHAR(09),dVr/fRydberg2eV
    end do  

    !Cleanup
    close(2)        
    close(3)
    deallocate(task)
        
end program    
       
       
       
       
