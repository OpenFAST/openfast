
!subroutine sc_init_obfuscator()
!
!   CALL RANDOM_SEED ( SIZE = 1 )
!   CALL RANDOM_SEED ( PUT=3459872 )      
!  
!end subroutine sc_init_obfuscator   
   
   
   
!=======================================================================
!SUBROUTINE sc_init (  ) BIND (C, NAME='sc_init')
!subroutine sc_init ( nTurbines, nInpGlobal, NumCtrl2SC, NumParamGlobal, ParamGlobal, NumParamTurbine, &
!                        ParamTurbine, NumStatesGlobal, NumStatesTurbine, NumSC2CtrlGlob, &
!                        NumSC2Ctrl, errStat, errMsg )  bind (C, NAME='sc_init')
subroutine sc_init ( nTurbines, nInpGlobal, NumCtrl2SC, NumParamGlobal,  NumParamTurbine, &
                         NumStatesGlobal, NumStatesTurbine, NumSC2CtrlGlob, &
                        NumSC2Ctrl, errStat, errMsg )  bind (C, NAME='sc_init')
!subroutine sc_init ( t, nTurbines, nInpGlobal, to_SCglob, NumCtrl2SC, to_SC, &
!                        nStatesGlobal, StatesGlob, nStatesTurbine, StatesTurbine, NumSC2CtrlGlob, from_SCglob, &
!                        NumSC2Ctrl, from_SC, errStat, errMsg )  bind (C, NAME='sc_calcOutputs')
         

   ! This DLL super controller is used to implement a ...
   
   ! Modified by B. Jonkman to conform to ISO C Bindings (standard Fortran 2003) and 
   ! compile with either gfortran or Intel Visual Fortran (IVF)
   ! DO NOT REMOVE or MODIFY LINES starting with "!DEC$" or "!GCC$"
   ! !DEC$ specifies attributes for IVF and !GCC$ specifies attributes for gfortran
   !
   ! Note that gfortran v5.x on Mac produces compiler errors with the DLLEXPORT attribute,
   ! so I've added the compiler directive IMPLICIT_DLLEXPORT.
   
   use, intrinsic :: ISO_C_Binding

   implicit                        none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: sc_init
!GCC$ ATTRIBUTES DLLEXPORT :: sc_init
#endif
   integer(C_INT),            intent(in   ) :: nTurbines         !< number of turbines connected to this supercontroller
   integer(C_INT),            intent(  out) :: nInpGlobal          !< number of global inputs to supercontroller
   integer(C_INT),            intent(  out) :: NumCtrl2SC          !< number of turbine controller outputs [inputs to supercontroller]
   integer(C_INT),            intent(  out) :: NumParamGlobal      !< number of global parameters
   integer(C_INT),            intent(  out) :: NumParamTurbine     !< number of parameters per turbine
   integer(C_INT),            intent(  out) :: NumStatesGlobal       !< number of global states
   integer(C_INT),            intent(  out) :: NumStatesTurbine      !< number of states per turbine
   integer(C_INT),            intent(  out) :: NumSC2CtrlGlob      !< number of global controller inputs [from supercontroller]
   integer(C_INT),            intent(  out) :: NumSC2Ctrl          !< number of turbine specific controller inputs [output from supercontroller]
   integer(C_INT),            intent(  out) :: errStat             !< error status code (uses NWTC_Library error codes)
   character(kind=C_CHAR),    intent(inout) :: errMsg          (*) !< Error Message from DLL to simulation code        
   
   !errMsg = TRANSFER( TRIM(avcMSG)//C_NULL_CHAR, avcMSG, SIZE(avcMSG) )
   errStat           = 0
   !errMsg            = ''
   
   nInpGlobal        = 0
   NumCtrl2SC        = 2 
   NumParamGlobal    = 5
   NumParamTurbine   = 4
   NumStatesGlobal   = 1
   NumStatesTurbine  = 2
   NumSC2CtrlGlob    = 2 
   NumSC2Ctrl        = 3 
    
  
   return
   
   end subroutine sc_init
subroutine sc_getInitData(nTurbines, NumParamGlobal, NumParamTurbine, ParamGlobal, ParamTurbine, &
        NumSC2CtrlGlob, from_SCglob, NumSC2Ctrl, from_SC,&
        & nStatesGlobal, StatesGlob, nStatesTurbine, StatesTurbine,&
        & errStat, errMsg )  bind (C, NAME='sc_getInitData')
use, intrinsic :: ISO_C_Binding

   implicit                        none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: sc_getInitData
!GCC$ ATTRIBUTES DLLEXPORT :: sc_getInitData
#endif
   integer(C_INT),            intent(in   ) :: nTurbines         !< number of turbines connected to this supercontroller
   integer(C_INT),            intent(in   ) :: NumParamGlobal      !< number of global parameters
   integer(C_INT),            intent(in   ) :: NumParamTurbine     !< number of parameters per turbine
   real(C_FLOAT),             intent(inout) :: ParamGlobal     (*) !< global parameters
   real(C_FLOAT),             intent(inout) :: ParamTurbine    (*) !< turbine-based parameters
   integer(C_INT),            intent(in   ) :: NumSC2CtrlGlob    !< number of global controller inputs [from supercontroller]
   real(C_FLOAT),             intent(inout) :: from_SCglob  (*)  !< global outputs of the super controller (to the turbine controller)
   integer(C_INT),            intent(in   ) :: NumSC2Ctrl        !< number of turbine specific controller inputs [output from supercontroller]
   real(C_FLOAT),             intent(inout) :: from_SC      (*)  !< turbine specific outputs of the super controller (to the turbine controller)
   integer(C_INT),            intent(in   ) :: nStatesGlobal     !< number of global states
   real(C_FLOAT),             intent(inout) :: StatesGlob   (*)  !< global states at time increment, n=0 (total of nStatesGlobal of these states)
   integer(C_INT),            intent(in   ) :: nStatesTurbine    !< number of states per turbine
   real(C_FLOAT),             intent(inout) :: StatesTurbine(*)  !< turbine-dependent states at time increment, n=0 (total of nTurbines*nStatesTurbine of these states)

   integer(C_INT),            intent(inout) :: errStat             !< error status code (uses NWTC_Library error codes)
   character(kind=C_CHAR),    intent(inout) :: errMsg          (*) !< Error Message from DLL to simulation code        
   integer                                  :: i,j
   real(C_FLOAT), allocatable               :: mask1(:)
   integer                                  :: seedVal(1), nSeeds
     
       ! Add a data obfuscator for your proprietary Parameter data
   
   
   
   !nSeeds     = 1
   !seedVal(1) = 3459872
   !call random_seed ( size = nSeeds  )
   !call random_seed ( put  = seedVal )     
   !allocate(mask1(NumParamGlobal), stat = errStat)
   !call random_number( mask1 )
   do i = 1, NumParamGlobal
      ParamGlobal(i) = real(0.6,C_FLOAT)  !real(i*mask1(i),C_FLOAT) 
   end do
   
   do j = 1, nTurbines
      do i = 1, NumParamTurbine
         ParamTurbine((j-1)*NumParamTurbine+i) = real((j-1)*NumParamTurbine+i,C_FLOAT)
      end do
   end do
   
   do i = 1, NumSC2CtrlGlob
      from_SCglob(i) = real(i,C_FLOAT)  !real(i*mask1(i),C_FLOAT) 
   end do
   
   do j = 1, nTurbines
      do i = 1, NumSC2Ctrl
         from_SC((j-1)*NumSC2Ctrl+i) = real((j-1)*NumSC2Ctrl+i,C_FLOAT)
      end do
   end do
     
   end subroutine sc_getInitData
!=======================================================================
!SUBROUTINE sc_calcOutputs (  ) BIND (C, NAME='sc_calcOutputs')                      
subroutine sc_calcOutputs ( t, nTurbines, nParamGlobal, paramGlobal, nParamTurbine, paramTurbine, nInpGlobal, to_SCglob, NumCtrl2SC, to_SC, &
                        nStatesGlobal, StatesGlob, nStatesTurbine, StatesTurbine, NumSC2CtrlGlob, from_SCglob, &
                        NumSC2Ctrl, from_SC, errStat, errMsg )  bind (C, NAME='sc_calcOutputs')
         

   ! This DLL super controller is used to implement a ...
   
   ! Modified by B. Jonkman to conform to ISO C Bindings (standard Fortran 2003) and 
   ! compile with either gfortran or Intel Visual Fortran (IVF)
   ! DO NOT REMOVE or MODIFY LINES starting with "!DEC$" or "!GCC$"
   ! !DEC$ specifies attributes for IVF and !GCC$ specifies attributes for gfortran
   !
   ! Note that gfortran v5.x on Mac produces compiler errors with the DLLEXPORT attribute,
   ! so I've added the compiler directive IMPLICIT_DLLEXPORT.
   
   use, intrinsic :: ISO_C_Binding

   implicit                        none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: sc_calcOutputs
!GCC$ ATTRIBUTES DLLEXPORT :: sc_calcOutputs
#endif
   
   real(C_DOUBLE),         INTENT(IN   ) :: t                 !< time (s)
   integer(C_INT),         intent(in   ) :: nTurbines         !< number of turbines connected to this supercontroller
   integer(C_INT),         intent(in   ) :: nParamGlobal        !< number of global parameters for the supercontroller
   real(C_FLOAT),          intent(in   ) :: paramGlobal    (*)  !< global parameters for the supercontroller
   integer(C_INT),         intent(in   ) :: nParamTurbine        !< number of turbine-based parameters for supercontroller
   real(C_FLOAT),          intent(in   ) :: paramTurbine    (*)  !< turbine-based parameters for the supercontroller
   integer(C_INT),         intent(in   ) :: nInpGlobal        !< number of global inputs to supercontroller
   real(C_FLOAT),          intent(in   ) :: to_SCglob    (*)  !< global inputs to the supercontroller
   integer(C_INT),         intent(in   ) :: NumCtrl2SC        !< number of turbine controller outputs [inputs to supercontroller]
   real(C_FLOAT),          intent(in   ) :: to_SC        (*)  !< inputs to the super controller (from the turbine controller)
   integer(C_INT),         intent(in   ) :: nStatesGlobal     !< number of global states
   real(C_FLOAT),          intent(in   ) :: StatesGlob   (*)  !< global states at time increment, n (total of nStatesGlobal of these states)
   integer(C_INT),         intent(in   ) :: nStatesTurbine    !< number of states per turbine
   real(C_FLOAT),          intent(in   ) :: StatesTurbine(*)  !< turbine-dependent states at time increment, n (total of nTurbines*nStatesTurbine of these states)
   integer(C_INT),         intent(in   ) :: NumSC2CtrlGlob    !< number of global controller inputs [from supercontroller]
   real(C_FLOAT),          intent(inout) :: from_SCglob  (*)  !< global outputs of the super controller (to the turbine controller)
   integer(C_INT),         intent(in   ) :: NumSC2Ctrl        !< number of turbine specific controller inputs [output from supercontroller]
   real(C_FLOAT),          intent(inout) :: from_SC      (*)  !< turbine specific outputs of the super controller (to the turbine controller)
   integer(C_INT),         intent(inout) :: errStat           !< error status code (uses NWTC_Library error codes)
   character(kind=C_CHAR), intent(inout) :: errMsg       (*)  !< Error Message from DLL to simulation code        
   integer                               :: i, j, c
   
   ! For this demo control we have:
   ! nInpGlobal        = 0
   ! NumCtrl2SC        = 2 
   ! NumParamGlobal    = 5
   ! NumParamTurbine   = 4
   ! NumStatesGlobal   = 1
   ! NumStatesTurbine  = 2
   ! NumSC2CtrlGlob    = 2 
   ! NumSC2Ctrl        = 3 
   
   !c = 1
   do j = 1, nTurbines  
      do i = 1, NumSC2Ctrl
         from_SC((j-1)*NumSC2Ctrl+i) = (j-1)*NumSC2Ctrl+i! StatesTurbine(c) + StatesTurbine(c+2)
         !from_SC((i-1)*NumSC2Ctrl+2) = StatesTurbine(c+1) + StatesTurbine(c+2)
         !c = c+3
      end do
   end do
   
   do i = 1, NumSC2CtrlGlob
      from_SCglob(i) = StatesGlob(1)
   end do
   
   !errMsg = TRANSFER( TRIM(avcMSG)//C_NULL_CHAR, avcMSG, SIZE(avcMSG) )
   return
end subroutine sc_calcOutputs

!=======================================================================
!SUBROUTINE sc_updateStates (  ) BIND (C, NAME='sc_updateStates')
subroutine sc_updateStates ( t, nTurbines, nParamGlobal, paramGlobal, nParamTurbine, paramTurbine, nInpGlobal, to_SCglob, NumCtrl2SC, to_SC, &
                        nStatesGlobal, StatesGlob, nStatesTurbine, StatesTurbine, errStat, errMsg )  bind (C, NAME='sc_updateStates')
         

   ! This DLL super controller is used to implement a ...
   
   ! Modified by B. Jonkman to conform to ISO C Bindings (standard Fortran 2003) and 
   ! compile with either gfortran or Intel Visual Fortran (IVF)
   ! DO NOT REMOVE or MODIFY LINES starting with "!DEC$" or "!GCC$"
   ! !DEC$ specifies attributes for IVF and !GCC$ specifies attributes for gfortran
   !
   ! Note that gfortran v5.x on Mac produces compiler errors with the DLLEXPORT attribute,
   ! so I've added the compiler directive IMPLICIT_DLLEXPORT.
   
   use, intrinsic :: ISO_C_Binding

   implicit                        none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: sc_updateStates
!GCC$ ATTRIBUTES DLLEXPORT :: sc_updateStates
#endif
   
   real(C_DOUBLE),         INTENT(IN   ) :: t                 !< time (s)
   integer(C_INT),         intent(in   ) :: nTurbines         !< number of turbines connected to this supercontroller
   integer(C_INT),         intent(in   ) :: nParamGlobal        !< number of global parameters for the supercontroller
   real(C_FLOAT),          intent(in   ) :: paramGlobal    (*)  !< global parameters for the supercontroller
   integer(C_INT),         intent(in   ) :: nParamTurbine        !< number of turbine-based parameters for supercontroller
   real(C_FLOAT),          intent(in   ) :: paramTurbine    (*)  !< turbine-based parameters for the supercontroller
   integer(C_INT),         intent(in   ) :: nInpGlobal        !< number of global inputs to supercontroller
   real(C_FLOAT),          intent(in   ) :: to_SCglob    (*)  !< global inputs to the supercontroller
   integer(C_INT),         intent(in   ) :: NumCtrl2SC        !< number of turbine controller outputs [inputs to supercontroller]
   real(C_FLOAT),          intent(in   ) :: to_SC        (*)  !< inputs to the super controller (from the turbine controller)
   integer(C_INT),         intent(in   ) :: nStatesGlobal     !< number of global states
   real(C_FLOAT),          intent(inout) :: StatesGlob   (*)  !< global states at time increment, n (total of nStatesGlobal of these states)
   integer(C_INT),         intent(in   ) :: nStatesTurbine    !< number of states per turbine
   real(C_FLOAT),          intent(inout) :: StatesTurbine(*)  !< turbine-dependent states at time increment, n (total of nTurbines*nStatesTurbine of these states)
   integer(C_INT),         intent(inout) :: errStat           !< error status code (uses NWTC_Library error codes)
   character(kind=C_CHAR), intent(inout) :: errMsg       (*)  !< Error Message from DLL to simulation code        
   integer                               :: i  
   real(C_FLOAT)                         :: sum
   ! Turbine-based inputs (one per turbine): to_SC
   ! 0 - Time
   ! 1 - GenTorque
   !
   ! Meaning of scOutputs
   ! 0 - Minimum Blade pitch
  
   ! Update the turbine-related states
   
   ! For this demo control we have:
   ! nInpGlobal        = 0
   ! NumCtrl2SC        = 2 
   ! NumParamGlobal    = 5
   ! NumParamTurbine   = 4
   ! NumStatesGlobal   = 1
   ! NumStatesTurbine  = 2
   ! NumSC2CtrlGlob    = 2 
   ! NumSC2Ctrl        = 3 
   sum = 0.0
   do i = 1, nTurbines*nStatesTurbine
   StatesTurbine(i) = i !paramGlobal(1)*to_SC(i)*paramTurbine(2*i-1) / paramTurbine(2*i) + (1-paramGlobal(1)*StatesTurbine(i))
   sum = sum + StatesTurbine(i)
   end do
   
   do i = 1,nStatesGlobal
      StatesGlob(i) = paramGlobal(2)*sum
   end do
   
   !double d2R = M_PI/180.0;
   ! Copy inputs into states first
   !for(int iTurb=0; iTurb < nTurbines; iTurb++) {
   !  for(int i=0; i < nScInputsTurbine; i++) {
   !    turbineStates_np1[iTurb][i] = sc_inputsTurbine[iTurb][i];
   !  }
   !}
   !
   !turbineStates_np1[0][nScInputsTurbine] = sc_inputsTurbine[0][0]/60.0 * 0.2 * d2R ;
   !turbineStates_np1[1][nScInputsTurbine] = sc_inputsTurbine[1][0]/60.0 * 0.45 * d2R ;
   !errMsg = TRANSFER( TRIM(avcMSG)//C_NULL_CHAR, avcMSG, SIZE(avcMSG) )
   
   return
end subroutine sc_updateStates

subroutine sc_end ( errStat, errMsg )  bind (C, NAME='sc_end')
         

   ! This DLL super controller is used to implement a ...
   
   ! Modified by B. Jonkman to conform to ISO C Bindings (standard Fortran 2003) and 
   ! compile with either gfortran or Intel Visual Fortran (IVF)
   ! DO NOT REMOVE or MODIFY LINES starting with "!DEC$" or "!GCC$"
   ! !DEC$ specifies attributes for IVF and !GCC$ specifies attributes for gfortran
   !
   ! Note that gfortran v5.x on Mac produces compiler errors with the DLLEXPORT attribute,
   ! so I've added the compiler directive IMPLICIT_DLLEXPORT.
   
   use, intrinsic :: ISO_C_Binding

   implicit                        none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: sc_end
!GCC$ ATTRIBUTES DLLEXPORT :: sc_end
#endif

   integer(C_INT),         intent(inout) :: errStat           !< error status code (uses NWTC_Library error codes)
   character(kind=C_CHAR), intent(inout) :: errMsg       (*)  !< Error Message from DLL to simulation code        
 

   
   return
end subroutine sc_end





