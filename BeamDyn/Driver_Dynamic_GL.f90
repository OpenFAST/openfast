!..................................................................................................................................
! LICENSING                                                                                                                         
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    Glue is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with Module2.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!  
!    ADD DESCRIPTION
!	
!    References:
!
!
!**********************************************************************************************************************************
PROGRAM MAIN

!   USE BeamDynDynamicGL
   USE BeamDynDynamicLSGL
   USE BeamDyn_Types

   USE NWTC_Library

   IMPLICIT NONE

   ! global glue-code-specific variables

   INTEGER(IntKi)                     :: ErrStat          ! Error status of the operation
   CHARACTER(1024)                    :: ErrMsg           ! Error message if ErrStat /= ErrID_None

   REAL(DbKi)                         :: dt_global        ! fixed/constant global time step

   REAL(DbKi)::t_initial
   REAL(DbKi)::t_final
   REAL(DbKi)::rhoinf,alfam,alfaf,gama,beta
   REAL(DbKi)::coef(9)
   REAL(DbKi)::StepEndTime

   INTEGER(IntKi)                     :: n_t_final        ! total number of time steps
   INTEGER(IntKi)                     :: n_t_global       ! global-loop time counter

   INTEGER(IntKi)                     :: pc_max           ! 1:explicit loose; 2:pc loose
   INTEGER(IntKi)                     :: pc               ! counter for pc iterations

   INTEGER(IntKi)                     :: BDyn_interp_order     ! order of interpolation/extrapolation
   
   INTEGER(IntKi),PARAMETER:: OutUnit = 10
   INTEGER(IntKi),PARAMETER:: QiDisUnit = 20
   INTEGER(IntKi),PARAMETER:: QiForUnit = 21
   INTEGER(IntKi),PARAMETER:: QiRMSu1Unit = 22
   INTEGER(IntKi),PARAMETER:: QiRMSu2Unit = 23
   INTEGER(IntKi),PARAMETER:: QiRMSu3Unit = 24

   ! BeamDyn Derived-types variables; see Registry_BeamDyn.txt for details

   TYPE(BDyn_InitInputType)           :: BDyn_InitInput
   TYPE(BDyn_ParameterType)           :: BDyn_Parameter
   TYPE(BDyn_ContinuousStateType)     :: BDyn_ContinuousState
   TYPE(BDyn_ContinuousStateType)     :: BDyn_ContinuousStateDeriv
   TYPE(BDyn_InitOutputType)          :: BDyn_InitOutput
   TYPE(BDyn_DiscreteStateType)       :: BDyn_DiscreteState
   TYPE(BDyn_ConstraintStateType)     :: BDyn_ConstraintState
   TYPE(BDyn_OtherStateType)          :: BDyn_OtherState

   TYPE(BDyn_InputType),Dimension(:),Allocatable  :: BDyn_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE           :: BDyn_InputTimes

   TYPE(BDyn_OutputType),Dimension(:),Allocatable  :: BDyn_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: BDyn_OutputTimes

   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops
   Integer(IntKi)                     :: j               ! counter for various loops
   REAL(DbKi):: DoubleTest
   REAL(ReKi):: SingleTest

   INTERFACE
       SUBROUTINE TiSchmComputeCoefficients(beta,gama,dt_global,alfam,alfaf,coef)
       USE BeamDyn_Types
       REAL(DbKi),INTENT(IN)::beta,gama,alfam,alfaf
       REAL(DbKi),INTENT(IN)::dt_global
       REAL(DbKi),INTENT(INOUT)::coef(:)
       END SUBROUTINE TiSchmComputeCoefficients
   END INTERFACE

   INTERFACE
      SUBROUTINE UpdateStructuralConfiguration(uuNf,vvNf,aaNf,xxNf,uuNi,vvNi,aaNi,xxNi)
      USE BeamDyn_Types
      REAL(ReKi),INTENT(IN)::uuNf(:),vvNf(:),aaNf(:),xxNf(:)
      REAL(ReKi),INTENT(INOUT)::uuNi(:),vvNi(:),aaNi(:),xxNi(:)
      END SUBROUTINE UpdateStructuralConfiguration
   END INTERFACE

    OPEN(unit = OutUnit, file = 'Dynamic.out', status = 'REPLACE',ACTION = 'WRITE')
    OPEN(unit = QiDisUnit, file = 'QiDisp.out', status = 'REPLACE',ACTION = 'WRITE')
    OPEN(unit = QiForUnit, file = 'QiForce.out', status = 'REPLACE',ACTION = 'WRITE')
    OPEN(unit = QiRMSu1Unit, file = 'QiRMSu1.out', status = 'REPLACE',ACTION = 'WRITE')
    OPEN(unit = QiRMSu2Unit, file = 'QiRMSu2.out', status = 'REPLACE',ACTION = 'WRITE')
    OPEN(unit = QiRMSu3Unit, file = 'QiRMSu3.out', status = 'REPLACE',ACTION = 'WRITE')

   DoubleTest = 1.
   SingleTest = 1.

   WRITE(*,*) "DoubleTest = ", DoubleTest
   WRITE(*,*) "SingleTest = ", SingleTest

   t_initial = 0.0D0
   t_final = 9.0D0
   
   dt_global = 1.0D-03
   
   n_t_final = ((t_final - t_initial) / dt_global)
   
   rhoinf = 1.0D0
   alfam = 0.0D0
   alfaf = 0.0D0
   gama = 0.0D0
   beta = 0.0D0
   
   coef = 0.0D0

   Allocate(BDyn_Input(1))
   Allocate(BDyn_InputTimes(1))

   Allocate(BDyn_Output(1))
   Allocate(BDyn_OutputTimes(1))


   CALL BDyn_Init( BDyn_InitInput        &
                   , BDyn_Input(1)         &
                   , BDyn_Parameter        &
                   , BDyn_ContinuousState  &
                   , BDyn_DiscreteState    &
                   , BDyn_ConstraintState  &
                   , BDyn_OtherState       &
                   , BDyn_Output(1)        &
                   , dt_global       &
                   , BDyn_InitOutput       &
                   , ErrStat               &
                   , ErrMsg )
                   
                   
   CALL TiSchmComputeParameters(rhoinf, alfam,alfaf,gama,beta)
   CALL TiSchmComputeCoefficients(beta,gama,dt_global,alfam,alfaf,coef)

   ! -------------------------------------------------------------------------
   ! Time-stepping loops
   ! -------------------------------------------------------------------------  


   DO n_t_global = 0,n_t_final
       WRITE(*,2000) n_t_global+1,(n_t_global)*dt_global 
!       WRITE(*,*) "***INITIAL TIME = ", (n_t_global)*dt_global
!       IF(n_t_global==2) STOP
       StepEndTime =  (n_t_global + 1) * dt_global
!       WRITE(*,*) "StepEndTime = ", StepEndTime
       CALL DynamicSolutionGL(BDyn_Parameter%uuN0,BDyn_OtherState%uuNi,BDyn_OtherState%vvNi,BDyn_OtherState%aaNi,&
                        &BDyn_OtherState%xxNi,BDyn_OtherState%uuNf,BDyn_OtherState%vvNf,BDyn_OtherState%aaNf,&
                        &BDyn_OtherState%xxNf,&
                        &BDyn_Parameter%Stif0,BDyn_Parameter%m00,BDyn_Parameter%mEta0,BDyn_Parameter%rho0,BDyn_Parameter%bc,&
                        &BDyn_Parameter%node_elem, BDyn_Parameter%dof_node,&
                        &BDyn_Parameter%elem_total, BDyn_Parameter%dof_total,BDyn_Parameter%node_total,&
                        &coef,BDyn_Parameter%niter,BDyn_Parameter%ngp,dt_global,StepEndTime)
                        
!       CALL ComputeRootForce(BDyn_Parameter%uuN0,BDyn_OtherState%uuNf,&
!                 &BDyn_Parameter%Stif0,BDyn_Parameter%node_elem,BDyn_Parameter%dof_node,&
!                 &BDyn_Parameter%ngp,BDyn_OtherState%RootForce)

       IF(n_t_global == 0) THEN
           WRITE(OutUnit,*) 'Initial Nodal Configurations (uuN0):'
           WRITE(OutUnit,*) '=========================================='
           WRITE(QiDisUnit,6000) 0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0
           WRITE(QiForUnit,6000) 0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0
           WRITE(QiRMSu1Unit,7000) 0.0D0,0.0D0
           WRITE(QiRMSu2Unit,7000) 0.0D0,0.0D0
           WRITE(QiRMSu3Unit,7000) 0.0D0,0.0D0
           DO i=1,BDyn_Parameter%node_total
           j = (i - 1) * BDyn_Parameter%dof_node
           WRITE(OutUnit,1000) i,BDyn_Parameter%uuN0(j+1),BDyn_Parameter%uuN0(j+2),BDyn_Parameter%uuN0(j+3),&
                                &BDyn_Parameter%uuN0(j+4),BDyn_Parameter%uuN0(j+5),BDyn_Parameter%uuN0(j+6)
           ENDDO
       ENDIF
       
       WRITE(OutUnit,3000) n_t_global+1
       WRITE(OutUnit,4000) (n_t_global)*dt_global
       WRITE(OutUnit,5000) (n_t_global + 1) * dt_global
                        
       CALL UpdateStructuralConfiguration(BDyn_OtherState%uuNf,BDyn_OtherState%vvNf,BDyn_OtherState%aaNf,BDyn_OtherState%xxNf,&
                                         &BDyn_OtherState%uuNi,BDyn_OtherState%vvNi,BDyn_OtherState%aaNi,BDyn_OtherState%xxNi)

!       DO i=BDyn_Parameter%dof_total-5,BDyn_Parameter%dof_total
!           WRITE(*,*) BDyn_OtherState%uuNf(i)
!       ENDDO
       j=BDyn_Parameter%dof_total
       WRITE(QiDisUnit,6000) (n_t_global+1)*dt_global,BDyn_OtherState%uuNf(j-5),BDyn_OtherState%uuNf(j-4),&
                            &BDyn_OtherState%uuNf(j-3),BDyn_OtherState%uuNf(j-2),BDyn_OtherState%uuNf(j-1),&
                            &BDyn_OtherState%uuNf(j)
!       WRITE(QiDisUnit,6000) (n_t_global+1)*dt_global,BDyn_OtherState%uuNf(1),BDyn_OtherState%uuNf(2),&
!                            &BDyn_OtherState%uuNf(3),BDyn_OtherState%uuNf(4),BDyn_OtherState%uuNf(5),&
!                            &BDyn_OtherState%uuNf(6)
       WRITE(QiForUnit,6000) (n_t_global+1)*dt_global,BDyn_OtherState%RootForce(1),BDyn_OtherState%RootForce(2),&
                            &BDyn_OtherState%RootForce(3),BDyn_OtherState%RootForce(4),&
                            &BDyn_OtherState%RootForce(5),BDyn_OtherState%RootForce(6)
       WRITE(QiRMSu1Unit,7000) (n_t_global+1)*dt_global,BDyn_OtherState%uuNf(j-5)
       WRITE(QiRMSu2Unit,7000) (n_t_global+1)*dt_global,BDyn_OtherState%uuNf(j-4)
       WRITE(QiRMSu3Unit,7000) (n_t_global+1)*dt_global,BDyn_OtherState%uuNf(j-3)
       WRITE(OutUnit,*) 'Nodal Displacements (uuNf):'
       WRITE(OutUnit,*) '=========================================='
       DO i=1,BDyn_Parameter%node_total
           j = (i - 1) * BDyn_Parameter%dof_node
           WRITE(OutUnit,1000) i,BDyn_OtherState%uuNf(j+1),BDyn_OtherState%uuNf(j+2),BDyn_OtherState%uuNf(j+3),&
                              &BDyn_OtherState%uuNf(j+4),BDyn_OtherState%uuNf(j+5),BDyn_OtherState%uuNf(j+6)
       ENDDO
       
       WRITE(OutUnit,*) 'Nodal Velocity (vvNf):'
       WRITE(OutUnit,*) '=========================================='
       DO i=1,BDyn_Parameter%node_total
           j = (i - 1) * BDyn_Parameter%dof_node
           WRITE(OutUnit,1000) i,BDyn_OtherState%vvNf(j+1),BDyn_OtherState%vvNf(j+2),BDyn_OtherState%vvNf(j+3),&
                              &BDyn_OtherState%vvNf(j+4),BDyn_OtherState%vvNf(j+5),BDyn_OtherState%vvNf(j+6)
       ENDDO
       
       WRITE(OutUnit,*) 'Nodal Acceleration (aaNf):'
       WRITE(OutUnit,*) '=========================================='
       DO i=1,BDyn_Parameter%node_total
           j = (i - 1) * BDyn_Parameter%dof_node
           WRITE(OutUnit,1000) i,BDyn_OtherState%aaNf(j+1),BDyn_OtherState%aaNf(j+2),BDyn_OtherState%aaNf(j+3),&
                              &BDyn_OtherState%aaNf(j+4),BDyn_OtherState%aaNf(j+5),BDyn_OtherState%aaNf(j+6)
       ENDDO


       WRITE(OutUnit,*) 'Root Force::'
       WRITE(OutUnit,*) '=========================================='
       WRITE(OutUnit,1000) 1,BDyn_OtherState%RootForce(1),BDyn_OtherState%RootForce(2),BDyn_OtherState%RootForce(3),&
                          &BDyn_OtherState%RootForce(4),BDyn_OtherState%RootForce(5),BDyn_OtherState%RootForce(6) 
   ENDDO    

   

   1000 FORMAT (' ',I5.2,6ES21.12)
   2000 FORMAT ('*TIME STEP NO:',I5.2,'       ','INITIAL TIME = ',ES12.5)
   3000 FORMAT ('TIME STE NO: ', I5.2)
   4000 FORMAT ('INITIAL TIME = ', ES12.5)
   5000 FORMAT ('TIME STEP END = ', ES12.5)
   6000 FORMAT (ES12.5,6ES21.12)
   7000 FORMAT (ES12.5,ES21.12)
   CLOSE (OutUnit)
   CLOSE (QiDisUnit)
   CLOSE (QiForUnit)
   CLOSE (QiRMSu1Unit)
   CLOSE (QiRMSu2Unit)
   CLOSE (QiRMSu3Unit)

   CALL BDyn_End( BDyn_Input(1), BDyn_Parameter, BDyn_ContinuousState, BDyn_DiscreteState, &
                    BDyn_ConstraintState, BDyn_OtherState, BDyn_Output(1), ErrStat, ErrMsg )

   do i = 2, BDyn_interp_order+1
      CALL BDyn_DestroyInput(BDyn_Input(i), ErrStat, ErrMsg )
      CALL BDyn_DestroyOutput(BDyn_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(BDyn_Input)
   DEALLOCATE(BDyn_Output)
   DEALLOCATE(BDyn_InputTimes)
   DEALLOCATE(BDyn_OutputTimes)


END PROGRAM MAIN


   SUBROUTINE TiSchmComputeParameters(rhoinf, alfam,alfaf,gama,beta)

   USE BeamDyn_Types

   REAL(DbKi),INTENT(IN)::rhoinf
   REAL(DbKi),INTENT(INOUT)::alfam, alfaf, gama, beta

   REAL(DbKi)::tr0


   tr0 = rhoinf + 1.0D0
   alfam = (2.0D0 * rhoinf - 1.0D0) / tr0
   alfaf = rhoinf / tr0
   gama = 0.5D0 - alfam + alfaf
   beta = 0.25 * (1.0D0 - alfam + alfaf) * (1.0D0 - alfam + alfaf)
   

   END SUBROUTINE TiSchmComputeParameters
   
   
   SUBROUTINE TiSchmComputeCoefficients(beta,gama,deltat,alfaM,alfaF,coef)

   USE BeamDyn_Types

   REAL(DbKi),INTENT(IN)::beta, gama, alfaM, alfaF

   REAL(DbKi),INTENT(IN)::deltat

   REAL(DbKi),INTENT(INOUT):: coef(:)

   REAL(DbKi)::oalfaM, tr0, tr1, tr2
   REAL(DbKi)::deltat2

   deltat2 = deltat * deltat
   oalfaM = 1.0D0 - alfaM
   tr0 =  alfaF / oalfaM
   tr1 = alfaM / oalfaM
   tr2 = (1.0D0 - alfaF) / oalfaM

   coef(1) = beta * tr0 * deltat2
   coef(2) = (0.5D0 - beta/oalfaM) * deltat2
   coef(3) = gama * tr0 * deltat
   coef(4) = (1.0D0 - gama / oalfaM) * deltat
   coef(5) = tr0
   coef(6) = -tr1
   coef(7) = gama * tr2 * deltat
   coef(8) = beta * tr2 * deltat2
   coef(9) = tr2 

   END SUBROUTINE TiSchmComputeCoefficients



   SUBROUTINE UpdateStructuralConfiguration(uuNf,vvNf,aaNf,xxNf,uuNi,vvNi,aaNi,xxNi)
   
   USE BeamDyn_Types

   REAL(ReKi),INTENT(IN)::uuNf(:),vvNf(:),aaNf(:),xxNf(:)

   REAL(ReKi),INTENT(INOUT)::uuNi(:),vvNi(:),aaNi(:),xxNi(:)
    
   uuNi = uuNf
   vvNi = vvNf
   aaNi = aaNf
   xxNi = xxNf      

   END SUBROUTINE UpdateStructuralConfiguration
