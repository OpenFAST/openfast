MODULE GridInterp

USE GridInterp_Types

IMPLICIT NONE

PRIVATE SetIndex
PRIVATE GetN1D
PRIVATE GetN1Ddx

PUBLIC GridInterp_SetParams
PUBLIC GridInterpSetup3D
PUBLIC GridInterpSetup4D
PUBLIC GridInterpSetupN

INTERFACE GridInterp3D
   MODULE PROCEDURE GridInterp3DR4
   MODULE PROCEDURE GridInterp3DR8
END INTERFACE

INTERFACE GridInterp3DVec
   MODULE PROCEDURE GridInterp3DVecR4
   MODULE PROCEDURE GridInterp3DVecR8
END INTERFACE

INTERFACE GridInterp3DVec6
   MODULE PROCEDURE GridInterp3DVec6R4
   MODULE PROCEDURE GridInterp3DVec6R8
END INTERFACE

INTERFACE GridInterp4D
   MODULE PROCEDURE GridInterp4DR4
   MODULE PROCEDURE GridInterp4DR8
END INTERFACE

INTERFACE GridInterp4DVec
   MODULE PROCEDURE GridInterp4DVecR4
   MODULE PROCEDURE GridInterp4DVecR8
END INTERFACE

INTERFACE GridInterp4DVec6
   MODULE PROCEDURE GridInterp4DVec6R4
   MODULE PROCEDURE GridInterp4DVec6R8
END INTERFACE

INTERFACE GridInterp4DVecN
   MODULE PROCEDURE GridInterp4DVecNR4
   MODULE PROCEDURE GridInterp4DVecNR8
END INTERFACE

INTERFACE GridInterpN
   MODULE PROCEDURE GridInterpNR4
   MODULE PROCEDURE GridInterpNR8
END INTERFACE

INTERFACE GridInterpS
   MODULE PROCEDURE GridInterpSR4
   MODULE PROCEDURE GridInterpSR8
END INTERFACE

CONTAINS

Subroutine GridInterp_SetParams(dim, n, delta, pZero, IsPeriodic, p, ErrStat, ErrMsg)

   Integer(IntKi),                 intent(in   )  :: dim
   Integer(IntKi),                 intent(in   )  :: n(:)
   Real(ReKi),                     intent(in   )  :: delta(:)
   Real(ReKi),                     intent(in   )  :: pZero(:)
   Logical,                        intent(in   )  :: IsPeriodic(:)
   Type(GridInterp_ParameterType), intent(inout)  :: p
   
   Integer(IntKi),                 INTENT(OUT)    :: ErrStat
   Character(*),                   INTENT(OUT)    :: ErrMsg

   Integer(IntKi)                                 :: i

   ErrStat = ErrID_None
   ErrMsg  = ""

   if (dim/=3 .and. dim/=4) then
      ErrStat = ErrID_Fatal
      ErrMsg  = 'GridInterp_Init: dim must be 3 or 4'
      return
   end if

   do i = 1,dim
      p%n(i)          = n(i)
      p%delta(i)      = delta(i)
      p%pZero(i)      = pZero(i)
      p%IsPeriodic(i) = IsPeriodic(i)
   end do

End Subroutine GridInterp_SetParams

Subroutine SetIndex(pIn,pZero,delta,nMax,IsPeriodic,Indx,isopc,Support,FirstWarn,ErrStat,ErrMsg)

   Real(ReKi),       intent(in   )  :: pIn
   Real(ReKi),       intent(in   )  :: pZero
   Real(ReKi),       intent(in   )  :: delta
   Integer(IntKi),   intent(in   )  :: nMax
   Logical,          intent(in   )  :: IsPeriodic
   Integer(IntKi),   intent(inout)  :: Indx(:)
   Real(ReKi),       intent(inout)  :: isopc
   Integer(IntKi),   intent(inout)  :: Support     ! = 0 for linear interpolation, = 1 for quadratic with one point to the left, = 2 for quadratic with one point to the right, = 3 for cubic
   Logical,          intent(inout)  :: FirstWarn
   Integer(IntKi),   intent(  out)  :: ErrStat
   Character(*),     intent(  out)  :: ErrMsg

   Real(ReKi)                       :: p, pMax
   Integer(IntKi)                   :: i

   ErrStat = ErrID_None
   ErrMsg  = ""

   isopc = 0.0_ReKi
   Indx  = 0_IntKi

   if ( nMax .EQ. 1_IntKi ) then ! Only one grid point, effectively ignore this dimension
      ! Construct a dummy linear interpolation for now
      Indx(1) =  0_IntKi
      Indx(2) =  0_IntKi
      Indx(3) =  0_IntKi
      Indx(4) =  0_IntKi
      isopc   = 0.5_ReKi
      Support = 0
      return
   end if

   ! Compute normalized coordinate
   p = (pIn-pZero) / delta

   if (isPeriodic) then

      ! Calculate normalized coordinate between 0 and 1
      isopc = p - floor(p,ReKi)

      ! Get the normalized coordinates of the two nearest nodes to the left and to the right
      Indx = floor( p, IntKi ) + [-1,0,1,2]

      ! Make sure the coordinates are not out of bound using periodicity
      do i = 1,4
         if ( Indx(i) < 0 ) then
            Indx(i) = Indx(i) - nMax*floor( Real(Indx(i),ReKi) / Real(nMax,ReKi) )
         end if
         Indx(i) = mod(Indx(i),nMax)
      end do

      ! Always use cubic interpolation for periodic dimensions
      support = 3 ! Cubic interpolation

   else

      pMax = Real(nMax-1_IntKi,ReKi)

      if (p<0) then

         p = 0.0_ReKi

         if (FirstWarn) then
            call SetErrStat(ErrID_Warn,'Position has been clamped to the grid boundary. Warning will not be repeated though condition may persist.',ErrStat,ErrMsg,'SetIndex')
            FirstWarn = .false.
         end if

      else if (p>pMax) then

         p = pMax

         if (FirstWarn) then ! don't warn if we are exactly at the boundary
            call SetErrStat(ErrID_Warn,'Position has been clamped to the grid boundary. Warning will not be repeated though condition may persist.',ErrStat,ErrMsg,'SetIndex')
            FirstWarn = .false.
         end if

      end if

      if ( EqualRealNos( p, pMax ) ) then
         ! Calculate normalized coordinate between 0 and 1
         isopc = 1.0_ReKi
         ! Get the normalized coordinates of the two nearest nodes to the left and to the right
         Indx = nMax + [-3,-2,-1,0]
      else
         ! Calculate normalized coordinate between 0 and 1
         isopc = p - floor(p,ReKi)
         ! Get the normalized coordinates of the two nearest nodes to the left and to the right
         Indx = floor( p, IntKi ) + [-1,0,1,2]
      end if

      ! Supported interpolation method
      if ( Indx(1) < 0_IntKi ) then
         Indx(1) = 0_IntKi
         if (Indx(4) > (nMax-1_IntKi) ) then
            support = 0 ! Linear interpolation
            Indx(4) = nMax-1_IntKi
         else
            support = 1 ! Quadratic interpolation with only one node to the left
         end if
      else if ( Indx(4) > (nMax-1_IntKi) ) then
         support = 2 ! Quadratic interpolation with only one node to the right
         Indx(4) = nMax-1_IntKi
      else
         support = 3 ! Cubic interpolation
      end if

   end if

End Subroutine SetIndex

Subroutine GetN1D(isopc, support, N1D)

   real(ReKi),           intent(in   )  :: isopc       ! isoparametric coordinates
   integer(IntKi),       intent(in   )  :: support
   real(ReKi),           intent(inout)  :: N1D(4)
   real(ReKi)                           :: isopc2,isopc3

   select case ( Support )
   case ( 3 )  ! Cubic interpolation

      isopc2 = isopc*isopc
      isopc3 = isopc*isopc2

      N1D(1) = -0.5_ReKi*isopc3 +          isopc2 - 0.5_ReKi*isopc
      N1D(2) =  1.5_ReKi*isopc3 - 2.5_ReKi*isopc2                  + 1.0_ReKi
      N1D(3) = -1.5_ReKi*isopc3 + 2.0_ReKi*isopc2 + 0.5_ReKi*isopc
      N1D(4) =  0.5_ReKi*isopc3 - 0.5_ReKi*isopc2

   case ( 1 )  ! Quadratic interpolation with only one node to the left

      isopc2 = isopc*isopc

      N1D(1) =                                     0.0_ReKi
      N1D(2) =  0.5_ReKi*isopc2 - 1.5_ReKi*isopc + 1.0_ReKi
      N1D(3) = -1.0_ReKi*isopc2 + 2.0_ReKi*isopc
      N1D(4) =  0.5_ReKi*isopc2 - 0.5_ReKi*isopc

   case ( 2 )  ! Quadratic interpolation with only one node to the right

      isopc2 = isopc*isopc

      N1D(1) =  0.5_ReKi*isopc2 - 0.5_ReKi*isopc
      N1D(2) = -1.0_ReKi*isopc2                  + 1.0_ReKi
      N1D(3) =  0.5_ReKi*isopc2 + 0.5_ReKi*isopc
      N1D(4) =                                     0.0_ReKi

   case default ! Support == 0 Linear interpolation

      N1D(1) =          0.0_ReKi
      N1D(2) = -isopc + 1.0_ReKi
      N1D(3) =  isopc
      N1D(4) =          0.0_ReKi

   end select

End Subroutine GetN1D

Subroutine GetN1Ddx(isopc, support, N1D)

   real(ReKi),           intent(in   )  :: isopc       ! isoparametric coordinates
   integer(IntKi),       intent(in   )  :: support
   real(ReKi),           intent(inout)  :: N1D(4)
   real(ReKi)                           :: isopc2

   select case ( Support )
   case ( 3 )  ! Cubic interpolation

      isopc2 = isopc*isopc

      N1D(1) = -1.5_ReKi*isopc2 + 2.0_ReKi*isopc - 0.5_ReKi
      N1D(2) =  4.5_ReKi*isopc2 - 5.0_ReKi*isopc
      N1D(3) = -4.5_ReKi*isopc2 + 4.0_ReKi*isopc + 0.5_ReKi
      N1D(4) =  1.5_ReKi*isopc2 -          isopc

   case ( 1 )  ! Quadratic interpolation with only one node to the left

      N1D(1) =                   0.0_ReKi
      N1D(2) =           isopc - 1.5_ReKi
      N1D(3) = -2.0_ReKi*isopc + 2.0_ReKi
      N1D(4) =           isopc - 0.5_ReKi

   case ( 2 )  ! Quadratic interpolation with only one node to the right

      N1D(1) =           isopc - 0.5_ReKi
      N1D(2) = -2.0_ReKi*isopc
      N1D(3) =           isopc + 0.5_ReKi
      N1D(4) =                   0.0_ReKi

   case default ! Support == 0 Linear interpolation

      N1D(1) =  0.0_ReKi
      N1D(2) = -1.0_ReKi
      N1D(3) =  1.0_ReKi
      N1D(4) =  0.0_ReKi

   end select

End Subroutine GetN1Ddx

Subroutine GridInterpSetup3D( position, p, m, ErrStat, ErrMsg )

   real(ReKi),                          intent(in   )  :: Position(3)       !< Array of 3 coordinates
   type(GridInterp_ParameterType),      intent(in   )  :: p                 !< Parameters
   type(GridInterp_MiscVarType),        intent(inout)  :: m                 !< MiscVars
   integer(IntKi),                      intent(  out)  :: ErrStat           !< Error status
   character(*),                        intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   character(*), parameter              :: RoutineName = 'GridInterpSetup3D'
   integer(IntKi)                       :: dim,i,j,k
   integer(IntKi)                       :: support
   real(ReKi)                           :: N1D(4,3)
   real(ReKi)                           :: isopc       ! isoparametric coordinates
   integer(IntKi)                       :: ErrStat2
   character(ErrMsgLen)                 :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ""

   do dim = 1,3
      call SetIndex(Position(dim), p%pZero(dim), p%delta(dim), p%n(dim), p%IsPeriodic(dim), m%Indx(:,dim), isopc, Support, m%FirstWarn_Clamp, ErrStat2, ErrMsg2)
         if (Failed()) return;
      call GetN1D(isopc, Support, N1D(:,dim))
   end do

   do k = 1,4
      do j = 1,4
         do i = 1,4
            m%N3D(i,j,k) = N1D(i,1)*N1D(j,2)*N1D(k,3)
         end do
      end do
   end do

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function

End Subroutine GridInterpSetup3D

Subroutine GridInterpSetup4D( position, p, m, ErrStat, ErrMsg )

   real(ReKi),                          intent(in   )  :: Position(4)       !< Array of 4 coordinates
   type(GridInterp_ParameterType),      intent(in   )  :: p                 !< Parameters
   type(GridInterp_MiscVarType),        intent(inout)  :: m                 !< MiscVars
   integer(IntKi),                      intent(  out)  :: ErrStat           !< Error status
   character(*),                        intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   character(*), parameter              :: RoutineName = 'GridInterpSetup4D'
   integer(IntKi)                       :: dim,i,j,k,l
   integer(IntKi)                       :: support
   real(ReKi)                           :: N1D(4,4)
   real(ReKi)                           :: isopc       ! isoparametric coordinates
   integer(IntKi)                       :: ErrStat2
   character(ErrMsgLen)                 :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ""

   do dim = 1,4
      call SetIndex(Position(dim), p%pZero(dim), p%delta(dim), p%n(dim), p%IsPeriodic(dim), m%Indx(:,dim), isopc, Support, m%FirstWarn_Clamp, ErrStat2, ErrMsg2)
         if (Failed()) return;
      call GetN1D(isopc, Support, N1D(:,dim))
   end do

   do l = 1,4
      do k = 1,4
         do j = 1,4
            do i = 1,4
               m%N4D(i,j,k,l) = N1D(i,1)*N1D(j,2)*N1D(k,3)*N1D(l,4)
            end do
         end do
      end do
   end do

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function

End Subroutine GridInterpSetup4D

Subroutine GridInterpSetupN( position, p, m, ErrStat, ErrMsg )

   real(ReKi),                          intent(in   )  :: Position(3)       !< Array of 3 coordinates
   type(GridInterp_ParameterType),      intent(in   )  :: p                 !< Parameters
   type(GridInterp_MiscVarType),        intent(inout)  :: m                 !< MiscVars
   integer(IntKi),                      intent(  out)  :: ErrStat           !< Error status
   character(*),                        intent(  out)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None

   character(*), parameter              :: RoutineName = 'GridInterpSetupN'
   integer(IntKi)                       :: dim,i,j,k
   integer(IntKi)                       :: support
   real(ReKi)                           :: N1D(4,3)
   real(ReKi)                           :: N1Ddx(4,2:3)
   real(ReKi)                           :: isopc       ! isoparametric coordinates
   integer(IntKi)                       :: ErrStat2
   character(ErrMsgLen)                 :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ""

   do dim = 1,3
      call SetIndex(Position(dim), p%pZero(dim), p%delta(dim), p%n(dim), p%IsPeriodic(dim), m%Indx(:,dim), isopc, Support, m%FirstWarn_Clamp, ErrStat2, ErrMsg2)
         if (Failed()) return;
      call GetN1D(isopc, Support, N1D(:,dim))
      if (dim>1) then
         call GetN1Ddx(isopc, Support, N1Ddx(:,dim))
      end if
   end do

   ! Need two sets of weights for d(.)/dx and d(.)/dy. Borrow m%N4D for this.
   do k = 1,4
      do j = 1,4
         do i = 1,4
            m%N4D(i,j,k,1) = N1D(i,1)*N1Ddx(j,2)*N1D  (k,3)
            m%N4D(i,j,k,2) = N1D(i,1)*N1D  (j,2)*N1Ddx(k,3)
         end do
      end do
   end do

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function

End Subroutine GridInterpSetupN


!=============================================================================================================
! INTERFACE GridInterp3D
!              - GridInterp3DR4
!              - GridInterp3DR8
!=============================================================================================================
function GridInterp3DR4( data, m )
   real(SiKi),                   intent(in   )  :: data(0:,0:,0:)   !< 3D grid of scalar data
   type(GridInterp_MiscVarType), intent(in   )  :: m                !< MiscVars

   character(*), parameter                      :: RoutineName = 'GridInterp3DR4'
   real(SiKi)                                   :: GridInterp3DR4
   integer(IntKi)                               :: i,j,k

   ! interpolate
   GridInterp3DR4 = 0.0_SiKi
   do k = 1,4
      do j = 1,4
         do i = 1,4
            GridInterp3DR4 = GridInterp3DR4 + m%N3D(i,j,k) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3) )
         end do
      end do
   end do

end function GridInterp3DR4

function GridInterp3DR8( data, m )
   real(DbKi),                   intent(in   )  :: data(0:,0:,0:)   !< 3D grid of scalar data
   type(GridInterp_MiscVarType), intent(in   )  :: m                !< MiscVars

   character(*), parameter                      :: RoutineName = 'GridInterp3DR8'
   real(DbKi)                                   :: GridInterp3DR8
   integer(IntKi)                               :: i,j,k

   ! interpolate
   GridInterp3DR8 = 0.0_DbKi
   do k = 1,4
      do j = 1,4
         do i = 1,4
            GridInterp3DR8 = GridInterp3DR8 + m%N3D(i,j,k) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3) )
         end do
      end do
   end do

end function GridInterp3DR8

!=============================================================================================================
! INTERFACE GridInterp3DVec
!              - GridInterp3DVecR4
!              - GridInterp3DVecR8
!=============================================================================================================
function GridInterp3DVecR4( data, m )
   real(SiKi),                   intent(in   )  :: data(0:,0:,0:,:)   !< 3D grid of vector data
   type(GridInterp_MiscVarType), intent(in   )  :: m                  !< MiscVars

   character(*),   parameter                    :: RoutineName = 'GridInterp3DVecR4'
   integer(IntKi), parameter                    :: vDim = 3
   integer(IntKi)                               :: i,j,k,vi
   real(SiKi)                                   :: GridInterp3DVecR4(vDim)

   ! interpolate
   GridInterp3DVecR4 = 0.0_SiKi
   do k = 1,4
      do j = 1,4
         do i = 1,4
            do vi = 1,vDim
               GridInterp3DVecR4(vi) = GridInterp3DVecR4(vi) + m%N3D(i,j,k) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3), vi )
            end do
         end do
      end do
   end do

end function GridInterp3DVecR4

function GridInterp3DVecR8( data, m )
   real(DbKi),                   intent(in   )  :: data(0:,0:,0:,:)   !< 3D grid of vector data
   type(GridInterp_MiscVarType), intent(in   )  :: m                  !< MiscVars

   character(*),   parameter                    :: RoutineName = 'GridInterp3DVecR8'
   integer(IntKi), parameter                    :: vDim = 3
   integer(IntKi)                               :: i,j,k,vi
   real(DbKi)                                   :: GridInterp3DVecR8(vDim)

   ! interpolate
   GridInterp3DVecR8 = 0.0_DbKi
   do k = 1,4
      do j = 1,4
         do i = 1,4
            do vi = 1,vDim
               GridInterp3DVecR8(vi) = GridInterp3DVecR8(vi) + m%N3D(i,j,k) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3), vi )
            end do
         end do
      end do
   end do

end function GridInterp3DVecR8

!=============================================================================================================
! INTERFACE GridInterp3DVec6
!              - GridInterp3DVec6R4
!              - GridInterp3DVec6R8
!=============================================================================================================
function GridInterp3DVec6R4( data, m )
   real(SiKi),                   intent(in   )  :: data(0:,0:,0:,:)   !< 3D grid of vector data
   type(GridInterp_MiscVarType), intent(in   )  :: m                  !< MiscVars

   character(*),   parameter                    :: RoutineName = 'GridInterp3DVec6R4'
   integer(IntKi), parameter                    :: vDim = 6
   integer(IntKi)                               :: i,j,k,vi
   real(SiKi)                                   :: GridInterp3DVec6R4(vDim)

   ! interpolate
   GridInterp3DVec6R4 = 0.0_SiKi
   do k = 1,4
      do j = 1,4
         do i = 1,4
            do vi = 1,vDim
               GridInterp3DVec6R4(vi) = GridInterp3DVec6R4(vi) + m%N3D(i,j,k) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3), vi )
            end do
         end do
      end do
   end do

end function GridInterp3DVec6R4

function GridInterp3DVec6R8( data, m )
   real(DbKi),                   intent(in   )  :: data(0:,0:,0:,:)   !< 3D grid of vector data
   type(GridInterp_MiscVarType), intent(in   )  :: m                  !< MiscVars

   character(*),   parameter                    :: RoutineName = 'GridInterp3DVec6R8'
   integer(IntKi), parameter                    :: vDim = 6
   integer(IntKi)                               :: i,j,k,vi
   real(DbKi)                                   :: GridInterp3DVec6R8(vDim)

   ! interpolate
   GridInterp3DVec6R8 = 0.0_DbKi
   do k = 1,4
      do j = 1,4
         do i = 1,4
            do vi = 1,vDim
               GridInterp3DVec6R8(vi) = GridInterp3DVec6R8(vi) + m%N3D(i,j,k) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3), vi )
            end do
         end do
      end do
   end do

end function GridInterp3DVec6R8

!=============================================================================================================
! INTERFACE GridInterp4D
!              - GridInterp4DR4
!              - GridInterp4DR8
!=============================================================================================================
function GridInterp4DR4( data, m )
   real(SiKi),                   intent(in   )  :: data(0:,0:,0:,0:)   !< 4D grid of scalar data
   type(GridInterp_MiscVarType), intent(in   )  :: m                   !< MiscVars

   character(*), parameter                      :: RoutineName = 'GridInterp4DR4'
   real(SiKi)                                   :: GridInterp4DR4
   integer(IntKi)                               :: i,j,k,l

   ! interpolate
   GridInterp4DR4 = 0.0_SiKi
   do l = 1,4
      do k = 1,4
         do j = 1,4
            do i = 1,4
               GridInterp4DR4 = GridInterp4DR4 + m%N4D(i,j,k,l) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3), m%Indx(l,4) )
            end do
         end do
      end do
   end do

end function GridInterp4DR4

function GridInterp4DR8( data, m )
   real(DbKi),                   intent(in   )  :: data(0:,0:,0:,0:)   !< 4D grid of scalar data
   type(GridInterp_MiscVarType), intent(in   )  :: m                   !< MiscVars

   character(*), parameter                      :: RoutineName = 'GridInterp4DR8'
   real(DbKi)                                   :: GridInterp4DR8
   integer(IntKi)                               :: i,j,k,l

   ! interpolate
   GridInterp4DR8 = 0.0_DbKi
   do l = 1,4
      do k = 1,4
         do j = 1,4
            do i = 1,4
               GridInterp4DR8 = GridInterp4DR8 + m%N4D(i,j,k,l) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3), m%Indx(l,4) )
            end do
         end do
      end do
   end do

end function GridInterp4DR8

!=============================================================================================================
! INTERFACE GridInterp4DVec
!              - GridInterp4DVecR4
!              - GridInterp4DVecR8
!=============================================================================================================
function GridInterp4DVecR4( data, m )
   real(SiKi),                   intent(in   )  :: data(0:,0:,0:,0:,:)   !< 4D grid of vector data
   type(GridInterp_MiscVarType), intent(in   )  :: m                     !< MiscVars

   character(*),   parameter                    :: RoutineName = 'GridInterp4DVecR4'
   integer(IntKi), parameter                    :: vDim = 3
   integer(IntKi)                               :: i,j,k,l,vi
   real(SiKi)                                   :: GridInterp4DVecR4(vDim)

   ! interpolate
   GridInterp4DVecR4 = 0.0_SiKi
   do l = 1,4
      do k = 1,4
         do j = 1,4
            do i = 1,4
               do vi = 1,vDim
                  GridInterp4DVecR4(vi) = GridInterp4DVecR4(vi) + m%N4D(i,j,k,l) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3), m%Indx(l,4), vi )
               end do
            end do
         end do
      end do
   end do

end function GridInterp4DVecR4

function GridInterp4DVecR8( data, m )
   real(DbKi),                   intent(in   )  :: data(0:,0:,0:,0:,:)   !< 4D grid of vector data
   type(GridInterp_MiscVarType), intent(in   )  :: m                     !< MiscVars

   character(*),   parameter                    :: RoutineName = 'GridInterp4DVecR8'
   integer(IntKi), parameter                    :: vDim = 3
   integer(IntKi)                               :: i,j,k,l,vi
   real(DbKi)                                   :: GridInterp4DVecR8(vDim)

   ! interpolate
   GridInterp4DVecR8 = 0.0_DbKi
   do l = 1,4
      do k = 1,4
         do j = 1,4
            do i = 1,4
               do vi = 1,vDim
                  GridInterp4DVecR8(vi) = GridInterp4DVecR8(vi) + m%N4D(i,j,k,l) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3), m%Indx(l,4), vi )
               end do
            end do
         end do
      end do
   end do

end function GridInterp4DVecR8

!=============================================================================================================
! INTERFACE GridInterp4DVec6
!              - GridInterp4DVec6R4
!              - GridInterp4DVec6R8
!=============================================================================================================
function GridInterp4DVec6R4( data, m )
   real(SiKi),                   intent(in   )  :: data(0:,0:,0:,0:,:)   !< 4D grid of vector data
   type(GridInterp_MiscVarType), intent(in   )  :: m                     !< MiscVars

   character(*),   parameter                    :: RoutineName = 'GridInterp4DVec6R4'
   integer(IntKi), parameter                    :: vDim = 6
   integer(IntKi)                               :: i,j,k,l,vi
   real(SiKi)                                   :: GridInterp4DVec6R4(vDim)

   ! interpolate
   GridInterp4DVec6R4 = 0.0_SiKi
   do l = 1,4
      do k = 1,4
         do j = 1,4
            do i = 1,4
               do vi = 1,vDim
                  GridInterp4DVec6R4(vi) = GridInterp4DVec6R4(vi) + m%N4D(i,j,k,l) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3), m%Indx(l,4), vi )
               end do
            end do
         end do
      end do
   end do

end function GridInterp4DVec6R4

function GridInterp4DVec6R8( data, m )
   real(DbKi),                   intent(in   )  :: data(0:,0:,0:,0:,:)   !< 4D grid of vector data
   type(GridInterp_MiscVarType), intent(in   )  :: m                     !< MiscVars

   character(*),   parameter                    :: RoutineName = 'GridInterp4DVec6R8'
   integer(IntKi), parameter                    :: vDim = 6
   integer(IntKi)                               :: i,j,k,l,vi
   real(DbKi)                                   :: GridInterp4DVec6R8(vDim)

   ! interpolate
   GridInterp4DVec6R8 = 0.0_DbKi
   do l = 1,4
      do k = 1,4
         do j = 1,4
            do i = 1,4
               do vi = 1,vDim
                  GridInterp4DVec6R8(vi) = GridInterp4DVec6R8(vi) + m%N4D(i,j,k,l) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3), m%Indx(l,4), vi )
               end do
            end do
         end do
      end do
   end do

end function GridInterp4DVec6R8

!=============================================================================================================
! INTERFACE GridInterp4DVec6
!              - GridInterp4DVec6R4
!              - GridInterp4DVec6R8
!=============================================================================================================
function GridInterp4DVecNR4( vDim, data, m )
   integer(IntKi),               intent(in   )  :: vDim
   real(SiKi),                   intent(in   )  :: data(0:,0:,0:,0:,:)   !< 4D grid of vector data
   type(GridInterp_MiscVarType), intent(in   )  :: m                     !< MiscVars

   character(*),   parameter                    :: RoutineName = 'GridInterp4DVecNR4'
   integer(IntKi)                               :: i,j,k,l,vi
   real(SiKi)                                   :: GridInterp4DVecNR4(vDim)

   ! interpolate
   GridInterp4DVecNR4 = 0.0_SiKi
   do l = 1,4
      do k = 1,4
         do j = 1,4
            do i = 1,4
               do vi = 1,vDim
                  GridInterp4DVecNR4(vi) = GridInterp4DVecNR4(vi) + m%N4D(i,j,k,l) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3), m%Indx(l,4), vi )
               end do
            end do
         end do
      end do
   end do

end function GridInterp4DVecNR4

function GridInterp4DVecNR8( vDim, data, m )
   integer(IntKi),               intent(in   )  :: vDim
   real(DbKi),                   intent(in   )  :: data(0:,0:,0:,0:,:)   !< 4D grid of vector data
   type(GridInterp_MiscVarType), intent(in   )  :: m                     !< MiscVars

   character(*),   parameter                    :: RoutineName = 'GridInterp4DVecNR8'
   integer(IntKi)                               :: i,j,k,l,vi
   real(DbKi)                                   :: GridInterp4DVecNR8(vDim)

   ! interpolate
   GridInterp4DVecNR8 = 0.0_DbKi
   do l = 1,4
      do k = 1,4
         do j = 1,4
            do i = 1,4
               do vi = 1,vDim
                  GridInterp4DVecNR8(vi) = GridInterp4DVecNR8(vi) + m%N4D(i,j,k,l) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3), m%Indx(l,4), vi )
               end do
            end do
         end do
      end do
   end do

end function GridInterp4DVecNR8

!=============================================================================================================
! INTERFACE GridInterpN
!              - GridInterpNR4
!              - GridInterpNR8
!=============================================================================================================
function GridInterpNR4( data, p, m )
   real(SiKi),                     intent(in   )  :: data(0:,0:,0:)   !< 3D grid of scalar data
   type(GridInterp_ParameterType), intent(in   )  :: p                !< Parameters
   type(GridInterp_MiscVarType),   intent(in   )  :: m                !< MiscVars

   character(*), parameter                        :: RoutineName = 'GridInterpNR4'
   real(SiKi)                                     :: GridInterpNR4(3)
   real(SiKi)                                     :: dZetadx, dZetady
   integer(IntKi)                                 :: i,j,k

   ! interpolate slope
   dZetadx = 0.0_SiKi
   dZetady = 0.0_SiKi
   do k = 1,4
      do j = 1,4
         do i = 1,4
            dZetadx = dZetadx + m%N4D(i,j,k,1) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3) )
            dZetady = dZetady + m%N4D(i,j,k,2) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3) )
         end do
      end do
   end do
   dZetadx = dZetadx / p%delta(2)
   dZetady = dZetady / p%delta(3)

   GridInterpNR4 = [-dZetadx,-dZetady,1.0_SiKi]
   GridInterpNR4 = GridInterpNR4 / TwoNorm(GridInterpNR4)

end function GridInterpNR4

function GridInterpNR8( data, p, m )
   real(DbKi),                     intent(in   )  :: data(0:,0:,0:)   !< 3D grid of scalar data
   type(GridInterp_ParameterType), intent(in   )  :: p                !< Parameters
   type(GridInterp_MiscVarType),   intent(in   )  :: m                !< MiscVars

   character(*), parameter                        :: RoutineName = 'GridInterpNR8'
   real(DbKi)                                     :: GridInterpNR8(3)
   real(DbKi)                                     :: dZetadx, dZetady
   integer(IntKi)                                 :: i,j,k

   ! interpolate slope
   dZetadx = 0.0_DbKi
   dZetady = 0.0_DbKi
   do k = 1,4
      do j = 1,4
         do i = 1,4
            dZetadx = dZetadx + m%N4D(i,j,k,1) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3) )
            dZetady = dZetady + m%N4D(i,j,k,2) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3) )
         end do
      end do
   end do
   dZetadx = dZetadx / p%delta(2)
   dZetady = dZetady / p%delta(3)

   GridInterpNR8 = (/-dZetadx,-dZetady,1.0_DbKi/)
   GridInterpNR8 = GridInterpNR8 / TwoNorm(GridInterpNR8)

end function GridInterpNR8

!=============================================================================================================
! INTERFACE GridInterpS
!              - GridInterpSR4
!              - GridInterpSR8
!=============================================================================================================
function GridInterpSR4( data, p, m )
   real(SiKi),                     intent(in   )  :: data(0:,0:,0:)   !< 3D grid of scalar data
   type(GridInterp_ParameterType), intent(in   )  :: p                !< Parameters
   type(GridInterp_MiscVarType),   intent(in   )  :: m                !< MiscVars

   character(*), parameter                        :: RoutineName = 'GridInterpSR4'
   real(SiKi)                                     :: GridInterpSR4(2)
   integer(IntKi)                                 :: i,j,k,dir

   ! interpolate slope
   GridInterpSR4 = 0.0_SiKi
   do k = 1,4
      do j = 1,4
         do i = 1,4
            do dir = 1,2
               GridInterpSR4(dir) = GridInterpSR4(dir) + m%N4D(i,j,k,dir) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3) )
            end do
         end do
      end do
   end do
   GridInterpSR4 = GridInterpSR4 / p%delta(2:3)

end function GridInterpSR4

function GridInterpSR8( data, p, m )
   real(DbKi),                     intent(in   )  :: data(0:,0:,0:)   !< 3D grid of scalar data
   type(GridInterp_ParameterType), intent(in   )  :: p                !< Parameters
   type(GridInterp_MiscVarType),   intent(in   )  :: m                !< MiscVars

   character(*), parameter                        :: RoutineName = 'GridInterpSR8'
   real(DbKi)                                     :: GridInterpSR8(2)
   integer(IntKi)                                 :: i,j,k,dir

   ! interpolate slope
   GridInterpSR8 = 0.0_DbKi
   do k = 1,4
      do j = 1,4
         do i = 1,4
            do dir = 1,2
               GridInterpSR8(dir) = GridInterpSR8(dir) + m%N4D(i,j,k,dir) * data( m%Indx(i,1), m%Indx(j,2), m%Indx(k,3) )
            end do
         end do
      end do
   end do
   GridInterpSR8 = GridInterpSR8 / p%delta(2:3)

end function GridInterpSR8

END MODULE GridInterp