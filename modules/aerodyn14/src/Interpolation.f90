!> 
module Interpolation
    use NWTC_Library
    implicit none

    interface interp_lin; module procedure interp_lin11, interp_lin10
    end interface

    private

    public:: interp_lin
    public:: interp_lin0

contains

    !--------------------------------------------------------------------
    ! --- Linear Interp
    !--------------------------------------------------------------------
    function interp_lin0(x,x0,x1,f0,f1)   ! Linear interpolation function                                     
        real(ReKi) ::interp_lin0
        real(ReKi),intent(in):: x,x0,x1,f0,f1
        real(ReKi) :: eps=1.0e-7
        if (abs(x0-x1)<eps) then    ! changed by tjul 16.08.2006 to avoid division by zero
            interp_lin0=f0
        else
            interp_lin0=(x-x1)/(x0-x1)*f0+(x-x0)/(x1-x0)*f1
        endif
    end function interp_lin0


    !> Linear interpolation between one dim arrays. ASSUMES xknown to have increasing values
    subroutine interp_lin11_incr(xknown,yknown,xnew,ynew)   
        real(ReKi),dimension(:),intent(in):: xknown,yknown
        real(ReKi),dimension(:),intent(in):: xnew
        real(ReKi),dimension(:),intent(out):: ynew
        integer i
        integer itmp
        integer:: nknown
        nknown=size(xknown)
        do i=1,size(xnew)
            itmp=minloc(abs(xnew(i)-xknown),dim=1)
            if (itmp==nknown) then
                if (xknown(itmp)>xnew(i)) then
                    ynew(i)=interp_lin0(xnew(i),xknown(itmp-1),xknown(itmp),yknown(itmp-1),yknown(itmp))
                else
                    ! The current x is above the max of xknown
                    ! extrapolation required, here fixed to upper bound
                    ynew(i)=yknown(nknown)
                endif
            elseif (xknown(itmp)<xnew(i)) then
                ! normal case, x between itmp and itmp+1
                ynew(i)=interp_lin0(xnew(i),xknown(itmp),xknown(itmp+1),yknown(itmp),yknown(itmp+1))
            elseif (itmp==1) then
                ! The current x is below the min of xknown
                ynew(i)=yknown(1)
            else
                ! normal case but inverted, x between itmp-1 and itmp
                ynew(i)=interp_lin0(xnew(i),xknown(itmp-1),xknown(itmp),yknown(itmp-1),yknown(itmp))
            endif
        enddo
    end subroutine

    !>
    subroutine interp_lin11(xknown,yknown,xnew,ynew)   ! Linear interpolation between one dim arrays
        real(ReKi),dimension(:),intent(in):: xknown,yknown
        real(ReKi),dimension(:),intent(in):: xnew
        real(ReKi),dimension(:),intent(out):: ynew
        !
        integer :: n
        n=size(xknown)
        if(n<2) then
            print*, 'interp_lin error, input argument should be at least of size 2'
        elseif(size(yknown)/=n) then
            print*, 'interp_lin error, dimension of two first inputs mismatch'
        elseif(size(xnew)/=size(ynew)) then
            print*, 'interp_lin error, dimension of two last inputs mismatch'
        else
            if(xknown(2)>xknown(1)) then
                ! increasing x values
                call interp_lin11_incr(xknown,yknown,xnew,ynew) 
            else
                ! decreasing x values
                call interp_lin11_incr(xknown(n:1:-1),yknown(n:1:-1),xnew,ynew) 
            endif
        endif

    end subroutine


    subroutine interp_lin10(xknown,yknown,xnew,ynew)   ! Linear interpolation between one dim array at scalar location
        real(ReKi),dimension(:),intent(in):: xknown,yknown
        real(ReKi),intent(in):: xnew
        real(ReKi),intent(out):: ynew
        integer itmp
        integer:: nknown
        nknown=size(xknown)
        if(size(xknown)/=size(yknown)) then
            print*, 'interp_lin error, dimension of two first inputs mismatch'
        else
            if(xknown(2)<xknown(1)) then
                print*,'interp_lin10, assumes increasing x values'
                print*,'x_in',xknown
                print*,'y_in',yknown
                print*,'x', xnew
                STOP
            endif
            itmp=minloc(abs(xnew-xknown),dim=1)
            if (itmp==nknown) then
                if (xknown(itmp)>xnew) then
                    ynew=interp_lin0(xnew,xknown(itmp-1),xknown(itmp),yknown(itmp-1),yknown(itmp))
                else
                    ! The current x is above the max of xknown
                    ! extrapolation required, here fixed to upper bound
                    ynew=yknown(nknown)
                endif
            elseif (xknown(itmp)<xnew) then
                ! normal case, x between itmp and itmp+1
                ynew=interp_lin0(xnew,xknown(itmp),xknown(itmp+1),yknown(itmp),yknown(itmp+1))
            elseif (itmp==1) then
                ! The current x is below the min of xknown
                ynew=yknown(1)
            else
                ! normal case but inverted, x between itmp-1 and itmp
                ynew=interp_lin0(xnew,xknown(itmp-1),xknown(itmp),yknown(itmp-1),yknown(itmp))
            endif
        endif
    end subroutine

end module Interpolation
