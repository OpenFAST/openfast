! A module for performing one-dimensional root finding via 
! 
!  brent's method (see sub_brent)
! 
! The function interface must be like
!     interface
!         function f(x, fcnArgs)
!         implicit none
!         real(ReKi), intent(in) :: x
!         real(ReKi) :: f
!         end function f
!     end interface
!
! Grey Gordon, 2013
! This code may be reproduced and modified provided it is not sold and the author's 
! name and this notice remain intact.
! 
module mod_root1dim
    use NWTC_Library
    implicit none
    integer, parameter, private :: SolveKi = ReKi
    real(SolveKi), parameter, private :: xtoler_def = 1d-6, toler_def = 1d-6, printmod_def = -1, maxiter_def = 100
contains

! Tests whether the two functions values bracket a root -- i.e., have different signs
function bracketsRoot(fa,fb) result(tf)
    logical :: tf
    real(SolveKi) :: fa,fb
    tf = sign(1.0_SolveKi,fb)/=sign(1.0_SolveKi,fa)
end function bracketsRoot



! Written using the FORTRAN program zero in "Algorithms for Minimization 
! without Derivatives" by Brent. 
! 
! Description: 
!       Returns a zero x of the function f in the given interval [a,b], to within a tolerance 6macheps|x| + 2t, 
! where macheps is the relative machine precision and t is a positive tolerance. The procedure assumes that 
! f(a) and f(b) have different signs.
subroutine sub_brent(x,f,a_in,b_in, toler_in,maxiter_in,fcnArgs,AFInfo,fa_in,fb_in,xtoler_in,printmod_in)
    use fminfcn
    
    implicit none 
    
    real(ReKi), intent(out) :: x !< solution
    interface
        function f(x, fcnArgs, AFInfo)
        use fminfcn
        implicit none
        real(ReKi),         intent(in)           :: x
        type(fmin_fcnArgs), intent(inout)        :: fcnArgs
        TYPE (AFInfoType),  INTENT(IN   )        :: AFInfo          ! The derived type for holding the constant parameters for this airfoil.
        real(ReKi) :: f
        end function f
    end interface
    
    real(ReKi), intent(in) :: a_in  !< lower bound of solution region
    real(ReKi), intent(in) :: b_in  !< upper bound of solution region
    
    type(fmin_fcnArgs), intent(inout) :: fcnArgs !< function arguments
    TYPE (AFInfoType),  INTENT(IN   ) :: AFInfo  !< The derived type for holding the constant parameters for this airfoil.
    real(ReKi), intent(in),  optional :: toler_in !< induction tolerance
    real(ReKi), intent(in),  optional :: fa_in !< starting value for f(a), if not present, will be evaluated
    real(ReKi), intent(in),  optional :: fb_in !< starting value for f(b), if not present, will be evaluated
    real(ReKi), intent(in),  optional :: xtoler_in !< 
    integer,    intent(in),  optional :: printmod_in !< print switch; otherwise uses default printmod_deff
    integer,    intent(in),  optional :: maxiter_in !< maximum number of iterations; otherwise uses default maxiter_def
    ! local
    real(SolveKi), parameter :: machep = epsilon(0.0_SolveKi)
    real(SolveKi) :: c,fa,fb,fc,toler,xtoler,e,d,m,p,q,tol,t,r,s
    real(ReKi) :: a,b
    integer :: maxiter,printmod,iter
    character(len=6) :: step

    ! Set of get parameters
    toler = 0.0_SolveKi; if (present(toler_in)) toler = toler_in ! Better to use custom toler here
    xtoler = xtoler_def; if (present(xtoler_in)) xtoler = xtoler_in
    maxiter = maxiter_def; if (present(maxiter_in)) maxiter = maxiter_in    
    printmod = printmod_def; if (present(printmod_in)) printmod = printmod_in

    ! Set the user chosen tolerance t to xtoler
    if (xtoler<0.0_SolveKi) then
        CALL WrScr('WARNING: xtoler must be positive. Resetting xtoler.')
        xtoler = 0.0_SolveKi
    end if
    t = xtoler
    
    ! Get initial bracket
    a=a_in
    b=b_in
    if (present(fa_in)) then
        fa = fa_in
    else
        fa = f(a, fcnArgs, AFInfo)
    end if
    if (present(fb_in)) then
        fb = fb_in
    else
        fb = f(b, fcnArgs, AFInfo)
    end if

    ! Test whether root is bracketed
    if (.not. bracketsRoot(fa,fb)) then
        if (abs(fa)<abs(fb)) then
            call WrScr( 'brent: WARNING: root is not bracketed, returning best endpoint a = '//trim(Num2Lstr(a))//' fa = '//trim(Num2Lstr(fa)) )
            x = a
        else
            call WrScr( 'brent: WARNING: root is not bracketed, returning best endpoint b = '//trim(Num2Lstr(b))//' fb = '//trim(Num2Lstr(fb)) )
            x = b
        end if
        return
    end if

    step = 'init'

    ! At any point in time, b is the best guess of the root, a is the previous value of b, 
    ! and the root is bracketed by b and c.
    do iter = 1,maxiter

        if (iter==1 .or. (fb>0.0_SolveKi .and. fc>0.0_SolveKi) .or. (fb<=0.0_SolveKi .and. fc<=0.0_SolveKi)) then
            c = a
            fc = fa
            e = b - a
            d = e
        end if

        ! If c is strictly better than b, swap b and c so b is the best guess. 
        if (abs(fc)<abs(fb)) then
            a = b
            b = c
            c  = a
            fa = fb
            fb = fc
            fc = fa
        end if

        ! Set the tolerance. Note: brent is very careful with these things, so don't deviate from this.
        tol = 2.0_SolveKi*machep*abs(b) + t

        ! Determine what half the length of the bracket [b,c] is
        m = 0.5_SolveKi*(c-b)

        ! If taking a bisection step would move the guess of the root less than tol, then return b the best guess.
        if ((abs(m)<=tol) .or. (fb==0.0_SolveKi)) then
            x = b
            return
        end if

        ! Display info
        if (printmod>0 .and. mod(iter,printmod)==0) write(*,"(A,i5,A,e13.5,A,e13.5,A,A)") &
          'brent: iter ',iter,' b',b,' fb',fb,' ',step

        ! If still here, then check whether need to do bisection or can do interpolation
        if ((abs(e)>=tol) .and. (abs(fa)>abs(fb))) then
            s = fb/fa
            if (a/=c) then
                ! Inverse quadratic interpolation
                q = fa/fc
                r = fb/fc
                p = s*(2.0_SolveKi*m*q*(q-r) - (b-a)*(r-1.0_SolveKi))
                q = (q-1.0_SolveKi)*(r-1.0_SolveKi)*(s-1.0_SolveKi)
                
                step = 'quad'
            else
                ! Linear interpolation
                p = 2.0_SolveKi*m*s
                q = 1.0_SolveKi-s

                step = 'linear'
            end if

            ! Ensure p is positive
            if (p<=0.0_SolveKi) then
                p = -p
            else
                q = -q
            end if

            s = e
            e = d
            if ((2.0_SolveKi*p>=3.0_SolveKi*m*q-abs(tol*q)) .or. &
                (p>=abs(0.5_SolveKi*s*q))) then
                ! Interpolation step failed to produce good step, bisect instead
                e = m
                d = m ! m is half the distance between b and c
                step = 'bisect'
            else
                !  Do interpolation step (either quadratic or linear)
                d = p/q
            end if
        else

            ! Do bisection step
            e = m 
            d = m
            
        end if

        ! Get new points. 
        !! Replace a (the old b) with b.
        a = b
        fa = fb

        !!! Increment b by d if that is greater than the tolerance. O/w, increment by tol.
        if (abs(d)<=tol) then
            ! m is .5*(c-b) with the bracket either [b,c] or [c,b]. 
            if (m > 0.0_SolveKi) then
                ! If m>0d0, then bracket is [b,c] so move towards c by tol
                b = b + tol
            else
                ! If m<=0d0, then bracket is [c,b] so move towards c by tol
                b = b - tol
            end if
        else
            b = b + d
        end if

        !!! Evaluate at the new point
        fb = f(b, fcnArgs, AFInfo)

        ! Check my custom tolerance 
        if (abs(fb)<toler) then
            x = b
            return
        end if
            
    end do

end subroutine sub_brent

end module mod_root1dim

! Contains test routines for mod_root1dim
module mod_root1dim_test
    use mod_root1dim
    implicit none
contains
    
    ! One dimensional test functions
    function f0(x)
        implicit none
        real(ReKi) :: f0
        real(ReKi), intent(in) :: x
        f0 = sin(x) + cos(x)
    end function 
    function f1(x)
        real(ReKi) :: f1
        real(ReKi), intent(in) :: x
        f1 = sign(sqrt(abs(x-5d0)),x-5d0)
    end function 
    function f2(x) 
        real(ReKi) :: f2
        real(ReKi), intent(in) :: x
        f2 = -10d0 + 2d0*x + x**2 
    end function 
    function f3(x)
        real(ReKi) :: f3
        real(ReKi), intent(in) :: x
        f3 = exp(x) - 2d0 
    end function 

!   ! type
!
!    ! Main routine
!    subroutine test_root1dim() 
!        implicit none
!        real(ReKi) :: a,b,x
!        logical :: tf(3)
!        
!        a = -10; b = 20; if (.not. all(passed(f0))) STOP 'test_root1dim: f0 test failed'
!        a = -100; b = 1000; if (.not. all(passed(f1))) STOP 'test_root1dim: f1 test failed'
!        a = -6; b = 2; if (.not. all(passed(f2))) STOP 'test_root1dim: f2 test failed' 
!        a = -100; b = 100; if (.not. all(passed(f3))) STOP 'test_root1dim: f3 test failed' 
!
!        print*, 'test_root1dim: passed'
!
!    contains 
!        function passed(f) result(tfout)
!            implicit none
!            logical :: tf(3), tfout(3)
!            real(ReKi) :: xstar(3)
!           
!            interface
!                  function f(x, fcnArgs)
!                 use fminfcn
!                 implicit none
!                 real(ReKi),            intent(in)           :: x
!                 type(fmin_fcnArgs), intent(in)              :: fcnArgs
!                 real(ReKi) :: f
!                 end function f
!            end interface
!            
!            tf = .true.
!            
!            call sub_bisect( xstar(1),f,a,b,xtoler_in=2d0*epsilon(0d0),printmod_in=1) 
!            call sub_brent(  xstar(2),f,a,b,xtoler_in=2d0*epsilon(0d0),printmod_in=1) 
!            
!            if (abs(f(xstar(1)))>1d-4) tf(1) = .false.
!            if (abs(f(xstar(2)))>1d-3) tf(2) = .false.
!            if (abs(f(xstar(3)))>1d-4) tf(3) = .false.
!
!!             print*,'xstar',xstar
!!             print*,'f(xstar)',f(xstar(1)),f(xstar(2)),f(xstar(3))
!!             print*,'tf',tf
!
!            tfout = tf
!
!        end function passed
!        
!    end subroutine test_root1dim

end module mod_root1dim_test
