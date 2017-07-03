@test
subroutine test_BD_CrvExtractCrv()

use pFUnit_mod
use BeamDyn_Subs
implicit none

REAL(BDKi) :: expected(3)
REAL(BDKi) :: input(3,3)
REAL(BDKi) :: output(3)

integer :: i

expected = (/-3.9999E+000, 0.000000000000E+000, 0.000000000000E+000/)

input = reshape( (/ 1.0                     ,   0.00000000000000000E+000,   0.00000000000000000E+000, &
                    0.00000000000000000E+000,  -9.99999917034452235E-001,   4.07346398941426172E-004, &
                    0.00000000000000000E+000,  -4.07346398941426172E-004,  -9.99999917034452235E-001 /), &
                    shape(input), order=(/2,1/) )

call BD_CrvExtractCrv(input, output)

@assertEqual(expected, output)!, 0.001)

end subroutine
