@test
subroutine test_BD_ComputeIniCoef()
    ! test branches
    ! - test the inputs/outputs from static_cantilever_beam
    ! - test randomly chosen integer position, no twist
    ! - test randomly chosen real-valued position, no twist
    ! - test randomly chosen real-valued position, with twist

    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------
    ! In BD_ComputeIniCoef(), the key point coordinates ([1 = x; 2 = y; 3 = z;
    ! 4 = -twist]) are interpolated via a cubic spline fit (the LaTeX below
    ! describes how this is done).
    ! This test verifies that the interpolation is occurring properly for a
    ! variety of key point coordinate inputs.
    ! --------------------------------------------------------------------------
    ! --------------------------------------------------------------------------

    use pFUnit_mod
    use BeamDyn
    use NWTC_Num
    use test_tools

    implicit none

    REAL(BDKi)      :: kp_coordinate(3, 4)
    INTEGER(IntKi)  :: kp_member
    REAL(BDKi)      :: SP_Coef(2, 4, 4)
    REAL(BDKi)      :: base_SP_Coef(2, 4, 4)

    integer(IntKi)  :: ErrStat ! Error status of the operation
    character(1024) :: ErrMsg  ! Error message if ErrStat /= ErrID_None

    character(1024) :: testname
    integer(IntKi)  :: accuracy
    real(BDKi)      :: tolerance

    ! initialize NWTC_Num constants
    call SetConstants()

    ! digits of desired accuracy
    accuracy = 16


    ! --------------------------------------------------------------------------
    testname = "test the inputs/outputs from static_cantilever_beam:"

    kp_member     = 3
    kp_coordinate = reshape((/ 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                               0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                               0.0000000000000000, 5.0000000000000000, 10.000000000000000,&
                               0.0000000000000000, 0.0000000000000000, 0.0000000000000000 /),&
                            (/ 3, 4 /))

    SP_Coef               = 0.0d0
    base_SP_Coef          = 0.0d0
    base_SP_Coef(1, 2, 3) = 1.0d0
    base_SP_Coef(2, 2, 3) = 1.0d0

    call BD_ComputeIniCoef(kp_member, kp_coordinate, SP_Coef, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, base_SP_Coef)
    @assertEqual(base_SP_Coef, SP_Coef, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "test randomly chosen integer position, no twist:"

    kp_member     = 3
    kp_coordinate = reshape((/ 2.0000000000000000, 3.0000000000000000, 4.0000000000000000,&
                               5.0000000000000000, 5.0000000000000000, 5.0000000000000000,&
                               0.0000000000000000, 5.0000000000000000, 10.000000000000000,&
                               0.0000000000000000, 0.0000000000000000, 0.0000000000000000 /),&
                            (/ 3, 4 /))

    SP_Coef               = 0.0d0
    base_SP_Coef          = 0.0d0
    base_SP_Coef(1, :, :) = reshape((/ 2.0000000000000000, 0.20000000000000001, 0.0000000000000000,&
                                       0.0000000000000000, 5.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                       1.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000 /),&
                                    (/ 4, 4 /))
    base_SP_Coef(2, :, :) = reshape((/ 2.0000000000000000, 0.200000000000000, 0.0000000000000000,&
                                       0.0000000000000000, 5.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                       1.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000 /),&
                                    (/ 4, 4 /))

    call BD_ComputeIniCoef(kp_member, kp_coordinate, SP_Coef, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, base_SP_Coef)
    @assertEqual(base_SP_Coef, SP_Coef, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "test randomly chosen real-valued position, no twist:"

    kp_member     = 3
    kp_coordinate = reshape((/ 6.5070816145726680, 9.9174285549003274, 3.0353426271406194,&
                               4.0676504339184724, 9.8600230042671247, 8.6152333976082973,&
                               1.2655438606634469, 1.9351441450244888, 6.2914376761408058,&
                               0.0000000000000000, 0.0000000000000000, 0.0000000000000000 /),&
                            (/ 3, 4 /))

    SP_Coef               = 0.0d0
    base_SP_Coef          = 0.0d0
    base_SP_Coef(1, :, :) = reshape((/ 1.5084746148222619, 0.77405984809861472, 3.7640457741166018,&
                                      -0.99141717936807128, -4.9422059009721213, 2.8665039656159514,&
                                       5.0407395052931054, -1.3276872923894700, 0.0000000000000000,&
                                       0.99999999999999989, 0.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                       0.0000000000000000 /),&
                                    (/ 4, 4 /))
    base_SP_Coef(2, :, :) = reshape((/ -6.7803428269260309, 13.623982366366723, -2.8762463980747857,&
                                        0.15238946147303528, -16.042434760593011, 20.074879199307517,&
                                       -3.8518152317462877, 0.20407710871097259, 0.0000000000000000,&
                                        0.99999999999999922, 0.0000000000000000, 0.0000000000000000,&
                                        0.0000000000000000, 0.0000000000000000, 0.0000000000000000,&
                                        0.0000000000000000 /),&
                                    (/ 4, 4 /))

    call BD_ComputeIniCoef(kp_member, kp_coordinate, SP_Coef, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, base_SP_Coef)
    @assertEqual(base_SP_Coef, SP_Coef, tolerance, testname)

    ! --------------------------------------------------------------------------
    testname = "test randomly chosen real-valued position, with twist:"

    kp_member     = 3
    kp_coordinate = reshape((/ 5.4865738082235779, 4.3652789706860627, 6.9069067544825052,&
                               1.0045718278258167, 6.8252852409960285, 8.1107452755498365,&
                               4.9035987218189865, 1.0004806979797973, 7.4265368279480926,&
                               -6.371698463301453, -8.748463941100916, -7.017037545817013 /),&
                            (/ 3, 4 /))

    SP_Coef               = 0.0d0
    base_SP_Coef          = 0.0d0
    base_SP_Coef(1, :, :) = reshape((/ 4.3153057788616724, -2.5435772421555584E-002, 8.0847145111863264E-002,&
                                      -5.4957695152967907E-003, 12.027690378189149, -6.3778859138083925,&
                                       1.2633336546839726, -8.5877993310141815E-002, 0.0000000000000000,&
                                       0.99999999999999944, 0.0000000000000000, 0.0000000000000000,&
                                     -10.102482767721879, 1.5898240522874485, -0.25358873848672225,&
                                       1.7238274233055084E-002 /),&
                                    (/ 4, 4 /))
    base_SP_Coef(2, :, :) = reshape((/ 4.3131449680939369, -1.8956454713281810E-002, 7.4370940502062083E-002,&
                                      -3.3380718462735709E-003, 11.993925117084293, -6.2766387997770119,&
                                       1.1621351865517835, -5.2161378108225374E-002, 0.0000000000000000,&
                                       0.99999999999999944, 0.0000000000000000, 0.0000000000000000,&
                                     -10.095705072901112, 1.5695007372016476, -0.23327518808356279,&
                                       1.0470344095679714E-002 /),&
                                    (/ 4, 4 /))

    call BD_ComputeIniCoef(kp_member, kp_coordinate, SP_Coef, ErrStat, ErrMsg)

    tolerance = AdjustTol(accuracy, base_SP_Coef)
    @assertEqual(base_SP_Coef, SP_Coef, tolerance, testname)

    ! --------------------------------------------------------------------------

end subroutine test_BD_ComputeIniCoef

! \noindent\LARGE{\textbf{Algorithm for \texttt{ComputeIniCoef()}}\normalsize\\

! \begin{itemize}
!    \item This subroutine computes the $z$-coordinate weights to interpolate the positions in $x$, $y$, $z$, and negative twist ($-t$)
!    \item We start with the \texttt{kp\textunderscore coordinate} matrix ($C$), of dimension \texttt{kp\textunderscore member} $\times$ 4, where the 4 columns are $(x, y, z, -t)$, with $-t$ being negative twist.
!    \item For a 3 key point example, we have 
!    \begin{equation}
!       C = 
!       \begin{pmatrix}
!          x_1 & y_1 & z_1 & -t_1\\
!          x_2 & y_2 & z_2 & -t_2\\
!          x_3 & y_3 & z_3 & -t_3
!       \end{pmatrix}
!    \end{equation}
!    \item From this, we build the coefficient matrix, $K$, based on the z-component of the key points, where $K$ has dimension $4 (\texttt{kp\textunderscore member} - 1) \times 4 (\texttt{kp\textunderscore member} - 1)$
!    \begin{equation}
!       K = 
!       \begin{pmatrix}
!          0 & 0 & 2 & 6z_1 & 0 & 0 & 0 & 0\\
!          1 & z_1 & z_1^2 & z_1^3 & 0 & 0 & 0 & 0\\
!          1 & z_2 & z_2^2 & z_2^3 & 0 & 0 & 0 & 0\\
!          0 & 1 & 2z_2 & 3z_2^2 & 0 & -1 & -2z_2 & -3z_2^2\\
!          0 & 0 & 2 & 6z_2 & 0 & 0 & -2 & -6z_2\\
!          0 & 0 & 0 & 0 & 1 & z_2 & z_2^2 & z_2^3\\
!          0 & 0 & 0 & 0 & 1 & z_3 & z_3^2 & z_3^3\\
!          0 & 0 & 0 & 0 & 0 & 0 & 2 & 6z_3
!       \end{pmatrix}
!    \end{equation}
!    \item Next, the coefficient matrices (for the first \texttt{kp\textunderscore member} - 1 key points), $S_j, \ j = 1, \dots, \texttt{kp\textunderscore member} - 1$, are solved for
!    \begin{itemize}
!       \item This is done sequentially for each column of $C$, i.e., each coordinate
!       \item The RHS, $b_i,\ i = x, y, z, -t$, is built for the relevant solve
!       \begin{equation}
!          b_i = 
!          \begin{pmatrix}
!             0\\
!             i_1\\
!             i_2\\
!             0\\
!             0\\
!             i_2\\
!             i_3\\
!             0
!          \end{pmatrix}
!       \end{equation}
!       \item The system is solved for $u_i$ in $K u_i = b_i,\ i = x, y, z, -t$
!       \item Thus the equations look like (using the $x$-interpolation for illustration)
!       \begin{align}
!          2 u_x^{(3)} + 6 z_1 u_x^{(4)} &= 0,\label{sys1}\\
!          u_x^{(1)} + z_1 u_x^{(2)} + z_1^2 u_x^{(3)} + z_1^3 u_x^{(4)} &= x_1,\label{sys2}\\
!          u_x^{(1)} + z_2 u_x^{(2)} + z_2^2 u_x^{(3)} + z_2^3 u_x^{(4)} &= x_2,\label{sys3}\\
!          u_x^{(2)} + 2 z_2 u_x^{(3)} + 3 z_2^2 u_x^{(4)} - u_x^{(6)} - 2 z_2 u_x^{(7)} - 3 z_2^2 u_x^{(8)} &= 0,\nonumber\\
!          \iff u_x^{(2)} + 2 z_2 u_x^{(3)} + 3 z_2^2 u_x^{(4)} &= u_x^{(6)} + 2 z_2 u_x^{(7)} + 3 z_2^2 u_x^{(8)},\label{sys4}\\
!          2 u_x^{(3)} + 6 z_2 u_x^{(4)} - 2 u_x^{(7)} - 6 z_2 u_x^{(8)} &= 0,\nonumber\\
!          \iff 2 u_x^{(3)} + 6 z_2 u_x^{(4)} &= 2 u_x^{(7)} + 6 z_2 u_x^{(8)},\label{sys5}\\
!          u_x^{(5)} + z_2 u_x^{(6)} + z_2^2 u_x^{(7)} + z_2^3 u_x^{(8)} &= x_2,\label{sys6}\\
!          u_x^{(5)} + z_3 u_x^{(6)} + z_3^2 u_x^{(7)} + z_3^3 u_x^{(8)} &= x_3,\label{sys7}\\
!          2 u_x^{(7)} + 6 z_3 u_x^{(8)} &= 0,\label{sys8}
!       \end{align}
!       where
!       \begin{itemize}
!          \item \eqref{sys1} imposes the zero second derivative at the first key point/boundary
!          \item \eqref{sys2} is equality at the first key point/boundary
!          \item \eqref{sys3} is equality at the second key point/interface
!          \item \eqref{sys4} is equality of first derivative at second key point/interface
!          \item \eqref{sys5} is equality of second derivative at second key point/interface
!          \item \eqref{sys6} is equality at the second key point/interface
!          \item \eqref{sys7} is equality at the third key point/boundary
!          \item \eqref{sys8} imposes the zero second derivative at the third key point/boundary
!       \end{itemize}
!       \item The coordinate-wise vectors are then stored in the relevant coefficient matrix, $S_j, \ j = 1, \dots, \text{\texttt{kp\textunderscore member} - 1}$, where the dimension of $S_j$ is $4 \times 4$ (order of spline coefficient $\times$ number of columns of \texttt{kp\textunderscore coordinate}, $x$, $y$, $z$, $-t$)
!       \begin{equation}
!          S_1 =
!          \begin{pmatrix}
!             u_x^{(1)} & u_y^{(1)} & u_z^{(1)} & u_{-t}^{(1)}\\
!             u_x^{(2)} & u_y^{(2)} & u_z^{(2)} & u_{-t}^{(2)}\\
!             u_x^{(3)} & u_y^{(3)} & u_z^{(3)} & u_{-t}^{(3)}\\
!             u_x^{(4)} & u_y^{(4)} & u_z^{(4)} & u_{-t}^{(4)}\\
!          \end{pmatrix}
!       \end{equation}
!       \begin{equation}
!          S_2 =
!          \begin{pmatrix}
!             u_x^{(5)} & u_y^{(5)} & u_z^{(5)} & u_{-t}^{(5)}\\
!             u_x^{(6)} & u_y^{(6)} & u_z^{(6)} & u_{-t}^{(6)}\\
!             u_x^{(7)} & u_y^{(7)} & u_z^{(7)} & u_{-t}^{(7)}\\
!             u_x^{(8)} & u_y^{(8)} & u_z^{(8)} & u_{-t}^{(8)}\\
!          \end{pmatrix}
!       \end{equation}
!    \end{itemize}
! \end{itemize}

