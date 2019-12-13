MODULE TNOConstants
    use ISO_FORTRAN_ENV
    REAL (kind=4),PARAMETER :: Cnuk = 5.5
    REAL (kind=4),PARAMETER :: kappa = 0.41
    REAL (kind=4),PARAMETER :: Cmu = 0.09
    REAL (kind=4),PARAMETER :: pi = 3.1415
    INTEGER (4),PARAMETER :: limit = 5000
    REAL (kind=4) :: omega ! NOTE: not a constant and used by function int1 and int2
END MODULE TNOConstants
!===========================================================
MODULE Atmosphere
    use ISO_FORTRAN_ENV
    !atmosphere constants
    REAL (kind=4) nu
    REAL (kind=4) co
    REAL (kind=4) rho
END MODULE Atmosphere
!===========================================================
MODULE Wavenumber
    use ISO_FORTRAN_ENV
    REAL (kind=4)           :: k
    REAL (kind=4)           :: k1
    REAL (kind=4)           :: k3
END MODULE Wavenumber
!===========================================================
MODULE BLParams
    use ISO_FORTRAN_ENV
    REAL (kind=4)           :: d99(2)
    REAL (kind=4)           :: Cf(2)
    REAL (kind=4)           :: edgevel(2)
END MODULE BLParams
!===========================================================
MODULE AirfoilParams
    use ISO_FORTRAN_ENV
    REAL(kind=4)          :: Mach
    LOGICAL               :: ISSUCTION
END MODULE AirfoilParams
!===========================================================
