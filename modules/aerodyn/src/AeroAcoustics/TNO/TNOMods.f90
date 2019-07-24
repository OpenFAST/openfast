MODULE TNOConstants
    use ISO_FORTRAN_ENV
    REAL (kind=4),PARAMETER :: Cnuk = 5.5
    REAL (kind=4),PARAMETER :: kappa = 0.41
    REAL (kind=4),PARAMETER :: Cmu = 0.09
    REAL (kind=4),PARAMETER :: pi = 3.1415
    INTEGER (4),PARAMETER :: limit = 5000
    INTEGER (4)             :: i_omega
    !!REAL (ReKi),ALLOCATABLE :: omega(:)
    REAL (kind=4) :: omega
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
    REAL (kind=4)           :: d_star(2)
    REAL (kind=4)           :: edgevel(2)
END MODULE BLParams
!===========================================================
MODULE AirfoilParams
    use ISO_FORTRAN_ENV
    CHARACTER*128         :: airfoil
    REAL(kind=4)          :: aofa,a_chord,Mach,Re
    REAL(kind=4)          :: xtrup, xtrlo
    LOGICAL               :: ISTRIPPED   
    LOGICAL               :: ISNACA   
    LOGICAL               :: ISSUCTION
END MODULE AirfoilParams
!===========================================================
MODULE Third_Octave_Bands
    use ISO_FORTRAN_ENV
    INTEGER (4),PARAMETER :: NumBands = 34
    REAL (kind=4),PARAMETER :: Third_Octave(NumBands) = (/10.,12.5,16.,20.,25.,31.5,40.,50.,63.,80.,	&
                                                        100.,125.,160.,200.,250.,315.,400.,500.,630.,800.,	& 
                                                        1000.,1250.,1600.,2000.,2500.,3150.,4000.,5000.,6300.,8000.,	&
                                                        10000.,12500.,16000.,20000./)
END MODULE Third_Octave_Bands
