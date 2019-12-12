FUNCTION Pressure(k1_in)

    USE Atmosphere,    only: rho
    USE AirfoilParams, only: ISSUCTION
    USE BLParams,      only: d99
    use ISO_FORTRAN_ENV
    USE Wavenumber, only: k, k1, k3 ! NOTE: that's all of it

    implicit none

    ! Variables
    REAL (kind=4):: a,b,answer
    REAL (kind=4):: abserr,resabs,resasc
    
    REAL (kind=4):: k1_in,Pressure

    REAL (kind=4), EXTERNAL :: int1

    a = 0. !1e-4*d99(1)
    IF (ISSUCTION)THEN
        b = d99(1)
    ELSE
        b = d99(2)
    ENDIF
    

    k1 = k1_in
    k3 = 0.
    k= sqrt(k1**2+k3**2)

    CALL qk61(int1,a,b,answer,abserr,resabs,resasc)
               
    Pressure = 4.*rho**2*k1**2./(k1**2.+k3**2.)*answer
               
RETURN
               
END FUNCTION Pressure
