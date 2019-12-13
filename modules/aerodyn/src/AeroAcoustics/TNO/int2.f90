FUNCTION int2(k1)

USE Atmosphere, only: co
USE TNOConstants, only: omega ! NOTE: Not a constant, set it as an argument!

implicit none

REAL (kind=4) :: int2
REAL (kind=4) :: k1
REAL (kind=4), EXTERNAL :: Pressure

int2 = omega/co/k1*Pressure(k1)

RETURN 

END FUNCTION int2
