*DECK AAAAAA
      SUBROUTINE AAAAAA (VER)
C***BEGIN PROLOGUE  AAAAAA
C***PURPOSE  SLATEC Common Mathematical Library disclaimer and version.
C***LIBRARY   SLATEC
C***CATEGORY  Z
C***TYPE      ALL (AAAAAA-A)
C***KEYWORDS  DISCLAIMER, DOCUMENTATION, VERSION
C***AUTHOR  SLATEC Common Mathematical Library Committee
C***DESCRIPTION
C
C   The SLATEC Common Mathematical Library is issued by the following
C
C           Air Force Weapons Laboratory, Albuquerque
C           Lawrence Livermore National Laboratory, Livermore
C           Los Alamos National Laboratory, Los Alamos
C           National Institute of Standards and Technology, Washington
C           National Energy Research Supercomputer Center, Livermore
C           Oak Ridge National Laboratory, Oak Ridge
C           Sandia National Laboratories, Albuquerque
C           Sandia National Laboratories, Livermore
C
C   All questions concerning the distribution of the library should be
C   directed to the NATIONAL ENERGY SOFTWARE CENTER, 9700 Cass Ave.,
C   Argonne, Illinois  60439, and not to the authors of the subprograms.
C
C                    * * * * * Notice * * * * *
C
C   This material was prepared as an account of work sponsored by the
C   United States Government.  Neither the United States, nor the
C   Department of Energy, nor the Department of Defense, nor any of
C   their employees, nor any of their contractors, subcontractors, or
C   their employees, makes any warranty, expressed or implied, or
C   assumes any legal liability or responsibility for the accuracy,
C   completeness, or usefulness of any information, apparatus, product,
C   or process disclosed, or represents that its use would not infringe
C   upon privately owned rights.
C
C *Usage:
C
C        CHARACTER * 16 VER
C
C        CALL AAAAAA (VER)
C
C *Arguments:
C
C     VER:OUT   will contain the version number of the SLATEC CML.
C
C *Description:
C
C   This routine contains the SLATEC Common Mathematical Library
C   disclaimer and can be used to return the library version number.
C
C***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
C                 and Lee Walton, Guide to the SLATEC Common Mathema-
C                 tical Library, April 10, 1990.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800424  DATE WRITTEN
C   890414  REVISION DATE from Version 3.2
C   890713  Routine modified to return version number.  (WRB)
C   900330  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   921215  Updated for Version 4.0.  (WRB)
C   930701  Updated for Version 4.1.  (WRB)
C***END PROLOGUE  AAAAAA
      CHARACTER * (*) VER
C***FIRST EXECUTABLE STATEMENT  AAAAAA
      VER = ' 4.1'
      RETURN
      END
