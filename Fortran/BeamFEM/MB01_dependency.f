* *******************************************************************
* COPYRIGHT (c) 1967 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of an HSL ARCHIVE
* Licence, see http://hsl.rl.ac.uk/archive/cou.html
*
* Please note that for an HSL ARCHIVE Licence:
*
* 1. The Package must not be copied for use by any other person.
*    Supply of any part of the library by the Licensee to a third party
*    shall be subject to prior written agreement between AEA
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
*######DATE 4 Oct 1992
C       Toolpack tool decs employed.
C       SAVE statement for COMMON FA01ED added.
C  EAT 21/6/93 EXTERNAL statement put in for block data on VAXs.
C
C
      DOUBLE PRECISION FUNCTION FA01AD(I)
      INTEGER I
      DOUBLE PRECISION R,S
      INTRINSIC DINT,MOD
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
      EXTERNAL FA01FD
      SAVE /FA01ED/
      R = GR*9228907D0/65536D0
      S = DINT(R)
      GL = MOD(S+GL*9228907D0,65536D0)
      GR = R - S
      IF (I.GE.0) FA01AD = (GL+GR)/65536D0
      IF (I.LT.0) FA01AD = (GL+GR)/32768D0 - 1.D0
      GR = GR*65536D0
      RETURN
      END
      SUBROUTINE FA01BD(MAX,NRAND)
      INTEGER MAX,NRAND
      DOUBLE PRECISION FA01AD
      EXTERNAL FA01AD
      INTRINSIC DBLE,INT
      NRAND = INT(FA01AD(1)*DBLE(MAX)) + 1
      RETURN
      END
      SUBROUTINE FA01CD(IL,IR)
      INTEGER IL,IR
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
      SAVE /FA01ED/
      IL = GL
      IR = GR
      RETURN
      END
      SUBROUTINE FA01DD(IL,IR)
      INTEGER IL,IR
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
      SAVE /FA01ED/
      GL = IL
      GR = IR
      RETURN
      END
      BLOCK DATA FA01FD
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
      SAVE /FA01ED/
      DATA GL/21845D0/
      DATA GR/21845D0/
      END

