C=======================================================================
C  FFTPACK_REALKIND
C
C  Returns the default REAL kind that this file -- and therefore
C  fftpack5.1.f -- was compiled with.
C
C  This file MUST be compiled with exactly the same real-size flags as
C  fftpack5.1.f.  See the FFTPACK_SOURCES block in
C  modules/nwtc-library/CMakeLists.txt and the RealKIND="realKIND4"
C  file configurations in vs-build/modules/NWTC-Library.vfproj.
C
C  FFTPACK 5.1 declares its arrays as bare REAL/COMPLEX, while the
C  NWTC_FFTPACK wrapper passes it explicitly kinded REAL(SiKi) and
C  COMPLEX(SiKi) buffers together with their element counts (LENSAV,
C  LENWRK).  A build that promotes the default REAL to 8 bytes in
C  fftpack5.1.f -- which DOUBLE_PRECISION does by default via
C  -fdefault-real-8 (GNU) or -real-size 64 (Intel) -- makes FFTPACK
C  write twice as many bytes as those buffers hold.  Nothing diagnoses
C  that at compile time, and the resulting heap/stack corruption
C  surfaces as a SIGBUS or SIGSEGV far away from the FFT call.
C
C  NWTC_FFTPACK calls this at initialization and aborts with a clear
C  message if the answer is not SiKi, so a build system that forgets
C  the flag fails loudly instead of corrupting memory.
C=======================================================================
      INTEGER FUNCTION FFTPACK_REALKIND ()
      FFTPACK_REALKIND = KIND(1.0)
      RETURN
      END
