!**********************************************************************************************************************************
!
!  PROGRAM:    Compare_LAPACK_MUMPS.f90
!
!  PURPOSE:  Compares a solve of a system of equations by LAPACK and MUMPS.
! 
!     

!**********************************************************************************************************************************

PROGRAM Compare_LAPACK_MUMPS

    USE NWTC_Library
    USE NWTC_LAPACK
    USE MUMPS
    IMPLICIT NONE   
    
   

    ! Variables
    
   
      ! File name
   CHARACTER(*),PARAMETER :: InputFilePath='C:\Users\bjonkman\Documents\DATA\DesignCodes\simulators\FAST\SVNdirectory\branches\FAST_MUMPS\CertTest\Test21'
   CHARACTER(1024)        :: InputFileRoot
   CHARACTER(1024)        :: InputFileName 
   CHARACTER(10)          :: TmpLine
   
   
   INTEGER(IntKi)    :: ErrStat
   CHARACTER(1024)   :: ErrMsg
      
   integer           :: i, j, K, iTest, i_rhs, N, NZ
   integer           :: UnIn, UnOut, UnOut2
   REAL              :: time_begin,  time_end
   REAL              :: time_begin2, time_end2
   REAL              :: init_time, init_time2
   REAL              :: total_time, total_time2
   
   real,   allocatable :: full_jac(:,:), rhs(:)
   integer,allocatable :: ipiv(:)
   
   type(SMUMPS_STRUC) :: mumps_par
   
      
   !...............................................................................................................................    
   ! Initialize the library
   !...............................................................................................................................    
   
   CALL NWTC_Init(  )
   
   
   UnOut  = -1
   UnOut2 = -1
   CALL GetNewUnit( UnOut )
   CALL OpenFOutFile( UnOut, 'LAPACK_MUMPS_timings.txt' )

   CALL GetNewUnit(   UnOut2 )
   CALL OpenFOutFile( UnOut2, 'LAPACK_MUMPS_table.txt' )
   
   
   !...............................................................................................................................    
   ! Initialize MUMPS
   !...............................................................................................................................    
   CALL CPU_TIME(time_begin)
   
      CALL NWTC_SMUMPS_Init( mumps_par, ErrStat, ErrMsg )   
   
   CALL CPU_TIME ( time_end )
   
   WRITE(UnOut,'(A,ES15.7,A,ES15.7,A)') 'Initialization: MUMPS = ', time_end - time_begin, ' seconds' !; LAPACK = ', time_end2 - time_begin2, ' seconds.'
   
   
   UnIn = -1   
   do iTest = 27,19,-1

      CALL WrScr( 'Test'//TRIM(num2lstr(iTest)) )
      
   WRITE(UnOut,'(A)' )  '-----------------------------------------------------------------'
      
      
      !...............................................................................................................................    
      ! Read the input files for matrix A
      !...............................................................................................................................    
   
      if (iTest==26) then
         InputFileRoot = TRIM(InputFilePath)//'/10/Test21'
      elseif (iTest==27) then
         InputFileRoot = TRIM(InputFilePath)//'/05/Test21'
      else
         InputFileRoot = InputFilePath//'/Test'//TRIM(num2lstr(iTest))
      end if
      
                  
      !......................
      CALL GetNewUnit( UnIn )
      InputFileName = TRIM(InputFileRoot)//'.0.Jacobian.Sparse'
      CALL OpenFInpFile( UnIn, TRIM(InputFileName) )
      
      READ(UnIn, *) N, NZ
      read(UnIn,*) tmpLine
      
         CALL NWTC_SMUMPS_CreateMatrix(mumps_par,n, ErrStat, ErrMsg, nz)      
         mumps_par%NZ = nz
         
      DO I = 1, mumps_par%NZ
         READ(UnIn,*) mumps_par%IRN(I),mumps_par%JCN(I),mumps_par%A(I)
      END DO
        
      CLOSE(UnIn)
        
      !......................            
      InputFileName = TRIM(InputFileRoot)//'.0.Jacobian'
      CALL OpenFInpFile( UnIn, TRIM(InputFileName) )

         ALLOCATE( full_jac ( mumps_par%N,  mumps_par%N ) )
         ALLOCATE( ipiv (     mumps_par%N               ) )
         ALLOCATE( rhs  (     mumps_par%N               ) )
         
      read(UnIn,*) tmpLine
      DO I = 1, mumps_par%N
        read(UnIn,*) full_jac(i,:)
      END DO
      CLOSE(UnIn)
      
                        
      !......................
      InputFileName = TRIM(InputFileRoot)//'.0.RHS'
      CALL OpenFInpFile( UnIn, TRIM(InputFileName) )
      !......................      
      
      !CALL WrScr( '  finished reading matrices' )
      
      !...............................................................................................................................    
      ! Factor matrix A
      !...............................................................................................................................    
      
         ! MUMPS
         
   CALL CPU_TIME(time_begin)
      CALL NWTC_SMUMPS_AnalyzeAndFactor(mumps_par, ErrStat, ErrMsg)
   CALL CPU_TIME ( time_end )
         if (ErrStat /= ErrID_None) THEN
            DEALLOCATE( full_jac )
            DEALLOCATE( ipiv     )
            DEALLOCATE( rhs      )            
            CALL WrScr(TRIM(ErrMsg))
            CALL NWTC_SMUMPS_DestroyMatrix(mumps_par)           
            CYCLE
         end if         

         
   init_time =  time_end  - time_begin       
         
         ! LAPACK
               
   CALL CPU_TIME(time_begin2)
      CALL LAPACK_getrf( M=mumps_par%N, N=mumps_par%N, A=FULL_JAC, IPIV=ipiv, ErrStat=ErrStat, ErrMsg=ErrMsg )      
   CALL CPU_TIME ( time_end2 )
         if (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
         
         
   init_time2 =  time_end2  - time_begin2       
         
   
         ! Compare results
   WRITE(UnOut,'( 2(A,ES15.7) )') 'Test'//trim(num2lstr(iTest))//' factorization times: MUMPS = ', init_time, &
                                                                            ' seconds; LAPACK = ', init_time2
   
   
      !CALL WrScr( '  finished factoring matrices' )
   
      !...............................................................................................................................    
      ! solve Ax=b
      !...............................................................................................................................    
      READ(UnIn,*,iostat=ErrStat) mumps_par%RHS
      
      i_rhs = 0
      total_time  = 0.
      total_time2 = 0.
      
      do while (ErrStat == 0  )    
         rhs = mumps_par%RHS
         i_rhs = i_rhs + 1
         !call wrscr('RHS'//trim(num2lstr(i_rhs)) )
         
         ! MUMPS
         
   CALL CPU_TIME(time_begin)
         CALL NWTC_SMUMPS_Solve(mumps_par, ErrStat, ErrMsg)
   CALL CPU_TIME ( time_end )         
   total_time = total_time + time_end  - time_begin
         if (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
         ! LAPACK
               
   CALL CPU_TIME(time_begin2)
         CALL LAPACK_getrs( TRANS='N', N=mumps_par%N, A=FULL_JAC, IPIV=ipiv, B=rhs, ErrStat=ErrStat, ErrMsg=ErrMsg )         
   CALL CPU_TIME ( time_end2 )
   total_time2 = total_time2 + time_end2  - time_begin2
         if (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   
   
   WRITE(UnOut,'( 3(A,ES15.7) )') 'Test'//trim(num2lstr(iTest))//' solve times: MUMPS = ', time_end  - time_begin, &
                                                                    ' seconds; LAPACK = ', time_end2 - time_begin2, &
                                                         ' seconds. Norm of difference = ', TwoNorm( rhs - mumps_par%rhs )
                  
         
         READ(UnIn,*,iostat=ErrStat) mumps_par%RHS

      end do !while there is a line in the file
      
      
   WRITE(UnOut,'( 2(A,ES15.7),A,I5 )') 'Test'//trim(num2lstr(iTest))//' solve times: MUMPS = ', total_time, &
                                                                         ' seconds; LAPACK = ', total_time2, &
                                                                  ' seconds. Number of solves = ', i_rhs
      
   WRITE(UnOut2, '( 3(I15,1x), 8(ES15.7,1x), I15, 2(1x,ES15.7) )') &
                        iTest, mumps_par%n, mumps_par%nz, 100.0*mumps_par%nz/(mumps_par%n**2), &
                        0.0, 0.0, 0.0,                                                    & !time savings columns
                        init_time, init_time2,                                            & !init time 
                        total_time/i_rhs, total_time2/i_rhs,                              & !solve time (calculated)
                        i_rhs, total_time, total_time2
   
      
      !...............................................................................................................................    
      ! Deallocate arrays and close open files
      !...............................................................................................................................    
   
      CLOSE(UnIn)
      
      DEALLOCATE( full_jac )
      DEALLOCATE( ipiv     )
      DEALLOCATE( rhs      )
         
   CALL CPU_TIME(time_begin)
      CALL NWTC_SMUMPS_DestroyMatrix(mumps_par)           
   CALL CPU_TIME ( time_end )
   !WRITE (UnOut,*) 'MUMPS exit time: ', time_end - time_begin, ' seconds'
   WRITE(UnOut,* )  !blank line
         
   end do ! next test case
   
   CLOSE( UnOut  )
   CLOSE( UnOut2 )
   CALL NWTC_SMUMPS_End(mumps_par)           

   
END PROGRAM Compare_LAPACK_MUMPS

