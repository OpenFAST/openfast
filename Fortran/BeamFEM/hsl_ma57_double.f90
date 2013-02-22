! *******************************************************************
! COPYRIGHT (c) 2000 Council for the Central Laboratory
!               of the Research Councils
! All rights reserved.
!
! None of the comments in this Copyright notice between the lines
! of asterisks shall be removed or altered in any way.
!
! This Package is intended for compilation without modification,
! so most of the embedded comments have been removed.
!
! ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
! SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
!
! Please note that for an ACADEMIC Licence:
!
! 1. The Packages may only be used for academic research or teaching
!    purposes by the Licensee, and must not be copied by the Licensee for
!    use by any other persons. Use of the Packages in any commercial
!    application shall be subject to prior written agreement between
!    Hyprotech UK Limited and the Licensee on suitable terms and
!    conditions, which will include financial conditions.
! 2. All information on the Package is provided to the Licensee on the
!    understanding that the details thereof are confidential.
! 3. All publications issued by the Licensee that include results obtained
!    with the help of one or more of the Packages shall acknowledge the
!    use of the Packages. The Licensee will notify the Numerical Analysis
!    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
! 4. The Packages may be modified by or on behalf of the Licensee
!    for such use in research applications but at no time shall such
!    Packages or modifications thereof become the property of the
!    Licensee. The Licensee shall make available free of charge to the
!    copyright holder for any purpose all information relating to
!    any modification.
! 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
!    direct or consequential loss or damage whatsoever arising out of
!    the use of Packages by the Licensee.
! *******************************************************************
!
! Original date 13 September 1999
! 01/11/00  Entries in IW initialized to zero in MA57O/OD to avoid copy
!           of unassigned variables by MA57E/ED.
!           AINPUT and IINPUT reset in call to MA57E/ED.
!           sinfo%more initialized to 0.
!           MA27 changed to MA57 in allocation error return.
!           calls to MC41 changed to MC71 calls.
! 19/12/01  Optional parameter added to enquire routine to return perturbation
!           factors%pivoting added to hold pivoting option
!           la increased when Schnabel-Eskow option invoked
!           sinfo%more removed.
! 10/05/02  factors of data type ma57_factors made intent(in) in solves
! 24/5/04 Statment functions in MA57U/UD replaced by in-line code.

! 12th July 2004 Version 1.0.0. Version numbering added.
! 20/07/04  Several changes incorporated for HSL 2004 code.
!           Removed unused INT,ABS from MA57U/UD
!           INFO(32), INFO(33), and INFO(34) added
!           INFO(32): number of zeros in the triangle of the factors
!           INFO(33): number of zeros in the rectangle of the factors
!           INFO(34): number of zero columns in rectangle of the factors
!           Static pivoting available (controlled by CNTL(4))

! 31st July 2004 Version 2.0.0 established at HSL 2004 release.
! 1 Nov 2004 Version 2.0.1. Defaults changed, since they are taken from MA57.
! 7 March 2005 Version 3.0.0.  Mostly a large set of changes inherited from
!           the MA57 (F77) code.  See the header for MA57A/AD for details
!           of these.  Changes within the F90 code include: change of
!           default for CONTROL%ORDERING, introduction of new component
!           AINFO%ORDERING, correction of error in deallocate statement
!           in solve1 and solve2, and changes to diagnostic output format.
! 14 June 2005  Version 3.1.0. The space needed for ma57bd now recomputed in
!           the factorize entry so that the user can change the value of
!           control%pivoting and control%scaling between calls to analyse
!           and factorize.  cond2, berr, berr2, error added to sinfo to
!           return condition number for category 2 equations, backward
!           error, and forward error for category 2 equations as in
!           Arioli, Demmel, and Duff (1989).
!           control%convergence added so that cntl(3) can be set for
!           ma57dd call.
!           Small change to write statement for failure in finalize.
! 28 February 2006. Version 3.2.0. To avoid memory leaks, ma57_initialize and
!        ma57_finalize leave pointer component null instead having zero size.
! 4 December 2006. Version 4.0.0. Pointer arrays replaced by allocatables.
! 20 September 2007.  Version 4.1.0. Code tidied to interface with
!           Version 3.2.0 of MA57 .. an added control parameter
!           rank_deficient allows the option of dropping blocks of
!           small entries during the factorization.

module hsl_ma57_double
   use hsl_zd11_double
   implicit none
   integer, parameter, private :: wp = kind(0.0d0)
   type ma57_factors
     private
      integer, allocatable :: keep(:)
      integer, allocatable :: iw(:)
      real(wp), allocatable :: val(:)
      integer :: n
      integer :: nrltot
      integer :: nirtot
      integer :: nrlnec
      integer :: nirnec
      integer :: pivoting
      integer :: scaling
      integer :: static
   end type ma57_factors
   type ma57_control
      real(wp) :: multiplier
      real(wp) :: reduce
      real(wp) :: u
      real(wp) :: static_tolerance
      real(wp) :: static_level
      real(wp) :: tolerance
      real(wp) :: convergence
      integer :: lp
      integer :: wp
      integer :: mp
      integer :: sp
      integer :: ldiag
      integer :: nemin
      integer :: factorblocking
      integer :: solveblocking
      integer :: la
      integer :: liw
      integer :: maxla
      integer :: maxliw
      integer :: pivoting
      integer :: thresh
      integer :: ordering
      integer :: scaling
      integer :: rank_deficient
   end type ma57_control
   type ma57_ainfo
      real(wp) :: opsa
      real(wp) :: opse
      integer :: flag
      integer :: more
      integer :: nsteps
      integer :: nrltot
      integer :: nirtot
      integer :: nrlnec
      integer :: nirnec
      integer :: nrladu
      integer :: niradu
      integer :: ncmpa
      integer :: oor
      integer :: dup
      integer :: maxfrt
      integer :: stat
   end type ma57_ainfo
   type ma57_finfo
      real(wp) :: opsa
      real(wp) :: opse
      real(wp) :: opsb
      real(wp) :: maxchange
      real(wp) :: smin
      real(wp) :: smax
      integer :: flag
      integer :: more
      integer :: maxfrt
      integer :: nebdu
      integer :: nrlbdu
      integer :: nirbdu
      integer :: nrltot
      integer :: nirtot
      integer :: nrlnec
      integer :: nirnec
      integer :: ncmpbr
      integer :: ncmpbi
      integer :: ntwo
      integer :: neig
      integer :: delay
      integer :: signc
      integer :: static
      integer :: modstep
      integer :: rank
      integer :: stat
   end type ma57_finfo
   type ma57_sinfo
      real(wp) :: cond
      real(wp) :: cond2
      real(wp) :: berr
      real(wp) :: berr2
      real(wp) :: error
      integer :: flag
      integer :: stat
   end type ma57_sinfo
   interface ma57_solve
      module procedure ma57_solve1,ma57_solve2
   end interface
   interface ma57_part_solve
      module procedure ma57_part_solve1,ma57_part_solve2
   end interface
contains
   subroutine ma57_initialize(factors,control)
      type(ma57_factors), intent(out), optional :: factors
      type(ma57_control), intent(out), optional :: control
      integer icntl(20),stat
      real(wp) cntl(5)
      if (present(factors)) then
        factors%n = 0
        deallocate(factors%keep,factors%val,factors%iw,stat=stat)
      end if
      if (present(control)) then
          call ma57id(cntl,icntl)
          control%u = cntl(1)
          control%tolerance = cntl(2)
          control%convergence = cntl(3)
          control%static_tolerance = cntl(4)
          control%static_level = cntl(5)
          control%lp = icntl(1)
          control%wp = icntl(2)
          control%mp = icntl(3)
          control%sp = icntl(4)
          control%ldiag = icntl(5)
          control%pivoting = icntl(7)
          control%ordering = icntl(6)
          control%scaling = icntl(15)
          control%factorblocking = icntl(11)
          control%nemin = icntl(12)
          control%solveblocking = icntl(13)
          control%thresh = icntl(14)
          control%rank_deficient = icntl(16)
          control%la = 0
          control%liw = 0
          control%maxla = huge(0)
          control%maxliw = huge(0)
          control%multiplier = 2.0
          control%reduce     = 2.0
      end if
    end subroutine ma57_initialize
   subroutine ma57_analyse(matrix,factors,control,ainfo,perm)
      type(zd11_type), intent(in) :: matrix
      type(ma57_factors), intent(inout) :: factors
      type(ma57_control), intent(in) :: control
      type(ma57_ainfo), intent(out) :: ainfo
      integer, intent(in), optional :: perm(matrix%n)
      integer, allocatable :: iw1(:)
      integer :: lkeep,n,ne,stat,icntl(20),info(40),rspace
      real(wp) rinfo(20)
      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%mp
      icntl(4) = control%sp
      icntl(5) = control%ldiag
      icntl(6) = control%ordering
      icntl(12)  = control%nemin
      icntl(14) = control%thresh
      icntl(7)  = control%pivoting
      icntl(15) = control%scaling
      n = matrix%n
      ne = matrix%ne
      stat = 0
      lkeep = 5*n+ne+max(n,ne)+42
      if(allocated(factors%keep)) then
         if(size(factors%keep)/=lkeep) then
            deallocate(factors%keep,stat=stat)
            if (stat/=0) go to 100
            allocate(factors%keep(lkeep),stat=stat)
            if (stat/=0) go to 100
         end if
      else
         allocate(factors%keep(lkeep),stat=stat)
         if (stat/=0) go to 100
      end if
      if (present(perm)) then
         factors%keep(1:n) = perm(1:n)
         icntl(6)=1
      end if
      allocate (iw1(5*n),stat=stat)
      if (stat/=0) go to 100
      call ma57ad(n,ne,matrix%row,matrix%col, &
               lkeep,factors%keep,iw1,icntl,info,rinfo)
      rspace = 0
      if (control%pivoting == 4) rspace = rspace + n + 5
      if (control%scaling == 1)  rspace = rspace + n
      factors%n = n
      factors%nrltot = info(9)  - rspace
      factors%nirtot = info(10)
      factors%nrlnec = info(11) - rspace
      factors%nirnec = info(12)
      ainfo%opsa   = rinfo(1)
      ainfo%opse   = rinfo(2)
      ainfo%flag   = info(1)
      if (info(1) == -18) ainfo%flag   = -10
      ainfo%more   = info(2)
      ainfo%oor    = info(3)
      ainfo%dup    = info(4)
      ainfo%nrladu = info(5)
      ainfo%niradu = info(6)
      ainfo%maxfrt = info(7)
      ainfo%nsteps  = info(8)
      ainfo%nrltot = info(9)
      ainfo%nirtot = info(10)
      ainfo%nrlnec = info(11)
      ainfo%nirnec = info(12)
      ainfo%ncmpa  = info(13)
      deallocate (iw1, stat=stat)
      if (stat/=0) go to 100
      return
  100 if (control%ldiag>0 .and. control%lp>0 ) &
          write (control%lp,'(/a/a,i5)') &
         'Error return from MA57_ANALYSE: flag = -3', &
         'Allocate or deallocate failed with STAT=',stat
       ainfo%flag = -3
       ainfo%stat = stat
   end subroutine ma57_analyse
   subroutine ma57_factorize(matrix,factors,control,finfo)
      type(zd11_type), intent(in) :: matrix
      type(ma57_factors), intent(inout) :: factors
      type(ma57_control), intent(in) :: control
      type(ma57_finfo), intent(out) :: finfo
      integer :: la,liw,lkeep,oldla,oldliw
      integer stat
      integer icntl(20),info(40),n,exla,expne
      real(wp) cntl(5),rinfo(20)
      integer, allocatable :: iwork(:)
      real(wp), allocatable :: temp(:)
      integer, allocatable :: itemp(:)
      n = matrix%n
      lkeep = 5*n+matrix%ne+max(n,matrix%ne)+42
      allocate (iwork(n),stat=stat)
      if (stat/=0) go to 100
      if(factors%n/=matrix%n) then
         if (control%ldiag>0 .and. control%lp>0 ) &
         write (control%lp,'(/a/a,i12,a,i12)') &
         'Error return from MA57_FACTORIZE: flag = -1', &
         'MATRIX%N has the value', &
         matrix%n,' instead of',factors%n
       finfo%flag = -1
       finfo%more = factors%n
       return
      end if
      cntl(1) = control%u
      cntl(2) = control%tolerance
      cntl(4) = control%static_tolerance
      cntl(5) =  control%static_level
      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%mp
      icntl(4) = control%sp
      icntl(5) = control%ldiag
      icntl(7) = control%pivoting
      factors%pivoting = control%pivoting
      factors%scaling = control%scaling
      icntl(8) = 1
      icntl(11) = control%factorblocking
      icntl(12) = control%nemin
      icntl(15) = control%scaling
      icntl(16) =  control%rank_deficient
      stat = 0
      expne = factors%keep(matrix%n+2)
      la = control%la
      if(la<factors%nrlnec)then
         la = 0
         if(allocated(factors%val))la = size(factors%val)
         if(la>control%reduce*factors%nrltot) la = factors%nrltot
         if(la<factors%nrlnec) la = factors%nrltot
      end if
         exla = 0
         if(control%pivoting == 4) exla =  factors%n + 5
         if(control%scaling == 1)  exla =  exla + factors%n
         la = max(la+exla,exla+expne+2)
         if(control%scaling == 1) la = max(la,exla+3*expne+3*factors%n+1)
      if(allocated(factors%val))then
         if(la/=size(factors%val))then
            deallocate(factors%val,stat=stat)
            if (stat/=0) go to 100
            allocate(factors%val(la),stat=stat)
            if (stat/=0) go to 100
         end if
      else
         allocate(factors%val(la),stat=stat)
         if (stat/=0) go to 100
      end if
      liw = control%liw
      if(liw<factors%nirnec)then
         liw = 0
         if(allocated(factors%iw))liw = size(factors%iw)
         if(liw>control%reduce*factors%nirnec) liw = factors%nirtot
         if(liw<factors%nirnec) liw = factors%nirtot
      end if
      if(control%scaling == 1) liw = max(liw,3*expne+5*factors%n+1)
      if(allocated(factors%iw))then
         if(liw/=size(factors%iw))then
            deallocate(factors%iw,stat=stat)
            if (stat/=0) go to 100
            allocate(factors%iw(liw),stat=stat)
            if (stat/=0) go to 100
         end if
      else
         allocate(factors%iw(liw),stat=stat)
         if (stat/=0) go to 100
      end if
      do
         call ma57bd(matrix%n,matrix%ne,matrix%val,factors%val, &
                  la,factors%iw,liw,lkeep,factors%keep,iwork,icntl,cntl, &
                  info,rinfo)
         finfo%flag = info(1)
         if (info(1)==11) then
            oldliw = liw
            liw = control%multiplier*liw
            if (liw>control%maxliw) then
               if (control%ldiag>0 .and. control%lp>0 ) &
                 write (control%lp,'(/a/a,i10)') &
                 'Error return from MA57_FACTORIZE: iflag = -8 ', &
                 'Main integer array needs to be bigger than', control%maxliw
               finfo%flag = -8
               return
            end if
            allocate (itemp(oldliw),stat=stat)
            if (stat/=0) go to 100
            itemp(1:oldliw) = factors%iw(1:oldliw)
            deallocate (factors%iw,stat=stat)
            if (stat/=0) go to 100
            allocate (factors%iw(liw),stat=stat)
            if (stat/=0) go to 100
            call ma57ed &
                  (matrix%n,1,factors%keep,factors%val,la,factors%val,la, &
                   itemp,oldliw,factors%iw,liw,info)
            deallocate (itemp,stat=stat)
            if (stat/=0) go to 100
         else if (info(1)==10) then
            oldla = la
            la = control%multiplier*la
            if (la>control%maxla) then
               if (control%ldiag>0 .and. control%lp>0 ) &
                 write (control%lp,'(/a/a,i10)') &
                 'Error return from MA57_FACTORIZE: flag = -7 ', &
                 'Main real array needs to be bigger than', control%maxla
               finfo%flag = -7
               return
            end if
            allocate (temp(oldla),stat=stat)
            if (stat/=0) go to 100
            temp(1:oldla) = factors%val(1:oldla)
            deallocate (factors%val,stat=stat)
            if (stat/=0) go to 100
            allocate (factors%val(la),stat=stat)
            if (stat/=0) go to 100
            call ma57ed &
                 (matrix%n,0,factors%keep,temp,oldla,factors%val,la, &
                 factors%iw,liw,factors%iw,liw,info)
            deallocate (temp,stat=stat)
            if (stat/=0) go to 100
         else
            exit
         end if
      end do
      deallocate (iwork,stat=stat)
      if (stat/=0) go to 100
      finfo%more = info(2)
      if (info(1)>=0) then
        finfo%nebdu  = info(14)
        finfo%nrlbdu = info(15)
        finfo%nirbdu = info(16)
        finfo%nrlnec = info(17)
        finfo%nirnec = info(18)
        finfo%nrltot = info(19)
        finfo%nirtot = info(20)
        finfo%maxfrt = info(21)
        finfo%ntwo   = info(22)
        finfo%delay  = info(23)
        finfo%neig   = info(24)
        finfo%rank   = info(25)
        finfo%signc  = info(26)
        finfo%static = info(35)
        factors%static = 0
        if (finfo%static > 0) factors%static = 1
        finfo%modstep= info(27)
        finfo%ncmpbr = info(28)
        finfo%ncmpbi = info(29)
        finfo%opsa   = rinfo(3)
        finfo%opse   = rinfo(4)
        finfo%opsb   = rinfo(5)
        finfo%smin   = rinfo(16)
        finfo%smax   = rinfo(17)
        if (finfo%modstep > 0) finfo%maxchange = rinfo(14)
      end if
      return
  100 if (control%ldiag>0 .and. control%lp>0 ) &
         write (control%lp,'(/a/a,i5)') &
         'Error return from MA57_ANALYSE: flag = -3', &
         'Allocate or deallocate failed with STAT=',stat
       finfo%flag = -3
       finfo%stat = stat
   end subroutine ma57_factorize
   subroutine ma57_solve2(matrix,factors,x,control,sinfo,rhs,iter,cond)
      type(zd11_type), intent(in) :: matrix
      type(ma57_factors), intent(in) :: factors
      real(wp), intent(inout) :: x(:,:)
      type(ma57_control), intent(in) :: control
      type(ma57_sinfo), intent(out) :: sinfo
      real(wp), optional, intent(in) :: rhs(:,:)
      integer, optional, intent(in) :: iter
      integer, optional, intent(in) :: cond
      integer icntl(20),info(40),job,stat
      real(wp) cntl(5),rinfo(20),zero
      integer i,lw,n,nrhs
      integer, allocatable :: iwork(:)
      real(wp), allocatable :: work(:),resid(:),start(:,:)
      parameter (zero=0.0d0)
      n = matrix%n
      nrhs = size(x,2)
      stat = 0
      sinfo%flag = 0
      lw = n*nrhs
      if (present(rhs))  lw = n
      if (present(iter)) lw = 3*n
      if (factors%static == 1) lw = 3*n
      if (present(cond)) lw = 4*n
      allocate (iwork(n),work(lw),resid(n),stat=stat)
      if (stat/=0) go to 100
      if (factors%static == 1 .and. .not. present(rhs)) &
          allocate(start(n,nrhs),stat=stat)
      if (stat/=0) go to 100
      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%mp
      icntl(4) = control%sp
      icntl(5) = control%ldiag
      icntl(13) = control%solveblocking
      icntl(15) = control%scaling
      cntl(3) = control%convergence
      if (present(iter)) then
        icntl(9)=100
        icntl(10)=0
        if (present(cond)) icntl(10)=1
        job = 2
        do i = 1,nrhs
        call ma57dd(job,matrix%n,matrix%ne,matrix%val,matrix%row,matrix%col, &
          factors%val,size(factors%val),factors%iw, &
          size(factors%iw),rhs(:,i), &
          x(:,i),resid,work,iwork,icntl,cntl,info,rinfo)
        enddo
        if (present(cond)) then
          sinfo%cond  = rinfo(11)
          sinfo%cond2 = rinfo(12)
          sinfo%berr  = rinfo(6)
          sinfo%berr2 = rinfo(7)
          sinfo%error = rinfo(13)
        endif
      else
        if(present(rhs)) then
          icntl(9) = 1
          icntl(10) = 0
          job = 2
          do i = 1,nrhs
          call ma57dd(job,matrix%n,matrix%ne,matrix%val,matrix%row, &
            matrix%col,factors%val,size(factors%val),factors%iw, &
            size(factors%iw),rhs(:,i), &
            x(:,i),resid,work,iwork,icntl,cntl,info,rinfo)
          enddo
        else
          if (factors%static == 1) then
            icntl(9) = 1
            icntl(10) = 0
            job = 2
            start = zero
            do i = 1,nrhs
            call ma57dd(job,matrix%n,matrix%ne,matrix%val,matrix%row, &
              matrix%col,factors%val,size(factors%val),factors%iw, &
              size(factors%iw),x(:,i), &
              start(:,i),resid,work,iwork,icntl,cntl,info,rinfo)
            enddo
            x = start
          else
            job=1
            call ma57cd(job,factors%n,factors%val,size(factors%val), &
              factors%iw,size(factors%iw),nrhs,x,size(x,1),   &
              work,nrhs*matrix%n,iwork,icntl,info)
          end if
        end if
      endif
      if (info(1) == -8 .or. info(1) == -14) sinfo%flag = -11
      deallocate (iwork,work,resid,stat=stat)
      if (factors%static == 1 .and. .not. present(rhs)) &
          deallocate (start,stat=stat)
      if (stat==0) return
  100 if (control%ldiag>0 .and. control%lp>0 )  &
          write (control%lp,'(/a/a,i5)') &
         'Error return from MA57_ANALYSE: flag = -3', &
         'Allocate or deallocate failed with STAT=',stat
      sinfo%flag = -3
      sinfo%stat = stat
   end subroutine ma57_solve2
   subroutine ma57_solve1(matrix,factors,x,control,sinfo,rhs,iter,cond)
      type(zd11_type), intent(in) :: matrix
      type(ma57_factors), intent(in) :: factors
      real(wp), intent(inout) :: x(:)
      type(ma57_control), intent(in) :: control
      type(ma57_sinfo), intent(out) :: sinfo
      real(wp), optional, intent(in) :: rhs(:)
      integer, optional, intent(in) :: iter
      integer, optional, intent(in) :: cond
      integer icntl(20),info(40),job,stat
      real(wp) cntl(5),rinfo(20),zero
      integer n,nrhs,lw
      integer, allocatable :: iwork(:)
      real(wp), allocatable :: work(:),resid(:),start(:)
      parameter (zero=0.0d0)
      n = matrix%n
      nrhs = 1
      stat = 0
      sinfo%flag = 0
      lw = n
      if (present(rhs))  lw = n
      if (present(iter)) lw = 3*n
      if (factors%static == 1) lw = 3*n
      if (present(cond)) lw = 4*n
      allocate (iwork(n),work(lw),resid(n),stat=stat)
      if (stat/=0) go to 100
      if (factors%static == 1 .and. .not. present(rhs)) &
          allocate(start(n),stat=stat)
      if (stat/=0) go to 100
      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%mp
      icntl(4) = control%sp
      icntl(5) = control%ldiag
      icntl(13) = control%solveblocking
      icntl(15) = control%scaling
      cntl(3) = control%convergence
      stat = 0
      if (present(iter)) then
        icntl(9)=100
        icntl(10)=0
        if (present(cond)) icntl(10)=1
        job = 2
        call ma57dd(job,matrix%n,matrix%ne,matrix%val,matrix%row, &
          matrix%col,factors%val,size(factors%val),factors%iw, &
          size(factors%iw),rhs, &
          x,resid,work,iwork,icntl,cntl,info,rinfo)
        if (present(cond)) then
          sinfo%cond  = rinfo(11)
          sinfo%cond2 = rinfo(12)
          sinfo%berr  = rinfo(6)
          sinfo%berr2 = rinfo(7)
          sinfo%error = rinfo(13)
        endif
      else
        if(present(rhs)) then
          icntl(9) = 1
          icntl(10) = 0
          job = 2
          call ma57dd(job,matrix%n,matrix%ne,matrix%val,matrix%row, &
            matrix%col,factors%val,size(factors%val),factors%iw, &
            size(factors%iw),rhs, &
            x,resid,work,iwork,icntl,cntl,info,rinfo)
        else
          if (factors%static == 1) then
            icntl(9) = 1
            icntl(10) = 0
            job = 2
            start = zero
            call ma57dd(job,matrix%n,matrix%ne,matrix%val,matrix%row, &
              matrix%col,factors%val,size(factors%val),factors%iw, &
              size(factors%iw),x, &
              start,resid,work,iwork,icntl,cntl,info,rinfo)
            x = start
          else
            job=1
            call ma57cd(job,factors%n,factors%val,size(factors%val), &
              factors%iw,size(factors%iw),nrhs,x,size(x,1),   &
              work,nrhs*matrix%n,iwork,icntl,info)
          end if
        end if
      endif
      deallocate (iwork,work,resid,stat=stat)
      if (factors%static == 1 .and. .not. present(rhs))  &
          deallocate (start,stat=stat)
      if (stat==0) return
  100 if (control%ldiag>0 .and. control%lp>0 )  &
         write (control%lp,'(/a/a,i5)')  &
         'Error return from MA57_ANALYSE: flag = -3', &
         'Allocate or deallocate failed with STAT=',stat
      sinfo%flag = -3
      sinfo%stat = stat
   end subroutine ma57_solve1
   subroutine ma57_finalize(factors,control,info)
      type(ma57_factors), intent(inout) :: factors
      type(ma57_control), intent(in) :: control
      integer, intent(out) :: info
      integer :: inf
      info = 0
      inf = 0
      if (allocated(factors%keep)) deallocate(factors%keep,stat=inf)
      if (inf/=0) info = inf
      if (allocated(factors%iw)) deallocate(factors%iw,stat=inf)
      if (inf/=0) info = inf
      if (allocated(factors%val)) deallocate(factors%val,stat=inf)
      if (inf/=0) info = inf
      if (info==0) return
      if (control%ldiag>0 .and. control%lp>0 ) &
         write (control%lp,'(/a/a,i5)') &
         'Error return from MA57_FINALIZE:', &
         'Deallocate failed with STAT=',info
    end subroutine ma57_finalize
   subroutine ma57_enquire(factors,perm,pivots,d,perturbation,scaling)
      type(ma57_factors), intent(in) :: factors
      integer, intent(out), optional :: perm(factors%n),pivots(factors%n)
      real(wp), intent(out), optional :: d(2,factors%n)
      real(wp), intent(out), optional :: perturbation(factors%n)
      real(wp), intent(out), optional :: scaling(factors%n)
      real(wp) one,zero
      parameter (one=1.0d0,zero=0.0d0)
      integer block
      integer i
      integer ka
      integer k2
      integer kd
      integer kp
      integer kw
      integer ncols
      integer nrows
      logical two
      if(present(perm)) then
        perm = factors%keep(1:factors%n)
      endif
      if (present(perturbation)) then
        if (factors%pivoting == 4) then
          if (factors%scaling == 1) then
            perturbation = factors%val(size(factors%val) -  2*factors%n : &
                                       size(factors%val) -  factors%n -1)
          else
            perturbation = factors%val(size(factors%val) -  factors%n : &
                                       size(factors%val) - 1)
          endif
        else
          perturbation = zero
        end if
      endif
      if (present(scaling)) then
        if (factors%scaling == 1) then
          scaling = factors%val(size(factors%val) - factors%n : &
                                size(factors%val) - 1)
        else
          scaling = one
        end if
      endif
      if(present(pivots).or.present(d)) then
        ka = 1
        k2 = factors%iw(1)
        kd = 0
        kp = 0
        kw = 4
        if(present(d)) d = 0
          do block = 1, abs(factors%iw(3))
            ncols = factors%iw(kw)
            nrows = factors%iw(kw+1)
            if(present(pivots)) then
              pivots(kp+1:kp+nrows) = factors%iw(kw+2:kw+nrows+1)
              kp = kp + nrows
            end if
            if(present(d)) then
            two = .false.
            do i = 1,nrows
              kd = kd + 1
              d(1,kd) = factors%val(ka)
              if(factors%iw(kw+1+i)<0) two = .not.two
              if (two) then
                d(2,kd) = factors%val(k2)
                k2 = k2 + 1
              endif
              ka = ka + nrows + 1 - i
            end do
            ka = ka + nrows*(ncols-nrows)
          end if
          kw = kw + ncols + 2
        end do
      endif
   end subroutine ma57_enquire
   subroutine ma57_alter_d(factors,d,info)
      type(ma57_factors), intent(inout) :: factors
      real(wp), intent(in) :: d(2,factors%n)
      integer, intent(out) :: info
      integer block
      integer i
      integer ka
      integer k2
      integer kd
      integer kw
      integer ncols
      integer nrows
      logical two
      info = 0
      ka = 1
      k2 = factors%iw(1)
      kd = 0
      kw = 4
      do block = 1, abs(factors%iw(3))
        ncols = factors%iw(kw)
        nrows = factors%iw(kw+1)
        two = .false.
        do i = 1,nrows
          kd = kd + 1
          factors%val(ka) = d(1,kd)
          if(factors%iw(kw+1+i)<0) two = .not.two
          if (two) then
            factors%val(k2) = d(2,kd)
            k2 = k2 + 1
          else
            if (d(2,kd) /= 0) info = kd
          end if
          ka = ka + nrows + 1 - i
        end do
        ka = ka + nrows*(ncols-nrows)
        kw = kw + ncols + 2
      end do
     end subroutine ma57_alter_d
   subroutine ma57_part_solve2(factors,control,part,x,info)
      type(ma57_factors), intent(in) :: factors
      type(ma57_control), intent(in) :: control
      character, intent(in) :: part
      real(wp), intent(inout) :: x(:,:)
      integer, intent(out) :: info
      integer icntl(20),n,nrhs,inf(40)
      integer, allocatable :: iwork(:)
      real(wp), allocatable :: work(:)
      n = factors%n
      nrhs = size(x,2)
      allocate (iwork(n),work(nrhs*n),stat=info)
      if (info/=0) go to 100
      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%mp
      icntl(4) = control%sp
      icntl(5) = control%ldiag
      icntl(13) = control%solveblocking
      icntl(15) = control%scaling
      info = 0
      if(part=='L')then
         call ma57cd(2,factors%n,factors%val,size(factors%val),factors%iw, &
           size(factors%iw),nrhs,x,size(x,1),   &
           work,nrhs*factors%n,iwork,icntl,inf)
      else if(part=='D') then
         call ma57cd(3,factors%n,factors%val,size(factors%val),factors%iw, &
           size(factors%iw),nrhs,x,size(x,1),   &
           work,nrhs*factors%n,iwork,icntl,inf)
      else if(part=='U') then
         call ma57cd(4,factors%n,factors%val,size(factors%val),factors%iw, &
           size(factors%iw),nrhs,x,size(x,1),   &
           work,nrhs*factors%n,iwork,icntl,inf)
      end if
      deallocate (iwork,work,stat=info)
      if (info==0) return
  100 if (control%ldiag>0 .and. control%lp>0 )  &
          write (control%lp,'(/a/a,i5)') &
         'Error return from MA57_ANALYSE: flag = -3', &
         'Allocate or deallocate failed with STAT=',info
end subroutine ma57_part_solve2
   subroutine ma57_part_solve1(factors,control,part,x,info)
      type(ma57_factors), intent(in) :: factors
      type(ma57_control), intent(in) :: control
      character, intent(in) :: part
      real(wp), intent(inout) :: x(:)
      integer, intent(out) :: info
      integer inf(40),icntl(20),n,nrhs
      integer, allocatable :: iwork(:)
      real(wp), allocatable :: work(:)
      n = factors%n
      nrhs = 1
      allocate (iwork(n),work(n),stat=info)
      if (info/=0) go to 100
      icntl(1) = control%lp
      icntl(2) = control%wp
      icntl(3) = control%mp
      icntl(4) = control%sp
      icntl(5) = control%ldiag
      icntl(13) = control%solveblocking
      icntl(15) = control%scaling
      info = 0
      if(part=='L')then
         call ma57cd(2,factors%n,factors%val,size(factors%val),factors%iw, &
           size(factors%iw),nrhs,x,size(x,1),   &
           work,nrhs*factors%n,iwork,icntl,inf)
      else if(part=='D') then
         call ma57cd(3,factors%n,factors%val,size(factors%val),factors%iw, &
           size(factors%iw),nrhs,x,size(x,1),   &
           work,nrhs*factors%n,iwork,icntl,inf)
      else if(part=='U') then
         call ma57cd(4,factors%n,factors%val,size(factors%val),factors%iw, &
           size(factors%iw),nrhs,x,size(x,1),   &
           work,nrhs*factors%n,iwork,icntl,inf)
      end if
      deallocate (iwork,work,stat=info)
      if (info==0) return
  100 if (control%ldiag>0 .and. control%lp>0 ) &
          write (control%lp,'(/a/a,i5)') &
         'Error return from MA57_ANALYSE: flag = -3', &
         'Allocate or deallocate failed with STAT=',info
end subroutine ma57_part_solve1
end module hsl_ma57_double

