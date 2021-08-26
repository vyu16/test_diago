program test_scalapack

   use MPI

   implicit none

   integer, parameter :: dp = selected_real_kind(15,300)

   character(10) :: arg1
   character(10) :: arg2

   integer :: nproc
   integer :: nprow
   integer :: npcol
   integer :: myid
   integer :: myprow
   integer :: mypcol
   integer :: comm
   integer :: ierr
   integer :: blk
   integer :: ctxt
   integer :: desc(9)
   integer :: nbasis
   integer :: nlrow
   integer :: nlcol
   integer :: ldm
   integer :: nrand
   integer :: lwork

   real(dp) :: t1
   real(dp) :: t2

   integer, allocatable :: seed(:)

   real(dp), allocatable :: mat(:,:)
   real(dp), allocatable :: tmp(:,:)
   real(dp), allocatable :: work(:)
   real(dp), allocatable :: evec(:,:)
   real(dp), allocatable :: eval(:)

   integer, external :: numroc

   real(dp), parameter :: zero = 0.0_dp
   real(dp), parameter :: one = 1.0_dp

   ! Initialize MPI
   call MPI_Init(ierr)
   comm = MPI_COMM_WORLD
   call MPI_Comm_size(comm,nproc,ierr)
   call MPI_Comm_rank(comm,myid,ierr)

   ! Read command line arguments
   if(command_argument_count() == 2) then
      call get_command_argument(1,arg1)
      call get_command_argument(2,arg2)

      read(arg1,*) nbasis
      if(nbasis <= 0) then
         nbasis = 1000
      end if

      read(arg2,*) blk
   else
      if(myid == 0) then
         write(6,"(2X,A)") "################################################"
         write(6,"(2X,A)") "##  Wrong number of command line arguments!!  ##"
         write(6,"(2X,A)") "##  Arg#1: Size of test matrix.               ##"
         write(6,"(2X,A)") "##  Arg#2: Block size                         ##"
         write(6,"(2X,A)") "################################################"
         flush(6)
         call MPI_Abort(comm,0,ierr)
         stop
      end if
   end if

   ! Set up square-like processor grid
   do npcol = nint(sqrt(real(nproc))),2,-1
      if(mod(nproc,npcol) == 0) exit
   end do
   nprow = nproc/npcol

   ! Set up BLACS
   ctxt = comm

   call BLACS_Gridinit(ctxt,"r",nprow,npcol)
   call BLACS_Gridinfo(ctxt,nprow,npcol,myprow,mypcol)

   nlrow = numroc(nbasis,blk,myprow,0,nprow)
   nlcol = numroc(nbasis,blk,mypcol,0,npcol)

   ldm = max(nlrow,1)

   call descinit(desc,nbasis,nbasis,blk,blk,0,0,ctxt,ldm,ierr)

   ! Generate a random matrix
   call random_seed(size=nrand)

   allocate(seed(nrand))
   allocate(mat(nlrow,nlcol))
   allocate(tmp(nlrow,nlcol))

   seed(:) = myid+1

   call random_seed(put=seed)
   call random_number(mat)

   ! Symmetrize test matrix
   tmp(:,:) = mat

   call pdtran(nbasis,nbasis,one,tmp,1,1,desc,one,mat,1,1,desc)

   deallocate(seed)
   deallocate(tmp)
   allocate(evec(nlrow,nlcol))
   allocate(eval(nbasis))
   allocate(work(1))

   call pdsyev("v","u",nbasis,mat,1,1,desc,eval,evec,1,1,desc,work,-1,ierr)

   lwork = ceiling(work(1))

   deallocate(work)
   allocate(work(lwork))

   if(myid == 0) then
      write(6,*)
      write(6,"(2X,A)") "Generated random matrix"
      write(6,"(2X,A,I10)") "| Size       :",nbasis
      write(6,"(2X,A,I10)") "| Block size :",blk
      write(6,"(2X,A,I10)") "| MPI tasks  :",nproc
      write(6,*)
      flush(6)
   end if

   t1 = MPI_Wtime()

   ! Solve
   call pdsyev("v","u",nbasis,mat,1,1,desc,eval,evec,1,1,desc,work,lwork,ierr)

   t2 = MPI_Wtime()

   if(ierr /= 0) then
      write(6,"(2X,A)") "Error"
      flush(6)
      call MPI_Abort(comm,0,ierr)
      stop
   end if

   if(myid == 0) then
      write(6,"(2X,A)") "ScaLAPACK finished"
      write(6,"(2X,A,F10.3,A)") "| Time  :",t2-t1,"s"
      flush(6)
   end if

   deallocate(mat)
   deallocate(eval)
   deallocate(evec)
   deallocate(work)

   call BLACS_Gridexit(ctxt)
   call BLACS_Exit(1)
   call MPI_Finalize(ierr)

end program
