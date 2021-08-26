program test_elpa

   use ELPA
   use MPI

   implicit none

   integer, parameter :: dp = selected_real_kind(15,300)

   character(10) :: arg1
   character(10) :: arg2
   character(10) :: arg3
   character(10) :: arg4

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
   integer :: method
   integer :: cpu
   integer :: nlrow
   integer :: nlcol
   integer :: ldm
   integer :: nrand

   real(dp) :: t1
   real(dp) :: t2

   integer, allocatable :: seed(:)

   real(dp), allocatable :: mat(:,:)
   real(dp), allocatable :: tmp(:,:)
   real(dp), allocatable :: evec(:,:)
   real(dp), allocatable :: eval(:)

   class(elpa_t), pointer :: eh

   integer, external :: numroc

   real(dp), parameter :: zero = 0.0_dp
   real(dp), parameter :: one = 1.0_dp

   ! Initialize MPI
   call MPI_Init(ierr)
   comm = MPI_COMM_WORLD
   call MPI_Comm_size(comm,nproc,ierr)
   call MPI_Comm_rank(comm,myid,ierr)

   ! Read command line arguments
   if(command_argument_count() == 4) then
      call get_command_argument(1,arg1)
      call get_command_argument(2,arg2)
      call get_command_argument(3,arg3)
      call get_command_argument(4,arg4)

      read(arg1,*) nbasis
      if(nbasis <= 0) then
         nbasis = 1000
      end if

      read(arg2,*) blk
      read(arg3,*) method
      read(arg4,*) cpu
   else
      if(myid == 0) then
         write(6,"(2X,A)") "################################################"
         write(6,"(2X,A)") "##  Wrong number of command line arguments!!  ##"
         write(6,"(2X,A)") "##  Arg#1: Size of test matrix.               ##"
         write(6,"(2X,A)") "##  Arg#2: Block size                         ##"
         write(6,"(2X,A)") "##  Arg#3: 1 = ELPA 1-stage                   ##"
         write(6,"(2X,A)") "##         2 = ELPA 2-stage                   ##"
         write(6,"(2X,A)") "##  Arg#4: 1 = CPU                            ##"
         write(6,"(2X,A)") "##         2 = GPU                            ##"
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

   deallocate(tmp)
   allocate(evec(nlrow,nlcol))
   allocate(eval(nbasis))

   if(myid == 0) then
      write(6,*)
      write(6,"(2X,A)") "Generated random matrix"
      write(6,"(2X,A,I10)") "| Size       :",nbasis
      write(6,"(2X,A,I10)") "| Block size :",blk
      write(6,"(2X,A,I10)") "| MPI tasks  :",nproc
      write(6,*)
      flush(6)
   end if

   ! Initialize ELPA
   ierr = elpa_init(20180525)
   eh => elpa_allocate(ierr)

   call eh%set("na",nbasis,ierr)
   call eh%set("nev",nbasis,ierr)
   call eh%set("nblk",blk,ierr)
   call eh%set("local_nrows",nlrow,ierr)
   call eh%set("local_ncols",nlcol,ierr)
   call eh%set("mpi_comm_parent",comm,ierr)
   call eh%set("process_row",myprow,ierr)
   call eh%set("process_col",mypcol,ierr)
   call eh%set("timings",1,ierr)

   ierr = eh%setup()

   if(ierr /= 0) then
      write(6,"(2X,A)") "Error: setup"
      flush(6)
      call MPI_Abort(comm,0,ierr)
      stop
   end if

   if(method == 1) then
      call eh%set("solver",ELPA_SOLVER_1STAGE,ierr)
   else
      call eh%set("solver",ELPA_SOLVER_2STAGE,ierr)
   end if

   if(cpu == 1) then
      call eh%set("gpu",0,ierr)
      call eh%set("real_kernel",ELPA_2STAGE_REAL_GENERIC,ierr)
   else
      call eh%set("gpu",1,ierr)
      call eh%set("real_kernel",ELPA_2STAGE_REAL_GPU,ierr)
   end if

   t1 = MPI_Wtime()

   ! Solve
   call eh%eigenvectors(mat,eval,evec,ierr)

   t2 = MPI_Wtime()

   if(ierr /= 0) then
      write(6,"(2X,A)") "Error: solve"
      flush(6)
      call MPI_Abort(comm,0,ierr)
      stop
   end if

   if(myid == 0) then
      write(6,"(2X,A)") "ELPA finished"
      write(6,"(2X,A,F10.3,A)") "| Time  :",t2-t1,"s"
      flush(6)

!      call eh%print_times()
   end if

   ! Finalize ELPA
   call elpa_deallocate(eh,ierr)

   nullify(eh)

   deallocate(mat)
   deallocate(eval)
   deallocate(evec)
   deallocate(seed)

   call BLACS_Gridexit(ctxt)
   call BLACS_Exit(1)
   call MPI_Finalize(ierr)

end program
