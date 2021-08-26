program test_cusolver

   use MPI
   use cublas_v2
   use cusolverdn

   implicit none

   integer, parameter :: dp = selected_real_kind(15,300)

   character(10) :: arg1

   integer :: ierr
   integer :: nbasis
   integer :: nrand
   integer :: lwork
   integer, device :: ierr_d

   real(dp) :: t1
   real(dp) :: t2

   integer, allocatable :: seed(:)

   real(dp), allocatable :: mat(:,:)
   real(dp), allocatable :: tmp(:,:)

   real(dp), device, allocatable :: mat_d(:,:)
   real(dp), device, allocatable :: eval_d(:)
   real(dp), device, allocatable :: work_d(:)

   type(cusolverDnHandle) :: cusolver_h

   real(dp), parameter :: zero = 0.0_dp
   real(dp), parameter :: one = 1.0_dp

   ! Read command line arguments
   if(command_argument_count() == 1) then
      call get_command_argument(1,arg1)

      read(arg1,*) nbasis
      if(nbasis <= 0) then
         nbasis = 1000
      end if
   else
      write(6,"(2X,A)") "################################################"
      write(6,"(2X,A)") "##  Wrong number of command line arguments!!  ##"
      write(6,"(2X,A)") "##  Arg#1: Size of test matrix.               ##"
      write(6,"(2X,A)") "################################################"
      flush(6)
      stop
   end if

   ! Generate a random matrix
   call random_seed(size=nrand)

   allocate(seed(nrand))
   allocate(mat(nbasis,nbasis))
   allocate(tmp(nbasis,nbasis))
   allocate(mat_d(nbasis,nbasis))
   allocate(eval_d(nbasis))

   seed = 12345678

   call random_seed(put=seed)
   call random_number(tmp)

   ! Symmetrize test matrix
   mat = tmp+transpose(tmp)
   mat_d = mat

   deallocate(seed)
   deallocate(mat)
   deallocate(tmp)

   ierr = cusolverDnCreate(cusolver_h)
   ierr = cusolverDnDsyevd_bufferSize(cusolver_h,CUSOLVER_EIG_MODE_VECTOR,&
   & CUBLAS_FILL_MODE_UPPER,nbasis,mat_d,nbasis,eval_d,lwork)

   allocate(work_d(lwork))

   write(6,*)
   write(6,"(2X,A)") "Generated random matrix"
   write(6,"(2X,A,I10)") "| Size  :",nbasis
   write(6,*)
   flush(6)

   t1 = MPI_Wtime()

   ! Solve
   ierr = cusolverDnDsyevd(cusolver_h,CUSOLVER_EIG_MODE_VECTOR,&
   & CUBLAS_FILL_MODE_UPPER,nbasis,mat_d,nbasis,eval_d,work_d,lwork,ierr_d)

   t2 = MPI_Wtime()

   if(ierr /= 0) then
      write(6,"(2X,A)") "Error"
      flush(6)
      stop
   end if

   write(6,"(2X,A)") "cuSOLVER finished"
   write(6,"(2X,A,F10.3,A)") "| Time  :",t2-t1,"s"
   flush(6)

   deallocate(mat_d)
   deallocate(eval_d)
   deallocate(work_d)

   ierr = cusolverDnDestroy(cusolver_h)

end program
