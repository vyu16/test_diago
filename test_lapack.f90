program test_lapack

   use MPI

   implicit none

   integer, parameter :: dp = selected_real_kind(15,300)

   character(10) :: arg1

   integer :: ierr
   integer :: nbasis
   integer :: nrand
   integer :: lwork

   real(dp) :: t1
   real(dp) :: t2

   integer, allocatable :: seed(:)

   real(dp), allocatable :: mat(:,:)
   real(dp), allocatable :: tmp(:,:)
   real(dp), allocatable :: work(:)
   real(dp), allocatable :: eval(:)

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

   seed(:) = 12345678

   call random_seed(put=seed)
   call random_number(tmp)

   ! Symmetrize test matrix
   mat(:,:) = tmp+transpose(tmp)

   deallocate(seed)
   deallocate(tmp)
   allocate(eval(nbasis))
   allocate(work(1))

   call dsyev("v","u",nbasis,mat,nbasis,eval,work,-1,ierr)

   lwork = ceiling(work(1))

   deallocate(work)
   allocate(work(lwork))

   write(6,*)
   write(6,"(2X,A)") "Generated random matrix"
   write(6,"(2X,A,I10)") "| Size  :",nbasis
   write(6,*)
   flush(6)

   t1 = MPI_Wtime()

   ! Solve
   call dsyev('v','u',nbasis,mat,nbasis,eval,work,lwork,ierr)

   t2 = MPI_Wtime()

   if(ierr /= 0) then
      write(6,"(2X,A)") "Error"
      flush(6)
      stop
   end if

   write(6,"(2X,A)") "LAPACK finished"
   write(6,"(2X,A,F10.3,A)") "| Time  :",t2-t1,"s"
   flush(6)

   deallocate(mat)
   deallocate(eval)
   deallocate(work)

end program
