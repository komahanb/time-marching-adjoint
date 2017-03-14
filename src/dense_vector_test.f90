#include "scalar.fpp"
!=====================================================================!
! Module to test the dense vector implementations
!=====================================================================!

module dense_vector_test

  use constants, only : WP, TINY
  use dense_vector_class
  
contains

  subroutine test_dense_vector(nvals)

    type(integer)      :: nvals
    type(dense_vector) :: A
    type(scalar)       :: val

    A = dense_vector(nvals)

    if (nvals .ne. A % get_size()) then
       print *, "Size mismatch"
    end if
    
    do i = 1, nvals

          ! Generate a random number
          call random_number(val)

          ! Insert the vector entries into the vector
          call A % add_entry(i, val)

          ! Get the vector entries from the vector
          if ( abs(val - A % get_entry(i)) .gt. TINY ) then
             print *, "DENSE VECTOR error", i, val, A % get_entry(i)
          end if

       end do

  end subroutine test_dense_vector

end module dense_vector_test
