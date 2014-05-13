module test
use utilities
implicit none
contains
subroutine test_utilities()
        integer, parameter :: n = 1000000
        integer, dimension(n) :: numbers
        integer, dimension(n/100) :: tofind, indices
        integer, dimension(4,n) :: numbers2
        integer, dimension(4,n/100)  :: tofind2
        integer :: i, j
        logical :: test
        do i = 1, n
                numbers(i) = irand()
        end do
        tofind(:) = numbers(1:n/100)
        call sort(numbers)
        test = .FALSE.
        do i = 2, n
        if ( numbers(i) < numbers(i-1) ) then
                test = .TRUE.
        end if
        end do
        if ( test ) then
                write(*,*) "Error sort: ", __FILE__, ": ", __LINE__
        end if
        call binsearch(numbers, tofind, indices)
        if ( any(numbers(indices) /= tofind ) ) then
                write(*,*) "Error binsearch: ", __FILE__, ": ", __LINE__
        end if
        do i = 1, n
        do j = 1, 4
                numbers2(j, i) = irand()
        end do
        end do
        tofind2(:,:) = numbers2(:,1:n/100)
        call sort(numbers2)
        test = .FALSE.
        do i = 2, n
        if ( arrayunder(numbers2(:,i), numbers2(:,i-1)) ) then
                test = .TRUE.
        end if
        end do
        if ( test ) then
                write(*,*) "Error sort: ", __FILE__, ": ", __LINE__
        end if
        call binsearch(numbers2, tofind2, indices)
        if ( any(numbers2(:,indices) /= tofind2 ) ) then
                write(*,*) "Error binsearch: ", __FILE__, ": ", __LINE__
        end if
                
end subroutine
subroutine test_concat()
        integer, dimension(:), allocatable :: a, b
        allocate(a(1000))
        a = 1
        call concat(a, b)

end subroutine        

subroutine test_sorts()
        integer, dimension(2,100) :: array, ori
        real*8, dimension(2,100) :: numbers
        integer, dimension(100) :: indices
        integer :: i, j, ind
        do i = 1, 100
        do j = 1, 2
        array(j,i) = irand()
        ori(j,i) = array(j,i)
        numbers(j,i) = real(array(j,i), kind=8)
        end do
        indices(i) = i 
        end do
        call sort(array, indices, numbers)
        do i = 1, 100
                call binsearch(array, ori(:,i), ind)
                if ( ind == 0) then
                        write(*,*) "error"
                end if
        end do
end subroutine

end module test
