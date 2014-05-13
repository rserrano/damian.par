module utilities

use iso_c_binding, only : c_ptr, c_f_pointer
implicit none
include 'mpif.h' 
! Numeric constants.
real*8, parameter :: u_pi   = 3.14159265358979323846D0
real*8, parameter :: u_tol  = 1.0D-9
real*8, parameter :: u_inf  = 1.0D0/0.0D0
real*8, parameter :: u_ninf = -1.0D0/0.0D0
real*8, parameter :: u_nan  = 0.0D0/0.0D0
!Declaration of list types.

! linked_list_1d: Every node has a single value.

type linked_list_1d
        ! Single value.
        integer :: val
        ! Pointer to next node.
        type(linked_list_1d), pointer :: next => null(), last => null()
end type linked_list_1d

! linked_list_r1d: real*8 Every node has a single value.

type linked_list_r1d
        ! Single value.
        real*8 :: val
        ! Pointer to next node.
        type(linked_list_r1d), pointer :: next => null(), last => null()
end type linked_list_r1d

type linked_list_r2d
        real*8, dimension(:), allocatable :: val
        type(linked_list_r2d), pointer :: next => null(), last => null()
end type linked_list_r2d

! linked_list_2d: Every node has many values. The values can have different
! sizes.

type linked_list_2d
        ! Allocated array of values.
        integer, allocatable, dimension(:) :: val
        ! Pointer to next node.
        type(linked_list_2d), pointer :: next => null(), last => null()
end type linked_list_2d

! linked_list_2d: Every node has many values. The values can have different
! sizes.

type linked_list_il2d
        ! Allocated array of values.
        integer, allocatable, dimension(:) :: val
        logical*1, allocatable, dimension(:) :: logi
        ! Pointer to next node.
        type(linked_list_il2d), pointer :: next => null(), last => null()
end type linked_list_il2d

type ll2dp
        type(linked_list_2d), pointer :: p => null()
end type

type ll1dp
        type(linked_list_1d), pointer :: p => null()
end type

! Elements to redistribute for each process.
type csr
        integer, dimension(:), allocatable :: idx
        integer, dimension(:), allocatable :: vals
end type

! INTERFACE
! reduce_csr: reduces a CSR format array.

interface reduce_csr
        procedure reduce_csrs, reduce_csrc
end interface reduce_csr

! INTERFACE
! concat_csr: concatenates two CSR format arrays.

interface concat_csr
        procedure concat_csrs, concat_csrc
end interface concat_csr

! INTERFACE
! toarray_csr: converts a linked list to a CSR array.

interface toarray_csr
        procedure toarray_csrs, toarray_csrc
end interface toarray_csr

! INTERFACE
! insertll: inserts an element to a linked list.

interface insertll
        procedure insertll_1d, insertll_2d, insertll_r1d, insertll_r2d, insertll_il2d
end interface insertll

! INTERFACE
! toarray: converts a linked list to an array
!          deleting the list in the process.

interface toarray
        procedure toarray_1d, toarray_2d, toarray_r1d, toarray_r2d
end interface toarray

! INTERFACE
! concat: concatenates two arrays in the first one,
!         deleting the second

interface concat
        procedure concati, concati2, concatr, concatr2
end interface concat

! INTERFACE
! reduce: reduces an array according to the mask.

interface reduce 
        procedure reducei, reducei2, reducer, reducer2, reducel
end interface reduce

! INTERFACE:
! sort: sorts by an index array.

interface sort
        procedure sorts, sorti, sortr, sorts2, sorti2, sorti2ir, sortr2
end interface sort

interface binsearch
        procedure binsearch1, binsearch2, binsearch1s, binsearch2s
end interface binsearch

interface get_unique
        procedure get_uniques, get_uniquec
end interface 

interface ctofort
        procedure ctofortr1, ctofortr2, ctofortc1, ctofortc2, ctoforti1,&
                  ctoforti2, ctofortl1, ctofortl2
end interface

interface
subroutine c_free(ptr) bind(C, name="free")
        use iso_c_binding, only : c_ptr
        type(c_ptr), value :: ptr
        end subroutine
end interface

contains
! Include templates for utilities.
#include "util_temps.f90"

! FUNCTION
! real[comp]: Compares two real numbers. (equal, under or over)
! parameters:
! a, b:     Real numbers to compare

function realequal(a, b)
        logical*1 :: realequal
        real*8, intent(in) :: a, b
        realequal = (abs(a-b) < u_tol)
end function

function realunder(a, b)
        logical*1 :: realunder
        real*8, intent(in) :: a, b
        realunder = ((a-b) < -u_tol)
end function

function realover(a, b)
        logical*1 :: realover
        real*8, intent(in) :: a, b
        realover = ((a-b) > u_tol)
end function

function realfloor(a)
        real*8 :: realfloor
        real*8, intent(in) :: a
        realfloor = floor(a+u_tol)
end function

function realceiling(a)
        real*8 :: realceiling
        real*8, intent(in) :: a
        realceiling = ceiling(a-u_tol)
end function

function arrayunder(a, b)
        logical*1 :: arrayunder
        integer, dimension(:), intent(in) :: a, b
        integer :: i
        do i = 1, min(size(a,1), size(b,1))
                if ( a(i) < b(i) ) then
                        arrayunder = .TRUE.
                        return
                else if ( b(i) < a(i) ) then
                        arrayunder = .FALSE.
                        return
                end if
        end do
        arrayunder = .FALSE.
end function

! SURBOUTINE:
! distance:     Distance to point.
! p1:           point 1.
! p2:           point 2.

function distance(p1, p2)
        real*8 :: distance
        real*8, intent(in), dimension(2) :: p1, p2
        distance = sqrt((p1(1)-p2(1))**2 + (p1(2)-p2(2))**2)
end function

! SUBROUTINE
! insertll_1d: Insert an integer to a linked_list_2d.
! PARAMETERS
! ll:          the linked_list_2d object.
! val:         the value to insert.

subroutine insertll_1d(ll, val)
        type(linked_list_1d), pointer, intent(inout) :: ll
        integer, intent(in) :: val
        type(linked_list_1d), pointer :: node
        if ( .not. associated(ll) ) then
                allocate(ll)
                ll%val = val
                ll%last => ll
                node => ll%last
                nullify(node%next)
        else
                node => ll%last
                allocate(node%next)
                node => node%next
                node%val = val
                nullify(node%next)
                ll%last => node
        end if
end subroutine

! SUBROUTINE
! insertll_r1d: Insert a real*8 to a linked_list_r1d.
! PARAMETERS
! ll:          the linked_list_r1d object.
! val:         the real value to insert.

subroutine insertll_r1d(ll, val)
        type(linked_list_r1d), pointer, intent(inout) :: ll
        real*8, intent(in) :: val
        type(linked_list_r1d), pointer :: node
        if ( .not. associated(ll) ) then
                allocate(ll)
                ll%val = val
                ll%last => ll
                node => ll%last
                nullify(node%next)
        else
                node => ll%last
                allocate(node%next)
                node => node%next
                node%val = val
                nullify(node%next)
                ll%last => node
        end if
end subroutine

! SUBROUTINE
! toarray_1d: converts a linked_list_1d into a 1d array
!             and deletes all the information in the original list.
! PARAMETERS
! ll:         the linked list.
! array:      a 1d fortran array of the values.

subroutine toarray_1d(ll, array)
        type(linked_list_1d), pointer, intent(inout) :: ll
        integer, dimension(:), allocatable, intent(out) :: array
        type(linked_list_1d), pointer :: node, prev
        integer :: numcols, i
        numcols = 0
        if ( .not. associated(ll) ) then
                return
        end if
        node => ll
        do while( associated(node) )
                numcols = numcols + 1
                node => node%next
        end do
        allocate( array(numcols) )
        node => ll
        i = 0
        do while( associated(node) )
                i = i + 1
                array(i) = node%val
                prev => node
                node => node%next
                deallocate(prev)
        end do
        nullify(ll)
end subroutine


! SUBROUTINE
! insertll_2d: Insert an array to a linked_list_2d.
! PARAMETERS
! ll:          the linked_list_2d object.
! val:         the value to insert.

subroutine insertll_2d(ll, val)
        type(linked_list_2d), pointer, intent(inout) :: ll
        integer, dimension(:), intent(in) :: val
        type(linked_list_2d), pointer :: node
        if ( size(val,1) == 0 ) then
                return
        end if
        if ( .not. associated(ll) ) then
                allocate(ll)
                allocate(ll%val(size(val,1)))
                ll%val = val
                ll%last => ll
                node => ll%last
                nullify(node%next)
        else
                node => ll%last
                allocate(node%next)
                node => node%next
                allocate(node%val(size(val,1)))
                node%val = val
                nullify(node%next)
                ll%last => node
        end if
end subroutine

subroutine insertll_r2d(ll, val)
        type(linked_list_r2d), pointer, intent(inout) :: ll
        real*8, dimension(:), intent(in) :: val
        type(linked_list_r2d), pointer :: node 
        if ( size(val,1) == 0 ) then
                return
        end if
        if ( .not. associated(ll) ) then
                allocate(ll)
                allocate(ll%val(size(val,1)))
                ll%val = val
                ll%last => ll
                node => ll%last
                nullify(node%next)
        else
                node => ll%last
                allocate(node%next)
                node => node%next
                allocate(node%val(size(val,1)))
                node%val = val
                nullify(node%next)
                ll%last => node
        end if
end subroutine



! SUBROUTINE
! insertll_il2d: Insert an integer and a logical to array
!                to a linked_list_il2d.
! PARAMETERS
! ll:            the linked_list_il2d object.
! val:           the value to insert.

subroutine insertll_il2d(ll, val, logi)
        type(linked_list_il2d), pointer, intent(inout) :: ll
        integer, dimension(:), intent(in) :: val
        logical*1, dimension(:), intent(in) :: logi
        type(linked_list_il2d), pointer :: node 
        if ( size(val,1) == 0 ) then
                write(*,*) "Error: ", __FILE__, ":", __LINE__
                return
        end if
        if ( size(logi,1) /= size(val,1) ) then
                write(*,*) "Error: ", __FILE__, ":", __LINE__
        end if
        if ( .not. associated(ll) ) then
                allocate(ll)
                allocate(ll%val(size(val,1)))
                allocate(ll%logi(size(logi,1)))
                ll%val  = val
                ll%logi = logi
                ll%last => ll
                node => ll%last
                nullify(node%next)
        else
                node => ll%last
                allocate(node%next)
                node => node%next
                allocate(node%val(size(val,1)))
                allocate(node%logi(size(logi,1)))
                node%val = val
                node%logi = logi
                nullify(node%next)
                ll%last => node
        end if
end subroutine


! SUBROUTINE
! toarray_2d: converts a linked_list_2d into a 2d array
!             and deletes all the information in the original list.
! PARAMETERS
! ll:         the linked list.
! array:      a 2d fortran array with every value in a column.

subroutine toarray_2d(ll, array)
        type(linked_list_2d), pointer, intent(inout) :: ll
        integer, dimension(:,:), allocatable, intent(out) :: array
        type(linked_list_2d), pointer :: node, prev
        integer :: numcols, collen, i
        numcols = 0
        collen = 0
        if ( .not. associated(ll) ) then
                return
        end if
        node => ll
        do while(associated(node))
                numcols = numcols + 1
                collen = max(size(node%val(:),1),collen)
                node => node%next
        end do
        allocate(array(collen,numcols))
        node => ll
        i = 0
        do while(associated(node))
                i = i + 1
                array(:,i) = node%val(:)
                deallocate(node%val)
                prev => node
                node => node%next
                deallocate(prev)
        end do
        nullify(ll)
end subroutine

subroutine toarray_r2d(ll, array)
        type(linked_list_r2d), pointer, intent(inout) :: ll
        real*8, dimension(:,:), allocatable, intent(out) :: array
        type(linked_list_r2d), pointer :: node, prev
        integer :: numcols, collen, i
        numcols = 0
        collen = 0
        if ( .not. associated(ll) ) then
                return
        end if
        node => ll
        do while(associated(node))
                numcols = numcols + 1
                collen = max(size(node%val(:),1),collen)
                node => node%next
        end do
        allocate(array(collen,numcols))
        node => ll
        i = 0
        do while(associated(node))
                i = i + 1
                array(:,i) = node%val(:)
                deallocate(node%val)
                prev => node
                node => node%next
                deallocate(prev)
        end do
        nullify(ll)
end subroutine


! SUBROUTINE
! toarray_r1d: converts a linked_list_r1d into a 1d real*8 array
!             and deletes all the information in the original list.
! PARAMETERS
! ll:         the linked list.
! array:      a real 1d fortran array of the values.

subroutine toarray_r1d(ll, array)
        type(linked_list_r1d), pointer, intent(inout) :: ll
        real*8, dimension(:), allocatable, intent(out) :: array
        type(linked_list_r1d), pointer :: node, prev
        integer :: numcols, i
        numcols = 0
        if ( .not. associated(ll) ) then
                return
        end if
        node => ll
        do while( associated(node) )
                numcols = numcols + 1
                node => node%next
        end do
        allocate( array(numcols) )
        node => ll
        i = 0
        do while( associated(node) )
                i = i + 1
                array(i) = node%val
                prev => node
                node => node%next
                deallocate(prev)
        end do
        nullify(ll)
end subroutine

! SUBROUTINE
! toarray_csr: converts a linked_list_2d into csr format 
!             and deletes all the information in the original list.
! PARAMETERS
! ll:         the linked list that has its information deleted.
! idx:        the indices of the csr.
! vals:       the values of the csr.

subroutine toarray_csrs(ll, idx, vals)
        type(linked_list_2d), pointer, intent(inout) :: ll
        integer, dimension(:), allocatable, intent(out) :: idx, vals
        type(linked_list_2d), pointer :: node, prev
        integer :: numnodes, numvals, i
        numnodes = 0
        numvals = 0
        if (.not. associated(ll)) then
                return
        end if
        node => ll
        do while(associated(node))
                numnodes = numnodes + 1
                ! The next are errors. They should never happen
                ! if the list was created with insert_ll_2d.
                if (.not. allocated(node%val)) then
                        return
                end if
                if (size(node%val,1) == 0) then
                        return
                end if
                numvals = numvals + size(node%val,1)
                node => node%next
        end do
        allocate(idx(numnodes+1), vals(numvals))
        node => ll
        i = 0
        idx(1) = 1
        do while(associated(node))
                i = i + 1
                idx(i+1) = idx(i) + size(node%val,1)
                vals(idx(i):idx(i+1)-1) = node%val(:)
                deallocate(node%val)
                prev => node
                node => node%next
                deallocate(prev)
        end do
        nullify(ll)
end subroutine

! SUBROUTINE
! toarray_csr: converts a linked_list_2d into csr format 
!             and deletes all the information in the original list.
! PARAMETERS
! ll:         the linked list that has its information deleted.
! idx:        the indices of the csr.
! vals:       the values of the csr.

subroutine toarray_csrc(ll, idx, vals, logi)
        type(linked_list_il2d), pointer, intent(inout) :: ll
        integer, dimension(:), allocatable, intent(out) :: idx, vals
        logical*1, dimension(:), allocatable, intent(out) :: logi
        type(linked_list_il2d), pointer :: node, prev
        integer :: numnodes, numvals, i
        numnodes = 0
        numvals = 0
        if (.not. associated(ll)) then
                return
        end if
        node => ll
        do while(associated(node))
                numnodes = numnodes + 1
                ! The next are errors. They should never happen
                ! if the list was created with insert_ll_2d.
                if (.not. allocated(node%val)) then
                        write(*,*) "Error: ", __FILE__, ":", __LINE__
                        return
                end if
                if (size(node%val,1) == 0) then
                        write(*,*) "Error: ", __FILE__, ":", __LINE__
                        return
                end if
                if (.not. allocated(node%logi)) then
                        write(*,*) "Error: ", __FILE__, ":", __LINE__
                        return
                end if
                if (size(node%logi,1) == 0) then
                        write(*,*) "Error: ", __FILE__, ":", __LINE__
                        return
                end if
                numvals = numvals + size(node%val,1)
                node => node%next
        end do
        allocate(idx(numnodes+1), vals(numvals), logi(numvals))
        node => ll
        i = 0
        idx(1) = 1
        do while(associated(node))
                i = i + 1
                idx(i+1) = idx(i) + size(node%val,1)
                vals(idx(i):idx(i+1)-1) = node%val(:)
                logi(idx(i):idx(i+1)-1) = node%logi(:)
                deallocate(node%val, node%logi)
                prev => node
                node => node%next
                deallocate(prev)
        end do
        nullify(ll)
end subroutine



! SUBROUTINE
! reduce_csr: returns a CSR format list where the elements in MASK 
!             are removed.
! PARAMETERS
! idx:        pointers of the csr.
! vals:       values of the csr.
! mask:       mask to filter.

subroutine reduce_csrs(idx, vals, mask)
        integer, dimension(:), allocatable, intent(inout) :: idx, vals
        logical*1, dimension(:), intent(in) :: mask
        integer, dimension(:), allocatable :: tidx
        logical*1, dimension(:), allocatable :: tmask
        integer :: i, j, nidx
        nidx = count(mask .EQV. .FALSE.)
        allocate(tidx(nidx+1))
        allocate(tmask(size(vals,1)))
        tmask = .TRUE.
        tidx(1) = 1
        j = 1
        do i = 1, size(mask,1)
                if ( .NOT. mask(i) ) then
                        tidx(j+1) = tidx(j) + (idx(i+1)-idx(i))
                        tmask(idx(i):idx(i+1)-1) = .FALSE.
                        j = j + 1
                end if
        end do
        deallocate(idx)
        allocate(idx(size(tidx,1)))
        idx = tidx
        deallocate(tidx)
        call reduce(vals, tmask)
        deallocate(tmask)
end subroutine

subroutine reduce_csrc(idx, vals, logi, mask)
        integer, dimension(:), allocatable, intent(inout) :: idx, vals
        logical*1, dimension(:), allocatable, intent(inout) :: logi
        logical*1, dimension(:), intent(in) :: mask
        integer, dimension(:), allocatable :: tidx
        logical*1, dimension(:), allocatable :: tmask
        integer :: i, j, nidx
        nidx = count(mask .EQV. .FALSE.)
#ifdef DEBUG
        write(*,*) "Count is: ",nidx , "Sizes are: idx: ", size(idx), " vals: ",&
                   size(vals), " logi: ", size(logi)
#endif
        allocate(tidx(nidx+1))
        allocate(tmask(size(vals,1)))
        tmask = .TRUE.
        tidx(1) = 1
        j = 1
        do i = 1, size(mask,1)
                if ( .NOT. mask(i) ) then
                        tidx(j+1) = tidx(j) + (idx(i+1)-idx(i))
                        tmask(idx(i):idx(i+1)-1) = .FALSE.
                        j = j + 1
                end if
        end do
        deallocate(idx)
        allocate(idx(size(tidx,1)))
        idx = tidx
        deallocate(tidx)
        call reduce(logi, tmask)
        call reduce(vals, tmask)
#ifdef DEBUG
        write(*,*) "Sizes are: idx: ", size(idx), " vals: ", size(vals), " logi: ", size(logi)
#endif

        deallocate(tmask)
end subroutine

! SUBROUTINE
! concat_csr: concatenates 2 csr format arrays into one.
! PARAMETERS
! idxa:       the indices array that is increased.
! valsa:      values of the csr that are increased.
! idxa:       the indices array that is concatenated and deleted afterwards.
! valsa:      the values of the array that are concatenated and then deleted.

subroutine concat_csrs(idxa, valsa, idxb, valsb)
        integer, dimension(:), allocatable, intent(inout) :: idxa, idxb, &
                                                             valsa, valsb
        integer, dimension(:), allocatable :: idxc, valsc
        allocate(idxc(size(idxa,1)+size(idxb,1)-1))
        allocate(valsc(size(valsa,1)+size(valsb,1)))
        idxc(1:size(idxa,1))=idxa(:)
        idxc(size(idxa,1)+1:)=idxb(2:)+(idxa(size(idxa,1))-1)
        valsc(1:size(valsa,1)) = valsa(:)
        valsc(size(valsa,1)+1:) = valsb(:)
        deallocate(idxa, idxb, valsa, valsb)
        allocate(idxa(size(idxc,1)), valsa(size(valsc,1)))
        idxa = idxc
        valsa = valsc
        deallocate(idxc, valsc)
end subroutine

subroutine concat_csrc(idxa, valsa, logia, idxb, valsb, logib)
        integer, dimension(:), allocatable, intent(inout) :: idxa, idxb, &
                                                             valsa, valsb
        logical*1, dimension(:), allocatable, intent(inout) :: logia, logib
        integer, dimension(:), allocatable :: idxc, valsc
        logical*1, dimension(:), allocatable :: logic
        
        allocate(idxc(size(idxa,1)+size(idxb,1)-1))
        allocate(valsc(size(valsa,1)+size(valsb,1)))
        allocate(logic(size(logia,1)+size(logib,1)))
#ifdef DEBUG
        write(*,*) "Sizes are: idx: ", size(idxa), size(idxb), " vals: ", &
                   size(valsa), size(valsb), " logi: ", size(logia), size(logib)
#endif
       

        idxc(1:size(idxa,1))=idxa(:)
        idxc(size(idxa,1)+1:)=idxb(2:)+(idxa(size(idxa,1))-1)
        valsc(1:size(valsa,1)) = valsa(:)
        valsc(size(valsa,1)+1:) = valsb(:)
        logic(1:size(logia,1)) = logia(:)
        logic(size(logia,1)+1:) = logib(:)
        deallocate(idxa)
        deallocate(idxb)
        deallocate(valsa)
        deallocate(valsb) 
        deallocate(logia)
        deallocate(logib)
        allocate(idxa(size(idxc,1)), valsa(size(valsc,1)))
        allocate(logia(size(logic,1)))
        idxa = idxc
        valsa = valsc
        logia = logic
#ifdef DEBUG
        write(*,*) "Sizes are: idxa: ", size(idxa), " valsa: ", &
                   size(valsa), " logia: ", size(logia)
#endif
        deallocate(idxc)
        deallocate(valsc)
        deallocate(logic)
end subroutine
! SUBROUTINE
! concat[x]:  concatenates 2 arrays into one.
! PARAMETERS
! arra:       first and returned array.
! arrb:       second array, gets added to the first and deleted.

subroutine concati(arra, arrb)
        integer, dimension(:), allocatable, intent(inout) :: arra, arrb
        integer, dimension(:), allocatable :: tmp
        allocate(tmp(size(arra,1)+size(arrb,1)))
        tmp(:size(arra,1)) = arra
        tmp(size(arra,1)+1:) = arrb
        deallocate(arra, arrb)
        allocate(arra(size(tmp,1)))
        arra = tmp
        deallocate(tmp)
end subroutine

subroutine concatr(arra, arrb)
        real*8, dimension(:), allocatable, intent(inout) :: arra, arrb
        real*8, dimension(:), allocatable :: tmp
        allocate(tmp(size(arra,1)+size(arrb,1)))
        tmp(:size(arra,1)) = arra
        tmp(size(arra,1)+1:) = arrb
        deallocate(arra, arrb)
        allocate(arra(size(tmp,1)))
        arra = tmp
        deallocate(tmp)
end subroutine

subroutine concati2(arra, arrb)
        integer, dimension(:,:), allocatable, intent(inout) :: arra, arrb
        integer, dimension(:,:), allocatable :: tmp
        if ( size(arra,1) /= size(arrb,1) ) then
                return
        end if
        allocate(tmp(size(arra,1),size(arra,2)+size(arrb,2)))
        tmp(:,:size(arra,2)) = arra(:,:)
        tmp(:,size(arra,2)+1:) = arrb(:,:)
        deallocate(arra, arrb)
        allocate(arra(size(tmp,1),size(tmp,2)))
        arra = tmp
        deallocate(tmp)
end subroutine

subroutine concatr2(arra, arrb)
        real*8, dimension(:,:), allocatable, intent(inout) :: arra, arrb
        real*8, dimension(:,:), allocatable :: tmp
        if ( size(arra,1) /= size(arrb,1) ) then
                return
        end if
        allocate(tmp(size(arra,1),size(arra,2)+size(arrb,2)))
        tmp(:,:size(arra,2)) = arra(:,:)
        tmp(:,size(arra,2)+1:) = arrb(:,:)
        deallocate(arra, arrb)
        allocate(arra(size(tmp,1),size(tmp,2)))
        arra = tmp
        deallocate(tmp)
end subroutine

! SUBROUTINE
! reduce[x]:  reduces the array according to a mask.
! PARAMETERS
! arra:       array to reduce
! mask:       true if the element is deleted.

subroutine reducei(arra, mask)
        integer, dimension(:), allocatable, intent(inout) :: arra
        logical*1, dimension(:), intent(in) :: mask        
        integer, dimension(:), allocatable :: tmp
        integer :: n
        n = count(mask .EQV. .FALSE.)
        allocate(tmp(n))
        tmp = pack(arra, mask .EQV. .FALSE.)
        deallocate(arra)
        allocate(arra(n))
        arra = tmp
        deallocate(tmp)
end subroutine

subroutine reducer(arra, mask)
        real*8, dimension(:), allocatable, intent(inout) :: arra
        logical*1, dimension(:), intent(in) :: mask        
        real*8, dimension(:), allocatable :: tmp
        integer :: n
        n = count(mask .EQV. .FALSE.)
        allocate(tmp(n))
        tmp = pack(arra, mask .EQV. .FALSE.)
        deallocate(arra)
        allocate(arra(n))
        arra = tmp
        deallocate(tmp)
end subroutine

subroutine reducel(arra, mask)
        logical*1, dimension(:), allocatable, intent(inout) :: arra
        logical*1, dimension(:), intent(in) :: mask        
        logical*1, dimension(:), allocatable :: tmp
        integer :: n
        n = count(mask .EQV. .FALSE.)
        allocate(tmp(n))
        tmp = pack(arra, mask .EQV. .FALSE.)
        deallocate(arra)
        allocate(arra(n))
        arra = tmp
        deallocate(tmp)
end subroutine

subroutine reducei2(arra, mask)
        integer, dimension(:,:), allocatable, intent(inout) :: arra
        logical*1, dimension(:), intent(in) :: mask
        integer, dimension(:,:), allocatable :: tmp
        integer :: n, i
        n = count(mask .EQV. .FALSE.)
        allocate(tmp(size(arra,1),n))
        do i = 1, size(arra,1)
                tmp(1,:) = pack(arra(1,:), mask .EQV. .FALSE.)
        end do
        deallocate(arra)
        allocate(arra(size(tmp,1),n))
        arra = tmp
        deallocate(tmp)
end subroutine

subroutine reducer2(arra, mask)
        real*8, dimension(:,:), allocatable, intent(inout) :: arra
        logical*1, dimension(:), intent(in) :: mask
        real*8, dimension(:,:), allocatable :: tmp
        integer :: n, i
        n = count(mask .EQV. .FALSE.)
        allocate(tmp(size(arra,1),n))
        do i = 1, size(arra,1)
                tmp(i,:) = pack(arra(i,:), mask .EQV. .FALSE.)
        end do
        deallocate(arra)
        allocate(arra(size(tmp,1),n))
        arra = tmp
        deallocate(tmp)
end subroutine

! SUBROUTINE
! sort[x]:        sorts a two-dimensional array using an array of  indices
! PARAMETERS
! ind:          indices to sort.
! right:        the array sorted by the indices.

recursive subroutine sorti(ind, right)
        integer, dimension(:), intent(inout) :: ind
        integer, dimension(:,:), intent(inout) :: right
        integer :: i, pivote, valpiv, siz
        siz = size(ind)
        if ( siz < 2 ) then
                return
        end if
        pivote = siz/2
        valpiv = ind(pivote)
        ind((/pivote, siz/)) = ind((/siz, pivote/))
        right(:,(/pivote, siz/)) = right(:,(/siz, pivote/))
        pivote = 1
        do i=1,siz
                if ( ind(i) < valpiv ) then
                        ind((/pivote, i/)) = ind((/i, pivote/))
                        right(:,(/pivote, i/)) = right(:,(/i, pivote/))
                        pivote = pivote + 1
                end if
        end do
        ind((/pivote, siz/)) = ind((/siz, pivote/))
        right(:,(/pivote, siz/)) = right(:,(/siz, pivote/))
        if ( pivote < siz ) then
                call sorti(ind((pivote+1):siz), right(:,(pivote+1):siz))
        end if
        if ( 1 < pivote ) then
                call sorti(ind(1:(pivote-1)), right(:,1:(pivote-1)))
        end if
end subroutine

recursive subroutine sorti2(ind, right)
        integer, dimension(:,:), intent(inout) :: ind
        integer, dimension(:,:), intent(inout) :: right
        integer, dimension(:), allocatable :: valpiv
        allocate(valpiv(size(ind,1)))
        call sorti2imp(ind, right, valpiv)
        deallocate(valpiv)
end subroutine

recursive subroutine sorti2imp(ind, right, valpiv)
        integer, dimension(:,:), intent(inout) :: ind
        integer, dimension(:,:), intent(inout) :: right
        integer, dimension(:), intent(inout) :: valpiv
        integer :: i, pivote, siz
        siz = size(ind,2)
        if ( siz < 2 ) then
                return
        end if
        pivote = siz/2
        valpiv = ind(:,pivote)
        ind(:, (/pivote, siz/)) = ind(:,(/siz, pivote/))
        right(:,(/pivote, siz/)) = right(:,(/siz, pivote/))
        pivote = 1
        do i=1,siz
                if ( arrayunder(ind(:,i), valpiv(:)) ) then
                        ind(:,(/pivote, i/)) = ind(:,(/i, pivote/))
                        right(:,(/pivote, i/)) = right(:,(/i, pivote/))
                        pivote = pivote + 1
                end if
        end do
        ind(:,(/pivote, siz/)) = ind(:,(/siz, pivote/))
        right(:,(/pivote, siz/)) = right(:,(/siz, pivote/))
        if ( pivote < siz ) then
                call sorti2imp(ind(:,(pivote+1):siz),&
                               right(:,(pivote+1):siz),valpiv)
        end if
        if ( 1 < pivote ) then
                call sorti2imp(ind(:,1:(pivote-1)),&
                               right(:,1:(pivote-1)),valpiv)
        end if
end subroutine

recursive subroutine sorti2ir(ind, righti, rightr)
        integer, dimension(:,:), intent(inout) :: ind
        integer, dimension(:),   intent(inout) :: righti
        real*8,  dimension(:,:), intent(inout) :: rightr
        integer :: i, pivote, siz
        integer, dimension(:), allocatable :: valpiv
        allocate(valpiv(size(ind,1)))
        siz = size(ind,2)
        if ( siz < 2 ) then
                return
        end if
        pivote = siz/2
        valpiv = ind(:,pivote)
        ind(:, (/pivote, siz/)) = ind(:,(/siz, pivote/))
        righti((/pivote, siz/)) = righti((/siz, pivote/))
        rightr(:,(/pivote, siz/)) = rightr(:,(/siz, pivote/))
        pivote = 1
        do i=1,siz
                if ( arrayunder(ind(:,i), valpiv(:)) ) then
                        ind(:,(/pivote, i/)) = ind(:,(/i, pivote/))
                        righti((/pivote, i/)) = righti((/i, pivote/))
                        rightr(:,(/pivote, i/)) = rightr(:,(/i, pivote/))
                        pivote = pivote + 1
                end if
        end do
        deallocate(valpiv)
        ind(:,(/pivote, siz/)) = ind(:,(/siz, pivote/))
        righti((/pivote, siz/)) = righti((/siz, pivote/))
        rightr(:,(/pivote, siz/)) = rightr(:,(/siz, pivote/))
        if ( pivote < siz ) then
                call sorti2ir(ind(:,(pivote+1):siz), righti((pivote+1):siz),&
                              rightr(:,(pivote+1):siz))
        end if
        if ( 1 < pivote ) then
                call sorti2ir(ind(:,1:(pivote-1)), righti(1:(pivote-1)),&
                              rightr(:,1:(pivote-1)))
        end if
end subroutine

subroutine sorts2(ind)
        integer, dimension(:,:), intent(inout) :: ind
        integer, dimension(:), allocatable :: valpiv
        allocate(valpiv(size(ind,1)))
        call sorts2imp(ind, valpiv)
        deallocate(valpiv)
end subroutine

recursive subroutine sorts2imp(ind, valpiv)
        integer, dimension(:,:), intent(inout) :: ind
        integer, dimension(:), intent(inout) :: valpiv
        integer :: i, pivote, siz
        siz = size(ind,2)
        if ( siz < 2 ) then
                return
        end if
        pivote = siz/2
        valpiv = ind(:,pivote)
        ind(:, (/pivote, siz/)) = ind(:,(/siz, pivote/))
        pivote = 1
        do i=1,siz
                if ( arrayunder(ind(:,i), valpiv(:)) ) then
                        ind(:,(/pivote, i/)) = ind(:,(/i, pivote/))
                        pivote = pivote + 1
                end if
        end do
        ind(:,(/pivote, siz/)) = ind(:,(/siz, pivote/))
        if ( pivote < siz ) then
                call sorts2imp(ind(:,(pivote+1):siz), valpiv)
        end if
        if ( 1 < pivote ) then
                call sorts2imp(ind(:,1:(pivote-1)), valpiv)
        end if

end subroutine

recursive subroutine sortr(ind, right)
        integer, dimension(:), intent(inout) :: ind
        real*8, dimension(:,:), intent(inout) :: right
        integer :: i, pivote, valpiv, siz
        siz = size(ind)
        if ( siz < 2 ) then
                return
        end if
        pivote = siz/2
        valpiv = ind(pivote)
        ind((/pivote, siz/)) = ind((/siz, pivote/))
        right(:,(/pivote, siz/)) = right(:,(/siz, pivote/))
        pivote = 1
        do i=1,siz
                if ( ind(i) < valpiv ) then
                        ind((/pivote, i/)) = ind((/i, pivote/))
                        right(:,(/pivote, i/)) = right(:,(/i, pivote/))
                        pivote = pivote + 1
                end if
        end do
        ind((/pivote, siz/)) = ind((/siz, pivote/))
        right(:,(/pivote, siz/)) = right(:,(/siz, pivote/))
        if ( pivote < siz ) then
                call sortr(ind((pivote+1):siz), right(:,(pivote+1):siz))
        end if
        if ( 1 < pivote ) then
                call sortr(ind(1:(pivote-1)), right(:,1:(pivote-1)))
        end if
end subroutine

recursive subroutine sortr2(ind, right)
        integer, dimension(:,:), intent(inout) :: ind
        real*8, dimension(:,:), intent(inout) :: right
        integer, dimension(:), allocatable :: valpiv
        allocate(valpiv(size(ind,1)))
        call sortr2imp(ind, right, valpiv)
        deallocate(valpiv)
end subroutine


recursive subroutine sortr2imp(ind, right, valpiv)
        integer, dimension(:,:), intent(inout) :: ind
        real*8, dimension(:,:), intent(inout) :: right
        integer, dimension(:), intent(inout) :: valpiv
        integer :: i, pivote, siz
        siz = size(ind,2)
        if ( siz < 2 ) then
                return
        end if
        pivote = siz/2
        valpiv = ind(:,pivote)
        ind(:, (/pivote, siz/)) = ind(:,(/siz, pivote/))
        right(:,(/pivote, siz/)) = right(:,(/siz, pivote/))
        pivote = 1
        do i=1,siz
                if ( arrayunder(ind(:,i), valpiv(:)) ) then
                        ind(:,(/pivote, i/)) = ind(:,(/i, pivote/))
                        right(:,(/pivote, i/)) = right(:,(/i, pivote/))
                        pivote = pivote + 1
                end if
        end do
        ind(:,(/pivote, siz/)) = ind(:,(/siz, pivote/))
        right(:,(/pivote, siz/)) = right(:,(/siz, pivote/))
        if ( pivote < siz ) then
                call sortr2imp(ind(:,(pivote+1):siz),  &
                               right(:,(pivote+1):siz),&
                               valpiv)
        end if
        if ( 1 < pivote ) then
                call sortr2imp(ind(:,1:(pivote-1)),  &
                            right(:,1:(pivote-1)),&
                            valpiv)
        end if
end subroutine


! SUBROUTINE
! sorts:        sorts a one dimensional integer array array.
! PARAMETERS
! ind:          array to sort.

recursive subroutine sorts(ind)
        integer, dimension(:), intent(inout) :: ind
        integer :: i, pivote, valpiv, siz
        siz = size(ind)
        if ( siz < 2 ) then
                return
        end if
        pivote = siz/2
        valpiv = ind(pivote)
        ind((/pivote, siz/)) = ind((/siz, pivote/))
        pivote = 1
        do i=1,siz
                if ( ind(i) < valpiv ) then
                        ind((/pivote, i/)) = ind((/i, pivote/))
                        pivote = pivote + 1
                end if
        end do
        ind((/pivote, siz/)) = ind((/siz, pivote/))
        if ( pivote < siz ) then
                call sorts(ind((pivote+1):siz))
        end if
        if ( 1 < pivote ) then
                call sorts(ind(1:(pivote-1)))
        end if
end subroutine

! SUBROUTINE
! binsearch[x]:  searches into a sorted list for a set of elements and returns the
!             indices. 
! PARAMETERS
! list:       Sorted list.
! values:     The values to search in the list. it's an array of any size.
! indices:    The returned indices for the values searched.
! found:      An array of logicals. True if the indices exist.

subroutine binsearch1(list, values, indices)
        implicit none
        integer, dimension(:), intent(in)    :: list
        integer, dimension(:), intent(in)    :: values
        integer, dimension(:), intent(inout) :: indices
        integer :: n, nval, i
        integer, allocatable, dimension(:) :: maxim, minim
        n = size(list,1)
        nval = size(values,1)
        allocate(maxim(nval), minim(nval))
        minim(:) = 1
        maxim(:) = n
        do
                do i = 1, nval
                        if ( (maxim(i)-minim(i)) <=1 ) then
                                if ( list(minim(i)) == values(i) ) then
                                        maxim(i) = minim(i)
                                else if ( list(maxim(i)) == values(i) ) then
                                        minim(i) = maxim(i)
                                else
                                        minim(i) = 0 
                                        maxim(i) = 0
                                end if
                        end if
                end do
                indices(:) = ( maxim(:) + minim(:) ) / 2
                if ( all(minim(:) == maxim(:)) ) then
                        exit
                end if
                do i = 1, nval
                        if ( indices(i) == 0 ) then
                                cycle
                        end if
                        if ( list(indices(i)) < values(i) ) then
                                minim(i) = indices(i) + 1
                        else if (  values(i) < list(indices(i)) ) then
                                maxim(i) = indices(i) - 1
                        else
                                minim(i) = indices(i)
                                maxim(i) = indices(i)
                        end if
                end do
        end do
        deallocate(maxim, minim)
end subroutine

subroutine binsearch2(list, values, indices)
        implicit none
        integer, dimension(:,:), intent(in)    :: list
        integer, dimension(:,:), intent(in)    :: values
        integer, dimension(:), intent(inout) :: indices
        integer :: n, nval, i
        integer, allocatable, dimension(:) :: maxim, minim
        n = size(list,2)
        nval = size(values,2)
        allocate(maxim(nval), minim(nval))
        minim(:) = 1
        maxim(:) = n
        do
                do i = 1, nval
                        if ( (maxim(i)-minim(i)) <=1 ) then
                                if ( all(list(:,minim(i)) == values(:,i)) ) then
                                        maxim(i) = minim(i)
                                else if ( all(list(:,maxim(i)) == values(:,i)) ) then
                                        minim(i) = maxim(i)
                                else
                                        minim(i) = 0 
                                        maxim(i) = 0
                                end if
                        end if
                end do
                indices(:) = ( maxim(:) + minim(:) ) / 2
                if ( all(minim(:) == maxim(:)) ) then
                        exit
                end if
                do i = 1, nval
                        if ( indices(i) == 0 ) then
                                cycle
                        end if
                        if ( arrayunder(list(:,indices(i)), values(:,i)) ) then
                                minim(i) = indices(i) + 1
                        else if ( arrayunder(values(:,i), list(:,indices(i))) ) then
                                maxim(i) = indices(i) - 1
                        else
                                minim(i) = indices(i)
                                maxim(i) = indices(i)
                        end if
                end do
        end do
        deallocate(maxim, minim)
end subroutine

subroutine binsearch2s(list, val, ind)
        implicit none
        integer, dimension(:,:), intent(in)    :: list
        integer, dimension(:), intent(in)    :: val
        integer, intent(out) :: ind
        integer :: n
        integer :: maxim, minim
        n = size(list,2)
        minim = 1
        maxim = n
        do
                if ( (maxim-minim) <=1 ) then
                        if ( all( list(:,minim) == val(:) ) ) then
                                ind = minim
                        else if ( all( list(:,maxim) == val(:) ) ) then
                                ind = maxim
                        else
                                ind = 0
                        end if
                        return
                end if
                ind = ( maxim + minim ) / 2
                if ( arrayunder(list(:,ind), val(:)) ) then
                        minim = ind + 1
                else if ( arrayunder(val(:), list(:,ind)) ) then
                        maxim = ind - 1
                else
                        return
                end if
        end do
end subroutine

subroutine binsearch1s(list, val, ind)
        integer, dimension(:), intent(in)    :: list
        integer, intent(in) :: val
        integer, intent(out) :: ind
        integer :: n
        integer :: maxim, minim
        n = size(list,1)
        minim = 1
        maxim = n
        do
                if ( (maxim-minim) <=1 ) then
                        if ( list(minim) == val )  then
                                ind = minim
                        else if ( list(maxim) == val )  then
                                ind = maxim
                        else
                                ind = 0
                        end if
                        return
                end if
                ind = ( maxim + minim ) / 2
                if ( list(ind) < val ) then
                        minim = ind + 1
                else if ( val < list(ind)) then
                        maxim = ind - 1
                else
                        return
                end if
        end do

end subroutine

subroutine reverse_indices(list, values)
        integer, dimension(:), intent(in)    :: list
        integer, dimension(:), intent(inout) :: values
        integer, dimension(:), allocatable   :: indices
        allocate(indices(size(values,1)))
        call binsearch(list, values, indices)
! #ifdef DEBUG
        if ( count( indices == 0 ) > 0 ) then
                write(*,*) "Not found: ", count(indices == 0)
        end if
! #endif
        values(:) = indices(:)
        deallocate(indices)
end subroutine

subroutine get_uniquec(list, ret)
        integer, dimension(:), intent(in) :: list
        integer, dimension(:), allocatable, intent(out) :: ret
        type(linked_list_1d), pointer :: lsnodes
        integer :: i
        nullify(lsnodes)
        if ( allocated(ret) ) then
                deallocate(ret)
        end if
        allocate(ret(size(list,1)))
        ret = list
        call sort(ret)
        if ( size(ret,1) >= 1 ) then
                call insertll_1d(lsnodes, ret(1))
        end if
        do i = 2, size(ret,1)
                if ( ret(i) /= ret(i-1) ) then
                        call insertll_1d(lsnodes, ret(i))
                end if
        end do
        deallocate(ret)
        call toarray_1d(lsnodes, ret)
end subroutine

subroutine get_uniques(list)
        integer, dimension(:), allocatable, intent(inout) :: list
        type(linked_list_1d), pointer :: lsnodes
        integer :: i
        nullify(lsnodes)
        call sort(list)
        if ( allocated(list) ) then
                call insertll_1d(lsnodes, list(1))
        end if
        do i = 2, size(list,1)
                if ( list(i) /= list(i-1) ) then
                        call insertll_1d(lsnodes, list(i))
                end if
        end do
        deallocate(list)
        call toarray(lsnodes, list)
end subroutine

! SUBROUTINE
! factorize:  factorizes a number into the two largest numbers that multiplied 
!             give as a result that number. Not the fastest algorithm.
!             if you have a number of 128 bits of processors, the least of 
!             your problems is factorization.
! PARAMETERS
! tofactor:   A positive integer to factorize.
! ret:        Two positive numbers such that: ret(1)*ret(2) = tofactor.

subroutine factorize(tofactor, ret)
        integer, intent(in) :: tofactor
        integer, dimension(2), intent(out) :: ret
        type(linked_list_1d), pointer :: ll 
        integer, dimension(:), allocatable :: factors
        integer :: num, div, i
        nullify(ll)
        num = tofactor
        div = 2
        do while( num /= 1 )
                if ( mod(num, div) == 0 ) then
                        call insertll_1d(ll, div)
                        num = num/div
                else
                        div = div + 1
                end if
        end do
        call toarray_1d(ll, factors)
        ret(1) = 1
        ret(2) = 1
        if ( allocated(factors) ) then
                do i = 1, size(factors,1)
                        ret(2-mod(i,2)) = ret(2-mod(i,2))*factors(i)
                end do
        deallocate(factors)
        end if
end subroutine

subroutine error_exfem(msg)
        character(len=*), intent(in) :: msg
        integer :: ierr
        write(*,*) msg
        call mpi_abort(MPI_COMM_WORLD, 1, ierr)
        call exit(1)
end subroutine

function isinf(r)
        real*8, intent(in) :: r
        logical*1 :: isinf
        if ( r*0.0D0 .ne. 0.0D0 ) then
                isinf = .true.
        else
                isinf = .false.
        end if
end function

end module
