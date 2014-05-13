define(`name_func',`$1_$2_$3$4')
define(`name_sele',`name_func(`share',`$1',`$2',ifelse($3,`1',,`2'))')
define(`share_func',`subroutine name_sele($1,$2,$4)(outt, innt)
type($1), dimension(:), intent(inout) :: outt, innt
integer, dimension(:), allocatable  :: inn, outn
integer, dimension(:), allocatable    :: reqs
integer, dimension(:,:), allocatable :: stats
integer, parameter :: TAG_DATA  = 1
integer :: i, reqc, ierr
logical :: flag
allocate(inn(c_size), outn(c_size))
outn = 0
inn = 0
do i = 1, c_size
if ( allocated(outt(i)%$2) ) then
ifelse(`1',`$4',`outn(i) = size(outt(i)%$2,1)',`outn(i) = size(outt(i)%$2,2)')
end if
end do
call share_num(outn, inn)
allocate(reqs(2*(c_size-1)))
allocate(stats(MPI_STATUS_SIZE, 2*(c_size-1))) 
reqs = 0
stats = 0
reqc = 0
do i = 0, (c_size-1)
if ( i /= c_rank ) then
if ( inn(i+1) > 0 ) then
ifelse(`1',`$4',`allocate(innt(i+1)%$2(inn(i+1)))',`allocate(innt(i+1)%$2($4,inn(i+1)))')
reqc = reqc + 1
call mpi_irecv(innt(i+1)%$2,&
               $4*inn(i+1),&  
               MPI_$3,i,TAG_DATA,&
               MPI_COMM_WORLD,&
               reqs(reqc),ierr)
call check_mpi(ierr,__LINE__)
end if
if ( outn(i+1) > 0 ) then
reqc = reqc + 1
call mpi_isend(outt(i+1)%$2,&
               $4*outn(i+1),&
               MPI_$3,i,TAG_DATA,&
               MPI_COMM_WORLD,&
               reqs(reqc),ierr)
call check_mpi(ierr,__LINE__)
end if
end if
end do
call mpi_waitall(reqc,reqs,stats,ierr)
call check_mpi(ierr,__LINE__)
do i = 1,reqc
call mpi_test_cancelled(stats(:,i),flag,ierr)
call check_mpi(ierr,__LINE__)
if ( flag ) then
        write(*,*) "Cancelled communication"
end if
end do
end subroutine')

share_func(`edges',`links',`INTEGER',`2')
share_func(`edges',`hanging',`LOGICAL1',`1')
share_func(`elemcsr', `idx', `INTEGER', `1')
share_func(`elemcsr', `vals', `INTEGER', `1')
share_func(`elemcsr', `logi', `LOGICAL1', `1')
share_func(`elemcsr', `prop', `INTEGER', `1')
share_func(`elemcsr', `gnode', `INTEGER', `1')
share_func(`elemcsr', `nodes', `REAL8', `2')
share_func(`elemedge', `idx', `INTEGER', `1')
share_func(`elemedge', `vals', `INTEGER', `1')
share_func(`elemedge', `links', `INTEGER', `2')
share_func(`division', `underedges', `INTEGER', `2')
share_func(`division', `overgnodes', `INTEGER', `1')
share_func(`division', `overedges', `INTEGER', `2')
share_func(`division', `overpoints', `REAL8', `2')
share_func(`shnodes', `hgnodes', `INTEGER', `1')
share_func(`shnodes', `annodes', `INTEGER', `1')

