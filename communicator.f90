module communicator

        use model
        use mesh_utils
        use io_utils
        use iso_c_binding

        implicit none

        type neighbor
                ! Shared nodes between communicators.
                integer :: reqsnd
                integer :: reqrcv
                integer, dimension(MPI_STATUS_SIZE) :: statsnd
                integer, dimension(MPI_STATUS_SIZE) :: statrcv
                integer, dimension(:), allocatable  :: locdofs
                integer, dimension(:), allocatable  :: gnodes
                real*8,  dimension(:), allocatable  :: outvalues
                real*8,  dimension(:), allocatable  :: invalues
                integer, dimension(:), allocatable  :: inhdofs
                integer, dimension(:), allocatable  :: outhdofs
                real*8,  dimension(:), allocatable  :: inhvalues
                real*8,  dimension(:), allocatable  :: outhvalues
                integer, dimension(:), allocatable  :: inadofs
                integer, dimension(:), allocatable  :: outadofs
                real*8, dimension(:), allocatable   :: inavalues
                real*8, dimension(:), allocatable   :: outavalues
        end type neighbor

        type neighborll
                type(linked_list_1d), pointer :: gnodes => null() 
        end type

        integer     :: c_rank
        integer     :: c_size
        integer, parameter :: MIN_ELEMS = 10 
        integer, dimension(:), allocatable :: c_neighlist
        type(c_ptr) :: c_comm_world
        type(neighbor), dimension(:), allocatable :: c_neighbors

        !Interfaces to access C functions from MPI and ParMETIS that we
        !need to use.
        interface
                function mpi_comm_f2c(fcomm) bind(C, name="MPI_Comm_f2c")
                        use iso_c_binding
                        type(c_ptr)    :: mpi_comm_f2c
                        integer(c_int) :: fcomm
                end function mpi_comm_f2c
        end interface
        ! parmetis mesh partition
        interface
                function parmetis_v3_partmeshkway(elmdist, eptr, eind,          &
                                                    elmwgt, wgtflag, numflag,   &
                                                    ncon, ncommonnodes, nparts, &
                                                    tpwgts, ubvec, options,     &
                                                    edgecut, part, comm)
                        use iso_c_binding
                        integer     :: parmetis_v3_partmeshkway
                        integer     :: elmdist(*)
                        integer     :: eptr(*)
                        integer     :: eind(*)
                        type(c_ptr) :: elmwgt(*)
                        integer     :: wgtflag(*)
                        integer     :: numflag(*)
                        integer     :: ncon(*)
                        integer     :: ncommonnodes(*)
                        integer     :: nparts(*)
                        real        :: tpwgts(*)
                        real        :: ubvec(*)
                        integer     :: options(*)
                        integer     :: edgecut(*)
                        integer     :: part(*)
                        type(c_ptr) :: comm
                end function
        end interface
        ! parmetis mesh dual
        interface
                function parmetis_v3_mesh2dual(elmdist, eptr, eind,                  &
                                               numflag, ncommonnodes,                &
                                               xadj, adjncy, comm)
                        use iso_c_binding
                        integer(c_int) :: parmetis_v3_mesh2dual
                        integer(c_int) :: elmdist(*)
                        integer(c_int) :: eptr(*)
                        integer(c_int) :: eind(*)
                        integer(c_int) :: numflag(*)
                        integer(c_int) :: ncommonnodes(*)
                        type(c_ptr)    :: xadj
                        type(c_ptr)    :: adjncy
                        type(c_ptr)    :: comm
                end function
        end interface

        ! metis default options (not parallel)

        interface 
                function metis_setdefaultoptions(options)
                        use iso_c_binding
                        integer(c_int) :: metis_setdefaultoptions
                        integer(c_int) :: options(*)
                end function
        end interface

        ! metis partition (not parallel)
        interface
                function metis_partmeshdual(ne, nn, eptr, eind, vwgt, vsize,  &
                                            ncommon, nparts, tpwgts, options, &
                                            objval, epart, npart)
                        use iso_c_binding
                        integer(c_int) :: metis_partmeshdual
                        integer(c_int) :: ne
                        integer(c_int) :: nn
                        integer(c_int) :: eptr(*)
                        integer(c_int) :: eind(*)
                        type(c_ptr), value :: vwgt
                        type(c_ptr), value :: vsize
                        integer(c_int) :: ncommon
                        integer(c_int) :: nparts
                        type(c_ptr), value :: tpwgts
                        integer(c_int) :: options(*)
                        integer(c_int) :: objval
                        integer(c_int) :: epart(*)
                        integer(c_int) :: npart(*)
                end function
        end interface



        ! frees metis pointer
        interface
                subroutine metis_free(ptr)
                        use iso_c_binding
                        type(c_ptr), value :: ptr
                end subroutine
        end interface

        contains

        ! imports the code generated with m4
#include "comm_temps.f90"


        ! subroutine to check errors in MPI
        subroutine check_mpi(ierr, line)

                integer :: ierr, ierrc, line, length
                character(len=MPI_MAX_ERROR_STRING) :: serror
                ierrc = ierr
                if ( ierr .NE. MPI_SUCCESS ) then
                        call mpi_error_string(ierrc, serror, length, ierr)
                        write(*,*) "Failed with error code: ", ierrc, "line: ", line
                        write(*,*) "Failed with error string: ", serror(1:length)
                        call mpi_abort(MPI_COMM_WORLD, ierrc, ierr)
                        call exit(1)
                end if
        end subroutine


        ! subroutine that starts the main parameters of mpi and initializes
        ! them accordingly
        subroutine init_mpi()
                integer :: ierr
                call mpi_init(ierr)
                call mpi_errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_comm_rank(MPI_COMM_WORLD, c_rank, ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_comm_size(MPI_COMM_WORLD, c_size, ierr)
                call check_mpi(ierr,__LINE__)
                ! C MPI World communicator needed by parmetis functions.
                c_comm_world = mpi_comm_f2c(MPI_COMM_WORLD)
        end subroutine

        ! Calls mpi_finalize to end the mpi routines.
        subroutine end_mpi()
                integer :: ierr
                call mpi_finalize(ierr)
        end subroutine

        ! Stops the program until gdb is used to debug it and from gdb the flags to
        ! continue are set.
        subroutine wait_debug()
                logical*1 :: st, ch
                character(len=MPI_MAX_PROCESSOR_NAME) :: hostname
                integer :: length, ierr, i
                st = .FALSE.
                call mpi_get_processor_name(hostname, length, ierr)
                write(*,*) "Rank = ", c_rank, ", Proc name = ", hostname(1:length), &
                           " Ready for GDB attachment, PID = ", getpid()
                do while(.NOT. st)
                        call sleep(5)
                        do i = 0, (c_size-1)
                                if ( c_rank == i ) then
                                        ch = st
                                end if
                                call mpi_bcast(ch,1,MPI_LOGICAL1,i,MPI_COMM_WORLD,ierr)
                                if ( ch .AND. logical(i /= c_rank, kind=1) ) then
                                        st = .TRUE.
                                end if
                        end do
                        call mpi_barrier(MPI_COMM_WORLD, ierr)
                end do
        end subroutine

        ! Broadcasts the main parameters needed by the simulation. to all nodes from the
        ! 0th node.
        subroutine bcast_problem()
                integer :: ierr, n
                ! Broadcasting initial mesh conditions.
                call mpi_bcast(m_blength,1,MPI_DOUBLE_PRECISION,&
                               0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(m_mlength,1,MPI_DOUBLE_PRECISION,&
                               0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(m_totcols,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(m_totrows,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(m_lastgnode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(m_totelems,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                ! Broadcasting time controls.
                call mpi_bcast(m_nsteps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(m_deltat,1,MPI_DOUBLE_PRECISION,&
                               0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(m_duration,1,MPI_DOUBLE_PRECISION,&
                               0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                ! Broadcasting materials.
                call mpi_bcast(m_nummat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                if ( c_rank /= 0 ) then
                        allocate(m_matconst(5,m_nummat))
                        allocate(m_probconst(2,m_nummat))
                end if
                call mpi_bcast(m_matconst,5*m_nummat,MPI_DOUBLE_PRECISION,&
                               0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(m_probconst,2*m_nummat,MPI_DOUBLE_PRECISION,&
                               0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                ! Broadcasting problem. 
                if ( c_rank == 0) then
                        n = size(m_pparams)
                end if
                call mpi_bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                if ( c_rank /= 0 ) then
                        allocate(m_pparams(n))
                end if
                call mpi_bcast(m_pparams,n,MPI_DOUBLE_PRECISION,&
                               0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(m_prob,6,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                ! Broadcasting geometry.
                n = size(m_gparams)
                call mpi_bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                if ( c_rank /= 0 ) then
                        allocate(m_gparams(n))
                end if
                call mpi_bcast(m_gparams,n,MPI_DOUBLE_PRECISION,&
                               0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(m_geom,6,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                ! Wait until all messages are sent.
                call mpi_barrier(MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)

        end subroutine

        subroutine bcast_environ()
                integer :: ierr
                ! Broadcasting io parameters.
                call mpi_bcast(io_interval,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(io_newsave,1,MPI_LOGICAL1,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(io_fnodlen,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(io_outdir,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(io_fncdlen,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(io_cvmdir,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(io_fnvdlen,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(io_visdir,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(io_fnindlen,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(io_indir,256,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)       

                ! Wait until all messages are sent.
                call mpi_barrier(MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
        end subroutine



        ! Initial divide operation. It finds 2 integer numbers that multiplied 
        ! give the number of nodes in MPI and divides the elements by the number
        ! of columns and rows in the starting model. (before dividing elements).

        subroutine global_initial_mesh()
                integer :: i, j, ind, prev, ierr
                integer :: fcol, lcol, frow, lrow
                integer, dimension(2) :: ppdim, npproc, plus
                integer, dimension(4) :: localranges
                integer, dimension(:), allocatable :: ranges
                if ( c_rank == 0 ) then
                        allocate(ranges(4*c_size))
                        call factorize(c_size, ppdim)
                        if ( m_totcols > m_totrows ) then
                                if ( ppdim(1) < ppdim(2) ) then
                                        ppdim((/1,2/)) = ppdim((/2,1/))
                                end if
                        else
                                if ( ppdim(1) > ppdim(2) ) then
                                        ppdim((/1,2/)) = ppdim((/2,1/))
                                end if
                        end if
                        npproc(1) = (m_totcols-1)/ppdim(1)
                        npproc(2) = (m_totrows-1)/ppdim(2)
                        plus(1)   = mod((m_totcols-1),ppdim(1))
                        plus(2)   = mod((m_totrows-1),ppdim(2))
                        do i = 1, ppdim(1)
                        do j = 1, ppdim(2)
                        ind = ((i-1)*ppdim(2)+j-1)*4
                        if ( i == 1 ) then
                                ranges(ind+1) = 1
                        else
                                prev = ((i-2)*ppdim(2)+j-1)*4
                                ranges(ind+1) = ranges(prev+2)
                        end if
                        if ( j == 1 ) then
                                ranges(ind+3) = 1
                        else
                                prev = ((i-1)*ppdim(2)+j-2)*4
                                ranges(ind+3) = ranges(prev+4)
                        end if
                        if ( i <= plus(1) ) then
                                ranges(ind+2) = ranges(ind+1)+npproc(1)+1
                        else
                                ranges(ind+2) = ranges(ind+1)+npproc(1)
                        end if
                        if ( j <= plus(2) ) then
                                ranges(ind+4) = ranges(ind+3)+npproc(2)+1
                        else
                                ranges(ind+4) = ranges(ind+3)+npproc(2)
                        end if
                        end do
                        end do
                end if
                call mpi_scatter(ranges,4,MPI_INTEGER,localranges, &
                                 4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                fcol = localranges(1)
                lcol = localranges(2)
                frow = localranges(3)
                lrow = localranges(4)
                if ( c_rank == 0 ) then
                        deallocate(ranges)
                end if
                call initial_mesh(frow, lrow, fcol, lcol)
        end subroutine

        ! Creates the element distribution array needed by parmetis for its operations.
        ! The array contains the number of elements in each node in a parmetis way.
        subroutine update_elmdist()
                integer, dimension(:), allocatable :: nelems
                integer :: i, ierr
                if (.NOT. allocated(m_elmdist)) then
                        allocate(m_elmdist(c_size+1))
                end if
                allocate(nelems(c_size))
                nelems(c_rank+1) = m_nelems
                m_elmdist(1) = 1
                do i = 1, c_size
                        call mpi_bcast(nelems(i),1,MPI_INTEGER, &
                                       i-1,MPI_COMM_WORLD,ierr)
                        call check_mpi(ierr,__LINE__)
                        m_elmdist(i+1) = m_elmdist(i) + nelems(i)
                end do
                call mpi_barrier(MPI_COMM_WORLD, ierr)
                call check_mpi(ierr, __LINE__)
                m_totelems = (m_elmdist(c_size+1)-1)
                deallocate(nelems)
        end subroutine

        subroutine prefix_sum_newnodes(newlocnodes)
                integer, intent(inout) :: newlocnodes
                integer, dimension(:), allocatable :: newnodes
                integer :: i, ierr, curr, next
                if ( c_rank == 0 ) then
                        allocate(newnodes(c_size))
                end if
                call mpi_gather(newlocnodes, 1, MPI_INTEGER,&
                                newnodes, 1, MPI_INTEGER,&
                                0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                if ( c_rank == 0 ) then
                        next = m_lastgnode 
                        do i = 1, c_size
                                curr = newnodes(i)
                                newnodes(i) = next
                                next = curr + next
                        end do
                        m_lastgnode = next
                end if
                call mpi_scatter(newnodes, 1, MPI_INTEGER,&
                                 newlocnodes, 1, MPI_INTEGER,&
                                 0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_bcast(m_lastgnode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_barrier(MPI_COMM_WORLD, ierr)
                call check_mpi(ierr,__LINE__)
                if ( c_rank == 0 ) then
                        deallocate(newnodes)
                end if
        end subroutine

        ! copies the array stored in the c pointer to a newly allocated
        ! fortran array.

        subroutine getmetisptr(caddr, larray, sh)
                type(c_ptr), intent(in) :: caddr
                integer, intent(in) :: sh
                integer, dimension(:), allocatable, intent(out) :: larray
                integer, pointer, dimension(:) :: fptr
                call c_f_pointer(caddr, fptr, (/sh/))
                allocate(larray(sh))
                larray(1:sh) = fptr(1:sh)
                call metis_free(caddr)
        end subroutine

        ! Use parmetis to find a partition for every mesh.
        subroutine partmesh(part, eind)
                integer, dimension(:), intent(inout) :: part, eind
                integer, dimension(1) :: wgtflag, numflag, ncon, ncmnnds, nparts, edgecut
                integer, dimension(4) :: options
                real, dimension(1) :: ubvec
                real, dimension(:), allocatable :: tpwgts
                integer :: ierr
                allocate(tpwgts(c_size))
                
                wgtflag = 0
                numflag = 1
                ncon    = 1
                ncmnnds = 2
                nparts  = c_size
                tpwgts  = 1.0/real(c_size, kind=4)
                ubvec   = 1.05
                options(1) = 0
                ierr = parmetis_v3_partmeshkway(m_elmdist, m_eptr, eind,       &
                                                C_NULL_PTR, wgtflag, numflag,  &
                                                ncon, ncmnnds, nparts,         &
                                                tpwgts, ubvec, options,        &
                                                edgecut, part, c_comm_world)
                part = part - 1
                deallocate(tpwgts)
        end subroutine

        subroutine singlepartmesh(epart)
                use iso_c_binding
                integer, dimension(:), intent(inout) :: epart
                integer, dimension(:), allocatable   :: npart
                integer :: ierr, objval
                integer, parameter :: METIS_NUMBERING  = 17
                integer, parameter :: METIS_NUM_OPTS = 40
                integer, dimension(METIS_NUM_OPTS) :: options
                
                ierr = metis_setdefaultoptions(options)
        !        options(METIS_NUMBERING) = 1
                allocate(npart(m_nnodes))
                m_eptr = m_eptr - 1
                m_eind = m_eind - 1
                ierr = metis_partmeshdual(m_nelems, m_nnodes, m_eptr, m_eind, &
                                          C_NULL_PTR, C_NULL_PTR, 2, 2, C_NULL_PTR, &
                                          options, objval, epart, npart)
                m_eptr = m_eptr + 1
                m_eind = m_eind + 1
                deallocate(npart)
        end subroutine

        ! Finds the dual graph of the mesh and stores it into fortran arrays with the
        ! parmetis format.
        subroutine find_dual(numcom, dxadj, dadjncy)
                integer :: numcom
                integer, dimension(:), allocatable :: dxadj, dadjncy
                integer, dimension(:), allocatable :: eind
                type(c_ptr) :: xadj, adjncy
                integer, dimension(1) :: ncmnnds, numflag
                integer :: ierr
                
                allocate(eind(size(m_eind,1)))
                eind(:) = m_gnode(m_eind)
                numflag = 1
                ncmnnds = numcom 
                ierr = parmetis_v3_mesh2dual(m_elmdist, m_eptr, eind, &
                                             numflag, ncmnnds, xadj,  &
                                             adjncy, c_comm_world)
                call getmetisptr(xadj, dxadj, m_nelems+1)
                call getmetisptr(adjncy, dadjncy, dxadj(m_nelems+1)-1)
                deallocate(eind)
        end subroutine

        subroutine find_dual_corresp(xadj, adjncy, corresp)
                integer, dimension(:), intent(in)  :: xadj, adjncy
                integer, dimension(:), allocatable, intent(out) :: corresp
                integer :: i, j, local, proc, ia, ea, ib, eb, inda, indb, glob, indc
                type(elemedgell), dimension(:), allocatable :: outll, inll
                type(elemedge), dimension(:), allocatable :: outedges, inedges
                
                allocate(corresp(size(adjncy,1)))
                allocate(outll(c_size), inll(c_size))
                allocate(outedges(c_size), inedges(c_size))
                corresp(:) = 0
                do i = 1, m_nelems
                        ia = m_eptr(i)
                        ea = m_eptr(i+1)-1
                        do j = xadj(i), (xadj(i+1)-1)
                                call local_element(adjncy(j), local, proc)
                                if ( proc == c_rank ) then
                                        ib = m_eptr(local)
                                        eb = m_eptr(local+1)-1
                                        call edge_connects(m_eind(ia:ea),&
                                                           m_eind(ib:eb),&
                                                           inda, indb)
                                        corresp(j) = inda
                                else
                                        call insertll(outll(proc+1)%elem, m_eind(ia:ea))
                                        call global_element(i, c_rank, glob)
                                        call insertll(outll(proc+1)%link, (/glob, adjncy(j)/))
                                end if
                        end do 
                end do
                do i = 1, c_size
                        call toarray_csr(outll(i)%elem, outedges(i)%idx,&
                                    outedges(i)%vals)
                        call toarray(outll(i)%link, outedges(i)%links)
                        if ( allocated(outedges(i)%vals) ) then
                                outedges(i)%vals = m_gnode(outedges(i)%vals)
                        end if
                end do
                call share_elemedge_idx(outedges, inedges)
                call share_elemedge_vals(outedges, inedges)
                call share_elemedge_links2(outedges, inedges)
                do i = 1, c_size
                        if ( .NOT. allocated(inedges(i)%idx) ) then
                                cycle
                        end if
                        do j = 1, size(inedges(i)%idx)-1
                                ib = inedges(i)%idx(j)
                                eb = inedges(i)%idx(j+1)-1
                                call local_element(inedges(i)%links(2,j), local, proc)
                                if ( proc /= c_rank ) then
                                        write(*,*) "Error: proc ", c_rank,      &
                                                   " received elem: ",          &
                                                   inedges(i)%links(2,j), proc
                                end if
                                ia = m_eptr(local)
                                ea = m_eptr(local+1)-1
                                do indc = xadj(local), xadj(local+1)-1
                                        if ( adjncy(indc) == inedges(i)%links(1,j) ) then
                                                exit
                                        end if
                                end do
                                if ( adjncy(indc) /=  inedges(i)%links(1,j) ) then
                                        write(*,*) "Error: proc ", c_rank,         &
                                                   "The edge does not correspond", &
                                                   inedges(i)%links(:,j)
                                end if
                                call edge_connects(m_gnode(m_eind(ia:ea)),&
                                                   inedges(i)%vals(ib:eb),&
                                                   inda, indb)
                                corresp(indc) = inda
                        end do
                end do
                if ( count(corresp == 0 ) > 0 ) then
                        write(*,*) "Some corresponding were not found. ", count(corresp == 0)
                end if
                do i = 1, c_size
                        if ( allocated(outedges(i)%idx) ) then
                                deallocate(outedges(i)%idx)
                                deallocate(outedges(i)%vals)
                                deallocate(outedges(i)%links)
                        end if
                        if ( allocated(inedges(i)%idx) ) then
                                deallocate(inedges(i)%idx)
                                deallocate(inedges(i)%vals)
                                deallocate(inedges(i)%links)
                        end if
                end do
                deallocate(inll, outll, outedges, inedges)
        end subroutine

        subroutine share_num(outnums, innums)
        integer, dimension(:), intent(in)  :: outnums
        integer, dimension(:), intent(inout) :: innums
        integer :: ierr, i
        do i = 0, (c_size-1)
                call mpi_scatter(outnums, 1, MPI_INTEGER,     &
                                 innums(i+1), 1, MPI_INTEGER, &
                                 i, MPI_COMM_WORLD, ierr)
        end do
end subroutine

subroutine share_elems(outelems, inelems)
        type(elemcsr), dimension(:) :: outelems, inelems
        call share_elemcsr_idx(outelems,inelems)
        call share_elemcsr_vals(outelems,inelems)
        call share_elemcsr_logi(outelems,inelems)
        call share_elemcsr_prop(outelems,inelems)
        call share_elemcsr_gnode(outelems,inelems)
        call share_elemcsr_nodes2(outelems,inelems)
end subroutine
        
! divide
subroutine select_global_divide(adjptr, adjncy, todivide)
        logical*1, dimension(:), intent(inout)     :: todivide
        integer, dimension(:), intent(in) :: adjptr, adjncy
        type(edges), dimension(:), allocatable :: outedges, inedges
        integer :: i
        
        ! We only need to find the nodes that share an edge.
        allocate(outedges(c_size), inedges(c_size))
        call select_local_divide(adjptr, adjncy, c_rank, c_size, outedges, todivide)
        call share_edges_links2(outedges, inedges)
        call check_neighbor_divide(c_rank, c_size, outedges, inedges)
        call share_edges_hanging(outedges, inedges)
        call update_local_divide(c_rank, c_size, outedges, inedges, todivide)
        do i = 1, c_size
                if ( allocated(outedges(i)%links) ) then
                        deallocate(outedges(i)%links)
                end if
                if ( allocated(outedges(i)%hanging) ) then
                        deallocate(outedges(i)%hanging)
                end if
                if ( allocated(inedges(i)%links) ) then
                        deallocate(inedges(i)%links)
                end if
                if ( allocated(inedges(i)%hanging) ) then
                        deallocate(inedges(i)%hanging)
                end if
        end do

        deallocate(outedges, inedges)
end subroutine

subroutine create_global_divide(adjptr, adjncy, corresp, todivide)
        logical*1, dimension(:), intent(inout)     :: todivide
        integer, dimension(:), intent(in) :: adjptr, adjncy, corresp
        real*8, dimension(:,:), allocatable :: centerpts, singlepts, nodes
        integer, dimension(:,:), allocatable :: corredges
        integer, dimension(:), allocatable :: correlems, gnodes
        logical*1, dimension(:), allocatable :: ishanging
        type(division), dimension(:), allocatable :: outdiv, indiv
        integer :: i, newnodes, s1, s2
        ! We only need to find the nodes that share an edge.
        allocate(outdiv(c_size), indiv(c_size))
        call select_underedges(adjptr, adjncy, todivide, c_size, c_rank, outdiv)
        call share_division_underedges2(outdiv, indiv)
        do i = 1, c_size
                if ( allocated(outdiv(i)%underedges) ) then
                        if ( (i-1) == c_rank ) then
                                allocate(indiv(i)%underedges(2,&
                                         size(outdiv(i)%underedges,2)))
                                indiv(i)%underedges = outdiv(i)%underedges
                        end if               
                        deallocate(outdiv(i)%underedges)
                end if
        end do
        call create_centernodes(todivide, centerpts, correlems)
        call create_singlenodes(adjptr, corresp, todivide,&
                                singlepts, corredges)
        call create_edgenodes(adjptr, adjncy, corresp, todivide, c_size,&
                              c_rank, outdiv, indiv, nodes)
        do i = 1, c_size
                if ( allocated(indiv(i)%underedges) ) then
                        deallocate(indiv(i)%underedges)
                end if
        end do
        if ( allocated(nodes) ) then
                s1 = size(nodes,2)
        else
                s1 = 0
        end if
        if ( allocated(centerpts) ) then
                s2 = s1 + size(centerpts,2)
                call concat(nodes, centerpts)
        else
                s2 = s1
        end if
        if ( allocated(singlepts) ) then
                newnodes = s2 + size(singlepts,2)
                call concat(nodes, singlepts)
        else
                newnodes = s2
        end if
        allocate(gnodes(newnodes))
        gnodes = (/ ( i, i=1, newnodes ) /)
        call prefix_sum_newnodes(newnodes)
        gnodes = gnodes + newnodes
        do i = 1, c_size
                if ( allocated(outdiv(i)%overgnodes) ) then
                        outdiv(i)%overgnodes = outdiv(i)%overgnodes + newnodes
                end if
        end do
        call share_division_overgnodes(outdiv, indiv)
        call share_division_overedges2(outdiv, indiv)
        call share_division_overpoints2(outdiv, indiv)
        call select_newhanging(todivide, c_size, c_rank, &
                               outdiv, indiv, ishanging)
        call final_division(adjptr, adjncy, corresp, nodes, gnodes,&
                            correlems, corredges, c_size, c_rank, s1, s2,&
                            todivide, ishanging, outdiv, indiv)
        deallocate(ishanging)
        deallocate(indiv, outdiv)
end subroutine

subroutine global_divide(further)
        logical*1, intent(out) :: further
        logical*1, allocatable, dimension(:) :: nfurther
        integer :: ierr
        logical*1, dimension(:), allocatable :: todivide
        integer, dimension(:), allocatable :: adjptr, adjncy
        integer, dimension(:), allocatable :: corresp
        allocate(todivide(m_nelems))
        call find_dual(2,adjptr, adjncy)
        call find_dual_corresp(adjptr, adjncy, corresp)
        call select_global_divide(adjptr, adjncy, todivide)
        
        ! Checking if we need to divide further.
        further = any(todivide)
        if ( c_rank == 0 ) then
                allocate(nfurther(c_size))
                nfurther = .FALSE.
        else
                allocate(nfurther(1)) ! so there is no shitty warning.
        end if
        call mpi_gather(further, 1, MPI_LOGICAL1, nfurther,&
                        1, MPI_LOGICAL1, 0, MPI_COMM_WORLD, ierr)
        call check_mpi(ierr, __LINE__)
        if ( c_rank == 0 ) then
                further = any(nfurther)
        end if
        deallocate(nfurther)
        call mpi_bcast(further,1,MPI_LOGICAL1,0,MPI_COMM_WORLD,ierr)
        call check_mpi(ierr, __LINE__)
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call check_mpi(ierr, __LINE__)
        if ( .NOT. further ) then
                return
        end if
        m_mlength = m_mlength / 2.0D0
        call create_global_divide(adjptr, adjncy, corresp, todivide)
        call update_elmdist()
        deallocate(todivide, adjptr, adjncy, corresp)
end subroutine

subroutine global_cut_geom()
        integer :: ierr
        call cut_geom()
        call update_elmdist()
        call assign_minimum()
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call check_mpi(ierr, __LINE__)
end subroutine

subroutine assign_minimum()
        integer, dimension(:), allocatable   :: part
        integer, dimension(:,:), allocatable :: nelems
        type(elemcsr), dimension(:), allocatable :: outelems, inelems
        integer, dimension(2) :: locelems
        integer :: i, ierr, undermin, currmin
        allocate(nelems(2, c_size))
        locelems(1) = m_nelems 
        locelems(2) = c_rank
        call mpi_gather(locelems, 2, MPI_INTEGER, nelems, 2, &
                        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call check_mpi(ierr, __LINE__)
        if ( c_rank == 0 ) then
                call sort(nelems)
        end if
        call mpi_bcast(nelems, 2*c_size, MPI_INTEGER, &
                       0, MPI_COMM_WORLD, ierr)
        call check_mpi(ierr, __LINE__)
        undermin = 0
        do i = 1, c_size
                if ( nelems(1, i) < MIN_ELEMS ) then
                        undermin = undermin + 1
                end if
        end do
        if ( undermin == 0 ) then
        ! Moving elements is not necessary.
                deallocate(nelems)
                return
        end if
        
        ! Moving elements is necessary.
        allocate(part(m_nelems))
        allocate(outelems(c_size), inelems(c_size))
        currmin = 0
        do i = c_size, 1, -1
        currmin = currmin + 1
        if ( nelems(2,i) == c_rank ) then
                if ( currmin <= undermin ) then
                        call singlepartmesh(part)
                        where ( part == 1 )
                                part = c_rank
                        elsewhere
                                part = nelems(2,currmin)
                        endwhere
                else
                        part = c_rank
                end if
        end if
        end do
        deallocate(nelems)
        call reduce_elems(part, c_rank, c_size, outelems)
        if ( size(m_ediv) /= (count(.not. m_ediv) + count(m_ediv)) ) then
                write(*,*) "Error: ", __FILE__, ":", __LINE__
        end if
        call share_elems(outelems, inelems)
        
        call add_elems(c_size, inelems)
        if ( size(m_ediv) /= (count(.not. m_ediv) + count(m_ediv)) ) then
                write(*,*) "Error: ", __FILE__, ":", __LINE__
        end if
        do i = 1,c_size
                if ( allocated(inelems(i)%idx) ) then
                        deallocate(inelems(i)%idx, inelems(i)%vals, inelems(i)%logi, inelems(i)%prop)
                end if
                if ( allocated(inelems(i)%gnode) ) then
                        deallocate(inelems(i)%gnode, inelems(i)%nodes)
                end if
                if ( allocated(outelems(i)%idx) ) then
                        deallocate(outelems(i)%idx, outelems(i)%vals, outelems(i)%logi, outelems(i)%prop)
                end if
                if ( allocated(outelems(i)%gnode) ) then
                        deallocate(outelems(i)%gnode, outelems(i)%nodes)
                end if
        end do
        call update_elmdist()
        deallocate(part, outelems, inelems)
end subroutine



! Redistributes the elements according to the results given by parmetis mesh
! partition algorithm.

subroutine redistribute()
        integer, dimension(:), allocatable :: part, eind
        type(elemcsr), dimension(:), allocatable :: outelems, inelems
        integer :: i
        allocate(eind(size(m_eind,1)))
        allocate(part(m_nelems))
        allocate(outelems(c_size), inelems(c_size))
        eind(:) = m_gnode(m_eind)
        call partmesh(part, eind)
        call reduce_elems(part, c_rank, c_size, outelems)
        if ( size(m_ediv) /= (count(.not. m_ediv) + count(m_ediv)) ) then
                write(*,*) "Error: ", __FILE__, ":", __LINE__
        end if
        call share_elems(outelems, inelems)
        
        call add_elems(c_size, inelems)
        if ( size(m_ediv) /= (count(.not. m_ediv) + count(m_ediv)) ) then
                write(*,*) "Error: ", __FILE__, ":", __LINE__
        end if
        do i = 1,c_size
                if ( allocated(inelems(i)%idx) ) then
                        deallocate(inelems(i)%idx, inelems(i)%vals, inelems(i)%logi, inelems(i)%prop)
                end if
                if ( allocated(inelems(i)%gnode) ) then
                        deallocate(inelems(i)%gnode, inelems(i)%nodes)
                end if
                if ( allocated(outelems(i)%idx) ) then
                        deallocate(outelems(i)%idx, outelems(i)%vals, outelems(i)%logi, outelems(i)%prop)
                end if
                if ( allocated(outelems(i)%gnode) ) then
                        deallocate(outelems(i)%gnode, outelems(i)%nodes)
                end if
        end do
        call update_elmdist()
        deallocate(eind, part, outelems, inelems)
end subroutine

subroutine init_exfem_comm()
        call init_shared_nodes()
        call init_shared_hanging()
        call init_shared_anchors()
end subroutine

subroutine init_shared_anchors()
        integer :: i, j, k, n
        type(shnodes), dimension(:), allocatable   :: outsh, insh
        type(shnodesll), dimension(:), allocatable :: outll
        integer, dimension(:), allocatable :: ghang, indic
        allocate(outll(c_size))
        allocate(outsh(c_size), insh(c_size))
        if ( allocated(m_hanging) ) then
                allocate(ghang(size(m_hanging,2)*2))
                ghang(1:size(m_hanging,2))  = m_gnode(m_hanging(2,:))
                ghang(size(m_hanging,2)+1:) = m_gnode(m_hanging(3,:))
                call get_unique(ghang)
                allocate(indic(size(ghang)))
                do i = 1, c_size
                if ( allocated(c_neighbors(i)%gnodes) ) then
                        call binsearch(c_neighbors(i)%gnodes, ghang, indic)
                        do j = 1, size(indic)
                        if ( indic(j) /= 0 ) then
                                call insertll(outll(i)%annodes, ghang(j))
                        end if
                        end do 
                        call toarray(outll(i)%annodes, outsh(i)%annodes)
                end if
                end do
        end if
        call share_shnodes_annodes(outsh, insh)
        do i = 1, c_size
        if ( allocated(outsh(i)%annodes) ) then
                call reverse_indices(m_gnode, outsh(i)%annodes)
                ! Locally, we have to keep track of the hanging nodes, so their
                ! value can be passed to the other nodes.
                ! TODO: This is probably VERY SLOW. but it only has to be done
                ! once.
                do j = 1, size(outsh(i)%annodes)
                do k = 1, size(m_hanging)
                if ( m_hanging(2,k) == outsh(i)%annodes(j) .OR. &
                     m_hanging(3,k) == outsh(i)%annodes(j) ) then
                        outsh(i)%annodes(j) = m_hanging(1,k)
                        exit
                end if
                end do
                end do
                n = 2*size(outsh(i)%annodes)
                allocate(c_neighbors(i)%outadofs(n))
                c_neighbors(i)%outadofs(2:n:2) = 2*outsh(i)%annodes
                c_neighbors(i)%outadofs(1:n-1:2) = &
                               c_neighbors(i)%outadofs(2:n:2) - 1
                allocate(c_neighbors(i)%outavalues(n))
        end if
        if ( allocated(insh(i)%annodes) ) then
                call reverse_indices(m_gnode, insh(i)%annodes)
                n = 2*size(insh(i)%annodes)
                allocate(c_neighbors(i)%inadofs(n))
                c_neighbors(i)%inadofs(2:n:2) = 2*insh(i)%annodes
                c_neighbors(i)%inadofs(1:n-1:2) = c_neighbors(i)%inadofs(2:n:2) - 1
                allocate(c_neighbors(i)%inavalues(n))
        end if
        end do
        if ( allocated(m_hanging) ) then
                deallocate(ghang, indic)
        end if
        deallocate(outsh, insh, outll)
end subroutine

subroutine init_shared_hanging()
        integer :: i, j, n
        type(shnodes), dimension(:), allocatable   :: outsh, insh
        type(shnodesll), dimension(:), allocatable :: outll
        integer, dimension(:), allocatable :: ghang, indic
        allocate(outll(c_size))
        allocate(outsh(c_size), insh(c_size))
        if ( allocated(m_hanging) ) then
                allocate(ghang(size(m_hanging,2)), indic(size(m_hanging,2)))
                ghang = m_gnode(m_hanging(1,:))
                do i = 1, c_size
                if ( allocated(c_neighbors(i)%gnodes) ) then
                        call binsearch(c_neighbors(i)%gnodes, ghang, indic)
                        do j = 1, size(indic)
                        if ( indic(j) /= 0 ) then
                                call insertll(outll(i)%hgnodes, ghang(j))
                        end if
                        end do 
                        call toarray(outll(i)%hgnodes, outsh(i)%hgnodes)
                end if
                end do
        end if
        call share_shnodes_hgnodes(outsh, insh)
        do i = 1, c_size
        if ( allocated(outsh(i)%hgnodes) ) then
                call reverse_indices(m_gnode, outsh(i)%hgnodes)
                n = 2*size(outsh(i)%hgnodes)
                allocate(c_neighbors(i)%outhdofs(n))
                c_neighbors(i)%outhdofs(2:n:2) = 2*outsh(i)%hgnodes
                c_neighbors(i)%outhdofs(1:n-1:2) = &
                               c_neighbors(i)%outhdofs(2:n:2) - 1
                allocate(c_neighbors(i)%outhvalues(n))
        end if
        if ( allocated(insh(i)%hgnodes) ) then
                call reverse_indices(m_gnode, insh(i)%hgnodes)
                n = 2*size(insh(i)%hgnodes)
                allocate(c_neighbors(i)%inhdofs(n))
                c_neighbors(i)%inhdofs(2:n:2) = 2*insh(i)%hgnodes
                c_neighbors(i)%inhdofs(1:n-1:2) = &
                               c_neighbors(i)%inhdofs(2:n:2) - 1
                allocate(c_neighbors(i)%inhvalues(n))
        end if
        end do
        if ( allocated(m_hanging) ) then
                deallocate(ghang, indic)
        end if
        deallocate(outsh, insh, outll)
end subroutine

subroutine init_shared_nodes()
        integer, dimension(:), allocatable :: adjptr, adjncy
        integer :: i, j, k, local, proc, ia, ea, ib, eb, global, n
        integer, dimension(8) :: nodes
        type(elemedgell), dimension(:), allocatable :: outll, inll
        type(elemedge), dimension(:), allocatable :: outedges, inedges
        type(neighborll), dimension(:), allocatable :: llneighs
        type(linked_list_1d), pointer :: llneighlist
        nullify(llneighlist)
        allocate(llneighs(c_size))
        allocate(outll(c_size), inll(c_size))
        allocate(outedges(c_size), inedges(c_size))

        ! GLOBAL allocation of information to make this step faster in the
        ! timesteps.

        allocate(c_neighbors(c_size))
        

        call find_dual(1, adjptr, adjncy)
        do i = 1, m_nelems 
                do j = adjptr(i), adjptr(i+1)-1
                ia = m_eptr(i)
                ea = m_eptr(i+1)-1
                call local_element(adjncy(j), local, proc)
                if ( proc /= c_rank ) then
                        call insertll(outll(proc+1)%elem, m_eind(ia:ea))
                        call global_element(i, c_rank, global)
                        call insertll(outll(proc+1)%link, (/global, adjncy(j)/))
                end if
                end do
        end do
        deallocate(adjncy, adjptr)
        do i = 1, c_size
                call toarray_csr(outll(i)%elem, outedges(i)%idx,&
                            outedges(i)%vals)
                call toarray(outll(i)%link, outedges(i)%links)
                if ( allocated(outedges(i)%vals) ) then
                        outedges(i)%vals = m_gnode(outedges(i)%vals)
                end if
        end do
        call share_elemedge_idx(outedges, inedges)
        call share_elemedge_vals(outedges, inedges)
        call share_elemedge_links2(outedges, inedges)
        do i = 1, c_size
                if ( .NOT. allocated(inedges(i)%idx) ) then
                        cycle
                end if
                do j = 1, size(inedges(i)%idx)-1
                        ib = inedges(i)%idx(j)
                        eb = inedges(i)%idx(j+1)-1
                        call local_element(inedges(i)%links(2,j), local, proc)
                        if ( proc /= c_rank ) then
                                write(*,*) "Error: proc ", c_rank,      &
                                           " received elem: ",          &
                                           inedges(i)%links(2,j), proc
                        end if
                        ia = m_eptr(local)
                        ea = m_eptr(local+1)-1
                        call nodes_shared(m_gnode(m_eind(ia:ea)),&
                                          inedges(i)%vals(ib:eb),&
                                          nodes)
                        do k = 1, 8
                                if ( nodes(k) /= 0 ) then
                                        call insertll(llneighs(i)%gnodes,&
                                                      nodes(k))
                                end if
                        end do
                end do
        end do
        do i = 1, c_size
                call toarray(llneighs(i)%gnodes, c_neighbors(i)%gnodes)
                if ( allocated(c_neighbors(i)%gnodes) ) then
                        call insertll(llneighlist, (i-1))
                        call get_unique(c_neighbors(i)%gnodes)
                        call reverse_indices(m_gnode, c_neighbors(i)%gnodes)
                        n = size(c_neighbors(i)%gnodes,1)
                        allocate(c_neighbors(i)%locdofs(2*n))
                        c_neighbors(i)%locdofs(2:2*n:2)   = 2*c_neighbors(i)%gnodes 
                        c_neighbors(i)%locdofs(1:2*n-1:2) = c_neighbors(i)%locdofs(2:2*n:2)-1
                        c_neighbors(i)%gnodes = m_gnode(c_neighbors(i)%gnodes)
                        allocate(c_neighbors(i)%invalues(2*n))
                        allocate(c_neighbors(i)%outvalues(2*n))
                end if
        end do
        call toarray(llneighlist, c_neighlist)
        do i = 1, c_size
                if ( allocated(outedges(i)%idx) ) then
                        deallocate(outedges(i)%idx)
                        deallocate(outedges(i)%vals)
                        deallocate(outedges(i)%links)
                end if
                if ( allocated(inedges(i)%idx) ) then
                        deallocate(inedges(i)%idx)
                        deallocate(inedges(i)%vals)
                        deallocate(inedges(i)%links)
                end if
        end do
        deallocate(inll, outll, outedges, inedges, llneighs)
end subroutine

subroutine calc_global_totnodes()
        integer, dimension(:), allocatable :: numnodes, undergnodes, cpgnodes
        integer :: i, ierr, locnodes
        if ( c_rank == 0 ) then
                allocate(numnodes(c_size))
        end if
        do i = 0, c_rank-1
        if ( allocated(c_neighbors(i+1)%gnodes) ) then
                if ( allocated(undergnodes) ) then
                        allocate(cpgnodes(size(c_neighbors(i+1)%gnodes)))
                        cpgnodes = c_neighbors(i+1)%gnodes
                        call concat(undergnodes, cpgnodes)
                else
                        allocate(undergnodes(size(c_neighbors(i+1)%gnodes)))
                        undergnodes = c_neighbors(i+1)%gnodes
                end if
        end if
        end do
        if ( allocated(undergnodes ) ) then
                call get_unique(undergnodes)
                locnodes = m_nnodes - size(undergnodes)
                deallocate(undergnodes)
        else
                locnodes = m_nnodes
        end if
        call mpi_gather(locnodes, 1, MPI_INTEGER, numnodes, 1,&
                        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if ( c_rank == 0 ) then
                m_totnodes = 0
                do i = 0, c_size
                        m_totnodes = m_totnodes + numnodes(i) 
                end do
                deallocate(numnodes)
        end if
        call mpi_bcast(m_totnodes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call check_mpi(ierr, __LINE__)
end subroutine

subroutine set_neighbor_displacements(disps)
        real*8,  dimension(:), intent(inout) :: disps
        integer, parameter :: TAG_DATA  = 1
        integer :: i, j, ierr
        logical :: flag
        if (.not. allocated(c_neighlist) ) then
        return
        end if
        do j = 1,size(c_neighlist)
        i = c_neighlist(j)       
        if ( allocated(c_neighbors(i+1)%outhdofs) ) then
                c_neighbors(i+1)%outhvalues = disps(c_neighbors(i+1)%outhdofs)
        end if
        end do
        do j = 1,size(c_neighlist)
        i = c_neighlist(j)
        if ( allocated(c_neighbors(i+1)%inhdofs) ) then
                call mpi_irecv(c_neighbors(i+1)%inhvalues,      &
                               size(c_neighbors(i+1)%inhvalues),& 
                               MPI_REAL8,i,TAG_DATA,            &
                               MPI_COMM_WORLD,                  &
                               c_neighbors(i+1)%reqrcv,ierr)
                call check_mpi(ierr,__LINE__)
        end if
        if ( allocated(c_neighbors(i+1)%outhdofs) ) then
                call mpi_isend(c_neighbors(i+1)%outhvalues,      &
                               size(c_neighbors(i+1)%outhvalues),&
                               MPI_REAL8,i,TAG_DATA,             &
                               MPI_COMM_WORLD,                   &
                               c_neighbors(i+1)%reqsnd,ierr)
                call check_mpi(ierr,__LINE__)
        end if
        end do
        do j = 1,size(c_neighlist)
        i = c_neighlist(j)
        if ( allocated(c_neighbors(i+1)%inhdofs) ) then
                call mpi_wait(c_neighbors(i+1)%reqrcv,c_neighbors(i+1)%statrcv,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_test_cancelled(c_neighbors(i+1)%statrcv,&
                                        flag,ierr)
                call check_mpi(ierr,__LINE__)
                if ( flag ) then
                        write(*,*) "Cancelled communication receive: ", c_rank, i
                end if
        end if
        if ( allocated(c_neighbors(i+1)%outhdofs) ) then
                call mpi_wait(c_neighbors(i+1)%reqsnd,c_neighbors(i+1)%statsnd,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_test_cancelled(c_neighbors(i+1)%statsnd,&
                                        flag,ierr)
                call check_mpi(ierr,__LINE__)
                if ( flag ) then
                        write(*,*) "Cancelled communication send: ", c_rank, i
                end if
        end if
        end do
        do j = 1,size(c_neighlist)
        i = c_neighlist(j)
        if ( allocated(c_neighbors(i+1)%inhdofs) ) then
                disps(c_neighbors(i+1)%inhdofs) = c_neighbors(i+1)%inhvalues
        end if
        end do
end subroutine

subroutine sum_neighbor_anchors(reactions)
        real*8,  dimension(:), intent(inout) :: reactions
        integer, parameter :: TAG_DATA  = 1
        integer :: i, j, ierr
        logical :: flag
        if (.not. allocated(c_neighlist) ) then
        return
        end if
        do j = 1,size(c_neighlist)
        i = c_neighlist(j)       
        if ( allocated(c_neighbors(i+1)%outadofs) ) then
                c_neighbors(i+1)%outavalues = reactions(c_neighbors(i+1)%outadofs)/2.0D0
        end if
        end do
        do j = 1,size(c_neighlist)
        i = c_neighlist(j)
        if ( allocated(c_neighbors(i+1)%inadofs) ) then
                call mpi_irecv(c_neighbors(i+1)%inavalues,      &
                               size(c_neighbors(i+1)%inavalues),& 
                               MPI_REAL8,i,TAG_DATA,            &
                               MPI_COMM_WORLD,                  &
                               c_neighbors(i+1)%reqrcv,ierr)
                call check_mpi(ierr,__LINE__)
        end if
        if ( allocated(c_neighbors(i+1)%outadofs) ) then
                call mpi_isend(c_neighbors(i+1)%outavalues,      &
                               size(c_neighbors(i+1)%outavalues),&
                               MPI_REAL8,i,TAG_DATA,             &
                               MPI_COMM_WORLD,                   &
                               c_neighbors(i+1)%reqsnd,ierr)
                call check_mpi(ierr,__LINE__)
        end if
        end do
        do j = 1,size(c_neighlist)
        i = c_neighlist(j)
        if ( allocated(c_neighbors(i+1)%inadofs) ) then
                call mpi_wait(c_neighbors(i+1)%reqrcv,c_neighbors(i+1)%statrcv,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_test_cancelled(c_neighbors(i+1)%statrcv,&
                                        flag,ierr)
                call check_mpi(ierr,__LINE__)
                if ( flag ) then
                        write(*,*) "Cancelled communication receive: ", c_rank, i
                end if
        end if
        if ( allocated(c_neighbors(i+1)%outadofs) ) then
                call mpi_wait(c_neighbors(i+1)%reqsnd,c_neighbors(i+1)%statsnd,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_test_cancelled(c_neighbors(i+1)%statsnd,&
                                        flag,ierr)
                call check_mpi(ierr,__LINE__)
                if ( flag ) then
                        write(*,*) "Cancelled communication send: ", c_rank, i
                end if
        end if
        end do
        do j = 1,size(c_neighlist)
        i = c_neighlist(j)
        if ( allocated(c_neighbors(i+1)%inadofs) ) then
                reactions(c_neighbors(i+1)%inadofs) = reactions(c_neighbors(i+1)%inadofs) + c_neighbors(i+1)%inavalues
        end if
        end do
end subroutine

subroutine sum_neighbor_reactions(reactions)
        real*8,  dimension(:), intent(inout) :: reactions
        integer, parameter :: TAG_DATA  = 1
        integer :: i, j, ierr
        logical :: flag
        if (.not. allocated(c_neighlist) ) then
                return
        end if
        do j = 1,size(c_neighlist)
        i = c_neighlist(j)
        if ( allocated(c_neighbors(i+1)%gnodes) ) then
                c_neighbors(i+1)%outvalues = reactions(c_neighbors(i+1)%locdofs)
        end if
        end do
        do j = 1,size(c_neighlist)
        i = c_neighlist(j)       
        if ( allocated(c_neighbors(i+1)%gnodes) ) then
                call mpi_irecv(c_neighbors(i+1)%invalues,      &
                               size(c_neighbors(i+1)%invalues),& 
                               MPI_REAL8,i,TAG_DATA,        &
                               MPI_COMM_WORLD,              &
                               c_neighbors(i+1)%reqrcv,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_isend(c_neighbors(i+1)%outvalues,      &
                               size(c_neighbors(i+1)%outvalues),&
                               MPI_REAL8,i,TAG_DATA,         &
                               MPI_COMM_WORLD,               &
                               c_neighbors(i+1)%reqsnd,ierr)
                call check_mpi(ierr,__LINE__)
        end if
        end do
        do j = 1,size(c_neighlist)
        i = c_neighlist(j)
        if ( allocated(c_neighbors(i+1)%gnodes) ) then
                call mpi_wait(c_neighbors(i+1)%reqsnd,c_neighbors(i+1)%statsnd,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_test_cancelled(c_neighbors(i+1)%statsnd,&
                                        flag,ierr)
                call check_mpi(ierr,__LINE__)
                if ( flag ) then
                        write(*,*) "Cancelled communication send: ", c_rank, i
                end if
                call mpi_wait(c_neighbors(i+1)%reqrcv,c_neighbors(i+1)%statrcv,ierr)
                call check_mpi(ierr,__LINE__)
                call mpi_test_cancelled(c_neighbors(i+1)%statrcv,&
                                        flag,ierr)
                call check_mpi(ierr,__LINE__)
                if ( flag ) then
                        write(*,*) "Cancelled communication receive: ", c_rank, i
                end if
        end if
        end do
        do j = 1,size(c_neighlist)
        i = c_neighlist(j)
        if ( allocated(c_neighbors(i+1)%gnodes) ) then
                reactions(c_neighbors(i+1)%locdofs) = reactions(c_neighbors(i+1)%locdofs) +&
                                                      c_neighbors(i+1)%invalues
        end if
        end do
end subroutine

subroutine bcast_visual_params(vtype, params)
        character(len=6), intent(inout) :: vtype
        real*8, dimension(5), intent(inout) :: params
        integer :: ierr
        
        call mpi_bcast(vtype,6,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
        call check_mpi(ierr,__LINE__)
        call mpi_bcast(params,5,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call check_mpi(ierr,__LINE__)
        
        ! Wait until all messages are sent.
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call check_mpi(ierr,__LINE__)


end subroutine

subroutine share_sheet(npoints, lhigh, ghigh, locfield, allfield)
        integer, intent(in) :: npoints
        integer, dimension(:), allocatable, intent(in) :: lhigh
        integer, dimension(:), intent(in)   :: ghigh
        real*8, dimension(:), intent(in)    :: locfield
        real*8, dimension(:), intent(inout) :: allfield
        integer :: i, cnt, ierr
        integer, dimension(MPI_STATUS_SIZE) :: st
        integer, parameter :: DATA_TAG = 1
        cnt = 1
        do i = 1, npoints
        if ( allocated(lhigh) .AND. size(lhigh) >= cnt ) then
        if ( lhigh(cnt) == i ) then
                if ( c_rank /= 0 ) then

                        call mpi_send(locfield(cnt*2-1:cnt*2), 2, MPI_REAL8, 0,&
                                      DATA_TAG, MPI_COMM_WORLD, ierr)
                        call check_mpi(ierr,__LINE__)
                else        
                        allfield(i*2-1:i*2) = locfield(cnt*2-1:cnt*2)
                end if
                cnt = cnt + 1
        end if
        end if
        if ( c_rank == 0 ) then
                if ( ghigh(i) > 0 ) then
                        call mpi_recv(allfield(i*2-1:i*2), 2, MPI_REAL8,ghigh(i),&
                                      DATA_TAG, MPI_COMM_WORLD, st, ierr)
                        call check_mpi(ierr,__LINE__)
                else if ( ghigh(i) < 0 ) then
                        allfield(i*2-1:i*2) = 0
                end if
        end if
        end do
        call mpi_barrier(MPI_COMM_WORLD, ierr)
        call check_mpi(ierr, __LINE__)
end subroutine

subroutine share_highest(points, elems, lhigh, ghigh)
        real*8, dimension(:,:)  :: points
        integer, dimension(:)   :: elems
        integer, dimension(:), allocatable :: lhigh
        integer, dimension(:)   :: ghigh
        real*8, dimension(:), allocatable :: allvals
        real*8 :: locval, minvalu
        integer :: i, j, ierr, chosen
        type(linked_list_1d), pointer :: lllhigh
        
        nullify(lllhigh)
        if ( c_rank == 0 ) then
                allocate(allvals(c_size))
        end if
        do i = 1, size(points, 2)
                if ( elems(i) /= 0 ) then
                        locval = points(2,i)
                else
                        locval = u_inf
                end if
                call mpi_gather(locval, 1, MPI_REAL8, allvals, 1, &
                                MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
                call check_mpi(ierr,__LINE__)
                if ( c_rank == 0 ) then
                        minvalu = u_inf
                        ghigh(i) = -1 
                        do j = 0, (c_size-1)
                        if ( realunder(allvals(j+1), minvalu ) ) then
                                minvalu = allvals(j+1)
                                ghigh(i) = j
                        end if
                        end do
                        chosen = ghigh(i)
                        if ( (minvalu - minvalu) .ne. 0 ) then
                                ghigh(i) = -1 
                        else
                                points(2,i) = minvalu
                        end if
                end if
                call mpi_bcast(chosen,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                call check_mpi(ierr,__LINE__)
                if ( chosen == c_rank ) then
                        call insertll(lllhigh, i)
                end if
        end do
        call toarray(lllhigh, lhigh)
end subroutine

subroutine share_visu_nodes(nodes, nredu, owner)
        integer,   dimension(:), intent(in), allocatable :: nodes
        logical*1, dimension(:), intent(inout), allocatable :: nredu
        integer,   dimension(:), intent(out) :: owner
        integer,   dimension(:), allocatable :: inn, innodes, outnodes, indices
        integer :: nnodes, outn, i, ierr
        integer, dimension(MPI_STATUS_SIZE) :: st
        integer, parameter :: DATA_TAG = 1 
        ! The code for 0 so the rest can be for the other nodes.
        nnodes = 0
        if ( c_rank == 0 ) then
                nnodes = size(owner)

                allocate( inn(c_size) )
                owner = -1
                if ( allocated(nodes) ) then
                        owner(nodes) = 0
                end if
        end if
        if ( allocated(nodes) ) then
                outn = size(nodes)
        else
                outn = 0
        end if
        call mpi_gather(outn, 1, MPI_INTEGER, inn, 1,      &
                        MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call check_mpi(ierr,__LINE__)
        
        do i = 1, (c_size-1) ! from 0 it would never stop. Well, may be.
        if ( c_rank == i .AND. outn > 0 ) then
                call mpi_send(nodes, outn, MPI_INTEGER, 0,      &
                              DATA_TAG, MPI_COMM_WORLD, ierr)
                call check_mpi(ierr,__LINE__)
        end if
        
        if ( c_rank == 0 .AND. inn(i+1) > 0 ) then
                allocate(innodes(inn(i+1)))
                call mpi_recv(innodes, inn(i+1), MPI_INTEGER, i, &
                              DATA_TAG, MPI_COMM_WORLD, st, ierr )
                call check_mpi(ierr,__LINE__)
                owner(innodes) = i
                deallocate(innodes)
        end if
        end do
        if ( c_rank == 0 ) then
                do i = 0, (c_size-1)
                        inn(i+1) = count( owner == i )
                end do
        end if
        call mpi_scatter(inn, 1, MPI_INTEGER, outn, 1,       &
                         MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call check_mpi(ierr,__LINE__)
        do i = 1, (c_size-1) ! from 0 it would never stop.
        if ( c_rank == i ) then
        if ( outn > 0 ) then
                allocate(innodes(outn)) 
                allocate(indices(size(nodes)))
                call mpi_recv(innodes, outn, MPI_INTEGER, 0,    &
                              DATA_TAG, MPI_COMM_WORLD, st, ierr)
                call check_mpi(ierr,__LINE__)
                call binsearch(innodes, nodes, indices)
                nredu = (indices == 0)
                deallocate(innodes, indices)
        else
                if ( allocated(nredu) ) then
                        nredu = .FALSE.
                end if
        end if
        end if
        if ( c_rank == 0 .AND. inn(i+1) > 0 ) then
                allocate(outnodes(inn(i+1)))
                outnodes = pack((/(i, i = 1, nnodes)/), owner == i)
                call mpi_send(outnodes, inn(i+1), MPI_INTEGER,     &
                              i, DATA_TAG, MPI_COMM_WORLD, ierr)
                call check_mpi(ierr,__LINE__)
                deallocate(outnodes)
        end if
        end do
        if ( c_rank == 0 ) then
                deallocate(inn)
                if ( allocated(nodes) ) then
                        nredu = (owner(nodes) /= 0)
                end if
        end if
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call check_mpi(ierr,__LINE__)
end subroutine

subroutine share_field(field, fielddofs, nown, allfield)
        real*8, dimension(:), intent(in), allocatable :: field
        integer, dimension(:), intent(in) :: fielddofs, nown
        real*8, dimension(:), intent(out) :: allfield
        integer :: i, n, ierr, arrb, arre
        integer, dimension(MPI_STATUS_SIZE) :: st
        integer, parameter :: DATA_TAG = 1
        arre = 0
        if ( c_rank == 0 .AND. allocated(field) ) then
                arrb = 1
                arre = size(field) 
                allfield(arrb:arre) = field
        end if
        
        do i = 1, (c_size-1)
        if ( c_rank == i ) then
                if ( allocated(field) ) then
                        n = size(field)
                else
                        n = 0
                end if
                if ( n > 0 ) then
                call mpi_send(field, n, MPI_REAL8, 0,         &
                              DATA_TAG, MPI_COMM_WORLD, ierr)
                call check_mpi(ierr,__LINE__)
                end if
        end if
        if ( c_rank == 0 ) then
                n = nown(i+1) 
                if ( n > 0 ) then
                        arrb = arre + 1
                        arre = arre + 2*n
                        call mpi_recv(allfield(arrb:arre), 2*n, MPI_REAL8,i,&
                                      DATA_TAG, MPI_COMM_WORLD, st, ierr)
                        call check_mpi(ierr,__LINE__)
                end if
        end if
        end do
        
        if ( c_rank == 0 ) then
                if ( arre < size(allfield) ) then
                        arrb = arre + 1
                        allfield(arrb:size(allfield)) = 0.0D0 
                end if
                allfield(fielddofs) = allfield
        end if
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call check_mpi(ierr,__LINE__)
end subroutine

subroutine save_io_ranks()
        character(len=256) :: filename
        integer :: fnlen, i, ierr
        fnlen = io_fnodlen + 10
        write(filename,"(A,A10)") io_outdir(1:io_fnodlen),"/ranks.bin"
        do i = 0, (c_size-1)
        if (c_rank == i ) then
                call remove_rank_file(fnlen, filename)
        end if
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call check_mpi(ierr,__LINE__)
        end do
        do i = 0, (c_size-1)
        if (c_rank == i ) then
                call save_rank_file(fnlen, filename, c_rank)
        end if
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call check_mpi(ierr,__LINE__)
        end do
end subroutine

subroutine load_io_ranks()
        character(len=256) :: filename
        integer :: fnlen, i, ierr
        fnlen = io_fnodlen + 10
        write(filename,"(A,A10)") io_outdir(1:io_fnodlen),"/ranks.bin"
        do i = 0, (c_size-1)
        if (c_rank == i ) then
                call restart_rank_file(fnlen, filename)
        end if
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call check_mpi(ierr,__LINE__)
        end do
        do i = 0, (c_size-1)
        if (c_rank == i ) then
                call load_rank_file(fnlen, filename, io_locprocs, io_rank)
        end if
        call mpi_barrier(MPI_COMM_WORLD,ierr)
        call check_mpi(ierr,__LINE__)
        end do
end subroutine

subroutine recalculate_deltat()
        real*8 :: mindt, eldt
        integer :: i, ierr
        if ( c_rank == 0 ) then
                write(*,*) "Original deltat: ", m_deltat
        end if
        mindt = m_deltat
        do i = 1, m_nelems
                call recalcstep(m_matconst(4, m_ematprop(i)),&
                                elength(i), eldt)
                mindt = min(mindt, eldt)
        end do
        call mpi_reduce(mindt, m_deltat, 1, MPI_REAL8, &
                        MPI_MIN, 0, MPI_COMM_WORLD, ierr)
        call check_mpi(ierr,__LINE__)
        if ( c_rank == 0 ) then
                write(*,*) "Updated deltat", m_deltat
                m_nsteps = int(m_duration/m_deltat)
        end if
        call mpi_bcast(m_nsteps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call check_mpi(ierr,__LINE__)
        call mpi_bcast(m_deltat,1,MPI_DOUBLE_PRECISION,&
                       0,MPI_COMM_WORLD,ierr)
        call check_mpi(ierr,__LINE__)
        if ( c_rank == 0 ) then
                io_interval = int(m_inttime/m_deltat)
                m_inttime = m_deltat*io_interval
        end if
        call mpi_bcast(m_inttime,1,MPI_DOUBLE_PRECISION,0, &
                       MPI_COMM_WORLD,ierr)
        call check_mpi(ierr,__LINE__)
        call mpi_bcast(io_interval,1,MPI_INTEGER, 0,MPI_COMM_WORLD,ierr)
        call check_mpi(ierr,__LINE__)
end subroutine

end module communicator

