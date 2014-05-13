module visualization
use explicit
use quadtree
implicit none

contains

subroutine init_visualization()
        if ( c_rank == 0 ) then
                call load_problem()
        end if
        call bcast_problem()
        if ( c_rank == 0 ) then
                call load_environ()
        end if
        call bcast_environ()
        if ( m_prob(1:4) == 'WAVE' ) then
                call init_wave()
        end if
        call load_tparams(c_rank, m_deltat, io_interval)
end subroutine

subroutine create_visualization()
        character(len=6)      :: vtype
        real*8, dimension(5)  :: params
        real*8, dimension(2)  :: origin, space
        integer, dimension(2) :: dims
        integer, dimension(:), allocatable  :: elems, owner, &
                                               selems, lhigh, ghigh
        real*8, dimension(:,:), allocatable :: points, spoints
        call load_mesh(io_rank, m_eptr, m_eind, m_ediv, m_nodes)
        m_nelems = size(m_eptr,1)-1
        m_nnodes = size(m_nodes,2)
        call create_quadtree()
        if ( c_rank == 0 ) then
                call visual_params(vtype, params) 
        end if
        call bcast_visual_params(vtype, params)
        
        if ( vtype(1:5) == 'SHEET' ) then 
                ! Visualization of the sheet.
                call find_highest(params(1), params(2), params(3),&
                                  points, elems, lhigh, ghigh)
                call visualize_sheet(points, elems, lhigh, ghigh)
        else if ( vtype(1:6) == 'ANIMAT' ) then
                ! Visualization of the solution.
                call find_local(params(1), params(2), params(3),    &
                                params(4), params(5), dims,         &
                                origin, space, points, elems, owner)
                call visualize_solution(dims, origin, space, points, elems, owner)
        end if
end subroutine

subroutine find_local(xmin, xmax, ymin, ymax, dist, dims,&
                      origin, space, lpoints, elems, owner)
        real*8, intent(inout) :: xmin, xmax, ymin, ymax, dist
        integer, dimension(:), allocatable, intent(inout) :: elems, owner
        real*8, dimension(:,:), allocatable, intent(inout) :: lpoints
        real*8, dimension(2), intent(out)  :: origin, space
        integer, dimension(2), intent(out) :: dims
        real*8, dimension(:,:), allocatable  :: points
        type(linked_list_1d), pointer        :: lllocal, llelems
        integer, dimension(:), allocatable   :: local
        logical*1, dimension(:), allocatable :: nredu
        integer :: i, elem
        call visu_points(xmin, xmax, ymin, ymax, dist, dims, origin, space, points)
        nullify(lllocal)
        nullify(llelems)
        do i = 1, size(points,2)
        call search_quadtree(points(:,i), elem)
        if ( elem /= 0 ) then
                call insertll(lllocal, i)
                call insertll(llelems, elem)
        end if
        end do
        call toarray(lllocal, local)
        call toarray(llelems, elems)
        if ( allocated(local) ) then
                allocate(nredu(size(local,1)))
        end if
        if ( c_rank == 0 ) then
                allocate(owner(size(points,2)))
        end if
        call share_visu_nodes(local, nredu, owner)
        if ( allocated(local) ) then
                call reduce(local, nredu)
                call reduce(elems, nredu)
                allocate(lpoints(2, size(elems,1)))
                lpoints = points(:,local)
                deallocate(local)
                deallocate(nredu)
                
        end if
        call find_error(owner)
        deallocate(points)
end subroutine

subroutine find_error(owner)
        integer, dimension(:), intent(in) :: owner
        integer, dimension(:), allocatable :: nfound
        type(linked_list_1d), pointer :: llnfound
        integer :: n, i, ierr
        if ( c_rank == 0 ) then
                nullify(llnfound)
                do i = 1, size(owner)
                if ( owner(i) == -1 ) then
                        call insertll(llnfound, i)
                end if
                end do
                call toarray(llnfound, nfound)
        end if
        if ( allocated(nfound) ) then
                n = size(nfound)
        else
                n = 0
        end if
        call mpi_bcast(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call check_mpi(ierr, __LINE__)
        if ( c_rank == 0 ) then

                if ( n == 0 ) then
                        write(*,*) "All points found in the quadtree."
                else
                        write(*,*) "Found: ", n, "points not in the quadtree."
                
                end if
        end if
end subroutine

subroutine visu_points(xmin, xmax, ymin, ymax, dist, dims, origin, space, points)
                real*8, intent(inout) :: xmin, xmax, ymin, ymax, dist
                real*8, dimension(:,:), intent(inout), allocatable :: points
                real*8, dimension(2), intent(out)  :: origin, space
                integer, dimension(2), intent(out) :: dims
                integer :: nx, ny, npts, i, j, ind
                nx = int((xmax-xmin)/dist)
                ny = int((ymax-ymin)/dist)
                xmax = xmin + real(nx, kind=8)*dist
                ymax = ymin + real(ny, kind=8)*dist
                npts = (nx+1)*(ny+1)
                allocate(points(2,npts))
                ind = 0
                do i = 0, nx
                do j = 0, ny
                        ind = ind + 1
                        points(:,ind) = (/xmin + real(i, kind=8)*dist,&
                                          ymin + real(j, kind=8)*dist/)
                end do
                end do
                dims = (/nx, ny/)
                origin = (/xmin, ymin/)
                space  = (/dist, dist/)
end subroutine


subroutine visu_field(ele, pt, time, field)
        integer, intent(in) :: ele
        real*8, dimension(2), intent(in)  :: pt
        real*8, intent(in) :: time
        real*8, dimension(2), intent(out) :: field
        
        integer, dimension(4)  :: nodes
        integer, dimension(8)  :: ndofs
        real*8, dimension(2,4) :: corners
        real*8, dimension(8)   :: dofs
        real*8, dimension(2)   :: outp
        real*8, dimension(2,3) :: disps
        real*8 :: outt

        call elemallnhnodes(ele, nodes)
        call elemallnhdofs(ele, ndofs)
        corners = m_nodes(:,nodes)
        dofs = e_u(ndofs)
        call calcfieldat(corners, dofs, pt, field)
        if ( inside_soil(pt) ) then
                call timetransform(m_wttoorigin, time, outt)
                call pointtransform(m_worigin,pt,outp)
                call planewave(m_prob,outt,outp(1),outp(2),disps)
                field = field + disps(:,1)
        end if
end subroutine


subroutine visu_material(pt, fld)
        real*8, dimension(2), intent(in)  :: pt
        real*8, dimension(2), intent(out) :: fld
        integer :: prop
        prop = material_prop(pt)
        fld = (/ real(prop, kind=8),  m_matconst(4, prop) /)
end subroutine

subroutine visualize_solution(dims, origin, space, points, elems, owner)
        real*8, dimension(2), intent(in)  :: origin, space
        integer, dimension(2), intent(in) :: dims
        real*8, dimension(:,:), intent(in), allocatable :: points
        integer, dimension(:), intent(in), allocatable :: elems, owner
        real*8, dimension(2)   :: field
        real*8, dimension(:), allocatable  :: allfield, locfield
        integer, dimension(:), allocatable :: fielddofs, nown, eptr, eind
        real*8, dimension(:,:), allocatable  :: pnodes
        integer :: i, j, st, step, nlocal, arrb, arre
        real*8 :: time
        if ( c_rank == 0 ) then
                call create_struct_mesh(dims, origin, space, &
                                        eptr, eind, pnodes)
                allocate(allfield(2*size(owner)))
                allocate(fielddofs(2*size(owner)))
                allocate(nown(c_size))
                fielddofs = 0
                arre = 0
                do i = 0, (c_size-1)
                nown(i+1) = count(owner == i)
                if (nown(i+1) > 0 ) then
                        arrb = arre + 1
                        arre = arre + 2*nown(i+1)
                        fielddofs(arrb:arre-1:2) =   &
                                  pack((/(2*j-1,j=1, &
                                  size(owner))/),    &
                                  owner == i)
                        fielddofs(arrb+1:arre:2) =  &
                                  pack((/(2*j,j=1,  &
                                  size(owner))/),   &
                                  owner == i)
                end if
                if ( any(owner == -1) ) then
                        arrb = arre + 1
                        fielddofs(arrb:2*size(owner)-1:2) =        &
                                  pack((/(2*j-1,j=1,size(owner))/),&
                                  owner == -1)
                        fielddofs(arrb+1:2*size(owner):2) =      &
                                  pack((/(2*j,j=1,size(owner))/),&
                                  owner == -1)
                end if
                end do
        end if
        if ( allocated(elems) ) then
                nlocal = size(elems)
        else
                nlocal = 0
        end if
        if ( nlocal > 0 ) then
                allocate(locfield(2*nlocal))
        end if
        
        if ( m_geom(1:3) == 'CVM' ) then
                call init_cvm()
        end if
        do i = 1, nlocal 
                call visu_material(points(:,i), field)
                locfield((/2*i-1, 2*i/)) = field
        end do
        call share_field(locfield, fielddofs, nown, allfield)
        if ( c_rank == 0 ) then
                call print_vtkun_mater(eptr, eind, pnodes, allfield)
        end if
        if ( m_geom(1:3) == 'CVM' ) then
                call end_cvm()
        end if
        
        st = 0
        step = 0
        call load_step(io_rank, step, e_u, st)
        do while ( st /= -1 )
                time = m_deltat*step*io_interval
                do i = 1, nlocal 
                        call visu_field(elems(i), points(:,i), time, field)
                        locfield((/2*i-1, 2*i/)) = field
                end do
                call share_field(locfield, fielddofs, nown, allfield)
                if ( c_rank == 0 ) then
                        call print_vtkun_field(step+1, eptr, eind, pnodes, allfield)
                end if
                step = step + 1
                call load_step(io_rank, step, e_u, st)
        end do
        if ( nlocal > 0 ) then
                deallocate(locfield)
        end if
        if ( c_rank == 0 ) then
                deallocate(allfield)
                deallocate(nown)
                deallocate(fielddofs)
        end if
end subroutine

subroutine find_highest(xmin, xmax, dist, points, elems, lhigh, ghigh)
        real*8 :: xmin, xmax, dist
        real*8, dimension(:,:), allocatable  :: points
        integer, dimension(:), allocatable   :: elems
        integer, dimension(:), allocatable   :: lhigh, ghigh
        integer :: i, npoints
        integer, dimension(4) :: nodes
        
        npoints = int((xmax-xmin)/dist)+1
        allocate(points(2, npoints))
        allocate(elems(npoints))
        do i = 1, npoints
                points(1, i) = xmin+(dist*real(i,kind=8))
                call highest_quadtree(points(1,i), points(2,i), elems(i))
        end do
        if ( c_rank == 0 ) then
                allocate(ghigh(npoints))
        end if
        call share_highest(points, elems, lhigh, ghigh)
end subroutine

subroutine visualize_sheet(points, elems, lhigh, ghigh)
        real*8, dimension(:,:), allocatable :: points
        integer, dimension(:), allocatable :: elems
        integer, dimension(:), allocatable :: lhigh, ghigh
        real*8, dimension(:), allocatable  :: locfield, allfield
        real*8, dimension(2) :: field
        integer :: i, step, st, nlocal
        real*8  :: time
        if ( c_rank == 0 ) then
                call save_sheet_points(points)
                allocate(allfield(size(points,2)*2))
        end if
        if ( allocated(lhigh) ) then
                allocate(locfield(2*size(lhigh)))
                nlocal = size(lhigh)
        else
                nlocal = 0
        end if
        st = 0
        step = 0
        call load_step(io_rank, step, e_u, st)
        do while ( st /= -1 )
                time = m_deltat*step*io_interval
                do i = 1, nlocal
                        call visu_field(elems(lhigh(i)), points(:,lhigh(i)), time, field)
                        locfield((/2*i-1, 2*i/)) = field
                end do
                call share_sheet(size(points,2), lhigh, ghigh, locfield, allfield)
                if ( c_rank == 0 ) then
                        call save_sheet_step(allfield)
                end if
                step = step + 1
                call load_step(io_rank, step, e_u, st)
        end do
        if ( c_rank == 0 ) then
                call save_sheet_close(m_inttime)
        end if
        if ( allocated(points) ) then
                deallocate(points)
                deallocate(elems)
        end if
        if ( allocated(lhigh) ) then
                deallocate(lhigh)
        end if
        if ( allocated(ghigh) ) then
                deallocate(ghigh)
        end if
end subroutine

end module

