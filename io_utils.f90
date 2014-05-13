module io_utils
use utilities
implicit none

integer :: io_interval
logical*1 :: io_newsave
character(len=256) :: io_outdir, io_cvmdir, io_visdir, io_indir
integer :: io_fnodlen, io_fncdlen, io_fnvdlen, io_fnindlen
integer :: io_rank
integer :: io_locprocs
real  :: io_tstart, io_tend


interface 
        subroutine remove_rank_file(nfname, filename)                &
                                    bind(C, name="remove_rank_file")
                use iso_c_binding, only : c_int, c_char
                integer (c_int)    :: nfname
                character (c_char) :: filename
        end subroutine
end interface

interface 
        subroutine load_rank_file(nfname, filename, nranks, myrank)  &
                                  bind(C, name="load_rank_file")
                use iso_c_binding, only : c_int, c_char
                integer (c_int)    :: nfname
                character (c_char) :: filename
                integer (c_int)    :: nranks
                integer (c_int)    :: myrank
        end subroutine
end interface

interface
        subroutine restart_rank_file(nfname, filename)                 &
                                     bind(C, name="restart_rank_file")
                use iso_c_binding, only : c_int, c_char
                integer (c_int)    :: nfname
                character (c_char) :: filename
        end subroutine
end interface

interface 
        subroutine save_rank_file(nfname, filename, myrank)      &
                                  bind(C, name="save_rank_file")
                use iso_c_binding, only : c_int, c_char
                integer (c_int)    :: nfname
                character (c_char) :: filename
                integer (c_int)    :: myrank
        end subroutine
end interface

interface
        subroutine load_environ_file(nfname, filename, flen, &
                                     fname, clen, cnam, vlen, vnam)   &
                                    bind(C, name="load_environ_file")
                use iso_c_binding, only : c_int, c_char
                integer (c_int)    :: nfname
                character(c_char)  :: filename
                integer (c_int)    :: flen
                character (c_char) :: fname
                integer (c_int)    :: clen
                character (c_char) :: cnam
                integer (c_int)    :: vlen
                character (c_char) :: vnam
        end subroutine
end interface

interface
        subroutine load_problem_file(nfname, filename, height, width, time, &
                                     inttime, chper, nperlen, sperlen,      &
                                      nmats,  cmats, vmats, prob,           &
                                     nppars, vppars, geom, ngpars, vgpars)  &
                                     bind(C, name="load_problem_file")
                use iso_c_binding, only : c_ptr, c_double, c_int, c_char
                integer(c_int)     :: nfname
                character(c_char)  :: filename
                real (c_double)    :: height
                real (c_double)    :: width
                real (c_double)    :: time
                real (c_double)    :: inttime
                real (c_double)    :: chper
                real (c_double)    :: nperlen 
                real (c_double)    :: sperlen
                integer (c_int)    :: nmats
                type(c_ptr)        :: cmats
                type(c_ptr)        :: vmats
                character (c_char) :: prob
                integer (c_int)    :: nppars
                type(c_ptr)        :: vppars
                character (c_char) :: geom
                integer (c_int)    :: ngpars
                type(c_ptr)        :: vgpars
        end subroutine
end interface

interface
        subroutine save_binmesh(nfname, filename,  nelems, eptr, eind,  &
                                ediv, nnodes, nodes) &
                                bind(C, name="save_binmesh")
                use iso_c_binding, only : c_int, c_char, c_double, c_bool
                implicit none
                integer   (c_int)    :: nfname
                character (c_char)   :: filename
                integer   (c_int)    :: nelems
                integer   (c_int)    :: eptr(*)
                integer   (c_int)    :: eind(*)
                logical   (c_bool)   :: ediv(*)
                integer   (c_int)    :: nnodes
                real      (c_double) :: nodes(*)
        end subroutine
end interface

interface
        subroutine load_binmesh(nfname, filename,  nelems, eptr, eind,  &
                                ediv, nnodes, nodes) &
                                bind(C, name="load_binmesh")
                use iso_c_binding, only : c_ptr, c_int, c_char, c_double
                implicit none
                integer   (c_int)    :: nfname
                character (c_char)   :: filename
                integer   (c_int)    :: nelems
                type(c_ptr)          :: eptr
                type(c_ptr)          :: eind
                type(c_ptr)          :: ediv
                integer   (c_int)    :: nnodes
                type(c_ptr)          :: nodes
        end subroutine
end interface

interface
        subroutine save_binmass(nfname, filename, ndofs, mass) &
                             bind(C, name="save_binmass")
                use iso_c_binding, only : c_int, c_char, c_double
                implicit none
                integer (c_int)    :: nfname
                character (c_char) :: filename
                integer (c_int)    :: ndofs
                real (c_double)    :: mass(*)
        end subroutine
end interface

interface
        subroutine load_binmass(nfname, filename, ndofs, mass) &
                             bind(C, name="load_binmass")
                use iso_c_binding, only : c_ptr, c_char, c_int, c_double
                implicit none
                integer (c_int)    :: nfname
                character (c_char) :: filename
                integer (c_int)    :: ndofs
                type (c_ptr)       :: mass
        end subroutine
end interface
interface
        subroutine save_binstep(nfname, filename, mode, step, ndofs, u) &
                                 bind(C, name="save_binstep")
                use iso_c_binding, only : c_int, c_double, c_char
                implicit none
                integer (c_int)    :: nfname
                character (c_char) :: filename
                character (c_char) :: mode
                integer (c_int)    :: step
                integer (c_int)    :: ndofs
                real (c_double)    :: u(*)
        end subroutine
end interface
interface
        subroutine load_binstep(nfname, filename, step, ndofs, u, st) &
                                bind(C, name="load_binstep")
                use iso_c_binding, only : c_int, c_char, c_ptr
                integer (c_int)    :: nfname
                character (c_char) :: filename
                integer (c_int)    :: step
                integer (c_int)    :: ndofs
                type(c_ptr)        :: u
                integer (c_int)    :: st
        end subroutine
end interface


interface
        subroutine save_bintparams(nfname, filename, deltat, iointerval) &
                                   bind(C, name="save_bintparams")
                use iso_c_binding, only : c_int, c_char, c_double
                integer (c_int)    :: nfname
                character (c_char) :: filename
                real(c_double)     :: deltat
                integer(c_int)     :: iointerval
        end subroutine
end interface


interface
        subroutine load_bintparams(nfname, filename, deltat, iointerval) &
                                   bind(C, name="load_bintparams")
                use iso_c_binding, only : c_int, c_char, c_double
                integer (c_int)    :: nfname
                character (c_char) :: filename
                real(c_double)     :: deltat
                integer(c_int)     :: iointerval
        end subroutine
end interface

interface
        subroutine visual_params(vtype, params) &
                        bind(C, name="visual_params")
                use iso_c_binding, only : c_double, c_char
                character (c_char) :: vtype
                real(c_double)     :: params(*)
        end subroutine
end interface

interface
        subroutine open_cvm_pipe(npath, path) &
                      bind(C, name="open_cvm_pipe")
                use iso_c_binding, only : c_int, c_char
                integer   (c_int)  :: npath
                character (c_char) :: path
        end subroutine
end interface

interface
        subroutine cvm_mat(point, mat) &
                        bind(C, name="cvm_mat")
                use iso_c_binding, only : c_int, c_double
                real (c_double) :: point(*)
                integer (c_int) :: mat
        end subroutine
end interface

interface
        subroutine close_cvm_pipe() &
                        bind(C, name="close_cvm_pipe")
        end subroutine
end interface

interface
        subroutine write_vtk_unst_field(fnlen, filename, nelems, eptr,     &
                                        eind, nnodes, nodes, nfldname,     &
                                        fldname, field)                    &
                                        bind(C, name="write_vtk_unst_field")
                use iso_c_binding, only : c_int, c_char, c_double
                integer(c_int)    :: fnlen
                character(c_char) :: filename
                integer(c_int)    :: nelems
                integer(c_int)    :: eptr(*)
                integer(c_int)    :: eind(*)
                integer(c_int)    :: nnodes
                real(c_double)    :: nodes(*)
                integer(c_int)    :: nfldname
                character(c_char) :: fldname
                real(c_double)    :: field(*)
        end subroutine
end interface

interface
        subroutine write_sheet_close(nfn, fname, inttime) &
                   bind(C, name="write_sheet_close")
                use iso_c_binding, only : c_int, c_char, c_double
                integer(c_int)    :: nfn
                character(c_char) :: fname
                real(c_double) :: inttime 
        end subroutine
end interface

interface
        subroutine write_sheet_step(nfn, fname, n, field) &
                   bind(C, name="write_sheet_step")
                use iso_c_binding, only : c_int, c_char, c_double
                integer(c_int)    :: nfn
                character(c_char) :: fname
                integer(c_int)    :: n
                real(c_double)    :: field(*)
        end subroutine
end interface

interface
        subroutine write_sheet_points(nfn, fname, n, points) &
                   bind(C, name="write_sheet_points")
                use iso_c_binding, only : c_int, c_char, c_double
                integer(c_int)    :: nfn
                character(c_char) :: fname
                integer(c_int)    :: n
                real(c_double)    :: points(*)
        end subroutine
end interface

contains

subroutine save_mesh(proc, eptr, eind, ediv, nodes)
        integer, dimension(:), intent(in) :: eptr, eind
        logical*1, dimension(:), intent(in) :: ediv
        real*8, dimension(:,:), intent(in) :: nodes
        integer, intent(in) :: proc
        character(len=256) :: filename
        integer :: fnlen
        fnlen = io_fnodlen + 13
        write(filename,"(A,A6,I3.3,A4)") io_outdir(1:io_fnodlen),"/mesh.",proc,".bin"
        call save_binmesh(fnlen, filename, size(eptr,1)-1, eptr, eind,&
                          ediv, size(nodes,2), nodes)
end subroutine

subroutine load_mesh(proc, eptr, eind, ediv, nodes)
        integer, intent(in) :: proc
        integer, dimension(:), allocatable, intent(inout)   :: eptr, eind
        logical*1, dimension(:), allocatable, intent(inout) :: ediv
        real*8, dimension(:,:), allocatable, intent(inout) :: nodes
        type(c_ptr) :: peptr, peind, pediv, pnodes
        integer :: nelems, nnodes
        character(len=256) :: filename
        integer :: fnlen
        fnlen = io_fnodlen + 13
        write(filename,"(A,A6,I3.3,A4)") io_outdir(1:io_fnodlen),"/mesh.",proc,".bin"
        call load_binmesh(fnlen, filename, nelems, peptr,  &
                          peind, pediv, nnodes, pnodes)
        if ( allocated(eptr) ) then
                deallocate(eptr)
        end if
        if ( allocated(eind) ) then
                deallocate(eind)
        end if
        if ( allocated(ediv) ) then
                deallocate(ediv)
        end if
        if ( allocated(nodes) ) then
                deallocate(nodes)
        end if
        call ctofort(peptr, eptr, (/nelems+1/))
        call ctofort(peind, eind, (/eptr(nelems+1)/))
        call ctofort(pediv, ediv, (/eptr(nelems+1)/))
        call ctofort(pnodes, nodes, (/2, nnodes/))
end subroutine

subroutine save_mass(proc, mass)
        integer, intent(in)  :: proc
        real*8, dimension(:) :: mass
        character(len=256) :: filename
        integer :: fnlen
        fnlen = io_fnodlen + 13
        write(filename,"(A,A6,I3.3,A4)") io_outdir(1:io_fnodlen),"/mass.",proc,".bin"
        call save_binmass(fnlen, filename, size(mass,1), mass)
end subroutine

subroutine load_mass(proc, mass)
        integer, intent(in) :: proc
        real*8, dimension(:), allocatable, intent(inout)  :: mass
        type(c_ptr) :: pmass
        integer :: ndofs
        character(len=256) :: filename
        integer :: fnlen
        fnlen = io_fnodlen + 13
        write(filename,"(A,A6,I3.3,A4)") io_outdir(1:io_fnodlen),"/mass.",proc,".bin"
        call load_binmass(fnlen, filename, ndofs, pmass)
        if ( allocated(mass) ) then
                deallocate(mass)
        end if
        call ctofort(pmass, mass, (/ndofs/))
end subroutine

subroutine save_tparams(proc, deltat, iointerval)
        integer, intent(in)  :: proc
        real*8 :: deltat
        integer :: iointerval
        character(len=256) :: filename
        integer :: fnlen
        fnlen = io_fnodlen + 13
        write(filename,"(A,A6,I3.3,A4)") io_outdir(1:io_fnodlen),"/tpar.",proc,".bin"
        call save_bintparams(fnlen, filename, deltat, iointerval)
end subroutine

subroutine load_tparams(proc, deltat, iointerval)
        integer, intent(in) :: proc
        real*8 :: deltat
        integer :: iointerval
        character(len=256) :: filename
        integer :: fnlen
        fnlen = io_fnodlen + 13
        write(filename,"(A,A6,I3.3,A4)") io_outdir(1:io_fnodlen),"/tpar.",proc,".bin"
        call load_bintparams(fnlen, filename, deltat, iointerval)
end subroutine

subroutine save_step(proc, step, u, timer)
        integer, intent(in) :: proc, step
        real*8, intent(in), dimension(:) :: u
        integer, intent(inout) :: timer
        character(len=256) :: filename
        character(len=27)  :: elstr
        integer :: fnlen
        fnlen = io_fnodlen + 13
        if ( mod(step, io_interval) == 0 ) then
                write(filename,"(A,A6,I3.3,A4)") &
                      io_outdir(1:io_fnodlen),"/step.",proc,".bin"
                if ( io_newsave ) then
                        call save_binstep(fnlen, filename, "n", step/io_interval, size(u,1), u)
                        io_newsave = .FALSE.
                else
                        call save_binstep(fnlen, filename, "o", step/io_interval, size(u,1), u)
                end if
                call end_timer(elstr, timer)
                if ( proc == 0 ) then
                write(*,'(A9,I6,A28)') "In step ", step, elstr
                end if
                call start_timer(timer)
        end if
end subroutine

subroutine load_step(proc, step, u, st)
        integer, intent(in) :: proc, step
        real*8, dimension(:), allocatable, intent(inout)  :: u
        integer, intent(out) :: st
        type(c_ptr) :: pu
        integer :: ndofs
        character(len=256) :: filename
        integer :: fnlen
        fnlen = io_fnodlen + 13
        write(filename,"(A,A6,I3.3,A4)") io_outdir(1:io_fnodlen),"/step.",proc,".bin"
        call load_binstep(fnlen, filename, step, ndofs, pu, st)
        if ( st == -1 ) then
                return
        end if
        if ( allocated(u) ) then
                deallocate(u)
        end if
        call ctofort(pu, u, (/ndofs/))
end subroutine

subroutine print_vtkun_locfield(proc, step, eptr, eind, nodes, u)
        integer, intent(in) :: proc, step
        integer, dimension(:), intent(in)  :: eptr, eind
        real*8, dimension(:), intent(in)   :: u
        real*8, dimension(:,:), intent(in) :: nodes
        character(len=256) :: filename
        integer :: fnlen
        fnlen = io_fnvdlen + 20 
        write(filename,"(A,A7,I3.3,A1,I5.5,A4)") & 
              io_visdir(1:io_fnvdlen),"/exfem.",proc,".",step,".vtu"
        call write_vtk_unst_field(fnlen, filename, size(eptr,1)-1,  &
                                  eptr, eind, size(nodes,2), nodes, &
                                  13, "displacements", u)
end subroutine

subroutine print_vtkun_field(step, eptr, eind, nodes, u)
        integer, intent(in) :: step
        integer, dimension(:), intent(in)  :: eptr, eind
        real*8, dimension(:), intent(in)   :: u
        real*8, dimension(:,:), intent(in) :: nodes
        character(len=256) :: filename
        integer :: fnlen
        fnlen = io_fnvdlen + 16
        write(filename,"(A,A7,I5.5,A4)") io_visdir(1:io_fnvdlen),"/exfem.",step,".vtu"
        call write_vtk_unst_field(fnlen, filename, size(eptr,1)-1, eptr,&
                                  eind, size(nodes,2), nodes, 13,       &
                                  "displacements", u)
end subroutine

subroutine print_vtkun_mater(eptr, eind, nodes, u)
        integer, dimension(:), intent(in)  :: eptr, eind
        real*8, dimension(:), intent(in)   :: u
        real*8, dimension(:,:), intent(in) :: nodes
        character(len=256) :: filename
        integer :: fnlen
        fnlen = io_fnvdlen + 10
        write(filename,"(A,A10)") io_visdir(1:io_fnvdlen),"/mater.vtu"
        call write_vtk_unst_field(fnlen, filename, size(eptr,1)-1, eptr,      &
                                  eind, size(nodes,2), nodes, 8, "material", u)
end subroutine

subroutine print_vtkst_field(step, dims, origin, space, u)
        integer, intent(in) :: step
        integer, dimension(2), intent(in) :: dims
        real*8, dimension(2), intent(in) :: origin, space
        real*8, dimension(:), intent(in) :: u
        character(len=24) :: filename
        write(filename,"(A13,I5.5,A4)") &
              "output/exfem.",step,".vtk"
        call savevtk_struct_field(filename, dims, origin, space, u)
end subroutine

subroutine save_sheet_points(points)
        real*8, dimension(:,:), intent(in) :: points
        character(len=256) :: filename
        integer :: fnlen
        fnlen = io_fnvdlen + 10
        write(filename,"(A,A10)") io_visdir(1:io_fnvdlen),"/sheet.txt"
        call write_sheet_points(fnlen, filename, 2*size(points,2), points)
end subroutine

subroutine save_sheet_step(field)
        real*8, dimension(:), intent(in) :: field
        character(len=256) :: filename
        integer :: fnlen
        fnlen = io_fnvdlen + 10
        write(filename,"(A,A10)") io_visdir(1:io_fnvdlen),"/sheet.txt"
        call write_sheet_step(fnlen, filename, size(field), field)
end subroutine

subroutine save_sheet_close(inttime)
        real*8 :: inttime
        character(len=256) :: filename
        integer :: fnlen
        fnlen = io_fnvdlen + 10
        write(filename,"(A,A10)") io_visdir(1:io_fnvdlen),"/sheet.txt"
        call write_sheet_close(fnlen, filename, inttime)
end subroutine

subroutine savevtk_mesh(filename, ptr, ind, points)
        integer, dimension(:), intent(in)   :: ptr, ind
        real*8, dimension(:,:), intent(in)  :: points
        character(len=*), intent(in) :: filename
        integer :: i, j
        integer, save :: cons = 0 
        
        cons = cons + 1
        open(5, file=filename)
        write(5, "(A)") "# vtk DataFile Version 2.0"
        write(5, "(A)") "Explicit FEM Simulation"
        write(5, "(A)") "ASCII"
        write(5, "(A)") "DATASET UNSTRUCTURED_GRID"
        write(5, "(A6,1X,I0,1X,A5)") "POINTS",size(points,2),"float"
        do i=1,size(points,2),4
                write(5, "(4(3(F12.3,1X),2X))") (/( (/points(:,j), 0.0D0/), &
                                                 j=i,min(i+3,size(points,2) ) )/)
        end do
        write(5, "(A5,1X,I0,1X,I0)") "CELLS",size(ptr,1)-1,(size(ind,1)+size(ptr,1)-1)
        do i=1,size(ptr,1)-1,3
                write(5, "(3(5(I0,1X),2X))") (/( (/(ptr(j+1)-ptr(j)), ind(ptr(j):(ptr(j+1)-1))-1/), &
                                               j=i,min(i+2,size(ptr,1)-1) )/)
        end do
        write(5, "(A10,1X,I0)") "CELL_TYPES", size(ptr,1)-1
        do i=1,size(ptr,1)-1,50
                write(5, "(50(I0,1X))") (/ (7, j=i,min(i+49,size(ptr,1)-1) ) /)
        end do
        
        write(5, "(A,1X,I0)") "CELL_DATA",size(ptr,1)-1
        write(5, "(A)") "SCALARS prop int 1"
        write(5, "(A)") "LOOKUP_TABLE default"
        do i=1,size(ptr,1)-1,50
                write(5, "(50(I0,1X))") (/ (0, j=i,min(i+49,size(ptr,1)-1) ) /)
        end do
        close(5)
end subroutine

subroutine savevtk_unst_prop(filename, ptr, ind, points, field)
        integer, dimension(:), intent(in)   :: ptr, ind
        real*8, dimension(:,:), intent(in)  :: points
        integer, dimension(:), intent(in)   :: field
        character(len=*), intent(in) :: filename
        integer :: i, j
        integer, save :: cons = 0 
        cons = cons + 1
        open(5, file=filename)
        write(5, "(A)") "# vtk DataFile Version 2.0"
        write(5, "(A)") "Explicit FEM Simulation"
        write(5, "(A)") "ASCII"
        write(5, "(A)") "DATASET UNSTRUCTURED_GRID"
        write(5, "(A6,1X,I0,1X,A5)") "POINTS",size(points,2),"float"
        do i=1,size(points,2),4
                write(5, "(4(3(F12.3,1X),2X))") (/( (/points(:,j), 0.0D0/), &
                                                 j=i,min(i+3,size(points,2) ) )/)
        end do
        write(5, "(A5,1X,I0,1X,I0)") "CELLS",size(ptr,1)-1,(size(ind,1)+size(ptr,1)-1)
        do i=1,size(ptr,1)-1,3
                write(5, "(3(5(I0,1X),2X))") (/( (/(ptr(j+1)-ptr(j)), ind(ptr(j):(ptr(j+1)-1))-1/), &
                                               j=i,min(i+2,size(ptr,1)-1) )/)
        end do
        write(5, "(A10,1X,I0)") "CELL_TYPES", size(ptr,1)-1
        do i=1,size(ptr,1)-1,50
                write(5, "(50(I0,1X))") (/ (7, j=i,min(i+49,size(ptr,1)-1) ) /)
        end do
        
        write(5, "(A,1X,I0)") "CELL_DATA",size(ptr,1)-1
        write(5, "(A)") "SCALARS prop int 1"
        write(5, "(A)") "LOOKUP_TABLE default"
        do i=1,size(ptr,1)-1,50
                write(5, "(50(I0,1X))") (/ (field(j), j=i,min(i+49,size(ptr,1)-1) ) /)
        end do
        close(5)
end subroutine

subroutine savevtk_unst_field2(filename, ptr, ind, points, field)
        integer, dimension(:), intent(in)   :: ptr, ind
        real*8, dimension(:,:), intent(in)  :: points
        real*8, dimension(:,:), intent(in)   :: field
        character(len=*), intent(in) :: filename
        integer :: i, j
        integer, save :: cons = 0 
        cons = cons + 1
        open(5, file=filename)
        write(5, "(A)") "# vtk DataFile Version 2.0"
        write(5, "(A)") "Explicit FEM Simulation"
        write(5, "(A)") "ASCII"
        write(5, "(A)") "DATASET UNSTRUCTURED_GRID"
        write(5, "(A6,1X,I0,1X,A5)") "POINTS",size(points,2),"float"
        do i=1,size(points,2),4
                write(5, "(4(3(F12.3,1X),2X))") (/( (/points(:,j), 0.0D0/), &
                                                 j=i,min(i+3,size(points,2) ) )/)
        end do
        write(5, "(A5,1X,I0,1X,I0)") "CELLS",size(ptr,1)-1,(size(ind,1)+size(ptr,1))
        do i=1,size(ptr,1)-1,3
                write(5, "(3(5(I0,1X),2X))") (/( (/(ptr(j+1)-ptr(j)), ind(ptr(j):(ptr(j+1)-1))-1/), &
                                               j=i,min(i+2,size(ptr,1)-1) )/)
        end do
        write(5, "(A10,1X,I0)") "CELL_TYPES", size(ptr,1)-1
        do i=1,size(ptr,1)-1,80
                write(5, "(80(I0,1X))") (/ (7, j=i,min(i+79,size(ptr,1)-1) ) /)
        end do
        write(5, "(A,1X,I0)") "POINT_DATA",size(points,2)
        write(5, "(A)") "VECTORS U float"
        do i=1,size(points,2),3
                write(5, "(4(3F12.2,2X))") (/ ( (/real(field(1,j), kind=4), real(field(2,j), kind=4), 0.0/), &
                                                  j=i,min(i+2,size(points,2)) ) /)
        end do
        close(5)
end subroutine

subroutine savevtk_unst_field(filename, ptr, ind, points, field)
        integer, dimension(:), intent(in)   :: ptr, ind
        real*8, dimension(:,:), intent(in)  :: points
        real*8, dimension(:), intent(in)   :: field
        character(len=*), intent(in) :: filename
        integer :: i, j
        open(5, file=filename)
        write(5, "(A)") "# vtk DataFile Version 2.0"
        write(5, "(A)") "Explicit FEM Simulation"
        write(5, "(A)") "ASCII"
        write(5, "(A)") "DATASET UNSTRUCTURED_GRID"
        write(5, "(A6,1X,I0,1X,A5)") "POINTS",size(points,2),"float"
        do i=1,size(points,2),4
                write(5, "(4(3(F12.3,1X),2X))") (/( (/points(:,j), 0.0D0/), &
                                                 j=i,min(i+3,size(points,2) ) )/)
        end do
        write(5, "(A5,1X,I0,1X,I0)") "CELLS",size(ptr,1)-1,(size(ind,1)+size(ptr,1))
        do i=1,size(ptr,1)-1,3
                write(5, "(3(5(I0,1X),2X))") (/( (/(ptr(j+1)-ptr(j)), ind(ptr(j):(ptr(j+1)-1))-1/), &
                                               j=i,min(i+2,size(ptr,1)-1) )/)
        end do
        write(5, "(A10,1X,I0)") "CELL_TYPES", size(ptr,1)-1
        do i=1,size(ptr,1)-1,80
                write(5, "(80(I0,1X))") (/ (7, j=i,min(i+79,size(ptr,1)-1) ) /)
        end do
        write(5, "(A,1X,I0)") "POINT_DATA",size(points,2)
        write(5, "(A)") "VECTORS U float"
        do i=1,size(field,1),8
                write(5, "(4(3(F12.6,1X),2X))") (/ ( (/field(j:j+1), 0.0D0/), &
                                                  j=i,min(i+6,size(field,1)),2 ) /)
        end do
        close(5)
end subroutine

subroutine savevtk_struct_field(filename, dims, origin, space, field)
        integer, dimension(2), intent(in) :: dims
        real*8, dimension(2), intent(in)  :: origin, space
        real*8, dimension(:), intent(in)  :: field
        character(len=*), intent(in) :: filename
        integer :: i, j
        open(5, file=filename)
        write(5, "(A)") "# vtk DataFile Version 2.0"
        write(5, "(A)") "Explicit FEM Simulation"
        write(5, "(A)") "ASCII"
        write(5, "(A)") "DATASET STRUCTURED_POINTS"
        write(5, "(A10,1X,I0,1X,I0,1X,I0)") "DIMENSIONS", dims(1), dims(2), 1
        write(5, "(A10,1X,F12.3,1X,F12.3,1X,F12.3)") "ORIGIN",origin(1),origin(2),0.0D0
        write(5, "(A10,1X,F12.3,1X,F12.3,1X,F12.3)") "SPACING",space(1),space(2),0.0D0
        write(5, "(A,1X,I0)") "POINT_DATA",(dims(1)*dims(2))
        write(5, "(A)") "VECTORS U float"
        do i=1,size(field,1),8
                write(5, "(4(3(F12.6,1X),2X))") (/ ( (/field(j:j+1), 0.0D0/), &
                                                  j=i,min(i+7,size(field,1)),2 ) /)
        end do
        close(5)
end subroutine

subroutine init_cvm()
        if ( io_fncdlen /= -1 ) then
                call open_cvm_pipe(io_fncdlen, io_cvmdir)
        else
                write(*,*) "No CVM model"
        end if        
end subroutine

subroutine end_cvm()
        if ( io_fncdlen /= -1 ) then
                call close_cvm_pipe()
        else
                write(*,*) "No CVM model"
        end if
end subroutine

subroutine load_environ()
        character(len=256) :: filename
        integer :: fnlen
        fnlen = 10
        call load_environ_file(fnlen, "environ.in", io_fnodlen, io_outdir, &
                               io_fncdlen, io_cvmdir, io_fnvdlen, io_visdir)
        write(io_outdir, '(A,A1,A)') io_outdir(1:io_fnodlen), "/", io_indir(1:io_fnindlen)
        write(io_visdir, '(A,A1,A)') io_visdir(1:io_fnvdlen), "/", io_indir(1:io_fnindlen)
        io_fnodlen = len_trim(io_outdir)
        io_fnvdlen = len_trim(io_visdir)
        io_newsave = .TRUE.
end subroutine

subroutine start_timer(timer)
        integer, intent(out) :: timer
        timer = time() 
end subroutine

subroutine end_timer(elstr, timer)
        integer, intent(in) :: timer
        character(len=27) :: elstr
        integer  :: tend
        tend = time()
        write(elstr,'(A14,1X,I5,A1)') "Ellapsed time:", (tend - timer), "s"
end subroutine

end module io_utils 

