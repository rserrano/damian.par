module mesh_utils 
use utilities
use model
implicit none

type elemll
        type(linked_list_il2d),  pointer :: elem   => null()
        type(linked_list_1d),  pointer :: prop   => null()
end type

type elemcsr
        integer, dimension(:),   allocatable :: idx
        integer, dimension(:),   allocatable :: vals
        logical*1, dimension(:),   allocatable :: logi
        integer, dimension(:),   allocatable :: prop
        integer, dimension(:),   allocatable :: gnode
        real*8,  dimension(:,:), allocatable :: nodes
end type

type elemedgell
       type(linked_list_2d), pointer :: elem => null()
       type(linked_list_2d), pointer :: link => null()
end type elemedgell

type elemedge
        integer, dimension(:),   allocatable :: idx
        integer, dimension(:),   allocatable :: vals
        integer, dimension(:,:), allocatable :: links
end type

type edges
        integer, dimension(:,:), allocatable :: links 
        logical*1, dimension(:),   allocatable :: hanging
end type

type division
        integer, dimension(:,:), allocatable :: underedges
        integer, dimension(:), allocatable   :: undergnodes
        integer, dimension(:,:), allocatable :: overedges
        integer, dimension(:), allocatable   :: overgnodes
        real*8, dimension(:,:), allocatable  :: overpoints
end type

type divisionll
        type(linked_list_2d), pointer  :: underedges => null()
        type(linked_list_2d), pointer  :: overedges  => null()
        type(linked_list_1d), pointer  :: overgnodes => null()
        type(linked_list_r2d), pointer :: overpoints => null()
end type

type shnodes
        integer, dimension(:), allocatable :: hgnodes
        integer, dimension(:), allocatable :: annodes
end type

type shnodesll
        type(linked_list_1d), pointer :: hgnodes => null()
        type(linked_list_1d), pointer :: annodes => null()
end type

contains

!Initializes the mesh
subroutine initial_mesh(frow, lrow, fcol, lcol)
        integer, intent(in) :: frow, lrow, fcol, lcol
        integer :: ind, i, j, il, jl
        real*8, dimension(2)   :: pt
        type(linked_list_2d), pointer :: lsabs
        nullify(lsabs)
        
        m_ncols  = lcol-fcol+1
        m_nrows  = lrow-frow+1
        m_nnodes = m_ncols*m_nrows
        m_nelems = (m_nrows-1)*(m_ncols-1)
        
        allocate(m_nodes(2, m_nnodes))
        allocate(m_eptr(m_nelems+1))
        allocate(m_eind(4*m_nelems))
        allocate(m_ediv(4*m_nelems))
        allocate(m_ematprop(m_nelems))
        allocate(m_gnode(m_nnodes))
        ind = 0
        do i=fcol,lcol
                do j=frow,lrow
                        ind = ind + 1
                        m_nodes(:,ind) = (/ real(i-1,kind=8)*m_blength, &
                                            real(j-1,kind=8)*m_blength /)
                        m_gnode(ind) = (i-1)*m_totrows+j
                end do
        end do
        ind=0
        il = 0
        m_eptr(1) = 1
        do i=fcol,lcol-1
                il = il + 1
                jl = 0
                do j=frow,lrow-1
                        jl = jl + 1
                        ind = ind+1
                        m_eptr(ind+1) = m_eptr(ind)+4
                        m_eind(m_eptr(ind):m_eptr(ind+1)-1) = (/(il-1)*m_nrows+jl,    &
                                                                il*m_nrows+jl,        &
                                                                il*m_nrows+(jl+1),    &
                                                                (il-1)*m_nrows+(jl+1)/)
                        m_ediv(m_eptr(ind):m_eptr(ind+1)-1) = (/.FALSE., .FALSE., &
                                                                .FALSE., .FALSE. /)
                        pt = (/(real(i-1,kind=8)+0.5D0)*m_blength, & 
                               (real(j-1,kind=8)+0.5D0)*m_blength /)
                        m_ematprop(ind) = material_prop(pt)
                end do
        end do
end subroutine

subroutine select_local_divide(adjptr, adjncy, crank, csize, outedges, todivide)
        integer, dimension(:), intent(in) :: adjptr, adjncy
        integer, intent(in)               :: crank, csize
        type(edges), dimension(:), intent(out) :: outedges
        logical*1, dimension(:), intent(inout)  :: todivide
        type(ll2dp), dimension(:), allocatable :: lledges
        integer :: i, j, lelem, gelem, proc
        todivide = .FALSE.
        allocate(lledges(csize))
        ! Find the elements that are going to be divided.
        do i = 1, m_nelems
                if ( needs_division(i) .AND. .NOT. hanging_element(i) ) then
                        todivide(i) = .TRUE.
                        do j = adjptr(i), (adjptr(i+1)-1)
                                call local_element(adjncy(j), lelem, proc)
                                if ( proc == crank ) then
                                        if ( hanging_element(lelem) ) then
                                                todivide(i) = .FALSE.
                                                exit
                                        end if
                                end if
                        end do
                        if ( todivide(i) ) then
                                do j = adjptr(i), (adjptr(i+1)-1)
                                        call local_element(adjncy(j), lelem, proc)
                                        if ( proc /= crank ) then
                                                call global_element(i, crank, gelem)
                                                call insertll_2d(lledges(proc+1)%p, &
                                                                 (/gelem, adjncy(j)/))
                                        end if
                                end do
                        end if
                end if
        end do
        do i = 0, csize-1
                if ( i /= crank ) then
                        call toarray_2d(lledges(i+1)%p, outedges(i+1)%links)
                end if
        end do
        deallocate(lledges)
end subroutine

subroutine check_neighbor_divide(crank, csize, outedges, inedges)
        integer, intent(in)               :: crank, csize
        type(edges), dimension(:), intent(inout) :: inedges
        type(edges), dimension(:), intent(inout) :: outedges
        integer :: i, j, lelem, proc
        do i = 0, csize-1
                if ( allocated(inedges(i+1)%links) ) then
                        allocate(outedges(i+1)%hanging(size(inedges(i+1)%links,2)))
                        outedges(i+1)%hanging = .FALSE.
                        do j = 1, size(inedges(i+1)%links,2)
                                call local_element(inedges(i+1)%links(2,j), lelem, proc)
                                if ( proc /= crank ) then
                                        write(*,*) "Received:", inedges(i+1)%links(2,j), &
                                                   "It's not local, local proc: ", crank,&
                                                   "proc destinated: ",proc
                                end if
                                if ( hanging_element(lelem) ) then
                                        outedges(i+1)%hanging(j) = .TRUE.
                                end if
                        end do
                end if
        end do
end subroutine

subroutine update_local_divide(crank, csize, outedges, inedges, todivide)
        integer, intent(in)                   :: crank, csize
        type(edges), dimension(:), intent(in) :: inedges
        type(edges), dimension(:), intent(in) :: outedges
        logical*1, dimension(:), intent(inout)  :: todivide
        integer :: i, j, lelem, proc
        do i = 0, csize-1
                if ( allocated(outedges(i+1)%links) ) then
                        do j = 1, size(outedges(i+1)%links,2) 
                                call local_element(outedges(i+1)%links(1,j), lelem, proc)
                                if ( proc /= crank ) then
                                        write(*,*) "Received:", outedges(i+1)%links(1,j), &
                                                   "It's not local, local proc: ", crank,&
                                                   "proc destinated: ",proc
                                end if
                                if ( inedges(i+1)%hanging(j) ) then
                                        todivide(lelem) = .FALSE.
                                end if
                        end do
                end if
        end do
end subroutine

subroutine create_centernodes(todivide, newnodes, correlems)
        logical*1, dimension(:), intent(in)  :: todivide
        real*8, dimension(:,:), allocatable, intent(out) :: newnodes
        integer, dimension(:), allocatable, intent(out)  :: correlems
        integer :: i
        integer, dimension(4) :: nds
        real*8, dimension(2,4) :: pts
        real*8, dimension(2) :: point
        type(linked_list_r2d), pointer :: llnewnodes
        type(linked_list_1d), pointer  :: llcorrelems
        nullify(llnewnodes)
        nullify(llcorrelems)
        do i = 1, m_nelems
                if ( todivide(i) ) then
                        ! First, create the middle node.
                        nds = m_eind(m_eptr(i):(m_eptr(i+1)-1))
                        pts(:,:) = m_nodes(:,nds)
                        point = (/sum(pts(1,:)), sum(pts(2,:))/)
                        point = point/4.0D0
                        call insertll(llnewnodes, point)
                        call insertll(llcorrelems, i)
                end if
        end do
        call toarray(llnewnodes, newnodes)
        ! write(*,*) newnodes
        call toarray(llcorrelems, correlems)
        ! write(*,*) correlems
end subroutine

subroutine create_singlenodes(adjptr, corresp, todivide, newnodes, corredges)
        integer, dimension(:), intent(in)  :: adjptr, corresp
        logical*1, dimension(:), intent(in)  :: todivide
        real*8, dimension(:,:), allocatable :: newnodes
        integer, dimension(:,:), allocatable :: corredges
        integer :: i, j
        integer, dimension(2) :: nds
        integer, dimension(4) :: reverse
        real*8, dimension(2,2) :: pts
        real*8, dimension(2) :: point
        type(linked_list_r2d), pointer :: llnewnodes
        type(linked_list_2d), pointer :: llcorredges
        nullify(llnewnodes)
        nullify(llcorredges)
        do i = 1, m_nelems
        if ( todivide(i) ) then
                reverse = 0
                do j = adjptr(i), adjptr(i+1)-1
                        reverse(corresp(j)) = (j-adjptr(i))+1
                end do
                do j = 1, 4
                if ( reverse(j) == 0 ) then
                        call elementnodes(i, j, 2, nds)
                        pts = m_nodes(:,nds)
                        point = (/sum(pts(1,:)), sum(pts(2,:))/)/2.0D0
                        call insertll(llnewnodes, point)
                        call insertll(llcorredges, (/i, j/))
                end if
                end do
        end if
        end do
        call toarray(llnewnodes, newnodes)
        ! write(*,*) newnodes
        call toarray(llcorredges, corredges)
        ! write(*,*) corredges
end subroutine

subroutine select_underedges(xadj, adjncy, todivide, csize, crank, outdiv)
        integer, dimension(:), intent(in)  :: xadj, adjncy
        logical*1, dimension(:), intent(in)  :: todivide
        integer, intent(in) :: csize, crank
        type(division), dimension(:), intent(inout) :: outdiv
        type(divisionll), dimension(:), allocatable :: divll
        integer :: i, j, global, local, proc
        allocate(divll(csize))
        do i = 1, m_nelems
                if ( todivide(i) ) then
                        call global_element(i, crank, global)
                        
                        do j = xadj(i), xadj(i+1)-1
                                if ( global > adjncy(j) ) then
                                        call local_element(adjncy(j),local,proc)
                                        call insertll(divll(proc+1)%underedges,&
                                                      (/global, adjncy(j)/))
                                end if
                        end do
                end if
        end do
        do i = 1, csize
                
                call toarray(divll(i)%underedges, outdiv(i)%underedges)
                if ( allocated(outdiv(i)%underedges) ) then
                        call sort(outdiv(i)%underedges)
                end if
        end do
        deallocate(divll)
end subroutine

subroutine create_edgenodes(xadj, adjncy, corresp, todivide, csize, &
                            crank, outdiv, indiv, uppoints)
        integer, dimension(:), intent(in)  :: xadj, adjncy, corresp
        logical*1, dimension(:), intent(in)  :: todivide
        integer, intent(in) :: csize, crank
        type(division), dimension(:), intent(inout) :: outdiv, indiv
        real*8, dimension(:,:), allocatable, intent(inout) :: uppoints
        type(divisionll), dimension(:), allocatable :: divll
        type(linked_list_2d), pointer :: llupedges
        type(linked_list_r2d), pointer :: lluppoints
        integer, dimension(:,:), allocatable :: upedges
        integer :: i, j, k, ind, uind, global, local, proc
        integer, dimension(2)  :: nds
        real*8, dimension(2,2) :: pts
        real*8, dimension(2)   :: point
        logical*1 :: found
        nullify(llupedges)
        nullify(lluppoints)
        allocate(divll(csize))
        do i = 1, csize
                if (allocated(indiv(i)%underedges)) then
                        allocate(indiv(i)%undergnodes(size(indiv(i)%underedges,2)))
                        indiv(i)%undergnodes = 0
                end if
        end do
        ind = 0
        do i = 1, m_nelems
        if ( todivide(i) ) then
                call global_element(i, crank, global)
                do j = xadj(i), xadj(i+1)-1
                if ( global < adjncy(j) ) then
                        call local_element(adjncy(j), local, proc)
                        ind = ind + 1
                        call insertll(divll(proc+1)%overgnodes,ind)
                        
                        if ( allocated(indiv(proc+1)%underedges) ) then
                                call binsearch(indiv(proc+1)%underedges,(/adjncy(j),global/),uind)
                        else
                                uind = 0
                        end if
                        
                        if ( uind /= 0 ) then
                                indiv(proc+1)%undergnodes(uind) = ind
                        end if
                        
                        call insertll(llupedges,(/global, adjncy(j)/))
                        call elementnodes(i, corresp(j), 2, nds)
                        pts = m_nodes(:,nds)
                        point = (/sum(pts(1,:)), sum(pts(2,:))/)/2.0D0
                        call insertll(lluppoints, point)
                end if
                end do
        end if
        end do
        
        do i = 1, csize
        if ( allocated(indiv(i)%underedges) ) then
                do j = 1, size(indiv(i)%underedges,2)
                if ( indiv(i)%undergnodes(j) == 0 ) then
                        ind = ind + 1
                        call insertll(divll(i)%overgnodes,ind)
                        indiv(i)%undergnodes(j) = ind
                        call insertll(llupedges,indiv(i)%underedges((/2,1/),j))
                        call local_element(indiv(i)%underedges(2,j), local, proc)
                        if ( proc /= crank ) then
                                write(*,*) "Error: ", __FILE__, ": ", __LINE__
                        end if
                        found = .FALSE.
                        do k = xadj(local), xadj(local+1)-1
                        if ( adjncy(k) == indiv(i)%underedges(1,j) ) then
                                found = .TRUE.
                                exit
                        end if
                        end do
                        if ( .NOT. found ) then
                                write(*,*) "Error: ", __FILE__, ": ", __LINE__
                        end if
                        call elementnodes(local, corresp(k), 2, nds)
                        pts = m_nodes(:,nds)
                        point = (/sum(pts(1,:)), sum(pts(2,:))/)/2.0D0
                        call insertll(lluppoints, point)
                end if
                end do
        end if
        end do
        call toarray(llupedges, upedges)
        call toarray(lluppoints, uppoints)
        if ( allocated(upedges) .and. size(upedges) /= size(uppoints) ) then
                write(*,*) "Error at: ", __FILE__, ":", __LINE__
        end if
        do i = 1, csize
        call toarray(divll(i)%overgnodes, outdiv(i)%overgnodes)
        if ( allocated( outdiv(i)%overgnodes ) ) then
                if ( any(outdiv(i)%overgnodes <= 0 ) .OR.&
                     any(outdiv(i)%overgnodes > size(upedges,2) ) ) then
                        write(*,*) "Error at: ", __FILE__, ":", __LINE__
                end if
                allocate(outdiv(i)%overedges(2,size(outdiv(i)%overgnodes,1)))
                allocate(outdiv(i)%overpoints(2,size(outdiv(i)%overgnodes,1)))
                outdiv(i)%overedges(:,:)  = upedges(:,outdiv(i)%overgnodes)
                outdiv(i)%overpoints(:,:) = uppoints(:,outdiv(i)%overgnodes)
                call sort(outdiv(i)%overedges, outdiv(i)%overgnodes, &
                          outdiv(i)%overpoints)
        end if
        end do
        if ( allocated(upedges) ) then
                deallocate(upedges)
        end if
        deallocate(divll)
end subroutine

subroutine select_newhanging(todivide, csize, crank, outdivs, indivs, newhanging)
        logical*1, dimension(:), intent(in) :: todivide
        integer, intent(in) :: csize, crank
        type(division), dimension(:), intent(in) :: outdivs, indivs
        logical*1, dimension(:), allocatable, intent(out) :: newhanging
        integer :: i, j, local, proc
        allocate(newhanging(m_nelems))
        newhanging = .FALSE.
        do i = 1, csize
        if ( (i-1) /= crank ) then
                if ( allocated(outdivs(i)%overgnodes) ) then
                        do j = 1, size(outdivs(i)%overgnodes)
                        call local_element(outdivs(i)%overedges(1,j),local,proc)
                        if ( .NOT. todivide(local) ) then
                                newhanging(local) = .TRUE.
                        end if
                        end do
                end if
                if ( allocated(indivs(i)%overgnodes) ) then
                        do j = 1, size(indivs(i)%overgnodes)
                        call local_element(indivs(i)%overedges(2,j),local,proc)
                        if ( .NOT. todivide(local) ) then
                                newhanging(local) = .TRUE.
                        end if
                        end do
                end if
        else
                if (allocated(outdivs(i)%overgnodes) ) then
                        do j = 1, size(outdivs(i)%overgnodes)
                        call local_element(outdivs(i)%overedges(1,j),local,proc)
                        if ( .NOT. todivide(local) ) then
                                newhanging(local) = .TRUE.
                        end if
                        call local_element(outdivs(i)%overedges(2,j),local,proc)
                        if ( .NOT. todivide(local) ) then
                                newhanging(local) = .TRUE.
                        end if
                        end do
                end if
        end if
        end do
end subroutine

subroutine final_division(adjptr, adjncy, corresp, nodes, gnode,&
                          elems, sinedges, csize, crank, s1, s2,&
                          todivide, ishanging, outdivs, indivs)
        integer, dimension(:), intent(in) :: adjptr, adjncy, corresp
        real*8, dimension(:,:), allocatable, intent(inout) :: nodes
        integer, dimension(:), allocatable,  intent(inout) :: gnode
        integer, dimension(:), intent(in) :: elems
        integer, dimension(:,:), intent(in) :: sinedges
        integer, intent(in) :: csize, crank, s1, s2
        logical*1, dimension(:), intent(in) :: todivide, ishanging
        type(division), dimension(:), intent(inout) :: outdivs, indivs
        type(linked_list_il2d), pointer :: llnewelems
        type(linked_list_1d), pointer   :: llmatprop
        integer, dimension(:), allocatable :: neptr, neind, matprop
        logical*1, dimension(:), allocatable :: nelogi
        integer :: i, j, local, global, proc, idx, nidx, cidx, cgno, nndsh, ridx
        integer, dimension(4)  :: reverse, toinsert
        integer, dimension(8)  :: hangelem
        logical*1, dimension(8)  :: hanglogi
        real*8, dimension(2,4) :: pelem
        real*8 :: diff
        m_eind = m_gnode(m_eind)
        nullify(llnewelems)
        nullify(llmatprop)
        do i = 1, m_nelems
        if ( ishanging(i) ) then
                toinsert = 0
                call global_element(i, crank, global)
                do j = adjptr(i), adjptr(i+1)-1
                call local_element(adjncy(j), local, proc)
                
                if ( crank == proc ) then
                        idx = 0
                        if ( global < adjncy(j) ) then
                                if ( allocated(outdivs(proc+1)%overedges) ) then
                                        call binsearch(outdivs(proc+1)%overedges,& 
                                                       (/global, adjncy(j)/),&
                                                       idx)
                                end if
                        else
                                if ( allocated(outdivs(proc+1)%overedges) ) then
                                        call binsearch(outdivs(proc+1)%overedges,&
                                                       (/adjncy(j), global/),&
                                                       idx)
                                end if
                        end if
                        if ( idx /= 0 ) then
                                if ( .not. allocated(outdivs(proc+1)%overgnodes) ) then
                                        write(*,*) "ERROR: ", __FILE__, __LINE__
                                end if
                                if ( size(outdivs(proc+1)%overgnodes) < idx ) then
                                        write(*,*) "ERROR: ", __FILE__, __LINE__
                                end if

                                nidx = outdivs(proc+1)%overgnodes(idx)
                        end if
                else
                        idx = 0
                        if ( global < adjncy(j) ) then
                                if ( allocated(outdivs(proc+1)%overedges) ) then
                                        call binsearch(outdivs(proc+1)%overedges,&
                                                       (/global, adjncy(j)/),&
                                                       idx)
                                end if
                                if ( idx /= 0 ) then
                                        nidx = outdivs(proc+1)%overgnodes(idx)
                                end if
                        else
                                if ( allocated(indivs(proc+1)%overedges) ) then
                                        call binsearch(indivs(proc+1)%overedges, &
                                                       (/adjncy(j), global/),&
                                                       idx)
                                end if
                                if ( idx /= 0 ) then
                                        nidx = indivs(proc+1)%overgnodes(idx)
                                end if
                        end if
                        
                end if 
                if ( idx /= 0 ) then
                        toinsert(corresp(j)) = nidx
                end if
                
                end do
                hangelem(1:7:2) = m_eind(m_eptr(i):m_eptr(i+1)-1)
                hangelem(2:8:2) = toinsert(1:4)
                hanglogi(1:7:2) = .FALSE.
                hanglogi(2:8:2) = .TRUE.
                nndsh = count(hangelem /= 0)
                if ( count(toinsert /= 0) == 0 ) then
                        write(*,*) "Error: ", __FILE__, ":", __LINE__
                end if
                hanglogi(1:nndsh) = pack(hanglogi, hangelem /= 0)
                hangelem(1:nndsh) = pack(hangelem, hangelem /= 0)
                if ( count(.not. hanglogi(1:nndsh) ) /= 4 ) then
                        write(*,*) "Error: ", __FILE__, ":", __LINE__, &
                                    count( .not. hanglogi(1:nndsh) ),  &
                                    count( hanglogi(1:nndsh) )
                end if
                call insertll(llnewelems, hangelem(1:nndsh), hanglogi(1:nndsh))
                call insertll(llmatprop, m_ematprop(i))
        end if
        if ( todivide(i) ) then

                reverse = 0
                call global_element(i, crank, global)
                
                do j = adjptr(i), adjptr(i+1)-1
                        reverse(corresp(j)) = (j-adjptr(i))+1
                end do
                toinsert = 0
                do j = 1, 4
                
                if ( reverse(j) == 0 ) then
                        call binsearch(sinedges, (/i, j/), idx)
                        if ( idx == 0 ) then
                                write(*,*) "Error: ", __FILE__, ":", __LINE__
                        end if
                        toinsert(j)    = gnode(s2+idx)
                else
                        ridx = adjptr(i)+reverse(j)-1
                        call local_element(adjncy(ridx), local, proc)
                        if ( crank == proc ) then
                                if ( global < adjncy(ridx) ) then
                                        call binsearch(outdivs(proc+1)%overedges,& 
                                                       (/global, adjncy(ridx)/),&
                                                       idx)
                                else
                                        call binsearch(outdivs(proc+1)%overedges,&
                                                       (/adjncy(ridx), global/),&
                                                       idx)
                                end if
                                if ( idx == 0 ) then
                                        write(*,*) "Error: ", __FILE__, ":", __LINE__
                                        write(*,*) global, adjncy(ridx)
                                        write(*,*) todivide(local),&
                                                   ishanging(local)
                                end if
                                nidx = outdivs(proc+1)%overgnodes(idx)
                        else
                                if ( global < adjncy(ridx) ) then
                                        call binsearch(outdivs(proc+1)%overedges,&
                                                       (/global, adjncy(ridx)/),&
                                                       idx)
                                        if ( idx == 0 ) then
                                                write(*,*) "Error: ", __FILE__,&
                                                            __LINE__
                                        end if
                                        nidx = outdivs(proc+1)%overgnodes(idx)
                                else
                                        call binsearch(indivs(proc+1)%overedges, &
                                                       (/adjncy(ridx), global/),&
                                                       idx)
                                        if ( idx == 0 ) then
                                                write(*,*) "Error: ",__FILE__,&
                                                           __LINE__
                                        end if
                                        nidx = indivs(proc+1)%overgnodes(idx)
                                end if
                        
                        end if
                        toinsert(j)    = nidx
                end if
                end do
                
                
                ! These points will be evaluated to know
                ! To which material the new element belongs.
                do j = m_eptr(i), m_eptr(i+1)-1
                        call binsearch(m_gnode, m_eind(j), idx)
                        pelem(:,(j-m_eptr(i))+1) = m_nodes(:,idx)
                end do
                diff = (pelem(2,3)-pelem(2,1))
                pelem(:,1) = pelem(:,1) + (/ diff,  diff/)
                pelem(:,2) = pelem(:,2) + (/-diff,  diff/)
                pelem(:,3) = pelem(:,3) + (/-diff, -diff/)
                pelem(:,4) = pelem(:,4) + (/ diff, -diff/)
                
                call binsearch(elems, i, cidx)
                cgno = gnode(s1+cidx)
                call insertll(llnewelems,&
                              (/m_eind(m_eptr(i)),toinsert(1),cgno,toinsert(4)/),&
                              logical((/.FALSE., .FALSE., .FALSE.,.FALSE./),kind=1))
                call insertll(llmatprop, material_prop(pelem(:,1)))
                call insertll(llnewelems,&
                              (/toinsert(1),m_eind(m_eptr(i)+1),toinsert(2),cgno/),&
                              logical((/.FALSE., .FALSE., .FALSE.,.FALSE./),kind=1))
                call insertll(llmatprop, material_prop(pelem(:,2)))
                call insertll(llnewelems,&
                              (/cgno,toinsert(2),m_eind(m_eptr(i)+2),toinsert(3)/),&
                              logical((/.FALSE., .FALSE., .FALSE.,.FALSE./),kind=1))
                call insertll(llmatprop, material_prop(pelem(:,3)))
                call insertll(llnewelems,&
                              (/toinsert(4),cgno,toinsert(3),m_eind(m_eptr(i)+3)/),&
                              logical((/.FALSE., .FALSE., .FALSE.,.FALSE./),kind=1))
                call insertll(llmatprop, material_prop(pelem(:,4)))
        end if
        end do
        call toarray_csr(llnewelems, neptr, neind, nelogi)
        call toarray(llmatprop, matprop)
        call reduce_csr(m_eptr, m_eind, m_ediv, ishanging .OR. todivide)
        call reduce(m_ematprop, ishanging .OR. todivide)
        if ( allocated(neptr) ) then
                call concat_csr(m_eptr, m_eind, m_ediv, neptr, neind, nelogi)
                call concat(m_ematprop, matprop)
        end if
        if ( allocated(nodes) ) then 
                call concat(m_nodes, nodes)
                call concat(m_gnode, gnode)
        end if
        do i = 1, csize
                if ( allocated(outdivs(i)%overgnodes) ) then
                        deallocate(outdivs(i)%overedges)
                        deallocate(outdivs(i)%overgnodes)
                        deallocate(outdivs(i)%overpoints)
                end if
                if ( allocated(indivs(i)%overgnodes) ) then
                        deallocate(indivs(i)%overedges)
                        call concat(m_nodes, indivs(i)%overpoints)
                        call concat(m_gnode, indivs(i)%overgnodes)
                end if
        end do
        call sort(m_gnode, m_nodes)
        call reverse_indices(m_gnode, m_eind)
        
        m_nelems = size(m_eptr,1)-1
        m_nnodes = size(m_nodes,2)
end subroutine

subroutine add_elems(csize, inelems)
        integer, intent(in) :: csize
        type(elemcsr), dimension(:), allocatable, intent(inout) :: inelems
        integer, dimension(:), allocatable :: indices
        integer :: i
        m_eind = m_gnode(m_eind)
        do i = 1, csize
                if ( allocated( inelems(i)%idx ) ) then
                       
                        if ( size(inelems(i)%logi) /= &
                            (count(.not. inelems(i)%logi) + &
                             count(inelems(i)%logi)) ) then
                                write(*,*) "Error: ", __FILE__, ":", __LINE__,&
                                           (i-1)
                        end if
                        call concat_csr(m_eptr, m_eind, m_ediv, inelems(i)%idx,&
                                        inelems(i)%vals, inelems(i)%logi)
                        call concat(m_ematprop, inelems(i)%prop)
                end if
                if ( allocated( inelems(i)%gnode ) ) then
                        allocate(indices(size(inelems(i)%gnode,1)))
                        call binsearch(m_gnode, inelems(i)%gnode, indices)
                        call reduce(inelems(i)%gnode,&
                                    logical(indices /= 0,kind=1))
                        call reduce(inelems(i)%nodes,&
                                    logical(indices /= 0,kind=1))
                        call concat(m_gnode, inelems(i)%gnode)
                        call concat(m_nodes, inelems(i)%nodes)
                        call sort(m_gnode, m_nodes)
                        deallocate(indices)
                end if
        end do
        call reverse_indices(m_gnode, m_eind)
        m_nelems = size(m_eptr,1)-1
        m_nnodes = size(m_nodes,2)
end subroutine

subroutine cut_geom()
        logical*1, dimension(:), allocatable :: mask
        logical*1, dimension(:), allocatable :: nmask
        integer, dimension(:), allocatable   :: tmpnodes
        integer :: i
        allocate(mask(m_nelems))
        do i = 1, m_nelems
                mask(i) = needs_deletion(i)
        end do
        ! Reduce the elements and indices saved in csr.
        call reduce_csr(m_eptr, m_eind, m_ediv, mask)
        call reduce(m_ematprop, mask)
        deallocate(mask)
        ! Reduce the number of nodes saved locally.
        allocate(nmask(m_nnodes))
        call get_unique(m_eind, tmpnodes)
        if ( .NOT. allocated(tmpnodes) ) then
                allocate(tmpnodes(0))
        end if
        m_eind = m_gnode(m_eind)
        nmask = .TRUE.
        nmask(tmpnodes) = .FALSE.
        call reduce(m_nodes, nmask)
        call reduce(m_gnode, nmask)
        call reverse_indices(m_gnode, m_eind)
        deallocate(nmask)
        deallocate(tmpnodes)
        ! update number of elements.
        m_nnodes = size(m_nodes,2)
        m_nelems = size(m_eptr,1)-1
end subroutine

subroutine reduce_elems(part, crank, csize, outelems)
        integer, dimension(:), intent(in) :: part
        integer, intent(in) :: crank, csize
        type(elemcsr), dimension(:), intent(out) :: outelems
        type(elemll), dimension(:), allocatable  :: llstored
        integer :: i, nnodes
        logical*1, dimension(:),   allocatable :: mask, nmask
        integer, dimension(:),   allocatable :: tmpnodes
        
        allocate(llstored(csize))
        allocate(mask(m_nelems))
        mask = .FALSE.
        do i = 1, m_nelems
                if (part(i)/=crank) then
                        mask(i) = .TRUE.
                        call insertll(llstored(part(i)+1)%elem,       &
                                      m_eind(m_eptr(i):m_eptr(i+1)-1),&
                                      m_ediv(m_eptr(i):m_eptr(i+1)-1))
                        if ( size(m_ediv(m_eptr(i):m_eptr(i+1)-1)) /= &
                             (count(.not. m_ediv(m_eptr(i):m_eptr(i+1)-1)) + &
                             count(m_ediv(m_eptr(i):m_eptr(i+1)-1)) ) ) then
                                write(*,*) "Error: ", __FILE__, ":", __LINE__
                        end if

                        call insertll(llstored(part(i)+1)%prop,   &
                                      m_ematprop(i))
                end if
        end do
        do i = 0, (csize-1)
                if ( crank /= i ) then
                        call toarray_csr(llstored(i+1)%elem,  &
                                         outelems(i+1)%idx,   &
                                         outelems(i+1)%vals,  &
                                         outelems(i+1)%logi)
                        if ( allocated(outelems(i+1)%logi) ) then
                                if ( size(outelems(i+1)%logi) /= &
                                    (count(.not. outelems(i+1)%logi) + &
                                    count(outelems(i+1)%logi)) ) then
                                        write(*,*) "Error: ", __FILE__, ":", __LINE__
                                end if
                        end if

                        call toarray(llstored(i+1)%prop, &
                                     outelems(i+1)%prop)
                        if ( allocated(outelems(i+1)%idx) ) then
                                call get_unique(outelems(i+1)%vals, tmpnodes)
                                !tmpnodes allocated
                                nnodes = size(tmpnodes,1)
                                allocate(outelems(i+1)%nodes(2,nnodes))
                                allocate(outelems(i+1)%gnode(nnodes))
                                outelems(i+1)%nodes(:,:) = m_nodes(:,tmpnodes)
                                outelems(i+1)%gnode(:)   = m_gnode(tmpnodes)
                                deallocate(tmpnodes)
                                outelems(i+1)%vals = m_gnode(outelems(i+1)%vals)
                        end if
                end if
        end do
#ifdef DEBUG
        write(*,*) "Local: ", crank, "Count: ", count(part==0),&
                   count(part==1), count(part==2), count(part==3),&
                   char(10), "sizes out: ", size(outelems(1)%idx), size(outelems(2)%idx), &
                   size(outelems(3)%idx), size(outelems(4)%idx)
#endif
        deallocate(llstored)
        ! Reduce the elements and indices saved in csr.
        call reduce_csr(m_eptr, m_eind, m_ediv, mask)
        call reduce(m_ematprop, mask)
        deallocate(mask)
        ! Reduce the number of nodes saved locally.
        allocate(nmask(m_nnodes))
        call get_unique(m_eind, tmpnodes)
        if ( .NOT. allocated(tmpnodes) ) then
                allocate(tmpnodes(0))
        end if
        
        if ( any(tmpnodes <= 0) .OR. any(tmpnodes > m_nnodes) ) then
                write(*,*) "Error at: ", __FILE__, ":", __LINE__
        end if
        m_eind = m_gnode(m_eind)
        nmask = .TRUE.
        nmask(tmpnodes) = .FALSE.
        call reduce(m_nodes, nmask)
        call reduce(m_gnode, nmask)
        call reverse_indices(m_gnode, m_eind)
        if ( any(m_eind == 0) ) then
                write(*,*) "Error at: ", __FILE__, ":", __LINE__
        end if
        deallocate(nmask)
        deallocate(tmpnodes)
        ! update number of elements.
        m_nnodes = size(m_nodes,2)
        m_nelems = size(m_eptr,1)-1
end subroutine

subroutine create_struct_mesh(dims, origin, space, eptr, eind, nodes)
        integer, dimension(2), intent(in) :: dims
        real*8, dimension(2), intent(in) :: origin, space
        integer, dimension(:), allocatable, intent(out)  :: eptr
        integer, dimension(:), allocatable, intent(out)  :: eind
        real*8, dimension(:,:), allocatable, intent(out) :: nodes
        integer :: nnodes, nelems, ind, i, j
        nnodes = (dims(1)+1)*(dims(2)+1)
        nelems = dims(1)*dims(2)
        
        allocate(nodes(2, nnodes))
        allocate(eptr(nelems+1))
        allocate(eind(4*nelems))
        ind = 0
        do i=0, dims(1)
                do j=0, dims(2)
                        ind = ind + 1
                        nodes(:,ind) = (/  origin(1)+real(i,kind=8)*space(1), &
                                           origin(2)+real(j,kind=8)*space(2) /)
                end do
        end do
        ind=0
        eptr(1) = 1
        do i=1,dims(1)
                do j=1,dims(2)
                        ind = ind+1
                        eptr(ind+1) = eptr(ind)+4
                        eind(eptr(ind):eptr(ind+1)-1) = (/(i-1)*(dims(2)+1)+j,      &
                                                          i*(dims(2)+1)+j,          &
                                                          i*(dims(2)+1)+(j+1),      &
                                                          (i-1)*(dims(2)+1)+(j+1)/)
                end do
        end do
end subroutine

end module mesh_utils

