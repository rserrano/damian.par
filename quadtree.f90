module quadtree
use utilities
use model
use communicator
implicit none
type qtree
        real*8, dimension(2) :: center
        type(qtree), pointer :: downleft  => null()
        type(qtree), pointer :: downright => null()
        type(qtree), pointer :: upright   => null()
        type(qtree), pointer :: upleft    => null()
        integer :: elem
        real*8, dimension(2,2) :: bbox
end type qtree
type(qtree), pointer :: q_elemtree => null()
contains

subroutine create_quadtree()
        integer, dimension(:), allocatable :: elems
        integer :: i
        i = 1
        if ( associated(q_elemtree) ) then
                write(*,*) "ERROR: the main quadtree is already allocated"
        end if
        allocate(q_elemtree)
        allocate(elems(m_nelems))
        elems = (/(i, i=1, m_nelems)/)
        call order_quadtree(q_elemtree, elems)
        call bbox_quadtree(q_elemtree)
end subroutine

recursive subroutine order_quadtree(tree, elems)
        type(qtree), pointer, intent(inout) :: tree
        integer, dimension(:), intent(inout)   :: elems
        real*8, dimension(2) :: mind, maxd, dist
        integer :: i, stpiv, endpiv, midpiv
        integer, dimension(4) :: nodes, enodes
        integer :: nelems
        real*8 :: dlength 
        nelems = size(elems,1)

        if ( nelems == 1 ) then
                tree%elem = elems(1)
                return ! Nothing to do here.
        else if ( size(elems,1) < 1 ) then
                write(*,*) "ERROR:", __FILE__, __LINE__
                return
        else
                tree%elem = 0
        end if
        ! Find maximum and minimum in x and y, and with their distances find
        ! center based on m_blength.
        mind = u_inf;
        maxd = u_ninf;
        do i = 1, nelems
                call elemallnhnodes(elems(i), nodes) 
                mind(1) = min(mind(1), m_nodes(1,nodes(1)))
                mind(2) = min(mind(2), m_nodes(2,nodes(1)))
                maxd(1) = max(maxd(1), m_nodes(1,nodes(3)))
                maxd(2) = max(maxd(2), m_nodes(2,nodes(3)))
        end do
        dist(:) = maxd(:) - mind(:)
        if ( realunder( m_blength, dist(1) ) .OR.&
             realunder( m_blength, dist(2) ) ) then
                mind(1) = realfloor(mind(1)/m_blength)*m_blength
                mind(2) = realfloor(mind(2)/m_blength)*m_blength
                maxd(1) = realceiling(maxd(1)/m_blength)*m_blength
                maxd(2) = realceiling(maxd(2)/m_blength)*m_blength
                dist(1) = realfloor(((maxd(1)-mind(1))/2.0D0)/m_blength)*m_blength
                dist(2) = realfloor(((maxd(2)-mind(2))/2.0D0)/m_blength)*m_blength
                tree%center = mind + dist
        else
                dlength = m_blength/2.0D0
                do while( .NOT. realunder(dlength, dist(1)) .AND. &
                          .NOT. realunder(dlength, dist(2)) )
                        dlength = dlength/2.0D0
                end do 
                mind(1) = realfloor(mind(1)/dlength)*dlength
                mind(2) = realfloor(mind(2)/dlength)*dlength
                maxd(1) = realceiling(maxd(1)/dlength)*dlength
                maxd(2) = realceiling(maxd(2)/dlength)*dlength
                dist(1) = realfloor(((maxd(1)-mind(1))/2.0D0)/dlength)*dlength
                dist(2) = realfloor(((maxd(2)-mind(2))/2.0D0)/dlength)*dlength
                tree%center = mind + dist
        end if
        
        ! With a pivotal like scheme, find the four divisions of the elements.
        ! First the x coordinate.
        
        stpiv = 1
        endpiv = nelems 
        do while ( stpiv <= endpiv )
                call elemallnhnodes(elems(stpiv), nodes)
                call elemallnhnodes(elems(endpiv), enodes)
                if ( realover(tree%center(1), m_nodes(1,nodes(1))) ) then
                        stpiv = stpiv + 1
                        cycle
                end if
                if ( realunder(tree%center(1), m_nodes(1,enodes(3))) ) then
                        endpiv = endpiv - 1
                        cycle
                end if
                elems((/stpiv, endpiv/)) = elems((/endpiv, stpiv/))
                stpiv = stpiv + 1
                endpiv = endpiv - 1
        end do
        midpiv = stpiv
        
        ! Second the y coordinate, first part.
        if ( midpiv > 1 ) then
        stpiv  = 1
        endpiv = midpiv - 1
        do while ( stpiv <= endpiv ) 
                call elemallnhnodes(elems(stpiv), nodes)
                call elemallnhnodes(elems(endpiv), enodes)
                if ( realover(tree%center(2), m_nodes(2,nodes(1))) ) then
                        stpiv = stpiv + 1
                        cycle
                end if
                if ( realunder(tree%center(2), m_nodes(2,enodes(3))) ) then
                        endpiv = endpiv - 1
                        cycle
                end if
                elems((/stpiv, endpiv/)) = elems((/endpiv, stpiv/))
                stpiv = stpiv + 1
                endpiv = endpiv - 1
        end do
        if ( stpiv > 1 ) then
                allocate(tree%downleft)
                call order_quadtree(tree%downleft, elems(1:stpiv-1))
        end if
        if ( stpiv < midpiv ) then
                allocate(tree%upleft)
                call order_quadtree(tree%upleft, elems(stpiv:midpiv-1))
        end if
        end if
        ! Third the y coordinate, second part.
        if ( midpiv <= nelems)  then
        stpiv = midpiv
        endpiv = nelems 
        do while ( stpiv <= endpiv )
                call elemallnhnodes(elems(stpiv), nodes)
                call elemallnhnodes(elems(endpiv), enodes)
                if ( realover(tree%center(2), m_nodes(2,nodes(1))) ) then
                        stpiv = stpiv + 1
                        cycle
                end if
                if ( realunder(tree%center(2), m_nodes(2,enodes(3))) ) then
                        endpiv = endpiv - 1
                        cycle
                end if
                elems((/stpiv, endpiv/)) = elems((/endpiv, stpiv/))
                stpiv = stpiv + 1
                endpiv = endpiv - 1
        end do
        if ( stpiv > midpiv ) then
                allocate(tree%downright)
                call order_quadtree(tree%downright, elems(midpiv:stpiv-1))
        end if
        if ( stpiv <= nelems ) then
                allocate(tree%upright)
                call order_quadtree(tree%upright, elems(stpiv:nelems))
        end if
        end if
        
end subroutine

recursive subroutine bbox_quadtree(tree)
        type(qtree), pointer, intent(in) :: tree
        integer, dimension(4)            :: nodes
        type(qtree), pointer :: dl, dr, ul, ur
        
        if ( tree%elem /= 0 ) then
                call elemallnhnodes(tree%elem, nodes)
                tree%bbox(1,1) = m_nodes(1, nodes(1))
                tree%bbox(2,1) = m_nodes(2, nodes(1))
                tree%bbox(1,2) = m_nodes(1, nodes(3))
                tree%bbox(2,2) = m_nodes(2, nodes(3))
        else
                dl => tree%downleft
                dr => tree%downright
                ul => tree%upleft
                ur => tree%upright
                
                tree%bbox(1,1) = u_inf
                tree%bbox(2,1) = u_inf
                tree%bbox(1,2) = u_ninf
                tree%bbox(2,2) = u_ninf
                
                if ( associated(dl) ) then
                        call bbox_quadtree(dl)
                        tree%bbox(1,1) = min(tree%bbox(1,1), dl%bbox(1,1))
                        tree%bbox(2,1) = min(tree%bbox(2,1), dl%bbox(2,1))
                        tree%bbox(1,2) = max(tree%bbox(1,2), dl%bbox(1,2))
                        tree%bbox(2,2) = max(tree%bbox(2,2), dl%bbox(2,2))
                end if
                if ( associated(dr) ) then
                        call bbox_quadtree(dr)
                        tree%bbox(1,1) = min(tree%bbox(1,1), dr%bbox(1,1))
                        tree%bbox(2,1) = min(tree%bbox(2,1), dr%bbox(2,1))
                        tree%bbox(1,2) = max(tree%bbox(1,2), dr%bbox(1,2))
                        tree%bbox(2,2) = max(tree%bbox(2,2), dr%bbox(2,2))
                end if
                if ( associated(ul) ) then
                        call bbox_quadtree(ul)
                        tree%bbox(1,1) = min(tree%bbox(1,1), ul%bbox(1,1))
                        tree%bbox(2,1) = min(tree%bbox(2,1), ul%bbox(2,1))
                        tree%bbox(1,2) = max(tree%bbox(1,2), ul%bbox(1,2))
                        tree%bbox(2,2) = max(tree%bbox(2,2), ul%bbox(2,2))
                end if
                if ( associated(ur) ) then
                        call bbox_quadtree(ur)
                        tree%bbox(1,1) = min(tree%bbox(1,1), ur%bbox(1,1))
                        tree%bbox(2,1) = min(tree%bbox(2,1), ur%bbox(2,1))
                        tree%bbox(1,2) = max(tree%bbox(1,2), ur%bbox(1,2))
                        tree%bbox(2,2) = max(tree%bbox(2,2), ur%bbox(2,2))
                end if
        end if
end subroutine

subroutine search_quadtree(point, elem)
        real*8, dimension(2), intent(in) :: point
        integer, intent(out) :: elem
        type(qtree), pointer :: tree
        logical*1, dimension(2) :: under
        tree => q_elemtree
        elem = 0
        do while( associated(tree) )
        if ( tree%elem /= 0 ) then
                if ( inside_element(tree%elem, point) ) then
                        elem = tree%elem
                end if
                exit
        end if
        under(1) = point(1) < tree%center(1)
        under(2) = point(2) < tree%center(2)
        if ( under(1) ) then
                if ( under(2) ) then
                        tree => tree%downleft
                else
                        tree => tree%upleft
                end if
        else
                if ( under(2) ) then
                        tree => tree%downright
                else
                        tree => tree%upright
                end if
        end if
        end do
end subroutine

recursive subroutine highest_recursive(xcoord, ycoord, elem, tree)
        real*8, intent(in)               :: xcoord
        real*8, intent(out)              :: ycoord
        integer, intent(inout)           :: elem
        logical*1                        :: under, down
        type(qtree), pointer, intent(in) :: tree

        if ( .not. associated(tree) ) then
                return
        end if
        if ( realover(xcoord, tree%bbox(1,2) )    .OR. &
             realunder(xcoord, tree%bbox(1,1) ) ) then
                return
        end if
        if ( tree%elem /= 0 ) then
                elem = tree%elem
                ycoord = tree%bbox(2,1)
                return
        end if
        if ( xcoord < tree%center(1) ) then
                call highest_recursive(xcoord, ycoord, elem, tree%downleft)
                if ( elem == 0 ) then
                        call highest_recursive(xcoord, ycoord, elem, tree%upleft)
                end if
        else
                call highest_recursive(xcoord, ycoord, elem, tree%downright)
                if ( elem == 0 ) then
                        call highest_recursive(xcoord, ycoord, elem, tree%upright)
                end if
        end if
end subroutine

subroutine highest_quadtree(xcoord, ycoord, elem)
        real*8, intent(in)               :: xcoord
        real*8, intent(out)              :: ycoord
        integer, intent(inout)           :: elem
        elem = 0
        call highest_recursive(xcoord, ycoord, elem, q_elemtree)
end subroutine


!subroutine highest_quadtree(xcoord, ycoord, elem)
!        real*8, intent(in)    :: xcoord
!        real*8, intent(out)   :: ycoord
!        integer, intent(out)  :: elem
!        logical*1             :: under, down
!        type(qtree), pointer  :: tree
!        integer, dimension(4) :: nodes
!        
!        if ( realequal(xcoord, 19400D0) ) then
!                if ( .not. realover(xcoord, q_elemtree%bbox(1,2) ) .and.&
!                     .not. realunder(xcoord, q_elemtree%bbox(1,1) ) ) then
!                        write(*,*) "Searching 19400", q_elemtree%bbox, c_rank 
!                end if
!                call wait_debug()
!        end if
!
!        tree => q_elemtree
!        elem = 0
!        do while( associated(tree) )
!        if ( realover(xcoord, tree%bbox(1,2) ) .OR.&
!             realunder(xcoord, tree%bbox(1,1) ) ) then
!                exit
!        end if
!        if ( tree%elem /= 0 ) then
!                elem = tree%elem
!                ycoord = tree%bbox(2,1)
!                exit
!        end if
!        if ( realequal(xcoord, tree%center(1)) ) then
!                write(*,*) "Search both", xcoord
!        end if
!                under = xcoord < tree%center(1)
!        if ( under ) then
!                down = .TRUE.
!                if ( .not. associated(tree%downleft) ) then
!                        down = .FALSE.
!                else if ( realover(xcoord, tree%downleft%bbox(1,2)) .OR.&
!                          realunder(xcoord, tree%downleft%bbox(1,1)) ) then
!                        down = .FALSE.
!                end if
!                if ( down ) then
!                        tree => tree%downleft
!                else
!                        tree => tree%upleft
!                end if
!        else
!                down = .TRUE.
!                if ( .not. associated(tree%downright) ) then
!                        down = .FALSE.
!                else if ( realover(xcoord, tree%downright%bbox(1,2)) .OR.&
!                          realunder(xcoord, tree%downright%bbox(1,1)) ) then
!                        down = .FALSE.
!                end if
!                if ( down ) then
!                        tree => tree%downright
!                else
!                        tree => tree%upright
!                end if
!        end if
!        end do
!end subroutine

recursive subroutine test_quadtree(nelems, tree, presnt)
        integer, intent(in) :: nelems
        type(qtree), pointer, intent(in) :: tree
        logical*1, dimension(:), allocatable, intent(inout) :: presnt
        
        if ( .not. allocated(presnt) ) then
                allocate(presnt(nelems))
                presnt = .FALSE.
        end if
        
        if ( associated(tree) ) then
                if ( tree%elem /= 0 ) then
                        if ( tree%elem <= nelems .AND. tree%elem > 0 ) then
                                if ( presnt(tree%elem) ) then
                                        write(*,*) "repeated element", tree%elem
                                else
                                        presnt(tree%elem) = .TRUE.
                                end if
                        else
                                write(*,*) "elem outside rank",tree%elem,nelems
                        end if
                end if
                call test_quadtree(nelems, tree%downleft, presnt)
                call test_quadtree(nelems, tree%upleft, presnt)
                call test_quadtree(nelems, tree%downright, presnt)
                call test_quadtree(nelems, tree%upright, presnt)
        end if
end subroutine

recursive subroutine delete_quadtree(tree)
        type(qtree), pointer, intent(inout) :: tree
        if ( associated(tree) ) then
                call delete_quadtree(tree%downleft)
                call delete_quadtree(tree%upleft)
                call delete_quadtree(tree%downright)
                call delete_quadtree(tree%upright)
                deallocate(tree)
                nullify(tree)
        end if
end subroutine

end module quadtree
