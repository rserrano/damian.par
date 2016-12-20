module model
use utilities
use formulation
use pwaves
use io_utils
implicit none
! Kinds of boundary conditions for nodes.
integer, parameter :: BC_NEUMANN   = 2, BC_DIRICLET = 1, &
                      BC_ABSORBING = 3, BC_SEISMIC  = 4
! Actual boundary conditions. Add a definition here.
! Newmann boundary conditions for nodes.
integer, parameter :: RC_ZRICKER = 1
! Diriclet boundary conditions for nodes.
integer, parameter :: DC_XFIXED = 1, DC_ZFIXED = 2,&
                      DC_FIXED = 3, DC_ZRICKER = 4

integer, parameter :: AB_NOABS = 0,  AB_ZDIR    = 1, AB_XDIR    = 2

integer, parameter :: SE_SVERL = 2, SE_SVERR = 4, SE_SHORI = 1, &
                      SE_CORNL = 6, SE_CORNR = 5
! Global information
integer :: m_nummat
integer :: m_lastgnode
integer :: m_totnodes
integer :: m_totelems
real*8  :: m_duration
real*8  :: m_deltat
real*8  :: m_blength
real*8  :: m_mlength
integer :: m_nsteps
integer :: m_totrows
integer :: m_totcols

! Constants of the material (mu, lambda, rho, vp, vs) x Num materials.
real*8, allocatable, dimension(:,:)   :: m_matconst
! Given a material and the characteristic frequency, constants to solve the
! problem for that material. (characteristic length, characteristic time step).
real*8, allocatable, dimension(:,:)   :: m_probconst

! Local information.
integer :: m_nnodes
integer :: m_nelems
integer :: m_nrows
integer :: m_ncols

! Stiffness matrix of the material. 8 x 8 x Num materials.
real*8, allocatable, dimension(:,:,:) :: m_kmat
real*8, dimension(8,8)                :: m_mmat

! Number of material property of the elements. Num elements.
integer, allocatable, dimension(:)    :: m_ematprop

! Location of the points in the grid 2 x Num Nodes.
real*8, allocatable, dimension(:,:)   :: m_nodes

! CSR format elements for distributed use.
integer, allocatable, dimension(:)    :: m_elmdist
integer, allocatable, dimension(:)    :: m_eptr
integer, allocatable, dimension(:)    :: m_eind
logical*1, allocatable, dimension(:)  :: m_ediv
integer, allocatable, dimension(:)    :: m_gnode
integer, allocatable, dimension(:,:)  :: m_hanging
! Absorbing boundaries.
integer, allocatable, dimension(:)    :: m_absorbing
! Diriclet and neumann boundary conditions.
integer, allocatable, dimension(:,:)  :: m_ndiriclet

! Parameters defining the geometry to tacle.
character(len=6) :: m_geom
real*8, allocatable, dimension(:)     :: m_gparams

! Parameters defining the boundary conditions to tacle.
character(len=6) :: m_prob
real*8, allocatable, dimension(:) :: m_pparams

! Constants specific to the problems

! Parameters of the incident plane wave.
real*8, dimension(2) :: m_worigin
real*8               :: m_wttoorigin
real*8               :: m_wangle
real*8               :: m_wamp
real*8               :: m_wfreq
real*8, dimension(2) :: m_wvelo
real*8, dimension(2) :: m_wlame
integer              :: m_wzero
integer              :: m_wrows
integer              :: m_wcols

integer, allocatable, dimension(:,:) :: m_seismic

! Parameters of the force or displacement at a point.
real*8 :: m_pointamp
real*8 :: m_pointsigma
real*8 :: m_pointstart

integer, allocatable, dimension(:) :: m_pointnodes


! Interval for saving the problem.
real*8 :: m_inttime

contains

! UTILITIES:

! FUNCTION
! elength:      returns the element length.
! PARAMETERS
! ele(in):      the intex of the element in m_eptr
! return:       the length of the element
function elength(ele)
        integer, intent(in)   :: ele
        real*8 :: elength
        integer, dimension(4) :: nodes
        call elemallnhnodes(ele, nodes)
        elength = (m_nodes(2,nodes(3)) - m_nodes(2,nodes(1)))
end function

! SUBROUTINE
! elemallnhdofs:  returns the DOFs of the non hanging 
!                 nodes of the element.
! PARAMETERS
! ele:          the local index of the element.
! dofs:         the dofs for the program.
subroutine elemallnhdofs(ele, dofs)
        integer, intent(in) :: ele
        integer, dimension(8), intent(out) :: dofs
        integer, dimension(4) :: nodes
        call elemallnhnodes(ele, nodes)
        dofs(1:7:2) = nodes*2-1
        dofs(2:8:2) = dofs(1:7:2)+1
end subroutine

! SUBROUTINE
! elementnodes: returns the number of nodes starting from until.
! PARAMETERS
! ele(in):      the index of the element in m_eptr
! st(in):       the node to start.
! n(in):        number of nodes.
! nodes(out):   nodes.

subroutine elementnodes(ele, st, n, nodes)
        integer, intent(in) :: ele, st, n
        integer, dimension(n), intent(out) :: nodes
        integer :: ind, i
        ind = m_eptr(ele)+st-1
        do i = 1, n
                if ( ind >= m_eptr(ele+1) ) then
                        ind = m_eptr(ele)
                end if
                nodes(i) = m_eind(ind)
                ind = ind + 1
        end do
end subroutine

! SUBROUTINE
! elemallnhnodes: returns the nodes that are not hanging nodes. (4)
! PARAMETERS
! ele(in):      the index of the element in m_eptr
! nodes(out):   nodes.

subroutine elemallnhnodes(ele, nodes)
        integer, intent(in)                :: ele
        integer, dimension(4), intent(out) :: nodes
        if ( count(.not. m_ediv(m_eptr(ele):m_eptr(ele+1)-1) ) /= 4 ) then
                write(*,*) "Error at: ", __FILE__, ":", __LINE__, ele,&
                           count(.not. m_ediv(m_eptr(ele):m_eptr(ele+1)-1)),&
                           count(m_ediv(m_eptr(ele):m_eptr(ele+1)-1)),&
                           m_eptr(ele+1)-m_eptr(ele)
        end if
        nodes = pack(m_eind(m_eptr(ele):m_eptr(ele+1)-1),&
                     .not. m_ediv(m_eptr(ele):m_eptr(ele+1)-1))
end subroutine

function inside_element(ele, point)
        logical*1 :: inside_element
        real*8, dimension(2), intent(in) :: point
        integer, intent(in) :: ele
        integer, dimension(4) :: nodes
        call elemallnhnodes(ele, nodes)
        inside_element = .NOT. realunder(point(1), m_nodes(1,nodes(1))) .AND. &
                         .NOT. realunder(point(2), m_nodes(2,nodes(1))) .AND. &
                         .NOT. realover(point(1),  m_nodes(1,nodes(3))) .AND. &
                         .NOT. realover(point(2),  m_nodes(2,nodes(3)))
end function


! FUNCTION 
! inside_wave:         returns true if a point is inside the place where
!                       the plane wave loads are applied.
! PARAMETERS
! pt(in):               point to be evaluated.
! return:               boolean true if inside the domain. False in the oposite
!                       case.

function inside_wave(pt)
        logical*1 :: inside_wave
        real*8, dimension(2) :: pt
        real*8 :: left, right, top
        if ( m_prob(1:4) /= 'WAVE' ) then
                inside_wave = .FALSE.
                return
        end if
        left  = real(m_wcols-1,kind=8)*m_blength
        right = real((m_totcols-m_wcols)-1,kind=8)*m_blength
        top   = real(m_wrows-1,kind=8)*m_blength
        if ( realover(pt(1), left)   .AND. &
             realunder(pt(1), right) .AND. &
             realunder(pt(2), top) ) then
                inside_wave = .TRUE.
        else
                inside_wave = .FALSE.
        end if
end function

function inside_soil(pt)
        logical*1 :: inside_soil
        real*8, dimension(2), intent(in) :: pt
        real*8 :: left, right, top
        if ( m_prob(1:4) /= 'WAVE' ) then
                inside_soil = .FALSE.
                return
        end if
        left  = real(m_wcols-1,kind=8)*m_blength
        right = real((m_totcols-m_wcols)-1,kind=8)*m_blength
        top   = real(m_wrows-1,kind=8)*m_blength
        if ( realunder(pt(1), left) .OR. &
             realover(pt(1), right) .OR. &
             realover(pt(2), top) ) then
                inside_soil = .TRUE.
        else
                inside_soil = .FALSE.
        end if
end function

function void_soil(pt)
        logical*1 :: void_soil
        real*8, dimension(2), intent(in) :: pt
        void_soil = realunder(pt(2), m_worigin(2))
end function

! FUNCTION
! needs_division:       Returns true if the geometry and wave velocities of 
!                       the elements dictate that it needs to be divided 
!                       further.
! 
! PARAMETERS:
! ele(in)               element to be evaluated.
! return:               boolean true if the element has to be divided further.

function needs_division(ele)
        logical*1 :: needs_division
        integer, intent(in)   :: ele
        integer :: i, mat
        integer, dimension(4) :: nodes
        real*8  :: leng, elen
        real*8, dimension(2)  :: center, pt
        elen = elength(ele)
        if ( realover(elen, m_mlength) ) then
                needs_division = .FALSE.
                return
        end if
        call elemallnhnodes(ele, nodes)
        center = (/sum(m_nodes(1, nodes)), sum(m_nodes(2, nodes))/)/4.0D0
        mat = material_prop(center)
        leng = m_probconst(1,mat)
        if ( leng < elen ) then
                needs_division = .TRUE.
                return
        end if
        if ( inside_soil(center) ) then
                needs_division = .FALSE.
                return
        end if
        leng = u_inf
        do i = -1,1,2
                pt = center
                pt(1) = pt(1) + real(i, kind=8)*elen
                mat = material_prop(pt)
                leng = min(leng, m_probconst(1,mat))
                pt = center
                pt(2) = pt(2) + real(i, kind=8)*elen
                mat = material_prop(pt)
                leng = min(leng, m_probconst(1,mat))
        end do
        if ( (leng*2.0D0) < elen ) then
                needs_division = .TRUE.
        else
                needs_division = .FALSE.
        end if
end function


! FUNCTION
! needs_deletion:       Returns true if the geometry dictates the element has to be
!                       cut from the model.
! 
! PARAMETERS:
! ele(in)               element to be evaluated.
! return:               boolean true if the element has to be divided further.

function needs_deletion(ele)
        logical*1 :: needs_deletion
        integer, intent(in)    :: ele
        integer, dimension(4)  :: nodes
        real*8, dimension(2)   :: point
        real*8, dimension(3)   :: trpt
        integer :: cvmprop, i
        real*8 :: elen
        call elemallnhnodes(ele, nodes)
        elen = elength(ele)
                if ( realover(elen, m_mlength) ) then
                needs_deletion = .FALSE.
                return
        end if
        point = (/sum(m_nodes(1, nodes)), sum(m_nodes(2, nodes))/)/4.0D0
        if ( inside_soil(point) ) then
                if ( void_soil(point) ) then
                        needs_deletion=.TRUE.
                else
                        needs_deletion=.FALSE.
                end if
                return
        end if

        needs_deletion = .FALSE.
        select case(m_geom)
        case('CANYON')
                needs_deletion = .TRUE.
                do i = 1, 4 
                point = m_nodes(:,nodes(i))
                if ( distance(point, m_gparams(1:2)) > m_gparams(3) ) then
                        needs_deletion = .FALSE.
                end if
                end do
        case('RECTAN')
                needs_deletion = .TRUE.
                do i = 1, 4 
                point = m_nodes(:,nodes(i))
                if ( abs(point(1)-m_gparams(1)) > m_gparams(3) .OR. &
                     abs(point(2)-m_gparams(2)) > m_gparams(4) ) then
                        needs_deletion = .FALSE.
                end if
                end do
        case('CVMPLX')
                needs_deletion = .TRUE.
                do i = 1, 4 
                point = m_nodes(:,nodes(i))
                trpt(1) = m_gparams(1)
                trpt(2) = m_gparams(2)+point(1)
                trpt(3) = m_gparams(3)-point(2)
                call cvm_mat(trpt, cvmprop)
                if ( (cvmprop + 1) /= 18 ) then
                        needs_deletion = .FALSE.
                        exit
                end if
                end do
        case('CVMPLY')
                needs_deletion = .TRUE.
                do i = 1, 4 
                point = m_nodes(:,nodes(i))
                trpt(1) = m_gparams(1)+point(1)
                trpt(2) = m_gparams(2)
                trpt(3) = m_gparams(3)-point(2)
                call cvm_mat(trpt, cvmprop)
                if ( (cvmprop + 1) /= 18 ) then
                        needs_deletion = .FALSE.
                        exit
                end if
                end do
        case('CVMPLZ')
                needs_deletion = .TRUE.
                do i = 1, 4 
                point = m_nodes(:,nodes(i))
                trpt(1) = m_gparams(1)+point(1)
                trpt(2) = m_gparams(2)+point(2)
                trpt(3) = m_gparams(3)
                call cvm_mat(trpt, cvmprop)
                if ( (cvmprop + 1) /= 18 ) then
                        needs_deletion = .FALSE.
                        exit
                end if
                end do
        case default
                needs_deletion = .FALSE.
        end select
end function

! FUNCTION
! material_prop
! 
! PARAMETERS:
! pt(in)                point of the material property.
! return:               integer with the material property given the point.

function material_prop(pt)
        integer :: material_prop
        real*8, dimension(2), intent(in) :: pt
        real*8, dimension(3) :: trpt
        integer :: cvmprop
        material_prop = 1
        if ( inside_soil(pt) ) then
                material_prop = int(m_pparams(5))
                return
        end if 
        select case(m_geom)
        case('VALLEY')
                if ( distance(pt, m_gparams(1:2)) < m_gparams(3) ) then
                        material_prop = 2
                else
                        material_prop = 1
                end if
        case('VARECT')
                if ( abs(pt(1)-m_gparams(1)) < m_gparams(3) .AND. &
                     abs(pt(2)-m_gparams(2)) < m_gparams(4) ) then
                        material_prop = 2
                else
                        material_prop = 1
                end if
        case('CVMPLX')
                trpt(1) = m_gparams(1)
                trpt(2) = m_gparams(2)+pt(1)
                trpt(3) = m_gparams(3)-pt(2)
                call cvm_mat(trpt, cvmprop)
                material_prop = cvmprop + 1
        case('CVMPLY')
                trpt(1) = m_gparams(1)+pt(1)
                trpt(2) = m_gparams(2)
                trpt(3) = m_gparams(3)-pt(2)
                call cvm_mat(trpt, cvmprop)
                material_prop = cvmprop + 1
        case('CVMPLZ')
                trpt(1) = m_gparams(1)+pt(1)
                trpt(2) = m_gparams(2)+pt(2)
                trpt(3) = m_gparams(3)
                call cvm_mat(trpt, cvmprop)
                material_prop = cvmprop + 1
        case('CVMPLN')
                ! General planes functionality, not implemented.
                material_prop = 1
        case default
                material_prop = 1
        end select
end function

! FUNCTION
! hanging_element:      returns true if the element has a hanging node.
! PARAMETERS
! ele:                  the element to be evaluated.
! return:               .TRUE. if the element has a hanging node. .FALSE.
!                       otherwise.

function hanging_element(ele)
        integer, intent(in) :: ele
        logical*1 :: hanging_element
        hanging_element = ((m_eptr(ele+1)-m_eptr(ele)) > 4 )
end function

! SUBROUTINE
! local_element:        returns true if a point is inside the finite where 
!                       loads are applied.
! PARAMETERS
! ele(in):              global element index.
! local(out):           local index in the processor.
! proc(out):            processor.
! TODO:                 make the method faster with binary search. (Is it really
!                       necessary?)

subroutine local_element(ele, local, proc)
        integer, intent(in)  :: ele
        integer, intent(out) :: local, proc
        integer :: i
        do i = 1, size(m_elmdist)-1
                if ( ele >= m_elmdist(i) .AND.   &
                     ele <  m_elmdist(i+1) ) then
                        local = (ele-m_elmdist(i))+1
                        proc  = (i-1)
                        exit
                end if
        end do
end subroutine

! SUBROUTINE
! edge_connects:        returns the index of the edge that connects both
!                       elemens. The edge is between 1 and the number of nodes
!                       that the element has.
! PARAMETERS:
! elea(in):             The nodes in element a.
! eleb(in):             The nodes in element b.
! inda(out):            The index of the edge in element a.
! indb(out):            The index of the edge in element b.

subroutine edge_connects(elea, eleb, inda, indb)
        integer, dimension(:), intent(in) :: elea, eleb
        integer, intent(out) :: inda, indb
        integer :: nxta, nxtb
        do inda = 1, size(elea)
                do indb = 1, size(eleb)
                        nxta = inda + 1
                        if ( nxta > size(elea) ) then
                                nxta = 1
                        end if
                        nxtb = indb + 1
                        if ( nxtb > size(eleb) ) then
                                nxtb = 1
                        end if
                        if ( elea(inda) == eleb(nxtb) .AND.  &
                             elea(nxta) == eleb(indb) ) then
                             return
                        end if
                end do
        end do
        inda = 0 
        indb = 0
end subroutine

! SUBROUTINE
! nodes_shared:        returns the nodes shared by both elements
! PARAMETERS:
! elea(in):             The nodes in element a.
! eleb(in):             The nodes in element b.
! nodes:                The nodes shared in the positions of element a.

subroutine nodes_shared(elea, eleb, nodes)
        integer, dimension(:), intent(in)  :: elea, eleb
        integer, dimension(:), intent(out) :: nodes
        integer :: i, j
        nodes = 0
        do i = 1, size(elea)
                do j = 1, size(eleb)
                        if ( elea(i) == eleb(j) ) then
                                nodes(i) = elea(i)
                                exit
                        end if
                end do
        end do
end subroutine


! SUBROUTINE
! global_element:       Returns true if a point is inside the finite where 
!                       loads are applied.
! PARAMETERS
! ele(in):              local element index.
! proc(in):             processor.
! global(out):          global index of the node.


subroutine global_element(ele, proc, global)
        integer, intent(in) :: ele, proc
        integer, intent(out) :: global
        global = (m_elmdist(proc+1)+ele-1)
end subroutine


! SUBROUTINE           
! getmassof:            Returns the mass matrix for the current elem lumped.
!
! PARAMETERS
! ele(in):              local element index.
! mass(out):            the mass matrix.

subroutine getmassof(ele, mass)
        integer, intent(in) :: ele
        real*8, intent(out) :: mass
#ifdef DEBUG
        if ( ele > m_nelems .OR. ele <= 0 ) then
                write(*,*) "Error: ", __FILE__, ":", __LINE__
        end if
        if ( m_ematprop(ele) > m_nummat .OR. m_ematprop(ele) <= 0 ) then
                write(*,*) "Error: ", __FILE__, ":", __LINE__ 
        end if 
#endif
        call calcm(elength(ele)/2.0D0,&
                   m_matconst(3,m_ematprop(ele)), mass)
end subroutine

! SUBROUTINE
! getdampingof:         Returns the value for the damping matrix lumped.
!
! PARAMETERS
! ele(in):              local element index.
! damping(out):         the damping lumped matrix.
subroutine getdampingof(ele, damping)
        integer, intent(in)               :: ele
        integer, dimension(4)             :: nodes
        real*8, dimension(8), intent(out) :: damping
        integer :: i, epr
        epr = m_ematprop(ele)
        damping(:) = 0.0D0
        call elemallnhnodes(ele, nodes)
        do i = 1, 4
                if (m_absorbing(nodes(i)) == AB_ZDIR) then
                        call lysmer(elength(ele),m_matconst(3,epr),&
                                    (/m_matconst(4,epr), m_matconst(5,epr)/),&
                                    logical(.TRUE., kind=1), damping(2*i-1:2*i))
                else if ( m_absorbing(nodes(i)) == AB_XDIR ) then
                        call lysmer(elength(ele),m_matconst(3,epr),&
                                    (/m_matconst(4,epr), m_matconst(5,epr)/),&
                                    logical(.FALSE., kind=1), damping(2*i-1:2*i))
                endif
        end do
end subroutine

! SUBROUTINE
! set_dispshanging:     Sets the displacements for the hanging nodes.
!
! PARAMETERS
! ele(in):              local element index.
! damping(out):         the damping lumped matrix.
subroutine set_dispshanging(u)
        real*8, dimension(:), intent(inout) :: u
        integer, dimension(6) :: dofs
        integer :: i
        if ( .NOT. allocated(m_hanging) ) then
                return
        end if
        do i = 1, size(m_hanging, 2)
        dofs(2:6:2) = 2*m_hanging(:,i)
        dofs(1:5:2) = dofs(2:6:2) - 1
        u(dofs(1:2)) = (u(dofs(3:4)) + u(dofs(5:6)))/2.0D0
        end do
end subroutine

! SUBROUTINE
! sum_reactanchors:     Sets the displacements for the hanging nodes.
!
! PARAMETERS
! ele(in):              local element index.
! damping(out):         the damping lumped matrix.
subroutine sum_reactanchors(u)
        real*8, dimension(:), intent(inout) :: u
        integer, dimension(6) :: dofs
        integer :: i
        real*8, dimension(2) :: avg
        if ( .NOT. allocated(m_hanging) ) then
                return
        end if
        do i = 1, size(m_hanging, 2)
        dofs(2:6:2) = 2*m_hanging(:,i)
        dofs(1:5:2) = dofs(2:6:2) - 1
        red = u(dofs(1:2))-(u(dofs(3:4)) + u(dofs(5:6)))/2.0D0
        u(dofs(3:4)) = u(dofs(3:4)) + red
        u(dofs(5:6)) = u(dofs(5:6)) + red
        end do
end subroutine


! SUBROUTINE
! setup_hanging:        Finds which are the hanging nodes in this processor.
!                       and initializes the model so they are considered by
!                       set_dispshanging.
subroutine setup_hanging()
        integer :: i, j, ind, n
        n = count(m_ediv)
        if ( n > 0 ) then
                allocate(m_hanging(3,n))
                ind = 0
                do i = 1, size(m_eptr)-1
                do j = m_eptr(i), m_eptr(i+1)-1
                if ( m_ediv(j) ) then
                        ind = ind + 1
                        m_hanging(1,ind) = m_eind(j)
                        m_hanging(2,ind) = m_eind(j-1)
                        if ( (j+1) < m_eptr(i+1) ) then
                                m_hanging(3,ind) = m_eind(j+1)
                        else
                                m_hanging(3,ind) = m_eind(m_eptr(i))
                        end if
                end if
                end do
                end do
        end if
end subroutine


! Subroutines specific to problems. Depending on the problem being simulated,
! they get called. They setup and execute boundary conditions.

! SUBROUTINE
! summed_reactions:    Sets the reactions to the problem depending on the
!                      parameters read from the problem file.
!
! PARAMETERS
! time(in):            The time for the reactions to be calculated.
! reactions(out):      Sums the reactions to the vector.

subroutine additive_reactions(time, reactions)
        real*8, intent(in)                  :: time
        real*8, dimension(:), intent(inout) :: reactions
        select case(m_prob)
        case('WAVEP ')
                call wave_reactions(time, reactions)
        case('WAVESV')
                call wave_reactions(time, reactions)
        case('WAVESH')
                ! TODO: Implement SH incomming wave.
        case('DISLOC')
                ! TODO: Implement dislocation.
        case default

        end select
end subroutine

subroutine absolute_reactions(time, reactions)
        real*8, intent(in)                  :: time
        real*8, dimension(:), intent(inout) :: reactions
        select case(m_prob)
        case('POINTN')
                call point_reactions(time, reactions)
        case default

        end select
end subroutine
 

        
! SUBROUTINE
! problem_displacements:        Sets the displacementss to the problem
!                               depending on the parameters read from the 
!                               problem file.
!
! PARAMETERS
! time(in):                     The time for the displacements to be calculated.
! displacements(out):           Sums the displacements to the vector.

subroutine problem_displacements(time, displacements)
        real*8, intent(in)                  :: time
        real*8, dimension(:), intent(inout) :: displacements
        select case(m_prob)
        case('POINTD')
                call point_displacements(time, displacements)
        case default
        end select
end subroutine

! SUBROUTINE
! wave_reactions:       At each iteration it sets the plane wave, P and SV
!                       reactions vector.
! PARAMETERS
! time(in)
! reactions(out)

subroutine wave_reactions(time, reactions)
        real*8, intent(in) :: time
        real*8, dimension(:), intent(inout) :: reactions
        real*8, dimension(8)    :: u, a, p
        real*8, dimension(2,3)  :: up
        integer, dimension(8)   :: edofs
        integer, dimension(4)   :: nodes
        logical*1, dimension(8) :: e
        integer :: ele, bck
        real*8, dimension(2)    :: outp
        real*8 :: outt, mass
        integer :: i, j

        if ( .NOT. allocated(m_seismic) ) then
                return
        end if
        call timetransform(m_wttoorigin, time, outt)
        
        do i = 1, size(m_seismic,2)
                ele = m_seismic(1,i)
                bck = m_seismic(2,i)
                call elemallnhnodes(ele, nodes)
                call elemallnhdofs(ele, edofs)
                call getmassof(ele, mass)
                do j = 1,4
                        call pointtransform(m_worigin,           &
                                            m_nodes(:,nodes(j)), &
                                            outp)
                        call planewave(m_prob,outt,outp(1),     &
                                       outp(2),up)
                        u(2*j-1:2*j) = up(1:2,1)
                        a(2*j-1:2*j) = up(1:2,3)
                end do
                if ( bck == SE_SVERL ) then
                        e = (/.TRUE., .TRUE., .FALSE., .FALSE.,&
                              .FALSE., .FALSE., .TRUE., .TRUE./)
                else if ( bck == SE_SVERR ) then
                        e = (/.FALSE., .FALSE., .TRUE., .TRUE.,&
                              .TRUE., .TRUE., .FALSE., .FALSE./)
                else if (bck == SE_SHORI ) then
                        e = (/.FALSE., .FALSE., .FALSE., .FALSE.,&
                              .TRUE., .TRUE., .TRUE., .TRUE./)
                else if (bck == SE_CORNL ) then        
                        e = (/.TRUE., .TRUE., .FALSE., .FALSE.,&
                              .TRUE., .TRUE., .TRUE., .TRUE.  /)
                else if (bck == SE_CORNR ) then
                        e = (/.FALSE., .FALSE., .TRUE., .TRUE.,&
                              .TRUE., .TRUE., .TRUE., .TRUE.  /)
                end if
                call multdrm(m_kmat(:,:,m_ematprop(ele)), u, e, p)
                reactions(edofs) = reactions(edofs)+p
                call multdrm(m_mmat, a, e, p)
                p = p*mass
                reactions(edofs) = reactions(edofs)+p
        end do
end subroutine

! SUBROUTINE
! point_reactions:      At each iteration, it sets the singular point reactions
!                       vector.
! PARAMETERS
! time(in):             The time for the reactions to be calculated.
! reactions(out):       Sums the reactions to the vector.

subroutine point_reactions(time, reactions)
        real*8, intent(in) :: time
        real*8, dimension(:), intent(inout) :: reactions
        real*8 :: val
        integer :: i, dof
        if ( .NOT. allocated(m_pointnodes) ) then
                return
        end if
        do i = 1, size(m_pointnodes)
                dof = 2*m_pointnodes(i)-1
                call ricker(m_pointamp , m_pointsigma, time-m_pointstart, val)
                reactions(dof+1) = reactions(dof+1)+val
        end do
end subroutine

! SUBROUTINE
! point_displacements:  At each iteration, it sets the singular point
!                       displacements vector.
! PARAMETERS
! time(in):             The time for the displacements to be calculated.
! reactions(out):       Sums the displacements to the vector.

subroutine point_displacements(time, displacements)
        real*8, intent(in) :: time
        real*8, dimension(:), intent(inout) :: displacements
        integer :: i, dof
        if ( .NOT. allocated(m_pointnodes) ) then
                return
        end if
        do i = 1, size(m_pointnodes)
                dof = 2*m_pointnodes(i)-1
                call ricker(m_pointamp, m_pointsigma, time-m_pointstart, &
                            displacements(dof+1))
        end do
end subroutine

! SUBROUTINE
! setup_problem_params: Starts the problem by setting the global parameters of
!                       it. 
!
! PARAMETERS
! height(in):           The height of the problem.
! width(in):            The width of the problem.
! duration(in):         Whole duration of the problem, (not reorganized with
!                       time step.
! chper(in):            Characteristic period. (1/freq).

subroutine setup_problem_params(height, width, duration, chper)
        real*8, intent(in)   :: height, width, duration, chper
        integer :: n, m, steps, i
        real*8  :: length, deltat
        allocate(m_probconst(2,m_nummat))
        length = 0.0D0
        deltat = u_inf
        do i = 1, m_nummat
                call explicitparams(chper,m_matconst(4,i),m_matconst(5,i),&
                                    m_probconst(1,i),m_probconst(2,i))
                if ( length < m_probconst(1,i) ) then
                        length = m_probconst(1,i)
                end if
                if ( deltat > m_probconst(2,i) ) then
                        deltat = m_probconst(2,i)
                end if
        end do
        ! It's more important to respect the vertical parameter than the.
        ! horizontal.
        m = ceiling(height/length)
        length = height/real(m, 8)
        n = ceiling(width/length)
        steps = ceiling(duration/deltat)
        deltat = duration/real(steps, 8)
        m_duration = deltat*steps
        m_blength = length
        m_mlength = m_blength
        m_totcols = n+1
        m_totrows = m+1
        m_nsteps = steps
        m_deltat = deltat
        m_lastgnode = m_totrows*m_totcols
        m_totelems = (m_totrows-1)*(m_totcols-1)
end subroutine


! Sets up the materials and their stiffness matrices.
subroutine setup_materials(nmats, cmats, vmats)
        integer, intent(inout) :: nmats
        character(len=1), dimension(:,:), intent(inout) :: cmats
        real*8, dimension(:,:), intent(inout) :: vmats
        integer, dimension(:,:), allocatable  :: imats
        integer :: i
        m_nummat = nmats
        allocate(m_matconst(5,m_nummat))
        allocate(imats(2,m_nummat))
        do i = 1, m_nummat
        imats(1,i) = c_to_p(cmats(1,i))
        imats(2,i) = c_to_p(cmats(2,i))
        if ( any(imats(:,i) == 0) ) then
                if ( imats(2,i) == f_nu ) then
                        imats((/2,1/),i) = imats((/1,2/),i)
                        cmats((/2,1/),i) = cmats((/1,2/),i)
                        vmats((/2,1/),i) = vmats((/1,2/),i)
                end if
                if ( imats(1,i) == f_nu ) then
                        imats(2,i) = c_to_v(cmats(2,i))
                        if ( imats(2,i) == 0 ) then
                        call error_exfem("Not recognized material")
                        else
                        call v_nu_to_v(imats(2,i), vmats(2,i),      &
                                       vmats(1,i), m_matconst(4:5,i))
                        m_matconst(3,i)   = vmats(3,i)
                        call v_to_e(m_matconst(4:5,i), m_matconst(3,i),    &
                                   (/f_shear, f_lambda/), m_matconst(1:2,i))
                        end if
                else if ( all(imats(:,i) == 0 ) ) then
                        imats(1,i) = c_to_v(cmats(1,i))
                        imats(2,i) = c_to_v(cmats(2,i))
                        if ( any(imats(:,i) == 0 ) ) then
                        call error_exfem("Not recognized material")
                        end if
                        m_matconst(3+imats(1,i),i) = vmats(1,i)
                        m_matconst(3+imats(2,i),i) = vmats(2,i)
                        m_matconst(3,i) = vmats(3,i)
                        call v_to_e(m_matconst(4:5,i), m_matconst(3,i),     &
                                   (/f_shear, f_lambda/), m_matconst(1:2,i))
                else
                        call error_exfem("Not recognized material")
                end if
        else
                call e_to_e(imats(:,i), (/f_shear, f_lambda/),  &
                            vmats(1:2,i), m_matconst(1:2,i))
                m_matconst(3,i)   = vmats(3,i)
                call e_to_v(m_matconst(1:2,i), m_matconst(3,i),    &
                           (/f_shear, f_lambda/), m_matconst(4:5,i))
        end if
        end do
        deallocate(imats)
end subroutine

! SUBROUTINE
! problem_setup:        This should be executed after the mesh
!                       is completely created. It depends on
!                       geometric information instead of node 
!                       number information, so it should work in
!                       every region regardless.
!                       If first set ups the shared parameters of the problem
!                       that depend on geometric information, absorbing 
!                       boundaries and hanging nodes.
!                       The problem that gets setted up depends on m_prob

subroutine problem_setup()
        ! The following have to be called regardless of the problem.
        call calcsqm(m_mmat) 
        call init_stiffness()
        call setup_hanging()
        call create_absorbingbc()
        select case(m_prob)
        case('WAVEP ')
                call setup_wave()
        case('WAVESV')
                call setup_wave()
                ! TODO: Implement SV incoming wave.
        case('WAVESH')
                ! TODO: Implement SH incomming wave.
        case('DISLOC')
                ! TODO: Implement dislocation.
        case('POINTN')
                call setup_point_source()
        case('POINTD')
                call setup_point_source()
        case default
        end select
end subroutine

subroutine init_wave()
        ! For a PWave the vector m_pparams is:                 
        ! (/ width of wave (aprox), height of wave, (aprox),   
        !    X coordinate of origin, Y coordinate of origin,   
        !    index of material const, time to origin,          
        !    angle, amplitude, frequency of ricker pulse /)
        m_wrows = ceiling(m_pparams(2)/m_blength)+1
        m_wcols = (m_totcols-(1+ceiling(m_pparams(1)/m_blength)))/2
        m_worigin = m_pparams(3:4)
        m_wvelo(1:2) = m_matconst(4:5,int(m_pparams(5)))
        m_wlame(1:2) = m_matconst(1:2,int(m_pparams(5)))
        m_wttoorigin = m_pparams(6)
        m_wangle     = m_pparams(7)
        m_wamp       = m_pparams(8)
        m_wfreq      = m_pparams(9)
        m_wzero      = floor(m_worigin(2)/m_blength)
        m_worigin(2) = real(m_wzero,kind=8)*m_blength
        call initplanewaves(m_prob, m_wamp, m_wvelo(1), m_wvelo(2), &
                            m_wangle, m_wfreq)
end subroutine


! SUBROUTINE:   Sets up the wave.
!
subroutine setup_wave()
        call wave_elems()
end subroutine

! SUBPROCEDURE
! wave_elems:          Selects the elements that will contribute to the plane
!                       wave reactions. Sets them up with their corresponding
!                       function in the wave.

subroutine wave_elems()
        integer :: ind
        real*8 :: left, right, top
        integer, dimension(4) :: nodes
        type(linked_list_2d), pointer :: lselems
        nullify(lselems)
        left  = real(m_wcols-1,kind=8)*m_blength
        right = real((m_totcols-m_wcols)-1,kind=8)*m_blength
        top   = real(m_wrows-1,kind=8)*m_blength
        do ind = 1, (size(m_eptr,1)-1)
                call elemallnhnodes(ind, nodes)
                if ( realequal(m_nodes(1,nodes(2)), left) .AND. &
                     realequal(m_nodes(2,nodes(2)), top) ) then
                        call insertll_2d(lselems, (/ind, SE_CORNL/))
                else if ( realequal(m_nodes(1,nodes(1)), right) .AND. &
                          realequal(m_nodes(2,nodes(1)), top) ) then
                        call insertll_2d(lselems, (/ind, SE_CORNR/))
                else if ( realequal(m_nodes(1,nodes(2)), left) .AND.  &
                          realunder(m_nodes(2,nodes(2)), top) ) then
                        call insertll_2d(lselems, (/ind, SE_SVERL/))
                else if ( realequal(m_nodes(1,nodes(1)), right) .AND. &
                          realunder(m_nodes(2,nodes(1)), top) ) then
                        call insertll_2d(lselems, (/ind, SE_SVERR/))
                else if ( realover(m_nodes(1,nodes(2)), left)   .AND. &
                          realunder(m_nodes(1,nodes(1)), right) .AND. &
                          realequal(m_nodes(2,nodes(1)), top) ) then
                        call insertll_2d(lselems, (/ind, SE_SHORI/))
                end if
        end do
        call toarray(lselems, m_seismic)
end subroutine

! SUBPROCEDURE
! setup_point_source:           Sets up the nodes for the punctual source
!                               and the parameters like amplitude, sigma
!                               
subroutine setup_point_source()
        integer :: i, j
        real*8  :: dist
        real*8, dimension(2) :: pt
        type(linked_list_1d), pointer :: lanodes
        character(len=40) :: hn
        integer :: pid, stat
        ! For a punctual source, the vector m_pparams is:                 
        ! (/ amplitude, sigma of ricker pulse, time to start ricker pulse,
        !    distance to points, number of points, X coord, Y coord, ... /)
        m_pointamp   = m_pparams(1)
        m_pointsigma = m_pparams(2)
        m_pointstart = m_pparams(3)
        dist         = m_pparams(4)
        nullify(lanodes)
        do j = 5, size(m_pparams), 2
        pt = m_pparams(j:j+1)
        do i = 1, m_nnodes
                if ( distance(pt, m_nodes(:,i)) < dist ) then
                        call insertll(lanodes, i)
                end if
        end do
        end do
        call toarray(lanodes, m_pointnodes)
        if ( allocated(m_pointnodes) ) then
                call get_unique(m_pointnodes)
        end if
end subroutine

subroutine create_absorbingbc()
        integer :: ind
        real*8 :: left, right, top
        left  = real(0,kind=8)*m_blength
        right = real(m_totcols-1,kind=8)*m_blength
        top   = real(m_totrows-1,kind=8)*m_blength
        allocate(m_absorbing(m_nnodes))
        m_absorbing = AB_NOABS
        do ind = 1, m_nnodes 
                if ( realequal(m_nodes(1,ind), left)  .OR.  &
                     realequal(m_nodes(1,ind), right) ) then
                        m_absorbing(ind) = AB_XDIR
                else if ( realequal(m_nodes(2,ind), top) ) then
                        m_absorbing(ind) = AB_ZDIR
                end if
        end do
end subroutine

! Sets up stiffness matrices.
subroutine init_stiffness()
        integer :: i
        allocate(m_kmat(8,8,m_nummat))
        do i=1,m_nummat
                call calck(m_matconst(1,i), m_matconst(2,i), &
                           m_kmat(:,:,i))
        end do
end subroutine


subroutine load_problem()
        real*8  :: height, width, time, chper
        integer :: nmats, nppars, ngpars
        type(c_ptr) :: cmptr, vmptr, vpptr, vgptr
        character(len=1), dimension(:,:), allocatable :: cmats
        real*8, dimension(:,:), allocatable           :: vmats
        character(len=256) :: filename
        integer :: fnlen
        fnlen = io_fnindlen + 11
        write(filename,"(A,A11)") io_indir(1:io_fnindlen),"/problem.in"
        call load_problem_file(fnlen, filename, height, width, time, &
                               m_inttime, chper, f_nperlen, f_sperlen,   &
                               nmats, cmptr, vmptr, m_prob, nppars,  &
                               vpptr, m_geom, ngpars, vgptr)
        call ctofort(cmptr, cmats, (/2, nmats/))
        call ctofort(vmptr, vmats, (/3, nmats/))
        call ctofort(vgptr, m_gparams, (/ ngpars /))
        call ctofort(vpptr, m_pparams, (/ nppars /))
        call setup_materials(nmats, cmats, vmats)
        call setup_problem_params(height, width, time, chper)
        deallocate(cmats, vmats)
end subroutine

! Sets the interval in which 
subroutine set_save_interval()
        write(*,*) m_inttime
        io_interval = int(m_inttime/m_deltat)
        m_inttime = m_deltat * real(io_interval, kind=8)
end subroutine

end module model
