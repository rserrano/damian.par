module explicit
        use model
        use formulation
        use pwaves
        use communicator
        use io_utils 
        implicit none
        real*8, dimension(:), allocatable  :: e_mfull
        real*8, dimension(4)               :: e_a
        real*8, dimension(:), allocatable  :: e_u
        real*8, dimension(:), allocatable  :: e_up
        real*8, dimension(:), allocatable  :: e_ro
        integer, dimension(:), allocatable :: e_watch
        contains
        
        subroutine init_model()
                
                logical*1 :: further
                character(len=27) :: elstr
                integer :: timer_inter, timer_total

                call start_timer(timer_total)
                if ( c_rank == 0 ) then
                        call load_problem()
                end if
                call bcast_problem()
                if ( c_rank == 0 ) then
                        call load_environ()
                end if
                call bcast_environ()
                if ( m_geom(1:3) == 'CVM' ) then
                        call init_cvm()
                end if
                if ( m_prob(1:4) == 'WAVE' ) then
                        call init_wave()
                end if
                call global_initial_mesh()
                
                call update_elmdist()
                further = logical(.TRUE., kind=1)
                call start_timer(timer_inter)
                do while(further)
                ! call wait_debug()
                call global_cut_geom()
                call redistribute()
                call global_divide(further)
                call end_timer(elstr, timer_inter)
                if ( c_rank == 0 ) then
                write(*,'(A,A)') "Iterated model. ", elstr
                end if
                call start_timer(timer_inter)
                end do
                call problem_setup()
                call init_exfem_comm()
                call calc_global_totnodes()
                if ( m_geom(1:3) == 'CVM' ) then
                        call end_cvm()
                end if
                call recalculate_deltat()
                call end_timer(elstr, timer_total)
                if ( c_rank == 0 ) then
                        write(*,'(A,A)') "Finished constructing model. ", elstr
                        write(*,*) "elements:  ", m_totelems
                        write(*,*) "dofs:      ", m_totnodes*2
                        write(*,*) "no steps:  ", m_nsteps
                        write(*,*) "time step: ", m_deltat
                end if
                call save_mesh(c_rank, m_eptr, m_eind, m_ediv, m_nodes)
                call save_tparams(c_rank, m_deltat, io_interval)
        end subroutine
        
        subroutine setupconstants()
                e_a(1) = 1.0D0/m_deltat**2.0D0
                e_a(2) = 1.0D0/(2.0D0*m_deltat)
                e_a(3) = 2.0D0*e_a(1)
                e_a(4) = 1.0D0/e_a(2)
        end subroutine
        
        subroutine init_explicit()
                implicit none
                integer :: i
                integer, dimension(8) :: dofs
                real*8 :: mass
                real*8, dimension(8) :: damping

                call setupconstants()
                allocate(e_u(2*m_nnodes))
                allocate(e_up(2*m_nnodes))
                allocate(e_mfull(2*m_nnodes))
                allocate(e_ro(2*m_nnodes))
                e_u(:) = 0.0D0
                e_up(:) = 0.0D0
                call problem_displacements(0.0D0, e_u)
                call problem_displacements(-m_deltat, e_up)
                e_mfull(:) = 0.0D0
                do i = 1, m_nelems
                        call elemallnhdofs(i, dofs)
                        call getmassof(i,mass)
                        
                        call getdampingof(i,damping)
                        e_mfull(dofs) = e_mfull(dofs) +  &
                                        e_a(2)*damping + &
                                        e_a(1)*mass
                end do
                call sum_neighbor_reactions(e_mfull)
                call sum_reactanchors(e_mfull)
                call sum_neighbor_anchors(e_mfull)
                call save_mass(c_rank, e_mfull)
        end subroutine
        
        subroutine run_explicit()
                implicit none
                integer :: i, j
                !logical :: debug
                real*8, dimension(8)    :: localu, localup, localr, damping
                real*8  :: mass
                integer, dimension(8)   :: dofs
                integer :: timer_step, timer_total
                character(len=27) :: timstr
                integer :: maxind, minind, k, l, ndind

                call start_timer(timer_total)
                call start_timer(timer_step)
                call save_step(c_rank, 0, e_u, timer_step)
                do i=1, m_nsteps
                        e_ro(:) = 0.0D0
                        !! Check what to do with punctual loads.
                        call additive_reactions(real(i,8)*m_deltat, e_ro) 
                        do j=1, m_nelems
                                call elemallnhdofs(j, dofs)
                                call getmassof(j,mass)
                                call getdampingof(j,damping)
                                localu = e_u(dofs)
                                localup = e_up(dofs)
                                localr = (e_a(3)*localu(:) - &
                                          e_a(1)*localup(:))*mass
                                localr = localr + e_a(2)*damping*localup
                                
                                call dgemv('N',8,8,-1.0D0,m_kmat(:,:,m_ematprop(j)),&
                                           8,localu,1,1.0D0,localr,1)
                                e_ro(dofs)=e_ro(dofs) + localr
                        end do
                        call sum_neighbor_reactions(e_ro)
                        call absolute_reactions(real(i,8)*m_deltat, e_ro)
                        call sum_reactanchors(e_ro)
                        call sum_neighbor_anchors(e_ro)
                        e_up(:) = e_u(:)
                        e_u(:)  = e_ro(:)/e_mfull(:)
                        
                        call problem_displacements(real(i,8)*m_deltat, e_u)
                        call set_dispshanging(e_u)
                        call set_neighbor_displacements(e_u)
                        call save_step(c_rank, i, e_u, timer_step)
                end do
                call end_timer(timstr, timer_total)
                if ( c_rank == 0 ) then
                        write(*,'(A,A)') "Solved exfem problem. ", timstr
                end if
        end subroutine

end module explicit

