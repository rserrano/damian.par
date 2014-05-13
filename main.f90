program exfem

        use explicit
        use visualization 
        implicit none
        integer :: nargs
        character(len=256) :: arg
        nargs = command_argument_count() 
        if ( nargs < 2 ) then
                
                call helpmsg()
        else
                call get_command_argument(1, arg)
                select case(arg)
                case("solve")
                        call get_command_argument(2, io_indir)
                        io_fnindlen = len_trim(io_indir)
                        call solve()
                case("visualize")
                        call get_command_argument(2, io_indir)
                        io_fnindlen = len_trim(io_indir)
                        call visualize()
                case("help")
                        call helpmsg()
                case default
                        call helpmsg()
                end select
        end if
contains
        subroutine visualize()
                call init_mpi()
                if ( c_rank == 0 ) then
                        write(*,*) "VISUALIZING MODEL", io_indir
                end if
                ! call wait_debug()
                call init_visualization()
                call load_io_ranks()
                call create_visualization()
                call end_mpi()
        end subroutine

        subroutine solve()
                call init_mpi()
                if ( c_rank == 0 ) then
                        write(*,*) "SOLVING MODEL: ", io_indir
                end if
                call init_model()
                call save_io_ranks()
                call init_explicit()
                ! call wait_debug()
                call run_explicit()
                call end_mpi()
        end subroutine

        subroutine helpmsg()
                call init_mpi()
                if ( c_rank == 0 ) then
                write(*,*) "Usage:      mpirun -np [N] ./exfem [option] [params]"
                write(*,*) "solve:      solves the current problem specified"
                write(*,*) "            in dir [param]."
                write(*,*) "visualize:  creates vtk files from the binary"
                write(*,*) "            solution for the range selected."
                end if
                call end_mpi()
        end subroutine
end program exfem
