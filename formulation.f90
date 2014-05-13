module formulation
use utilities
implicit none

integer, parameter :: f_young  = 1, &
                      f_nu     = 2, &
                      f_shear  = 3, &
                      f_bulk   = 4, &
                      f_lambda = 5, &
                      f_vp     = 1, &
                      f_vs     = 2 
! Parameters of the discretization.

real*8 :: f_nperlen
real*8 :: f_sperlen

! precalculated values
contains

subroutine explicitparams(chperiod, alpha, beta, length, step )
        real*8, intent(in)  :: chperiod, alpha, beta
        real*8, intent(out) :: length, step
        length = (chperiod*beta)/f_nperlen
        step   = (length/alpha)/f_sperlen
end subroutine

subroutine recalcstep(alpha, length, step)
        real*8, intent(in)  :: alpha, length
        real*8, intent(out) :: step
        step   = (length/alpha)/f_sperlen
end subroutine

subroutine pointtransform(origin, inp, outp)
        real*8, dimension(2), intent(in)  :: origin
        real*8, dimension(2), intent(in)  :: inp
        real*8, dimension(2), intent(out) :: outp
        outp(1) = inp(1)-origin(1)
        outp(2) = inp(2)-origin(2)
end subroutine

subroutine timefromlower(height, width, e, alpha, time)
        real*8, intent(in) :: height, width, e, alpha
        real*8, intent(out) :: time
        time = (height+(width/2)*tan(e))/alpha
end subroutine

subroutine timetransform(start, inti, outti)
        real*8, intent(in) :: start, inti
        real*8, intent(out) :: outti
        outti = inti-start
end subroutine

subroutine ricker(amp, sigma, t, r)
        real*8, intent(in) :: amp, sigma, t
        real*8, intent(out) :: r
        
        r = exp(-t**2/(4*sigma))*(t**2-2*sigma)/(8*sqrt(u_pi)*sigma**(2.5D0))
        r = r*amp
end subroutine

! Calculates a mechanical parameter from velocity and density.
subroutine v_to_e(param, rho, pto, res)
        real*8, dimension(2), intent(in) :: param
        real*8, intent(in) :: rho
        integer, dimension(:), intent(in) :: pto
        real*8, dimension(:), intent(out) :: res
        real*8, dimension(2) :: lame
        lame(1) = rho*param(f_vs)**2
        lame(2) = rho*param(f_vp)**2-2*lame(1)
        call e_to_e((/f_shear, f_lambda/), pto, lame, res)
end subroutine
! Calculates the velocity from a mechanical parameter from velocity and
! density.
subroutine e_to_v(param, rho, pfrom, res)
        real*8, dimension(2), intent(in) :: param
        real*8, intent(in) :: rho
        integer, dimension(2), intent(in) :: pfrom
        real*8, dimension(2), intent(out) :: res
        real*8, dimension(2) :: lame
        call e_to_e(pfrom, (/f_shear, f_lambda/), param, lame)
        res(f_vp) = sqrt((lame(2)+2*lame(1))/rho)
        res(f_vs) = sqrt(lame(1)/rho)
end subroutine

subroutine v_nu_to_v(pfrom, vin, nu, res)
        real*8, intent(in)  :: vin, nu
        integer, intent(in) :: pfrom
        real*8, dimension(2), intent(out) :: res
        real*8 :: rel
        rel = sqrt(2.0D0*(nu-1.0D0)/(2.0D0*nu-1.0D0))
        res(pfrom) = vin
        if ( pfrom == f_vs ) then
                res(f_vp) = rel*res(f_vs)
        else if ( pfrom == f_vp ) then
                res(f_vs) = res(f_vp)/rel
        end if
end subroutine

! Returns the numeric parameter from the character used to describe it.
function c_to_p(ch)
        integer :: c_to_p
        character(len=*), intent(in) :: ch
        c_to_p = 0
        select case(ch)
        case("E")
                c_to_p = f_young
        case("v")
                c_to_p = f_nu
        case("G")
                c_to_p = f_shear
        case("K")
                c_to_p = f_bulk
        case("l")
                c_to_p = f_lambda
        case default
                c_to_p = 0
        end select
end function

function c_to_v(ch)
        integer :: c_to_v
        character(len=*), intent(in) :: ch
        c_to_v = 0
        select case(ch)
        case("p")
                c_to_v = f_vp
        case("s")
                c_to_v = f_vs
        case default
                c_to_v = 0
        end select
end function

! Calculates a mechanical parameter from another mechanical parameter.
subroutine e_to_e(pfrom, pto, param, res)
        integer, dimension(2), intent(in)  :: pfrom
        real*8, dimension(2), intent(in)   :: param
        integer, dimension(:), intent(in)  :: pto
        real*8, dimension(:), intent(out)  :: res
        logical*1, dimension(5)            :: found
        real*8, dimension(5)               :: from
        integer                            :: i,j
        from(:) = 0.0D0
        found(:) = .FALSE.
        do i=1,2
                from(pfrom(i)) = param(i)
                found(pfrom(i)) = .TRUE.
        end do
        
        do i = 1,3
        if(.NOT. any(.NOT. found(pto(:)))) then
                exit
        end if
        do j = 1,5
                if (found(j)) then
                        cycle
                end if
                select case(j)
                case(f_young)
                if (found(f_shear) .AND. found(f_nu)) then
                found(f_young) = .TRUE.
                from(f_young) = 2*(1+from(f_nu))*from(f_shear)
                else if (found(f_bulk) .AND. found(f_nu)) then
                found(f_young) = .TRUE.
                from(f_young) = 3*(1-2*from(f_nu))*from(f_bulk)
                else if (found(f_lambda) .AND. found(f_nu)) then
                found(f_young) = .TRUE.
                from(f_young) = (1-2*from(f_nu))*(1+from(f_nu))
                from(f_young) = from(f_young)*from(f_lambda)/from(f_nu)
                end if
                case(f_nu)
                if (found(f_young) .AND. found(f_shear)) then
                found(f_nu) = .TRUE.
                from(f_nu) = from(f_young)-2*from(f_shear)
                from(f_nu) = from(f_nu)/(2*from(f_shear))
                else if (found(f_young) .AND. found(f_bulk)) then
                found(f_nu) = .TRUE.
                from(f_nu) = 3*from(f_bulk)-from(f_young)
                from(f_nu) = from(f_nu)/(6*from(f_bulk))
                else if (found(f_bulk) .AND. found(f_shear)) then
                found(f_nu) = .TRUE.
                from(f_nu) = 3*from(f_bulk)-2*from(f_shear)
                from(f_nu) = from(f_nu)/(6*from(f_bulk)+2*from(f_shear))
                else if (found(f_shear) .AND. found(f_lambda)) then
                found(f_nu) = .TRUE.
                from(f_nu) = 2*(from(f_lambda)+from(f_shear))
                from(f_nu) = from(f_lambda)/from(f_nu)
                else if (found(f_bulk) .AND. found(f_lambda)) then
                found(f_nu) = .TRUE.
                from(f_nu) = 3*from(f_bulk)-from(f_lambda)
                from(f_nu) = from(f_lambda)/from(f_nu)
                end if
                case(f_shear)
                if (found(f_young) .AND. found(f_nu)) then
                found(f_shear) = .TRUE.
                from(f_shear) = from(f_young)/(2*(from(f_nu)+1))
                end if
                case(f_bulk)
                if (found(f_young) .AND. found(f_nu)) then 
                found(f_bulk) = .TRUE.
                from(f_bulk) = from(f_young)/(3*(1-2*from(f_nu)))
                end if
                case(f_lambda)
                if ( found(f_young) .AND. found(f_nu)) then
                found(f_lambda) = .TRUE.
                from(f_lambda) = from(f_young)*from(f_nu)
                from(f_lambda) = from(f_lambda)/((1+from(f_nu))*(1-2*from(f_nu)))
                end if
                end select
        end do
        end do
        if ( any(.NOT. found( pto(:) ) ) ) then
                call error_exfem("Conversión no válida.")
                
        end if
        res(:) = from( pto(:) )
end subroutine

subroutine lysmer(scl, rho, v, iszdir, res)
        real*8,               intent(in)  :: scl, rho
        real*8, dimension(2), intent(in)  :: v
        logical*1,              intent(in)  :: iszdir
        real*8, dimension(2), intent(out) :: res
        res(1) = v(1)*rho*scl
        res(2) = v(2)*rho*scl
        if ( iszdir ) then
                res((/1,2/)) = res((/2,1/))
        endif
end subroutine

! Calculates a body force at the boundary.
subroutine calcbf(scl, rho, g, bf)
        real*8, intent(in)                :: rho, scl
        real*8, dimension(2), intent(in)  :: g
        real*8, dimension(8), intent(out) :: bf
        bf(1:7:2) = g(1)*(scl**2)*rho
        bf(2:8:2) = g(2)*(scl**2)*rho
end subroutine
! Calculates a surface force.
subroutine calcsf(scl, f, sf)
        real*8, intent(in) :: scl
        real*8, dimension(2), intent(in)  :: f
        real*8, dimension(4), intent(out) :: sf
        sf = (/ f(1)*scl, &
                f(2)*scl, &
                f(1)*scl, &
                f(2)*scl /)
end subroutine
! Calculates the local mass matrix.
subroutine calcm(scl, rho, m)
        real*8, intent(out) :: m
        real*8, intent(in)  :: scl, rho
        m = (scl**2.0D0)*rho
end subroutine

subroutine calcsqm(m)
        real*8, intent(out), dimension(8,8) :: m
        m = reshape( (/4.0D0,0.0D0,2.0D0,0.0D0,1.0D0,0.0D0,2.0D0,0.0D0, &
                       0.0D0,4.0D0,0.0D0,2.0D0,0.0D0,1.0D0,0.0D0,2.0D0, &
                       2.0D0,0.0D0,4.0D0,0.0D0,2.0D0,0.0D0,1.0D0,0.0D0, &
                       0.0D0,2.0D0,0.0D0,4.0D0,0.0D0,2.0D0,0.0D0,1.0D0, &
                       1.0D0,0.0D0,2.0D0,0.0D0,4.0D0,0.0D0,2.0D0,0.0D0, &
                       0.0D0,1.0D0,0.0D0,2.0D0,0.0D0,4.0D0,0.0D0,2.0D0, &
                       2.0D0,0.0D0,1.0D0,0.0D0,2.0D0,0.0D0,4.0D0,0.0D0, &
                       0.0D0,2.0D0,0.0D0,1.0D0,0.0D0,2.0D0,0.0D0,4.0D0  &
                      /), shape(m) )
        m = m/9.0D0 
end subroutine

! Inverts the mass matrix.
subroutine invertmfull(mfull)
        real*8, dimension(:), intent(inout) :: mfull
        mfull(:) = 1.0D0 / mfull(:)
end subroutine
! Calculates the local stiffness matrix.
subroutine calck(m, l, k)
        real*8, intent(in)                  :: m, l
        real*8, dimension(8,8), intent(out) :: k
        k = reshape( (/                                &
                (3.0D0*m+l)/3.0D0,(m+l)/4.0D0,         &
                -(3.0D0*m+2.0D0*l)/6.0D0,-(m-l)/4.0D0, &
                -(3.0D0*m+l)/6.0D0,-(m+l)/4.0D0,       &
                l/6.0D0,(m-l)/4.0D0,                   &
                (m+l)/4.0D0,(3.0D0*m+l)/3.0D0,         &
                (m-l)/4.0D0,l/6.0D0,                   &
                -(m+l)/4.0D0,-(3.0D0*m+l)/6.0D0,       &
                -(m-l)/4.0D0,-(3.0D0*m+2.0D0*l)/6.0D0, &
                -(3.0D0*m+2.0D0*l)/6.0D0,(m-l)/4.0D0,  &
                (3.0D0*m+l)/3.0D0,-(m+l)/4.0D0,        &
                l/6.0D0,-(m-l)/4.0D0,                  &
                -(3.0D0*m+l)/6.0D0,(m+l)/4.0D0,        &
                -(m-l)/4.0D0,l/6.0D0,                  &
                -(m+l)/4.0D0,(3.0D0*m+l)/3.0D0,        &
                (m-l)/4.0D0,-(3.0D0*m+2.0D0*l)/6.0D0,  &
                (m+l)/4.0D0,-(3.0D0*m+l)/6.0D0,        &
                -(3.0D0*m+l)/6.0D0,-(m+l)/4.0D0,       &
                l/6.0D0,(m-l)/4.0D0,                   &
                (3.0D0*m+l)/3.0D0,(m+l)/4.0D0,         &
                -(3.0D0*m+2.0D0*l)/6.0D0,-(m-l)/4.0D0, &
                -(m+l)/4.0D0,-(3.0D0*m+l)/6.0D0,       &
                -(m-l)/4.0D0,-(3.0D0*m+2.0D0*l)/6.0D0, &
                (m+l)/4.0D0,(3.0D0*m+l)/3.0D0,         &
                (m-l)/4.0D0,l/6.0D0,                   &
                l/6.0D0,-(m-l)/4.0D0,                  &
                -(3.0D0*m+l)/6.0D0,(m+l)/4.0D0,        &
                -(3.0D0*m+2.0D0*l)/6.0D0,(m-l)/4.0D0,  &
                (3.0D0*m+l)/3.0D0,-(m+l)/4.0D0,        &
                (m-l)/4.0D0,-(3.0D0*m+2.0D0*l)/6.0D0,  &
                (m+l)/4.0D0,-(3.0D0*m+l)/6.0D0,        &
                -(m-l)/4.0D0,l/6.0D0,                  &
                -(m+l)/4.0D0,(3.0D0*m+l)/3.0D0         &
        /), shape(k) )
end subroutine

subroutine calcfieldat(nodes, dofs, pt, field)
        real*8, dimension(2,4), intent(in) :: nodes
        real*8, dimension(8), intent(in)   :: dofs
        real*8, dimension(2), intent(in)   :: pt
        real*8, dimension(2), intent(out)  :: field
        real*8, dimension(2)   :: param
        real*8, dimension(2,8) :: nmat
        call shapeparams(nodes, pt, param)
        call calcnmat(param, nmat)
        call calcfield(nmat, dofs, field)
end subroutine

subroutine calcfield(nmat, dofs, field)
        real*8, dimension(2,8), intent(in) :: nmat
        real*8, dimension(8), intent(in)   :: dofs
        real*8, dimension(2), intent(out)  :: field
        field = matmul(nmat, dofs)
end subroutine

subroutine shapeparams(nodes, point, param)
        real*8, dimension(2,4), intent(in) :: nodes
        real*8, dimension(2), intent(in)   :: point
        real*8, dimension(2), intent(out)  :: param
        param(1) = (point(1)-nodes(1,1))/(nodes(1,3)-nodes(1,1))
        param(2) = (point(2)-nodes(2,1))/(nodes(2,3)-nodes(2,1))
        param = 2.0D0*param - 1.0D0
        param = max((/-1.0D0,-1.0D0/), param)
        param = min((/ 1.0D0, 1.0D0/), param)
end subroutine

subroutine calcnmat(param, nmat)
        real*8, dimension(2), intent(in)    :: param
        real*8, dimension(2,8), intent(out) :: nmat
        real*8, dimension(4) :: inter
        call shapefuncs(param(1), param(2), inter)
        nmat(:,:)     = 0.0D0
        nmat(1,1:7:2) = inter(:)
        nmat(2,2:8:2) = inter(:)
end subroutine

! Calculates the shape functions, (not used).
subroutine shapefuncs(r, s, inter)
        real*8, intent(in) :: r,s
        real*8, dimension(4), intent(out) :: inter
        inter = (/ (1-r)*(1-s)/4, &
                   (1+r)*(1-s)/4, &
                   (1+r)*(1+s)/4, &
                   (1-r)*(1+s)/4 /)
end subroutine
! Calculates the derivatives of the shape functions, (not used).
subroutine derivs(r, s, deriv)
        real*8, intent(in) :: r,s
        real*8, dimension(4,2), intent(out) :: deriv
        deriv =  reshape( (/ (s-1)/4,   &
                             (1-s)/4,   &
                             (s+1)/4,   &
                             (-1-s)/4,  &
                             (r-1)/4,   &
                             (-1-r)/4,  &
                             (r+1)/4,   &
                             (1-r)/4 /), shape(deriv))
end subroutine

subroutine multdrm(k, u, e, p)
        real*8, dimension(8,8), intent(in)  :: k
        real*8, dimension(8), intent(in)    :: u
        logical*1, dimension(8), intent(in) :: e
        real*8, dimension(8), intent(out)   :: p
        integer :: i, j
        p = 0.0D0
        do i = 1,8
        do j = 1,8
        if ( e(j) .AND. .NOT. e(i) ) then
                p(i) = p(i)-k(i,j)*u(j)
                p(j) = p(j)+k(j,i)*u(i)
        end if
        end do
        end do
end subroutine


end module formulation

