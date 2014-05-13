module pwaves 
use utilities
implicit none

! precalculated values


real*8, parameter      :: f_wangzero = 0.0001
real*8                 :: f_walpha, f_wbeta, f_wsigma, f_wfreq
real*8                 :: f_wsine, f_wcose, f_wsinf, f_wcosf
real*8                 :: f_wa, f_wa1, f_wb, f_wb1
real*8, dimension(4)   :: f_wcons
real*8, dimension(2,3) :: f_wdv, f_wui
contains

subroutine initplanewaves(pwave, mag, alpha, beta, angle, freq)
        character(len=6), intent(in) :: pwave
        real*8, intent(in) :: mag, alpha, beta, angle, freq
        real*8 :: a, b, e, f, sin2e, sin2f, cos2f, cos2e, fcr
        integer :: i
        f_walpha = alpha
        f_wbeta  = beta
        f_wfreq  = freq
        if ( pwave == 'WAVEP ' ) then
                e = angle
                a = mag
                f_wsine = sin(e)
                f_wsinf = f_wsine*(beta/alpha)
                f = asin(f_wsinf)
                f_wcose = cos(e)
                f_wcosf = cos(f)
                if ( f_wcose < f_wangzero ) then
                        ! Grazing incidence. Not implemented.
                        write(*,*) "Grazing incidence. Not implemented."
                        call exit(1)
                else if ( abs(f_wsine) < f_wangzero ) then
                        f_wa = a
                        f_wa1 = -f_wa
                        f_wb1 = 0.0D0
                        f_wcose = 1
                        f_wsine = 0
                else
                        sin2e = sin(2*e)
                        sin2f = sin(2*f)
                        cos2f = cos(2*f)
                        f_wa  = a
                        f_wa1 = (sin2e*sin2f-(alpha/beta)**2*(cos2f**2))
                        f_wa1 = f_wa1/(sin2e*sin2f+(alpha/beta)**2*(cos2f**2))
                        f_wa1 = f_wa1*a
                        f_wb1 = -(2*(alpha/beta)*sin2e*cos2f)
                        f_wb1 = f_wb1/(sin2e*sin2f+(alpha/beta)**2*(cos2f**2))
                        f_wb1 = f_wb1*a
                end if
                f_wui(1,1) =  f_wsine
                f_wui(2,1) = -f_wcose
                f_wui(1,2) =  f_wsine
                f_wui(2,2) =  f_wcose
                f_wui(1,3) = -f_wcosf
                f_wui(2,3) =  f_wsinf
        else if ( pwave == 'WAVESV' ) then
                f = angle
                b = mag
                f_wsinf = sin(f)
                fcr = asin(beta/alpha)
                f_wsine = f_wsinf*(alpha/beta)
                e = asin(f_wsine)
                f_wcose = cos(e)
                f_wcosf = cos(f)
                
                if ( abs(f_wsinf) < f_wangzero ) then
                        f_wb = b
                        f_wb1 = -f_wb
                        f_wa1 = 0.0D0
                        f_wcosf = 1
                        f_wsinf = 0
                else if ( (fcr-f) > f_wangzero ) then
                        sin2e = sin(2*e)
                        sin2f = sin(2*f)
                        cos2f = cos(2*f)
                        f_wb  = b
                        f_wb1 = (sin2e*sin2f-(alpha/beta)**2*(cos2f**2))
                        f_wb1 = f_wb1/(sin2e*sin2f+(alpha/beta)**2*(cos2f**2))
                        f_wb1 = f_wb1*b
                        f_wa1 = (alpha/beta)*sin(4.0D0*f)
                        f_wa1 = f_wa1/(sin2e*sin2f+(alpha/beta)**2*(cos2f**2))
                        f_wa1 = f_wa1*b
                else
                        write(*,*) "Inhomogeneous waves. Not implemented"
                        call exit(1)
                end if
                f_wui(1,1) =  f_wcosf
                f_wui(2,1) =  f_wsinf
                f_wui(1,2) =  f_wsine
                f_wui(2,2) =  f_wcose
                f_wui(1,3) = -f_wcosf
                f_wui(2,3) =  f_wsinf
        end if
end subroutine

subroutine planewave(pwave, t, x1, x3, u)
        character(len=6), intent(in) :: pwave
        real*8, intent(in) :: t, x1, x3
        real*8, dimension(2,3), intent(out) :: u
        if ( pwave == 'WAVEP ' ) then
                call wavep(t, x1, x3, u)
        else if ( pwave == 'WAVESV' ) then
                call wavesv(t, x1, x3, u)
        else
                u = 0
        end if
end subroutine

subroutine wavepold(t, x1, x3, u)
        real*8, intent(in) :: t, x1, x3
        real*8, dimension(2,3), intent(out) :: u
        real*8, dimension(2,3) :: ip, rp, rs
        real*8, dimension(3) :: v, ev2
        v(1) = (x1*f_wsine-x3*f_wcose)/f_walpha
        v(2) = (x1*f_wsine+x3*f_wcose)/f_walpha
        v(3) = (x1*f_wsinf+x3*f_wcosf)/f_wbeta
        
        ev2  = exp(-(t-v)**2/(4*f_wsigma))
        
        ip(:,1) = -ev2(1)*(2*f_wsigma-(t-v(1))**2)/f_wcons(2)
        ip(:,2) = ev2(1)*(t-v(1))*(6*f_wsigma-(t-v(1))**2)/f_wcons(3)
        ip(:,3) = ev2(1)*((t-v(1))**4-12.0D0*f_wsigma*(t-v(1))**2+12.0D0*f_wsigma**2)/f_wcons(4)
        
        ip(1,:) = ip(1,:)*f_wui(1,1)*f_wa
        ip(2,:) = ip(2,:)*f_wui(2,1)*f_wa
        
        rp(:,1) = -ev2(2)*(2*f_wsigma-(t-v(2))**2)/f_wcons(2)
        rp(:,2) = ev2(2)*(t-v(2))*(6*f_wsigma-(t-v(2))**2)/f_wcons(3)
        rp(:,3) = ev2(2)*((t-v(2))**4-12.0D0*f_wsigma*(t-v(2))**2+12.0D0*f_wsigma**2)/f_wcons(4)
        
        rp(1,:) = rp(1,:)*f_wui(1,2)*f_wa1
        rp(2,:) = rp(2,:)*f_wui(2,2)*f_wa1
        
        rs(:,1) = -ev2(3)*(2*f_wsigma-(t-v(3))**2)/f_wcons(2)
        rs(:,2) = ev2(3)*(t-v(3))*(6*f_wsigma-(t-v(3))**2)/f_wcons(3)
        rs(:,3) = ev2(3)*((t-v(3))**4-12.0D0*f_wsigma*(t-v(3))**2+12.0D0*f_wsigma**2)/f_wcons(4)
        
        rs(1,:) = rs(1,:)*f_wui(1,3)*f_wb1
        rs(2,:) = rs(2,:)*f_wui(2,3)*f_wb1
        
        u = ip + rp + rs
end subroutine

subroutine wavep(t, x1, x3, u)
        real*8, intent(in) :: t, x1, x3
        real*8, dimension(2,3), intent(out) :: u
        real*8, dimension(2,3) :: ip, rp, rs
        real*8, dimension(3) :: v, v2
        v(1) = (x1*f_wsine-x3*f_wcose)/f_walpha
        v(2) = (x1*f_wsine+x3*f_wcose)/f_walpha
        v(3) = (x1*f_wsinf+x3*f_wcosf)/f_wbeta
        
        v2   = (f_wfreq*u_pi*(t-v))**2
      
        ip(:,1) = (1.0D0-2.0D0*v2(1))/exp(v2(1))
        ip(:,2) = (2.0D0*v2(1)-3.0D0)*(f_wfreq*u_pi)**2*(t-v(1))/exp(v2(1))
        ip(:,3) = -2.0D0*(4.0D0*v2(1)**2-12.0D0*v2(1)+3)*(f_wfreq*u_pi)**2/exp(v2(1))
        
        rp(:,1) = (1.0D0-2.0D0*v2(2))/exp(v2(2))
        rp(:,2) = (2.0D0*v2(2)-3.0D0)*(f_wfreq*u_pi)**2*(t-v(2))/exp(v2(2))
        rp(:,3) = -2.0D0*(4.0D0*v2(2)**2-12.0D0*v2(2)+3)*(f_wfreq*u_pi)**2/exp(v2(2))

        rs(:,1) = (1.0D0-2.0D0*v2(3))/exp(v2(3))
        rs(:,2) = (2.0D0*v2(3)-3.0D0)*(f_wfreq*u_pi)**2*(t-v(3))/exp(v2(3))
        rs(:,3) = -2.0D0*(4.0D0*v2(3)**2-12.0D0*v2(3)+3)*(f_wfreq*u_pi)**2/exp(v2(3))
        
        ip(1,:) = ip(1,:)*f_wui(1,1)*f_wa
        ip(2,:) = ip(2,:)*f_wui(2,1)*f_wa
        
        rp(1,:) = rp(1,:)*f_wui(1,2)*f_wa1
        rp(2,:) = rp(2,:)*f_wui(2,2)*f_wa1
        
        rs(1,:) = rs(1,:)*f_wui(1,3)*f_wb1
        rs(2,:) = rs(2,:)*f_wui(2,3)*f_wb1
        
        u = ip + rp + rs
end subroutine

subroutine wavesv(t, x1, x3, u)
        real*8, intent(in) :: t, x1, x3
        real*8, dimension(2,3), intent(out) :: u
        real*8, dimension(2,3) :: is, rp, rs
        real*8, dimension(3) :: v, v2
        
        v(1) = (x1*f_wsinf-x3*f_wcosf)/f_wbeta
        v(2) = (x1*f_wsine+x3*f_wcose)/f_walpha
        v(3) = (x1*f_wsinf+x3*f_wcosf)/f_wbeta
        
        v2   = (f_wfreq*u_pi*(t-v))**2
      
        is(:,1) = (1.0D0-2.0D0*v2(1))/exp(v2(1))
        is(:,2) = (2.0D0*v2(1)-3.0D0)*(f_wfreq*u_pi)**2*(t-v(1))/exp(v2(1))
        is(:,3) = -2.0D0*(4.0D0*v2(1)**2-12.0D0*v2(1)+3.0D0)*(f_wfreq*u_pi)**2/exp(v2(1))
        
        rp(:,1) = (1.0D0-2.0D0*v2(2))/exp(v2(2))
        rp(:,2) = (2.0D0*v2(2)-3.0D0)*(f_wfreq*u_pi)**2*(t-v(2))/exp(v2(2))
        rp(:,3) = -2.0D0*(4.0D0*v2(2)**2-12.0D0*v2(2)+3.0D0)*(f_wfreq*u_pi)**2/exp(v2(2))

        rs(:,1) = (1.0D0-2.0D0*v2(3))/exp(v2(3))
        rs(:,2) = (2.0D0*v2(3)-3.0D0)*(f_wfreq*u_pi)**2*(t-v(3))/exp(v2(3))
        rs(:,3) = -2.0D0*(4.0D0*v2(3)**2-12.0D0*v2(3)+3.0D0)*(f_wfreq*u_pi)**2/exp(v2(3))
        
        is(1,:) = is(1,:)*f_wui(1,1)*f_wb
        is(2,:) = is(2,:)*f_wui(2,1)*f_wb
        
        rp(1,:) = rp(1,:)*f_wui(1,2)*f_wa1
        rp(2,:) = rp(2,:)*f_wui(2,2)*f_wa1
        
        rs(1,:) = rs(1,:)*f_wui(1,3)*f_wb1
        rs(2,:) = rs(2,:)*f_wui(2,3)*f_wb1
        
        u = is + rp + rs
end subroutine

end module pwaves

