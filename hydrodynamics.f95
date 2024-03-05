module hydrodynamics

USE global_variables
!USE phys_prop
USE tools
    implicit none

 real*8 :: A_r = 3.D0                  !Turbulence parameter
! real*8 :: B_ = 1.86578208D0    ! Turbulence Parameter Pall Ring   2.212797246D0  ! Turbulence Parameter Envipac  1.86578208D0
 real*8 :: Kappa = 0.41D0              !Kharman constant
 real*8 :: C_mu_turb                   !Turbulence parameter
 real*8 :: u_tau                       !Turbulence parameter
 real*8 :: sigma                           ! surface tension
 real*8 :: rho_g, rho_l
 real*8 :: mu_l, mu_g
 real*8 :: q_l, q_g
 real*8 :: alpha=pi/4
 real*8 :: eps_l,eps_g                 ! pressure drop term of liquid and gas side
 real*8 :: C_l1,C_l2,C_g2,C_g1!,C_turb1,C_turb2,C_turb3,C_turb4              ! integration constants of liquid and gas flow

 logical :: turbulent=.true.                       ! if true turbulent calculations will be performed else laminar
 logical, dimension (3) :: vconst=(/.true.,.true.,.true./)           ! (film,drop,jet) if true calculation will be performed with constant velocity profiles
 logical :: shift=.false., vgmin=.false.!.true.


    contains

    subroutine hydro_film(eta_liq, eta_gas, rho_g, rho_l, alpha, mu_l, mu_g, R_c, q_g, q_l,delta, ul, ug, int_ul,int_ug)
        implicit none

 real*8, intent(in) :: eta_liq(nnl), eta_gas(nng)
 real*8, intent(in) :: rho_g, rho_l
 real*8, intent(in) :: mu_l, mu_g
 real*8, intent(in) :: q_l, q_g
 real*8, intent(in) :: R_c
 real*8, intent(in) :: alpha                       ! gravity flow angle
 real*8 :: r(NNl+NNg)                 ! radial coordinates

 real*8 :: delta                       ! liquid film thickness calculated by Brent's algorithm
 real*8 :: eps_l,eps_g                 ! pressure drop term of liquid and gas side
 real*8 :: C_l1,C_l2,C_g2!,C_turb1,C_turb2,C_turb3,C_turb4              ! integration constants of liquid and gas flow
 real*8 :: ul(NNl),ug(NNg)             ! liquid and gas velocity profiles
 real*8 :: int_ul,int_ug               ! integral mean velocities
 integer :: i                             !Laufvariable
 real*8 :: minimum=1D-5                    ! on input: domain, so that F(min)<0 & F(max)>0
 real*8 :: maximum                    ! on output: max is the best approximation to zero of F
! real*8 :: errabs=0D0                  ! absolute error of function result with estimated zero
 real*8 :: errrel=1D-24                 ! accepted change between two successive approximations
 INTEGER :: MAXFN=1000                 ! on input: max iterations, on output: needed iterations



maximum=1D-3
minimum=1D-5
call brent(mass_balance_film,minimum,maximum,errrel,MAXFN)

delta=maximum

eps_l=-g*rho_l*sin(alpha)

eps_g=(8*mu_l*q_l+((-4*log(R_c-delta)+4*log(R_c)+3)*delta**4+(16*R_c*log(R_c-delta)-16*R_c*log(R_c)-12*R_c)*delta**3+(-24*R_c**2*&
log(R_c-delta)+24*R_c**2*log(R_c)+14*R_c**2)*delta**2+(16*R_c**3*log(R_c-delta)-16*R_c**3*log(R_c)-4*R_c**3)*delta-4*R_c**4*&
log(R_c-delta)+4*R_c**4*log(R_c))*eps_l*pi)/(((4*log(R_c-delta)-4*log(R_c)-2)*delta**4+(-16*R_c*log(R_c-delta)+16*R_c*log(R_c)+&
8*R_c)*delta**3+(24*R_c**2*log(R_c-delta)-24*R_c**2*log(R_c)-10*R_c**2)*delta**2+(-16*R_c**3*log(R_c-delta)+16*R_c**3*log(R_c)+&
4*R_c**3)*delta+4*R_c**4*log(R_c-delta)-4*R_c**4*log(R_c))*pi)


!eps_g=-eps_l+2*mu_l/(R_c-delta)**2*(-q_l/(2*pi)+eps_l/(16*mu_l)*(R_c**4-(R_c-delta)**4)&
!    -eps_l/(8*mu_l)*R_c**2*(R_c**2-(R_c-delta)**2))/((R_c**2/2*log(R_c)-R_c**2/4)-((R_c-delta)&
!    **2/2*log(R_c-delta)-(R_c-delta)**2/4)-log(R_c)/2*(R_c**2-(R_c-delta)**2))


C_l1=(eps_l-eps_g)/2*(R_c-delta)**2
C_l2=1/mu_l*(eps_l/4*R_c**2-C_l1*log(R_c))

!((eps_l-eps_g)/2*(R_c-delta)^2)/mu_l*log(R_c-delta)+(eps_g/(4*mu_g)-eps_l/(4*mu_l))*(R_c-delta)^2+(1/mu_l*(eps_l/4*R_c^2-((eps_l-eps_g)/2*(R_c-delta)^2)*log(R_c)))


do i=1,NNl

    r(i)=eta_liq(i)*delta+R_c-delta
    ul(i)=C_l1/mu_l*log(r(i))-eps_l/(4*mu_l)*r(i)**2+C_l2

end do

int_ul=2*pi*(((delta**4-4*R_c*delta**3+6*R_c**2*delta**2-4*R_c**3*delta+R_c**4)*eps_l)/(16*mu_l)-(R_c**4*eps_l)/(16*mu_l)-&
(C_l1*log(R_c-delta)*delta**2-2*C_l1*R_c*log(R_c-delta)*delta+C_l1*R_c**2*log(R_c-delta))/(2*mu_l)+(C_l1*delta**2-2*C_l1*R_c*&
delta+C_l1*R_c**2)/(4*mu_l)+(C_l1*R_c**2*log(R_c))/(2*mu_l)-(C_l1*R_c**2)/(4*mu_l)-(C_l2*delta**2-2*C_l2*R_c*delta+C_l2*R_c**2)&
/2+(C_l2*R_c**2)/2)


call integration(r(1:NNl)*ul(1:nnl),r(1:NNl),int_ul)
!call integration(ul(1:nnl),r(1:NNl),int_ul)

int_ul=C_l1/mu_l*((R_c*log(R_c)-R_c)-((R_c-delta)*log(R_c-delta)-(R_c-delta)))&
        -eps_l/(12*mu_l)*(R_c**3-(R_c-delta)**3)+C_l2*(R_c-(R_c-delta))

if (turbulent) then
    u_tau = sqrt(1.D0/(8.D0*200.D0**0.25D0)*((q_g/(pi*(R_c-delta)**2))**1.75D0)*(mu_g/rho_g)**0.25D0*(R_c-delta)**(-0.25D0))
    C_mu_turb = B_*rho_g*(R_c-delta)*u_tau*Kappa/(2*A_r)
    C_g2=C_l1/mu_l*log(R_c-delta)-eps_l/(4D0*mu_l)*(R_c-delta)**2D0+C_l2-(((R_c-delta)**2*eps_g*log(abs((4*C_mu_turb*(R_c-delta)&
        **2-(R_c-delta)**2*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))-C_mu_turb*(R_c-delta)**2)/(4*C_mu_turb*(R_c-delta)**2+(R_c-delta)&
        **2*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))-C_mu_turb*(R_c-delta)**2))))/(4*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))))

if (vgmin) then
    C_g2=-3D-1-((R_c-delta)**2*eps_g*log(abs((4*C_mu_turb*(R_c-delta)**2-(R_c-delta)**2*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))-&
                C_mu_turb*(R_c-delta)**2)/(4*C_mu_turb*(R_c-delta)**2+(R_c-delta)**2*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))-&
                C_mu_turb*(R_c-delta)**2))))/(4*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb)))
end if


    do i=1,NNg

        r(i)=eta_gas(i)*(R_c-delta)
        ug(i)=((R_c-delta)**2*eps_g*log(abs((4*C_mu_turb*r(i)**2-(R_c-delta)**2*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))-&
                C_mu_turb*(R_c-delta)**2)/(4*C_mu_turb*r(i)**2+(R_c-delta)**2*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))-&
                C_mu_turb*(R_c-delta)**2))))/(4*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb)))+C_g2

!        (R_c-delta)**2*eps_g/(2*sqrt(C_mu_turb*(9*C_mu_turb+8*mu_g)))*datanh(sqrt(C_mu_turb)&
!                *(4*r(i)**2-(R_c-delta)**2)/((R_c-delta)**2*sqrt(9*C_mu_turb+8*mu_g)))+C_g2

    end do




else

    C_g2=C_l1/mu_l*log(R_c-delta)+(eps_g/(4*mu_g)-eps_l/(4*mu_l))*(R_c-delta)**2+C_l2

    ug(1)=-eps_g/(4*mu_g)*(R_c-delta)**2+C_g2

    do i=1,NNg

        r(i)=eta_gas(i)*(R_c-delta) !
        ug(i)=-eps_g/(4*mu_g)*r(i)**2+C_g2

    end do
    int_ug=(12*C_g2*(R_c-delta)*mu_g-(R_c-delta)**3*eps_g)/(12*mu_g)!maxval(abs(ug))
end if

if (vconst(1)) then
    ug=q_g/(pi*(R_c-delta)**2)
    ul=-q_l/(pi*(R_c**2-(R_c-delta)**2))
end if

!call integration(2D0*pi*ul*(eta_liq*delta+R_c-delta),(eta_liq*delta+R_c-delta),int_ul)

!call integration(2D0*pi*ug*(eta_gas*(R_c-delta)),(eta_gas*(R_c-delta)),int_ug)

call integration(2*pi*r(1:nng)*ug,r(1:NNg),int_ug)
int_ug=(C_g2/2D0*(R_c-delta)**2D0-(R_c-delta)**4D0*eps_g/(16*mu_g))*2D0*pi
!int_ug=1d0-int_ug/q_g

int_ul=maxval(abs(ul))!abs(int_ul)  abs(sum(ul)/size(ul))!
int_ug=maxval(abs(ug))!abs(int_ug)  abs(sum(ug)/size(ug))!
ul(:)=ul(:)/abs(int_ul)
if (shift) ug=ug+abs(minval(ug))
ug(:)=ug(:)/abs(int_ug)

    end subroutine hydro_film


    subroutine hydro_jet(eta_liq, eta_gas, rho_g, rho_l, mu_l, mu_g, R_c, q_g, q_l,delta, ul, ug, int_ul,int_ug)
        implicit none

 !input variables
! integer, intent(in) :: NNg
 !integer, intent(in) :: NNl
 real*8, intent(in) :: eta_liq(:), eta_gas(:)
 real*8, intent(in) :: rho_g, rho_l
 real*8 :: mu_l, mu_g
 real*8, intent(in) :: q_l, q_g
 real*8, intent(in) :: R_c
 real*8 :: r(NNl+NNg)                 ! radial coordinates

 real*8, intent(out) :: delta                       ! liquid film thickness calculated by Brent's algorithm
 real*8 :: eps_l,eps_g                 ! pressure drop term of liquid and gas side
 real*8 :: C_l1,C_l2,C_g2,C_g1!,C_turb1,C_turb2,C_turb3,C_turb4              ! integration constants of liquid and gas flow
 real*8, intent(out) :: ul(NNl),ug(NNg)             ! liquid and gas velocity profiles
 real*8, intent(out) :: int_ul,int_ug               ! integral mean velocities
 integer :: i                             !Laufvariable
 real*8 :: minimum=1D-4                    ! on input: domain, so that F(min)<0 & F(max)>0
 real*8 :: maximum                    ! on output: max is the best approximation to zero of F
 real*8 :: errrel=1D-24                 ! accepted change between two succesive approximations
 INTEGER :: MAXFN=1000                 ! on input: max iterations, on output: needed iterations

maximum=5D-3
minimum=1D-4
call brent(mass_balance_jet,minimum,maximum,errrel,MAXFN)

delta=maximum

eps_l=-g*rho_l

eps_g=-(8*mu_g*mu_l*q_l+((4*log(R_c)*delta**4-4*delta**4*log(delta))*eps_l*mu_l+delta**4*eps_l*mu_g)*pi)/((4*delta**4*log(delta)+&
        (-4*log(R_c)-2)*delta**4+2*R_c**2*delta**2)*mu_l*pi)



C_g1=(-eps_l+eps_g)/2*(delta)**2
C_g2=eps_g*R_c**2/(4*mu_g)-C_g1/mu_g*log(R_c)
C_l1=0
C_l2=C_g1/mu_g*log(delta)-eps_g/(4*mu_g)*(delta)**2+C_g2+eps_l/(4*mu_l)*(delta)**2


if (turbulent) then
    u_tau = sqrt(1.D0/(8.D0*200.D0**0.25D0)*((q_g/(pi*(R_c-delta)**2))**1.75D0)*(mu_g/rho_g)**0.25D0*(R_c-delta)**(-0.25D0))
    C_mu_turb =B_*rho_g*(R_c-delta)*u_tau*Kappa/(2*A_r)

    eps_g=((-64*C_mu_turb*R_c**2*mu_g**2-136*C_mu_turb**2*R_c*mu_g-72*C_mu_turb**3)*mu_l*q_l+(((-8*C_mu_turb*R_c**2*delta**4*eps_l&
    *mu_g-9*C_mu_turb**2*R_c*delta**4*eps_l)*log(abs(-R_c**5*mu_g+2*C_mu_turb*delta**4-C_mu_turb*R_c**2*delta**2-C_mu_turb*R_c**4)&
    )+(8*C_mu_turb*R_c**2*delta**4*eps_l*mu_g+9*C_mu_turb**2*R_c*delta**4*eps_l)*log(abs(-R_c**5*mu_g-C_mu_turb*R_c**4-C_mu_turb*&
    R_c**2*R_c**2+2*C_mu_turb*R_c**4))+(32*C_mu_turb*R_c**2*delta**4*log(delta)-32*C_mu_turb*log(R_c)*R_c**2*delta**4)*eps_l*mu_g+(&
    36*C_mu_turb**2*R_c*delta**4*log(delta)-36*C_mu_turb**2*log(R_c)*R_c*delta**4)*eps_l)*mu_l-8*C_mu_turb*R_c**2*delta**4*eps_l*&
    mu_g**2-17*C_mu_turb**2*R_c*delta**4*eps_l*mu_g-9*C_mu_turb**3*delta**4*eps_l)*pi+sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb**2)*&
    (C_mu_turb*R_c*delta**4*eps_l*log(abs(-(sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb**2)*(2*R_c**7*mu_g-4*C_mu_turb*R_c**2*delta**4+2&
    *C_mu_turb*R_c**4*delta**2+2*C_mu_turb*R_c**6)+(8*C_mu_turb*R_c**5*delta**2-2*C_mu_turb*R_c**7)*mu_g-16*C_mu_turb**2*delta**6+&
    12*C_mu_turb**2*R_c**2*delta**4+6*C_mu_turb**2*R_c**4*delta**2-2*C_mu_turb**2*R_c**6)/(sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb&
    **2)*(2*R_c**7*mu_g+12*C_mu_turb*R_c**2*delta**4-6*C_mu_turb*R_c**4*delta**2+3*C_mu_turb*R_c**6)+(24*C_mu_turb*R_c**5*delta**2&
    -6*C_mu_turb*R_c**7)*mu_g+16*C_mu_turb**2*delta**6-12*C_mu_turb**2*R_c**2*delta**4+30*C_mu_turb**2*R_c**4*delta**2-7*C_mu_turb&
    **2*R_c**6)))-C_mu_turb*R_c*delta**4*eps_l*log(abs(-(sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb**2)*(2*R_c**7*mu_g+2*C_mu_turb*R_c&
    **6+2*C_mu_turb*R_c**2*R_c**4-4*C_mu_turb*R_c**4*R_c**2)+(8*C_mu_turb*R_c**2*R_c**5-2*C_mu_turb*R_c**7)*mu_g-2*C_mu_turb**2*&
    R_c**6+6*C_mu_turb**2*R_c**2*R_c**4+12*C_mu_turb**2*R_c**4*R_c**2-16*C_mu_turb**2*R_c**6)/(sqrt(8*C_mu_turb*R_c*mu_g+9*&
    C_mu_turb**2)*(2*R_c**7*mu_g+3*C_mu_turb*R_c**6-6*C_mu_turb*R_c**2*R_c**4+12*C_mu_turb*R_c**4*R_c**2)+(24*C_mu_turb*R_c**2*R_c&
    **5-6*C_mu_turb*R_c**7)*mu_g-7*C_mu_turb**2*R_c**6+30*C_mu_turb**2*R_c**2*R_c**4-12*C_mu_turb**2*R_c**4*R_c**2+16*C_mu_turb**2&
    *R_c**6))))*mu_l*pi)/(sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb**2)*((2*R_c**4*delta**2*mu_g+C_mu_turb*R_c*delta**4+2*C_mu_turb*&
    R_c**3*delta**2)*log(abs(-(sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb**2)*(2*R_c**7*mu_g-4*C_mu_turb*R_c**2*delta**4+2*C_mu_turb*&
    R_c**4*delta**2+2*C_mu_turb*R_c**6)+(8*C_mu_turb*R_c**5*delta**2-2*C_mu_turb*R_c**7)*mu_g-16*C_mu_turb**2*delta**6+12*&
    C_mu_turb**2*R_c**2*delta**4+6*C_mu_turb**2*R_c**4*delta**2-2*C_mu_turb**2*R_c**6)/(sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb**2)*&
    (2*R_c**7*mu_g+12*C_mu_turb*R_c**2*delta**4-6*C_mu_turb*R_c**4*delta**2+3*C_mu_turb*R_c**6)+(24*C_mu_turb*R_c**5*delta**2-6*&
    C_mu_turb*R_c**7)*mu_g+16*C_mu_turb**2*delta**6-12*C_mu_turb**2*R_c**2*delta**4+30*C_mu_turb**2*R_c**4*delta**2-7*C_mu_turb**&
    2*R_c**6)))+(-2*R_c**4*delta**2*mu_g-C_mu_turb*R_c*delta**4-2*C_mu_turb*R_c**3*delta**2)*log(abs(-(sqrt(8*C_mu_turb*R_c*mu_g+9&
    *C_mu_turb**2)*(2*R_c**7*mu_g+2*C_mu_turb*R_c**6+2*C_mu_turb*R_c**2*R_c**4-4*C_mu_turb*R_c**4*R_c**2)+(8*C_mu_turb*R_c**2*&
    R_c**5-2*C_mu_turb*R_c**7)*mu_g-2*C_mu_turb**2*R_c**6+6*C_mu_turb**2*R_c**2*R_c**4+12*C_mu_turb**2*R_c**4*R_c**2-16*&
    C_mu_turb**2*R_c**6)/(sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb**2)*(2*R_c**7*mu_g+3*C_mu_turb*R_c**6-6*C_mu_turb*R_c**2*R_c**4+12&
    *C_mu_turb*R_c**4*R_c**2)+(24*C_mu_turb*R_c**2*R_c**5-6*C_mu_turb*R_c**7)*mu_g-7*C_mu_turb**2*R_c**6+30*C_mu_turb**2*R_c**2*&
    R_c**4-12*C_mu_turb**2*R_c**4*R_c**2+16*C_mu_turb**2*R_c**6))))*mu_l*pi+((-8*C_mu_turb*R_c**2*delta**4*mu_g-9*C_mu_turb**2*R_c&
    *delta**4)*log(abs(-R_c**5*mu_g+2*C_mu_turb*delta**4-C_mu_turb*R_c**2*delta**2-C_mu_turb*R_c**4))+(8*C_mu_turb*R_c**2*delta**4&
    *mu_g+9*C_mu_turb**2*R_c*delta**4)*log(abs(-R_c**5*mu_g-C_mu_turb*R_c**4-C_mu_turb*R_c**2*R_c**2+2*C_mu_turb*R_c**4))+(32*&
    C_mu_turb*R_c**2*delta**4*log(delta)-32*C_mu_turb*log(R_c)*R_c**2*delta**4)*mu_g+36*C_mu_turb**2*R_c*delta**4*log(delta)-36*&
    C_mu_turb**2*log(R_c)*R_c*delta**4)*mu_l*pi)

    C_g1=(-eps_l+eps_g)/2*(delta)**2

    C_g2=-(C_g1*((C_mu_turb*R_c*log(abs((4*C_mu_turb*R_c**2-R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2)/(4*&
    C_mu_turb*R_c**2+R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2))))/((4*R_c*mu_g+4*C_mu_turb)*sqrt(C_mu_turb&
    *(8*R_c*mu_g+9*C_mu_turb)))-(R_c*log(abs(2*C_mu_turb*R_c**4-C_mu_turb*R_c**2*R_c**2-R_c**5*mu_g-C_mu_turb*R_c**4)))/(4*R_c*&
    mu_g+4*C_mu_turb)+(2*R_c*log(R_c))/(2*R_c*mu_g+2*C_mu_turb))+(R_c**3*eps_g*log(abs((4*C_mu_turb*R_c**2-R_c**2*sqrt(C_mu_turb*&
    (8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2)/(4*C_mu_turb*R_c**2+R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*&
    R_c**2))))/(4*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))))

    c_l2=C_g1*((C_mu_turb*R_c*log(abs((4*C_mu_turb*delta**2-R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2)/(4*&
    C_mu_turb*delta**2+R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2))))/((4*R_c*mu_g+4*C_mu_turb)*sqrt(&
    C_mu_turb*(8*R_c*mu_g+9*C_mu_turb)))-(R_c*log(abs(2*C_mu_turb*delta**4-C_mu_turb*R_c**2*delta**2-R_c**5*mu_g-C_mu_turb*R_c**4)&
    ))/(4*R_c*mu_g+4*C_mu_turb)+(2*R_c*log(delta))/(2*R_c*mu_g+2*C_mu_turb))+(R_c**3*eps_g*log(abs((4*C_mu_turb*delta**2-R_c**2*&
    sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2)/(4*C_mu_turb*delta**2+R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb&
    ))-C_mu_turb*R_c**2))))/(4*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb)))+C_g2+eps_l/(4*mu_l)*(delta)**2



    do i=1,NNg

        r(i)=eta_gas(i)*(R_c-delta)+delta
        ug(i)=C_g1*((C_mu_turb*R_c*log(abs((4*C_mu_turb*r(i)**2-R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2)/&
        (4*C_mu_turb*r(i)**2+R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2))))/((4*R_c*mu_g+4*C_mu_turb)*sqrt&
        (C_mu_turb*(8*R_c*mu_g+9*C_mu_turb)))-(R_c*log(abs(2*C_mu_turb*r(i)**4-C_mu_turb*R_c**2*r(i)**2-R_c**5*mu_g-C_mu_turb*R_c&
        **4)))/(4*R_c*mu_g+4*C_mu_turb)+(2*R_c*log(r(i)))/(2*R_c*mu_g+2*C_mu_turb))+(R_c**3*eps_g*log(abs((4*C_mu_turb*r(i)**2-R_c&
        **2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2)/(4*C_mu_turb*r(i)**2+R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*&
        C_mu_turb))-C_mu_turb*R_c**2))))/(4*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb)))+C_g2

    end do

call integration(ug,r(1:nng),int_ug)


else


    do i=1,NNg

        r(i)=eta_gas(i)*(R_c-delta)+delta !
        ug(i)=C_g1/mu_g*log(r(i))-eps_g/(4*mu_g)*r(i)**2+C_g2

    end do
    int_ug=-((12*C_g2*delta-12*C_g2*R_c)*mu_g+(R_c**3-delta**3)*eps_g+12*C_g1*delta*log(delta)-12*C_g1*delta-12*C_g1*R_c*&
                log(R_c)+12*C_g1*R_c)/(12*mu_g)
end if

do i=1,NNl

    r(i)=eta_liq(i)*delta
    ul(i)=-eps_l/(4*mu_l)*r(i)**2+C_l2

end do


if (vconst(3)) then
!!!! const velocity profiles !!!
    ug=q_g/(pi*(R_c**2-delta**2))
    ul=-q_l/(pi*delta**2)
end if


call integration(2*pi*r(1:nnl)*ul,r(1:nnl),int_ul)
int_ul=maxval(abs(ul))  !abs(sum(ul)/size(ul))
int_ug=maxval(abs(ug))  !abs(sum(ug)/size(ug))
ul(:)=ul(:)/abs(int_ul)
if (shift) ug=ug+abs(minval(ug))
ug(:)=ug(:)/abs(int_ug)

    end subroutine hydro_jet



    subroutine hydro_drop(eta_liq, eta_gas, rho_g, rho_l, sigma, mu_l, mu_g, R_c, q_g, q_l,delta, ul, ug, int_ul,int_ug)
        implicit none

 real*8, intent(in) :: eta_liq(:), eta_gas(:)
 real*8, intent(in) :: rho_g, rho_l
 real*8, intent(in) :: mu_l, mu_g
 real*8, intent(in) :: q_l, q_g
 real*8, intent(in) :: R_c
 real*8, intent(in) :: sigma                       ! surface tension
 real*8 :: r(NNl+NNg)                 ! radial coordinates
 real*8 :: d_Drop

 real*8, intent(out) :: delta                       ! liquid film thickness calculated by Brent's algorithm
 real*8 :: eps_g                 ! pressure drop term of liquid and gas side
 real*8 :: C_g1,C_g2!,C_turb1,C_turb2,C_turb3,C_turb4              ! integration constants of liquid and gas flow
 real*8, intent(out) :: ul(NNl),ug(NNg)             ! liquid and gas velocity profiles
 real*8, intent(out) :: int_ul,int_ug               ! integral mean velocities
 integer :: i                             !Laufvariable
 real*8 :: minimum=1D-1                    ! on input: domain, so that F(min)<0 & F(max)>0
 real*8 :: maximum=5D7                    ! on output: max is the best approximation to zero of F
 real*8 :: errrel=1D-24                 ! accepted change between two succesive approximations
 INTEGER :: MAXFN=1000                 ! on input: max iterations, on output: needed iterations

d_Drop=sqrt(sigma/((rho_l-rho_g)*9.81))
ul=-8D-1*5.25D-1*(2*R_c/d_Drop)**(25D-2)*sqrt((rho_l-rho_g)*g*d_Drop/rho_g)+q_g/(pi*(R_c**2))!-0.075*(2*R_c/d_Drop)**(25D-2)*sqrt((rho_l-rho_g)/rho_g)
delta=q_l*6D0/(d_Drop*abs(ul(1))*pi*2D0)

eps_g=-((8*delta**2*log(delta)-8*log(R_c-delta)*delta**2+(16*R_c*log(R_c-delta)-16*R_c*log(R_c)-8*R_c)*delta-8*R_c**2*&
        log(R_c-delta)+8*R_c**2*log(R_c)+4*R_c**2)*mu_g*pi*ul(1)+(8*log(delta)-8*log(R_c))*mu_g*q_g)/(((2*delta**4-4*R_c*delta&
        **3+4*R_c**2*delta**2-R_c**4)*log(delta)-2*log(R_c-delta)*delta**4+(4*R_c*log(R_c-delta)-2*R_c)*delta**3+(R_c**2-4*R_c&
        **2*log(R_c))*delta**2+(-4*R_c**3*log(R_c-delta)+4*R_c**3*log(R_c)+2*R_c**3)*delta+2*R_c**4*log(R_c-delta)-R_c**4*&
        log(R_c)-R_c**4)*pi)

int_ul=ul(1)


if (turbulent) then
    minimum=1D-1
    maximum=5D7
    call brent(mass_balance_drop,minimum,maximum,errrel,MAXFN)
    eps_g=maximum
    u_tau = sqrt(1.D0/(8.D0*200.D0**0.25D0)*((q_g/(pi*(R_c-delta)**2))**1.75D0)*(mu_g/rho_g)**0.25D0*(R_c-delta)**(-0.25D0))
    C_mu_turb = B_*rho_g*(R_c-delta)*u_tau*Kappa/(2*A_r)*R_c**(-5D0)

    C_g1=-((32*C_mu_turb*mu_g**2+68*C_mu_turb**2*R_c**4*mu_g+36*C_mu_turb**3*R_c**8)*ul(1)+sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c&
    **4)*(-eps_g*mu_g*log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*&
    C_mu_turb*delta**2-C_mu_turb*R_c**2)+(4*C_mu_turb*delta**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb*delta**2&
    -C_mu_turb*R_c**2)-(C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb*delta**2-C_mu_turb*R_c**2)))-&
    C_mu_turb*R_c**4*eps_g*log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*&
    C_mu_turb*delta**2-C_mu_turb*R_c**2)+(4*C_mu_turb*delta**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb*delta**2&
    -C_mu_turb*R_c**2)-(C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb*delta**2-C_mu_turb*R_c**2))))+&
    sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)*(eps_g*mu_g+C_mu_turb*R_c**4*eps_g)*log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2&
    *R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)-(C_mu_turb*R_c**2)/(sqrt(8*&
    C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)+(4*C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*&
    C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*C_mu_turb*R_c**2))))/(-C_mu_turb*R_c**2*sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)*&
    log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb*delta**2-&
    C_mu_turb*R_c**2)+(4*C_mu_turb*delta**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb*delta**2-C_mu_turb*R_c**2)-&
    (C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb*delta**2-C_mu_turb*R_c**2)))+C_mu_turb*R_c**2*&
    sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)*log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*&
    C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)-(C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-&
    C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)+(4*C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*&
    C_mu_turb*R_c**2)))+(8*C_mu_turb*log(abs(-mu_g+2*C_mu_turb*delta**4-C_mu_turb*R_c**2*delta**2-C_mu_turb*R_c**4))-8*C_mu_turb*&
    log(abs(-mu_g-C_mu_turb*R_c**4-C_mu_turb*R_c**2*R_c**2+2*C_mu_turb*R_c**4))-32*C_mu_turb*log(delta)+32*C_mu_turb*log(R_c))*&
    mu_g+9*C_mu_turb**2*R_c**4*log(abs(-mu_g+2*C_mu_turb*delta**4-C_mu_turb*R_c**2*delta**2-C_mu_turb*R_c**4))-9*C_mu_turb**2*R_c&
    **4*log(abs(-mu_g-C_mu_turb*R_c**4-C_mu_turb*R_c**2*R_c**2+2*C_mu_turb*R_c**4))-36*C_mu_turb**2*R_c**4*log(delta)+36*C_mu_turb&
    **2*log(R_c)*R_c**4)

    C_g2=-(eps_g*mu_g*log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-&
    C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)-(C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*&
    C_mu_turb*R_c**2)+(4*C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)))+&
    C_mu_turb*R_c**4*eps_g*log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-&
    C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)-(C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*&
    C_mu_turb*R_c**2)+(4*C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)))+&
    C_g1*C_mu_turb*R_c**2*log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-&
    C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)-(C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*&
    C_mu_turb*R_c**2)+(4*C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)))-&
    C_g1*log(abs(-mu_g-C_mu_turb*R_c**4-C_mu_turb*R_c**2*R_c**2+2*C_mu_turb*R_c**4))*sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+&
    4*C_g1*log(R_c)*sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4))/(4*mu_g*sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb&
    *R_c**4*sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4))



    do i=1,NNg

        r(i)=eta_gas(i)*(R_c-delta)+delta
        ug(i)=C_g1*((C_mu_turb*R_c**2*log(abs((4*C_mu_turb*r(i)**2-sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb*R_c**4))-C_mu_turb*R_c**2)&
        /(4*C_mu_turb*r(i)**2+sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb*R_c**4))-C_mu_turb*R_c**2))))/((4*mu_g+4*C_mu_turb*R_c**4)*&
        sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb*R_c**4)))-log(abs(2*C_mu_turb*r(i)**4-C_mu_turb*R_c**2*r(i)**2-mu_g-C_mu_turb*R_c&
        **4))/(4*mu_g+4*C_mu_turb*R_c**4)+(2*log(r(i)))/(2*mu_g+2*C_mu_turb*R_c**4))+(eps_g*log(abs((4*C_mu_turb*r(i)**2-&
        sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb*R_c**4))-C_mu_turb*R_c**2)/(4*C_mu_turb*r(i)**2+sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb*&
        R_c**4))-C_mu_turb*R_c**2))))/(4*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb*R_c**4)))+C_g2

    end do


else

    C_g1=(ul(1)+eps_g/(4*mu_g)*(delta**2-R_c**2))*mu_g/log(delta)/(1-log(R_c)/log(delta))
    C_g2=eps_g*R_c**2/(4*mu_g)-C_g1/mu_g*log(R_c)
    do i=1,NNg

        r(i)=eta_gas(i)*(R_c-delta)+delta !
        ug(i)=C_g1/mu_g*log(r(i))-eps_g/(4*mu_g)*r(i)**2+C_g2

    end do
    int_ug=-((12*C_g2*delta-12*C_g2*R_c)*mu_g+(R_c**3-delta**3)*eps_g+12*C_g1*delta*log(delta)-12*C_g1*delta-12*C_g1*R_c*&
                log(R_c)+12*C_g1*R_c)/(12*mu_g)
end if


if (vconst(2)) then
!!!! const velocity profiles !!!
    ug=q_g/(pi*(R_c**2-delta**2))
endif

int_ul=maxval(abs(ul))  !abs(sum(ul)/size(ul))
int_ug=maxval(abs(ug))  !abs(sum(ug)/size(ug))
ul(:)=ul(:)/abs(int_ul)
if (shift) ug=ug+abs(minval(ug))
ug(:)=ug(:)/abs(int_ug)

    end subroutine hydro_drop




real*8 FUNCTION mass_balance_film(X)
IMPLICIT NONE
real*8 :: int_ug,int_ul, u_tau              !integral of velocity. in mass_balance: 2*pi*int(ug*r)dr
real*8 :: X, eps_l,eps_g,C_l1,C_l2,C_g2,C_mu_turb             ! Variable for Brent's algorithm, equals delta
real*8 :: r_t(50001),ug_t(50001)
integer :: l
    eps_l=-g*rho_l*sin(alpha)
    eps_g=(8*mu_l*q_l+((-4*log(R_c-X)+4*log(R_c)+3)*X**4+(16*R_c*log(R_c-X)-16*R_c*log(R_c)-12*R_c)*X**3+(-24*R_c**2*&
log(R_c-X)+24*R_c**2*log(R_c)+14*R_c**2)*X**2+(16*R_c**3*log(R_c-X)-16*R_c**3*log(R_c)-4*R_c**3)*X-4*R_c**4*&
log(R_c-X)+4*R_c**4*log(R_c))*eps_l*pi)/(((4*log(R_c-X)-4*log(R_c)-2)*X**4+(-16*R_c*log(R_c-X)+16*R_c*log(R_c)+&
8*R_c)*X**3+(24*R_c**2*log(R_c-X)-24*R_c**2*log(R_c)-10*R_c**2)*X**2+(-16*R_c**3*log(R_c-X)+16*R_c**3*log(R_c)+&
4*R_c**3)*X+4*R_c**4*log(R_c-X)-4*R_c**4*log(R_c))*pi)

!    -eps_l+2D0*mu_l/(R_c-X)**2D0*(-q_l/(2D0*pi)+eps_l/(16D0*mu_l)*(R_c**4D0-(R_c-X)**4D0)&
!            -eps_l/(8D0*mu_l)*R_c**2D0*(R_c**2D0-(R_c-X)**2D0))/((R_c**2D0/2D0*log(R_c)-R_c**2D0/4D0)-((R_c-X)&
!            **2D0/2D0*log(R_c-X)-(R_c-X)**2D0/4D0)-log(R_c)/2D0*(R_c**2D0-(R_c-X)**2D0))

    C_l1=(eps_l-eps_g)/2D0*(R_c-X)**2D0
    C_l2=1D0/mu_l*(eps_l/4D0*R_c**2D0-C_l1*log(R_c))



    IF (turbulent) THEN !turbulent calculations following

        u_tau = sqrt(1.D0/(8.D0*200.D0**0.25D0)*((q_g/(pi*(R_c-X)**2))**1.75D0)*(mu_g/rho_g)**0.25D0*(R_c-X)**(-0.25D0))

        C_mu_turb = B_*rho_g*(R_c-X)*u_tau*Kappa/(2D0*A_r)

        C_g2=C_l1/mu_l*log(R_c-X)-eps_l/(4D0*mu_l)*(R_c-X)**2D0+C_l2-(((R_c-X)**2*eps_g*log(abs((4*C_mu_turb*(R_c-X)**2-(R_c-X)**2&
                *sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))-&
                C_mu_turb*(R_c-X)**2)/(4*C_mu_turb*(R_c-X)**2+(R_c-X)**2*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))-&
                C_mu_turb*(R_c-X)**2))))/(4*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))))

if (vgmin) then
        C_g2=-3D-1-((R_c-X)**2*eps_g*log(abs((4*C_mu_turb*(R_c-X)**2-(R_c-X)**2*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))-&
               C_mu_turb*(R_c-X)**2)/(4*C_mu_turb*(R_c-X)**2+(R_c-X)**2*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))-&
               C_mu_turb*(R_c-X)**2))))/(4*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb)))
end if


    r_t(1)=0
    ug_t(1)=((R_c-X)**2*eps_g*log(abs((4*C_mu_turb*r_t(1)**2-(R_c-X)**2*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))-&
        C_mu_turb*(R_c-X)**2)/(4*C_mu_turb*r_t(1)**2+(R_c-X)**2*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))-&
        C_mu_turb*(R_c-X)**2))))/(4*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb)))+C_g2

    do l=2,50001

        r_t(l)=(dfloat(l)-1)/50000D0*(R_c-X)

!        ug_t(l)=(R_c-X)**2D0*eps_g/(2D0*sqrt(C_mu_turb*(9D0*C_mu_turb+8D0*mu_g)))*datanh(abs(sqrt(C_mu_turb)&
!                *(4D0*r_t(l)**2D0-(R_c-X)**2D0)/((R_c-X)**2D0*sqrt(9D0*C_mu_turb+8D0*mu_g))))+C_g2

        ug_t(l)=((R_c-X)**2*eps_g*log(abs((4*C_mu_turb*r_t(l)**2-(R_c-X)**2*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))-&
                C_mu_turb*(R_c-X)**2)/(4*C_mu_turb*r_t(l)**2+(R_c-X)**2*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb))-&
                C_mu_turb*(R_c-X)**2))))/(4*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb)))+C_g2

    end do


    call integration(2D0*pi*ug_t*r_t,r_t,int_ug)

    mass_balance_film=1D0-int_ug/q_g
    ELSE !laminar calculation following

    C_g2=C_l1/mu_l*log(R_c-X)+(eps_g/(4*mu_g)-eps_l/(4*mu_l))*(R_c-X)**2+C_l2

!    int_ug=(C_g2/2D0*(R_c-X)**2D0-(R_c-X)**4D0*eps_g/(16*mu_g))*2D0*pi
    int_ug=(((8*C_g2*x**2-16*C_g2*R_c*x+8*C_g2*R_c**2)*mu_g+(-x**4+4*R_c*x**3-6*R_c**2*x**2+4*R_c**3*x-R_c**4)*eps_g)*pi)/(8*mu_g)

    int_ul=2*pi*(((X**4-4*R_c*X**3+6*R_c**2*X**2-4*R_c**3*X+R_c**4)*eps_l)/(16*mu_l)-(R_c**4*eps_l)/(16*mu_l)-&
(C_l1*log(R_c-X)*X**2-2*C_l1*R_c*log(R_c-X)*X+C_l1*R_c**2*log(R_c-X))/(2*mu_l)+(C_l1*X**2-2*C_l1*R_c*&
X+C_l1*R_c**2)/(4*mu_l)+(C_l1*R_c**2*log(R_c))/(2*mu_l)-(C_l1*R_c**2)/(4*mu_l)-(C_l2*X**2-2*C_l2*R_c*X+C_l2*R_c**2)&
/2+(C_l2*R_c**2)/2)

    mass_balance_film=1D0-int_ug/q_g!int_ul/q_l+1D0

    END IF !End of calculations





END FUNCTION mass_balance_film

real*8 FUNCTION mass_balance_jet(X)
IMPLICIT NONE
real*8 :: int_ug, u_tau              !integral of velocity. in mass_balance: 2*pi*int(ug*r)dr
real*8 :: X, eps_l,eps_g,C_l1,C_l2,C_g1,C_g2,C_mu_turb             ! Variable for Brent's algorithm, equals delta
real*8 :: r(1000),ug(1000)
integer :: i

eps_l=-g*rho_l

eps_g=-(8*mu_g*mu_l*q_l+((4*log(R_c)*X**4-4*X**4*log(X))*eps_l*mu_l+X**4*eps_l*mu_g)*pi)/((4*X**4*log(X)+&
        (-4*log(R_c)-2)*X**4+2*R_c**2*X**2)*mu_l*pi)!-(8*mu_g*mu_l*q_g+((4*log(R_c)*X**4-4*X**4*log(X))*eps_l*mu_l+X**4*eps_l*mu_g)*pi)/((4*X**4*log(X)+(-4*log(R_c)-2)*X**4+2*&
            !R_c**2*X**2)*mu_l*pi)


C_g1=(-eps_l+eps_g)/2*(X)**2
C_g2=eps_g*R_c**2/(4*mu_g)-C_g1/mu_g*log(R_c)
C_l1=0
C_l2=C_g1/mu_g*log(X)-eps_g/(4*mu_g)*(X)**2+C_g2+eps_l/(4*mu_l)*(X)**2


if (turbulent) then
    u_tau = sqrt(1.D0/(8.D0*200.D0**0.25D0)*((q_g/(pi*(R_c-X)**2))**1.75D0)*(mu_g/rho_g)**0.25D0*(R_c-X)**(-0.25D0))
    C_mu_turb = B_*rho_g*(R_c-X)*u_tau*Kappa/(2*A_r)

    eps_g=((-64*C_mu_turb*R_c**2*mu_g**2-136*C_mu_turb**2*R_c*mu_g-72*C_mu_turb**3)*mu_l*q_l+(((-8*C_mu_turb*R_c**2*X**4*eps_l&
    *mu_g-9*C_mu_turb**2*R_c*X**4*eps_l)*log(abs(-R_c**5*mu_g+2*C_mu_turb*X**4-C_mu_turb*R_c**2*X**2-C_mu_turb*R_c**4)&
    )+(8*C_mu_turb*R_c**2*X**4*eps_l*mu_g+9*C_mu_turb**2*R_c*X**4*eps_l)*log(abs(-R_c**5*mu_g-C_mu_turb*R_c**4-C_mu_turb*&
    R_c**2*R_c**2+2*C_mu_turb*R_c**4))+(32*C_mu_turb*R_c**2*X**4*log(X)-32*C_mu_turb*log(R_c)*R_c**2*X**4)*eps_l*mu_g+(&
    36*C_mu_turb**2*R_c*X**4*log(X)-36*C_mu_turb**2*log(R_c)*R_c*X**4)*eps_l)*mu_l-8*C_mu_turb*R_c**2*X**4*eps_l*&
    mu_g**2-17*C_mu_turb**2*R_c*X**4*eps_l*mu_g-9*C_mu_turb**3*X**4*eps_l)*pi+sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb**2)*&
    (C_mu_turb*R_c*X**4*eps_l*log(abs(-(sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb**2)*(2*R_c**7*mu_g-4*C_mu_turb*R_c**2*X**4+2&
    *C_mu_turb*R_c**4*X**2+2*C_mu_turb*R_c**6)+(8*C_mu_turb*R_c**5*X**2-2*C_mu_turb*R_c**7)*mu_g-16*C_mu_turb**2*X**6+&
    12*C_mu_turb**2*R_c**2*X**4+6*C_mu_turb**2*R_c**4*X**2-2*C_mu_turb**2*R_c**6)/(sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb&
    **2)*(2*R_c**7*mu_g+12*C_mu_turb*R_c**2*X**4-6*C_mu_turb*R_c**4*X**2+3*C_mu_turb*R_c**6)+(24*C_mu_turb*R_c**5*X**2&
    -6*C_mu_turb*R_c**7)*mu_g+16*C_mu_turb**2*X**6-12*C_mu_turb**2*R_c**2*X**4+30*C_mu_turb**2*R_c**4*X**2-7*C_mu_turb&
    **2*R_c**6)))-C_mu_turb*R_c*X**4*eps_l*log(abs(-(sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb**2)*(2*R_c**7*mu_g+2*C_mu_turb*R_c&
    **6+2*C_mu_turb*R_c**2*R_c**4-4*C_mu_turb*R_c**4*R_c**2)+(8*C_mu_turb*R_c**2*R_c**5-2*C_mu_turb*R_c**7)*mu_g-2*C_mu_turb**2*&
    R_c**6+6*C_mu_turb**2*R_c**2*R_c**4+12*C_mu_turb**2*R_c**4*R_c**2-16*C_mu_turb**2*R_c**6)/(sqrt(8*C_mu_turb*R_c*mu_g+9*&
    C_mu_turb**2)*(2*R_c**7*mu_g+3*C_mu_turb*R_c**6-6*C_mu_turb*R_c**2*R_c**4+12*C_mu_turb*R_c**4*R_c**2)+(24*C_mu_turb*R_c**2*R_c&
    **5-6*C_mu_turb*R_c**7)*mu_g-7*C_mu_turb**2*R_c**6+30*C_mu_turb**2*R_c**2*R_c**4-12*C_mu_turb**2*R_c**4*R_c**2+16*C_mu_turb**2&
    *R_c**6))))*mu_l*pi)/(sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb**2)*((2*R_c**4*X**2*mu_g+C_mu_turb*R_c*X**4+2*C_mu_turb*&
    R_c**3*X**2)*log(abs(-(sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb**2)*(2*R_c**7*mu_g-4*C_mu_turb*R_c**2*X**4+2*C_mu_turb*&
    R_c**4*X**2+2*C_mu_turb*R_c**6)+(8*C_mu_turb*R_c**5*X**2-2*C_mu_turb*R_c**7)*mu_g-16*C_mu_turb**2*X**6+12*&
    C_mu_turb**2*R_c**2*X**4+6*C_mu_turb**2*R_c**4*X**2-2*C_mu_turb**2*R_c**6)/(sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb**2)*&
    (2*R_c**7*mu_g+12*C_mu_turb*R_c**2*X**4-6*C_mu_turb*R_c**4*X**2+3*C_mu_turb*R_c**6)+(24*C_mu_turb*R_c**5*X**2-6*&
    C_mu_turb*R_c**7)*mu_g+16*C_mu_turb**2*X**6-12*C_mu_turb**2*R_c**2*X**4+30*C_mu_turb**2*R_c**4*X**2-7*C_mu_turb**&
    2*R_c**6)))+(-2*R_c**4*X**2*mu_g-C_mu_turb*R_c*X**4-2*C_mu_turb*R_c**3*X**2)*log(abs(-(sqrt(8*C_mu_turb*R_c*mu_g+9&
    *C_mu_turb**2)*(2*R_c**7*mu_g+2*C_mu_turb*R_c**6+2*C_mu_turb*R_c**2*R_c**4-4*C_mu_turb*R_c**4*R_c**2)+(8*C_mu_turb*R_c**2*&
    R_c**5-2*C_mu_turb*R_c**7)*mu_g-2*C_mu_turb**2*R_c**6+6*C_mu_turb**2*R_c**2*R_c**4+12*C_mu_turb**2*R_c**4*R_c**2-16*&
    C_mu_turb**2*R_c**6)/(sqrt(8*C_mu_turb*R_c*mu_g+9*C_mu_turb**2)*(2*R_c**7*mu_g+3*C_mu_turb*R_c**6-6*C_mu_turb*R_c**2*R_c**4+12&
    *C_mu_turb*R_c**4*R_c**2)+(24*C_mu_turb*R_c**2*R_c**5-6*C_mu_turb*R_c**7)*mu_g-7*C_mu_turb**2*R_c**6+30*C_mu_turb**2*R_c**2*&
    R_c**4-12*C_mu_turb**2*R_c**4*R_c**2+16*C_mu_turb**2*R_c**6))))*mu_l*pi+((-8*C_mu_turb*R_c**2*X**4*mu_g-9*C_mu_turb**2*R_c&
    *X**4)*log(abs(-R_c**5*mu_g+2*C_mu_turb*X**4-C_mu_turb*R_c**2*X**2-C_mu_turb*R_c**4))+(8*C_mu_turb*R_c**2*X**4&
    *mu_g+9*C_mu_turb**2*R_c*X**4)*log(abs(-R_c**5*mu_g-C_mu_turb*R_c**4-C_mu_turb*R_c**2*R_c**2+2*C_mu_turb*R_c**4))+(32*&
    C_mu_turb*R_c**2*X**4*log(X)-32*C_mu_turb*log(R_c)*R_c**2*X**4)*mu_g+36*C_mu_turb**2*R_c*X**4*log(X)-36*&
    C_mu_turb**2*log(R_c)*R_c*X**4)*mu_l*pi)

    C_g1=(-eps_l+eps_g)/2*(X)**2

    C_g2=-(C_g1*((C_mu_turb*R_c*log(abs((4*C_mu_turb*R_c**2-R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2)/(4*&
    C_mu_turb*R_c**2+R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2))))/((4*R_c*mu_g+4*C_mu_turb)*sqrt(C_mu_turb&
    *(8*R_c*mu_g+9*C_mu_turb)))-(R_c*log(abs(2*C_mu_turb*R_c**4-C_mu_turb*R_c**2*R_c**2-R_c**5*mu_g-C_mu_turb*R_c**4)))/(4*R_c*&
    mu_g+4*C_mu_turb)+(2*R_c*log(R_c))/(2*R_c*mu_g+2*C_mu_turb))+(R_c**3*eps_g*log(abs((4*C_mu_turb*R_c**2-R_c**2*sqrt(C_mu_turb*&
    (8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2)/(4*C_mu_turb*R_c**2+R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*&
    R_c**2))))/(4*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))))

    c_l2=C_g1*((C_mu_turb*R_c*log(abs((4*C_mu_turb*X**2-R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2)/(4*&
    C_mu_turb*X**2+R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2))))/((4*R_c*mu_g+4*C_mu_turb)*sqrt(&
    C_mu_turb*(8*R_c*mu_g+9*C_mu_turb)))-(R_c*log(abs(2*C_mu_turb*X**4-C_mu_turb*R_c**2*X**2-R_c**5*mu_g-C_mu_turb*R_c**4)&
    ))/(4*R_c*mu_g+4*C_mu_turb)+(2*R_c*log(X))/(2*R_c*mu_g+2*C_mu_turb))+(R_c**3*eps_g*log(abs((4*C_mu_turb*X**2-R_c**2*&
    sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2)/(4*C_mu_turb*X**2+R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb&
    ))-C_mu_turb*R_c**2))))/(4*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb)))+C_g2+eps_l/(4*mu_l)*(X)**2



    do i=1,1000

        r(i)=(1D3-i)/1D3*(R_c-X)+X
        ug(i)=C_g1*((C_mu_turb*R_c*log(abs((4*C_mu_turb*r(i)**2-R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2)/&
        (4*C_mu_turb*r(i)**2+R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2))))/((4*R_c*mu_g+4*C_mu_turb)*sqrt&
        (C_mu_turb*(8*R_c*mu_g+9*C_mu_turb)))-(R_c*log(abs(2*C_mu_turb*r(i)**4-C_mu_turb*R_c**2*r(i)**2-R_c**5*mu_g-C_mu_turb*R_c&
        **4)))/(4*R_c*mu_g+4*C_mu_turb)+(2*R_c*log(r(i)))/(2*R_c*mu_g+2*C_mu_turb))+(R_c**3*eps_g*log(abs((4*C_mu_turb*r(i)**2-R_c&
        **2*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb))-C_mu_turb*R_c**2)/(4*C_mu_turb*r(i)**2+R_c**2*sqrt(C_mu_turb*(8*R_c*mu_g+9*&
        C_mu_turb))-C_mu_turb*R_c**2))))/(4*sqrt(C_mu_turb*(8*R_c*mu_g+9*C_mu_turb)))+C_g2

    end do


else


    do i=1,1000

        r(i)=(1D3-i)/1D3*(R_c-X)+X !
        ug(i)=C_g1/mu_g*log(r(i))-eps_g/(4*mu_g)*r(i)**2+C_g2

    end do

end if

call integration(2*pi*ug*r,r,int_ug)


mass_balance_jet=1-int_ug/q_g

END FUNCTION mass_balance_jet


real*8 FUNCTION mass_balance_drop(X)
IMPLICIT NONE
real*8 :: int_ug, int_ul, ul, delta, d_Drop              !integral of velocity. in mass_balance: 2*pi*int(ug*r)dr
real*8 :: X,C_g1,C_g2,C_mu_turb, u_tau             ! Variable for Brent's algorithm, equals delta
real*8 :: r(3000),ug(3000)
integer :: i
d_Drop=sqrt(sigma/((rho_l-rho_g)*9.81))
ul=-0.075*(R_c/d_Drop)**(25D-2)*sqrt((rho_l-rho_g)/rho_g)
delta=(q_l/(d_Drop/6D0)/abs(ul)/pi)!5D-1*(q_l/(pi/6*d_Drop**3)*pi*d_Drop**2/(pi*0.1))

!int_ul=ul*delta

!ul=ul/int_ul

    u_tau = sqrt(1.D0/(8.D0*200.D0**0.25D0)*((q_g/(pi*(R_c-delta)**2))**1.75D0)*(mu_g/rho_g)**0.25D0*(R_c-delta)**(-0.25D0))
    C_mu_turb = B_*rho_g*(R_c-delta)*u_tau*Kappa/(2*A_r)*R_c**(-5D0)

    C_g1=-((32*C_mu_turb*mu_g**2+68*C_mu_turb**2*R_c**4*mu_g+36*C_mu_turb**3*R_c**8)*ul+sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c&
    **4)*(-X*mu_g*log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*&
    C_mu_turb*delta**2-C_mu_turb*R_c**2)+(4*C_mu_turb*delta**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb*delta**2&
    -C_mu_turb*R_c**2)-(C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb*delta**2-C_mu_turb*R_c**2)))-&
    C_mu_turb*R_c**4*X*log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*&
    C_mu_turb*delta**2-C_mu_turb*R_c**2)+(4*C_mu_turb*delta**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb*delta**2&
    -C_mu_turb*R_c**2)-(C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb*delta**2-C_mu_turb*R_c**2))))+&
    sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)*(X*mu_g+C_mu_turb*R_c**4*X)*log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2&
    *R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)-(C_mu_turb*R_c**2)/(sqrt(8*&
    C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)+(4*C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*&
    C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*C_mu_turb*R_c**2))))/(-C_mu_turb*R_c**2*sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)*&
    log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb*delta**2-&
    C_mu_turb*R_c**2)+(4*C_mu_turb*delta**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb*delta**2-C_mu_turb*R_c**2)-&
    (C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb*delta**2-C_mu_turb*R_c**2)))+C_mu_turb*R_c**2*&
    sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)*log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*&
    C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)-(C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-&
    C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)+(4*C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*&
    C_mu_turb*R_c**2)))+(8*C_mu_turb*log(abs(-mu_g+2*C_mu_turb*delta**4-C_mu_turb*R_c**2*delta**2-C_mu_turb*R_c**4))-8*C_mu_turb*&
    log(abs(-mu_g-C_mu_turb*R_c**4-C_mu_turb*R_c**2*R_c**2+2*C_mu_turb*R_c**4))-32*C_mu_turb*log(delta)+32*C_mu_turb*log(R_c))*&
    mu_g+9*C_mu_turb**2*R_c**4*log(abs(-mu_g+2*C_mu_turb*delta**4-C_mu_turb*R_c**2*delta**2-C_mu_turb*R_c**4))-9*C_mu_turb**2*R_c&
    **4*log(abs(-mu_g-C_mu_turb*R_c**4-C_mu_turb*R_c**2*R_c**2+2*C_mu_turb*R_c**4))-36*C_mu_turb**2*R_c**4*log(delta)+36*C_mu_turb&
    **2*log(R_c)*R_c**4)

    C_g2=-(X*mu_g*log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-&
    C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)-(C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*&
    C_mu_turb*R_c**2)+(4*C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)))+&
    C_mu_turb*R_c**4*X*log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-&
    C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)-(C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*&
    C_mu_turb*R_c**2)+(4*C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)))+&
    C_g1*C_mu_turb*R_c**2*log(abs(-sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-&
    C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)-(C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*&
    C_mu_turb*R_c**2)+(4*C_mu_turb*R_c**2)/(sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)-C_mu_turb*R_c**2+4*C_mu_turb*R_c**2)))-&
    C_g1*log(abs(-mu_g-C_mu_turb*R_c**4-C_mu_turb*R_c**2*R_c**2+2*C_mu_turb*R_c**4))*sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+&
    4*C_g1*log(R_c)*sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4))/(4*mu_g*sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4)+4*C_mu_turb&
    *R_c**4*sqrt(8*C_mu_turb*mu_g+9*C_mu_turb**2*R_c**4))

    do i=1,3000

        r(i)=(3D3-i)/3D3*(R_c-delta)+delta
        ug(i)=C_g1*((C_mu_turb*R_c**2*log(abs((4*C_mu_turb*r(i)**2-sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb*R_c**4))-C_mu_turb*R_c**2)&
        /(4*C_mu_turb*r(i)**2+sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb*R_c**4))-C_mu_turb*R_c**2))))/((4*mu_g+4*C_mu_turb*R_c**4)*&
        sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb*R_c**4)))-log(abs(2*C_mu_turb*r(i)**4-C_mu_turb*R_c**2*r(i)**2-mu_g-C_mu_turb*R_c&
        **4))/(4*mu_g+4*C_mu_turb*R_c**4)+(2*log(r(i)))/(2*mu_g+2*C_mu_turb*R_c**4))+(X*log(abs((4*C_mu_turb*r(i)**2-&
        sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb*R_c**4))-C_mu_turb*R_c**2)/(4*C_mu_turb*r(i)**2+sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb*&
        R_c**4))-C_mu_turb*R_c**2))))/(4*sqrt(C_mu_turb*(8*mu_g+9*C_mu_turb*R_c**4)))+C_g2

    end do

call integration(2*pi*ug*r,r,int_ug)


mass_balance_drop=1-int_ug/q_g

END FUNCTION mass_balance_drop

end module hydrodynamics
