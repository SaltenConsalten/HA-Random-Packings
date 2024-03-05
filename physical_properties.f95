module phys_prop

USE global_variables

    implicit none
    !integer, private :: ncomp

    !integer, public :: NNg, MM, NNl, NN
!    real*8, allocatable :: Dist_coef(:)            ! Distribution coefficient

! common /mesh/ ncomp, nnl, nng, nn, mm

    contains

    subroutine physical_properties(ncomp,comp_order,T_0_l,T_0_g,p_0_g,x_0,y_0,rho_l, rho_g,conc_l,&
                        conc_g,D_l,D_g,lambda_l, lambda_g,cp_l,cp_g,mu_l, mu_g,Dist_coef,sigma,H_lv)

        implicit none
        !######################input####################
        integer, intent(in) :: ncomp
        real*8, intent(in) :: T_0_l,T_0_g
        real*8, intent(in) :: p_0_g
        real*8, intent(in) :: x_0(ncomp),y_0(ncomp)
        character*14, intent(in) :: comp_order(ncomp)


        real*8 :: M(ncomp)                              ! Array of Molecular weights
        real*8 :: xsi(ncomp)
        real*8 :: mu_g_par(ncomp,5), mu_l_par(ncomp,5)
        real*8 :: psi(ncomp)
        real*8 :: cp_id(ncomp)                          ! [J/(mol K)] ideal heat capacities for each component

        logical :: condensable(ncomp)        ! if true, gas phase component (=vapour) is
                                            ! condensable at ambient conditions


    real*8 :: V_mol(ncomp)                ! [cm3/mol] molar volume of component i at its normal boiling temperature

    integer :: z_ion(ncomp)                ! [-] charge of ion i
    real*8 :: cond(ncomp)
    real*8 :: H_par(ncomp,4)            ! [-] array with parameters for calculation of the Henry-coefficient
    real*8 :: cp_gas_par(ncomp,4)        ! parameter for calculating the heat capacity of pure gases
    real*8 :: T_crit(ncomp)                ! [K] critical temperature
    real*8 :: acentric(ncomp)            ! [-] Pitzer acentric factor listed in [1]
    real*8 :: dH_excess(ncomp)            ! excess enthalpiy/ absorption heat
    real*8 :: lambda_gas_par(ncomp,4)    ! parameter for calculating the heat conductivity of pure gases
    real*8 :: lambda_liq_par(ncomp,3)    ! parameter for calculating the heat conductivity of pure liquids


        !######################output####################
        real*8, intent(out) :: rho_l, rho_g             ! Densities
        real*8, intent(out) :: conc_l, conc_g           ! [kmol/m³] molar densities
        real*8, intent(out) :: D_l(ncomp), D_g(ncomp)   ! Diffusion coefficients
        real*8, intent(out) :: lambda_l, lambda_g       ! Thermal conductivities
        real*8, intent(out) :: mu_l, mu_g               ! viscosities
        real*8, intent(out) :: Dist_coef(ncomp)         ! distribution coefficient (Henry)
        real*8, intent(out) :: sigma                    ! surface tension
        real*8, intent(out) :: H_lv(ncomp)              ! Heat of Evaporation
        real*8, intent(out) :: cp_l, cp_g                ! heat capacity of gas and liquid

        call DB_read(comp_order,condensable,M,xsi,V_mol,psi,&
                            z_ion,cond,mu_g_par,mu_l_par,H_par,&
                            cp_gas_par,T_crit,acentric,dH_excess,&
                            lambda_gas_par,lambda_liq_par,ncomp)
        call liquid_density(ncomp,T_0_l,x_0(:),M(:),rho_l,conc_l)
        call gas_density(ncomp,T_0_g,p_0_g,y_0,M,rho_g,conc_g)
        call gas_DiffCoef(p_0_g,T_0_g,y_0,M,xsi,z_ion,D_g,ncomp)
        call liquid_ThermConduct(ncomp,x_0,T_0_l,lambda_liq_par,V_mol,lambda_l)
        call gas_ThermConduct(ncomp,T_0_g,y_0,lambda_gas_par,lambda_g)
        call gas_cp(ncomp,T_0_g,y_0,cp_gas_par,cp_id,cp_g)
        call liquid_cp(ncomp,comp_order,x_0,acentric,T_crit,T_0_l,&
                        cp_id,cp_l)
        call liquid_viscosity(T_0_l,mu_l)
        call gas_viscosity(ncomp,T_0_g,y_0,z_ion,mu_g_par,M,mu_g)
        call Distribution_Coef(ncomp,T_0_l,p_0_g,Dist_coef)
        call surface_tension(ncomp,x_0,T_0_l,M,sigma)
        call Heat_Evap(ncomp,condensable,acentric,&
                        T_crit,(T_0_g+T_0_l)/2,dH_excess,H_lv)
        call liquid_DiffCoef(T_0_l,x_0,M,psi,V_mol,z_ion,cond&
                            ,mu_l_par,D_l,ncomp)


    end subroutine physical_properties








    subroutine liquid_density(ncomp,T,x,M,rho_liq,conc_l)


        !################################################################################################!
!
!    SUBROUTINE "liquid density_DB"
!
!  Purpose
!    It is used to calculate the density of the liquid mixture.
!
!  Use
!
!    - The equation is based on experimental data [Gme53]; valid from 10° to 30°C
!
!    - The liquid density is calculated from the mass fraction of sodium hydroxide and the temperature
!
!    - It is assumed that sodium hydroxide is dissociated totally. So instead of taking concentration of
!
!      sodium hydroxide you can also take the concentration of hydroxide ions. But you can only do it in
!
!      that way if you set the start concentration of sodium hydroxide to zero (the physical properties are
!
!      calculated with the start concentrations).
!
!    - hydroxide-ion must be component 4
!
!
!  Record of revisions
!    Date        Programmer            Description of change
!    ----        ----------            ---------------------
!    2008/08/01    Anna Janzen          (Original code)
!   2016/31/05  A. Salten           code transferred
!
!################################################################################################!

        implicit none

    integer, intent(in) :: ncomp
    real*8, intent(in):: T                          ! [K] temperature
    real*8, intent(in) :: x(:)                   ! [-] mole fraction
    real*8, intent(in) :: M(:)                   ! [kg/kmol] molecular weights
    real*8 :: w_NaOH                     ! [-] mass fraction of sodium hydroxide
    real*8 :: M_NaOH                     ! [kg/kmol] molecular weight of sodium hydroxide
    integer :: i                         ! [-] control variable
    real*8, intent(out) :: rho_liq                    ! [kg/m³] liquid density
    real*8, intent(out) :: conc_l           ! [kmol/m³] molar density

    M_NaOH=3.9997D1
    w_NaOH=0.D0

    if(ncomp>=4) then
        do i=1,ncomp
            w_NaOH=w_NaOH+x(i)*M(i)
        end do

        w_NaOH=x(4)*M(4)/w_NaOH
    end if

    rho_liq=1072.2D0*w_NaOH-0.41D0*T+1122.D0

    conc_l=rho_liq/dot_product(x,M)

    end subroutine liquid_density

    subroutine gas_density(ncomp,t,p,y,M,rho_gas,conc_g)

        !###################################################
        !
        !   calculates the gas density by use of ideal gas law
        !
        !
        !####################################################
        implicit none

        integer, intent(in) :: ncomp
        real*8, intent(in) :: T,p               !Temperature [K], pressure [Pa], mean molecular weight [kg/kmol]
        real*8 :: M_m
        real*8, intent(in) :: y(:), M(:)    !mole fraction of gaseous components [-], molecular weight [kg/kmol]
        real*8 :: R=8.314472D0          !Universal gas constant [J/(mol*K)]
        real*8, intent(out) :: rho_gas               !gas density [kg/m³]
        real*8, intent(out) :: conc_g           ! [kmol/m³] molar density
        integer :: i

        M_m=0
        do i=1,ncomp
            M_m=M_m+y(i)*M(i)
        end do

        conc_g=p/(R*T)*1D-3!rho_gas/dot_product(y,M)
        rho_gas=conc_g*M_m!1.D-3*p*M_m/(R*T)

    end subroutine gas_density

    !################################################################################################!
!
!    SUBROUTINE "surface tension_DB"
!
!  Purpose
!    It is used to calculate the surface tension
!
!  Use
!
!    - The equation is valid from 20° to 170°C and to a sodium hydroxide mass fraction of 0.7 [Fel69]
!
!    - It is assumed that sodium hydroxide is dissociated totally. So instead of taking concentration of
!
!      sodium hydroxide you can also take the concentration of hydroxide ions. But you can only do it in
!
!      that way if you set the start concentration of sodium hydroxide to zero (the physical properties are
!
!      calculated with the start concentrations).
!
!    - hydroxide-ion must be component 4
!
!
!  Record of revisions
!    Date       Programmer             Description of change
!    ----       ----------              ---------------------
!    2007/06/11  Manuel Braß            (Original code)
!    2008/08/01  Anna Janzen            M(ncomp) is input variable
!   2016/31/05  A. Salten             code transferred/ minor improvements
!
!################################################################################################!


    subroutine surface_tension(ncomp,x,T,M,sigma)


    implicit none


!### input variables ########################################################################

    integer, intent(in) :: ncomp
    real*8, intent(in) :: x(:)                  ![-]mole fraction of component i
    real*8, intent(in) :: T                         ![K] temperature
    real*8, intent(in) :: M(:)                  !array of the molecular weights of the components

!    ### internal variables #####################################################################

    real*8 :: w_NaOH                    ! [-] mass fraction of sulfur dioxide
    real*8 :: M_NaOH                     ! [kg/kmol] molecular weight of sodium hydroxide
    integer :: i,k                      ! control variables
    real*8 :: teta                      ! [°C] temperature
    real*8 :: a(5,5)                    ! coefficient matrix
    real*8 :: sigma_1(5)            ! parameter for calculation of sigma

!    ### output variables #######################################################################

    real*8, intent(out) :: sigma                     ! [N/m] surface tension

!    ### definition of the coefficients, constants and parameters

    M_NaOH=39.997D0

    a(1,1)=+75.4534D0
    a(1,2)=-0.135601D0
    a(1,3)=-0.363901D-03
    a(1,4)=+0.572853D-06
    a(1,5)=-0.385295D-09
    a(2,1)=-1.34722D0
    a(2,2)=-0.757382D0
    a(2,3)=+0.0105068D0
    a(2,4)=-0.402082D-04
    a(2,5)=+0.535398D-07
    a(3,1)=+440.876D0
    a(3,2)=+7.92225D0
    a(3,3)=-0.0937907D0
    a(3,4)=+0.357642D-03
    a(3,5)=-0.474984D-06
    a(4,1)=-522.231D0
    a(4,2)=-23.042D0
    a(4,3)=+0.253931D0
    a(4,4)=-0.964081D-03
    a(4,5)=+0.127614D-05
    a(5,1)=+11.6158D0
    a(5,2)=+19.5414D0
    a(5,3)=-0.205044D0
    a(5,4)=+0.775193D-03
    a(5,5)=-0.102274D-05

!    ### calculations ##########################################################################

    if (ncomp>=4) then

    teta=T-273.15D0

    w_NaOH=0.D0
    do i=1,ncomp
        w_NaOH=w_NaOH+x(i)*M(i)
    end do

    w_NaOH=(x(4)*M_NaOH)/w_NaOH

    sigma_1=0.D0
    sigma=0.D0
    do i=1,5
        do k=1,5
            sigma_1(i)=sigma_1(i)+a(i,k)*teta**(k-1.D0)
        end do
        sigma=sigma+sigma_1(i)*w_NaOH**(i-1.D0)
    end do

    sigma=sigma*1.D-03

    else

    sigma=73.8D-03

    end if

    end subroutine surface_tension

    subroutine gas_viscosity(ncomp,T,y,z,eta_gas_par,M,eta_gas)
!################################################################################################!
!
!    SUBROUTINE "gas viscosity_DB"
!
!  Purpose
!    It is used to calculate the viscosity of the gas mixture.
!
!  Use
!
!    - The viscosity of the gas mixture is calculated from Wilke's method (1950) [Pol01]
!
!    - The viscosity of the pure components are calculated from paramters and equation of [VDI]
!
!    - An approximation from Herning and Zipperer (1936) is used to calculate phi [Pol01]
!
!  Record of revisions
!    Date       Programmer             Description of change
!    ----       ----------              ---------------------
!    2007/06/12  Manuel Braß            (Original code)
!    2008/08/01  Anna Janzen            transformed into general type
!   2016/02/06  A. Salten               code transferred
!
!################################################################################################!
    implicit none



    integer, intent(in) :: ncomp
    real*8, intent(in) :: T                                ! [K] temperature
    real*8, intent(in) :: y(:)                         ! [-] mole fraction of component i,j
    integer, intent(in) :: z(:)                        ! [-] charge of ion i
    real*8, intent(in) :: eta_gas_par(:,:)             ! Parameters for viscosity calculation of pure components in the gas-phase
    real*8, intent(in) :: M(:)                         ! array which contains the molecular weights of the components


    real*8 :: phi(ncomp,ncomp)                 ! [-] square root of the quotient of the molecular weights of component j and i
    real*8 :: eta_pure(ncomp)                  ! [Pa*s] viscosity of the pure components
    real*8 :: A(ncomp)                         ! parameter for calculating gas viscosity
    integer :: i,j                             ! [-] control variables


    real*8, intent(out) :: eta_gas                          ! [Pa*s] viscosity of the gas mixture



    forall (i=1:ncomp,j=1:ncomp)
        phi(i,j)=(M(j)/M(i))**0.5D0             ! for calculation of phi
    endforall

    A=0.D0                                       ! calculation of the sum in equation (70)
    do i=1,ncomp                                !             "
        do j=1,ncomp                               !             "
            if (z(i)==0.and.z(j)==0) then
                A(i) = A(i)+y(j)*phi(i,j)                 !             "
            endif
        end do                                     !             "
    end do                                      !

    eta_pure=0                              !calculation of the dynamic viscosity of the pure non ionic components
        do i=1,ncomp
            eta_pure(i)=eta_gas_par(i,1)+eta_gas_par(i,2)*T+&
                        eta_gas_par(i,3)*T**2+eta_gas_par(i,4)*T**3+&
                        eta_gas_par(i,5)*T**4
        end do

        eta_gas=0                           ! [Pa*s] calculation of the viscosity of the gas mixture
    do i=1,ncomp
        if (z(i)==0) then
            eta_gas=eta_gas+(y(i)*eta_pure(i))/A(i)
        endif
!        eta_gas=eta_gas*1D-10
    end do


    end subroutine gas_viscosity

    subroutine liquid_viscosity(T,mu_l)
!################################################################################################!
!
!    SUBROUTINE "liquid viscosity_DB"
!
!  Purpose
!    It is used to calculate the viscosity of the liquid mixture.
!   From VDI Wärmeatlas Parametres for Water
!
!
!
!  Record of revisions
!    Date       Programmer             Description of change
!    ----       ----------              ---------------------
!    2007/06/12  Manuel Braß            (Original code)
!    2008/08/01  Anna Janzen            deleted do-loop in the calculations
!   2016/02/06  A. Salten              replaced code
!################################################################################################!

    implicit none
    real*8, intent(in) :: T         ![K] Temperature
    real*8, intent(out) :: mu_l   ![Pa*s] liquid viscosity


    mu_l=exp(-22.968D0+3275.89D0/T+1.7637D-02*T+6.93D-07*T**2-1.2933D-08*T**3)

    end subroutine liquid_viscosity





    !################################################################################################!
!
!    SUBROUTINE "phase equilibria_DB"
!
!  Purpose
!    It is used to calculate the phase equilibria.
!
!  Use
!
!    - A Correlation is used to calculate the Henry coefficient for water [Bre04] valid from 0° to 150°C
!
!    - A Correlation is used to calculate the Henry coefficient for sulfur dioxide [Edw78] valid from 0° to 100°C
!
!     - The phase equilibria are calculated from Henry's law [Bre04]
!
!  Record of revisions
!    Date       Programmer             Description of change
!    ----       ----------              ---------------------
!    2007/06/13  Manuel Braß            (Original code)
!
!################################################################################################!

    subroutine Distribution_Coef(ncomp,T,p,K_F)

    implicit none

!    ### input variables ########################################################################

    integer, intent(in):: NCOMP
    real*8, intent(in) :: T                        ![K] temperature
    real*8, intent(in) :: p                     ![Pa] total pressure
    integer :: i                    ! [-] control variable
    real*8 :: gamma(ncomp)         ! [-] activity coefficient of component i
    real*8 :: M_w                   ! [kg/mol] molecular weight of water


!     ### output variables #######################################################################

    real*8, intent(out) :: K_F(ncomp)           ! [-] equilibrium constant at the liquid-gas interface


!     ### calculations ###########################################################################
    gamma=1
    M_w=0.018                         ![kg/mol] Molmasse Wasser

    !calculation of the Henry coefficients

    K_F=0

    if (ncomp>=4) then
     do i=4,ncomp
      K_F(i)=1.D-005
     end do
    endif


    !calculation of the phase-equelibria

    K_F(1)=EXP(66.41843682D0-7362.7D0/T-9.D0*DLOG(T)+0.006952D0*T)/M_w& !611.21D0/p*EXP((17.5132D0-T/234.5D0)*((T-273.15D0)/(T-16.01D0)))!
            *101325.D0/(8.3145D0*T*55.5D3)

    K_F(2)=(3.3D-4)**(-1D0)*DEXP(24D2*(1/298.15D0-1/T))/(8.3145D0*T)!EXP(170.7126D0-8477.711D0/T-21.95743D0*DLOG(T)&
          !+0.005780748D0*T)

    K_F(3)=(6.4D-6)**(-1D0)*DEXP(16D2*(1/298.15D0-1/T))/(8.3145D0*T) !1/EXP(-181.587D0+8632.13D0/T+24.7981D0*DLOG(T))*101325.D0


    !do i=1,ncomp

        !K_F(i)=1.01325D05*gamma(i)*H(i)/(p*x(1)*M_w)
        !K_F(i)=gamma(i)*K_F(i)/p

    !end do

    end subroutine Distribution_Coef



    subroutine liquid_DiffCoef(T,x,M,psi,V_mol,z,lamda&
                            ,eta_liq_par,D_l,NCOMP)

    implicit none

!     ### input variables ########################################################################

    integer :: NCOMP
    real*8 :: T                        ! [K] temperature
    real*8 :: x(:)                 ! [-] mole fraction
    real*8 :: M(:)                 ! [kg/kmol] molecular weight of component j
    real*8 :: psi(:)               ! [-] association factor of solvent j
    real*8 :: V_mol(:)             ! [cm3/s] molar volume of component i at its normal boiling temperature
    integer :: z(:)                 ! [-] charge of ion i
    real*8 :: lamda(:)             ! [cm2/(ohm*mol)] equivalent conductance of ion i
    real*8 :: eta_liq_par(:,:)       ! Parameters for viscosity calculation of pure components in the liquid-phase

!    ### internal variables #####################################################################


    real*8 :: eta_pure_liq (ncomp)     ! [Pa*s] liquid viscosity of component j
    real*8 :: D_bin_inf(ncomp,ncomp)   ! [m2/s] binary diffusion coefficient of non ionic components at infinite dilution
    real*8 :: D_bin_conc(ncomp,ncomp)  ! [m2/s] binary diffusion coefficient at concentrated dilution
    real*8 :: D_ion_conc(ncomp)  ! [m2/s] diffusion coefficient of an ion in an electrolyte mixture [Asp00]
    real*8 :: sum_x                       ! Parameter for calclation of D_ion_conc
    integer :: i,j                      ! [-] control variables
    real*8, parameter :: R= 8.314472D0                       ! [J/(mol*K)] universal gas constant
    real*8, parameter :: F=96500.D0                        ! [C/mol] Faraday number
    real*8 :: D_eff_liq_sum(ncomp)       ! [m2/s] parameter for calculation of D_l
      !cmya
    real*8 :: w_NaOH                     ! [-] mass fraction of sodium hydroxide
    real*8 :: M_NaOH                     ! [kg/kmol] molecular weight of sodium hydroxide
      !cmya
!    ### output variables #######################################################################

    real*8, intent(out) :: D_l(ncomp)         ! [m2/s] effective diffusion coefficient of component i

!    ### calculations ###########################################################################

    !dynamic viscosity of the pure non ionic components:

    eta_pure_liq=0.D0
        do i=1,ncomp
        eta_pure_liq(i)=exp(eta_liq_par(i,1)+&
                    eta_liq_par(i,2)/T+eta_liq_par(i,3)*T+eta_liq_par(i,4)*T**2&
                    +eta_liq_par(i,5)*T**3)
        end do

    !diffusion coefficients of the ionic components

    sum_x=0.D0
        do i=1,ncomp
            if(z(i)==0) then
            sum_x=sum_x+x(i)
            end if
        end do

    D_ion_conc=0.D0
        do i=1,ncomp
            if(z(i)/=0) then
     D_ion_conc(i)=lamda(i)*(R*T/(dabs(DBLE(z(i)))*F**2))*sum_x*1.D-04
            end if
        end do

    !Diffusion Coefficients of non ionic components at infinite dilution:

    D_bin_inf=0
    do i=1,ncomp
     do j=1,ncomp
      if (i.ne.j) then
        if(z(i).eq.0 .and. z(j).eq.0) then
       D_bin_inf(i,j)=7.4D-15*dsqrt(psi(j)*M(j))*T/(eta_pure_liq(j)*&
                         (V_mol(i)**0.6))
        end if
      end if
     end do
    end do

    !binary diffusion coefficients at concentrated dilution:

    D_bin_conc=0.D0
    do i=1,ncomp
     do j=1,ncomp
      if (i.ne.j) then
        if(z(i)==0.and.z(j)==0) then
        D_bin_conc(i,j)=(D_bin_inf(i,j)**((1+x(j)-x(i))/2.D0))*&
                           (D_bin_inf(j,i)**((1+x(i)-x(j))/2.D0))
        else if(z(i)==0.and.z(j)/=0) then
        D_bin_conc(i,j)=D_ion_conc(j)
        else if(z(i)/=0.and.z(j)==0) then
        D_bin_conc(i,j)=D_ion_conc(i)
        else
        D_bin_conc(i,j)=0.5D0*(D_ion_conc(i)+D_ion_conc(j))
        endif
      end if
     end do
    end do

    !effective Diffusion coefficients for all components in the mixture:

    D_eff_liq_sum=0.D0
    do i=1,ncomp
     do j=1,ncomp
      if (i.ne.j) then
       D_eff_liq_sum(i)=D_eff_liq_sum(i)+(x(j)/D_bin_conc(i,j))
      end if
     end do
    end do

    do i=1,ncomp
     D_l(i)=(1-x(i))/D_eff_liq_sum(i)
    end do
!mya liquid_diffusion_coefficient_DB

    M_NaOH=39.997D0

    w_NaOH=0.D0

    if (ncomp>=4) then

     do i=1,ncomp
      w_NaOH=w_NaOH+x(i)*M(i)
     end do

      w_NaOH=(x(4)*M_NaOH)/w_NaOH

    endif

    end subroutine liquid_DiffCoef




    subroutine gas_DiffCoef(p,T,y,M,xsi,z,D_g,ncomp)
    implicit none

!    ### input variables ########################################################################

    integer :: NCOMP
    real*8 :: T,p                      ! [K,Pa] temperature, pressure
    real*8 :: M(NCOMP)                 ! [kg/kmol] molecular weight of pure component
    real*8 :: xsi(NCOMP)               ! [-] Fuller et al. volume parameter
    integer :: z(NCOMP)                   ! [-] charge of ion i
    real*8 :: y(NCOMP)                   ! [-] mole fraction of the component

!    ### internal variables #####################################################################

    real*8 :: M_bin(NCOMP,NCOMP)         ! [kg/kmol] molecular weight of binary mixture
    real*8 :: D_bin_gas(NCOMP,NCOMP)     ! [m2/s] binary diffusion coefficient
    integer :: i,j                       ! [-] control variables

! ### output variables #######################################################################

 real*8 :: D_g(NCOMP)         ! [m2/s] effective diffusion coefficient

!    ### calculations ###########################################################################


    do i=1,ncomp
     do j=1,ncomp
      if (i.ne.j) then
      M_bin(i,j)=2.D0*((1.D0/M(i))+(1.D0/M(j)))**(-1)
          end if
     end do
    end do

    D_g=0.D0
    D_bin_gas=0.D0
    do i=1,ncomp
     do j=1,ncomp
      if (i.ne.j) then
      if (z(i)==0.and.z(j)==0) then
       D_bin_gas(i,j) = (10.D0*0.00143D0*(T**1.75D0))/&![m2/s] calculation of the binary diffusion coefficients
                            (p*(M_bin(i,j)**0.5D0)*((xsi(i)**(1.D0/3.D0))+&
                            (xsi(j)**(1.D0/3.D0)))**2)
         D_g(i)=D_g(i)+(y(j)/D_bin_gas(i,j))                  ![m2/s] calculation of the sum in equation of effective diffusion coefficients
      end if
      end if
     end do
    end do

    do i=1,ncomp                                                        !calculation of the effective diffusion coefficients [m²/s]
        if(z(i)==0) then
        D_g(i)=(1.D0-y(i))/D_g(i)
        else
        D_g(i)=1.D-20
        endif
    end do

    end subroutine gas_DiffCoef


    subroutine liquid_ThermConduct(ncomp,x,temp,lambda_liq_par,V_mol,lambda_liq_mix)

    implicit none

!    ### input variables ######################################################

    integer :: ncomp,i,j
    real*8  :: lambda_liq_par(:,:)    ![J/(s m K)] parameter for calculating the heat conductivity of pure liquids
    real*8  :: temp             ! [K] temperature
    real*8  :: x(:)                ! [mol/mol] liquid phase 1 molar fraction
    real*8  :: V_mol(:)            ! [cm3/mol] molar volume of component i at its normal boiling temperature
    real*8    :: lambda_bin(ncomp,ncomp)    ! [J/(s m K)] binary conductivity
    real*8    :: phi(ncomp)                    ! [-] superficial volume fraction

!   ### internal variables ###################################################

    real*8  :: temp_max                ! [K]
    real*8  :: temp_min                ! [K]
    real*8  :: lambda_liq(ncomp)         ! [J/(s m K)] conductivity of pure liq

!    ### output variables #####################################################

    real*8  :: lambda_liq_mix                ! [J/(s m K)] conductivity of liquid mixture

!    ##########################################################################
!    ### calculations #########################################################

    temp_max = 623.d0

    temp_min = 64.d0


    if ((temp.le.temp_min).or.(temp.ge.temp_max)) then ! validity check


     write (*,*) "An error occured in subroutine COND_LIQ. The&
     &              temperature is not valid"

     lambda_liq = 0.d0

    else

     lambda_liq = lambda_liq_par(:,1)+temp*lambda_liq_par(:,2)+temp**2*lambda_liq_par(:,3)

    endif

! eq. 10-12.18
    do i=1,ncomp
     do j=1,ncomp
      lambda_bin(i,j)=2.d0/(1.d0/lambda_liq(i)+1.d0/lambda_liq(j))
     enddo
    enddo

!    # eq. 10-12.19 with "molar volume = molar mass/density"

      phi = x*V_mol/dot_product(x,V_mol)

!    # eq. 10-12.17
    lambda_liq_mix = 0.d0
    do i=1,ncomp
     do j=1,ncomp
      lambda_liq_mix=lambda_liq_mix+phi(i)*phi(j)*lambda_bin(i,j)
     enddo
    enddo


    end subroutine liquid_ThermConduct



!#############################################################################!
!
!    SUBROUTINE COND_GAS, calculates the heat conductivity of pure gases
!
! Reference
!    [1]    table 10-3
!
!#############################################################################!

    subroutine gas_ThermConduct(ncomp,temp,y,lambda_gas_par,lambda_g)
    implicit none

!    ### input variables ######################################################

    integer :: ncomp                ! total number of components
    real*8  :: temp                    ! [K] temperature
    real*8 :: y(:)
    real*8  :: lambda_gas_par(:,:)    ! parameter for calculating the heat conductivity of pure gases

!   ### internal variables ###################################################

    real*8  :: temp_max                ! [K]
    real*8  :: temp_min                ! [K] sollte die Siedetemperatur der komponente sein
    real*8  :: lambda_gas(ncomp)         ! [J/(s m K)] conductivity of pure components

!    ### output variables #####################################################

    real*8 :: lambda_g              ! Conductivity of gas mixture

!    ##########################################################################
!    ### calculations #########################################################

    temp_max = 1070.d0

    temp_min = 115.d0

    if ((temp.le.temp_min).or.(temp.ge.temp_max)) then ! validity check

        write (*,*) "An error occured in subroutine COND_GAS. The temperature is not valid"

        lambda_gas = 0.d0

    else

        lambda_gas=lambda_gas_par(:,1)+temp*lambda_gas_par(:,2)+temp**2*&
                    lambda_gas_par(:,3)+temp**3*lambda_gas_par(:,4)

    endif

    lambda_g=dot_product(lambda_gas,y)


    endsubroutine gas_ThermConduct



!#############################################################################!
!
!    SUBROUTINE CP_ID_GAS, calculates heat capacities of pure ideal gases
!
! Reference
!    [1]    p. 191/201, 657, 668
!
!#############################################################################!

    subroutine gas_cp(ncomp,temp,y,cp_gas_par,cp_id,cp_g)
    implicit none

!    ### input variables ######################################################
    integer, intent(in) :: ncomp
    real*8  :: temp             ! [K] temperature
    real*8 :: y(:)
    real*8 :: cp_gas_par(:,:)        ! parameter for calculating the heat capacity of pure gases

!    ### output variables #####################################################

    real*8 :: cp_id(ncomp)         ! [J/(mol K)] heat capacity of the ideal gascomponent
    real*8, intent(out) :: cp_g     ! [J/(kmol K)]
!    ##########################################################################
!    ### calculations #########################################################


     cp_id(:) = cp_gas_par(:,1)+temp*cp_gas_par(:,2)+temp**2*&
              cp_gas_par(:,3)+temp**3*cp_gas_par(:,4)


                cp_g= 1E3*dot_product(cp_id,y)


    end subroutine gas_cp


!#############################################################################!
!
!    SUBROUTINE CP_LIQUID, calculates heat capacities of pure liquids
!
! Require subroutines
!    cp_id_gas
!
! Reference
!    [1]    p. 140, pp. 658, [2]
!
!#############################################################################!

    subroutine liquid_cp(ncomp,comp_order,x,acentric,temp_crit,temp,&
                        cp_id_g,cp_l)
    implicit none

!    ### input variables ######################################################

    integer :: i            ! component index of pure component
    integer :: ncomp
    character*14 :: comp_order(:)
    real*8    :: acentric(:)! [-] Pitzer acentric factor listed in [1]
    real*8, intent(in) :: x(:)
    real*8    :: temp_crit(:)    ! [K] critical temperature listed in [1]
    real*8  :: temp             ! [K] temperature


!   ### internal variables ###################################################

    real*8  :: r_id_gas = 8.314D0 !  [J/(mol K)] gas constant
    real*8  :: temp_redu            ! [-] reduced temperature = temp/temp_crit
    real*8    :: nist_a                ! [J/(mol K)]  nist coefficient A from [2]
    real*8    :: nist_b                ! [J/(mol K2)] nist coefficient B from [2]
    real*8    :: nist_c                ! [J/(mol K3)] nist coefficient C from [2]
    real*8    :: nist_d                ! [J/(mol K4)] nist coefficient D from [2]
    real*8    :: nist_e                ! [J K/mol]    nist coefficient E from [2]
    real*8  :: temp_max                ! [K] validity range of function [1]
    real*8  :: temp_min     = 0.d0        ! [K]
    real*8  :: temp2_min = 298.0d0     ! [K] validity range of function [2]
    real*8  :: temp2_max = 500.0d0  ! [K]
    real*8  :: temp2_extr = 20.0d0  ! [K] range of extrapolation
    real*8  :: cp_liq(ncomp)         ! [J/(mol K)] heat capacity of the ideal liquidcomponent

!    ### output variables #####################################################

    real*8 :: cp_l                      ! [J/(kmol K)] liquid heat capacity
    real*8 :: cp_id_g(ncomp)            ! [J/(mol K)] heat capacity of the ideal gas component

!    ##########################################################################
!    ### initialisation #######################################################


    nist_a = -203.6060d0
    nist_b = 1523.290d0
    nist_c = -3196.413d0
    nist_d = 2474.455d0
    nist_e = 3.855326



!    ##########################################################################
!    ### calculations #########################################################



    DO i=1,ncomp

    temp_max = temp_crit(i)

    temp_redu = temp/temp_crit(i)


    if (comp_order(i).eq."H2O") then            ! using [2] for water

       if((temp.le.temp2_min).or.(temp.ge.temp2_max)) then ! validity check

         if ((temp.le.(temp2_min-temp2_extr)).or.(temp.ge.(temp2_max+&
                 temp2_extr))) then


         write (*,*) "An error occured in subroutine CP_LIQUID. The&
                        & temperature is not valid"

         cp_liq = 0.d0

         RETURN

         else ! validity check II

         write (*,*) "WARNING from subroutine CP_LIQUID: temperature&
                        & is slightly above or below validity limit. "

         !NO RETURN

         endif

       endif

       cp_liq(i) = nist_a+nist_b*temp/1000.d0+nist_c*(temp/1000.d0)*&
                (temp/1000.d0)+nist_d*(temp/1000.d0)*(temp/1000.d0)*&
                (temp/1000.d0)+nist_e/((temp/1000.d0)*(temp/1000.d0))

    else                        ! using [1] for other components

       if((temp.le.temp_min)) then ! validity check


        write (*,*) "An error occured in subroutine CP_LIQUID.&
        & The temperature is not valid"

        cp_liq = 0.d0

        RETURN

        elseif    (temp.ge.temp_max) then ! validity check II

      write (*,*) "WARNING from subroutine CP_LIQUID: for component",comp_order(i),&
                    & "temperature is above upper validity limit. Reduced Temperature is set to 0.9."

        temp_redu = 0.9

        !NO RETURN

       endif


     cp_liq(i) = cp_id_g(i)+r_id_gas*(1.45d0+0.45d0/(1.0d0-temp_redu)+&
      0.25D0*acentric(i)*(17.11d0+25.2d0*((1.0d0-temp_redu)**(1./3.))*&
       1.d0/temp_redu+1.742d0/(1.0d0-temp_redu)))

    endif

    END DO

    cp_l= 1E3*dot_product(cp_liq,x)


    end subroutine liquid_cp


  !#############################################################################!
!
!    SUBROUTINE HEAT_EVAP_SOLU, calculates the heat of evaporation and solution
!   in water of pure components
!
! Reference
!    [1]    eq. 7-9.4
!
! Verbesserung und Erweiterung auf polare Komponenten
! kann anhand folgender Quellen durchgeführt werden
! - Halm, R.L. and L.I. Stiel: Aiche Journal, Vol 13, p. 351-355, 1967
! - [1] S. 222 ff.
!
! Zur Berechnung von Lösungsenthalpien wurde nicht viel gefunden. Daher wird
! die Lösungsenthalpie für CO2 nach
! - Axel Schönbucher: Thermische Verfahrenstechnik, Springer, 2002. S. 537 ff.
! abgeschätzt und für N2 =0 gesetzt
!
!#############################################################################!

    subroutine Heat_Evap(ncomp,condensable,acentric,&
                        temp_crit,temp,h_solu,dH)
    implicit none

!    ### input variables ######################################################

    integer :: ncomp
    logical :: condensable(ncomp)        ! if true, gas phase component (=vapour) is
                                ! condensable at ambient conditions
    real*8    :: acentric(ncomp)            ! [-] Pitzer acentric factor listed in [1]
    real*8    :: temp_crit(ncomp)        ! [K] critical temperature listed in [1]
    real*8  :: temp             ! [K] temperature
    real*8  :: h_solu(ncomp)            ! [J/mol] heat of solution in water at ambient conditions

!   ### internal variables ###################################################

    real*8, parameter  :: r_id_gas = 8.314D0 !  [J/(mol K)] gas constant
    real*8  :: temp_redu(ncomp)            ! [-] reduced temperature = temp/temp_crit
    real*8  :: temp_min(ncomp)                ! [K]
    real*8, parameter  :: temp_extr = 250.d0    ! [K] range of extrapolation
    real*8  :: dHLV(ncomp)                    ! [J/mol] heat of evaporation
    integer :: i



!    ### output variables #####################################################

    real*8  :: dH(ncomp)        ! [J/kmol] Enthalpy of the component i, sum of heat of solution and evaporation

!    ##########################################################################
!    ### calculations #########################################################

    temp_min = 0.6d0 * temp_crit

    temp_redu = temp/temp_crit
do i=1,ncomp

    dHLV(i) = 0.d0
    dH(i) = 0.D0


       if (acentric(i)*temp_crit(i).le.(1.0d-10)) then ! error handler

        write (*,*) "An error occured in subroutine Heat_Evap.&
                         &One or more paramters equal 0"

        RETURN

       elseif((temp.le.temp_min(i)).or.(temp.ge.temp_crit(i))) then ! validity check

      if ((temp.le.(temp_min(i)-temp_extr)).or.(temp.ge.(temp_crit(i)))) then

         write (*,*) "WARNING from subroutine Heat_Evap. Heat of evaporation was set to 0"

         RETURN

        else ! validity check II

         write (*,*) "WARNING from subroutine Heat_Evap: temperature is slightly above or below validity limit."

         !NO RETURN2

        endif

       endif

    if (condensable(i)) then

    dHLV(i)=r_id_gas*temp_crit(i)*(7.08d0*((1.d0-temp_redu(i))**0.354)&
             +10.95d0*acentric(i)*((1.d0-temp_redu(i))**0.456))

    endif

    dH(i)=(dHLV(i)+H_SOLU(i))*1.D03
enddo


    end subroutine Heat_Evap






    !################################################################################################!
!
!    SUBROUTINE "DB_read"
!
!  Purpose
!    It is used to read the properties of the components needed in the main subroutine from a csv-file.
!    The properties are first read as a big Matrix, then split into Vectors and Matritzen.
!
!
! List of references
!   [1]    Reid, R.C., R. Prober and B.E. Poling: The Properties of Gases and
!        Liquids. 4 ed.; Mc-Graw-Hill: New York, 1987.
!   [2] Chase, M.W., Jr., NIST-JANAF Themochemical Tables, Fourth Edition,
!       J. Phys. Chem. Ref. Data, Monograph 9, 1998, 1-1951.
!        http://webbook.nist.gov/cgi/inchi/InChI%3D1S/H2O/h1H2
!
!
!
!  Record of revisions
!    Date        Programmer             Description of change
!    ---------   ----------              ---------------------
!    2008/08/18   Anna Janzen            (Original code)
!
!################################################################################################!


    subroutine DB_read(comp_order,condensable,M_weight,xsi,V_mol,psi,&
                            z_ion,lamda,eta_gas_par,eta_liq_par,H_par,&
                            cp_gas_par,T_crit,acentric,dH_excess,&
                            lambda_gas_par,lambda_liq_par,ncomp_ncc)

    implicit none

!    ### input variables ########################################################################

    integer :: ncomp_ncc
    character*14 :: comp_order(NCOMP_NCC)  ! Vector with names of the the components in the order like they are used in the main subroutine/Aspen Plus

!     ### internal variables #####################################################################

    logical:: fileexists                        !variable needed to see if a Database-csv-file is in the folder or not
    integer:: i,j                                ! count variables
    integer:: total_set(2)                        !total number of rows and columns in the DB-file
    integer:: n,m                                !n-> number of rows, m->number of columns with data in the Database-csv-file
    real*8, dimension(:,:),allocatable :: Mat_all, Mat_need    !matrices with all columns of DB
    character*14,dimension(:), allocatable:: comp_name            !vector with names of the components in the DB

!    ### output variables #######################################################################

    logical :: condensable(ncomp_ncc)        ! if true, gas phase component (=vapour) is
                                            ! condensable at ambient conditions
    real*8 :: M_weight(NCOMP_NCC)           ! [kg/kmol] molecular weight of component i
    real*8 :: xsi(NCOMP_NCC)                ! [-] Fuller et al. volume parameter
    real*8 :: V_mol(NCOMP_NCC)                ! [cm3/mol] molar volume of component i at its normal boiling temperature
    real*8 :: psi(NCOMP_NCC)                ! [-] association factor of solvent j
    integer :: z_ion(NCOMP_NCC)                ! [-] charge of ion i
    real*8 :: lamda(NCOMP_NCC)                ! [cm2/(ohm*mol)] equivalent conductance of ion i
    real*8 :: eta_gas_par(NCOMP_NCC,5)        ! Parameters for viscosity calculation of pure components in the gas-phase
    real*8 :: eta_liq_par(NCOMP_NCC,5)        ! Parameters for viscosity calculation of pure components in the liquid-phase
    real*8 :: H_par(NCOMP_NCC,4)            ! [-] array with parameters for calculation of the Henry-coefficient
    real*8 :: cp_gas_par(NCOMP_NCC,4)        ! parameter for calculating the heat capacity of pure gases
    real*8 :: T_crit(NCOMP_NCC)                ! [K] critical temperature
    real*8 :: acentric(NCOMP_NCC)            ! [-] Pitzer acentric factor listed in [1]
    real*8 :: dH_excess(NCOMP_NCC)            ! excess enthalpiy/ absorption heat
    real*8 :: lambda_gas_par(NCOMP_NCC,4)    ! parameter for calculating the heat conductivity of pure gases
    real*8 :: lambda_liq_par(NCOMP_NCC,3)    ! parameter for calculating the heat conductivity of pure liquids

!     ### definition of the columns in the Database file: ###############################

    !1: Name of component i
    !2: M_weight
    !3: xsi
    !4: V_mol
    !5: psi
    !6: z_ion
    !7: lamda
    !8-12: eta_gas_par
    !13-17: eta_liq_par
    !18-21: H_par


!     ################### Read counter ##################################################

      fileexists = .false.
    inquire (file='initialisation/DB_HT.csv',exist=fileexists)
      IF (fileexists) THEN
        do i=1,2
            open(13,file='initialisation/DB_HT.csv')
            read(13,*)    total_set(i)
        end do
    else

        write(*,*)  "The physical parameter file: initialisation/DB_HT.csv is missing"
    return

    ENDIF



    n=total_set(1)-2
    m=total_set(2)-1

    allocate(Mat_all(n,m))
    allocate(Mat_need(NCOMP_NCC,m))
    allocate(comp_name(n))


    Mat_all=0.D0
    do i=1,n
        read(13,*) comp_name(i),Mat_all(i,1:m)
    end do

    close(13)

!    ##################extraction of the parameters needed in the main subroutine########

    Mat_need=0.D0
    do i=1,n
        do j=1,NCOMP_NCC
         if (comp_name(i).eq.comp_order(j)) then
         Mat_need(j,:)=Mat_all(i,:)
         endif
        end do
    end do

    do i=1,NCOMP_NCC

        if (Mat_need(i,1)==1.D0) then
            condensable(i)=.true.
        else
            condensable(i)=.false.
        endif

        M_weight(i)=Mat_need(i,2)
        xsi(i)=Mat_need(i,3)
        V_mol(i)=Mat_need(i,4)
        psi(i)=Mat_need(i,5)
        z_ion(i)=int(Mat_need(i,6))
        lamda(i)=Mat_need(i,7)
        eta_gas_par(i,1:5)=Mat_need(i,8:12)
        eta_liq_par(i,1:5)=Mat_need(i,13:17)
        H_par(i,1:4)=Mat_need(i,18:21)
        cp_gas_par(i,1:4)=Mat_need(i,22:25)
        T_crit(i)=Mat_need(i,26)
        acentric(i)=Mat_need(i,27)
        dH_excess(i)=Mat_need(i,28)
        lambda_gas_par(i,1:4)=Mat_need(i,29:32)
        lambda_liq_par(i,1:3)=Mat_need(i,33:35)
    end do

      end subroutine DB_read



end module
