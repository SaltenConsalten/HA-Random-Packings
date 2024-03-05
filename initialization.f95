module initialization
USE global_variables
use omp_lib
    implicit none

    contains

    subroutine flow_parametres(feedl,feedg,mu_l, conc_l, conc_g,sigma, R_c, FVOLL, FVOLG, q_l, q_g, zl_film,zl_drop,&
                                zg,H,l_mix,verh,ae,rho_l)!,rho_g)
        ! receive geometric parametres from file and calculate packing parametres
        implicit none

        real*8, intent(in) :: feedg,feedl       !Feeds of gas and Liquid (reboiler and condenser)
        real*8, intent(in) :: mu_l              !liquid viscosity
        real*8, intent(in) :: conc_l,rho_l             !liquid density
        real*8, intent(in) :: conc_g!,rho_g             !gas density
        real*8, intent(in) :: sigma
        real*8, intent(out) :: zl_film,zl_drop, zg           ! length of undisturbed fluid/gas flow
        real*8, intent(out) :: R_c               !channel radius
        real*8, intent(out) :: FVOLL            !liquid volume flux
        real*8, intent(out) :: FVOLG            !gas volume flux
        real*8, intent(out) :: q_l              !liquid flux per wetted channel
        real*8, intent(out) :: q_g              !gas flux per channel
        real*8, intent(out) ::verh(3)           !volume fraction of film, drop and jet flow


        character(20) :: Pack_type
!        integer :: num       !number of channels
!        integer :: k, i        !control variable
        real*8 :: NC, NC_w  ! number of channels / number of wetted channels
        real*8 :: uls, ugs  ! superficial liquid/gas velocity
        real*8 :: ae_korrel, ae        ! wetted area from korrelation, recalculated because jets and drops have smaller surface than films
        real*8 :: D_Pac     ! packing diametre [m]
        real*8 :: A, winkel         ! packing specific area [m²/m³]
        real*8 :: B_0       ! basis length[m]
        real*8 :: H         ! packing height [m]
        real*8 :: LFK, L_u       ! characteristic length of packing [m], flow length circumference
        real*8 :: Phi       ! ratio of void area to geometrical surface area [m²/m²]
        real*8 :: Ga        ! average position angle [°]
        real*8 :: Eps       ! void fraction of the packing
        real*8 :: Re, Fr, We            !Reynolds, Froude and Weber number
        real*8 ::pi = 4*atan(1D0)
        real*8 :: psi, l_mix
        real*8 :: C_ae, ae_err=1          !liquid flow depending factor for ae transfer


        logical :: Isotherm         ! if true only isothermal calculations will be performed



        NAMELIST /GlobRunPar/ Pack_Type,D_pac,A,B_0,H,LFK,PHI,GA,Eps,Isotherm


        !execution section
        !write(*,*) 'Achtung Überprüfen der Winkelberechnung Unterschiede Dissertation / Implementierung'
        !calculate length of undisturbed liquid flow




open (2,file='initialisation/globalparameters.dat',DELIM='APOSTROPHE')

      read (2,GlobRunPar)

      close(2)

      R_c=2D0*eps/a


if (Pack_type.eq.'Pall Ring 50') then

    zl_Drop=-1D0/log(Phi)*LFK/(6.3D-1)

    Ga=Ga/360*2*pi

    if(B_0/cos(pi/2E0-Ga).gt.dsqrt(B_0**2+LFK**2)) then        !prüfung, ob die planare mischungslänge größer als der halbe umfang des FK ist
        zl_Film=dsqrt(B_0**2+LFK**2)
        B_0=zl_film*cos(pi/2E0-Ga)                              !Basislängenkorrektur
    else
        zl_Film=B_0/cos(pi/2E0-Ga)     !planar mixing length/projection of Packing
    end if

        l_mix = zl_Film*cos(Ga)


    if(zl_film*sin(pi/2E0-Ga).gt.LFK)then      !berechnung der Projektion in umfangsrichtung
        L_u=LFK
    else
        L_u=zl_film*sin(pi/2E0-Ga)
    end if

    psi=pi-2E0*acos(L_u/LFK)               !Überströmwinkel

    zl_film=dsqrt(B_0**2+(psi*LFK/2E0)**2)

    zg=LFK*sin(pi/4D0)
else
    zl_Drop=-1D0/log(Phi)*LFK/(6.95D-1)
    Ga=Ga/360*2*pi
    zl_film=pi/4D0*LFK/3D0
    zg=LFK
    l_mix=LFK/6D0
endif

        FVOLL=FEEDL/conc_l

        FVOLG=FEEDG/conc_g

        ugs=FVOLG/(pi/4D0*D_Pac**2)

        uls=FVOLL/(pi/4D0*D_Pac**2)

        NC=Eps*D_pac**2/(R_c**2*4)
!##########################################################################
 !       R_c=R_c/2D0     !Parametervariation Aufpassen
!##########################################################################

        q_g=FVOLG/NC


        Re=rho_l*uls/(A*mu_l)

        Fr=uls**2D0*A/(9.81D0)
        We=rho_l*uls**2D0/(a*sigma)


        if (Pack_type.eq.'Pall Ring 50') then

            B_=2.212797246D0

            verh(2)=1D-2*(13.4439818263535D0*(uls*1D3)**(1.5D0)*EXP(-1.02494133067777D0*(uls*1D3))+16.1240627089786D0)
!            verh(1)=(1-verh(2))*9D-1
!            verh(3)=(1-verh(2))*1D-1
            verh(1)=(1-verh(2))
            verh(3)=verh(2)*7D-1
            verh(2)=verh(2)*3D-1

        else            !Envipac3

            B_=1.86578208D0

            verh(2)=1D-2*(17.6303082620858D0*(uls*1D3)**4D0*EXP(-2.45486576309387D0*(uls*1D3))+11.3339202952846D0)

            if(q_g.gt.7D-4) then
                 C_ae=2D-1
                verh(2)=((verh(2)-((6.0993D-4*(C_ae+verh(2))-3.4548D-3*verh(2))/(6.0993D-4-3.4548D-3)))/6.0993D-4)*q_g&
                        +((6.0993D-4*(C_ae+verh(2))-3.4548D-3*verh(2))/(6.0993D-4-3.4548D-3))

            endif

!            verh(1)=(1-verh(2))*9D-1
!            verh(3)=(1-verh(2))*1D-1
            verh(1)=(1-verh(2))
            verh(3)=verh(2)*5D-1
            verh(2)=verh(2)*5D-1
        end if

        ae_korrel=1.045D0*Re**4.1D-2*We**1.33D-1*(sigma/29D-3)**(-1.82D-1)    !0.509598D0!3.42D0*(Fr)**(-1D0/6D0)*(We)**(5D-1)

        ae=ae_korrel

        do while (abs(ae_err).gt.1D-6)

            !ae=1D0
            NC_w=ae*NC

            q_l=FVOLL/NC_w


            C_ae=1/(verh(1)+verh(2)*1.1867D-1*q_l**2.338D-1+verh(3)*6.033D3*q_l**1.0079D0)

            ae_err=(ae-C_ae*ae_korrel)/ae

            ae=C_ae*ae_korrel

        end do
        q_l=FVOLL/NC
!        winkel=1.02D2/3.6E2*2*pi
!        ae=0.76*(eps**-6D-1)*((Fr*We)**15D-2)*(138*LFK**1.1D0)*(Re**-2D-1)*(1-0.93*cos(winkel))**-1D0
        !calculation of number of (wetted) channels


    end subroutine flow_parametres




    subroutine input_streams(ncomp,feed_1_m,feed_2_m,T1IN_m,T2IN_m,p1ini_m,p2INi_m,x_0_m,y_0_m,comp_order)

integer :: ncomp                              ! number of components
real*8 :: feed_1_m,feed_2_m                   ! [kmol/s]
real*8, allocatable :: T1IN_m(:), T2IN_m(:)   ! [K]
real*8 :: p1INi_m,p2INi_m  ! [Pa]
real*8, allocatable :: x_0_m(:),y_0_m(:)      ! [kmol/kmol]
character*14, allocatable :: comp_order(:)
logical :: l=.false.

namelist /ManInput/ feed_1_m,feed_2_m,T1IN_m,T2IN_m,p1ini_m,&
                    p2ini_m,x_0_m,y_0_m,comp_order

    inquire (file='initialisation/ManualInputStreams.dat', exist=l)
    open(unit=32,file='initialisation/ManualInputStreams.dat', DELIM='APOSTROPHE', access='sequential')
    read(32,*) ncomp

    allocate(x_0_m(ncomp),y_0_m(ncomp), comp_order(ncomp))
    allocate(T1IN_m(ncomp),T2IN_m(ncomp))

    read (unit=32,NML=ManInput)
    close(unit=32)
 !   write(*,*) 'file closed'
!call sleep(5)
    IF((sum(x_0_m).NE.1).AND.(sum(y_0_m).NE.1)) then

       write(*,*) 'The molar input concentrations do not sum up to one'

    END IF

    end subroutine input_streams



subroutine grid(NNl, NNg, MM, eta_liq, eta_gas)
        ! open gridparameter extrahiere NNg, NNl, a,b, MM

        implicit none

        integer, intent(out) :: nnl     !number of liquid side nodes
        integer, intent(out) :: nng     !number of gas side nodes
        integer, intent(out) :: MM      ! number of nodes per liquid mixing length
        integer :: i                    !control variable
        real*8 :: a,b               !grid parameter
        real*8, allocatable, intent(out) :: eta_liq(:), eta_gas(:) !dimensionless radial coordinates
        namelist /GRPGrid/  NNl,NNg,MM,a,b

        open(32,file='initialisation/gridparameters.dat',DELIM='APOSTROPHE', access='sequential')

            read(32,NML=GRPGrid)

        close(32)

        allocate(eta_gas(NNg),eta_liq(NNl))

        do i=1,NNg
            eta_gas(i) = (((10.D0)**a+b)*(dble(i)/dble(NNg+1))**a)/((10.D0*dble(i)/dble(NNg+1))**a+b)
        end do

        !a=2D0
        !b=1D1
        do i=1,NNl
            eta_liq(nnl+1-i) = 1-(((10.D0)**a+b)*(dble(i-1)/dble(NNl))**a)/((10.D0*dble(i-1)/dble(nnl))**a+b)
        end do

    end subroutine




    subroutine interpol(r11,r12,r21,r22,r31,r32,H,x_int)
implicit none

integer :: zx, zr
real*8,intent(in) :: r11(NNl),r12(NNg),r21(NNl),r22(NNg),r31(NNl),r32(NNg)
real*8, intent(in) :: H
integer :: NNl_old, NNg_old, MM_old, NN_old
integer :: i,j,k,l
real*8 :: dr(3)
real*8, allocatable :: r_old(:,:), x_old(:)
real*8 :: x(MM*nzl), r_new(3,NN)
real*8, allocatable :: c(:,:,:,:)
real*8 :: a,b
logical :: v=.false.
real*8, allocatable, intent(out) :: x_int(:,:,:,:)

allocate (x_int(3,ncomp,NN,MM*nzl))

inquire(file='output/Grid-old.dat', exist=v)
if (v) then
    open(12, file='output/Grid-old.dat')
        read(12,*) NNl_old
        read(12,*) NNg_old
        read(12,*) MM_old
    close(12)
else
    x_int=0
    write(*,*) 'No interpolation is performed. Standard initialisation is used.'
    goto 100
endif

NN=NNl+NNg
NN_old=NNl_old+NNg_old

r_new(1,1:NNl)=r11
r_new(1,NNl+1:NN)=r12
r_new(2,1:NNl)=r21
r_new(2,NNl+1:NN)=r22
r_new(3,1:NNl)=r31
r_new(3,NNl+1:NN)=r32
do i=1,MM*nzl
    x(i)=H-H/dble(MM*nzl)*i
end do


dr(1)=dble(NNl_old)/dble(NNl)
dr(2)=dble(NNg_old)/dble(NNg)
dr(3)=dble(MM_old)/dble(MM)
allocate (c(3,ncomp,0:NN_old+1,0:MM_old*nzl+1),r_old(3,0:NN_old+1),x_old(0:MM_old*nzl+1))
c=0
x_old=0
x_int=0
r_old=0


open(21, file="output/film/radial-profiles.csv")
do i=1,ncomp
    read(21,'(A14)') comp_order(i)
    read(21,*) r_old(1,1:NN_old)
    do j=1,mm_old*nzl
        read(21,*) c(1,i,1:NN_old,j)
    end do
end do
close(21)


open(22, file="output/jet/radial-profiles.csv")
do i=1,ncomp
    read(22,'(A14)') comp_order(i)
    read(22,*) r_old(2,1:NN_old)
    do j=1,mm_old*nzl
        read(22,*) c(2,i,1:NN_old,j)
    end do
end do
close(22)

open(23, file="output/drop/radial-profiles.csv")
do i=1,ncomp
    read(23,'(A14)') comp_order(i)
    read(23,*) r_old(3,1:NN_old)
    do j=1,mm_old*nzl
        read(23,*) c(3,i,1:NN_old,j)
    end do
end do
close(23)


open(24, file="output/profilesgasvap.csv")
read(24,*)
do i=1,MM_old*nzl
    read(24,*) x_old(i)
end do
close(24)

if(((nnl.eq.nnl_old).and.(nng.eq.nng_old)).and.(mm.eq.mm_old)) then
    x_int=c(:,:,1:NN_old,1:MM_old*nzl)
    goto 100
end if

!|\1234567890/|
!|1\12345678/0|
!|21\123456/89|
!|321\1234/678|
!|4321\12/4567|
!|54321\/23456|
!|65432/\12345|
!|7654/21\1234|
!|876/4321\123|
!|98/654321\12|
!|0/87654321\1|
!|/0987654321\|



!call omp_set_num_threads(3)
!$omp parallel private(i,j,k,l,a,b,zr,zx) shared(r_new,x,x_int)
!$omp do
do i=1,3    !film, jet, drop
    do j=1,ncomp    !components

        do k=1,MM*NZL    ! axial nodes
            zx=int(dble(k)*dr(3))
            if (zx.eq.0) then
                a=0
                zx=1
            elseif (zx.eq.MM_old*Nzl) then
                a=0
            else
                a=(x(k)-x_old(zx))/(x_old(zx+1)-x_old(zx))
            end if

            do l=1,NNl
                zr=int(dble(l)*dr(1))
                if ((zr.eq.0))then
                    b=0
                    zr=1
                elseif (zr.eq.NNl_old) then
                    b=0
                else
                    b=(r_new(i,l)-r_old(i,zr))/(r_old(i,zr+1)-r_old(i,zr))
                end if
                x_int(i,j,l,k)=a*b*c(i,j,zr+1,zx+1)+(b-a*b)*c(i,j,zr+1,zx)+(a-a*b)*c(i,j,zr,zx+1)+(1+a*b-a-b)*c(i,j,zr,zx)

                if (isnan(x_int(i,j,l,k))) then
                    call  sleep (1)
                end if

            end do

            do l=1,NNg
                zr=int(dble(l)*dr(2))
                if ((zr.eq.0))then
                    b=0
                    zr=1
                elseif (zr.eq.NNg_old) then
                    b=0
                else
                    b=(r_new(i,nnl+l)-r_old(i,NNl_old+zr))/(r_old(i,NNl_old+zr+1)-r_old(i,NNl_old+zr))
                end if
                x_int(i,j,l+NNl,k)=a*b*c(i,j,zr+1+NNl_old,zx+1)+(b-a*b)*c(i,j,zr+1+NNl_old,zx)+(a-a*b)*c(i,j,zr+NNl_old,zx+1)&
                                    +(1+a*b-a-b)*c(i,j,zr+NNl_old,zx)

                if (isnan(x_int(i,j,l+NNl,k))) then
                    call  sleep (1)
                end if
            end do

        end do

    end do

end do
!$omp end do
!$omp end parallel

100 continue

end subroutine interpol


    subroutine pattern_distribution(holdup,verh)

        implicit none

        real*8, intent(in) :: holdup(3)
        real*8, intent(inout) :: verh(3)
        real*8 :: NC(4)
        real*8 :: ae_calc

        NC(1)=(holdup(1)*verh(2)/(holdup(2)*verh(1))+holdup(1)*verh(3)/(holdup(3)*verh(1))+1D0)**(-1D0)
        NC(2)=NC(1)*holdup(1)/holdup(2)*verh(2)/verh(1)
        NC(3)=NC(1)*holdup(1)/holdup(3)*verh(3)/verh(1)

        NC(4)=sum(NC(1:3))

      !  NC(4)=NC(3)*holdup(3)/(sum(NC(1:3)*holdup))

        verh=NC(1:3)

        ae_calc=NC(1)*(1-delta_film/R_c)+NC(2)*delta_drop/R_c+NC(3)*delta_jet/R_c

        ae_calc=ae_calc

    end subroutine pattern_distribution


end module initialization
