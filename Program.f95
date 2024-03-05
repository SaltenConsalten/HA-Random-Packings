!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                            !!!
!!! CODE für die STUDIENARBEIT !!!
!!!     Frankenstein Code      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program HydrodynamicAnalogy
!**********************************************************************
!   This program contains the hydrodynamic analogy model as it was proposed by Shilkin in 2006.
!   For stability reasons the original TDMA algorithm is replaced by a combination of Cholesky decomposition and
!   the Conjugate Gradient method.
!**********************************************************************


USE omp_lib
USE matrix
USE Algorithms
USE initialization
USE hydrodynamics
USE phys_prop
USE tools
USE global_variables
USE num_solution_mod
!use zbren_int

  implicit none

             !number of components
!  character*14, allocatable :: comp_order(:)       ! Array with component names
                 !number of radial nodes in liquid phase
 ! integer :: NNg                    !number of radial nodes in gaseous phase
 ! integer :: NN         !number of radial nodes
 ! integer :: MM                     !number of axial nodes per liquid mixing length
  real*8 :: zl_film,zl_drop,l_mix,zg                     !liquid mixing length
!  real*8 :: delta_film, delta_drop, delta_jet                       ! liquid film width
!  real*8 :: H                           ! packing height
  REAL*8, ALLOCATABLE :: eta_gas(:),eta_liq(:), e_liq_f(:),e_liq(:), e_gas_f(:),e_gas(:)
  real*8, allocatable :: Conc_film(:,:,:),Conc_drop(:,:,:),Conc_jet(:,:,:), conc(:,:,:)
  real*8, allocatable :: temp_film(:,:,:),temp_drop(:,:,:),temp_jet(:,:,:), Temp(:,:,:)
  real*8, allocatable :: Conc_old_film(:,:,:),Conc_old_drop(:,:,:),Conc_old_jet(:,:,:), Conc_old(:,:,:,:)
  real*8, allocatable :: temp_old_film(:,:,:),temp_old_drop(:,:,:),temp_old_jet(:,:,:), Temp_old(:,:,:,:)
  real*8, allocatable :: dissflux(:,:,:)        !The difference between the output concentration of the nth element to the n+1th element (Ncc,#of liquid mixing points,gas and liquid side)
  real*8, allocatable :: Pr_l(:), Pr_g(:)                          ! Prandtl number
  real*8, allocatable :: Sc_l(:), Sc_g(:)       ! Schmidt number
  real*8, allocatable :: Dist_Coef(:)           ! Henry coefficient of components
  real*8 :: conc_l, conc_g                      ! [kmol/m³] molar densities
  real*8 :: Re_l, Re_g                                  ! Reynolds number
  real*8 :: lambda_l, lambda_g                  ! thermal conductivity
  real*8 :: cp_l, cp_g                          ! heat capacities of liquid and gas
  real*8, allocatable :: H_lv(:)                                ! heat of evaporation + heat of solution
  real*8, allocatable :: D_l(:), D_g(:)            ! diffusion coefficient
  real*8, allocatable :: Ug_film(:),Ug_jet(:),Ug_drop(:),Ul_film(:),Ul_jet(:),Ul_drop(:)
  real*8, allocatable :: turb_ga_film(:,:,:), turb_ga_jet(:,:,:), turb_ga_drop(:,:,:)
  integer :: nzl_all, Nzg                   !number of liquid mixing points
  logical :: recalc                         ! declares if the cholesky decomposition has to be recalculated
!  real*8, allocatable :: C0l(:),C0g(:)         !boundary values for liquid/gas inlets
  real*8 :: dx(2,3) , dh(2,3)                            !axial increment
!  real*8 :: dx(3,2) , dh(3,2)                            !axial increment
  integer :: i,j,k,z,w                                  !control variables
  real*8, allocatable :: T_0_l(:),T_0_g(:)                      ! input temperatures
  real*8 :: T_0lm, T_0gm, T_1lm, T_1gm
  real*8, allocatable :: t_diff(:,:)
  real*8 :: p_0_l,p_0_g                     ! input pressures
  real*8 :: feedl,feedg, FVOLL, FVOLG        ! molar feeds and volume flux
  real*8, allocatable :: IC_l(:,:,:,:),IC_g(:,:,:,:)                        ! integral mean concentrations
  real*8 :: int_ul_film,int_ul_jet,int_ul_drop,int_ug_film,int_ug_jet,int_ug_drop                    ! integral mean velocities
  real*8, allocatable :: x_0(:),y_0(:)                    ! input concentrations
  real*8, allocatable :: rms_c(:,:), rms_t(:,:), Mim(:)               ! root mean square errors of concentration field, molar imbalances
!  real*8, allocatable :: grad_l(:,:,:),grad_g(:,:,:)
  real*8 :: grad_l, grad_g, grad_bp, grad_bp2!, grad_bm, grad_bm2
  real*8 :: grad_l_film, grad_g_film, grad_bp_f, grad_bp2_f!, grad_bm_f, grad_bm2_f
  real*8 :: A_dg_j,A_dg_d,A_dg_f,v_d,v_J,v_F
!  real*8 :: flux_l, flux_g
!  real*8 :: n_in, n_out, diff_flux, balance
!  real*8, allocatable :: incr_x(:)
  logical :: calc=.false.
  logical, dimension (2) :: note=(/.true.,.false./)                   ! log. var. to decide which solution is written
  real*8, allocatable :: val_vec(:,:,:), meanmat(:,:,:), M(:,:,:,:), L(:,:,:,:)
  integer, allocatable :: row(:,:,:), col_ind(:,:,:)
  real*8 :: verh(3)
  real*8, allocatable :: k_reac(:,:)
  real*8, allocatable :: n_A(:,:,:), h_r(:,:,:,:)
  real*8 :: dh_r
  real*8 :: start=0, finish=0
  logical, dimension(2) :: Mnew=(/.false.,.false./)
  real*8, allocatable :: Test(:,:,:)
  real*8 :: flux_test(3)

!  logical :: isocalc=.false.                ! Isotherm calc, only Concentration profiles are calculated
!  logical :: isotherm=.false.               ! only temperature profiles will be calculated
!  logical, dimension(3) :: tcalc=(/.false.,.true.,.true./)                   ! when true, then temp. calc. otherwise conc. calc.

 call control()


  contains


 subroutine control()

implicit none


call input_streams(ncomp,feedl,feedg,T_0_l,T_0_g,p_0_l,p_0_g,x_0,y_0,comp_order)
allocate(D_l(ncomp),D_g(ncomp),H_lv(ncomp),Dist_Coef(ncomp),Sc_l(ncomp),Sc_g(ncomp),Pr_l(ncomp),Pr_g(ncomp))

T_0lm=dble(sum(T_0_l)/ncomp)
T_0gm=dble(sum(T_0_g)/ncomp)
call physical_properties(ncomp,comp_order,T_0lm,T_0gm,p_0_g,x_0,y_0, rho_l,&
                        rho_g,conc_l,conc_g,D_l, D_g,lambda_l, lambda_g,&
                        cp_l,cp_g,mu_l, mu_g,Dist_coef,sigma,H_lv)


forall (w=1:ncomp)
Pr_g(w)=(mu_g*cp_g*(conc_g/rho_g))/lambda_g
Pr_l(w)=(mu_l*cp_l*(conc_l/rho_l))/lambda_l
endforall


forall (w=1:ncomp)
Sc_l(w)=mu_l/(rho_l*D_l(w))
Sc_g(w)=mu_g/(rho_g*D_g(w))
endforall



call flow_parametres(feedl,feedg,mu_l, conc_l, conc_g,sigma, R_c, FVOLL, FVOLG, q_l, q_g, zl_film,zl_Drop,&
                    zg,H,l_mix,verh,ae,rho_l)!,rho_g)
Nzl=nint(H/l_mix)
nzl_all=nint(H/zl_drop)
Nzg=int(H/zg)

call grid(NNl, NNg, MM, eta_liq, eta_gas)

NN=NNl+NNg

allocate(Ul_film(NNl),Ul_jet(NNl),Ul_drop(NNl),Ug_film(nng),Ug_jet(NNg),Ug_drop(NNg),turb_ga_film(NNg+1,ncomp,2)&!,incr_x(MM)&
            ,turb_ga_jet(NNg+1,ncomp,2),turb_ga_drop(NNg+1,ncomp,2),k_reac(3,ncomp))
allocate(Test(nzl*mm,2,3))


k_reac=0D0
!if (reac) k_reac(2,:)=3898.771311D0

do k=1,ncomp
    call kreaction(t_0_l(k),k_reac(1,k))
    k_reac(2,k)=k_reac(1,k)
    k_reac(3,k)=0D0
enddo


!Allocate memory for Matrix operations
if(NN*MM.lt.NMmax) then
    write(*,*) 'Cholesky-Decomposition is used'
    write(*,*)
        allocate(M(NN*MM,NN*MM,3,2),val_vec(5*NN*MM-2*(NN+(MM-1)),3,2),col_ind(NN*MM+1,3,2),&
        row(5*NN*MM-2*(NN+(MM-1)),3,2), meanmat(3,3,2))
  else
    write(*,*) 'ICCG Method is used'
    write(*,*)
        threshold=.true.
        num_it=int(dfloat(nn*mm))
        allocate(M(NN*MM,11,3,2),L(NN*MM+2*NN,11,3,2),val_vec(5*NN*MM-2*(NN+(MM-1)),3,2),&
        col_ind(NN*MM+1,3,2),row(5*NN*MM-2*(NN+(MM-1)),3,2),meanmat(3,3,2))

    L=0D0
end if

M=0D0
val_vec=0D0
row=0D0
col_ind=0D0
meanmat=0D0


    call hydro_film(1-eta_liq,1-eta_gas, rho_g, rho_l, alpha, mu_l, mu_g, R_c, q_g, q_l,delta_film, ul_film,&
                ug_film, int_ul_film,int_ug_film)
    call hydro_drop(eta_liq,eta_gas, rho_g, rho_l, sigma, mu_l, mu_g, R_c, q_g, q_l,delta_drop, ul_drop&
                    ,ug_drop, int_ul_drop,int_ug_drop)
    call hydro_jet(eta_liq,eta_gas, rho_g, rho_l, mu_l, mu_g, R_c, q_g, q_l,delta_jet, ul_jet, ug_jet&
                    ,int_ul_jet,int_ug_jet)




    call trubulence_parameter(R_c, zl_film, rho_l,rho_g, mu_l,mu_g, dx,dh, q_g,re_l,re_g, eta_gas,&
                                turb_ga_film(:,:,1),turb_ga_drop(:,:,1),turb_ga_jet(:,:,1), Sc_g(:),&
                                int_ul_film,int_ul_jet,int_ul_drop,int_ug_film,int_ug_jet,int_ug_drop)

    call trubulence_parameter(R_c, zl_film, rho_l,rho_g, mu_l,mu_g, dx,dh, q_g,re_l,re_g, eta_gas,&
                                turb_ga_film(:,:,2),turb_ga_drop(:,:,2),turb_ga_jet(:,:,2), Pr_g(:),&
                                int_ul_film,int_ul_jet,int_ul_drop,int_ug_film,int_ug_jet,int_ug_drop)



open(23, file='output/check.csv', status='replace')


write(23,*) 'film'
write(23,"(2000(e13.7,:,','))") (1-eta_liq(:))*delta_film+R_c-delta_film,(1-eta_gas(:))*(R_c-delta_film)
write(23,"(2000(e13.7,:,','))") ul_film(:)*int_ul_film,Ug_film(:)*int_ug_film
write(23,"(2000(e13.7,:,','))") ul_film(:)*int_ul_film,turb_ga_film(:,2,1)

write(23,*) 'drop'
write(23,"(2000(e13.7,:,','))") eta_liq(:)*delta_drop,eta_gas(:)*(R_c-delta_drop)+delta_drop
write(23,"(2000(e13.7,:,','))") ul_drop(:)*int_ul_drop,ug_drop(:)*int_ug_drop
write(23,"(2000(e13.7,:,','))") ul_drop(:)*int_ul_drop,turb_ga_drop(:,2,1)

write(23,*) 'jet'
write(23,"(2000(e13.7,:,','))") eta_liq(:)*delta_jet,eta_gas(:)*(R_c-delta_jet)+delta_jet
write(23,"(2000(e13.7,:,','))") ul_jet(:)*int_ul_jet,Ug_jet(:)*int_ug_jet
write(23,"(2000(e13.7,:,','))") ul_jet(:)*int_ul_jet,turb_ga_jet(:,2,1)


close(23)


!allocate(grad_g(ncomp,nzl,MM),grad_l(ncomp,nzl,MM))
allocate(rms_c(ncomp,3), rms_t(ncomp,3), Mim(ncomp))!, C0l(NNl),C0g(NNg))
allocate (dissflux(ncomp,Nzl,2))
dissflux=0D0
allocate (Conc_film(NN*MM,Nzl,ncomp),Conc_drop(NN*MM,Nzl,ncomp),Conc_jet(NN*MM,Nzl,ncomp),Temp(NN*MM,Nzl,ncomp),&
            temp_old_film(NN*MM,Nzl,ncomp),temp_old_drop(NN*MM,Nzl,ncomp),temp_old_jet(NN*MM,Nzl,ncomp),&
            temp_film(NN*MM,Nzl,ncomp),temp_drop(NN*MM,Nzl,ncomp),temp_jet(NN*MM,Nzl,ncomp),&
            Conc_old_film(NN*MM,Nzl,ncomp),Conc_old_drop(NN*MM,Nzl,ncomp),Conc_old_jet(NN*MM,Nzl,ncomp),Temp_old(3,ncomp,NN*MM,Nzl)&
            ,IC_l(nzl*MM,3,ncomp,2),IC_g(nzl*MM,3,ncomp,2),Conc_old(3,ncomp,NN,MM*nzl),Conc(NN*MM,Nzl,ncomp),&
            n_A(Nzl*MM,3,ncomp), h_r(NN*MM,Nzl,3,ncomp))
allocate (e_liq_f(0:NNl),e_liq(0:NNl),e_gas_f(0:NNg+1),e_gas(0:NNg+1),t_diff(3,4))

!#######################################################################################
!       Initialisation
!#######################################################################################

    h_r=0
    n_A=0

    call interpol((1-eta_liq(:))*delta_film+R_c-delta_film,(1-eta_gas(:))*(R_c-delta_film),&
                    eta_liq(:)*delta_jet,eta_gas(:)*(R_c-delta_jet)+delta_jet,eta_liq(:)*delta_drop,&
                    eta_gas(:)*(R_c-delta_drop)+delta_drop,H,Conc_old)!,note(1))
    IC_l=0
    IC_g=0

    Conc_old_drop(:,:,:)=0D0
    Conc_old_film(:,:,:)=0D0
    Conc_old_jet(:,:,:)=0D0
    do k=1,ncomp
    do i=1,nzl
        do j=1,MM
            Conc_drop(nn*(j-1)+1:nn*(j-1)+nnl,i,k)=conc_old(3,k,1:NNl,MM*(i-1)+j)/conc_l
            Conc_film(nn*(j-1)+1:nn*(j-1)+nnl,i,k)=Conc_old(1,k,1:NNl,MM*(i-1)+j)/conc_l
            Conc_jet(nn*(j-1)+1:nn*(j-1)+nnl,i,k)=conc_old(2,k,1:NNl,MM*(i-1)+j)/conc_l
            Conc_drop(nn*(j-1)+NNl+1:NN*j,i,k)=conc_old(3,k,nnl+1:NN,MM*(i-1)+j)/conc_g
            Conc_film(nn*(j-1)+NNl+1:NN*j,i,k)=Conc_old(1,k,nnl+1:NN,MM*(i-1)+j)/conc_g
            Conc_jet(nn*(j-1)+nnl+1:NN*j,i,k)=conc_old(2,k,nnl+1:NN,MM*(i-1)+j)/conc_g
        end do
    end do
    end do
    deallocate (Conc_old)


!    call interpol((1-eta_liq(:))*delta_film+R_c-delta_film,(1-eta_gas(:))*(R_c-delta_film),&
!                    eta_liq(:)*delta_jet,eta_gas(:)*(R_c-delta_jet)+delta_jet,eta_liq(:)*delta_drop,&
!                    eta_gas(:)*(R_c-delta_drop)+delta_drop,H,temp_old,note(2))
!    IC_l=0
!    IC_g=0
!
    temp_old_drop(:,:,:)=0D0
    temp_old_film(:,:,:)=0D0
    temp_old_jet(:,:,:)=0D0
!    do k=1,ncomp
!    do i=1,nzl
!        do j=1,MM
!    temp_drop(k,nn*(j-1)+1:NNl*j,i)=Temp_old(3,k,1:NNl,MM*(i-1)+j)/T_0_l
!    temp_film(k,nn*(j-1)+1:NNl*j,i)=temp_old(1,k,1:NNl,MM*(i-1)+j)/T_0_l
!    temp_jet(k,nn*(j-1)+1:NNl*j,i)=temp_old(2,k,1:NNl,MM*(i-1)+j)/T_0_l
!    temp_drop(k,nn*(j-1)+NNl+1:NN*j,i)=temp_old(3,k,nnl+1:NN,MM*(i-1)+j)/T_0_g
!    temp_film(k,nn*(j-1)+NNl+1:NN*j,i)=temp_old(1,k,nnl+1:NN,MM*(i-1)+j)/T_0_g
!    temp_jet(k,nn*(j-1)+nnl+1:NN*j,i)=temp_old(2,k,nnl+1:NN,MM*(i-1)+j)/T_0_g
!        end do
!    end do
!    end do
!
    deallocate (temp_old)


!goto 100
!Conc_drop=0
!Conc_film=0
!Conc_jet=0
!###############################################################################################
    !Start Iterations
!###############################################################################################

rms_c=1D0
rms_t=1D0
Mim=1D0

do z=1,3
    t_diff(z,1)=T_0_l(k)!(IC_l(nzl*mm,z,k,1)-IC_l(nzl*mm,z,k,2))
    t_diff(z,2)=T_0_g(k)!(IC_g(1,z,k,1)-IC_g(1,z,k,2))
end do

!    lambda_l=mu_l/(Pr_l(2)*rho_l)
!    lambda_g=mu_g/(Pr_g(2)*rho_g)

call pattern_distribution((/pi*(R_c**2-(R_c-delta_film)**2)*h,q_l*h/(abs(ul_drop(1)*int_ul_drop)),pi*delta_jet**2*h/),verh)!,R_c)
ae=1

open(666, file='output/std-out.csv', status='replace')
call cpu_time(start)
do k=2,ncomp-1
    write(*,*) k, "of ", ncomp, "components are calculated"

    j=1

    do while((DMAX1(maxval(abs(rms_t(k,:))),maxval(abs(rms_c(k,:)))).gt.1D-6).and.(j.lt.itmax))
        write(*,*)
        write(*,*) 'outer loop iteration: ', j
        write(*,*) 'the current RMS for concentration are', abs(RMS_C(k,:))
        write(666,"(2000(e13.7,:,','))") abs(RMS_C(k,:))
        write(*,*) 'Imbalance', Mim(:)
            Conc_old_film(:,:,:)=Conc_film(:,:,:)
            Conc_old_drop(:,:,:)=Conc_drop(:,:,:)
            Conc_old_jet(:,:,:)=Conc_jet(:,:,:)

        write(*,*) 'the current RMS for temperature are', abs(RMS_t(k,:))
        write(666,"(2000(e13.7,:,','))") abs(RMS_t(k,:))

        if(maxval(abs(rms_t(k,:))).lt.1D-4) then
!            if ((mod(j,10).eq.0).or.tcalc(3)) then
                do i=1,nzl
                    do z=1,MM
                        call mean_int(temp_film((z-1)*NN+1:(z-1)*NN+NNl,i,k)*t_0_l(k),&!*conc_l,&
                                    temp_film((z-1)*NN+NNl+1:z*NN,i,k)*t_0_g(k),&!*conc_g,&
                                    ul_film,ug_film,(1-eta_liq)*delta_film+R_c-delta_film,(1-eta_gas)*(R_c-delta_film),&
                                    IC_l(((i-1)*MM)+z,1,k,2),IC_g(((i-1)*MM)+z,1,k,2))
                        call mean_int(temp_drop((z-1)*NN+1:(z-1)*NN+NNl,i,k)*t_0_l(k),&!*conc_l,&
                                    temp_drop((z-1)*NN+NNl+1:z*NN,i,k)*t_0_g(k),&!*conc_g,&
                                    ul_drop,ug_drop,eta_liq*delta_drop,eta_gas*(R_c-delta_drop)+delta_drop,&
                                    IC_l(((i-1)*MM)+z,2,k,2),IC_g(((i-1)*MM)+z,2,k,2))
                        call mean_int(temp_jet((z-1)*NN+1:(z-1)*NN+NNl,i,k)*t_0_l(k),&!*conc_l,&
                                    temp_jet((z-1)*NN+NNl+1:z*NN,i,k)*t_0_g(k),&!*conc_g,&
                                    ul_jet,ug_jet,eta_liq*delta_jet,eta_gas*(R_c-delta_jet)+delta_jet,&
                                    IC_l(((i-1)*MM)+z,3,k,2),IC_g(((i-1)*MM)+z,3,k,2))


                        call mean_int(temp_old_film((z-1)*NN+1:(z-1)*NN+NNl,i,k)*t_0_l(k),&!*conc_l,&
                                    temp_old_film((z-1)*NN+NNl+1:z*NN,i,k)*t_0_g(k),&!*conc_g,&
                                    ul_film,ug_film,(1-eta_liq)*delta_film+R_c-delta_film,(1-eta_gas)*(R_c-delta_film),&
                                    IC_l(((i-1)*MM)+z,1,k,1),IC_g(((i-1)*MM)+z,1,k,1))
                        call mean_int(temp_old_drop((z-1)*NN+1:(z-1)*NN+NNl,i,k)*t_0_l(k),&!*conc_l,&
                                    temp_old_drop((z-1)*NN+NNl+1:z*NN,i,k)*t_0_g(k),&!*conc_g,&
                                    ul_drop,ug_drop,eta_liq*delta_drop,eta_gas*(R_c-delta_drop)+delta_drop,&
                                    IC_l(((i-1)*MM)+z,2,k,1),IC_g(((i-1)*MM)+z,2,k,1))
                        call mean_int(temp_old_jet((z-1)*NN+1:(z-1)*NN+NNl,i,k)*t_0_l(k),&!*conc_l,&
                                    temp_old_jet((z-1)*NN+NNl+1:z*NN,i,k)*t_0_g(k),&!*conc_g,&
                                    ul_jet,ug_jet,eta_liq*delta_jet,eta_gas*(R_c-delta_jet)+delta_jet,&
                                    IC_l(((i-1)*MM)+z,3,k,1),IC_g(((i-1)*MM)+z,3,k,1))

                    end do
                end do

                do z=1,3
                    t_diff(z,3)=t_diff(z,1)-IC_l(nzl*mm,z,k,2)
                    t_diff(z,4)=t_diff(z,2)-IC_g(1,z,k,2)
                end do

        if ((maxval(abs(t_diff(:,3:4))).gt.0.5D0).or.tcalc(3)) then

                do z=1,3
                    t_diff(z,1)=IC_l(nzl*mm,z,k,1)
                    t_diff(z,2)=IC_g(1,z,k,1)
                end do



!                T_1lm=verh(1)*(IC_l(2,k,1,1)+(((IC_l(2,k,1,1)-IC_g(2,k,1,1))-(IC_l(2,k,mm*nzl,1)-IC_g(2,k,mm*nzl,1)))/&
!                        log((IC_l(2,k,1,1)-IC_g(2,k,1,1))/(IC_l(2,k,mm*nzl,1)-IC_g(2,k,mm*nzl,1)))))+&
!                      verh(2)*(IC_l(2,k,1,2)+(((IC_l(2,k,1,2)-IC_g(2,k,1,2))-(IC_l(2,k,mm*nzl,2)-IC_g(2,k,mm*nzl,2)))/&
!                        log((IC_l(2,k,1,2)-IC_g(2,k,1,2))/(IC_l(2,k,mm*nzl,2)-IC_g(2,k,mm*nzl,2)))))+&
!                      verh(3)*(IC_l(2,k,1,3)+(((IC_l(2,k,1,3)-IC_g(2,k,1,3))-(IC_l(2,k,mm*nzl,3)-IC_g(2,k,mm*nzl,3)))/&
!                        log((IC_l(2,k,1,3)-IC_g(2,k,1,3))/(IC_l(2,k,mm*nzl,3)-IC_g(2,k,mm*nzl,3)))))
!
!                T_1gm=verh(1)*(IC_g(2,k,mm*nzl,1)+(((IC_l(2,k,1,1)-IC_g(2,k,1,1))-(IC_l(2,k,mm*nzl,1)-IC_g(2,k,mm*nzl,1)))/&
!                        log((IC_l(2,k,1,1)-IC_g(2,k,1,1))/(IC_l(2,k,mm*nzl,1)-IC_g(2,k,mm*nzl,1)))))+&
!                      verh(2)*(IC_g(2,k,mm*nzl,2)+(((IC_l(2,k,1,2)-IC_g(2,k,1,2))-(IC_l(2,k,mm*nzl,2)-IC_g(2,k,mm*nzl,2)))/&
!                        log((IC_l(2,k,1,2)-IC_g(2,k,1,2))/(IC_l(2,k,mm*nzl,2)-IC_g(2,k,mm*nzl,2)))))+&
!                      verh(3)*(IC_g(2,k,mm*nzl,3)+(((IC_l(2,k,1,3)-IC_g(2,k,1,3))-(IC_l(2,k,mm*nzl,3)-IC_g(2,k,mm*nzl,3)))/&
!                        log((IC_l(2,k,1,3)-IC_g(2,k,1,3))/(IC_l(2,k,mm*nzl,3)-IC_g(2,k,mm*nzl,3)))))

                T_1lm=(verh(1)*(IC_l(1,1,k,2)+IC_l(mm*nzl,1,k,2))/2D0&
                        +verh(2)*(IC_l(1,2,k,2)+IC_l(mm*nzl,2,k,2))/2D0&
                        +verh(3)*(IC_l(1,3,k,2)+IC_l(mm*nzl,3,k,2))/2D0)!(verh(1)*IC_l(2,k,nzl,1)+verh(2)*IC_l(2,k,nzl,2)+verh(3)*IC_l(2,k,nzl,3))
                T_1gm=(verh(1)*(IC_g(1,1,k,2)+IC_g(mm*nzl,1,k,2))/2D0&
                        +verh(2)*(IC_g(1,2,k,2)+IC_g(mm*nzl,2,k,2))/2D0&
                        +verh(3)*(IC_g(1,3,k,2)+IC_g(mm*nzl,3,k,2))/2D0)!(verh(1)*IC_g(2,k,nzl,1)+verh(2)*IC_g(2,k,nzl,2)+verh(3)*IC_g(2,k,nzl,3))

                        write(*,*)
                        write(*,*) 'neue mittlere Temp: ', T_1lm, T_1gm
!                        write(*,*)
!                do z=1,3
!                    t_diff(z,1,1,1)=IC_l()
!                enddo

                call physical_properties(ncomp,comp_order,T_1lm,T_1gm,p_0_g,x_0,y_0, rho_l,&
                                        rho_g,conc_l,conc_g,D_l, D_g,lambda_l, lambda_g,&
                                        cp_l,cp_g,mu_l, mu_g,Dist_coef,sigma,H_lv)


                call kreaction(T_1lm,k_reac(1,k))
                k_reac(2,k)=k_reac(1,k)
                k_reac(3,k)=0D0


                forall (w=1:ncomp)
                Pr_g(w)=mu_g*cp_g*conc_g/(rho_g*lambda_g)
                Pr_l(w)=mu_l*cp_l*conc_l/(rho_l*lambda_l)
                endforall


                forall (w=1:ncomp)
                Sc_l(w)=mu_l/(rho_l*D_l(w))
                Sc_g(w)=mu_g/(rho_g*D_g(w))
                endforall


                call hydro_film(1-eta_liq,1-eta_gas, rho_g, rho_l, alpha, mu_l, mu_g, R_c, q_g, q_l,delta_film, ul_film,&
                                ug_film, int_ul_film,int_ug_film)
                call hydro_drop(eta_liq,eta_gas, rho_g, rho_l, sigma, mu_l, mu_g, R_c, q_g, q_l,delta_drop, ul_drop&
                                ,ug_drop, int_ul_drop,int_ug_drop)
                call hydro_jet(eta_liq,eta_gas, rho_g, rho_l, mu_l, mu_g, R_c, q_g, q_l,delta_jet, ul_jet, ug_jet&
                                ,int_ul_jet,int_ug_jet)


                call trubulence_parameter(R_c, zl_film, rho_l,rho_g, mu_l,mu_g, dx,dh, q_g,re_l,re_g, eta_gas,&
                                turb_ga_film(:,:,1),turb_ga_drop(:,:,1),turb_ga_jet(:,:,1), Sc_g(:),&
                                int_ul_film,int_ul_jet,int_ul_drop,int_ug_film,int_ug_jet,int_ug_drop)
                call trubulence_parameter(R_c, zl_film, rho_l,rho_g, mu_l,mu_g, dx,dh, q_g,re_l,re_g, eta_gas,&
                                turb_ga_film(:,:,2),turb_ga_drop(:,:,2),turb_ga_jet(:,:,2), Pr_g(:),&
                                int_ul_film,int_ul_jet,int_ul_drop,int_ug_film,int_ug_jet,int_ug_drop)

                    tcalc(3)=.false.    !erste Neuberechnung der Properties erzwingen
!                    tcalc(4)=.true.     !wenn true matrix neu aufstellen
                    Mnew=.true.
            endif
        endif



!            temp_old(:,:,:)=temp(:,:,:)
            temp_old_film(:,:,:)=temp_film(:,:,:)
            temp_old_drop(:,:,:)=temp_drop(:,:,:)
            temp_old_jet(:,:,:)=temp_jet(:,:,:)

            if(maxval(rms_c(k,:)).lt.1D-4) then
                calc=.false.
            end if
            !Algorithmus

!    !! calculation of concentration profiles

!        if(maxval(abs(rms_c(k,:))).gt.1D-6) then

            call calculation(j,eta_gas(:),eta_liq(:), conc_film(:,:,k),conc_drop(:,:,k),conc_jet(:,:,k), R_c,&
                        conc_old_film(:,:,k),conc_old_drop(:,:,k),conc_old_jet(:,:,k),Sc_l(k), Dist_coef(k),conc_l,conc_g,&
                        D_l(k),D_g(k),delta_film,delta_drop,delta_jet,Ug_film,Ug_drop,Ug_jet,Ul_film,Ul_drop,Ul_jet,&
                        turb_ga_film(:,k,1),turb_ga_drop(:,k,1),turb_ga_jet(:,k,1),nzl_all,Nzg, dx,dh, IC_l(:,:,k,1),&
                        IC_g(:,:,k,1),x_0(k),y_0(k),verh,k_reac(:,k),val_vec(:,:,1),meanmat(:,:,1),row(:,:,1),&
                        col_ind(:,:,1),M(:,:,:,1),L(:,:,:,1),n_A(:,:,k),h_r(:,:,:,k),tcalc(1),Mnew(1))

    !!!  calculate RMS and imbalances

                    call rms_err(eta_gas, eta_liq, conc_film, conc_drop, conc_jet, conc, conc_l,conc_g,IC_l(:,:,:,1),&
                            IC_g(:,:,:,1),rms_c, conc_old_film, conc_old_drop, conc_old_jet, ug_drop, ul_drop,&!
!                            ug_film, ug_jet,ul_film, ul_jet,
                            k, fvoll, fvolg, tcalc(1), x_0, y_0, Mim)

!        endif

                    e_liq_f(1:nnl)=1-eta_liq(:)
                    e_liq_f(0)=1
                    e_gas_f(1:nng)=1-eta_gas(:)
                    e_gas_f(0)=1
                    e_gas_f(nng+1)=0
                    e_liq(1:nnl)=eta_liq(:)
                    e_liq(0)=0
                    e_gas(1:nng)=eta_gas(:)
                    e_gas(0)=0
                    e_gas(nng+1)=1

                    grad_l_film=((e_liq_f(NNl-1)-e_liq_f(NNl))**2-(e_liq_f(NNl-2)-e_liq_f(NNl))**2)/&
                            ((e_liq_f(NNl-2)-e_liq_f(NNl-1))*(e_liq_f(NNl-2)-e_liq_f(NNl))*(e_liq_f(NNl-1)-e_liq_f(NNl)))
                    grad_l=((e_liq(NNl-1)-e_liq(NNl))**2-(e_liq(NNl-2)-e_liq(NNl))**2)/&
                            ((e_liq(NNl-2)-e_liq(NNl-1))*(e_liq(NNl-2)-e_liq(NNl))*(e_liq(NNl-1)-e_liq(NNl)))

                    grad_g=((e_gas(1)-e_gas(0))**2-(e_gas(2)-e_gas(0))**2)/&!
                            ((e_gas(2)-e_gas(1))*(e_gas(2)-e_gas(0))*(e_gas(1)-e_gas(0)))
                    grad_g_film=((e_gas_f(1)-e_gas_f(0))**2-(e_gas_f(2)-e_gas_f(0))**2)/&!
                                ((e_gas_f(2)-e_gas_f(1))*(e_gas_f(2)-e_gas_f(0))*(e_gas_f(1)-e_gas_f(0)))

                    grad_bp=(e_gas(2)-e_gas(0))/((e_gas(2)-e_gas(1))*(e_gas(1)-e_gas(0)))
                    grad_bp_f=(e_gas_f(2)-e_gas_f(0))/((e_gas_f(2)-e_gas_f(1))*(e_gas_f(1)-e_gas_f(0)))

                    grad_bp2=(e_gas(1)-e_gas(0))/((e_gas(2)-e_gas(1))*(e_gas(2)-e_gas(0)))
                    grad_bp2_f=(e_gas_f(1)-e_gas_f(0))/((e_gas_f(2)-e_gas_f(1))*(e_gas_f(2)-e_gas_f(0)))

!                    grad_bm=(e_liq(NNl-2)-e_liq(NNl))/((e_liq(NNl-2)-e_liq(NNl-1))*(e_liq(NNl-1)-e_liq(NNl)))
!                    grad_bm_f=(e_liq_f(NNl-2)-e_liq_f(NNl))/((e_liq_f(NNl-2)-e_liq_f(NNl-1))*(e_liq_f(NNl-1)-e_liq_f(NNl)))
!
!                    grad_bm2=(e_liq(NNl-1)-e_liq(NNl))/((e_liq(NNl-2)-e_liq(NNl-1))*(e_liq(NNl-2)-e_liq(NNl)))
!                    grad_bm2_f=(e_liq_f(NNl-1)-e_liq_f(NNl))/((e_liq_f(NNl-2)-e_liq_f(NNl-1))*(e_liq_f(NNl-2)-e_liq_f(NNl)))
!
!write(*,*) grad_L, grad_l_film
!write(*,*) grad_g, grad_g_film

                    A_dg_j=2*pi*dh(2,2)*delta_jet
                    v_j=(delta_jet**2-(delta_jet*(eta_liq(nnl-1)))**2)*pi*dh(2,2)
                    A_dg_d=2*pi*dh(2,2)*delta_drop
                    v_d=(delta_drop**2-(delta_drop*(eta_liq(nnl-1)))**2)*pi*dh(2,2)
                    A_dg_f=-2*pi*dh(2,2)*(R_c-delta_film)
                    v_f=pi*((R_c-delta_film+delta_film*(1D0-e_liq(nnl-1)))**2-(R_c-delta_film)**2)*dh(2,2)
                    dh_R=-108.68

                    do i=1,nzl
                        do z=1,MM
                            n_A(((i-1)*MM)+z,1,k)=A_dg_f*h_lv(k)*D_g(k)/(R_c-delta_film)*&
                                                    (conc_film((z-1)*NN+NNl,i,k)*conc_l*grad_g_film*Dist_coef(k))&
                                                    +A_dg_f*h_lv(k)*D_g(k)/(R_c-delta_film)*(&
                                                    +conc_film((z-1)*nn+nnl+1,i,k)*conc_g*grad_bp_f&
                                                    -conc_film((z-1)*nn+nnl+2,i,k)*conc_g*grad_bp2_f)&
                                                    -dh_r*k_reac(k,2)*conc_film((z-1)*NN+nnl,i,k)*conc_l*v_F



                            n_A(((i-1)*MM)+z,2,k)=A_dg_d*h_lv(k)*D_g(k)/(R_c-delta_drop)*&
                                                    (conc_drop((z-1)*NN+NNl,i,k)*conc_l*grad_g*Dist_coef(k))&
                                                    +A_dg_d*h_lv(k)*D_g(k)/(R_c-delta_drop)*(&
                                                    +conc_drop((z-1)*nn+nnl+1,i,k)*conc_g*grad_bp&
                                                    -conc_drop((z-1)*nn+nnl+2,i,k)*conc_g*grad_bp2)&
                                                    +dh_r*k_reac(k,2)*Conc_drop((z-1)*NN+nnl,i,k)*conc_l*v_D

                            n_A(((i-1)*MM)+z,3,k)=A_dg_j*h_lv(k)*D_g(k)/(R_c-delta_jet)*&
                                                    (conc_jet((z-1)*NN+NNl,i,k)*conc_l*grad_g*Dist_coef(k))&
                                                    +A_dg_j*h_lv(k)*D_g(k)/(R_c-delta_jet)*(&
                                                    +conc_jet((z-1)*nn+nnl+1,i,k)*conc_g*grad_bp&
                                                    -conc_jet((z-1)*nn+nnl+2,i,k)*conc_g*grad_bp2)&
                                                    +dh_r*k_reac(k,2)*conc_jet((z-1)*NN+nnl,i,k)*conc_l*v_J


                        end do
                    end do
!
!

                   !*1D3                                ![J/kmol]

                    do i=1,nzl
                        do z=1,MM
                            do w=1,NNl

                                h_r((z-1)*NN+w,i,1,k)=dh_r*k_reac(k,2)*conc_film((z-1)*NN+w,i,k)*conc_l*delta_film**2/&
                                                    (cp_l*mu_l*T_0_l(k))
!
                                h_r((z-1)*NN+w,i,2,k)=dh_r*k_reac(k,2)*conc_drop((z-1)*NN+w,i,k)*conc_l*delta_drop**2/&
                                                    (cp_l*mu_l*T_0_l(k))
!
                                h_r((z-1)*NN+w,i,3,k)=dh_r*k_reac(k,2)*conc_jet((z-1)*NN+w,i,k)*conc_l*delta_jet**2/&
                                                    (cp_l*mu_l*T_0_l(k))
!
                            end do
                        end do
                    end do

!        if(maxval(abs(rms_c(k,:))).lt.1D-6) then

    !! calculation of temperatur profiles
h_r(:,:,:,:)=0
n_A(:,:,:)=0

            call calculation(j,eta_gas,eta_liq, temp_film(:,:,k),temp_drop(:,:,k),temp_jet(:,:,k), R_c,&
                        temp_old_film(:,:,k),temp_old_drop(:,:,k),temp_old_jet(:,:,k),Pr_l(k), 1D0, T_0_l(k),T_0_g(k),&
                        lambda_l,lambda_g ,delta_film,delta_drop,delta_jet,Ug_film,Ug_drop,Ug_jet,Ul_film,Ul_drop,Ul_jet,&
                        turb_ga_film(:,k,2),turb_ga_drop(:,k,2),turb_ga_jet(:,k,2),nzl_all,Nzg, dx,dh,IC_l(:,:,k,2),&
                        IC_g(:,:,k,2), T_0_l(k)/T_0_l(k),T_0_g(k)/T_0_g(k),verh,k_reac(3,k),val_vec(:,:,2),meanmat(:,:,2),&
                        row(:,:,2),col_ind(:,:,2),M(:,:,:,2),L(:,:,:,2),n_A(:,:,k),h_r(:,:,:,k),tcalc(2),Mnew(2))

    !!!   calculate RMS and imbalances

                    call rms_err(eta_gas, eta_liq, temp_film,temp_drop,temp_jet,temp, conc_l,conc_g, IC_l(:,:,:,2),&
                            IC_g(:,:,:,2), rms_t, temp_old_film, temp_old_drop, temp_old_jet, ug_drop, ul_drop,&!
!                            ug_film, ug_jet,ul_film, ul_jet,
                            k, fvoll, fvolg, tcalc(2), t_0_l, t_0_g, Mim)


            num_it=num_it+int(dfloat(nn*mm/20))

do i=1,nzl
do z=1,mm
!####################################################################################################################################
!#####################                       FILM                                                            ########################
!####################################################################################################################################

Test((i-1)*mm+z,1,1)=T_0_l(1)*lambda_l/(delta_film)*(-&
Temp_film(nn*(z-1)+nnl-2,i,2)*((e_liq_f(NNl-1)-e_liq_f(NNl))/((e_liq_f(NNl-2)-e_liq_f(NNl-1))*(e_liq_f(NNl-2)-e_liq_f(NNl))))+&
Temp_film(nn*(z-1)+nnl-1,i,2)*((e_liq_f(NNl-2)-e_liq_f(NNl))/((e_liq_f(NNl-2)-e_liq_f(NNl-1))*(e_liq_f(NNl-1)-e_liq_f(NNl))))+&
Temp_film(nn*(z-1)+nnl,i,2)*(((e_liq_f(NNl-1)-e_liq_f(NNl))**2-(e_liq_f(NNl-2)-e_liq_f(NNl))**2)/((e_liq_f(NNl-2)-e_liq_f(NNl-1))*&
                                    (e_liq_f(NNl-2)-e_liq_f(NNl))*(e_liq_f(NNl-1)-e_liq_f(NNl)))))

Test((i-1)*mm+z,2,1)=T_0_g(1)*lambda_g/(R_c-delta_film)*(-&
Temp_film(nn*(z-1)+nnl+2,i,2)*((e_gas_f(1)-e_gas_f(0))/((e_gas_f(2)-e_gas_f(1))*(e_gas_f(2)-e_gas_f(0))))+&
Temp_film(nn*(z-1)+nnl+1,i,2)*((e_gas_f(2)-e_gas_f(0))/((e_gas_f(2)-e_gas_f(1))*(e_gas_f(1)-e_gas_f(0))))+T_0_l(1)/T_0_g(1)*&
Temp_film(nn*(z-1)+nnl,i,2)*(((e_gas_f(1)-e_gas_f(0))**2-(e_gas_f(2)-e_gas_f(0))**2)/&
            ((e_gas_f(2)-e_gas_f(1))*(e_gas_f(2)-e_gas_f(0))*(e_gas_f(1)-e_gas_f(0)))))

!####################################################################################################################################
!#####################                        JET                                                            ########################
!####################################################################################################################################
Test((i-1)*mm+z,1,2)=T_0_l(1)*lambda_l/(delta_jet)*(-&
Temp_Jet(nn*(z-1)+nnl-2,i,2)*((e_liq(NNl-1)-e_liq(NNl))/((e_liq(NNl-2)-e_liq(NNl-1))*(e_liq(NNl-2)-e_liq(NNl))))+&
Temp_Jet(nn*(z-1)+nnl-1,i,2)*((e_liq(NNl-2)-e_liq(NNl))/((e_liq(NNl-2)-e_liq(NNl-1))*(e_liq(NNl-1)-e_liq(NNl))))+&
Temp_Jet(nn*(z-1)+nnl,i,2)*(((e_liq(NNl-1)-e_liq(NNl))**2-(e_liq(NNl-2)-e_liq(NNl))**2)/((e_liq(NNl-2)-e_liq(NNl-1))*&
                                    (e_liq(NNl-2)-e_liq(NNl))*(e_liq(NNl-1)-e_liq(NNl)))))

Test((i-1)*mm+z,2,2)=T_0_g(1)*lambda_g/(R_c-delta_jet)*(-&
Temp_Jet(nn*(z-1)+nnl+2,i,2)*((e_gas(1)-e_gas(0))/((e_gas(2)-e_gas(1))*(e_gas(2)-e_gas(0))))+&
Temp_Jet(nn*(z-1)+nnl+1,i,2)*((e_gas(2)-e_gas(0))/((e_gas(2)-e_gas(1))*(e_gas(1)-e_gas(0))))+T_0_l(1)/T_0_g(1)*&
Temp_Jet(nn*(z-1)+nnl,i,2)*(((e_gas(1)-e_gas(0))**2-(e_gas(2)-e_gas(0))**2)/&
            ((e_gas(2)-e_gas(1))*(e_gas(2)-e_gas(0))*(e_gas(1)-e_gas(0)))))
!####################################################################################################################################
!#####################                        DROP                                                           ########################
!####################################################################################################################################
Test((i-1)*mm+z,1,3)=T_0_l(1)*lambda_l/(delta_drop)*(-&
Temp_Drop(nn*(z-1)+nnl-2,i,2)*((e_liq(NNl-1)-e_liq(NNl))/((e_liq(NNl-2)-e_liq(NNl-1))*(e_liq(NNl-2)-e_liq(NNl))))+&
Temp_Drop(nn*(z-1)+nnl-1,i,2)*((e_liq(NNl-2)-e_liq(NNl))/((e_liq(NNl-2)-e_liq(NNl-1))*(e_liq(NNl-1)-e_liq(NNl))))+&
Temp_Drop(nn*(z-1)+nnl,i,2)*(((e_liq(NNl-1)-e_liq(NNl))**2-(e_liq(NNl-2)-e_liq(NNl))**2)/((e_liq(NNl-2)-e_liq(NNl-1))*&
                                    (e_liq(NNl-2)-e_liq(NNl))*(e_liq(NNl-1)-e_liq(NNl)))))

Test((i-1)*mm+z,2,3)=T_0_g(1)*lambda_g/(R_c-delta_drop)*(-&
Temp_Drop(nn*(z-1)+nnl+2,i,2)*((e_gas(1)-e_gas(0))/((e_gas(2)-e_gas(1))*(e_gas(2)-e_gas(0))))+&
Temp_Drop(nn*(z-1)+nnl+1,i,2)*((e_gas(2)-e_gas(0))/((e_gas(2)-e_gas(1))*(e_gas(1)-e_gas(0))))+T_0_l(1)/T_0_g(1)*&
Temp_Drop(nn*(z-1)+nnl,i,2)*(((e_gas(1)-e_gas(0))**2-(e_gas(2)-e_gas(0))**2)/&
            ((e_gas(2)-e_gas(1))*(e_gas(2)-e_gas(0))*(e_gas(1)-e_gas(0)))))
end do
end do
flux_test(1)=a_dg_f*abs(sum(test(:,1,1)))!-a_dg_f*abs(sum(test(:,2,1))))/(a_dg_f*abs(sum(test(:,1,1))))
flux_test(2)=a_dg_j*abs(sum(test(:,1,2)))!-a_dg_j*abs(sum(test(:,2,2))))/(a_dg_j*abs(sum(test(:,1,2))))
flux_test(3)=a_dg_d*abs(sum(test(:,1,3)))!-a_dg_d*abs(sum(test(:,2,3))))/(a_dg_d*abs(sum(test(:,1,3))))
            j=j

            call mean_int(T_0_l(2)*Temp_film((MM-1)*NN+1:(MM-1)*NN+NNl,nzl,2),T_0_g(2)*Temp_film(NNL+1:NN,1,2),ul_film,Ug_film,&
                    (1-eta_liq)*delta_film+R_c-delta_film,(1-eta_gas)*(R_c-delta_film),IC_l(nzl,1,2,2),IC_g(1,1,2,2))


            mim(1)=abs(q_l*cp_l*conc_l*(T_0_l(2)-IC_l(nzl,1,2,2)))/flux_test(1)!abs(q_g*cp_g*conc_g*(T_0_g(2)-IC_g(1,1,2,2)))/flux_test(1)!
!            1-(abs(q_l*cp_l*conc_l/rho_l*T_0_l(2))+abs(q_g*cp_g*conc_g/rho_g*T_0_g(2)))/&
!                    (abs(q_l*cp_l*conc_l/rho_l*IC_l(nzl,1,2,2))+abs(q_g*cp_g*conc_g/rho_g*IC_g(1,1,2,2)))


            call mean_int(T_0_l(2)*Temp_drop((MM-1)*NN+1:(MM-1)*NN+NNl,nzl,2),T_0_g(2)*Temp_drop(NNL+1:NN,1,2),ul_drop,ug_drop,&
                    eta_liq*delta_drop,eta_gas*(R_c-delta_drop)+delta_drop,IC_l(nzl,2,2,2),IC_g(1,2,2,2))
            mim(2)=abs(q_l*cp_l*conc_l*(T_0_l(2)-IC_l(nzl,2,2,2)))/flux_test(2)!abs(q_g*cp_g*conc_g*(T_0_g(2)-IC_g(1,2,2,2)))/flux_test(2)!
!            1-(abs(q_l*cp_l*conc_l/rho_l*T_0_l(2))+abs(q_g*cp_g*conc_g/rho_g*T_0_g(2)))/&
!                    (abs(q_l*cp_l*conc_l/rho_l*IC_l(nzl,2,2,2))+abs(q_g*cp_g*conc_g/rho_g*IC_g(1,2,2,2)))


            call mean_int(T_0_l(2)*Temp_jet((MM-1)*NN+1:(MM-1)*NN+NNl,nzl,2),T_0_g(2)*Temp_jet(NNL+1:NN,1,2),Ul_jet,ug_jet,&
                    eta_liq*delta_jet,eta_gas*(R_c-delta_jet)+delta_jet,IC_l(nzl,3,2,2),IC_g(1,3,2,2))
            mim(3)=abs(q_l*cp_l*conc_l*(T_0_l(2)-IC_l(nzl,3,2,2)))/flux_test(3)!abs(q_g*cp_g*conc_g*(T_0_g(2)-IC_g(1,3,2,2)))/flux_test(3)!
!            1-(abs(q_l*cp_l*conc_l/rho_l*T_0_l(2))+abs(q_g*cp_g*conc_g/rho_g*T_0_g(2)))/&
!                    (abs(q_l*cp_l*conc_l/rho_l*IC_l(nzl,3,2,2))+abs(q_g*cp_g*conc_g/rho_g*IC_g(1,3,2,2)))




!do i=nzl-1,nzl-1
!do z=1,mm
!!####################################################################################################################################
!!#####################                       FILM                                                            ########################
!!####################################################################################################################################
!
!Test((i-1)*mm+z,1,1)=conc_l*D_l(2)/(delta_film)*(-&
!conc_film(nn*(z-1)+nnl-2,i,2)*((e_liq_f(NNl-1)-e_liq_f(NNl))/((e_liq_f(NNl-2)-e_liq_f(NNl-1))*(e_liq_f(NNl-2)-e_liq_f(NNl))))+&
!conc_film(nn*(z-1)+nnl-1,i,2)*((e_liq_f(NNl-2)-e_liq_f(NNl))/((e_liq_f(NNl-2)-e_liq_f(NNl-1))*(e_liq_f(NNl-1)-e_liq_f(NNl))))+&
!conc_film(nn*(z-1)+nnl,i,2)*(((e_liq_f(NNl-1)-e_liq_f(NNl))**2-(e_liq_f(NNl-2)-e_liq_f(NNl))**2)/((e_liq_f(NNl-2)-e_liq_f(NNl-1))*&
!                                    (e_liq_f(NNl-2)-e_liq_f(NNl))*(e_liq_f(NNl-1)-e_liq_f(NNl)))))
!
!Test((i-1)*mm+z,2,1)=(conc_g)*D_g(2)/(R_c-delta_film)*(-&
!conc_film(nn*(z-1)+nnl+2,i,2)*((e_gas_f(1)-e_gas_f(0))/((e_gas_f(2)-e_gas_f(1))*(e_gas_f(2)-e_gas_f(0))))+&
!conc_film(nn*(z-1)+nnl+1,i,2)*((e_gas_f(2)-e_gas_f(0))/((e_gas_f(2)-e_gas_f(1))*(e_gas_f(1)-e_gas_f(0))))+&
!conc_l/(conc_g)*dist_coef(2)*&
!conc_film(nn*(z-1)+nnl,i,2)*(((e_gas_f(1)-e_gas_f(0))**2-(e_gas_f(2)-e_gas_f(0))**2)/&
!            ((e_gas_f(2)-e_gas_f(1))*(e_gas_f(2)-e_gas_f(0))*(e_gas_f(1)-e_gas_f(0)))))
!
!!####################################################################################################################################
!!#####################                        JET                                                            ########################
!!####################################################################################################################################
!Test((i-1)*mm+z,1,2)=conc_l*D_l(2)/(delta_jet)*(-&
!conc_Jet(nn*(z-1)+nnl-2,i,2)*((e_liq(NNl-1)-e_liq(NNl))/((e_liq(NNl-2)-e_liq(NNl-1))*(e_liq(NNl-2)-e_liq(NNl))))+&
!conc_Jet(nn*(z-1)+nnl-1,i,2)*((e_liq(NNl-2)-e_liq(NNl))/((e_liq(NNl-2)-e_liq(NNl-1))*(e_liq(NNl-1)-e_liq(NNl))))+&
!conc_Jet(nn*(z-1)+nnl,i,2)*(((e_liq(NNl-1)-e_liq(NNl))**2-(e_liq(NNl-2)-e_liq(NNl))**2)/((e_liq(NNl-2)-e_liq(NNl-1))*&
!                                    (e_liq(NNl-2)-e_liq(NNl))*(e_liq(NNl-1)-e_liq(NNl)))))
!
!Test((i-1)*mm+z,2,2)=(conc_g)*D_g(2)/(R_c-delta_jet)*(-&
!conc_Jet(nn*(z-1)+nnl+2,i,2)*((e_gas(1)-e_gas(0))/((e_gas(2)-e_gas(1))*(e_gas(2)-e_gas(0))))+&
!conc_Jet(nn*(z-1)+nnl+1,i,2)*((e_gas(2)-e_gas(0))/((e_gas(2)-e_gas(1))*(e_gas(1)-e_gas(0))))+&
!conc_l/(conc_g)*dist_coef(2)*&
!conc_Jet(nn*(z-1)+nnl,i,2)*(((e_gas(1)-e_gas(0))**2-(e_gas(2)-e_gas(0))**2)/&
!            ((e_gas(2)-e_gas(1))*(e_gas(2)-e_gas(0))*(e_gas(1)-e_gas(0)))))
!!####################################################################################################################################
!!#####################                        DROP                                                           ########################
!!####################################################################################################################################
!Test((i-1)*mm+z,1,3)=conc_l*D_l(2)/(delta_drop)*(-&
!conc_Drop(nn*(z-1)+nnl-2,i,2)*((e_liq(NNl-1)-e_liq(NNl))/((e_liq(NNl-2)-e_liq(NNl-1))*(e_liq(NNl-2)-e_liq(NNl))))+&
!conc_Drop(nn*(z-1)+nnl-1,i,2)*((e_liq(NNl-2)-e_liq(NNl))/((e_liq(NNl-2)-e_liq(NNl-1))*(e_liq(NNl-1)-e_liq(NNl))))+&
!conc_Drop(nn*(z-1)+nnl,i,2)*(((e_liq(NNl-1)-e_liq(NNl))**2-(e_liq(NNl-2)-e_liq(NNl))**2)/((e_liq(NNl-2)-e_liq(NNl-1))*&
!                                    (e_liq(NNl-2)-e_liq(NNl))*(e_liq(NNl-1)-e_liq(NNl)))))
!
!Test((i-1)*mm+z,2,3)=(conc_g)*D_g(2)/(R_c-delta_drop)*(-&
!conc_Drop(nn*(z-1)+nnl+2,i,2)*((e_gas(1)-e_gas(0))/((e_gas(2)-e_gas(1))*(e_gas(2)-e_gas(0))))+&
!conc_Drop(nn*(z-1)+nnl+1,i,2)*((e_gas(2)-e_gas(0))/((e_gas(2)-e_gas(1))*(e_gas(1)-e_gas(0))))+&
!conc_l/(conc_g)*dist_coef(2)*&
!conc_Drop(nn*(z-1)+nnl,i,2)*(((e_gas(1)-e_gas(0))**2-(e_gas(2)-e_gas(0))**2)/&
!            ((e_gas(2)-e_gas(1))*(e_gas(2)-e_gas(0))*(e_gas(1)-e_gas(0)))))
!end do
!end do
!flux_test(1)=a_dg_f*abs(sum(test(:,2,1)))!-a_dg_f*abs(sum(test(:,2,1))))/(a_dg_f*abs(sum(test(:,1,1))))
!flux_test(2)=a_dg_j*abs(sum(test(:,2,2)))!-a_dg_j*abs(sum(test(:,2,2))))/(a_dg_j*abs(sum(test(:,1,2))))
!flux_test(3)=a_dg_d*abs(sum(test(:,2,3)))!-a_dg_d*abs(sum(test(:,2,3))))/(a_dg_d*abs(sum(test(:,1,3))))
!            j=j
!
!            call mean_int(conc_l*conc_film((MM-1)*NN+1:(MM-1)*NN+NNl,nzl-2,2),conc_g*conc_film(NNL+1:NN,nzl,2),ul_film,Ug_film,&
!                    (1-eta_liq)*delta_film+R_c-delta_film,(1-eta_gas)*(R_c-delta_film),IC_l(nzl,1,2,2),IC_g(1,1,2,2))
!            call mean_int(conc_l*conc_film(1:NNl,nzl,2),conc_g*conc_film((mm-1)*nn+NNL+1:mm*NN,nzl-2,2),ul_film,Ug_film,&
!                    (1-eta_liq)*delta_film+R_c-delta_film,(1-eta_gas)*(R_c-delta_film),IC_l(1,1,2,2),IC_g(nzl,1,2,2))
!
!
!            mim(1)=abs(q_l*(IC_l(1,1,2,2)-IC_l(nzl,1,2,2)))/flux_test(1)!abs(q_g*(IC_g(nzl,1,2,2)-IC_g(1,1,2,2)))/flux_test(1)!
!!            1-(abs(q_l*cp_l*conc_l/rho_l*conc_l)+abs(q_g*cp_g*conc_g/rho_g*conc_g))/&
!!                    (abs(q_l*cp_l*conc_l/rho_l*IC_l(nzl,1,2,2))+abs(q_g*cp_g*conc_g/rho_g*IC_g(1,1,2,2)))
!
!
!            call mean_int(conc_l*conc_drop((MM-1)*NN+1:(MM-1)*NN+NNl,nzl-2,2),conc_g*conc_drop(NNL+1:NN,nzl,2),ul_drop,ug_drop,&
!                    eta_liq*delta_drop,eta_gas*(R_c-delta_drop)+delta_drop,IC_l(nzl,2,2,2),IC_g(1,2,2,2))
!            call mean_int(conc_l*conc_drop(1:NNl,nzl,2),conc_g*conc_drop((mm-1)*nn+NNL+1:mm*NN,nzl-2,2),ul_drop,ug_drop,&
!                    eta_liq*delta_drop,eta_gas*(R_c-delta_drop)+delta_drop,IC_l(1,2,2,2),IC_g(nzl,2,2,2))
!            mim(2)=abs(q_l*(IC_l(1,2,2,2)-IC_l(nzl,2,2,2)))/flux_test(2)!abs(q_g*(IC_g(nzl,2,2,2)-IC_g(1,2,2,2)))/flux_test(2)!
!!            1-(abs(q_l*cp_l*conc_l/rho_l*conc_l)+abs(q_g*cp_g*conc_g/rho_g*conc_g))/&
!!                    (abs(q_l*cp_l*conc_l/rho_l*IC_l(nzl,2,2,2))+abs(q_g*cp_g*conc_g/rho_g*IC_g(1,2,2,2)))
!
!
!            call mean_int(conc_l*conc_jet((MM-1)*NN+1:(MM-1)*NN+NNl,nzl-2,2),conc_g*conc_jet(NNL+1:NN,nzl,2),Ul_jet,ug_jet,&
!                    eta_liq*delta_jet,eta_gas*(R_c-delta_jet)+delta_jet,IC_l(nzl,3,2,2),IC_g(1,3,2,2))
!            call mean_int(conc_l*conc_jet(1:NNl,nzl,2),conc_g*conc_jet((mm-1)*nn+NNL+1:mm*NN,nzl-2,2),Ul_jet,ug_jet,&
!                    eta_liq*delta_jet,eta_gas*(R_c-delta_jet)+delta_jet,IC_l(1,3,2,2),IC_g(nzl,3,2,2))
!            mim(3)=abs(q_l*(IC_l(1,3,2,2)-IC_l(nzl,3,2,2)))/flux_test(3)!abs(q_g*(IC_g(nzl,3,2,2)-IC_g(1,3,2,2)))/flux_test(3)!
!!            1-(abs(q_l*cp_l*conc_l/rho_l*T_0_l(2))+abs(q_g*cp_g*conc_g/rho_g*conc_g))/&
!!                    (abs(q_l*cp_l*conc_l/rho_l*IC_l(nzl,3,2,2))+abs(q_g*cp_g*conc_g/rho_g*IC_g(1,3,2,2)))
!









            if(mod(j,10).eq.0) then

                    call writing(tcalc(1),eta_liq,delta_film,delta_drop,delta_jet,eta_gas,R_c, conc_film, conc_drop, conc_jet,&
                                conc_l,conc_g, x_0,y_0 ,ul_film,ug_film,ul_drop,ug_drop,ul_jet,ug_jet,&
                                IC_l(:,:,:,1),IC_g(:,:,:,1),H,verh,rms_c,Mim)

                    call writing(tcalc(2),eta_liq,delta_film,delta_drop,delta_jet,eta_gas,R_c,temp_film, temp_drop, temp_jet,&
                                conc_l,conc_g, T_0_l,T_0_g, ul_film,ug_film,ul_drop,ug_drop,ul_jet,ug_jet,&
                                IC_l(:,:,:,2),IC_g(:,:,:,2),H,verh,rms_c,Mim)


            end if

            j=j+1
        end do
!    endif

end do

call cpu_time(finish)

write(*,*) 'iterations took' , finish-start , 'seconds'
write(666,*) 'iterations took' , finish-start , 'seconds'
close(666)

write(*,*) 'Start writing output. Close all open files in output directory and press Enter to continue.'
!read(*,*)

    call writing(tcalc(1),eta_liq,delta_film,delta_drop,delta_jet,eta_gas,R_c, conc_film, conc_drop, conc_jet,&
                conc_l,conc_g, x_0,y_0 ,ul_film,ug_film,ul_drop,ug_drop,ul_jet,ug_jet,IC_l(:,:,:,1),IC_g(:,:,:,1),H,verh,rms_c,Mim)

    call writing(tcalc(2),eta_liq,delta_film,delta_drop,delta_jet,eta_gas,R_c,temp_film, temp_drop, temp_jet,&
                conc_l,conc_g, T_0_l,T_0_g, ul_film,ug_film,ul_drop,ug_drop,ul_jet,ug_jet,IC_l(:,:,:,2),IC_g(:,:,:,2),&
                H,verh,rms_c,Mim)


end subroutine control

end program



