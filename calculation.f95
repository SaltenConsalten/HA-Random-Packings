module num_solution_mod

USE omp_lib
USE global_variables
USE phys_prop
USE Algorithms
USE tools
USE hydrodynamics

implicit none

contains


subroutine trubulence_parameter(R_c, zl_film, rho_l,rho_g, mu_l,mu_g, dx,dh, q_g,re_l, re_g, &
                                eta_gas,turb_ga_film,turb_ga_drop,turb_ga_jet,SP_g,&
                                int_ul_film,int_ul_jet,int_ul_drop,int_ug_film,int_ug_jet,int_ug_drop)

    implicit none

    real*8, intent(in) :: rho_l,rho_g, mu_l,mu_g           ! Density and dyn. viscosity
    real*8, intent(in) :: q_g
    real*8, intent(in) :: R_c, zl_film
    real*8 :: eta_gas(nng)
    real*8, intent(inout) :: Re_l,Re_g                                    ! Reynolds number
    real*8, intent(in):: SP_g(ncomp)                                  ! Schmidt (concentration calc) or Prandtl (temp calc) number
!    real*8 :: dx(3,2), dh(3,2)
    real*8 :: dx(2,3), dh(2,3)
    real*8, intent(in) :: int_ul_film,int_ul_jet,int_ul_drop,int_ug_film,int_ug_jet,int_ug_drop
    real*8, intent(inout) :: turb_ga_film(NNg+1,ncomp),turb_ga_drop(NNg+1,ncomp),turb_ga_jet(NNg+1,ncomp)
    integer :: k,i




turb_ga_film(:,:)=0D0
turb_ga_drop(:,:)=0D0
turb_ga_jet(:,:)=0D0
!SP_g(:)=1D0
dh(:,:)=H/dfloat(nzl*mm)
dx(:,:)=-dh


!######################################################################################################
          !film
!######################################################################################################

    Re_l=4*delta_film*int_ul_film*rho_l/mu_l
    Re_g=2*(R_c-delta_film)*int_ug_film*rho_g/mu_g

    dx(1,1)=4*dx(1,1)/(Re_l*delta_film)
    dx(2,1)=2*dx(2,1)/(Re_g*(R_c-delta_film))


    if (turbulent) then

        eta_gas=1-eta_gas
        u_tau = sqrt(1.D0/(8.D0*200.D0**0.25D0)*((q_g/(pi*(R_c-delta_film)**2))**1.75D0)*(mu_g/rho_g)**0.25D0*&
                (R_c-delta_film)**(-0.25D0))
        C_mu_turb = B_*rho_g*(R_c-delta_film)*u_tau*Kappa/(2*A_r)

    do k=1,ncomp
        !Turbulent Gamma calculation
        do i=1,(NNg+1)
            if (i.eq.1) then
                turb_ga_film(i,k)=(SP_g(k)**(-1)+C_mu_turb/mu_g*(1-((eta_gas(1)+1)/2)**2)*&
                                    (1+2*((eta_gas(1)+1)/2)**2))
            elseif (i.eq.NNg+1) then
                turb_ga_film(i,k)=(SP_g(k)**(-1)+C_mu_turb/mu_g*(1-((eta_gas(NNg))/2)**2)*&
                                    (1+2*((eta_gas(NNg))/2)**2))
            else
                turb_ga_film(i,k)=(SP_g(k)**(-1)+C_mu_turb/mu_g*(1-((eta_gas(i-1)+eta_gas(i))/2)**2)*&
                                    (1+2*((eta_gas(i-1)+eta_gas(i))/2)**2))
            end if

        end do
    end do

        eta_gas=1-eta_gas

    else
        forall (k=1:ncomp) turb_ga_film(:,k)=1D0/SP_g(k)
    endif

    !######################################################################################################
              !Drop
    !######################################################################################################

    Re_l=4*delta_drop*int_ul_drop*rho_l/mu_l
    Re_g=2*(R_c-delta_drop)*int_ug_drop*rho_g/mu_g

    dx(1,2)=4*dx(1,2)/(Re_l*delta_drop)
    dx(2,2)=2*dx(2,2)/(Re_g*(R_c-delta_drop))

    if (turbulent) then

        u_tau = sqrt(1.D0/(8.D0*200.D0**0.25D0)*((q_g/(pi*(R_c-delta_drop)**2))**1.75D0)*(mu_g/rho_g)**0.25D0*&
                (R_c-delta_drop)**(-0.25D0))
        C_mu_turb = B_*rho_g*(R_c-delta_drop)*u_tau*Kappa/(2*A_r)
    do k=1,ncomp
        !Turbulent Gamma calculation

        do i=1,(NNg+1)
            if (i.eq.1) then
                turb_ga_drop(i,k)=(SP_g(k)**(-1)+C_mu_turb/mu_g*(1-((1+eta_gas(1))/2)**2)*&
                                    (1+2*((1+eta_gas(1))/2)**2))
            elseif (i.eq.NNg+1) then
                turb_ga_drop(i,k)=(SP_g(k)**(-1)+C_mu_turb/mu_g*(1-(eta_gas(NNg)&
                                    /2)**2)*(1+2*((eta_gas(NNg))/2)**2))
            else
                turb_ga_drop(i,k)=(SP_g(k)**(-1)+C_mu_turb/mu_g*(1-&
                                    ((eta_gas(i-1)+eta_gas(i))/2)**2)*(1+2*((eta_gas(i-1)+eta_gas(i))/2)**2))
            end if

        end do
    enddo

    else

        forall (k=1:ncomp) turb_ga_drop(:,k)=1D0/SP_g(k)

    end if

    !######################################################################################################
              !Jet
    !######################################################################################################

    Re_l=4*delta_jet*int_ul_jet*rho_l/mu_l
    Re_g=2*(R_c-delta_jet)*int_ug_jet*rho_g/mu_g

    dx(1,3)=4*dx(1,3)/(Re_l*delta_jet)
    dx(2,3)=2*dx(2,3)/(Re_g*(R_c-delta_jet))

    if (turbulent) then

        u_tau = sqrt(1.D0/(8.D0*200.D0**0.25D0)*((q_g/(pi*(R_c-delta_jet)**2))**1.75D0)*(mu_g/rho_g)**0.25D0*&
                (R_c-delta_jet)**(-0.25D0))
        C_mu_turb = B_*rho_g*(R_c-delta_jet)*u_tau*Kappa/(2*A_r)
    do k=1,ncomp
        !Turbulent Gamma calculation

        do i=1,(NNg+1)
            if (i.eq.1) then
                turb_ga_jet(i,k)=(SP_g(k)**(-1)+C_mu_turb/mu_g*(1-((1+eta_gas(1))/2)**2)*&
                                    (1+2*((1+eta_gas(1))/2)**2))
            elseif (i.eq.NNg+1) then
                turb_ga_jet(i,k)=(SP_g(k)**(-1)+C_mu_turb/mu_g*(1-(eta_gas(NNg)&
                                    /2)**2)*(1+2*((eta_gas(NNg))/2)**2))
            else
                turb_ga_jet(i,k)=(SP_g(k)**(-1)+C_mu_turb/mu_g*(1-&
                                    ((eta_gas(i-1)+eta_gas(i))/2)**2)*(1+2*((eta_gas(i-1)+eta_gas(i))/2)**2))
            end if

        end do
    enddo

    else

        forall (k=1:ncomp) turb_ga_jet(:,k)=1D0/SP_g(k)

    end if

end subroutine



subroutine calculation(j,eta_gas,eta_liq, conc_film,conc_drop,conc_jet,R_c,conc_old_film,conc_old_drop,&
                    conc_old_jet,Sc_l,Dist_coef,conc_l,conc_g, D_l,D_g,delta_film,delta_drop,delta_jet,&
                    Ug_film,Ug_drop,Ug_jet, Ul_film,Ul_drop,Ul_jet,turb_ga_film,turb_ga_drop,turb_ga_jet,&
                    nzl_all,Nzg, dx,dh, IC_l,IC_g, x_0,y_0, verh, k_reac,val_vec,meanmat,row,col_ind,M,L,&
                    n_A,h_r,tcalc,Mnew)

implicit none


  real*8, intent(in) :: eta_gas(nng),eta_liq(nnl)
  real*8 :: conc_film(NN*MM,Nzl),conc_drop(NN*MM,Nzl),conc_jet(NN*MM,Nzl)!, conc(NN*MM,Nzl)
  real*8 :: conc_old_film(NN*MM,Nzl),conc_old_drop(NN*MM,Nzl),conc_old_jet(NN*MM,Nzl)!, Temp_old(:,:,:),conc_old(:,:,:,:)
  real*8, intent(in) :: conc_l, conc_g , x_0,y_0                      ! [kmol/m³] molar densities
  real*8, intent(in) :: D_l, D_g            ! diffusion coefficient
  real*8, intent(in) :: R_c, Sc_l
  real*8, intent(in) :: Dist_coef           ! Henry coefficient of components
  real*8, intent(in) :: delta_film,delta_drop,delta_jet
  real*8, intent(in) :: Ug_film(nng),Ug_jet(nng),Ug_drop(nng),Ul_film(nnl),Ul_jet(nnl),Ul_drop(nnl)
  real*8, intent(in) :: turb_ga_film(nng+1), turb_ga_jet(nng+1), turb_ga_drop(nng+1)
  integer, intent(in) :: nzl_all, Nzg                   !number of liquid mixing points
  logical :: recalc                         ! declares if the cholesky decomposition has to be recalculated
!  real*8, intent(in) :: dx(3,2) , dh(3,2)                            !axial increment
  real*8, intent(in) :: dx(2,3) , dh(2,3)                            !axial increment
  integer :: i,j,z                              !control variables
  real*8, intent(in) :: n_A(nzl*mm,3), h_r(nn*mm,nzl,3)
  real*8  :: IC_l(nzl*mm,3),IC_g(nzl*mm,3)                        ! integral mean concentrations
!  logical :: calc=.false.
  logical, intent(in) :: tcalc
  logical :: Mnew
!  real*8 :: flux_l, flux_g
  real*8 :: val_vec(5*NN*MM-2*(NN+(MM-1)),3), meanmat(3,3), M(:,:,:), L(:,:,:)
  integer :: row(5*NN*MM-2*(NN+(MM-1)),3), col_ind(NN*MM+1,3)
  real*8, intent(in) :: verh(3)
  real*8, intent(in) :: k_reac(3)
!  real*8 :: start=0, finish=0


        !$omp parallel num_threads(3)
         !$omp sections
          !$omp section

!######################################################################################################
!          !film
!######################################################################################################

        call solve(Nzl,Nzg,j,R_c,delta_film,Sc_l,turb_ga_film(:),D_g,D_l,Dist_coef,k_reac(:),dx(:,1),dh(:,1),&
                1-eta_gas(:),1-eta_liq(:),ul_film(:),ug_film(:),x_0,conc_l,y_0,conc_g,conc_old_film(:,:),&
                conc_film(:,:),val_vec(:,1),meanmat(:,1),row(:,1),col_ind(:,1),M(:,:,1),L(:,:,1),&
                n_A(:,2),h_r(:,:,2),tcalc,Mnew)

          !$omp section
!######################################################################################################
!          !drop
!######################################################################################################

        call solve(Nzl,Nzg,j,R_c,delta_drop,Sc_l,turb_ga_drop(:),D_g,D_l,Dist_coef,k_reac(:),dx(:,2),dh(:,1),&
                eta_gas(:),eta_liq(:),ul_drop(:),ug_drop(:),x_0,conc_l,y_0,conc_g,conc_old_drop(:,:),&
                conc_drop(:,:),val_vec(:,2),meanmat(:,2),row(:,2),col_ind(:,2),M(:,:,2),L(:,:,2),&
                n_A(:,2),h_r(:,:,2),tcalc,Mnew)

          !$omp section
!!!######################################################################################################
!!!          !jet
!!!######################################################################################################

        call solve(Nzl,Nzg,j,R_c,delta_jet,Sc_l,turb_ga_jet(:),D_g,D_l,Dist_coef,k_reac(:),dx(:,3),dh(:,1),&
                eta_gas(:),eta_liq(:),ul_jet(:),ug_jet(:),x_0,conc_l,y_0,conc_g,conc_old_jet(:,:),&
                conc_jet(:,:),val_vec(:,3),meanmat(:,3),row(:,3),col_ind(:,3),M(:,:,3),L(:,:,3),&
                n_A(:,2),h_r(:,:,2),tcalc,Mnew)

          !$omp end sections
         !$omp barrier
        !$omp end parallel

if (Mnew) then
    Mnew=.false.
end if

!!flux_g=y_0(k)
!do z=Nzl-1,2,-1!2,Nzl-1!
!!flux_g=1
!    call mean_int(Conc_film((MM-1)*NN+1:(MM-1)*NN+NNl,z-1),conc_film(NNL+1:NN,z+1),ul_film,Ug_film,&        !film mixing
!                    (1-eta_liq)*delta_film+R_c-delta_film,(1-eta_gas)*(R_c-delta_film),IC_l(z,1),IC_g(z,1))
!
!    Conc_film((MM-1)*NN+1:(MM-1)*NN+NNl,z-1)=IC_l(z,1)
!
!    if ((mod(z,nint(dble(nzl)/dble(nzl_all))).eq.0)) then                                                          !overall liquid mixing
!
!        call mean_int(Conc_drop((MM-1)*NN+1:(MM-1)*NN+NNl,z-1),Conc_drop(NNL+1:NN,z+1),ul_drop,ug_drop,&
!                    eta_liq*delta_drop,eta_gas*(R_c-delta_drop)+delta_drop,IC_l(z,2),IC_g(z,2))
!
!        call mean_int(Conc_jet((MM-1)*NN+1:(MM-1)*NN+NNl,z-1),Conc_jet(NNL+1:NN,z+1),Ul_jet,ug_jet,&
!                    eta_liq*delta_jet,eta_gas*(R_c-delta_jet)+delta_jet,IC_l(z,3),IC_g(z,3))
!
!        Conc_film((MM-1)*NN+1:(MM-1)*NN+NNl,z-1)=verh(1)*IC_l(z,1)+verh(2)*IC_l(z,2)+verh(3)*IC_l(z,3)
!        Conc_drop((MM-1)*NN+1:(MM-1)*NN+NNl,z-1)=verh(1)*IC_l(z,1)+verh(2)*IC_l(z,2)+verh(3)*IC_l(z,3)
!        Conc_jet((MM-1)*NN+1:(MM-1)*NN+NNl,z-1)=verh(1)*IC_l(z,1)+verh(2)*IC_l(z,2)+verh(3)*IC_l(z,3)
!
!
!
!
!    endif
!
!   if (mod(z,int(dble(nzl)/dble(nzg))).eq.0)then                                                           !gas mixing
!
!        call mean_int(Conc_drop((MM-1)*NN+1:(MM-1)*NN+NNl,z-1),Conc_drop(NNL+1:NN,z+1),ul_drop,ug_drop,&
!                    eta_liq*delta_drop,eta_gas*(R_c-delta_drop)+delta_drop,IC_l(z,2),IC_g(z,2))
!
!        call mean_int(Conc_jet((MM-1)*NN+1:(MM-1)*NN+NNl,z-1),Conc_jet(NNL+1:NN,z+1),Ul_jet,ug_jet,&
!                    eta_liq*delta_jet,eta_gas*(R_c-delta_jet)+delta_jet,IC_l(z,3),IC_g(z,3))
!
!            Conc_film(NNL+1:NN,z+1)=ae*(verh(1)*IC_g(z,1)+verh(2)*IC_g(z,2)+verh(3)*IC_g(z,3))!+(1-ae)*flux_g
!            Conc_drop(NNL+1:NN,z+1)=ae*(verh(1)*IC_g(z,1)+verh(2)*IC_g(z,2)+verh(3)*IC_g(z,3))!+(1-ae)*flux_g
!            Conc_jet(NNL+1:NN,z+1)=ae*(verh(1)*IC_g(z,1)+verh(2)*IC_g(z,2)+verh(3)*IC_g(z,3))!+(1-ae)*flux_g
!
!
!           ! flux_g=ae*(verh(1)*IC_g(k,z,1)+verh(2)*IC_g(k,z,2)+verh(3)*IC_g(k,z,3))+(1-ae)*flux_g
!
!
!    endif
!
!    if (z.eq.2) then
!
!        call mean_int(Conc_film((MM-1)*NN+1:(MM-1)*NN+NNl,z-1),conc_film(NNL+1:NN,1),ul_film,Ug_film,&        ! mixing at the end of column (ae)
!                    (1-eta_liq)*delta_film+R_c-delta_film,(1-eta_gas)*(R_c-delta_film),IC_l(z,1),IC_g(z,1))
!
!        call mean_int(Conc_drop((MM-1)*NN+1:(MM-1)*NN+NNl,z-1),Conc_drop(NNL+1:NN,1),ul_drop,ug_drop,&
!                    eta_liq*delta_drop,eta_gas*(R_c-delta_drop)+delta_drop,IC_l(z,2),IC_g(z,2))
!
!        call mean_int(Conc_jet((MM-1)*NN+1:(MM-1)*NN+NNl,z-1),Conc_jet(NNL+1:NN,1),Ul_jet,ug_jet,&
!                    eta_liq*delta_jet,eta_gas*(R_c-delta_jet)+delta_jet,IC_l(z,3),IC_g(z,3))
!
!            Conc_film(NNL+1:NN,1)=ae*(verh(1)*IC_g(z,1)+verh(2)*IC_g(z,2)+verh(3)*IC_g(z,3))!+(1-ae)*flux_g
!            Conc_drop(NNL+1:NN,1)=ae*(verh(1)*IC_g(z,1)+verh(2)*IC_g(z,2)+verh(3)*IC_g(z,3))!+(1-ae)*flux_g
!            Conc_jet(NNL+1:NN,1)=ae*(verh(1)*IC_g(z,1)+verh(2)*IC_g(z,2)+verh(3)*IC_g(z,3))!+(1-ae)*flux_g
!
!    end if
!
!end do


end subroutine



subroutine rms_err(eta_gas, eta_liq, trans_film, trans_drop, trans_jet, trans,conc_l, conc_g, IC_l, IC_g, rms,&
                trans_old_film, trans_old_drop, trans_old_jet, ug_drop, ul_drop,&!, ug_film, ug_jet, ul_film, ul_jet,&
                k, fvoll, fvolg, tcalc, transl_0, transg_0, Mim)

    implicit none

  real*8, intent(in) :: eta_gas(nng),eta_liq(nnl)
  real*8, intent(in) :: trans_film(NN*MM,Nzl,ncomp),trans_drop(NN*MM,Nzl,ncomp),trans_jet(NN*MM,Nzl,ncomp), trans(NN*MM,Nzl,ncomp)
  real*8, intent(in) :: trans_old_film(NN*MM,Nzl,ncomp),trans_old_drop(NN*MM,Nzl,ncomp),trans_old_jet(NN*MM,Nzl,ncomp)!, trans_old(:,:,:,:)
  real*8, intent(in) :: conc_l, conc_g                      ! [kmol/m³] molar densities
!  real*8, intent(in) :: R_c, H
  real*8, intent(in) :: ug_drop(nng),ul_drop(nnl)!,ul_film(nnl),ul_jet(nnl),ug_film(nng),ug_jet(nng)
!  integer :: nzl_all, Nzg                   !number of liquid mixing points
  integer :: k                              !control variables
  real*8 :: IC_l(nzl*MM,3,ncomp),IC_g(nzl*MM,3,ncomp)                        ! integral mean concentrations
!  real*8 :: transl_0(:),transg_0(:)                    ! input concentrations or temprature
  real*8, intent(in) :: fvoll,fvolg
  real*8 :: rms(ncomp,3),Mim(ncomp)                        ! root mean square errors of concentration field/ temperatur field
  real*8, intent(in) :: transl_0(ncomp),transg_0(ncomp)                ! molar imbalances
!  real*8 :: rms_t(2)=0                              ! root mean square error of temperature field
!  real*8, intent(in) :: k_reac(ncomp,2)
!  real*8 :: test(3)
  logical, intent(in) :: tcalc
!  integer :: thread
!  thread=omp_get_thread_num()


!###############################################################################
!   calculate RMS and imbalances
!###############################################################################

        RMS(k,1)=dsqrt(sum((trans_old_film(:,:,k)-trans_film(:,:,k))**2)/size(trans_film(:,:,k)))/(sum(trans_film(:,:,k))&
                            /size(trans_film(:,:,k)))

        RMS(k,2)=dsqrt(sum((trans_old_drop(:,:,k)-trans_drop(:,:,k))**2)/size(trans_drop(:,:,k)))/(sum(trans_drop(:,:,k))&
                            /size(trans_drop(:,:,k)))

        RMS(k,3)=dsqrt(sum((trans_old_jet(:,:,k)-trans_jet(:,:,k))**2)/size(trans_drop(:,:,k)))/(sum(trans_jet(:,:,k))&
                            /size(trans_drop(:,:,k)))

!        test(1)=sum(trans_film(k,:,:))
!        test(2)=sum(trans_drop(k,:,:))
!        test(3)=sum(trans_jet(k,:,:))
!        if ((isnan(rms(k,1))).or.(rms(k,1).eq.0)) then
!                RMS(k,1)=1
!            endif
!        if ((isnan(rms(k,2))).or.(rms(k,2).eq.0)) then
!                RMS(k,2)=1
!            endif
!        if ((isnan(rms(k,3))).or.(rms(k,3).eq.0)) then
!                RMS(k,3)=1
!            endif

        call mean_int(trans((MM-1)*NN+1:(MM-1)*NN+NNl,Nzl,k),trans(NNl+1:NN,1,k),ul_drop,Ug_drop,&
                    eta_liq*delta_drop+R_c-delta_drop,eta_gas*(R_c-delta_drop),IC_l(Nzl,1,k),IC_g(1,1,k))

        if (.not.tcalc) then
            Mim(k)=(IC_l(Nzl,1,k)*FVOLL+IC_g(1,1,k)*FVOLG)
            Mim(k)=1-(IC_l(Nzl,1,k)*transl_0(k)*conc_l*FVOLL+IC_g(1,1,k)*transg_0(k)*conc_g*FVOLG)&
                    /(transl_0(k)*conc_l*fvoll+transg_0(k)*conc_g*FVOLG)
        end if

    end subroutine


subroutine kreaction(temp,k_reac)

    implicit none
    real*8, intent(in) :: temp
    real*8, intent(out) :: k_reac
    real*8 ::  k, k_0
    real*8 :: EA, COH, beta

    EA=54971D0
    COH=1.177
    k_0=3.27869D13

    beta=2.83D-4*temp**2-1.7367D-1*temp+26.809D0
    k=k_0*exp(-EA/(8.3145*temp))
    k=k*exp(beta*1D0)

    k_reac=k*COH

    end subroutine

end module
