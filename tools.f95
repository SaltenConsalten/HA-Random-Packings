module tools


USE global_variables
    contains
!#############################################################
!this module contains some subroutines needed in the control program
!#############################################################

    subroutine integration(f,x,intf,t)

!###############################################################################
! This subroutine contains the trapezoidal rule for numerical integration of a
! function f over the variable x.
!###############################################################################

        implicit none

        real*8, intent(in) :: f(:), x(:)    !f = function to be integrated, x = integration variable
        integer :: i,n
        real*8, intent(out) :: intf         !integration value
        integer, optional :: t

    if (present(t)) then

        write (*,*) size(x), size(f)
    end if

        n=size(x)

        intf=0D0 !x(1)*f(1)          !zero gradient at x=0 integration formula results in x(1)*f(1)
        do i=1,n-1
            intf=intf+(x(i+1)-x(i))*(f(i)+f(i+1))/2D0
        end do

        if (x(1).gt.x(n))   intf=-intf


    end subroutine integration

    subroutine restore_balance(calc,c_l_in, c_l_out, c_g_in, c_g_out, ul, ug, rl, rg,int_l_in,int_g_in, dissflux_l,&
                                dissflux_g)
!#############################################################################
! due to conflicts in the inlet boundary and inter phase boundary condition
! there is a lag of mass/heat continuity. This subroutine adds the dissipated
! mass/heat to the inlet boundary.
!#############################################################################
        implicit none

        real*8, intent (in) :: c_l_in(NNl), c_l_out(NNl), c_g_in(NNg), c_g_out(NNg)      !in and out coming concentration profiles for nth and n+1st stage (liquid and gas)
        real*8, intent(in) :: ul(nnl), ug(nng)                                !liquid and gas velocity profiles
        real*8, intent(in) :: rl(nnl), rg(nng)                                !radial coordinates for mentioned profiles
       ! real*8, intent(in) :: rms
        real*8 :: int_l_in,int_l_out,int_g_in,int_g_out             !mean conc/temp values
        real*8 :: nom_l_in,nom_l_out,nom_g_in,nom_g_out             !nominator of conc/temp values (Shilkin 6.14)
        real*8 :: denom_l,denom_g, test                                   !denominator of conc/temp values (Shilkin 6.14)
        real*8, intent(inout) :: dissflux_l, dissflux_g               ! dissipated mass/heat
        logical, intent(in) :: calc!=.false.

        call integration(c_l_in*ul*rl,rl,nom_l_in)
        call integration(c_l_out*ul*rl,rl,nom_l_out)
        call integration(c_g_in*ug*rg,rg,nom_g_in)
        call integration(c_g_out*ug*rg,rg,nom_g_out)
        call integration(ul*rl,rl,denom_l)
        call integration(ug*rg,rg,denom_g)

        int_l_in=nom_l_in/denom_l
        int_l_out=nom_l_out/denom_l
        int_g_in=nom_g_in/denom_g
        int_g_out=nom_g_out/denom_g



        if(calc) then
            dissflux_l=dissflux_l+1D-1*(int_l_in-int_l_out)!1D-1*
            dissflux_g=dissflux_g+1D-2*(int_g_in-int_g_out)!1D-2*

        else

            dissflux_l=0!dissflux_l/abs(dissflux_l)*5E-2*int_l_in
            dissflux_g=0
        end if

        test=int_l_in-int_l_out
        test=test

    end subroutine

    subroutine mean_int(c_l, c_g, ul, ug, rl, rg,int_l,int_g)

        implicit none

        real*8, intent (in) :: c_l(NNl), c_g(NNg)   !in and out coming concentration profiles for nth and n+1st stage (liquid and gas)
        real*8, intent(in) :: ul(NNl), ug(NNg)                                !liquid and gas velocity profiles
        real*8, intent(in) :: rl(NNl), rg(NNg)                                !radial coordinates for mentioned profiles
        real*8 :: int_l,int_g             !mean conc/temp values
        real*8 :: nom_l,nom_g           !nominator of conc/temp values (Shilkin 6.14)
        real*8 :: denom_l,denom_g                                   !denominator of conc/temp values (Shilkin 6.14)

        call integration(c_l*ul*rl,rl,nom_l)
        call integration(c_g*ug*rg,rg,nom_g)
        call integration(ul*rl,rl,denom_l)
        call integration(ug*rg,rg,denom_g)

        int_l=abs(nom_l/denom_l)
        int_g=abs(nom_g/denom_g)

    end subroutine


 subroutine brent(f,a,b,t,maxiter)

        implicit none

        real*8, intent(in) :: t
        integer, intent(in) :: maxiter
        real*8, intent(inout) :: a, b
        integer :: iter
        real*8 :: fa, fb, fc, c, d, e, tol, m, p, q, s, r


        interface brent_algorithm
        real*8 function f(x)
            real*8 :: x
        endfunction
        endinterface brent_algorithm

        fa=f(a)
        fb=f(b)

        if (fa*fb.gt.0) then
            write(*,*) 'f(a) und f(b) sollten unterschiedliche Vorzeichen haben'
        endif

        c=a
        fc=fa   !Zu Beginn ist c = a

        c=a
        fc=fa
        d=b-a
        e=d

        iter=0

        do while (iter.lt.maxiter)
            iter=iter+1

            if (fb*fc.gt.0) then
                c=a
                fc=fa
                d=b-a
                e=d
            endif

            if (abs(fc).lt.abs(fb)) then
                a=b
                b=c
                c=a
                fa=fb
                fb=fc
                fc=fa
            endif

            tol=2*epsilon(a)*abs(b)+t
            m=(c-b)/2                !Toleranz

            if ((abs(m)>tol).and.(abs(fb)>0)) then !Verfahren muss noch durchgeführt werden

                if ((abs(e)<tol).or.(abs(fa)<=abs(fb))) then
                    d=m
                    e=m
                else
                    s=fb/fa
                    if (a.eq.c) then
                        p=2*m*s
                        q=1-s
                    else
                        q=fa/fc
                        r=fb/fc
                        p=s*(2*m*q*(q-r)-(b-a)*(r-1))
                        q=(q-1)*(r-1)*(s-1)
                    endif
                    if (p.gt.0) then
                        q=-q
                    else
                        p=-p
                    endif
                    s=e
                    e=d
                    if ((2*p.lt.3*m*q-abs(tol*q)).and.(p.lt.abs(s*q/2))) then
                        d=p/q
                    else
                        d=m
                        e=m
                    endif
                endif

                a=b
                fa=fb

                if (abs(d).gt.tol) then
                    b=b+d
                else
                    if (m.gt.0) then
                        b=b+tol
                    else
                    b=b-tol
                endif
            endif
            else
                return
            endif
            fb=f(b)
        enddo
    end subroutine brent

subroutine writing(tcalc,eta_liq,delta_film,delta_drop,delta_jet,eta_gas,R_c,trans_film,trans_drop,trans_jet&
                ,conc_l,conc_g,transl_0,transg_0,ul_film,ug_film,ul_drop,ug_drop,ul_jet,ug_jet,IC_l,IC_g,H&
                ,verh, rms_c, Mim)
implicit none

    real*8, intent(in) :: Ug_film(NNg),Ug_jet(NNg),Ug_drop(NNg),Ul_film(NNl),Ul_jet(NNl),Ul_drop(NNl)
    real*8, intent(in) :: conc_l, conc_g
    real*8, intent(in) :: trans_film(NN*MM,Nzl,ncomp),trans_drop(NN*MM,Nzl,ncomp),trans_jet(NN*MM,Nzl,ncomp)
    real*8, intent(in) :: IC_l(nzl*MM,3,ncomp),IC_g(nzl*MM,3,ncomp)
    real*8, intent(in) :: transl_0(ncomp),transg_0(ncomp)
    real*8, intent(in) :: eta_gas(NNg),eta_liq(NNl)
    real*8, intent(in) :: delta_film, delta_drop, delta_jet
    real*8, intent(in) :: H, R_c
    real*8, intent(in) :: verh(3)
    real*8, intent(in) :: rms_c(ncomp,3),Mim(ncomp)
    integer :: k,i,j
    logical, intent(in) :: tcalc

!##############################################################################################
!   if note is true, then write down the concentration results, else the temperatur results
!##############################################################################################

if(.not.tcalc)then

    open(16, file='output/film/radial-profiles.csv', status='replace')

    open(26, file='output/drop/radial-profiles.csv', status='replace')

    open(36, file='output/jet/radial-profiles.csv', status='replace')

    do k=1,ncomp
    write(16,*) comp_order(k)
    write(16,"(2010(e13.7,:,','))") (1-eta_liq(:))*delta_film+R_c-delta_film,(1-eta_gas(:))*(R_c-delta_film)

    write(26,*) comp_order(k)
    write(26,"(2010(e13.7,:,','))") eta_liq(:)*delta_drop,eta_gas(:)*(R_c-delta_drop)+delta_drop

    write(36,*) comp_order(k)
    write(36,"(2010(e13.7,:,','))") eta_liq(:)*delta_jet,eta_gas(:)*(R_c-delta_jet)+delta_jet

    do i=1,nzl
            do j=1,MM
                write(16,"(2010(e13.7,:,','))") trans_film((j-1)*NN+1:(j-1)*NN+NNl,i,k)*conc_l,&
                                                trans_film((j-1)*NN+NNl+1:j*NN,i,k)*conc_g
                write(26,"(2010(e13.7,:,','))") trans_drop((j-1)*NN+1:(j-1)*NN+NNl,i,k)*conc_l,&
                                                trans_drop((j-1)*NN+NNl+1:j*NN,i,k)*conc_g
                write(36,"(2010(e13.7,:,','))") trans_jet((j-1)*NN+1:(j-1)*NN+NNl,i,k)*conc_l,&
                                                trans_jet((j-1)*NN+NNl+1:j*NN,i,k)*conc_g
            end do
    end do
    end do

    close(16)
    close(26)
    close(36)


    do k=1,ncomp
        do i=1,nzl
            do j=1,MM
                call mean_int(trans_film((j-1)*NN+1:(j-1)*NN+NNl,i,k)*conc_l,&
                trans_film((j-1)*NN+NNl+1:j*NN,i,k)*conc_g,&
                ul_film,ug_film,(1-eta_liq)*delta_film+R_c-delta_film,(1-eta_gas)*(R_c-delta_film),IC_l(((i-1)*MM)+j,1,k),&
                IC_g(((i-1)*MM)+j,1,k))
                call mean_int(trans_drop((j-1)*NN+1:(j-1)*NN+NNl,i,k)*conc_l,&
                trans_drop((j-1)*NN+NNl+1:j*NN,i,k)*conc_g,&
                ul_drop,ug_drop,eta_liq*delta_drop,eta_gas*(R_c-delta_drop)+delta_drop,IC_l(((i-1)*MM)+j,2,k),&
                IC_g(((i-1)*MM)+j,2,k))
                call mean_int(trans_jet((j-1)*NN+1:(j-1)*NN+NNl,i,k)*conc_l,&
                trans_jet((j-1)*NN+NNl+1:j*NN,i,k)*conc_g,&
                ul_jet,ug_jet,eta_liq*delta_jet,eta_gas*(R_c-delta_jet)+delta_jet,IC_l(((i-1)*MM)+j,3,k),&
                IC_g(((i-1)*MM)+j,3,k))

            end do
        end do
    end do



    open(13, file='output/film/profilesgasvap.csv', status='replace')  !gas
    open(14, file='output/film/profilesliquid.csv', status='replace')  !liquid

    open(23, file='output/drop/profilesgasvap.csv', status='replace')  !gas
    open(24, file='output/drop/profilesliquid.csv', status='replace')  !liquid

    open(33, file='output/jet/profilesgasvap.csv', status='replace')  !gas
    open(34, file='output/jet/profilesliquid.csv', status='replace')  !liquid

    open(43, file='output/profilesgasvap.csv', status='replace')  !gas
    open(45, file='output/profilesliquid.csv', status='replace')  !liquid

    write(13,'(12A24)') "Height", comp_order
    write(14,'(12A24)') "Height", comp_order

    write(23,'(12A24)') "Height", comp_order
    write(24,'(12A24)') "Height", comp_order

    write(33,'(12A24)') "Height", comp_order
    write(34,'(12A24)') "Height", comp_order

    write(43,'(12A24)') "Height", comp_order
    write(45,'(12A24)') "Height", comp_order

    do i=1,nzl
            do j=1,MM
                write(13,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)), IC_g((((i-1)*MM)+j),1,:)!/sum(IC_g(:,(((i-1)*MM)+j),1)),&
                                                    !dissflux(:,i,2)       !gas concentrations
                write(14,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)), IC_l((((i-1)*MM)+j),1,:)!/sum(IC_l(:,(((i-1)*MM)+j),1)),&
                                               ! dissflux(:,i,1)       !liquid concentrations
                write(23,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)), IC_g((((i-1)*MM)+j),2,:)!/sum(IC_g(:,(((i-1)*MM)+j),2)),&
                                                !   dissflux(:,i,2)       !gas concentratpythonions
                write(24,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)), IC_l((((i-1)*MM)+j),2,:)!/sum(IC_l(:,(((i-1)*MM)+j),2)),&
                                              !  dissflux(:,i,1)       !liquid concentrations
                write(33,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)), IC_g((((i-1)*MM)+j),3,:)!/sum(IC_g(:,(((i-1)*MM)+j),3)),&
                                                  !  dissflux(:,i,2)       !gas concentrations
                write(34,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)), IC_l((((i-1)*MM)+j),3,:)!/sum(IC_l(:,(((i-1)*MM)+j),3)),&
                                                !dissflux(:,i,1)       !liquid concentrations
                write(43,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)),&
                (verh(1)*IC_g((((i-1)*MM)+j),1,:)&!/sum(IC_g(:,(((i-1)*MM)+j),1))
                +verh(2)*IC_g((((i-1)*MM)+j),2,:)&!/sum(IC_g(:,(((i-1)*MM)+j),2))
                +verh(3)*IC_g((((i-1)*MM)+j),3,:))!/sum(IC_g(:,(((i-1)*MM)+j),3)))*conc_g

                write(45,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)), verh(1)*IC_l((((i-1)*MM)+j),1,:)&!/sum(IC_l(:,(((i-1)*MM)+j),1))
                +verh(2)*IC_l((((i-1)*MM)+j),2,:)&!/sum(IC_l(:,(((i-1)*MM)+j),2))
                +verh(3)*IC_l((((i-1)*MM)+j),3,:)!/sum(IC_l(:,(((i-1)*MM)+j),3))
            end do
        end do


    close(13) !gas
    close(14) !liquid

    close(23) !gas
    close(24) !liquid

    close(33) !gas
    close(34) !liquid

    close(43) !gas
    close(45) !liquid


    open(15, file='output/film/Imbalances.csv', status='replace')
    open(25, file='output/drop/Imbalances.csv', status='replace')
    open(35, file='output/jet/Imbalances.csv', status='replace')

    write(15,*) 'Molar imbalances'
    write(25,*) 'Molar imbalances'
    write(35,*) 'Molar imbalances'
    do i=1,ncomp

        write(15,*) comp_order(i), Mim(i), rms_c(i,1)
        write(25,*) comp_order(i), Mim(i), rms_c(i,2)
        write(35,*) comp_order(i), Mim(i), rms_c(i,3)
    end do

    close(15)
    close(25)
    close(35)

!    open(12, file='output/Grid-old.dat', status='replace')
!
!        write(12,*) NNl
!        write(12,*) NNg
!        write(12,*) MM
!    close(12)

!####################################################################################
!   Temperatur Results
!####################################################################################

else

    open(16, file='output/film/radial-temp_profiles.csv', status='replace')

    open(26, file='output/drop/radial-temp_profiles.csv', status='replace')

    open(36, file='output/jet/radial-temp_profiles.csv', status='replace')

    do k=1,ncomp
    write(16,*) comp_order(k)
    write(16,"(2010(e13.7,:,','))") (1-eta_liq(:))*delta_film+R_c-delta_film,(1-eta_gas(:))*(R_c-delta_film)

    write(26,*) comp_order(k)
    write(26,"(2010(e13.7,:,','))") eta_liq(:)*delta_drop,eta_gas(:)*(R_c-delta_drop)+delta_drop

    write(36,*) comp_order(k)
    write(36,"(2010(e13.7,:,','))") eta_liq(:)*delta_jet,eta_gas(:)*(R_c-delta_jet)+delta_jet

    do i=1,nzl
            do j=1,MM
                write(16,"(2010(e13.7,:,','))") trans_film((j-1)*NN+1:(j-1)*NN+NNl,i,k)*transl_0(k),&!*conc_l,&
                                                trans_film((j-1)*NN+NNl+1:j*NN,i,k)*transg_0(k)!*conc_g
                write(26,"(2010(e13.7,:,','))") trans_drop((j-1)*NN+1:(j-1)*NN+NNl,i,k)*transl_0(k),&!*conc_l,&
                                                trans_drop((j-1)*NN+NNl+1:j*NN,i,k)*transg_0(k)!*conc_g
                write(36,"(2010(e13.7,:,','))") trans_jet((j-1)*NN+1:(j-1)*NN+NNl,i,k)*transl_0(k),&!*conc_l,&
                                                trans_jet((j-1)*NN+NNl+1:j*NN,i,k)*transg_0(k)!*conc_g
            end do
    end do
    end do


    close(16)
    close(26)
    close(36)

    do k=1,ncomp
        do i=1,nzl
            do j=1,MM
                call mean_int(trans_film((j-1)*NN+1:(j-1)*NN+NNl,i,k)*transl_0(k),&!*conc_l,&
                            trans_film((j-1)*NN+NNl+1:j*NN,i,k)*transg_0(k),&!*conc_g,&
                            ul_film,ug_film,(1-eta_liq)*delta_film+R_c-delta_film,(1-eta_gas)*(R_c-delta_film),&
                            IC_l(((i-1)*MM)+j,1,k),IC_g(((i-1)*MM)+j,1,k))
                call mean_int(trans_drop((j-1)*NN+1:(j-1)*NN+NNl,i,k)*transl_0(k),&!*conc_l,&
                            trans_drop((j-1)*NN+NNl+1:j*NN,i,k)*transg_0(k),&!*conc_g,&
                            ul_drop,ug_drop,eta_liq*delta_drop,eta_gas*(R_c-delta_drop)+delta_drop,&
                            IC_l(((i-1)*MM)+j,2,k),IC_g(((i-1)*MM)+j,2,k))
                call mean_int(trans_jet((j-1)*NN+1:(j-1)*NN+NNl,i,k)*transl_0(k),&!*conc_l,&
                            trans_jet((j-1)*NN+NNl+1:j*NN,i,k)*transg_0(k),&!*conc_g,&
                            ul_jet,ug_jet,eta_liq*delta_jet,eta_gas*(R_c-delta_jet)+delta_jet,&
                            IC_l(((i-1)*MM)+j,3,k),IC_g(((i-1)*MM)+j,3,k))

            end do
        end do
    end do


    open(13, file='output/film/tempprofilesgasvap.csv', status='replace')  !gas
    open(14, file='output/film/tempprofilesliquid.csv', status='replace')  !liquid

    open(23, file='output/drop/tempprofilesgasvap.csv', status='replace')  !gas
    open(24, file='output/drop/tempprofilesliquid.csv', status='replace')  !liquid

    open(33, file='output/jet/tempprofilesgasvap.csv', status='replace')  !gas
    open(34, file='output/jet/tempprofilesliquid.csv', status='replace')  !liquid

    open(43, file='output/tempprofilesgasvap.csv', status='replace')  !gas
    open(44, file='output/tempprofilesliquid.csv', status='replace')  !liquid

    write(13,'(12A24)') "Height", comp_order
    write(14,'(12A24)') "Height", comp_order

    write(23,'(12A24)') "Height", comp_order
    write(24,'(12A24)') "Height", comp_order

    write(33,'(12A24)') "Height", comp_order
    write(34,'(12A24)') "Height", comp_order

    write(43,'(12A24)') "Height", comp_order
    write(44,'(12A24)') "Height", comp_order

    do i=1,nzl
            do j=1,MM
                write(13,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)), IC_g((((i-1)*MM)+j),1,:)!/sum(IC_g(:,(((i-1)*MM)+j),1)),&
                                                    !dissflux(:,i,2)       !gas concentrations
                write(14,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)), IC_l((((i-1)*MM)+j),1,:)!/sum(IC_l(:,(((i-1)*MM)+j),1)),&
                                               ! dissflux(:,i,1)       !liquid concentrations
                write(23,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)), IC_g((((i-1)*MM)+j),2,:)!/sum(IC_g(:,(((i-1)*MM)+j),2)),&
                                                !   dissflux(:,i,2)       !gas concentratpythonions
                write(24,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)), IC_l((((i-1)*MM)+j),2,:)!/sum(IC_l(:,(((i-1)*MM)+j),2)),&
                                              !  dissflux(:,i,1)       !liquid concentrations
                write(33,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)), IC_g((((i-1)*MM)+j),3,:)!/sum(IC_g(:,(((i-1)*MM)+j),3)),&
                                                  !  dissflux(:,i,2)       !gas concentrations
                write(34,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)), IC_l((((i-1)*MM)+j),3,:)!/sum(IC_l(:,(((i-1)*MM)+j),3)),&
                                                !dissflux(:,i,1)       !liquid concentrations
                write(43,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)),&
                (verh(1)*IC_g((((i-1)*MM)+j),1,:)&!/sum(IC_g(:,(((i-1)*MM)+j),1))
                +verh(2)*IC_g((((i-1)*MM)+j),2,:)&!/sum(IC_g(:,(((i-1)*MM)+j),2))
                +verh(3)*IC_g((((i-1)*MM)+j),3,:))!/sum(IC_g(:,(((i-1)*MM)+j),3)))*conc_g

                write(44,"(2010(e13.7,:,','))") H-(H/(Nzl*MM)*(((i-1)*MM)+j)), verh(1)*IC_l((((i-1)*MM)+j),1,:)&!/sum(IC_l(:,(((i-1)*MM)+j),1))
                +verh(2)*IC_l((((i-1)*MM)+j),2,:)&!/sum(IC_l(:,(((i-1)*MM)+j),2))
                +verh(3)*IC_l((((i-1)*MM)+j),3,:)!/sum(IC_l(:,(((i-1)*MM)+j),3))
            end do
        end do


    close(13) !gas
    close(14) !liquid

    close(23) !gas
    close(24) !liquid

    close(33) !gas
    close(34) !liquid

    close(43) !gas
    close(44) !liquid

!        open(12, file='output/Grid-old.dat', status='replace')
!
!            write(12,*) NNl
!            write(12,*) NNg
!            write(12,*) MM
!        close(12)

end if

 end subroutine



end module tools





