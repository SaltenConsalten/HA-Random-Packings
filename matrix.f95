module matrix
USE global_variables
USE omp_lib
    implicit none

!******************************************************************************************************************
!   This module contains two different ways to set up the coefficient matrix of the model and the right hand-side.
!   Which way is chosen depends on the size of matrix (and therefore the storage demand)
!   and the solution algorithm.
!   The matrix is afterwards multiplied with its transpose to obtain a positive definite and symmetric matrix which
!   is demanded by the Cholesky method as well as by the Conjugate Gradient method.
!******************************************************************************************************************



    contains

    subroutine right_hand_side(NNg, NNl, MM, Ug,dx,c0g,Ul,c0l, b_int,meanmat, val_vec, row, col_ind, n_A,h_r,tcalc)

        integer,intent(in) :: NNg, NNl, MM
        real*8, intent(in) :: Ug(nng), Ul(NNl), dx(2), c0g(NNg), c0l(NNl)

        integer :: i,j,z
        real*8, allocatable :: b(:),Cp(:),Cm(:),c_c(:)
        real*8, intent(in) :: meanmat(3),val_vec(5*NN*MM-2*(NN+(MM-1)))
        integer, intent(in) :: row(5*NN*MM-2*(NN+(MM-1))), col_ind(NN*MM+1)
        real*8, allocatable, intent(out) :: b_int(:)
        real*8, intent(in) :: n_A(mm), h_r(nn*mm)
        real*8 :: nAm(mm), hRm(nn*mm)
        logical, intent(in) :: tcalc

!        real*8, dimension(2) :: deb_out

NN=NNG+NNl

allocate(b(NN*MM),Cp(NN*MM),Cm(NN*MM),b_int(NN*MM),c_c(NN*MM))
b=0
Cp=0
Cm=0
b_int=0
nAm=0
hRm=0


!##############################################
!Set up RHS
!##############################################
  do i=1,MM
    !liquid side
        do j=1,NNl-1

            Cp((i-1)*NN+j)=-Ul(j)/(2*dx(1)*meanmat(1))
            Cm((i-1)*NN+j)=Ul(j)/(2*dx(1)*meanmat(1))

        end do

!gas side of the domain
        do J=1,NNg

            Cp((i-1)*NN+NNl+j)=-Ug(j)/(2*dx(2)*meanmat(2))
            Cm((i-1)*NN+NNl+j)=Ug(j)/(2*dx(2)*meanmat(2))

        end do


 end do
!upside boundary condition

    do i=1,NNl
        b(i)=-Cm(i)*c0l(i)
    end do

 !   do i=NNl+1,NN
  !      Cm(i)=0
  !  end do

!outlet boundary condition zero gradient, nothing changes in outlet
   ! do i=(MM-1)*NN+1,(MM-1)*NN+NNl-1
   !     Cp(i)=0
   ! end do

    do i=1,NNg
        b((MM-1)*NN+NNl+i)=-Cp((MM-1)*NN+NNl+i)*C0g(i)
    end do


!!! Verdampfungsenthalpie   &   Reaktionswärme   !!!

    if(tcalc)then
        do i=1,MM
            nAm(i)=n_A(i)/meanmat(3)
            b((i-1)*NN+NNl)=b((i-1)*NN+NNl)-nAm(i)
        end do

        do i=1,MM
            do z=1,NNl-1
                hRm((i-1)*NN+z)=h_r((i-1)*NN+z)/meanmat(1)
                b((i-1)*NN+z)=b((i-1)*NN+z)-hRm((i-1)*NN+z)
            end do
        end do
    endif



!#####################################################
! Calculate M'*b=b_int
!#####################################################


do i=1,NN*MM
        c_c=0
        do z=col_ind(i)+1,col_ind(i+1)
            c_c(row(z))=val_vec(z)
        end do
        b_int(i)=dot_product(c_c,b)
enddo






    end subroutine right_hand_side



    subroutine FullMatrix(NNg, NNl, MM,R_c,delta, Ug,t_gam,D_g,k_reac,dx,dh,eta_gas,c0g,Ul,Sc_l,D_l,eta_liq,c0l,K,x,y,&
                            b_int,M,recalc,meanmat, val_vec, row, col_ind, n_A,h_r,tcalc)
!****************************************************************
!   Here, the full coefficient matrix is set up. This necessary for the Cholesky decomposition.
!   Use ONLY in case of rather small matrices to prevent excessive usage of RAM
!************************************************************************************************
    implicit none


    integer, intent(in) :: NNg, NNl, MM         !grid nodes
    real*8, intent(in) :: eta_gas(nng), eta_liq(nnl), dx(2), dh(2)      !dimensionless radial/axial coordinates
    real*8, intent(in) :: Ug(nng),Ul(nnl)                 !velocity profiles
    real*8, intent(in) :: t_gam(nng+1),k_reac(3)
    real*8, intent(in) ::  Sc_l,D_g, D_l, R_c,delta               ! gas turbulent sc or pr / liquid sc or pr number
    real*8, intent(in) :: n_A(mm),h_r(nn*mm)           ! temp conductivity, Enthalpie of vap
    Real*8, intent(in) :: c0l(nnl),c0g(nng)
    real*8, intent(in) :: K,x,y             !distribution coefficient
    real*8 :: b((NNl+NNg)*MM)
    real*8, intent(inout) :: val_vec(5*NN*MM-2*(NN+(MM-1)))
    integer,intent(inout) :: row(5*NN*MM-2*(NN+(MM-1))),col_ind(NN*MM+1)
    real*8, intent(inout) :: meanmat(3),M(NN*MM,NN*MM)

    integer :: i,z, j                           !control variables
    real*8, allocatable, intent(out) :: b_int(:)                !coefficient matrix and right hand side

    real*8 :: c_c((NNl+NNg)*MM)
    real*8 :: c_r((NNl+NNg)*MM)
    logical, intent(out) :: recalc

    logical, intent(in) :: tcalc
    !Set up positive definite symmetric Matrix and right hand side for cholesky algorithm

    NN=NNl+NNg
    allocate(b_int((NNl+NNg)*MM))
    M=0
    b=0
    b_int=0
    call initial(NNg, NNl, MM, R_c,delta, Ug,t_gam,D_g,k_reac(:),dx,dh,eta_gas,c0g,Ul,Sc_l,D_l,eta_liq,c0l,K,x,y, val_vec,&
                row,col_ind,b,meanmat, n_A,h_r,tcalc)


    do i=1,NN*MM
        c_c=0
        do z=col_ind(i)+1,col_ind(i+1)
            c_c(row(z))=val_vec(z)
        end do
        b_int(i)=dot_product(c_c,b)
        M(i,i)=dot_product(c_c,c_c)

        if (i.lt.NN*MM-(1)+1) then
        j=i+1
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,j)=dot_product(c_c,c_r)
        M(j,i)=M(i,j)

        end if

        if (i.lt.NN*MM-(2)+1) then
        j=i+2
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,j)=dot_product(c_c,c_r)
        M(j,i)=M(i,j)

        end if
        if (i.lt.NN*MM-(3)+1) then
        j=i+3
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,j)=dot_product(c_c,c_r)
        M(j,i)=M(i,j)

        end if

        if (i.lt.NN*MM-(4)+1) then


        j=i+4
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,j)=dot_product(c_c,c_r)
        M(j,i)=M(i,j)

        end if

        if (i.lt.NN*MM-(NN-2)+1) then
        j=i+(NN-2)
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,j)=dot_product(c_c,c_r)
        M(j,i)=M(i,j)

        end if

        if (i.lt.NN*MM-(NN-1)+1) then
        j=i+(NN-1)
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do
        M(i,j)=dot_product(c_c,c_r)
        M(j,i)=M(i,j)
        end if


        if (i.lt.NN*MM-NN+1) then
        j=i+NN
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do
        M(i,j)=dot_product(c_c,c_r)
        M(j,i)=M(i,j)
        end if



        if (i.lt.NN*MM-(NN+1)+1) then
        j=i+NN+1
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,j)=dot_product(c_c,c_r)
        M(j,i)=M(i,j)
        end if

        if (i.lt.NN*MM-(NN+2)+1) then
        j=i+NN+2;
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,j)=dot_product(c_c,c_r)
        M(j,i)=M(i,j)

        end if

        if (i.lt.NN*MM-NN*2+1) then
        j=i+NN*2
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,j)=dot_product(c_c,c_r)
        M(j,i)=M(i,j)
        end if
    end do

    recalc=.true.


    end subroutine


   subroutine CompMatrix(NNg, NNl, MM, R_c,delta, Ug,t_gam,D_g,k_reac,dx,dh,eta_gas,c0g,Ul,Sc_l,D_l,eta_liq,c0l, K,x,y, b_int,M,&
                            meanmat, val_vec, row, col_ind,n_a,h_r,tcalc)
!*********************************************************************************************
!   Here, the coefficient matrix is calculated in Compressed Sparse Columns (CSC). Because the Matrix is a sparse band matrix
!   only the entries in the bands are stored. This is benefitial for the RAM usage. It also saves computational time
!   for the Conjugate Gradient algorithm.
!*****************************************************************************************************************************

        implicit none


    integer, intent(in) :: NNg, NNl, MM         !grid nodes
    real*8, intent(in) :: eta_gas(NNg), eta_liq(nnl), dx(2), dh(2)      !dimensionless radial/axial coordinates
    real*8, intent(in) :: Ug(nng),Ul(nnl)                 !velocity profiles
    real*8, intent(in) :: t_gam(nng+1), Sc_l,D_g,D_l, R_c,delta,k_reac(2)               ! gas turbulent sc or pr / liquid sc or pr number
    real*8, optional, intent(in) :: n_A(:),h_r(:)
    Real*8, intent(in) :: c0l(NNl),c0g(NNg)
    real*8, intent(in) :: K,x,y             !distribution coefficient
    real*8 :: b((NNl+NNg)*MM)
    real*8 :: val_vec(5*NN*MM-2*(NN+(MM-1)))
    integer :: row(5*NN*MM-2*(NN+(MM-1))), col_ind(NN*MM+1)
    real*8 :: meanmat(3), M(NN*MM,11)
    integer :: NN                               !number of radial nodes
    integer :: i,z, j                           !control variables
    real*8, allocatable, intent(out) :: b_int(:)                !coefficient matrix and right hand side
    logical, intent(in) :: tcalc

    real*8, allocatable :: c_c(:)
    real*8, allocatable :: c_r(:)
    !Set up positive definite symmetric Matrix and right hand side for cholesky algorithm

    NN=NNl+NNg
    allocate(b_int(NN*MM),c_r(NN*MM),c_c(NN*MM))
    M=0
    b=0
    b_int=0
    call initial(NNg, NNl, MM, R_c,delta, Ug,t_gam,D_g,k_reac,dx,dh,eta_gas,c0g,Ul,Sc_l,D_l,eta_liq,c0l,K,x,y, val_vec,&
                row,col_ind,b,meanmat, n_a,h_r,tcalc)

!Set up coefficients for conjugate gradient method
do i=1,NN*MM
        c_c=0
        do z=col_ind(i)+1,col_ind(i+1)
            c_c(row(z))=val_vec(z)
        end do
        b_int(i)=dot_product(c_c,b)
        M(i,1)=dot_product(c_c,c_c)

        if (i.lt.NN*MM-(1)+1) then
        j=i+1
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,2)=dot_product(c_c,c_r)

        end if

        if (i.lt.NN*MM-(2)+1) then
        j=i+2
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,3)=dot_product(c_c,c_r)

        end if
        if (i.lt.NN*MM-(3)+1) then
        j=i+3
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,4)=dot_product(c_c,c_r)

        end if

        if (i.lt.NN*MM-(4)+1) then


        j=i+4
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,5)=dot_product(c_c,c_r)

        end if

        if (i.lt.NN*MM-(NN-2)+1) then
        j=i+(NN-2)
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,6)=dot_product(c_c,c_r)

        end if

        if (i.lt.NN*MM-(NN-1)+1) then
        j=i+(NN-1)
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do
        M(i,7)=dot_product(c_c,c_r)

        end if


        if (i.lt.NN*MM-NN+1) then
        j=i+NN
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do
        M(i,8)=dot_product(c_c,c_r)

        end if



        if (i.lt.NN*MM-(NN+1)+1) then
        j=i+NN+1
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,9)=dot_product(c_c,c_r)

        end if

        if (i.lt.NN*MM-(NN+2)+1) then
        j=i+NN+2;
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,10)=dot_product(c_c,c_r)

        end if

        if (i.lt.NN*MM-NN*2+1) then
        j=i+NN*2
        c_r=0

        do z=col_ind(j)+1,col_ind(j+1)
               c_r(row(z))=val_vec(z)
        end do

        M(i,11)=dot_product(c_c,c_r)

        end if
    end do


    end subroutine CompMatrix





    subroutine Initial(NNg, NNl, MM, R_c,delta, Ug,t_gam,D_g,k_reac,dx,dh,eta_gas,c0g,Ul,Sc_l,D_l,eta_liq,c0l, K,x,y, val_vec,&
                        row, col_ind,b,meanmat, n_A,h_r,tcalc)

        implicit none
        integer, intent(in) :: NNg, NNl, MM         !grid nodes
        real*8, intent(in) :: eta_gas(nng), eta_liq(nnl), dx(2), dh(2)      !dimensionless radial/axial coordinates
        real*8, intent(in) :: K ,x,y             !distribution coefficient
        real*8 :: Ug(nng),Ul(nnl)                 !velocity profiles
        real*8 :: t_gam(NNg+1)
        real*8, intent(in) :: Sc_l,D_g,D_l, R_c,delta               ! gas turbulent sc or pr / liquid sc or pr number
        real*8, intent(in) :: c0l(NNl)
        real*8, intent(in) :: c0g(NNg)
        real*8, intent(out) :: meanmat(3)
        real*8, allocatable :: A(:), Bp(:), Bm(:), Cp(:), Cm(:), Bp2(:), Bm2(:) !Matrix diagonals
        real*8, intent(out) :: b(NN*MM)         !right hand side
        real*8 :: val_vec(5*NN*MM-2*(NN+(MM-1)))
        real*8, dimension(0:nnl) :: e_liq
        real*8, dimension(0:nng+1) :: e_gas
        real*8 :: k_reac(3),a_dg,a_dl, v!,a_k!, k_r(2)
        integer :: row(5*NN*MM-2*(NN+(MM-1))), col_ind(NN*MM+1)
        integer :: i, j, z                           !control variables
        real*8, intent(in) ::  n_A(mm), h_r(nn*mm)                  ! thermal conductivity, Entalpie of Vap
        real*8 :: nAm(mm), hRm(nn*mm)
        logical, intent(in) :: tcalc


        e_liq(1:nnl)=eta_liq
        e_gas(1:nng)=eta_gas
        if (tcalc) then
            k_reac=0D0

                if(eta_liq(nnl).gt.eta_liq(1)) then         !! Drop & Jet
                    e_liq(0)=0
                    e_gas(0)=0
                    e_gas(NNg+1)=1

                    A_dg=2*pi*dh(2)*delta
                    A_dl=-2*pi*dh(1)*(delta*eta_liq(nnl-1))
                    !A_k=(delta**2-(delta*(1-eta_liq(nnl-1)))**2)*pi
                    v=0D0!(delta**2-(delta*(eta_liq(nnl-1)))**2)*pi*dh(1)
!
!                    A_dl=-A_dl!/A_dg
!        !            A_k=A_k!/A_dg
!                    v=v/A_dg
!                    A_dg=1

                else                    !! Film
                    e_liq(0)=1
                    e_gas(0)=1
                    e_gas(NNg+1)=0

                    A_dg=-2*pi*dh(2)*(R_c-delta)
                    A_dl=2*pi*dh(1)*((R_c-delta)+delta*e_liq(nnl-1))
        !            A_k=pi*((R_c-delta+delta*e_liq(nnl-1))**2-(R_c-delta)**2)
                    v=0D0!pi*((R_c-delta+delta*e_liq(nnl-1))**2-(R_c-delta)**2)*dh(1)

        !            A_dl=A_dl!/A_dg
        !            A_k=A_k!/A_dg
                    !v=v/A_dg
        !            A_dg=-1

                end if

            else
                k_reac=0D0
                k_reac(2)=k_reac(1)*delta**2/(Sc_l*D_l)

                if(eta_liq(nnl).gt.eta_liq(1)) then         !! Drop & Jet
                    e_liq(0)=0
                    e_gas(0)=0
                    e_gas(NNg+1)=1

                    A_dg=2*pi*dh(2)*delta
                    A_dl=2*pi*dh(1)*(delta*eta_liq(nnl-1))
        !            A_k=(delta**2-(delta*(1-eta_liq(nnl-1)))**2)*pi
                    v=(delta**2-(delta*(eta_liq(nnl-1)))**2)*pi*dh(1)

                    A_dl=-A_dl/A_dg
!                    A_k=A_k!/A_dg
                    v=v/A_dg
                    A_dg=1

                else                    !! Film
                    e_liq(0)=1
                    e_gas(0)=1
                    e_gas(NNg+1)=0

                    A_dg=2*pi*dh(2)*(R_c-delta)
                    A_dl=2*pi*dh(1)*((R_c-delta)+delta*e_liq(nnl-1))
        !            A_k=pi*((R_c-delta+delta*e_liq(nnl-1))**2-(R_c-delta)**2)
                    v=pi*((R_c-delta+delta*e_liq(nnl-1))**2-(R_c-delta)**2)*dh(1)

                    A_dl=A_dl/A_dg
        !            A_k=A_k!/A_dg
                    v=v/A_dg
                    A_dg=-1

                end if
        end if

!        if(eta_liq(nnl).gt.eta_liq(1)) then         !! Drop & Jet
!            e_liq(0)=0
!            e_gas(0)=0
!            e_gas(NNg+1)=1
!
!            A_dg=2*pi*dh(2)*delta
!            A_dl=-2*pi*dh(1)*(delta*eta_liq(nnl-1))
!!            A_dg=2*pi*dh(2)*delta
!!            A_dl=2*pi*dh(1)*delta
!!            A_k=(delta**2-(delta*(1-eta_liq(nnl-1)))**2)*pi
!            v=(delta**2-(delta*(eta_liq(nnl-1)))**2)*pi*dh(1)
!
!!            A_dl=-A_dl!/A_dg
!!!            A_k=A_k!/A_dg
!            v=v/A_dg
!!            A_dg=1
!
!        else                    !! Film
!            e_liq(0)=1
!            e_gas(0)=1
!            e_gas(NNg+1)=0
!
!            A_dg=-2*pi*dh(2)*(R_c-delta)
!            A_dl=2*pi*dh(1)*((R_c-delta)+delta*e_liq(nnl-1))
!!            A_dg=2*pi*dh(2)*(R_c-delta)
!!            A_dl=2*pi*dh(1)*(R_c-delta)
!!            A_k=pi*((R_c-delta+delta*e_liq(nnl-1))**2-(R_c-delta)**2)
!            v=pi*((R_c-delta+delta*e_liq(nnl-1))**2-(R_c-delta)**2)*dh(1)
!
!!            A_dl=A_dl!/A_dg
!!            A_k=A_k!/A_dg
!            v=v/A_dg
!!            A_dg=-1

!        end if


        NN=NNg+NNl


        !k_reac(1)=0
        allocate(Cp(NN*MM),Cm(NN*MM),Bp2(NN*MM),Bm2(NN*MM),A(NN*MM),Bp(NN*MM),Bm(NN*MM))

        A=0
        Bm=0
        Bp=0
        Cp=0
        Cm=0
        Bp2=0
        Bm2=0
        row=0
        col_ind=0
        val_vec=0
        nAm=0
        hRm=0




!!! liquid side of the domain
    do i=1,MM
        do j=1,NNl-1

            A((i-1)*NN+j)=Sc_l**(-1)/(2*e_liq(j))*(-(e_liq(j)+e_liq(j+1))/((e_liq(j+1)-e_liq(j))*(e_liq(j+1)-&
                        e_liq(j-1))/2)-(e_liq(j-1)+e_liq(j))/((e_liq(j)-e_liq(j-1))*(e_liq(j+1)-e_liq(j-1))/2))-k_reac(2)
            Bp((i-1)*NN+j)=Sc_l**(-1)/(2*e_liq(j))*(e_liq(j)+e_liq(j+1))/((e_liq(j+1)-e_liq(j))*(e_liq(j+1)-&
                        e_liq(j-1))/2)
            Bm((i-1)*NN+j)=Sc_l**(-1)/(2*e_liq(j))*(e_liq(j-1)+e_liq(j))/((e_liq(j)-e_liq(j-1))*(e_liq(j+1)-&
                            e_liq(j-1))/2)
            Cp((i-1)*NN+j)=-Ul(j)/(2*dx(1))
            Cm((i-1)*NN+j)=Ul(j)/(2*dx(1))

        end do

!!!!!!! Interface liquid side

        Bp((i-1)*NN+NNl)=A_dg*y*D_g/(R_c-delta)*(e_gas(2)-e_gas(0))/((e_gas(2)-e_gas(1))*(e_gas(1)-e_gas(0)))     !
        Bp2((i-1)*NN+NNl)=-A_dg*y*D_g/(R_c-delta)*(e_gas(1)-e_gas(0))/((e_gas(2)-e_gas(1))*(e_gas(2)-e_gas(0))) !
        Bm((i-1)*NN+NNl)=A_dl*x*D_l/delta*(e_liq(NNl-2)-e_liq(NNl))/((e_liq(NNl-2)-e_liq(NNl-1))*(e_liq(NNl-1)-e_liq(NNl)))!


        Bm2((i-1)*NN+NNl)=-A_dl*x*D_l/delta*(e_liq(NNl-1)-e_liq(NNl))/((e_liq(NNl-2)-e_liq(NNl-1))*(e_liq(NNl-2)-&!
        e_liq(NNl)))

        A((i-1)*NN+NNl)=A_dg*k*x*D_g/(R_c-delta)*((e_gas(1)-e_gas(0))**2-(e_gas(2)-e_gas(0))**2)/&!
                                                ((e_gas(2)-e_gas(1))*(e_gas(2)-e_gas(0))*(e_gas(1)-e_gas(0)))&
                        +A_dl*x*D_l/delta*((e_liq(NNl-1)-e_liq(NNl))**2-(e_liq(NNl-2)-e_liq(NNl))**2)/&
                        ((e_liq(NNl-2)-e_liq(NNl-1))*(e_liq(NNl-2)-e_liq(NNl))*(e_liq(NNl-1)-e_liq(NNl)))-x*k_reac(1)*v

        Cm((i-1)*NN+NNl)=0!x*Ul((NNl))*a_k
        Cp((i-1)*NN+NNl)=0!-x*Ul((NNl))*a_k


!!!!!!! inteface gas side
        A((i-1)*NN+NNl+1)=t_gam(2)/(2*e_gas(1))*(-(e_gas(1)+e_gas(2))/((e_gas(2)-e_gas(1))*(e_gas(2)-&
                        e_gas(0))/2))-t_gam(1)/(2*e_gas(1))*((e_gas(0)+e_gas(1))/((e_gas(1)-e_gas(0))*(e_gas(2)-e_gas(0))/2))
        Bp((i-1)*NN+NNl+1)=t_gam(2)/(2*e_gas(1))*(e_gas(1)+e_gas(2))/((e_gas(2)-e_gas(1))*(e_gas(2)-&
                        e_gas(0))/2)
        Bm((i-1)*NN+NNl+1)=x/y*k*t_gam(1)/(2*e_gas(1))*(e_gas(0)+e_gas(1))/((e_gas(1)-e_gas(0))*(e_gas(2)-e_gas(0))/2)!K*x/y*
        Cp((i-1)*NN+NNl+1)=-Ug(1)/(2*dx(2))
        Cm((i-1)*NN+NNl+1)=Ug(1)/(2*dx(2))
        !b((i-1)*NN+NNl+1)=-Bm((i-1)*NN+NNl+1)
        !Bm((i-1)*NN+NNl+1)=0


!!! gas side of the domain
        do j=2,NNg

            A((i-1)*NN+NNl+j)=t_gam(j+1)/(2*e_gas(j))*(-(e_gas(j)+e_gas(j+1))/((e_gas(j+1)-e_gas(j))*(e_gas(j+1)-e_gas(j-1))/2))&
                            -t_gam(j)/(2*e_gas(j))*((e_gas(j-1)+e_gas(j))/((e_gas(j)-e_gas(j-1))*(e_gas(j+1)-e_gas(j-1))/2))
            Bp((i-1)*NN+NNl+j)=t_gam(j+1)/(2*e_gas(j))*(e_gas(j)+e_gas(j+1))/((e_gas(j+1)-e_gas(j))*(e_gas(j+1)-&
                        e_gas(j-1))/2)
            Bm((i-1)*NN+NNl+j)=t_gam(j)/(2*e_gas(j))*(e_gas(j-1)+e_gas(j))/((e_gas(j)-e_gas(j-1))*(e_gas(j+1)-&
                            e_gas(j-1))/2)
            Cp((i-1)*NN+NNl+j)=-Ug(j)/(2*dx(2))
            Cm((i-1)*NN+NNl+j)=Ug(j)/(2*dx(2))

        end do


    end do

    !###############################################################################
    ! start preconditioning
    !###############################################################################


!if (tcalc) then
        meanmat=0
!else
!        meanmat=1D0
!end if

    do i=1,MM

        if (meanmat(1).lt.maxval(max(abs(A((i-1)*NN+1:(i-1)*NN+NNl-1)),abs(Bm((i-1)*NN+1:(i-1)*NN+NNl-1)),abs(Bp((i-1)*NN+1:(i-1)*&
                NN+NNl-1)),abs(Cm((i-1)*NN+1:(i-1)*NN+NNl-1)),abs(Cp((i-1)*NN+1:(i-1)*NN+NNl-1))))) then

            meanmat(1)=maxval(max(abs(A((i-1)*NN+1:(i-1)*NN+NNl-1)),abs(Bm((i-1)*NN+1:(i-1)*NN+NNl-1)),abs(Bp((i-1)*NN+1:(i-1)*&
                NN+NNl-1)),abs(Cm((i-1)*NN+1:(i-1)*NN+NNl-1)),abs(Cp((i-1)*NN+1:(i-1)*NN+NNl-1))))

        end if
        if (meanmat(2).lt.maxval(max(abs(A((i-1)*NN+NNl+1:i*NN)),abs(Bm((i-1)*NN+NNl+1:i*NN)),abs(Bp((i-1)*NN+NNl+1:i*NN)),abs(Cm(&
                    (i-1)*NN+NNl+1:i*NN)),abs(Cp((i-1)*NN+NNl+1:i*NN))))) then

            meanmat(2)=maxval(max(abs(A((i-1)*NN+NNl+1:i*NN)),abs(Bm((i-1)*NN+NNl+1:i*NN)),abs(Bp((i-1)*NN+NNl+1:i*NN)),abs(Cm(&
                    (i-1)*NN+NNl+1:i*NN)),abs(Cp((i-1)*NN+NNl+1:i*NN))))

        end if

        if (meanmat(3).lt.max(abs(A((i-1)*NN+NNl)),abs(Bm((i-1)*NN+NNl)),abs(Bp((i-1)*NN+NNl)),abs(Bm2((i-1)*NN+NNl)),abs(&
                    Bp2((i-1)*NN+NNl)))) then

            meanmat(3)=max(abs(A((i-1)*NN+NNl)),abs(Bm((i-1)*NN+NNl)),abs(Bp((i-1)*NN+NNl)),abs(Bm2((i-1)*NN+NNl)),abs(&
                    Bp2((i-1)*NN+NNl)))
        end if

    end do

    do i=1,MM

        A((i-1)*NN+1:(i-1)*NN+NNl-1)=A((i-1)*NN+1:(i-1)*NN+NNl-1)/meanmat(1)
        Bp((i-1)*NN+1:(i-1)*NN+NNl-1)=Bp((i-1)*NN+1:(i-1)*NN+NNl-1)/meanmat(1)
        Bm((i-1)*NN+1:(i-1)*NN+NNl-1)=Bm((i-1)*NN+1:(i-1)*NN+NNl-1)/meanmat(1)
        Cp((i-1)*NN+1:(i-1)*NN+NNl-1)=Cp((i-1)*NN+1:(i-1)*NN+NNl-1)/meanmat(1)
        Cm((i-1)*NN+1:(i-1)*NN+NNl-1)=Cm((i-1)*NN+1:(i-1)*NN+NNl-1)/meanmat(1)

        A((i-1)*NN+NNl+1:i*NN)=A((i-1)*NN+NNl+1:i*NN)/meanmat(2)
        Bp((i-1)*NN+NNl+1:i*NN)=Bp((i-1)*NN+NNl+1:i*NN)/meanmat(2)
        Bm((i-1)*NN+NNl+1:i*NN)=Bm((i-1)*NN+NNl+1:i*NN)/meanmat(2)
        Cp((i-1)*NN+NNl+1:i*NN)=Cp((i-1)*NN+NNl+1:i*NN)/meanmat(2)
        Cm((i-1)*NN+NNl+1:i*NN)=Cm((i-1)*NN+NNl+1:i*NN)/meanmat(2)

        A((i-1)*NN+NNl)=A((i-1)*NN+NNl)/meanmat(3)
        Bp((i-1)*NN+NNl)=Bp((i-1)*NN+NNl)/meanmat(3)
        Bm((i-1)*NN+NNl)=Bm((i-1)*NN+NNl)/meanmat(3)
        Bp2((i-1)*NN+NNl)=Bp2((i-1)*NN+NNl)/meanmat(3)
        Bm2((i-1)*NN+NNl)=Bm2((i-1)*NN+NNl)/meanmat(3)

    end do





    !###############################################################################
    ! end preconditioning
    !###############################################################################


!upside boundary condition

    do i=1,NNl
        b(i)=-Cm(i)*c0l(i)
        Cm(i)=0
    end do

    do i=NNl+1,NN
        A(i)=A(i)+Cm(i)
        Cm(i)=0
    end do

!Wall boundary condition: zero gradient
    do i=0,MM-1
        A(i*NN+1)=A(i*NN+1)+Bm(i*NN+1)
        Bm(i*NN+1)=0
    end do

!outlet boundary condition zero gradient, nothing changes in outlet
    do i=(MM-1)*NN+1,(MM-1)*NN+NNl-1
        A(i)=A(i)+Cp(i)
        Cp(i)=0
    end do

    do i=(MM-1)*NN+NNl+1,(MM)*NN
        b(i)=-Cp(i)*C0g(i-((MM-1)*NN+NNl+1)+1)
        Cp(i)=0
    end do


!Pipe centre boundary condition C(NN)=K*C2
    do i=1,MM
        A(i*NN)=A(i*NN)+Bp(i*NN)
        Bp(i*NN)=0
    end do

    if(tcalc)then
        do i=1,MM
            nAm(i)=n_A(i)/meanmat(3)
            b((i-1)*NN+NNl)=b((i-1)*NN+NNl)-nAm(i)
        end do

        do i=1,MM
            do z=1,NNl-1
                hRm((i-1)*NN+z)=h_r((i-1)*NN+z)/meanmat(1)
                b((i-1)*NN+z)=b((i-1)*NN+z)-hRm((i-1)*NN+z)
            end do
        end do
    endif



    do i=1,NN*MM

        if ((i.le.NN).or.(i.gt.NN*(MM-1))) then
            if ((i.eq.1).or.(i.eq.NN*MM).or.(i.eq.NN).or.(i.eq.NN*(MM-1)+1).or.(i.eq.NNL).or.(i.eq.NN*(MM-1)+NNL)) then

                col_ind(i+1)=col_ind(i)+3
            elseif ((i.eq.NNL-2).or.(i.eq.NNL+2).or.(i.eq.NN*(MM-1)+NNL-2).or.(i.eq.NN*(MM-1)+NNL+2)) then
                col_ind(i+1)=col_ind(i)+5
            else
                col_ind(i+1)=col_ind(i)+4
            end if
        else
            if ((mod(i,NN).eq.0).or.(mod((i-1),NN).eq.0)) then
                col_ind(i+1)=col_ind(i)+4
            elseif ((mod(i-(NNL-2),NN).eq.0).or.(mod(i-(NNL+2),NN).eq.0)) then

            col_ind(i+1)=col_ind(i)+6
            elseif (mod(i-NNL,NN).eq.0) then
                col_ind(i+1)=col_ind(i)+3
            else
                col_ind(i+1)=col_ind(i)+5
            end if
        end if
    end do



    val_vec(1)=A(1)
    row(1)=1
    val_vec(2)=Bm(2)
    row(2)=2
    val_vec(3)=Cm(NN+1)
    row(3)=NN+1
    val_vec(size(val_vec))=A(NN*MM)
    row(size(row))=NN*MM
    val_vec(size(val_vec)-1)=Bp(NN*MM-1)
    row(size(row)-1)=NN*MM-1
    val_vec(size(val_vec)-2)=Cp(NN*(MM-1))
    row(size(row)-2)=NN*(MM-1)

    do i=2,NN*MM-1

        z=col_ind(i+1)-col_ind(i)

        if (z.eq.3) then
            if (Bm2(i).ne.0) then
        val_vec(col_ind(i)+1)=Bp(i-1)
        val_vec(col_ind(i)+2)=A(i)
        val_vec(col_ind(i)+3)=Bm(i+1)

        row(col_ind(i)+1)=i-1
        row(col_ind(i)+2)=i
        row(col_ind(i)+3)=i+1

            elseif (i.eq.NN) then

        val_vec(col_ind(i)+1)=Bp(i-1)
        val_vec(col_ind(i)+2)=A(i)
        val_vec(col_ind(i)+3)=Cm(i+NN)

        row(col_ind(i)+1)=i-1
        row(col_ind(i)+2)=i
        row(col_ind(i)+3)=i+NN
            else
        val_vec(col_ind(i)+1)=Cp(i-NN)
        val_vec(col_ind(i)+2)=A(i)
        val_vec(col_ind(i)+3)=Bm(i+1)

        row(col_ind(i)+1)=i-NN
        row(col_ind(i)+2)=i
        row(col_ind(i)+3)=i+1
            end if


        elseif (z.eq.4) then
            if (i.le.NN) then
        val_vec(col_ind(i)+1)=Bp(i-1)
        val_vec(col_ind(i)+2)=A(i)
        val_vec(col_ind(i)+3)=Bm(i+1)
        val_vec(col_ind(i)+4)=Cm(i+NN)

        row(col_ind(i)+1)=i-1
        row(col_ind(i)+2)=i
        row(col_ind(i)+3)=i+1
        row(col_ind(i)+4)=i+NN
            elseif (i.gt.NN*(MM-1)) then !mod(i,NNL+NNG)==0 || mod(i-1,NNL+NNG)==0

        val_vec(col_ind(i)+1)=Cp(i-NN)
        val_vec(col_ind(i)+2)=Bp(i-1)
        val_vec(col_ind(i)+3)=A(i)
        val_vec(col_ind(i)+4)=Bm(i+1)

        row(col_ind(i)+1)=i-NN
        row(col_ind(i)+2)=i-1
        row(col_ind(i)+3)=i
        row(col_ind(i)+4)=i+1
            elseif (Bm(i+1).eq.0) then

        val_vec(col_ind(i)+1)=Cp(i-NN)
        val_vec(col_ind(i)+2)=Bp(i-1)
        val_vec(col_ind(i)+3)=A(i)
        val_vec(col_ind(i)+4)=Cm(i+NN)

        row(col_ind(i)+1)=i-NN
        row(col_ind(i)+2)=i-1
        row(col_ind(i)+3)=i
        row(col_ind(i)+4)=i+NN
            else

        val_vec(col_ind(i)+1)=Cp(i-NN)
        val_vec(col_ind(i)+2)=A(i)
        val_vec(col_ind(i)+3)=Bm(i+1)
        val_vec(col_ind(i)+4)=Cm(i+NN)

        row(col_ind(i)+1)=i-NN
        row(col_ind(i)+2)=i
        row(col_ind(i)+3)=i+1
        row(col_ind(i)+4)=i+NN
            end if

         elseif (z.eq.5) then
             if (i.eq.NNL-2) then
        val_vec(col_ind(i)+1)=Bp(i-1)
        val_vec(col_ind(i)+2)=A(i)
        val_vec(col_ind(i)+3)=Bm(i+1)
        val_vec(col_ind(i)+4)=Bm2(i+2)
        val_vec(col_ind(i)+5)=Cm(i+NN)

        row(col_ind(i)+1)=i-1
        row(col_ind(i)+2)=i
        row(col_ind(i)+3)=i+1
        row(col_ind(i)+4)=i+2
        row(col_ind(i)+5)=i+NN

             elseif (i.eq.NNL+2) then
        val_vec(col_ind(i)+1)=Bp2(i-2)
        val_vec(col_ind(i)+2)=Bp(i-1)
        val_vec(col_ind(i)+3)=A(i)
        val_vec(col_ind(i)+4)=Bm(i+1)
        val_vec(col_ind(i)+5)=Cm(i+NN)

        row(col_ind(i)+1)=i-2
        row(col_ind(i)+2)=i-1
        row(col_ind(i)+3)=i
        row(col_ind(i)+4)=i+1
        row(col_ind(i)+5)=i+NN
             elseif (i.eq.NN*(MM-1)+NNL-2) then


        val_vec(col_ind(i)+1)=Cp(i-NN)
        val_vec(col_ind(i)+2)=Bp(i-1)
        val_vec(col_ind(i)+3)=A(i)
        val_vec(col_ind(i)+4)=Bm(i+1)
        val_vec(col_ind(i)+5)=Bm2(i+2)

        row(col_ind(i)+1)=i-NN
        row(col_ind(i)+2)=i-1
        row(col_ind(i)+3)=i
        row(col_ind(i)+4)=i+1
        row(col_ind(i)+5)=i+2

               elseif (i.eq.NN*(MM-1)+NNL+2)  then
        val_vec(col_ind(i)+1)=Cp(i-NN)
        val_vec(col_ind(i)+2)=Bp2(i-2)
        val_vec(col_ind(i)+3)=Bp(i-1)
        val_vec(col_ind(i)+4)=A(i)
        val_vec(col_ind(i)+5)=Bm(i+1)

        row(col_ind(i)+1)=i-NN
        row(col_ind(i)+2)=i-2
        row(col_ind(i)+3)=i-1
        row(col_ind(i)+4)=i
        row(col_ind(i)+5)=i+1

             else

        val_vec(col_ind(i)+1)=Cp(i-NN)
        val_vec(col_ind(i)+2)=Bp(i-1)
        val_vec(col_ind(i)+3)=A(i)
        val_vec(col_ind(i)+4)=Bm(i+1)
        val_vec(col_ind(i)+5)=Cm(i+NN)

        row(col_ind(i)+1)=i-NN
        row(col_ind(i)+2)=i-1
        row(col_ind(i)+3)=i
        row(col_ind(i)+4)=i+1
        row(col_ind(i)+5)=i+NN

             end if
         elseif (z.eq.6) then
             if (mod(i+2-NNL,NN).eq.0) then
        val_vec(col_ind(i)+1)=Cp(i-NN)
        val_vec(col_ind(i)+2)=Bp(i-1)
        val_vec(col_ind(i)+3)=A(i)
        val_vec(col_ind(i)+4)=Bm(i+1)
        val_vec(col_ind(i)+5)=Bm2(i+2)
        val_vec(col_ind(i)+6)=Cm(i+NN)

        row(col_ind(i)+1)=i-NN
        row(col_ind(i)+2)=i-1
        row(col_ind(i)+3)=i
        row(col_ind(i)+4)=i+1
        row(col_ind(i)+5)=i+2
        row(col_ind(i)+6)=i+NN

             else
        val_vec(col_ind(i)+1)=Cp(i-NN)
        val_vec(col_ind(i)+2)=Bp2(i-2)
        val_vec(col_ind(i)+3)=Bp(i-1)
        val_vec(col_ind(i)+4)=A(i)
        val_vec(col_ind(i)+5)=Bm(i+1)
        val_vec(col_ind(i)+6)=Cm(i+NN)

        row(col_ind(i)+1)=i-NN
        row(col_ind(i)+2)=i-2
        row(col_ind(i)+3)=i-1
        row(col_ind(i)+4)=i
        row(col_ind(i)+5)=i+1
        row(col_ind(i)+6)=i+NN
             end if
         end if
    end do

    deallocate(A,Bp,Bm,Cp,Cm,Bp2,Bm2)

    end subroutine Initial

    subroutine precond(NN,MM,M,L)
implicit none

integer, intent(in) :: NN, MM
real*8, intent (in) :: M(:,:)
real*8 :: row_M(2*NN+NN*MM), row_N(2*NN+NN*MM)
!real*8 :: tmp_dot, tmp_L, tmp
real*8 :: alpha=1!(3)
integer :: i,k=0,j!, ll
real*8, intent (inout) :: L(NN*MM+2*NN,11)

row_n=0
row_m=0
L(:,:)=0

10 L(2*NN+1:NN*MM+2*NN,:)=M

forall (i=2*NN+1:nn*mm+2*nn) L(i,1)=L(i,1)*alpha





do i=2*NN+1,NN*MM+2*NN
    do j=2*NN+1,i


if(j.eq.i) then
row_M=0
row_M(i)=L(i,1)
row_M(i-1)=L(i-1,2)
row_M(i-2)=L(i-2,3)
row_M(i-3)=L(i-3,4)
row_M(i-4)=L(i-4,5)
row_M(i-(NN-2))=L(i-(NN-2),6)
row_M(i-(NN-1))=L(i-(NN-1),7)
row_M(i-NN)=L(i-(NN),8)
row_M(i-(NN+1))=L(i-(NN+1),9)
row_M(i-(NN+2))=L(i-(NN+2),10)
row_M(i-2*NN)=L(i-2*NN,11)

    L(i,1)=dsqrt(L(i,1)-dot_product(row_M(2*nn+1:i-1),row_M(2*nn+1:i-1)))
elseif (((j+1).eq.i)) then

row_M=0

row_M(i)=L(i,1)
row_M(i-1)=L(i-1,2)
row_M(i-2)=L(i-2,3)
row_M(i-3)=L(i-3,4)
row_M(i-4)=L(i-4,5)
row_M(i-(NN-2))=L(i-(NN-2),6)
row_M(i-(NN-1))=L(i-(NN-1),7)
row_M(i-NN)=L(i-(NN),8)
row_M(i-(NN+1))=L(i-(NN+1),9)
row_M(i-(NN+2))=L(i-(NN+2),10)
row_M(i-2*NN)=L(i-2*NN,11)
        row_N=0
    row_N(j)=L(j,1)
    row_N(j-1)=L(j-1,2)
    row_N(j-2)=L(j-2,3)
    row_N(j-3)=L(j-3,4)
    row_N(j-4)=L(j-4,5)
    row_N(j-(NN-2))=L(j-(NN-2),6)
    row_N(j-(NN-1))=L(j-(NN-1),7)
    row_N(j-NN)=L(j-(NN),8)
    row_N(j-(NN+1))=L(j-(NN+1),9)
    row_N(j-(NN+2))=L(j-(NN+2),10)
    row_N(j-2*NN)=L(j-2*NN,11)
    L(j,2)=1D0/L(j,1)*(L(j,2)-dot_product(row_M(2*nn+1:j-1),row_N(2*nn+1:j-1)))
elseif ((j+2).eq.i) then

row_M=0

row_M(i)=L(i,1)
row_M(i-1)=L(i-1,2)
row_M(i-2)=L(i-2,3)
row_M(i-3)=L(i-3,4)
row_M(i-4)=L(i-4,5)
row_M(i-(NN-2))=L(i-(NN-2),6)
row_M(i-(NN-1))=L(i-(NN-1),7)
row_M(i-NN)=L(i-(NN),8)
row_M(i-(NN+1))=L(i-(NN+1),9)
row_M(i-(NN+2))=L(i-(NN+2),10)
row_M(i-2*NN)=L(i-2*NN,11)
        row_N=0
    row_N(j)=L(j,1)
    row_N(j-1)=L(j-1,2)
    row_N(j-2)=L(j-2,3)
    row_N(j-3)=L(j-3,4)
    row_N(j-4)=L(j-4,5)
    row_N(j-(NN-2))=L(j-(NN-2),6)
    row_N(j-(NN-1))=L(j-(NN-1),7)
    row_N(j-NN)=L(j-(NN),8)
    row_N(j-(NN+1))=L(j-(NN+1),9)
    row_N(j-(NN+2))=L(j-(NN+2),10)
    row_N(j-2*NN)=L(j-2*NN,11)
    L(j,3)=1D0/L(j,1)*(L(j,3)-dot_product(row_M(2*nn+1:j-1),row_N(2*nn+1:j-1)))
elseif ((j+3).eq.i) then

row_M=0

row_M(i)=L(i,1)
row_M(i-1)=L(i-1,2)
row_M(i-2)=L(i-2,3)
row_M(i-3)=L(i-3,4)
row_M(i-4)=L(i-4,5)
row_M(i-(NN-2))=L(i-(NN-2),6)
row_M(i-(NN-1))=L(i-(NN-1),7)
row_M(i-NN)=L(i-(NN),8)
row_M(i-(NN+1))=L(i-(NN+1),9)
row_M(i-(NN+2))=L(i-(NN+2),10)
row_M(i-2*NN)=L(i-2*NN,11)
        row_N=0
    row_N(j)=L(j,1)
    row_N(j-1)=L(j-1,2)
    row_N(j-2)=L(j-2,3)
    row_N(j-3)=L(j-3,4)
    row_N(j-4)=L(j-4,5)
    row_N(j-(NN-2))=L(j-(NN-2),6)
    row_N(j-(NN-1))=L(j-(NN-1),7)
    row_N(j-NN)=L(j-(NN),8)
    row_N(j-(NN+1))=L(j-(NN+1),9)
    row_N(j-(NN+2))=L(j-(NN+2),10)
    row_N(j-2*NN)=L(j-2*NN,11)
    L(j,4)=1D0/L(j,1)*(L(j,4)-dot_product(row_M(2*nn+1:j-1),row_N(2*nn+1:j-1)))
elseif ((j+4).eq.i) then

row_M=0

row_M(i)=L(i,1)
row_M(i-1)=L(i-1,2)
row_M(i-2)=L(i-2,3)
row_M(i-3)=L(i-3,4)
row_M(i-4)=L(i-4,5)
row_M(i-(NN-2))=L(i-(NN-2),6)
row_M(i-(NN-1))=L(i-(NN-1),7)
row_M(i-NN)=L(i-(NN),8)
row_M(i-(NN+1))=L(i-(NN+1),9)
row_M(i-(NN+2))=L(i-(NN+2),10)
row_M(i-2*NN)=L(i-2*NN,11)
        row_N=0
    row_N(j)=L(j,1)
    row_N(j-1)=L(j-1,2)
    row_N(j-2)=L(j-2,3)
    row_N(j-3)=L(j-3,4)
    row_N(j-4)=L(j-4,5)
    row_N(j-(NN-2))=L(j-(NN-2),6)
    row_N(j-(NN-1))=L(j-(NN-1),7)
    row_N(j-NN)=L(j-(NN),8)
    row_N(j-(NN+1))=L(j-(NN+1),9)
    row_N(j-(NN+2))=L(j-(NN+2),10)
    row_N(j-2*NN)=L(j-2*NN,11)
    L(j,5)=1D0/L(j,1)*(L(j,5)-dot_product(row_M(2*nn+1:j-1),row_N(2*nn+1:j-1)))
elseif ((j+NN-2).eq.i) then

row_M=0

row_M(i)=L(i,1)
row_M(i-1)=L(i-1,2)
row_M(i-2)=L(i-2,3)
row_M(i-3)=L(i-3,4)
row_M(i-4)=L(i-4,5)
row_M(i-(NN-2))=L(i-(NN-2),6)
row_M(i-(NN-1))=L(i-(NN-1),7)
row_M(i-NN)=L(i-(NN),8)
row_M(i-(NN+1))=L(i-(NN+1),9)
row_M(i-(NN+2))=L(i-(NN+2),10)
row_M(i-2*NN)=L(i-2*NN,11)
        row_N=0
    row_N(j)=L(j,1)
    row_N(j-1)=L(j-1,2)
    row_N(j-2)=L(j-2,3)
    row_N(j-3)=L(j-3,4)
    row_N(j-4)=L(j-4,5)
    row_N(j-(NN-2))=L(j-(NN-2),6)
    row_N(j-(NN-1))=L(j-(NN-1),7)
    row_N(j-NN)=L(j-(NN),8)
    row_N(j-(NN+1))=L(j-(NN+1),9)
    row_N(j-(NN+2))=L(j-(NN+2),10)
    row_N(j-2*NN)=L(j-2*NN,11)
    L(j,6)=1D0/L(j,1)*(L(j,6)-dot_product(row_M(2*nn+1:j-1),row_N(2*nn+1:j-1)))
elseif ((j+NN-1).eq.i) then

row_M=0

row_M(i)=L(i,1)
row_M(i-1)=L(i-1,2)
row_M(i-2)=L(i-2,3)
row_M(i-3)=L(i-3,4)
row_M(i-4)=L(i-4,5)
row_M(i-(NN-2))=L(i-(NN-2),6)
row_M(i-(NN-1))=L(i-(NN-1),7)
row_M(i-NN)=L(i-(NN),8)
row_M(i-(NN+1))=L(i-(NN+1),9)
row_M(i-(NN+2))=L(i-(NN+2),10)
row_M(i-2*NN)=L(i-2*NN,11)
        row_N=0
    row_N(j)=L(j,1)
    row_N(j-1)=L(j-1,2)
    row_N(j-2)=L(j-2,3)
    row_N(j-3)=L(j-3,4)
    row_N(j-4)=L(j-4,5)
    row_N(j-(NN-2))=L(j-(NN-2),6)
    row_N(j-(NN-1))=L(j-(NN-1),7)
    row_N(j-NN)=L(j-(NN),8)
    row_N(j-(NN+1))=L(j-(NN+1),9)
    row_N(j-(NN+2))=L(j-(NN+2),10)
    row_N(j-2*NN)=L(j-2*NN,11)
    L(j,7)=1D0/L(j,1)*(L(j,7)-dot_product(row_M(2*nn+1:j-1),row_N(2*nn+1:j-1)))
elseif ((j+NN).eq.i) then

row_M=0

row_M(i)=L(i,1)
row_M(i-1)=L(i-1,2)
row_M(i-2)=L(i-2,3)
row_M(i-3)=L(i-3,4)
row_M(i-4)=L(i-4,5)
row_M(i-(NN-2))=L(i-(NN-2),6)
row_M(i-(NN-1))=L(i-(NN-1),7)
row_M(i-NN)=L(i-(NN),8)
row_M(i-(NN+1))=L(i-(NN+1),9)
row_M(i-(NN+2))=L(i-(NN+2),10)
row_M(i-2*NN)=L(i-2*NN,11)
        row_N=0
    row_N(j)=L(j,1)
    row_N(j-1)=L(j-1,2)
    row_N(j-2)=L(j-2,3)
    row_N(j-3)=L(j-3,4)
    row_N(j-4)=L(j-4,5)
    row_N(j-(NN-2))=L(j-(NN-2),6)
    row_N(j-(NN-1))=L(j-(NN-1),7)
    row_N(j-NN)=L(j-(NN),8)
    row_N(j-(NN+1))=L(j-(NN+1),9)
    row_N(j-(NN+2))=L(j-(NN+2),10)
    row_N(j-2*NN)=L(j-2*NN,11)
    L(j,8)=1D0/L(j,1)*(L(j,8)-dot_product(row_M(2*nn+1:j-1),row_N(2*nn+1:j-1)))
elseif ((j+NN+1).eq.i) then

row_M=0

row_M(i)=L(i,1)
row_M(i-1)=L(i-1,2)
row_M(i-2)=L(i-2,3)
row_M(i-3)=L(i-3,4)
row_M(i-4)=L(i-4,5)
row_M(i-(NN-2))=L(i-(NN-2),6)
row_M(i-(NN-1))=L(i-(NN-1),7)
row_M(i-NN)=L(i-(NN),8)
row_M(i-(NN+1))=L(i-(NN+1),9)
row_M(i-(NN+2))=L(i-(NN+2),10)
row_M(i-2*NN)=L(i-2*NN,11)
        row_N=0
    row_N(j)=L(j,1)
    row_N(j-1)=L(j-1,2)
    row_N(j-2)=L(j-2,3)
    row_N(j-3)=L(j-3,4)
    row_N(j-4)=L(j-4,5)
    row_N(j-(NN-2))=L(j-(NN-2),6)
    row_N(j-(NN-1))=L(j-(NN-1),7)
    row_N(j-NN)=L(j-(NN),8)
    row_N(j-(NN+1))=L(j-(NN+1),9)
    row_N(j-(NN+2))=L(j-(NN+2),10)
    row_N(j-2*NN)=L(j-2*NN,11)
    L(j,9)=1D0/L(j,1)*(L(j,9)-dot_product(row_M(2*nn+1:j-1),row_N(2*nn+1:j-1)))

elseif ((j+NN+2).eq.i) then

row_M=0

row_M(i)=L(i,1)
row_M(i-1)=L(i-1,2)
row_M(i-2)=L(i-2,3)
row_M(i-3)=L(i-3,4)
row_M(i-4)=L(i-4,5)
row_M(i-(NN-2))=L(i-(NN-2),6)
row_M(i-(NN-1))=L(i-(NN-1),7)
row_M(i-NN)=L(i-(NN),8)
row_M(i-(NN+1))=L(i-(NN+1),9)
row_M(i-(NN+2))=L(i-(NN+2),10)
row_M(i-2*NN)=L(i-2*NN,11)
        row_N=0
    row_N(j)=L(j,1)
    row_N(j-1)=L(j-1,2)
    row_N(j-2)=L(j-2,3)
    row_N(j-3)=L(j-3,4)
    row_N(j-4)=L(j-4,5)
    row_N(j-(NN-2))=L(j-(NN-2),6)
    row_N(j-(NN-1))=L(j-(NN-1),7)
    row_N(j-NN)=L(j-(NN),8)
    row_N(j-(NN+1))=L(j-(NN+1),9)
    row_N(j-(NN+2))=L(j-(NN+2),10)
    row_N(j-2*NN)=L(j-2*NN,11)
    L(j,10)=1D0/L(j,1)*(L(j,10)-dot_product(row_M(2*nn+1:j-1),row_N(2*nn+1:j-1)))
elseif ((j+2*NN).eq.i) then

row_M=0

row_M(i)=L(i,1)
row_M(i-1)=L(i-1,2)
row_M(i-2)=L(i-2,3)
row_M(i-3)=L(i-3,4)
row_M(i-4)=L(i-4,5)
row_M(i-(NN-2))=L(i-(NN-2),6)
row_M(i-(NN-1))=L(i-(NN-1),7)
row_M(i-NN)=L(i-(NN),8)
row_M(i-(NN+1))=L(i-(NN+1),9)
row_M(i-(NN+2))=L(i-(NN+2),10)
row_M(i-2*NN)=L(i-2*NN,11)
        row_N=0
    row_N(j)=L(j,1)
    row_N(j-1)=L(j-1,2)
    row_N(j-2)=L(j-2,3)
    row_N(j-3)=L(j-3,4)
    row_N(j-4)=L(j-4,5)
    row_N(j-(NN-2))=L(j-(NN-2),6)
    row_N(j-(NN-1))=L(j-(NN-1),7)
    row_N(j-NN)=L(j-(NN),8)
    row_N(j-(NN+1))=L(j-(NN+1),9)
    row_N(j-(NN+2))=L(j-(NN+2),10)
    row_N(j-2*NN)=L(j-2*NN,11)
    L(j,11)=1D0/L(j,1)*(L(j,11)-dot_product(row_M(2*nn+1:j-1),row_N(2*nn+1:j-1)))
end if


end do
enddo



do j=1,11
do i=1,nn*mm+2*nn
    if (isnan(L(i,j))) then
        alpha=alpha*1.15D0
        k=k+1
        goto 10
    end if
end do
enddo

!L=0
!L(:,1)=1

end subroutine precond

end module
