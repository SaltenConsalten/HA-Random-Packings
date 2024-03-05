module Algorithms
   USE global_variables
   USE omp_lib
   USE tools
   USE matrix
    !use liblapack
    !use librefblas
EXTERNAL         DPOTRF, DPOTRS
!    implicit none
real*8 :: eps=1.0D-8



    contains

    subroutine Cholesky(M, b, NN, MM,recalc)
!*************************************************************
!   Here, the Cholesky decomposition of the coefficient matrix is calculated.
!   Afterwards the SLE is solved.
!****************************************************************************
implicit none

    real*8, intent(inout) :: M(NN*MM,NN*MM)                         !coefficient matrix
    real*8, intent(inout) :: b(NN*MM)                        !right hand side(s)
    character(1) :: UPLO = 'U'      !store either the upper 'U' triangle of M or lower 'L'
    integer, INTENT(IN) :: NN
    integer, INTENT(IN) :: MM
    integer :: LDM, LDB                 !Leading dimensions of matrix M and right hand side b
    integer :: N                !order of M, in this case N=LDM
    integer :: NRHS           !Number of right hand sides, equals number of components -1
    integer :: err_inf                  !Error output of dpotrs. if 0 no error occurred
    logical, intent(inout) :: recalc

LDM = NN*MM
LDB = NN*MM
N = LDM
NRHS = 1

call openblas_set_num_threads(1)

        if(recalc) then

            call dpotrf(uplo, N, M, LDM, err_inf)

            recalc=.false.

            if(err_inf.eq.0) then

                call dpotrs(UPLO, N, NRHS, M, LDM, b, LDB, err_inf)

            end if

        else

            call dpotrs(UPLO, N, NRHS, M, LDM, b, LDB, err_inf)

        endif

    end subroutine Cholesky


    subroutine ConGrad(M,b,NN,MM,sol_old)
!****************************************************************************************
!   Here, SLE is solved by the Conjugate Gradient algorithm
!****************************************************************************************
implicit none

integer, intent(in) :: NN, MM
real*8, intent (inout) :: M(:,:)
real*8, intent (inout) :: b(NN*MM)
integer :: i,z,j
real*8, intent(in), optional :: sol_old(NN*MM)
real*8 :: A(NN*MM+2*NN,11)
real*8 :: x(NN*MM+4*NN), x_old(NN*MM+4*NN)
real*8 :: y(NN*MM+4*NN)
real*8 :: r(NN*MM+4*NN)
real*8 :: p(NN*MM+4*NN)
real*8 :: c(21),row_M(21)
real*8 :: eps, gam, alpha, beta, gam2


if(present(sol_old))then
    x(2*NN+1:NN*MM+2*NN)=sol_old
else
    x=1
end if

A=0
y=0
r=0
p=0
gam2=1D0




A(2*NN+1:NN*MM+2*NN,:)=M

do i=2*NN+1,NN*MM+2*NN
row_M=0
c=0
row_M(1)=A(i-2*NN,11)
row_M(2)=A(i-(NN+2),10)
row_M(3)=A(i-(NN+1),9)
row_M(4)=A(i-(NN),8)
row_M(5)=A(i-(NN-1),7)
row_M(6)=A(i-(NN-2),6)
row_M(7)=A(i-4,5)
row_M(8)=A(i-3,4)
row_M(9)=A(i-2,3)
row_M(10)=A(i-1,2)
row_M(11)=A(i,1)
row_M(12)=A(i,2)
row_M(13)=A(i,3)
row_M(14)=A(i,4)
row_M(15)=A(i,5)
row_M(16)=A(i,6)
row_M(17)=A(i,7)
row_M(18)=A(i,8)
row_M(19)=A(i,9)
row_M(20)=A(i,10)
row_M(21)=A(i,11)

c(1)=x(i-2*NN)
c(2)=x(i-(NN+2))
c(3)=x(i-(NN+1))
c(4)=x(i-(NN))
c(5)=x(i-(NN-1))
c(6)=x(i-(NN-2))
c(7)=x(i-4)
c(8)=x(i-3)
c(9)=x(i-2)
c(10)=x(i-1)
c(11)=x(i)
c(12)=x(i+1)
c(13)=x(i+2)
c(14)=x(i+3)
c(15)=x(i+4)
c(16)=x(i+(NN-2))
c(17)=x(i+(NN-1))
c(18)=x(i+NN)
c(19)=x(i+(NN+1))
c(20)=x(i+(NN+2))
c(21)=x(i+(NN*2))

z=i-2*NN
r(i)=b(z)-dot_product(row_M,c)
end do
p=r
gam=sum(r**2)
gam=norm2(r)**2
eps=1E-24
j=0
do while (gam2.gt.eps)
!    call omp_set_num_threads(3) !firstprivate(A) x,y,p,c,row_M,
! !$omp parallel shared(A,p,y) private(c,row_M,i)
! !$omp do
    do i=2*NN+1,NN*MM+2*NN
row_M=0
c=0
row_M(1)=A(i-2*NN,11)
row_M(2)=A(i-(NN+2),10)
row_M(3)=A(i-(NN+1),9)
row_M(4)=A(i-(NN),8)
row_M(5)=A(i-(NN-1),7)
row_M(6)=A(i-(NN-2),6)
row_M(7)=A(i-4,5)
row_M(8)=A(i-3,4)
row_M(9)=A(i-2,3)
row_M(10)=A(i-1,2)
row_M(11)=A(i,1)
row_M(12)=A(i,2)
row_M(13)=A(i,3)
row_M(14)=A(i,4)
row_M(15)=A(i,5)
row_M(16)=A(i,6)
row_M(17)=A(i,7)
row_M(18)=A(i,8)
row_M(19)=A(i,9)
row_M(20)=A(i,10)
row_M(21)=A(i,11)

c(1)=p(i-2*NN)
c(2)=p(i-(NN+2))
c(3)=p(i-(NN+1))
c(4)=p(i-(NN))
c(5)=p(i-(NN-1))
c(6)=p(i-(NN-2))
c(7)=p(i-4)
c(8)=p(i-3)
c(9)=p(i-2)
c(10)=p(i-1)
c(11)=p(i)
c(12)=p(i+1)
c(13)=p(i+2)
c(14)=p(i+3)
c(15)=p(i+4)
c(16)=p(i+(NN-2))
c(17)=p(i+(NN-1))
c(18)=p(i+NN)
c(19)=p(i+(NN+1))
c(20)=p(i+(NN+2))
c(21)=p(i+(NN*2))

            y(i)=dot_product(row_M,c)
    end do
! !$omp end do
! !$omp end parallel
        alpha=gam/dot_product(y(2*NN+1:NN*MM+2*NN),p(2*NN+1:NN*MM+2*NN))

    gam2=dsqrt(sum((x_old(2*NN+1:NN*MM+2*NN)-x(2*NN+1:NN*MM+2*NN))**2)/size(x(2*NN+1:NN*MM+2*NN)))/(sum(abs(x(2*NN+1:NN*MM+2*NN)))&
                            /size(x(2*NN+1:NN*MM+2*NN)))

        x_old(2*NN+1:NN*MM+2*NN)=x(2*NN+1:NN*MM+2*NN)


        x(2*NN+1:NN*MM+2*NN)=x(2*NN+1:NN*MM+2*NN)+alpha*p(2*NN+1:NN*MM+2*NN)
        r(2*NN+1:NN*MM+2*NN)=r(2*NN+1:NN*MM+2*NN)-alpha*y(2*NN+1:NN*MM+2*NN)
        beta=norm2(r)**2/gam
        gam=norm2(r)**2
        p(2*NN+1:NN*MM+2*NN)=r(2*NN+1:NN*MM+2*NN)+beta*p(2*NN+1:NN*MM+2*NN)


        if (j>1D6) then
            write(*,*) 'CG iteration exceeds maximum, subroutine will stop'
            !call sleep(10)
            return
        endif
        j=j+1
    enddo
b=x(2*NN+1:NN*MM+2*NN)

    end subroutine ConGrad

subroutine ICCG(M,L,b,NN,MM,sol_old)
!****************************************************************************************
!   Here, SLE is solved by the Conjugate Gradient algorithm
!****************************************************************************************
implicit none

integer, intent(in) :: NN, MM
real*8, intent (in) :: L(NN*MM+2*NN,11)
real*8, intent (inout) :: M(NN*MM,11)
real*8, intent (inout) :: b(NN*MM)
integer :: i,k,j!, w
real*8, intent(in), optional :: sol_old(NN*MM)
real*8 :: A(NN*MM+2*NN,11)
real*8 :: x(NN*MM+4*NN)
real*8 :: y(NN*MM+4*NN)
real*8 :: r(NN*MM+4*NN), x_old(NN*MM+4*NN)
real*8 :: z(NN*MM+4*NN)
real*8 :: v(NN*MM+4*NN)
real*8 :: p(NN*MM+4*NN)
real*8 :: c(21),row_M(21)
real*8 :: gam, alpha, beta, gam2
!######################################################################
!
!GLS Initialisieren
!
!######################################################################
gam2=1
x_old=0
x_old(2*NN+1:NN*MM+2*NN)=2
x=0
x(2*NN+1:NN*MM+2*NN)=1
if(present(sol_old))then
    x(2*NN+1:NN*MM+2*NN)=sol_old
end if

A=0
y=0
z=0
v=0
r=0
p=0


A(2*NN+1:NN*MM+2*NN,:)=M

    do i=2*NN+1,NN*MM+2*NN
    row_M=0
    c=0
    row_M(1)=A(i-2*NN,11)
    row_M(2)=A(i-(NN+2),10)
    row_M(3)=A(i-(NN+1),9)
    row_M(4)=A(i-(NN),8)
    row_M(5)=A(i-(NN-1),7)
    row_M(6)=A(i-(NN-2),6)
    row_M(7)=A(i-4,5)
    row_M(8)=A(i-3,4)
    row_M(9)=A(i-2,3)
    row_M(10)=A(i-1,2)
    row_M(11)=A(i,1)
    row_M(12)=A(i,2)
    row_M(13)=A(i,3)
    row_M(14)=A(i,4)
    row_M(15)=A(i,5)
    row_M(16)=A(i,6)
    row_M(17)=A(i,7)
    row_M(18)=A(i,8)
    row_M(19)=A(i,9)
    row_M(20)=A(i,10)
    row_M(21)=A(i,11)

    c(1)=x(i-2*NN)
    c(2)=x(i-(NN+2))
    c(3)=x(i-(NN+1))
    c(4)=x(i-(NN))
    c(5)=x(i-(NN-1))
    c(6)=x(i-(NN-2))
    c(7)=x(i-4)
    c(8)=x(i-3)
    c(9)=x(i-2)
    c(10)=x(i-1)
    c(11)=x(i)
    c(12)=x(i+1)
    c(13)=x(i+2)
    c(14)=x(i+3)
    c(15)=x(i+4)
    c(16)=x(i+(NN-2))
    c(17)=x(i+(NN-1))
    c(18)=x(i+NN)
    c(19)=x(i+(NN+1))
    c(20)=x(i+(NN+2))
    c(21)=x(i+(NN*2))

    k=i-2*NN
    r(i)=b(k)-dot_product(row_M,c)
    end do
!#####################################################################################
!	z=r
!   z initalisieren, Berechnung der Inversen durch Vorwärts-, Rückwärtseinsetzen
!#####################################################################################
    do i=2*NN+1,NN*MM+2*NN
        row_M=0
        c=0

        row_M(1)=L(i-2*NN,11)
        row_M(2)=L(i-(NN+2),10)
        row_M(3)=L(i-(NN+1),9)
        row_M(4)=L(i-(NN),8)
        row_M(5)=L(i-(NN-1),7)
        row_M(6)=L(i-(NN-2),6)
        row_M(7)=L(i-4,5)
        row_M(8)=L(i-3,4)
        row_M(9)=L(i-2,3)
        row_M(10)=L(i-1,2)
        row_M(11)=L(i,1)

        c(1)=v(i-2*NN)
        c(2)=v(i-(NN+2))
        c(3)=v(i-(NN+1))
        c(4)=v(i-(NN))
        c(5)=v(i-(NN-1))
        c(6)=v(i-(NN-2))
        c(7)=v(i-4)
        c(8)=v(i-3)
        c(9)=v(i-2)
        c(10)=v(i-1)
        c(11)=v(i)

        if ((i-NN*2).eq.1) then
            v(i)=r(i)/L(i,1)
        else
            v(i)=(1/L(i,1))*(r(i)-dot_product(row_M(1:10),c(1:10)))
        end if
        if (isnan(v(i))) then
            call sleep(1)
        end if

    end do

    do i=NN*MM+2*NN,2*NN+1,-1
        row_M=0
        c=0

        row_M(11)=L(i,1)
        row_M(12)=L(i,2)
        row_M(13)=L(i,3)
        row_M(14)=L(i,4)
        row_M(15)=L(i,5)
        row_M(16)=L(i,6)
        row_M(17)=L(i,7)
        row_M(18)=L(i,8)
        row_M(19)=L(i,9)
        row_M(20)=L(i,10)
        row_M(21)=L(i,11)

        c(11)=z(i)
        c(12)=z(i+1)
        c(13)=z(i+2)
        c(14)=z(i+3)
        c(15)=z(i+4)
        c(16)=z(i+(NN-2))
        c(17)=z(i+(NN-1))
        c(18)=z(i+NN)
        c(19)=z(i+(NN+1))
        c(20)=z(i+(NN+2))
        c(21)=z(i+(NN*2))

        if (i.eq.nn*mm+2*nn) then
                z(i)=v(i)/L(i,1)
            else
                z(i)=(1/L(i,1))*(v(i)-dot_product(row_M(12:21),c(12:21)))
            end if

    end do

    p=z
    gam=norm2(r)**2
    j=0

!#####################################################################################
!
!       Start CG Verfahren
!
!#####################################################################################

    do while ((gam2.gt.eps).and.(j.lt.num_it))!(gam.gt.1E-200)
    !$omp parallel
    !$omp do
        do i=2*NN+1,NN*MM+2*NN
            row_M=0
            c=0
            row_M(1)=A(i-2*NN,11)
            row_M(2)=A(i-(NN+2),10)
            row_M(3)=A(i-(NN+1),9)
            row_M(4)=A(i-(NN),8)
            row_M(5)=A(i-(NN-1),7)
            row_M(6)=A(i-(NN-2),6)
            row_M(7)=A(i-4,5)
            row_M(8)=A(i-3,4)
            row_M(9)=A(i-2,3)
            row_M(10)=A(i-1,2)
            row_M(11)=A(i,1)
            row_M(12)=A(i,2)
            row_M(13)=A(i,3)
            row_M(14)=A(i,4)
            row_M(15)=A(i,5)
            row_M(16)=A(i,6)
            row_M(17)=A(i,7)
            row_M(18)=A(i,8)
            row_M(19)=A(i,9)
            row_M(20)=A(i,10)
            row_M(21)=A(i,11)

            c(1)=p(i-2*NN)
            c(2)=p(i-(NN+2))
            c(3)=p(i-(NN+1))
            c(4)=p(i-(NN))
            c(5)=p(i-(NN-1))
            c(6)=p(i-(NN-2))
            c(7)=p(i-4)
            c(8)=p(i-3)
            c(9)=p(i-2)
            c(10)=p(i-1)
            c(11)=p(i)
            c(12)=p(i+1)
            c(13)=p(i+2)
            c(14)=p(i+3)
            c(15)=p(i+4)
            c(16)=p(i+(NN-2))
            c(17)=p(i+(NN-1))
            c(18)=p(i+NN)
            c(19)=p(i+(NN+1))
            c(20)=p(i+(NN+2))
            c(21)=p(i+(NN*2))

            y(i)=dot_product(row_M,c)
    end do
 !$omp end do
 !$omp end parallel
    gam=dot_product(z(2*NN+1:NN*MM+2*NN),r(2*NN+1:NN*MM+2*NN))

    alpha=gam/dot_product(y(2*NN+1:NN*MM+2*NN),p(2*NN+1:NN*MM+2*NN))
    gam2=dsqrt(sum((x_old(2*NN+1:NN*MM+2*NN)-x(2*NN+1:NN*MM+2*NN))**2)/size(x(2*NN+1:NN*MM+2*NN)))/(sum(abs(x(2*NN+1:NN*MM+2*NN)))&
                            /size(x(2*NN+1:NN*MM+2*NN)))

    x_old(2*NN+1:NN*MM+2*NN)=x(2*NN+1:NN*MM+2*NN)

    x(2*NN+1:NN*MM+2*NN)=x(2*NN+1:NN*MM+2*NN)+alpha*p(2*NN+1:NN*MM+2*NN)

    r(2*NN+1:NN*MM+2*NN)=r(2*NN+1:NN*MM+2*NN)-alpha*y(2*NN+1:NN*MM+2*NN)

!#######################################################################################
! 			z_k=C⁻¹*r_k
! 	lösen durch vorwärts-, rückwärtseinsetzen
!			Lv=r, v=Rz
!#######################################################################################

     do i=2*NN+1,NN*MM+2*NN
        row_M=0
        c=0

        row_M(1)=L(i-2*NN,11)
        row_M(2)=L(i-(NN+2),10)
        row_M(3)=L(i-(NN+1),9)
        row_M(4)=L(i-(NN),8)
        row_M(5)=L(i-(NN-1),7)
        row_M(6)=L(i-(NN-2),6)
        row_M(7)=L(i-4,5)
        row_M(8)=L(i-3,4)
        row_M(9)=L(i-2,3)
        row_M(10)=L(i-1,2)
        row_M(11)=L(i,1)

        c(1)=v(i-2*NN)
        c(2)=v(i-(NN+2))
        c(3)=v(i-(NN+1))
        c(4)=v(i-(NN))
        c(5)=v(i-(NN-1))
        c(6)=v(i-(NN-2))
        c(7)=v(i-4)
        c(8)=v(i-3)
        c(9)=v(i-2)
        c(10)=v(i-1)
        c(11)=v(i)

        if ((i-NN*2).eq.1) then
            v(i)=r(i)/L(i,1)
        else
            v(i)=(1/L(i,1))*(r(i)-dot_product(row_M(1:10),c(1:10)))
        end if
    end do

    do i=NN*MM+2*NN,2*NN+1,-1
        row_M=0
        c=0

        row_M(11)=L(i,1)
        row_M(12)=L(i,2)
        row_M(13)=L(i,3)
        row_M(14)=L(i,4)
        row_M(15)=L(i,5)
        row_M(16)=L(i,6)
        row_M(17)=L(i,7)
        row_M(18)=L(i,8)
        row_M(19)=L(i,9)
        row_M(20)=L(i,10)
        row_M(21)=L(i,11)

        c(11)=z(i)
        c(12)=z(i+1)
        c(13)=z(i+2)
        c(14)=z(i+3)
        c(15)=z(i+4)
        c(16)=z(i+(NN-2))
        c(17)=z(i+(NN-1))
        c(18)=z(i+NN)
        c(19)=z(i+(NN+1))
        c(20)=z(i+(NN+2))
        c(21)=z(i+(NN*2))

        if (i.eq.nn*mm+2*nn) then
                z(i)=v(i)/L(i,1)
            else
                z(i)=(1/L(i,1))*(v(i)-dot_product(row_M(12:21),c(12:21)))
            end if
    end do

    beta=dot_product(z(2*NN+1:NN*MM+2*NN),r(2*NN+1:NN*MM+2*NN))/gam
    p(2*NN+1:NN*MM+2*NN)=z(2*NN+1:NN*MM+2*NN)+beta*p(2*NN+1:NN*MM+2*NN)
!   eps=dsqrt(sum((x_old-x(2*NN+1:NN*MM+2*NN))**2)/size(x_old))/(sum(x(2*NN+1:NN*MM+2*NN))/size(x_old))!

        if (j>1D6) then
            write(*,*) 'ICCG iteration exceeds maximum, subroutine will stop'
            !call sleep(10)
            return
        endif
        j=j+1
    enddo

write(*,*) 'ICCG iteration: ', j
b=x(2*NN+1:NN*MM+2*NN)


end subroutine ICCG

    subroutine solve(Nzl,Nzg,j,R_c,delta,Sc_l,turb_ga,D_g,D_l,Dist_coef,k_reac,dx,dh,eta_gas,eta_liq,ul,ug,x_0,conc_l,y_0,conc_g,&
                    trans_old,trans,val_vec,meanmat,row,col_ind,M,L, n_A,h_R,tcalc,Mnew)
        implicit none

        integer, intent(in) :: Nzl,Nzg,j
        real*8, intent(in) :: R_c, delta, Sc_l, D_g, D_l, dx(2),dh(2)
        real*8, intent(in) :: n_A(nzl*mm), h_r(nn*mm,nzl)
        real*8, intent(in) :: conc_l,conc_g,Dist_coef,k_reac(3),x_0,y_0
        real*8, intent(in) :: turb_ga(nng+1),eta_gas(nng), eta_liq(nnl), ul(nnl), ug(nng), trans_old(nn*mm,Nzl)
        logical, intent(in):: tcalc
        real*8, intent(inout) :: trans(nn*mm,Nzl)
        integer :: i
!        integer :: thread
!        real*8, allocatable :: IC_l(:),IC_g(:), dissflux(:,:)
        real*8 :: IC_l(Nzl), IC_g(Nzl)!, dissflux(Nzl,2)
        real*8, intent(inout) :: M(:,:),L(:,:)
!        real*8, allocatable :: M(:,:), L(:,:)
        real*8, allocatable :: b(:)
        logical :: recalc=.false., Mnew
        real*8, intent(inout) :: val_vec(5*NN*MM-2*(NN+(MM-1))), meanmat(3)
        real*8, allocatable :: c0g(:),c0l(:), r_l(:),r_g(:)!c0g(NNg),c0l(NNl), r_l(nnl),r_g(nng)
        integer, intent(inout) :: row(5*NN*MM-2*(NN+(MM-1))), col_ind(NN*MM+1)
!        integer, allocatable :: row(:), col_ind(:)
!        real*8, allocatable :: val_vec(:), meanmat(:)
!        real*8 :: mvt

        logical :: precondition=.true.
!        logical :: restorediff=.false.


!        allocate(IC_g(Nzl),IC_l(Nzl), dissflux(Nzl,2))
        allocate(c0g(NNg), c0l(NNl), r_l(nnl), r_g(nng))





       if(eta_liq(nnl).gt.eta_liq(1)) then
            r_l=eta_liq*delta
            r_g=eta_gas*(R_c-delta)+delta

        else
            r_l=eta_liq*delta+(R_c-delta)
            r_g=eta_gas*(R_c-delta)

        end if


        if (.not.threshold) then
            do i=1,Nzl
                if(maxval(trans_old).eq.0)then                          !field initialization

                    if(i.eq.1) then
                        C0l(1:NNl)=x_0
!                        C0g(1:NNg)=(dble(i)/dble(Nzl))*y_0
                        if(tcalc)then
                            C0g(1:NNg)=y_0-(dble(Nzl-i)/dble(Nzl))*(conc_g-conc_l)/conc_g
                        else
                            C0g(1:NNg)=dble(i)/dble(Nzl)*y_0
                        endif
                        call FullMatrix(NNg,NNl,MM,R_c,delta,ug,turb_ga(:),D_g,k_reac(:),dx,dh,eta_gas,c0g,ul,&
                                        Sc_l,D_l,eta_liq,c0l,Dist_coef,conc_l,conc_g,b,M,recalc,meanmat,val_vec,row,col_ind,&
                                        n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
                        call Cholesky(M,b,NN,MM,recalc)
                    elseif(i.eq.Nzl) then
                        c0g(1:NNg)=y_0
                        C0l(1:NNl)=trans(NN*(MM-1)+1:NN*(MM-1)+NNl,i-1)
                        call mean_int(c0l,c0g,ul,ug,r_l,r_g,IC_l(i),IC_g(i))
                        c0l=IC_l(i)
                        call right_hand_side(NNg, NNl, MM, Ug,dx,c0g,Ul,c0l, b,meanmat, val_vec, row, col_ind,&
                                            n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
                        call Cholesky(M,b,NN,MM,recalc)
                    else
!                        C0g(1:NNg)=(dble(i)/dble(Nzl))*y_0
                        if(tcalc)then
                            C0g(1:NNg)=y_0-(dble(Nzl-i)/dble(Nzl))*(conc_g-conc_l)/conc_g
                        else
                            C0g(1:NNg)=dble(i)/dble(Nzl)*y_0
                        endif
                        C0l(1:NNl)=trans(NN*(MM-1)+1:NN*(MM-1)+NNl,i-1)
                        call mean_int(c0l,c0g,ul,ug,r_l,r_g,IC_l(i),IC_g(i))
                        c0l=IC_l(i)
                        call right_hand_side(NNg, NNl, MM, Ug,dx,c0g,Ul,c0l, b,meanmat, val_vec, row, col_ind,&
                                            n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
                        call Cholesky(M,b,NN,MM,recalc)

                    end if


                else                                !Iterations

                    if(i.eq.1) then
                        C0l(1:NNl)=x_0

                        c0l=c0l
                        C0g(1:NNg)=trans_old(NNL+1:NN,i+1)
                        if (j.eq.1) then
                            call FullMatrix(NNg,NNl,MM,R_c,delta,ug,turb_ga(:),D_g,k_reac(:),dx,dh,eta_gas,c0g,ul,&
                                        Sc_l,D_l,eta_liq,c0l,Dist_coef,conc_l,conc_g,b,M,recalc,meanmat, val_vec,row,col_ind,&
                                        n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
                        elseif (Mnew) then
                            call FullMatrix(NNg,NNl,MM,R_c,delta,ug,turb_ga(:),D_g,k_reac(:),dx,dh,eta_gas,c0g,ul,&
                                        Sc_l,D_l,eta_liq,c0l,Dist_coef,conc_l,conc_g,b,M,recalc,meanmat,val_vec,row,col_ind,&
                                        n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
!                            Mnew=.false.
                        end if
                        call right_hand_side(NNg, NNl, MM, Ug,dx,c0g,Ul,c0l, b,meanmat, val_vec, row, col_ind,&
                                            n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
                        call Cholesky(M,b,NN,MM,recalc)


                    elseif(i.eq.Nzl) then
                        c0g(1:NNg)=y_0
                        C0l(1:NNl)=trans_old((MM-1)*NN+1:(MM-1)*NN+NNl,i-1)
                        call right_hand_side(NNg, NNl, MM, Ug,dx,c0g,Ul,c0l, b,meanmat, val_vec, row, col_ind,&
                                            n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
                        call Cholesky(M,b,NN,MM,recalc)

                    else

                        c0l=trans_old((MM-1)*NN+1:(MM-1)*NN+NNl,i-1)
                        c0g=trans_old(NNL+1:NN,i+1)
                        call right_hand_side(NNg, NNl, MM, Ug,dx,c0g,Ul,c0l, b,meanmat, val_vec, row, col_ind,&
                                            n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
                        call Cholesky(M,b,NN,MM,recalc)
                    endif


                endif

                      trans(:,i)=b
!                      if(minval(b).lt.0) then
!                        b=b
!                      end if

            enddo

        else
!#################################################################
!code for conjugate gradient method
!#################################################################

            do i=1,nzl
                if(maxval(trans_old).eq.0)then                          !field initialization

                    if(i.eq.1) then
                        C0l(1:NNl)=x_0
!                        C0g(1:NNg)=dble(i)/dble(Nzl)*y_0
                        if(tcalc)then
                            C0g(1:NNg)=y_0-(dble(Nzl-i)/dble(Nzl))*(conc_g-conc_l)/conc_g
                        else
                            C0g(1:NNg)=dble(i)/dble(Nzl)*y_0
                        endif
                        call CompMatrix(NNg,NNl,MM,R_c,delta,ug,turb_ga(:),D_g,k_reac(:),dx,dh,eta_gas,c0g,ul,Sc_l,D_l,&
                                eta_liq,c0l,Dist_coef,conc_l,conc_g,b,M,meanmat, val_vec, row, col_ind,&
                                n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
                        if (precondition) then
                            call precond(NN,MM,M,L)
                            call ICCG(M,L,b,NN,MM)
                        else
                            call ConGrad(M,b,NN,MM)
                        endif

                    elseif(i.eq.Nzl) then
                        c0g(1:NNg)=y_0
                        C0l(1:NNl)=trans(NN*(MM-1)+1:NN*(MM-1)+NNl,i-1)
                        call mean_int(c0l,c0g,ul,ug,r_l,r_g,IC_l(i),IC_g(i))
                        c0l=IC_l(i)
                        call right_hand_side(NNg, NNl, MM, Ug,dx,c0g,Ul,c0l, b,meanmat, val_vec, row, col_ind,&
                                            n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
                        if (precondition) then
                            call ICCG(M,L,b,NN,MM)
                        else
                            call ConGrad(M,b,NN,MM)
                        endif
                    else
!                        C0g(1:NNg)=dble(i)/dble(Nzl)*y_0
                        if(tcalc)then
                            C0g(1:NNg)=y_0-(dble(Nzl-i)/dble(Nzl))*(conc_g-conc_l)/conc_g
                        else
                            C0g(1:NNg)=dble(i)/dble(Nzl)*y_0
                        endif
                        C0l(1:NNl)=trans(NN*(MM-1)+1:NN*(MM-1)+NNl,i-1)
                        call mean_int(c0l,c0g,ul,ug,r_l,r_g,IC_l(i),IC_g(i))
                        c0l=IC_l(i)
                        call right_hand_side(NNg, NNl, MM, Ug,dx,c0g,Ul,c0l, b,meanmat, val_vec, row, col_ind,&
                                            n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
                        if (precondition) then
                            call ICCG(M,L,b,NN,MM)
                        else
                            call ConGrad(M,b,NN,MM)
                        endif

                    end if


                else                                !Iterations

                    if(i.eq.1) then
                        C0l(1:NNl)=x_0

                        c0l=c0l
                        C0g(1:NNg)=trans_old(NNL+1:NN,i+1)
                        if ((j.eq.1).or.Mnew) then
                            call CompMatrix(NNg,NNl,MM,R_c,delta,ug,turb_ga(:),D_g,k_reac(:),dx,dh,eta_gas,c0g,ul,Sc_l,D_l,&
                                            eta_liq,c0l,Dist_coef,conc_l,conc_g,b,M,meanmat, val_vec, row, col_ind,&
                                            n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
                            if (precondition) then
                                call precond(NN,MM,M,L)
                            end if
                        end if
                        call right_hand_side(NNg, NNl, MM, Ug,dx,c0g,Ul,c0l, b,meanmat, val_vec, row, col_ind,&
                                            n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
                        if (precondition) then
                            call ICCG(M,L,b,NN,MM,trans(:,i))
                        else
                            call ConGrad(M,b,NN,MM,trans(:,i))
                        endif


                    elseif(i.eq.Nzl) then
                        c0g(1:NNg)=y_0
                        C0l(1:NNl)=trans_old((MM-1)*NN+1:(MM-1)*NN+NNl,i-1)
                        c0g(1:NNg)=c0g(1:NNg)
                        call right_hand_side(NNg, NNl, MM, Ug,dx,c0g,Ul,c0l, b,meanmat, val_vec, row, col_ind,&
                                            n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
                        if (precondition) then
                            call ICCG(M,L,b,NN,MM,trans(:,i))
                        else
                            call ConGrad(M,b,NN,MM,trans(:,i))
                        endif

                    else

                        c0l=trans_old((MM-1)*NN+1:(MM-1)*NN+NNl,i-1)
                        c0g=trans_old(NNL+1:NN,i+1)
                        call right_hand_side(NNg, NNl, MM, Ug,dx,c0g,Ul,c0l, b,meanmat, val_vec, row, col_ind,&
                                            n_A((i-1)*mm+1:(i-1)*mm+MM),h_r(:,i),tcalc)
                        if (precondition) then
                            call ICCG(M,L,b,NN,MM,trans(:,i))
                        else
                            call ConGrad(M,b,NN,MM,trans(:,i))
                        endif
                    endif


                endif

                    trans(:,i)=b

            enddo

        endif


    end subroutine

end module
