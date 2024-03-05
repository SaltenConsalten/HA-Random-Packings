module global_variables
    implicit none
    integer, public :: ncomp
    integer, public :: NNl, NNg, NN, MM, nzl
    logical :: threshold=.false.
    real*8 :: ae                                ! wetted area
    real*8 :: g=9.81                            ! gravity acceleration
    real*8, PARAMETER :: pi = 4D0*datan(1D0)    ! pi
    real*8 :: delta_film, delta_drop, delta_jet                       ! liquid film width
    real*8 :: H!, dH                           ! packing height, discrete height
    real*8 :: R_c           !channel radius
    character*14, allocatable :: comp_order(:)
    integer :: num_it
    integer :: itmax=1000, NMmax=1.5D4
    logical, dimension(3) :: tcalc=(/.false.,.true.,.true./)              !(/.false.,.true.,.true./)                   ! when true, then temp. calc. otherwise conc. calc.
    real*8 :: B_
!    logical, dimension(2) :: Mnew=(/.false.,.false./)
end module
