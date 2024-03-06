
bool use_stiff = .True.;

void mini_ch_dlsode(double T_in, double P_in, double t_end, double * VMR, char * network){

    int  ncall;
    double P_cgs;

    // Time controls
    double t_begin, t_now, t_old;

    // DLSODE variables
    double  y[n_sp], y_old[n_sp];
    double *rwork, *rtol, *atol;
    integer, allocatable, dimension(:) :: iwork
    integer :: itol, itask, istate, iopt, mf
    integer :: rworkdim, iworkdim

    !! Find the number density of the atmosphere
    P_cgs = P_in * 10.0_dp   ! Convert pascal to dyne cm-2
    nd_atm = P_cgs/(kb*T_in)  ! Find initial number density [cm-3] of atmosphere

    allocate(Keq(n_reac), re_f(n_reac), re_r(n_reac))

    ! First find the reverse reaction coefficents (Keq)
    call reverse_reactions(T_in)
    ! Find the forward, backward and net reaction rates
    call reaction_rates(T_in, P_cgs, nd_atm)

    !! Find initial number density of all species from VMR
    y(:) = nd_atm * VMR(:)

    ! -----------------------------------------
    ! ***  parameters for the DLSODE solver  ***
    ! -----------------------------------------

    itask = 1
    istate = 1
    iopt = 1

    ! Method flag
    if (use_stiff .eqv. .True.) then
      ! Problem is stiff (usual)
      ! mf = 21 - full jacobian matrix with jacobian save
      ! mf = 22 - internal calculated jacobian
      mf = 21
      rworkdim = 22 +  9*n_sp + n_sp**2
      iworkdim = 20 + n_sp
      allocate(rtol(n_sp), atol(n_sp), rwork(rworkdim), iwork(iworkdim))

      itol = 4
      rtol(:) = 1.0e-3_dp           ! Relative tolerances for each scalar
      atol(:) = 1.0e-99_dp               ! Absolute tolerance for each scalar (floor value)

      rwork(:) = 0.0_dp
      iwork(:) = 0

      rwork(1) = 0.0_dp               ! Critical T value (don't integrate past time here)
      rwork(5) = 0.0_dp              ! Initial starting timestep (start low, will adapt in DVODE)
      rwork(6) = 0.0_dp       ! Maximum timestep

      iwork(5) = 0               ! Max order required
      iwork(6) = 100000               ! Max number of internal steps
      iwork(7) = 1                ! Number of error messages

    else
      ! Problem is not too stiff (not typical)
      ! mf = 11 - full jacobian matrix with jacobian save
      mf = 11
      rworkdim = 22 + 16*n_sp + 2*n_sp**2
      iworkdim = 30 + n_sp
      allocate(rtol(n_sp), atol(n_sp), rwork(rworkdim), iwork(iworkdim))
      itol = 4
      rtol(:) = 1.0e-3_dp
      atol(:) = 1.0e-99_dp

      rwork(1) = t_end

      rwork(5:10) = 0.0_dp
      iwork(5:10) = 0

      rwork(5) = 1.0e-99_dp
      rwork(6) = t_end
      iwork(6) = 100000
    end if

    t_begin = 0.0_dp
    t_now = t_begin

    ! Set the printing flag
    ! 0 = no printing, 1 = printing
    call xsetf(1)

    ncall = 0

    do while (t_now < t_end)

      y_old(:) = y(:)
      t_old = t_now

      select case(network)
      case('HO')
        call DLSODE (RHS_update, n_sp, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_HO, mf)
      case('CHO')
        call DLSODE (RHS_update, n_sp, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_CHO, mf)
      case('NCHO')
        call DLSODE (RHS_update, n_sp, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_NCHO, mf)
      case default
        print*, 'Invalid network provided: ', trim(network)
        stop
      end select

      ! call check_con(n_sp,y(:),y_old(:),t_now,t_old,con)
      ! if (con .eqv. .True.) then
      !   exit
      ! end if

      ncall = ncall + 1

      if (mod(ncall,50) == 0) then
        istate = 1
      else  if (istate == -1) then
        istate = 2
      else if (istate < -1) then
        print*, 'dlsode: ', istate
        exit
      end if

    end do

    VMR(:) = y(:)/nd_atm

    deallocate(Keq, re_r, re_f, rtol, atol, rwork, iwork)

  end subroutine mini_ch_dlsode

  subroutine RHS_update(NEQ, time, y, f, rpar, ipar)
    implicit none

    integer, intent(in) ::  NEQ
    real(dp), intent(inout) :: time
    real(dp), dimension(NEQ), intent(inout) :: y
    real(dp), dimension(NEQ), intent(inout) :: f
    real(dp), intent(inout) :: rpar
    integer, intent(inout) :: ipar

    integer :: i, k, j
    real(dp) :: msum, msum2, frate, rrate
    real(dp), dimension(n_reac) :: net_pr, net_re
    real(dp), dimension(NEQ) :: f_pr, f_re, t_pr, t_re
    real(dp), dimension(NEQ) :: c_pr, c_re

    ! Calculate the rate of change of number density for all species [cm-3/s] this is the f vector
    f_pr(:) = 0.0_dp
    c_pr(:) = 0.0_dp
    f_re(:) = 0.0_dp
    c_re(:) = 0.0_dp

    ! Loop through reactions add rates to the f array
    do i = 1, n_reac
      ! Do the forward and backward flux calculation for each speices in the reaction

      ! Find number density multiple for reactants in reaction
      msum = y(re(i)%gi_re(1))
      do k = 2, re(i)%n_re
         msum = msum * y(re(i)%gi_re(k))
      end do

      ! Find number density multiple for products in reaction
      msum2 = y(re(i)%gi_pr(1))
      do k = 2, re(i)%n_pr
         msum2 = msum2 * y(re(i)%gi_pr(k))
      end do

      if (re(i)%re_t == 3) then
        ! Mutliply both msum and msum2 by atmosphere nd if neutral body involved
        msum = msum * nd_atm
        msum2 = msum2 * nd_atm
      end if

      frate = msum * re_f(i)
      rrate = msum2 * re_r(i)

      net_pr(i) = frate - rrate
      net_re(i) = -net_pr(i)

      !! Perform the Kahan-Babushka-Neumaier compensation sum algorithm
      ! This is slightly slower than peicewise addition for small timesteps, but faster for larger timesteps, 
      ! and more general (should work for all networks)

      !! Add the product rates
      do j = 1, re(i)%n_pr
        t_pr(re(i)%gi_pr(j)) = f_pr(re(i)%gi_pr(j)) + net_pr(i)
        if (abs(f_pr(re(i)%gi_pr(j))) >= abs(net_pr(i))) then
          c_pr(re(i)%gi_pr(j)) = c_pr(re(i)%gi_pr(j)) + (f_pr(re(i)%gi_pr(j)) - t_pr(re(i)%gi_pr(j))) + net_pr(i)
        else
          c_pr(re(i)%gi_pr(j)) = c_pr(re(i)%gi_pr(j)) + (net_pr(i) - t_pr(re(i)%gi_pr(j))) + f_pr(re(i)%gi_pr(j))
        end if
        f_pr(re(i)%gi_pr(j)) = t_pr(re(i)%gi_pr(j))
      end do
      f_pr(re(i)%gi_pr(:)) =  f_pr(re(i)%gi_pr(:)) + c_pr(re(i)%gi_pr(:))

      !! Add the reactant rates
      do j = 1, re(i)%n_re
        t_re(re(i)%gi_re(j)) = f_re(re(i)%gi_re(j)) + net_re(i)
        if (abs(f_re(re(i)%gi_re(j))) >= abs(net_re(i))) then
          c_re(re(i)%gi_re(j)) = c_re(re(i)%gi_re(j)) + (f_re(re(i)%gi_re(j)) - t_re(re(i)%gi_re(j))) + net_re(i)
        else
          c_re(re(i)%gi_re(j)) = c_re(re(i)%gi_re(j)) + (net_re(i) - t_re(re(i)%gi_re(j))) + f_re(re(i)%gi_re(j))
        end if
        f_re(re(i)%gi_re(j)) = t_re(re(i)%gi_re(j))
      end do
      f_re(re(i)%gi_re(:)) =  f_re(re(i)%gi_re(:)) + c_re(re(i)%gi_re(:))

    end do

    !! Sum product and reactant rates to get net rate for species
    f(:) = f_pr(:) + f_re(:)
 
  end subroutine RHS_update

  subroutine jac_dummy (NEQ, X, Y, ML, MU, PD, NROWPD)
    integer, intent(in) :: NEQ, ML, MU, NROWPD
    real(dp), intent(in) :: X
    real(dp), dimension(NEQ), intent(in) :: Y
    real(dp), dimension(NROWPD, NEQ), intent(inout) :: PD
  end subroutine jac_dummy

void jac_NCHO(N, X, Y, ML, MU, DFY, NROWPD){

    dfy[0][0] = -f[0]*y[1] - f[1]*y[4] - r[2]*y[3];
    dfy[0][1] = -f[0]*y[0] + f[2]*y[6];
    dfy[0][2] = r[0]*y[3];
    dfy[0][3] = r[0]*y[2] + re_r[1]*y[5] - re_r[2]*y[0];
    dfy[0][4] = -f[1]*y[0];
    dfy[0][5] = r[1]*y[3];
    dfy[0][6] = f[2]*y[1];
    dfy[0][7] = 0.0;
    dfy[0][8] = 0.0;
    dfy[0][9] = 0.0;
    dfy[0][10] = 0.0;
    dfy[0][11] = 0.0;

    dfy[1][0] = -f[0]*y[1] + r[2]*y[3];
    dfy[1][1] = -nd_atm*r[3] - 9.0*r[4]*y[1]*y[1]*y[4] - 
            9.0*r[5]*y[1]*y[1]*y[8] - 9.0*r[7]*y[10]*y[1]*y[1] - 
            9.0*r[8]*y[11]*y[1]*y[1] - f[0]*y[0] - f[2]*y[6];
    dfy[1][2] = r[0]*y[3] + 3.0*f[4]*y[7];
    dfy[1][3] = 2.0*nd_atm*f[3]*y[3] + r[0]*y[2] + r[2]*y[0];
    dfy[1][4] = -3.0*r[4]*y[2]*y[2]*y[2];
    dfy[1][5] = 0.0;
    dfy[1][6] = -f[2]*y[1];

    dfy[1][7] = 6.0*f[5]*y[7] + 3.0*f[8]*y[9] + 3.0*f[4]*y[2];
    dfy[1][8] = -3.0*r[5]*y[1]*y[1]*y[1];
    dfy[1][9] = 6.0*f[7]*y[9] + 3.0*f[8]*y[7];
    dfy[1][10] = -3.0*r[7]*y[1]*y[1]*y[1];
    dfy[1][11] = -3.0*r[8]*y[1]*y[1]*y[1];

    dfy[2][0] = f[0]*y[1];
    dfy[2][1] = 3.0*r[4]*y[1]*y[1]*y[4] + f[0]*y[0];
    dfy[2][2] = -r[6]*y[8] - re_r[9]*y[11] - r[0]*y[3] - f[4]*y[7];
    dfy[2][3] = -r[0]*y[2];
    dfy[2][4] = r[4]*y[1]*y[1]*y[1] + f[6]*y[7] + f[9]*y[9];
    dfy[2][5] = 0.0;
    dfy[2][6] = 0.0;
    dfy[2][7] = f[6]*y[4] - f[4]*y[2];
    dfy[2][8] = -re_r[6]*y[2];
    dfy[2][9] = f[9]*y[4];
    dfy[2][10] = 0.0;
    dfy[2][11] = -r[9]*y[2];

    dfy[3][0] = f[0]*y[1] + f[1]*y[4] - r[2]*y[3];
    dfy[3][1] = 2.0*nd_atm*r[3] + f[0]*y[0] + f[2]*y[6];
    dfy[3][2] = -r[0]*y[3];
    dfy[3][3] = -4.0*nd_atm*f[3]*y[3] - re_r[0]*y[2] - 
            r[1]*y[5] - r[2]*y[0];
    dfy[3][4] = f[1]*y[0];
    dfy[3][5] = -r[1]*y[3];
    dfy[3][6] = f[2]*y[1];
    dfy[3][7] = 0.0;
    dfy[3][8] = 0.0;
    dfy[3][9] = 0.0;
    dfy[3][10] = 0.0;
    dfy[3][11] = 0.0;

    dfy[4][0] = -f[1]*y[4];
    dfy[4][1] = -3.0*r[4]*y[1]*y[1]*y[4];
    dfy[4][2] = r[6]*y[8] + r[9]*y[11] + f[4]*y[7];
    dfy[4][3] = r[1]*y[5];
    dfy[4][4] = -r[4]*y[1]*y[1]*y[1] - f[6]*y[7] - f[9]*y[9] - f[1]*y[0];
    dfy[4][5] = r[1]*y[3];
    dfy[4][6] = 0.0;
    dfy[4][7] = -f[6]*y[4] + f[4]*y[2];
    dfy[4][8] = r[6]*y[2];
    dfy[4][9] = -f[9]*y[4];
    dfy[4][10] = 0.0;
    dfy[4][11] = r[9]*y[2];

    dfy[5][0] = f[1]*y[4];
    dfy[5][1] = 0.0;
    dfy[5][2] = 0.0;
    dfy[5][3] = -r[1]*y[5];
    dfy[5][4] = f[1]*y[0];
    dfy[5][5] = -r[1]*y[3];
    dfy[5][6] = 0.0;
    dfy[5][7] = 0.0;
    dfy[5][8] = 0.0;
    dfy[5][9] = 0.0;
    dfy[5][10] = 0.0;
    dfy[5][11] = 0.0;

dfy[6][0] = r[2]*y[5];
dfy[6][1] = -r[2]*y[6];
dfy[6][2] = 0.0;
dfy[6][3] = r[2]*y[0];
dfy[6][4] = 0.0;
dfy[6][5] = 0.0;
dfy[6][6] = -f[2]*y[1];
dfy[6][7] = 0.0;
dfy[6][8] = 0.0;
dfy[6][9] = 0.0;
dfy[6][10] = 0.0;
dfy[6][11] = 0.0;

dfy[7][0] = 0.0;
dfy[7][1] = 3.0*r[4]*y[1]*y[1]*y[4] + 6.0*r[5]*y[1]*y[2]*y[8] + 
            3.0*r[8]*y[1]*y[1]*y[11];
dfy[7][2] = r[6]*y[8] - f[4]*y[7];
dfy[7][3] = 0.0;
dfy[7][4] = r[4]*y[1]*y[1]*y[1] - f[6]*y[7];
dfy[7][5] = 0.0;
dfy[7][6] = 0.0;
dfy[7][7] = -4.0*f[5]*y[7] - f[6]*y[4] - f[8]*y[9] - f[4]*y[2];
dfy[7][8] = 2.0*r[5]*y[1]*y[1]*y[1] + r[6]*y[2];
dfy[7][9] = -f[8]*y[7];
dfy[7][10] = 0.0;
dfy[7][11] = r[8]*y[1]*y[1]*y[1];

dfy[8][0] = 0.0;
dfy[8][1] = -3.0*r[5]*y[1]*y[1]*y[8];
dfy[8][2] = -r[6]*y[8];
dfy[8][3] = 0.0;
dfy[8][4] = f[6]*y[7];
dfy[8][5] = 0.0;
dfy[8][6] = 0.0;
dfy[8][7] = 2.0*f[5]*y[7] + f[6]*y[4];
dfy[8][8] = -r[5]*y[1]*y[1] - r[6]*y[2];
dfy[8][9] = 0.0;
dfy[8][10] = 0.0;
dfy[8][11] = 0.0;

dfy[9][0] = 0.0;
dfy[9][1] = 6.0*r[7]*y[10]*y[1]*y[1] + 3.0*r[8]*y[11]*y[1]*y[1];
dfy[9][2] = r[9]*y[11];
dfy[9][3] = 0.0;
dfy[9][4] = -f[9]*y[9];
dfy[9][5] = 0.0;
dfy[9][6] = 0.0;
dfy[9][7] = -f[9]*y[9];
dfy[9][8] = 0.0;
dfy[9][9] = -4.0*f[7]*y[9] - r[8]*y[7] - f[9]*y[4];
dfy[9][10] = 2.0*r[7]*y[1]*y[1]*y[1];
dfy[9][12] = r[8]*y[1]*y[1]*y[1] + r[9]*y[2];

dfy[10][0] = 0.0;
dfy[10][1] = -3.0*r[7]*y[1]*y[1]*y[10];
dfy[10][2] = 0.0;
dfy[10][3] = 0.0;
dfy[10][4] = 0.0;
dfy[10][5] = 0.0;
dfy[10][6] = 0.0;
dfy[10][7] = 0.0;
dfy[10][8] = 0.0;
dfy[10][9] = 2.0*f[7]*y[9];
dfy[10][10] = -r[7]*y[1]*y[1]*y[1];
dfy[10][11] = 0.0;

dfy[11][0] = 0.0;
dfy[11][1] = -3.0*r[8]*y[11]*y[1]*y[1];
dfy[11][2] = -r[9]*y[11];
dfy[11][3] = 0.0;
dfy[11][4] = f[9]*y[9];
dfy[11][5] = 0.0;
dfy[11][6] = 0.0;
dfy[11][7] = f[8]*y[9];
dfy[11][8] = 0.0;
dfy[11][9] = f[8]*y[7] + f[9]*y[4];
dfy[11][10] = 0.0;
dfy[11][11] = -r[8]*y[1]*y[1]*y[1] - r[9]*y[2];
}

void jac_CHO(N, X, Y, ML, MU, DFY, NROWPD){

    dfy[0][0] = -f[0]*y[1] - f[1]*y[4] - r[2]*y[3]
    dfy[0][1] = -f[0]*y[0] + f[2]*y[6]
    dfy[0][2] = r[0]*y[3]
    dfy[0][3] = r[0]*y[2] + r[1]*y[5] - r[2]*y[0]
    dfy[0][4] = -f[1]*y[0]
    dfy[0][5] = r[1]*y[3]
    dfy[0][6] = f[2]*y[1]
    dfy[0][7] = 0.0
    dfy[0][8] = 0.0
    
    dfy[1][0] = -f[0]*y[1] + r[2]*y[3]
    dfy[1][1] = -nd_atm*r[3] - 9.0*r[4]*y[1]*y[1]*y[4] - 
      9.0*r[5]*y[1]*y[1]*y[8] - f[0]*y[0] - f[2]*y[6]
    dfy[1][2] = r[0]*y[3] + 3.0*f[4]*y[7]
    dfy[1][3] = 2.0*nd_atm*f[3]*y[3] + r[0]*y[2] + r[2]*y[0]
    dfy[1][4] = -3.0*r[4]*y[1]*y[1]*y[1]
    dfy[1][5] = 0.0
    dfy[1][6] = -f[2]*y[1]
    dfy[1][7] = 6.0*f[5]*y[7] + 3.0*f[4]*y[2]
    dfy[1][8] = -3.0*r[5]*y[1]*y[1]*y[1]
    
    dfy[2][0] = f[0]*y[1]
    dfy[2][1] = 3.0*r[4]*y[1]*y[1]*y[4] + f[0]*y[0]
    dfy[2][2] = -r[6]*y[8] - r[0]*y[3] - f[4]*y[7]
    dfy[2][3] = -r[0]*y[2]
    dfy[2][4] = r[4]*y[1]*y[1]*y[1] + f[6]*y[7]
    dfy[2][5] = 0.0
    dfy[2][6] = 0.0
    dfy[2][7] = f[6]*y[4] - f[4]*y[2]
    dfy[2][8] = -r[6]*y[2]
    
    dfy[3][0] = f[0]*y[1] + f[1]*y[6] - r[2]*y[3]
    dfy[3][1] = 2.0*nd_atm*r[3] + f[0]*y[0] + f[2]*y[6]
    dfy[3][2] = -r[0]*y[3]
    dfy[3][3] = -4.0*nd_atm*f[3]*y[3] - r[0]*y[2] - r[1]*y[5] - r[2]*y[0]
    dfy[3][4] = f[1]*y[0]
    dfy[3][5] = -r[1]*y[3]
    dfy[3][6] = f[2]*y[1]
    dfy[3][7] = 0.0
    dfy[3][8] = 0.0
   
    dfy[4][0] = -f[1]*y[4]
    dfy[4][1] = -3.0*r[4]*y[1]*y[1]*y[4]
    dfy[4][2] = r[6]*y[8] + f[4]*y[7];
    dfy[4][3] = r[1]*y[5];
    dfy[4][4] = -3.0*r[4]*y[2]*y[2]*y[5] - f[6]*y[7] - f[1]*y[6];
    dfy[4][5] = r[2]*y[5];
    dfy[4][6] = 0.0;
    dfy[4][7] = -f[6]*y[4] + f[4]*y[2];
    dfy[4][8] = r[6]*y[2];
   
    dfy[5][0] = f[1]*y[4];
    dfy[5][1] = 0.0;
    dfy[5][2] = 0.0;
    dfy[5][3] = -r[1]*y[5];
    dfy[5][4] = f[1]*y[0];
    dfy[5][5] = -r[1]*y[3];
    dfy[5][6] = 0.0;
    dfy[5][7] = 0.0;
    dfy[5][8] = 0.0;

    dfy[5][0] = r[2]*y[3];  
    dfy[5][1] = 0.0;
    dfy[5][2] = 0.0;
    dfy[5][3] = -r[1]*y[5];
    dfy[5][4] = f[1]*y[0];
    dfy[5][5] = -r[1]*y[3];
    dfy[5][6] = 0.0;
    dfy[5][7] = 0.0;
    dfy[5][8] = 0.0;
  
    dfy[6][0] = r[2]*y[3];
    dfy[6][1] = -f[2]*y[6];
    dfy[6][2] = 0.0;
    dfy[6][3] = r[2]*y[0];
    dfy[6][4] = 0.0;
    dfy[6][5] = 0.0;
    dfy[6][6] = -f[2]*y[1];
    dfy[6][7] = 0.0;
    dfy[6][8] = 0.0;
  
    dfy[7][0] = 0.0;
    dfy[7][1] = 3.0*r[4]*y[1]*y[1]*y[4] + 6.0*r[5]*y[2]*y[1]*y[8];
    dfy[7][2] = r[6]*y[8] - f[4]*y[7];
    dfy[7][3] = 0.0;
    dfy[7][4] = r[4]*y[1]*y[1]*y[1] - f[6]*y[7];
    dfy[7][5] = 0.0;
    dfy[7][6] = 0.0;
    dfy[7][7] = -4.0*f[5]*y[7] - f[6]*y[4] - f[4]*y[2];
    dfy[7][8] = 2.0*r[5]*y[1]*y[1]*y[1] + r[6]*y[2];
  
    dfy[8][0] = 0.0;
    dfy[8][1] = -3.0*r[5]*y[1]*y[1]*y[8];
    dfy[8][2] = -r[6]*y[8];
    dfy[8][3] = 0.0;
    dfy[8][4] = f[6]*y[7];
    dfy[8][5] = 0.0;
    dfy[8][6] = 0.0;
    dfy[8][7] = 2.0*f[5]*y[7] + f[6]*y[4];
    dfy[8][8] = -r[5]*y[1]*y[1]*y[1] - r[6]*y[2];
}

void jac_HO(double *N, double *X, Y, ML, MU, DFY, NROWPD)

    dfy[0][0] = -f[0]*y[1] - r[1]*y[3];
    dfy[0][1] = -f[0]*y[0] + f[1]*y[4];
    dfy[0][2] = r[0]*y[3];
    dfy[0][3] = r[0]*y[2] - r[1]*y[0];
    dfy[0][4] = f[1]*y[1];
   
    dfy[1][0] = -f[0]*y[1] + r[1]*y[3];
    dfy[1][1] = -nd_atm*r[2] - f[0]*y[0] - f[1]*y[4];
    dfy[1][2] = r[0]*y[3];
    dfy[1][3] = 2.0*nd_atm*f[2]*y[3] + r[0]*y[2] + r[1]*y[0];
    dfy[1][4] = -f[1]*y[1];
  
    dfy[2][0] = f[0]*y[1];
    dfy[2][1] = f[0]*y[0];
    dfy[2][2] = -r[0]*y[3];
    dfy[2][3] = -r[0]*y[2];
    dfy[2][4] = 0.0;

    dfy[3][0] = f[0]*y[1] - r[1]*y[3];
    dfy[3][1] = 2.0*nd_atm*r[2] + f[0]*y[0] + f[1]*y[4];
    dfy[3][2] = -r[0]*y[3];
    dfy[3][3] = -4.0*nd_atm*f[2]*y[3] - r[0]*y[2] - r[1]*y[0];
    dfy[3][4] = f[1]*y[1];

    dfy[4][0] = r[1]*y[3];
    dfy[4][1] = -f[1]*y[4];
    dfy[4][2] = 0.0;
    dfy[4][3] = r[1]*y[0];
    dfy[4][4] = -f[1]*y[1];
}
