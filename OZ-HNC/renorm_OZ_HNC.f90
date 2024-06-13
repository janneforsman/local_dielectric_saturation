program renorm_OZ
    implicit none
    integer(kind=4) :: i, restart, nn
    integer(kind=4), parameter :: N = 262144 !2**18 points
    real(kind=8), parameter :: pi = acos(-1.d0) 
    real(kind=8), parameter :: kb = 1.380649D-23
    real(kind=8), parameter :: el_charge = 1.602176634D-19
    real(kind=8), parameter :: permittivity = 8.85418782D-12
    real(kind=8), parameter :: Na = 6.02214076D23

    real(kind=8) :: sigma_p, sigma_n, sigma_pn, delr, delk, rho_p, rho_n, rho_tot, mixing, cost, kappa, T, beta
    real(kind=8) :: eps_bulk, z_p, z_n, lb_bulk, iteration_counter, int1, int2, int3, int4
    real(kind=8) :: int5, int6, phi, mu_p, mu_n, kappa_d

    real(kind=8), dimension(N) :: r, k
    real(kind=8), dimension(N) :: us_pp, us_pn, us_nn
    real(kind=8), dimension(N) :: gr_pp, gr_pn, gr_nn
    real(kind=8), dimension(N) :: hr_pp, hr_pn, hr_nn
    real(kind=8), dimension(N) :: cr_pp, cr_pn, cr_nn
    real(kind=8), dimension(N) :: cr_0pp, cr_0pn, cr_0nn, ck_0pp, ck_0pn, ck_0nn
    real(kind=8), dimension(N) :: taur_pp, taur_pn, taur_nn, tauk_pp, tauk_pn, tauk_nn
    real(kind=8), dimension(N) :: taur_0pp, taur_0pn, taur_0nn, tauk_0pp, tauk_0pn, tauk_0nn
    real(kind=8), dimension(N) :: qr_pp, qr_pn, qr_nn, qk_pp, qk_pn, qk_nn
    real(kind=8), dimension(N) :: Br_pp, Br_pn, Br_nn
    real(kind=8), dimension(N) :: phi_pp, phi_pn, phi_nn

    open(unit=10, file='OZ_input.txt', action='read', form='formatted')
    read(10,*)
    read(10,*) delr
    read(10,*)
    read(10,*) sigma_p, sigma_n
    read(10,*)
    read(10,*) T
    read(10,*)
    read(10,*) eps_bulk
    read(10,*)
    read(10,*) rho_p, rho_n
    read(10,*)
    read(10,*) restart
    read(10,*)
    read(10,*) mixing
    close(unit=10)

    nn = N
    beta = 1.d0 / (kb * T)
    z_p = 1.d0
    z_n = -1.d0
    sigma_pn = (sigma_p + sigma_n) / 2.d0
    rho_p = (rho_p * 1D3 * Na) / 1D30
    rho_n = (rho_n * 1D3 * Na) / 1D30
    rho_tot = rho_p + rho_n
    lb_bulk = (1D10 * beta * el_charge ** 2.d0)/ (4.d0 * pi * permittivity * eps_bulk)
    kappa_d = dsqrt(4.d0 * pi * lb_bulk * (rho_p * z_p ** 2.d0 + rho_n * z_n ** 2.d0))
    delk = pi / (delr * float(nn))
    cost = 100.d0
    iteration_counter = 0

    write(*,*) rho_p
    write(*,*) 'kappa_debye for the system is: ', kappa_d

    kappa = kappa_d
    !r,k creation and potential definitions
    do i = 1, N
        r(i) = float(i) * delr
        if (r(i) < sigma_p) then
            us_pp(i) = 1.0d100
        else
            us_pp(i) = 0.0
        end if
        if (r(i) < sigma_n) then
            us_nn(i) = 1.0d100
        else
            us_nn(i) = 0.0
        end if
        if (r(i) < sigma_pn) then
            us_pn(i) = 1.0d100
        else
            us_pn(i) = 0.0
        end if
    end do
    k = delk * r / delr

    !Bridge function setup (HNC)
    Br_pp = 0.0
    Br_pn = 0.0
    Br_nn = 0.0

    !q function definitions
    qr_pp = -lb_bulk * z_p * z_p * dexp(-kappa * r) / r
    qr_pn = -lb_bulk * z_p * z_n * dexp(-kappa * r) / r
    qr_nn = -lb_bulk * z_n * z_n * dexp(-kappa * r) / r
    qk_pp = -4.d0 * pi * lb_bulk * z_p * z_p / (k ** 2.d0 + kappa ** 2.d0)
    qk_pn = -4.d0 * pi * lb_bulk * z_p * z_n / (k ** 2.d0 + kappa ** 2.d0)
    qk_nn = -4.d0 * pi * lb_bulk * z_n * z_n / (k ** 2.d0 + kappa ** 2.d0)

    !total (phi) potential definition
    phi_pp = -lb_bulk * z_p * z_p / r
    phi_pn = -lb_bulk * z_p * z_n / r
    phi_nn = -lb_bulk * z_n * z_n / r

    !Define tau functions
    if (restart == 0) then
        taur_0pp = 0.0d0
        taur_0pn = 0.0d0
        taur_0nn = 0.0d0
        else
        open(555,file='tau.txt')
        do i = 1,N
           read(555,*) taur_0pp(i), taur_0pn(i), taur_0nn(i)
        end do
        close(555)
    end if

    !Preparation for FFT
    taur_0pp = taur_0pp * r 
    taur_0pn = taur_0pn * r
    taur_0nn = taur_0nn * r

    !FFT to k-space
    call FFT(1, taur_0pp, tauk_0pp, nn, delr)
    call FFT(1, taur_0pn, tauk_0pn, nn, delr)
    call FFT(1, taur_0nn, tauk_0nn, nn, delr)

    !Post-FFT cleanup
    taur_0pp = taur_0pp / r 
    taur_0pn = taur_0pn / r
    taur_0nn = taur_0nn / r
    tauk_0pp = tauk_0pp / k
    tauk_0pn = tauk_0pn / k
    tauk_0nn = tauk_0nn / k

    !Closure relation (HNC)
    cr_0pp = dexp(-us_pp + taur_0pp + qr_pp + Br_pp) - 1.d0 - taur_0pp - qr_pp
    cr_0pn = dexp(-us_pn + taur_0pn + qr_pn + Br_pn) - 1.d0 - taur_0pn - qr_pn
    cr_0nn = dexp(-us_nn + taur_0nn + qr_nn + Br_nn) - 1.d0 - taur_0nn - qr_nn

    !Preparation for FFT
    cr_0pp = cr_0pp * r 
    cr_0pn = cr_0pn * r
    cr_0nn = cr_0nn * r

    !FFT to k-space
    call FFT(1, cr_0pp, ck_0pp, nn, delr)
    call FFT(1, cr_0pn, ck_0pn, nn, delr)
    call FFT(1, cr_0nn, ck_0nn, nn, delr)

    !Post-FFT cleanup
    cr_0pp = cr_0pp / r 
    cr_0pn = cr_0pn / r
    cr_0nn = cr_0nn / r
    ck_0pp = ck_0pp / k
    ck_0pn = ck_0pn / k
    ck_0nn = ck_0nn / k

    do while ((cost > 1.d-6) .and. (iteration_counter < 1D5))
        iteration_counter = iteration_counter + 1.d0

        call tau_comp(nn, ck_0pp, ck_0pn, ck_0nn, qk_pp, qk_pn, qk_nn, rho_p, rho_n, tauk_pp, tauk_pn, tauk_nn)

        !Preparation for iFFT
        tauk_pp = tauk_pp * k
        tauk_pn = tauk_pn * k
        tauk_nn = tauk_nn * k

        !iFFT to r-space
        call FFT(-1, tauk_pp, taur_pp, nn, delk)
        call FFT(-1, tauk_pn, taur_pn, nn, delk)
        call FFT(-1, tauk_nn, taur_nn, nn, delk)

        !Post-FFT cleanup
        taur_pp = taur_pp / r 
        taur_pn = taur_pn / r
        taur_nn = taur_nn / r
        tauk_pp = tauk_pp / k
        tauk_pn = tauk_pn / k
        tauk_nn = tauk_nn / k

        cost = 0.d0
        do i = 1,nn
            cost = cost + (taur_pp(i) - taur_0pp(i))**2.0d0 + (taur_pn(i) - taur_0pn(i))**2.0d0 + (taur_nn(i) - taur_0nn(i))**2.0d0
        end do
        cost = sqrt(cost) * delr
        write(*,'(a9,i8,e15.5)') 'Cycle # =',int(iteration_counter), cost

        !Picard mixing
        taur_0pp = (1.0d0 - mixing) * taur_pp + mixing * taur_0pp
        taur_0pn = (1.0d0 - mixing) * taur_pn + mixing * taur_0pn
        taur_0nn = (1.0d0 - mixing) * taur_nn + mixing * taur_0nn

        !Closure relation (HNC)
        cr_0pp = dexp(-us_pp + taur_0pp + qr_pp + Br_pp) - 1.d0 - taur_0pp - qr_pp
        cr_0pn = dexp(-us_pn + taur_0pn + qr_pn + Br_pn) - 1.d0 - taur_0pn - qr_pn
        cr_0nn = dexp(-us_nn + taur_0nn + qr_nn + Br_nn) - 1.d0 - taur_0nn - qr_nn

        !Preparation for FFT
        cr_0pp = cr_0pp * r 
        cr_0pn = cr_0pn * r
        cr_0nn = cr_0nn * r

        !FFT to k-space
        call FFT(1, cr_0pp, ck_0pp, nn, delr)
        call FFT(1, cr_0pn, ck_0pn, nn, delr)
        call FFT(1, cr_0nn, ck_0nn, nn, delr)

        !Post-FFT cleanup
        cr_0pp = cr_0pp / r 
        cr_0pn = cr_0pn / r
        cr_0nn = cr_0nn / r
        ck_0pp = ck_0pp / k
        ck_0pn = ck_0pn / k
        ck_0nn = ck_0nn / k
    end do

    !Total correlation function
    hr_pp = taur_0pp + cr_0pp + qr_pp
    hr_pn = taur_0pn + cr_0pn + qr_pn
    hr_nn = taur_0nn + cr_0nn + qr_nn

    !Direct correlation function
    cr_pp = cr_0pp + phi_pp
    cr_pn = cr_0pn + phi_pn
    cr_nn = cr_0nn + phi_nn

    !RDF
    gr_pp = hr_pp + 1.0d0
    gr_pn = hr_pn + 1.0d0
    gr_nn = hr_nn + 1.0d0

    !Saving the RDF
    do i = 1, N
        if (gr_pp(i) < 1.0d-5) gr_pp(i) = 0.0d0
        if (gr_pn(i) < 1.0d-5) gr_pn(i) = 0.0d0
        if (gr_nn(i) < 1.0d-5) gr_nn(i) = 0.0d0
    end do

    open(10,file='rdf_' // real_to_string(sigma_p) // 'A_' // real_to_string(rho_p) //'.txt')     
    do i = 1, N
        write(10,'(7e20.8)') r(i), gr_pp(i), gr_pn(i), gr_nn(i)
    end do
    close (10)

    open(556,file='tau.txt')     
    do i = 1, N
        write(556,'(3e20.10)') taur_0pp(i), taur_0pn(i), taur_0nn(i)
    end do
    close (556)

    !Chemical potential via Hansen-Vieillefosse-Belloni equation
    call simpson((N - 1) / 2, r, cr_pp * r ** 2.d0, int1)
    call simpson((N - 1) / 2, r, cr_pn * r ** 2.d0, int2)
    call simpson((N - 1) / 2, r, cr_nn * r ** 2.d0, int3)
    call simpson((N - 1) / 2, r, hr_pp * (hr_pp - cr_pp) * r ** 2.d0, int4)
    call simpson((N - 1) / 2, r, hr_pn * (hr_pn - cr_pn) * r ** 2.d0, int5)
    call simpson((N - 1) / 2, r, hr_nn * (hr_nn - cr_nn) * r ** 2.d0, int6)

    phi = 0.5d0
    mu_p = 4.d0 * pi * (-rho_p*int1 - rho_n*int2 + 0.5d0*(rho_p*int4 + rho_n*int5))
    mu_n = 4.d0 * pi * (-rho_n*int3 - rho_p*int2 + 0.5d0*(rho_n*int6 + rho_p*int5))

    write(*,*) 'mu_p = ', mu_p
    write(*,*) 'mu_n = ', mu_n

contains

! *******************************************************************

    SUBROUTINE FFT(KAM,A,B,N,dDR)
        !
        !     FAST FOURIER ROUTINE  - IN ONE DIMENSION
        !
        !       B(K)=KOEF*DR*SUMA( SIN(I*K*PI/N) * A(I) )
        !
        !               N = 2 ** NM
        !               DK*DR = PI/N
        !               KOEF  = 4*PI    DIRECT TRANS.  KAM > 0
        !                     = 2/PI**2 INVERS. TRANS. KAM < 0
        !
        !               !!!!!!! A(N) MUST BE ZERO !!!!!!!
        !
        implicit none        
        integer(4):: I,J,K,L,L1,L2,N,NM,KAM
        real(8), dimension(N):: A,B
        real(8), dimension(2*N):: AA,BB
        real(8):: P1,P2,P3,P4,P5,P6,P7,P8,COEF,dDR
        
                NM = nint(log(float(N))/log(2.0d0))
        
        !        N=2**NM
        
                AA(1)=0.D0
                AA(N+1)=0.D0
                BB(1)=0.D0
                BB(N+1)=0.D0
                DO 10 I=2,N
                AA(I)=A(I-1)
                AA(N+I)=0.D0
                BB(I)=0.D0
        10      BB(N+I)=0.D0
                NM=NM+1
                N=N+N
                J=1
                L=N-1
                DO 40 I=1,L
                IF (I.GE.J) GO TO 20
                P1=AA(J)
                AA(J)=AA(I)
                AA(I)=P1
        20      K=N/2
        30      IF (K.GE.J) GO TO 40
                J=J-K
                K=K/2
                GO TO 30
        40      J=J+K
                DO 60 L=1,NM
                L1=2**L
                L2=L1/2
                P1=1.D0
                P2=0.D0
                P3=dCOS(PI/L2)
                P4=-dSIN(PI/L2)
                DO 60 J=1,L2
                DO 50 I=J,N,L1
                P7=AA(I+L2)
                P8=BB(I+L2)
                P5=P7*P1-P8*P2
                P6=P7*P2+P8*P1
                AA(I+L2)=AA(I)-P5
                AA(I)=AA(I)+P5
                BB(I+L2)=BB(I)-P6
        50      BB(I)=BB(I)+P6
                P5=P1
                P1=P1*P3-P2*P4
        60      P2=P5*P4+P2*P3
                NM=NM-1
                N=N/2
                COEF=-4.D0*PI*dDR
                IF(KAM.LT.0) COEF=-0.5D0/PI/PI*dDR
                DO 70 I=2,N
        70      B(I-1)=COEF*BB(I)
                B(N)=0.0D0
                RETURN
                END

    subroutine tau_comp(nop,ck11,ck12,ck22,qk11,qk12,qk22,density_1,density_2,tauk11,tauk12,tauk22)
        implicit none
        integer(4):: nop
        real(8), dimension(nop):: ck11,ck12,ck22
        real(8), dimension(nop):: qk11,qk12,qk22
        real(8), dimension(nop):: tauk11,tauk12,tauk22
        real(8), dimension(nop):: rtaur11,rtaur12,rtaur22
        real(8):: density_1,density_2
        
        !rho*tau*rho = Vc * [ I - Vc ]^-1 * V - rho * c * rho
        
        rtaur11 = (density_1**2.0d0*(ck11*density_1*(ck11 + 2.0d0*qk11 + qk11*(ck11 + qk11)*density_1) + &
             & (2.0d0*ck12*qk12*(1.0d0 + (ck11 + qk11)*density_1) + ck12**2.0d0*(1.0d0 + qk11*density_1)*(1.0d0 + & 
             & (ck11 + qk11)*density_1)  - ck11*ck22*density_1*(ck11 + 2.0d0*qk11 + qk11*(ck11 + qk11)*density_1))*density_2 + &
             & (ck12**2.0d0*(1.0d0 + (ck11 + qk11)*density_1)*(qk22 - qk12**2.0d0*density_1 + qk11*qk22*density_1) + &
             & ck22*(-ck11*qk22*density_1*(ck11 + 2.0d0*qk11 + qk11*(ck11 + qk11)*density_1) + &
             & qk12**2.0d0*(1 + ck11*density_1*(1.0d0 + (ck11 + qk11)*density_1))))*density_2**2.0d0))/&
             & (1.0d0 + ck11*density_1*(-1.0d0 + ck22*density_2*(1.0d0 + qk22*density_2 - qk12**2.0d0*density_1*density_2) + &
             & qk11*density_1*(-1.0d0 + ck22*density_2*(1.0d0 + qk22*density_2))) - density_2*(ck22 + ck22*qk22*density_2 + &
             & ck12*density_1*(2.0d0*qk12 + ck12*(1.0d0 + qk11*density_1 + (qk22 - &
             & qk12**2.0d0*density_1 + qk11*qk22*density_1)*density_2))))
        
        rtaur12 = (density_1*density_2*(-ck12*qk11*density_1 - (ck12*ck22 + ck22*qk12 + ck12*qk22 + &
             & ck12*(ck12**2.0d0 + 3.0d0*ck12*qk12 + qk12**2.0d0 + qk11*qk22)*density_1 + &
             & ck12**2.0d0*qk11*(ck12 + qk12)*density_1**2.0d0)*density_2 - (ck12 + qk12)*(ck22*qk22 + &
             & ck12**2.0d0*density_1*(qk22 - qk12**2.0d0*density_1 + qk11*qk22*density_1))*density_2**2.0d0 - &
             & ck11*(ck12 + qk12)*density_1*(1.0d0 + ck22*density_2*(-1.0d0 - qk22*density_2 + qk12**2.0d0*density_1*density_2) - &
             & qk11*density_1*(-1.0d0 + ck22*density_2*(1.0d0 + qk22*density_2)))))/(-1.0d0 + &
             & ck11*density_1*(1.0d0 + ck22*density_2*(-1.0d0 - qk22*density_2 + qk12**2.0d0*density_1*density_2) - &
             & qk11*density_1*(-1.0d0 + ck22*density_2*(1.0d0 + qk22*density_2))) + density_2*(ck22 + ck22*qk22*density_2 + &
             & ck12*density_1*(2.0d0*qk12 + ck12*(1.0d0 + qk11*density_1 + &
             & (qk22 - qk12**2.0d0*density_1 + qk11*qk22*density_1)*density_2))))
        
        rtaur22 = -(density_2**2.0d0*(2.0d0*ck12*qk12*density_1*(1.0d0 + (ck22 + qk22)*density_2) + &
             & ck22*density_2*(ck22 + 2.0d0*qk22 + qk22*(ck22 + qk22)*density_2) + &
             & ck12**2.0d0*density_1*(1.0d0 + (ck22 + qk22)*density_2)*(1.0d0 + qk11*density_1 + &
             & (qk22 - qk12**2.0d0*density_1 + qk11*qk22*density_1)*density_2) + &
             & ck11*density_1*(-ck22*(1.0d0 + qk11*density_1)*density_2*(ck22 + 2.0d0*qk22 + &
             & qk22*(ck22 + qk22)*density_2) + qk12**2.0d0*density_1*(1.0d0 + ck22*density_2*(1.0d0 + (ck22 + qk22)*density_2)))))/&
             & (-1.0d0 + ck11*density_1*(1.0d0 + ck22*density_2*(-1.0d0 - qk22*density_2 + qk12**2.0d0*density_1*density_2) - &
             & qk11*density_1*(-1.0d0 + ck22*density_2*(1.0d0 + qk22*density_2))) + density_2*(ck22 + ck22*qk22*density_2 + &
             & ck12*density_1*(2.0d0*qk12 + ck12*(1.0d0 + qk11*density_1 + (qk22 - &
             & qk12**2.0d0*density_1 + qk11*qk22*density_1)*density_2))))
        
        
        !final tau calculation
        tauk11 = rtaur11*density_2**2.0d0/((density_1*density_2)**2.0d0)
        tauk12 = rtaur12*density_1*density_2/((density_1*density_2)**2.0d0)
        tauk22 = rtaur22*density_1**2.0d0/((density_1*density_2)**2.0d0)
        
    end subroutine tau_comp

    subroutine simpson(M,x_space,y_space,int)
        implicit none
        integer :: ii, i, N, M
        real(kind=8) :: int, simp
        real(kind=8), dimension(M) :: x_space, y_space
        
        int = 0.0d0
        do ii = 1, M-2, 2
           simp = ((x_space(ii+2)-x_space(ii))/6.0d0)*(y_space(ii) + 4.0d0*y_space(ii+1) + y_space(ii+2))
           int = int + simp
        end do
        
    end subroutine simpson

    pure function real_to_string(x) result(s)
    real(8), intent(in) :: x
    character (len=:), allocatable :: s
    character (len=100) :: str
    character (len=*), parameter :: fmt = "(e0.3)" ! choose the desired format
    write (str,fmt) x ! convert x to fixed-length string
    s = trim(str) ! create a string with no trailing spaces
    end function real_to_string
          
end program renorm_OZ