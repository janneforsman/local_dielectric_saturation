program ionic_fluid
    implicit none  

    integer :: NoP, NoPh, max_iterations, RDF_sample, RDF_points, max_subsimulations, initialization
    integer :: move_att
    real(kind=8) :: L, sigma, T, beta, lb, step, dR, R_max
    real(kind=8), parameter :: kb = 1.380649D-23
    real(kind=8), parameter :: el_charge = 1.602176634D-19
    real(kind=8), parameter :: permittivity = 8.85418782D-12
    real(kind=8), parameter :: pi = acos(-1.d0) 
    integer :: i, j, k, m, HS_flag, rand_particle_index, rdf_index, number_of_RDF_samples, rdf_local, widom_count
    integer :: widom_frequency, widom_counter, accepted_move, rej, charge_rej, charge_acc, clust_rej, clust_acc, clust_att
    real(kind=8) :: u_new, u_old, rand_particle, tx, ty, tz, rand_x, rand_y, rand_z, dx, dy, dz, r_squared
    real(kind=8) :: L_half, interaction_energy, deltaU, c_interaction_energy, xi, rel_deltaU, dist_scalar
    real(kind=8) :: total_V, ref_pp, ref_pn, ref_nn, r, r_upper, r_lower, dV, norm_ref_pp, norm_ref_pn, norm_ref_nn
    real(kind=8) :: Np, Nn, Q_diff, exp_Q_diff_sum, mu_ex_uncorr, K_widom, U_correction, net_charge, mu_correction, widom_net_charge
    real(kind=8) :: mu_ex, C_constant, sigma2, iL_half, eps_bulk, eps_contact, eps_delta, r_dist, lb_ir
    character(len=50) :: coordinate_input, RDF_output, coordinate_output, energy_output
    real(kind=8) :: cttx, ctty, cttz, net_u_charge, rdm_decision
    real(kind=8) :: frac_clust_move, R_clust_max, clust_step, cluster_energy_old, cluster_energy_new, delta_cluster

    real(kind=8), dimension(15000) :: xxx, yyy, zzz, charges
    real(kind=8), dimension(:), allocatable :: energy_data
    real(kind=8), dimension(15000, 4) :: coordinates
    real(kind=8), dimension(15000) :: gRpn, gRpp, gRnn, phi

    !Input file read
    open(unit=10, file='input_veps.txt', action='read', form='formatted')
    read(10,*)
    read(10,*) L
    read(10,*)
    read(10,*) Np, Nn
    read(10,*)
    read(10,*) sigma
    read(10,*)
    read(10,*) T
    read(10,*)
    read(10,*) eps_bulk, eps_contact, eps_delta
    read(10,*)
    read(10,*) step
    read(10,*)
    read(10,*) dR
    read(10,*)
    read(10,*) max_iterations
    read(10,*)
    read(10,*) max_subsimulations
    read(10,*)
    read(10,*) widom_frequency
    read(10,*)
    read(10,*) RDF_sample
    read(10,*)
    read(10,*) initialization
    read(10,*)
    read(10,*) frac_clust_move
    read(10,*)
    read(10,*) R_clust_max, clust_step
    read(10,*)
    read(10,*) coordinate_input
    read(10,*)
    read(10,*) RDF_output
    read(10,*)
    read(10,*) coordinate_output
    read(10,*) 
    read(10,*) energy_output
    close(unit=10)

    allocate(energy_data(max_subsimulations+1))

    NoP = int(Np) + int(Nn)
    NoPh = NoP / 2
    L_half = L / 2.d0
    beta = 1.d0 / (kb * T)
    R_max = L_half
    RDF_points = int(R_max / dR)
    !dr_intra = 0.01d0 * (sigma / 2.d0)
    !intra_RDF_points = int(sigma / dr_intra)
    sigma2 = sigma ** 2

    !Bulk lb
    lb = (1D10 * beta * el_charge ** 2.d0)/ (4.d0 * pi * permittivity * eps_bulk)

    accepted_move = 0
    charge_rej = 0
    charge_acc = 0
    clust_rej = 0
    clust_acc = 0
    rej = 0
    iL_half =  1.d0 / L_half
    widom_count = 0
    exp_Q_diff_sum = 0.d0
    K_widom = 2.d0 * (6.d0 * dlog(2.d0 + dsqrt(3.d0)) - pi)
    net_charge = Np - Nn 
    net_u_charge = 0.d0
    !We remove one cation, hence a net charge of -1
    widom_net_charge = net_charge - 1.d0
    C_constant = (lb * K_widom) / (8.d0 * L)
    U_correction = - C_constant * net_charge ** 2
    mu_correction = (C_constant *  widom_net_charge ** 2) + U_correction
    write(*,*) 'N+, N-: ', Np, Nn
    write(*,*) 'Net charge, reverse Widom charge: ', net_charge, widom_net_charge
    write(*,*) 'sigma: ', sigma
    write(*,*) 'L: ', L
    write(*,*) 'L_b (bulk): ', lb
    write(*,*) 'Epsilon (bulk), epsilon (contact): ', eps_bulk, eps_contact
    write(*,*) 'U correction: ', U_correction
    write(*,*) 'mu_ex correction: ', mu_correction
    write(*,*) 'Maximum trial configurations per macrostep: ', max_iterations
    write(*,*) 'Number of macrosteps: ', max_subsimulations

    !RDF initialization
    do i = 1, RDF_points 
        gRpn(i) = 0.d0
        gRpp(i) = 0.d0
        gRnn(i) = 0.d0
    end do
    number_of_RDF_samples = 0

    !Coordinate initialization
    if (initialization == 1) then
        write (*,*)
        write (*,*) '##################  Initial random configuration generation  ##################'
        call initialize(Np, Nn, L, xxx, yyy, zzz, charges, sigma)
        write (*,*) '###########  Generation of initial random configuration successful! ###########'
        write (*,*)
    else
        !Coordinate reading
        open(unit=10, file=coordinate_input, action='read', form="formatted")
        read(10, *)
        i = 1
        do while (i <= NoP) 
            read(10, *) coordinates(i, :)
            xxx(i) = coordinates(i, 1)
            yyy(i) = coordinates(i, 2)
            zzz(i) = coordinates(i, 3)
            charges(i) = coordinates(i, 4)
            i = i + 1
        end do
        close (unit=10)
        write(*,*) 'Coordinates read from file'
    end if

    !RNG setup
    call random_seed()

    !Check for HS overlap in initial configuration
    HS_flag = 0
    do i = 1, NoP
        do j = i + 1, NoP
            dx = dabs(xxx(j) - xxx(i))
            dx = dx - dint(dx * iL_half) * L
            dy = dabs(yyy(j) - yyy(i))
            dy = dy - dint(dy * iL_half) * L
            dz = dabs(zzz(j) - zzz(i))
            dz = dz - dint(dz * iL_half) * L
            
            r_squared = dx * dx + dy * dy + dz * dz

            if (r_squared .le. sigma2) then
                HS_flag = HS_flag + 1 
            end if
        end do
    end do

    !Initial interaction energy
    interaction_energy = 0.d0
    do k = 1, NoP
        do j = k + 1, NoP
            dx = dabs(xxx(j) - xxx(k))
            dx = dx - dint(dx * iL_half) * L
            dy = dabs(yyy(j) - yyy(k))
            dy = dy - dint(dy * iL_half) * L
            dz = dabs(zzz(j) - zzz(k))                
            dz = dz - dint(dz * iL_half) * L
        
            r_squared = dx * dx + dy * dy + dz * dz
            r_dist = dsqrt(r_squared)
            lb_ir = lb_vareps(r_dist) / r_dist
            interaction_energy = interaction_energy + lb_ir * (charges(k) * charges(j))
        end do
    end do

    write (*,*) 'Number of HS overlaps in initial configuration: ', HS_flag
    write (*,*) 'Initial interaction energy = ', interaction_energy + U_correction
    
    call flush
    !Main MC loop with m sub-simulations (i runs across MC configurations)
    do m = 1, max_subsimulations
        rdf_local = 0
        widom_counter = 0
        move_att = 0
        clust_att = 0
        do i = 1, max_iterations
        
            call random_number(rdm_decision)

            !Cluster move
            if ((rdm_decision) .le. (frac_clust_move)) then
                clust_att = clust_att + 1
                call cluster_move
            else
            !Normal move
            u_new = 0.d0
            u_old = 0.d0
            deltaU = 0.d0
            move_att = move_att + 1
            !Random trial move
            call random_number(rand_x)
            call random_number(rand_y)
            call random_number(rand_z)
            call random_number(rand_particle)
            rand_x = (rand_x - 0.5d0) * step
            rand_y = (rand_y - 0.5d0) * step
            rand_z = (rand_z - 0.5d0) * step
            rand_particle_index = 1 + nint(rand_particle * dfloat((NoP-1)))

            cttx = xxx(rand_particle_index)
            ctty = yyy(rand_particle_index)
            cttz = zzz(rand_particle_index)

            tx = xxx(rand_particle_index) + rand_x
            ty = yyy(rand_particle_index) + rand_y
            tz = zzz(rand_particle_index) + rand_z
            
            !PBC
            if (tx.gt.L_half) then
                tx = tx - L
            else if (tx.lt.-L_half) then
                tx = tx + L
            end if
                
            if (ty.gt.L_half) then
                ty = ty - L
            else if (ty.lt.-L_half) then
                ty = ty + L
            end if
                
            if (tz.gt.L_half) then
                tz = tz - L
            else if (tz.lt.-L_half) then
                tz = tz + L
            end if

            !New energy calculation
            do j = 1, NoP
                if (j == rand_particle_index) then
                    cycle
                end if
                dx = dabs(tx - xxx(j))
                dx = dx - dint(dx * iL_half) * L
                dy = dabs(ty - yyy(j))
                dy = dy - dint(dy * iL_half) * L
                dz = dabs(tz - zzz(j))
                dz = dz - dint(dz * iL_half) * L
        
                r_squared = dx * dx + dy * dy + dz * dz
                if (r_squared .le. sigma2) go to 89

                r_dist = dsqrt(r_squared)
                lb_ir = lb_vareps(r_dist) / r_dist
                
                u_new = u_new + lb_ir * (charges(rand_particle_index) * charges(j))
            end do

            !Old energy calculation
            do j = 1, NoP
                if (j == rand_particle_index) then
                    cycle
                end if
                dx = dabs(cttx - xxx(j))
                dx = dx - dint(dx * iL_half) * L
                dy = dabs(ctty - yyy(j))
                dy = dy - dint(dy * iL_half) * L
                dz = dabs(cttz - zzz(j))
                dz = dz - dint(dz * iL_half) * L
        
                r_squared = dx * dx + dy * dy + dz * dz
                r_dist = dsqrt(r_squared)
                lb_ir = lb_vareps(r_dist) / r_dist
                
                u_old = u_old + lb_ir * (charges(rand_particle_index) * charges(j))
            end do

            !Metropolis criteria
            deltaU = u_new - u_old
            call random_number(xi)
            if (xi < dexp(-1. * deltaU)) then
                xxx(rand_particle_index) = tx
                yyy(rand_particle_index) = ty
                zzz(rand_particle_index) = tz
                ! c_xxx(rand_particle_index) = ctx
                ! c_yyy(rand_particle_index) = cty
                ! c_zzz(rand_particle_index) = ctz
                interaction_energy = interaction_energy + deltaU
                accepted_move = accepted_move + 1
            else 
                go to 89
            end if
go to 90
89          rej = rej+1
90          continue
            end if

            !Widom sample
           if (mod(widom_counter, widom_frequency) == 0) then
                widom_count = widom_count + 1
                !Cations are in the first half of particles (currently)
                rand_particle_index = 1 + nint(rand_particle * (Np-1.))

                Q_diff = 0.d0
                do j = 1, NoP
                    if (j == rand_particle_index) then
                        cycle
                    end if
                    dx = dabs(xxx(rand_particle_index) - xxx(j))
                    dx = dx - dint(dx * iL_half) * L
                    dy = dabs(yyy(rand_particle_index) - yyy(j))
                    dy = dy - dint(dy * iL_half) * L
                    dz = dabs(zzz(rand_particle_index) - zzz(j))
                    dz = dz - dint(dz * iL_half) * L
            
                    r_squared = dx * dx + dy * dy + dz * dz
                    r_dist = dsqrt(r_squared)
                    lb_ir = lb_vareps(r_dist) / r_dist

                    Q_diff = Q_diff + lb_ir * (charges(rand_particle_index) * charges(j))

                end do
                exp_Q_diff_sum = exp_Q_diff_sum + dexp(Q_diff)
            end if 
            widom_counter = widom_counter + 1

            !RDF histogram count
            if (mod(rdf_local,RDF_sample) == 0) then
                number_of_RDF_samples = number_of_RDF_samples + 1
                do j = 1, NoP
                    tx = xxx(j)
                    ty = yyy(j)
                    tz = zzz(j)
                    ! ctx = c_xxx(j)
                    ! cty = c_yyy(j)
                    ! ctz = c_zzz(j)
!                    charge_distance = dist(ctx, cty, ctz, tx, ty, tz)
!                    if (charges(j) == -1.d0) then
!                        rdf_index = int(charge_distance / dr_intra) + 1
!                        gRintra(rdf_index) = gRintra(rdf_index) + 1.d0
!                   end if
                    do k = j + 1, NoP
                        dist_scalar = dist(tx, ty, tz, xxx(k), yyy(k), zzz(k))
                        if (dist_scalar < R_max) then
                            rdf_index = int(dist_scalar / dR) + 1
                            if (charges(j) == 1.d0) then
                                if (charges(k) == 1.d0) then
                                    gRpp(rdf_index) = gRpp(rdf_index) + 1.d0
                                else if (charges(k) == -1.d0) then
                                    gRpn(rdf_index) = gRpn(rdf_index) + 1.d0
                                end if
                            else if (charges(j) == -1.d0) then
                                if (charges(k) == -1.d0) then
                                    gRnn(rdf_index) = gRnn(rdf_index) + 1.d0
                                end if
                            end if
                        end if 
                    end do
                end do
            end if 
            rdf_local = rdf_local + 1
        end do

        !Sub-simulation statistics
        write (*,*)
        write (*,*) m, 'out of', max_subsimulations, 'macrosteps completed'

        c_interaction_energy = 0.d0
        do k = 1, NoP
            do j = k + 1, NoP
                dx = dabs(xxx(j) - xxx(k))
                dx = dx - dint(dx * iL_half) * L
                dy = dabs(yyy(j) - yyy(k))
                dy = dy - dint(dy * iL_half) * L
                dz = dabs(zzz(j) - zzz(k))
                dz = dz - dint(dz * iL_half) * L
                
                r_squared = dx * dx + dy * dy + dz * dz
                r_dist = dsqrt(r_squared)
                lb_ir = lb_vareps(r_dist) / r_dist
                  
                c_interaction_energy = c_interaction_energy + lb_ir * (charges(k) * charges(j))

            end do
        end do

!        c_charge_energy = sum(charge_energy_tracker)

        rel_deltaU = dabs(interaction_energy - c_interaction_energy) / interaction_energy
!        rel_charge_deltaU = dabs(charge_energy - c_charge_energy) / charge_energy
        write (*,*) 'Relative energy check is: ', rel_deltaU
!        write (*,*) 'Relative charge energy check is: ', rel_charge_deltaU
        write (*,*) 'Energy is: ', c_interaction_energy + U_correction
        write (*,*) 'Energy per particle is: ', (c_interaction_energy + U_correction) / float(NoP)
!        write(*,*)  'Charge energy is: ', c_charge_energy
        interaction_energy = c_interaction_energy
!        charge_energy = c_charge_energy
        write (*,*) 'Move acceptance percentage:', (accepted_move / dfloat(move_att)) * 100.d0
        write (*,*) 'Move rejection percentage:', (rej / dfloat(move_att)) * 100.d0
        write (*,*) 'Cluster move acceptance percentage:', (dfloat(clust_acc) / dfloat(clust_att)) * 100.d0
        write (*,*) 'Cluster move  rejection percentage:', (dfloat(clust_rej) / dfloat(clust_att)) * 100.d0
!        write (*,*) 'Charge perturbation acceptance percentage:', (charge_acc / dfloat(max_iterations)) * 100.d0
!        write (*,*) 'Charge perturbation rejection percentage:', (charge_rej / dfloat(max_iterations)) * 100.d0
        accepted_move = 0.d0
        move_att = 0
        rej = 0.d0
        charge_acc = 0.d0
        charge_rej = 0.d0
        clust_att = 0
        clust_acc = 0.0
        clust_rej = 0.0
        energy_data(m) = c_interaction_energy

        !Widom for the simulation
        mu_ex_uncorr = dlog(exp_Q_diff_sum / widom_count)
        mu_ex = mu_ex_uncorr + mu_correction
        write (*,*) 'mu_ex_uncorr:', mu_ex_uncorr
        write (*,*) 'mu_ex:', mu_ex

        call flush
    end do

    write (*,*)

    !Check for HS overlap in final simulation
    HS_flag = 0
    do i = 1, NoP
        do j = i + 1, NoP
            dx = dabs(xxx(j) - xxx(i))
            dx = dx - dint(dx * iL_half) * L
            dy = dabs(yyy(j) - yyy(i))
            dy = dy - dint(dy * iL_half) * L
            dz = dabs(zzz(j) - zzz(i))
            dz = dz - dint(dz * iL_half) * L
            
            r_squared = dx * dx + dy * dy + dz * dz

            if (r_squared .le. sigma2) then
                HS_flag = HS_flag + 1 
            end if
        end do
    end do
    if (HS_flag > 0) then
        write (*,*) 'There are', HS_flag, 'instances of HS overlap in the final configuration'
    else
        write (*,*) 'There is no HS overlap in the final configuration'
    end if

    !Final RDF calculation
    open(unit=10, file=RDF_output, status="replace", action="write", form="formatted")
    write (*,*) 'Number of RDF samples: ', number_of_RDF_samples
    total_V = L ** 3
    ref_pp = 0.5d0 * Np * (Np - 1.d0) * dfloat(number_of_RDF_samples) / total_V
    ref_pn =  Np * Nn * dfloat(number_of_RDF_samples)  / total_V
    ref_nn = 0.5d0 * Nn * (Nn - 1.d0) * dfloat(number_of_RDF_samples)  / total_V
    r = -dR * 0.5d0
    do i = 1, RDF_points
        r = r + dR
        r_upper = r + 0.5d0 * dR
        r_lower = r - 0.5d0 * dR
        dV = 4.d0 * pi * (r_upper ** 3 - r_lower ** 3) / 3.d0
        norm_ref_pp = ref_pp * dV
        norm_ref_pn = ref_pn * dV
        norm_ref_nn = ref_nn * dV
        gRpp(i) = gRpp(i) / norm_ref_pp
        gRpn(i) = gRpn(i) / norm_ref_pn
        gRnn(i) = gRnn(i) / norm_ref_nn
        write(10, '(4(1X, E16.8))')  r, gRpp(i), gRpn(i), gRnn(i)
    end do
    close(10)
    write(*,*) "RDFs written to file!"

    open(unit=10, file=coordinate_output, status="replace", action="write", form="formatted")
    write(10, *) "x, y, z, charge, chrage_x, charge_y, charge_z, charge_distance"
    do i = 1, NoP
        write(10, *) xxx(i), yyy(i), zzz(i), charges(i)
    end do
    close(10)
    write(*,*) "Coordinates written to file!"

    open(unit=10, file=energy_output, status="replace", action="write", form="formatted")
    write(10, *) "energy"
    do i = 1, max_subsimulations
        write(10, *) energy_data(i)
    end do
    close(10)
    write(*,*) "Energy data written to file!"

    !Potential visualisation
    open(unit=10, file="potential_plot.txt", status="replace", action="write", form="formatted")
    write(10, *) "r, phi(r), kappa(r)"
    r = -dR * 0.5d0
    do i=1, RDF_points
        r = r + dR
        if (r > sigma) then
            phi(i) = lb_vareps(r) / r
        end if
        write(10, *) r, phi(i)
    end do
    close(10)
    write(*,*) "Potential data written to file!"

contains   

    function dist(xx, yy, zz, x_coord, y_coord, z_coord)
        real(kind=8) :: delx, dely, delz, xx, yy, zz, x_coord, y_coord, z_coord, dist
        delx = abs(xx - x_coord)
        delx = delx-dint(delx * iL_half) * L
        dely = abs(yy - y_coord)
        dely = dely-dint(dely* iL_half) * L
        delz = abs(zz - z_coord)
        delz = delz-dint(delz* iL_half) * L
        dist = dsqrt(delx * delx + dely * dely + delz * delz)
        return
    end function

    function lb_vareps(r_value)
        real(kind=8) :: r_value, current_epsilon, lb_vareps

        if (r_value < sigma) then
            current_epsilon = eps_contact
        else if (r_value > (sigma + eps_delta)) then
            current_epsilon = eps_bulk
        else
            current_epsilon = eps_contact + (eps_bulk - eps_contact) * (r_value - sigma) / eps_delta !To je shady
        end if

        lb_vareps = (1D10 * beta * el_charge ** 2.d0)/ (4.d0 * pi * permittivity * current_epsilon)
        return
    end function

    subroutine initialize(Nop_p, Nop_n, length, x, y, z, c, diameter)
        real(8) ::  Nop_p, Nop_n, length, HS_step, u_new_HS
        real(8) , dimension(15000) :: x, y, z, c
        real(8) ::  x_rdm, y_rdm, z_rdm, rdm_particle, vec_x, vec_y, vec_z, half_length
        real(8) :: delx, dely, delz, r2, diameter
        integer :: Nop_tot, rdm_particle_index, HS_max_iterations

        HS_step = 20.d0
        HS_max_iterations = 1E6
        half_length = length / 2.d0

        Nop_tot = int(Nop_p) + int(Nop_n)

        do k = 1, int(Nop_p)
            call random_number(x(k))
            call random_number(y(k))
            call random_number(z(k))
            x(k) = (x(k) - 0.5d0) * length
            y(k) = (y(k) - 0.5d0) * length
            z(k) = (z(k) - 0.5d0) * length
            c(k) = 1.d0
        end do
    
        do k = int(Nop_p) + 1, Nop_tot + 1
            call random_number(x(k))
            call random_number(y(k))
            call random_number(z(k))
            x(k) = (x(k) - 0.5d0) * length
            y(k) = (y(k) - 0.5d0) * length
            z(k) = (z(k) - 0.5d0) * length
            c(k) = -1.d0
        end do

14        do i = 1, HS_max_iterations
            u_new_HS = 0.d0
    
13          call random_number(x_rdm)
            call random_number(y_rdm)
            call random_number(z_rdm)
            call random_number(rdm_particle)
            x_rdm = (x_rdm - 0.5d0) * HS_step
            y_rdm = (y_rdm - 0.5d0) * HS_step
            z_rdm = (z_rdm - 0.5d0) * HS_step
            rdm_particle_index = 1 + nint(rdm_particle * float(Nop_tot-1))
    
            vec_x = x(rdm_particle_index) + x_rdm
            vec_y = y(rdm_particle_index) + y_rdm
            vec_z = z(rdm_particle_index) + z_rdm
            
            if (vec_x.gt.half_length) then
                vec_x = vec_x - length
            else if (vec_x.lt.-half_length) then
                vec_x = vec_x + length
            end if
                
            if (vec_y.gt.half_length) then
                vec_y = vec_y - length
            else if (vec_y.lt.-half_length) then
                vec_y = vec_y + length
            end if
                
            if (vec_z.gt.half_length) then
                vec_z = vec_z - length
            else if (vec_z.lt.-half_length) then
                vec_z = vec_z + length
            end if
    
            do j = 1, Nop_tot
                if (j == rdm_particle_index) then
                    cycle
                end if
                delx = dabs(vec_x - x(j))
                delx = delx - dint(delx / half_length) * length
                dely = dabs(vec_y - y(j))
                dely = dely - dint(dely / half_length) * length
                delz = dabs(vec_z - z(j))
                delz = delz - dint(delz / half_length) * length
        
                r2 = delx * delx + dely * dely + delz * delz
    
                if (r2 <= (diameter ** 2)) go to 13
            end do
            
            x(rdm_particle_index) = vec_x
            y(rdm_particle_index) = vec_y
            z(rdm_particle_index) = vec_z
                
            if (modulo(i, HS_max_iterations / 10) == 0) then
                write(*, *) "Percentage completed: ", float(i) * 100. / HS_max_iterations, "%"
                HS_flag = 0
                do j = 1, Nop_tot
                    do k = j + 1, Nop_tot
                        delx = dabs(x(j) - x(k))
                        delx = delx - dint(delx / half_length) * length
                        dely = dabs(y(j) - y(k))
                        dely = dely - dint(dely / half_length) * length
                        delz = dabs(z(j) - z(k))
                        delz = delz - dint(delz / half_length) * length
                        
                        r2 = delx * delx + dely * dely + delz * delz
    
                        if (r2 <= (diameter ** 2)) then
                            HS_flag = HS_flag + 1 
                        end if
                    end do
                end do
                print *, 'Number of HS overalps in subsimulation:', HS_flag
            end if
        end do

        if (HS_flag > 0) then
            print *, 'HS overalps still present, restarting HS loop:' 
            go to 14
        end if 
    end subroutine

    subroutine cluster_move
        implicit None
        real(8) :: x_rdm, y_rdm, z_rdm
        real(8) :: rdm_particle, Nop_tot, rdm_radius, vec_x, vec_y, vec_z, old_x, old_y, old_z
        real(8) :: t_clust_x, t_clust_y, t_clust_z
        integer :: rdm_particle_index, index, f, p, g, index_counter
        real(8) , dimension(:), allocatable :: clust_x, clust_y, clust_z
        integer, dimension(:), allocatable :: Indexes, cluster_tracker

        Nop_tot = float(NoP)
        allocate(cluster_tracker(NoP+1))
        allocate(clust_x(NoP+1))
        allocate(clust_y(NoP+1))
        allocate(clust_z(NoP+1))
        allocate(Indexes(NoP+1))

        !Random displacement of the central particle
        call random_number(x_rdm)
        call random_number(y_rdm)
        call random_number(z_rdm)
        call random_number(rdm_particle)
        call random_number(rdm_radius)
        x_rdm = (x_rdm - 0.5d0) * clust_step
        y_rdm = (y_rdm - 0.5d0) * clust_step
        z_rdm = (z_rdm - 0.5d0) * clust_step
        rdm_radius = sigma + rdm_radius * (R_clust_max - sigma)
        rdm_particle_index = 1 + nint(rdm_particle * (Nop_tot-1))

        clust_x = 0.0
        clust_y = 0.0
        clust_z = 0.0

        cluster_tracker = 0
        Indexes = 0

        vec_x = xxx(rdm_particle_index) + x_rdm
        vec_y = yyy(rdm_particle_index) + x_rdm
        vec_z = zzz(rdm_particle_index) + x_rdm
        old_x = xxx(rdm_particle_index)
        old_y = yyy(rdm_particle_index)
        old_z = zzz(rdm_particle_index)

        !PBC for cluster centre
        if (vec_x.gt.L_half) then
            vec_x = vec_x - L
        else if (vec_x.lt.-L_half) then
            vec_x = vec_x + L
        end if
            
        if (vec_y.gt.L_half) then
            vec_y = vec_y - L
        else if (vec_y.lt.-L_half) then
            vec_y = vec_y + L
        end if
            
        if (vec_z.gt.L_half) then
            vec_z = vec_z - L
        else if (vec_z.lt.-L_half) then
            vec_z = vec_z + L
        end if
        
        !Idenitfy particles in the vicinity of the cluster centre at the new position
        do f = 1, NoP
            if (f == rdm_particle_index) then
                cycle
            end if
            dx = dabs(vec_x - xxx(f))
            dx = dx - dint(dx * iL_half) * L
            dy = dabs(vec_y - yyy(f))
            dy = dy - dint(dy * iL_half) * L
            dz = dabs(vec_z - zzz(f))
            dz = dz - dint(dz * iL_half) * L
    
            r_squared = dx * dx + dy * dy + dz * dz
            if (r_squared .le. (rdm_radius ** 2.0)) go to 29
        end do

        !Identify particles in the cluster at the original position
        do f = 1, NoP
            dx = dabs(old_x - xxx(f))
            dx = dx - dint(dx * iL_half) * L
            dy = dabs(old_y - yyy(f))
            dy = dy - dint(dy * iL_half) * L
            dz = dabs(old_z - zzz(f))
            dz = dz - dint(dz * iL_half) * L
            r_squared = dx * dx + dy * dy + dz * dz
            if (r_squared .le. rdm_radius ** 2.0) then
                cluster_tracker(f) = 1
            end if
        end do

        index_counter = 1
        do f = 1, NoP
            if (cluster_tracker(f) == 1) then
                t_clust_x = xxx(f) + x_rdm
                t_clust_y = yyy(f) + y_rdm
                t_clust_z = zzz(f) + z_rdm

                if (t_clust_x.gt.L_half) then
                    t_clust_x = t_clust_x - L
                else if (t_clust_x.lt.-L_half) then
                    t_clust_x = t_clust_x + L
                end if
                    
                if (t_clust_y.gt.L_half) then
                    t_clust_y = t_clust_y - L
                else if (t_clust_y.lt.-L_half) then
                    t_clust_y = t_clust_y + L
                end if
                    
                if (t_clust_z.gt.L_half) then
                    t_clust_z = t_clust_z - L
                else if (t_clust_z.lt.-L_half) then
                    t_clust_z = t_clust_z + L
                end if
                clust_x(f) = t_clust_x
                clust_y(f) = t_clust_y
                clust_z(f) = t_clust_z
                
                Indexes(index_counter) = f
                index_counter = index_counter + 1
            end if
        end do

        cluster_energy_old = 0.0
        cluster_energy_new = 0.0
        delta_cluster = 0.0

        !Calculate cluster energy in original position
        do p = 1, index_counter - 1
            index = Indexes(p)
            do f = 1, NoP
                do g = 1, index_counter - 1
                    if (f == Indexes(g)) then
                        cycle
                    end if
                    dx = dabs(xxx(index) - xxx(f))
                    dx = dx - dint(dx * iL_half) * L
                    dy = dabs(yyy(index) - yyy(f))
                    dy = dy - dint(dy * iL_half) * L
                    dz = dabs(zzz(index) - zzz(f))
                    dz = dz - dint(dz * iL_half) * L
        
                    r_squared = dx * dx + dy * dy + dz * dz
                    r_dist = dsqrt(r_squared)
                    lb_ir = lb_vareps(r_dist) / r_dist
                        
                    cluster_energy_old = cluster_energy_old + lb_ir * (charges(f) * charges(index))
                end do
            end do
        end do

        !Calculate cluster energy in new position
        do p = 1, index_counter - 1
            index = Indexes(p)
            do f = 1, NoP
                do g = 1, index_counter - 1
                    if (f == Indexes(g)) then
                        cycle
                    end if
                    dx = dabs(clust_x(index) - xxx(f))
                    dx = dx - dint(dx * iL_half) * L
                    dy = dabs(clust_y(index) - yyy(f))
                    dy = dy - dint(dy * iL_half) * L
                    dz = dabs(clust_z(index) - zzz(f))
                    dz = dz - dint(dz * iL_half) * L
        
                    r_squared = dx * dx + dy * dy + dz * dz
                    if (r_squared .le. sigma2) go to 29
                    r_dist = dsqrt(r_squared)
                    lb_ir = lb_vareps(r_dist) / r_dist
                        
                    cluster_energy_new = cluster_energy_new + lb_ir * (charges(f) * charges(index))
                end do
            end do
        end do

        delta_cluster = cluster_energy_new - cluster_energy_old 

        call random_number(xi)
        if (xi < dexp(-1. * delta_cluster)) then
            interaction_energy = interaction_energy + delta_cluster
            clust_acc = clust_acc + 1
            do p = 1, index_counter - 1
                index = Indexes(p)
                xxx(index) = clust_x(index)
                yyy(index) = clust_y(index)
                zzz(index) = clust_z(index)
            end do
        else
            go to 29
        end if

go to 30
29      clust_rej = clust_rej + 1
30      continue
    end subroutine

end program ionic_fluid
