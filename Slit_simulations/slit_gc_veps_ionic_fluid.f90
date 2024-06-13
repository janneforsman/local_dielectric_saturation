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
    real(kind=8) :: ref_pp, r, r_upper, r_lower, dV, norm_ref_pp
    real(kind=8) :: Np, Nn, Q_diff, exp_Q_diff_sum, K_widom, U_correction, net_charge, mu_correction, widom_net_charge
    real(kind=8) :: C_constant, sigma2, iL_half, eps_bulk, eps_contact, eps_delta, r_dist, lb_ir
    character(len=50) :: coordinate_input, RDF_output, coordinate_output, energy_output
    real(kind=8) :: cttx, ctty, cttz, net_u_charge
    real(kind=8) :: frac_clust_move, R_clust_max, clust_step, cluster_energy_old, cluster_energy_new, delta_cluster
    integer :: z_dist_int
    real(kind=8) :: slit_h, half_slit_h, single_ion_ext_pot, rj, dif_z, half_dR, z_low, z_high, u_ext_new
    real(kind=8) :: frac_GC_move, chem_pot, volume, rdm_decision, u_ext_c, rel_deltaU_ext, u_ext
    real(kind=8) :: u_ext_old, delta_u_ext
    integer :: creation_rej, creation_acc, GC_cre_att, GC_dest_att, dest_rej, dest_acc

    real(kind=8), dimension(15000) :: xxx, yyy, zzz, charges, ext_pot
    real(kind=8), dimension(:), allocatable :: energy_data
    real(kind=8), dimension(15000, 4) :: coordinates
    real(kind=8), dimension(15000, 3) :: n_z
    real(kind=8), dimension(15000) :: nz_p, nz_n, phi, z_axis, charge_density

    !Input file read
    open(unit=10, file='input_GC_slit_veps.txt', action='read', form='formatted')
    read(10,*)
    read(10,*) L
    read(10,*)
    read(10,*) slit_h
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
    read(10,*) frac_GC_move, chem_pot
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

    NoPh = NoP / 2
    half_dR = dR / 2.d0
    L_half = L / 2.d0
    beta = 1.d0 / (kb * T)
    R_max = slit_h
    RDF_points = int(slit_h / dR)
    sigma2 = sigma ** 2
    volume = L * L * (slit_h - sigma)

    !Bulk lb
    lb = (1D10 * beta * el_charge ** 2.d0)/ (4.d0 * pi * permittivity * eps_bulk)

    accepted_move = 0
    charge_rej = 0
    charge_acc = 0
    clust_rej = 0
    clust_acc = 0
    u_ext = 0.0
    rej = 0
    GC_cre_att = 0
    GC_dest_att = 0
    creation_acc = 0
    creation_rej = 0
    dest_acc = 0
    dest_rej = 0
    iL_half =  1.d0 / L_half
    half_slit_h = slit_h / 2.d0
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
    
    number_of_RDF_samples = 0

    !n(z) initialization
    if (initialization == 1) then
        do i = 1, RDF_points 
            nz_p(i) = 0.d0
            nz_n(i) = 0.d0
            charge_density(i) = 0.d0
        end do
    else
        open(unit=162, file=RDF_output, action='read')
        i = 1
        do while (i <= RDF_points) 
            read(162, *) n_z(i, :)
            z_axis(i) = n_z(i, 1)
            nz_p(i) = n_z(i, 2)
            nz_n(i) = n_z(i, 3)
            charge_density(i) = nz_p(i) - nz_n(i)
            i = i + 1
        end do
        close (unit=162)
        write(*,*) 'n(z) read from file'

        open(unit=123, file = 'GC_NoP.txt', action='read')
        read(123,*)
        read(123,*) Np, Nn
        close (unit=123)
    end if

    !Charge density simetrization
    j = RDF_points+1
    do i = 1,RDF_points/2
        j = j-1
        charge_density(i) = 0.5d0*(charge_density(i)+charge_density(j))
    end do
    j = RDF_points/2+1
    do i = RDF_points/2+1,RDF_points
        j = j-1
        charge_density(i) = charge_density(j)
    end do

    NoP = int(Np) + int(Nn)

    !Coordinate initialization
    if (initialization == 1) then
        write (*,*)
        write (*,*) '##################  Initial random configuration generation  ##################'
        call initialize(Np, Nn, L, slit_h, xxx, yyy, zzz, charges, sigma)
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
            dz = dz - dint(dz * half_slit_h) * (slit_h - sigma)
            
            r_squared = dx * dx + dy * dy + dz * dz

            if (r_squared .le. sigma2) then
                HS_flag = HS_flag + 1 
            end if
        end do
    end do

    !Initial external potential vector
    open(unit=123, file='ext_pot.txt', status="replace", action="write", form="formatted")
    r = -dR * 0.5d0 - half_slit_h
    do i = 1, RDF_points
        r = r + dR
        single_ion_ext_pot = 0.d0
        rj = -dR * 0.5d0 - half_slit_h
        do j = 1, RDF_points
            rj = rj + dR
            dif_z = dabs(r - rj)
            single_ion_ext_pot = single_ion_ext_pot + charge_density(j) * (-2.d0 * pi * dif_z - Phiw(dif_z))
        end do
        ext_pot(i) = single_ion_ext_pot * dR * lb
        write(123, '(4(1X, E16.8))') r, ext_pot(i), charge_density(i)
    end do
    close(123)

    !Initial interaction and external energy
    interaction_energy = 0.d0
    do k = 1, NoP
        do j = k + 1, NoP
            dx = dabs(xxx(j) - xxx(k))
            dx = dx - dint(dx * iL_half) * L
            dy = dabs(yyy(j) - yyy(k))
            dy = dy - dint(dy * iL_half) * L
            dz = dabs(zzz(j) - zzz(k))                
        
            r_squared = dx * dx + dy * dy + dz * dz
            r_dist = dsqrt(r_squared)
            lb_ir = lb_vareps(r_dist) / r_dist
            interaction_energy = interaction_energy + lb_ir * (charges(k) * charges(j))
        end do

        !External potential contribution
        z_dist_int = int((zzz(k) + half_slit_h + half_dR) / dR)
        z_low = dfloat(z_dist_int) * dR - half_dR - half_slit_h
        z_high = z_low + dR

        if (charges(k) == 1.d0) then
            u_ext=u_ext+((z_high - zzz(k))*ext_pot(z_dist_int) + (-z_low + zzz(k))*ext_pot(z_dist_int + 1)) / dR
        else
            u_ext=u_ext-((z_high - zzz(k))*ext_pot(z_dist_int) + (-z_low + zzz(k))*ext_pot(z_dist_int + 1)) / dR
        end if
    end do

    write(*,*) 'N+, N-: ', Np, Nn
    write(*,*) 'Net charge, reverse Widom charge: ', net_charge, widom_net_charge
    write(*,*) 'sigma: ', sigma
    write(*,*) 'L: ', L
    write(*,*) 'H: ', slit_h
    write(*,*) 'L_b (bulk): ', lb
    write(*,*) 'Epsilon (bulk), epsilon (contact): ', eps_bulk, eps_contact
    write(*,*) 'U correction: ', U_correction
    write(*,*) 'mu_ex correction: ', mu_correction
    write(*,*) 'Maximum trial configurations per macrostep: ', max_iterations
    write(*,*) 'Number of macrosteps: ', max_subsimulations
    write (*,*) 'Number of HS overlaps in initial configuration: ', HS_flag
    write (*,*) 'Initial interaction energy = ', interaction_energy + U_correction
    write (*,*) 'Initial external potential energy = ', u_ext
    write (*,*) 'Initial total energy = ', interaction_energy + U_correction + u_ext
    
    call flush
    !Main MC loop with m sub-simulations (i runs across MC configurations)
    do m = 1, max_subsimulations
        rdf_local = 0
        widom_counter = 0
        move_att = 0
        do i = 1, max_iterations

            call random_number(rdm_decision)

            !GC move
            if (rdm_decision .le. frac_GC_move) then
                call random_number(xi)
                if (xi .le. 0.5) then
                    call GC_creation_move
                    GC_cre_att = GC_cre_att + 1
                else
                    call GC_destruction_move
                    GC_dest_att = GC_dest_att + 1
                end if
            end if

            !Cluster move
            if ((rdm_decision + frac_GC_move) .le. (frac_clust_move + frac_GC_move)) then
                clust_att = clust_att + 1
                call cluster_move
            else
                !Translation move
                u_new = 0.d0
                u_old = 0.d0
                deltaU = 0.d0
                u_ext_new = 0.0
                u_ext_old = 0.0
                delta_u_ext = 0.0
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

                dz = abs(tz)
                if (dz > (half_slit_h - sigma/2.0)) go to 89

                !External potential contribution
                z_dist_int = int((tz + half_slit_h + half_dR) / dR)
                z_low = dfloat(z_dist_int) * dR - half_dR - half_slit_h
                z_high = z_low + dR

                if (charges(rand_particle_index) == 1.d0) then
                    u_ext_new = ((z_high-tz)*ext_pot(z_dist_int)+(-z_low+tz)*ext_pot(z_dist_int+1))/dR
                else
                    u_ext_new = -((z_high-tz)*ext_pot(z_dist_int)+(-z_low+tz)*ext_pot(z_dist_int+1))/dR
                end if

                z_dist_int = int((cttz + half_slit_h + half_dR) / dR)
                z_low = dfloat(z_dist_int) * dR - half_dR - half_slit_h
                z_high = z_low + dR

                if (charges(rand_particle_index) == 1.d0) then
                    u_ext_old=((z_high-cttz)*ext_pot(z_dist_int)+(-z_low+cttz)*ext_pot(z_dist_int+1))/dR
                else
                    u_ext_old=-((z_high-cttz)*ext_pot(z_dist_int)+(-z_low+cttz)*ext_pot(z_dist_int+1))/dR
                end if

                delta_u_ext = u_ext_new - u_ext_old

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
            
                    r_squared = dx * dx + dy * dy + dz * dz
                    r_dist = dsqrt(r_squared)
                    lb_ir = lb_vareps(r_dist) / r_dist
                    
                    u_old = u_old + lb_ir * (charges(rand_particle_index) * charges(j))
                end do

                !Metropolis criteria
                deltaU = u_new - u_old
                call random_number(xi)
                if (xi < dexp(-1. * (deltaU + delta_u_ext))) then
                    xxx(rand_particle_index) = tx
                    yyy(rand_particle_index) = ty
                    zzz(rand_particle_index) = tz
                    interaction_energy = interaction_energy + deltaU
                    u_ext = u_ext + delta_u_ext
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
                    dist_scalar = tz + half_slit_h
                        if (dist_scalar < R_max) then
                            rdf_index = int(dist_scalar / dR) + 1
                            if (charges(j) == 1.d0) then
                                nz_p(rdf_index) = nz_p(rdf_index) + 1.d0
                            else if (charges(j) == -1.d0) then
                                nz_n(rdf_index) = nz_n(rdf_index) + 1.d0
                            end if
                        end if 
                end do
            end if 
            rdf_local = rdf_local + 1
        end do

        !Sub-simulation statistics
        write (*,*)
        write (*,*)
        write (*,*) m, 'out of', max_subsimulations, 'macrosteps completed'

        c_interaction_energy = 0.d0
        u_ext_c = 0.0
        do k = 1, NoP
            do j = k + 1, NoP
                dx = dabs(xxx(j) - xxx(k))
                dx = dx - dint(dx * iL_half) * L
                dy = dabs(yyy(j) - yyy(k))
                dy = dy - dint(dy * iL_half) * L
                dz = dabs(zzz(j) - zzz(k))
                
                r_squared = dx * dx + dy * dy + dz * dz
                r_dist = dsqrt(r_squared)
                lb_ir = lb_vareps(r_dist) / r_dist
                  
                c_interaction_energy = c_interaction_energy + lb_ir * (charges(k) * charges(j))
            end do
            !External potential contribution
            z_dist_int = int((zzz(k) + half_slit_h + half_dR) / dR)
            z_low = dfloat(z_dist_int) * dR - half_dR - half_slit_h
            z_high = z_low + dR

            if (charges(k) == 1.d0) then
                u_ext_c=u_ext_c+((z_high-zzz(k))*ext_pot(z_dist_int)+(-z_low+zzz(k))*ext_pot(z_dist_int+1))/dR
            else
                u_ext_c=u_ext_c-((z_high-zzz(k))*ext_pot(z_dist_int)+(-z_low+zzz(k))*ext_pot(z_dist_int+1))/dR
            end if
        end do

        rel_deltaU = dabs(interaction_energy - c_interaction_energy) / interaction_energy
        rel_deltaU_ext = dabs(u_ext - u_ext_c) / u_ext
        write (*,*)
        write (*,*) 'Relative interaction energy check is: ', rel_deltaU
        write (*,*) 'Relative external field energy check is: ', rel_deltaU_ext
        write (*,*) 'Interaction energy is: ', c_interaction_energy + U_correction
        write (*,*) 'Interaction energy per particle is: ', (c_interaction_energy + U_correction) / float(NoP)
        write (*,*) 'External field energy is: ', u_ext
        interaction_energy = c_interaction_energy
        u_ext = u_ext_c
        write (*,*) 'Move acceptance percentage:', (accepted_move / dfloat(move_att)) * 100.d0
        write (*,*) 'Move rejection percentage:', (rej / dfloat(move_att)) * 100.d0
        write (*,*) 'Cluster move acceptance percentage:', (dfloat(clust_acc) / dfloat(clust_att)) * 100.d0
        write (*,*) 'Cluster move rejection percentage:', (dfloat(clust_rej) / dfloat(clust_att)) * 100.d0
        write (*,*) 'GC insertion acceptance percentage:', (dfloat(creation_acc) / dfloat(GC_cre_att)) * 100.d0
        write (*,*) 'GC insertion rejection percentage:', (dfloat(creation_rej) / dfloat(GC_cre_att)) * 100.d0
        write (*,*) 'GC deletion acceptance percentage:', (dfloat(dest_acc) / dfloat(GC_dest_att)) * 100.d0
        write (*,*) 'GC deletion rejection percentage:', (dfloat(dest_rej) / dfloat(GC_dest_att)) * 100.d0
        write (*,*) 'Current number of particles:', NoP
        write (*,*) 'Current number of cations, anions:', int(Np), int(Nn)
        accepted_move = 0
        move_att = 0
        rej = 0
        clust_att = 0
        clust_acc = 0.0
        clust_rej = 0.0
        GC_cre_att = 0
        GC_dest_att = 0
        creation_acc = 0
        creation_rej = 0
        dest_acc = 0
        dest_rej = 0
        
        energy_data(m) = c_interaction_energy

        !Widom for the simulation
        !mu_ex_uncorr = dlog(exp_Q_diff_sum / widom_count)
        !mu_ex = mu_ex_uncorr + mu_correction
        !write (*,*) 'mu_ex_uncorr:', mu_ex_uncorr
        !write (*,*) 'mu_ex:', mu_ex

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
            dz = dz - dint(dz * half_slit_h) * slit_h
            
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

 !Final n(z) calculation
    open(unit=10, file=RDF_output, status="replace", action="write", form="formatted")
    write (*,*) 'Number of RDF samples: ', number_of_RDF_samples
    ref_pp =  dfloat(number_of_RDF_samples) 
    r = -dR * 0.5d0  - half_slit_h
    do i = 1, RDF_points
        r = r + dR
        r_upper = r + 0.5d0 * dR
        r_lower = r - 0.5d0 * dR
        dV = L ** 2 * dR
        norm_ref_pp = ref_pp * dV
        nz_p(i) = nz_p(i) / norm_ref_pp
        nz_n(i) = nz_n(i) / norm_ref_pp
        write(10, '(3(1X, E16.8))') r, nz_p(i), nz_n(i)
    end do
    close(10)
    write(*,*) "N(z) written to file!"

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

    !GC NoP
    open(unit=123, file="GC_NoP.txt", status="replace", action="write", form="formatted")
    write(123, *) "Number of cations, anions for \mu = ", chem_pot
    write(123, *) Np, Nn
    close(123)
    write(*,*) "Number of particles written to file!"

contains   

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

    subroutine initialize(Nop_p, Nop_n, length, slit_distance, x, y, z, c, diameter)
        real(8) ::  Nop_p, Nop_n, length, slit_distance, HS_step, u_new_HS
        real(8) , dimension(15000) :: x, y, z, c
        real(8) ::  x_rdm, y_rdm, z_rdm, rdm_particle, vec_x, vec_y, vec_z, half_length, half_slit_distance
        real(8) :: delx, dely, delz, r2, diameter
        integer :: Nop_tot, rdm_particle_index, HS_max_iterations

        HS_step = 20.d0
        HS_max_iterations = 1E5
        half_length = length / 2.d0
        half_slit_distance = slit_distance / 2.d0

        Nop_tot = int(Nop_p) + int(Nop_n)

        do k = 1, int(Nop_p)
            call random_number(x(k))
            call random_number(y(k))
            call random_number(z(k))
            x(k) = (x(k) - 0.5d0) * length
            y(k) = (y(k) - 0.5d0) * length
            z(k) = (z(k) - 0.5d0) * (slit_distance - diameter)
            c(k) = 1.d0
        end do
    
        do k = int(Nop_p) + 1, Nop_tot + 1
            call random_number(x(k))
            call random_number(y(k))
            call random_number(z(k))
            x(k) = (x(k) - 0.5d0) * length
            y(k) = (y(k) - 0.5d0) * length
            z(k) = (z(k) - 0.5d0) * (slit_distance - diameter)
            c(k) = -1.d0
        end do

        do i = 1, HS_max_iterations
            u_new_HS = 0.d0
    
13          call random_number(x_rdm)
            call random_number(y_rdm)
            call random_number(z_rdm)
            call random_number(rdm_particle)
            x_rdm = (x_rdm - 0.5d0) * HS_step
            y_rdm = (y_rdm - 0.5d0) * HS_step
            z_rdm = (z_rdm - 0.5d0) * HS_step
            rdm_particle_index = 1 + nint(rdm_particle * Nop_tot)
    
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
                
            if (vec_z.gt.half_slit_distance) then
                vec_z = vec_z - slit_distance
            else if (vec_z.lt.-half_slit_distance) then
                vec_z = vec_z + slit_distance
            end if

            if (abs(vec_z) > ((slit_distance - diameter)/ 2.d0)) go to 13
    
            do j = 1, Nop_tot
                if (j == rdm_particle_index) then
                    cycle
                end if
                delx = dabs(vec_x - x(j))
                delx = delx - dint(delx / half_length) * length
                dely = dabs(vec_y - y(j))
                dely = dely - dint(dely / half_length) * length
                delz = dabs(vec_z - z(j))
                delz = delz - dint(delz / half_slit_distance) * slit_distance
        
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
                        delz = delz - dint(delz / half_slit_distance) * slit_distance
                        
                        r2 = delx * delx + dely * dely + delz * delz
    
                        if (r2 <= (diameter ** 2)) then
                            HS_flag = HS_flag + 1 
                        end if
                    end do
                end do
                print *, 'Number of HS overalps in subsimulation:', HS_flag
                HS_flag = 0
            end if
        end do
    end subroutine

    subroutine cluster_move
        implicit None
        real(8) :: x_rdm, y_rdm, z_rdm
        real(8) :: rdm_particle, Nop_tot, rdm_radius, vec_x, vec_y, vec_z, old_x, old_y, old_z
        real(8) :: t_clust_x, t_clust_y, t_clust_z, u_ext_cluster_new, u_ext_cluster_old
        integer :: rdm_particle_index, index, f, p, g, index_counter
        real(8) , dimension(:), allocatable :: clust_x, clust_y, clust_z
        integer, dimension(:), allocatable :: Indexes, cluster_tracker

        Nop_tot = float(NoP)
        allocate(cluster_tracker(NoP))
        allocate(clust_x(NoP))
        allocate(clust_y(NoP))
        allocate(clust_z(NoP))
        allocate(Indexes(NoP))

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
        rdm_particle_index = 1 + nint(rdm_particle * dfloat((NoP-1)))

        clust_x = 0.0
        clust_y = 0.0
        clust_z = 0.0

        cluster_tracker = 0
        Indexes = 0

        vec_x = xxx(rdm_particle_index) + x_rdm
        vec_y = yyy(rdm_particle_index) + y_rdm
        vec_z = zzz(rdm_particle_index) + z_rdm
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

                if (abs(t_clust_z) > ((slit_h - sigma)/ 2.d0)) go to 29
                    
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
        u_ext_cluster_new = 0.0
        u_ext_cluster_old = 0.0
        delta_u_ext = 0.0

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
        
                    r_squared = dx * dx + dy * dy + dz * dz
                    r_dist = dsqrt(r_squared)
                    lb_ir = lb_vareps(r_dist) / r_dist
                        
                    cluster_energy_old = cluster_energy_old + lb_ir * (charges(f) * charges(index))
                end do
            end do

            !External potential contribution
            z_dist_int = int((zzz(index) + half_slit_h + half_dR) / dR)
            z_low = dfloat(z_dist_int) * dR - half_dR - half_slit_h
            z_high = z_low + dR

            if (charges(index) == 1.0) then
u_ext_cluster_old=u_ext_cluster_old+((z_high-zzz(index))*ext_pot(z_dist_int)+(-z_low+zzz(index))*ext_pot(z_dist_int+1))/dR
            else
u_ext_cluster_old=u_ext_cluster_old-((z_high-zzz(index))*ext_pot(z_dist_int)+(-z_low+zzz(index))*ext_pot(z_dist_int+1))/dR
            end if
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
        
                    r_squared = dx * dx + dy * dy + dz * dz
                    if (r_squared .le. sigma2) go to 29
                    r_dist = dsqrt(r_squared)
                    lb_ir = lb_vareps(r_dist) / r_dist
                        
                    cluster_energy_new = cluster_energy_new + lb_ir * (charges(f) * charges(index))
                end do
            end do

            !External potential contribution
            z_dist_int = int((clust_z(index) + half_slit_h + half_dR) / dR)
            z_low = dfloat(z_dist_int) * dR - half_dR - half_slit_h
            z_high = z_low + dR

            if (charges(index) == 1.d0) then
u_ext_cluster_new=u_ext_cluster_new+((z_high-clust_z(index))*ext_pot(z_dist_int)+(-z_low+clust_z(index))*ext_pot(z_dist_int+1))/dR
            else
u_ext_cluster_new=u_ext_cluster_new-((z_high-clust_z(index))*ext_pot(z_dist_int)+(-z_low+clust_z(index))*ext_pot(z_dist_int+1))/dR
            end if
        end do

        !Computing the differences in energy
        delta_cluster = cluster_energy_new - cluster_energy_old
        delta_u_ext = u_ext_cluster_new - u_ext_cluster_old

        !Metropolis energy cgeck
        call random_number(xi)
        if (xi < dexp(-1.0 * (delta_cluster + delta_u_ext))) then
            interaction_energy = interaction_energy + delta_cluster
            u_ext = u_ext + delta_u_ext
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

    function Phiw(z)
        real(8) :: z, Phiw, z2, L_half_2
        z2 = z*z
        L_half_2 = L_half*L_half
        Phiw = 8.d0*L_half*dlog((dsqrt(2.d0*L_half_2+z2)+L_half)/(dsqrt(L_half_2+z2)))
        Phiw = Phiw -2.d0*z*(dasin((L_half_2**2-z2*z2-2.d0*L_half_2*z2)/(L_half_2+z2)**2) + pi / 2.d0)
        return
    end  

    subroutine GC_creation_move
        implicit None
        real(8) :: xp_rdm, yp_rdm, zp_rdm, xn_rdm, yn_rdm, zn_rdm, creation_energy, u_ext_creation_p
        real(8) :: u_ext_creation_n, u_ext_creation, GC_cre_factor, probability
        integer :: f

        creation_energy = 0.0
        u_ext_creation_p = 0.0
        u_ext_creation_n = 0.0

        !Random insertion of ion pair
        call random_number(xp_rdm)
        call random_number(yp_rdm)
        call random_number(zp_rdm)
        call random_number(xn_rdm)
        call random_number(yn_rdm)
        call random_number(zn_rdm)
        xp_rdm = (xp_rdm - 0.5d0) * L
        yp_rdm = (yp_rdm - 0.5d0) * L
        zp_rdm = (zp_rdm - 0.5d0) * (slit_h - sigma)
        xn_rdm = (xn_rdm - 0.5d0) * L
        yn_rdm = (yn_rdm - 0.5d0) * L
        zn_rdm = (zn_rdm - 0.5d0) * (slit_h - sigma)

        !New energy calculation
        do f = 1, NoP
            !Cation
            dx = dabs(xp_rdm - xxx(f))
            dx = dx - dint(dx * iL_half) * L
            dy = dabs(yp_rdm - yyy(f))
            dy = dy - dint(dy * iL_half) * L
            dz = dabs(zp_rdm - zzz(f))
            r_squared = dx * dx + dy * dy + dz * dz
            if (r_squared .le. sigma2) go to 49
            r_dist = dsqrt(r_squared)
            lb_ir = lb_vareps(r_dist) / r_dist
            creation_energy = creation_energy + lb_ir * charges(f)

            !Anion
            dx = dabs(xn_rdm - xxx(f))
            dx = dx - dint(dx * iL_half) * L
            dy = dabs(yn_rdm - yyy(f))
            dy = dy - dint(dy * iL_half) * L
            dz = dabs(zn_rdm - zzz(f))
            r_squared = dx * dx + dy * dy + dz * dz
            if (r_squared .le. sigma2) go to 49
            r_dist = dsqrt(r_squared)
            lb_ir = lb_vareps(r_dist) / r_dist
            creation_energy = creation_energy - lb_ir * charges(f)
        end do

        !Ion pair interaction energy
        dx = dabs(xp_rdm - xn_rdm)
        dx = dx - dint(dx * iL_half) * L
        dy = dabs(yp_rdm - yn_rdm)
        dy = dy - dint(dy * iL_half) * L
        dz = dabs(zp_rdm - zn_rdm)
        r_squared = dx * dx + dy * dy + dz * dz
        if (r_squared .le. sigma2) go to 49
        r_dist = dsqrt(r_squared)
        lb_ir = lb_vareps(r_dist) / r_dist
        creation_energy = creation_energy - lb_ir

        !External potential contribution
        !Cation
        z_dist_int = int((zp_rdm + half_slit_h + half_dR) / dR)
        z_low = dfloat(z_dist_int) * dR - half_dR - half_slit_h
        z_high = z_low + dR
        u_ext_creation_p = ((z_high-zp_rdm)*ext_pot(z_dist_int)+(-z_low+zp_rdm)*ext_pot(z_dist_int+1))/dR

        !Anion
        z_dist_int = int((zn_rdm + half_slit_h + half_dR) / dR)
        z_low = dfloat(z_dist_int) * dR - half_dR - half_slit_h
        z_high = z_low + dR
        u_ext_creation_n = -((z_high-zn_rdm)*ext_pot(z_dist_int)+(-z_low+zn_rdm)*ext_pot(z_dist_int+1))/dR

        u_ext_creation = u_ext_creation_p + u_ext_creation_n

        delta_u_ext = u_ext_creation

        !Metropolis energy check
        GC_cre_factor=(chem_pot-creation_energy-delta_u_ext)+dlog(volume**2.0/((Np+1.0)*(Nn+1.0)))
        if (GC_cre_factor .lt. -100.0) go to 49
        probability = dexp(GC_cre_factor)
        call random_number(xi)
        if (xi .lt. probability)  then 
            creation_acc = creation_acc + 1
         !   write(*,*) interaction_energy, creation_energy
            interaction_energy = interaction_energy + creation_energy
            u_ext = u_ext + u_ext_creation
            xxx(NoP+1) = xp_rdm
            yyy(NoP+1) = yp_rdm
            zzz(NoP+1) = zp_rdm
            charges(NoP+1) = 1.0
            xxx(NoP+2) = xn_rdm
            yyy(NoP+2) = yn_rdm
            zzz(NoP+2) = zn_rdm
            charges(NoP+2) = -1.0
            NoP = NoP + 2
            Np = Np + 1.0
            Nn = Nn + 1.0
          !  write(*,*) interaction_energy, creation_energy
        else
            go to 49
        end if

go to 50
49      creation_rej = creation_rej + 1
50      continue
    end subroutine

    subroutine GC_destruction_move
        implicit None
        real(8) :: rdm_particle_p, rdm_particle_n, dest_energy, u_ext_destruction_p
        real(8) :: u_ext_destruction_n, u_ext_destruction, GC_dest_factor, probability
        integer :: f, rdm_p_index, rdm_n_index

        dest_energy = 0.0
        u_ext_destruction_p = 0.0
        u_ext_destruction_n = 0.0

        !Random deletion of ion pair
158        call random_number(rdm_particle_p)
159        call random_number(rdm_particle_n)

        rdm_p_index = 1 + nint(rdm_particle_p * dfloat((NoP-1)))
        if (charges(rdm_p_index) .ne. 1.0) go to 158
        rdm_n_index = 1 + nint(rdm_particle_n * dfloat((NoP-1)))
        if (charges(rdm_n_index) .ne. -1.0) go to 159

        if (abs(zzz(rdm_p_index)) > ((slit_h - sigma)/ 2.0)) go to 158
        if (abs(zzz(rdm_n_index)) > ((slit_h - sigma)/ 2.0)) go to 159

        !New energy calculation
        do f = 1, NoP
            if (f .eq. rdm_p_index) then
                cycle
            end if
            if (f .eq. rdm_n_index) then
                cycle
            end if

            !Cation
            dx = dabs(xxx(rdm_p_index) - xxx(f))
            dx = dx - dint(dx * iL_half) * L
            dy = dabs(yyy(rdm_p_index)  - yyy(f))
            dy = dy - dint(dy * iL_half) * L
            dz = dabs(zzz(rdm_p_index)  - zzz(f))
            r_squared = dx * dx + dy * dy + dz * dz
            r_dist = dsqrt(r_squared)
            lb_ir = lb_vareps(r_dist) / r_dist
            dest_energy = dest_energy + lb_ir * charges(f)

            !Anion
            dx = dabs(xxx(rdm_n_index) - xxx(f))
            dx = dx - dint(dx * iL_half) * L
            dy = dabs(yyy(rdm_n_index) - yyy(f))
            dy = dy - dint(dy * iL_half) * L
            dz = dabs(zzz(rdm_n_index) - zzz(f))
            r_squared = dx * dx + dy * dy + dz * dz
            r_dist = dsqrt(r_squared)
            lb_ir = lb_vareps(r_dist) / r_dist
            dest_energy = dest_energy - lb_ir * charges(f)
        end do

        !Ion pair interaction energy
        dx = dabs(xxx(rdm_p_index) - xxx(rdm_n_index))
        dx = dx - dint(dx * iL_half) * L
        dy = dabs(yyy(rdm_p_index) - yyy(rdm_n_index))
        dy = dy - dint(dy * iL_half) * L
        dz = dabs(zzz(rdm_p_index) - zzz(rdm_n_index))
        r_squared = dx * dx + dy * dy + dz * dz
        r_dist = dsqrt(r_squared)
        lb_ir = lb_vareps(r_dist) / r_dist
        dest_energy = dest_energy - lb_ir

        !External potential contribution
        !Cation
        z_dist_int = int((zzz(rdm_p_index) + half_slit_h + half_dR) / dR)
        z_low = dfloat(z_dist_int) * dR - half_dR - half_slit_h
        z_high = z_low + dR
u_ext_destruction_p = ((z_high-zzz(rdm_p_index))*ext_pot(z_dist_int)+(-z_low+zzz(rdm_p_index))*ext_pot(z_dist_int+1))/dR

        !Anion
        z_dist_int = int((zzz(rdm_n_index) + half_slit_h + half_dR) / dR)
        z_low = dfloat(z_dist_int) * dR - half_dR - half_slit_h
        z_high = z_low + dR
u_ext_destruction_n = -((z_high-zzz(rdm_n_index))*ext_pot(z_dist_int)+(-z_low+zzz(rdm_n_index))*ext_pot(z_dist_int+1))/dR

        u_ext_destruction = u_ext_destruction_p + u_ext_destruction_n
        delta_u_ext = u_ext_destruction

        !Metropolis energy check
        GC_dest_factor=(-chem_pot+dest_energy+delta_u_ext)+dlog((Np*Nn)/(volume**2.0))
        if (GC_dest_factor .lt. -100.0) go to 59
        probability = dexp(GC_dest_factor)
        call random_number(xi)
        if (xi .lt. probability) then 
            if ((Np-1.0) .le. 1) go to 59 
            if ((Nn-1.0) .le. 1) go to 59 
            dest_acc = dest_acc + 1
            interaction_energy = interaction_energy - dest_energy
            u_ext = u_ext - u_ext_destruction

            ! write(*,*) probability, xi
            ! write(*,*) rdm_p_index, rdm_n_index, dest_energy, u_ext_destruction, GC_dest_factor
            ! stop
            
           ! write(*,*) rdm_p_index, rdm_n_index
            if (rdm_p_index < rdm_n_index) then
                !Cation deletion
                do f = 1, NoP-1
                    if (f < rdm_p_index) then
                        cycle
                    end if
                    xxx(f) = xxx(f + 1)
                    yyy(f) = yyy(f + 1)
                    zzz(f) = zzz(f + 1)
                    charges(f) = charges(f + 1)
                end do
                !Anion deletion
                do f = 1, NoP-2
                    if (f < (rdm_n_index-1)) then
                        cycle
                    end if
                    xxx(f) = xxx(f + 1)
                    yyy(f) = yyy(f + 1)
                    zzz(f) = zzz(f + 1)
                    charges(f) = charges(f + 1)
                end do
            else
                !Anion deletion
                do f = 1, NoP-1
                    if (f < rdm_n_index) then
                        cycle
                    end if
                    xxx(f) = xxx(f + 1)
                    yyy(f) = yyy(f + 1)
                    zzz(f) = zzz(f + 1)
                    charges(f) = charges(f + 1)
                end do
                !Cation deletion
                do f = 1, NoP-2
                    if (f < (rdm_p_index-1)) then
                        cycle
                    end if
                    xxx(f) = xxx(f + 1)
                    yyy(f) = yyy(f + 1)
                    zzz(f) = zzz(f + 1)
                    charges(f) = charges(f + 1)
                end do
            end if
            NoP = NoP - 2
            Np = Np - 1.0
            Nn = Nn - 1.0
        else
            go to 59
        end if
go to 60
59      dest_rej = dest_rej + 1
60      continue
    end subroutine

end program ionic_fluid
