function phi = PHI_order1(U, params)
    % PHI = RHS

    % Parameters
    N = params.N;
    dx = params.dx;
    Nb_species = params.Nb_species;

    % Indices
    % U = [species; G; A_half; xi; q_half]

    offset_species = 0;
    offset_G       = Nb_species*N;
    offset_A       = offset_G + N;
    offset_xi      = offset_A + (N-1);
    offset_q       = offset_xi + N;


    % Extract variables

    U_species = reshape(U(offset_species + (1:Nb_species*N)), Nb_species, N);

    G = U(offset_G + (1:N));

    A_half = U(offset_A + (1:N-1));

    xi = U(offset_xi + (1:N));

    q_half = U(offset_q + (1:N-1));


    % Interface quantities
    G_half = 0.5 .* (G(1:end-1) + G(2:end));

    xi_half = 0.5 .* (xi(1:end-1) + xi(2:end));

    U_species_half = 0.5 .* (U_species(:,1:end-1) + U_species(:,2:end));

    % F = A G at interfaces
    F_half = A_half .* G_half;

    % Nodal F
    F = InterfaceToNode(F_half);

    % Physical densities at nodes
    u_species = U_species ./ F';


    % Strain and stress at interfaces

    epsilon_half = params.EpsilonLawA(A_half, G_half, xi_half, U_species_half, params);

    % sigma = (q + b1*epsilon)/a1
    sigma_half = (q_half + params.visco.b1 .* epsilon_half) ./ params.visco.a1;

    sigma_node = InterfaceToNode(sigma_half);


    % Species equations
    reaction = params.Reaction(u_species,params);

    phi_U_species = F' .* reaction + DiffusionAhalf(U_species,F,F_half,params);

    phi_U_species = phi_U_species(:);


    % Growth equation
    growth_rate = params.S(u_species,sigma_node,params);

    phi_G = growth_rate .* G;


    % Geometric compatibility
    % 0 = A_half*G_half - 1 - diff(xi)/dx

    phi_A = F_half - 1 - diff(xi)./dx;


    % Mechanical equilibrium
    % 0 = sigma_X - force_density*xi

    xi_total = params.xi_ref + xi;

    force_density = params.s_vec(:)' * u_species;
    force_density = force_density';

    phi_xi = zeros(N,1);

    phi_xi(2:N-1) = (sigma_half(2:end) - sigma_half(1:end-1))./dx - force_density(2:N-1) .* xi_total(2:N-1);

    % Boundary conditions
    switch params.mechanics_bc

        case "fixed-fixed"

            phi_xi(1) = xi_total(1);
            phi_xi(N) = xi_total(N);

        case "free-free"

            if all(abs(params.s_vec) == 0)

                % Gauge condition plus zero stress
                phi_xi(1) = xi_total(1);
                phi_xi(N) = sigma_half(end);

            else

                phi_xi(1) = sigma_half(1);
                phi_xi(N) = sigma_half(end);

            end

        case "spring-spring"

            phi_xi(1) = sigma_half(1) - params.ks .* xi_total(1);

            phi_xi(N) = sigma_half(end) + params.ks .* xi_total(N);

        case "traction-traction"

            phi_xi(1) = sigma_half(1) - params.traction_left;

            phi_xi(N) = sigma_half(end) - params.traction_right;
        
        case "contact with a rigid wall"

            x_left  = params.X(1)   + xi_total(1);
            x_right = params.X(end) + xi_total(end);

            penetration_left  = max(0, params.wall_left - x_left);
            penetration_right = max(0, x_right - params.wall_right);

            phi_xi(1) = sigma_half(1) + params.kwall * penetration_left;

            phi_xi(N) = sigma_half(end) + params.kwall * penetration_right;

    end


    % q evolution
    % q_t = -(a0/a1)q + (b0 - a0*b1/a1)epsilon

    phi_q = -(params.visco.a0/params.visco.a1) .* q_half + (params.visco.b0 - params.visco.a0*params.visco.b1/params.visco.a1) .* epsilon_half;


    phi = [phi_U_species; phi_G; phi_A; phi_xi; phi_q];

end
