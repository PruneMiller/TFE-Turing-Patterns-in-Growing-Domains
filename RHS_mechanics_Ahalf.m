function [rhs_A,rhs_xi] = RHS_mechanics_Ahalf(U_species,G,A_half,xi,params)

    N = params.N;
    dx = params.dx;
    ks = params.ks;
    s_vec = params.s_vec(:);

    % Total displacement
    xi_total = params.xi_ref + xi;

    % Interface quantities

    G_half = 0.5 .* (G(1:end-1) + G(2:end));

    xi_half = 0.5 .* (xi(1:end-1) + xi(2:end));

    U_species_half = 0.5 .* (U_species(:,1:end-1) + U_species(:,2:end));

    % F = A G
    F_half = A_half .* G_half;

    % Elastic strain at interfaces
    epsilon_half = params.EpsilonLawA(A_half, G_half, xi_half, U_species_half, params);

    % Stress at interfaces
    sigma_half = params.SigmaLaw(epsilon_half, F_half, G_half, xi_half, U_species_half, params);

   
    % Nodal F and physical densities

    F_node = zeros(N,1);

    F_node(1) = F_half(1);
    F_node(N) = F_half(end);

    F_node(2:N-1) = 0.5 .* (F_half(1:end-1) + F_half(2:end));

    u_species = U_species ./ F_node';

    force_density = s_vec' * u_species;
    % force_density = 0.01 + s_vec' * u_species;
    force_density = force_density';



    % Geometric compatibility: N-1 equations

    rhs_A = F_half - 1 - diff(xi)./dx;


    % Mechanical equilibrium: N equations

    rhs_xi = zeros(N,1);

    rhs_xi(2:N-1) = (sigma_half(2:end) - sigma_half(1:end-1))./dx - force_density(2:N-1) .* xi_total(2:N-1);


    % Boundary conditions

    switch params.mechanics_bc

        case "fixed-fixed"

            rhs_xi(1) = xi_total(1);
            rhs_xi(N) = xi_total(N);

        case "free-free"

            if all(abs(params.s_vec) < 1e-14)

                % No body-force anchoring: remove rigid translation with one gauge condition
                rhs_xi(1) = xi_total(1);
                rhs_xi(N) = sigma_half(end);

            else

                % Body-force term removes translation invariance
                rhs_xi(1) = sigma_half(1);
                rhs_xi(N) = sigma_half(end);
                % rhs_xi(1) = sigma_half(1) - 0.008;
                % rhs_xi(N) = sigma_half(end) - 0.008;

            end

        case "spring-spring"

            rhs_xi(1) = sigma_half(1) - ks*xi_total(1);

            rhs_xi(N) = sigma_half(end) + ks*xi_total(N);

        case "traction-traction"

            rhs_xi(1) = sigma_half(1) - params.traction_left;

            rhs_xi(N) = sigma_half(end) - params.traction_right;

        case "contact with a rigid wall penalty"

             % Physical positions of the tissue boundaries
            x_left  = params.X(1)   + xi_total(1);
            x_right = params.X(end) + xi_total(end);

            % Signed penetrations:
            % z < 0 : tissue is inside the admissible region
            % z > 0 : tissue has crossed the wall
            z_left  = params.wall_left - x_left;
            z_right = x_right - params.wall_right;


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fonction de pénalisation n°1 : p(z) = k max(0,z)

            % penetration_left  = max(0,z_left);
            % penetration_right = max(0,z_right);
            % 
            % pressure_left  = params.kwall .* penetration_left;
            % pressure_right = params.kwall .* penetration_right;

            % ne marche pas car max n'est pas différentiable en 0
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fonction de pénalisation n°2 : p(z) = k ​max(0,z)²

            % penetration_left  = max(0,z_left);
            % penetration_right = max(0,z_right);
            % 
            % pressure_left  = params.kwall .* penetration_left.^3;
            % pressure_right = params.kwall .* penetration_right.^3;

            % --> marche jusqu'à t = 163
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fonction de pénalisation n°3 : p(z) = k/2 (z + sqrt(z² + eps²))

            eps_penalty = params.penalty_eps;
            % 
            % % Smooth approximation of max(0,z):
            % %
            % % positive_part(z)
            % % = 1/2 * (z + sqrt(z^2 + eps^2))
            % %
            % % Its derivative at z = 0 is 1/2, so the contact
            % % contribution is visible in the Newton Jacobian.
            % % 
            % penetration_left = 0.5 .* (z_left + sqrt(z_left.^2 + eps_penalty.^2));
            % 
            % penetration_right = 0.5 .* (z_right + sqrt(z_right.^2 + eps_penalty.^2));
            % 
            % % Penalty pressures
            % pressure_left = params.kwall .* penetration_left;
            % pressure_right = params.kwall .* penetration_right;

            % --> marche jusqu'à t = 160
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fonction de pénalisation n°4 : p(z) = k/2 (z + sqrt(z² + delta (1-z²/c0²)²eps²)) pour |z|<c0 et p(z) = k ​max(0,z)² pour |z|>c0

            penetration_left = LocalSmoothPositivePart(z_left, params.penalty_eps, params.penalty_c0, params.penalty_delta);

            penetration_right = LocalSmoothPositivePart(z_right, params.penalty_eps, params.penalty_c0, params.penalty_delta);

            pressure_left  = params.kwall .* penetration_left;
            pressure_right = params.kwall .* penetration_right;

            % --> marche  jusqu'à t = 160 avec : params.penalty_eps = 1e-1;
            params.penalty_c0  = 0.5;
            params.penalty_delta = 0.5;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



            % Boundary mechanical equilibrium
            rhs_xi(1) = sigma_half(1) + pressure_left;

            rhs_xi(N) = sigma_half(end) + pressure_right;


    end

end