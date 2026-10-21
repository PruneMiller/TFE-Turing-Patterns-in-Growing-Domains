function phi = PHI_order0(U,params)

    N = params.N;
    Nb_species = params.Nb_species;

    % -------------------------------------------------------------
    % Extract variables
    % -------------------------------------------------------------

    U_species = reshape(U(1:Nb_species*N), Nb_species, N);

    G = U(Nb_species*N + (1:N));

    A_half = U(Nb_species*N + N + (1:N-1));

    xi = U(Nb_species*N + 2*N - 1 + (1:N));
    
    
    % -------------------------------------------------------------
    % Reconstruct F
    % -------------------------------------------------------------

    G_half = 0.5 .* (G(1:end-1) + G(2:end));

    F_half = A_half .* G_half; 

    F = zeros(N,1);

    F(1) = F_half(1);
    F(N) = F_half(end);

    F(2:N-1) = 0.5 .* (F_half(1:end-1) + F_half(2:end));


    % Physical species densities
    u_species = U_species ./ F';

    % -------------------------------------------------------------
    % Nodal strain and stress
    % -------------------------------------------------------------
    
    xi_half = 0.5 .* (xi(1:end-1) + xi(2:end));

    U_species_half = 0.5 .* (U_species(:,1:end-1) + U_species(:,2:end));

    epsilon_half = params.EpsilonLawA(A_half, G_half, xi_half, U_species_half, params);

    sigma_half = params.SigmaLaw(epsilon_half, F_half, G_half, xi_half, U_species_half, params);

    sigma_node = InterfaceToNode(sigma_half);

    % -------------------------------------------------------------
    % Species equations
    % -------------------------------------------------------------

     reaction = params.Reaction(u_species,params);

    phi_U_species = F' .* reaction + DiffusionAhalf(U_species, F, F_half, params);

    phi_U_species = phi_U_species(:);

    % -------------------------------------------------------------
    % Growth equation
    % -------------------------------------------------------------

    S = params.S(u_species, sigma_node, params);

    phi_G = S .* G;

    % -------------------------------------------------------------
    % Algebraic mechanical equations
    % -------------------------------------------------------------

    [phi_A,phi_xi] = RHS_mechanics_Ahalf(U_species, G, A_half, xi, params);

    phi = [phi_U_species; phi_G; phi_A; phi_xi];

end