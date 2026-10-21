%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% General parameters
params.Linit = 5;  
params.N = 8000;                         % number of points
N = params.N;
params.E = 0.05;                        % Young's modulus 
Tfinal = 200;   


% Number of species
params.Nb_species = 2;


% Diffusion parameters
params.D1 = 1;                         
params.D2 = 40;                    
params.D3 = 10; 
params.chi01 = 0.01;                     % chemotaxis strength
params.chi02 = 0.0;
params.k1 = 1;                          % chemotaxis saturation parameter

% params.D = @(u,params) [ params.D1,     0,              -params.chi01*params.k1*u(1)/(params.k1 + u(3));
%                                  0,     params.D2,      -params.chi02*params.k1*u(2)/(params.k1 + u(3));
%                                  0,     0,               params.D3];

params.D = @(u,params) [ params.D1,     -params.chi01*params.k1*u(1)/(params.k1 + u(2));
                                 0,     params.D2];

% params.D = @(u,params) [params.D1];


% Penalty parameters
% params.wall_left  = -20;
% params.wall_right = 25;
params.wall_left  = -200;
params.wall_right = 205;

% Penalty stiffness
params.kwall = 0.8;


% Reaction parameters
params.a = 0.01;                        % reaction parameter
params.b = 1.1;                         % reaction parameter
params.rho1 = 0.05;
params.rho2 = 0.1;
params.beta = 0.1;
params.r = 0.1;
params.K = 1;

% params.Reaction = @(u,params) [
%     params.a - u(1,:) + u(1,:).^2 .* u(2,:);
%     params.b - u(1,:).^2 .* u(2,:);
%     params.rho1*u(1,:) + params.rho2*u(2,:) - params.beta*u(3,:)
% ];

params.Reaction = @(u,params) [
    params.a - u(1,:) + u(1,:).^2 .* u(2,:);
    params.b - u(1,:).^2 .* u(2,:)
    ];

% params.Reaction = @(u,params) [params.r .* u(1,:) .* (1 - u(1,:)./params.K)];


% Body force density parameter
params.s_vec = [0.0001, 0.00001];          % one traction coefficient s for each species


% Growth law
params.S0     = 0.049082;
params.sigmac = 0.05;
params.S = @(u, sigma, params) params.S0 .* ( (sigma >= 0) + ((sigma < 0) & (sigma > -params.sigmac)) .* (1 - abs(sigma) ./ params.sigmac) );

% Stress law
params.eta = 1;  
params.E1 = params.E/2;
params.E2 = params.E/2;
params.visco.a2 = 0;
params.visco.a1 = 0;
params.visco.a0 = 1;
params.visco.b2 = 0;
params.visco.b1 = 0;
params.visco.b0 = params.E;

% params.Sigma = @(F,G,params) params.E .* (F./G - 1);
params.EpsilonLaw = @(F,G,xi,U_species,params) F./G - 1;
params.EpsilonLawA = @(A,G,xi,U_species,params) A - 1;
params.SigmaLaw = @(epsilon,F,G,xi,U_species,params) (params.visco.b0/params.visco.a0).*epsilon;        % Only used for order 0

% Boundary condition : can be "fixed-fixed", "free-free", "spring-spring" or "traction-traction"
params.mechanics_bc = "free-free";      
% params.mechanics_bc = "fixed-fixed"; 
% params.mechanics_bc = "spring-spring"; 
% params.mechanics_bc = "contact with a rigid wall penalty";
params.ks = 0.02;                          % spring constant
% params.mechanics_bc = "traction-traction"; 
params.traction_left  = 2e-5;
params.traction_right = 2e-5;


% Initial conditions
U_species0 = zeros(params.Nb_species, N);

u1_star = params.a + params.b;
u2_star = params.b/(params.a + params.b)^2;
u3_star = params.b/(params.a + params.b)^2;

% u1_star = params.K;

eps0 = 1e-3;
U_species0(1,:) = u1_star + eps0*randn(1,N);
U_species0(2,:) = u2_star + eps0*randn(1,N);
% U_species0(3,:) = u3_star + eps0*randn(1,N);


% Plot parameters
s_str = strjoin(compose('%.2e', params.s_vec), ', ');

% line1 = sprintf(['L_0=%.2f, N=%d, D_1=%.2f, D_2=%.2f, D_3=%.2f  ' ...
%     'E=%.3f, s=[%s]'], ...
%     params.Linit, params.N, params.D1, params.D2, params.D3, ...
%     params.E, s_str);
% 
% line2 = sprintf('a=%.3f, b=%.3f, \\chi_{01}=%.3f, \\chi_{02}=%.3f, k_{G}=%.3f, \\sigma_{0}=%.5f', ...
%         params.a, params.b, params.chi01, params.chi02, params.kG, params.sigma_c);
 
% line2 = sprintf('a=%.3f, b=%.3f, \\chi_{01}=%.3f, \\chi_{02}=%.3f, k_{s}=%.3f', ...
%         params.a, params.b, params.chi01, params.chi02,  params.ks);


line1 = sprintf(['L_0=%.2f, N=%d, D_1=%.2f, D_2=%.2f ' ...
    'E=%.3f, s=[%s]'], ...
    params.Linit, params.N, params.D1, params.D2, ...
    params.E, s_str);

line2 = sprintf('a=%.3f, b=%.3f, \\chi_{0}=%.3f, k_{G}=%.3f, \\sigma_{0}=%.5f, k_{s}=%.3f', ...
        params.a, params.b, params.chi01,  params.kG, params.sigma_c,  params.ks);

% line2 = sprintf('a=%.3f, b=%.3f, k_{s}=%.3f', ...
%     params.a, params.b,  params.ks);

% line2 = sprintf('a=%.3f, b=%.3f, \\chi_{01}=%.3f, \\chi_{02}=%.3f, k_1=%.2f', ...
%     params.a, params.b, params.chi01, params.chi02, params.k1);

% line1 = sprintf(['L_0=%.2f, N=%d, D_1=%.2f, ' ...
%     'E=%.3f, s=[%s]'], ...
%     params.Linit, params.N, params.D1, ...
%     params.E, s_str);

% line2 = sprintf('r=%.3f, K=%.3f, k_G=%.3f', ...
%     params.r, params.K, params.kG);

% line2 = sprintf('r=%.3f, K=%.3f, k_G=%.3f, \\sigma_0=%.3f', ...
%     params.r, params.K, params.kG, params.sigma0);

% line2 = sprintf('r=%.3f, K=%.3f, \\sigma_0=%.3f', ...
%     params.r, params.K, params.sigma0);

% line2 = sprintf('r=%.3f, K=%.3f, \\sigma_0=%.3f, k_{s}=%.3f', ...
%     params.r, params.K, params.sigma0, params.ks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



params.L0 = params.Linit;                        % initial length of the domain
L0 = params.L0;

L_initial = params.L0;


% Space grid
X = linspace(0, L0, N)';
dx = L0/(N-1);

params.X = X;
params.dx = dx;

% Width of the smooth transition at first contact
params.penalty_eps = 0.1;
% params.penalty_eps = 1e-2;
params.penalty_c0  = 0.5;
params.penalty_delta = 0.5;

% Order of the model
if params.visco.a2 == 0 && params.visco.b2 == 0

    if params.visco.a1 == 0 && params.visco.b1 == 0
        params.visco.modelOrder = 0;   % algebraic law
    elseif params.visco.a1 ~= 0 
        params.visco.modelOrder = 11;   % q = a1*sigma - b1*epsilon
    else
        % a1 = 0 and b1 ~= 0:
        % Kelvin-Voigt / purely viscous type
        params.visco.modelOrder = 12;

    end

else
    params.visco.modelOrder = 2;       % sigma, epsilon, y, z
end


% Initialization
G0 = ones(N,1);
A_half0 = ones(N-1,1);
xi0 = zeros(N,1);
F0 = ones(N,1);
epsilon0 = params.EpsilonLaw(F0, G0, xi0, U_species0, params);

if params.visco.modelOrder == 0

    % a0*sigma = b0*epsilon
    if params.visco.a0 == 0
        error('a0 must be nonzero.');
    end
    U0 = [U_species0(:); G0; A_half0; xi0];

elseif params.visco.modelOrder == 11

    % q = a1*sigma - b1*epsilon
    G_half0 = 0.5 .* (G0(1:end-1) + G0(2:end));
    xi_half0 = 0.5 .* (xi0(1:end-1) + xi0(2:end));
    U_species_half0 = 0.5 .* (U_species0(:,1:end-1) + U_species0(:,2:end));
    epsilon_half0 = params.EpsilonLawA(A_half0, G_half0, xi_half0, U_species_half0, params);
    sigma_half0 = zeros(N-1,1);
    q_half0 = params.visco.a1 .* sigma_half0 - params.visco.b1 .* epsilon_half0;
    U0 = [U_species0(:); G0; A_half0; xi0; q_half0];

elseif params.visco.modelOrder == 12

    G_half0 = 0.5 .* (G0(1:end-1) + G0(2:end));
    xi_half0 = 0.5 .* (xi0(1:end-1) + xi0(2:end));
    U_species_half0 = 0.5 .* (U_species0(:,1:end-1) + U_species0(:,2:end));

    epsilon_half0 = params.EpsilonLawA(A_half0, G_half0, xi_half0, U_species_half0, params);

    sigma_half0 = zeros(N-1,1);
    U0 = [U_species0(:); G0; A_half0; xi0; sigma_half0; epsilon_half0 ];

elseif params.visco.modelOrder == 2

    G_half0 = 0.5 .* (G0(1:end-1) + G0(2:end));

    xi_half0 = 0.5 .* (xi0(1:end-1) + xi0(2:end));

    U_species_half0 = 0.5 .* (U_species0(:,1:end-1) + U_species0(:,2:end));

    epsilon_half0 = params.EpsilonLawA(A_half0, G_half0, xi_half0, U_species_half0, params);

    % Relaxed initial mechanical state
    sigma_half0 = zeros(N-1,1);

    % y = sigma_t and z = epsilon_t
    y_half0 = zeros(N-1,1);
    z_half0 = zeros(N-1,1);

    U0 = [U_species0(:); G0; A_half0; xi0; sigma_half0; epsilon_half0; y_half0; z_half0];

end


% Building of the mass matrix M
switch params.visco.modelOrder

    case 0
        % U = [species; G; A; xi]
    
        M = spdiags([ones(params.Nb_species*N,1); ...   % species
            ones(N,1); ...                    % G
            zeros(N-1,1); ...                   % A
            zeros(N,1)], ...                  % xi
            0, params.Nb_species*N + 3*N - 1, params.Nb_species*N + 3*N - 1);

    case 11
        % U = [species; G; A; xi; q]
        M = spdiags([ones(params.Nb_species*N,1);     % species
            ones(N,1);            % G
            zeros(N-1,1);           % A
            zeros(N,1);             % xi
            ones(N-1,1)             % q
            ], 0, params.Nb_species*N + 4*N - 2, params.Nb_species*N + 4*N -2);

    case 12
        % U = [species; G; A; xi; sigma; epsilon]
        M = spdiags([ones(params.Nb_species*N,1);     % species
            ones(N,1);            % G
            zeros(N-1,1);         % A
            zeros(N,1);           % xi
            zeros(N-1,1);             % sigma
            params.visco.b1*ones(N-1,1)       % b1*epsilon_t
            ], 0, params.Nb_species*N + 5*N - 3, params.Nb_species*N + 5*N - 3);

    case 2
        % U = [species; G; A; xi; sigma; epsilon; y; z]
        M = spdiags([ones(params.Nb_species*N,1);...            % Cells
            ones(N,1); ...                                      % G
            zeros(N-1,1); ...                                     % A
            zeros(N,1); ...                                     % \xi
            ones(N-1,1); ...                                      % sigma
            ones(N-1,1); ...                                      % epsilon 
            zeros(N-1,1); ...                                     % epsilon = E(...)
            -params.visco.b2 * ones(N-1,1)], ...                  % z_t
            0, ...                                                                  % Diagonal
            params.Nb_species*N + 7*N - 5, ...                                          % Size
            params.Nb_species*N + 7*N - 5);                                             % Size
        
        M(params.Nb_species*N + 6*N - 4 + (1:N-1), params.Nb_species*N + 5*N - 3 + (1:N-1)) = params.visco.a2*speye(N-1);         % Non - diagonal term
end

% Solver options
JPattern = jacobian_sparsity_3modes(params.N, params.Nb_species, params.visco.modelOrder);
options = odeset('Mass', M, 'RelTol', 1e-6, 'AbsTol', 1e-8, 'JPattern', JPattern,'MaxStep',0.1);
% options = odeset('Mass', M, 'RelTol', 1e-3, 'AbsTol', 1e-5);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resolution with penalty contact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.xi_ref = zeros(N,1);

switch params.visco.modelOrder

    case 0
        System = @(t,U) PHI_order0(U,params);

    case 11
        System = @(t,U) PHI_order1(U,params);

    case 12
        System = @(t,U) PHI_order12(U,params);

    case 2
        System = @(t,U) PHI_order2(U,params);
end

% outputDt = 0.1;
outputDt = 0.5;

Tspan = 0:outputDt:Tfinal;

if Tspan(end) < Tfinal
    Tspan = [Tspan,Tfinal];
end


[t,U] = ode15s(System, Tspan, U0, options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U_species = U(:,1:params.Nb_species*N);

G = U(:, params.Nb_species*N + (1:N));

A_half = U(:,params.Nb_species*N + N + (1:N-1));

xi = U(:, params.Nb_species*N + 2*N - 1 + (1:N));

% Without remeshing, xi is already the total displacement
xi_total = xi;

% Deformation gradient:
% F = dx/dX = 1 + dxi/dX
Nt = length(t);

F = zeros(Nt,N);
F_half_all = zeros(Nt,N-1);


for k = 1:Nt

    G_half_k = 0.5 .* (G(k,1:end-1)' + G(k,2:end)');

    F_half_k = A_half(k,:)' .* G_half_k;

    F_k = zeros(N,1);

    F_k(1) = F_half_k(1);
    F_k(N) = F_half_k(end);

    F_k(2:N-1) = 0.5 .* (F_half_k(1:end-1) + F_half_k(2:end));

    F(k,:) = F_k';
    F_half_all(k,:) = F_half_k';

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check monotonicity of the physical mapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The physical mapping is x(X,t) = X + xi(X,t)
% Its derivative is dx/dX = F so the mapping is strictly increasing if F > 0

minFhalf_each_time = min(F_half_all,[],2);

first_nonmonotone_idx = find(minFhalf_each_time <= 0,1,'first');
if ~isempty(first_nonmonotone_idx)

    [minF_value,interval_idx] = min( ...
        F_half_all(first_nonmonotone_idx,:));

    error(['The physical mapping is no longer monotone. ' ...
        'First detected at t = %.6e between nodes %d and %d, ' ...
        'with F_{i+1/2} = %.6e.'], ...
        t(first_nonmonotone_idx), ...
        interval_idx, ...
        interval_idx+1, ...
        minF_value);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stress and strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch params.visco.modelOrder

    case 0

    epsilon = zeros(Nt,N);
    sigma   = zeros(Nt,N);

    epsilon_half_all = zeros(Nt,N-1);
    sigma_half_all   = zeros(Nt,N-1);

    for k = 1:Nt

        G_half_k = 0.5 .* (G(k,1:end-1)' + G(k,2:end)');

        xi_half_k = 0.5 .* (xi(k,1:end-1)' + xi(k,2:end)');

        U_species_k = reshape(U_species(k,:), params.Nb_species, N);

        U_species_half_k = 0.5 .* (U_species_k(:,1:end-1) + U_species_k(:,2:end));

        A_half_k = A_half(k,:)';

        F_half_k = A_half_k .* G_half_k;

        epsilon_half_k = params.EpsilonLawA(A_half_k, G_half_k, xi_half_k, U_species_half_k, params);

        sigma_half_k = params.SigmaLaw(epsilon_half_k, F_half_k, G_half_k, xi_half_k, U_species_half_k, params);

        epsilon_half_all(k,:) = epsilon_half_k';
        sigma_half_all(k,:)   = sigma_half_k';

        epsilon(k,:) = InterfaceToNode(epsilon_half_k)';

        sigma(k,:) = InterfaceToNode(sigma_half_k)';

    end

    q = [];
    y = [];
    z = [];

    case 11
    
        offset_A = params.Nb_species*N + N;
        offset_xi = offset_A + (N-1);
        offset_q = offset_xi + N;
    
        q_half_all = U(:,offset_q + (1:N-1));
    
        epsilon_half_all = zeros(Nt,N-1);
        sigma_half_all = zeros(Nt,N-1);
    
        epsilon = zeros(Nt,N);
        sigma = zeros(Nt,N);
        q = zeros(Nt,N);
    
        for k = 1:Nt
    
            U_species_k = reshape(U_species(k,:), params.Nb_species, N);
    
            G_half_k = 0.5 .* (G(k,1:end-1)' + G(k,2:end)');
    
            xi_half_k = 0.5 .* (xi(k,1:end-1)' + xi(k,2:end)');
    
            U_species_half_k = 0.5 .* (U_species_k(:,1:end-1) + U_species_k(:,2:end));
    
            A_half_k = A_half(k,:)';
    
            epsilon_half_k = params.EpsilonLawA(A_half_k, G_half_k, xi_half_k, U_species_half_k, params);
    
            q_half_k = q_half_all(k,:)';
    
            sigma_half_k = (q_half_k + params.visco.b1 .* epsilon_half_k) ./ params.visco.a1;
    
            epsilon_half_all(k,:) = epsilon_half_k';
            sigma_half_all(k,:) = sigma_half_k';
    
            epsilon(k,:) = InterfaceToNode(epsilon_half_k)';
            sigma(k,:) = InterfaceToNode(sigma_half_k)';
            q(k,:) = InterfaceToNode(q_half_k)';
    
        end
    
        y = [];
        z = [];

    case 12

        offset_A = params.Nb_species*N + N;
        offset_xi = offset_A + (N-1);
        offset_sigma = offset_xi + N;
        offset_epsilon = offset_sigma + (N-1);
    
        sigma_half_all = U(:,offset_sigma + (1:N-1));
    
        epsilon_half_all = U(:,offset_epsilon + (1:N-1));
    
        sigma = zeros(Nt,N);
        epsilon = zeros(Nt,N);

    for k = 1:Nt

        sigma(k,:) = InterfaceToNode(sigma_half_all(k,:)')';

        epsilon(k,:) = InterfaceToNode(epsilon_half_all(k,:)')';

    end

    q = [];
    y = [];
    z = [];

    case 2

        offset_A = params.Nb_species*N + N;
        offset_xi = offset_A + (N-1);
        offset_sigma = offset_xi + N;
        offset_epsilon = offset_sigma + (N-1);
        offset_y = offset_epsilon + (N-1);
        offset_z = offset_y + (N-1);
    
        sigma_half_all = U(:,offset_sigma + (1:N-1));
    
        epsilon_half_all = U(:,offset_epsilon + (1:N-1));
    
        y_half_all = U(:,offset_y + (1:N-1));
    
        z_half_all = U(:,offset_z + (1:N-1));
    
        sigma = zeros(Nt,N);
        epsilon = zeros(Nt,N);
        y = zeros(Nt,N);
        z = zeros(Nt,N);
    
        for k = 1:Nt
    
            sigma(k,:) = InterfaceToNode(sigma_half_all(k,:)')';
    
            epsilon(k,:) = InterfaceToNode(epsilon_half_all(k,:)')';
    
            y(k,:) = InterfaceToNode(y_half_all(k,:)')';
    
            z(k,:) = InterfaceToNode(z_half_all(k,:)')';
    
        end
    
        q = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruct the physical domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Physical position:
% x(X,t) = X + xi(X,t)
x_phys = X' + xi;



% Physical domain length
L_phys = x_phys(:,end) - x_phys(:,1);

% Centre the physical domain at each time for plotting
x_center = 0.5 .* (x_phys(:,1) + x_phys(:,end));
x_plot = x_phys - x_center;

% Initial centred material coordinates
X_initial_centered = X - 0.5*(X(1) + X(end));

% Final displacement measured between centred configurations
displacement_total_final = x_plot(end,:)' - X_initial_centered;


%% Reconstruct physical species densities

Nt = length(t);
u_plot = zeros(Nt, params.Nb_species, N);

for k = 1:Nt
    U_species_k = reshape(U_species(k,:), params.Nb_species, N);
    u_plot(k,:,:) = U_species_k ./ F(k,:);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Réduction de la quantité de donnée envoyée à surf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Downsampling for graphics 
space_step = 1;
time_step  = 2;

idx_x = 1:space_step:N;
idx_t = 1:time_step:length(t);

% Make sure the final point is included
if idx_x(end) ~= N
    idx_x = [idx_x,N];
end

if idx_t(end) ~= length(t)
    idx_t = [idx_t,length(t)];
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(t,L_phys,'LineWidth',2);
xlabel('t');
ylabel('L_{phys}(t)');
title("Physical domain length with boundary conditions " + params.mechanics_bc);
subtitle({line1; line2}, 'Interpreter','tex');


for i = 1:params.Nb_species

    ui = squeeze(u_plot(:,i,:));

    figure;
    % surf(x_plot,t,ui,'EdgeColor','none');
    surf(x_plot(idx_t,idx_x), t(idx_t), ui(idx_t,idx_x), 'EdgeColor','none');
    view(2);
    colorbar;
    xlabel('$x$', 'Interpreter','latex', 'FontSize',30);
    ylabel('time', 'FontName','Cambria', 'FontSize',30);
    ylim([0 Tfinal])
    % title("Species " + i + " density with boundary conditions " + params.mechanics_bc);
    % subtitle({line1; line2}, 'Interpreter','tex');
    ax = gca;
    ax.FontName = 'Cambria';
    ax.FontSize = 30;

    figure;
    plot(x_plot(end,:),ui(end,:),'LineWidth',2);
    xlabel('x');
    ylabel("u_" + i);
    title("Final species " + i + " density with boundary conditions " + params.mechanics_bc);
    subtitle({line1; line2}, 'Interpreter','tex');

end


figure;
% surf(x_plot,t,G,'EdgeColor','none');
surf(x_plot(idx_t,idx_x), t(idx_t), G(idx_t,idx_x), 'EdgeColor','none');
view(2);
colorbar;
xlabel('x');
ylabel('t');
ylim([0 Tfinal])
title("Growth G with boundary conditions " + params.mechanics_bc);
subtitle({line1; line2}, 'Interpreter','tex');

xi_center = 0.5 .* (xi_total(:,1) + xi_total(:,end));
xi_total_centered = xi_total - xi_center;
figure;
% surf(x_plot,t,xi_total_centered,'EdgeColor','none');
surf(x_plot(idx_t,idx_x), t(idx_t), xi_total_centered(idx_t,idx_x), 'EdgeColor','none');
view(2);
colorbar;
xlabel('x');
ylabel('t');
ylim([0 Tfinal])
title("Displacement \xi with boundary conditions " + params.mechanics_bc);
subtitle({line1; line2}, 'Interpreter','tex');


figure;
plot(x_plot(end,:),G(end,:),'LineWidth',2);
xlabel('x');
ylabel('G');
title("Final growth with boundary conditions " + params.mechanics_bc);
subtitle({line1; line2}, 'Interpreter','tex');


figure;
plot(X_initial_centered, displacement_total_final,'LineWidth',2);
xlabel('initial position X');
ylabel('total displacement');
title("Final displacement from initial position with boundary conditions " + params.mechanics_bc);
subtitle({line1; line2}, 'Interpreter','tex');


figure;
% surf(x_plot,t,sigma,'EdgeColor','none');
surf(x_plot(idx_t,idx_x), t(idx_t), sigma(idx_t,idx_x), 'EdgeColor','none');
view(2);
colorbar;
xlabel('$x$', 'Interpreter','latex', 'FontSize',30);
ylabel('time', 'FontName','Cambria', 'FontSize',30);
ylim([0 Tfinal])
% title("Mechanical stress \sigma with boundary conditions " + params.mechanics_bc);
% subtitle({line1; line2}, 'Interpreter','tex');
ax = gca;
ax.FontName = 'Cambria';
ax.FontSize = 30;

