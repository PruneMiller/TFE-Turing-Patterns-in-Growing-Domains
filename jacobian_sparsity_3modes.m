function S = jacobian_sparsity_3modes(N,Nb_species,modelOrder)
% JACOBIAN_SPARSITY_3MODES
% Sparsity pattern of dPHI/dU.
%
% modelOrder = 0:
%   U = [species; G; A_half; xi]
%
%   species : Nb_species*N variables
%   G       : N variables
%   A_half  : N-1 variables
%   xi      : N variables
%
% modelOrder = 11:
%   U = [species; G; xi; q]
%
% modelOrder = 12:
%   U = [species; G; xi; sigma; epsilon]
%
% modelOrder = 2:
%   U = [species; G; xi; sigma; epsilon; y; z]


    % Species are stored node by node because U_species(:) gives:
    %
    % [u1(X1),...,un(X1),u1(X2),...,un(X2),...]

    idx_species = @(i) (i-1)*Nb_species + (1:Nb_species);

    idx_G = @(i) Nb_species*N + i;


    switch modelOrder

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Model order 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        case 0

            % State:
            %
            % U = [species; G; A_half; xi]

            Nspecies = Nb_species*N;

            offset_G = Nspecies;

            offset_A = offset_G + N;

            offset_xi = offset_A + (N-1);

            Ntot = offset_xi + N;

            idx_A = @(k) offset_A + k;

            idx_xi = @(i) offset_xi + i;


            % Conservative estimate of the number of nonzero entries

            S = spalloc(Ntot, Ntot, 30*(Nb_species+4)^2*N);



            % Species and growth equations

            for i = 1:N

                rows_species = idx_species(i);

                row_G = idx_G(i);


                % Neighbouring nodal indices

                neigh1 = max(i-1,1):min(i+1,N);

                neigh2 = max(i-2,1):min(i+2,N);



                % Neighbouring interface indices
                %
                % Interface k lies between nodes k and k+1.
                %
                % A nodal F near node i depends on interfaces i-1 and i.
                % Diffusion may extend this dependence by one additional
                % node on each side.

                A_neigh_diff = max(i-2,1):min(i+1,N-1);

                A_neigh_local = max(i-1,1):min(i,N-1);


                % Species equations
                %
                % phi_species = F*Reaction(u) + Diffusion(U_species,F)
                %
                % The diffusion stencil involves species at i-1:i+1.
                %
                % Since
                %
                %     F_half = A_half .* G_half,
                %
                % the nodal F values entering nearby fluxes depend on:
                %
                %     G at approximately i-2:i+2,
                %     A_half at approximately i-2:i+1

                for j = neigh1

                    S(rows_species, idx_species(j)) = 1;

                end


                for j = neigh2

                    S(rows_species, idx_G(j)) = 1;

                end


                for k = A_neigh_diff

                    S(rows_species, idx_A(k)) = 1;

                end


                % The species equations do not depend directly on xi
                % because F is now reconstructed from A_half and G.


                % Growth equation
                %
                % phi_G = S(u,sigma)*G
                %
                % At node i, u and sigma may depend on:
                %
                % - local species;
                % - nearby G through F_node;
                % - interfaces A_{i-1/2} and A_{i+1/2};
                % - local xi if EpsilonLaw or SigmaLaw depends explicitly
                %   on xi.

                S(row_G, idx_species(i)) = 1;


                for j = neigh1

                    S(row_G, idx_G(j)) = 1;

                end


                for k = A_neigh_local

                    S(row_G, idx_A(k)) = 1;

                end


                S(row_G, idx_xi(i)) = 1;

            end


            % Compatibility equations
            %
            % For k = 1,...,N-1:
            % A_half(k)*G_half(k) - 1 - (xi(k+1)-xi(k))/dx = 0
            %
            % The compatibility rows occupy the same positions as A_half.


            for k = 1:N-1

                row_A = idx_A(k);


                % Dependence on A_half(k)

                S(row_A, idx_A(k)) = 1;


                % G_half(k) = 0.5*(G(k)+G(k+1))

                S(row_A, idx_G(k)) = 1;

                S(row_A, idx_G(k+1)) = 1;


                % Geometric compatibility involves xi(k) and xi(k+1)

                S(row_A, idx_xi(k)) = 1;

                S(row_A, idx_xi(k+1)) = 1;

            end


            % Mechanical-equilibrium equations
            %
            % Interior equation:
            %
            % (sigma_{i+1/2}-sigma_{i-1/2})/dx
            % - force_density_i*xi_i
            % = 0.
            %
            % The stress at the two adjacent interfaces can depend on:
            %
            % - species at i-1:i+1;
            % - G at i-1:i+1;
            % - A_half at interfaces i-1 and i;
            % - xi at i-1:i+1.
            %
            % The same conservative pattern is also used for boundary rows.

            for i = 1:N

                row_xi = idx_xi(i);


                nodal_neigh = max(i-1,1):min(i+1,N);

                interface_neigh = max(i-1,1):min(i,N-1);


                for j = nodal_neigh

                    S(row_xi, idx_species(j)) = 1;

                    S(row_xi, idx_G(j)) = 1;

                    S(row_xi, idx_xi(j)) = 1;

                end


                for k = interface_neigh

                    S(row_xi, idx_A(k)) = 1;

                end

            end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Model order 11
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case 11
        
            % State:
            % U = [species; G; A_half; xi; q_half]
        
            Nspecies = Nb_species*N;
        
            offset_G  = Nspecies;
            offset_A  = offset_G + N;
            offset_xi = offset_A + (N-1);
            offset_q  = offset_xi + N;
        
            Ntot = offset_q + (N-1);
        
            idx_A  = @(k) offset_A + k;
            idx_xi = @(i) offset_xi + i;
            idx_q  = @(k) offset_q + k;
        
            S = spalloc(Ntot, Ntot, 40*(Nb_species+6)^2*N);
        

            % Species and growth equations
        
            for i = 1:N
        
                rows_species = idx_species(i);
                row_G = idx_G(i);
        
                neigh1 = max(i-1,1):min(i+1,N);
                neigh2 = max(i-2,1):min(i+2,N);
        
                A_neigh_diff = max(i-2,1):min(i+1,N-1);
                interface_neigh = max(i-1,1):min(i,N-1);
        
                % Species equations
                for j = neigh1
                    S(rows_species,idx_species(j)) = 1;
                end
        
                for j = neigh2
                    S(rows_species,idx_G(j)) = 1;
                end
        
                for k = A_neigh_diff
                    S(rows_species,idx_A(k)) = 1;
                end
        
                % Growth equation
                S(row_G,idx_species(i)) = 1;
        
                for j = neigh1
                    S(row_G,idx_G(j)) = 1;
                end
        
                for k = interface_neigh
                    S(row_G,idx_A(k)) = 1;
                    S(row_G,idx_q(k)) = 1;
                end
        
                S(row_G,idx_xi(i)) = 1;
        
            end
        

            % Compatibility equations
        
            for k = 1:N-1
        
                row_A = idx_A(k);
        
                S(row_A,idx_A(k)) = 1;
        
                S(row_A,idx_G(k)) = 1;
                S(row_A,idx_G(k+1)) = 1;
        
                S(row_A,idx_xi(k)) = 1;
                S(row_A,idx_xi(k+1)) = 1;
        
            end
        
            % Mechanical equilibrium equations
        
            for i = 1:N
        
                row_xi = idx_xi(i);
        
                nodal_neigh = max(i-1,1):min(i+1,N);
                interface_neigh = max(i-1,1):min(i,N-1);
        
                for j = nodal_neigh
                    S(row_xi,idx_species(j)) = 1;
                    S(row_xi,idx_G(j)) = 1;
                    S(row_xi,idx_xi(j)) = 1;
                end
        
                for k = interface_neigh
                    S(row_xi,idx_A(k)) = 1;
                    S(row_xi,idx_q(k)) = 1;
                end
        
            end
        
            % q equations
        
            for k = 1:N-1
        
                row_q = idx_q(k);
        
                S(row_q,idx_q(k)) = 1;
                S(row_q,idx_A(k)) = 1;
        
                S(row_q,idx_G(k)) = 1;
                S(row_q,idx_G(k+1)) = 1;
        
                S(row_q,idx_xi(k)) = 1;
                S(row_q,idx_xi(k+1)) = 1;
        
                S(row_q,idx_species(k)) = 1;
                S(row_q,idx_species(k+1)) = 1;
        
            end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Model order 12
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        case 12
        
            % State:
            % U = [species; G; A_half; xi; sigma_half; epsilon_half]
        
            Nspecies = Nb_species*N;
        
            offset_G = Nspecies;
            offset_A = offset_G + N;
            offset_xi = offset_A + (N-1);
            offset_sigma = offset_xi + N;
            offset_epsilon = offset_sigma + (N-1);
        
            Ntot = offset_epsilon + (N-1);
        
            idx_A = @(k) offset_A + k;
            idx_xi = @(i) offset_xi + i;
            idx_sigma = @(k) offset_sigma + k;
            idx_epsilon = @(k) offset_epsilon + k;
        
            S = spalloc(Ntot, Ntot, 40*(Nb_species+6)^2*N);
        

            % Species and growth equations
        
            for i = 1:N
        
                rows_species = idx_species(i);
                row_G = idx_G(i);
        
                neigh1 = max(i-1,1):min(i+1,N);
                neigh2 = max(i-2,1):min(i+2,N);
        
                A_neigh_diff = max(i-2,1):min(i+1,N-1);
                interface_neigh = max(i-1,1):min(i,N-1);
        
                % Species equations
                for j = neigh1
                    S(rows_species,idx_species(j)) = 1;
                end
        
                for j = neigh2
                    S(rows_species,idx_G(j)) = 1;
                end
        
                for k = A_neigh_diff
                    S(rows_species,idx_A(k)) = 1;
                end
        
                % Growth equation
                S(row_G,idx_species(i)) = 1;
        
                for j = neigh1
                    S(row_G,idx_G(j)) = 1;
                end
        
                for k = interface_neigh
                    S(row_G,idx_A(k)) = 1;
                    S(row_G,idx_sigma(k)) = 1;
                end
        
            end
        

            % Compatibility equations

            for k = 1:N-1
        
                row_A = idx_A(k);
        
                S(row_A,idx_A(k)) = 1;
        
                S(row_A,idx_G(k)) = 1;
                S(row_A,idx_G(k+1)) = 1;
        
                S(row_A,idx_xi(k)) = 1;
                S(row_A,idx_xi(k+1)) = 1;
        
            end
        

            % Mechanical equilibrium equations
        
            for i = 1:N
        
                row_xi = idx_xi(i);
        
                nodal_neigh = max(i-1,1):min(i+1,N);
                interface_neigh = max(i-1,1):min(i,N-1);
        
                for j = nodal_neigh
                    S(row_xi,idx_species(j)) = 1;
                    S(row_xi,idx_G(j)) = 1;
                    S(row_xi,idx_xi(j)) = 1;
                end
        
                for k = interface_neigh
                    S(row_xi,idx_A(k)) = 1;
                    S(row_xi,idx_sigma(k)) = 1;
                end
        
            end
        

            % Strain constraint and constitutive equation
        
            for k = 1:N-1
        
                row_sigma = idx_sigma(k);
                row_epsilon = idx_epsilon(k);
        
                % epsilon_half - EpsilonLawA(...) = 0
                S(row_sigma,idx_epsilon(k)) = 1;
                S(row_sigma,idx_A(k)) = 1;
        
                S(row_sigma,idx_G(k)) = 1;
                S(row_sigma,idx_G(k+1)) = 1;
        
                S(row_sigma,idx_xi(k)) = 1;
                S(row_sigma,idx_xi(k+1)) = 1;
        
                S(row_sigma,idx_species(k)) = 1;
                S(row_sigma,idx_species(k+1)) = 1;
        
                % b1*epsilon_t = a0*sigma - b0*epsilon
                S(row_epsilon,idx_sigma(k)) = 1;
                S(row_epsilon,idx_epsilon(k)) = 1;
        
            end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Model order 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        case 2
        
            % State:
            % U = [species; G; A_half; xi;
            %      sigma_half; epsilon_half; y_half; z_half]
        
            Nspecies = Nb_species*N;
        
            offset_G       = Nspecies;
            offset_A       = offset_G + N;
            offset_xi      = offset_A + (N-1);
            offset_sigma   = offset_xi + N;
            offset_epsilon = offset_sigma + (N-1);
            offset_y       = offset_epsilon + (N-1);
            offset_z       = offset_y + (N-1);
        
            Ntot = offset_z + (N-1);
        
            idx_A       = @(k) offset_A + k;
            idx_xi      = @(i) offset_xi + i;
            idx_sigma   = @(k) offset_sigma + k;
            idx_epsilon = @(k) offset_epsilon + k;
            idx_y       = @(k) offset_y + k;
            idx_z       = @(k) offset_z + k;
        
            S = spalloc(Ntot, Ntot, 50*(Nb_species+8)^2*N);
        

            % Species and growth equations
        
            for i = 1:N
        
                rows_species = idx_species(i);
                row_G = idx_G(i);
        
                neigh1 = max(i-1,1):min(i+1,N);
                neigh2 = max(i-2,1):min(i+2,N);
        
                A_neigh_diff = max(i-2,1):min(i+1,N-1);
        
                interface_neigh = max(i-1,1):min(i,N-1);
        
                % Species
                for j = neigh1
                    S(rows_species,idx_species(j)) = 1;
                end
        
                for j = neigh2
                    S(rows_species,idx_G(j)) = 1;
                end
        
                for k = A_neigh_diff
                    S(rows_species,idx_A(k)) = 1;
                end
        
                % Growth
                S(row_G,idx_species(i)) = 1;
        
                for j = neigh1
                    S(row_G,idx_G(j)) = 1;
                end
        
                for k = interface_neigh
                    S(row_G,idx_A(k)) = 1;
                    S(row_G,idx_sigma(k)) = 1;
                end
        
            end
        

            % Compatibility
        
            for k = 1:N-1
        
                row_A = idx_A(k);
        
                S(row_A,idx_A(k)) = 1;
        
                S(row_A,idx_G(k)) = 1;
                S(row_A,idx_G(k+1)) = 1;
        
                S(row_A,idx_xi(k)) = 1;
                S(row_A,idx_xi(k+1)) = 1;
        
            end
        

            % Mechanical equilibrium
        
            for i = 1:N
        
                row_xi = idx_xi(i);
        
                nodal_neigh = max(i-1,1):min(i+1,N);
        
                interface_neigh = max(i-1,1):min(i,N-1);
        
                for j = nodal_neigh
                    S(row_xi,idx_species(j)) = 1;
                    S(row_xi,idx_G(j)) = 1;
                    S(row_xi,idx_xi(j)) = 1;
                end
        
                for k = interface_neigh
                    S(row_xi,idx_A(k)) = 1;
                    S(row_xi,idx_sigma(k)) = 1;
                end
        
            end
        

            % Interface equations
        
            for k = 1:N-1
        
                row_sigma   = idx_sigma(k);
                row_epsilon = idx_epsilon(k);
                row_y       = idx_y(k);
                row_z       = idx_z(k);
        
                % sigma_t = y
                S(row_sigma,idx_y(k)) = 1;
        
                % epsilon_t = z
                S(row_epsilon,idx_z(k)) = 1;
        
                % epsilon - EpsilonLawA(...) = 0
                S(row_y,idx_epsilon(k)) = 1;
                S(row_y,idx_A(k)) = 1;
        
                S(row_y,idx_G(k)) = 1;
                S(row_y,idx_G(k+1)) = 1;
        
                S(row_y,idx_xi(k)) = 1;
                S(row_y,idx_xi(k+1)) = 1;
        
                S(row_y,idx_species(k)) = 1;
                S(row_y,idx_species(k+1)) = 1;
        
                % Constitutive equation
                S(row_z,idx_sigma(k)) = 1;
                S(row_z,idx_epsilon(k)) = 1;
                S(row_z,idx_y(k)) = 1;
                S(row_z,idx_z(k)) = 1;
        
            end

        otherwise

            error('Unknown viscoelastic model order.');

    end


    % Convert every stored nonzero value to exactly one

    S = spones(S);

end