
% make the mesh
clear all;

addpath('./util');

% r values to loop over
r_vals = [1]';
% k values to loop over
k_vals = linspace(6,50,20)';
% initiating and storaging for h
init_h = 1;
h_iters = 12;
h_storage = zeros(h_iters,size(r_vals,1));

% Domain radius loop
for r_idx = 1:size(r_vals,1)
    % Init circular mesh
    %mesh = discmesh(r_vals(r_idx),init_h);
    
    % Mesh density loop
    for h_idx=1:h_iters
        %mesh = refine_tri(mesh);
        
        mesh = discmesh(r_vals(r_idx),init_h*(0.7^h_idx));

        h_storage(h_idx,r_idx) = longest_edge(mesh);

        % K value loop 
        for k_idx = 1:size(k_vals,1)
            k=k_vals(k_idx);
            
            % exact u and the derivatives of u
            uexact = @(x,y)(exp(-1i*k*x));
            uexact_x = @(x,y)(-1i*k*exp(-1i*k*x));
            uexact_y = @(x,y)(0);
            
            % Load function 
            linf = @(V,dV,gX)(0*heaviside(1-gX{1}.^2-gX{2}.^2).*V);
            b_linf = @(V,gX)(1*k*1i*(ones(size(gX{1}))-gX{1}./sqrt(gX{1}.^2+gX{2}.^2)).*exp(-1i*k*gX{1}).*V);
            % System
            bilin = @(U,V,dU,dV,gX)( dU{1}.*dV{1} + dU{2}.*dV{2} - (k^2)*U.*V);
            % System edge
            b_bilin = @(U,V,gX)(1i*k*U.*V);
        

            [K,b] = combined_assembly(mesh,bilin,linf,b_bilin,b_linf);

            % FEM solution

            x = full(K\b);

            %errors(h_idx,k_idx) = L2_error(mesh, x, uexact);
            errors(h_idx,k_idx) = total_error(mesh, x, uexact, uexact_x, uexact_y);
            %if (k_idx == size(k_vals,1))
                %tri = delaunay(mesh.p(1,:)', mesh.p(2,:)');
                %trisurf(tri, mesh.p(1,:)', mesh.p(2,:)', real(x));
                %xlabel('X'); ylabel('Y');
                %pause;
            %end
        end
    end
    plot_error(h_storage(:,r_idx),k_vals,errors,r_vals(r_idx));
end

