
% make the mesh
clear all;

addpath('./util');

% r values to loop over
r_vals = [10]';
% k values to loop over
k_vals = [0.1,2]';
% initiating and storaging for h
init_h = 1;
h_iters = 2;
h_storage = zeros(h_iters,size(r_vals,1));

% Domain radius loop
for r_idx = 1:size(r_vals,1)
    % Init circular mesh
    mesh = discmesh(r_vals(r_idx),init_h);
    
    % Mesh density loop
    for h_idx=1:h_iters
        mesh = refine_tri(mesh);
        
        % Evaluate the length of the longest edge in the mesh 
        px = mesh.p(1,:);    % 1xNpoint - vector
        py = mesh.p(2,:);    % 1xNpoint - vector
        ex = px(mesh.edges); % 2xNedges - matrix.
        ey = py(mesh.edges); % 2xNedges - matrix.
        
        % pythagoras
        len = (ex(1,:)-ex(2,:) ).^ 2  + (ey(1,:) -ey(2,:)).^2;
        h_storage(h_idx,r_idx) = sqrt(max(len));
        
        % K value loop 
        for k_idx = 1:size(k_vals,1)
            k=k_vals(k_idx);
            
            % Load function
            linf = @(V,dV,gX)(0*heaviside(1-gX{1}.^2-gX{2}.^2).*V);
            b_linf = @(V,gX)(k*1i*(ones(size(gX{1}))-gX{1}./r_vals(r_idx)).*exp(-1i*k*gX{1}));
            % System
            bilin = @(U,V,dU,dV,gX)(- dU{1}.*dV{1} - dU{2}.*dV{2} + (k^2)*U.*V);
            % System edge
            b_bilin = @(U,V,gX)(-1i*k*U.*V);
        
            [K,b] = combined_assembly(mesh,bilin,linf,b_bilin,b_linf);

            % FEM solution
            x = full(K\b);

            % errors(k_idx) = H1_error(mesh, x, uexact_x, uexact_y);
            %if (k_idx == size(k_vals,1))
                tri = delaunay(mesh.p(1,:)', mesh.p(2,:)');
                trisurf(tri, mesh.p(1,:)', mesh.p(2,:)', real(x));
                xlabel('X'); ylabel('Y');
                pause;
            %end
        end
        % plot_error(k_vals, error)
    end
end

