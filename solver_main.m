
% make the mesh
clear all;

addpath('./util');

% r values to loop over
r_vals = [10,20]';
% k values to loop over
k_vals = [1,2]';
% initiating and storaging for h
init_h = 1
h_iters = 8;
h_storage = zeros(h_iters,size(r_vals,1));

% Domain radius loop
for r_idx = 1:size(r_vals,1)
    % Init circular mesh
    mesh = discmesh(r_vals(r_idx),1);
    
    % Mesh density loop
    for h_idx=1:h_iters
        mesh = refine_tri(mesh);
        
        % Evaluate the length of the longest edge in the mesh 
        px = mesh.p(1,:); % 1xNpoint - vector
        py = mesh.p(2,:); % 1xNpoint - vector
        ex = px(mesh.edges); % 2xNedges - matrix.
        ey = py(mesh.edges); % 2xNedges - matrix.
        
        % pythagoras
        len = (ex(1,:)-ex(2,:) ).^ 2  + (ey(1,:) -ey(2,:)).^2;
        h_storage(h_idx,r_idx) = sqrt(max(len));
        
        % K value loop 
        for k = k_vals
            % System
            bilin = @(U,V,dU,dV,gX)(dU{1}.*dV{1} + dU{2}.*dV{2}-(k^2)*U.*V);
            % System edge
            bilin_edge = @(U,V,dU,dV,gX)(-1i*k*U.*V);
        
            [K,b] = simple_assembly(mesh,bilin,linf);
            % TODO: System needs an additional path integral over the boundary edge
            %K = K + ;
        end    
    end
end

