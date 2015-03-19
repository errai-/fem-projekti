
% make the mesh
clear all;

% TODO: initial circular mesh
%mesh = ;

% Storages
h = zeros(8,1);

k = 1;

% System
bilin = @(U,V,dU,dV,gX)(dU{1}.*dV{1} + dU{2}.*dV{2}-(k^2)*U.*V);
% System edge
bilin_edge = @(U,V,dU,dV,gX)(-1i*k*U.*V);
% Load (could be almost anyting)
linf = @(V,dV,gX)( (heaviside(1-gX{1}.^2-gX{2}.^2)).*V );

for h_idx=1:8
    mesh = refine_tri(mesh);
    
    [K,b] = simple_assembly(mesh,bilin,linf);
    % TODO: System needs an additional path integral over the boundary edge
    %K = K + ;
    
    % Evaluate the length of the longest edge in the mesh 
    px = mesh.p(1,:); % 1xNpoint - vector
    py = mesh.p(2,:); % 1xNpoint - vector

    ex = px(mesh.edges); % 2xNedges - matrix.
    ey = py(mesh.edges); % 2xNedges - matrix.

    % evaluate the length
    len = (ex(1,:)-ex(2,:) ).^ 2  + (ey(1,:) -ey(2,:)).^2;

    h(idx) = sqrt(max(len));
end


