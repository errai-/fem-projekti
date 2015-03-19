
% make the mesh
clear all;

mesh = make_rect_mesh(1);

% Storages
h = zeros(8,1);

% System example
bilin = @(U,V,dU,dV,gX)(dU{1}.*dV{1} + dU{2}.*dV{2});
lect_bilin = @(U,V,dU,dV,gX)(bilin(U,V,dU,dV,gX)+epsilon*U.*V );

% Load example
linf = @(V,dV,gX)( cos(pi*gX{1}).*cos(pi*gX{2}).*V );

for i=1:8
    % TODO: replace with a round mesh
    mesh = refine_tri(mesh);
    
    [K,b] = simple_assembly(mesh,bilin,linf);
    % TODO: System needs an additional path integral over the boundary edge
    %K = K + 
    
    % Evaluate the length of the longest edge in the mesh 
    px = mesh.p(1,:); % 1xNpoint - vector
    py = mesh.p(2,:); % 1xNpoint - vector

    ex = px(mesh.edges); % 2xNedges - matrix.
    ey = py(mesh.edges); % 2xNedges - matrix.

    % evaluate the length
    len = (ex(1,:)-ex(2,:) ).^ 2  + (ey(1,:) -ey(2,:)).^2;

    h(i) = sqrt(max(len));
end


