
% make the mesh
clear all;

addpath('./util');

% r value
r_val = 10;
% k values to loop over
k = 1;
% initiating and storaging for h
h_val = 1;
% Init circular mesh
mesh = discmesh(r_val,h_val);
for i=1:2,
    mesh=refine_tri(mesh);
end

% Load function
linf = @(V,dV,gX)( 0*heaviside(1-gX{1}.^2-gX{2}.^2).*V);
b_linf = @(V,gX)(k*1i*(ones(size(gX{1}))-gX{1}./r_val).*exp(-1i*k*gX{1}));
% System
bilin = @(U,V,dU,dV,gX)( -dU{1}.*dV{1} - dU{2}.*dV{2} + (k^2)*U.*V);
% System edge
b_bilin = @(U,V,gX)(-1i*k*U.*V);

[K,b] = combined_assembly(mesh,bilin,linf,b_bilin,b_linf);

% FEM solution
x = full(K\b);

tri = delaunay(mesh.p(1,:)', mesh.p(2,:)');
trisurf(tri, mesh.p(1,:)', mesh.p(2,:)', real(x));
xlabel('X'); ylabel('Y');

