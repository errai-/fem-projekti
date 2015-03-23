
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

% Load function
linf = @(V,dV,gX)(0*heaviside(1-gX{1}.^2-gX{2}.^2).*V);
b_linf = @(V,gX)(exp(-1i*k*gX{2}).*V);
%b_linf = @(V,gX)((-0.5)*exp(-1i*k*rssq([gX{1}; gX{2}])./rssq([gX{1}; gX{2}]).^3));
% System
bilin = @(U,V,dU,dV,gX)(dU{1}.*dV{1} + dU{2}.*dV{2}-(k^2)*U.*V);
% System edge
b_bilin = @(U,V,gX)(1i*k*U.*V);

[K,b] = combined_assembly(mesh,bilin,linf,b_bilin,b_linf);

% FEM solution
x = full(K\b);

tri = delaunay(mesh.p(1,:)', mesh.p(2,:)');
trisurf(tri, mesh.p(1,:)', mesh.p(2,:)', abs(x));
xlabel('X'); ylabel('Y');
