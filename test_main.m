
% make the mesh
clear all;

addpath('./util');

% r value
r_val = 30;
% k values to loop over
k = 10;
% initiating and storaging for h
h_val = 0.5;
% Init circular mesh
mesh = discmesh(r_val,h_val);
% Evaluate the length of the longest edge in the mesh 
px = mesh.p(1,:);    % 1xNpoint - vector
py = mesh.p(2,:);    % 1xNpoint - vector
ex = px(mesh.edges); % 2xNedges - matrix.
ey = py(mesh.edges); % 2xNedges - matrix.

% pythagoras
len = (ex(1,:)-ex(2,:) ).^ 2  + (ey(1,:) -ey(2,:)).^2;
h_val = sqrt(max(len));
display(h_val*k);

% Load function
linf = @(V,dV,gX)( 0*heaviside(1-gX{1}.^2-gX{2}.^2).*V);
b_linf = @(V,gX)(1*k*1i*(ones(size(gX{1}))-gX{1}./r_val).*exp(-1i*k*gX{1}).*V);
% System
bilin = @(U,V,dU,dV,gX)(dU{1}.*dV{1} + dU{2}.*dV{2} - (k^2)*U.*V);
% System edge
b_bilin = @(U,V,gX)(1i*k*U.*V);

[K,b] = combined_assembly(mesh,bilin,linf,b_bilin,b_linf);

% FEM solution
x = full(K\b);

tri = delaunay(mesh.p(1,:)', mesh.p(2,:)');
trisurf(tri, mesh.p(1,:)', mesh.p(2,:)', real(x),'EdgeColor','none','LineStyle','none','FaceLighting','phong');
xlabel('X'); ylabel('Y');
set(gca,'XLim',[-30 30],'YLim',[-30 30]); colorbar;

figure;
trisurf(tri, mesh.p(1,:)', mesh.p(2,:)', abs(x).^2,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
xlabel('X'); ylabel('Y');
%set(gca,'XLim',[-20 20],'YLim',[-20 20],'ZLim',[0,0.6]); colorbar;

