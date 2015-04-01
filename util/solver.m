function [x,err] = solver(mesh,k)
% exact u and the derivatives of u
uexact = @(x,y)(exp(-1i*k*x));
uexact_x = @(x,y)(-1i*k*exp(-1i*k*x));
uexact_y = @(x,y)(0); 
% Load function 
linf = @(V,dV,gX)(0*gX{1});
b_linf = @(V,gX)(1*k*1i*(ones(size(gX{1}))-gX{1}./sqrt(gX{1}.^2+gX{2}.^2)).*exp(-1i*k*gX{1}).*V);
% System
bilin = @(U,V,dU,dV,gX)( dU{1}.*dV{1} + dU{2}.*dV{2} - (k^2)*U.*V);
% System edge
b_bilin = @(U,V,gX)(1i*k*U.*V);

[K,b] = combined_assembly(mesh,bilin,linf,b_bilin,b_linf);
% FEM solution
x = full(K\b);
%errors(h_idx,k_idx) = L2_error(mesh, x, uexact);
err = total_error(mesh, x, uexact, uexact_x, uexact_y);
end
