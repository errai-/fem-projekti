function x = solver_interference(mesh,k)

% Load function 
linf = @(V,dV,gX)(100*(heaviside(0.05-(gX{1}-0.2).^2-gX{2}.^2)+heaviside(0.05-(gX{1}+0.2).^2-gX{2}.^2)).*V);
b_linf = @(V,gX)(0*k*1i*(ones(size(gX{1}))-gX{1}./sqrt(gX{1}.^2+gX{2}.^2)).*exp(-1i*k*gX{1}).*V);
% System
bilin = @(U,V,dU,dV,gX)( dU{1}.*dV{1} + dU{2}.*dV{2} - (k^2)*U.*V);
% System edge
b_bilin = @(U,V,gX)(1i*k*U.*V);

[K,b] = combined_assembly(mesh,bilin,linf,b_bilin,b_linf);
% FEM solution
x = full(K\b);
end
