% Assembles the system and load for a weak Helmholtz problem

function [K,F] = simple_assembly(mesh,bilin,linf,b_bilin,b_linf)

[K,F] = simple_assembly(mesh,bilin,linf);
[Ktmp,Ftmp] = boundary_assembly(mesh,b_bilin,b_linf);
K = K+Ktmp;
F = F+Ftmp;

end