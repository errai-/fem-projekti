% Assembly routine for P1-elements
%
% The bilinear form is given as a structure bilin :
%
% bilin = function specifying the bilinear form. The input for
%                  this function is :
%
%                  bilin(U,V,gX)
%
%                  - U,V  basisfunction values
%
%                  - global coordinate gX in a cell array
%
% linf = function specifying the load functional. The input for
%        this function is :
%
%                   linf(V,gX).
%
% ALL INPUT WILL BE IN MATRIX FORM.
%


function [K,F] = boundary_assembly(mesh,bilin,linf)
% initialize
[X,W] = gaussint(5,0,1);
X=X';

% Boundary edges are those which are part of only one triangle
% boundary_edges is a 2-by-N matrix with each column having indices to the endpoints of an edge
boundary_edges = mesh.edges(:,xor(mesh.e2t(1,:) == 0, mesh.e2t(2,:) == 0));

% B and A are N-by-2 matrices
B = mesh.p(:,boundary_edges(1,:))';
A = mesh.p(:,boundary_edges(2,:))' - B;

Ne = size(boundary_edges,2);
for i=1:2
    gX{i} = bsxfun(@plus,A(:,i)*X,B(:,i));
end

L{1} = 1-X(1,:);
L{2} = X(1,:);

% define arrays for matrix entries.

ff = zeros(2,Ne);
iind = [];
jind = [];
kk = zeros(4,Ne);
mind = 1;

for i=1:2
    Li = repmat(L{i},Ne,1);
    ff(i,:) = linf(Li,gX)*W.*rssq(A,2);
    for j=1:2
        
        Lj = repmat(L{j},Ne,1);
               
        % Keep track of indices
        iind = [iind i];
        jind = [jind j];
        
        kk(mind,:) = bilin(Lj,Li,gX)*W.*rssq(A,2);
        mind = mind+1;
    end
end
K = sparse(boundary_edges(iind,:),boundary_edges(jind,:),kk,size(mesh.p,2),size(mesh.p,2));
F = sparse(boundary_edges,ones(size(boundary_edges)),ff,size(mesh.p,2),1);

end
