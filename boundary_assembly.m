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

% TODO: is it ok to assume that only the second column has zeros?
bedge_ind = find( mesh.e2t(2,:) == 0);
boundary_edges = mesh.edges(:,bedge_ind);

B = mesh.p(:,boundary_edges(1,:))';
A = mesh.p(:,boundary_edges(2,:))' - B;

Nip = size(X,2);
Ne = size(boundary_edges,2);
for i=1:2
    gX{i} = bsxfun(@plus,A(:,i)*X,B(:,i));
end

L{1} = 1-X(1,:);
L{2} = X(1,:);

dL{1} = [ -ones(1,Nip); zeros(1,Nip) ];
dL{2} = [  ones(1,Nip); zeros(1,Nip) ];


% define arrays for matrix entries.
ffind = [];
ff = zeros(2,Ne);

iind = [];
jind = [];
kk = zeros(4,Ne);

mind = 1;

for i=1:2
    
    Li = repmat(L{i},Ne,1);
    
    %dLi{1} = bPx*dL{i};
    %dLi{2} = bPy*dL{i};
    
    ff(i,:) = linf(Li,gX)*W.*rssq(A,2);
    
    for j=1:2
        
        Lj = repmat(L{j},Ne,1);
        
        %dLj{1} = bPx*dL{j};
        %dLj{2} = bPy*dL{j};
               
        % Keep track of indeces : can be eliminated ??  !! ??
        iind = [iind i] ; jind = [jind j];
        
        % SIMPLIFIED HERE !!!!
        kk(mind,:) = bilin(Lj,Li,gX)*W.*rssq(A,2);
        mind = mind+1;
    end
end
%K = sparse(mesh.t(iind,:),mesh.t(jind,:),kk,Ndof,Ndof);
K = sparse(size(mesh.p,2),size(mesh.p,2));
%F = sparse(mesh.t,ones(size(mesh.t)),ff,Ndof,1);
F = sparse(size(mesh.p,2), 1);

