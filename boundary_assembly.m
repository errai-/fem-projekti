% Assembly routine for P1-elements
%
% The bilinear form is given as a structure bilin :
%
% bilin = function specifying the bilinear form. The input for
%                  this function is :
%
%                  bilin(U,V,dU,dV,gX)
%
%                  - U,V  basisfunction values
%
%                   - dU,dV = cell array of function derivatives
%                     dU{i} = xi - derivative
%
%                  - global coordinate gX in a cell array
%
% linf = function specifying the load functional. The input for
%        this function is :
%
%                   linf(V,dV,gX).
%
% ALL INPUT WILL BE IN MATRIX FORM.
%


function [K,F] = boundary_assembly(mesh,bilin,linf)
%{
Ndof = size(mesh.p,2);

% initialize
[Ax,Ay,bx,by,detA,Px,Py] = affine_tri(mesh);

[X,W] = gaussint(5,0,1);
X=X';
X=[X; zeros(size(X))];

% TODO: is it ok to assume that only the second column has zeros?
bedge_ind = find( mesh.e2t(2,:) == 0);
boundary_edges = mesh.edges(:,bedge_ind);
% indices of the respective triangles:
bind = mesh.e2t(1,bedge_ind);
bAx = Ax(bind,:);
bAy = Ax(bind,:);
bbx = bx(bind);
bby = by(bind);
bdetA = detA(bind);
bPx = Px(bind);
bPy = Py(bind);

Nip = size(X,2);
Nt = size(mesh.t,2);
gX{1} = bsxfun(@plus,bAx*X,bbx);
gX{2} = bsxfun(@plus,bAy*X,bby);

L{1} = 1-X(1,:);
L{2} = X(1,:);

dL{1} = [ -ones(1,Nip); zeros(1,Nip) ];
dL{2} = [  ones(1,Nip); zeros(1,Nip) ];


% define arrays for matrix entries.
ffind = [];
ff = zeros(3,Nt);

iind = [];
jind = [];
kk = zeros(4,Nt);

mind = 1;

for i=1:2
    
    Li = repmat(L{i},Nt,1);
    
    dLi{1} = bPx*dL{i};
    dLi{2} = bPy*dL{i};
    
    ff(i,:) = linf(Li,dLi,gX)*W.*abs(bdetA);
    

    
    for j=1:2
        
        Lj = repmat(L{j},Nt,1);
        
        
        dLj{1} = bPx*dL{j};
        dLj{2} = bPy*dL{j};
               
        % Keep track of indeces : can be eliminated ??  !! ??
        iind = [iind i] ; jind = [jind j];
        
        % SIMPLIFIED HERE !!!!
        kk(mind,:) = bilin(Lj,Li,dLj,dLi,gX)*W.*abs(bdetA);
        mind = mind+1;
    end
end
%}
%K = sparse(mesh.t(iind,:),mesh.t(jind,:),kk,Ndof,Ndof);
K = sparse(size(mesh.p,2),size(mesh.p,2));
%F = sparse(mesh.t,ones(size(mesh.t)),ff,Ndof,1);
F = sparse(size(mesh.p,2), 1);

