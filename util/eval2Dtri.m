% Antti Hannukainen
%
% Evaluate value of P1 - element at mapped points X
%
%

function [U,dU,gX] = eval2Dtri(mesh,u,X)

Ndof = size(mesh.p,2);

% initialize
[Ax,Ay,bx,by,detA,Px,Py] = affine_tri(mesh);


gX{1} = bsxfun(@plus,Ax*X,bx);
gX{2} = bsxfun(@plus,Ay*X,by);

L{1} = 1-X(1,:)-X(2,:);
L{2} = X(1,:);
L{3} = X(2,:);

Nip = size(X,2);
Nt = size(mesh.t,2);

dL{1} = [ -ones(1,Nip); -ones(1,Nip) ];
dL{2} = [  ones(1,Nip); zeros(1,Nip) ];
dL{3} = [ zeros(1,Nip);  ones(1,Nip) ];

for i=1:3
    
    Li = repmat(L{i},Nt,1);
    
    dLi{1} = Px*dL{i};
    dLi{2} = Py*dL{i};

    if(i==1)
      U = bsxfun(@times,Li,u(mesh.t(i,:) ));
    else
      U = U + bsxfun(@times,Li,u(mesh.t(i,:) ));
    end
   
    
    if(i==1)
      dU{1} = bsxfun(@times,dLi{1},u(mesh.t(i,:) ));
      dU{2} = bsxfun(@times,dLi{2},u(mesh.t(i,:) ));

    else
      dU{1} = dU{1} + bsxfun(@times,dLi{1},u(mesh.t(i,:) ));
      dU{2} = dU{2} + bsxfun(@times,dLi{2},u(mesh.t(i,:) ));
      
    end

    
end