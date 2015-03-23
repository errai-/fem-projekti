%
% Evaluate the finite element function at the image of a given reference
% point
%
% x = vector of coefficients
% X = coordinates of one point on reference element

function [U,dU,gX] = myeval(mesh,x,X)

U = zeros(size(mesh.t,2),1);
dU = { zeros(size(mesh.t,2),1), zeros(size(mesh.t,2),1)};

% evaluate basisfunctions
L{1} = 1 - X(1,:) - X(2,:);
L{2} = X(1,:);
L{3} = X(2,:);

gradL{1} = [ -1 ; -1];
gradL{2} = [ 1 ; 0];
gradL{3} = [ 0 ; 1];

% Loop over all elements in the mesh
for i=1:size(mesh.t,2)
   
    [A,b] = affine_map(mesh,i);
    % loop over local basisfunctions
    abu = A*X + b;
    gX{1}(i,1) = abu(1);
    gX{2}(i,1) = abu(2);
    
    for j=1:3
       
        % generate the value of the basisfunction !!
        
        Lj = L{j};
        dLj = inv(A)'* ; 
        
        % sum the values. First for the function
        
        U(i) = U(i) + x(mesh.t(j,i))*
            
        % then for the derivative
        % this is the x-derivative
        dU{1}(i) = dU{1}(i) + x(mesh.t(j,i))*
        
        % this is the y-derivative 
        dU{2}(i) = dU{2}(i) + x(mesh.t(j,i))*
    
    end
end

    
    
end