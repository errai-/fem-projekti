% 
% total error
%

function y = total_error(mesh,x,uexact,uexact_x,uexact_y)

% sum L2 and H1 errors

y = L2_error(mesh, x, uexact) + H1_error(mesh, x, uexact_x, uexact_y);

%y = y/size(mesh.p,2);
% take square root ?
y=sqrt(y);
