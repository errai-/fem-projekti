% 
% evaluate the L2-error
%

function y = L2_error(mesh,x,uexact)

[X, W] = inttri(3);

% evaluate values at every element and every mapped integration point
[U,dU,gX] = eval2Dtri(mesh,x,X);

UE = uexact(gX{1}, gX{2});

% evaluate the integral

val = (U - UE).^2;

[Ax, Ay, bx, by, detA] = affine_tri(mesh);

% compute the integral

y = sum(val*W.*abs(detA));
