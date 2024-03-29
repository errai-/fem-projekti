%
% evaluate the H1-error
%
function y = H1_error(mesh, x, uexact_x, uexact_y)

[X, W] = inttri(3);

% evaluate values at every element and every mapped integration point
[U,dU,gX] = eval2Dtri(mesh,x,X);

UEX = uexact_x(gX{1}, gX{2});
UEY = uexact_y(gX{1}, gX{2});

% evaluate the integral

val = abs(UEX - dU{1}).^2 + abs(UEY - dU{2}).^2;
[Ax, Ay, bx, by, detA] = affine_tri(mesh);

% compute
y = sum(val*W.*abs(detA));