% Plot the H1-error
function plot_error(mesh, k_vec, x, uexact_x, uexact_y)

error = H1_error(mesh, x, uexact_x, uexact_y);

figure;
plot(k_vec, error);