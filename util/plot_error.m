% Plot the H1-error
function plot_error(h_vec, k_vec, error, r)

figure;
loglog(h_vec, error);
hold on;
title(['H1-error of a plane wave in a disk of radius ' num2str(r)]);
xlabel('h');
ylabel('H1 error');
for i=1:size(k_vec)
    asd{i} = num2str(k_vec(i)); 
end
legend(strcat('k = ', asd));
end