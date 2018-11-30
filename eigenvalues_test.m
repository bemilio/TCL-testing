clear all
close all

k=0;
for a = linspace(0, 1, 100)
    F = [ a, 0, 0, 1;
        1-a, a, 0, 0;
        0, 1-a, a, 0;
        0, 0, 1-a, 0 ];
    [V, D, W] = eig(F);
    [eigens, indexes] = sort(diag(D));
    A_on_left_principal_eigenvector = V(:, indexes(4));
    A_on_left_principal_eigenvector = A_on_left_principal_eigenvector/sum(A_on_left_principal_eigenvector);
    k = k+1;
    eigenvalue_norm(1, k) = norm(eigens(end-1));
    eigenvalue_real(1, k) = real(eigens(end-1));
    eigenvalue_imag(1, k) = imag(eigens(end-1));
    eigenvalue_norm(2, k) = norm(eigens(end-2));
    eigenvalue_real(2, k) = real(eigens(end-2));
    eigenvalue_imag(2, k) = imag(eigens(end-2));
    eigenvalue_norm(3, k) = norm(eigens(end-3));
    eigenvalue_real(3, k) = real(eigens(end-3));
    eigenvalue_imag(3, k) = imag(eigens(end-3));
end

plot(linspace(0,1, 100), eigenvalue_norm(1, :));
hold on 
plot(linspace(0,1, 100), eigenvalue_norm(2, :));
plot(linspace(0,1, 100), eigenvalue_norm(3, :));
title("norm");
figure
plot(linspace(0,1, 100), eigenvalue_real(1, :));
hold on 
plot(linspace(0,1, 100), eigenvalue_real(2, :));
plot(linspace(0,1, 100), eigenvalue_real(3, :));
title("real");
figure
plot(linspace(0,1, 100), eigenvalue_imag(1, :));
hold on 
plot(linspace(0,1, 100), eigenvalue_imag(2, :));
plot(linspace(0,1, 100), eigenvalue_imag(3, :));
title("imag");

% figure
% scatter(second_eigenvalue_real, second_eigenvalue_imaginary)