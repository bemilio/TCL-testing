close all
clear all
% size of the matrix A
r = 200;
% transition probabilities
p_stay = 0.2;
p_go = 1-p_stay;

timesteps = 150;

epsilon = eps(1);

%------END PARAMETERS------%

% grid indices
%  x_{1} ------------------ x_{r/2}
%  |                         |
%  |                         |
%  x_{r/2+1}----------------x_{r}

s = [ones(r/2, 1);
     zeros(r/2, 1) ];

%A is the matrix greek Pi!
A=zeros(r);
for i = 1:r
    if i == r/2
        A(i,i) = p_stay;
        A(r, i) = p_go;
    elseif i == r/2+1
        A(i,i) = p_stay;
        A(1, i) = p_go;
    else
        if i<=r/2
            A(i,i) = p_stay;
            A(i+1, i) = p_go;
        else
            A(i,i) = p_stay;
            A(i-1, i) = p_go;
        end
    end
end



for i = 1:r
    if i == r/2
        A_on(i,i) = p_stay;
        A_on(r, i) = p_go;
    elseif i == r/2+1
        A_on(i,i) = p_stay;
        A_on(1, i) = p_go;
    elseif i== r
        A_on(i, i) = p_stay;
        A_on(i-1,i) = p_go;
    elseif i ==1
        A_on(i, i) = p_stay;
        A_on(i+1,i) = p_go;
    else
        if i<r/2
            A_on(i,i) = p_stay;
            
            A_on(i+1, i) = p_go;
        else
            A_on(i-r/2, i) = 1;
        end
    end
end


for i = 1:r
    if i == r/2
        A_off(i,i) = p_stay;
        A_off(r, i) = p_go;
    elseif i == r/2+1
        A_off(i,i) = p_stay;
        A_off(1, i) = p_go;
    elseif i== r
        A_off(i, i) = p_stay;
        A_off(i-1,i) = p_go;
    elseif i ==1
        A_off(i, i) = p_stay;
        A_off(i+1,i) = p_go;
    else
        if i>r/2 +1
            A_off(i,i) = p_stay;
            A_off(i-1, i) = p_go;
        else
            A_off(i+r/2, i) = 1;
        end
    end
end


P_on = A_on';
[V, D, W] = eig(A_on);
[eigens, indexes] = sort(diag(D));
A_on_left_principal_eigenvector = V(:, indexes(r));
A_on_left_principal_eigenvector = A_on_left_principal_eigenvector/sum(A_on_left_principal_eigenvector);


iterations = 1000;
norm_fac = 1/iterations;

%% check if the average infinite cost converges to the eigenvector
P_pow = P_on;
P_star = norm_fac * (eye(r) + P_pow);
for k = 1:iterations
    P_pow = P_on * P_pow;
    P_star = P_star + norm_fac* P_pow;
end
P_star(find(P_star < 10^(-3))) = 0;
error = A_on_left_principal_eigenvector' - P_star(1,:);

% h = ((eye(r)/(eye(r) - P_on + P_star)) - P_star) * s;

%equations set to solve for the mean and bias costs of an MDP under WA
B = [ eye(r) - P_on, ones(r,1);
    zeros(1, r-1), 1, 0];
x = B\[s; 0];
h = x(1:r);
lambda = x(r+1);

power_expected = A_on_left_principal_eigenvector'*s;
mean_power_error = power_expected - lambda

P = A';
for i=1:r
    h_next_it(i) = max([s(i)-lambda + P_on(i, :)*h; s(i)-lambda + P(i,:)*h] );
end

improvement = h_next_it' - h;
distance_from_optimality = sum(improvement)


%% Checking whether phi tends to h plus a constant
A_inf = A_on_left_principal_eigenvector * ones(1, r);
sum2 = zeros(r,r);
for k=0:iterations
    sum2 = sum2 +  (A_on^k - A_inf);
end

check_if_h = s'*sum2;
check_if_h' - h

    
R = [  eye(r) - A_on', ones(r, 1);
        A_on_left_principal_eigenvector', 0 ];
det(R)
check_if_h_lambda = R\[s; 0];
check_if_h_lambda - [check_if_h'; lambda]

%% Simulation for phi_top under various distributions

skewed_right_distribution = zeros(r,1);
for i=1:r
    if i<=r/2
        skewed_right_distribution(i) = i;
    else
        skewed_right_distribution(i) = r-i;
    end
end
skewed_right_distribution = skewed_right_distribution/sum(skewed_right_distribution);

skewed_left_distribution = zeros(r,1);
for i=1:r
    if i<=r/2
        skewed_left_distribution(i) = r-i;
    else
        skewed_left_distribution(i) = i;
    end
end

skewed_left_distribution = skewed_left_distribution/sum(skewed_left_distribution);


initial_conditions = [ 
    ones(r,1)'/r; 
    skewed_right_distribution'; 
    skewed_left_distribution';
    [zeros(r-1, 1); 1]';
    [zeros(r/2-1, 1); 1; zeros(r/2 , 1)]';
    [zeros(r/2, 1); 1; zeros(r/2 - 1, 1)]'
    ]';
labels = [
    "uniform";
    "skewed to the right";
    "skewed to the left";
    "bottom right delta"
    "top right delta";
    "bottom left delta"
    ];

figure();
plots_label = [];
colors = lines;
for i=1:size(initial_conditions, 2)
    sum2 = 0;
    A_power = eye(r);
    for k=0:iterations
        power(k+1) = s' * A_power * initial_conditions(:,i);
        sum2 = sum2 + s'* (A_on^k - A_inf) * initial_conditions(:,i);
        A_power = A_power * A_on;
    end
    phi_top = sum2;
    plots_label(i) = plot(1:iterations+1, power, 'Color', colors(i, :));
    text(2+ i*2 , power(2+ i*2), ['\leftarrow', num2str(phi_top)], 'Color', colors(i, :));
    hold on;   
end
mytitle = ['Power over time and values of $\bar{\phi}$ with $\Pi_{stay} =$ ', num2str(p_stay)];
title(mytitle, 'Interpreter', 'latex');

legend(plots_label, labels)
%for i= 1: length(plots_label)

%end

%%  Plot h for each initial state
I = eye(r);
phi_mat = zeros(2, r/2);
for i=1:r
    sum2 = 0;
    A_power = eye(r);
    probability_mass = I(:,i);
    for k=0:iterations
        power(k+1) = s' * A_power * probability_mass;
        sum2 = sum2 + s'* (A_on^k - A_inf) * probability_mass;
        A_power = A_power * A_on;
    end
    phi_mat(floor((i-1)/(r/2)) + 1, rem(i-1, r/2)+1)  = sum2;
end 
figure();
imagesc(phi_mat);
title(['Value of $\bar{\phi}$ for each initial state with $\Pi_{stay} =$ ', num2str(p_stay)],'Interpreter', 'latex');
hold on
for i = 1:size(phi_mat, 1)
    for j = 1:size(phi_mat, 2)
        text(j, i, num2str(phi_mat(i, j), 2));
    end
end

figure();

bar(1:r, [phi_mat(1,:), phi_mat(2, :)]);
xlabel("state");
ylabel('$\bar{\phi}$', 'Interpreter', 'latex')
title(['Value of $\bar{\phi}$ for each initial state with $\Pi_{stay} =$ ', num2str(p_stay)],'Interpreter', 'latex');

%% Simulation for phi_bottom under various distributions

[V, D, W] = eig(A_off);
[eigens, indexes] = sort(diag(D));
A_off_left_principal_eigenvector = V(:, indexes(r));
A_off_left_principal_eigenvector = A_off_left_principal_eigenvector/sum(A_off_left_principal_eigenvector);
A_off_inf = A_off_left_principal_eigenvector * ones(1, r);


figure();
plots_label = [];
colors = lines;
for i=1:size(initial_conditions, 2)
    sum2 = 0;
    A_power = eye(r);
    for k=0:iterations
        power(k+1) = s' * A_power * initial_conditions(:,i);
        sum2 = sum2 + s'* (A_power - A_off_inf) * initial_conditions(:,i);
        A_power = A_power * A_off;
    end
    phi_top = sum2;
    plots_label(i) = plot(1:iterations+1, power, 'Color', colors(i, :));
    text(2+ i*2 , power(2+ i*2), ['\leftarrow', num2str(phi_top)], 'Color', colors(i, :));
    hold on;   
end
mytitle = ['Power over time and values of $\underline{\phi}$ with $\Pi_{stay} =$ ', num2str(p_stay)];
title(mytitle, 'Interpreter', 'latex');

legend(plots_label, labels)

%%  Plot phi_bottom for each initial state
I = eye(r);
phi_mat_bot = zeros(2, r/2);
for i=1:r
    sum2 = 0;
    A_power = eye(r);
    initial_conditions = I(:,i);
    for k=0:iterations
        power(k+1) = s' * A_power * initial_conditions;
        sum2 = sum2 + s'* (A_off_inf -A_power) * initial_conditions;
        A_power = A_power * A_off;
    end
    phi_mat_bot(floor((i-1)/(r/2)) + 1, rem(i-1, r/2)+1)  = sum2;
end 
figure();
imagesc(phi_mat_bot);
title(['Value of $\underline{\phi}$ for each initial state with $\Pi_{stay} =$ ', num2str(p_stay)],'Interpreter', 'latex');
hold on
for i = 1:size(phi_mat_bot, 1)
    for j = 1:size(phi_mat_bot, 2)
        text(j, i, num2str(phi_mat_bot(i, j), 2));
    end
end

figure();

bar(1:r, [phi_mat_bot(1,:), phi_mat_bot(2, :)]);
xlabel("state");
ylabel('$\bar{\phi}$', 'Interpreter', 'latex')
title(['Value of $\underline{\phi}$ for each initial state with $\Pi_{stay} =$ ', num2str(p_stay)],'Interpreter', 'latex');


phi_tot = phi_mat + phi_mat_bot; 
figure;
imagesc(phi_tot);
hold on
for i = 1:size(phi_tot, 1)
    for j = 1:size(phi_tot, 2)
        text(j, i, num2str(phi_tot(i, j), 2));
    end
end
title(['Value of $\underline{\phi}+ \bar{\phi}$ for each initial state with $\Pi_{stay} =$ ', num2str(p_stay)],'Interpreter', 'latex');