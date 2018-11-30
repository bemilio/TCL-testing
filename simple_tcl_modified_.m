clf
cla
clc
clear variables
close all


%User parameters
%chices: 
% 1 - sinusoid
% 2 - constant
reference = 2;


%% Building the model
% size of the matrix A
r = 20;
% transition probabilities
p_stay = 0.2;
p_go = 0.8;

% grid indices
%  x_{1} ------------------ x_{r/2}
%  |                         |
%  |                         |
%  x_{r/2+1}----------------x_{r}

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
            A(i+1,i) = p_go;
        else
            A(i,i) = p_stay;
            A(i-1,i) = p_go;
        end
    end
end


% grid indices
%  x_{1} ------------------ x_{r/2}
%  |                         |
%  |                         |
%  x_{r/2+1}----------------x_{r}

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
s = [ones(r/2, 1);
     zeros(r/2, 1) ];

[V, D, W] = eig(A_on);
[eigens, indexes] = sort(diag(D));
A_on_left_principal_eigenvector = V(:, indexes(r));
A_on_left_principal_eigenvector = A_on_left_principal_eigenvector/sum(A_on_left_principal_eigenvector);
lambda = A_on_left_principal_eigenvector'*s;

%B has value 1 if the corresponding transition prob. can be optimized
B=zeros(r);
for i = 1:r
    if  i>1 && i<r/2
        
        B(i+r/2, i) = 1 ;
%         B(i,i) = 1;
%         B(i,i+1) = 1;
    end
    if  i>r/2+1 && i<r
        
        B(i-r/2, i) = 1;
%         B(i,i) = 1;
%         B(i, i-1) = 1;
    end
end
% eigenvalues and eigenvectors for steady state
[eigvects,eigvals] = eig(A);
eigvals = diag(eigvals);

%% Simulation

P = A;

n = size(P,1);

% number of timesteps
N = 100;
one = ones(n,1);
on = one;
on(r/2+1:end)=zeros(r/2,1);
% reference to be tracked
if reference == 1
    p_star(1)= 0.5;
    for k = 1:N
        p_star(k+1) = 0.5+0.2*sin(0.1*k);
    end
elseif reference == 2
    p_star = ones(N+1, 1);
end

% initial condition
% random_init = rand(n,1);
epsilon = 1*10^(-7);
rho0 = [epsilon * ones(r-1, 1); 1 - (r-1)*epsilon];


Constraints = [];
Objective = 0;

% Define variables
M = sdpvar(n,n,N,'full');
rho = sdpvar(n,N+1);
sdpvar alpha;

% Define constraints and objective
rho(:,1) = rho0;
for k = 1:N
    for i=1:n
        for j = 1:n
            if B(i,j)==0
                if P(i,j) == 0
                    %if B=0 and P=0 the controller can't create the
                    %non-existing edge of the MC
                    M(i,j,k)=0;
                else
                    %if B=0 and P!=0 the transition is constrained
                    
                    if (j==1|| j==r/2 || j==r/2+1 || j==r) 
                        %case when i == {1; r/2; r/2+1; r} (incontrollable
                        %states)
                        Constraints = [Constraints; M(i,j,k)==P(i,j)*rho(j,k)];
                    else
                        %P_controlled(i_off,j) = P(i_off,j)(1-P_controlled(i_on,j))
                        index_complementary = find(B(:,j)>0);
                        Constraints = [Constraints; M(i,j,k) == P(i,j) * (rho(j, k) - M(index_complementary, j, k))];
                    end 
                end
            end
        end
    end 
    
    Constraints = [Constraints;
        M(:,:,k)>=zeros(n,n);
        rho(:,k+1)>=epsilon*ones(n,1);
        M(:,:,k)*one==rho(:,k+1);
        one'*M(:,:,k)==rho(:,k)';
        one'*rho(:,k+1)==1;
        ];
    
    
    Objective = Objective + abs(on'*rho(:,k+1));% on'*rho(:,k) - lambda;  %- alpha*p_star(k+1));%+alpha* sum(sum(abs(M(:,:,k)-P*diag(rho(:,k)))));
    
end



% Set some options for YALMIP and solver
options = sdpsettings('verbose',0,'solver','mosek');


M_opt = optimizer(Constraints,Objective,options,alpha,M);
rho_opt = optimizer(Constraints,Objective,options,alpha,rho);


%% Simulate with different alpha parameters

alpha_val = 0; % linspace(0,1, 21);
% alpha_val = 5;
counter =0;
palette = parula;
for l=1:numel(alpha_val)
    
    alpha = alpha_val(l);
    
    disp('l = ');disp(l)
    rho = rho_opt(alpha);
    M = M_opt(alpha);
    
    tracking(l) = 0;
    disconfort(l) = 0;
    
    
    for k = 1:N
        
        disconfort(l) = disconfort(l) + sum(sum(abs(M(:,:,k)-P*diag(rho(:,k)))));
        tracking(l) = tracking(l) + abs(on'*rho(:,k+1) - p_star(k+1));
        
    end
    
    figure(1)
    jump_factor = floor( size(palette,1)/numel(alpha_val));
    c = palette(jump_factor * (counter+1), :);
    plot(1:N+1,on'*rho, 'Color', c);
    hold on
    %plot(1:N+1, alpha*p_star, ':', 'Color', c);
    xlabel('time')
    ylabel('power')
    
    counter = counter + 1;
    
%     if (counter==5)
%     
%         figure(3)
%         M_test = M(:, :,10);
%         imagesc(M_test);
% 
%         title(' M matrix at step 5');
%     end
end

for i=1:N
    pi_opt(:,:,i)=  M (:,:, i)/ diag(rho(:, i));
end 



figure(1)
plot(1:N+1,p_star,'r')
xlabel('time')
ylabel('power')
legend('Normalized power','Reference')



figure(2)
plot(disconfort,tracking,'-x')
xlabel('total disconfort')
ylabel('tracking error')
legend('Normalized power','Reference')



% figure
% for k = 1:N
% 
%     Pi(:,:,k) = M(:,:,k)/diag(rho(:,k));
%     U(:,:,k) = Pi(:,:,k)-P;
% 
% 
%     subplot(2,1,1)
%     plot(1:n/2,rho(1:n/2,k))
%     hold on
%     subplot(2,1,2)
%     plot(1:n/2,rho(n/2+1:end,k))
%     hold on
% 

% end






