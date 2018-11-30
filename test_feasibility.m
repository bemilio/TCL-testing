clf
cla
clc
clear variables
% close all


%User parameters
%chices: 
% 1 - sinusoid
% 2 - constant
reference = 2;


%% Building the model
% size of the matrix A
r = 14;
% number of timesteps
N = 30;
% transition probabilities
p_stay = 0.2;
p_go = 1-p_stay;

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

one = ones(n,1);
on = one;
on(r/2+1:end)=zeros(r/2,1);

epsilon = 1*10^(-7);

for rho0_iterations= 1:1
    P_track = zeros(N,1);
    Constraints = [];
    Objective = 0;

%     rho0 = rand(r);
    rho0 = ones(r,1);
    rho0 = rho0/sum(rho0);
    counter =0;
    palette = parula;
    [ top_limit, bottom_limit ] = calculateTube(rho0, N, A_on, A_off, on);
    for reference_iterations = 1:1
        for i = 1:N
            P_track(i) = rand * (top_limit(i) - bottom_limit(i)) + bottom_limit(i);
        end


        % Define variables
        M = sdpvar(n,n,N,'full');
        rho = sdpvar(n,N+1);
%         P_track = sdpvar(N, 1);

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
                %power tracking 
                on'*rho(:,k+1) <= P_track(k)+epsilon;
                on'*rho(:,k+1) >= P_track(k)-epsilon;
                ];
            Objective = 0;

        end
        
        % Set some options for YALMIP and solver
        options = sdpsettings('verbose',0,'solver','mosek');
        % M_opt = optimizer(Constraints,Objective,options,[ rho0' , P_track' ],M);
%          rho_opt = optimizer(Constraints,Objective,options,[ rho0' , P_track' ],rho);
        %% Simulate with different parameters
        error_flag = optimize(Constraints,Objective,options);
       
%         [ rho, error_flag ] = rho_opt([initial_cond', P_track']);
        if(rho0_iterations == 1 &&  reference_iterations == 1)
            rho_opt = value(rho);
            figure(1)
            plot(1:N+1,on'*rho_opt, 'k');
            hold on
            plot(2:N+1, P_track, '--m');
            plot(2:N+1, top_limit, '--r');
            plot(2:N+1, bottom_limit, '--g');
            xlabel('time')
            ylabel('power')
        end
%         M = M_opt([rho_0', P_track']);
        counter = counter + 1;
        if(error_flag.problem == 1)
            disp("problem is infeasible!!!");
        end
        if(error_flag.problem > 1)
            disp("something else happened");
            error_flag
        end
    end
end

% for i=1:N
%     pi_opt(:,:,i)=  M (:,:, i)/ diag(rho(:, i));
% end 



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






