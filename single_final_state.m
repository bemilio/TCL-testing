clf
cla
clc
% clear variables
% close all


%User parameters
%chices: 
% 1 - sinusoid
% 2 - constant
reference = 2;


%% Building the model
% size of the matrix A
r = 12;
% number of timesteps
N = 40;
% transition probabilities
p_stay = 0.4;
p_go = 1-p_stay;

% grid indices
%  x_{1} ------------------ x_{r/2}
%  |                         |
%  |                         |
%  x_{r/2+1}----------------x_{r}

A=zeros(r);
for i = 1:r
    if i == r/2
        A(i,i) = 1 - p_go/2;
        A(r, i) = p_go/2;
    elseif i == r/2+1
        A(i,i) = 1 - p_go/2;
        A(1, i) = p_go/2;
    else
        if i<r/2
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
        A_on(i,i) = A(i,i);
        A_on(r, i) = A(r,i); %maybe modify?
    elseif i == r/2+1
        A_on(i,i) = A(i,i);
        A_on(1, i) = A(1,i); %maybe modify?
    else
        if i<r/2
            A_on(i,i) = p_stay;
            A_on(i+1, i) = p_go;
        else
            A_on(i-(r/2 + 1), i) = 1;
        end
    end
end


for i = 1:r
    if i == r/2
        A_off(i,i) = A(i,i);
        A_off(r, i) = A(r,i);
    elseif i == r/2+1
        A_off(i,i) = A(i,i);
        A_off(1, i) = A(1,i);
    else
        if i>r/2 +1
            A_off(i,i) = p_stay;
            A_off(i-1, i) = p_go;
        else
            A_off(i+(r/2 + 1), i) = 1;
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
    if  i>=1 && i<r/2
        B(i+(r/2 +1), i) = 1 ;
    end
    if  i>r/2+1 && i<=r
        B(i-(r/2 +1), i) = 1;
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
on(r/2) = 0.5;
on(r/2 + 1) = 0.5;

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



initial_conditions = [ 
%     ones(1,r)/r;
    [epsilon*ones(r-1, 1); 1-(r-1)*epsilon]';
%     [epsilon*ones(r/2-1, 1); 1-(r-1)*epsilon; epsilon*ones(r/2 , 1)]';
%     [epsilon*ones(r/2, 1); 1-(r-1)*epsilon; epsilon*ones(r/2 - 1, 1)]'
    ]';
labels = [
%     "uniform"
    "bottom right delta"
%     "top right delta";
%     "bottom left delta"
    ];

    

for rho0_iterations= 1:size(initial_conditions, 2)
%     P_track = zeros(N,1);
    Constraints = [];
    Objective = 0;
    rho0 = initial_conditions(:,rho0_iterations);
    rho0 = rho0/sum(rho0);
    counter =0;
    palette = parula;
    [ top_limit, bottom_limit ] = calculateTube(rho0, N, A_on, A_off, on);
%     [ top_limit_LP, bottom_limit_LP ] = calculateTubeSingleFinalStateLP(rho0, N, P, B, on); 
    for reference_iterations = 1:50
        for k = 1:N
%             P_track(k) = ((top_limit(k)+bottom_limit(k)) /2 ) + ((top_limit(k)-bottom_limit(k))/2)*cos(pi * reference_iterations*k);
%             P_track(k) = ((top_limit(k)+bottom_limit(k)) /2 )  + (top_limit(k)-bottom_limit(k)) * (rand-0.5); 
        end
%         for k = 1:N/2
%             P_track(k) = bottom_limit(k);
%         end
%         for k = N/2 +1:N
%             P_track(k) = top_limit(k);
%         end
        Constraints = [];
        Objective = 0;
        rho0 = initial_conditions(:,rho0_iterations);
        rho0 = rho0/sum(rho0);
        % Define variables
        M = sdpvar(n,n,N,'full');
        rho = sdpvar(n,N+1, 'full');
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
                            if (j==1 + n/2 || j==n/2) 
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
%                 on'*rho(:,k+1) <= P_track(k)+r*epsilon;
%                 on'*rho(:,k+1) >= P_track(k)-r*epsilon;
                ];
            Objective = Objective +(P_track(k) - on'*rho(:,k+1))^2;

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
            figure()
            plot(1:N+1,on'*rho_opt, 'k');
            hold on
            plot(2:N+1, P_track, '--m');
            plot(2:N+1, top_limit, '--r');
            plot(2:N+1, bottom_limit, '--g');
%             plot(2:N+1, top_limit_LP(2:end), 'b');
%             plot(2:N+1, bottom_limit_LP(2:end), 'c');
            xlabel('time')
            ylabel('power')
        end
%         M = M_opt([rho_0', P_track']);
        counter = counter + 1;
        if(error_flag.problem == 1)
            figure,
            rho_opt = value(rho);
            plot(1:N+1,on'*rho_opt, 'k');
            hold on
            plot(2:N+1, P_track, '--m');
            plot(2:N+1, top_limit, '--r');
            plot(2:N+1, bottom_limit, '--g');
            xlabel('time')
            ylabel('power')
            disp("problem is infeasible!!!");
         end
        if(error_flag.problem > 1)
            disp("something else happened");
            error_flag
        end
        rho_val = value(rho);
        rho_to_check = rho_val(r/2, :);
        temp = zeros(r, N);
        temp(:, 1) = rho0;
        for i=1:N
            temp(:, i+1) = (A_on^i) * rho0;
            %this has to be bigger than the other
        end
        rho_to_check_2 = temp(r/2,:);
        x = find(rho_to_check > rho_to_check_2);
        if max(x)>0
            disp("tour hyp. is invalid");
        end
    end
end

% for i=1:N
%     pi_opt(:,:,i)=  M (:,:, i)/ diag(rho(:, i));
% end 



% figure(1)
% plot(1:N+1,p_star,'r')
% xlabel('time')
% ylabel('power')
% legend('Normalized power','Reference')


% 
% figure(2)
% plot(disconfort,tracking,'-x')
% xlabel('total disconfort')
% ylabel('tracking error')
% legend('Normalized power','Reference')
% 


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






