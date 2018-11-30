function [top_limit,bottom_limit] = calculateTubeSingleFinalStateLP(rho0, time_window, P, B, on)
%CALCULATE_TUBE Summary of this function goes here
%   Detailed explanation goes here

epsilon = 1*10^(-7);
r = size(P,1);
one = ones(r,1);

expect_power_right = 1/2;
expect_power_left = 1/2;

%in the model with only 1 extreme state, s is the vector that gives the
%maximum power obtainable at each state
s_top = one;
s_top(r/2 - 1) = P(r/2-1, r/2-1) +  P(r/2, r/2-1) * expect_power_right;
s_top(r/2) = P(r/2, r/2) * expect_power_right;
s_top(r/2 + 1) = P(r/2+1, r/2 + 1) * expect_power_left + P(1, r/2 +1);

s_bottom = zeros(r,1);
%...

top_limit =zeros(time_window+1,1);
top_limit(1) = 1; % nonsense bound at time 0
top_limit(2) = s_top'* rho0;

bottom_limit(1) = 0; % nonsense bound at time 0
bottom_limit(2) = s_bottom'* rho0;


for N=1:time_window-1
    Constraints = [];
    Objective_top = 0;
    M = sdpvar(r,r,N,'full');
    rho = sdpvar(r,N+1, 'full');
    rho(:,1) = rho0;
    
    for k = 1:N
        for i=1:r
            for j = 1:r
                if B(i,j)==0
                    if P(i,j) == 0
                        %if B=0 and P=0 the controller can't create the
                        %non-existing edge of the MC
                        M(i,j,k)=0;
                    else
                        %if B=0 and P!=0 the transition is constrained
                        if (j==1 + r/2 || j==r/2) 
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
            M(:,:,k)>=zeros(r,r);
            rho(:,k+1)>=epsilon*ones(r,1);
            M(:,:,k)*one==rho(:,k+1);
            one'*M(:,:,k)==rho(:,k)';
            one'*rho(:,k+1)==1;
%             on'*rho(:,k+1)<= top_limit(k+1);
%             on'*rho(:,k+1)>= bottom_limit(k+1);
            ];
    end
    Objective_top = s_top' * rho(:, N+1); % minimum power obtainable at N+2 if I turn everything on at N+1 is given by <s, rho(t+1)>
    Objective_bottom =  s_bottom' * rho(:, N+1);
    options = sdpsettings('verbose',0,'solver','mosek');
    errors_top = optimize(Constraints, Objective_top,options);
    if(errors_top.problem ==1)
        disp("problem is infeasible!")
    end

    top_limit(N+2) = value(Objective_top);
    
    options = sdpsettings('verbose',0,'solver','mosek');
%     errors_bot = optimize(Constraints, Objective_bottom,options);
%     if(errors_bot.problem ==1)
%         disp("problem is infeasible!")
%     end

    bottom_limit(N+2) = 0; % value(Objective_bottom);
end

top_limit

end

