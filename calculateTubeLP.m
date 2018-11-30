function [top_limit,bottom_limit] = calculateTubeLP(rho0, time_window, P, B)
%CALCULATE_TUBE Summary of this function goes here
%   Detailed explanation goes here

epsilon = 1*10^(-7);
r = size(P,1);
one = ones(r,1);
on = ones(r, 1);
on(r/2+1:end)=zeros(r/2,1);

rho_top(1) = 1; % nonsense bound at time 0
rho_top(2) = P(r,r/2) * rho0(r/2) +  rho0(r);

rho_bot(1) = 0; % nonsense bound at time 0
rho_bot(2) = P(1,r/2 + 1) * rho0(r/2 +1) +  rho0(1);


for N=1:time_window-1
    Constraints = [];
    Objective = 0;
    M = sdpvar(r,r,N,'full');
    rho = sdpvar(r,N+1, 'full');
    rho(:,1) = rho0;
    if(N==r-1)
        disp("odjwqodj")
    end
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
                        if (j==1 || j==r) 
                            %case when i == {1; r} (incontrollable
                            %states)
                            Constraints = [Constraints; M(i,j,k)==P(i,j)*rho(j,k)];
                        elseif (j==1 + r/2 || j==r/2)
                            if(i ~= j)
                                Constraints = [Constraints; M(i,j,k)>=P(i,j)*rho(j,k)]; %if I am too warm / too cold, I can only increase this probability
                            end
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
            on'*rho(:,k+1)<= 1-rho_top(k+1);
            on'*rho(:,k+1)>= rho_bot(k+1);
            ];
    end
    Objective_top = (-P(r,r/2) * rho(r/2, N+1)) -  (rho(r, N+1));
    Objective_bottom = (-P(1,r/2 +1) * rho(r/2 +1, N+1)) -  (rho(1, N+1));
    options = sdpsettings('verbose',0,'solver','mosek');
    errors_top = optimize(Constraints, Objective_top,options);
    if(errors_top.problem ==1)
        disp("problem is infeasible!")
    end
%     Power_rho_top = on'*value(rho);
%     Constraints_top = [Constraints; 
%         on'*rho == Power_rho_top];
%     Objective_top_2 = (P(r,r/2) * rho(r/2, N+1)) +  (rho(r, N+1)); 
%     errors_top_2 = optimize(Constraints, Objective_top_2,options);
%     if(errors_top_2.problem ==1)
%         disp("problem is infeasible!")
%     end
%     Constraints_top = [];
%     rho_top(N+2) = value(Objective_top_2);
    rho_top(N+2) = -value(Objective_top);
    
    options = sdpsettings('verbose',0,'solver','mosek');
    errors_bot = optimize(Constraints, Objective_bottom,options);
    if(errors_bot.problem ==1)
        disp("problem is infeasible!")
    end
%     Power_rho_bottom = on'*value(rho);
%     Constraints_bottom = [Constraints; 
%         on'*rho == Power_rho_bottom];
%     Objective_bottom_2 = (P(1,r/2 +1) * rho(r/2 +1, N+1)) + (rho(1, N+1));
%     errors_bot_2 = optimize(Constraints, Objective_bottom_2,options);
%     if(errors_bot_2.problem ==1)
%         disp("problem is infeasible!")
%     end
%     Constraints_bottom = [];
%     rho_bot(N+2) = value(Objective_bottom_2);
    rho_bot(N+2) = - value(Objective_bottom);
end


% rho_top =optimizer(Constraints,Objective_top,options,index_parameter,Objective_top);
% rho_bot =optimizer(Constraints,Objective_bot,options,index_parameter,Objective_bot);

for i=1:time_window+1
    top_limit(i) =  1- rho_top(i);
end
for i=1:time_window+1
    bottom_limit(i) = rho_bot(i);
end

end

