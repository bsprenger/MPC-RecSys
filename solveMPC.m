function [state_results,input_results,cost_results] = solveMPC(A,B,Lambda,x0,T,iters)
% solveModelFree Simulates the evolution of a system under a MPC-based recommendation system.
%
% Inputs:
%   A : System dynamics matrix
%   B : Input matrix
%   Lambda : Bias matrix
%   x0 : Initial state vector
%   T : Prediction horizon for MPC
%   iters : Number of iterations to simulate
%
% Outputs:
%   state_results : A matrix where each column represents the state vector of the system at each iteration.
%   input_results : A row vector where each element represents the computed control input at each iteration.
%   cost_results : A row vector where each element represents the cost associated with the state and control input at each iteration.

%% Solve for optimal steady-state
num_users = size(A,1);

% x = (I - A)^(-1)(B*u + Lambda*x0);

% cost (x-u*ones)'*(x-u*ones)
% x'x - 2u*ones'*x + u*ones'*ones*u
% 0 = - ones'*x + ones'*ones*u
% 0 = -ones'* (I - A)^(-1)(B*u + Lambda*x0) + ones'*ones*u
% 0 = -ones'*(I-A)^(-1)*Bu - ones'*(I-A)^(-1)*Lambda*x0 + ones'*ones*u
% ones'*(I-A)^(-1)*Lambda*x0 = (ones'*ones - ones'*(I-A)^(-1)*B)u
% u = (ones'*(ones - (I-A)^(-1)*B))^(-1)*ones'*(I-A)^(-1)*Lambda*x0

u_s = (ones(num_users,1)'*(ones(num_users,1) - (eye(num_users) - A)\B))\...
    (ones(num_users,1)'*((eye(num_users) - A)\(Lambda*x0)));

% disp("x_s: ")
% disp(u_s);
x_s = (eye(num_users)-A)\(B*u_s + Lambda*x0);
% disp(x_s)

% double check that u is in [0,1] here, max/min it

% then calculate x? or leave as u?

% X = S_x*x0 + S_u*U + C
% we want to set X(end chunk) = x_s
% Choose selection matrix M_end
% x_s = M_end*S_x*x0 + M_end*S_u*U + M_end*C
% x_s - M_end*C - M_end*S_x*x0 = M_end*S_u*U simple equality constraint
%% Construct QP
    

    % Unroll the Friedkin-Johnsen dynamics over the horizon to create a big QP.
    % X = S_x*x0 + S_u*U + C, where X is the vector of the entire
    % state evolution, and U is the vector of all control inputs,
    % and C is a constant vector.
    % In this way, we remove X from our optimization and optimize only over U.
    S_x = zeros(size(A, 1) * (T + 1), size(A, 1));
    S_u = zeros(size(A, 1) * (T + 1), size(B, 2) * (T+1));
    C = zeros(size(A, 1) * (T+1), 1); % addition of constant vector is necessary to capture effect of Lambda * x0
    
    % Construct S_x matrix
    for i = 1:T+1
        S_x(((i - 1) * num_users + 1):(i * num_users), :) = A^(i-1);
    end
    
    % Construct S_u matrix
    for i = 1:T
        for j = 1:i
            S_u(((i - 1) * num_users + 1 + num_users):(i * num_users + num_users), ((j - 1) * size(B, 2) + 1):(j * size(B, 2))) = A^(i-j)*B;
        end
    end
    
    % Construct C matrix
    current_c = zeros(num_users, 1);
    for i = 2:T+1
        current_c = A*current_c + Lambda*x0;
        C(((i - 1) * size(A, 1) + 1):(i * size(A, 1)), 1) = current_c;
    end
    
    M = zeros(num_users,num_users*(T+1));
    % disp(size(M))
    M(:,end-num_users+1:end) = eye(num_users);
    % disp(M)
    %% Solve MPC Recommendation System
    % initial condition
    x_t = x0;
    
    % create variables to store results
    state_results = zeros(num_users,iters+1);
    state_results(:,1) = x_t;
    input_results = zeros(1,iters+1);
    cost_results = zeros(1,iters+1);
    
    for i = 1:iters
        % Construct QP: see Matlab's quadprog
        % for notation
        H = 2*(S_u - kron(eye(T+1),ones(num_users,1)))'*(S_u - kron(eye(T+1),ones(num_users,1)));
        H = (H + H') / 2; % correct minor numerical asymmetry errors
        f = 2*(S_u - kron(eye(T+1),ones(num_users,1)))'*(S_x*x_t + C);

        Aeq = M*S_u;
        beq = x_s - M*C - M*S_x*x0;

        % constraints
        lb = zeros(T+1,1);
        ub = ones(T+1,1);
        
        % solve for optimal inputs over horizon
        options = optimoptions('quadprog', 'Display', 'off');
        u = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
        % disp("iter:")
        % disp(i)
        % disp("U: ")
        % disp(u)
        % disp("done u")
        % calculate the cost for the FIRST optimal input
        cost_results(1,i) = (x_t - ones(num_users,1)*u(1))'*(x_t - ones(num_users,1)*u(1));
    
        % move system forward one step
        x_t = A*x_t + B*u(1) + Lambda*x0;
    
        % record results
        state_results(:,i+1) = x_t;
        input_results(1,i) = u(1);
    end
    
    % get the last input and store
    H = 2*(S_u - kron(eye(T+1),ones(num_users,1)))'*(S_u - kron(eye(T+1),ones(num_users,1)));
    H = (H + H') / 2; % correct minor numerical asymmetry errors
    f = 2*(S_u - kron(eye(T+1),ones(num_users,1)))'*(S_x*x_t + C);
    
    lb = zeros(T+1,1);
    ub = ones(T+1,1);
    
    u = quadprog(H,f,[],[],[],[],lb,ub,[],options);
    input_results(1,end) = u(1);
    % disp("x_end:")
    % disp(state_results(:,end));
    cost_results(1,end) = (x_t - ones(num_users,1)*u(1))'*(x_t - ones(num_users,1)*u(1));
end