%% Monte Carlo Simulation for Model-free vs. Model-based Approaches
clear
close all
clc

%% Parameters
num_trials = 1000; % Number of Monte Carlo trials
T = 20; % Time horizon
iters = 50; % Number of iterations

%% Initialize storage structures
results = struct();

%% Run trials
for trial=1:num_trials
    disp(trial);
    %% Generate dynamics with random number of users
    num_users = randi([10, 20]);
    
    W = generateRowStochasticMatrix(num_users,num_users+1);
    Lambda = diag(rand([num_users 1]));
    
    A = (eye(num_users) - Lambda)*W(:,1:end-1);
    B = (eye(num_users) - Lambda)*W(:,end);
    
    x0 = rand([num_users 1]);
    
    %% Solve Recommendation Systems
    [mpc_state,mpc_input,mpc_cost] = solveMPC(A,B,Lambda,x0,T,iters);
    [mf_state,mf_input,mf_cost] = solveModelFree(A,B,Lambda,x0,iters);

    controlled_steady_state_mpc = mpc_state(:, end);
    controlled_steady_state_mf = mf_state(:, end);
    
    %% Calculate uncontrolled steady-state
    % remove rec sys column
    W_new = W(:,1:end-1);
    
    % calculate the sum of each row without the rec sys
    W_new_sums = sum(W_new,2);
    
    % re-scale each row so that it is row-stochastic
    for i=1:num_users
        W_new(i,:) = W_new(i,:)/W_new_sums(i);
    end
    uncontrolled_steady_state = (eye(num_users) - (eye(num_users) - Lambda)*W_new)\Lambda;
    
    %% Calculate statistics
    % Total cost of each approach over the horizon
    total_mpc_cost = sum(mpc_cost);
    total_mf_cost = sum(mf_cost);
    total_improvement_mpc_over_mf = total_mf_cost - total_mpc_cost;
    total_percentage_improvement_mpc_over_mf = ((total_mf_cost - total_mpc_cost) / total_mf_cost) * 100;
    
    % steady-state cost
    steady_state_mpc_cost = mpc_cost(end);
    steady_state_mf_cost = mf_cost(end);
    % positive value means MPC performs better (lower cost)
    steady_state_improvement_mpc_over_mf = steady_state_mf_cost - steady_state_mpc_cost;
    steady_state_percentage_improvement_mpc_over_mf = ((steady_state_mf_cost - steady_state_mpc_cost) / steady_state_mf_cost) * 100;
    
    % final opinion shift relative to the uncontrolled system 
    opinion_shift_mpc = norm(controlled_steady_state_mpc - uncontrolled_steady_state);
    opinion_shift_mf = norm(controlled_steady_state_mf - uncontrolled_steady_state);
    percentage_shift_mpc_over_mf = ((opinion_shift_mpc - opinion_shift_mf) / opinion_shift_mf) * 100;
    
    %% Save results
    results(trial).mpc_state = mpc_state;
    results(trial).mpc_input = mpc_input;
    results(trial).mpc_cost = mpc_cost;
    results(trial).mf_state = mf_state;
    results(trial).mf_input = mf_input;
    results(trial).mf_cost = mf_cost;
    results(trial).A = A;
    results(trial).B = B;
    results(trial).Lambda = Lambda;
    results(trial).x0 = x0;
    results(trial).total_mpc_cost = total_mpc_cost;
    results(trial).total_mf_cost = total_mf_cost;
    results(trial).total_improvement_mpc_over_mf = total_improvement_mpc_over_mf;
    results(trial).total_percentage_improvement_mpc_over_mf = total_percentage_improvement_mpc_over_mf;
    results(trial).steady_state_mpc_cost = steady_state_mpc_cost;
    results(trial).steady_state_mf_cost = steady_state_mf_cost;
    results(trial).steady_state_improvement_mpc_over_mf = steady_state_improvement_mpc_over_mf;
    results(trial).steady_state_percentage_improvement_mpc_over_mf = steady_state_percentage_improvement_mpc_over_mf;
    results(trial).opinion_shift_mpc = opinion_shift_mpc;
    results(trial).opinion_shift_mf = opinion_shift_mf;
    results(trial).percentage_shift_mpc_over_mf = percentage_shift_mpc_over_mf;

end

save('monte_carlo_data.mat','results');