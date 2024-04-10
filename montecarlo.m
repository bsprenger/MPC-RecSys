%% Monte Carlo Simulation for Model-free vs. Model-based Approaches
clear
close all
clc

%% Parameters
num_trials = 1000; % Number of Monte Carlo trials
T = 50; % Time horizon
iters = 50; % Number of iterations
num_users = 20;
convergence_threshold = 1e-3;

%% Initialize storage structures
results = struct();

%% Run trials
for trial=1:num_trials
    fprintf("Trial %i\n",trial);

    % Determine sparsity level based on the current trial number
    if trial <= 250
        sparsity_factor = 1.00; % 100% connections
    elseif trial <= 500
        sparsity_factor = 0.75; % 75% connections
    elseif trial <= 750
        sparsity_factor = 0.50; % 50% connections
    else
        sparsity_factor = 0.25; % 25% connections
    end

    %% Generate dynamics
    max_connections = num_users * (num_users + 1);
    num_connections = round(max_connections * sparsity_factor);

    % Initialize W
    W = [];
    
    % Check for at least one non-zero entry in the last column
    hasNonZeroInLastColumn = false; % Initial condition
    
    while ~hasNonZeroInLastColumn
        W = generateSparseRowStochasticMatrix(num_users, num_users + 1, num_connections);
        
        % Check if the last column has at least one non-zero entry
        if any(W(:, end) ~= 0)
            hasNonZeroInLastColumn = true;
        end
    end

    Lambda = diag(rand([num_users 1]));
    
    A = (eye(num_users) - Lambda)*W(:,1:end-1);
    B = (eye(num_users) - Lambda)*W(:,end);
    
    x0 = rand([num_users 1]);

    %% Solve Recommendation Systems
    [mpc_state,mpc_input,mpc_cost] = solveMPC(A,B,Lambda,x0,T,iters);
    [mf_state,mf_input,mf_cost] = solveModelFree(A,B,Lambda,x0,iters);

    ss_mpc = mpc_state(:, end);
    ss_mf = mf_state(:, end);

    mpc_converge_step = find_convergence(mpc_state,convergence_threshold);
    mf_converge_step = find_convergence(mf_state,convergence_threshold);

    %% Calculate uncontrolled steady-state
    % remove rec sys column
    W_new = W(:,1:end-1);
    
    % calculate the sum of each row without the rec sys
    W_new_sums = sum(W_new,2);
    
    % re-scale each row so that it is row-stochastic
    for i=1:num_users
        W_new(i,:) = W_new(i,:)/W_new_sums(i);
    end
    ss_uncontrolled = (eye(num_users) - (eye(num_users) - Lambda)*W_new)\Lambda*x0;
    
    uncontrolled_state = solveUncontrolled(W_new,Lambda,x0,iters);
    uncontrolled_converge_step = find_convergence(uncontrolled_state,convergence_threshold);

    %% Calculate statistics
    % Transient cost performance comparison
    transient_mpc_cost = sum(mpc_cost);
    transient_mf_cost = sum(mf_cost);
    transient_pct_mpc_improvement = ((transient_mf_cost - transient_mpc_cost) / transient_mf_cost) * 100;
    
    % Steady-state cost performance comparison
    ss_mpc_cost = mpc_cost(end);
    ss_mf_cost = mf_cost(end);
    % positive value means MPC performs better (lower cost)
    ss_pct_mpc_improvement = ((ss_mf_cost - ss_mpc_cost) / ss_mf_cost) * 100;

    % Percentage shift per user, compared to uncontrolled
    pct_shift_mpc = abs(((ss_mpc - ss_uncontrolled)./ss_uncontrolled) * 100);
    pct_shift_mf = abs(((ss_mf - ss_uncontrolled)./ss_uncontrolled) * 100);
    
    
    %% Save results
    results(trial).transient_pct_mpc_improvement = transient_pct_mpc_improvement;
    results(trial).ss_pct_mpc_improvement = ss_pct_mpc_improvement;
    results(trial).pct_shift_mpc = pct_shift_mpc;
    results(trial).pct_shift_mf = pct_shift_mf;

    results(trial).mpc_converge_step = mpc_converge_step;
    results(trial).mf_converge_step = mf_converge_step;
    results(trial).uncontrolled_converge_step = uncontrolled_converge_step;
end

save('monte_carlo_data.mat','results');

%% Functions

function convergence_index = find_convergence(matrix, epsilon)
    % Compute the differences between consecutive columns
    column_diffs = diff(matrix, 1, 2);
    
    % Compute the norms of these differences
    norms = sqrt(sum(column_diffs.^2, 1));
    
    % Find the index where the norm is less than epsilon
    convergence_index = find(norms < epsilon, 1) + 1; % +1 to correspond to column of matrix, not column of norms
    
    % If convergence is not found, set the index to NaN
    if isempty(convergence_index)
        convergence_index = NaN;
    end
end