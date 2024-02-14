%% calculate_statistics.m
% Quick and dirty script to calculate the important statistics
% for the paper from Monte Carlo simulations

clear
close all
clc


%% Get dataset
load('monte_carlo_data.mat');


%% Calculate statistics
% Sims were chunked into batches of 250 based on network sparsity
disp('-------------------------------------------------------------------')
disp('100% FILLED:');
calculateStatistics(results(1:250))

disp('-------------------------------------------------------------------')
disp('75% FILLED:');
calculateStatistics(results(251:500))

disp('-------------------------------------------------------------------')
disp('50% FILLED:');
calculateStatistics(results(501:750))

disp('-------------------------------------------------------------------')
disp('25% FILLED:');
calculateStatistics(results(751:1000))
disp('-------------------------------------------------------------------')


%% Helper functions

function calculateStatistics(results)
    % Takes the struct array as input. Check Monte Carlo sim for details
    % Prints important stats to console

    % Initialize counters and accumulators
    num_trials = length(results);
    num_negative_transient_improvement = 0; % check if MF ever outperforms MPC in transient
    num_negative_steady_state_improvement = 0; % check if MF ever outperforms MPC in steady-state
    avg_ss_pct_mpc_improvement = 0; % 
    pct_shift_mpc = [];
    pct_shift_mf = [];
    
    for trial = 1:num_trials
        % Check for negative values in transient improvement MPC over MF
        % These would indicate MF outperforms MPC in transient
        if results(trial).transient_pct_mpc_improvement < 0
            num_negative_transient_improvement = num_negative_transient_improvement + 1;
        end
        
        % Check for negative values in steady state improvement MPC over MF
        % These would indicate MF outperforms MPC in transient
        if results(trial).ss_pct_mpc_improvement < 0
            num_negative_steady_state_improvement = num_negative_steady_state_improvement + 1;
        end
        
        % Accumulate values for average calculations
        avg_ss_pct_mpc_improvement = avg_ss_pct_mpc_improvement + results(trial).ss_pct_mpc_improvement;
        pct_shift_mpc = [pct_shift_mpc;results(trial).pct_shift_mpc];
        pct_shift_mf = [pct_shift_mf;results(trial).pct_shift_mf];
    end

    % Calculate averages and maxes
    avg_ss_pct_mpc_improvement = avg_ss_pct_mpc_improvement / num_trials;
    avg_pct_shift_mpc = mean(pct_shift_mpc);
    avg_pct_shift_mf = mean(pct_shift_mf);
    max_pct_shift_mpc = max(pct_shift_mpc);
    max_pct_shift_mf = max(pct_shift_mf);
    
    % Display results
    fprintf('Number of trials with negative total improvement: %d\n', num_negative_transient_improvement);
    fprintf('Number of trials with negative steady state improvement: %d\n', num_negative_steady_state_improvement);
    fprintf('Average steady state percentage improvement MPC over MF: %.5f%%\n', avg_ss_pct_mpc_improvement);
    fprintf('Average opinion shift MPC: %.5f%%\n', avg_pct_shift_mpc);
    fprintf('Average opinion shift MF: %.5f%%\n', avg_pct_shift_mf);
    fprintf('Max opinion shift MPC: %.5f%%\n', max_pct_shift_mpc);
    fprintf('Max opinion shift MF: %.5f%%\n', max_pct_shift_mf);
end
