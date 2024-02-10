clear
close all
clc

load('monte_carlo_data.mat');

% Initialize counters and accumulators
num_trials = length(results);
num_negative_total_improvement = 0;
num_negative_steady_state_improvement = 0;
avg_steady_state_percentage_improvement_mpc_over_mf = 0;
avg_opinion_shift_mpc = 0;
avg_opinion_shift_mf = 0;
avg_percentage_shift_mpc_over_mf = 0;

for trial = 1:num_trials
    % Check for negative values in total improvement MPC over MF
    if results(trial).total_improvement_mpc_over_mf < 0
        num_negative_total_improvement = num_negative_total_improvement + 1;
    end
    
    % Assuming steady_state_improvement_mpc_over_mf is calculated and stored
    % Check for negative values in steady state improvement MPC over MF
    if results(trial).steady_state_improvement_mpc_over_mf < 0
        num_negative_steady_state_improvement = num_negative_steady_state_improvement + 1;
    end
    
    % Accumulate values for average calculations
    avg_steady_state_percentage_improvement_mpc_over_mf = avg_steady_state_percentage_improvement_mpc_over_mf + results(trial).steady_state_percentage_improvement_mpc_over_mf;
    avg_opinion_shift_mpc = avg_opinion_shift_mpc + results(trial).opinion_shift_mpc;
    avg_opinion_shift_mf = avg_opinion_shift_mf + results(trial).opinion_shift_mf;
    avg_percentage_shift_mpc_over_mf = avg_percentage_shift_mpc_over_mf + results(trial).percentage_shift_mpc_over_mf;
end

% Calculate averages
avg_steady_state_percentage_improvement_mpc_over_mf = avg_steady_state_percentage_improvement_mpc_over_mf / num_trials;
avg_opinion_shift_mpc = avg_opinion_shift_mpc / num_trials;
avg_opinion_shift_mf = avg_opinion_shift_mf / num_trials;
avg_percentage_shift_mpc_over_mf = avg_percentage_shift_mpc_over_mf / num_trials;

% Display results
fprintf('Number of trials with negative total improvement: %d\n', num_negative_total_improvement);
fprintf('Number of trials with negative steady state improvement: %d\n', num_negative_steady_state_improvement);
fprintf('Average steady state percentage improvement MPC over MF: %.2f%%\n', avg_steady_state_percentage_improvement_mpc_over_mf);
fprintf('Average opinion shift MPC: %.2f\n', avg_opinion_shift_mpc);
fprintf('Average opinion shift MPC: %.2f\n', avg_opinion_shift_mf);
fprintf('Average percentage shift MPC over MF: %.2f%%\n', avg_percentage_shift_mpc_over_mf);
