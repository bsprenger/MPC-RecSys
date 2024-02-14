%% base_example.m
% This script is a simple example of how the recommendation systems work
% for a given network. The MPC, Model-Free, and uncontrolled systems are
% simulated and plotted.

clear
clc
close all


%% Define social network
% The network with the influence of the recommendation system is based on 
% an extension of Friedkin-Johnsen opinion dynamics - see paper for more
% details.

num_users = 20;

sparsity_factor = 0.5; % pct of possible connections in network
max_connections = num_users * (num_users + 1);
connections = round(0.5*max_connections);

% Network adjacency matrix
W = generateSparseRowStochasticMatrix(num_users,num_users+1,connections);

% Biases
Lambda = diag(rand([num_users 1]));
    
% Dynamics matrices
A = (eye(num_users) - Lambda)*W(:,1:end-1);
B = (eye(num_users) - Lambda)*W(:,end);

% Initial opinion
x0 = rand([num_users 1]);


%% Simulate the system with MPC and Model-Free Rec. Systems
T = 20; % horizon for MPC
iters = 50; % number of simulation steps

[mpc_state,mpc_input,mpc_cost] = solveMPC(A,B,Lambda,x0,T,iters);
[mf_state,mf_input,mf_cost] = solveModelFree(A,B,Lambda,x0,iters);


%% Simulate network without recommendation system
% We can remove the column of the adjacency matrix associated with the
% recommendation system and rescale so that it is still row-stochastic.
% Then, we simulate the uncontrolled state evolution

% remove recommendation system's column from adj. matrix
W_new = W(:,1:end-1);

% calculate the sum of each row without the rec sys
W_new_sums = sum(W_new,2);

% re-scale each row so that it is row-stochastic
for i=1:num_users
    W_new(i,:) = W_new(i,:)/W_new_sums(i);
end

% initialize results variable
uncontrolled_state = zeros(num_users,iters+1);

% same initial condition
x_t = x0;

% calculate the state evolution
uncontrolled_state(:,1) = x_t;
for i=1:iters
    % Standard Friedkin-Johnsen model update step
    x_t = (eye(num_users)-Lambda)*W_new*x_t + Lambda*x0;
    
    % Log data
    uncontrolled_state(:,i+1) = x_t;
end


%% Plots
% configure plots to conform to IEEE standards for journal articles
figure_configuration_IEEE_standard

figure;
hold on;

% It may be necessary to set font to Times/Times New Roman, if it is not
% the default in Matlab. Examples below of some ways to do this:
% set(gca, 'FontName', 'Times')
% set(groot,'defaultAxesFontName','Times New Roman')
% set(gca, 'FontName', 'Times New Roman')

% Plot the state evolution of MPC and MF for the first user only
% This is to make it easier to generate the legend
h_mpc_state = plot(0:iters, mpc_state(1,:), 'Color', [1 0 0 1], 'LineWidth', 1.5);
h_mf_state = plot(0:iters, mf_state(1,:), 'Color', [0 0 1 1], 'LineWidth', 1.5);
h_uncontrolled_state = plot(0:iters, uncontrolled_state(1,:), 'Color', [0.6 0.6 0.6 1], 'LineStyle', '--', 'LineWidth', 1.5);

% plot remaining
plot(0:iters,mpc_state(2:,:)','Color',[1 0 0 1],'LineWidth',1.5);
plot(0:iters,mf_state(2:,:)','Color', [0 0 1 1],'LineWidth',1.5);
plot(0:iters,uncontrolled_state(2:,:)','Color', [0.6 0.6 0.6 1],'LineStyle','--','LineWidth',1.5);

% plot recommendations
h_mpc_input = plot(0:iters,mpc_input,'Color',[1 0 0 1],'LineStyle','-');
h_mf_input = plot(0:iters,mf_input,'Color',[0 0 1 1],'LineStyle','-');

xlabel('Update Step (t)')
ylabel('Opinion')

% Generate legend
legendEntries = {'User Opinions (MPC)', 'User Opinions (Model-Free)', 'User Opinions (Uncontrolled)', 'Recommendation (MPC)', 'Recommendation (Model-Free)'};
legend([h_mpc_state, h_mf_state, h_uncontrolled_state, h_mpc_input, h_mf_input], legendEntries, 'Location', 'northeast');

% plot the cost at each step
figure;
hold on;

plot(0:iters,mpc_cost,'Color',[1 0 0 1])
plot(0:iters,mf_cost,'Color',[0 0 1 1]);

legend({'MPC','Model-Free'},'Location','northeast')
xlabel("Update Step (t)")
ylabel("Step Cost")