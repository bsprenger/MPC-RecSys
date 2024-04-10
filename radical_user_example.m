%% radical_user_example.m
% This script is an example of a situation in which the two recommendation
% system approaches can vary significantly in both performance and in their
% influence on users.
% The scenario is of one radical user with a polar opinion who is unwilling
% to be influenced by other users or even the recommendation system. Note
% how both recommendation systems steer users towards this radical opinion,
% but especially the MPC, which strategically makes use of its knowledge of
% the dynamics to manipulate users to a more convenient opinion to cater
% to.

clear
clc
close all


%% Define social network
% Network is set up to be particularly vulnerable to being influenced by
% the recommendation system (low biases), except for one radical user (5th
% row) who is unwilling to be influenced by anyone.
% The network with the influence of the recommendation system is based on 
% an extension of Friedkin-Johnsen opinion dynamics - see paper for more
% details.

% Network adjacency matrix
W = [0, 0.041, 0, 0.397, 0, 0, 0.562;
    0, 0.191, 0, 0, 0, 0.011, 0.798;
    0, 0, 0, 0, 0, 0.224, 0.776;
    1, 0, 0, 0, 0, 0, 0;
    0, 0, 0.472, 0.171, 0, 0, 0.357;
    0, 0, 0, 1, 0, 0, 0];

% Biases
Lambda = diag([0.011, 0.001, 0.092, 0.064, 1.000, 0.055]);

% Dynamics matrices
A = (eye(6) - Lambda)*W(1:6,1:6);
B = (eye(6) - Lambda)*W(1:6,7);

% Initial opinion. Note the radical user (#5) with the extreme opinion
x0 = [0.67; 0.74; 0.83; 0.68; 0.; 0.59];


%% Simulate the system with MPC and Model-Free Rec. Systems
T = 50; % horizon for MPC
iters = 50; % number of simulation steps

[mpc_state,mpc_input,mpc_cost] = solveMPC(A,B,Lambda,x0,T,iters);
[mf_state,mf_input,mf_cost] = solveModelFree(A,B,Lambda,x0,iters);


%% Simulate network without recommendation system
% We can remove the column of the adjacency matrix associated with the
% recommendation system and rescale so that it is still row-stochastic.
% Then, we calculate the steady-state opinion from the classic
% Friedkin-Johnsen model steady-state solution.

% remove recommendation system's column from adj. matrix
W_new = W(:,1:end-1);
    
% calculate the sum of each row without the rec sys
W_new_sums = sum(W_new,2);

% re-scale each row so that it is row-stochastic
for i=1:6
    W_new(i,:) = W_new(i,:)/W_new_sums(i);
end

% calculate steady-state solution
ss_uncontrolled = (eye(6) - (eye(6) - Lambda)*W_new)\Lambda*x0;


%% Calculate statistics
% Transient cost performance comparison
transient_mpc_cost = sum(mpc_cost);
transient_mf_cost = sum(mf_cost);
% positive value of below means MPC performs better (lower cost)
transient_pct_mpc_improvement = ((transient_mf_cost - transient_mpc_cost) / transient_mf_cost) * 100;

% Steady-state cost performance comparison
ss_mpc_cost = mpc_cost(end);
ss_mf_cost = mf_cost(end);
% positive value of below means MPC performs better (lower cost)
ss_pct_mpc_improvement = ((ss_mf_cost - ss_mpc_cost) / ss_mf_cost) * 100;

% Percentage opinion shift per user, compared to uncontrolled
pct_shift_mpc = abs(((mpc_state(:, end) - ss_uncontrolled)./ss_uncontrolled) * 100);
pct_shift_mf = abs(((mf_state(:, end) - ss_uncontrolled)./ss_uncontrolled) * 100);
avg_pct_shift_mpc = mean(rmmissing(pct_shift_mpc));
avg_pct_shift_mf = mean(rmmissing(pct_shift_mf));


%% Plots
% configure plots to conform to IEEE standards for journal articles
figure_configuration_IEEE_standard

figure;
hold on;

% It may be necessary to set font to Times/Times New Roman, if it is not
% the default in Matlab. Examples below of some ways to do this:
% set(gca, 'FontName', 'Times')
% set(groot,'defaultAxesFontName','Times New Roman')

% plot the state evolution of MPC and MF for first user only
% This is to create legend easier
h_mpc_state = plot(0:iters,mpc_state(1,:),'Color',[1 0 0 0.7],'LineWidth',1.2);
h_mf_state = plot(0:iters,mf_state(1,:),'Color', [0 0 1 0.7],'LineWidth',1.2);

% plot remaining
plot(0:iters,mpc_state(2:end,:)','Color',[1 0 0 0.7],'LineWidth',1.2);
plot(0:iters,mf_state(2:end,:)','Color', [0 0 1 0.7],'LineWidth',1.2);

% highlight radical user opinion
h_radical = yline(0,'Color',[0 0.7 0],'LineWidth',4);

% plot recommendations (control input)
h_mpc_input = plot(0:iters,mpc_input,'Color',[1 0 0 1],'LineStyle','-');
h_mf_input = plot(0:iters,mf_input,'Color',[0 0 1 1],'LineStyle','-');

% update font name if other than Times is desired
xlabel('Update Step (t)')
ylabel('Opinion')
ylim([-0.05 1])

legendEntries = {'User Opinions (MPC)', 'User Opinions (Model-Free)', 'Radical User', 'Recommendation (MPC)', 'Recommendation (Model-Free)'};
legend([h_mpc_state, h_mf_state, h_radical, h_mpc_input, h_mf_input], legendEntries, 'Location', 'northeast');

% plot the cost at each step
figure;
hold on
plot(0:iters,mpc_cost,'Color',[1 0 0 1])
plot(0:iters,mf_cost,'Color',[0 0 1 1]);
legend({'MPC','Model-Free'},'Location','northeast')
xlabel("Update Step (t)")
ylabel("Step Cost")