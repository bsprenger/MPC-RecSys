%% base_example.m
% This script is a simple example of how the recommendation systems work
% for a given network. The MPC, Model-Free, and uncontrolled systems are
% simulated and plotted.

clear
clc
close all

%% Initialize the problem
% create the adjacency matrix - column 7 corresponds to users' trust in rec. sys.
W = [0, 0.354, 0, 0, 0.132, 0, 0.514;
    0, 0, 0.395, 0.605, 0, 0, 0;
    0.333, 0.282, 0.084, 0, 0.071, 0.23, 0;
    0.3, 0, 0, 0.372, 0, 0.328, 0;
    0, 0, 1, 0, 0, 0, 0;
    0, 0, 0.102, 0, 0, 0.069, 0.829];

% biases
Lambda = diag([0.333, 0.295, 0.152, 0.005, 0.319, 0.108]);

% construct A and B (see the paper for derivation)
A = (eye(6) - Lambda)*W(1:6,1:6);
B = (eye(6) - Lambda)*W(1:6,7);

% initial opinions
x0 = [0.79; 0.68; 0.11; 0.1; 0.92; 0.02];

T = 10; % horizon
iters = 30; % number of simulation steps

%% Solve Recommendation Systems
[mpc_state,mpc_input,mpc_cost] = solveMPC(A,B,Lambda,x0,T,iters);
[mf_state,mf_input,mf_cost] = solveModelFree(A,B,Lambda,x0,iters);

%% Calculate uncontrolled (open-loop) state evolution
% First, we remove the recommendation system from the graph and re-scale so
% that matrix is row-stochastic

% remove rec sys column
W_new = W(:,1:6);

% calculate the sum of each row without the rec sys
W_new_sums = sum(W_new,2);

% re-scale each row so that it is row-stochastic
for i=1:6
    W_new(i,:) = W_new(i,:)/W_new_sums(i);
end

% initialize results variable
uncontrolled_results = zeros(6,iters+1);

% initial condition
x_t = x0;

% calculate the state evolution
uncontrolled_results(:,1) = x_t;
for i=1:iters
    x_t = (eye(6)-Lambda)*W_new*x_t + Lambda*x0;
    uncontrolled_results(:,i+1) = x_t;
end


%% Plots
figure_configuration_IEEE_standard

figure;
% For IEEE:
set(gca, 'FontName', 'Times')
set(groot,'defaultAxesFontName','Times New Roman')

% plot the states
plot(0:iters,mpc_state','Color',[1 0 0 1],'LineWidth',1.5);
set(gca, 'FontName', 'Times New Roman')
hold on;
plot(0:iters,mf_state','Color', [0 0 1 1],'LineWidth',1.5);
plot(0:iters,uncontrolled_results','Color', [0.6 0.6 0.6 1],'LineStyle','--','LineWidth',1.5)

% plot recommendations
plot(0:iters,mpc_input,'Color',[1 0 0 1],'LineStyle','-')
plot(0:iters,mf_input,'Color',[0 0 1 1],'LineStyle','-');

xlabel('Update Step (t)','FontName','Times New Roman')
ylabel('Opinion')
legend({'User Opinions (MPC)','','','','','','User Opinions (Model-Free)',...
    '','','','','','User Opinions (Uncontrolled)','','','','','',...
    'Recommendation (MPC)','Recommendation (Model-Free)'},...
    'Location','northeast')

% plot the cost at each step
figure;
plot(0:iters,mpc_cost,'Color',[1 0 0 1])
hold on
plot(0:iters,mf_cost,'Color',[0 0 1 1]);
legend({'MPC','Model-Free'},'Location','northeast')
xlabel("Update Step (t)")
ylabel("Step Cost")