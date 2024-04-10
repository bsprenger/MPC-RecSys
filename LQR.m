clear
close all
clc

%% Initialize the problem
% create the adjacency matrix - column 7 corresponds to users' trust in rec. sys.
W = [0, 0.041, 0, 0.397, 0, 0, 0.562;
    0, 0.191, 0, 0, 0, 0.011, 0.798;
    0, 0, 0, 0, 0, 0.224, 0.776;
    1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 1, 0, 0, 0];

% biases
Lambda = diag([0.011, 0.001, 0.092, 0.064, 1.000, 0.055]);

% construct A and B (see the paper for derivation)
A = (eye(6) - Lambda)*W(1:6,1:6);
B = (eye(6) - Lambda)*W(1:6,7);

% initial opinions
x0 = [0.67; 0.74; 0.83; 0.68; 0.; 0.59];

%% Reformulate dynamics
A_tilde = [A, eye(6);
           zeros(6), eye(6)];

B_tilde = [B;zeros(6,1)];

x0_tilde = [x0;Lambda*x0];


%% Setup LQR
M = eye(12);
N = 6;
S = -[ones(6,1);zeros(6,1)];

T = 30;

% Calculate Riccati recursion
P_matrices = zeros(12,12,T+1);
P_matrices(:,:,T+1) = zeros(12,12);
for k=T:-1:1
    P_matrices(:,:,k) = M + A_tilde'*P_matrices(:,:,k+1)*A_tilde - (S' + B_tilde'*P_matrices(:,:,k+1)*A_tilde)'*((N+B_tilde'*P_matrices(:,:,k+1)*B_tilde)\(S' + B_tilde'*P_matrices(:,:,k+1)*A_tilde));
end

% Calculate gains from 0 to T-1
K_matrices = zeros(1,12,T);
for k=1:T
    K_matrices(:,:,k) = -(N + B_tilde'*P_matrices(:,:,k+1)*B_tilde)\(S' + B_tilde'*P_matrices(:,:,k+1)*A_tilde);
end

%% Simulate system
% Initialize state
x = x0_tilde;

% Initialize arrays to store state and control input history
x_history = zeros(12, T+1);
u_history = zeros(1, T);

% Set initial state in history
x_history(:,1) = x;

% Simulate system
for k = 1:T
    % Compute control input using the current state and the gain matrix
    u = K_matrices(:,:,k) * x;
    u_history(:,k) = u;

    % Update state using system dynamics
    x = A_tilde * x + B_tilde * u;
    
    % Store the state
    x_history(:,k+1) = x;
end

plot(0:T, x_history(1:6,:),'b');
xlabel('Time Step');
ylabel('Opinions');
title("Control of Opinion Dynamics with LQR")

