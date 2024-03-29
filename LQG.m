clear
close all
clc

%% Initialize the problem
% create the adjacency matrix - column 7 corresponds to users' trust in rec. sys.
W = [0, 0.354, 0, 0, 0.132, 0, 0.514;
    0, 0, 0.395, 0.605, 0, 0, 0;
    0.333, 0.282, 0.084, 0, 0.071, 0.23, 0;
    0.3, 0, 0, 0.372, 0, 0.328, 0;
    0, 0, 1, 0, 0, 0, 0;
    0, 0, 0.102, 0, 0, 0.069, 0.829];

% assume known biases for now - for Kalman filter you always assume the
% dynamics are known, so since this is part of the dynamics, this is a fair
% assumption
Lambda = diag([0.333, 0.295, 0.152, 0.005, 0.319, 0.108]);

% construct A and B (see the paper for derivation)
A = (eye(6) - Lambda)*W(1:6,1:6);
B = (eye(6) - Lambda)*W(1:6,7);

% initial opinions
mu0 = [0.79; 0.68; 0.11; 0.1; 0.92; 0.02]; % mean
R0 = 0.01*eye(6); % variance matrix

%% Reformulate dynamics
A_tilde = [A, eye(6);
           zeros(6), eye(6)];

B_tilde = [B;zeros(6,1)];

C_tilde = [eye(6),zeros(6)];

mu0_tilde = [mu0;Lambda*mu0];
R0_tilde = [eye(6);Lambda]*R0*[eye(6),Lambda];

%% Setup noise variables
% noise is zero-mean by assumption
Sigma = 0.00005*eye(12); % dynamics noise covariance
Gamma = 0.01*eye(6); % measurement noise covariance, must be > 0
Upsilon = zeros(12,6); % covariance between dynamics and measurement

%% Setup LQG
% horizon
T = 30;

% cost matrices
M = eye(12);
N = 6;
S = -[ones(6,1);zeros(6,1)];

% ------------ Controller Matrices ------------ %
% from 1 to T i.e. of length T
P_matrices = zeros(12,12,T);
P_matrices(:,:,T) = zeros(12,12);
for k=T-1:-1:1
    P_matrices(:,:,k) = M + A_tilde'*P_matrices(:,:,k+1)*A_tilde - (S' + B_tilde'*P_matrices(:,:,k+1)*A_tilde)'*((N+B_tilde'*P_matrices(:,:,k+1)*B_tilde)\(S' + B_tilde'*P_matrices(:,:,k+1)*A_tilde));
end

% from 0 to T-1 i.e. of length T
K_matrices = zeros(1,12,T);
for k=1:T
    K_matrices(:,:,k) = -(N + B_tilde'*P_matrices(:,:,k)*B_tilde)\(S' + B_tilde'*P_matrices(:,:,k)*A_tilde);
end

% ------------ Estimator Matrices ------------ %
% from 0 to T-1 i.e. of length T
R_matrices = zeros(12,12,T);
R_matrices(:,:,1) = R0_tilde;
for k=2:T
    R_matrices(:,:,k) = Sigma + A_tilde*R_matrices(:,:,k-1)*A_tilde' - ...
    (Upsilon + A_tilde*R_matrices(:,:,k-1)*C_tilde')*...
    ((Gamma + C_tilde*R_matrices(:,:,k-1)*C_tilde')\(Upsilon + A_tilde*R_matrices(:,:,k-1)*C_tilde')');
end

% H from 1 to T i.e. of length T
H_matrices = zeros(12,6,T);
for k=1:T
    H_matrices(:,:,k) = (Upsilon + A_tilde*R_matrices(:,:,k)*C_tilde')/(Gamma + C_tilde*R_matrices(:,:,k)*C_tilde');
end

%% Simulate system
% Initialize state and estimation vectors
x = zeros(12,T+1); % States from 0 to T i.e. of length T+1
rng(3);
x(:,1) = mvnrnd(mu0_tilde,R_matrices(:,:,1),1); % Initial state 
y = zeros(6,T); % measurements from 1 to T i.e. length T
u = zeros(1,T); % Control input from 0 to T-1 i.e. length T
mu_hat = zeros(12,T+1); % State estimate from 0 to T i.e. of length T+1
mu_hat(:,1) = mu0_tilde; % Initial state estimate

% Simulation loop
for k = 1:T
    % Compute the control input based on the state estimate
    u(k) = K_matrices(:,:,k)*mu_hat(:,k);

    % Generate the system noise and measurement noise
    v = sqrt(Sigma)*[randn(6,1);zeros(6,1)]; % System noise
    w = sqrt(Gamma)*randn(6,1); % Measurement noise

    % Obtain the measurement
    y(:,k) = C_tilde*x(:,k) + w;

    % Update the state estimate
    mu_hat(:,k+1) = A_tilde*mu_hat(:,k) + B_tilde*u(k) + H_matrices(:,:,k)*(y(:,k) - C_tilde*mu_hat(:,k));
    
    % Update the actual state
    x(:,k+1) = A_tilde*x(:,k) + B_tilde*u(k) + v;
end

% Extract the actual states and estimates for plotting
actual_states = x(1:6,:);
estimated_states = mu_hat(1:6,:);

% Plot the actual and estimated states
time = 0:T;
figure;
plot(time, actual_states(1:6,:),'b',time, estimated_states(1:6,:), 'r--');
xlabel('Time');
ylabel('Opinion');
title("Control of Opinion Dynamics with LQG")
legend('Actual Opinions','','','','','','Estimated Opinions')
