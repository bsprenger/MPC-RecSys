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
Sigma = 0.02*eye(12);
Gamma = 0.01*eye(6); % must be > 0
Upsilon = zeros(12,6);

%% Setup LQG
T = 30;

% from 0 to T
R_matrices = zeros(12,12,T+1);
R_matrices(:,:,1) = R0_tilde;
for k=2:T+1
    R_matrices(:,:,k) = Sigma + A_tilde*R_matrices(:,:,k-1)*A_tilde' - ...
    (Upsilon + A_tilde*R_matrices(:,:,k-1)*C_tilde')*...
    ((Gamma + C_tilde*R_matrices(:,:,k-1)*C_tilde')\(Upsilon + A_tilde*R_matrices(:,:,k-1)*C_tilde')');
end

