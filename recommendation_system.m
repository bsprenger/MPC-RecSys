clear
clc
close all

W = [0, 0.041, 0, 0.397, 0, 0, 0.562;
    0, 0.191, 0, 0, 0, 0.011, 0.798;
    0, 0, 0, 0, 0, 0.224, 0.776;
    1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 1, 0, 0, 0];

Lambda = diag([0.011, 0.001, 0.092, 0.064, 1.000, 0.055]);

A = (eye(6) - Lambda)*W(1:6,1:6); % check this
B = (eye(6) - Lambda)*W(1:6,7);

x0 = [0.67; 0.74; 0.83; 0.68; 0.; 0.59];

T = 10;
iters = 30;

% Initialize the matrices
S_x = zeros(size(A, 1) * (T + 1), size(A, 1));
S_u = zeros(size(A, 1) * (T + 1), size(B, 2) * (T+1));
C = zeros(size(A, 1) * (T+1), 1);

% Construct S_x matrix
for i = 1:T+1
    S_x(((i - 1) * size(A, 1) + 1):(i * size(A, 1)), :) = A^(i-1);
end

% Construct S_u matrix
for i = 1:T
    for j = 1:i
        S_u(((i - 1) * size(A, 1) + 1 + 6):(i * size(A, 1) + 6), ((j - 1) * size(B, 2) + 1):(j * size(B, 2))) = A^(i-j)*B;
    end
end

% Construct C matrix
current_c = zeros(size(A, 1), 1);
for i = 2:T+1
    current_c = A*current_c + Lambda*x0;
    C(((i - 1) * size(A, 1) + 1):(i * size(A, 1)), 1) = current_c;
end

%% Solve QP
x_t = x0;
state_results = zeros(6,iters+1);
state_results(:,1) = x_t;
input_results = zeros(1,iters);
for i = 1:iters
    H = 2*(S_u - kron(eye(T+1),ones(6,1)))'*(S_u - kron(eye(T+1),ones(6,1)));
    f = 2*(S_u - kron(eye(T+1),ones(6,1)))'*(S_x*x_t + C);
    
    lb = zeros(T+1,1);
    ub = ones(T+1,1);
    
    u = quadprog(H,f,[],[],[],[],lb,ub);
    x_t = A*x_t + B*u(1) + Lambda*x0;

    state_results(:,i+1) = x_t;
    input_results(1,i) = u(1);
end

figure_configuration_IEEE_standard
plot(0:iters,state_results')

figure
plot(input_results)