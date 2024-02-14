function [state_results,input_results,cost_results] = solveModelFree(A,B,Lambda,x0,iters)
% solveModelFree Simulates the evolution of a system under a model-free recommendation system.
%
% Inputs:
%   A : The system dynamics matrix
%   B : The control input matrix
%   Lambda : Bias matrix
%   x0 : The initial state vector
%   iters : The number of iterations for which the system evolution is computed.
%
% Outputs:
%   state_results : A matrix where each column represents the state vector of the system at each iteration.
%   input_results : A row vector where each element represents the computed control input at each iteration.
%   cost_results : A row vector where each element represents the cost associated with the state and control input at each iteration.

    num_users = size(A,1);
    x_t = x0; % initial state

    % initialize output variables
    state_results = zeros(num_users,iters+1);
    state_results(:,1) = x_t;
    input_results = zeros(1,iters+1);
    cost_results = zeros(1,iters+1);
    
    for i=1:iters
        % Calculate input
        u = mean(x_t);

        % Calculate cost
        cost_results(1,i) = (x_t - ones(num_users,1)*u)'*(x_t - ones(num_users,1)*u);

        % Calculate state update
        x_t = A*x_t + B*u(1) + Lambda*x0;

        % Record results
        state_results(:,i+1) = x_t;
        input_results(1,i) = u;
    end

    u = mean(x_t);
    input_results(1,end) = u;
    cost_results(1,end) = (x_t - ones(num_users,1)*u)'*(x_t - ones(num_users,1)*u);
end