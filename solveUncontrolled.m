function state_results = solveUncontrolled(W,Lambda,x0,iters)
% solveUncontrolled Simulates the evolution of a Friedkin-Johnsen system
% (without a recommendation system).
%
% Inputs:
%   W : The system dynamics matrix
%   Lambda : Bias matrix
%   x0 : The initial state vector
%   iters : The number of iterations for which the system evolution is computed.
%
% Outputs:
%   state_results : A matrix where each column represents the state vector of the system at each iteration.

    num_users = size(W,1);
    x_t = x0; % initial state

    % Preallocate memory for state_results matrix
    state_results = zeros(num_users, iters+1);
    state_results(:,1) = x_t;
    
    % Iterate over the specified number of iterations
    for i = 1:iters
        % Calculate state update
        x_t = (eye(num_users)-Lambda)*W*x_t + Lambda*x0;

        % Record results
        state_results(:,i+1) = x_t;
    end
end

