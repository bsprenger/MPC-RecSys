function [state_results,input_results,cost_results] = solveModelFree(A,B,Lambda,x0,iters)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    num_users = size(A,1);
    x_t = x0;
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