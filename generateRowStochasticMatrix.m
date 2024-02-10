function R = generateRowStochasticMatrix(m, n)
    % Generates an m x n row stochastic matrix
    %
    % Parameters:
    % m - the number of rows
    % n - the number of columns
    %
    % Returns:
    % R - an m x n row stochastic matrix

    % Check if inputs are valid
    if m <= 0 || n <= 0 || floor(m) ~= m || floor(n) ~= n
        error('Inputs must be positive integers');
    end
    
    % Initialize the matrix with random numbers
    R = rand(m, n);
    
    % Normalize each row to sum to 1
    rowSums = sum(R, 2); % Calculate the sum of elements in each row
    R = R ./ rowSums; % Divide each element by its row sum to normalize
end
