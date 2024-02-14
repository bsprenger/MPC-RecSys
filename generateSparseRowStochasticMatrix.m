function R = generateSparseRowStochasticMatrix(m, n, nonZeroEntries)
    % Generates an m x n sparse row stochastic matrix with a specified number of non-zero entries
    %
    % Parameters:
    % m - the number of rows
    % n - the number of columns
    % nonZeroEntries - total number of non-zero entries in the matrix
    %
    % Returns:
    % R - an m x n sparse row stochastic matrix

    % Check if inputs are valid
    if m <= 0 || n <= 0 || floor(m) ~= m || floor(n) ~= n || nonZeroEntries < m || nonZeroEntries > m*n
        error('Invalid inputs. Ensure positive integers, and nonZeroEntries is between m and m*n.');
    end
    
    % Initialize the matrix with zeros
    R = zeros(m, n);
    
    % Ensure at least one non-zero entry per row for row stochastic property
    for i = 1:m
        colIndex = randi(n-1); % Random column index for the non-zero entry -> can't have only value being rec sys (must trust another user too)
        R(i, colIndex) = rand(); % Assign a random value
    end
    
    % Fill the rest of the non-zero entries randomly across the matrix
    additionalEntries = nonZeroEntries - m; % Subtract the already filled entries
    while additionalEntries > 0
        rowIndex = randi(m);
        colIndex = randi(n);
        if R(rowIndex, colIndex) == 0 % Check if the position is already filled
            R(rowIndex, colIndex) = rand(); % Assign a random value
            additionalEntries = additionalEntries - 1;
        end
    end
    
    % Normalize each row to sum to 1
    rowSums = sum(R, 2); % Calculate the sum of elements in each row
    R = R ./ rowSums; % Divide each element by its row sum to normalize
end