function trimmedMatrix = trim3DMatrix(matrix, desiredSize)
    % Get the current size of the matrix
    currentSize = size(matrix);

    if length(currentSize) == 2
        desiredSize = desiredSize(1:2);
    end
    
    % Determine the amount to trim from each dimension
    trimSize = (currentSize - desiredSize) / 2;
    
    % Calculate the indices for trimming
    startIndices = trimSize + 1;
    endIndices = currentSize - trimSize;
    
    % Trim the matrix
    if length(currentSize) == 2
        trimmedMatrix = matrix(startIndices(1):endIndices(1), startIndices(2):endIndices(2));
    else
        trimmedMatrix = matrix(startIndices(1):endIndices(1), startIndices(2):endIndices(2), startIndices(3):endIndices(3));
    end
end