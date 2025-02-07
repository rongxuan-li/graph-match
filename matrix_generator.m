clear, close all

% Define the size of the matrices
m = 50;

% Create a .mat file to save all the matrices
fileName = 'matricesnoise_50.mat';

% Loop to create and save 50 matrices
for k = 1:50
    % Generate matrix A (you can replace this with any matrix generation logic)
    A = creatingA(m);  
    C = (1/20)*creatingA(m); 

    % Create random permutation matrix P
    idx = randperm(m);
    I = eye(m);
    P = I(:, idx);
    
    % Compute matrix B
    B = P' * A * P+C;
    
    % Save each matrix as a separate variable in the .mat file
    % The variable names will be A1, A2, ..., A50, and also save P and B
    eval(sprintf('A%d = A;', k));
    eval(sprintf('P%d = P;', k));
    eval(sprintf('B%d = B;', k));
end

% Save all matrices and permutations in the file
save(fileName);
