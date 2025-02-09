
% LP_norm regularization slover for solving graph matching problem

%===============================================================

%Denote:
%f=1/2*||A*X-X*B||_F^2
%F=1/2*||A*X-X*B||_F^2+lp_norm(X + epsilon, p)^p;
%grad_F: Gradient of Function F

%Input:
%P: Random Generated n*n Permutation Matrix
%A: Random Generated n*n Matrix
%B: Random Generated n*n Matrix
%epsilon_init: initial epsilon, epsilon_0 in paper
%lambda_init: initial lambda, lambda_0 in paper
%tol: Stopping condition (Sparsity test), lp_norm(X, p)^p / n - 1 tolorence

%Output:
%X_best: Optimal X obtained
%f_best: Optimal function value obtained
%out.time_tracking: time (t) tracking vector
%out.objerror_tracking: objective_error_tracking vector at time (t) 
%out.solerror_tracking: solution_error_tracking vector at time (t) 

%===============================================================

function [X_best, f_best, out] = lp_norm_solver(P, A, B, p, epsilon_init, lambda_init, tol)

    max_iter = 1000; % Maximum number of iterations
    epsilon = epsilon_init; % initial epsilon
    lambda = lambda_init; % initial lambda
    c_hat = 1; % Regularization parameter increment
    n = size(B,1);
    X = 1/n*ones(n,n); % initial X
    lambdamax=1e6; % safeguards 

    X_best = X; % The best X obtained
    f_best = inf; % f(X_best)

    % Time tracking
    time_tracking = [];
    objerror_tracking = [];
    solerror_tracking = [];
    tic; % Start timer

    % Iterate until a sparse solution obtained
    for k = 1:max_iter
        h = lp_norm(X, p)^p / n - 1;
        if h <= tol
            break;
        end
        
        % Step 3: Set X_k_0
        X_k_0 = X;
        
        % Step 4:Find KKT point
        [X, f_val] = find_KKT_point(A, B, X_k_0, lambda, epsilon, p);
        
        % Update the best solution found
        if f_val < f_best
            f_best = f_val;
            X_best = X;
            lambda = min(lambda + c_hat, lambdamax);
        else        
            % Update parametersm lambda and epsilon
            lambda = min(lambda + c_hat, lambdamax);
            epsilon = max(0.9 * epsilon, 1e-3);
        end

        % Store the error and time at each iteration
        elapsed_time = toc;  % Get elapsed time
        time_tracking = [time_tracking, elapsed_time];
        objective_error = abs(norm(A * X - X * B, 'fro')^2-norm(A * P - P * B, 'fro')^2);
        sol_error = norm(P - X, 'fro');
        objerror_tracking = [objerror_tracking, objective_error];
        solerror_tracking = [solerror_tracking, sol_error];
        
    end

    % Output tracking data
    out.time_tracking = time_tracking;
    out.objerror_tracking = objerror_tracking;
    out.solerror_tracking = solerror_tracking;
end

%%
%compute an approximate KKT point

function [X, f_val] = find_KKT_point(A, B, X_init, lambda, epsilon, p)
    % Parameters for the projected gradient method
    max_iter = 500;
    alpha_init = 1e-4;
    tol_x = 1e-2;
    tol_f = 1e-2;
    theta = 0.5;
    n = size(A, 1);
    X = X_init;
    eta = 0.85;
    alpha = alpha_init;
    Ci = objective_F(A, B, X, lambda, epsilon, p);
    Qi = 1;
    sigam = 1e-4;

    for i = 1:max_iter
        % Compute gradient
        grad_F = gradient_F(A, B, X, lambda, epsilon, p);
        
        % Projection onto the set D_n
        D = dualBB_projection(X - alpha * grad_F) - X;
        
        % Line search
        for j = 0:100
            if objective_F(A, B, X + theta^j * D, lambda, epsilon, p) <= ...
                    Ci + sigam * theta^j * grad_F(:)' * D(:)
                break;
            end
        end

        F_val = objective_F(A, B, X, lambda, epsilon, p); 
        Ci = (eta * Qi * Ci + F_val) / (eta * Qi + 1);
        Qi = eta * Qi + 1;
        
        % Update X
        X0 = X;
        X_new = X + theta^j * D;
        X = X_new;
        
        % Check convergence
        if (norm(X - X0, 'fro') / sqrt(n) < tol_x && ...
                abs(objective_F(A, B, X, lambda, epsilon, p) - objective_F(A, B, X0, lambda, epsilon, p)) ...
                / (1 + abs(objective_F(A, B, X0, lambda, epsilon, p))) < tol_f)
            break;
        end
    end

    f_val = objective_f(A, B, X);
end

%%
%f=1/2*||A*X-X*B||^2
%F=1/2*||A*X-X*B||^2+lp_norm(X + epsilon, p)^p;
%grad_F: Gradient of Function F

function grad_F = gradient_F(A, B, X, lambda, epsilon, p)
    grad_F = A' * (A * X - X * B) - (A * X - X * B) * B' + lambda * p * (X + epsilon).^(p-1);
end

function F_val = objective_F(A, B, X, lambda, epsilon, p)
    F_val = 0.5 * norm(A * X - X * B, 'fro')^2 + lambda * lp_norm(X + epsilon, p)^p;
end

function f_val = objective_f(A, B, X)
    f_val = 0.5 * norm(A * X - X * B, 'fro')^2;
end

function lp_norm_value = lp_norm(X, p)
    lp_norm_value = sum(abs(X(:)).^p)^(1/p);
end

%%
%projected gradient method with a nonmonotone line search and Barzilai-Borwein step sizes

function D = dualBB_projection(C)
    tol = 1e-2;
    maxIter = 5000;
    n = size(C, 1);
    e = ones(n, 1);
    [y, z] = dualBB(C, maxIter, tol);
    D = projection_nonnegative(C + y * e' + e * z');
end

function [y_opt, z_opt] = dualBB(C, max_iter, tol)
    n = size(C, 1);
    y = zeros(n, 1);
    z = zeros(n, 1);
    e = ones(n, 1);
    
    % Initial learning rate (not used directly but as a fallback)
    alpha = 0.01;
    
    % Previous gradients and iterates for BB step size calculation
    grad_y_prev = zeros(n, 1);
    grad_z_prev = zeros(n, 1);
    y_prev = y;
    z_prev = z;
    
    for iter = 1:max_iter
        % Calculate the matrix and its projection
        M = C + y * e' + e * z';
        P_plus = max(M, 0);
        
        % Compute the gradients
        grad_y = P_plus * e - e;
        grad_z = P_plus' * e - e;
        
        % Compute the BB step size for y and z updates
        if iter > 1
            s_y = y - y_prev;
            s_z = z - z_prev;
            y_prev = y;
            z_prev = z;
            
            y_diff = grad_y - grad_y_prev;
            z_diff = grad_z - grad_z_prev;
            
            alpha_y = (s_y' * s_y) / (s_y' * y_diff);
            alpha_z = (s_z' * s_z) / (s_z' * z_diff);
            
            % Use the average of alpha_y and alpha_z as the step size
            alpha = 0.5 * (alpha_y + alpha_z);
        end
        
        % Update y and z
        y = y - alpha * grad_y;
        z = z - alpha * grad_z;
        
        % Store current gradients for the next iteration
        grad_y_prev = grad_y;
        grad_z_prev = grad_z;
        
        % Stopping criterion based on the magnitude of the gradient
        if norm([grad_y; grad_z]) < tol
            break;
        end
    end
    
    y_opt = y;
    z_opt = z;
end

function P_plus = projection_nonnegative(M)
    P_plus = max(M, 0);
end



