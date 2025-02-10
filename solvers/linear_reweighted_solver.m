
% Linear reweighted regularization slover for solving graph matching problem

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



function [X,f_best,out]= linear_reweighted_solver(P, A, B, epsilon_init, lambda_init,tol1)
    max_iter1 = 1000; 
    epsilon = epsilon_init;
    lambda = lambda_init;
    c_hat = 1; 
    n = size(B,1);
    X_n = 1/n*ones(n,n);
    X_n_i = X_n;
    lambdamax=1e6;

    time_tracking = [];
    objerror_tracking = [];
    solerror_tracking = [];
    tic; 

for k = 1:max_iter1
        h= 1-sum(max(X_n))/n;       
        if 1 <= tol1
            break;
        end
        
        % X_n: Subproblem initial solution
        % X_n_i: Subproblem iterative solution
        [X, f_val] = find_min(A, B, X_n, X_n_i, lambda, epsilon);
        
        f_best = f_val;
        lambda = min(lambda + c_hat,lambdamax);
        epsilon = max(0.9 * epsilon, 1e-3); 

        X_n = X;
        X_n_i = X;

        elapsed_time = toc; 
        objective_error = abs(norm(A * X - X * B, 'fro')^2-norm(A * P - P * B, 'fro')^2);
        sol_error = norm(P - X, 'fro');
        time_tracking = [time_tracking, elapsed_time];
        objerror_tracking = [objerror_tracking, objective_error];
        solerror_tracking = [solerror_tracking, sol_error];
end

    out.time_tracking = time_tracking;
    out.objerror_tracking = objerror_tracking;
    out.solerror_tracking = solerror_tracking;
   
end

function [X_new, f_val] = find_min(A, B, X_n, X_n_i, lambda, epsilon)
    max_iter = 2000;
    alpha_init = 1e-4;
    tol_x = 1e-2;
    tol_f = 1e-2;
    theta = 0.5;
    m = size(A, 1);
    eta=0.85;
    alpha = alpha_init;
    Ci = objective_F(A, B, X_n, X_n_i, lambda, epsilon);
    Qi = 1;
    sigam=1e-4;

    for i = 1:max_iter
  
        grad_F = gradient_F(A, B, X_n, X_n_i, lambda, epsilon);

        D = dualBB_projection(X_n_i - alpha * grad_F)-X_n_i;
        
        for j = 0:100
            if objective_F(A, B, X_n, X_n_i + theta^j * D, lambda, epsilon) <= ...
            Ci+ sigam * theta^j * grad_F(:)' * D(:)
                break;
            end
        end
       
        F_val = objective_F(A, B, X_n,X_n_i, lambda, epsilon); 
        Ci = (eta * Qi * Ci + F_val) / (eta * Qi + 1);
        Qi = eta * Qi + 1;
         
        X0=X_n_i;
        X_new = X_n_i + theta^j * D;
               
        if (norm(X_new - X0, 'fro') / sqrt(m) < tol_x && ...
                abs(objective_F(A, B, X_n, X_new, lambda, epsilon) - objective_F(A, B, X_n, X0, lambda, epsilon)) ...
                / (1 + abs(objective_F(A, B,X_n, X0, lambda, epsilon))) < tol_f)
            break;
        end

        X_n_i= X_n_i + theta^j * D;
        
       
    end
    
    f_val = objective_f(A, B, X_new);
end




function grad_F = gradient_F(A, B, X0, X,lambda, epsilon)
    grad_F = A' * (A * X - X * B) - (A * X - X * B) * B' + lambda * 1 ./ (X0 + epsilon);
end

function F_val = objective_F(A, B, X0, X, lambda, epsilon)
    F_val = 0.5 * norm(A * X - X * B, 'fro')^2 + lambda * reweighted(X0,X,epsilon);
end


function f_val = objective_f(A, B, X)
    f_val = 0.5 * norm(A * X - X * B, 'fro')^2;
end


function reweighted_value = reweighted(X0,X,epsilon)
     reweighted_value = sum(sum(1./(X0+epsilon).*X));
end

function D = dualBB_projection(C)
tol=1e-1;
maxIter=5000;
m = size(C, 1);
e = ones(m,1);
[y, z] = dualBB(C, maxIter,tol);
D = projection_nonnegative(C + y * e' + e * z');
end


function [y_opt, z_opt] = dualBB(C, max_iter, tol)
  
    m = size(C, 1);
    y = zeros(m, 1);
    z = zeros(m, 1);
    e = ones(m, 1);
    
    % Initial learning rate (not used directly but as a fallback)
    alpha = 0.01;
    
    % Previous gradients and iterates for BB step size calculation
    grad_y_prev = zeros(m, 1);
    grad_z_prev = zeros(m, 1);
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

end



