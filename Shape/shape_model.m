function [r, s] = shape_model(B, s_FBG, kappa_FBG, L, insertion_case, weights)

    % Optimization parameters
    [params0, lb, ub] = InsertionCase(insertion_case, L);
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter');
    
    % Cost function handle
    CostFunctionHandle = @(optimvars) CostFunction(optimvars, B, s_FBG, kappa_FBG, L, insertion_case, weights);
    
    % Run optimization
    [params_opt, ~] = fmincon(CostFunctionHandle, params0, [], [], [], [], lb, ub, [], options);
    
    % Optimized parameters % USELESS FOR NOW
    % [omega_init, kappa_c, s_star] = OptimParam(insertion_case, params_opt);
    
    % Compute needle shape and orientation along the shaft
    % Solve Euler-Poincaré equations with optimized parameters
    [kappa_model_s, R_s] = EulerPoincareSolver(params_opt, B, L, insertion_case);

    % 3D coordinates
    s = kappa_model_s.s;
    r = coordinates(R_s, kappa_model_s.s);

end