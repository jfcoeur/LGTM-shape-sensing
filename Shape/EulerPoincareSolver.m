function [kappa_model_s, R_s] = EulerPoincareSolver(optimvars, B, L, insertion_case)

% Euler-Poincaré Equations Solver

% Initilization
N = 25; % Number of points
s_vals = linspace(0, L, N)'; % Longitudinal positions
kappa = zeros(N,3); % Curvature values
R_s = zeros(3,3,N+1); % Rotation matrices at position s
R = eye(3); % Initial rotation matrix assuming initial rotation is identity
R_s(:,:,1) = R;

% Extract parameters
kappa0 = optimvars(1:3); % Inital curvature
kappa_params = optimvars(4:end); % Intrinsic curvatures

% Set up ODE solver options
options = odeset('RelTol',1e-6,'AbsTol',1e-8);

% Solve the Euler-Poincaré equations using an ODE solver
[~, kappa] = ode45(@(s, k) EulerPoincareODE(s, k, kappa_params, B, L, insertion_case), s_vals, kappa0, options);
s_vals(end+1) = 2*s_vals(end) - s_vals(end-1);

% Compute rotation matrices along the needle
for i = 2:length(s_vals)

    ds = s_vals(i) - s_vals(i-1);
    omega = kappa(i-1,:)';

    % Update rotation matrix using exponential map
    R = R * expm(skew(omega * ds));
    R_s(:,:,i) = R;

end

% Output curvature and rotation matrices
kappa_model_s.s = s_vals(1:end-1);
kappa_model_s.kappa = kappa;

end % function EulerPoincareSolver
