function r_GT = GTvalidation(k, alpha, s_vals)

% Initialization
s_vals(end+1) = 2*s_vals(end) - s_vals(end-1);
N = length(s_vals);
kappa = ones(N-1,3); % Curvature initialization

% Curvatures
kappa = [k*sind(alpha), k*cosd(alpha), 0].*kappa;

% Rotation matrix
R = eye(3); % Initial rotation matrix assuming initial rotation is identity
r_GT = zeros(3, N);
e3 = [0; 0; 1];  % Initial tangent vector

% Compute rotation matrices along the needle
for i = 2:N

    ds = s_vals(i) - s_vals(i-1);
    omega = kappa(i-1,:)';

    % Update rotation matrix using exponential map    
    R = R * expm(skew(omega * ds));
    t = R * e3;
    r_GT(:,i) = r_GT(:,i-1) + t * ds;

end

r_GT = transpose(r_GT);

end