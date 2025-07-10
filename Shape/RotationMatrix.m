function R_s = RotationMatrix(s_vals, kappa)

kappa(:,[1 2]) = kappa(:,[2 1]);
N = length(s_vals);
R_s = zeros(3,3,N); % Rotation matrices at position s
R = eye(3); % Initial rotation matrix assuming initial rotation is identity
R_s(:,:,1) = R;

% Compute rotation matrices along the needle
for i = 2:N

    ds = s_vals(i) - s_vals(i-1);
    omega = kappa(i-1,:)';

    % Update rotation matrix using exponential map
    R = R * expm(skew(omega * ds));
    R_s(:,:,i) = R;

end

end