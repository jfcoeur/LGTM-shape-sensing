function [exp_aligned, R] = tipRegistration(truth, exp)
    
    % SVD-based approach to compute best rotation R that aligns in a least squares sense
    M = exp' * truth;       
    [U,~,V] = svd(M);
    R = V * U';
    
    % Enforce a proper rotation (no reflection) by checking determinant:
    if det(R) < 0
        % If det(R) is -1, we have a reflection. Fix it by flipping one axis.
        V(:,end) = -V(:,end);
        R = V * U';
    end
    
    exp_aligned = (R * exp')';

end % function tipRegistration