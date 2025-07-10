function [s_FBG, kappa_FBG] = FBGexperiments(directory, C, filename, idx)

% Load data
data = load(directory + filename).results;

% Strain data
sig1 = data.fiber1{1};
sig2 = data.fiber2{1};
sig3 = data.fiber3{1};

if length(idx) == 4
    % FBG positions
    % s_FBG = load(directory + "s_FBG.mat").s_FBG;
    s_FBG = [0.0199, 0.0548, 0.0897, 0.1096];
    
    % Strain data
    sig1 = sig1(:, idx, :);
    sig2 = sig2(:, idx, :);
    sig3 = sig3(:, idx, :);
else
    % FBG positions
    s_FBG = data.x{1};
    s_FBG = s_FBG(end-24:end); %%%
    s_FBG = s_FBG - s_FBG(1);
    
    % Strain data
    sig1 = sig1(:, end-24:end, :); %%%
    sig2 = sig2(:, end-24:end, :); %%%
    sig3 = sig3(:, end-24:end, :); %%%
end

sig = cat(4, sig1, sig3, sig2); %%%%
sig = permute(sig,[1 2 4 3]);

% Temperature compensation
sig = Tcomp(sig);

% FBG experimental curvature
kappa_FBG = [];
for j = 1:size(sig,4)

    sig_temp = sig(:,:,:,j);

    for i = 1:size(sig_temp,2)
    
        AA = reshape(sig_temp(:,i,:),[size(sig_temp,1),3]);
        kappa_FBG_temp(:,:,i) = AA*C(:,:,i);
    
    end

    temp = kappa_FBG_temp;    
    temp(:,3,:) = zeros(size(temp(:,1,:)));
    temp = permute(temp, [1 3 2]);
    kappa_FBG = cat(4, kappa_FBG, temp);

end

end % function FBGexperiments