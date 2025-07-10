function [s_FBG, kappa_FBG] = FBGvalidation(directory, C, idx)

% Load data
calib = load(directory + "Calibration A_sig.mat").results;
valid = load(directory + "Validation A_sig.mat").results;

% Strain data
[sig1c, sig2c, sig3c] = ExtractSig(calib, 6, idx);
[sig1v, sig2v, sig3v] = ExtractSig(valid, 5, idx);
sig1 = [sig1c; sig1v];
sig2 = [sig2c; sig2v];
sig3 = [sig3c; sig3v];
sig = cat(3, sig1, sig2, sig3);

% Temperature compensation
sig = Tcomp(sig);

% FBG positions
if length(idx) == 4
    s_FBG = [0.0199, 0.0548, 0.0897, 0.1096];
else
    s_FBG = calib.x{1};
    s_FBG = s_FBG(end-24:end); %%%
    s_FBG = s_FBG - s_FBG(1);
end

% FBG experimental curvature
for i = 1:size(sig,2)

    AA = reshape(sig(:,i,:),[size(sig,1),3]);
    kappa_FBG(:,:,i) = AA*C(:,:,i);

end

kappa_FBG(:,3,:) = zeros(size(kappa_FBG(:,1,:)));
kappa_FBG = permute(kappa_FBG, [1 3 2]);


%%%%%%%%
function [sig1, sig3, sig2] = ExtractSig(data, N, idx) %%%

    sig1 = [];
    sig2 = [];
    sig3 = [];

    for j = 1:N
    
        field = fieldnames(data.fiber1{1});
    
        data1 = data.fiber1{1}.(field{j});
        sig1 = [sig1; data1(:,idx)];
        data2 = data.fiber2{1}.(field{j});
        sig2 = [sig2; data2(:,idx)];
        data3 = data.fiber3{1}.(field{j});
        sig3 = [sig3; data3(:,idx)];
    
    end

    if length(idx) == 30

        sig1 = sig1(:, end-24:end, :); %%%
        sig2 = sig2(:, end-24:end, :); %%%
        sig3 = sig3(:, end-24:end, :); %%%

    end

end
%%%%%%%%

end % function FBGvalidation