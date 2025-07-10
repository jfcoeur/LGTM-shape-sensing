function [C, weights] = Crealdata(directory, idx, kcalib, kvalid, angcalib, angvalid)

% Strain measurements
calib = load(directory + "Calibration A_sig.mat").results;
valid = load(directory + "Validation A_sig.mat").results;

% Exclude validation curvatures
exclude = [0.8, 1.25, 3.125];
kvalid = kvalid(~ismember(kvalid, exclude));
valid.fiber1{1} = rmfield(valid.fiber1{1}, {'x0_8', 'x1_25', 'x3_125'});
valid.fiber2{1} = rmfield(valid.fiber2{1}, {'x0_8', 'x1_25', 'x3_125'});
valid.fiber3{1} = rmfield(valid.fiber3{1}, {'x0_8', 'x1_25', 'x3_125'});

% Curvatures
kx = [kcalib*cosd(angcalib); kvalid*cosd(angvalid)]';
ky = [kcalib*sind(angcalib); kvalid*sind(angvalid)]';
% kx = (kcalib*cosd(angcalib))';
% ky = (kcalib*sind(angcalib))';
k = [ky(:), kx(:)];

% Strain
[sig1c, sig2c, sig3c] = ExtractSig(calib, kcalib, idx);
[sig1v, sig2v, sig3v] = ExtractSig(valid, kvalid, idx);
sig1 = [sig1c; sig1v];
sig2 = [sig2c; sig2v];
sig3 = [sig3c; sig3v];

C = [];
err = [];
MSE = [];
for i = 1:size(sig1,2)
    
    % Temperature compensation
    AAi = [sig1(:,i) sig2(:,i) sig3(:,i)];
    AAi = AAi - mean(AAi,2);
    AA(:,:,i) = AAi;

    % Linear mapping
    CAAi = pinv(AAi)*k;
    C = cat(3, C, CAAi);

    % Fitting validation
    KAAi = AAi*CAAi;
    KAA(:,:,i) = KAAi;

    % Error
    erri = [mean(abs(k-KAAi), "all"), std(abs(k-KAAi), 0, "all")];
    err = [err; erri];
 
    % MSE
    MSEi = mean((k - KAAi).^2, "all");
    MSE = [MSE MSEi];    

end

% Reliability weights for each sensing area
weights = 1./MSE;
weights = weights/sum(weights); % Normalize so sum = 1

% Weighted shape
k_weight = 0;
for j = 1:length(weights)
    KAAiw = weights(j)*KAA(:,:,j);
    k_weight = k_weight + KAAiw;
end

% Weighted error
errw = [mean(abs(k-k_weight), "all"), std(abs(k-k_weight), 0, "all")];


%%%%%%%%
function [sig1, sig3, sig2] = ExtractSig(data, k, idx) %%%%

    sig1 = [];
    sig2 = [];
    sig3 = [];

    for i = 1:length(k)
    
        field = fieldnames(data.fiber1{1});
    
        data1 = data.fiber1{1}.(field{i});
        sig1 = [sig1; data1(:,idx)];
        data2 = data.fiber2{1}.(field{i});
        sig2 = [sig2; data2(:,idx)];
        data3 = data.fiber3{1}.(field{i});
        sig3 = [sig3; data3(:,idx)];
    
    end

    if length(idx) == 30

        sig1 = sig1(:, end-24:end, :); %%%
        sig2 = sig2(:, end-24:end, :); %%%
        sig3 = sig3(:, end-24:end, :); %%%

    end
    
end
%%%%%%%%

end % function Crealdata