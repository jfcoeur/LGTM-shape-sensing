clearvars
clc
close all

calib = load("D:\Distributed needles\Hopkins\Processed data\strain0_sig.mat").results;
valid = load("D:\Distributed needles\Hopkins\Processed data\strain90_sig.mat").results;
idx = [10 17 24 28];
kcalib = [0.5; 1.6; 2; 2.5; 3.2; 4];
kvalid = [0.25; 0.8; 1; 1.25; 3.125];
ang_calib = 5;
ang_valid = ang_calib + 9;
k = [kcalib, zeros(size(kcalib)) ; zeros(size(kvalid)), kvalid];
k_plot = [kcalib ; kvalid];

% Strain
[sig1c, sig2c, sig3c] = ExtractSig(calib, kcalib, ang_calib, idx);
[sig1v, sig2v, sig3v] = ExtractSig(valid, kvalid, ang_valid, idx);

% Wavelength shift
k_idx = 6;
wvl_shift = [sig1c(k_idx,:); sig2c(k_idx,:); sig3c(k_idx,:)];

% Temperature compensation
[AA1c, AA2c, AA3c, AA4c] = Tcomp(sig1c, sig2c, sig3c);
[AA1v, AA2v, AA3v, AA4v] = Tcomp(sig1v, sig2v, sig3v);

% Combine data of sensing areas
AA1 = [AA1c; AA1v];
AA2 = [AA2c; AA2v];
AA3 = [AA3c; AA3v];
AA4 = [AA4c; AA4v];

% Linear mapping
CAA1 = pinv(AA1)*k;
CAA2 = pinv(AA2)*k;
CAA3 = pinv(AA3)*k;
CAA4 = pinv(AA4)*k;

% Fitting validation
KAA1 = AA1*CAA1;
KAA2 = AA2*CAA2;
KAA3 = AA3*CAA3;
KAA4 = AA4*CAA4;

% Errors
MSE1 = mean((k - KAA1).^2, "all");
MSE2 = mean((k - KAA2).^2, "all");
MSE3 = mean((k - KAA3).^2, "all");
MSE4 = mean((k - KAA4).^2, "all");
MSE = [MSE1 MSE2 MSE3 MSE4];

% Reliability weights for each sensing area
weights = 1./MSE;
weights = weights/sum(weights); % Normalize so sum = 1

% Weighted shape
KAA1w = weights(1)*KAA1;
KAA2w = weights(2)*KAA2;
KAA3w = weights(3)*KAA3;
KAA4w = weights(4)*KAA4;
k_weight = KAA1w + KAA2w + KAA3w + KAA4w;
k_weight_plot = [k_weight(1:6, 1) ; k_weight(7:end, 2)];

% Plotting before weigthing
% i = 1;
% i = plotting(k_plot, [KAA1(1:6, 1) ; KAA1(7:end, 2)], i);
% i = plotting(k_plot, [KAA2(1:6, 1) ; KAA2(7:end, 2)], i);
% i = plotting(k_plot, [KAA3(1:6, 1) ; KAA3(7:end, 2)], i);
% i = plotting(k_plot, [KAA4(1:6, 1) ; KAA4(7:end, 2)], i);

% Error table
err1 = [mean(abs(k-KAA1), "all"), std(abs(k-KAA1), 0, "all")];
err2 = [mean(abs(k-KAA2), "all"), std(abs(k-KAA2), 0, "all")];
err3 = [mean(abs(k-KAA3), "all"), std(abs(k-KAA3), 0, "all")];
err4 = [mean(abs(k-KAA4), "all"), std(abs(k-KAA4), 0, "all")];
err = [err1; err2; err3; err4];

% Ploting after weigthing
% [k_plot, idx] = sort(k_plot);
% k_weight_plot = k_weight_plot(idx);
% plot(k_plot,k_weight_plot,'o--', 'LineWidth',2)
% title("Weighted curvature")
% xlabel("Real curvature [m^{-1}]")
% ylabel("Exp. curvature [m^{-1}")
% set(gca,"FontSize",15)
% for j = 1:length(k_plot)
%     text(k_plot(j), k_weight_plot(j), sprintf('(%.2f, %.2f)', k_plot(j), k_weight_plot(j)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
% end
errw = [mean(abs(k-k_weight), "all"), std(abs(k-k_weight), 0, "all")];


function [sig1, sig2, sig3] = ExtractSig(data, k, ang, idx)

    for i = 1:length(k)
    
        field = fieldnames(data.fiber1{1});
    
        data1 = data.fiber1{1}.(field{i});
        sig1(i,:) = data1(ang,idx);
        data2 = data.fiber2{1}.(field{i});
        sig2(i,:) = data2(ang,idx);
        data3 = data.fiber3{1}.(field{i});
        sig3(i,:) = data3(ang,idx);
    
    end

end

function [AA1, AA2, AA3, AA4] = Tcomp(sig1, sig2, sig3)

    AA1 = [sig1(:,1) sig2(:,1) sig3(:,1)];
    AA2 = [sig1(:,2) sig2(:,2) sig3(:,2)];
    AA3 = [sig1(:,3) sig2(:,3) sig3(:,3)];
    AA4 = [sig1(:,4) sig2(:,4) sig3(:,4)];
    
    AA1 = AA1 - mean(AA1,2);
    AA2 = AA2 - mean(AA2,2);
    AA3 = AA3 - mean(AA3,2);
    AA4 = AA4 - mean(AA4,2);

end

function [i] = plotting(k, KAA, i)

    [k, idx] = sort(k);
    KAA = KAA(idx);

    subplot(2,2,i)
    plot(k,KAA,'o--', 'LineWidth',2)
    
    tit = "AA" + string(i);
    title(tit)
    xlabel("Real curvature [m^{-1}]")
    ylabel("Exp. curvature [m^{-1}")
    set(gca,"FontSize",15)
    
    for j = 1:length(k)
        text(k(j), KAA(j), sprintf('(%.2f, %.2f)', k(j), KAA(j)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    end

    i = i + 1;

end