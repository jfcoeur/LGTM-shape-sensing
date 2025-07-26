%% Clear workspace
% clearvars; clc; close all;

% Add all MATLAB functions to path
needle = "Hopkins 3in1";
[datapath, GTpath] = directories(needle);


%% Physical constants
% Needle diameter [m]
d = 1.27e-3; 

% Stiffness matrix
B = stiffmatrix(d);

% Insertion length [m]
L = 0.12;

% Number of sensing areas and their index on the distributed needle
% N = 4; % m
% FBGidx = [10 17 24 28]; % AA positions [10 17 24 28] - 5 = [5 12 19 23]
N = 30; % g
FBGidx = 1:N; % AA positions

% Calibration/Validation curvatures and orientations
kcalib = [0.5; 1.6; 2; 2.5; 3.2; 4];
kvalid = [0.25; 0.8; 1; 1.25; 3.125];
angcalib = 0:20:340;
angvalid = 10:20:350;
k = [kcalib; kvalid]';
ang = [angcalib; angvalid];


%% Calibration

% Simulations
% [C, weights] = Csimuldata(N);

% Real data (distributed sensor)
[C, weights] = Crealdata(datapath, FBGidx, kcalib, kvalid, angcalib, angvalid);


%% Validation

% Compute and plot shape errors
% validation(datapath, C, FBGidx, k, ang, B, L, weights)


%% Experiments

% Simulations
% [s_FBG, kappa_FBG] = FBGsimuldata(C,N,L);

% Real data (distributed sensor)
if needle == "Hopkins"

    filenames = ["Soft_sig.mat", "Hard_sig.mat", "2-layer_sig.mat", "Meat_sig.mat"];
    fieldnames = ["soft", "hard", "layers", "meat"];
    expnames = ["Soft gel", "Hard gel", "2-layer gel", "Meat"];
    plotnames = ["Soft gel", "Hard gel", "2-layer gel", "Ex vivo tissue"];
    insertion_case = ["single_layer_C", "single_layer_C", "double_layer_C", "single_layer_C"];
    filetype = ["*.csv", "*.csv", "*.csv", "*.xlsx"];
    soft_idx = 1:12; hard_idx = [1,3:14]; layer_idx = [1,5:8,10:18]; meat_idx = 2:11;
    GTidx = {soft_idx, hard_idx, layer_idx, meat_idx};
    Box2Needle_frame = load(GTpath + "Meat" + filesep + "M(4,4,ref,res).mat").M;
    
    % Compute and plot shape errors
    [rmse, elapsedTime] = experiments(GTpath, filenames, fieldnames, expnames, insertion_case, filetype, FBGidx, Box2Needle_frame, datapath, C, GTidx, B, L, weights, plotnames);

end


%% Pig insertions
if needle == "Hopkins 3in1"

    % Real data (distributed sensor)
    filenames = "Pig_apnea_sig.mat";
    fieldnames = "pig";
    expnames = "Pig insertions";
    plotnames = "Pig model";
    insertion_case = "single_layer_C";
    filetype = "*.dat";
    pig_idx = [1:7, 9:10];
    GTidx = {pig_idx};
    Box2Needle_frame = eye(4);
    Box2Needle_frame = repmat(Box2Needle_frame, [1, 1, 4, 4]);
    
    % Compute and plot shape errors
    [rmse.(fieldnames), elapsedTime.(fieldnames)] = experiments(GTpath, filenames, fieldnames, expnames, insertion_case, filetype, FBGidx, Box2Needle_frame, datapath, C, GTidx, B, L, weights, plotnames);

end


%% Funtions    
function [datapath, GTpath] = directories(needle)
    if ispc
        directory = "C:\Users\jfcoe\OneDrive - Johns Hopkins\Documents\MATLAB\LGTM shape sensing";
        datapath = "D:\Distributed needles\" + needle + "\Processed data\Discrete sensing\";
        GTpath = "D:\Distributed needles\" + needle + "\Ground truths\";
    else
        directory = "/home/jfcoeur/OneDrive/Documents/MATLAB/LGTM shape sensing";
        datapath = "/media/jfcoeur/T7 Shield/Distributed needles/" + needle + "/Processed data/Discrete sensing/";
        GTpath = "/media/jfcoeur/T7 Shield/Distributed needles/" + needle + "/Ground truths/";
    end
    addpath(genpath(directory));
end