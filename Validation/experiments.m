function [out, elapsedTime] = experiments(GTpath, filenames, fieldnames, expnames, insertion_case, filetype, FBGidx, Box2Needle_frame, datapath, C, GTidx, B, L, weights, plotnames)

elapsedTime = struct("soft", [], "hard", [], "layers", [], "meat", [], "pig", []); % Ajout

for i = 1:length(expnames) % experiments

    rmse = [];
    directory = fullfile(GTpath, expnames(i));

    % Get FBG curvature data
    [s_FBG, kappa_FBG] = FBGexperiments(datapath, C, filenames(i), FBGidx);

    for ref = 1:size(kappa_FBG,4)

        % Initialization
        g = 1;        
        Mr = [];
        MxyzGT = [];
        r_cell = {};
        xyzGT_cell = {};
    
        for j = GTidx{i} % trials
    
            kexp = kappa_FBG(j,:,:,ref);
            kexp = reshape(kexp,[size(kexp,2),3]);
    
            % Compute shape using model
            timerVal = tic; % Ajout
            [r, ~] = shape_model(B, s_FBG, kexp, L, insertion_case(i), weights);
            temp = elapsedTime.(fieldnames(i)); % Ajout
            temp(ref,g) = toc(timerVal); % Ajout
            elapsedTime.(fieldnames(i)) = temp; % Ajout
    
            % Ground truths
            if filetype(i) == "*.xlsx" || filetype(i) == "*.dat"
                index = i;
            else
                index = g;
            end
            
            [xyzGT, r] = FBGGroundTruths(directory, filetype(i), index, 1000*r, Box2Needle_frame(:,:,ref,3));            
            xyzGT = xyzGT - xyzGT(end, :);
            r = r - r(end, :);

            % Shorten one shape
            if length(r) > length(xyzGT)
                r = r(end-length(xyzGT)+1:end, :);
            else
                xyzGT = xyzGT(end-length(r)+1:end, :);
            end 
    
            % For CT registration
            Mr = [Mr; r];
            MxyzGT = [MxyzGT; xyzGT];
            r_cell{j} = r;
            xyzGT_cell{j} = xyzGT;       
    
            % Align with GT and compute error
            [r_aligned, ~] = tipRegistration(xyzGT, r);
            [rmse(g,ref), err_tip(g)] = errors(xyzGT, r_aligned);
    
            g = g + 1;
    
        end
    
        % CT registration
        [~, R] = tipRegistration(MxyzGT, Mr);
    
        % Compute error with CT
        if i == 4
            g = 1;
            for j = GTidx{i} % trials
                r = r_cell{j};
                xyzGT = xyzGT_cell{j};
                r_aligned = (R * r')';
                [rmse(g,ref), err_tip(g)] = errors(xyzGT, r_aligned);
                g = g + 1;
            end
        end        

    end % r = ref

    [~,idxmin] = min(mean(rmse,1));
    rmse = rmse(:,idxmin);
    out.(fieldnames(i)) = rmse;
    errorplot(categorical(plotnames(i)), rmse, "Scenario", "Experimental Results")

end

end