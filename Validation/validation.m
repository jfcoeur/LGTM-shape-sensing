function validation(datapath, C, idx, k, ang, B, L, weights)
    
    % Real data (distributed sensor)
    [s_FBG, kappa_FBG] = FBGvalidation(datapath, C, idx);

    n = 1;
    alpha = ang(1,:);

    t = 1;

    for i = k % curvature
    
        if i == k(7)
            alpha = ang(2,:);
        end
    
        rmse = [];
        for j = 1:length(alpha) % angle
    
            kappa = reshape(kappa_FBG(n,:,:),[size(kappa_FBG,2),3]);
            [r, s] = shape_model(B, s_FBG, kappa, L, 'single_layer_C', weights);
            r_GT = GTvalidation(i, alpha(j), s);
            [rmse(j), ~] = errors(r_GT, r);
            n = n + 1;

            temp(t,j) = 1000*rmse(j);
        end
        
        t = t + 1;

        errorplot(i, 1000*rmse, "Curvature [m^{-1}]", "Calibration Shape-Reconstruction Error")
    
    end

    [ksorted,I] = sort(k);
    ksmall = ksorted < 1.5;
    klarge = ksorted > 1.5;
    temp = temp(I,:);
    small = temp(ksmall,:);
    large = temp(klarge,:);

    column = [mean(temp,"all"), std(temp, 0, "all");
        mean(small,"all"), std(small, 0, "all");
        mean(large,"all"), std(large, 0, "all")];
end