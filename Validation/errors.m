function [rmse, err_tip] = errors(truth, exp)

    subs = truth - exp;    
    err2 = sum(subs.^2, 2);
    rmse = sqrt(sum(err2)/length(err2));

    err = sqrt(err2);
    err_tip = err(1);

end