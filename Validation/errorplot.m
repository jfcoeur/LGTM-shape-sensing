function [] = errorplot(exp, rmse, x_label, tit)

rmse_ave = mean(rmse);
rmse_std = std(rmse);

errorbar(exp, rmse_ave, rmse_std, '_', 'LineWidth', 3, 'MarkerSize', 35, 'Color', 'g')
hold on

xlabel(x_label, 'FontSize', 40)
ylabel('RMSE [mm]', 'FontSize', 40)
title(tit, 'FontSize', 45)
set(gca, 'FontSize', 40)

end