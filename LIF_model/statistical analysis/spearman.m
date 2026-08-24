

clear
clc


data = readmatrix('dot_texture_firerate_rmse2.csv');
% data = readmatrix('dot_texture_tempvar_rmse2.csv');
% data = readmatrix('single_dot_respwidth2_dw_4_1.csv');
% data = readmatrix('single_dot_respwidth2_dw_2_6.csv');
% data = readmatrix('single_dot_respwidth2_dw_1_6.csv');

[rho, pValue] = corr(data(:,1), data(:,2), 'Type', 'Spearman');

fprintf('Spearman rho = %.3f\n', rho);
fprintf('p = %.4f\n', pValue);


