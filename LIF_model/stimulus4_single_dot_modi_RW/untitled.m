


% x = [1; 2; 3; 4; 5];
% y = [2; 1; 4; 3; 5];

[rho, pValue] = corr(data(:,1), data(:,2), 'Type', 'Spearman');

fprintf('Spearman rho = %.3f\n', rho);
fprintf('p = %.4f\n', pValue);