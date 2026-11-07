
%% master_fminsearch_multistart.m
% Revised version: T_ref is fixed at 2.1 ms and beta, K_d and V_thr are optimised (plt_main_FST2.m,
% plt_main_ISI2.m). The 150 sets in param.mat were obtained with NSGA-II in Python (python/optimize.py;
% population 150, 100 generations, seed 1) instead of gamultiobj, which returns a different front
% (by default at most 35% of the population).
clear; clc;

%% --- Set up paths & data ---
modelNum    = '61-1';
mainFolder  = pwd;
modelFolder = fullfile(mainFolder, ['finger', modelNum]);
cd(modelFolder);

disp('  read dataMain...');
data        = xlsread('f61-1.xlsx','Sheet1');
actual_data = xlsread('f61-1.xlsx','Sheet3');
cd(mainFolder);

save('dataForOptimize.mat','data');
save('actual_dataForOptimize.mat','actual_data');

%% 1) 探索範囲
lb = [-2.5e-6, 0.1e-7, 30.0]; % beta, K_d, V_thr (v1.0 also optimised T_ref in [0, 1] ms)
ub = [ 2.5e-6, 2.5e-7, 35.0];
nVar = 3;

%% 2) 目的関数ハンドル
multiObj = @(x) [
        plt_main_FST2(x)          ; ...
        plt_main_ISI2(x)            % ← maximise FSS
        ]';

%% 3) GA オプション & 実行
opts = optimoptions('gamultiobj',...
        'PopulationSize' , 150 , ...
        'MaxGenerations' , 100 , ...
        'PlotFcn'        , {@gaplotpareto} , ...   % リアルタイム Pareto 図
        'Display'        , 'iter');

[popSol , popCost] = gamultiobj(multiObj , nVar , ...
                                [] , [] , [] , [] , lb , ub , opts);



paramTemp = [popSol(:,1:2), 2.1*ones(size(popSol,1),1), popSol(:,3)]; % columns of param.mat: beta, K_d, T_ref, V_thr
popCostTemp = popCost;
save('paramTemp.mat','paramTemp');
save('popCostTemp.mat','popCostTemp');



