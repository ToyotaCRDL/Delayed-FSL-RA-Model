
%% master_fminsearch_multistart.m
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
lb = [-2.5e-6, 0.1e-7, 0.0, 30.0];
ub = [ 2.5e-6, 2.5e-7, 1.0, 35.0];
nVar = 4;

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



paramTemp = popSol;
popCostTemp = popCost;
save('paramTemp.mat','paramTemp');
save('popCostTemp.mat','popCostTemp');



