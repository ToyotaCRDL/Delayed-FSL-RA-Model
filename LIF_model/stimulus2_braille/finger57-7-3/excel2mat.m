


clear
clc


mainFolder = pwd;

modelNum = extractAfter(mainFolder, 'finger');

disp('  read nodeNum...')
nodeNum = xlsread('node_num.xlsx','Sheet1');

disp('  read dataMain...')
dataMain = xlsread(strcat('f',modelNum,'.xlsx'),'sheet1');

%***SED/sed_data_proc.m/dataSED
disp('  read dataSED...')
dataSED = xlsread(strcat('f',modelNum,'.xlsx'),'sheet2');

disp('  read iniYcoordCenter...')
iniYcoordCenter = xlsread(strcat('f',modelNum,'.xlsx'),'sheet3');

disp('  read dataTime...')
dataTime = xlsread(strcat('f',modelNum,'.xlsx'),'sheet4');


save ('nodeNum.mat', 'nodeNum')
save ('dataMain.mat', 'dataMain')
save ('dataSED.mat', 'dataSED')
save ('iniYcoordCenter.mat', 'iniYcoordCenter')
save ('dataTime.mat', 'dataTime')



