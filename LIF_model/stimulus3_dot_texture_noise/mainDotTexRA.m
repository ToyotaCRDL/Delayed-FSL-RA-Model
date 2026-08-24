

clear
clc




%**********************************************
global R;R = 66.11e6;
global C;C = 1.01e-6;
global tau;tau = 66.77; % tau=R*C
global M;M = 1.0; % number of merkels
%**********************************************


load param.mat
load expFireRate.mat
load expTempVar.mat
dotSpace = [1.3; 2.4; 3.2; 4.3; 5.2; 6.2];


for mlp2 = 1:size(param,1)

receptorName = 'RA';
paramSet='setA';
global beta;beta  = param(mlp2,1);
global lamda;lamda = param(mlp2,2);
global v_thr;v_thr = param(mlp2,4);
global ref;ref = param(mlp2,3);
global u_thr;u_thr = 1.1;

figure(11); clf


for mlp = 1:6

figure(1); clf

clearvars -except R C tau M beta alpha lamda v_thr ref u_thr receptorName paramSet mlp mlp2 param ...
                  vec_allFiringRate_1D_SEP vec_allFiringRateHistogram vec_tempoVari vec_allFiringRateHistogramArray vec_tempoVariArray ...
                  expFireRate expTempVar dotSpace setFiringRateR2 setTempoVariR2 setFiringRateRMSE setTempoVariRMSE

%*****************************************************************************
if     mlp==1; modelNum='58-28' ; %dot spacing=6.2mm %**********
elseif mlp==2; modelNum='58-27'; %dot spacing=5.2mm %**********
elseif mlp==3; modelNum='58-26'; %dot spacing=4.3mm %**********
elseif mlp==4; modelNum='58-25' ; %dot spacing=3.2mm %**********
elseif mlp==5; modelNum='58-24'; %dot spacing=2.4mm %**********
elseif mlp==6; modelNum='58-23'; %dot spacing=1.3mm %**********
end
%*****************************************************************************


mainFolder = pwd;
modelFolder = strcat('finger',modelNum);


%*****************************************************************************
disp('  read coordinate...')
dataCoordY = xlsread('data_coord.xlsx','sheet1'); % ドットFEMの座標系は指の幅方向がX軸，指の長軸方向がY軸で左に正をとるので，XY入れ替えてX正負反転
dataCoordX = -1 * (xlsread('data_coord.xlsx','sheet2')); % sheet1=FEM X座標 / sheet2=FEM Y座標 
%*****************************************************************************

cd (modelFolder);

disp('  read nodeNum...')
load nodeNum.mat
numNodeY = size(nodeNum,1);

disp('  read dataMain...')
load dataMain.mat

disp('  read dataSED...')
load dataSED.mat

disp('  read iniYcoordCenter...')
load iniYcoordCenter.mat

disp('  read dataTime...')
load dataTime.mat

cd (mainFolder);


for i=1:size(dataTime,1)
    dataSEDWhole{i,1} = reshape(dataSED(i,:),[size(dataCoordX,1),size(dataCoordY,1)]);
end

j=1;
for i=5:-0.2:-5
    resamplePre(j,1) = i;
    j=j+1;
end
[dataCoordResampX,dataCoordResampY] = meshgrid(resamplePre,resamplePre);
dataCoordResampX = -1*dataCoordResampX;
dataCoordResampX(:,1) = dataCoordX(1,1);
dataCoordResampX(:,size(dataCoordResampX,1)) = dataCoordX(1,size(dataCoordX,1));
dataCoordResampY(1,:) = dataCoordY(1,1);
dataCoordResampY(size(dataCoordResampY,1),:) = dataCoordY(size(dataCoordY,1),1);

for i=1:size(dataTime,1)
    dataSEDWholeResamp{i,1} = griddata(dataCoordX, dataCoordY, dataSEDWhole{i,1}, dataCoordResampX, dataCoordResampY);
end

for i=1:size(dataTime,1)
    for j=1:( size(dataSEDWholeResamp{1,1},1)*size(dataSEDWholeResamp{1,1},2) )
        dataSEDPost1(i,:) = reshape(dataSEDWholeResamp{i,1},1,[]);
    end
end


dataSED = dataSEDPost1;
iniYcoordCenter = dataCoordResampY(:,1);

clear dataMain
dataMain(:,1) = dataSED(:,1);
dataMain(:,2) = dataTime(:,1)*1000;
dataMain(1:size(dataSED,1),3) = (dataTime(2,1)-dataTime(1,1))*1000;
dataMain(:,4) = dataTime(:,1);


jj=1;
kk=1;
%***start main loop
for ii=1276:1326+1

disp(['  mlp2 = ',num2str(mlp2),',  mlp = ',num2str(mlp),'  calculate ',num2str(ii),'th node out of ',num2str(size(dataSED,2)),' nodes'])
              
dataMain(:,1) = dataSED(:,ii);%***SED
dataMain(:,4) = dataTime(:,ii);%***yposi or time


[firingRate,spikeTime,spikePosition,spikeMarker,time,membPotential,time2,cur2,dataSED2] = spikeRA(dataMain);

allSpikeTime{ii,1} = spikeTime;
allFiringRate{ii,1} = firingRate;


if ( 1275<ii && ii<1327 ) 
    clear spikePositionYcoord err
    spikePositionYcoord(length(spikeTime),1)=0;
    for nn=1:length(spikeTime)
        spikePositionYcoord(nn,1) = iniYcoordCenter(ii-1275,1);
        err(nn,1)=0.075;
    end
    
    figure(1)
    er=errorbar(spikeTime/1000,spikePositionYcoord,err,'CapSize',0,'color',[50/255 50/255 50/255]);
    er.LineStyle = 'none';  
    if strcmp(modelNum,'58-27')
        xlim([1.15 1.65]) %**********
    else
        xlim([1.2 1.7]) %**********
    end
    ylim([-6 6])
    hold on

    %*********************************************************display multiple SEPs
    figure(11)
    subplot(1,6,mlp)
    er=errorbar(spikeTime/1000,spikePositionYcoord,err,'CapSize',0,'color',[50/255 50/255 50/255]);
    er.LineStyle = 'none';  
    if strcmp(modelNum,'58-27')
        xlim([1.15 1.65]) %**********
    else
        xlim([1.2 1.7]) %**********
    end
    ylim([-6 6])
    hold on
    %*********************************************************

end


end % end main loop


%**********************************************************************************output each spikeTime
% maxLength = max(cellfun(@length, allSpikeTime));
% outSpikeTimeTemp = zeros(51, maxLength);
% for i=1:51
%     for j=1:size(allSpikeTime{i+1275,1},1)
%         outSpikeTimeTemp(i,j) = allSpikeTime{i+1275,1}(j,1);
%     end
% end
% 
% cd ('result1')
% 
% if     mlp==1; ds='6-2'; %dot spacing=6.2mm %**********
% elseif mlp==2; ds='5-2'; %dot spacing=5.2mm %**********
% elseif mlp==3; ds='4-3'; %dot spacing=4.3mm %**********
% elseif mlp==4; ds='3-2'; %dot spacing=3.2mm %**********
% elseif mlp==5; ds='2-4'; %dot spacing=2.4mm %**********
% elseif mlp==6; ds='1-3'; %dot spacing=1.3mm %**********
% end
% % filename = strcat('dot_texture_sep_param_',num2str(mlp2),'_ds_',ds,'_noise_2_5.csv');
% filename = strcat('dot_texture_sep_param_',num2str(mlp2),'_ds_',ds,'_noise_5_0.csv');
% 
% if exist(filename, 'file') == 2
%     delete(filename);
% end
% 
% writematrix(outSpikeTimeTemp, filename);
% 
% cd (mainFolder);
%**********************************************************************************


%*********************************************************Temporal Variation
if strcmp(modelNum,'58-27')
    dataTime2=[1.15:0.01:1.65]'; %**********
else
    dataTime2=[1.2:0.01:1.7]'; %**********
end
for ii=1276:1326
    for i=1:size(dataTime2,1)-1
        numSpikesRA=0;
        counts=0; index=0;
        for j=1:size(allSpikeTime{ii,1},1)
            if dataTime2(i,1)<=allSpikeTime{ii,1}(j,1)/1000 && allSpikeTime{ii,1}(j,1)/1000<=dataTime2(i+1,1) && index==0
                counts = (allSpikeTime{ii,1}(j,1)/1000-dataTime2(i,1)) / (allSpikeTime{ii,1}(j,1)/1000-allSpikeTime{ii,1}(j-1,1)/1000);
                index=1;
            elseif dataTime2(i,1)<=allSpikeTime{ii,1}(j,1)/1000 && allSpikeTime{ii,1}(j,1)/1000<=dataTime2(i+1,1) && index==1
                counts = counts + 1;
            elseif dataTime2(i+1,1)<=allSpikeTime{ii,1}(j,1)/1000 && index==1
                counts = counts + (dataTime2(i+1,1)-allSpikeTime{ii,1}(j-1,1)/1000) / (allSpikeTime{ii,1}(j,1)/1000-allSpikeTime{ii,1}(j-1,1)/1000);
                index=0;
            end
        end
        allFiringRateHistogram(i,ii-1275) = counts/0.01;
    end    
end
tempoVari=0;
for i=1:size(allFiringRateHistogram,2)
% for i=26:26
    tempoVariPre=0;
    for j=1:size(allFiringRateHistogram,1)-13
        tempoVari_tamp = abs(allFiringRateHistogram(j+13,i)-allFiringRateHistogram(j,i));
        tempoVariPre = tempoVariPre + tempoVari_tamp/(size(allFiringRateHistogram,1)-13);
    end
    tempoVari = tempoVari + tempoVariPre/(size(allFiringRateHistogram,2));
end
%*********************************************************


vec_allFiringRateHistogram(7-mlp,1) = mean(allFiringRateHistogram,'all');
vec_tempoVari(7-mlp,1) = mean(tempoVari,'all');


end % end mlp


vec_allFiringRateHistogramArray(:,mlp2) = vec_allFiringRateHistogram;
vec_tempoVariArray(:,mlp2) = vec_tempoVari;


% 実測値と予測値のベクトルを定義
predicted = vec_allFiringRateHistogram;   % 予測値
actual = expFireRate;      % 実測値
% RMSEの計算
errors = actual - predicted;
squaredErrors = errors .^ 2;
meanSquaredError = mean(squaredErrors);
rmse = sqrt(meanSquaredError);
disp(['  FireRate_RMSE = ',num2str(rmse)])
setFiringRateRMSE(mlp2,1) = rmse;


% 実測値と予測値のベクトルを定義
predicted = vec_tempoVari;   % 予測値
actual = expTempVar;      % 実測値
% RMSEの計算
errors = actual - predicted;
squaredErrors = errors .^ 2;
meanSquaredError = mean(squaredErrors);
rmse = sqrt(meanSquaredError);
disp(['  TempVar_RMSE = ',num2str(rmse)])
setTempoVariRMSE(mlp2,1) = rmse;


end % end mlp2




