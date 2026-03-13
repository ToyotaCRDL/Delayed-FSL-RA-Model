

clear
clc


%**********************************************
global R;   R   = 66.11e6;
global C;   C   = 1.01e-6;
global tau; tau = 66.77; % tau=R*C
global M;   M   = 1.0; % number of merkels
%**********************************************


load param.mat

for mlp2 = 1:size(param,1)
% for mlp2 = 53:53

receptorName = 'RA';
paramSet='setA';
global beta;  beta  = param(mlp2,1);
global lamda; lamda = param(mlp2,2);
global v_thr; v_thr = param(mlp2,4);
global ref;   ref   = param(mlp2,3);
global u_thr; u_thr = 1.1;



figure(1); clf

clearvars -except R C tau M beta alpha lamda v_thr ref u_thr receptorName paramSet mlp mlp2 param setRespWidth

%*****************************************************************************
modelNum='57-7-3';
%*****************************************************************************

mainFolder = pwd;
modelFolder = strcat('finger',modelNum);

%*****************************************************************************
disp('  read coordinate...')
if strcmp(modelNum,'57-7-3')
    dataCoordX = xlsread('data_coord.xlsx','sheet1'); % 点字FEMの座標系
    dataCoordY = xlsread('data_coord.xlsx','sheet2');
elseif  strcmp(modelNum,'70-2') || strcmp(modelNum,'70-3') || strcmp(modelNum,'70-4')
    dataCoordY = xlsread('data_coord.xlsx','sheet1'); % ドットFEMの座標系は指の幅方向がX軸，指の長軸方向がY軸で左に正をとるので，XY入れ替えてX正負反転
    dataCoordX = -1 * (xlsread('data_coord.xlsx','sheet2')); % sheet1=FEM X座標 / sheet2=FEM Y座標
end
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

disp(['  mlp2 = ',num2str(mlp2),',  calculate ',num2str(ii),'th node out of ',num2str(size(dataSED,2)),' nodes'])

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
    if strcmp(modelNum,'57-7-3')
        xlim([1.075 1.3])
    elseif  strcmp(modelNum,'70-2') || strcmp(modelNum,'70-3') || strcmp(modelNum,'70-4')
        xlim([1.0 2.0])
    end
    ylim([-6 6])
    hold on
end


for i=1:size(spikeTime,1)
    if i==size(spikeTime,1)
        allFiringRate2{ii,1}(i,1) = firingRate(i,1);
    else
        allFiringRate2{ii,1}(i,1) = ( firingRate(i,1) + firingRate(i+1,1) )*0.5;
    end
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
% cd ('result2')
% 
% filename = strcat('braille_sep_param_',num2str(mlp2),'.csv');
% 
% if exist(filename, 'file') == 2
%     delete(filename);
% end
% 
% writematrix(outSpikeTimeTemp, filename);
% 
% cd (mainFolder);
%**********************************************************************************


end % end mlp2




