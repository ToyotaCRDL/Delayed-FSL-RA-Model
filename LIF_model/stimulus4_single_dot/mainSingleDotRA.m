

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

receptorName = 'RA';
paramSet='setA';
global beta;  beta  = param(mlp2,1);
global lamda; lamda = param(mlp2,2);
global v_thr; v_thr = param(mlp2,4);
global ref;   ref   = param(mlp2,3);
global u_thr; u_thr = 1.1;


for mlp = 1:3

figure(1); clf
figure(2); clf

clearvars -except R C tau M beta alpha lamda v_thr ref u_thr receptorName paramSet mlp mlp2 param setRespWidth spikeWidth

%*****************************************************************************
if     mlp==1; modelNum='70-2';
elseif mlp==2; modelNum='70-3';
elseif mlp==3; modelNum='70-4';
end
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
dataCoordResampY(1,:) = dataCoordY(1,size(dataCoordY,2));
dataCoordResampY(size(dataCoordResampY,1),:) = dataCoordY(1,1);

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

spikeTime(1,:) = []; % ************************************************************
firingRate(1,:) = [];
firingRate(1,:) = 0;
spikeMarker(1,:) = [];

allSpikeTime{ii,1} = spikeTime;
allFiringRate{ii,1} = firingRate;


if ( 1275<ii && ii<1327 ) && ~isempty(spikeTime)
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
    if i==1
        allFiringRate2{ii,1}(i,1) = 0;
    elseif i==size(spikeTime,1)
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
% cd ('result1')
% 
% if     mlp==1; dw='4-1';
% elseif mlp==2; dw='2-6';
% elseif mlp==3; dw='1-6';
% end
% filename = strcat('single_dot_sep_param_',num2str(mlp2),'_dw_',dw,'.csv');
% 
% if exist(filename, 'file') == 2
%     delete(filename);
% end
% 
% writematrix(outSpikeTimeTemp, filename);
% 
% cd (mainFolder);
%**********************************************************************************


ii = 1301; % ************************************************************
for j=1:size(allSpikeTime{ii,1},1)
    if     strcmp(modelNum,'70-2') && ( 1245.0 <= allSpikeTime{ii,1}(j,1) && allSpikeTime{ii,1}(j,1) <= 1795.0 ) % within the stimulus period
        spikeWidthTemp(j,1) = allSpikeTime{ii,1}(j,1);
    elseif strcmp(modelNum,'70-3') && ( 1295.0 <= allSpikeTime{ii,1}(j,1) && allSpikeTime{ii,1}(j,1) <= 1750.0 )
        spikeWidthTemp(j,1) = allSpikeTime{ii,1}(j,1);
    elseif strcmp(modelNum,'70-4') && ( 1317.5 <= allSpikeTime{ii,1}(j,1) && allSpikeTime{ii,1}(j,1) <= 1727.5 )
        spikeWidthTemp(j,1) = allSpikeTime{ii,1}(j,1);
    end
end

spikeWidth(mlp2,mlp) = ( spikeWidthTemp(end,1) - spikeWidthTemp(find(spikeWidthTemp~=0,1,'first'),1) ) / 1000.0 * 10.0; % /1000=ms→s，*10(m/s)=速度をかけて時間スケールを距離スケールに

for i=0:50
    firingRateBinTemp=0; count=0;
    for ii=1298:1304
        for j=1:size(allSpikeTime{ii,1},1)
            if i*20+1000 <= allSpikeTime{ii,1}(j,1) && allSpikeTime{ii,1}(j,1) < (i+1)*20+1000
                if     strcmp(modelNum,'70-2') && ( allSpikeTime{ii,1}(j,1) <= 1245.0 || 1795.0 <= allSpikeTime{ii,1}(j,1) ) % outside the stimulus period

                elseif strcmp(modelNum,'70-3') && ( allSpikeTime{ii,1}(j,1) <= 1295.0 || 1750.0 <= allSpikeTime{ii,1}(j,1) )

                elseif strcmp(modelNum,'70-4') && ( allSpikeTime{ii,1}(j,1) <= 1317.5 || 1727.5 <= allSpikeTime{ii,1}(j,1) )

                else
                    firingRateBinTemp = firingRateBinTemp + allFiringRate2{ii,1}(j,1);
                    count = count + 1;
                end
            end
        end
    end
    if count==0
        firingRateBin(i+1,1) = 0;
    else
        firingRateBin(i+1,1) = firingRateBinTemp / count;
    end            
    firingRateBinTime(i+1,1) = i*20+1000;
end

pDataSt = find(firingRateBin~=0,1,'first') - 1;
pDataEn = find(firingRateBin~=0,1,'last') + 1;
p = polyfit(firingRateBinTime(pDataSt:pDataEn,1),firingRateBin(pDataSt:pDataEn,1),8);

figure(2); clf
handleBar = bar(firingRateBinTime/1000,firingRateBin, 'Barwidth', 1);
handleBar.FaceColor = 'flat';
for i=1:51; handleBar.CData(i,:) = [164/255 217/255 234/255]; end
hold on
plot(firingRateBinTime(pDataSt:pDataEn,1)/1000,polyval(p,firingRateBinTime(pDataSt:pDataEn,1)),'Color',[38/255 38/255 38/255])

% 使用しているパラメータだと#48からmlp==1の発火が片側1つになり、応答幅が求められなくなるのでNaNを挿入
if mlp2 > 47 && mlp==1
    setRespWidth(mlp2,mlp) = NaN;
else

ylim([0, max(firingRateBin)*1.2])
xlim([1 2])

y_target = max(polyval(p,firingRateBinTime(pDataSt:pDataEn,1)),[],'all') * 0.25;
y_target2(1:51,1) = y_target;
figure(2); plot((1:0.02:2),y_target2,'--','Color',[38/255 38/255 38/255]);

p_adj = p;
p_adj(end) = p_adj(end) - y_target;
x_values = roots(p_adj);
x_values_real = x_values(imag(x_values)==0);
x_values_real = sort( x_values_real( firingRateBinTime(pDataSt,1) < x_values_real & x_values_real < firingRateBinTime(pDataEn,1) ) , 'descend' );
respWidth = (x_values_real(1) - x_values_real(end)) / 1000.0 * 10.0; % /1000=ms→s，*10(m/s)=速度をかけて時間スケールを距離スケールに

setRespWidth(mlp2,mlp) = respWidth;

end


end % end mlp


end % end mlp2




