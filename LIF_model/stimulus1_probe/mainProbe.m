

clear
clf

close all

global R;R = 66.11e6;
global C;C = 1.01e-6;
global tau;tau = 66.77; % tau=R*C
global M;M = 1.0; % number of merkels
%**********************************************
global beta;
global lamda;
global v_thr;
global ref;
global u_thr;

%*****************************************************************************
modelNum='61-1'; %dot spacing=3.2mm
%*****************************************************************************


mainFolder = pwd;
modelFolder = strcat('finger',modelNum);

cd (modelFolder);

disp('  read dataMain...')
data          = xlsread(strcat('f',modelNum,'.xlsx'),'sheet1');
dataDisp      = xlsread(strcat('f',modelNum,'.xlsx'),'sheet2');
actual_data   = xlsread(strcat('f',modelNum,'.xlsx'),'sheet3');

cd (mainFolder);


red1 = [165/255  15/255  21/255];
red2 = [222/255  45/255  38/255];
red3 = [252/255 146/255 114/255];
red4 = [252/255 187/255 160/255];
redPallete{1,1}=red1;
redPallete{2,1}=red2;
redPallete{3,1}=red3;
redPallete{4,1}=red4;



%*****************************************************************************
load param.mat
%*****************************************************************************

for mlp = 1:size(param,1)
% for mlp = 53:53

beta  = param(mlp,1);
lamda = param(mlp,2);
v_thr = param(mlp,4);
ref   = param(mlp,3);
u_thr = 1.1;

start = 0.00;
finish = 200.0;
point = 401;

figure(1); clf
disp(['  mlp = ',num2str(mlp)])

for d_count = 1:6

    x = data(:,7)*1000;
    v = data(:,d_count);
    F = griddedInterpolant(x,v,'spline');
    xq = linspace(start,finish,point);
    vq = F(xq);

    dataMain(:,1) = vq; %SED
    dataMain(:,2) = xq; %Time
    dataMain(1:point,3) = 0.50; %dt
    dataMain(:,4) = linspace(start,finish/1000,point);
    
    jj=1;
    kk=1;

    for ii=1:1
        
        disp(['  calculate ',num2str(ii),'th node out of ',num2str(size(data,7)),' nodes'])
        [firingRate,spikeTime,spikePosition,spikeMarker,time,membPotential,time2,cur2,dataSED2] = spike(dataMain);

        adata = actual_data(:,d_count);
        is_nan = isnan(adata);
        adata(is_nan) = [];
        
        adata_without_zero = adata;
        indices = adata_without_zero==0;
        adata_without_zero(indices) = [];

        maker = ones(length(adata_without_zero));

        model_firingRate = [];
        model_firingRate(1,1) = 0.0;
        for firetime = 2:length(adata)
            model_firingRate(firetime, 1) = 1/(adata(firetime, 1) * 0.001 - adata(firetime - 1, 1) * 0.001); 
        end
        
        
        figure(1)
        subplot(2,6,d_count ) %%%%%
        stem(adata_without_zero,maker*1.2,'LineWidth',0.8,'Marker','none', 'Color', [50/255 50/255 50/255])
        hold on
        stem(adata_without_zero,maker*-1.2,'LineWidth',0.8,'Marker','none','Color', [50/255 50/255 50/255])
        xlim([0, 100])
        ylim([-2 2])
        yticks([])
        ax = gca;
        ax.TickLength = [0 0];
        ax.FontSize = 6;
        ax.FontName = 'Arial';

        subplot(2,6,d_count + 6) %%%%%
        stem(spikeTime,spikeMarker*1.2,'LineWidth',0.8,'Marker','none', 'Color', red4) % *********************************************************
        hold on
        stem(spikeTime,spikeMarker*-1.2,'LineWidth',0.8,'Marker','none','Color', red4)
        xlim([0, 100])
        ylim([-2 2])
        yticks([])
        ax = gca;
        ax.TickLength = [0 0];
        ax.FontSize = 6;
        ax.FontName = 'Arial';
        fig(2000) = gcf;
        fig(2000).PaperPosition = [0 0 20.6 2.5];
        saveas(gcf,'spikeTrain1.svg')
        

        
        setDiffFST(d_count,mlp) = spikeTime(2,1) - adata(2,1);

        setDiffMeanFirRate(d_count,mlp) = mean(firingRate(3:end,1),"all") - mean(model_firingRate(3:end,1),"all");


        indices = spikeTime==0;
        spikeTime(indices) = [];
        
        indices = adata==0;
        adata(indices) = [];


        %***************************************************************************output each spikeTime
        % cd ('result2')
        % 
        % if     d_count==1; ps='9-5';
        % elseif d_count==2; ps='14-8';
        % elseif d_count==3; ps='20-6';
        % elseif d_count==4; ps='34-0';
        % elseif d_count==5; ps='46-5';
        % elseif d_count==6; ps='60-6';
        % end
        % filename = strcat('probe_spiketrain_param_',num2str(mlp),'_ps_',ps,'.csv');
        % 
        % if exist(filename, 'file') == 2
        %     delete(filename);
        % end
        % 
        % writematrix(spikeTime, filename);
        % 
        % cd (mainFolder);
        %***************************************************************************


    end

end

end


%***************************************************************************
% save('setDiffFST.mat','setDiffFST');
% save('setDiffMeanFirRate.mat','setDiffMeanFirRate');
%***************************************************************************




