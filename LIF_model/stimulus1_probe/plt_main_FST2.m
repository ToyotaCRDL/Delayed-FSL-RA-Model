function sum_diff_FST = plt_main_FST(send_data)
global R;R = 66.11e6;
global C;C = 1.01e-6;
global tau;tau = 66.77; % tau=R*C
global M;M = 1.0; % number of merkels
global u_thr;u_thr = 1.1;
%**********************************************
global beta;
global lamda;
global ref;
global v_thr;


beta=send_data(1);
lamda=send_data(2);
ref=send_data(3);
v_thr=send_data(4);

load dataForOptimize.mat data
load actual_dataForOptimize.mat actual_data

sum_diff_FST = 0;

start = 0.00;
finish = 200.0;
point = 401;

for d_count = 1:6
% clear dataMain
    x = data(:,7)*1000;
    v = data(:,d_count);
    F = griddedInterpolant(x,v,'spline');
    xq = linspace(start,finish,point);
    vq = F(xq);

    dataMain(:,1) = vq; %SED
    dataMain(:,2) = xq; %Time
    dataMain(1:point,3) = 0.50; %dt
    dataMain(:,4) = linspace(start,finish/1000,point);


    [firingRate,spikeTime,spikePosition,spikeMarker,time,membPotential,time2,cur2,dataSED2] = spike(dataMain);
        
    spikeTime = spikeTime(~isnan(spikeTime) & spikeTime>0);

    %--- 3) Extract first spike timing (model)
    if isempty(spikeTime)
        fst_model(d_count,1) = 300;
    else
        fst_model(d_count,1) = spikeTime(1);
    end
    
    %--- 4) Extract first spike timing (actual)
    adata = actual_data(:,d_count);
    adata = adata(~isnan(adata) & adata>0);
    if isempty(adata)
        fst_actual(d_count,1) = 300;
    else
        fst_actual(d_count,1) = adata(1);
    end
    

end


    RSS = sum((fst_actual - fst_model).^2);
    TSS = sum(fst_actual.^2);
    sum_diff_FST = - (1 - RSS/TSS);


end


