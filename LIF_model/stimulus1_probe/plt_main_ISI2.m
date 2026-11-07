function sum_FSS = plt_main_ISI(send_data)
    global R;R = 66.11e6;
    global C;C = 1.01e-6;
    global tau;tau = 66.77; % tau=R*C
    global M;M = 1.0; % number of merkels
    global u_thr; theta = 0.13; dt = 0.0125; u_thr = theta * dt; % revised version; v1.0: u_thr = 1.1, dt = 0.5 ms
    %**********************************************
    global beta;
    global lamda;
    global ref;
    global v_thr;

    
    beta=send_data(1);
    lamda=send_data(2);
    ref=2.1; % T_ref fixed at 2.1 ms (revised version)
    v_thr=send_data(3);
    
    load dataForOptimize.mat data
    load actual_dataForOptimize.mat actual_data
    
    sum_diff_ISI = 0;
    sum_FSS = 0;
    FSS=0;
    
    start = 0.00;
    finish = 200.0;
    point = round((finish-start)/dt) + 1; % v1.0: 401
    
    for d_count = 1:6
    % clear dataMain
        x = data(:,7)*1000;
        v = data(:,d_count);
        F = griddedInterpolant(x,v,'spline');
        xq = linspace(start,finish,point);
        vq = F(xq);
    
        dataMain(:,1) = vq; %SED
        dataMain(:,2) = xq; %Time
        dataMain(1:point,3) = dt; %dt
        dataMain(:,4) = linspace(start,finish/1000,point);
    
    
        [firingRate,spikeTime,spikePosition,spikeMarker,time,membPotential,time2,cur2,dataSED2] = spike(dataMain);
            
        tm = spikeTime(~isnan(spikeTime)&spikeTime>0);
        ISI_m = diff(tm);  % [ms]
        % pred = 1 ./ ISI_m;
        fs(d_count,1) = 1/mean(ISI_m);

        % 実測データから ISI
        ta    = actual_data(:,d_count);
        ta    = ta(~isnan(ta)&ta>0);
        ISI_a = diff(ta);  % [ms]
        fo(d_count,1) = 1/mean(ISI_a);

    end


    RSS = sum((fo - fs).^2);
    TSS = sum(fo.^2);
    sum_FSS = - (1 - RSS/TSS);

end
    