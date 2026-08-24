% Spike function for RA1
function [firing_rate,t_spike,p_spike,spike_marker,time,v,time2,cur2,data_SED2] = new_spike(data_main)

    %**************************%
    %*****  model_v7-4  *****%
    %**************************%
    
    global beta;
    global alpha;
    global lamda;
    global M;
    
    global v_thr;
    global u_thr;
    
    global ref;
    
    cur(1,1) = 0.0;
    j=0;
    for i=1:length(data_main(:,1))-1
        if (data_main(i+1,1)-data_main(i,1)) > u_thr 
            cur(i+1,1) = beta + M * ( lamda*((data_main(i+1,1)-data_main(i,1))/data_main(i,3)));
        elseif (data_main(i+1,1)-data_main(i,1)) < -u_thr
            cur(i+1,1) = beta + M * ( lamda*((data_main(i,1)-data_main(i+1,1))/data_main(i,3)));
        else
            cur(i+1,1) = 0;
        end
    end

    cur = cur + 2.5e-6*randn(size(cur));
    % cur = cur + 5.0e-6*randn(size(cur));
    
    spike_time_last = -ref;
    tn = 0.0;
    vn = 0.0;
    time(1,1) = tn;
    time2(1,1) = tn;
    v(1,1) = vn;
    cur2(1,1) = cur(1,1);
    data_SED2(1,1) = data_main(1,1);
    t_spike(1,1) = 0.0;
    firing_rate(1,1) = 0.0;
    spike_marker(1,1) = 0.0;
    %*************************************
    p_spike(1,1) = -10;
    %*************************************
    
    i = 1;
    j = 1;
    while tn < data_main(length(data_main(:,2)),2)-1
        dv = runge_kutta4(tn,vn,cur,i,data_main);
        if tn >= (spike_time_last + ref)
            vn = vn + dv;
        end
        tn = data_main(i+1,2);
        
        if vn > v_thr
            vn = 0.0;
            t_spike(j+1,1) = tn;
            p_spike(j,1) = data_main(i+1,4);
            spike_marker(j+1,1) = 1.0; 
            firing_rate(j+1,1) = 1.0/((t_spike(j+1,1)-t_spike(j,1))*0.001);
            v(i+j,1) = v_thr*2.0;
            v(i+j+1,1) = 0.0;
            time2(i+j,1) = tn;
            time2(i+j+1,1) = tn;
            j = j+1;
    
            spike_time_last = tn;
        else
            v(i+j,1) = vn;
            time2(i+j,1) = tn;
        end
        
        time(i+1,1) = tn;
        cur2(i+1,1) = cur(i+1,1);
        data_SED2(i+1,1) = data_main(i+1,1);
        
        i = i+1;
        
    end
    
    
    end
    
    