function [firing_rate,t_spike,p_spike,spike_marker,time,v,time2,cur2,data_SED2] = spike(data_main)

%**************************%
%*****  model_v7-4  *****%
%**************************%

global beta;
global alpha;
global lamda;
global M;

global v_thr;

cur(1,1) = 0.0;
ttau = 0.0;
curp = 0.0;
for i=1:length(data_main(:,1))-1%*****'f54-1.xlsx','sheet1'/ (A—ñs”)-1
    if (data_main(i+1,1)-data_main(i,1)) > 0
        cur(i+1,1) = beta + M * ( alpha*data_main(i+1,1) + lamda*((data_main(i+1,1)-data_main(i,1))/data_main(i,3)) );
        ttau = data_main(i+1,2);
        curp = cur(i+1,1);
    elseif (data_main(i+1,1)-data_main(i,1)) == 0
        cur(i+1,1) = (beta + alpha*data_main(i+1,1)) + (curp-(beta + alpha*data_main(i+1,1)))*exp(-(data_main(i+1,2)-ttau)/100);
    elseif (data_main(i+1,1)-data_main(i,1)) < 0
        cur(i+1,1) = (beta + alpha*data_main(i+1,1));
        ttau = data_main(i+1,2);
        curp = cur(i+1,1);
    end
end


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
while tn < data_main(length(data_main(:,2)),2)-1%*****'f54-1.xlsx','sheet1'/ (B—ñÅIs‚ÌŠÔ)-1
    
    dv = runge_kutta4(tn,vn,cur,i,data_main);
    vn = vn + dv;
    tn = data_main(i+1,2);
    
    if vn > v_thr
        vn = 0.0;
        t_spike(j+1,1) = tn;
        p_spike(j,1) = data_main(i+1,4);
        spike_marker(j+1,1) = 1.0; 
        firing_rate(j+1,1) = 1.0/(t_spike(j+1,1)*0.001-t_spike(j,1)*0.001);
        v(i+j,1) = v_thr*2.0;
        v(i+j+1,1) = 0.0;
        time2(i+j,1) = tn;
        time2(i+j+1,1) = tn;
        j = j+1;
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

