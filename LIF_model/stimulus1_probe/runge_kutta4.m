function dv = runge_kutta4(tn,vn,cur,i,data)

% global data;

k1 = func(tn,vn,cur,i);
k2 = func(tn+(data(i,3)/2.0),vn+(data(i,3)/2.0)*k1,cur,i);
k3 = func(tn+(data(i,3)/2.0),vn+(data(i,3)/2.0)*k2,cur,i);
k4 = func(tn+data(i,3),vn+data(i,3)*k3,cur,i);

dv = (data(i,3)/6.0)*( k1 + 2.0*k2 + 2.0*k3 + k4 );

end

