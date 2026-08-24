function a = func(tn,vn,cur,i)

global R;
global C;
global tau;

a = (1/tau)*(R*cur(i,1) - vn);

end

