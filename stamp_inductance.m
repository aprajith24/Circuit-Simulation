function [time_M,time_I] = stamp_inductance(time_M,time_I,D,VOL,h,i)
global L_N1_ L_N2_ L_VALUE_ Y_ I_;
n1 = D(L_N1_);
n2 = D(L_N2_);
value = D(L_VALUE_);
geq= h/(2*value);
if (n1>0) && (n2>0),
    ieq=(h*(VOL(n1,1)-VOL(n2,1))/(2*value))+i;
elseif (n1<0)
    ieq=(h*(-VOL(n2,1))/(2*value))+i;
elseif (n2<0)
    ieq=(h*(VOL(n1,1))/(2*value))+i;
end

cond = [Y_ geq n1 n2];
[time_M] = stamp_conductance(time_M,cond);

curr = [I_ ieq n1 n2];
[time_I] = stamp_ind_isource(time_I,curr);