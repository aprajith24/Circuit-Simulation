function [time_M,time_I] = stamp_capacitance(time_M,time_I,D,VOL,h,i)
global C_N1_ C_N2_ C_VALUE_ Y_ I_;
n1 = D(C_N1_);
n2 = D(C_N2_);
value = D(C_VALUE_);
geq= 2*value/h;
ieq=0;
if (n1>0) && (n2>0),
    ieq=(2*value*(VOL(n1,1)-VOL(n2,1))/h)+i;
elseif (n1<0) && (n2>0),
    ieq=(2*value*(-VOL(n2,1))/h)+i;    
elseif (n2<0) && (n1>0),
    ieq=(2*value*(VOL(n1,1))/h)+i;
end    

cond = [Y_ geq n1 n2];
[time_M] = stamp_conductance(time_M,cond);

curr = [I_ ieq n2 n1];
[time_I] = stamp_ind_isource(time_I,curr);