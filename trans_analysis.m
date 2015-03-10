function [M,I,VOL,VOL_new, CL_I, time_M, time_I]= trans_analysis(LINELEM)
M = [0];
I=[0;0];
VOL=[0];
CL_I =[0];
[M,I,CL_I] = matrices(M,I,LINELEM,CL_I);
VOL = inv(M) * I;

[CL_I] = initialize_CL_I(VOL,CL_I,LINELEM);
h = 1.000000000000000e-11;
stop = 3.000000000000000e-09;
time_now = 0;
TIME = [0];
VOLTAGE1 = [VOL(2,1)];
VOLTAGE2 = [VOL(3,1)];

for a= h:h:100*h
    time_now = time_now+h
    time_M = [0];
    time_I = [0;0];
    VOL_new = [0];
    [time_M,time_I] = time_matrices(time_M,time_I,VOL,LINELEM,h,CL_I,time_now);
    VOL_new = inv(time_M) * time_I;
    [CL_I] = update_CL_I(CL_I,VOL,VOL_new,LINELEM,h);
    VOL = VOL_new;
    TIME = [TIME,time_now];
    VOLTAGE1 = [VOLTAGE1,VOL(20,1)];
    VOLTAGE2 = [VOLTAGE2,VOL(3,1)];

end
plot(TIME,VOLTAGE1);