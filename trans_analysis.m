function [M,I,VOL,CL_I,time_M,time_I]= trans_analysis(LINELEM,NLNELEM,INFO,PLOTNV)
M = [0];
I=[0;0];
VOL=[0];
CL_I =[0];
NL_I= [0];
e=1e-2;
[M,I,CL_I] = matrices(M,I,LINELEM,CL_I);
if(size(NLNELEM,1)>0),
    [VOL,NL_I] = nl_matrices(M,I,NLNELEM,NL_I,e);
end
if(size(NLNELEM,1)==0),
    VOL = inv(M) * I;
end

[CL_I] = initialize_CL_I(VOL,CL_I,LINELEM);
h = INFO(2,1);
stop = INFO(3,1);
time_now = 0;
TIME = (0);
no_of_plots = size(PLOTNV,1);
for plots = 1:no_of_plots
    VOLTAGE{plots} = VOL(PLOTNV(plots,1),1);
end

for a= h:h:stop
    time_now = time_now+h
    time_M = (0);
    time_I = [0;0];
    VOL_new = (0);
    [time_M,time_I] = time_matrices(time_M,time_I,VOL,LINELEM,h,CL_I,time_now);
    if(size(NLNELEM,1)>0),
        [VOL_new] = time_nl_matrices(time_M,time_I,VOL,NLNELEM,h,NL_I,e);
    end
    if(size(NLNELEM,1)==0),
        VOL_new = inv(time_M) * time_I;
    end
    [CL_I] = update_CL_I(CL_I,VOL,VOL_new,LINELEM,h);
    if(size(NLNELEM,1)>0),
        [NL_I] = update_NL_I(NL_I,VOL,VOL_new,NLNELEM,h);
    end
    VOL = VOL_new;
    TIME = [TIME,time_now];
    for plots = 1:no_of_plots
        VOLTAGE{plots} = [VOLTAGE{plots},VOL(PLOTNV(plots,1),1)];
    end
end
for plots = 1:no_of_plots
    figure(plots);
    plot(TIME,VOLTAGE{plots});
    str = sprintf('Plot of Node %d',PLOTNV(plots,1));
    title(str);
    xlabel('Time'); % x-axis label
    ylabel('Voltage'); % y-axis label
end
