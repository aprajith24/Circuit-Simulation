function[time_M,time_I] = time_matrices(time_M,time_I,VOL,LINELEM,h,CL_I,time_now)
global TYPE_ Y_ R_ I_ V_ L_ C_;
no_of_rows=size(LINELEM,1);
count = 0;
for n = 1:no_of_rows
    if (LINELEM(n,TYPE_) == Y_),
        [time_M]= stamp_comductance(time_M,LINELEM(n,:));
    end
    if (LINELEM(n,TYPE_) == R_),
        [time_M]= stamp_resistance(time_M,LINELEM(n,:));
    end
    if (LINELEM(n,TYPE_) == I_),
        [time_I] = stamp_ind_isource(time_I,LINELEM(n,:));
    end
    if (LINELEM(n,TYPE_) == C_),
        count = count +1;
        i = CL_I(count,1);
        [time_M,time_I]= stamp_capacitance(time_M,time_I,LINELEM(n,:),VOL,h,i);
    end
    if (LINELEM(n,TYPE_) == L_),
        count = count +1;
        i = CL_I(count,1);
        [time_M,time_I] = stamp_inductance(time_M,time_I,LINELEM(n,:),VOL,h,i);
    end
end

for n = 1:no_of_rows
    if (LINELEM(n,TYPE_) == V_),
        D = LINELEM(n,:);
        [vin] = voltage_pwl(D,time_now);
        D(1,2) = vin;
        [time_M,time_I,new_row] = stamp_ind_vsource(time_M,time_I,D);
    end
end


