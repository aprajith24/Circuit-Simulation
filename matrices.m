function[M,I,CL_I] = matrices(M,I,LINELEM,CL_I)
global TYPE_ Y_ R_ I_ V_ L_ L_N1_ L_N2_ C_;
no_of_rows=size(LINELEM,1);
count = 0;
for n = 1:no_of_rows
    if (LINELEM(n,TYPE_) == Y_),
        [M]= stamp_conductance(M,LINELEM(n,:));
    end
    if (LINELEM(n,TYPE_) == R_),
        [M]= stamp_resistance(M,LINELEM(n,:));
    end
    if (LINELEM(n,TYPE_) == I_),
        [I] = stamp_ind_isource(I,LINELEM(n,:));
    end
end

for n = 1:no_of_rows
    if (LINELEM(n,TYPE_) == V_),
        [M,I,new_row] = stamp_ind_vsource(M,I,LINELEM(n,:));
    end
    if (LINELEM(n,TYPE_) == L_),
        count = count +1;
        TEMP = [V_ 0 LINELEM(n,L_N1_) LINELEM(n,L_N2_)];
        [M,I,new_row] = stamp_ind_vsource(M,I,TEMP);
        CL_I(count,2) = new_row;
    end
    if (LINELEM(n,TYPE_) == C_),
        count = count +1;
    end
end


