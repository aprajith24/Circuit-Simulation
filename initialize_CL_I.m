function [CL_I] = initialize_CL_I(VOL,CL_I,LINELEM)
global C_ C_N1_ C_N2_ C_VALUE_ L_ L_N1_ L_N2_ L_VALUE_ TYPE_;
no_of_rows = size(LINELEM,1);
i = 0;
for n = 1:no_of_rows
    if (LINELEM(n,TYPE_) == C_),
        n1 = LINELEM(n,C_N1_);
        n2 = LINELEM(n,C_N2_);
        value = LINELEM(n,C_VALUE_);
        i = i+1;
        CL_I(i,1) = 0;
    end
    if (LINELEM(n,TYPE_) == L_),
        n1 = LINELEM(n,L_N1_);
        n2 = LINELEM(n,L_N2_);
        value = LINELEM(n,L_VALUE_);
        i = i+1;
        vol_row=CL_I(i,2);
        CL_I(i,1) = VOL(vol_row,1);
    end
end
    