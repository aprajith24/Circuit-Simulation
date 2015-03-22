function [CL_I] = initialize_CL_I(VOL,CL_I,LINELEM)
global C_ L_ TYPE_;
no_of_rows = size(LINELEM,1);
i = 0;
for n = 1:no_of_rows
    if (LINELEM(n,TYPE_) == C_),
        i = i+1;
        CL_I(i,1) = 0;
    end
    if (LINELEM(n,TYPE_) == L_),
        i = i+1;
        vol_row=CL_I(i,2);
        CL_I(i,1) = VOL(vol_row,1);
    end
end
    