function[CL_I] = update_CL_I(CL_I,VOL,VOL_new,LINELEM,h)

global C_N1_ C_N2_ C_VALUE_ C_ L_ L_VALUE_ L_N1_ L_N2_ TYPE_;
no_of_rows=size(LINELEM,1);
count = 0;

for n = 1:no_of_rows
    if (LINELEM(n,TYPE_) == C_),
        count = count +1;
        n1 = LINELEM(n,C_N1_);
        n2 = LINELEM(n,C_N2_);
        value = LINELEM(n,C_VALUE_);
        i_old = CL_I(count,1);
        if (n1>0) && (n2>0),
            CL_I(count,1) = (2*value*(VOL_new(n1,1)-VOL_new(n2,1)-(VOL(n1,1)-VOL(n2,1)))/h) - i_old;
        elseif(n1<0),
            CL_I(count,1) = (2*value*(-VOL_new(n2,1)-(-VOL(n2,1)))/h) - i_old;
        elseif(n2<0),
            CL_I(count,1) = (2*value*(VOL_new(n1,1)-(VOL(n1,1)))/h) - i_old;
        end
    end
    if (LINELEM(n,TYPE_) == L_),
        count = count +1;
        n1 = LINELEM(n,L_N1_);
        n2 = LINELEM(n,L_N2_);
        value = LINELEM(n,L_VALUE_);
        i_old = CL_I(count,1);
        if (n1>0) && (n2>0),
            CL_I(count,1) = (h*(VOL_new(n1,1)-VOL_new(n2,1)+(VOL(n1,1)-VOL(n2,1)))/(2*value)) + i_old;
        elseif(n1<0),
            CL_I(count,1) = (h*(-VOL_new(n2,1)+(-VOL(n2,1)))/(2*value)) + i_old;
        elseif(n2<0),
            CL_I(count,1) = (h*(VOL_new(n1,1)+(VOL(n1,1)))/(2*value)) + i_old;
        end
    end
end