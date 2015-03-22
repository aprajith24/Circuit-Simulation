function [new_M,new_I] = stamp_pmos(old_M,old_I,D,gm,gsd,i,VOL)
global G_ Y_ M_ND_ M_NG_ M_NS_ I_
d = D(1,M_ND_);
g = D(1,M_NG_);
s = D(1,M_NS_);
if(g>0 && s>0),
    vsg_old = VOL(s,1) - VOL(g,1);
end
if(s>0 && g<0),
    vsg_old = VOL(s,1);
end
if(s<0 && g>0),
    vsg_old = -VOL(g,1);
end            
if(s>0 && d>0),
    vsd_old = VOL(s,1) - VOL(d,1);
end
if(s>0 && d<0),
    vsd_old = VOL(s,1);
end
if(s<0 && d>0),
    vsd_old = -VOL(d,1);
end
i_src = i - (gm*vsg_old) - (gsd*vsd_old);

VOL_SRC=[G_ gm s d s g];
[new_M] = stamp_vctrl_isource(old_M,VOL_SRC);
old_M = new_M;
COND=[Y_ gsd s d];
[new_M] = stamp_conductance(old_M,COND);
I_SRC = [I_ i_src s d];
[new_I] = stamp_ind_isource(old_I,I_SRC);
