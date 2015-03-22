function [new_M,new_I] = stamp_nmos(old_M,old_I,D,gm,gds,i,VOL)
global G_ Y_ M_ND_ M_NG_ M_NS_ I_
d = D(1,M_ND_);
g = D(1,M_NG_);
s = D(1,M_NS_);
if(g>0 && s>0),
    vgs_old = VOL(g,1) - VOL(s,1);
end
if(g>0 && s<0),
    vgs_old = VOL(g,1);
end
if(g<0 && s>0),
    vgs_old = -VOL(s,1);
end            
if(d>0 && s>0),
    vds_old = VOL(d,1) - VOL(s,1);
end
if(d>0 && s<0),
    vds_old = VOL(d,1);
end
if(d<0 && s>0),
    vds_old = -VOL(s,1);
end
i_src = i - (gm*vgs_old) - (gds*vds_old);

VOL_SRC=[G_ gm d s g s];
[new_M] = stamp_vctrl_isource(old_M,VOL_SRC);
old_M = new_M;
COND=[Y_ gds d s];
[new_M] = stamp_conductance(old_M,COND);
I_SRC = [I_ i_src d s];
[new_I] = stamp_ind_isource(old_I,I_SRC);
