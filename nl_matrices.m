function[VOL_out,NL_I] = nl_matrices(M,I,NLNELEM,NL_I,e)
global TYPE_ M_ M_ND_ M_NG_ M_NS_ M_TYPE_ M_VT_ M_MU_ M_W_ M_L_ M_COX_ M_LAMBDA_ Y_
no_of_rows=size(NLNELEM,1);
count = 1;
temp_M = M;
temp_I = I;
VOL=[0];
gmin=1e-12;
conv=0;
cnt =0;

for n = 1:no_of_rows
    if (NLNELEM(n,TYPE_) == M_),
        if(NLNELEM(n,M_ND_)>0),
            VOL(NLNELEM(n,M_ND_),1)=0;
        end
        if(NLNELEM(n,M_NG_)>0),
            VOL(NLNELEM(n,M_NG_),1)=0;
        end
        if(NLNELEM(n,M_NS_)>0),
            VOL(NLNELEM(n,M_NS_),1)=0;
        end
        NL_I(count,1)=0;
        count = count+1;
        NL_I(count,1)=0;
        count = count+1;
        NL_I(count,1)=0;
        count = count+1;
        NL_I(count,1)=0;
        count = count+1;
    end
end

while(conv == 0),
    M=temp_M;
    I=temp_I;
    for n = 1:no_of_rows
    if (NLNELEM(n,TYPE_) == M_),
        d = NLNELEM(n,M_ND_);
        g = NLNELEM(n,M_NG_);
        s = NLNELEM(n,M_NS_);
        w = NLNELEM(n,M_W_);
        l = NLNELEM(n,M_L_);
        mu = NLNELEM(n,M_MU_);
        cox = NLNELEM(n,M_COX_);
        vt = NLNELEM(n,M_VT_);
        lam = NLNELEM(n,M_LAMBDA_);
        
        COND_D = [Y_ gmin d -1];
        [M] = stamp_conductance(M,COND_D);
        COND_G = [Y_ gmin g -1];
        [M] = stamp_conductance(M,COND_G);
        COND_S = [Y_ gmin s -1];
        [M] = stamp_conductance(M,COND_S);

        if(NLNELEM(n,M_TYPE_) == 0),
            if(g>0 && s>0),
                vsg = VOL(s,1) - VOL(g,1);
            end
            if(s>0 && g<0),
                vsg = VOL(s,1);
            end
            if(s<0 && g>0),
                vsg = -VOL(g,1);
            end            
            if(s>0 && d>0),
                vsd = VOL(s,1) - VOL(d,1);
            end
            if(s>0 && d<0),
                vsd = VOL(s,1);
            end
            if(s<0 && d>0),
                vsd = -VOL(d,1);
            end
            if (vsg <= -vt),
                gm = 0;
                gsd = 0;
                i=0;
            end
            if((vsg > -vt) && (vsd <= vsg + vt)),
                syms v_sg
                syms v_sd
                gm  = diff(mu*cox*(w/l)*((v_sg + vt)*vsd - (1/2)*(vsd^2))*(1 +lam*vsd),v_sg);
                gm  = vpa(subs(gm,v_sg,vsg));
                gsd = diff(mu*cox*(w/l)*((vsg + vt)*v_sd - (1/2)*(v_sd^2))*(1 +lam*v_sd),v_sd);
                gsd = vpa(subs(gsd,v_sd,vsd));
                i= mu*cox*(w/l)*((vsg + vt)*vsd - (1/2)*(vsd^2))*(1 +lam*vsd);
            end
            if((vsg > -vt) && (vsd > vsg + vt)),
                syms v_sg
                syms v_sd
                gm  = diff((1/2)*mu*cox*(w/l)*((v_sg + vt)^2)*(1 +lam*vsd),v_sg);
                gm  = vpa(subs(gm,v_sg,vsg));
                gsd = diff((1/2)*mu*cox*(w/l)*((vsg + vt)^2)*(1 +lam*v_sd),v_sd);
                gsd = vpa(subs(gsd,v_sd,vsd));
                i= (1/2)*mu*cox*(w/l)*((vsg + vt)^2)*(1 +lam*vsd);
            end
            [M,I] = stamp_pmos(M,I,NLNELEM(n,:),gm,gsd,i,VOL);
        elseif(NLNELEM(n,M_TYPE_) == 1),
            if(g>0 && s>0),
                vgs = VOL(g,1) - VOL(s,1);
            end
            if(g>0 && s<0),
                vgs = VOL(g,1);
            end
            if(g<0 && s>0),
                vgs = -VOL(s,1);
            end            
            if(d>0 && s>0),
                vds = VOL(d,1) - VOL(s,1);
            end
            if(d>0 && s<0),
                vds = VOL(d,1);
            end
            if(d<0 && s>0),
                vds = -VOL(s,1);
            end
            if (vgs <= vt),
                gm = 0;
                gds = 0;
                i=0;
            end
            if((vgs > vt) && (vds <= vgs - vt)),
                syms v_gs
                syms v_ds
                gm  = diff(mu*cox*(w/l)*((v_gs - vt)*vds - (1/2)*(vds^2))*(1 +lam*vds),v_gs);
                gm  = vpa(subs(gm,v_gs,vgs));
                gds = diff(mu*cox*(w/l)*((vgs - vt)*v_ds - (1/2)*(v_ds^2))*(1 +lam*v_ds),v_ds);
                gds = vpa(subs(gds,v_ds,vds));
                i= mu*cox*(w/l)*((vgs - vt)*vds - (1/2)*(vds^2))*(1 +lam*vds);
            end
            if((vgs > vt) && (vds > vgs - vt)),
                syms v_gs
                syms v_ds
                gm  = diff((1/2)*mu*cox*(w/l)*((v_gs - vt)^2)*(1 +lam*vds),v_gs);
                gm  = vpa(subs(gm,v_gs,vgs));
                gds = diff((1/2)*mu*cox*(w/l)*((vgs - vt)^2)*(1 +lam*v_ds),v_ds);
                gds = vpa(subs(gds,v_ds,vds));
                i= (1/2)*mu*cox*(w/l)*((vgs - vt)^2)*(1 +lam*vds);
            end            
            [M,I] = stamp_nmos(M,I,NLNELEM(n,:),gm,gds,i,VOL);
        end
    end
    end
    VOL_out = inv(M) * I;
    [conv] = compare_matrices(VOL_out,VOL,e);
    VOL=VOL_out;
    cnt = cnt +1;
end