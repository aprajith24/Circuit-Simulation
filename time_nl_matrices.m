function[VOL_new] = time_nl_matrices(time_M,time_I,VOL,NLNELEM,h,NL_I,e)
global TYPE_ M_ M_ND_ M_NG_ M_NS_ M_TYPE_ M_VT_ M_MU_ M_W_ M_L_ M_COX_ M_LAMBDA_ M_CJ0_ C_
no_of_rows=size(NLNELEM,1);
temp_M = time_M;
temp_I = time_I;
conv=0;
VOL_pre=VOL

while(conv == 0),
    time_M=temp_M;
    time_I=temp_I;
    upd=0;
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
        cj0 = NLNELEM(n,M_CJ0_);
        cgs = (1/2)*cox*w*l;
        cgd = (1/2)*cox*w*l;
        cd = cj0;
        cs = cj0;
 
        upd = upd +1;
        curr = NL_I(upd,1);
        C_GS = [C_ cgs g s];
        [time_M,time_I]= stamp_capacitance(time_M,time_I,C_GS,VOL,h,curr);
          
        upd = upd +1;
        curr = NL_I(upd,1);
        C_GD = [C_ cgd g d];
        [time_M,time_I]= stamp_capacitance(time_M,time_I,C_GD,VOL,h,curr);
        
        upd = upd +1;
        curr = NL_I(upd,1);
        C_D = [C_ cd d -1];
        [time_M,time_I]= stamp_capacitance(time_M,time_I,C_D,VOL,h,curr);

        upd = upd +1;
        curr = NL_I(upd,1);
        C_S = [C_ cs s -1];
        [time_M,time_I]= stamp_capacitance(time_M,time_I,C_S,VOL,h,curr);
        
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
            [time_M,time_I] = stamp_pmos(time_M,time_I,NLNELEM(n,:),gm,gsd,i,VOL);
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
            [time_M,time_I] = stamp_nmos(time_M,time_I,NLNELEM(n,:),gm,gds,i,VOL);
        end
    end
    end
    VOL_new = inv(time_M) * time_I
    [conv] = compare_matrices(VOL_new,VOL,e);
    VOL=VOL_new
end
VOL=VOL_pre;