function[NL_I] = update_NL_I(NL_I,VOL,VOL_new,NLNELEM,h)
global TYPE_ M_ M_ND_ M_NG_ M_NS_ M_W_ M_L_ M_COX_ M_CJ0_
no_of_rows=size(NLNELEM,1);
upd = 0;

for n = 1:no_of_rows
    if (NLNELEM(n,TYPE_) == M_),
        d = NLNELEM(n,M_ND_);
        g = NLNELEM(n,M_NG_);
        s = NLNELEM(n,M_NS_);
        w = NLNELEM(n,M_W_);
        l = NLNELEM(n,M_L_);
        cox = NLNELEM(n,M_COX_);
        cj0 = NLNELEM(n,M_CJ0_);
        
        cgs = (1/2)*cox*w*l;
        cgd = (1/2)*cox*w*l;
        cd = cj0;
        cs = cj0;
        
        upd = upd+1;
        i_old = NL_I(upd,1);
        if (g>0) && (s>0),
            NL_I(upd,1) = (2*cgs*(VOL_new(g,1)-VOL_new(s,1)-(VOL(g,1)-VOL(s,1)))/h) - i_old;
        elseif(g<0),
            NL_I(upd,1) = (2*cgs*(-VOL_new(s,1)-(-VOL(s,1)))/h) - i_old;
        elseif(s<0),
            NL_I(upd,1) = (2*cgs*(VOL_new(g,1)-(VOL(g,1)))/h) - i_old;
        end
        
        upd = upd+1;
        i_old = NL_I(upd,1);
        if (g>0) && (d>0),
            NL_I(upd,1) = (2*cgd*(VOL_new(g,1)-VOL_new(d,1)-(VOL(g,1)-VOL(d,1)))/h) - i_old;
        elseif(g<0),
            NL_I(upd,1) = (2*cgd*(-VOL_new(d,1)-(-VOL(d,1)))/h) - i_old;
        elseif(d<0),
            NL_I(upd,1) = (2*cgd*(VOL_new(g,1)-(VOL(g,1)))/h) - i_old;
        end

        upd = upd+1;
        i_old = NL_I(upd,1);
        if(d>0),
            NL_I(upd,1) = (2*cd*(VOL_new(d,1)-(VOL(d,1)))/h) - i_old;
        end

        upd = upd+1;
        i_old = NL_I(upd,1);
        if(s>0),
            NL_I(upd,1) = (2*cs*(VOL_new(s,1)-(VOL(s,1)))/h) - i_old;
        end
    end
end