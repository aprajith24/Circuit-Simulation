function [new_M] = stamp_vctrl_isource(old_M,D);
%STAMP_IND_ISOURCE : stamps entries corresponding to an independent current source.
%
% syntax : [new_I,new_row] = stamp_ind_isource(old_I,D)
%
% new_I,old_I,are the new and old current matrices (right hand side)
% D is the data vector corresponding to the source
global G_N1_ G_N2_ G_CN1_ G_CN2_ G_VALUE_
new_M=old_M;
length_M=length(old_M);

n1 = D(G_N1_);
n2 = D(G_N2_);
cn1 = D(G_CN1_);
cn2 = D(G_CN2_);
value=D(G_VALUE_);

if n1>length_M, new_M(n1,n1)=0;end;
if n2>length_M, new_M(n2,n2)=0;end;
if cn1>length_M, new_M(cn1,cn1)=0;end;
if cn2>length_M, new_M(cn2,cn2)=0;end;


if (n1>0) | (cn1>0),
 new_M(n1,cn1) = new_M(n1,cn1) + value;
end
if (n1>0) | (cn2>0),
 new_M(n1,cn2) = new_M(n1,cn2) - value;
end
if (n2>0) | (cn1>0),
 new_M(n2,cn1) = new_M(n2,cn1) - value;
end
if (n2>0) | (cn2>0),
 new_M(n2,cn2) = new_M(n2,cn2) + value;
end 