function [new_I] = stamp_ind_isource(old_I,D);
%STAMP_IND_ISOURCE : stamps entries corresponding to an independent current source.
%
% syntax : [new_I,new_row] = stamp_ind_isource(old_I,D)
%
% new_I,old_I,are the new and old current matrices (right hand side)
% D is the data vector corresponding to the source
global I_N1_ I_N2_ I_ I_VALUE_
new_I=old_I;
length_I=length(old_I);

n1 = D(I_N1_);
n2 = D(I_N2_);
if n1>length_I, new_I(n1)=0;end;
if n2>length_I, new_I(n2)=0;end;

if n1>0,new_I(n1)=-D(I_VALUE_);end
if n2>0,new_I(n2)=D(I_VALUE_);end
