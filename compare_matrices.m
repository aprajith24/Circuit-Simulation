function[conv] = compare_matrices(VOL_a,VOL_b,e)
no_of_rows=size(VOL_a,1);
conv=1
for i = 1:no_of_rows
    i=i
    y=abs(VOL_a(i,1) - VOL_b(i,1))
    if(abs(VOL_a(i,1) - VOL_b(i,1)) > e),
        conv=0;
        break;
    end
end
