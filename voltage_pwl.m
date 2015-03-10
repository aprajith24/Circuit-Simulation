function [vin] = voltage_pwl(D,time_now);
columns = size(D,2);
if(columns<7),
    vin = D(1,2);
else
    for n = 7:2:columns
        if(time_now <= D(1,n)),
            if(n>7),
                slope = (D(1,n+1) - D(1,n-1))/(D(1,n) - D(1,n-2));
            else
                slope = (D(1,n+1) - D(1,2))/D(1,n);
            end
            vin = D(1,n+1) + slope*(time_now - D(1,n));
            break;
        end
    end
    if(time_now > D(1,columns-1)),
        vin = D(1,columns);
    end
end

    
    

