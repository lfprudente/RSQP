function [hc,hcs,flag] = sevalhc(n,x,ind,sc,scaling)

[hc,flag] = evalhc(n,x,ind);

if ( scaling == true && flag == 1 )
    hcs = hc * sc.c(ind);
else
    hcs = hc;
end