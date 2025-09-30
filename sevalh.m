function [h,hs,flag] = sevalh(n,x,sc,scaling)

[h,flag] = evalh(n,x);

if ( scaling == true && flag == 1 )
    hs = h * sc.f;
else
    hs = h;
end