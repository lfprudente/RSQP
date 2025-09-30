function [merit] = eval_merit(x,n,m,equatn,lambda,rho,sc,scaling,LStype)

global evalfun

[~,fs,~] = sevalf(n,x,sc,scaling);
evalfun = evalfun + 1;

cs = [];
for i = 1:m
    [~,cs(i),~] = sevalc(n,x,i,sc,scaling);
    evalfun = evalfun + 1;
end
cs = cs';

merit = fs;
if LStype == 1
    for i = 1:m
        if ( equatn(i) == true || lambda(i) + rho * cs(i) > 0 )
            merit = merit + cs(i) * ( lambda(i) + 0.5 * rho * cs(i) );
        else
            merit = merit - 0.5 * lambda(i)^2 / rho;
        end
    end
elseif LStype == 2
    for i = 1:m
        if equatn(i)
            merit = merit + rho * abs(cs(i));
        else
            merit = merit + rho * max(0,cs(i));
        end
    end
end