function [nal,store] = evalnal(x,store,n,m,equatn,rho,lambar,sc,scaling)

global evalfun evalgrad nevalnal 

nevalnal = nevalnal + 1;

[~,gs,~] = sevalg(n,x,sc,scaling);
evalgrad = evalgrad + 1;

if ~isfield(store, 'c')
    cs = [];
    for i = 1:m
        [~,cs(i),~] = sevalc(n,x,i,sc,scaling);    evalfun = evalfun + 1;
    end
    cs = cs';

    store.c = cs;       
end

cs = store.c;

haveJc = isfield(store, 'Jc');
if haveJc
    Jc = store.Jc;
else
    Jc = zeros(m, n);
end

nal = [];

nal = gs;
for i = 1:m
    if ( equatn(i) == true || lambar(i) + rho * cs(i) > 0 )
        if ~haveJc
            [~,ncs,~] = sevalnc(n,x,i,sc,scaling);    evalgrad = evalgrad + 1;
            Jc(i,:) = ncs(:).';
        else
            ncs = Jc(i,:).';
        end

        nal = nal + ( lambar(i) + rho * cs(i) ) * ncs;
    end
end

if ~haveJc
    store.Jc = Jc;
end

[nal] = reshapevector(nal);