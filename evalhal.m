function [eh,store] = evalhal(x,u,store,n,m,equatn,rho,lambar,sc,scaling)

global evalgrad evalhess nevalhal 

nevalhal = nevalhal + 1;

% Ensure u is an ambient vector (packed)
if isstruct(u)
    u_vec = reshapevector(u);   % struct -> vetor
else
    u_vec = u(:);               % já é vetor, apenas força coluna
end

% ----- Objective Hessian-vector -----
[~, H,flag] = sevalh(n,x,sc,scaling);
if flag == 1 
    evalhess = evalhess + 1;
    % Start with the objective contribution
    eh = H * u_vec;
else
    eh = zeros(n,1);
end

% ----- Cache constraint values c_i(x) in store (as in evalal/evalnal) -----
if ~isfield(store, 'c')
  cs = zeros(m,1);
  for i = 1:m
      [~, cs(i), ~] = sevalc(n, x, i, sc, scaling);
  end
  store.c = cs;
end
cs = store.c;

haveJc = isfield(store, 'Jc');
if haveJc
    Jc = store.Jc;          % m-by-n
else
    Jc = zeros(m, n);       % prealoca
end

% ----- Constraint contributions (only when active) -----
for i = 1:m
    active_i = (equatn(i) == true) || (lambar(i) + rho * cs(i) > 0);
    if ~active_i
      continue;
    end
    
    % Euclidean gradient ∇c_i(x)
    if ~haveJc
        [~,nci,~] = sevalnc(n,x,i,sc,scaling);    evalgrad = evalgrad + 1;
        Jc(i, :) = nci(:).'; 
    else
        nci = Jc(i, :).';
    end
    
    % Hessian-vector ∇^2 c_i(x)[u]
    [~, Hci,flag] = sevalhc(n, x, i, sc, scaling);
    if flag == 1 
        evalhess = evalhess + 1;
        Hci_u = Hci * u_vec;
    else
        Hci_u = zeros(n,1);
    end
    
    % Add (λ_i + ρ c_i) ∇^2 c_i(x)[u]  +  ρ 〈∇c_i(x), u〉 ∇c_i(x)
    eh = eh + (lambar(i) + rho * cs(i)) * Hci_u ...
           + rho * (nci.' * u_vec) * nci;
end

if ~haveJc
    store.Jc = Jc;
end

eh = reshapevector(eh);  