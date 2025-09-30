function H = make_SPD(H)

% MAKE_SPD  SPD regularization for a symmetric-ish matrix.

global num_eig num_chol

d_eig = 2000;

% (0) Symmetrize.
H = (H + H.')/2;
d = size(H, 1);

 % Try raw Cholesky first: if already SPD, done.
if min(diag(H)) > 1e-12
    num_chol = num_chol + 1;
    
    [~,p] = chol(H);
    if p == 0, return; end
end

 % If the dimension less than or equal to d_eig, compute the eigenvalue
 % decomposition and correct the negative eigenvalues
if d <= d_eig 
    num_eig = num_eig + 1;

    [V, D] = eig(H,'vector');
    D = max(D, 1e-6);            % clip
    H = V*diag(D)*V.';
    H = (H + H.')/2;
    return
end

% s     = norm(H, inf);
% delta = max(1e-12*s, 100*eps(s));  % positive safety relative to scale
% absH  = abs(H);
% rad   = sum(absH,2) - abs(diag(H));
% need  = rad - diag(H) + delta;
% tau   = max(0, max(need));
% if tau > 0
%     H = H + tau*eye(d);
% end



% Parameters
diag_floor  = 1e-12;   % desired minimal diagonal
eps_shift   = 1e-12;   % tiny extra to escape zero
growth      = 10;      % geometric growth if needed
retry_maxit = 10;      % max iterations in retry loop
tau_0       = 1;      % initial shift for retry loop

% One-shot diagonal lift if min(diag) < floor.
md = min(diag(H));
if md < diag_floor
    tau = (diag_floor - md) + eps_shift;
    H = H + tau*eye(d);
end

% Try raw Cholesky. If it passes, we're done.
[~, p] = chol(H);
num_chol = num_chol + 1;

if p ~= 0
    % Short retry loop with geometric growth
    H0 = H;
    tau = tau_0;
    for it = 1:retry_maxit
        H = H0 + tau*eye(d);

        [~, p] = chol(H);
        num_chol = num_chol + 1;
        if p == 0
            break;
        end
        % increase diagonal shift
        tau = tau*growth;
    end

    % Last resort: scaled identity (very robust)
    if p ~= 0
        gamma = max(mean(abs(diag(H0))), 1e-6);  % scale heuristic
        H = gamma * eye(d);
    end
end