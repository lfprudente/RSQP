function [Htilde, gtilde, Aeq, beq, A, b, lb, ub, d_from_y, mult_recover] = ...
  build_subproblem(M, x, gf_vec, gl_vec, HL, c_vals, Jc, equatn)

% BUILD_SUBPROBLEM  Reduced QP data (and map y->d)
%
% Inputs:
%   M         : Manopt manifold object
%   x         : point on the manifold (struct or array, per problem)
%   gf_vec    : Euclidean gradient of the objective, vectorized (n×1)
%   gl_vec    : Euclidean gradient of the Lagrangian (needed for curvature)
%   HL        : Euclidean Hessian of the Lagrangian (n×n), symmetric
%   c_vals    : m×1 constraint values c(x)
%   Jc        : m×n Euclidean Jacobian; row i is (∇c_i)^T
%   equatn    : m×1 logical, true = equality, false = inequality
%
% Outputs:
%   Htilde, gtilde, Aeq, beq, A, b : reduced QP data in tangent coords
%   d_from_y    : function handle mapping y (reduced coords) -> d (tangent)
%   mult_recover: function handle used to return the Lagrange multipliers 
%                 from 'quadprog' in the ordering expected by 'evalcc'


global problemID

if problemID == 1
    % ---------- Problem 1: Sphere S^{n-1} (fast Householder ) ----------

    n = numel(x); 
    d = n - 1;

    % Householder vector v with H*e1 = x, H = I - 2 v v^T
    e1 = zeros(n,1); 
    e1(1) = 1;
    if norm(x - e1) < 1e-14
        v = zeros(n,1); 
        v(1) = 1;     % H = I
    else
        v = (e1 - x); 
        v = v / norm(v);
    end

    % Linear term (OBJECTIVE): g̃ = (H^T ∇f)_{2:n}
    w_f   = gf_vec(:) - 2*v*(v.'*gf_vec(:));
    gtilde = w_f(2:end);

    % Quadratic term (LAGRANGIAN): H̃ = B22 - α I, α = x' ∇L
    alpha = x.' * gl_vec(:);
    HL    = (HL + HL.')/2;
    t1    = HL * v;                 % n×1
    s     = v.' * t1;               % scalar
    v2    = v(2:end);
    t12   = t1(2:end);
    HL22  = HL(2:end, 2:end);
    B22   = HL22 - 2*(t12 * v2.') - 2*(v2 * t12.') + 4*s*(v2 * v2.');
    Htilde = (B22 - alpha * eye(d));
    Htilde = (Htilde + Htilde.')/2;

    % Constraints: rows = (H^T ∇c_i)_{2:n}^T, rhs = -c_i(x)
    idxI = find(~equatn(:));

    A = zeros(numel(idxI), d); 
    b = zeros(numel(idxI),1);
    for r = 1:numel(idxI)
        i   = idxI(r);
        gi  = Jc(i,:).';
        wi  = gi - 2*v*(v.'*gi);
        A(r,:) = wi(2:end).';
        b(r)   = - c_vals(i);
    end

    % No equalities in Problem 1
    Aeq = []; beq = [];

    % No bounds in Problem 1
    lb = []; ub = [];

    % This returns lambda (m×1) in the SAME ordering as evalcc uses.
    mult_recover = @(Lagmult) recover_multipliers_p1(Lagmult);

    % Map back y -> ambient tangent direction: d = H * [0; y]
    d_from_y = @(y) apply_householder_to_z(v, y);

    return
end

% ==================================================================

if problemID == 2

     % ---------- Problem 2: (S^1)^N × R^N × R (fast, block-wise) ----------
    
    if ~isstruct(x)
        x = reshapevector(x);
    end

    % Dimensions and fixed index layout
    N     = numel(x.uv);                     % number of circles
    m     = numel(c_vals);                   % total inequalities
    nck   = m - (3*N + 1);                   % number of pairwise constraints
    idx_u = 1:N;                             % u_1,...,u_N
    idx_v = N+(1:N);                         % v_1,...,v_N
    idx_s = 2*N+(1:N);                       % s_1,...,s_N
    idx_r = 3*N+1;                           % r
    d     = 2*N + 1;                         % reduced dim: [theta; y_s; y_r]

    % Unit tangents T_i = J * uv_i
    J  = [0 -1; 1 0];
    T  = zeros(2,N); 
    UV = zeros(2,N);
    for i = 1:N
        ui = x.uv{i}(:); 
        nui = norm(ui); 
        if nui==0 
            ui=[1;0]; 
            nui=1; 
        end
        ui = ui / nui; 
        UV(:,i) = ui; 
        T(:,i) = J*ui;
    end
    t1 = T(1,:); 
    t2 = T(2,:);

    % Reduced gradient (OBJECTIVE)
    Guv_f   = [gf_vec(idx_u).'; gf_vec(idx_v).'];
    g_theta = sum(T .* Guv_f, 1).';
    g_s     = gf_vec(idx_s);
    g_r     = gf_vec(idx_r);
    gtilde  = [g_theta; g_s; g_r];

    % Reduced Hessian (LAGRANGIAN)
    Huu = HL(idx_u, idx_u); 
    Huv = HL(idx_u, idx_v);
    Hvu = HL(idx_v, idx_u); 
    Hvv = HL(idx_v, idx_v);

    Hus = HL(idx_u, idx_s); 
    Hvs = HL(idx_v, idx_s);
    Hur = HL(idx_u, idx_r); 
    Hvr = HL(idx_v, idx_r);

    Hss = HL(idx_s, idx_s); 
    Hsr = HL(idx_s, idx_r); 
    Hrr = HL(idx_r, idx_r);

    D1 = diag(t1); 
    D2 = diag(t2);
    B  = D1*Huu*D1 + D1*Huv*D2 + D2*Hvu*D1 + D2*Hvv*D2;

    Guv_L = [gl_vec(idx_u).'; gl_vec(idx_v).'];
    alpha = sum(UV .* Guv_L, 1).';
    Htheta_theta = (B + B.')/2 - diag(alpha);

    Htheta_s = D1*Hus + D2*Hvs;
    Htheta_r = t1(:).*Hur + t2(:).*Hvr;

    Htilde = zeros(d,d);
    Htilde(1:N,1:N)         = Htheta_theta;
    Htilde(1:N,N+(1:N))     = Htheta_s;
    Htilde(1:N,2*N+1)       = Htheta_r;
    Htilde(N+(1:N),1:N)     = Htheta_s.';
    Htilde(N+(1:N),N+(1:N)) = (Hss + Hss.')/2;
    Htilde(N+(1:N),2*N+1)   = Hsr;
    Htilde(2*N+1,1:N)       = Htheta_r.';
    Htilde(2*N+1,N+(1:N))   = Hsr.';
    Htilde(2*N+1,2*N+1)     = Hrr;
    Htilde = (Htilde + Htilde.')/2;

    % --------------------- Box bounds in reduced coords ---------------------
    % From your evalcc: s_i ∈ [0,1], r ≥ 0  →  in y: 0 - s ≤ y_s ≤ 1 - s and 0 - r ≤ y_r
    lb = -inf(d,1);  ub = +inf(d,1);
    lb(N+(1:N)) = 0 - x.s(:);
    ub(N+(1:N)) = 1 - x.s(:);
    lb(2*N+1)   = 0 - x.r;
    % ub(2*N+1) stays +Inf unless you add an upper bound on r

    % --------------------- Linearized NON-box constraints -------------------
    % Using fixed index structure:
    %   Block 1: ind = 1 : N                       (nonlinear, keep)
    %   Block 2: ind = N+1 : N+nck                 (pairwise, keep)
    %   Block 3: ind = N+nck + (1:N)               (s lower)   -> moved to bounds
    %   Block 4: ind = 2N+nck + (1:N)              (s upper)   -> moved to bounds
    %   Block 5: ind = 3N+nck + 1                  (r lower)   -> moved to bounds

    keep_idx = [ (1:N), (N+1 : N+nck) ];
    A = zeros(numel(keep_idx), d);
    b = zeros(numel(keep_idx), 1);

    row = 0;
    for ii = keep_idx
        row = row + 1;

        % Euclidean gradient pieces for constraint ii
        gcU = Jc(ii, idx_u).';                    % N×1 (u-part)
        gcV = Jc(ii, idx_v).';                    % N×1 (v-part)
        row_theta = (t1 .* gcU.' + t2 .* gcV.');  % 1×N: projection onto theta
        row_s     = Jc(ii, idx_s);                % 1×N
        row_r     = Jc(ii, idx_r);                % scalar

        A(row,:)  = [row_theta, row_s, row_r];
        b(row)    = - c_vals(ii);
    end

    % No equalities in Problem 2
    Aeq = []; beq = [];

     % This returns lambda (m×1) in the SAME ordering as evalcc uses.
    mult_recover = @(Lagmult) recover_multipliers_p2(Lagmult, N, nck);
    
    % ---- Map back y -> tangent struct on product manifold ----
    d_from_y = @(y) y_to_tangent_struct_prod(M, x, T, y);

    return
end

% ==================================================================

if problemID == 3

    if ~isstruct(x)
        x = reshapevector(x);
    end

    % Dimensions
    d = M.dim();               % should be 5 (3 on SPD(2) + 2 on R^2)
    n = numel(gf_vec);         % 6 variables: [A11 A21 A12 A22 b1 b2]
    HL = (HL + HL.')/2;        % symmetrize for hygiene

    % === 1) Build an orthonormal basis Q of T_xM ===
    Q = make_onb_spd2_r2_fast(x);

    % === 2) Reduced gradient (objective): use Riemannian grad of f ===
    egrad_f = reshapevector(gf_vec);      % ambient struct (A,b)
    rgrad_f = M.egrad2rgrad(x, egrad_f);  % tangent struct
    gtilde  = zeros(d,1);
    for k = 1:d
        gtilde(k) = M.inner(x, Q{k}, rgrad_f);
    end

    % === 3) Reduced Hessian (Lagrangian): use Riemannian Hessian ===
    egrad_L = reshapevector(gl_vec);
    ehessL  = @(eta) reshapevector( HL * reshapevector(eta) );
    Htilde  = zeros(d,d);
    for j = 1:d
        qj  = Q{j};
        rhj = M.ehess2rhess(x, egrad_L, ehessL(qj), qj);  % tangent struct
        for i = 1:d
            Htilde(i,j) = M.inner(x, Q{i}, rhj);
        end
    end
    Htilde = (Htilde + Htilde.')/2;

    % === 4) Constraints: all inequalities A*y <= b (no equalities here) ===
    m   = size(Jc,1);            % m = 4
    idxE = find(equatn(:));      % empty in this problem
    idxI = find(~equatn(:));     % 1:4

    Aeq = zeros(numel(idxE), d); beq = zeros(numel(idxE), 1); % empty
    A   = zeros(numel(idxI), d); b   = zeros(numel(idxI), 1);

    for rr = 1:numel(idxI)
        i        = idxI(rr);
        egrad_ci = reshapevector(Jc(i,:).');     % ambient struct
        rgrad_ci = M.egrad2rgrad(x, egrad_ci);   % tangent struct
        for k = 1:d
            A(rr,k) = M.inner(x, Q{k}, rgrad_ci);
        end
        b(rr) = - c_vals(i);
    end

    % No bounds in reduced coords
    lb = []; ub = [];

    % === 5) Map back: y -> tangent direction d ===
    d_from_y = @(y) lincomb_many(M, x, Q, y);

    % === 6) Recover lambda: all constraints were sent as inequalities ===
    mult_recover = @(Lagmult) recover_multipliers_p3(Lagmult);
end

end

% ==================================================================
% Helpers (private)
% ==================================================================

% ==================================================================
% Problem 1
% ==================================================================

% ---------- Apply Householder to [0; y]  ----------
function d = apply_householder_to_z(v, y)
    % Apply d = H*[0;y] with H = I - 2 v v^T, without forming H explicitly.
    z = [0; y(:)];
    d = z - 2*v*(v.'*z);
end

function lambda = recover_multipliers_p1(Lagmult)
    lambda(:,1) = max(0, Lagmult.ineqlin(:));
end

% ==================================================================
% Problem 2
% ==================================================================

% ---------- Map y -> tangent struct  ----------
function d = y_to_tangent_struct_prod(M, x, T, y)
    % y = [theta(1..N), s(1..N), r]
    if ~isstruct(x)
        x = reshapevector(x); 
    end
    N = numel(x.uv);
    d = M.zerovec(x);
    y = y(:);
    y_theta = y(1:N);
    for i = 1:N
        d.uv{i} = T(:,i) * y_theta(i);
    end
    d.s = y(N+(1:N));
    d.r = y(2*N+1);
end

function lambda = recover_multipliers_p2(Lagmult, N_loc, nck_loc)
    % Sizes
    d_loc = 2*N_loc + 1;             % length of y
    m_loc = nck_loc + 3*N_loc + 1;   % total inequalities in evalcc
    lambda = zeros(m_loc,1);

    % 1&2) rows that remained as A*y<=b:
    lamI = Lagmult.ineqlin(:);
    if numel(lamI) ~= (N_loc + nck_loc)
        error('mult_recover(P2): ineqlin length %d, expected %d.', ...
               numel(lamI), N_loc + nck_loc);
    end
    lambda(1:N_loc)           = lamI(1:N_loc);             % block 1
    lambda(N_loc+1:N_loc+nck_loc) = lamI(N_loc+1:end);     % block 2

    % bounds: lower/upper may be empty depending on solver activity
    lamLower = zeros(d_loc,1);
    lamUpper = zeros(d_loc,1);
    if isfield(Lagmult,'lower') && ~isempty(Lagmult.lower)
        lamLower = Lagmult.lower(:);
    end
    if isfield(Lagmult,'upper') && ~isempty(Lagmult.upper)
        lamUpper = Lagmult.upper(:);
    end

    % 3) s-lower → lower bounds at reduced positions N+(1:N)
    lambda(N_loc+nck_loc + (1:N_loc))   = lamLower(N_loc+(1:N_loc));
    % 4) s-upper → upper bounds at reduced positions N+(1:N)
    lambda(2*N_loc+nck_loc + (1:N_loc)) = lamUpper(N_loc+(1:N_loc));
    % 5) r-lower → lower bound at reduced position 2N+1
    lambda(3*N_loc + nck_loc + 1)       = lamLower(2*N_loc+1);
end

% ==================================================================
% Problem 3
% ==================================================================

function Q = make_onb_spd2_r2_fast(x)
% Build an explicit orthonormal basis at (A,b) for SPD(2) x R^2.
% SPD(2) part: Q_i^{(A)} = A^{1/2} * \hat Q_i * A^{1/2}, where
%   { \hat Q_i } is Frobenius-orthonormal in Sym(2):
%     hatQ1 = [1 0; 0 0]
%     hatQ2 = [0 1; 1 0] / sqrt(2)
%     hatQ3 = [0 0; 0 1]
% These yield an orthonormal basis under the affine-invariant metric.

    % Compute A^{1/2} robustly
    A = (x.A + x.A.')/2;
    [U,D] = eig(A);
    d = max(diag(D), 1e-14);
    Ah = U*diag(sqrt(d))*U.';      % A^{1/2}

    hatQ1 = [1 0; 0 0];
    hatQ2 = [0 1; 1 0] / sqrt(2);
    hatQ3 = [0 0; 0 1];

    Q = cell(5,1);
    % SPD(2) block (3 dirs)
    q = struct('A', [], 'b', []);
    q.A = Ah * hatQ1 * Ah; q.b = [0;0]; Q{1} = q;
    q.A = Ah * hatQ2 * Ah; q.b = [0;0]; Q{2} = q;
    q.A = Ah * hatQ3 * Ah; q.b = [0;0]; Q{3} = q;
    % Euclidean R^2 block (2 dirs)
    q.A = zeros(2); q.b = [1;0]; Q{4} = q;
    q.A = zeros(2); q.b = [0;1]; Q{5} = q;
end

function d = lincomb_many(M, x, Q, y)
% Combine cell ONB Q with weights y into a tangent direction.
    d = M.zerovec(x);
    for i = 1:numel(Q)
        if y(i) ~= 0
            d = M.lincomb(x, 1.0, d, y(i), Q{i});
        end
    end
end

function lambda = recover_multipliers_p3(Lagmult)
    lambda(:,1) = max(0, Lagmult.ineqlin(:));
end