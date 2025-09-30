function [x,lambda,f,outiter,csupn,snorm,nlmannorm,time,evalfun,evalgrad,evalhess,num_eig,num_chol,sqpinfo] = RSQP( ...
    n, m, x, equatn, lambda, scaling, problem, LStype, epsopt, epsfeas, epscompl)

% RSQP: Riemannian SQP with reduced QP on the tangent space.
% LStype = 1: PHR merit function  
% LStype = 2: l1  merit function  

% Global counters

global evalfun evalgrad evalhess num_eig num_chol

% Parameters
maxoutit  = 100;
maxinit   = 1000;
eps       = 0.5;
LSparam   = 1e-4;
stpmin    = 1e-15;

E = (equatn == true);
I = ~E;

% Scaling factors (use Riemannian norms)
if scaling
    % Objective
    egrad_f_vec = evalg(n, x);                  % Euclidean gradient (vectorized)
    egrad_f     = reshapevector(egrad_f_vec);   % -> ambient struct (shape of M)
    rgrad_f     = problem.M.egrad2rgrad(x, egrad_f);
    sc.f        = 1 / max(1, problem.M.norm(x, rgrad_f));

    % Constraints
    sc.c = ones(1, m);
    for i = 1:m
        egrad_ci_vec = evalnc(n, x, i);           % Euclidean grad of c_i (vectorized)
        egrad_ci     = reshapevector(egrad_ci_vec);
        rgrad_ci     = problem.M.egrad2rgrad(x, egrad_ci);
        sc.c(i)      = 1 / max(1, problem.M.norm(x, rgrad_ci));
    end
else
    sc.f = 1;
    sc.c = ones(1, m);
end

% Print header
fprintf('----------------------------------------------------------------------\n')
fprintf('   RSQP algorithm for nonlinear optimization on Riemannian manifolds  \n')
fprintf('----------------------------------------------------------------------\n')
fprintf('Number of variables                : %i \n', n)
fprintf('Number of constraints              : %i \n\n', m)
fprintf('Optimality tolerance               : %.0e \n', epsopt)
fprintf('Feasibility tolerance              : %.0e \n', epsfeas)
fprintf('Complementarity tolerance          : %.0e \n\n', epscompl)
if LStype == 1
    fprintf('Line search merit function         : PHR\n\n', epscompl)
elseif LStype == 2
    fprintf('Line search merit function         : l1\n\n', epscompl)
end
if scaling
    fprintf('Objective function scale factor    : %.0e \n', sc.f)
    fprintf('Smallest constraints scale factor  : %.0e \n', min(sc.c))
end

% Start timing
tic;

% Initialization of counters and flags
outiter  = 0;
evalfun  = 0;
evalgrad = 0;
evalhess = 0;
num_eig  = 0;
num_chol = 0;
ISerror  = 0;

% Initial objective and gradient (scaled)
[f,fs,~] = sevalf(n,x,sc,scaling);    evalfun  = evalfun  + 1;
[~,gs,~] = sevalg(n,x,sc,scaling);    evalgrad = evalgrad + 1;

% Constraints and Jacobian (scaled)
c   = zeros(m,1);
cs  = zeros(m,1);
Jcs = zeros(m,n);
for i = 1:m
    [c(i), cs(i), ~] = sevalc(n,x,i,sc,scaling);    evalfun  = evalfun  + 1;
    [~, Jcs(i,:), ~] = sevalnc(n,x,i,sc,scaling);   evalgrad = evalgrad + 1;
end

% Violations and norms
csupn  = max( norm(abs(c(E)),Inf), norm(max(c(I),0),Inf) );
cssupn = max( norm(abs(cs(E)),Inf), norm(max(cs(I),0),Inf) );
snorm  = norm(min(-cs(I), lambda(I)), Inf);

% Euclidean Lagrangian gradient (vectorized), then to manifold struct
nl_vec = gs + Jcs' * lambda;
nl     = reshapevector(nl_vec);                 % ambient struct
nlman  = problem.M.egrad2rgrad(x, nl);          % tangent struct
nlmannorm = problem.M.norm(x, nlman);

% ==================================================================
% Main loop
% ==================================================================

while true

    % ==================================================================
    % Print information of this iteration
    % ==================================================================

    if scaling
        if mod(outiter,10) == 0
            fprintf('\n')
            fprintf('%3s   %-9s %-6s  %-6s  %+9s %-6s %-4s    %-5s  %-5s  %-5s \n',...
                'out','objective','infeas','scaled','scaled','comple','norm','norm','line','inner')
            fprintf('%3s   %-9s %-6s  %-6s  %-6s %-6s %-5s  %-5s %-6s %-5s\n',...
                'ite','function ','ibilty','obj-funct','infeas','mentar','graLag','direct','search','stop')
        end

        if ( outiter == 0 )
            fprintf('%3i  %10.3e %6.0e %10.3e  %6.0e %6.0e %6.0e    -      -      - \n',...
                outiter, f, csupn, fs, cssupn, snorm, nlmannorm)
        else
            fprintf('%3i  %10.3e %6.0e %10.3e  %6.0e %6.0e %6.0e  %6.0e   %i     %2i\n',...
                outiter, f, csupn, fs, cssupn, snorm, nlmannorm, dsupn, flagLS, exitflag)
        end
    else
        if mod(outiter,10) == 0
            fprintf('\n')
            fprintf('%3s   %-9s %-6s  %-6s %-6s %-4s   %-5s   %-5s  %-5s\n',...
                'out','objective','infeas','infeas','comple','norm','norm','line','inner')
            fprintf('%3s   %-9s %-6s  %-6s %-6s %-5s %-6s  %-6s %-5s\n',...
                'ite','function ','ibilty','+compl','mentar','graLag','direct','search','stop ')
        end

        if ( outiter == 0 )
            fprintf('%3i  %10.3e %6.0e  %6.0e %6.0e %6.0e    -      -      - \n',...
                outiter, f, csupn, cssupn, snorm, nlmannorm)
        else
            fprintf('%3i  %10.3e %6.0e  %6.0e %6.0e %6.0e %6.0e    %i     %2i\n',...
                outiter, f, csupn, cssupn, snorm, nlmannorm, dsupn, flagLS, exitflag)
        end
    end

    % ==================================================================
    % Test stopping criteria
    % ==================================================================

    if (csupn <= epsfeas) && (snorm <= epscompl) && (nlmannorm <= epsopt)
        sqpinfo = 1;
        time = toc;
        fprintf('\nSolution was found.\n')
        fprintf('Number of function evaluations   : %i\n',evalfun)
        fprintf('Number of gradient evaluations   : %i\n',evalgrad)
        fprintf('Number of Hessian evaluations    : %i\n',evalhess)
        fprintf('CPU time(s)                      : %.1f \n',time)
        return
    end

    if outiter >= maxoutit
        sqpinfo = 2;
        time = toc;
        fprintf('\nMaximum of iterations reached.\n')
        fprintf('Number of function evaluations   : %i\n',evalfun)
        fprintf('Number of gradient evaluations   : %i\n',evalgrad)
        fprintf('Number of Hessian evaluations    : %i\n',evalhess)
        fprintf('CPU time(s)                      : %.1f \n',time)
        return
    end

    if ISerror >= 3
        sqpinfo = 3;
        time = toc;
        fprintf('\nThe subproblem cannot be solved to sufficient accuracy.\n')
        fprintf('Number of function evaluations   : %i\n',evalfun)
        fprintf('Number of gradient evaluations   : %i\n',evalgrad)
        fprintf('Number of Hessian evaluations    : %i\n',evalhess)
        fprintf('CPU time(s)                      : %.1f \n',time)
        return
    end

    % ==================================================================
    % Iteration
    % ==================================================================

    outiter = outiter + 1;

    % Build Euclidean Hessian of the Lagrangian
    [~, H_obj, flag] = sevalh(n,x,sc,scaling);         
    if flag == 1 
        HL = H_obj;
        evalhess = evalhess + 1;
    else
        HL = zeros(n,n);
    end
    
    for i_h = 1:m
        [~, Hc_i, flag] = sevalhc(n,x,i_h,sc,scaling); 
        if flag == 1
            HL = HL + lambda(i_h) * Hc_i;
            evalhess = evalhess + 1;
        end
    end
    HL = (HL + HL.')/2;  % symmetrize numerically

    % ==================================================================
    % Subproblem
    % ==================================================================

    % Reduced QP in tangent coordinates
    [Htilde, gtilde, Aeq, beq, A, b, lb, ub, d_from_y, mult_recover] = ...
        build_subproblem(problem.M, x, gs, nl_vec, HL, cs, Jcs, equatn);

    % Regularize Htilde to SPD
    Htilde = make_SPD(Htilde);

    % Solve reduced QP
    optsQP = optimoptions('quadprog','Display','off','Algorithm','interior-point-convex','MaxIterations',maxinit);
    [y,~,exitflag,~,Lagmult] = quadprog(Htilde, gtilde, A, b, Aeq, beq, lb, ub, zeros(numel(gtilde),1), optsQP);

    if exitflag ~= 1
        ISerror = ISerror + 1;
    else
        ISerror = 0;
    end

    % Tangent step and its norm
    d = d_from_y(y);               % tangent struct
    dsupn = problem.M.norm(x, d);  % Riemannian norm    

    % Update the Lagrange multipliers
    if exitflag == 0 || exitflag == 1 || exitflag == 2
        lambda = mult_recover(Lagmult);
    else
        lambda(:) = 0;
    end

    % ==================================================================
    % Line search
    % ==================================================================

    % Line search with exponential map
    if outiter == 1
        rho = comprhoini(cs,fs,E,I,LStype);
    end

    % Update penalty parameter
    v = max( max(abs(Lagmult.eqlin)), max(Lagmult.ineqlin) );
    if rho < v
        rho = v + eps;
    end

    % Merit at current iterate
    merit = fs;
    if LStype == 1

        % PHR merit function  
        for i = 1:m
            if equatn(i) || (lambda(i) + rho*cs(i) > 0)
                merit = merit + cs(i) * (lambda(i) + 0.5*rho*cs(i));
            else
                merit = merit - 0.5 * (lambda(i)^2) / rho;
            end
        end

    elseif LStype == 2

        % l1 merit function
        for i = 1:m
            if equatn(i)
                merit = merit + rho * abs(cs(i));
            else
                merit = merit + rho * max(0,cs(i));
            end
        end
    end

    % Directional model in the tangent space:
    % ⟨H d, d⟩_R = y' * Htilde * y   (because the tangent basis is orthonormal)
    ftol = LSparam * dot( Htilde * y, y); % β * ⟨H d, d⟩_R  (>= 0)

    % Backtracking with exponential map
    stp = 1;
    while true
        % Trial point via the manifold's exponential
        xtrial = problem.M.exp(x, d, stp);

        % Merit at trial
        merit_trial = eval_merit(xtrial, n, m, equatn, lambda, rho, sc, scaling, LStype);

        % Armijo-type test:
        if merit_trial <= merit - ftol * stp
            flagLS = 1;    % success
            break;
        end
        stp = 0.5 * stp;
        if stp <= stpmin
            flagLS = 2;    % failed backtracking
            break;
        end
    end

    % Accept step
    x = xtrial;

    % ==================================================================
    % Prepare next iteration
    % ==================================================================

    % Objective and grad (scaled)
    [f,fs,~] = sevalf(n,x,sc,scaling);    evalfun  = evalfun  + 1;
    [~,gs,~] = sevalg(n,x,sc,scaling);    evalgrad = evalgrad + 1;

    % Constraints and Jacobian (scaled)
    for i = 1:m
        [c(i), cs(i), ~] = sevalc(n,x,i,sc,scaling);    evalfun  = evalfun  + 1;
        [~, Jcs(i,:), ~] = sevalnc(n,x,i,sc,scaling);   evalgrad = evalgrad + 1;
    end

    % Violations
    csupn  = max( norm(abs(c(E)),Inf), norm(max(c(I),0),Inf) );
    cssupn = max( norm(abs(cs(E)),Inf), norm(max(cs(I),0),Inf) );
    snorm  = norm(min(-cs(I), lambda(I)), Inf);

    % Lagrangian grad (Euclidean), then Riemannian and its norm
    nl_vec = gs + Jcs.' * lambda;
    nl     = reshapevector(nl_vec);
    nlman  = problem.M.egrad2rgrad(x, nl);
    nlmannorm = problem.M.norm(x, nlman);

    % ==================================================================
    % Iterate
    % ==================================================================
end

% ==================================================================
% End of main loop
% ==================================================================
end