clear all

global problemID pairs a b nballs A Points nPoints LabelPoints emin

problemID = 1;

% =================================================================

if ( problemID == 1 )

    dim_set  = [10, 20, 50, 200, 500, 1000];          % Dimension of "the Cov Matrix"
    snrset   = [0.05, 0.1, 0.25, 0.5, 1.0, 2.0];      % Signal Strength
    deltaset = [0.1, 0.3, 0.7, 0.9];                  % Sparsity
   
    for dim = dim_set

        rng(2025*dim);
        
        for snr = snrset
            
            for delta = deltaset

                % Set up data

                T = dim;
                samplesize = floor(delta*dim);
                S = randsample(dim, samplesize);
                v = zeros(dim,1);
                v(S) = 1/sqrt(samplesize);
                A = sqrt(snr) * v * (v.');
                B = randn(dim)/sqrt(T);
                for ii = 1: dim
                    B(ii,ii) = randn * 2/sqrt(T);
                end
                A = A+B;
                A = (A+A')/2;      

                % Number of variables
                
                n = dim;

                % Constraints

                m = dim;

                E = [];
                I = [];
                equatn = [];
                lambda = [];

                equatn = false(m,1);
                lambda = zeros(m,1);
        
                E = (equatn==true);
                I = (equatn==false);

                M = spherefactory(dim);
                problem.M = M;

                % Set initial guess

                x = problem.M.rand();

                % Checking derivatives?
                
                checkder = false;
                
                % Scale the problem?
                
                scaling = true;

                if ( checkder )
                    l(1:n) = - Inf;
                    u(1:n) =   Inf;
                    checkd(n,m,x0,l,u)
                end
    
                % Set the feseabilty, optimality, and complementarity tolerances
                
                epsfeas   = 10^(-4);
                epsopt    = 10^(-4);
                epscompl  = 10^(-4);

                x0 = x;
                lambda0 = lambda;

                for Alg = 1:4

                    fprintf('dim     = %i\n',dim)
                    fprintf('snr     = %f\n',snr)
                    fprintf('delta   = %f\n\n',delta)

                    % Call the solver

                    % Alg = 1: RSQP-PHR
                    % Alg = 2: RSQP-l1
                    % Alg = 3: AL-TR
                    % Alg = 4: AL-BFGS

                    if Alg == 1 || Alg == 2

                        if Alg == 1, LStype = 1;end
                        if Alg == 2, LStype = 2;end

                        [x,lambda,f,outiter,csupn,snorm,nlmannorm,time,evalfun,evalgrad,evalhess,num_eig,num_chol,info] = ...
                            RSQP(n,m,x0,equatn,lambda0,scaling,problem,LStype,epsopt,epsfeas,epscompl);
                    end

                    if Alg == 3 || Alg == 4

                        if Alg == 3, InSolver = 1;end
                        if Alg == 4, InSolver = 2;end
                        
                        [x,lambda,f,outiter,csupn,snorm,nlmannorm,time,evalfun,evalgrad,evalhess,info] = ...
                            auglag(n,m,x0,equatn,lambda0,scaling,problem,epsopt,epsfeas,epscompl,InSolver);
                    end
                end
            end
        end
    end
end

% ==================================================================

if ( problemID == 2 )

    balls = [2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100];

    for j = 1:length(balls)
    
        % Set the number of balls to be packed
    
        nballs = balls(j);         
    
        rng(2025*nballs);
    
        % Print solution?
        
        print = false;
    
        % Problem data
            
        a = 2;
        b = 1;
    
        % Number of variables
    
        n = 3 * nballs + 1;
    
        % Constraints
    
        if ( nballs == 1 )
            m = 3 * nballs + 1;
        else
            m = nchoosek(nballs,2) + 3 * nballs + 1;
        end
    
        equatn = [];
        lambda = [];
        E = [];
        I = [];
    
        equatn = false(m,1);
        lambda = zeros(m,1);
        
        E = (equatn==true);
        I = (equatn==false);
    
        sphere = spherefactory(2);
        manifold = productmanifold(struct('uv',powermanifold(sphere, nballs),'s', euclideanfactory(nballs),'r', euclideanfactory(1)));
        problem.M = manifold;
    
        if ( nballs == 1 )
            pairs = 0;
        else
            pairs = nchoosek(1:nballs, 2);
        end
    
        % Checking derivatives?
        
        checkder = false;
        
        % Scale the problem?
        
        scaling = true;
    
        % Set initial guess
    
        x = [];

        for i = 1:nballs
            x0 = 2 * rand(2,1) - 1;
            x0 = x0 / norm(x0);
            x.uv{i}(:,1) = x0;
        end
        x.uv = x.uv';
        x.s = rand(nballs,1);
        x.r = rand;
    
        if ( checkder )
            xa = cell2mat(x.uv);
            xb = [];
            xb(1:nballs,1) = xa(1:2:end,1);
            xb(nballs+1:2*nballs,1) = xa(2:2:end,1);
            xz = [xb; x.s; x.r];
            checkd(n,m,xz,-Inf(n,1),Inf(n,1))
        end               

        % Set the feseabilty, optimality, and complementarity tolerances
                    
        epsfeas   = 10^(-4);
        epsopt    = 10^(-4);
        epscompl  = 10^(-4);
    
        % Call the solver
    
        x0 = x;
        lambda0 = lambda;

        for Alg = 1:4

            fprintf('\n\nnballs = %i \n',nballs)

            % Call the solver

            % Alg = 1: RSQP-PHR
            % Alg = 2: RSQP-l1
            % Alg = 3: AL-TR
            % Alg = 4: AL-BFGS

            if Alg == 1 || Alg == 2

                if Alg == 1, LStype = 1;end
                if Alg == 2, LStype = 2;end

                [x,lambda,f,outiter,csupn,snorm,nlmannorm,time,evalfun,evalgrad,evalhess,num_eig,num_chol,info] = ...
                    RSQP(n,m,x0,equatn,lambda0,scaling,problem,LStype,epsopt,epsfeas,epscompl);
            end

            if Alg == 3 || Alg == 4

                if Alg == 3, InSolver = 1;end
                if Alg == 4, InSolver = 2;end
                
                [x,lambda,f,outiter,csupn,snorm,nlmannorm,time,evalfun,evalgrad,evalhess,info] = ...
                    auglag(n,m,x0,equatn,lambda0,scaling,problem,epsopt,epsfeas,epscompl,InSolver);
            end

            if ( print ) printfigure(x), end;
        end 
    end
end

% ==================================================================

if ( problemID == 3 )

    % fileID = fopen('Problem3.dat','w');
    
    for Ptype = 1:4

        rng(2025*Ptype);
    
        % Number of variables
        
        n = 6;
        
        % Constraints  

        m = 1;
        
        equatn = false(m,1);
        lambda = zeros(m,1);
        
        E = (equatn==true);
        I = (equatn==false);
    
        manifold = productmanifold(struct('A',sympositivedefinitefactory(2),'b', euclideanfactory(2)));
        problem.M = manifold;
        
        % Checking derivatives?
        
        checkder = false;
        
        % Scale the problem?
        
        scaling = false;
    
        % Print solution?

        print = false;
    
        % Set the problem data

        nPoints = 10000;
        Points  = -10 + 20 * rand(2,nPoints);
    
        for i = 1:nPoints
            if (  ( Ptype == 1 && norm(Points(:,i)) <= 7 ) || ...
                  ( Ptype == 2 && norm(Points(:,i),Inf) <= 7 ) || ...
                  ( Ptype == 3 && abs(Points(1,i)) <= 7 && abs(Points(2,i)) <= 3.5 ) || ...
                  ( Ptype == 4 && 2 * Points(1,i) - Points(2,i) <= 7 && - Points(1,i) + 2 * Points(2,i) <= 7 && - Points(1,i) - Points(2,i) <= 7 ) )
                LabelPoints(i) = 1;
            else
                LabelPoints(i) = 0;
            end
        end
    
        % Set initial guess    

        x = problem.M.rand();

        emin = 0.95;
  
        if ( checkder )
            xz = [x.A(:);x.b(:)];
            l(1:n) = - Inf;
            u(1:n) =   Inf;
            checkd(n,m,xz,l,u)
        end

        % Set the feseabilty, optimality, and complementarity tolerances
                    
        epsfeas   = 10^(-4);
        epsopt    = 10^(-4);
        epscompl  = 10^(-4);

        % Call the solver

        x0 = x;
        lambda0 = lambda;


        for Alg = 1:4

            fprintf('\n\nPtype   = %i \n',Ptype)

            % Call the solver

            % Alg = 1: RSQP-PHR
            % Alg = 2: RSQP-l1
            % Alg = 3: AL-TR
            % Alg = 4: AL-BFGS

            if Alg == 1 || Alg == 2

                if Alg == 1, LStype = 1;end
                if Alg == 2, LStype = 2;end

                [x,lambda,f,outiter,csupn,snorm,nlmannorm,time,evalfun,evalgrad,evalhess,num_eig,num_chol,info] = ...
                    RSQP(n,m,x0,equatn,lambda0,scaling,problem,LStype,epsopt,epsfeas,epscompl);
            end

            if Alg == 3 || Alg == 4

                if Alg == 3, InSolver = 1;end
                if Alg == 4, InSolver = 2;end
                
                [x,lambda,f,outiter,csupn,snorm,nlmannorm,time,evalfun,evalgrad,evalhess,info] = ...
                    auglag(n,m,x0,equatn,lambda0,scaling,problem,epsopt,epsfeas,epscompl,InSolver);
            end

            if ( print ) printfigure(x); end
        end 
 
    end
end