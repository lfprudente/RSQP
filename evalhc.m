function [hc,flag] = evalhc(n,x,ind)

global problemID pairs nballs a b N rankY

% ==================================================================

if ( problemID == 1 )

    flag = 0;

    hc = zeros(n);

    return
end

% ==================================================================

if ( problemID == 2 )

    flag = 1;

    if (~isstruct(x))
        x = reshapevector(x);
    end

    % Conveniences
    if (nballs == 1)
        nck = 0;
    else
        nck = nchoosek(nballs, 2);
    end
    cte = (b / a)^2;

    % Indices base
    % u_i -> i
    % v_i -> N + i
    % s_i -> 2N + i
    % r   -> 3N + 1
    N  = nballs;
    Pr = 3*N + 1;

    % Inicialize the Hessian
    hc = zeros(n);

    % --------- Case 1: 1 <= ind <= N ------------------------------------
    % c_i = r^2 - b^2 (s_i - 1)^2 ( cte*u_i^2 + v_i^2 )
    if (ind >= 1 && ind <= N)

        i  = ind;
        Pui = i;
        Pvi = N + i;
        Psi = 2*N + i;

        ui = x.uv{i}(1);
        vi = x.uv{i}(2);
        si = x.s(i);
        ri = x.r;

        % Second derivatives (all other second derivatives are zero):
        % ∂²c/∂u_i²
        hc(Pui,Pui) = -2 * b^2 * cte * (si - 1)^2;

        % ∂²c/∂v_i²
        hc(Pvi,Pvi) = -2 * b^2 * (si - 1)^2;

        % ∂²c/∂s_i²
        hc(Psi,Psi) = -2 * b^2 * ( cte*ui^2 + vi^2 );

        % ∂²c/∂r²
        hc(Pr,Pr) = 2;

        % Mixed:
        % ∂²c/∂u_i∂s_i = ∂²c/∂s_i∂u_i
        val = -4 * b^2 * cte * ui * (si - 1);
        hc(Pui,Psi) = val;
        hc(Psi,Pui) = val;

        % ∂²c/∂v_i∂s_i = ∂²c/∂s_i∂v_i
        val = -4 * b^2 * vi * (si - 1);
        hc(Pvi,Psi) = val;
        hc(Psi,Pvi) = val;

        % (crossed with r or between u_i and v_i are 0)

    % --------- Case 2: N+1 <= ind <= N + nck -----------------------------
    % c_ij = 4 r^2
    %        - a^2 * ( alpha_i*u_i - alpha_j*u_j )^2
    %        - b^2 * ( s_i*v_i - s_j*v_j )^2
    % where alpha_i = 1 + cte*(s_i - 1)
    elseif (ind >= N + 1 && ind <= N + nck)

        i = pairs(ind - N, 1);
        j = pairs(ind - N, 2);

        Pui = i;      Pvi = N + i;    Psi = 2*N + i;
        Puj = j;      Pvj = N + j;    Psj = 2*N + j;

        ui = x.uv{i}(1);  vi = x.uv{i}(2);  si = x.s(i);
        uj = x.uv{j}(1);  vj = x.uv{j}(2);  sj = x.s(j);
        r  = x.r;

        alpha_i = 1 + cte*(si - 1);
        alpha_j = 1 + cte*(sj - 1);

        X = alpha_i*ui - alpha_j*uj;   % termo "x"
        Y = si*vi - sj*vj;             % termo "y"

        % ----- block (u,u) -----
        % ∂²c/∂u_i² = -2 a^2 alpha_i^2
        hc(Pui,Pui) = -2 * a^2 * alpha_i^2;

        % ∂²c/∂u_j² = -2 a^2 alpha_j^2
        hc(Puj,Puj) = -2 * a^2 * alpha_j^2;

        % ∂²c/∂u_i∂u_j = ∂²c/∂u_j∂u_i = +2 a^2 alpha_i alpha_j
        val =  2 * a^2 * alpha_i * alpha_j;
        hc(Pui,Puj) = val; hc(Puj,Pui) = val;

        % ----- block (v,v) -----
        % ∂²c/∂v_i² = -2 b^2 s_i^2
        hc(Pvi,Pvi) = -2 * b^2 * si^2;

        % ∂²c/∂v_j² = -2 b^2 s_j^2
        hc(Pvj,Pvj) = -2 * b^2 * sj^2;

        % ∂²c/∂v_i∂v_j = ∂²c/∂v_j∂v_i = +2 b^2 s_i s_j
        val =  2 * b^2 * si * sj;
        hc(Pvi,Pvj) = val; hc(Pvj,Pvi) = val;

        % ----- block (s,s) -----
        % ∂²c/∂s_i² = -2 a^2 cte^2 u_i^2 - 2 b^2 v_i^2
        hc(Psi,Psi) = -2 * a^2 * cte^2 * ui^2 - 2 * b^2 * vi^2;

        % ∂²c/∂s_j² = -2 a^2 cte^2 u_j^2 - 2 b^2 v_j^2
        hc(Psj,Psj) = -2 * a^2 * cte^2 * uj^2 - 2 * b^2 * vj^2;

        % ∂²c/∂s_i∂s_j = ∂²c/∂s_j∂s_i = +2 a^2 cte^2 u_i u_j + 2 b^2 v_i v_j
        val =  2 * a^2 * cte^2 * ui * uj + 2 * b^2 * vi * vj;
        hc(Psi,Psj) = val; hc(Psj,Psi) = val;

        % ----- mixed u-s -----
        % ∂²c/∂u_i∂s_i = -2 a^2 cte (alpha_i*u_i + X)
        val = -2 * a^2 * cte * (alpha_i*ui + X);
        hc(Pui,Psi) = val; hc(Psi,Pui) = val;

        % ∂²c/∂u_j∂s_j =  2 a^2 cte (-alpha_j*u_j + X)
        val =  2 * a^2 * cte * (-alpha_j*uj + X);
        hc(Puj,Psj) = val; hc(Psj,Puj) = val;

        % ∂²c/∂u_i∂s_j = +2 a^2 cte * alpha_i * u_j
        val =  2 * a^2 * cte * alpha_i * uj;
        hc(Pui,Psj) = val; hc(Psj,Pui) = val;

        % ∂²c/∂u_j∂s_i = +2 a^2 cte * alpha_j * u_i
        val =  2 * a^2 * cte * alpha_j * ui;
        hc(Puj,Psi) = val; hc(Psi,Puj) = val;

        % ----- mixed v-s -----
        % ∂²c/∂v_i∂s_i = -2 b^2 (vi*si + Y)
        val = -2 * b^2 * (vi*si + Y);
        hc(Pvi,Psi) = val; hc(Psi,Pvi) = val;

        % ∂²c/∂v_j∂s_j =  2 b^2 (-vj*sj + Y)
        val =  2 * b^2 * (-vj*sj + Y);
        hc(Pvj,Psj) = val; hc(Psj,Pvj) = val;

        % ∂²c/∂v_i∂s_j = +2 b^2 * si * v_j
        val =  2 * b^2 * si * vj;
        hc(Pvi,Psj) = val; hc(Psj,Pvi) = val;

        % ∂²c/∂v_j∂s_i = +2 b^2 * sj * v_i
        val =  2 * b^2 * sj * vi;
        hc(Pvj,Psi) = val; hc(Psi,Pvj) = val;

        % ----- r -----
        % ∂²c/∂r² = 8
        hc(Pr,Pr) = 8;

        % (all other second derivatives are zero)

    % --------- Caso 3: N + nck + 1  ...  2N + nck  -----------------------
    % c = - s_i  -> Null Hessian
    elseif (ind >= N + nck + 1 && ind <= 2*N + nck)
        % hc is already zeroed
        flag = 0;

    % --------- Caso 4: 2N + nck + 1 ... 3N + nck  ------------------------
    % c = s_i - 1 -> Null Hessian
    elseif (ind >= 2*N + nck + 1 && ind <= 3*N + nck)
        % hc is already zeroed
        flag = 0;

    % --------- Case 5: last (r >= 0): c = - r -> Null Hessian ---------
    else
        % hc is already zeroed
        flag = 0;
    end

    % Symmetrizes for numerical safety
    hc = (hc + hc.')/2;

    return
end

% ==================================================================

if ( problemID == 3 )

    flag = 0;

    if ~isstruct(x)
        x = reshapevector(x);
    end

    hc = zeros(n,n);  

    if ind == 1

        A  = (x.A + x.A')/2;
    
        % Derivatives w.r.t. (s,p) (phi = e^2):
        s  = trace(A);
        s12= 0.5*(A(1,2)+A(2,1));
        p  = A(1,1)*A(2,2) - s12^2;
    
        Delta = sqrt(max(s*s - 4*p, 0));
        h     = s + Delta;
    
        % phi e derivadas de 1a ordem
        phi   = 2*Delta / max(h,1e-16);
        f     = 1 / max(Delta*h*h, 1e-16);          % f = 1/(Delta*(s+Delta)^2)
    
        phi_s = 8*p * f;
        phi_p = -4*s * f;
    
        % f_s e f_p
        epsd  = 1e-16;
        Delta2 = max(Delta*Delta, epsd);
        fs = - ( (s/Delta2) + 2/h + (2*s)/(Delta*h) ) * f;
        fp =   ( 2/Delta2 + 4/(Delta*h) ) * f;
    
        % 2o ordem derivatives
        phi_ss = 8*p * fs;           % d/ds (8 p f) = 8 p f_s
        phi_sp = 8*f + 8*p * fp;     % d/dp (8 p f) = 8 f + 8 p f_p
        phi_pp = -4*s * fp;          % d/dp (-4 s f) = -4 s f_p
    
        % ds/dA e dp/dA (simetrized) in order [A11 A21 A12 A22]
        gs = [1; 0; 0; 1];
        s12 = 0.5*(A(1,2)+A(2,1));
        gp = [A(2,2); -s12; -s12; A(1,1)];
    
        % Hessiana of p simetrized
        Hp = [ 0   0    0   1;
               0  -1/2 -1/2 0;
               0  -1/2 -1/2 0;
               1   0    0   0 ];
    
        % He_A(g5)
        He =  phi_ss*(gs*gs.') + phi_sp*(gs*gp.' + gp*gs.') + phi_pp*(gp*gp.') + phi_p*Hp;
        He = - (He + He.')/2;               
    
        hc(1:4,1:4) = He;   
    else
        flag = 1;
        hc(n,n) = NaN;
    end
    
    return
end