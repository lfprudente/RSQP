function checkd(n,m,xini,l,u)

macheps12 = 10^(-8);

rng(123456);

for i = 1:n
    if ( l(i) < xini(i) && xini(i) < u(i) )
        x(i) = xini(i) + macheps12 * ( 2.0d0 * rand - 1.0d0 ) * max( 1.0d0, abs( xini(i) ) );
    elseif ( xini(i) == l(i) )
        x(i) = xini(i) + macheps12 * rand * max( 1.0d0, abs( xini(i) ) );
    else
        x(i) = xini(i) - macheps12 * rand * max( 1.0d0, abs( xini(i) ) );
    end
    x(i) = max( l(i), min( x(i), u(i) ) );
end
x = x';

fprintf(' Derivatives will be tested at the perturbed initial guess:\n')
for i = 1:n
    fprintf(' x( %6i ) = %15.8f \n',i,x(i))
end

% CHECK OBJECTIVE FUNCTION GRADIENT

fprintf('\n Would you like to check subroutine evalg?\n')
prompt = " Type Y(es), N(o), A(bort checking) or S(kip checking this subroutine): ";
answer = input(prompt,"s");
if isempty(answer)
    answer = 'S';
end

if ( answer == 'A' || answer == 'a' )
    return
elseif ( answer == 'Y' || answer == 'y' )
    checkg(n,x)
end

% CHECK CONSTRAINT GRADIENTS

% if ( m == 0 ) 
%     fprintf('\n')
%     fprintf('Press any key to continue.\n')
%     fprintf('\n\n')
%     pause
%     return
% end

if ( m > 0 ) 
    fprintf('\n Would you like to check subroutine evalnc?\n')
    prompt = " Type Y(es), N(o), A(bort checking) or S(kip checking this subroutine): ";
    answer = input(prompt,"s");
    if isempty(answer)
        answer = 'S';
    end
    
    if ( answer == 'A' || answer == 'a' )
        return
    elseif ( answer == 'Y' || answer == 'y' )
        for i = 1:m
            fprintf('\n Check gradient of constraint %5i ?\n',i)
            prompt = " Type Y(es), N(o), A(bort checking) or S(kip checking this subroutine): ";
            answerc = input(prompt,"s");
            if isempty(answer)
                answerc = 'S';
            end
            
            if ( answerc == 'A' || answerc == 'a' )
                return
            elseif ( answerc == 'Y' || answerc == 'y' )
                checknc(n,x,i)
            end
        end
    end
end

% CHECK HESSIAN OF THE OBJECTIVE FUNCTION

fprintf('\n Would you like to check subroutine evalh?\n')
prompt = " Type Y(es), N(o), A(bort checking) or S(kip checking this subroutine): ";
answer = input(prompt,"s");
if isempty(answer)
    answer = 'S';
end

if ( answer == 'A' || answer == 'a' )
    return
elseif ( answer == 'Y' || answer == 'y' )
    checkh(n,x)
end

% CHECK HESSIAN OF THE CONSTRAINTS

if ( m > 0 ) 
    fprintf('\n Would you like to check subroutine evalhc?\n')
    prompt = " Type Y(es), N(o), A(bort checking) or S(kip checking this subroutine): ";
    answer = input(prompt,"s");
    if isempty(answer)
        answer = 'S';
    end
    
    if ( answer == 'A' || answer == 'a' )
        return
    elseif ( answer == 'Y' || answer == 'y' )
        for i = 1:m
            fprintf('\n Check the Hessian of constraint %5i ?\n',i)
            prompt = " Type Y(es), N(o), A(bort checking) or S(kip checking this subroutine): ";
            answerc = input(prompt,"s");
            if isempty(answer)
                answerc = 'S';
            end
            
            if ( answerc == 'A' || answerc == 'a' )
                return
            elseif ( answerc == 'Y' || answerc == 'y' )
                checkhc(n,x,i)
            end
        end
    end
end

fprintf('\n')
fprintf('Press any key to continue.\n')
fprintf('\n\n')
pause

end

% ==================================================================

function checkg(n,x)

[g,~] = evalg(n,x);

fprintf('\n')
fprintf(' Gradient vector of the objective function.\n')
fprintf(' Index             evalg  Central diff (two different steps)    Absolute error\n')


maxerr = 0.0d0;

eps = 10^(-16/3);

for i = 1:n
	 tmp  = x(i);

	 step1 = eps * max( abs( tmp ), 1.0d0 );

	 x(i) = tmp + step1;
     [fplus,~] = evalf(n,x);
	 
	 x(i) = tmp - step1;
     [fminus,~] = evalf(n,x);

	 gdiff1 = ( fplus - fminus ) / ( 2.0d0 * step1 );

	 step2 = eps * max( abs( tmp ), 1.0d-03 );

	 x(i) = tmp + step2;
     [fplus,~] = evalf(n,x);

	 x(i) = tmp - step2;
	 [fminus,~] = evalf(n,x);
	 
	 x(i) = tmp;

	 gdiff2 = ( fplus - fminus ) / ( 2.0d0 * step2 );

	 tmp = min( abs( g(i) - gdiff1 ), abs( g(i) - gdiff2 ) );
     fprintf(' %5i   %15.8f   %15.8f   %15.8f   %15.8f\n',i,g(i),gdiff1,gdiff2,tmp)

	 maxerr = max( maxerr, tmp );

end

fprintf('\n Maximum absolute error = %15.8f \n',maxerr)

end

% ==================================================================

function checknc(n,x,ind)

[g,~] = evalnc(n,x,ind);

fprintf('\n')
fprintf(' Gradient vector of the objective function.\n')
fprintf(' Index             evalg  Central diff (two different steps)    Absolute error\n')


maxerr = 0.0d0;

eps = 10^(-16/3);

for i = 1:n
	 tmp  = x(i);

	 step1 = eps * max( abs( tmp ), 1.0d0 );

	 x(i) = tmp + step1;
     [fplus,~] = evalcc(n,x,ind);
	 
	 x(i) = tmp - step1;
     [fminus,~] = evalcc(n,x,ind);

	 gdiff1 = ( fplus - fminus ) / ( 2.0d0 * step1 );

	 step2 = eps * max( abs( tmp ), 1.0d-03 );

	 x(i) = tmp + step2;
     [fplus,~] = evalcc(n,x,ind);

	 x(i) = tmp - step2;
     [fminus,~] = evalcc(n,x,ind);
	 
	 x(i) = tmp;

	 gdiff2 = ( fplus - fminus ) / ( 2.0d0 * step2 );

	 tmp = min( abs( g(i) - gdiff1 ), abs( g(i) - gdiff2 ) );
     fprintf(' %5i   %15.8f   %15.8f   %15.8f   %15.8f\n',i,g(i),gdiff1,gdiff2,tmp)

	 maxerr = max( maxerr, tmp );

end

fprintf('\n Maximum absolute error = %15.8f \n',maxerr)

end

% ==================================================================

function checkh(n,x)

eps = 1.0d-8;

% Compute the gradient of the objective function at x

[g,~] = evalg(n,x);

% Compute the Hessian of the objective function at x

[h,~] = evalh(n,x);

% Test column by column

fprintf('\n Hessian matrix of the objective function column by column.\n')

maxerr = 0.0d0;

for j = 1:n

    tmp = x(j);
    
    step1 = eps * max( abs( tmp ), 1.0d0 );
    
    x(j) = tmp + step1;

    [gplus1,~] = evalg(n,x);
    
    step2 = eps * max( abs( tmp ), 1.0d-03 );
    
    x(j) = tmp + step2;

    [gplus2,~] = evalg(n,x);
    
    x(j) = tmp;

    fprintf('\n Column: %6i\n',j)
    
    maxcoe(j) = 0.0d0;
    
    nullcol = true;
    
    for i = 1:n
        if ( i >= j )
            elem = h(i,j);
        else
            elem = h(j,i);
        end
        hdiff1 = ( gplus1(i) - g(i) ) / step1;
        hdiff2 = ( gplus2(i) - g(i) ) / step2;
        tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) );
        if ( elem ~= 0.0d0 || hdiff1 ~= 0.0d0 ||hdiff2 ~= 0.0d0 )
            if ( nullcol )
                nullcol = false;
                fprintf('\n Index             evalh   Incr. Quoc. (two different steps)    Absolute error\n')
            end
            fprintf(' %5i   %15.8f   %15.8f   %15.8f   %15.8f\n',i,elem,hdiff1,hdiff2,tmp)
        end
        maxcoe(j) = max( maxcoe(j), tmp );
    end
    
    maxerr = max( maxerr, maxcoe(j) );
    
    if ( nullcol )
        fprintf(' All the elements of this column are null.\n')
    else
        fprintf(' Maximum absolute error = %15.8f\n',maxcoe(j))
    end

end

for j = 1:n
    fprintf(' Column %6i Maximum absolute error = %15.8f\n',j,maxcoe(j))
end
fprintf(' Overall maximum absolute error = %15.8f\n',maxerr)

end

% ==================================================================

function checkhc(n,x,ind)

eps = 1.0d-8;

% Compute the gradient of the objective function at x

[g,~] = evalnc(n,x,ind);

% Compute the Hessian of the objective function at x

[h,~] = evalhc(n,x,ind);

% Test column by column

fprintf('\n Hessian matrix of constraint %5i column by column.\n',ind)

maxerr = 0.0d0;

for j = 1:n

    tmp = x(j);
    
    step1 = eps * max( abs( tmp ), 1.0d0 );
    
    x(j) = tmp + step1;

    [gplus1,~] = evalnc(n,x,ind);
    
    step2 = eps * max( abs( tmp ), 1.0d-03 );
    
    x(j) = tmp + step2;

    [gplus2,~] = evalnc(n,x,ind);
    
    x(j) = tmp;

    fprintf('\n Column: %6i\n',j)
    
    maxcoe(j) = 0.0d0;
    
    nullcol = true;
    
    for i = 1:n
        if ( i >= j )
            elem = h(i,j);
        else
            elem = h(j,i);
        end
        hdiff1 = ( gplus1(i) - g(i) ) / step1;
        hdiff2 = ( gplus2(i) - g(i) ) / step2;
        tmp = min( abs( elem - hdiff1 ), abs( elem - hdiff2 ) );
        if ( elem ~= 0.0d0 || hdiff1 ~= 0.0d0 || hdiff2 ~= 0.0d0 )
            if ( nullcol )
                nullcol = false;
                fprintf('\n Index             evalh   Incr. Quoc. (two different steps)    Absolute error\n')
            end
            fprintf(' %5i   %15.8f   %15.8f   %15.8f   %15.8f\n',i,elem,hdiff1,hdiff2,tmp)
        end
        maxcoe(j) = max( maxcoe(j), tmp );
    end
    
    maxerr = max( maxerr, maxcoe(j) );
    
    if ( nullcol )
        fprintf(' All the elements of this column are null.\n')
    else
        fprintf(' Maximum absolute error = %15.8f\n',maxcoe(j))
    end

end

for j = 1:n
    fprintf(' Column %6i Maximum absolute error = %15.8f\n',j,maxcoe(j))
end
fprintf(' Overall maximum absolute error = %15.8f\n',maxerr)

end