function [rhoini] = comprhoini(c,f,E,I,LStype)

rhoinimin = 10^(-8);
rhoinimax = 10^8;

if LStype == 1
    sumc = 0.5 * ( sum(c(E).^2) + sum(max(c(I),0).^2) );
elseif LStype == 2
    sumc = sum(abs(c(E))) + sum(max(c(I),0));
end
rhoini = 10 * max( 1, abs( f ) ) / max( 1, sumc );

rhoini = max( rhoinimin, min( rhoini, rhoinimax ) );