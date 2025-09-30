function [v] = reshapevector(v)

global problemID nballs

if ( problemID == 2 )

    if ( isstruct(v) )
        % Transform v in a n-dimensional vector
        va = cell2mat(v.uv);
        vb(1:nballs,1) = va(1:2:end,1);
        vb(nballs+1:2*nballs,1) = va(2:2:end,1);
        v = [vb; v.s; v.r];
        
    else
        % Trasform v in the struct of the Manifold
        xa(:,1) = v(1:2*nballs);
        xa = reshape(xa, nballs, [])';
        xb.uv = mat2cell(xa, 2, ones(1, nballs))';
        xb.s = v(2*nballs+1:3*nballs);
        xb.r = v(end);
        v = [];
        v = xb; 
    end

    return
end

% ==================================================================

if ( problemID == 3 )

    if ( isstruct(v) )
        % Transform v in a n-dimensional vector
        v = [v.A(:);v.b(:)];
    else
        % Trasform v in the struct of the Manifold
        A = reshape(v(1:4), [2, 2]);
        b = reshape(v(5:6), [2, 1]);
    
        v = [];
        v.A = A;
        v.b = b;
    end

    return
end