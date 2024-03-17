function dS = AttractorNetworkUpdateS( t, S, x, A)

    m = height(A)-1;
    n = length(S)/m;

    S = reshape(S, n, m );
    S = [S ones(n,1)];

    dS = zeros(n,m);
    x = [x;1];
    
    for i=1:n
        si = S(i,:)';
        d =  - 5*(4*si'*A*si*A*si) - 0.5*mean(1./(1+S*si).*(S-repmat(si',n,1)))'  -  (si-x);
        dS(i,:)= d(1:end-1)';
    end


    dS = dS(:);

end