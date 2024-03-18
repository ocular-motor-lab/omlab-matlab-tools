function dS = AttractorNetworkUpdateS( t, S, x, A)

    m = height(A)-1;
    n = length(S)/m;

    S = reshape(S, n, m );
    S = [S ones(n,1)];

    dS = zeros(n,m);
    x = [x;1];
    
    for i=1:n
        si = S(i,:)';
        dist = S-repmat(si',n,1);
        d =  - (4*si'*A*si*A*si) - 1*mean(1./(1+sum(dist.*dist,2)).*(dist))'  -  (si-x);
        dS(i,:)= d(1:end-1)';
    end


    dS = dS(:);

end