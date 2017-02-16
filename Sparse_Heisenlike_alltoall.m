function [ H ] = Sparse_Heisenlike_alltoall( N, co )

n = 2^N;
I = (1:n)';
base = n - I;
nonzn = 2*N*(N-1)/2;
J = zeros(n,nonzn);
Val = zeros(n,nonzn);
idxj = 0;
bitval = 2.^((N-1):-1:0);
for i = 1:N-1
    for j = i+1:N
        ibit = 2*bitget(base,N+1-i)-1;
        jbit = 2*bitget(base,N+1-j)-1;
        idxj = idxj + 1;
        J(:,idxj) = I + ibit * bitval(i) + jbit * bitval(j);
        Val(:,idxj) = co(i,j,1) - co(i,j,2) * ibit .* jbit;
        idxj = idxj + 1;
        J(:,idxj) = I;
        Val(:,idxj) = co(i,j,3) * ibit .* jbit;
    end
end
H = sparse( repmat(I,nonzn,1),J(:),Val(:) );

end