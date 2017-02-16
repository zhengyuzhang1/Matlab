function [ H ] = Sparse_Heisenlike_local( N, co )

n = 2^N;
I = (1:n)';
base = n - I;
nonzn = size(co,1);
J = zeros(n,2*nonzn);
Val = zeros(n,2*nonzn);
basei = (N-co(:,1))';
basej = (N-co(:,2))';
bitvali = repmat(2.^basei,n,1);
bitvalj = repmat(2.^basej,n,1);
ibit = 2 * bsxfun(@bitget,base,basei+1) - 1;
jbit = 2 * bsxfun(@bitget,base,basej+1) - 1;
J(:,1:nonzn) = repmat(I,1,nonzn) + ibit .* bitvali + jbit .* bitvalj;
Val(:,1:nonzn) = repmat(co(:,3)',n,1) - repmat(co(:,4)',n,1) .* ibit .* jbit;
J(:,nonzn+1:end) = repmat(I,1,nonzn);
Val(:,nonzn+1:end) = repmat(co(:,5)',n,1) .* ibit .* jbit;
H = sparse( repmat(I,2*nonzn,1),J(:),Val(:) );

end