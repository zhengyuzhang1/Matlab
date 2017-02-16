function [state,up,H] = mblexact(ndt,step,A,up)

n = size(A,1); %number of spins
%initial state
if nargin < 4 
    up = randi([1,n],1);
end
state = zeros(n,step);
state(up,1) = 1;

H = diag(sum(sum(A(:,:,end))) / 2 * ones(n,1));

for i = 1:n
    H(i,i) = H(i,i) - 2 * sum(A(i,:,end));
    for j = i + 1:n
        H(i,j) = A(i,j,1) + A(i,j,4);
        H(j,i) = conj(H(i,j));
    end
end        
H = H - mean(diag(H)) * eye(n);
%random Hermitian
% for i = 1:n
%     H(i,i) = 1000 * rand(1);
%     for j = i + 1:n
%         H(i,j) = 1000 * (rand(1) + 1i * rand(1));
%         H(j,i) = conj(H(i,j));
%     end
% end

dt = ndt / max(abs(H(up,:)));
U = expm(-1i * dt * H);

for t = 2:step
    state(:,t) = U * state(:,t-1);
end
       