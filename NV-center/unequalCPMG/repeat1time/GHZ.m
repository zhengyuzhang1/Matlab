% w1 = w1([1;2;4;5],:);
pulse = taominall{1,6,30}([1;2;4;5],1:end-1);
U0 = eye(16);
U1 = eye(16);
for i = 1:4
    [~, u0, u1] = unequaltimeCnot(w,w1,pulse(i,:),i,1);
    U0 = U0 * kron(u0(:,:,1),kron(u0(:,:,2),kron(u0(:,:,3),u0(:,:,4))));
    U1 = U1 * kron(u1(:,:,1),kron(u1(:,:,2),kron(u1(:,:,3),u1(:,:,4))));
end
U = blkdiag(U0,U1);
s0 = kron(1/sqrt(2)*[1;1],sparse(1,1,1,16,1));
sT = sparse([1,32],[1,1],1/sqrt(2),32,1);
fid = norm(sT' * U * s0);