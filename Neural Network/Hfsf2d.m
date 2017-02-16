function [ H ] = Hfsf2d( N, gamma, lambda )
% Hamiltonian for free spinless fermion model with periodic
% boundary condition for N * N sites in 2D

sx = [0 1;1 0]; sy = [0 -1i;1i 0]; sz = [1 0;0 -1];
dim = 2^(N*N);
H = zeros(dim);
for i = 1:N*N
    H = H - lambda * kron(kron(eye(2^(i-1)),eye(2)+sz),eye(2^(N*N-i)));
end
for i = 1:N
    for j = 1:N
        i1 = (i-1)*N+j;
        i2 = mod(i,N)*N+j;
        a = min(i1,i2);
        b = max(i1,i2);
        hzstring = 1;
        for k = 1:b-a-1
            hzstring = kron(hzstring,sz);
        end
        if a == i2
            beta = -1;
        else
            beta = 1;
        end
        H = H + (-1)^abs(a-b-1) * ...
            ( 0.5*(1-beta*gamma)*kron(kron(kron(kron(eye(2^(a-1)),sx),hzstring),sx),eye(2^(N*N-b))) ...
            + 0.5*(1+beta*gamma)*kron(kron(kron(kron(eye(2^(a-1)),sy),hzstring),sy),eye(2^(N*N-b))) );
        i2 = (i-1)*N+mod(j,N)+1;
        a = min(i1,i2);
        b = max(i1,i2);
        hzstring = 1;
        for k = 1:b-a-1
            hzstring = kron(hzstring,sz);
        end
        if a == i2
            beta = -1;
        else
            beta = 1;
        end
        H = H + (-1)^abs(a-b-1) * ...
            ( 0.5*(1-beta*gamma)*kron(kron(kron(kron(eye(2^(a-1)),sx),hzstring),sx),eye(2^(N*N-b))) ...
            + 0.5*(1+beta*gamma)*kron(kron(kron(kron(eye(2^(a-1)),sy),hzstring),sy),eye(2^(N*N-b))) );
    end
end

