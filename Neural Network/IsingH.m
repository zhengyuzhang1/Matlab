function [ H ] = IsingH( N, h )
% Hamiltonian for 1D Ising model with periodic
% boundary condition for N sites

sx = [0 1;1 0]; sy = [0 -1i;1i 0]; sz = [1 0;0 -1];
dim = 2^N;
H = zeros(dim);
for i = 1:N-1
    H = H - kron(kron(kron(eye(2^(i-1)),sz),sz),eye(2^(N-i-1)));
end
H = H - kron(sz,kron(eye(2^(N-2)),sz));
for i = 1:N
   H = H - h * kron(kron(eye(2^(i-1)),sx),eye(2^(N-i)));
end
end