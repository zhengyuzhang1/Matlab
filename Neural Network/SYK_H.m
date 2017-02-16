function [ H ] = SYK_H( N, J )
% Hamiltonian for anti-ferromagnetic Heisenberg model with periodic
% boundary condition for N sites

sx = [0 1;1 0]; sy = [0 -1i;1i 0]; sz = [1 0;0 -1];
dim = 2^N;
H = zeros(dim);
for i = 1:N-1
    for j = i+1:N
    H = H + J(i,j)*( kron(kron(kron(kron(eye(2^(i-1)),sx),eye(2^(j-i-1))),sx),eye(2^(N-j))) ...
        + kron(kron(kron(kron(eye(2^(i-1)),sy),eye(2^(j-i-1))),sy),eye(2^(N-j))) ...
        + kron(kron(kron(kron(eye(2^(i-1)),sz),eye(2^(j-i-1))),sz),eye(2^(N-j))) );

    end
end

end