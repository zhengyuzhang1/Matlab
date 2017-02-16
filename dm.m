function f = dm(statevector)
%density matrix of purestate with statevector a column vector

norms = statevector / norm(statevector);
f = norms * norms';
end