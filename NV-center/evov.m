function sn=evov(H, s0, deltat, n)
% direct integration to get the state sn=U(deltat*n)*s0, deltat being the
% step in time, n the total steps, H is hamiltonian and s0 initial state
dU=expm(-1i*deltat*H);
sn=dU^n*s0;
end