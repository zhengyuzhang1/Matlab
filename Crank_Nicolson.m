function [finalstate,states] = Crank_Nicolson(Hamiltonian, t, steps, state0)

dt = t / steps;
n = length(state0);
states = zeros(length(state0),length(steps));
states(:,1) = state0;
opts.SYM = true;
for i = 2:steps
    states(:,i) = linsolve(eye(n) + 1i * dt / 2 * Hamiltonian(i*dt), (eye(n) - 1i * dt / 2 * Hamiltonian(i*dt)) * states(:,i-1), opts);
end
finalstate = states(:,end);
end