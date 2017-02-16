w = 2 * pi * 1.47e9;
Rabi = 2 * pi * 30e6 / sqrt(2);
A = 2 * pi * 1e6 * [13.7; 6.5];
delta = 0.5 * sort([A(1) + A(2), -A(1)-A(2), A(1)-A(2), A(2)-A(1)]);
s0 = [1;0];
tspan = [0; pi / Rabi];
options = odeset('RelTol',1e-4);
rabi = [Rabi; Rabi];
wrf = [w + (delta(1)+delta(2))/2; w + (delta(3)+delta(4))/2];
phi = [0;0];
f = 0;
sf = zeros(2,4);
for i = 1:length(delta)
    [T, Y] = ode45(@(t,y) rfH(t,y,w+delta(i),rabi,wrf,phi),tspan,s0,options);
    f = f + Y(end);
    sf(:,i) = transpose(Y(end,:));
end
f = norm(f)/4;
    