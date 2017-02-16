function f = searchrf(wrf, trf, phi, delta)

w = 2 * pi * 1.47e9;
nrf = length(wrf);
Rabi = 2 * pi * 30e6 / sqrt(nrf) * ones(nrf,1);
s0 = [1;0];
options = odeset('RelTol',1e-4);
f = 0;
for i = 1:length(delta)
    [~, Y] = ode45(@(t,y) rfH(t,y,w+delta(i),Rabi,wrf,phi),[0 trf],s0,options);
    f = f + Y(end);
end
f = 1-norm(f)^2/16;
end
    