function fid = dephasingNMR(w,delta,gamma2,t)

fun = @(t,y) dephasingNMRdy(t,y,w,delta,gamma2);
options = odeset('RelTol',1e-5);
[~,Y] = ode45(fun, [0 t], [0;0;1],options);
% if abs(Y(end))>1
%     Y(end) = sign(Y(end));
% end
fid = 0.5 * (1 - Y(end));

end