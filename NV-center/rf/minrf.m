w = 2 * pi * 1.47e9;
A = 2 * pi * 1e6 * [13.7,6.5];
maxRabi = 2 * pi * 30e6;
dw0 = (abs(A(1)+A(2))+abs(A(1)-A(2)))/2; 
delta = pi * 1e6 * [-20.2, -7.2, 7.2, 20.2];
%wrf = [w-x(1); w+x(1)];
%wrf = [w-0.5*(delta(3)+delta(4)); w+0.5*(delta(3)+delta(4))];
%%
x01 = pi / maxRabi;
fun1 = @(x) searchrf(w,x,0,delta);
[x1, xval1, exitflag1] = fminsearch(fun1,x01);
%%
x02 = pi / maxRabi / sqrt(2);
fun2 = @(x) searchrf([w-dw0;w+dw0],x,[0;0],delta);
[x2, xval2, exitflag2] = fminsearch(fun2,x02);
%%
x03 = [dw0,pi / maxRabi / sqrt(2)];
fun3 = @(x) searchrf([w-x(1);w+x(1)],x(2),[0;0],delta);
[x3, xval3, exitflag3] = fminsearch(fun3,x03);
%%
x04 = [-dw0,dw0,pi / maxRabi / sqrt(2)];
fun4 = @(x) searchrf([w+x(1);w+x(2)],x(3),[0;0],delta);
[x4, xval4, exitflag4] = fminsearch(fun4,x04);
%%
% x05 = [dw0, pi / maxRabi / sqrt(2)];
% fun5 = @(x) searchrf([w-x(1);w+x(1)],x(2),delta);
% options = optimoptions(@fmincon,'TolX',1e-10);
% [x5, xval5, exitflag5] = fmincon(fun5,x05,[],[],[],[],[dw0/2,0],[2*dw0,1e-6],[],options);
%%
x05 = [delta(3),delta(4), pi / maxRabi / 2 *1.2];
fun5 = @(x) searchrf([w-x(1);w-x(2);w+x(1);w+x(2)],x(3),delta);
options = optimoptions(@fmincon,'TolX',1e-10);
[x5, xval5, exitflag5] = fmincon(fun5,x05,[],[],[],[],[delta(3)/2,delta(4)/2,0],[2*delta(3),2*delta(4),1e-6],[],options);