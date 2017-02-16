M = zeros(1000,3);
for i = 1:40
    filename = ['dipoledynamic_',num2str(i),'.mat'];
    load(filename,'mx','my','mz');
    M(:,1) = M(:,1) + mx';
    M(:,2) = M(:,2) + my';
    M(:,3) = M(:,3) + mz';
end
M = M / 40 * 2;
save dipoledynamic
%% 
fun = @(x,xdata) exp(-x(1) * xdata) + x(2) ./ xdata;
xdata = (0:0.001:0.999)';
ydata = M(:,1);
x0 = [2,1];
x1 = lsqcurvefit(fun,x0,xdata,ydata);
figure;plot(xdata,[M,fun(x1,xdata)]);
%%
figure;plot(xdata,[M,fun(x1,xdata),fun(x2,xdata)]);
legend('mx','my','mz',['exp(-',num2str(x1),'x)'],['exp(-',num2str(x2),'x)']);
xlabel('t*A');
ylabel('P');