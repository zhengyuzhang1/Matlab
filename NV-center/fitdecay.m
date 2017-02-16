figure;plot(t,data);hold on;
fun1=@(a,xdata) exp(-a*(32*xdata));
pa1=lsqcurvefit(fun1,0,t(DataIndex)*1e9,data(DataIndex));
plot(t,fun1(pa1,t*1e9));
fun2=@(a,xdata) exp(-a*(32*xdata).^3);
pa2=lsqcurvefit(fun2,0,t(DataIndex)*1e3,data(DataIndex));
plot(t,fun2(pa2,t*1e3));
text(0.6e-5,1,['exp(-',num2str(1e9*pa2,4),'t^3)']);
text(0.9e-5,0.8,['exp(-',num2str(1e9*pa1,4),'t)']);