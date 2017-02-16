function infid = optimGauss(gaussP,delta,gamma2)
%gaussP = [gaussA, gaussSigma, gaussT];
w = @(t) gaussP(1) / sqrt(2*pi) / gaussP(2) * exp(-(t-gaussP(3)/2)^2/2/gaussP(2)^2);
infid = 1/7 * (1 - dephasingNMR(w,0,gamma2,gaussP(3)) ...
    + 2 * (dephasingNMR(w,delta,gamma2,gaussP(3)) + dephasingNMR(w,1.5*delta,gamma2,gaussP(3)) + dephasingNMR(w,2*delta,gamma2,gaussP(3))));
end