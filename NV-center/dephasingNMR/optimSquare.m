function infid = optimSquare(P,delta,gamma2)
%gaussP = [gaussA, gaussSigma, gaussT];
w = @(t) P(1);
infid = 0.5 * (-dephasingNMR(w,0,gamma2,P(2)) + 1 + dephasingNMR(w,delta,gamma2,P(2)));
end