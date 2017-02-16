function dy = dephasingNMRdy(t,y,w,delta,gamma2)
%gamma2 = 1/T2
dy = zeros(3,1);
dy(1) = delta * y(2) - y(1) * t * gamma2;
dy(2) = -w(t) * y(3) - delta * y(1) - y(2) * t * gamma2;
dy(3) = w(t) * y(2);
end