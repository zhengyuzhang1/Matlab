function fid = aveUGateFid(U0, U1)

d = size(U0,1);
X = diag(ones(d-1,1),-1); X(1,end) = 1;
Z = diag(exp(2i*pi/d*(0:d-1)));

fid = d^2;
for i = 0:d-1
    for j = 0:d-1
        Ub = X^i * Z^j;
        fid = fid + trace(U0*Ub'*U0'*U1*Ub*U1');
    end
end
fid = fid / d^2 / (d+1);

end