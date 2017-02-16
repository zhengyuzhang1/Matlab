delta=0:300;
t=0.1;
sigma=(1:50)';

a=sqrt(pi/2)*sigma./erf(t/2/sqrt(2)*sigma);

basicstep=500;

sz=[1,0;0,-1];
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
s0=[0;1];

U1=cell(length(sigma),length(delta));
w0=1000;

for i=1:length(sigma)
    for j=1:length(delta)
        U1{i,j}=eye(2);
        if t*sigma(i)<6
            dt=linspace(-t/2,t/2,basicstep);
        else
            dt=-t/2:6/sigma(i)/basicstep:t/2;
        end
        for k=1:length(dt)-1
            H=w0/2*sz+a(i)*exp(-dt(k)^2/2*sigma(i)^2)*cos((w0+delta(j))*dt(k))*sx;
            U1{i,j}=expm(-1i*(dt(k+1)-dt(k))*H)*U1{i,j};
        end
    end
end
save U1