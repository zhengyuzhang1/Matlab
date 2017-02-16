delta=0:2:500;
t=0.1;
sigma=1./t'*(8);

a=sqrt(pi/2)*sigma./erf(repmat(t',1,size(sigma,2))/2/sqrt(2).*sigma);

basicstep=300;

sz=[1,0;0,-1];
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
s0=[0;1];

soptmp=zeros(length(t),size(sigma,2),length(delta));
w0=300;

for i=1:length(t)
    for j=1:size(sigma,2)
        for k=1:length(delta)
            s0=[0;1];
            if t(i)*sigma(i,j)<6
                dt=linspace(-t(i)/2,t(i)/2,basicstep);
            else
                dt=-t(i)/2:6/sigma(i,j)/basicstep:t(i)/2;
            end
            for l=1:length(dt)-1
                H=w0/2*sz+a(i,j)*exp(-dt(l)^2/2*sigma(i,j)^2)*cos((w0+delta(k))*dt(l))*sx;
                s0=expm(-1i*(dt(l+1)-dt(l))*H)*s0;
            end
            soptmp(i,j,k)=norm(s0(1))^2;
        end
    end
end