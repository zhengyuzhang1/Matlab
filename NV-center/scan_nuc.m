delta=300;
%t=0.1:0.01:0.3;
t=0.1;
sigma=1./t'*linspace(5,10,20);

a=sqrt(pi/8)*sigma./erf(repmat(t',1,size(sigma,2))/2/sqrt(2).*sigma);

basicstep=100;

sz=[1,0;0,-1];
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
s0=[0;1];

sop=zeros(size(sigma));
w0=1000;

for i=1:length(t)
    for j=1:size(sigma,2)
        s0=[0;1];
        if t(i)*sigma(i,j)<6
            dt=linspace(-t(i)/2,t(i)/2,basicstep);
        else
            dt=-t(i)/2:6/sigma(i,j)/basicstep:t(i)/2;
%             tmp=linspace(3/sigma(i,j),t(i)/2,50);
%             tmp=tmp(2:end);
%             dt=[fliplr(-tmp),linspace(-3/sigma(i,j),3/sigma(i,j),basicstep),tmp];
        end
        for l=1:length(dt)-1
            H=w0/2*sz+a(i,j)*exp(-dt(l)^2/2*sigma(i,j)^2)*2*cos((w0+delta)*dt(l))*sx;
            s0=expm(-1i*(dt(l+1)-dt(l))*H)*s0;
        end
        sop(i,j)=norm(s0(1))^2;
    end
end
%  [popt,sigmaopt]=min(sop,[],2);
%  figure;plot(t*delta,popt);
%  figure;plot(t*delta,sigmaopt);
%%
figure;
plot(diff(sop,2,2)');
%%
figure;plot(t*sigma(1:100),sop(1,1:100));