t1=t(7);
delta1=0:2:500;
sigma1=sigma(7,2);

a=sqrt(pi/8)*sigma1./erf(repmat(t1',1,size(sigma1,2))/2/sqrt(2).*sigma1);
dt=min(t1)/300;
tn=round(t1/dt);

sz=[1,0;0,-1];
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
s0=[0;1];

sn=zeros(length(t1),length(sigma1),length(delta1));
w0=1000;

for i=1:length(t1)
    for j=1:length(sigma1)
        for k=1:length(delta1)
            s0=[0;1];
            for l=0:tn(i)-1
                H=w0/2*sz+a(i,j)*exp(-(l*dt-t1(i)/2)^2/2*sigma1(j)^2)*2*cos((w0+delta1(k))*dt*l)*sx;
                s0=expm(-1i*dt*H)*s0;
            end
            sn(i,j,k)=norm(s0(1))^2;
        end
    end
end
figure;plot(reshape(sn(1,:,:),1,251)');
%% plot truncated gauss pulse shape
t1=t1(10);
sigma1=sigma1(10,1);
a=sqrt(pi/8)*sigma1/erf(t1/2/sqrt(2)*sigma1);
%figure;plot(linspace(-t(i)/2,t(i)/2,tn(i)+1),sqrt(pi/2)/sigma(j)*exp(-linspace(-t(i)/2,t(i)/2,tn(i)+1).^2/2/sigma(j)^2));
figure;plot(linspace(-t1/2,t1/2,100),a*exp(-linspace(-t1/2,t1/2,100).^2/2*sigma1^2));
%%
t1=t;
delta1=0;
sigma1=sigma;
alpha=linspace(0,7*pi/16,8);

a=sqrt(pi/8)*sigma1./erf(repmat(t1',1,size(sigma1,2))/2/sqrt(2).*sigma1);
basicstep = 300;

sz=[1,0;0,-1];
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
s0=[0;1];

sn=zeros(length(t),length(sigma1),length(alpha));
w0=1000;

for i = 1:length(t1)
    for j = 1:size(sigma1,2)
        for k = 1:length(alpha)
            s0 = [0;1];
            if t1(i) * sigma1(i,j) < 6
                dt = linspace(-t1(i)/2,t1(i)/2,basicstep);
            else
                dt = -t1(i)/2 : 6/sigma1(i,j)/basicstep : t1(i)/2;
            end
            for l = 1:length(dt)-1
                 H = w0 / 2 * sz + a(i,j) * exp(-dt(l)^2 / 2 * sigma1(i,j)^2) * 2 * cos((w0 + delta1) * dt(l)) * (cos(alpha(k)) * sx + sin(alpha(k)) * sz);
                 s0 = expm(-1i * (dt(l+1) - dt(l)) * H) * s0;
            end
            sn(i,j,k)=norm(s0(1))^2;
        end
    end
end
%%
i = 1;
j = 1;
k = 4;
figure;plot(t,reshape(sn(:,j,k),1,length(t)));
figure;plot(sigma(i,:),reshape(sn(i,:,k),1,size(sigma,2)));
figure;plot(alpha,reshape(sn(i,j,:),1,length(alpha)));