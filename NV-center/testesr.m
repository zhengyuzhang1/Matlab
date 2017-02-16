omiga=20;
t=pi/2/omiga;
A=Aall(2:2);
nuc=length(A);
n=1000;
delta=linspace(0,5000,n);
%sn=zeros(1,n);
%sneff=zeros(1,n);
Aeff=zeros(1,2^nuc);
for i=1:2^nuc
    Aeff(i)=(2*de2bi(i-1,nuc)-1)*A/2;
    %tmpd=Aeff(i)-delta;
    %rabif=sqrt(tmpd.^2+omiga^2);
    %sn=sn+omiga^2./rabif.^2.*sin(rabif.*t/2).^2;
    %sneff=sneff+4*omiga^2./rabif.^6.*(tmpd.^2+omiga^2.*cos(pi/4/omiga.*rabif)).^2.*sin(pi/4/omiga.*rabif).^2;
end
%sn=1-sn;%/2^nuc;
%sneff=1-sneff/2^nuc;

sz=[1,0;0,-1];
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
s0=[0;1];
sn8=zeros(size(delta));
for i=1:n
    for j=1:2^nuc
        H=(Aeff(j)/2-delta(i))*sz+omiga*sx;
        U=expm(-1i*t/16*H);
        sf=(U*sx*U)^8*s0;
        sn8(i)=sn8(i)+norm(sf(1))^2;
    end
end
sn8=1-sn8/2^nuc;

% snde=zeros(1,n);
% tn=100;
% dt=t/4/tn;
% for i=1:n
%     for j=2:2%2^nuc
%         s0=[0;1];
%         for k=1:tn
%             H=(Aeff(j)+w0)*sz+omiga*2*cos((w0+delta(i))*dt*k)*sx;
%             s0=expm(-1i*dt*H)*s0;
%         end
%         s0=2*sx*s0;
%         for k=1:2*tn
%             H=(Aeff(j)+w0)*sz+omiga*2*cos((w0+delta(i))*(t/4+dt*k))*sx;
%             s0=expm(-1i*dt*H)*s0;
%         end
%         s0=2*sx*s0;
%         for k=1:tn
%             H=(Aeff(j)+w0)*sz+omiga*2*cos((w0+delta(i))*(3*t/4+dt*k))*sx;
%             s0=expm(-1i*dt*H)*s0;
%         end
%         snde(i)=snde(i)+norm(s0(1))^2;
%     end
% end
% figure;
% plot(delta,1-snde);
figure;
plot(delta,sn8,[A(1)/4,A(1)/4],[0,1])
