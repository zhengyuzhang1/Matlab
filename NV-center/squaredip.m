A=Aall(2:11);
width=80;
nuc=length(A);
n=2001;
delta=linspace(0,2000,n);

sdip=zeros(1,n);
Aeff=zeros(1,2^nuc);
for i=1:2^nuc
    Aeff(i)=(2*de2bi(i-1,nuc)-1)*A/2;
    sdip=sdip+double(abs(delta-Aeff(i))<width/2);
end
sdip=1-sdip/2^nuc;
figure;
plot(delta,sdip)