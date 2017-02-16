function [statef,fidelity]=ghz(w,azx,pulse,phase)
non=size(azx,1);
stata0=[ones(1,non);zeros(1,non)];
stata1=[zeros(1,non);ones(1,non)];
a=1;b=1;
for i=1:non
    gate0=eye(2);gate1=eye(2);
    for j=1:non
        gate0=gate0*(U0(w,azx(i,:),pulse(j,3))^pulse(j,4))*(U0(w,azx(i,:),pulse(j,1))^pulse(j,2));
        gate1=gate1*(U1(w,azx(i,:),pulse(j,3))^pulse(j,4))*(U1(w,azx(i,:),pulse(j,1))^pulse(j,2));
        aa=abs(gate0(1));bb=abs(gate1(2));
    end
    a=a*abs(gate0(1));
    b=b*abs(gate1(2));
end
statef=1;
fidelity=(a+b)/2;
end
