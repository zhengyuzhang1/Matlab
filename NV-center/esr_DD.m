function sn=esr_DD(nD, A, omiga, bound)
%nD is number of Pi(x) DD pulses, A are coupling coefficients of nucleis, omiga is microwave intensity
t=pi/2/omiga;
nuc=length(A);
delta=(-bound:2*bound);

sz=[1,0;0,-1];
sx=[0,1;1,0];
s0=[0;1];
sn_origin=zeros(size(delta));
for i=1:length(delta)
    H=delta(i)*sz+omiga*sx;
    U=expm(-1i*t/2/nD*H);
    sf=(U*sx*U)^nD*s0;
    sn_origin(i)=norm(sf(1))^2;
end
Aeff=zeros(1,2^nuc);
sn=zeros(1,bound+1);
for i=1:2^nuc
    Aeff(i)=round((2*de2bi(i-1,nuc)-1)*A/2/2); %last /2 because we used Pauli matrices sigma(x) in hamiltonian, not 1/2*sigma(x)
    sn=sn+sn_origin(bound-Aeff(i)+1:2*bound-Aeff(i)+1);
end
sn=1-sn/2^nuc;
end