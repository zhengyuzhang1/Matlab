sgmz=[1,0;0,-1];
sgmx=[0,1;1,0];
sgmy=[0,-1i;1i,0];
npipulse=16;
fliprate=0;
x0=ones(size(t));
for nuc=1:size(sample,1)
    H(:,:,1)=w*sgmz/2;
    H(:,:,2)=(sample(nuc,1)+w)*sgmz/2+sample(nuc,2)*sgmx/2;
    x1nuc=zeros(size(t));
    for i=1:numel(t)
        flipt=zeros(10000,1);
        qpipulse=zeros(size(flipt));
        T=-log(rand())/fliprate;
        flipn=0;
        while T<2*npipulse*t(i)
            flipn=flipn+1;
            flipt(flipn)=T;
            T=T-log(rand())/fliprate;
        end
        for j=flipn+1:flipn+npipulse
            flipt(j)=(2*j-1)*t(i);
            qpipulse(j)=1;
        end
        [flipt(1:flipn+npipulse),sortindx]=sort(flipt(1:flipn+npipulse));
        qpipulse(1:flipn+npipulse)=qpipulse(sortindx);
        flipt(2:flipn+npipulse+1)=flipt(1:flipn+npipulse);
        flipt(1)=0;
        flipt(flipn+npipulse+2)=2*npipulse*t(i);
        qpipulse(2:flipn+npipulse+1)=qpipulse(1:flipn+npipulse);
        U0=eye(2);
        U1=eye(2);
        qoddpipulse=0; %if have applied odd number of pi pulses
        for j=2:flipn+npipulse+2
            U0=expm(-1i*(flipt(j)-flipt(j-1))*H(:,:,1+qoddpipulse))*U0;
            U1=expm(-1i*(flipt(j)-flipt(j-1))*H(:,:,2-qoddpipulse))*U1;
            if qpipulse(j)
                qoddpipulse=1-qoddpipulse;
            else
                U0=sgmx*U0;
                U1=sgmx*U1;
            end
        end
        x1nuc(i)=real(trace(U0*U1'))/2;
    end
    x0=x0.*x1nuc;
    disp(nuc);
end
