function [meanmx,meanmy,meanmz]=mean_field_evo(w,dt,step,A,clust,eachclustn)
clustn=size(clust,1);
clustsize=size(clust,2);
H0=zeros(2^clustsize,2^clustsize,clustn)+0j;
sx=zeros(2^clustsize,2^clustsize,clustsize);
sy=zeros(2^clustsize,2^clustsize,clustsize);
sz=zeros(2^clustsize,2^clustsize,clustsize);
Sx=1/2*[0,1;1,0];
Sy=1/2*[0,-1j;1j,0];
Sz=1/2*[1,0;0,-1];
for i=1:clustsize
    sx(:,:,i)=kron(eye(2^(i-1)),kron(Sx,eye(2^(clustsize-i))));
    sy(:,:,i)=kron(eye(2^(i-1)),kron(Sy,eye(2^(clustsize-i))));
    sz(:,:,i)=kron(eye(2^(i-1)),kron(Sz,eye(2^(clustsize-i))));
end
sxforclust=zeros(2^clustsize,2^clustsize,clustn);
syforclust=zeros(2^clustsize,2^clustsize,clustn);
szforclust=zeros(2^clustsize,2^clustsize,clustn);
for i=1:clustn
    for j=1:eachclustn(i)
        szforclust(:,:,i)=szforclust(:,:,i)+sz(:,:,j);
        sxforclust(:,:,i)=sxforclust(:,:,i)+sx(:,:,j);
        syforclust(:,:,i)=syforclust(:,:,i)+sy(:,:,j);
        H0(:,:,i)=H0(:,:,i)-sz(:,:,j)*w;
        for k=j+1:eachclustn(i)
            x=clust(i,j);
            y=clust(i,k);
            axx=A(x,y,1);
            axy=A(x,y,2);
            axz=A(x,y,3);
            ayy=A(x,y,4);
            ayz=A(x,y,5);
            azz=A(x,y,6);
            H0(:,:,i)=H0(:,:,i)+(axx*sx(:,:,j)*sx(:,:,k)+ayy*sy(:,:,j)*sy(:,:,k)+azz*sz(:,:,j)*sz(:,:,k)...
                +axy/2*(sx(:,:,j)*sy(:,:,k)+sy(:,:,j)*sx(:,:,k))+ayz/2*(sy(:,:,j)*sz(:,:,k)+sz(:,:,j)*sy(:,:,k))+axz/2*(sz(:,:,j)*sx(:,:,k)+sx(:,:,j)*sz(:,:,k)));
        end
    end
end
H=zeros(size(H0))+0j;

state=zeros(2^clustsize,clustn,step);

%initial state: all X
% for i=1:clustn
%     tmpstate=1;
%     for j=1:eachclustn(i)
%         tmpstate=kron(tmpstate,1/sqrt(2)*[1;1]);
%     end
%     for j=eachclustn(i)+1:clustsize
%         tmpstate=kron(tmpstate,[1;0]);
%     end
%     state(:,i,1)=tmpstate;
% end

%initial state: random one up, all others down
rand1=randi([1,clustn],1);
rand2=randi([1,eachclustn(i)],1);
for i=[1:rand1-1,rand1+1:clustn]
    tmpstate=1;
    for j=1:eachclustn(i)
        tmpstate=kron(tmpstate,[0;1]);
    end
    for j=eachclustn(i)+1:clustsize
        tmpstate=kron(tmpstate,[1;0]);
    end
    state(:,i,1)=tmpstate;
end
tmpstate=1;
for j=1:rand2-1
    tmpstate=kron(tmpstate,[0;1]);
end
tmpstate=kron(tmpstate,[1;0]);
for j=rand2+1:eachclustn(rand1)
    tmpstate=kron(tmpstate,[0;1]);
end
state(:,rand1,1)=tmpstate;
    
    
means=zeros([size(clust),3]);
mx=zeros(clustn,step);
my=zeros(clustn,step);
mz=zeros(clustn,step);
for t=2:step
    for i=1:clustn
        for j=1:eachclustn(i)
            means(i,j,1)=state(:,i,t-1)'*sx(:,:,j)*state(:,i,t-1);
            means(i,j,2)=state(:,i,t-1)'*sy(:,:,j)*state(:,i,t-1);
            means(i,j,3)=state(:,i,t-1)'*sz(:,:,j)*state(:,i,t-1);
        end
    end
    for i=1:clustn
        for j=1:eachclustn(i)
            x=clust(i,j);
            ax=0;ay=0;az=0;
            for k=[1:i-1,i+1:clustn]
                for l=1:eachclustn(k)
                    y=clust(k,l);
                    ax=A(x,y,1)*means(k,l,1)+A(x,y,2)*means(k,l,2)+A(x,y,3)*means(k,l,3);
                    ay=A(x,y,2)*means(k,l,1)+A(x,y,4)*means(k,l,2)+A(x,y,5)*means(k,l,3);
                    az=A(x,y,3)*means(k,l,1)+A(x,y,5)*means(k,l,2)+A(x,y,6)*means(k,l,3);
                end
            end
            H(:,:,i)=H0(:,:,i)+ax*sx(:,:,j)+ay*sy(:,:,j)+az*sz(:,:,j);
        end
        state(:,i,t)=expm(-1j*dt*H(:,:,i))*state(:,i,t-1);
        mx(i,t)=state(:,i,t)'*sxforclust(:,:,i)*state(:,i,t);
        my(i,t)=state(:,i,t)'*syforclust(:,:,i)*state(:,i,t);
        mz(i,t)=state(:,i,t)'*szforclust(:,:,i)*state(:,i,t);
    end
end
meanmx=real(sum(mx,1))/max(clust(:));
meanmy=real(sum(my,1))/max(clust(:));
meanmz=real(sum(mz,1))/max(clust(:));
end
       