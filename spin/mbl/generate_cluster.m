function [clust,eachclustn]=generate_cluster(indx,clustsize)
clustn=0;
N=max(indx(:));
clust=zeros(N,clustsize);
eachclustn=zeros(N,1);
inclust=zeros(N,1);
if clustsize==1
    clust=(1:N)';
    eachclustn=ones(N,1);
else
    for i=1:size(indx,1)
        x=indx(i,1);
        y=indx(i,2);
        clx=inclust(x);
        cly=inclust(y);
        if ~clx && ~cly
            clustn=clustn+1;
            clust(clustn,1:2)=[x;y];
            inclust(x)=clustn;
            inclust(y)=clustn;
            eachclustn(clustn)=2;
        elseif clx && ~cly
            if eachclustn(clx)<clustsize
                eachclustn(clx)=eachclustn(clx)+1;
                clust(clx,eachclustn(clx))=y;
                inclust(y)=clx;
            else
                clustn=clustn+1;
                clust(clustn,1)=y;
                inclust(y)=clustn;
                eachclustn(clustn)=1;
            end
        elseif ~clx && cly
            if eachclustn(cly)<clustsize
                eachclustn(cly)=eachclustn(cly)+1;
                clust(cly,eachclustn(cly))=x;
                inclust(x)=cly;
            else
                clustn=clustn+1;
                clust(clustn,1)=x;
                inclust(x)=clustn;
                eachclustn(clustn)=1;
            end
        elseif clx~=cly && eachclustn(clx)+eachclustn(cly)<=clustsize
            clust(clx,eachclustn(clx)+1:eachclustn(clx)+eachclustn(cly))=clust(cly,1:eachclustn(cly));
            inclust(clust(cly,1:eachclustn(cly)))=clx;
            eachclustn(clx)=eachclustn(clx)+eachclustn(cly);
            eachclustn(cly)=0;
            clust(cly,:)=0;
        end
    end
end
clust(~any(clust,2),:)=[];
eachclustn=eachclustn(eachclustn~=0);
end