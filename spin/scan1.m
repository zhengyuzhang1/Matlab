nN=[200];
nclustsize=[5];
nw=[0.1];
ndt=[0.1];
mxA=zeros(size(nN));
A=cell(length(nN),1);
r=cell(length(nN),1);
mz=cell(length(nN),length(nclustsize),length(nw),length(ndt));
for i=1:length(nN)
    [sortdr,indx,A{i},r{i}]=generate_spin_sample(nN(i),0.01);
    mxA(i)=max(A{i}(:));
    for j=1:length(nclustsize)
        [clust,eachclustn]=generate_cluster(indx,nclustsize(j));
        for k=1:length(nw)
            for l=1:length(ndt)
                mz{i,j,k,l}=mean_field_evo(nw(k)*mxA(i),ndt(l)/max([mxA,nw(k)*mxA(i)]),round(10000/ndt(l)),A{i},clust,eachclustn);
                save('scan1.mat')
            end
        end
    end
end
