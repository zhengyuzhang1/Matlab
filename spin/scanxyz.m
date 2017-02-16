nN=[200,400,600];
nclustsize=[4,6,8];
nw=[0.1,1,10];
ndt=[1,0.1];
mxA=zeros(size(nN));
A=cell(length(nN),1);
r=cell(length(nN),1);
mx=cell(length(nN),length(nclustsize),length(nw),length(ndt));
my=cell(length(nN),length(nclustsize),length(nw),length(ndt));
mz=cell(length(nN),length(nclustsize),length(nw),length(ndt));
for i=1:1%length(nN)
    [sortdr,indx,A{i},r{i}]=generate_spin_sample(nN(i),0.01);
    mxA(i)=max(A{i}(:));
    for j=1:1%length(nclustsize)
        [clust,eachclustn]=generate_cluster(indx,nclustsize(j));
        for k=2:2%length(nw)
            for l=1:1%length(ndt)
                [mx{i,j,k,l},my{i,j,k,l},mz{i,j,k,l}]=mean_field_evo(nw(k)*mxA(i),ndt(l)/max([mxA(i),nw(k)*mxA(i)]),round(1000/ndt(l)),A{i},clust,eachclustn);
                %save('scan.mat')
            end
        end
    end
end
