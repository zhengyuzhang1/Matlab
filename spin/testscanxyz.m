nN=[200,400,600];
nclustsize=[4,6,8];
ndt=[1,0.1];
mxA=zeros(size(nN));
A=cell(length(nN),1);
r=cell(length(nN),1);
mx=cell(length(nN),length(nclustsize),length(ndt));
my=cell(length(nN),length(nclustsize),length(ndt));
mz=cell(length(nN),length(nclustsize),length(ndt));
for i=1:1%length(nN)
    [sortdr,indx,A{i},r{i}]=generate_spin_sample(nN(i),0.01);
    mxA(i)=max(A{i}(:));
    for j=1:1%length(nclustsize)
        [clust,eachclustn]=generate_cluster(indx,nclustsize(j));
        for l=1:1%length(ndt)
            [mx{i,j,l},my{i,j,l},mz{i,j,l}]=mean_field_evo(0,ndt(l)/mxA(i),round(1000/ndt(l)),A{i},clust,eachclustn);
            %save('scan.mat')
        end
    end
end