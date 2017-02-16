filename = mfilename; %Filename is Name+ _ +number.

nN = 200;
nclustsize = 6;
nw = 100000;
ndt = 0.001;
[sortdr,indx,A,r] = generate_spin_sample(nN,0.01);
At=sqrt(A(:,:,1).^2+A(:,:,4).^2+A(:,:,6).^2);
meanA = mean(At(:));
[clust,eachclustn] = generate_cluster(indx,nclustsize);
mz = mbl(nw,ndt/meanA,10000,A,clust,eachclustn);
save(filename)