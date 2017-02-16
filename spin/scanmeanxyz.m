nN = 200;
nclustsize = 6;
nw = 1000;
ndt = 0.001;
[sortdr,indx,A,r] = generate_spin_sample(nN,0.01);
At=sqrt(A(:,:,1).^2+A(:,:,4).^2+A(:,:,6).^2);
meanA = mean(At(:));
[clust,eachclustn] = generate_cluster(indx,nclustsize);
[mx,my,mz] = mean_field_evo(nw,ndt/meanA,1000,A,clust,eachclustn);
save('scanmeanxyz.mat')