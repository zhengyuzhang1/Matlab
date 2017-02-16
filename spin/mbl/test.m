[maxmz,x0] = max(mz,[],2);
% figure;plot(mz(:,x0(1)));
figure;plot(mz(:,x0(1)));
figure;plot(sum(mz,2));
%% 
initialclust = rem(find(clust==x0(1)),size(clust,1));
[~,ind] = sort(clustr(:,initialclust));
clust = clust(ind,:);
cc = clust';
cc = cc(cc>0);
%% 
figure;plot(mz(1000,cc))
%% 
clustr = zeros(size(clust,1));
for i = 1:size(clust,1)
    for j = i + 1:size(clust,1)
        clustr(i,j) = norm(mean(r(clust(i,1:eachclustn(i)),:)) - mean(r(clust(j,1:eachclustn(j)),:)));  
        clustr(j,i) = clustr(i,j);
    end
end