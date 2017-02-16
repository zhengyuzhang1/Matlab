N = 9; M = 9; ntry = 3;
idx = [1:40];
prefile = 'fsf.o22400496-';
n = ntry * length(idx);
r = zeros(n,1);
X = zeros(N+M+N*M,n);
for i = 1:length(idx)
    fileID = fopen([prefile,int2str(idx(i))],'r');
    A = fscanf( fileID, '%f',((N+M+N*M)*2+1)*ntry);
    A = reshape(A,(N+M+N*M)*2+1,ntry);
    X(:,(i-1)*ntry+1:(i-1)*ntry+ntry) = A(2:1+(end-1)/2,:) + A(1+(end-1)/2+1:end,:)*1i;
    r((i-1)*ntry+1:(i-1)*ntry+ntry) = A(1,:)';
    fclose(fileID);
end
figure;plot(r);
%%
[minr,ind] = min(r);
ind = ind +1+2;
gstate = RBMprob(X(:,ind),N,M);
H = Hfsf2d(3,1,1);
[V,D] = eig(H);
overlap = abs(gstate' * V).^2;
figure;plot(overlap);