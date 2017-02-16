N = 4;
M = 4;
fileID = fopen(['SYK_exact',int2str(N),'_',int2str(M),'.txt'],'r');
A = fscanf( fileID, '%f', N*(N-1)/2 );
J = zeros(N);
k = 0;
for i = 1:N-1
    for j = i+1:N
        k = k + 1;
        J(i,j) = A(k);
    end
end
H = SYK_H(N, J);
[VV,DD]=eig(H);
E = min(diag(DD))/N;
B = fscanf( fileID, '%f');
fclose(fileID);
X = B(1:end/2)+B(end/2+1:end)*1i;
arbm = X(1:N);
crbm = X(N+1:N+M);
wrbm = reshape(X(N+M+1:end),M,N)';
% figure;plot(abs(abs(C' * V).^2));
% fclose(fileID);
% %%
% N = 6; M=64;
% fileID = fopen(['SYK',int2str(N),'_',int2str(M),'.txt'],'r');
% A = fscanf( fileID, '%f' );
% fclose(fileID);
% B = reshape(A(1:end/2), M+N+M*N, M+N+M*N) + reshape(A(end/2+1:end), M+N+M*N, M+N+M*N) * 1i;