N = 200;
a = rand(N);
tic 
shiftindex = 0:N-1;
[m,n]=size(a);
S=full(sparse(mod(shiftindex,m)+1,1:n,1,m,n));
a_new=ifft(fft(a).*fft(S),'symmetric');
toc
tic
for i = 2:N
    a(:,i) = circshift(a(:,i),i-1);
end
toc
max(a(:)-a_new(:))