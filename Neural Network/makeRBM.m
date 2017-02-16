n = 4;
N = 2^n;
V = rand(N,1) + rand(N,1)*1i;
p = 0.5;%ratio of non-zero states
V(rand(N,1)>p) = 0;
V = VV(:,1);
V(abs(V)<1e-6) = 0;
V = V / sum(abs(V));
M = sum(abs(V)>0);

a = zeros(M,1);
a(:) = 50;
lambda = zeros(M,1);
lambda(1) = 10;
w = zeros(n,M);
c = zeros(M,1);
v = zeros(n,M);
[~, indx] = sort(abs(V),'descend');
for i = 1:M
    v(:,i) = flip(de2bi(indx(M+1-i)-1,n));
    w(:,i) = a(i) * (v(:,i) - 1/2);
end
c(1) = -w(:,1)' * v(:,1) + lambda(1);

p1 = (1 + exp(lambda(1))) * 2^(-n) / (1 + exp(lambda(1))*2^(-n));
p2 = (1-p1)/(2^n-1);
for i = 2:M
%     display(RBMprob(w(:,1:i-1),c(1:i-1)));
    lambda(i) = log(V(indx(M+1-i))/V(indx(M+1-i+1))*p1/p2-1);
    temp = p1;
    p1 = (1+exp(lambda(i)))*p2 / (1+exp(lambda(i))*p2);
    p2 = p2 / (1+exp(lambda(i))*p2);
    c(i) = -w(:,i)' * v(:,i) + lambda(i);
end  

%%
vm = zeros(N,1);
for i = 1:N
    state = flip(de2bi(i-1,n)');
    alpha = 0;
    theta = c + w' * state;
    vm(i) = exp(Phi_log_bin(alpha,theta));
end
vm = vm / sum(abs(vm));
eps = abs(V'*vm)/sqrt(V'*V)/sqrt(vm'*vm);
%%
ap = 1/4 * sum(w,2);
cp = c/2 + 1/4 * sum(w,1)';
wp = w/4;
vmpx = zeros(N,1);
for i = 1:N
    state = flip(de2bi(i-1,n)');
    state(state==0) = -1;
    alpha = ap' * state;
    theta = cp + wp' * state;
    vmpx(i) = Phi_log(alpha,theta);
end
vmpx = vmpx - min(real(vmpx));
vmp = exp(vmpx);
vmp = vmp / sum(abs(vmp));
epsp = abs(V'*vmp)/sqrt(V'*V)/sqrt(vmp'*vmp);
