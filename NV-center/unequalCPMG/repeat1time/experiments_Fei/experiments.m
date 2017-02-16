%global parameter
w = 495*2*pi*0.00107; %unit is MRad
rng('shuffle');
rngState = rng;
rng(rngState.Seed+uint32(feature('getpid')));

%control paratemer
maxoptn = 1000;
n = 2;
tn = 15;
maxtotalt = 40; %time unit is ms
minsinglet = 0;

azx = 2*pi* [0.07919, 0.10282; -0.07089,0.05594];
w1 = [sqrt((w+azx(:,1)).^2+azx(:,2).^2),(w+azx(:,1))./sqrt((w+azx(:,1)).^2+azx(:,2).^2),...
    azx(:,2)./sqrt((w+azx(:,1)).^2+azx(:,2).^2)];

fun = cell(n,4);
X = [0 0 0 1 ; 1 0 0 0 ; 0 1 0 0 ; 0 0 1 0];
Z = [1 0 0 0 ; 0 1i 0 0 ; 0 0 -1 0 ; 0 0 0 -1i];
Ub = zeros(4,4,16);
for i = 0:3
    for j = 0:3
        Ub(:,:,1+4*i+j) = X^i * Z^j;
    end
end
Ub = Ub(:,:,2:end);

gate = zeros(4,4,4);
gate(:,:,1) = 1/sqrt(2) * blkdiag([1 -1i ; -1i 1], [1 1i ; 1i 1]); %conditional X_pi/2
gate(:,:,2) = 1/sqrt(2) * blkdiag([1 -1i ; -1i 1], [1 -1i ; -1i 1]); %unconditional X_pi/2
gate(:,:,3) = 1/sqrt(2) * blkdiag([1-1i 0 ; 0 1+1i], [1-1i 0 ; 0 1+1i]); %unconditional Z_pi/2
gate(:,:,4) = 1/sqrt(2) * blkdiag([1+1i 0 ; 0 1-1i], [1+1i 0 ; 0 1-1i]); %unconditional Z_-pi/2

for i = 1:n
    Ut = repmat(eye(4),1,1,n);
    for j = 1:4     
        Ut(:,:,i) = gate(:,:,j);
        fun{i,j} = @(x) infidGate(w,w1,x,Ut,Ub);
    end
end
fmin = n*ones(n,4); %we need 4 gates for each nucleus
taomin = zeros(n,tn,4);

Aeq = ones(1,tn);
Aeq(2:2:end) = -1;
beq = 0;
A = ones(1,tn);
b = maxtotalt;
lb = minsinglet * ones(1,tn);
opts = optimoptions(@fmincon,'Display','off','MaxFunctionEvaluations',5000);

runtimes = 0;
while max(fmin(:))>0.001 && runtimes<maxoptn
    t0 = rand(1) * maxtotalt/2 * diff([[0,sort(rand(1,(tn-1)/2)),1];[0,sort(rand(1,(tn-3)/2)),1,0]],1,2);
    t0 = t0(1:end-1);
    for j = 1:n
        for k = 1:4
            [tao,f,existflag] = fmincon(fun{j,k},t0,A,b,Aeq,beq,lb,[],[],opts);
            if f<fmin(j,k) && existflag==1
                fmin(j,k) = f;
                taomin(j,:,k) = tao;
            end
        end    
    end
    runtimes = runtimes + 1;
end

save(mfilename);
