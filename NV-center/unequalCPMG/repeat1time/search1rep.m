%global parameter
w = 505*2*pi*0.00107; %unit is MHz
rng('shuffle');
rngState = rng;
rng(rngState.Seed+uint32(feature('getpid')));

%control paratemer
maxoptn = 500;
nrange = [5];
tnrange = 5:2:21;
maxtotalt = 100; %time unit is ms
minsinglet = 0;

u0 = cell(length(nrange),length(tnrange));
u1 = cell(size(u0));
axes0 = cell(size(u0));
axes1 = cell(size(u0));
angle0 = cell(size(u0));
angle1 = cell(size(u0));

fun = cell(max(nrange),1);
fmin = cell(length(nrange),length(tnrange));
taomin = cell(size(fmin));

load('nucSample.mat','azx');
% azx = cell(length(nrange),1);
w1 = cell(length(nrange),1);
for nidx = 1:length(nrange)
    n = nrange(nidx);
%     azx{nidx} = simNVsam(1,[],n);
%     while min(azx{nidx}(:))<1 || min(min(abs(diff(azx{nidx},1))))<1
%         azx{nidx} = simNVsam(1,[],n);
%     end
    w1{nidx} = [sqrt((w+azx{nidx}(:,1)).^2+azx{nidx}(:,2).^2),...
        (w+azx{nidx}(:,1))./sqrt((w+azx{nidx}(:,1)).^2+azx{nidx}(:,2).^2),...
        azx{nidx}(:,2)./sqrt((w+azx{nidx}(:,1)).^2+azx{nidx}(:,2).^2)];
    for i = 1:max(nrange)
        fun{i} = @(x) unequaltimeCnot(w,w1{nidx},x,i,1);
    end
    for tnidx = 1:length(tnrange)
        tn = tnrange(tnidx);
        fmin{nidx,tnidx} = n*ones(n,1);
        taomin{nidx,tnidx} = zeros(n,tn);

        Aeq = ones(1,tn);
        Aeq(2:2:end) = -1;
        beq = 0;
        A = ones(1,tn);
        b = maxtotalt;
        lb = minsinglet * ones(1,tn);
        opts = optimoptions(@fmincon,'Display','off','MaxFunctionEvaluations',5000);
                        
        runtimes = 0;
        while max(fmin{nidx,tnidx})>0.001 && runtimes<maxoptn
            t0 = rand(1) * maxtotalt/2 * diff([[0,sort(rand(1,(tn-1)/2)),1];[0,sort(rand(1,(tn-3)/2)),1,0]],1,2);
            t0 = t0(1:end-1);
            for j = 1:n
                [tao,f,existflag] = fmincon(fun{j},t0,A,b,Aeq,beq,lb,[],[],opts);
                if f<fmin{nidx,tnidx}(j) && existflag==1
                    fmin{nidx,tnidx}(j) = f;
                    taomin{nidx,tnidx}(j,:) = tao;
                end
            end
            runtimes = runtimes + 1;
        end
        
        u0{nidx,tnidx} = zeros(2,2,n,n);%third index is rotating qubit, fouth index is target qubit
        u1{nidx,tnidx} = zeros(2,2,n,n);
        axes0{nidx,tnidx} = zeros(3,n,n);
        axes1{nidx,tnidx} = zeros(3,n,n);
        angle0{nidx,tnidx} = zeros(n,n);
        angle1{nidx,tnidx} = zeros(n,n);
        for k = 1:n
            [~,b,c]=fun{k}(taomin{nidx,tnidx}(k,:));
            [~,u0{nidx,tnidx}(:,:,:,k),u1{nidx,tnidx}(:,:,:,k)] = fun{k}(taomin{nidx,tnidx}(k,:));
            for j = 1:n
                [axes0{nidx,tnidx}(:,j,k),angle0{nidx,tnidx}(j,k)] = axisangle(u0{nidx,tnidx}(:,:,j,k));
                [axes1{nidx,tnidx}(:,j,k),angle1{nidx,tnidx}(j,k)] = axisangle(u1{nidx,tnidx}(:,:,j,k));
            end
        end

    end
end
save(mfilename);
