%global parameter
w = 505*2*pi*0.00107e6;
rng('shuffle');

%control paratemer
maxoptn = 1000;
nrange = [5];
tnrange = [5];
maxtotalt = 1e-4;
minsinglet = 1e-7;

u0 = cell(length(nrange),length(tnrange));
u1 = cell(size(u0));
axes0 = cell(size(u0));
axes1 = cell(size(u0));
angle0 = cell(size(u0));
angle1 = cell(size(u0));

fun = cell(max(nrange),1);
fmin = cell(length(nrange),length(tnrange));
taomin = cell(size(fmin));

load('result5.mat','azx');
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
        fun{i} = @(x) unequaltimeCnot(w,w1{nidx},[x(1:end-1),sum(x(2:2:end-1))-sum(x(1:2:end-1))],i,x(end));
    end
    for tnidx = 1:length(tnrange)
        tn = tnrange(tnidx);
        fmin{nidx,tnidx} = n*ones(n,1);
        taomin{nidx,tnidx} = zeros(n,tn);

        A = ones(1,tn);
        A(2:2:end) = -1;
        A(end) = 0;
        b = -minsinglet;
        nonlc = @(x) nonlcon(x,maxtotalt);
        opts = gaoptimset('display','off','TolCon',minsinglet);
                
        runtimes = 0;
        while max(fmin{nidx,tnidx})>0.001 && runtimes<maxoptn
            for j = 1:n
                [tao,f,existflag] = ga(fun{j},tn,A,b,[],[],[minsinglet*ones(1,tn-1),1],[maxtotalt*ones(1,tn-1),floor(maxtotalt/minsinglet/tn)],nonlc,tn,opts);
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
