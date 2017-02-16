%global parameter
w = 505*2*pi*0.00107e6;

%control paratemer
rng('shuffle');
maxoptn = 80;
nrange = [2,4,6,8,10];
tnrange = [3,5,7,9,11];
ub = 1e-5;

fun = cell(max(nrange),1);
fmin = cell(length(nrange),length(tnrange));
taomin = cell(size(fmin));

u0 = cell(length(nrange),length(tnrange));
u1 = cell(size(u0));
axes0 = cell(size(u0));
axes1 = cell(size(u0));
angle0 = cell(size(u0));
angle1 = cell(size(u0));

azx = cell(length(nrange),1);
for nidx = 1:length(nrange)
    n = nrange(nidx);
    azx{nidx} = simNVsam(1,[],n);
    while ~isempty(find(azx{nidx}(:,2)==0,1))
        azx{nidx} = simNVsam(1,[],n);
    end
    for i = 1:max(nrange)
        fun{i} = @(x) unequaltimeCPMG(w,azx{nidx},x,i);
    end
    for tnidx = 1:length(tnrange)
        tn = tnrange(tnidx);
        fmin{nidx,tnidx} = n*ones(n,1);
        taomin{nidx,tnidx} = zeros(n,tn);

        Aeq = ones(1,tn);
        Aeq(1:2:end) = -1;
        beq = 0;
        opts = optimoptions(@fmincon,'display','off');
                
        i = 0;
        while max(fmin{nidx,tnidx})>0.01 && i<maxoptn
            tao0 =  ub * rand(1,tn);
            tao0(end) = sum(tao0(2:2:end-1))-sum(tao0(1:2:end-2));
            if tao0(end) < 0
                tao0(end) = -tao0(end);
                tmp = tao0(1:2:end-2);
                tao0(1:2:end-2) = tao0(2:2:end-1);
                tao0(2:2:end-1) = tmp;
            end
            for j = 1:n
                [tao,f] = fmincon(fun{j},tao0,[],[],Aeq,beq,zeros(1,tn),ub*ones(1,tn),[],opts);
                if f<fmin{nidx,tnidx}(j)
                    fmin{nidx,tnidx}(j) = f;
                    taomin{nidx,tnidx}(j,:) = tao;
                end
            end
            i = i + 1;
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
        disp([nidx,tnidx]);
    end
end
