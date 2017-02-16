nof = 30;
n = 2;
tn = 15;
taominall = cell(nof,1);
fminall = zeros(n,4,nof);
for k = 1:nof
    taominall{k} = zeros(n,tn,4);
end
for i = 1:nof
    load(['experiments_',num2str(i),'.mat'],'fmin','taomin');
    if i == 1
        fminall(:,:,1) = fmin;
        taominall{1} = taomin;
    else
        for k = 1:n
            for l = 1:4
                if fmin(k,l)<fminall(k,l,i-1)
                    fminall(k,l,i) = fmin(k,l);
                    taominall{i}(k,:,l) = taomin(k,:,l);
                else
                    fminall(k,l,i) = fminall(k,l,i-1);
                    taominall{i}(k,:,l) = taominall{i-1}(k,:,l);
                end
            end
        end
    end                        
end
%% total gate time
for k = 1:nof
    taominall{k}(:,end+1,:) = sum(taominall{k},2);
end
%% experimental fideltiy of preparing 01+10 state
w = 495*2*pi*0.00107;
azx = 2*pi* [0.07919, 0.10282; -0.07089,0.05594];
w1 = [sqrt((w+azx(:,1)).^2+azx(:,2).^2),(w+azx(:,1))./sqrt((w+azx(:,1)).^2+azx(:,2).^2),...
    azx(:,2)./sqrt((w+azx(:,1)).^2+azx(:,2).^2)];

expT = taominall{end};
rou0 = kron(0.5*eye(2),kron(0.5*eye(2),0.5*eye(2)));

rou = kron([1 0 ; 0 0], partrace(rou0, [2 2 2], [2 3]));
U = kron(rSpin([0 1 0]*pi/2), eye(4));
U = UGate(w, w1, expT(2,1:end-1,1))*U;
U = kron(rSpin([1 0 0]*pi/2),eye(4))*UGate(w, w1, expT(2,1:end-1,3))*U;
U = UGate(w, w1, expT(2,1:end-1,1))*U;
rou = U * rou * U';
rou = kron([1 0 ; 0 0], partrace(rou, [2 2 2], [2 3]));
U = kron(rSpin([1 0 0]*pi/2), eye(4));
U = UGate(w, w1, expT(2,1:end-1,1))*U;
U = UGate(w, w1, expT(1,1:end-1,1))*U;
U = kron(rSpin([0 1 0]*pi), eye(4))*kron(rSpin([1 0 0]*pi/2), eye(4))*U;
U = UGate(w, w1, expT(1,1:end-1,4))*U;
U = UGate(w, w1, expT(1,1:end-1,1))*U;
U = kron(rSpin([0 1 0]*pi/2), eye(4))*kron(rSpin([1 0 0]*pi/2), eye(4))*U;
U = UGate(w, w1, expT(1,1:end-1,2))*U;
U = UGate(w, w1, expT(1,1:end-1,3))*U;
U = UGate(w, w1, expT(1,1:end-1,1))*U;
rou = U * rou * U';
rou = kron([1 0 ; 0 0], partrace(rou, [2 2 2], [2 3]));
phi = 1/2 * [1; 1i; 1i; 1];
fid = sqrt(phi'*rou(1:4,1:4)*phi);
totaltime = 3*expT(1,end,1)+3*expT(2,end,1)+expT(2,end,3)+expT(1,end,2)+expT(1,end,3)+expT(1,end,4);