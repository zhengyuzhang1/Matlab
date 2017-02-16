filename = mfilename; %Filename is Name+ _ +number.

nN = 200;
ndt = 1e5;
step = 10000;
[~,~,A,r] = generate_spin_sample(nN,0.01);
[state,up,H] = mblexact(ndt,step,A);
fliprate = A(up,:,1)+A(up,:,4);
[maxA,nearnb] = max(fliprate);
t = ndt / max(H(up,:)) * (0:step-1);
p = abs(state).^2;
figure;plot(t*maxA,p(up,:));
xlabel('time*Hii');
ylabel('pi');
figure;plot(mean(p,2));
xlabel('spins');
ylabel('<p>');
figure;plot(fliprate);
xlabel('spins');
ylabel('Hij');
c = mean(p,2) / max(mean(p,2)) * 10;
figure;scatter(r(:,1),r(:,2),[],c);
%caz = diag(H) / max(diag(H)) * 10;
%figure;scatter(r(:,1),r(:,2),[],caz);
figure;plot(diag(H));
xlabel('spins');
ylabel('Hii');
%% 
fig2pdf(1,[5,5],'pi_time_1');
fig2pdf(2,[5,5],'meanp_1');
fig2pdf(3,[5,5],'hop_1');
fig2pdf(4,[5,5],'meanp2D_1');
%% 
[a,b]=sort(mean(p,2),'descend');
[c,d]=sort(abs(fliprate'),'descend');
[e,f]=sort(abs(H(up,up)-diag(H)));
%% 
hoppingl = 10;
hopping = zeros(size(A,1),size(A,1),hoppingl);
%% 
distance = zeros(size(A,1));
for i = 1:size(distance,1)
    for j = i + 1:size(distance,1)
        distance(i,j) = norm(r(i,:) - r(j,:));
        distance(j,i) = distance(i,j);
    end
end
