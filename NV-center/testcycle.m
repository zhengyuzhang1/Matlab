%load experimentdata1.mat;
srchdata=data;
nuc=[];
update=true;
while update
    update=false;
    for i=1:size(nuc,1)+1    
        %srchdata(:,i)=data./prod([M(w,n,nuc(1:i-1,:),t');M(w,n,nuc(i+1:end,:),t')],1)';
        srchdata(:,i)=rmerrdata(data./prod([M(w,n,nuc(1:i-1,:),t');M(w,n,nuc(i+1:end,:),t')],1)');
        if min(srchdata(:,i))>0.8
            continue;
        end
        azx=srchhypc(w,n,t,srchdata(:,i),data,nuc,i);
        azx=azx(1,:);
        if isempty(nuc)
            minchnuc=1;
        else
            [minchnuc,loc]=min(abs((repmat(azx(1),size(nuc,1),1)-nuc(:,1))./nuc(:,1)));
        end
        if minchnuc>0.05
           loc=size(nuc,1)+1;
        end
        if minchnuc<0.05 || i==size(nuc,1)+1
            nuc(loc,:)=azx;
        end
        if minchnuc>0.005
            update=1;
        end
    end        
    %save(['expt',num2str(size(nuc,1)),'complete.mat'])
end
