load puresample.mat;
srchdata=data;
nuc=[];
update=true;
while size(nuc,1)<15 && update
    update=false;
    for i=1:size(nuc,1)+1    
        i
        nuc
        srchdata(:,i)=rmerrdata(data./prod([M(w,n,nuc(1:i-1,:),t);M(w,n,nuc(i+1:end,:),t)],1)');
        azx=srchhypcall(w,n,srchdata(:,i),data,t0,deltat,nuc);
        azx=azx(1,:);
        if isempty(nuc)
            minchnuc=1;
        else
            [minchnuc,loc]=min(abs((repmat(azx(1),size(nuc,1),1)-nuc(:,1))./nuc(:,1)));
        end
        if minchnuc>0.05
            nuc(end+1,:)=azx;
        else
            nuc(loc,:)=azx;
        end
        if minchnuc>0.005
            update=1;
        end
    end        
    save(['testcycleall',num2str(size(nuc,1)),'complete.mat']);
    nuc
end
