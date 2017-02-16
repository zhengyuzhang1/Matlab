load puresample.mat;
srchdata=data;
nuc=[];
update=true;
while size(nuc,1)<10 && update
    update=false;
    for i=1:size(nuc,1)+1            
        srchdata(:,i)=rmerrdata(data(:,1)./prod([M(w,n,nuc(1:i-1,:),t);M(w,n,nuc(i+1:end,:),t)],1)');
        azx=srchhypc(w,n,srchdata(:,i),data,t0,deltat,nuc);
        azx=azx(1,:);
        if isempty(nuc)
            minchnuc=1;
        else
            [minchnuc,loc]=min(max(abs((repmat(azx,size(nuc,1),1)-nuc)./nuc),[],2));
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
    save(['testcycle',num2str(size(nuc,1)),'complete.mat']);
    nuc
end
