function f = srchhypc(w, n, t, data, data0, nuc,i)
[dps,locs]=findpeaks(-data);
dps=-dps;
alldips=[locs,dps];
mid=4;lp=4;rp=4;mindv=dps(2);
gooddip=false;
gtu=true(size(data));
method=1;
if i==size(nuc,1)+1
    fig=figure(i);plot(t,data);
    method=input( 'Which method to use:');
end
if method==1
    while ~gooddip
        while ~gtu(mid) || ((find(dps==mindv)>1 && 1-mindv<0.5*(1-dps(find(dps==mindv)-1))) || (find(dps==mindv)<size(dps,1) && (1-mindv)<0.5*(1-dps(find(dps==mindv)+1))) || (data(mid)>0 && (data(lp)<0.95 || data(rp)<0.95)) || mid-lp<2 || rp-mid<2 || abs(data(mid)-data(lp))/abs(data(mid)-data(rp))<0.5 || abs(data(mid)-data(lp))/abs(data(mid)-data(rp))>2 || (abs(data(mid)-data(mid-2))<0.015 && abs(data(mid)-data(mid+2))<0.015) || min(abs(data(mid)-[data(mid+3),data(mid-3)]))<0.02  || isnan(abs(data(mid)-data(lp))/abs(data(mid)-data(rp))))
           [mindv,mindloc]=min(alldips(:,2));
            mid=alldips(mindloc,1);
            lp=mid;rp=mid;
            while lp>1 && data(lp-1)>data(lp)
                lp=lp-1;
            end
            while rp<length(data) && data(rp+1)>data(rp)
                rp=rp+1;
            end
            alldips(mindloc,:)=[];
        end
        dips=[locs(1-dps>0.3*(1-mindv)),dps(1-dps>0.3*(1-mindv))];
        [prds,gtu(mid)]=srchprd(w,t,data,find(dips(:,2)==mindv),dips(:,1),5);
        gooddip=gtu(mid);
    end
elseif method==2
    dcm_obj=datacursormode(fig);
    set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on');
    c_info=getCursorInfo(dcm_obj);
    prds(:,4)=input('Enter prds:');
    for c_index=1:length(c_info)
        prds(c_index,1)=c_info(c_index).DataIndex;
        prds(c_index,2)=prds(c_index,1);
        prds(c_index,3)=prds(c_index,1);
        while  prds(c_index,2)>1 && data( prds(c_index,2)-1)>data( prds(c_index,2))
             prds(c_index,2)= prds(c_index,2)-1;
        end
        while prds(c_index,3)<length(data) && data(prds(c_index,3)+1)>data(prds(c_index,3))
            prds(c_index,3)=prds(c_index,3)+1;
        end
    end
end
prds=[sortrows(prds(:,1:3),1),sort(prds(:,4))];
azx=fitdip(w,n,t,data,prds);
azx=azx(1:2);
%{
gos=1;
range=abs(azx./[10,3]);
data1=data0;
if azx(2)<5e5
    range=abs(azx).*[0.1,0.3];
    data1=data0;
    [~,gos]=fsrch(w,n,azx,range,range/10,data1,t);
end
if ~isempty(nuc) && gos~=0
    [minchnuc,loc]=min(max(abs((repmat(azx,size(nuc,1),1)-nuc)./nuc),[],2));
else
    minchnuc=1;
end
if minchnuc<0.05
    range=3*abs(nuc(loc,:)-azx);
end
step=min(range/2,max(min(2000*pi*[1,1],range/30),10*2*pi*[1,1]));
azx0=azx;
while gos>0 && min(step)>2*pi
    [azx,gos]=fsrch(w,n,azx,range,step,data1,t);
    range=range/2;
    step=step/2;
end
if gos>10
    azx=azx0;
end
%}
%[azx,gofs]=sfsrch(w,n,azx,range,step,data,t0,deltat,tn);
gos=0;
f=[azx;gos,0];    
end