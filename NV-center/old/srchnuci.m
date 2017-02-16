function [newnuc,newdata]=srchnuci(w,n,data,nuc,t0,deltat,i)
azx=srchhypc(w,n,data(:,i),data(:,1),t0,deltat);
newnuc=nuc;
newnuc(i,:)=azx(1,:);
newdata=data;
t=t0:deltat:t0+(length(data)-1)*deltat;
newdata(:,i+1)=rmerrdata(data(:,1)./prod(M(w,n,newnuc(1:i,:),t),1)');
end