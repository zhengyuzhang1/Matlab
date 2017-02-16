function sample=simNVsam(n,filter,p)
%p is the C13 abundance, usually 1.1%. n is the number of largest nuclei kept.
%Coupling unit is MRad.
%filter = [ifNoDuplicate,ifNoZero],default value is [0,0]

rng('shuffle');
rngState = rng;
rng(rngState.Seed+uint32(feature('getpid')));

if nargin < 2 || isempty(filter)
    filter = [0,0];
end
if nargin < 3 || isempty(p)
    p = 1.1e-2;
end    
l=10;
[xgrid,ygrid,zgrid]=meshgrid(-l:l,-l:l,-l:l);
xgrid=bsxfun(@plus,[0,0.25,0,0.25,0.5,0.75,0.5,0.75],xgrid(:));
ygrid=bsxfun(@plus,[0,0.25,0.5,0.75,0,0.25,0.5,0.75],ygrid(:));
zgrid=bsxfun(@plus,[0,0.25,0.5,0.75,0.5,0.75,0,0.25],zgrid(:));
lat=[xgrid(:),ygrid(:),zgrid(:)];
lat([l+1,3*l+2],:)=[];
mask=rand(8*(2*l+1)^3-2,1)<p;
lat=lat(mask,:);
nz=1/sqrt(3)*[1;1;1];
nx=normc([1;-1;-1]-[1,-1,-1]*nz*nz);
ny=cross(nz,nx);
lat=lat*[nx,ny,nz];
azfun=@(unitrz,r) (3*unitrz.^2-1)./realpow(r,3);
axfun=@(unitrz,r) 3*abs(unitrz).*realsqrt(1-realpow(unitrz,2))./realpow(r,3);
%coefficient
miu0=4e-7*pi;
gamma_e=1.76e11;
gamma_n=6.73e7;
a0=0.35668e-9;
h=1.054e-34;
coef=miu0/4/pi*gamma_e*gamma_n*h/a0^3;
r=sqrt(sum(lat.^2,2));
unitrz=lat(:,3)./r;
sample=coef*[azfun(unitrz,r),axfun(unitrz,r)];
sample(:,3)=sqrt(sample(:,1).^2+sample(:,2).^2);
[sample,order]=sortrows(sample,-3);
lat=lat(order,:);
sample=sortrows(sample,-3);
sample=sample * 1e-6;
%apply filter
if filter(1)
    sample(any(abs(diff(sample,1,1))<1e-3,2),:) = [];
end
if filter(2)
    sample(any(abs(sample)<1e-3,2),:) = [];
end
sample = sample(1:n,:);