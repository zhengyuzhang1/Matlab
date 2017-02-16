function f = partrace( m,v,n )
%partrace computes the partial trace of matrix m. v is the list of
%dimensions of subsystems and n is the list of subsystems of interest.
if size(m,1)~=size(m,2) || size(m,1)~=prod(v)
    error('Dimensions not match')
end
if numel(v)==1
    f=m;
    return;
end
dof=prod(v(n)); %dimension of output matrix
f=zeros(dof);
toff=1:numel(v);
toff(n)=[]; %subsystems to trace off
ntraceoff=numel(toff); %number of subsystems to trace off
nob=prod(v(toff)); %number of basis vectors to trace off
traceoffq=false(size(v));
traceoffq(toff)=true; %list of which subsystems should be traced off
ndgridcell=cell(ntraceoff,1);
for i=1:ntraceoff
    ndgridcell{i}=1:v(toff(i));
end
[ndgridcell{1:ntraceoff}]=ndgrid(ndgridcell{:});
for i=1:nob
    basis=1;
    for j=1:numel(v)
       if traceoffq(j)
           tmp=ndgridcell{find(toff==j)};
           basis=kron(basis,sparse(tmp(i),1,1,v(j),1));
       else
           basis=kron(basis,eye(v(j)));
       end
    end
    f=f+basis'*m*basis;
end
end

