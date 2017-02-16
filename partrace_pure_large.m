function Xpt = partrace_pure_large( phi,sys,dim )
%partrace_pure_large computes the partial trace of pure state phi which is 
%so large that the full density matrix cannot be stored in memory. 
%dim is the list of dimensions of subsystems and sys is the list of subsystems to trace on.

n_all = length(dim);
n_sys = length(sys);
dim_all = prod(dim);
dim_sys = prod(dim(sys));
dim_rem = dim_all / dim_sys;
Xpt = zeros(dim_rem);

perm = [setdiff(1:n_all,sys),sys];
phi_p = PermuteSystems(phi,perm,dim);

for i = 1:dim_rem
    for j = i:dim_rem
        Xpt(i,j) = phi_p((i-1)*dim_sys+1:i*dim_sys)' * phi_p((j-1)*dim_sys+1:j*dim_sys);
    end
end
Xpt = Xpt + triu(Xpt,1)';        

end

