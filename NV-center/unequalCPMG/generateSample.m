nrange = [5,7,9];
azx = cell(length(nrange),1);
for nidx = 1:length(nrange)
    n = nrange(nidx);
    azx{nidx} = simNVsam(n,[1,1]);
    while azx{nidx}(1,3)>5*azx{nidx}(2,3)
        azx{nidx} = simNVsam(n,[1,1]);
    end
end
save('nucSample.mat','azx');
    