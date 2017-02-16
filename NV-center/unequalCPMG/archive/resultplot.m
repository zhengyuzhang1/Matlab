avf = zeros(length(nrange),length(tnrange));
for nidx = 1:length(nrange)
    for tnidx = 1:length(tnrange)
        avf(nidx,tnidx) = mean(fmin{nidx,tnidx});
    end
end
figure;surface(nrange,tnrange,avf);