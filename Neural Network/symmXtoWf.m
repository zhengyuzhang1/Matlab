a = finalX(1,1) * ones(4,1);
b = finalX(2,1) * ones(4,1);
W = TransInv(finalX(3:end,1));
X = [a;b;reshape(W,[],1)];
ss = RBMprob(X,4,4);