
Data = scan_result;
Loop = length(Data) / 2;
Sig=[];
Ref=[];
Ratio=[];
Norm_Ratio = [];
Pointsize = 5;
ini_precess = 100;
for i = 1:Loop
   Sig(i) = Data(2*i-1);
   Ref(i) = Data(2*i);
   
   Ratio(i) = Sig(i)/Ref(i);
end
max_r = 1.16;
min_r = 0.975;

for i = 1:Loop
   
   
   Norm_Ratio(i) = (Ratio(i) - min_r)/(max_r - min_r);
end

figure
subplot(3,1,1);
x = 1:Loop;
time = Pointsize*(x) + ini_precess;
plot(time, Norm_Ratio(x));
%figure
% subplot(3,1,2);
% plot(time, Sig(x));
%figure
% subplot(3,1,3)
% plot(time, Ref(x));