% Copy files

for i=1:40
    copyfile('dipoledynamic.m',['dipoledynamic','_',num2str(i),'.m']);
end
