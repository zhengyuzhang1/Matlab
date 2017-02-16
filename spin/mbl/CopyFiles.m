% Copy files

for i=1:40
    copyfile('scanmbl.m',['scanmbl','_',num2str(i),'.m']);
end
