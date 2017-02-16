function copybatch(filename,n)

dotplace = strfind(filename,'.');
for i=1:n
    copyfile(filename,[filename(1:dotplace-1),'_',num2str(i),'.m']);
end
