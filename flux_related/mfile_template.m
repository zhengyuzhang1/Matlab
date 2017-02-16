filename = mfilename; %Filename is Name+ _ +number.
filenameInd = strfind(filename,'_');
FileIndex = str2double(filename(filenameInd+1:end)); %extract the number from the file.


save(filename)