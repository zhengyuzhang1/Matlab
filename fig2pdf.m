function fig2pdf(handle,size,filename)
%save figure with handle 'handle' to a margin-free pdf file with filename
% 'filename' and size 'size'
set(handle,'PaperSize',size);
set(handle,'PaperPositionMode','manual');
set(handle,'PaperPosition',[0 0 size]);
print(handle,'-dpdf',filename);
end