function [x,y]=fig2var(handle)
figure_info=findall(handle,'type','line');
x=get(figure_info(1,:),'xdata');
y=get(figure_info(1,:),'ydata');
end
