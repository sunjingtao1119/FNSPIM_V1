close all
load T.mat
load T3.mat
figure(1)
spy(T,'k')
title('Matrix T','FontName','Times','FontWeight','bold','FontSize',12)
set(gca,'FontName','Times','FontWeight','bold','FontSize',12)
pos=get(gcf,'Position');
pos(3)=pos(4); 
set(gcf,'Position',pos)
figure(2)
spy(T3,'k')
title('Matrix T3 ','FontName','Times','FontWeight','bold','FontSize',12)
set(gca,'FontName','Times','FontWeight','bold','FontSize',12)
pos(3)=pos(4); 
set(gcf,'Position',pos)
