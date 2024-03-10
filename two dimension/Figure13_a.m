PRK_4_MS_figure13_a
figure(1)
t=500*dt;
load MSPU500.mat
trisurf(Top, p(:,1), p(:,2), u4)
colorbar
colormap(gca,[seiscolor])
clim([-1e-8,1e-8])
shading interp
ylabel('Depth(m)')
xlabel('Distance(m)')
str="t="+num2str(t)+"s"+"(SPRK4)";
title(str,'FontName','Times','FontWeight','bold','FontSize',12)
%    colormap gray
set(gca,'FontName','Times','FontWeight','bold','FontSize',12)
 view([0,-90]);
 xlim([0,lex])
ylim([0,ley])
 dx=lex*5/1000;
 dy=ley*5/1000;
 set(gcf,'unit','centimeters','position',[10 5 dx dy]);
% figure(2)
% t=500*dt;
% load U500.mat; 
% trisurf(Top, p(:,1), p(:,2), U)
% colorbar
% colormap(gca,[seiscolor])
% clim([-1e-8,1e-8])
% shading interp
% xlabel('Depth(m)')
% ylabel('Distance(m)')
% str="t="+num2str(t)+"s"+"(ASPIM)";
% title(str,'FontName','Times','FontWeight','bold','FontSize',12)
% %    colormap gray
% view([90,90])
% set(gca,'FontName','Times','FontWeight','bold','FontSize',12)
% set(gcf,'Position',[100,100,500,400]);