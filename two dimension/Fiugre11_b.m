lay2_PRK_Figure_11
t=100*dt;
load PU100.mat; 
trisurf(Top, p(:,1), p(:,2), U)
colorbar
colormap(gca,[seiscolor])
clim([-1e-8,1e-8])
shading interp
xlabel('Depth(m)')
ylabel('Distance(m)')
str="t="+num2str(t)+"s";
title(str,'FontName','Times','FontWeight','bold','FontSize',12)
%    colormap gray
view([90,90])
set(gca,'FontName','Times','FontWeight','bold','FontSize',12)
set(gcf,'Position',[100,100,500,400]);
figure(2)
t=500*dt;
load PU500.mat; 
trisurf(Top, p(:,1), p(:,2), U)
colorbar
colormap(gca,[seiscolor])
clim([-1e-8,1e-8])
shading interp
xlabel('Depth(m)')
ylabel('Distance(m)')
str="t="+num2str(t)+"s";
title(str,'FontName','Times','FontWeight','bold','FontSize',12)
%    colormap gray
view([90,90])
set(gca,'FontName','Times','FontWeight','bold','FontSize',12)
set(gcf,'Position',[100,100,500,400]);