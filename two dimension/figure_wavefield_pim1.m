clear;clc
%% 生成单元
seiscolor=load('colorbar.txt');
lex=6000;   % Length in x direction
ley=6000;   % Length in y direction
nex=200;   % The number of units in x direction
ney=200;   % The number of units in x direction
x=linspace(0,lex,nex);
y=linspace(0,ley,ney);
[x,y]=meshgrid(x,y);
dex=lex/nex;
dey=ley/ney;
X=x(:);
Y=y(:);
Point=[X,Y];
DT = delaunayTriangulation(Point);
Top=DT.ConnectivityList;
p=DT.Points;
ppx=p(:,1);
V=3000.*ones(size(ppx));
index=find(ppx<800);
V(index)=1500;
F = freeBoundary(DT);
b=unique(F(:)); % 边界节点
%% 参数设置
nt = 1500;
v=2000;
ept=dex/v;
dt=0.5*ept;
figure(1)
t=100*dt;
load U100.mat; 
trisurf(Top, p(:,1), p(:,2), U)
colorbar
colormap(gca,[seiscolor])
clim([-1e-8,1e-8])
shading interp
xlabel('Depth(m)')
ylabel('Distance(m)')
str="t="+num2str(t)+"s"+"(ASPIM)";
title(str,'FontName','Times','FontWeight','bold','FontSize',12)
%    colormap gray
view([90,90])
set(gca,'FontName','Times','FontWeight','bold','FontSize',12)
set(gcf,'Position',[100,100,500,400]);
figure(2)
t=500*dt;
load U500.mat; 
trisurf(Top, p(:,1), p(:,2), U)
colorbar
colormap(gca,[seiscolor])
clim([-1e-8,1e-8])
shading interp
xlabel('Depth(m)')
ylabel('Distance(m)')
str="t="+num2str(t)+"s"+"(ASPIM)";
title(str,'FontName','Times','FontWeight','bold','FontSize',12)
%    colormap gray
view([90,90])
set(gca,'FontName','Times','FontWeight','bold','FontSize',12)
set(gcf,'Position',[100,100,500,400]);