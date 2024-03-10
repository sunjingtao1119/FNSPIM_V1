clear;clc
%% color spectrum of seismic 
seiscolor=load('colorbar.txt');

%% Generate triangular mesh
lex=4600;   % Length in x direction
ley=2500;   % Length in y direction
nex=460;   % The number of units in x direction
ney=250;   % The number of units in x direction
x=linspace(0,lex,nex);
y=linspace(0,ley,ney);
[x,y]=meshgrid(x,y);
dex=lex/nex;
dey=ley/ney;
X=x(:);
Y=y(:);
P=[X,Y];
DT = delaunayTriangulation(P); 
Top=DT.ConnectivityList; % Top matrix
p=DT.Points; % nodal coordinate
ppx=p(:,1);
%% velocity model
seismic=read_segy_file('MODEL_P-WAVE_VELOCITY_1.25m.segy');
data=seismic.traces;
data=data(1:2500,4000:8600);
imagesc(data)
[m,n]=size(data);
x=linspace(1,n,nex);
x=fix(x);
y=linspace(1,m,ney);
y=fix(y);
V=data(y,x);
V=V(:);
F = freeBoundary(DT);
b=unique(F(:));          % boundary node
nt = 3000;
v=max(V);
de=min([dex,dey]);
ept=dex/v;
dt=0.4*ept;
T = (1:nt)*dt;
number_of_notes = length(ppx);
fmain =40;
%% 参数设置

figure(1)
t=500*dt;
load MSPU500.mat; 
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
% set(gcf,'Position',[100,100,500,400]);
% figure(2)
% t=500*dt;
% load PU500.mat; 
% trisurf(Top, p(:,1), p(:,2), u4)
% colorbar
% colormap(gca,[seiscolor])
% clim([-1e-8,1e-8])
% shading interp
% xlabel('Depth(m)')
% ylabel('Distance(m)')
% str="t="+num2str(t)+"s"+"(SPRK4)";
% title(str,'FontName','Times','FontWeight','bold','FontSize',12)
% %    colormap gray
% view([90,90])
% set(gca,'FontName','Times','FontWeight','bold','FontSize',12)
% set(gcf,'Position',[100,100,500,400]);