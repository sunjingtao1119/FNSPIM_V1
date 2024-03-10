clear;clc
%% color spectrum of seismic 
seiscolor=load('colorbar.txt');

%% Generate triangular mesh
lex=2000;   % Length in x direction
ley=2000;   % Length in y direction
nex=200;   % The number of units in x direction
ney=200;   % The number of units in x direction
x=linspace(0,lex,nex);
y=linspace(0,ley,ney);
[x,y]=meshgrid(x,y);
dex=lex/nex;
dey=ley/ney;
X=x(:);
Y=y(:);
P=[X,Y];
DT = delaunayTriangulation(P);
F = freeBoundary(DT);
triplot(DT)
hold on
plot(x(F),y(F),'-r','LineWidth',2)
