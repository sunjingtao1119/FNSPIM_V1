clear;clc
%% color spectrum of seismic 
seiscolor=load('colorbar.txt');
%% Generate triangular mesh
lex=6000;   % Length in x direction
ley=6000;   % Length in y direction
nex=400;   % The number of units in x direction
ney=400;   % The number of units in x direction
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
V=3000.*ones(size(ppx)); %velocity
% index=find(ppx<3150);
% V(index)=1500;
F = freeBoundary(DT);
b=unique(F(:));          % boundary node
%% Parameter Setting
nt = 1500;
v=3000;
ept=dex/v;
dt=0.5*ept;
T = (1:nt)*dt;
number_of_notes = length(ppx);
fmain =40;
%% Setting the source
s_t = (1-2*pi^2*fmain*(T-0.2).^2).*exp(-fmain*pi^2*(T-0.2).^2);%雷克子波
fricker = @(x) (1-2*pi^2*fmain*(x-0.2).^2).*exp(-fmain*pi^2*(x-0.2).^2);
%%
dictK= TStiffness(p,Top,V); %  stiffness matrix
dictM= TMass(p,Top); %  mass matrix
K=dic2sp2(dictK,number_of_notes,number_of_notes); % stiffness matrix
clear dictK;
Minv=dic2sp2(dictM,number_of_notes,number_of_notes,1); %  mass matrix
clear dictM;
Fr = freeBoundary(DT);
boundarynodes=unique(Fr(:));    
U = zeros(number_of_notes,nt);
rec=100; y=linspace(1000,5000,100);
s_x=6000;s_y=6000;
[source_x,m] = find(abs(p(:,1))>=s_x/2-dex&abs(p(:,1))<=s_x/2+dex&abs(p(:,2))>=s_y/2-dey&abs(p(:,2))<=s_y/2+dey);
% [source_x,m] = find(abs(p(:,1))>=le/2-de&abs(p(:,1))<=le/2+de&abs(p(:,2))>=le/2-de&abs(p(:,2))<=le/2+de);
% rx=2500;
% for i=1:100
%     [rec,m] = find(abs(p(:,1))>=rx-dex&abs(p(:,1))<=rx+dex&abs(p(:,2))>=y(i)-dey&abs(p(:,2))<=y(i)+dey);
%     indrec(i,:)=rec; % 接受点位置
% end

rx=1500;yx=3000;
[rec1,m] = find(abs(p(:,1))>=rx-dex&abs(p(:,1))<=rx+dex&abs(p(:,2))>=yx-dey&abs(p(:,2))<=yx+dey);
%% The fourth-order R-K coefficient
e1=0.13883725894365473;
e2=0.46958619250378464;
e3=0.751399209882663;
e4=-0.3598226613301023;
d1=0.3726518368174738;
d2=0.41264784985125225;
d3=-0.04864313400799411;
d4=0.26334344733926796;

L=Minv*K;
isrc=source_x(1);
u = zeros(number_of_notes,1);
uold=u;
vold=u;
f=zeros(number_of_notes,1);
tic
% isrc=ne/2;

nt=500;
% disp=zeros(nt,100);
disp=zeros(nt,1);
for i = 1:nt
f(isrc)= fricker((i-1)*dt+d1*dt);   
u1=uold+e1*dt.*vold;
v1=vold+d1*dt.*(-L*u1+Minv*f);
f(isrc)= fricker((i-1)*dt+(d1+d2)*dt);    
u2=u1+e2*dt.*v1;
v2=v1+d2*dt.*(-L*u2++Minv*f);
f(isrc)= fricker((i-1)*dt+(d1+d2+d3)*dt);
u3=u2+e3*dt.*v2;
v3=v2+d3*dt.*(-L*u3+Minv*f);
f(isrc)= fricker((i-1)*dt+(d1+d2+d3+d4)*dt);
u4=u3+e4*dt.*v3;
v4=v3+d4*dt.*(-L*u4+Minv*f);
uold=u4;
vold=v4; 
u4(boundarynodes) = 0;
% disp(i,:)=u4(indrec(:,1));
 if i==10
 save PU101.mat u4
 end
 if i==500
 save PU5001.mat u4
 end
 disp(i)=u4(rec1(1));
trisurf(Top, p(:,1), p(:,2), u4)
     tsr=num2str(i);
     title(tsr)
% plot3(p(indrec(:,1),1),p(indrec(:,1),2),u4(indrec(:,1)),"*");
colorbar
clim([-1e-8,1e-8]);

colormap(gca,[seiscolor]);
 shading interp
 view([90,90])
 drawnow() 
end
 toc
%  save Prkdisp.mat disp
%  trisurf(Top, p(:,1), p(:,2), u4)
% colorbar
% clim([-1e-8,1e-8]);
% colormap(gca,[seiscolor]);
%  shading interp
%  view([90,90]);
%  drawnow() 
