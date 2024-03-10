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
% triplot(DT);
% hold on
% plot(p(nex/2,:),'o')
%% Parameter Setting
nt = 3000;
v=max(V);
de=min([dex,dey]);
ept=dex/v;
dt=0.4*ept;
T = (1:nt)*dt;
number_of_notes = length(ppx);
fmain =40;
%% Setting the source
s_t = (1-2*pi^2*fmain*(T-0.2).^2).*exp(-fmain*pi^2*(T-0.2).^2);%雷克子波
fricker = @(x) (1-2*pi^2*fmain*(x-0.2).^2).*exp(-fmain*pi^2*(x-0.2).^2);
[rec,m] = find(abs(p(:,2))==0);
b=setdiff(b,rec);
%%
dictK= TStiffness(p,Top,V); %  stiffness matrix
dictM= TMass(p,Top); %  mass matrix
K=dic2sp2(dictK,number_of_notes,number_of_notes);      % stiffness matrix
% boundary condition treatment
K(b,:)=0;
K(:,b)=0;
K(b,b)=1;
%
clear dictK;
Minv=dic2sp2(dictM,number_of_notes,number_of_notes,1); %  mass matrix
clear dictM;
Fr = freeBoundary(DT);
boundarynodes=unique(Fr(:));    
U = zeros(number_of_notes,nt);
% [source_x,m] = find(abs(p(:,1))>=lex/2-de&abs(p(:,1))<=lex/2+de&abs(p(:,2))>=ley/2-de&abs(p(:,2))<=ley/2+de);
[source_x,m] = find(abs(p(:,2))==0&abs(p(:,1))>=lex/2-dey&abs(p(:,1))<=lex/2+dey);
[rec,m] = find(abs(p(:,1))>=1000-dex&abs(p(:,1))<=1000+dex&abs(p(:,2))>=50-dey&abs(p(:,2))<=50+dey);
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
% isrc=nex/2;
u = zeros(number_of_notes,1);
uold=u;
vold=u;
f=zeros(number_of_notes,1);
disp=zeros(nt,1);
tic
for i = 1:500
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
% u4(boundarynodes)=0;
% disp(i,:)=u4(rec(1));
% trisurf(Top, p(:,1), p(:,2), u4)
% xlim([0,lex])
% ylim([0,ley])
% colorbar
% clim([-1e-8,1e-8]);
% colormap(gca,[seiscolor]);
%  shading interp
%  view([0,-90]);
%  dx=lex*5/1000;
%  dy=ley*5/1000;
%  set(gcf,'unit','centimeters','position',[10 5 dx dy]);
%  drawnow() 
end
save MSPU500.mat u4

