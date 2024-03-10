clear;clc
%% 生成单元
seiscolor=load('colorbar.txt');
le=1000;  % size 2000*2000
ne=100;   % 
[x,y]=meshgrid(linspace(0,le,ne));
de=le/ne;
X=x(:);
Y=y(:);
P=[X,Y];
DT = delaunayTriangulation(P);
F = freeBoundary(DT);
b=unique(F(:)); % 边界节点
%生成结构化三角网格
p=DT.Points;
p=p';
t=DT.ConnectivityList;
t=t';
num_nodes = length(p);
num_elements = length(t);
t_s = sprintf('nodes=%d  elements=%d',num_nodes,num_elements);
title(t_s)
%% 参数设置
nt = 1500;
v = 1500; 
ept=de/v;
dt=0.5*ept;
T = (1:nt)*dt;
number_of_notes = length(p);
fmain =40;
%% 设置震源
s_t = (1-2*pi^2*fmain*(T-0.2).^2).*exp(-fmain*pi^2*(T-0.2).^2);%雷克子波
fricker = @(x) (1-2*pi^2*fmain*(x-0.2).^2).*exp(-fmain*pi^2*(x-0.2).^2);

%% 调用函数计算矩阵
% FEM = FEM_2D_func();
% S = FEM.assemble_matrix_2D(p, t);%刚度矩阵
Top=t';
pp=p';
tic
dictK= TStiffness(pp,Top); % 
toc
tic
dictM= TMass(pp,Top); % 
toc
tic
K=dic2sp2(dictK,number_of_notes,number_of_notes);
toc
clear dictK;
tic
Minv=dic2sp2(dictM,number_of_notes,number_of_notes,1);
toc
clear dictM;
Fr = freeBoundary(DT);
boundarynodes=unique(Fr(:)); % 边界节点
U = zeros(number_of_notes,nt);
[m,source_x] = find(abs(p(1,:))>=le/2-de&abs(p(1,:))<=le/2+de&abs(p(2,:))>=le/2-de&abs(p(2,:))<=le/2+de);

%% 利用递推关系求波场值
e1=0.13883725894365473;
e2=0.46958619250378464;
e3=0.751399209882663;
e4=-0.3598226613301023;
d1=0.3726518368174738;
d2=0.41264784985125225;
d3=-0.04864313400799411;
d4=0.26334344733926796;
% Minv=diag(1./diag(M));
% Minv=sparse(Minv);
% clear M
L=Minv*K.*v^2;
% L=sparse(L);
isrc=source_x(1);
u = zeros(number_of_notes,1);
uold=u;
vold=u;
f=zeros(number_of_notes,1);
for i = 1:500
tic
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
  toc
trisurf(t(1: 3, :)', p(1, :)', p(2, :)', u4)
colorbar
clim([-1e-8,1e-8]);
colormap(gca,[seiscolor]);
 shading interp
 view([90,90]);
 drawnow() 
 toc
end
