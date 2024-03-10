clear,clc;
seismic=read_segy_file('MODEL_P-WAVE_VELOCITY_1.25m.segy');
data=seismic.traces;
data=data(1:2500,4000:8600);
imagesc(data)
% [m,n ]=size(data);
% x=linspace(1,m,200);
% x=fix(x);
% y=linspace(1,n,200);
% y=fix(y);
% newdata=data(x,y);
imagesc(data)
%
set(gca,'XAxisLocation','top')
ylabel("Distace(m)")

xlabel("Depth(m)")

c = colorbar;
c.Label.String = 'Velocity (m/s)';
set(gca,'FontName','Times','FontWeight','bold','FontSize',12)
set(gcf,'Position',[100,100,460,250]);