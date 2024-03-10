%%
disp("NSPIM time")
tic
NSPIM4
toc
%%
disp("FNSPIM_10 time")
tic
FNSPIM_10
toc
%%
disp("FNSPIM_15 time")
tic
FNSPIM_15
toc
%%
disp("FNSPIM_25 time")
tic
FNSPIM_25
toc
clear;
%

figure1 = figure(1);
figure2 = figure(2);
ax1=axes('Parent',figure1);
hold(ax1,'on');
ax2=axes('Parent',figure2);
hold(ax2,'on');
%%
load xg
load u100_10
plot(ax1,xg,u,'-.r',LineWidth=1.5)
load u200k_10
plot(ax2,xg,u,'-.r',LineWidth=1.5)

%%
load u100_15
plot(ax1,xg,u,'m',LineWidth=1.5)
load u200k_15
plot(ax2,xg,u,'m',LineWidth=1.5)
%%
load xg
load u100_25
plot(ax1,xg,u,'g',LineWidth=1.5)
load u200k_25
plot(ax2,xg,u,'g',LineWidth=1.5)
%%
load u100
plot(ax1,xg,u,'--b',LineWidth=1.5)
load u200k
plot(ax2,xg,u,'--b',LineWidth=1.5)

ylim(ax1,[-1,1])
ylabel(ax1,'Amplitude','FontSize',15)
xlabel(ax1,'x/(m)','FontSize',15)
set(ax1,'FontSize',12);


ylim(ax2,[-1,1])
ylabel(ax2,'Amplitude','FontSize',15)
xlabel(ax2,'x/(m)','FontSize',15)
set(ax2,'FontSize',12);

