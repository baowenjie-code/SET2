clc;clear;close all

%% 信号产生
fs = 200;
N = 400;
t = (0:N-1)/fs;
x1 = (1*exp(-0.5*(t-1).^2)).*cos(2*pi*(30*t+10*t.^2+1.5*sin(4*pi*t)));
ifreq = 30+20*t+4*pi*1.5*cos(4*pi*t);
%信号
x = x1;
load('noise.mat');
xn = x1+noise;
%窗参数
s = 0.018;

%% 时域图
figure;
plot(t,x,'b','linewidth',0.5);
ylim([-2 2]); xlim([0 2]);
xlabel('Time(s)','FontSize',8);
ylabel('Amplitude','FontSize',8);
set(gca,'FontSize',8);
set(gca,'Position',[.2 .28 .75 .65]);

figure;
plot(t,ifreq,'b','linewidth',0.5);
ylim([0 100]); xlim([0 2]);
xlabel('Time(s)','FontSize',8);
ylabel('Frequency(Hz)','FontSize',8);
set(gca,'FontSize',8);
set(gca,'Position',[.25 .28 .7 .65]);

%% SET1
tf_type ='SET1';
gamma = 0.0000;
%0)计算TFR
[~,SET1,~,~,Rep,~,q,t,f] = TFM(x,fs,s,tf_type,gamma);

figure('color',[1 1 1]);
imagesc(t,f,abs(SET1'));
xlim([0 2]);
axis xy
xlabel('Time(s)','FontSize',8);
ylabel('Amplitude','FontSize',8);
set(gca,'FontSize',8);
set(gca,'Position',[.2 .28 .75 .65]);

%1）重构
[ExtractTFR,~] = ExtractOneRidge2SubTFR(SET1,fs,s,'F',tf_type,Rep,q,0.1);
[Reconstruction,~] = ITFM(ExtractTFR,fs,s,'F',tf_type);
%2) 重构指标
xr = Reconstruction';

figure
plot(t,x,'b',t,x-xr,'k','linewidth',0.5);
ylim([-2 3]); xlim([0 2]);
xlabel('Time(s)','FontSize',8);
ylabel('Amplitude','FontSize',8);
set(gca,'FontSize',8);
set(gca,'Position',[.2 .28 .75 .65]);


%% SET2
tf_type ='SET2';
gamma = 0.0000;
%0)计算TFR
[~,SET2,~,~,Rep,~,q,~,~] = TFM(x,fs,s,tf_type,gamma);

figure('color',[1 1 1]);
imagesc(t,f,abs(SET2'));
xlim([0 2]);
axis xy
xlabel('Time(s)','FontSize',8);
ylabel('Amplitude','FontSize',8);
set(gca,'FontSize',8);
set(gca,'Position',[.2 .28 .75 .65]);
%1）重构
[ExtractTFR,~] = ExtractOneRidge2SubTFR(SET2,fs,s,'F',tf_type,Rep,q,0.1);
[Reconstruction,~] = ITFM(ExtractTFR,fs,s,'F',tf_type);
%2) 重构指标
xr = Reconstruction';
figure
plot(t,x,'b',t,x-xr,'k','linewidth',0.5);
ylim([-2 3]); xlim([0 2]);
xlabel('Time(s)','FontSize',8);
ylabel('Amplitude','FontSize',8);
set(gca,'FontSize',8);
set(gca,'Position',[.2 .28 .75 .65]);


