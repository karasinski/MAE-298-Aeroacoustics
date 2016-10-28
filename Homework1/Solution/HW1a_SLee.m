%=======================================================================
%                       EAE 298 HW#1
%         Read in a boom signal and compute SPL
%         Date : 09/30/2016  
%         Author: Prof. Seongkyu Lee
%=======================================================================
clear all
close all


%% Read in a file
filename='Boom_F1B2_6.wav';
[yy,fs]=audioread(filename);
dt=1/fs;
N=size(yy,1);
T=N*dt;
df=1/T;
times=(((1:N)-1)*dt)';
f=(((1:N)-1)*df)';

%% Applying the calibration factor
M=-116; %pascal/voltage ratio
y=yy*M;
peak=max(y)

%% Compute power spectral density and sound pressure level
[Sxx,Gxx,N,df,f]=PowerSpectralDensity(times,y);
SPL=20*log10(sqrt(abs(Gxx*df))/(20*10^-6));

%% Convert narrow band to one-thrid octave band and octave band spectrum and OASAPL
[SPL13,fcentre,SPLO,fcentreO,OASPL]=BroadbandSpectrum(SPL,f);
OASPL

%% Audioplay
fs=44100;
pvolts=yy/max(abs(yy));
pvolts=[pvolts];
ppobj=audioplayer(pvolts,fs);
play(ppobj);


%% Plot results
% Plot time history of a boom
figure(1)
hh=plot(times,y);
set(hh, 'linewidth', 2)
xlabel('Time [s]','fontsize', 16);
ylabel('Pascals [Pa]', 'fontsize', 16);
title('N wave shape of a sonic boom', 'fontsize', 18);
set(gcf,'Position',[100 100 900 600]);
set(gca,'ylim',[-100 100],'ytick',-100:20:100, 'fontsize', 18)

%Plot sound pressure level (narrow band)
figure(2)
h=plot(f(1:N/2+1),abs(Gxx));
set(h, 'linewidth',2);
xlabel('frequency [Hz]', 'fontsize', 16);
ylabel('Gxx (Pa^2/Hz)', 'fontsize', 16);
title('Gxx', 'fontsize', 18);
set(gcf,'Position',[100 100 900 600]);
set(gca,'fontsize', 18)
%set(gca, 'xlim', [0 200], 'xtick',0:200:1000, 'fontsize', 14)

%Plot sound pressure level (narrow band)
figure(3)
h=plot(f(1:N/2+1),SPL);
set(h, 'linewidth',2);
xlabel('frequency [Hz]', 'fontsize', 16);
ylabel('SPL[dB]', 'fontsize', 16);
title('Narrow band SPL', 'fontsize', 18);
set(gcf,'Position',[100 100 900 600]);
set(gca,'fontsize', 18)
%set(gca, 'xlim', [0 200], 'xtick',0:200:1000, 'fontsize', 14)

%Plot sound pressure level (one-thrid octave band)
figure(4)
h=semilogx(fcentre,SPL13);
set(h, 'linewidth',2);
xlabel('frequency [Hz]', 'fontsize', 16);
ylabel('SPL_{1/3} [dB]', 'fontsize', 16);
title('One-thrid octave band spectrum', 'fontsize', 18);
set(gcf,'Position',[100 100 900 600]);
set(gca,'fontsize', 18)

%Plot sound pressure level (one-thrid octave band)
figure(5)
h=semilogx(fcentreO,SPLO);
set(h, 'linewidth',2);
xlabel('frequency [Hz]', 'fontsize', 16);
ylabel('SPL_{oct} [dB]', 'fontsize', 16);
title('Octave band spectrum', 'fontsize', 18);
set(gcf,'Position',[100 100 900 600]);
set(gca,'fontsize', 18)












