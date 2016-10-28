function [Sxx,Gxx,N,df,f]=PowerSpectralDensity(times,volts)
%==========================================================================
%This function calculates power spectral density
%inputs : times : time history
%         volts : voltage for a total time histroy
%outputs : Sxx   : double-sided power spectral density
%         Gxx   : single-sided power spectral density         
%         N     : number of points
%         df    : df
%         f     : f
%==========================================================================
N=size(times,1);
dt=times(2)-times(1);
T=N*dt;

df=1/T;
fmax=1/dt;
f=linspace(0,fmax-df,N)';

LS=fft(volts)*dt;
Sxx=1/T.*LS.*conj(LS);

if mod(N,2)==0
 Gxx(2:N/2,1)=2.*Sxx(2:N/2);
 Gxx(1,1)=Sxx(1);
 Gxx(N/2+1,1)=Sxx(N/2+1);
else
 Gxx(1,1)=Sxx(1);
 Gxx(2:floor(N/2)+1,1)=2.*Sxx(2:floor(N/2)+1);
end
