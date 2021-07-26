function [noise,t]=europa_noise(name,npts,dt,comp,unit,flag)
% This function generates random noise using acceleration spectrum of VBB
% seismometer. The spectrum is calculated for March 2019. 
% 
% The variable 'comp' specifies the channel of the seismometer;
% 'unit' represents displacement (dis) or acceleration (acc);
% 'flag' decides weather to apply a filter.

% AGM adjusted May 4th. use matrix of inputted noise rather than read in VBB data 
%previous method read in period and amplitude 

%MHU_acc=csvread('MHU.02.noise.csv',1);
%MHV_acc=csvread('MHV.02.noise.csv',1);
%MHW_acc=csvread('MHW.02.noise.csv',1);
%BZC_acc=csvread('BZC.58.noise.
MHU=load([name 'E.mat']);
MHU_acc=MHU.median';
MHV=load([name 'N.mat']);
MHV_acc=MHV.median';
MHW=load([name 'Z.mat']);
MHW_acc=MHW.median';
fs=1./MHW_acc(:,1);

%MHUa(:,2)=diff(MHU(:,2),2);
%MHVa(:,2)=diff(MHV(:,2),2);
%MHWa(:,2)=diff(MHW(:,2),2);
%[MHU_acc(:,2),MHU_acc(:,1)]=pwelch(MHUa(:,2),800,775,[],fs,'psd');
%[MHV_acc(:,2),MHV_acc(:,1)]=pwelch(MHVa(:,2),800,775,[],fs,'psd');
%[MHW_acc(:,2),MHW_acc(:,1)]=pwelch(MHWa(:,2),800,775,[],fs,'psd');
%MHU_acc(1,:)=[]; % removes 0 freq
%MHV_acc(1,:)=[];
%MHW_acc(1,:)=[];

%MHU_acc(:,2)=fft(MHUa(:,2))./(sqrt(length(MHUa(:,2))./2./lh(MHU,'DELTA')));
%MHV_acc(:,2)=fft(MHVa(:,2))./(sqrt(length(MHVa(:,2))./2./lh(MHV,'DELTA')));
%MHW_acc(:,2)=fft(MHWa(:,2))./(sqrt(length(MHWa(:,2))./2./lh(MHW,'DELTA')));
%f=fftfreq(length(MHUa),lh(MHU,'DELTA'));
%MHU_acc(:,1)=1./f;
%%MHV_acc(:,1)=1./f;
%MHW_acc(:,1)=1./f;
%MHW_acc(:,1)=1./MHW_acc(:,1); % convert frequency to period
%MHV_acc(:,1)=1./MHV_acc(:,1); % convert frequency to period
%MHU_acc(:,1)=1./MHU_acc(:,1); % convert frequency to period
% Convert Acceleration to Displacement
w=2*pi./MHU_acc(:,1);
MHU_dis(:,1)=MHU_acc(:,1);
MHU_dis(:,2)=MHU_acc(:,2)./w.^2;
w=2*pi./MHV_acc(:,1);
MHV_dis(:,1)=MHV_acc(:,1);
MHV_dis(:,2)=MHV_acc(:,2)./w.^2;
w=2*pi./MHW_acc(:,1);
MHW_dis(:,1)=MHW_acc(:,1);
MHW_dis(:,2)=MHW_acc(:,2)./w.^2;
%w=2*pi./BZC_acc(:,1);
%BZC_dis(:,1)=BZC_acc(:,1);
%BZC_dis(:,2)=BZC_acc(:,2)./w.^2;

% Choose Component
switch unit
    case {'acc','ACC','acceleration'}
        switch comp
            case {'MHU','U','u'}
                 ps=MHU_acc;
            case {'MHV','V','v'}
                 ps=MHV_acc;
            case {'MHW','W','w'}
                 ps=MHW_acc;
            case {'BZC','C','c'}
                 ps=BZC_acc;
            case {'mean','M','m'}
                 % Average Model
                 ps(:,1)=MHU_acc(:,1);
                 ps(:,2)=(MHU_acc(:,2)+MHV_acc(:,2)+MHW_acc(:,2))/3;
        end        
    case {'dis','DIS','displacement'}
        switch comp
            case {'MHU','U','u'}
                 ps=MHU_dis;
            case {'MHV','V','v'}
                 ps=MHV_dis;
            case {'MHW','W','w'}
                 ps=MHW_dis;
            case {'BZC','C','c'}
                 ps=BZC_dis;
            case {'mean','M','m'}
                 % Average Model
                 ps(:,1)=MHU_dis(:,1);
                 ps(:,2)=(MHU_dis(:,2)+MHV_dis(:,2)+MHW_dis(:,2))/3;
        end
end
   
% Generate Noise 
ps(:,2)=10.^(ps(:,2)./20);
freq = 1./ps(:,1);
freq_int=abs(fftfreq(npts,dt)); %1/2*1/dt= NyQuist, period/2
%freq_int=logspace(log10(1/(npts*dt)),log10(0.5/dt),npts);
np_int=interp1(freq,ps(:,2),freq_int,'linear','extrap');
N=floor((length(freq_int)-1)/2);
N2=length(freq_int);
alph=rand(1,N)*2*pi;
phases=cos(alph)+1j*sin(alph);
np_int(2:N+1)=np_int(2:N+1).*phases;
np_int(N2:-1:N2-N+1)=conj(np_int(2:N+1));
nt=real(ifft(np_int*sqrt(npts/2/dt)));
nt=nt-mean(nt);

% Filter Noise and Apply window function
if  flag
    Fs=1/dt;
    Fnyq=0.5*Fs; % Nyquist Frequency
    % Tukey Windloow Function
    w=tukeywin(npts,0.05);
    nt=transpose(nt).*w;
    t=transpose(0:dt:(npts-1)*dt);
    % Bandpass filter
    band=[1/50,1/1];
    [B,A]=butter(4,band/Fnyq,'bandpass'); 
    noise=filtfilt(B,A,nt);
else
    noise = nt;
    m1=mean(abs(MHU_acc(:,2)+MHV_acc(:,2)+MHW_acc(:,2))/3); % median value of original noise signal
    m2=mean(abs(noise));
    %noise=noise*m1/m2;
    
    t=transpose(0:dt:(npts-1)*dt);
end

end
