%Europa_detection_comparison
% Created by AGM on 27 April 2021
% Uses noise models from Panning. PlanetProfile creates structure models,
% Axisem/Instaseis generates waveforms. 
clear
%close all

addpath('/Users/marusiak/Documents/MATLAB/saclab'); % add required sac code

% set desired model parameters
noise_strength='pref'; % set desired strength A=lease noise, D=most noise, pref=preferred
ice_thickness=20; % either 5 or 20 (km)
event_depth='deep'; % either surface or deep 

cat=1; % catalog, either 0, 2 3 or 3. Most only have 0 or 1

% set necessary scales for different magnitudes makes sure denominator
% matches python code that generate sac files
Mw=[2.5,3,3.5,4,4.5 5, 5.5 6];
scale=10.^(3/2.*Mw+9.1)/(3.98e13); % scales seismograms to desired magnitude


% get records without noise
for k=1:179 % distances
 
    record_name=strcat('/Users/marusiak/Documents/Gattaca/Europa_constant_Tc/Europa_',...
        num2str(ice_thickness),'km_',event_depth,'/',num2str(k,'%03d'),'_MX');
    
    if k==10
    record_E_10=rsac([record_name 'E.SAC']); % saves waveform at distance of 10 degrees
    record_N_10=rsac([record_name 'N.SAC']);
    record_Z_10=rsac([record_name 'Z.SAC']);
    end
    recordE=rsac([record_name 'E.SAC']); % creates array of seismogram data in displacement
    recordN=rsac([record_name 'N.SAC']);
    recordZ=rsac([record_name 'Z.SAC']);
    record_E(k,:)=recordE(:,2);
    record_N(k,:)=recordN(:,2);
    record_Z(k,:)=recordZ(:,2);
end

%read in noise models
noise_name=strcat('noise_records/ice',num2str(ice_thickness),'.',noise_strength,'_cat',num2str(cat),'.MX');
noise_E=rsac([noise_name 'E']);
noise_N=rsac([noise_name 'N']);
noise_Z=rsac([noise_name 'Z']);

% plot seismograms w/out noise for sanity check
figure(1)
ylimit=max([diff(record_Z_10(:,2),2)' diff(record_N_10(:,2),2)' diff(record_E_10(:,2),2)']);
subplot(3,1,1)
plot(record_E_10(1:end-2,1),diff(record_E_10(:,2),2))
ylim([-ylimit ylimit])
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title('Radial')
subplot(3,1,2)
plot(record_N_10(1:end-2,1),diff(record_N_10(:,2),2))
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title('Transverse')
 ylim([-ylimit ylimit])
subplot(3,1,3)
plot(record_Z_10(1:end-2,1),diff(record_Z_10(:,2),2))
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title('Vertical')
 ylim([-ylimit ylimit])
 
 %method 1: grab random noise segment from original noise time series
l1=length(record_E_10(:,2));
l2=length(noise_E(:,2));
 % resample to match sampling rates
 r = randi([1 l2-l1],1,1); % grab random spot in noise model;
 delr=lh(record_E_10,'DELTA');
 deln=lh(noise_E,'DELTA');
 noise_secE=noise_E(r:r+l1*ceil(delr/deln),:); % assuming sampling rate for noise is about double record; 
 noise_secE(:,1)=noise_secE(:,1)-noise_secE(1,1);
 noise_resampE=interp1(noise_secE(:,1),noise_secE(:,2),record_E_10(:,1));
 noise_secN=noise_N(r:r+l1*2,:); % assuming sampling rate for noise is about double record; 
 noise_secN(:,1)=noise_secN(:,1)-noise_secN(1,1);
 noise_resampN=interp1(noise_secN(:,1),noise_secN(:,2),record_N_10(:,1));
 noise_secZ=noise_Z(r:r+l1*2,:); % assuming sampling rate for noise is about double record; 
 noise_secZ(:,1)=noise_secZ(:,1)-noise_secZ(1,1);
 noise_resampZ=interp1(noise_secZ(:,1),noise_secZ(:,2),record_Z_10(:,1));
 %

figure(2)
% plot record with noise
subplot(3,4,1)
plot(record_E_10(1:end-2,1),diff(noise_resampE,2))
xlabel('Time (s)')
ylabel('noise amplitude')
xlim([00 1000])
subplot(3,4,2)
plot(record_E_10(1:end-2,1),diff(record_E_10(:,2),2)*scale(1)+diff(noise_resampE,2))
xlim([00 1000])
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title(['Radial Mag', num2str(Mw(1)),' with noise'])
xlim([00 1000])
subplot(3,4,3)
plot(record_E_10(1:end-2,1),diff(record_E_10(:,2),2)*scale(2)+diff(noise_resampE,2))
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title(['Radial Mag', num2str(Mw(2)),' with noise'])
xlim([00 1000])
subplot(3,4,4)
plot(record_E_10(1:end-2,1),diff(record_E_10(:,2),2)*scale(3)+diff(noise_resampE,2))
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title(['Radial Mag', num2str(Mw(3)),' with noise'])
xlim([00 1000])
subplot(3,4,5)
plot(record_N_10(1:end-2,1),diff(noise_resampN,2))
xlabel('Time (s)')
ylabel('noise amplitude')
xlim([00 1000])
subplot(3,4,6)
plot(record_N_10(1:end-2,1),diff(record_N_10(:,2),2)*scale(1)+diff(noise_resampN,2))
title('Background Noise')
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title(['Transverse Mag', num2str(Mw(1)),' with noise'])
xlim([00 1000])
subplot(3,4,7)
plot(record_N_10(1:end-2,1),diff(record_N_10(:,2),2)*scale(2)+diff(noise_resampN,2))
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title(['Transverse Mag', num2str(Mw(2)),' with noise'])
xlim([00 1000])
subplot(3,4,8)
plot(record_N_10(1:end-2,1),diff(record_N_10(:,2),2)*scale(3)+diff(noise_resampN,2))
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title(['TransverseMag', num2str(Mw(3)),' with noise'])
xlim([00 1000])
subplot(3,4,9)
plot(record_Z_10(1:end-2,1),diff(noise_resampZ,2))
xlabel('Time (s)')
ylabel('noise amplitude')
xlim([00 1000])
subplot(3,4,10)
plot(record_Z_10(1:end-2,1),diff(record_Z_10(:,2),2)*scale(1)+diff(noise_resampZ,2))
title('Background Noise')
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title(['Vertical Mag', num2str(Mw(1)),' with noise'])
xlim([00 1000])
subplot(3,4,11)
plot(record_Z_10(1:end-2,1),diff(record_Z_10(:,2),2)*scale(2)+diff(noise_resampZ,2))
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title(['Vertical Mag', num2str(Mw(2)),' with noise'])
xlim([00 1000])
subplot(3,4,12)
plot(record_Z_10(1:end-2,1),diff(record_Z_10(:,2),2)*scale(3)+diff(noise_resampZ,2))
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title(['Vertical Mag', num2str(Mw(3)),' with noise'])
xlim([00 1000])

%noise method 2 using FFT and spectral domain to calculate median noise
%levels

dt=lh(record_E_10,'DELTA');
Fs=1/diff(record_E_10(1:2,1)); % Sample Rate 20 Hz
Fnyq=0.5*Fs; % Nyquist Frequency

% apply very broad filter 
band1=1/1; % Low Pass 1 s
band2=1/100;% High Pass 100s
[B1,A1]=butter(10,band1/Fnyq,'low');
[B2,A2]=butter(5,band2/Fnyq,'high'); 
detrended=detrend(diff(record_E_10(:,2),2));
tapered=detrended.*tukeywin(length(record_E_10(:,2))-2,0.95);
    syn_all=filtfilt(B1,A1,tapered);
    syn_E=filtfilt(B2,A2,syn_all);
detrended=detrend(diff(record_N_10(:,2),2));
tapered=detrended.*tukeywin(length(record_N_10(:,2))-2,0.95);
    syn_all=filtfilt(B1,A1,tapered);
    syn_N=filtfilt(B2,A2,syn_all);
detrended=detrend(diff(record_Z_10(:,2),2));
tapered=detrended.*tukeywin(length(record_Z_10(:,2))-2,0.95);
    syn_all=filtfilt(B1,A1,tapered);
    syn_Z=filtfilt(B2,A2,syn_all);
    
 npts=length(syn_all);
 
 % Load noise files created using python code and grab noise time series
 % (nt)
noise_name=strcat('noise_records/ice',num2str(ice_thickness),'_',noise_strength,'_cat',num2str(cat),'_MX');
nt=transpose(europa_noise(noise_name,npts,dt,'m','acc',0));  



figure(3)
% plot record with noise
subplot(3,4,1)
plot(record_E_10(1:end-2,1),nt)
xlabel('Time (s)')
ylabel('Acc amplitude')
title('Background noise')

subplot(3,4,2)
plot(record_E_10(1:end-2,1),diff(record_E_10(:,2),2)*scale(1)+nt)

xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title(['Radial Mw ', num2str(Mw(1)),' with noise'])
xlim([00 1000])
subplot(3,4,3)
plot(record_E_10(1:end-2,1),diff(record_E_10(:,2),2)*scale(2)+nt)
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title(['Radial Mw ', num2str(Mw(2)),' with noise'])
xlim([00 1000])
subplot(3,4,4)
plot(record_E_10(1:end-2,1),diff(record_E_10(:,2),2)*scale(3)+nt)
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title(['Radial Mw ', num2str(Mw(3)),' with noise'])
xlim([00 1000])
subplot(3,4,5)
plot(record_N_10(1:end-2,1),nt)
xlabel('Time (s)')
ylabel('Acc amplitude')
title('Noise')
xlim([00 1000])
subplot(3,4,6)
plot(record_N_10(1:end-2,1),diff(record_N_10(:,2),2)*scale(1)+nt)
title('Background Noise')
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
xlim([00 1000])
title(['Transverse Mw ', num2str(Mw(1)),' with noise'])
subplot(3,4,7)
plot(record_N_10(1:end-2,1),diff(record_N_10(:,2),2)*scale(2)+nt)
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
title(['Transverse Mw ', num2str(Mw(2)),' with noise'])
xlim([00 1000])
subplot(3,4,8)
plot(record_N_10(1:end-2,1),diff(record_N_10(:,2),2)*scale(3)+nt)
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
xlim([00 1000])
title(['Transverse Mw ', num2str(Mw(3)),' with noise'])
subplot(3,4,9)
plot(record_Z_10(1:end-2,1),nt)
xlim([00 1000])
xlabel('Time (s)')
ylabel('Acc amplitude')
title('Background Noise')
subplot(3,4,10)
plot(record_Z_10(1:end-2,1),diff(record_Z_10(:,2),2)*scale(1)+nt)
title('Background Noise')
xlabel('Time (s)')
xlim([00 1000])
ylabel('Acc. Amp. (m/s^2)')
title(['Vertical Mw ', num2str(Mw(1)),' with noise'])
subplot(3,4,11)
plot(record_Z_10(1:end-2,1),diff(record_Z_10(:,2),2)*scale(2)+nt)
xlabel('Time (s)')
xlim([00 1000])
ylabel('Acc. Amp. (m/s^2)')
title(['Vertical Mw ', num2str(Mw(2)),' with noise'])
subplot(3,4,12)
plot(record_Z_10(1:end-2,1),diff(record_Z_10(:,2),2)*scale(3)+nt)
xlabel('Time (s)')
ylabel('Acc. Amp. (m/s^2)')
xlim([00 1000])
title(['Vertical Mw ', num2str(Mw(3)),' with noise'])

%%
figure(4)
jj=3; % set desited magnitude from scaled array 
% plot of spectrograms with noise 
%freqrange=logspace(-1,3,1000);
subplot(3,1,1)
spectrogram((diff(record_E_10(:,2),2).*scale(jj)+nt),50,40,logspace(-1.5,1.5,100),1/dt,'yaxis');
set(gca, 'YScale', 'log')

%caxis([-156 -150)
title(['Radial Mw ', num2str(Mw(jj)),' with noise'])
subplot(3,1,2)
spectrogram((diff(record_N_10(:,2),2).*scale(jj)+nt),50,40,logspace(-1.5,1.5,100),1/dt,'yaxis');
set(gca, 'YScale', 'log')
%caxis([-157 -151])
title(['Transverse Mw ', num2str(Mw(jj)),' with noise'])
subplot(3,1,3)
spectrogram((diff(record_Z_10(:,2),2).*scale(jj)+nt),50,40,logspace(-1.5,1.5,100),1/dt,'yaxis');
set(gca, 'YScale', 'log')
title(['Vertical Mw ', num2str(Mw(jj)),' with noise'])
%caxis([-156 -151])
addpath('/Users/marusiak/Documents/MATLAB/seizmo/')
model=strcat('/Users/marusiak/Documents/GitHub/PlanetProfile/Europa',...
    '/EuropaProfile_Seawater_35WtPct_Zb',num2str(ice_thickness),'km.tvel');
%%
jj=2; % scale for desired magnitude 
for h=1:179
    %times=tauptime('mod',model,'deg',h);
    %for k=1:length(times);
        seismoE(h,:)=diff(record_E(h,:),2)'.*scale(jj)+nt;
        seismoN(h,:)=diff(record_N(h,:),2)'.*scale(jj)+nt;
        seismoZ(h,:)=diff(record_Z(h,:),2)'.*scale(jj)+nt;
        
       % SNRE(h,:)=log10(abs(seismoE(h,:)./(abs(nt')))); %signal to noise ratio;
        %SNRN(h,:)=log10(abs(seismoN(h,:)./(abs(nt'))));
        %SNRZ(h,:)=log10(abs(seismoZ(h,:)./(abs(nt'))));
         
        SNRE(h,:)=(abs(seismoE(h,:)./(abs(nt')))); %signal to noise ratio;
        SNRN(h,:)=(abs(seismoN(h,:)./(abs(nt'))));
        SNRZ(h,:)=(abs(seismoZ(h,:)./(abs(nt'))));
end
figure(5) % plot distance vs time vs signal to noise ratios
colormap(flipud(gray))
imagesc(1:1:179,1:dt:1200,SNRE')
h=colorbar;
caxis([0 max(max(SNRE))])
%caxis([1 100])
ylabel(h, 'SNR')
set(gca,'Ydir','normal')
xlabel('Distance (deg)')
ylabel('Time (s)')

title([ event_depth ' Radial Mag ' num2str(Mw(jj))])
figure(6)
colormap(flipud(gray))
imagesc(1:1:179,1:dt:1200,SNRN')
h = colorbar;
ylabel(h, 'SNR')
caxis([0 max(max(SNRN))])
%caxis([1 7])
set(gca,'Ydir','normal')
xlabel('Distance (deg')
ylabel('Time (s)')

title([ event_depth ' Transverse Mag ' num2str(Mw(jj))])
figure(7)
colormap(flipud(gray))
imagesc(1:1:179,1:dt:1200,SNRZ')
h = colorbar;
ylabel(h, 'SNR')
caxis([0 max(max(SNRZ))])
%caxis([1 10])
set(gca,'Ydir','normal')
xlabel('Distance (deg')
ylabel('Time (s)')
title([event_depth ' Vertical Mag ' num2str(Mw(jj))])

%%
figure(11)
% plot periodogram vs instrument capabilities

% load instrument responses from Panning
STS2=load('/Users/marusiak/Documents/GitHub/EuropaNoise/noise_STS2.txt');
geophone=load('/Users/marusiak/Documents/GitHub/EuropaNoise/noise_10HZ_geophone.txt');
TC=load('/Users/marusiak/Documents/GitHub/EuropaNoise/noise_Trillium_Compact.txt');
SP=load('/Users/marusiak/Documents/GitHub/EuropaNoise/noise_SP_Imperial.txt');

plot(1./STS2(:,1),20.*log10(STS2(:,2)),'b')
set(gca, 'XScale', 'log')
xlabel('Period (s)')
ylabel('Power (m^2/s^4/Hz)(dB)')
hold on
plot(1./SP(:,1),20.*log10(SP(:,2)),'r')
plot(1./geophone(:,1),20.*log10(geophone(:,2)),'m')
plot(1./TC(:,1),20.*log10(TC(:,2)),'g')

[pxx f]=periodogram((diff(record_Z_10(:,2),2).*scale(jj)+nt),[],[],Fs,'power');
pxx(1)=[];
f(1)=[];
plot(1./f,10*log10(pxx),'c--')
%plot(1./f,10.*log10(smooth(pxx,10)),'b')

legend('STS2','SP','10 Hz Geophone','TC',['Magnitude ' num2str(Mw(jj))])
    xlim([min(1./f) max(1./f)])
    
    figure(12)
    STS2_resamp=interp1(1./STS2(:,1),20*log10(STS2(:,2)),1./f(3:end),'method','extrap');
    TC_resamp=interp1(1./TC(:,1),20*log10(TC(:,2)),1./f(3:end),'method','extrap');
    for ii=1:179
       [pxx f]=periodogram(seismoZ(ii,:),[],[],Fs,'power');
       pxx(1:3)=[];
        f(1:3)=[];
        over_inst(ii,:)=smooth(10*log10(pxx),2)-TC_resamp;
    end
%imagesc(1:1:179,log10(1./f)',over_inst)
contourf(1:1:179,log10(1./f)',over_inst',[0:1:25],'LineStyle','none')
h = colorbar;
ylabel(h, 'dB over TC')
caxis([0 20])
set(gca,'Ydir','normal')
xlabel('Distance (deg)')
ylabel('Period (s)')
title([num2str(ice_thickness) ' km shell, ' num2str(event_depth) ', Vertical Comp., Mw ' num2str(Mw(jj))])

%% Compare deep and shallow events
ice_thickness=20;
for ff=1:179
    
     record_name=strcat('/Users/marusiak/Documents/Gattaca/Europa_constant_Tc/Europa_',num2str(ice_thickness),...
        'km_surface/',num2str(ff,'%03d'),'_MX');
    
    recordEsurf=rsac([record_name 'E.SAC']);
    recordNsurf=rsac([record_name 'N.SAC']);
    recordZsurf=rsac([record_name 'Z.SAC']);
    
     record_name=strcat('/Users/marusiak/Documents/Gattaca/Europa_constant_Tc/Europa_',num2str(ice_thickness),...
        'km_deep/',num2str(ff,'%03d'),'_MX');
    
    recordEdeep=rsac([record_name 'E.SAC']);
    recordNdeep=rsac([record_name 'N.SAC']);
    recordZdeep=rsac([record_name 'Z.SAC']);
    
    
    record_E_ratio(ff)=rms(diff(recordEdeep(:,2),2))./rms(diff(recordEsurf(:,2),2));
    record_N_ratio(ff)=rms(diff(recordNdeep(:,2),2))./rms(diff(recordNsurf(:,2),2));
    record_Z_ratio(ff)=rms(diff(recordZdeep(:,2),2))./rms(diff(recordZsurf(:,2),2));
end
figure(20)
hold on
plot(1:179,smooth(record_E_ratio),'b--',1:179,smooth(record_N_ratio),'r--',...
    1:179,smooth(record_Z_ratio),'k--','LineWidth',2)
ylabel('RMS Ratio Deep/Surface')
xlabel('Distance (deg)')
title([num2str(ice_thickness) 'km thick'])

legend('Radial','Transverse','Vertical')
