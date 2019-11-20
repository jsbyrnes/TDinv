%%set up variables

clear all
close all
clc

re = 6371; %radius of the earth

dt = .05;%in seconds
len = 330;%in seconds. Note that final trace will be shifted when generated, and again when deconvolved
name = 'darwin\_test';%note that if the name is too long to program won't work right
pf = 'f';
modeconv = 'y';
phase = 'S';
fc = .15; %in Hz

sn = .1;%signal to noise ratio

halfspace = 1;%convolve with clipped waveform or not

filename = 'CULLED.Iguana7_mb55.13469.P';

depth = 300;%in km

high = .1;%in Hz
low = 0.03;%in Hz

rp = .11; %in s/km

cd '/research/users/jbyrnes/Reflectivity_Synthetics/matlabscripts'

%% make vp and thicknesss USE THIS FOR CUSTOM MODEL
%a depth/thickness of zero terminates the code(extends layer to infinity)

%Word of advice for S synthetics - if you want to look at precurser P
%arrivals, put a layer of very small velocity change far deeper than what
%you want to actually look at. This won't affect the result but will make
%it much easier to read

vp = [5.0 5.8 7.1 6.4];
thickness = [10 15 40 0.0];

% drive icmod by making a c script

writeicmod_5(vp, thickness, name, 2);

!chmod a+x icmodinput.csh
!./icmodinput.csh

%% load external velocity model - pick one, comment out the others
% program will not currently work if length(z) > 145!!!

%z = 1:8:depth;
%[vp, vs] = write_iaspei(z);

%[vp, z] = write_darwinsmodel(depth, name, 31, 31);

plot(vp, z);

%% run respknt to make sac format synthetics
%Uses a delta function for a source, so this is basically the Green's
%function

writerespknt(name, phase, dt, len, rp, pf, modeconv);

!chmod a+x respkntinput.csh
!./respkntinput.csh

zname = [name '_sp.z'];
rname = [name '_sp.r'];

[headerz, Zcomp] = sac2mat(zname, 'le');
[headerr, Rcomp] = sac2mat(rname, 'le');

n = length(Zcomp);

%find the max in the vericle channel to make time vector.

if phase == 'P'
    
    [amp arrival] = max(Zcomp);

elseif phase == 'S'
    
    [amp arrival] = max(Rcomp);
    
end
 
% Zcomp = Zcomp/amp;
% Rcomp = Rcomp/amp;

t = (arrival - 1)*dt*-1:dt:(n-arrival)*.05;

% for debugging/interest reasons, you can plot out some of the unfiltered traces
% some high frequency sinc noise
figure
plot(t, Zcomp);
xlabel('Seconds');
title( [name ' - Zcomp']);
figure
plot(t, Rcomp);
xlabel('Seconds');
title( [name ' - Rcomp']);

%% Convolution
%If you want to simulate a band limited source, convolve with a band
%limited source function and produce a more realistic trace, then
%deconvolved for some more accrute testing. Alternatively you could save in
%med_bow format and run it through too
%remember that the time series can not be on top of each other at first
%point. The traces will shift by the zero offset of the source.

%load/make source
%must be same length as the trace
% source = zeros(n,1);
% 
% source(n-100:n) = 1;

load('source45966.mat');

extra = n - length(source);

source = [source; zeros(extra, 1)];

S = fft(source);
Gz = fft(Zcomp);
Gr = fft(Rcomp);

Tz = Gz.*S;
Tr = Gr.*S;

Zcomp = real(ifft(Tz)); %should be real, but take the real part to avoid round off errors
Rcomp = real(ifft(Tr));

%redefine t = 0, avoids complications related to shifting
%This is useful because when you do real reciever functions, you don't have
%the same control so the max is the easiest wasy to set the zero
%good to renormalize too, can't hurt
if phase == 'P'
    
    [amp arrival] = max(Zcomp);

elseif phase == 'S'
    
    [amp arrival] = max(Rcomp);
    
end

Zcomp = Zcomp/amp;
Rcomp = Rcomp/amp;

t = (arrival - 1)*dt*-1:dt:(n-arrival)*.05;

figure(60)
plot(t, source);
title('Source function');
figure(61)
plot(t, Zcomp)
title('Zcomp post convolution');
xlabel('Time');
figure(62)
plot(t, Rcomp);
title('Rcomp post convolution');
xlabel('Time');

%% Add Noise if so desired
%adds real seismic noise, taken from the galapagos data set
%vector is 20,000 unit long. Make noise by taking a random starting point
%up to whatever comes after that
load('noise.mat')

startingpoint = round(rand(1)*(20000-n));
Zcomp = Zcomp + sn*noise(startingpoint:startingpoint + n - 1)';
startingpoint = rand(1)*(20000-n);
Rcomp = Rcomp + sn*noise(startingpoint:startingpoint + n - 1)';

figure(101)
plot(t, Zcomp)
title('Zcomp with noise');
xlabel('Time');
figure(102)
plot(t, Rcomp);
title('Rcomp with noise');
xlabel('Time');

%% Filtering and rotation into P-SV

Zcomp = bandpassfilt(Zcomp, dt, high, low);
Rcomp = bandpassfilt(Rcomp, dt, high, low);

figure(21)
plot(t, Zcomp);
title( [name ' - Zcomp']);
xlabel('Seconds');

figure(22)
plot(t, Rcomp);
title( [name ' - Rcomp']);
xlabel('Seconds');

% do the rotation to P-SV(velocity and rp known perfectly, but will leave ~.05 on the daughter channel)
if phase == 'P'
    
    [Pcomp, SVcomp] = ZR2PSV(Zcomp, Rcomp, rp, vp(1));
    
elseif phase == 'S'
    
    [SVcomp, Pcomp] = ZR2PSV(Rcomp, Zcomp, rp, vp(1)/1.73);

end

if phase == 'P'
    
    [amp arrival] = max(Pcomp);

elseif phase == 'S'
    
    [amp arrival] = max(SVcomp);
    
end

SVcomp = SVcomp/amp;
Pcomp = Pcomp/amp;
    
figure(23)
plot(t, Pcomp);
xlabel('Seconds');
title( [name ' - Pcomp']);

figure(24)
plot(t, SVcomp);
xlabel('Seconds');
title( [name ' - SVcomp']);

%% ETMTC deconvolution

deconshift = arrival*.05;

if phase == 'P'
   
    [RF, ~] = mtmdeconSP(Pcomp, SVcomp, length(SVcomp), 200, 3, 2.5, fc, dt, deconshift);

elseif phase == 'S'
    
    [RF, ~] = mtmdeconSP(SVcomp, Pcomp, length(SVcomp), 200, 7, 4, fc, dt, deconshift);
    
end

figure(30)
plot(t, RF);
xlabel('Seconds')
title('Multitaper Reciever Function')

%% load a file for comparasion if you want
%You don't have to, the code's feelings won't be hurt
%time is set adhoc
load('RF.Iguana7_mb55.P.14.MPLSQ.0.03-0.25Hz.mat')
t = -7.5:.05:42.45;
t = t';
figure(50)
plot(t, RF_D.C2)

%% fit the synthetics into a med_bow format file
% a lot of stuff is hardwired here, so be careful

load([filename '.mat']);

[num , ~] = size(D);

%if its from Med_bow all these should be the same length(or Medbow gets
%mad)
len = length(D(1).C1);

%clip out the arrival and converted phase train
%hardwired bits, so make sure it doesn't remove stuff you need
if phase == 'P'
    
    dataV = Pcomp(arrival - 100:end);
    dataSV = SVcomp(arrival - 100:end);

elseif phase == 'S'
    
    dataV = wrev(Pcomp(1:arrival + 40));
    dataSV = wrev(SVcomp(1:arrival + 40));

end 

dVlen = length(dataV);
dSVlen = length(dataSV);

for i = 1:num
   
    D(i).C1 = [zeros(300,1); dataV; zeros(len - (300 + dVlen), 1)];
    D(i).C2 = [zeros(300,1); dataSV; zeros(len - (300 + dSVlen), 1)];
    D(i).C3 = zeros(len, 1);
end

savename = [filename '_' name '.mat'];

save(savename, 'D', 'P', 'SC', 'T', 'TW', 'flip_time', 'is_post_cull');













