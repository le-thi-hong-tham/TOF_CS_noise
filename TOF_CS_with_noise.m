clear all;
close all;
% clc;

%% Create signal 
%Frequency components
f=1000; %kHz = 1Mhz frequency of light wave
f1=1000; %kHz
f2=1250; %kHz

T=1/f; %duty cycle of light wave
fs=50; % rate of camera sensor
Ts=1/fs; %duty cycle of camera sensor
k= 2; %sparse level per cycles
Nc1 = 100 % number sample per cycles
Nc2 = 80;
M = 200; % number tranfers - measuments
N = Nc1 * f1/fs; % length of signal

% Light speed
c=3e8; %m

% Maximum distance range
d_max= c/(2*abs(f1-f2)*10^3);

% The true object location
d = 540; %m 

%generate signal reference
refsig1 = zeros(N,1);
ref1= zeros(Nc1,1); % signal tranfers per cycle
ref1(1,1)= 1;
l=1;
refsig2 = zeros(N,1);
ref2= zeros(Nc2,1); % signal tranfers per cycle
ref2(l,1)= 2;
% Time of the light wave flighting from the object to the imaging sensor
shiftime = 360;

for i= 1:f1/fs
    refsig1((i-1)*Nc1+1:i*Nc1) = ref1(:,1);
end
for i= 1:f2/fs
    refsig2((i-1)*Nc2+1:i*Nc2) = ref2(:,1);
end

refsig = refsig1+refsig2;
objsig = circshift(refsig,shiftime);

% for i= 1:N
%     t(i) = (i-1)/(f2*Nc2);
%    
% end

figure(1);
plot(refsig)
hold on
plot(objsig);
ylim([-0.2 4]);
xlabel('sample');
ylabel('Amplitude');
title('Generated signal');
legend('ref','obj')
% 
%% Encoded signal
% generate encode signal
Phi = randi([0 1],N,N);
y=Phi*refsig;
y1=Phi*objsig;
figure(2);
plot(y)
hold on
plot(y1);
xlabel('Samples');
ylabel('Intensity');
title('Encoded signal');
legend('ref','obj')

%% Capture signal
d= 401; % capture start location
for i=1:M
   position(i,1) = d + 1;
   while(position(i) > N)
       position(i) = position(i)-N;
   end
   d = d +2;
end

%% Adding some measurement noise.
SNR=50; % Signal-to-noise ratio (dB) 
n=awgn(y1,SNR,'measured');
outputref = zeros(M,1);
outputobj = zeros(M,1);

%% Making random measurements
A=zeros(M,N);
for i= 1:M
    outputref(i) = y(position(i));
    A(i,:) = Phi(position(i),:);
    outputobj(i) = n(position(i));
end

figure(3)
plot(outputref);
hold on
plot(outputobj);
xlabel('sample');
ylabel('Intensity');
title(sprintf('Measurement vector with noise at SNR=%d dB', SNR));
legend('ref','obj')
% 
% %Adding some measurement noise.
% SNR=40;
% e=createNoise(outputobj,SNR);
% 
% 
% %Measurement vector with noise.
% outputobj = outputobj + e ;
% figure(4);
% plot(outputref);
% hold on
% plot(outputobj);
% title(sprintf('Measurement vector with noise at SNR=%d dB', SNR));
% % 
% 
%% Recoved signal
cvx_begin
    variable xp_ref(N);
    minimize (norm(xp_ref,1));
    subject to
    A*xp_ref==outputref;
cvx_end
% 

cvx_begin
    variable xp_obj(N);
    minimize (norm(xp_obj,1));
    subject to
    A*xp_obj==outputobj; 

%     norm(A*obj-outputobj,2) <= eps
%     minimize (norm(A*obj-outputobj,2)+0.01*norm(obj,1));
cvx_end
% 
%% Compute error recovered
diff_ref = refsig - xp_ref;
recovery_error_ref = norm(diff_ref) / norm(refsig);
fprintf('recovery error: %.3f%%\n', recovery_error_ref*100);
% 
diff_obj = objsig - xp_obj;
recovery_error_obj = norm(diff_obj) / norm(objsig);
fprintf('recovery error: %.3f%%\n', recovery_error_obj*100);
% 
% 
figure(5)
plot(xp_ref)
hold on 
plot(xp_obj)
xlabel('sample');
ylabel('Amplitude');
title('Recovered signal');
legend('ref','obj');


%% Calculate phase difference (in time domain)
MaxObjLoc = 0;
MaxRefLoc = 0;
Cycle = 400;
RcvObjCyc = zeros(1,Cycle);
RcvRefCyc = zeros(1,Cycle);

% 
for i = 1:Cycle
    RcvObjCyc(i) = xp_obj(i);
    RcvRefCyc(i) = xp_ref(i);
end
% 
for i = 1:length(RcvObjCyc)
    if(RcvObjCyc(i) == max(RcvObjCyc))
        MaxObjLoc = i;
        break;
    end
end

for i = 1:length(RcvRefCyc)
    if(RcvRefCyc(i) == max(RcvRefCyc))
        MaxRefLoc = i;
        break;
    end
end
% 
LocDif = abs(MaxRefLoc - MaxObjLoc);
PDS = ((2*pi*LocDif)/Cycle);

%% Distance calculation
Distance = (c/(2*abs(f1-f2)*10^3))*(PDS/(2*pi));
fprintf('Measured Distance = %.2fm\n',Distance)

