clear all;
close all;
% clc;

f=1000; %kHz = 1Mhz frequency of light wave
f1=1000;
f2=1250;
T=1/f; %duty cycle of light wave
fs=50; % rate of camera sensor
Ts=1/fs; %duty cycle of camera sensor
k= 1; %sparse level per cycles
Nc1 = 100 % number sample per cycles
Nc2 = 80;
M = Nc1 * (f1+f2)/500; % number tranfers - measuments
N = Nc1 * f1/fs; % length of signal

%generate signal reference
refsig1 = zeros(N,1);
ref1= zeros(Nc1,1); % signal tranfers per cycle
ref1(1,1)= 1;

refsig2 = zeros(N,1);
ref2= zeros(Nc2,1); % signal tranfers per cycle
ref2(50,1)= 2;
% time of the light wave flighting from the object to the imaging sensor
% shiftime = round(N/M);% plot(refsig);
shiftime = 18;

objsig1 = zeros(N,1);
obj1= zeros(Nc1,1); % signal tranfers per cycle
obj1(1+shiftime: 1+shiftime)= 1;

objsig2 = zeros(N,1);
obj2= zeros(Nc2,1); % signal tranfers per cycle
obj2(50+shiftime: 50+shiftime)= 2;

for i= 1:f1/fs
    refsig1((i-1)*Nc1+1:i*Nc1) = ref1(:,1);
    objsig1((i-1)*Nc1+1:i*Nc1) = obj1(:,1);
end
for i= 1:f2/fs
    refsig2((i-1)*Nc2+1:i*Nc2) = ref2(:,1);
    objsig2((i-1)*Nc2+1:i*Nc2) = obj2(:,1);
end

refsig = refsig1+refsig2;
objsig = objsig1+objsig2;


figure(1);
plot(refsig)
hold on
plot(objsig);
ylim([-0.2 4]);
xlabel('ms');
ylabel('Amplitude');
title('Reconstructed signal');
legend('ref','obj')


% generate encode signal
Phi = randi([0 1],N,N);
y=Phi*refsig;
y1=Phi*objsig;
% figure(2);
% plot(y)
% hold on
% plot(y1);

for i=1:M
   position(i,1) = (i-1) *shiftime+1;
   while(position(i) > N)
       position(i) = position(i)-N;
   end
end

%%Adding some measurement noise.
SNR=40;
n=awgn(y1,SNR,'measured');
outputref = zeros(M,1);
outputobj = zeros(M,1);

%Making random measurements
A=zeros(M,N);
for i=1 : M
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
%Compute error recovered
% diff_ref = refsig - xp_ref;
% recovery_error_ref = norm(diff_ref) / norm(refsig);
% fprintf('recovery error: %0.4f\n', recovery_error_ref);
% 
diff = objsig - xp_obj;
recovery_error = norm(diff) / norm(objsig);
fprintf('recovery error: %0.4f\n', recovery_error);
% 
% 
figure(5)
plot(xp_ref)
hold on 
plot(xp_obj)
xlabel('sample');
ylabel('Amplitude');
title('Reconstructed signal');
legend('ref','obj');

t1 = zeros(2,1);
dem=0;
for j = 1: length(xp_ref)
        if(round(xp_ref(j),1) == 2)
            t1(2) = (j-1)/(f1*Nc1);
            dem=dem+1;
        end
        if(round(xp_ref(j),1) == 1)
            t1(1) = (j-1)/(f2*Nc2);
            dem=dem+1;
        end
        if(dem ==2)
            break;
        end
        
        
        
end

t2 = zeros(2,1);
dem2=0;
for j = 1: length(xp_obj)
        if(round(xp_obj(j),1) == 2)
            t2(2) = (j-1)/(f1*10^3*Nc1);
            dem2=dem2+1;
        end
        if(round(xp_obj(j),1) == 1)
            t2(1) = (j-1)/(f2*10^3*Nc2);
            dem2=dem2+1;
        end
        if(dem2 ==2)
            break;
        end        
end

time_delay = t2(1)-t1(1);


LightSpeed = 3*10^8;
Distance = LightSpeed*time_delay;
fprintf('Measured Distance = %.2fm\n',Distance)

