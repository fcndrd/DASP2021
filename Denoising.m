% DASP Single Channel
clear;
clear sound;
rng(1)
[x,fs]=audioread('clean_speech.wav');
[noise,fs]=audioread('babble_noise.wav');
fplot=0; %plot flag

SNR=10^(10/10); %linear SNR
%noise=randn(length(x),1); % Gaussian Noise
noise=noise(1:length(x));
comp=0;  %More complex algorithm flag (bias compensation and optimized \alpha)

SNRr=norm(x).^2/norm(noise).^2;
noise=noise*sqrt(SNRr/SNR);
% NP=norm(noise).^2
% SP=norm(x).^2
y=x+noise;
% sound(y,fs);
%% framing


K=512;  %512 samples = 32ms with fs=16Khz
overlap=0.5;
Y=zeros(K,floor(length(y)/K));
i=1;
for l=1:K*(1-overlap):length(y)-K
    Y(:,i)=(fft(y(l:l+K-1).*(hann(K)))).';
    i=i+1;
end
L=size(Y,2);
K=K/2+1;
Y=Y(1:K,:);

%% SIGNAL PSD ESTIMATE
Pyy=zeros(K,L);
Pyy=abs(Y).^2;
Pyyb=Pyy;
M=32; % bartlett average length 0.2ms @16Khz

%Bartlett averaging
    %Transient
for i=1:M
    Pyyb(:,i)=1/i*sum(Pyy(:,1:i),2);
end
for i=M:L %Bartlett averaging
Pyyb(:,i)=1/M*sum(Pyy(:,(i-M+1):i),2);
end

%Exponential smoothing
Pyyex=Pyy;
alfa=0.85;

if alfa~=1   
for i=2:L
    Pyyex(:,i)=alfa*Pyyex(:,i-1)+(1-alfa)*Pyy(:,i);       
end
end   

if fplot   
figure;
hold on
plot(Pyyex(100,:));
plot(Pyy(100,:));
legend('Pyy smoothed','Pyy')
end

%% Minimum Statistics (No bias correction)
    

    D=100; %100*32/2*ms =1600 ms
    Pnn=zeros(K,L);

    
    for i=1:D
    Pnn(:,i)=min(Pyyex(:,1:i),[],2);
    end
    for i=D:L
    Pnn(:,i)=min(Pyyex(:,i-D+1:i),[],2);
    end

    
    
%Plotting
if fplot
figure;
hold on;
plot(Pyyex(100,:));
plot(Pnn(100,:));
legend('Pyy','Pnn')
end

%% Complete Mimimum Statistics

if comp==1
Pnn=noise_est_min_stat(Pyy,D);
end
%% GAIN CALCULATION
G=zeros(K,L);

G=1-(Pnn./Pyyb);
Yout=G.*Y;

Yout=[Yout; conj(flip(Yout(2:end-1,:)))]; %Reconstruct full Y through simmetry

%% IFFT
S=ifft(Yout);

%% Overlap and add
s=zeros(length(y),1);
i=1;
K=2*K-2;
for l=1:K*(1-overlap):length(y)-K
    s(l:l+K-1)=s(l:l+K-1)+S(:,i);
    i=i+1;
end

%% Output and Evaluation
clear sound;

% sound(s,fs); 
MSEpre=mean((y - x).^2);
MSEpost=mean((s - x).^2);
STOIpre=stoi(x,y,fs);
STOIpost=stoi(x,s,fs);


