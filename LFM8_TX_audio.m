

% TheGmr140 YouTube Channel
% Your Source for Real Radio Projects

% 8 LFM transmitter
% read text file
% create I/Q wave file for transmission

close all
clear all

% define TX message burst
numchar = 256;  % limit text message to 1024 char, one char = 8bits
maxbits = numchar * 8;  %max number of bits in message burst
bitspersymbol = 3;
Nsymbols = round(maxbits/bitspersymbol);

fs = 4e3;  %sample rate of file
rg = 1/fs;  %sample period of file





btime = 30;   % guard for symbol



% create preamble chirp waveform
pw = 500*rg;   % preamble time ,  500 range gates
bw = fs;  % bandwidth of chirp ,  
t = [rg:rg:pw];
t = t - pw/2;
slope = bw / (pw);
preamble = exp(1i*pi*slope*t.^2);


% create 8 LFM signal waveforms
fftsize = 128;   %keep power of 2
pwc = fftsize * rg;  %pw of lfm signal
ntones = 8;  %keep power of 2
bw1 = ([1:ntones]/ntones) * fs;
slope = bw1 / pwc;  % 1x8 different slopes

% make 8 lfm signals
t = [rg:rg:pwc];
t = t - pwc/2;
bb = fftsize + btime;
sww = zeros(ntones,bb);   % ntones x symbolsize
for k = 1:ntones
    xa = exp(-1i*pi*slope(k)*t.^2);
    sww(k,:) = [xa zeros(1,btime)];   %ntones x fftsize+btime
end



[~,stime] = size( sww );    %ntones x symbolsize






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% read in text file
[filename, pathname] = uigetfile('*.txt','Pick a Text File','c:\8LFM');

pathname = [pathname filename];
fid = fopen(pathname);
g = fread(fid);  % read in text file, reads in 8bit words
b = length(g);
sprintf('The max number of char in burst: %d',numchar)
sprintf('The number of char in file:  %d',b)
if b >= numchar
    g = g(1:numchar);
else
    x = numchar - b;
    g = [g ; zeros(x,1) ]; 
end
g = dec2bin(g,8);  % convert 8bit words to binary string
g = g';
g = reshape(g,1,[]);  %  '101011010101010101010........
z = floor(  length(g) / bitspersymbol);
g = g(1:(z*bitspersymbol));
fclose(fid);
g = reshape(g,bitspersymbol,[]);
g = g.';              % 100
                      % 100
                      % 110
[numbits,~] = size(g);



% define matrix for tx waveforms
iqdata = zeros(numbits,stime);

% build tx matrix for 8 level fsk waveforms
    for k = 1:numbits

            if g(k,:) == '000'
                iqdata(k,:) = sww(1,:);
            end

            if g(k,:) == '001'
                iqdata(k,:) = sww(2,:);
            end

            if g(k,:) == '010'
                iqdata(k,:) = sww(3,:);
            end

            if g(k,:) == '011'
                iqdata(k,:) = sww(4,:);
            end
            
            if g(k,:) == '100'
                iqdata(k,:) = sww(5,:);
            end
            
            
            if g(k,:) == '101'
                iqdata(k,:) = sww(6,:);
            end
            
            
            if g(k,:) == '110'
                iqdata(k,:) = sww(7,:);
            end
            
            if g(k,:) == '111'
                iqdata(k,:) = sww(8,:);
            end

    end



iqdata = reshape(iqdata.',1,[]);    %  1xN complex vector


message = [preamble iqdata ];

x = 5*fs;
message = [  zeros(1,x) message zeros(1,x)];
hlpf = fir1(24,0.92);
message = filter(hlpf,1,message);

message = message/max(abs(message));
% nS  = length(message);
% noise = 10^(-60/20) * (randn(1,nS) + 1i*randn(1,nS));
% 
% message = 10^(-50/20)*message + noise;
m = resample(message,2,1);
n = length(m);
fs2 = fs * 2;
t = [1:n]/fs2;
f4 = fs2/4;
message = m .* exp(1i*2*pi*f4*t);
message = real(message) + imag(message);
message = real(message);
message = message / max(message);


message = [real(message)'];

% strname = sprintf('8LFM_%dKHZ_IQ.wav',fs/1e3);
% x = ['c:\8LFM\'  strname ];
% 



audiowrite('c:\temp\audio8lfm256char.wav',message,fs2);













        