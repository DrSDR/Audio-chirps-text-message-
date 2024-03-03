

% TheGmr140 Youtube Channel
% your source for Real World Radio Projects
% read in wave IQ file from SDR radio 
% be sure sdr is tuned to center of TX freq
% decode wave file 8 lfm to text message
% text message is 256 char message burst
% text message signal will contain a chirp, and 256 char




%clear and close all the shit
close all
clear all

% debug plots enable
% set to one for debug plots
plotdebug = 0;


% define TX message burst
numchar = 256;         % limit text message to 256 char, one char = 8bits
maxbits = numchar * 8; % max number of bits in message burst
nbitspersymbol = 3;
bitpairs = round( maxbits / nbitspersymbol );
bits = zeros(bitpairs,nbitspersymbol);



fs = 4e3;         %sample rate of file
fs2 = fs*2;
f4 = fs2 / 4; 
rg = 1/fs;        %sample period of file
fftsize = 128;



btime = 30;

symsize = fftsize + btime;
nsamples = symsize * bitpairs;




% create preamble chirp waveform
pw = 500*rg;   % preamble time duration, 1/2 sec
bw = fs;  % bandwidth of chirp ,
t = [rg:rg:pw];
t = t - pw/2;
slope = bw / (pw);
preamble = exp(1i*pi*slope*t.^2);
hdet = conj( preamble(end:-1:1));  % chirp detect filter 



% create 8 lfm signal waveforms

pwc = fftsize * rg;
ntones = 8;  %keep power of 2
bw1 = ([1:ntones]/ntones) * fs;
slope = bw1 / pwc;  % 1x8 different slopes

% create 8 lfm signals
t = [rg:rg:pwc];
t = t - pwc/2;

sww = zeros(ntones,fftsize);   % ntones x symbolsize
for k = 1:ntones
    xa = exp(-1i*pi*slope(k)*t.^2);
    xa = conj( xa(end:-1:1) );
    sww(k,:) = xa; 
end






% read the  iq wave file
% [filename, pathname, filterindex] = uigetfile('*.*','Pick a FSK IQ wave file','c:\FSK4Level\');
[filename, pathname, filterindex] = uigetfile('*.*','Pick a FSK IQ wave file','c:\8LFM');

pathname = [pathname filename];
[message,fswave] = audioread(pathname);
[audiosamples,nch] = size(message);
% if nch == 2
%     message = message(:,1) + 1i*message(:,2);
%     message = message.';
% %     message = message / max(message);
% else
%     message = message';
% %     message = message / max(message);
% end
message = message(:,1);
message = message.';
m = hilbert(message);
n = length(m);
tn = [1:n]/fswave;
message = m .* exp(-1i*2*pi*f4*tn);



if fswave ~= fs
    
    x = gcd(fswave,fs);
    a = fs/x;
    b = fswave/x;
    message = resample(message,a,b);
end


% figure(1234)
% plot(abs(message))



% find preamble chirp signal
x = filter(hdet,1,message);
[~,imax] = max(abs(x));

if plotdebug == 1
    figure(1);
    whitebg(1,'k');
    plot(abs(x));
    title('Chirp preamble detection');
end

figure(1863)
r1 = imax;
r2 = imax + length(t);
plot(abs(x(r1:r2)));



istart = (imax + 1)  ;
iend = istart + nsamples ;
iqdata = message(istart:iend);
clear message
nloops = floor( length(iqdata) / symsize );




x1 = 1;
x2 = symsize;

% figure(200)
% plot( abs(iqdata(1:4*symsize)));
iqplot = zeros(nloops,fftsize);

sdet = zeros(1,ntones);

nstop = length(iqdata);


for k = 1:nloops
    
    if x1 >= nstop || x2 >= nstop
        break
    end
    
    iqk = iqdata(x1:x2);
    
    if plotdebug == 1
        figd = figure(99);
        set(figd,'Position',[50 50 1080 720]);
        
        whitebg(figd,'k');
        subplot(2,1,1)
        plot(abs(iqk),'g');
        title('FFT Data Frame');
        subplot(2,1,2)
        plot(real(iqk),'c')
        
        pause(0.001)
    end
    
    
    
    maxsig = zeros(1,ntones);
    maxindex = zeros(1,ntones);
    
    % cross corr each lfm across symbol data
    % find max output and resulting index sample
    for kk = 1:ntones
        bx = sww(kk,:);
        xdet = filter(bx,1,iqk);
        
%         if plotdebug == 1
%             figure(56)
%             set(56,'Position',[100,50,1080,720]);
%             whitebg(56,'k');
%             plot(abs(xdet))
%             title('Chirp compression')
%             axis([0 300 0 10])
%             pause(0.01)
%         end
        
        [ms,mi] = max(abs(xdet));
        maxsig(kk) = ms;
        maxindex(kk) = mi;
    end
    
    
   
        
    
    
    [maxvalbit,imaxbit] = max( maxsig );
    indexoffset = maxindex(imaxbit);
    
    if imaxbit == 1
        bits(k,:) = [0 0 0];
    end
    
    if imaxbit == 2
        bits(k,:) = [0 0 1];
    end
    
    if imaxbit == 3
        bits(k,:) = [0 1 0];
    end
    
    if imaxbit == 4
        bits(k,:) = [0 1 1];
    end
    
    
    if imaxbit == 5
        bits(k,:) = [1 0 0];
    end
    
    if imaxbit == 6
        bits(k,:) = [1 0 1];
    end
    
    if imaxbit == 7
        bits(k,:) = [1 1 0];
    end
    
    if imaxbit == 8
        bits(k,:) = [1 1 1];
    end
    
    
    

    
    x1 = x1 + indexoffset + btime;
    x2 = x1 + symsize;
    
    
    
    
end




       
bits = reshape(bits.',1,[]);  %bits 101010111101110111111111......


n8bits = floor(length(bits) / 8);
bits = bits(1: (n8bits * 8));


bits = reshape(bits,8,[]);
bits = bits';
w = [7 6 5 4 3 2 1 0];

w = 2.^(w);
w = w';
bits = bits*w;
bits = bits';
clc;
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
char(bits)
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

h = msgbox(char(bits),'replace');
% set(h,'FontSize',20);

% msgbox ( char(bits) );
clear bits
