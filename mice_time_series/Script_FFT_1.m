L = 10000000; % length of the input signal to be used for FFT.
Trt = 2000000;  %this is the cutoff transient time   
n = 2^nextpow2(L-Trt); % 

Fs = 1000; % sampling frequency.
load Mat2;% you can load Mat2 and Mat 3 similarly.
Mat_temp = Mat2(1:L);
Y=fft(Mat2(Trt:n));
f = Fs*(1:(n/2))/n;
P2=abs(Y/n);
P1 = P2(1:n/2+1);
P1(2:end-1) = 2*P1(2:end-1);


%[xmin xmax ymin ymax]

%plot(f,P1(2:n/2+1));
%plot(f,P1(3:n/2));
%plot(f(2:n/2+1),P1(2:n/2+1));
plot(f(1:1000000),P1(2:1000001))

M2 = P1;
save M2;

