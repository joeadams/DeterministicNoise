function [xcor,dt,acor] = autocorrelation(x,windowperiods,rangeperiods,freq,samplingfreq)

%inputs
%x - signal
%windowperiods - length of window in number of periods
%rangeperiods - range of lags in number of periods
%freq - frequency being examined
%samplingfreq - sampling frequency of signal

%outputs
%xcor - lags
%dt - lag scaling factor
%acor - autocorrelation

fs = 1/samplingfreq;                 % sampling rate
ACstep=samplingfreq;              % autocorrelation function is not not evaluated 
                         % with the delay equal to the sampling time but  
                         % with ACstep, which  must be greater or equal 1/fs; 
                         % it is useful in order to prevent too long
  	                     % computation
window = ceil(windowperiods*1/freq);

Delay_min = -ceil(rangeperiods*1/freq);         % Range of the delays  for which we plot 
Delay_max = ceil(rangeperiods*1/freq);          % the autocorrelation function

%%%%%%%%%%%%%%%%%%%
%close all

dt = 1/fs;
DN = floor(ACstep*fs);
Ds_min = floor(Delay_min*fs/DN);
Ds_max = floor(Delay_max*fs/DN);
WindowLength = floor(window*fs/DN);

if size(x,2)>1
    x=x';
end
x = (x-min(x))./(max(x)-min(x));

N=length(x);                        % number of samples 
 
% Finding the mid/point of the time series
 
if N-floor(N/2)*2==0
    N=N-1;
end
N=floor(N/2);
 
% Setting the window to 1/2 of the data length
 
M=N/2;
M=floor(M);
 
iter=floor((N-WindowLength/2)/DN);
if Ds_min <= -iter
	Ds_min = -iter;
end
if Ds_max >= iter
	Ds_max = iter;
end

acor=[];
xcor=[];

for i=Ds_min:Ds_max
    i1= N+1-WindowLength/2+i*DN;
    i2=i1+WindowLength;
    x1=x(i1:i2);
    x0=x(N+1-WindowLength/2:N+1+WindowLength/2);
    aa=x1'*x0/sqrt(x1'*x1)/sqrt(x0'*x0);
    acor=[acor; aa];
    xcor=[xcor; i*DN];
end

end

