function [xR0,PSD_aver,cond2,p2,beta2] = PSD(t,xMean,Fs,freqmin,freqmax)

    %inputs
    %t - time vector corresponding to signal
    %xMean - signal vector
    %Fs - sampling frequency of xMean
    %freqmin - minimum frequency to be included in PSD
    %freqmax - maximum frequency to be included in PSD
    
    %outputs
    %xR0 - frequency spectrum
    %PSD_aver - power corresponding to each frequency
    %cond2 - vector of included frequencies
    %p2 - polynomial of line of best fit
    %beta2 - gradient of line of best fit

    n=2^(nextpow2(length(t))+1);
    xdft = fft(xMean,n);
    xdft = xdft(1:ceil(n/2)+1);
    psdx = (1/(Fs*n)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    Freq = 0:Fs/n:Fs/2;
    xR0=log10(Freq);
    xR0(1)=NaN;
    yR0=log10(psdx);
    yR0(1)=NaN;
    PSD_aver = yR0;
    log_fX0=log10(freqmin);
    log_fX=log10(freqmax);
    cond2=(xR0>log_fX0 & xR0<log_fX);
    p2=polyfit(xR0(cond2),PSD_aver(cond2),1);
    beta2=-p2(1);
    
end