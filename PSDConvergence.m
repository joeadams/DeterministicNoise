function [B,C] = PSDConvergence(x,Fs,freqmin,freqmax,period)

%inputs
%x - signal
%Fs - sampling frequency of signal
%freqmin - minimum frequency to be included in PSD
%freqmax - maximum frequency to be included in PSD
%period - period being examined

%outputs
%B - cell of PSD variables at each signal length
%C - cell of PSD line of best fit gradients at each signal length

    Tfinal = [1:1:10,20:10:100,200:100:1000,2000:1000:10000]*period;
    B=cell(length(Tfinal));
    xT=cell(length(Tfinal));
    for i=1:length(Tfinal)
        xT{i} = x(1:Tfinal(i)*Fs+1);
    end
    clear x
    
    for i=1:length(Tfinal)
        t = 0:1/Fs:Tfinal(i);
	    n=2^(nextpow2(length(t))+1);

        xdft = fft(xT{i},n);
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
    
        A = {xR0,PSD_aver,cond2,p2,beta2};
        B{i}=A;
	    C{i}=beta2;
    end

end