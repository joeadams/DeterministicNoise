function [xMean,t] = integration(h,tend,N,fmin,fmax,famp,fmod,K)

    %inputs
    %h - integration step
    %tend - final time of integration
    %N - number of oscillators
    %fmin - minimum natural frequency of ensemble
    %fmax - maximum natural frequency of ensemble
    %famp - modulation amplitude of oscillators
    %fmod - modulation frequency of oscillators
    %K - coupling strength of network
    
    %output
    %xMean - integrated series
    %t - corresponding time

    t = 0:h:tend;
    x=zeros(N,1);
    fnat=zeros(N,1);
    for i=1:N
        x(i,1) = (i-1)*2*pi/(N-1);
        fnat(i) = fmin+(i-1)*(fmax-fmin)/(N-1)+rand(1,1);
    end
    n = length(t)-1;
    x_dot =@(fnat,famp,t,fmod,x,K,N) 2*pi*(fnat+fnat*famp*sin(2*pi*t*fmod)+K/N*sum(sin(x.*ones(1,N)-(ones(N,1).*x')))');
    for i = 1:n
         
          k11 = x_dot(fnat,famp,t(i),fmod,x(:,i),K,N);
          k12 = x_dot(fnat,famp,t(i)+.5*h,fmod,x(:,i)+.5*h*k11,K,N);
          k13 = x_dot(fnat,famp,t(i)+.5*h,fmod,x(:,i)+.5*h*k12,K,N);
          k14 = x_dot(fnat,famp,t(i)+h,fmod,x(:,i)+h*k13,K,N);
          x(:,i+1) = x(:,i)+((k11+2*k12+2*k13+k14)/6)*h;
          
    end
    xMean = 1/N*sum(sin(x));
end
