function [result, fft_result]=mtmdeconSP(T1, T2, n, ntap, nmtw, NW)
%T1 is the source,T2 is decon time series,n is length of ts,ntap is window length
%nmtw is the number of tapers, NW is the time bandwidth number (N) ntap time (W) frequency resolution

ts1(1,:) = T1;
ts2(1,:) = T2;

[mm,nn]=size(ts1);
nseg=length(0:ntap/4:n-ntap);
[E,V]=dpss(ntap,NW,nmtw);
tap1=zeros(nseg,n,nmtw);%tap1s=tap1;

for j=1:length(V)  %%%loop around each taper
    for i=1:nseg   %%%loop aroung each time delay
        tap1(i,[1:ntap]+ntap*(i-1)/4,j)=E(:,j);
    end
end

ts1b = zeros(n,nmtw);
ts2b = zeros(n,nmtw);

for i=1:mm
    [tts1,jk]=meshgrid(ts1(i,:),1:nseg);
    [tts2,jk]=meshgrid(ts2(i,:),1:nseg);
    
    for k=1:length(V)
        
        ts1b(:,k)=sum(fft(tap1(:,:,k).*tts1*V(k),[],2),1);
        
        ts2b(:,k)=sum(fft(tap1(:,:,k).*tts2*V(k),[],2),1);
        
        %Nois(:,k)=sum(fft(tap1(nseg,:,k).*tts2(nseg, :)*V(k),[],2),1);
        
        %Nois(:,k)=(fft(tap1(1,:,k).*tts2(1,:)*V(k),[],2));
        
    end
    
    size(ts1b);
%   top=sum(ts2b,2);
%   bottom=sum(ts1b,2);

    top=sum(ts2b.*conj(ts1b),2);
    bottom=sum(ts1b.*conj(ts1b),2);
    %S=0;
    
    %top = cosfiltfreq(nn, top, 1.5, .05);
    %bottom = cosfiltfreq(nn, bottom, 1.5, .05);

    %S=sum(Nois.*conj(Nois),2);
    size(top);
    size(bottom);
    
    %ts1b= abs(fft(tts1.*tap1*V(1),[],2));
    %ts2b=fft(tts2.*tap1*V(1),[],2).*conj(fft(tts1.*tap1s*V(1),[],2))+fft(tts2.*tap2*V(2),[],2).*conj(fft(tts1.*tap2s*V(2),[],2))+fft(tts2.*tap3*V(3),[],2).*conj(fft(tts1.*tap3s*V(3),[],2))+fft(tts2.*tap4*V(4),[],2).*conj(fft(tts1.*tap4s*V(4),[],2))+fft(tts2.*tap5*V(5),[],2).*conj(fft(tts1.*tap5s*V(5),[],2))+fft(tts2.*tap6*V(6),[],2).*conj(fft(tts1.*tap6s*V(6),[],2))+fft(tts2.*tap7*V(7),[],2).*conj(fft(tts1.*tap7s*V(7),[],2));;
    %top=sum((ts2b),1);
    %bottom=sum(ts1b,1);
    
    fft_result = top./(bottom);
    
 %   fft_result = cosfiltfreq(nn, fft_result, 2, .05);
    
    result = real(ifft(fft_result));%+S
end