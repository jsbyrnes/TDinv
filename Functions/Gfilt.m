function [A] = Gfilt(Fo,Fa,n,sps)
% [A] = Gfilt(Fa,Fa,n,sps)
% A is amplitude spectrum of Gussian filter
% n is number of points in time series
% Fo is freq where amplitude is max
% Fa is freq width where amplitude is down by half
% sps is samples per sec (ie, 1/(time interval))
% --- if sps left off, Fa is in freq points
%
% --- to filter time series X, ifft(Gfilt(fft(X)))
%set the optional argument to 1 if using a gpu
if nargin==4
  FreqPerPoint = sps/(n+1); % freq step
  Fa=Fa/FreqPerPoint;
  Fo=Fo/FreqPerPoint;
end

Fc=Fa/sqrt(log(2));  % Fc is at Amp = 1/(e^2)
n2= ceil(n/2);    
A = zeros(n,1);
nn = ((0:n2)-Fo)/Fc;
A(1:n2+1) = exp(-nn.^2);
A(n:-1:n2+1) = A(2:n2+1);
