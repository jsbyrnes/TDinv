function [yhat] = Kolmogorov(x,input,output)
%
% Minimum phase reconstruction from amplitude spectrum
% using spectral facorization (Claerbout '92)
%
% This function works along the column space of Signal
%
% Input and Output are flags for the form of the data
%
% Input options are 'freq' and 'time'
%
% Output options are 'time' and 'freq'
%
% Default for output is 'time', 'freq' input

if isempty(output)
    output='time';
end

if isempty(input)
    input='freq';
end

n = size(x,1);

switch lower(input)
    case 'time'
        spec=abs(fft(x));
    case 'freq'
        spec=abs(x);
    otherwise
        error('Not a valid input type, see documetation')
end
        

% if strcmpi(input, 'time')
%     xhat = real(ifft(log(abs(fft(x)))));
% else if strcmpi(input, 'freq')
%         xhat = real(ifft(log(abs((x)))));
%     else
%         error('Not a valid input type, see documetation')
%     end
% end

% check for zeros
if any(spec(:) == 0)
    fprintf('***Zero in the spectrum detected! Using eps... (Kolmogorov.m)')
    
    ipz = spec == 0;
    spec(ipz) = eps;
end

xhat = real(ifft(log(spec)));

odd = fix(rem(n,2));
wn = [1; 2*ones((n+odd)/2-1,1) ; ones(1-rem(n,2),1); zeros((n+odd)/2-1,1)];
wn = repmat(wn,1,size(x,2));
yhat = zeros(size(x));
yhat(:) =((exp(fft(wn.*xhat))));

if strcmpi(output,'time')
    yhat=real(ifft(yhat));
end
