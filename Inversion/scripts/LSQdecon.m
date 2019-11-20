function  [G S0] = LSQdecon(C1,C2,C3,Source,len, dt, shift, damping_factor)
%LSQdecon Does deconvolution for a set of waveforms. Format is [trace
%nevents]
npts=len;%/dt + 1;
[~, nevt] = size(C1);

filt_norm = 1;

Time=(0:npts-1)*dt;
df=1/(npts*dt);
Freq=Time./dt.*df;

% fft
FFT.C1=fft(C1,[],1);
FFT.C2=fft(C2,[],1);
FFT.C3=fft(C3,[],1);

% make sure there are no zeros
FFT.C1(FFT.C1==0)=eps;
FFT.C2(FFT.C2==0)=eps;
FFT.C3(FFT.C3==0)=eps;

if isempty(Source) % just use C1
    
    for qq = 1:nevt
        
        Source(qq).spec_est       = abs(FFT.C1(:, qq));
        Source(qq).min_phase_filt = Kolmogorov(FFT.C1(:, qq),'freq','freq')./FFT.C1(:, qq);
        
        if filt_norm % normalize
            Source(qq).min_phase_filt = Source(qq).min_phase_filt ./ abs(Source(qq).min_phase_filt);
        end
        
    end
    
elseif length(Source) ~= nevt
    
    error('Wrong number of events or Sources in LSQDecon');
       
end

% invert for amplitude spectrum
%================================
% create A and B matrix
B=zeros(nevt*4,npts,class(C1));

A=repmat([eye(3);[0,0,0]],nevt,1);
A=[A,zeros(nevt*4,nevt)];

for nn=1:nevt
        
    B(nn*4-3,:)=log(abs(FFT.C1(:,nn))).';
    B(nn*4-2,:)=log(abs(FFT.C2(:,nn))).';
    B(nn*4-1,:)=log(abs(FFT.C3(:,nn))).';
    %B(nn*4,:)=log(localsource(nn).spec_est).';
    B(nn*4,:)=log(Source(nn).spec_est).';
    
    A((nn-1)*4+1:(nn)*4,3+nn)=[1 1 1 1]';
end

% now invert!!
% using least squares

theta = eye(size(A'*A))*damping_factor;

Ainv=(A'*A + theta)\A'; % see mldivide for algorithm

M=zeros(nevt+3,npts,class(C1));
for nn=1:npts
    M(:,nn)=Ainv*B(:,nn);
end

% amp spec
G.FFT_C1=exp(M(1,:)).';
G.FFT_C2=exp(M(2,:)).';
G.FFT_C3=exp(M(3,:)).';

S0=exp(M(4:end,:)).'; % retrieve LSQR source spec estimate

% calculate and average all-pass filters

source_filt=[Source.min_phase_filt];

F.C2=(FFT.C2.*source_filt)./Kolmogorov(FFT.C2,'Freq','Freq');
F.C3=(FFT.C3.*source_filt)./Kolmogorov(FFT.C3,'Freq','Freq');

% create psedo all pass filters
F.C2=mean(F.C2,2);
F.C3=mean(F.C3,2);

%*********************************
% semi-normalize filters
if filt_norm == .5
    F.C2=F.C2./max(abs(F.C2));
    F.C3=F.C3./max(abs(F.C3));
end
%*********************************

% normalize filters
if filt_norm == 1
    F.C2=F.C2./abs(F.C2);
    F.C3=F.C3./abs(F.C3);
end

% it is possible to get NaNs by zero division or Kolmogorov
if any(isnan(F.C2))
   ip=find(isnan(F.C2));
   F.C2(ip)=ones(size(ip));
   warning(sprintf('NaN found in F.C2 at site : %d\n',ip)) 
end

if any(isnan(F.C3))
   ip=find(isnan(F.C3));
   F.C3(ip)=ones(size(ip));
   % warning(sprintf('NaN found in F.C3 at site : %d\n',ip))   
end

% time shift in sec

f_shift = exp(sqrt(-1).*Freq.*2*pi*(shift))';
%f_shift = exp(sqrt(-1).*Freq.*2*pi*shift)';

% reconstruct greens
one=ones(size(G.FFT_C1),class(C1)); % make sure C1 stays in the same class

G.FFT_C1=Kolmogorov(G.FFT_C1,'Freq','Freq').*one .*f_shift;
G.FFT_C2=Kolmogorov(G.FFT_C2,'Freq','Freq').*F.C2.*f_shift;
G.FFT_C3=Kolmogorov(G.FFT_C3,'Freq','Freq').*F.C3.*f_shift;

%G.FFT_C2=Kolmogorov(G.FFT_C2,'Freq','Freq').*F.C2;
%G.FFT_C3=Kolmogorov(G.FFT_C3,'Freq','Freq').*F.C3;

G.C1=real(ifft(G.FFT_C1));
G.C2=real(ifft(G.FFT_C2));
G.C3=real(ifft(G.FFT_C3));