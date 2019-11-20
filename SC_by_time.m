%%%%%%%%%%%%%
%Get data from the IRIS DMC
%Do SPAC and HV measurements for the inversion
%%
clear
clc
close all

addpath('./HVfunctions/');
addpath('./PVrfunctions/');
addpath('./wfTools/')
%%
%%%%%%%%%%%
%Parameters

%%%%%%%%%%%%%
%For the data
start_date    = '2011-02-01 00:00:00';
end_date      = '2011-04-01 00:00:00';
download_inc  = 1;%in days
window_length = 60*60;%in seconds
network       = 'TA';
station       = { '440A' }; %LRL 

eq_time = datenum('2019-07-04 16:33:49');

%work in late 2019 in SoCal
%{'ADO','ARV','BAK','BAR','BBR','BC3','CHF','CIA','CWC','DAN','DEC','DGR','DJJ','EDW2','FMP','FUR','GLA','GMR','GRA','GSC','HEC','IKP','IRM','ISA','LRL','MLAC','MPM','MUR','MWC','NEE2','OSI','PASC','PASC','PDM','PHL','PLM','RPV','RRX','SBC','SCI2','SDD','SHO','SLA','SMM','SVD','TIN','TUQ','USC','VES','VOG','VTV'}

vertical   = 'BHZ';
horz_1     = 'BHN';
horz_2     = 'BHE';

location   = '*';

%%%%%%%%%
%for correlation
onebit         = 1;
nsmooth_denom  = 20;
nsmooth_SC     = 3;%daily stack
low            = 1;
high           = 10;

time_keep = 20;%length of SC to keep
ccwin     = 3;%length over which you correlate
strech_factor = (-25:0.1:25)/100;

ZZ = [];
ZN = [];
ZE = [];

%%
%%Loop over the downloading inc 

disp([ 'Starting on ' start_date ' and working till ' end_date ]);
day_vec = datenum(start_date):(datenum(end_date) - 1);

SC = [];

for i = 1:length(day_vec) %assumes that the first day downloaded correctly
    
    disp([ 'Working on ' datestr(day_vec(i), 2) ])
    
    Z  = irisFetch.Traces(network, station, location, vertical, ...
        datestr(day_vec(i), 31), datestr(day_vec(i) + 1, 31));
    HN = irisFetch.Traces(network, station, location, horz_1, ...
        datestr(day_vec(i), 31), datestr(day_vec(i) + 1, 31));
    HE = irisFetch.Traces(network, station, location, horz_2, ...
        datestr(day_vec(i), 31), datestr(day_vec(i) + 1, 31));
    
    if length(Z) == 1 && length(HN) == 1 && length(HE) == 1
    
        %now check that you actually got stuff on each channel
        z_list  = unique({Z.station});
        h1_list = unique({HN.station});
        h2_list = unique({HE.station});

        kill = zeros(length(z_list), 1);
        for k = 1:length(z_list)

            if ~any(strcmp(h1_list, z_list(k))) || ~any(strcmp(h2_list, z_list(k)))

                kill(k) = 1;

            end

        end

        Z      = Z(~kill);%convets to logical
        z_list = z_list(~kill);

        for k = 1:length(h1_list)

            if ~any(strcmp(z_list, h1_list(k))) || ~any(strcmp(h2_list, h1_list(k)))

                kill(k) = 1;

            end

        end

        for k = 1:length(h2_list)

            if ~any(strcmp(z_list, h2_list(k))) || ~any(strcmp(h1_list, h2_list(k)))

                kill(k) = 1;

            end

        end

        %now check if you need to remove extra samples - not uncommon
        %or kill data that is totally the wrong ass length
        max_length = 0;

        for k = 1:length(z_list)

            max_length = max([ max_length length(Z(k).data)]);
            max_length = max([ max_length length(HN(k).data)]);
            max_length = max([ max_length length(HE(k).data)]);

        end

        %if off by one, clip. If off by more than one, kill
        for k = 1:length(z_list)

            if length(Z(k).data) == (max_length)

                Z(k).data = Z(k).data(1:(max_length - 1));

            elseif length(Z(k).data) < (max_length - 1)

                kill(k) = 1;

            end

            if length(HN(k).data) == (max_length)

                HN(k).data = HN(k).data(1:(max_length - 1));

            elseif length(HN(k).data) < (max_length - 1)

                kill(k) = 1;

            end

            if length(HE(k).data) == (max_length)

                HE(k).data = HE(k).data(1:(max_length - 1));

            elseif length(HE(k).data) < (max_length - 1)

                kill(k) = 1;

            end

        end
        
    else
        
        kill = 1;
        
    end

    %get length of output
    seg_length = window_length*Z.sampleRate;
    win_start     = 1:seg_length:( length(Z.data) - seg_length);
    
    if ~any(kill)
               
        [ SCtmp  ] = SC_function( Z, HN, HE, window_length, onebit, ...
            nsmooth_denom, low, high, time_keep);
               
        [q,w,e,r] = size(SCtmp);
        
        if isempty(SC)
        
            SC = SCtmp;
        
        else
            
            SC(:, end + 1:end + length(win_start), :, :) = SCtmp;
            
        end
            
    elseif any(kill) && i == 1
        
        error('Starting day invalid');
        
    elseif any(kill) && i > 1
        
        SC(:, end + 1:end + length(win_start), :, :) = nan(q, w, e, r);
        
    end
    
end

%save([ experiment_name '-data']);

%%
%%%%%%%%
%do raw correlations

[q,~,~,~] = size(SC);

for k = 1:q
   
    for i = 1:3
        
        for ii = 1:3
    
            SCsmooth(k, :, i, ii) = downsample(smoothdata(SC(k, :, i, ii), 2, 'movmean', ...
                nsmooth_SC), nsmooth_SC);

        end
        
    end
    
end

[q,w,e,r] = size(SCsmooth);

for k = 1:w
   
    for i = 1:3
        
        for ii = 1:3
    
            SCsmooth(:, k, i, ii) = SCsmooth(:, k, i, ii)/rms(SCsmooth(:, k, i, ii));

        end
        
    end
    
end

start_date = datenum(start_date);
dateaxis   = linspace(datenum(start_date),datenum(end_date), w); 

t     = linspace(0, time_keep, time_keep*Z.sampleRate);

% %%%%%%%%
% %now calculate the changes in seismic velocity
% preeventZN = nanmean(ZNsmooth(:, eq_time>dateaxis), 2);
% preeventZE = nanmean(ZEsmooth(:, eq_time>dateaxis), 2);
% 
% figure(2)
% hold on
% plot(t, preeventZN);
% plot(t, preeventZE);
% legend('Preevent ZN', 'Preevent ZE');
% 
% for k = 1:length(dateaxis)
%    
%     %first, align up up to two sample misalignments - not related to
%     %streching
%     [ xcZN, lags ] = xcorr(preeventZN, ZNsmooth(:, k), 2);
%     [~, ind] = max(xcZN);
%     ZNsmooth(:, k) = circshift(ZNsmooth(:, k), lags(ind));
%     
%     [ xcZE, lags ] = xcorr(preeventZE, ZEsmooth(:, k), 2);
%     [~, ind] = max(xcZE);
%     ZEsmooth(:, k) = circshift(ZEsmooth(:, k), lags(ind));
%     
%     for kk = 1:length(strech_factor)
%        
%         ZNstr = interp1(t, ZNsmooth(:, k), linspace(0, max(t)*(1 + strech_factor(kk)), length(t)));
%         ZEstr = interp1(t, ZEsmooth(:, k), linspace(0, max(t)*(1 + strech_factor(kk)), length(t)));
%         
%         ccZNtmp = corrcoef(preeventZN(ccwin > t), ZNstr(ccwin > t));
%         ccZEtmp = corrcoef(preeventZE(ccwin > t), ZEstr(ccwin > t));
% 
%         ccZNsearch(kk) = ccZNtmp(1,2);
%         ccZEsearch(kk) = ccZEtmp(1,2);
%         
%     end
%         
%     [ ccZN(k), ind ] = max(ccZNsearch);
%     epZN(k) = strech_factor(ind);
%     
%     [ ccZE(k), ind ] = max(ccZEsearch);
%     epZE(k) = strech_factor(ind);
%     
%     ep(k)  = (epZN(k)*(ccZN(k)^2) + epZE(k)*(ccZE(k)^2))/(ccZN(k)^2 + ccZE(k)^2);
%     ccM(k) = ((ccZN(k)^3) + (ccZE(k)^3))/(ccZN(k)^2 + ccZE(k)^2);
%     
% end

figure(30)
subplot(3, 3, 1)
imagesc(dateaxis, t, SCsmooth(:, :, 1, 1))
datetick('x','mm-dd HH:MM')
xlim([ min(dateaxis) max(dateaxis) ])
%caxis([ -max(abs(ZNsmooth(:))) max(abs(ZNsmooth(:))) ])

subplot(3, 3, 2)
imagesc(dateaxis, t, SCsmooth(:, :, 1, 2))
datetick('x','mm-dd HH:MM')
xlim([ min(dateaxis) max(dateaxis) ])

subplot(3, 3, 3)
imagesc(dateaxis, t, SCsmooth(:, :, 1, 3))
datetick('x','mm-dd HH:MM')
xlim([ min(dateaxis) max(dateaxis) ])

subplot(3, 3, 4)
imagesc(dateaxis, t, SCsmooth(:, :, 2, 1))
datetick('x','mm-dd HH:MM')
xlim([ min(dateaxis) max(dateaxis) ])

subplot(3, 3, 5)
imagesc(dateaxis, t, SCsmooth(:, :, 2, 2))
datetick('x','mm-dd HH:MM')
xlim([ min(dateaxis) max(dateaxis) ])

subplot(3, 3, 6)
imagesc(dateaxis, t, SCsmooth(:, :, 2, 3))
datetick('x','mm-dd HH:MM')
xlim([ min(dateaxis) max(dateaxis) ])

subplot(3, 3, 7)
imagesc(dateaxis, t, SCsmooth(:, :, 3, 1))
datetick('x','mm-dd HH:MM')
xlim([ min(dateaxis) max(dateaxis) ])

subplot(3, 3, 8)
imagesc(dateaxis, t, SCsmooth(:, :, 3, 2))
datetick('x','mm-dd HH:MM')
xlim([ min(dateaxis) max(dateaxis) ])

subplot(3, 3, 9)
imagesc(dateaxis, t, SCsmooth(:, :, 3, 3))
datetick('x','mm-dd HH:MM')
xlim([ min(dateaxis) max(dateaxis) ])
% figure(50)
% scatter(dateaxis, 100*ep, 20, ccM, 'filled', 'MarkerEdgeColor', 'k')
% grid on
% colormap(flipud(hot))
% h = colorbar; h.Label.String = 'cc';
% datetick('x','mm-dd-yy')
% xlabel('Date')
% ylabel('Change in velocity, %');
%xlim([ 737601 737637 ])
%hold on
plot([datenum('2011-03-11 17:00:00') datenum('2019-07-04 00:00:00')], [ -50 50 ], 'k--')
plot([datenum('2011-03-11 03:00:00') datenum('2019-07-06 00:00:00')], [ -50 50 ], 'k--')
%ylim([ -20 2])

save([network '.' station{:} '.SC.mat' ]);
