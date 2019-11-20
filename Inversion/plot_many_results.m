%plot the results of a chain

clear, close all
addpath('./scripts')
addpath('./mat_disperse')

vsaxis        = 0:10:3000;
curveaxis     = -10:0.1:10;
vpvsaxis      = 1.5:0.01:4;
HVaxis        = 0.2:0.1:3;
ZJ0axis       = -1:0.1:1;

log_plot_models = 1;

if log_plot_models

    vs_plot   = -4:0.25:-2;
    vpvs_plot = -1:0.25:1;
    curve_plot = -4:0.25:-2;
    
else
   
    vs_plot    = linspace(0.001, 0.01, 10);%-3:0.25:0;
    curve_plot = linspace(0.001, 0.01, 10);%-3:0.25:0;
    vpvs_plot  = linspace(0.01, 0.5, 10);%-3:0.25:0;

end

%name = './4Joe/IRISH_test3v1inversion1';
name = './TD_synBuriedLayerinversion';

inv = define_search;

modelfiles = dir([name '*']);
allmodels = [];

for k = 1:length(modelfiles)
   
    disp( ['Loading file ' modelfiles(k).name ', please wait']);
    load(modelfiles(k).name)

    for kk = 1:length(modelhist)
    
        chain = modelhist{kk};%remove starting position

        chain = chain(logical([chain.converged]));
        
        allmodels = [ allmodels chain ];
    
    end
    
end

%%%%%%%%%%%%%%%%%
%temp line to fix the noise on ZJ0
%ZJ0.value                = normrnd(ZJ0.value, ZJ0.error);
%ZJ0.value(ZJ0.value > 1) = 1;

allmodels = allmodels(:);

fits      = [allmodels.nfit];
%allmodels(fits > 1) = [];

figure(100)
histogram(([allmodels.nfit]))
xlabel('Fit'); ylabel('# of models')

figure(101);

disp('Interpolating for the best models');
for k = 1:length(allmodels)
   
    allmodels_int(k) = interpolate_model(allmodels(k),inv);
    
    n_vs(k)   = allmodels(k).vs.n;
    n_vpvs(k) = allmodels(k).vpvs.n;
    
end

n = inv.min_knots:inv.max_knots;
prior_pdf = 1./(n*(inv.max_knots - inv.min_knots));
prior_pdf = prior_pdf/sum(prior_pdf);
subplot(211)
hold on
histogram(n_vs, 'Normalization', 'Probability')
plot(n, prior_pdf, 'k', 'Linewidth', 3)
legend('Final Probability', 'Starting Probability')
xlabel('# of nodes for Vs model');
ylabel('Probability')

subplot(212)
hold on
histogram(n_vpvs, 'Normalization', 'Probability')
plot(n, prior_pdf, 'k', 'Linewidth', 3)
legend('Final Probability', 'Starting Probability')
xlabel('# of nodes for VpVs model');
ylabel('Probability')

figure(1000)
disp('Plotting Vs...')
z = allmodels_int(1).interp.z;

for k = 1:length(z)
   
    vsset    = [];
    curveset = [];
    
    for i = 1:length(allmodels_int)
        
        vsset    = [ vsset allmodels_int(i).interp.vs(k) ];
        
    end
    
    if log_plot_models
    
        pdfvs(k,:)    = log10( ksdensity( vsset, vsaxis) );
    
    else
       
        pdfvs(k,:)    = ksdensity( vsset, vsaxis);
    
    end
        
    [amp,ind] = max(pdfvs(k,:));
    %pdfvs(k,:) = pdfvs(k, :)/amp;
    vsmode(k) = vsaxis(ind);
    vsmean(k) = mean(vsset);
    
end

[VS,Z] = meshgrid(vsaxis,z);
subplot(121)
contourf(VS, Z, pdfvs, vs_plot)
xlabel('Vs, m/s');
ylabel('Depth, m');
h = colorbar;

if log_plot_models

    h.Label.String = 'Log_1_0(probability density)';

else
   
    h.Label.String = 'Probability Density';
    
end

set(gca, 'YDir', 'reverse')
colormap(flipud(pink))
drawnow

disp('Plotting VpVs...')

for k = 1:length(z)
   
    vpvsset = [];

    for i = 1:length(allmodels_int)
        
        vpvsset = [ vpvsset allmodels_int(i).interp.vpvs(k) ];
        
    end
    
    if log_plot_models
    
        pdfvpvs(k,:) = log10( ksdensity( vpvsset, vpvsaxis) );
        
    else
       
        pdfvpvs(k,:) = ksdensity( vpvsset, vpvsaxis);
        
    end
    
    [amp,ind] = max(pdfvpvs(k,:));
    %pdfvpvs(k,:) = pdfvpvs(k,:)/amp;
    vpvsmode(k) = vpvsaxis(ind);
    vpvsmean(k) = mean(vpvsset);
    
end

[VPVS,Z] = meshgrid(vpvsaxis,z);
subplot(122)
contourf(VPVS, Z, pdfvpvs, vpvs_plot)
xlabel('VpVs ratio');
ylabel('Depth, m');
h = colorbar;

if log_plot_models

    h.Label.String = 'Log_1_0(probability density)';

else
   
    h.Label.String = 'Probability Density';
    
end
    
set(gca, 'YDir', 'reverse')
colormap(flipud(pink))

[ inverse_parameters ] = define_search( );
[~, ~, model] = syn_data_linear(inverse_parameters);

subplot(121)
hold on
plot(vsmean, model.interp.z, 'r--', 'LineWidth', 2)
plot(vsmode, model.interp.z, 'k--', 'LineWidth', 2)
subplot(122)
hold on
plot(vpvsmean, model.interp.z, 'r--', 'LineWidth', 2)
plot(vpvsmode, model.interp.z, 'k--', 'LineWidth', 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now plot the predicted data

figure(1002)
HVmax  = -Inf(size(HV.value));
HVmin  = Inf(size(HV.value));
ZJ0max = -Inf(size(ZJ0.value));
ZJ0min = Inf(size(ZJ0.value));

for i = 1:length(allmodels_int)
    
    for j = 1:length(HV.frequency)
       
        HVmax(j) = max([ HVmax(j) allmodels_int(i).HV(j) ]);
        HVmin(j) = min([ HVmin(j) allmodels_int(i).HV(j) ]);
        ZJ0max(j) = max([ ZJ0max(j) allmodels_int(i).ZJ0(j) ]);
        ZJ0min(j) = min([ ZJ0min(j) allmodels_int(i).ZJ0(j) ]);
        
    end
    
end

subplot(211)
hold on
plot(HV.frequency, HVmax, 'k--');
plot(HV.frequency, HVmin, 'k--');
errorbar(HV.frequency, HV.value, HV.error, 'ko');
xlabel('Frequency, Hz');
ylabel('ZR(\omega)');
xlim([0 20])
subplot(212)
hold on
plot(ZJ0.frequency, ZJ0max, 'k--');
plot(ZJ0.frequency, ZJ0min, 'k--');
errorbar(ZJ0.frequency, ZJ0.value, ZJ0.error, 'ko');
xlabel('Frequency, Hz');
ylabel('CZ(\omega)');
xlim([0 20])
ylim([-0.75 1])


