%This script analyzes the first derivative of phase velocity and ellipicity
%ratios, and their robustness against different models. 
clear
close all

addpath(genpath('./scripts/'));
addpath('mat_disperse');

dvs   = 10;
dvpvs = 0.01;
drho  = 0.1;

f = 1:20;%[1 5 10 15];
model.vs.z       = [ 0   10  25  35  49.9 50   100];
model.vpvs.z     = [ 0   10  25  35  49.9 50   100];
%model.vs.value   = [ 600 600 600 300 300  800 800];
%model.vpvs.value = [ 1.7 1.7 1.7 3.5 3.5  1.7  1.7];
model.vs.value   = [ 1500 1500 1500 1500 1500 1500 1500];
model.vpvs.value = [ 1.7  1.7  1.7  1.7  1.7  1.7  1.7];
%model.knots.rho  = nafedrake_rho(model.knots.vs.*model.knots.vpvs/1000);

model.vs.n   = length(model.vs.z);
model.vpvs.n = length(model.vpvs.z);

inverse_parameters = define_search;
model              = interpolate_model(model, inverse_parameters);

% model.interp.vs   = normrnd(model.interp.vs, 10);
% model.interp.vpvs = normrnd(model.interp.vpvs, 0.1);
% model.knots.rho  = nafedrake_rho(model.vs.values*model.vpvs.value/1000);

[PVr.value,~,HV.value,~,~,~,ur,uy] = mat_disperse(diff(model.interp.z),model.interp.rho...
    ,model.interp.vpvs.*model.interp.vs,model.interp.vs,f,1e3);
PVr.value     = PVr.value';
HV.value      = (abs(ur)./abs(uy))';
HV.error      = 0.05;
HV.frequency  = f;
PVr.error     = 100;
PVr.frequency = f;

figure(1)
subplot(231)
plot(model.interp.vs, model.interp.z)
set(gca, 'YDir', 'reverse')
xlabel('Vs, m/s')
ylabel('Depth, m')
subplot(234)
plot(model.interp.vpvs, model.interp.z)
set(gca, 'YDir', 'reverse')
xlabel('VpVs ratio');
ylabel('Depth, m')

%subplot(313)
%plot(model.interp.rho, model.interp.z)
%set(gca, 'YDir', 'reverse')


%%%%%%%%%%%%%%%%%%%
%%%%Vs
for k = 1:length(model.interp.z)
   
    waitbar(k/length(model.interp.z))
    
    %perturb vs
    modelP = model;
    modelN = model;
    modelP.interp.vs(k) = model.interp.vs(k) + dvs;
    modelN.interp.vs(k) = model.interp.vs(k) - dvs;
    
    PVrP = PVr;
    PVrN = PVr;
    HVP  = HV;
    HVN  = HV;
    
    [PVrP.value,~,HVP.value,~,~,~,urP,uyP] = mat_disperse(diff(modelP.interp.z),modelP.interp.rho...
        ,modelP.interp.vpvs.*modelP.interp.vs,modelP.interp.vs,f,1e3);
    PVrP.value     = PVrP.value';
    %HVP.value      = (abs(urP)./abs(uyP))';
    [PVrN.value,~,HVN.value,~,~,~,urN,uyN] = mat_disperse(diff(modelN.interp.z),modelN.interp.rho...
        ,modelN.interp.vpvs.*modelN.interp.vs,modelN.interp.vs,f,1e3);
    PVrN.value     = PVrN.value';
    %HVN.value      = (abs(urN)./abs(uyN))';
    
    for i = 1:length(f)
    
        Vs1stPV(k, i) = (PVrP.value(i) - PVrN.value(i))/(2*dvs);
        %Vs2ndPV(k, i) = (PVrP.value(i) - 2*PVr.value(i) + PVrN.value(i))/(dvs*dvs);
        
        Vs1stHV(k, i) = (HVP.value(i) - HVN.value(i))/(2*dvs);
        %Vs2ndHV(k, i) = (HVP.value(i) - 2*HV.value(i) + HVN.value(i))/(dvs*dvs);
    
    end
    
end

subplot(232)
plot(Vs1stPV(1:end-1, :), model.interp.z(1:end-1))
set(gca, 'YDir', 'reverse')
xlabel('\delta(PV)/\delta(Vs)')
ylabel('Depth, m')

%subplot(212)
%plot(Vs2ndPV(1:end-1, :), model.interp.z(1:end-1))
%set(gca, 'YDir', 'reverse')
%xlabel('\delta^2(PV)/\delta(Vs)^2')
legend('1 Hz', '5 Hz', '10 Hz', '15 Hz', 'Location', 'SouthEast')
%figure(3)
subplot(233)
plot(Vs1stHV(1:end-1, :), model.interp.z(1:end-1))
set(gca, 'YDir', 'reverse')
xlabel('\delta(ZR)/\delta(Vs)')
ylabel('Depth, m')

%legend('1 Hz', '5 Hz', '10 Hz', '15 Hz')
%subplot(212)
%plot(Vs2ndHV(1:end-1, :), model.interp.z(1:end-1))
%set(gca, 'YDir', 'reverse')
%xlabel('\delta^2(HV)/\delta(Vs)^2')

%%%%%%%%%%%%%%%%%%%
%%%%VpVs
for k = 1:length(model.interp.z)
   
    waitbar(k/length(model.interp.z))
    
    %perturb vs
    modelP = model;
    modelN = model;
    modelP.interp.vpvs(k) = model.interp.vpvs(k) + dvpvs;
    modelN.interp.vpvs(k) = model.interp.vpvs(k) - dvpvs;
    
    PVrP = PVr;
    PVrN = PVr;
    HVP  = HV;
    HVN  = HV;
    
    [PVrP.value,~,HVP.value,~,~,~,urP,uyP] = mat_disperse(diff(modelP.interp.z),modelP.interp.rho...
        ,modelP.interp.vpvs.*modelP.interp.vs,modelP.interp.vs,f,1e3);
    PVrP.value     = PVrP.value';
    %HVP.value      = (abs(urP)./abs(uyP))';
    [PVrN.value,~,HVN.value,~,~,~,urN,uyN] = mat_disperse(diff(modelN.interp.z),modelN.interp.rho...
        ,modelN.interp.vpvs.*modelN.interp.vs,modelN.interp.vs,f,1e3);
    PVrN.value     = PVrN.value';
    %HVN.value      = (abs(urN)./abs(uyN))';
    
    for i = 1:length(f)
    
        VpVs1stPV(k, i) = (PVrP.value(i) - PVrN.value(i))/(2*dvpvs);
        %VpVs2ndPV(k, i) = (PVrP.value(i) - 2*PVr.value(i) + PVrN.value(i))/(dvpvs*dvpvs);
        
        VpVs1stHV(k, i) = (HVP.value(i) - HVN.value(i))/(2*dvpvs);
        %VpVs2ndHV(k, i) = (HVP.value(i) - 2*HV.value(i) + HVN.value(i))/(dvpvs*dvpvs);
    
    end
    
end

subplot(235)
plot(VpVs1stPV(1:end-1, :), model.interp.z(1:end-1))
set(gca, 'YDir', 'reverse')
xlabel('\delta(PV)/\delta(VpVs)')
ylabel('Depth, m')

% subplot(212)
% plot(VpVs2ndPV(1:end-1, :), model.interp.z(1:end-1))
% set(gca, 'YDir', 'reverse')
% xlabel('\delta^2(PV)/\delta(VpVs)^2')

subplot(236)
plot(VpVs1stHV(1:end-1, :), model.interp.z(1:end-1))
set(gca, 'YDir', 'reverse')
xlabel('\delta(ZR)/\delta(VpVs)')
ylabel('Depth, m')

% subplot(212)
% plot(VpVs2ndHV(1:end-1, :), model.interp.z(1:end-1))
% set(gca, 'YDir', 'reverse')
% xlabel('\delta^2(HV)/\delta(VpVs)^2')

% %%%%%%%%%%%%%%%%%%%
% %%%%Rho
% for k = 1:length(model.interp.z)
%    
%     waitbar(k/length(model.interp.z))
%     
%     %perturb vs
%     modelP = model;
%     modelN = model;
%     modelP.interp.rho(k) = model.interp.rho(k) + drho;
%     modelN.interp.rho(k) = model.interp.rho(k) - drho;
%     
%     PVrP = PVr;
%     PVrN = PVr;
%     HVP  = HV;
%     HVN  = HV;
%     
%     [PVrP.value,~,~,~,~,~,urP,uyP] = mat_disperse(diff(modelP.interp.z),modelP.interp.rho...
%         ,modelP.interp.vpvs.*modelP.interp.vs,modelP.interp.vs,f,1e3);
%     PVrP.value     = PVrP.value';
%     HVP.value      = (abs(urP)./abs(uyP))';
%     [PVrN.value,~,~,~,~,~,urN,uyN] = mat_disperse(diff(modelN.interp.z),modelN.interp.rho...
%         ,modelN.interp.vpvs.*modelN.interp.vs,modelN.interp.vs,f,1e3);
%     PVrN.value     = PVrN.value';
%     HVN.value      = (abs(urN)./abs(uyN))';
%     
%     for i = 1:length(f)
%     
%         Rho1stPV(k, i) = (PVrP.value(i) - PVrN.value(i))/(2*drho);
%         Rho2ndPV(k, i) = (PVrP.value(i) - 2*PVr.value(i) + PVrN.value(i))/(drho*drho);
%         
%         Rho1stHV(k, i) = (HVP.value(i) - HVN.value(i))/(2*drho);
%         Rho2ndHV(k, i) = (HVP.value(i) - 2*HV.value(i) + HVN.value(i))/(drho*drho);
%     
%     end
%     
% end
% 
% figure(6)
% subplot(211)
% plot(Rho1stPV(1:end-1, :), model.interp.z(1:end-1))
% set(gca, 'YDir', 'reverse')
% xlabel('\delta(PV)/\delta(Rho)')
% subplot(212)
% plot(Rho2ndPV(1:end-1, :), model.interp.z(1:end-1))
% set(gca, 'YDir', 'reverse')
% xlabel('\delta^2(PV)/\delta(Rho)^2')
% 
% figure(7)
% subplot(211)
% plot(Rho1stHV(1:end-1, :), model.interp.z(1:end-1))
% set(gca, 'YDir', 'reverse')
% xlabel('\delta(HV)/\delta(Rho)')
% subplot(212)
% plot(Rho2ndHV(1:end-1, :), model.interp.z(1:end-1))
% set(gca, 'YDir', 'reverse')
% xlabel('\delta^2(HV)/\delta(Rho)^2')
% 
% %%%%%% now do an experiment
% %if I take the first and second derviative, how much does that help me?
% 
% VpVsVector = 1.5:0.01:3.5;
% 
% find = 2;
% 
% for i = 1:length(VpVsVector)
% 
%     modelN.knots.z    = [ 0   10  25  35            49.9           50   100];
%     modelN.knots.vs   = [ 500 800 800 600           600            1500 1500];
%     modelN.knots.vpvs = [ 1.7 1.7 1.7 VpVsVector(i) VpVsVector(i)  1.7  1.7];
%     modelN.knots.rho  = nafedrake_rho(model.knots.vs.*model.knots.vpvs/1000);
%     modelN            = interpolate_model(modelN, inverse_parameters);
% 
%     PV_kernel1(i) = PVr.value(find);
%     PV_kernel2(i) = PVr.value(find);
%     
%     for k = 1:length(model.interp.z)
% 
%         PV_kernel1(i) = PV_kernel1(i) + VpVs1stPV(k, find)*(modelN.interp.vpvs(k) - model.interp.vpvs(k));
%         PV_kernel2(i) = PV_kernel2(i) + VpVs1stPV(k, find)*(modelN.interp.vpvs(k) - model.interp.vpvs(k)) ...
%             + VpVs2ndPV(k, find)*(modelN.interp.vpvs(k) - model.interp.vpvs(k)).^2/2;
% 
%     end
% 
%     HV_kernel1(i) = HV.value(find);
%     HV_kernel2(i) = HV.value(find);
%     
%     for k = 1:length(model.interp.z)
% 
%         HV_kernel1(i) = HV_kernel1(i) + VpVs1stHV(k, find)*(modelN.interp.vpvs(k) - model.interp.vpvs(k));
%         HV_kernel2(i) = HV_kernel1(i) + VpVs1stHV(k, find)*(modelN.interp.vpvs(k) - model.interp.vpvs(k)) ...
%             + VpVs2ndHV(k, find)*(modelN.interp.vpvs(k) - model.interp.vpvs(k)).^2/2;
% 
%     end
%         
%     %get the real value
%     [PV_real(i),~,~,~,~,~,urN,uyN] = mat_disperse(diff(modelN.interp.z),modelN.interp.rho...
%         ,modelN.interp.vpvs.*modelN.interp.vs,modelN.interp.vs,f(find),1e3);
%     HV_real(i)      = (abs(urN)./abs(uyN))';
%     
% end
% 
% figure(8)
% subplot(221)
% hold on
% plot(VpVsVector, PV_kernel1)
% plot(VpVsVector, PV_kernel2)
% plot(VpVsVector, PV_real, 'k');
% xlabel('VpVs varied in layer');
% ylabel('Phase velocity');
% subplot(222)
% hold on
% plot(VpVsVector, abs(PV_kernel1 - PV_real))
% plot(VpVsVector, abs(PV_kernel2 - PV_real))
% xlabel('VpVs varied in layer');
% ylabel('Error in Phase velocity');
% 
% subplot(223)
% hold on
% plot(VpVsVector, HV_kernel1)
% plot(VpVsVector, HV_kernel2)
% plot(VpVsVector, HV_real, 'k');
% xlabel('VpVs varied in layer');
% ylabel('HV ratio');
% subplot(224)
% hold on
% plot(VpVsVector, abs(HV_kernel1 - HV_real))
% plot(VpVsVector, abs(HV_kernel2 - HV_real))
% xlabel('VpVs varied in layer');
% ylabel('Error in HV ratio');
% 
% %%%%%% now do an experiment
% %if I take the first and second derviative, how much does that help me?
% 
% VsVector = 200:10:1000;
% 
% PV_kernel1 = [];
% HV_kernel1 = [];
% PV_kernel2 = [];
% HV_kernel2 = [];
% HV_real   = [];
% PV_real   = [];
% find = 3;
% 
% for i = 1:length(VsVector)
% 
%     modelN.knots.z    = [ 0           10  25  35  49.9 50   100];
%     modelN.knots.vs   = [ VsVector(i) 800 800 600 600  1500 1500];
%     modelN.knots.vpvs = [ 1.7         1.7 1.7 2.5 2.5  1.7  1.7];
%     modelN.knots.rho  = nafedrake_rho(model.knots.vs.*model.knots.vpvs/1000);
%     modelN            = interpolate_model(modelN, inverse_parameters);
% 
%     PV_kernel1(i) = PVr.value(find);
%     PV_kernel2(i) = PVr.value(find);
%     
%     for k = 1:length(model.interp.z)
% 
%         PV_kernel1(i) = PV_kernel1(i) + Vs1stPV(k, find)*(modelN.interp.vs(k) - model.interp.vs(k));
%         PV_kernel2(i) = PV_kernel2(i) + Vs1stPV(k, find)*(modelN.interp.vs(k) - model.interp.vs(k)) ...
%             + Vs2ndPV(k, find)*(modelN.interp.vs(k) - model.interp.vs(k)).^2/2;
% 
%     end
% 
%     HV_kernel1(i) = HV.value(find);
%     HV_kernel2(i) = HV.value(find);
%     
%     for k = 1:length(model.interp.z)
% 
%         HV_kernel1(i) = HV_kernel1(i) + Vs1stHV(k, find)*(modelN.interp.vs(k) - model.interp.vs(k));
%         HV_kernel2(i) = HV_kernel2(i) + Vs1stHV(k, find)*(modelN.interp.vs(k) - model.interp.vs(k)) ...
%             + Vs2ndHV(k, find)*(modelN.interp.vs(k) - model.interp.vs(k)).^2/2;
% 
%     end
%         
%     %get the real value
%     [PV_real(i),~,~,~,~,~,urN,uyN] = mat_disperse(diff(modelN.interp.z),modelN.interp.rho...
%         ,modelN.interp.vpvs.*modelN.interp.vs,modelN.interp.vs,f(find),1e3);
%     HV_real(i)      = (abs(urN)./abs(uyN))';
%     
% end
% 
% figure(9)
% subplot(221)
% hold on
% plot(VsVector, PV_kernel1)
% plot(VsVector, PV_kernel2)
% plot(VsVector, PV_real, 'k');
% xlabel('Vs varied in layer');
% ylabel('Phase velocity');
% subplot(222)
% hold on
% plot(VsVector, abs(PV_kernel1 - PV_real))
% plot(VsVector, abs(PV_kernel2 - PV_real))
% xlabel('Vs varied in layer');
% ylabel('Error in Phase velocity');
% 
% subplot(223)
% hold on
% plot(VsVector, HV_kernel1)
% plot(VsVector, HV_kernel2)
% plot(VsVector, HV_real, 'k');
% xlabel('Vs varied in layer');
% ylabel('HV ratio');
% subplot(224)
% hold on
% plot(VsVector, abs(HV_kernel1 - HV_real))
% plot(VsVector, abs(HV_kernel2 - HV_real))
% xlabel('Vs varied in layer');
% ylabel('Error in HV ratio');
% 
% figure
% subplot(211)
% plot(Vs1stPV(1:end-1, :), model.interp.z(1:end-1))
% set(gca, 'YDir', 'reverse')
% xlabel('\delta(PV)/\delta(Vs)')
% ylabel('Depth (m)');
% legend('1  Hz', '5 Hz', '10 Hz', '15 Hz', 'Location', 'SouthEast')
% subplot(212)
% plot(Vs1stHV(1:end-1, :), model.interp.z(1:end-1))
% set(gca, 'YDir', 'reverse')
% xlabel('\delta(HV)/\delta(Vs)')
% ylabel('Depth (m)');
% legend('1  Hz', '5 Hz', '10 Hz', '15 Hz', 'Location', 'SouthWest')
