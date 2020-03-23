function fitParam = fit_AzimAniso(theData,func,tf_demean,tf_dc)
% theData = [xdata, ydata, weights]

%%%% Things to add %%%%
% 1) Curve fit statistic (R2?)
% 2) Alternative error calculations for fit parameters (and when data error
%       is not given/ known)
%           a) F-test
%           b) standard error calculation
% 3) Non-linear error estimates

% Checks
if length(theData(1,:)) < 2
    error('theData must contain at least two columns, xdata and ydata')
end

% Demean residuals
if tf_demean
    theData(:,2) = theData(:,2) - mean(theData(:,2));
    tf_dc = 0;
end
% Fit DC-offset
if tf_dc
    A = ones(length(theData(:,1)),1);
else
    A = [];
end
% Range over which fitted function is calculated
xrange = [min(theData(:,1)),max(theData(:,1))];
xrange(1) = xrange(1) - 0.25*abs(xrange(1));
xrange(2) = xrange(2) + 0.25*abs(xrange(2));

% Build system of equations
switch func
    % Elliptical/ 2theta fit
    case '2theta'
        fitParam.header = {'Y';'A';'B';'X1';'alpha'};
        A = [A, cosd(2*theData(:,1)), sind(2*theData(:,1))];
        if length(theData(1,:)) == 3
            [fitParam.value,fitParam.error] = ...
                lscov(A,theData(:,2),theData(:,3));
        else
            fitParam.value = lsqr(A,theData(:,2));
            fitParam.error = zeros(size(fitParam.value));
        end
        
        if length(fitParam.value) ~= 3;
            fitParam.value = [0; fitParam.value];
            fitParam.error = [0; fitParam.error];
        end
        
        % Calculate rotation and amplitude
        A  = fitParam.value(2); B  = fitParam.value(3); 
        alpha = 0.5*atan2(B,A);
        X = A/cos(2*alpha);
        
        % Calculate errors in rotation and amplitude
        dA = fitParam.error(2); dB = fitParam.error(3);
        % Error Propogation Method
        dalpha = 0.5*sqrt((((1/((2*A)+((B^2)/A)))^2)*(dB^2))+(((-B/((4*A^2)+(B^2)))^2)*(dA^2)));
        dXa = (A^2)*sqrt(2*((dA/A)^2)); dXb = (B^2)*sqrt(2*((dB/B)^2)); dX3 = sqrt((dXa^2)+(dXb^2));
        dX     = (((A^2)+(B^2))^-0.5)*dX3;
        
        % Alternative error estimation (does not work if rotation is near
        % pi/2).
%         dalpha = std([0.5*atan2(B,A); 0.5*atan2(B+dB,A+dA); 0.5*atan2(B-dB,A-dA);...
%             0.5*atan2(B+dB,A-dA); 0.5*atan2(B-dB,A+dA); 0.5*atan2(B+dB,A);...
%             0.5*atan2(B-dB,A); 0.5*atan2(B,A+dA); 0.5*atan2(B,A-dA)]);
%         dX     = std([A/cos(2*alpha); (A+dA)/cos(2*(alpha+dalpha)); ...
%             (A-dA)/cos(2*(alpha-dalpha)); (A+dA)/cos(2*(alpha-dalpha)); (A-dA)/cos(2*(alpha+dalpha));...
%             (A)/cos(2*(alpha+dalpha)); (A)/cos(2*(alpha-dalpha)); (A+dA)/cos(2*alpha); (A-dA)/cos(2*alpha)]);
        
        fitParam.value = [fitParam.value; X; rad2deg(alpha)];
        fitParam.error = [fitParam.error; dX; rad2deg(dalpha)];
        % Define best-fit curve
        fitParam.x  = linspace(xrange(1),xrange(2),1000*diff(xrange));
        fitParam.y  = fitParam.value(1) + X*cosd(2*(fitParam.x - rad2deg(alpha)));
        
    % 2theta + 4theta fit    
    case '2,4theta'
        fitParam.header = {'Y';'A';'B';'C';'D';'X1';'X2';'alpha';'beta'};
        A = [A, cosd(2*theData(:,1)), sind(2*theData(:,1)), ...
            cosd(4*theData(:,1)), sind(4*theData(:,1))];
        if length(theData(1,:)) == 3
            [fitParam.value,fitParam.error] = ...
                lscov(A,theData(:,2),theData(:,3));
        else
            fitParam.value = lsqr(A,theData(:,2));
            fitParam.error = zeros(size(fitParam.value));
        end
        
        if length(fitParam.value) ~= 5;
            fitParam.value = [0; fitParam.value];
            fitParam.error = [0; fitParam.error];
        end
        
        % Calculate rotations and amplitudes
        A = fitParam.value(2); B = fitParam.value(3); C = fitParam.value(4); D = fitParam.value(5);
        alpha = 0.5*atan2(B,A);
        beta  = 0.25*atan2(D,C);
        X1    = A/cos(2*alpha);
        X2    = C/cos(4*beta);
        
        % Calculate errors in rotations and amplitudes
        dA = fitParam.error(2); dB = fitParam.error(3); dC = fitParam.error(4); dD = fitParam.error(5);
        
        % Error Propogation
        dalpha    = 0.5*sqrt((((1/((2*A)+((B^2)/A)))^2)*(dB^2))+(((-B/((4*A^2)+(B^2)))^2)*(dA^2)));
        dbeta     = 0.25*sqrt((((1/((2*C)+((D^2)/C)))^2)*(dD^2))+(((-D/((4*C^2)+(D^2)))^2)*(dC^2)));
        % 2-theta amplitude error
        dX1a = (A^2)*sqrt(2*((dA/A)^2));
        dX1b = (B^2)*sqrt(2*((dB/B)^2));
        dX1c = sqrt((dX1a^2)+(dX1b^2));
        dX1  = (((A^2)+(B^2))^-0.5)*dX1c;
        % 4-theta amplitude error
        dX2a = (C^2)*sqrt(2*((dC/C)^2));
        dX2b = (D^2)*sqrt(2*((dD/D)^2));
        dX2c = sqrt((dX2a^2)+(dX2b^2)); 
        dX2  = (((C^2)+(D^2))^-0.5)*dX2c;
        
        fitParam.value = [fitParam.value; X1; X2; rad2deg(alpha); rad2deg(beta)];
        fitParam.error = [fitParam.error; dX1; dX2; rad2deg(dalpha); rad2deg(dbeta)];
        % Define best-it curve
        fitParam.x     = linspace(xrange(1),xrange(2),1000*diff(xrange));
        fitParam.y     = fitParam.value(1) + X1*cosd(2*(fitParam.x - rad2deg(alpha))) ...
                            + X2*cosd(4*(fitParam.x - rad2deg(beta)));
        
    % Thomsen parameters--non-linear fit    
    case 'thomsen'
        x               = lsqcurvefit(@Thomsen,[0,0,0,0],theData(:,1),theData(:,2));
        fitParam.header = {'ao','delta','alpha','epsilon'};
        fitParam.value  = x(:);
        fitParam.error  = nan(size(fitParam.value));
        fitParam.x      = linspace(xrange(1),xrange(2),1000*diff(xrange))';
        fitParam.y      = x(1) + x(2)*((sind(fitParam.x-x(3))).^2).*((cosd(fitParam.x-x(3))).^2)...
                            + x(4)*((sind(fitParam.x-x(3))).^4);
                        
    case '1,2theta'
        fitParam.header = {'Y';'A';'B';'C';'D';'X1';'X2';'alpha';'beta'};
        A = [A, cosd(theData(:,1)), sind(theData(:,1)), ...
            cosd(2*theData(:,1)), sind(2*theData(:,1))];
        if length(theData(1,:)) == 3
            [fitParam.value,fitParam.error] = ...
                lscov(A,theData(:,2),theData(:,3));
        else
            fitParam.value = lsqr(A,theData(:,2));
            fitParam.error = zeros(size(fitParam.value));
        end
        
        if length(fitParam.value) ~= 5;
            fitParam.value = [0; fitParam.value];
            fitParam.error = [0; fitParam.error];
        end
        
        % Calculate rotations and amplitudes
        A = fitParam.value(2); B = fitParam.value(3); C = fitParam.value(4); D = fitParam.value(5);
        alpha = atan2(B,A);
        beta  = 0.5*atan2(D,C);
        X1    = A/cos(alpha);
        X2    = C/cos(2*beta);
        
        % Calculate errors in rotations and amplitudes
        dA = fitParam.error(2); dB = fitParam.error(3); dC = fitParam.error(4); dD = fitParam.error(5);
        
        % Error Propogation
        dalpha    = sqrt((((1/((2*A)+((B^2)/A)))^2)*(dB^2))+(((-B/((4*A^2)+(B^2)))^2)*(dA^2)));
        dbeta     = 0.5*sqrt((((1/((2*C)+((D^2)/C)))^2)*(dD^2))+(((-D/((4*C^2)+(D^2)))^2)*(dC^2)));
        % 1-theta amplitude error
        dX1a = (A^2)*sqrt(2*((dA/A)^2));
        dX1b = (B^2)*sqrt(2*((dB/B)^2));
        dX1c = sqrt((dX1a^2)+(dX1b^2));
        dX1  = (((A^2)+(B^2))^-0.5)*dX1c;
        % 2-theta amplitude error
        dX2a = (C^2)*sqrt(2*((dC/C)^2));
        dX2b = (D^2)*sqrt(2*((dD/D)^2));
        dX2c = sqrt((dX2a^2)+(dX2b^2)); 
        dX2  = (((C^2)+(D^2))^-0.5)*dX2c;
        
        fitParam.value = [fitParam.value; X1; X2; rad2deg(alpha); rad2deg(beta)];
        fitParam.error = [fitParam.error; dX1; dX2; rad2deg(dalpha); rad2deg(dbeta)];
        % Define best-it curve
        fitParam.x     = linspace(xrange(1),xrange(2),1000*diff(xrange));
        fitParam.y     = fitParam.value(1) + X1*cosd((fitParam.x - rad2deg(alpha))) ...
                            + X2*cosd(2*(fitParam.x - rad2deg(beta)));
        
        
end
% fitParam.table = table(fitParam.value,fitParam.error,'RowNames',fitParam.header);
% display(fitParam.table);