function [ G ] = drive_mtm(C1, C2, C3 , Model, rf_shift)
%DRIVE_MTM Do multi taper decon
%   C1 is the source

    if [~ 1] ~= size(C1)

        C1 = mean(C1,2);
        C2 = mean(C2,2);
        C3 = mean(C3,2);

    end

    n = length(C1);

    S = C1;
    
    C1 = timeshift(C1, n, rf_shift);
    C2 = timeshift(C2, n, rf_shift);
    C3 = timeshift(C3, n, rf_shift);

    [G.C1, G.FFT_C1] = mtmdeconSP(S, C1, n, Model.ntap, Model.nmtw, Model.NW);
    [G.C2, G.FFT_C2] = mtmdeconSP(S, C2, n, Model.ntap, Model.nmtw, Model.NW);
    [G.C3, G.FFT_C3] = mtmdeconSP(S, C3, n, Model.ntap, Model.nmtw, Model.NW);

end