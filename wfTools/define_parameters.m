function [Parameters] = define_parameters( )

    Parameters.run_name        = 'LandSlideTest_shortday';
    Parameters.data_dir        = './Data/';

    %directory for sac files
    Parameters.sac_directory   = 'C:\Users\Joseph Byrnes\Dropbox\Landslide_data_forJoe\Collected_short\';
    Parameters.experiment_name = 'IrishPark';
    Parameters.channels        = { 'EHZ' 'EHN' 'EHE' };%make it Z first then horizontals 
    Parameters.correlations    = { 'ZZ', 'RR', 'TT', 'ZR', 'ZT', 'RT'};
    
    %%%%%%%%%
    %for spac
    Parameters.normalize        = 0;
    Parameters.complex_flag     = 0;
    Parameters.integrate_flag   = 0;
    Parameters.onebit           = 0;
    Parameters.freq_range       = 0.5:10:40.5;
    Parameters.df               = 0.25*ones(size(Parameters.freq_range));
    Parameters.segment_length   = 100;

    %%%%%%%%%%%%%%
    %for ZR ratios
    %parameters to control the search for rayleigh waves.
    Parameters.baz_step           = 1; %when testing for a rayleigh wave, the code looks at test back azimuths, searching with this increment.
    Parameters.central_f          = 0.25:5:10.25; %central frequencies to look at. In Hz.
    Parameters.halfwidth          = 0.05*ones(size(Parameters.central_f));%half width of the Gaussian filter in Hz
    Parameters.phase_range        = 3; %90 degrees +/- this phase shift to for a Rayleigh wave
    Parameters.TR_max             = 5; %maximum ratio of transverse to radial energy for the test back azimuth. If the transverse is huge the rayleigh wave might not be very good. The results are not very sensitive to this threshold.
    Parameters.hitlength_cycles   = 1.5;
    Parameters.max_hits           = Inf;
    
end

