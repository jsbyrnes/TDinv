function [Parameters] = define_parameters( )

    Parameters.run_name        = 'DiamondArray';

    %directory for sac files
    Parameters.file_type       = 'SAC';
    Parameters.data_directory  = 'C:\Users\Joseph Byrnes\Dropbox\Landslide_data_forJoe\DiamondArray';%  '.\IrishParkData\' I:\Data\IrishPark\SAC_Data\1000_resp  C:\Users\Joseph Byrnes\Dropbox\Landslide_data_forJoe\Collected_short  I:\Data\IrishPark\SAC_Data\1000   C:\Users\Joseph Byrnes\Dropbox\Landslide_data_forJoe\DiamondArray L:\UMN\2001  
    Parameters.channels        = { 'EHZ' 'EHN' 'EHE' };%make it Z first then north equiv, then east equiv 
    Parameters.sign            = [ -1  1   1   ];%nodes have positive down for z
    Parameters.correlations    = { 'ZZ' 'RR' 'TT' 'ZR' 'ZT', 'RT' }; %, 'RR', 'TT', 'ZR', 'ZT', 'RT'
    
    %Parameters.time_window      = [ 737032.737110185 737032.759287269 ];%in datenum. Use [ -Inf Inf ] for all data
    Parameters.time_window       = [ 736497.546082949 736497.80046083 ];
    Parameters.downsample        = 5;%Will downsample to central_f*Parameters.downsample for ZR
    Parameters.bootstrap_samples = 1e4;%make it large if the dataset is small and vice versa
    
    %%%%%%%%%
    %for spac
    Parameters.normalize        = 0;
    Parameters.make_complex     = 0;
    Parameters.integrate        = 0;
    Parameters.onebit           = 0;
    Parameters.freq_range       = 1:0.5:30;
    Parameters.df               = 0.05*ones(size(Parameters.freq_range));
    Parameters.segment_length   = 50;
    Parameters.filter_type      = 'Gaussian';%Gaussian or Butterworth
    Parameters.time_cull        = 'EndPoints'; %Normally 'EndPoints'. Need 'Intersection' for irish park because of a clock issue
    
    %%%%%%%%%%%%%%
    %for ZR ratios
    %parameters to control the search for rayleigh waves.
    Parameters.baz_step           = 1; %when testing for a rayleigh wave, the code looks at test back azimuths, searching with this increment.
    Parameters.central_f          = 0.5:0.5:30; %central frequencies to look at. In Hz.
    Parameters.halfwidth          = 0.05*ones(size(Parameters.central_f));%half width of the Gaussian filter in Hz
    Parameters.phase_range        = 3; %90 degrees +/- this phase shift to for a Rayleigh wave
    Parameters.TR_max             = 5; %maximum ratio of transverse to radial energy for the test back azimuth. If the transverse is huge the rayleigh wave might not be very good. The results are not very sensitive to this threshold.
    Parameters.hitlength_cycles   = 1.5;
    Parameters.max_hits           = Inf;
    Parameters.sections           = 1;
    
end

