function [Parameters] = define_parameters( )

    Parameters.run_name        = 'Alluvium';

    %directory for sac files
    Parameters.file_type       = 'miniseed';%SAC, miniseed, or IRIS (for iris fetch)
    Parameters.data_directory  = '.\Alluvium\'; %I:\Data\IrishPark\SAC_Data\1000 '.\IrishParkData\' I:\Data\IrishPark\SAC_Data\1000_resp  C:\Users\Joseph Byrnes\Dropbox\Landslide_data_forJoe\Collected_short  I:\Data\IrishPark\SAC_Data\1000   C:\Users\Joseph Byrnes\Dropbox\Landslide_data_forJoe\DiamondArray L:\UMN\2001  
    Parameters.channels        = { 'HHZ' 'HHN' 'HHE' };%make it Z first then north equiv, then east equiv
    Parameters.sign            = [ 1   1   1 ];%nodes have positive down for z
    Parameters.correlations    = { 'ZZ', 'RT' };
    
    Parameters.time_window       = [  737998.264727902  737998.381874859 ];%in datenum. Use [ -Inf Inf ] for all data
    Parameters.download_sections = 15/60/24; %inc for download
    Parameters.downsample        = 5;%during ZR, factor of the central f
    Parameters.bootstrap_samples = 1e2;%make it large if the dataset is small and vice versa
        
    %%%%%%%%
    %parameters only used with IRIS fetch
    
%     stations = irisFetch.Stations('STATION', '*', '*', '*', '*', ...
%         'boxcoordinates', [ 33 38 -121 -114], 'StartBefore', datestr(Parameters.time_window(1), 31), 'EndAfter', datestr(Parameters.time_window(2), 31));
    
%     Parameters.station  = { 'CAG5' 'CAG4' 'CAG3' 'CAF5' 'CAF4' 'CAF3' 'CAH5' 'CAH4' 'CAH3'};%{ stations.StationCode }; %
%     Parameters.network  = { 'YB'   'YB'   'YB'   'YB'   'YB'   'YB'   'YB'   'YB'   'YB'};%can only enter one right now
%     Parameters.location = { '*' };
    
    %%%%%%%%%
    %for spac
    Parameters.normalize        = 1;
    Parameters.make_complex     = 0;
    Parameters.integrate        = 0;
    Parameters.onebit           = 0;
    Parameters.spac_freq        = logspace(0, 1.5, 25);
    Parameters.df               = Parameters.spac_freq/100; %0.01*ones(size(Parameters.freq_range));
    Parameters.segment_length   = 50;%in cycles
    Parameters.filter_type      = 'Gaussian';%Gaussian or Butterworth. Gaussian much faster
    Parameters.time_cull        = 'EndPoints'; %Normally 'EndPoints'. Need 'Intersection' for irish park because of a clock issue
    
    %%%%%%%%%%%%%%
    %for ZR ratios
    %parameters to control the search for rayleigh waves.
    Parameters.baz_step           = 1; %when testing for a rayleigh wave, the code looks at test back azimuths, searching with this increment.
    Parameters.central_f          = logspace(log10(15), log10(20), 20); %central frequencies to look at. In Hz.
    Parameters.halfwidth          = 0.05*Parameters.central_f;%0.01*ones(size(Parameters.central_f));%half width of the Gaussian filter in Hz
    Parameters.phase_range        = 10; %90 degrees +/- this phase shift to for a Rayleigh wave
    Parameters.TR_max             = 5; %maximum ratio of transverse to radial energy for the test back azimuth. If the transverse is huge the rayleigh wave might not be very good. The results are not very sensitive to this threshold.
    Parameters.hitlength_cycles   = 3;
    Parameters.max_hits           = Inf;
    Parameters.sections           = 1;
    
    %%%%%%%%%%%%
    %For SC
    Parameters.window_length = 60;
    Parameters.onebit        = 1;
    Parameters.nsmooth       = 1000;
    Parameters.low           = 1;
    Parameters.high          = 500;
    Parameters.time_keep     = 0.5;
    
    %%%%%%%%%%%%%
    %for wavedec
    Parameters.section_length    = 30;%in cycles
    Parameters.IC                = 'BIC';
    Parameters.maxwaves          = 30;
    Parameters.max_velocity      = 0.2;
    Parameters.iter              = 3;
    Parameters.writeOn           = 5;
    Parameters.taper             = 0;
    Parameters.stalllimit        = 20;
    Parameters.batch             = 10;
    
    Parameters.limits.kappa      = [ 1/1500 1/50 ];%m/s
    Parameters.limits.phi        = [-pi pi];
    Parameters.limits.psi        = [-pi pi];
    Parameters.limits.xi         = [-pi pi];
    
end

