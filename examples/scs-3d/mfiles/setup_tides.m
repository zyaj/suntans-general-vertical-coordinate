addpath '/Users/fringer/suntans/mfiles';
addpath '/Users/fringer/suntans/mfiles/suntides';
addpath '/Users/fringer/suntans/examples/tides/packages/m_map'

% Corners of simulation domain
MinLati = 18;
MaxLati = 23;
MinLong = 115;
MaxLong = 124;

% Boundaries of projection are larger than domain by amount Wb;
Wb = 1.1;

% Load the m_map files to convert from x,y to lon,lat using UTM
% SCS
m_proj('Transverse Mercator',...
       'lon',[MinLong-Wb MaxLong+Wb],...
       'lat',[MinLati-Wb MaxLati+Wb]);

% Location of xy points of suntans tidal boundaries
tideprefix = './data/tidexy.dat';

% Output location of tidal component files for input into
% suntans.
outputprefix = './data/tidecomponents.dat';

% Location of OTIS data; This example is for the Southern Pacific Ocean
tidespath = './pac';

% Year of simulation (so that phases are relative to start of year)
year = 2007;

% Number of processors
numprocs = 1;

% Data is in lon/lat (not xy projection)
% XY = [];
% Otherwise specify conversion function to convert from xy to lon/lat
XY = 'xy2ll';

% Verbose output
VERBOSE = 1;

% Initialize tidal files
initialize_tides(tidespath)