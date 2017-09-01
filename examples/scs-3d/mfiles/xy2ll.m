%
% XY2LL Convert x,y coordinates to lon,lat.
%   Input must be in meters relative to origin of SUNTANS simulation.
%   which is at (lon,lat)=(115,18).
%
function [lon,lat]=xy2ll(x_relative,y_relative)

  % The path of m_map
  addpath('c:\Users\Oliver Fringer\Documents\matlab\m_map1.4\m_map');

  % Corners of simulation domain
  MinLati = 18;
  MaxLati = 23;
  MinLong = 115;
  MaxLong = 124;

  % Boundaries of projection are larger than domain by amount Wb;
  Wb = 1.1;
  
  % Set up the projection method. 
  m_proj('Transverse Mercator',...
         'lon',[MinLong-Wb MaxLong+Wb],...
         'lat',[MinLati-Wb MaxLati+Wb]);

  % Radius of earth in m
  R = 6378e3; 
  
  % Origin of domain is southwest corner
  [x0,y0] = m_ll2xy(MinLong,MinLati);
  
  % Absolute coordinates for conversion with m_map
  x_absolute = x0 + x_relative/R;
  y_absolute = y0 + y_relative/R;

  [lon,lat] = m_xy2ll(x_absolute,y_absolute);
  lon = lon - 360;