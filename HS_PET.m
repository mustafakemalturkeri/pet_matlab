function [ ET0 ] = HS_PET( Tmax, Tmin, lat, timestep, calcdate, Tmean, C0)
% This function calculates Penman-Monteith (FAO56) ET0.
% Author: Mustafa Kemal Turkeri
% Written on 4/12/2022
%
%MIT License
%
%Copyright (c) 2022 Mustafa Kemal Türkeri
%
%Permission is hereby granted, free of charge, to any person obtaining a copy
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:
%
%The above copyright notice and this permission notice shall be included in all
%copies or substantial portions of the Software.
%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%SOFTWARE.
%
% INPUTS:
% Tmax: Max daily/montly temperature (°C)
% Tmin: Min daily/montly temperature (°C)
% lat: latitude of station (°)
% timestep: timestep, 'daily' or 'monthly'
% calcdate: date of calculation (first day of calculated month if montly
%
% Optional inputs:
% Tmean is the mean monthly temperature: Enter for only monthly calculations
% C0: This parameter is for optimisation only, use default value 0.0023 for
% original method.
%
% OUTPUTS:
% ET0: Reference evapotranspiration value based on Hargreaves-Samani
% method. (mm/day)
%
% Example:
% ET0 = HS_PET(25,13,36,'daily',datetime('01-Jan-2022'));

% Reference Document:
% Hargreaves, G.H. and Samani, Z.A. (1985) Reference Crop Evapotranspiration
% from Temperature. Applied Engineering in Agriculture, 1, 96-99.
% http://dx.doi.org/10.13031/2013.26773

if (nargin < 5)
	error('Not enough inputs!');
elseif	(nargin < 6 && strcmp(timestep,'monthly'))
	error('Not enought input arguments for monthly calculations!');
end

if (nargin < 7)
	C0 = 0.0023;
end

%% Air Temperature
if strcmp(timestep,'daily')
    Tmean = (Tmax+Tmin)/2; 																						%Tmean, °C, daily only!
end

%% Radiation
if strcmp(timestep,'daily')
    J = day(calcdate,'dayofyear');
elseif strcmp(timestep,'monthly')
    mon=month(calcdate);
    mons=eomday(year(calcdate),1:mon);
    J = sum(mons(1:end-1))+mons(end)/2;
end
dr = 1+0.033*cos(2*pi()/365*J); 																				%inverse relative distance Earth-Sun 
sdec = 0.409*sin(2*pi()/365*J-1.39); 																			%solar declination (Equation 3.17) [rad].
latrad=lat/180*pi(); 																							%latitude in rad
sunsethourangle = acos(-tan(latrad)*tan(sdec)); 																%solar declination (Equation 3.17) [rad].
Ra = 24*60/pi()*0.0820*dr*(sunsethourangle*sin(latrad)*sin(sdec)+cos(latrad)*cos(sdec)*sin(sunsethourangle)); 	%extraterrestrial radiation [MJ m-2 day-1], 


%% ET0 - Hargreaves-Samani Calculation:
ET0 = C0*Ra*(Tmax-Tmin)^(0.5)*(Tmean+17.8);

end