function [ ET0 ] = PM_PET( elev, Tmax, Tmin, RHmean, lat, Rs, uz, timestep, calcdate, Tmaxold, Tminold)
% This function calculates Penman-Monteith (FAO56) ET0.
% Author: Mustafa Kemal Turkeri
% Written on 8/4/2020
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
% elev: elevation above sea level (m)
% Tmax: Max daily/montly temperature (°C)
% Tmin: Min daily/montly temperature (°C)
% RHmean: Mean Relative Humidity (%)
% lat: latitude of station (°)
% Rs: Solar Radiation (MJ m-2 day-1)
% uz: Wind speed at elevation 2 meters (m s-1)
% timestep: timestep, 'daily' or 'monthly'
% calcdate: date of calculation (first day of calculated month if montly
%           calculation expected ie 5/1/2015 (1 May 2015 for month May) 
% Tmaxold (only required if monthly calculation expected): Max montly
%           temperature of previous month (°C)
% Tminold (only required if monthly calculation expected): Min montly
%           temperature of previous month (°C)
%
% OUTPUTS:
% ET0: Reference evapotranspiration value based on FAO56 method.
%
% Example:
% ET0 = PM_PET( 1050, 25, 13, 75, 36, 25, 3, 'daily', datetime('25-Dec-1989'));
%
% Reference Document:
% Allen, Richard & Pereira, L. & Raes, D. & Smith, M.. (1998). FAO 
% Irrigation and drainage paper No. 56. Rome: Food and Agriculture 
% Organization of the United Nations. 56. 26-40. 
%
% Based on http://www.ipcinfo.org/fileadmin/user_upload/faowater/docs/ReferenceManualV32.pdf

%% Atmospheric Parameters
P = 101.3*((293-0.0065*elev)/293)^(5.26); 																		%atmospheric pressure, kPa

%% Air Temperature
Tmean = (Tmax+Tmin)/2; 																							%Tmean, °C

lambda = 2.501-0.002361*Tmean;																	
psychr = (1.013*10^(-3)*P)/(0.622*lambda);																		%Psychrometric Constant


%% Air Humidity
e0Tmax = 0.6108*exp((17.27*Tmax)/(Tmax+237.3)); 																%max air humidity (calculated) 
e0Tmin = 0.6108*exp((17.27*Tmin)/(Tmin+237.3));																	%min air humidity (calculated) 
e0Tmean = 0.6108*exp((17.27*Tmean)/(Tmean+237.3));																%mean air humidity (calculated) 
es = (e0Tmax+e0Tmin)/2; 																						%saturation vapour pressure [kPa], 
delt = (4098*(0.6108*exp((17.27*Tmean)/(Tmean+237.3)))/(Tmean+237.3)^2); 										%slope vapour pressure curve [kPa °C-1] at Tmean 
ea = e0Tmean*RHmean/100; 																						%Smith, 1992 eq.

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
N = 24/pi()*sunsethourangle; 																					%daylight hours
Rso = (0.75+2*10^(-5)*elev)*Ra; 																				%clear-sky solar radiation [MJ m-2 day-1], 
albedo = 0.23; 																									%standard for grass FAO 56
Rns = (1-albedo)*Rs; 																							%net solar or shortwave radiation [MJ m-2 day-1], 
Rnl = 4.903*10^(-9)*(((Tmax+273.16)^4+(Tmin+273.16)^4)/2)*(0.34-0.14*sqrt(ea))*(1.35*Rs/Rso-0.35); 				%net outgoing longwave radiation [MJ m-2 day-1], 
Rn = Rns - Rnl; 																								%The net radiation (Rn)

%% Wind speed
% Adjust it if necessary to reduce 2 m standard height:
%uelev=elevation of gauge;
%u2 = uz*(4.87/(ln(67.8*uelev-5.42));

u2=uz; 																											%wind speed at 2 m above ground surface [m s-1], i.e no adjustment


%% ET0 - FAO 56 Penman-Monteith Calculation:
if strcmp(timestep,'daily')
    G = 0; 																										%soil heat flux density [MJ m-2 day-1], if it is equal to 24 hours it is zero.
elseif strcmp(timestep,'monthly')
    Tmeanold=(Tmaxold+Tminold)/2;
    G=0.14*(Tmean-Tmeanold); 																					%soil heat flux density [MJ m-2 day-1]
end
ET0 = (0.408*delt*(Rn-G)+psychr*900/(Tmean+273)*u2*(es-ea))/(delt+psychr*(1+0.34*u2));

end

