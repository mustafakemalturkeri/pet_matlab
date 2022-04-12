function [PET] = Thornthwaite_PET(T,NL)
% This program is based on a code written by Zekai Sen on 16 July 2018
% Modified by Mustafa Kemal Turkeri for use in MARRMoT.
% Tested in both Matlab and Octave
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
% Remark: If the temperature is below 0, then PET calculation gives
% a complex number. Here in this code, only the real part of the 
% complex number is taken and used.
%
% The code only works for monthly calculations!
%
% Only valid for 0-60 latitudes of northern hemisphere!
%
% T     : Monthly air temperature in Celsius 1x12 values
% NL    : North latitude (0 < NL <60)
%
%
% Reference:
% Thornthwaite CW (1948) An approach toward a rational classification of
% climate. Geographical Review, 38, 55-94.
%
%
N=length(T); % Number of months
NY=floor(N/12);     % Number of years
K=[0    1   1   1   1   1   1   1   1   1   1   1   1
    10  0.97    0.98    1   1.03    1.05    1.06    1.05    1.04    1.02    0.99    0.97    0.96
    20  0.92    0.98    1   1.05    1.09    1.11    1.1 1.07    1.02    0.98    0.93    0.91
    30  0.87    0.98    1   1.07    1.14    1.17    1.16    1.11    1.03    0.96    0.89    0.85
    40  0.8 0.89    0.99    1.1 1.2 1.25    1.23    1.15    1.04    0.93    0.83    0.78
    50 0.71 0.84 0.98 1.14 1.28 1.36 1.33 1.21 1.06 0.9 0.76 0.68
    60  0.54    0.67    0.97    1.19    1.33    1.56    1.55    1.33    1.07    0.84    0.58    0.48
    ]; % K coefficient in the Thornthwaite PET calculation depending on north latitude and months
for m=1:NY
    TT=T(12*(m-1)+1:12*m);
    I=(TT/5).^1.514;
    J=sum(real(I));
    c=0.000000675*J^3-0.0000771*J^2+0.01792*J+0.49239;
    MNL=K(:,1)/NL; % Middle value of latitude
    for i=1:7
        if MNL(i,1) == 1
            PETT=16*K(i,2:13).*(10*TT'./J).^c;
        else
        end
    end
    for i=1:6
        if MNL(i) < 1 && MNL(i+1) > 1
            LNL=K(i,1);
            UNL=K(i+1,1);
            R=(NL-LNL)/(UNL-LNL);
            CK=K(i,2:13)+R*(K(i+1,2:13)-K(i,2:13)); % Calculation of K
            PETT=16*CK(1,1:12)'.*(10*TT/J).^c; %modified for octave
        else
        end
    end
    PETM(12*(m-1)+1:12*m)=real(PETT); % Monthly PET
end
k=0;
NM=12*NY;
for i=1:NM
    k=k+1;
    PET(k)=PETM(1,i);
end
for i=1:N
    if PET(i) < 0
        PET(i)=0;
    else
    end
end
end