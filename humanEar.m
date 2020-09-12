function Hd = humanEar
%HUMANEAR Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.7 and Signal Processing Toolbox 8.3.
% Generated on: 08-Sep-2020 23:01:54

% Equiripple Bandpass filter designed using the FIRPM function.

% All frequency values are in Hz.
Fs = 44100;  % Sampling Frequency

Fstop1 = 20;              % First Stopband Frequency
Fpass1 = 2420;            % First Passband Frequency
Fpass2 = 17600;           % Second Passband Frequency
Fstop2 = 20000;           % Second Stopband Frequency
Dstop1 = 0.001;           % First Stopband Attenuation
Dpass  = 0.057501127785;  % Passband Ripple
Dstop2 = 0.0001;          % Second Stopband Attenuation
dens   = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
                          0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

% [EOF]