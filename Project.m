clear, clc;
%----------------------preprocess------------------------------------------
[x, fs] = audioread('.\Samples\3_8_4_5.wav');
time = (1:length(x))*1/fs;                                                  %time in seconds 1_5_7_6_2_9

Hd = humanEar(fs); filteredSignal = filter(Hd, x);                          %adjust bandpass filter to signal
figure("Name", 'Original signal vs Filtered signal'); plot(time, x);
hold on; plot(time, filteredSignal), xlabel('Time'), ylabel('Amplitude');
hold on;

Segments = cut_signal(filteredSignal(:,1), 2048, 2048);                     %split signal to specific windows of 2048 length with no overlap

%---------------------short time energy------------------------------------
E = short_time_energy(filteredSignal(:,1), 2048, 2048);                     %calculate energy for each window
DE = 10*log10(E);                                                           %convert energy to decibel
E_peaks = points_of_interest(DE);                                           %find out speech points in energy
DWE = detect_window_digits(E_peaks);                                        %get windows' indexes of speech points according to energy
figure("Name", 'Short time energy in decibel perceptible by humans');       %energy in decibel
plot(95 + 10*log10(E)); hold on;
figure("Name", 'Short time energy with peaks detected in regions of speech');
plot(DE, '-r'); hold on; plot(E_peaks, 'linestyle', 'none', 'marker','*'); 
xlabel('Window identifier'), ylabel('Energy in decibel'); hold on;

%--------------------zero crossing rate------------------------------------
ZCR = zcr(filteredSignal(:,1), 2048, 2048);                                 %calculate zero crossing rate for each window
DZCR = 10*log10(ZCR);                                                       %convert zero crossing rate to decibel
ZCR_peaks = points_of_interest(DZCR);                                       %find out speech points in zero crossing rate
DWZ = detect_window_digits(ZCR_peaks);                                      %get windows' indexes of speech points according to zero crossing rate
figure("Name", 'Zero crossing rate with peaks detected in regions of speech');
plot(DZCR, '-b'); hold on; plot(ZCR_peaks, 'linestyle', 'none', 'marker','*'); 
xlabel('Window identifier'), ylabel('Zero crossing rate'); hold on;

%-------------------segmentation in digits---------------------------------
Digits = get_digits(DWE, DWZ, Segments);                                    %get separated digits of input signal
%sound(Digits{1,1},fs);                                                     uncomment if you want to sound the first extracted digit
%audiowrite('./one.wav',Digits{1,1},fs);                                    uncomment if you want to write first extracted digit to wav file

%------------------template matching using svm-----------------------------
fprintf("Creating model...\n");
fprintf("Waiting for results...\n");

[iso, labels, max_len_iso] = template_digits();                             %load template words, labels according to filenames and calculate
ISO = [];                                                                   %it will be used for training data with template words
for i=1:length(iso)                                                         %for every template word           
    iel = length(iso{i,1});                                                 %get template word's length
    skip = false;
    temp = [];
    if(iel < max_len_iso)
        temp = zeros(1, max_len_iso - iel);
    else
        skip = true;
    end
    for j=1:iel
        ISO(i,j) = iso{i,1}(1,j);
    end
    if (skip)
        continue;
    else
        ISO(i,j+1:max_len_iso) = temp;
    end
end
results = [];
for m=1:length(Digits)
    DIGIT = [];
    temp = [];
    el = length(Digits{1,m});
    DIGIT(1,1:el) = reshape(Digits{1,m}(:,1),1,el);
    if(el < max_len_iso)
        temp = zeros(1, max_len_iso - el);
        DIGIT(1,el+1:max_len_iso) = temp;
    elseif(el > max_len_iso)
        d = el - max_len_iso;
        [irow, ~ ] = size(ISO);
        temp = zeros(irow, d);
        ISO = horzcat(ISO,temp);
    end
    model = fitcecoc(ISO, labels);
    [result, score] = predict(model,DIGIT);
    results(m) = str2double(result);
end

fprintf("\nRecognized digits:");
disp(results);
  
function Hd = humanEar(fs)                                                  %bandpass filter 20-20000Hz
    Fs = fs;  % Sampling Frequency

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
end

function Segments = cut_signal(X, N, L)                                     %signal, window, overlap (arguments)
    Segments = [];
    m = 0;
    while m * L + N-1 + 1 <= length(X)                                      %while current window has not reached end of signal (is not the last one)
        Segments = [ Segments X(m*L+1:m*L+N-1+1)];
        m = m + 1;
    end
end

function E = short_time_energy(X, N, L)                                     %signal, segment, overlap (arguments)
    m=0;
    E=[];
    while m * L + N-1 + 1 <= length(X)                                      %while current window has not reached end of signal (is not the last one)
        E = [ E sum( X(m*L+1:m*L+N-1+1).^2)/N];
        m = m + 1;
    end
end

function ZCR = zcr(X, N, L)                                                 %signal, segment, overlap (arguments)
    ZCR = [];
    m = 0;
    while m * L + N-1 + 1 <= length(X)
        odd = X(m*L+1:2:m*L+N-1+1);
        even = X(m*L+2:2:m*L+N-1+1);
        ZCR = [ ZCR sum( abs(sgn(even) - sgn(odd)))/(2*(N-1))];
        m = m + 1;
    end
end

function s = sgn(x)                                                         %return signature (-1 or 1) for every input parameter
    s = 1*(x>=0) + (-1)*(x<0);
end

function P = points_of_interest(Y)                                          %input in decibel, detect speech points in windows
    P = [];
    for i=1:length(Y)-1
        if( floor(abs(minus(Y(i), Y(i+1)))) > 0)
            P(i) = Y(i);
        else
            P(i) = NaN;
        end
    end
end

function D = detect_window_digits(peaks)
    D = {};                                                                 %return cell array, #rows = #digits, each row contains window identifiers per digit
    indexes = [];                                                           %column positions of ZCR_peaks when value exists
    j = 1;
    for i=1:length(peaks)
        if(not(isnan(peaks(i))))
            indexes(j) = i;
            j = j + 1;
        end
    end
    %--------------------------------------------------------------
    diff = [];                                                              %difference between ZCR peaks to detect voiced areas
    for i=1:length(indexes)-1
        diff(i)=indexes(i+1)-indexes(i);
    end
    %--------------------------------------------------------------
    j = 1;
    D{1,1}{1} = indexes(1);
    k = 2;
    for i=2:length(indexes)
        if(diff(i-1)>round(max(diff)/2))
            j = j + 1;
            k = 1;
        end 
        D{j,1}{k} = indexes(i);
        k = k + 1;
    end
end

function Digits = get_digits(DWE, DWZ, Segments)                            %return extracted digits
    digits = {};
    [rows, columns] = size(DWE);
    for i=1:rows
        for j=1:columns
            for k=1:length(DWE{i,j})
                if (k==1)
                    digits{i,j}{k} = DWZ{i,j}{k};
                    first = DWZ{i,j}{k};
                    p = 1;
                end
                if (DWE{i,j}{k} > first)
                    p = p + 1;
                    digits{i,j}{p} = DWE{i,j}{k};
                end
            end
        end
        first = 0;
    end
    [r, c] = size(digits);
    Digits={r};
    for i=1:r
        for j=1:c
            last = length(digits{i,j});
            digit = Segments(:,digits{i,j}{1}:digits{i,j}{last});
            [rows_d, columns_d] = size(digit);
            Digits{i} = reshape(digit,rows_d*columns_d,1);
        end
    end
end

function [iso, labels, max_len_iso] = template_digits()                     %load template words and necessary info
    iso = {9};
    labels = {9};
    max_len_iso = 0;                                                        %max length of signal values of all template words
    for i=1:9
        filename = sprintf('%i.wav',i);
        file = fullfile('.\IsolatedDigits',filename);
        x = audioread(file);
        [rows, columns] = size(x(:,1));
        x = reshape(x(:,1), 1, rows*columns);
        [~,name,~] = fileparts(file);
        iso{i,:} = x;
        labels{i,:} = name;
        if (max_len_iso < rows )
            max_len_iso = rows;
        end
    end
end