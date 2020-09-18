clear, clc;

[x, fs] = audioread('.\Samples\6_numbers.wav');
time = (1:length(x))*1/fs;  %time in seconds

Hd = humanEar(); filteredSignal = filter(Hd, x);
figure("Name", 'Original signal vs Filtered signal'); plot(time, x);
hold on; plot(time, filteredSignal), xlabel('Time'), ylabel('Amplitude');
hold on;
Segments = cut_signal(filteredSignal(:,1), 2048, 2048);

E = short_time_energy(filteredSignal(:,1), 2048, 2048);
DE = 10*log10(E); %energy in decibel
E_peaks = points_of_interest(DE);
DWE = detect_window_digits(E_peaks);
%figure, plot(DE, '-r'); hold on; plot(E_peaks, 'linestyle', 'none', 'marker','*'); hold on;

figure("Name", 'Short time energy in decibel perceptible by humans');
plot(95 + 10*log10(E));                          %energy in decibel
hold on;

ZCR = zcr(filteredSignal(:,1), 2048, 2048);
DZCR = 10*log10(ZCR);   %zero crossing rate in decibel
ZCR_peaks = points_of_interest(DZCR);
DWZ = detect_window_digits(ZCR_peaks);
%figure, plot(DZCR, '-b'); hold on; plot(ZCR_peaks, 'linestyle', 'none', 'marker','*'); hold on;

Digits = get_digits(DWE, DWZ, Segments);
%sound(Digits{1,6},fs); to sound last digit, the 6th

[iso, labels] = template_digits();
temp = [];
el = length(Digits{1,1});
ISO = [9,el];
for i=1:9
    iel = length(iso{i,1});
    if(el>iel)
        %temp = zeros(1, el - iel);
        temp = NaN(1, el - iel);
    end
    for j=1:length(iso{i,1})
        ISO(i,j) = iso{i,1}(1,j);
    end
    ISO(i,j+1:el) = temp; %horzcat(ISO(i,:), temp);
    temp = [];
end

Digit = reshape(Digits{1,1}(:,1),1,el);
model = fitcecoc(ISO, labels);
[result, score] = predict(model,Digit);


function Segments = cut_signal(X, N, L) %signal, window, overlap
    Segments = [];
    m = 0;
    while m * L + N-1 + 1 <= length(X) %while current window has not reached end of signal (is not the last one)
        Segments = [ Segments X(m*L+1:m*L+N-1+1)];
        m = m + 1;
    end
end

function E = short_time_energy(X, N, L) %signal, segment, overlap
    m=0;
    E=[];
    while m * L + N-1 + 1 <= length(X) %while current window has not reached end of signal (is not the last one)
        E = [ E sum( X(m*L+1:m*L+N-1+1).^2)/N];
        m = m + 1;
    end
end

function ZCR = zcr(X, N, L)
    ZCR = [];
    m = 0;
    while m * L + N-1 + 1 <= length(X)
        odd = X(m*L+1:2:m*L+N-1+1);
        even = X(m*L+2:2:m*L+N-1+1);
        ZCR = [ ZCR sum( abs(sgn(even) - sgn(odd)))/(2*(N-1))];
        m = m + 1;
    end
end

function s = sgn(x)
    s = 1*(x>=0) + (-1)*(x<0);
end

function P = points_of_interest(Y) %input in decibel
    P = [];
    for i=1:length(Y)-1
        if( floor(abs(minus(Y(i), Y(i+1)))) > 0)
            P(i) = Y(i);
        else
            P(i) = NaN;
        end
    end
end

function D = detect_window_digits(ZCR_peaks)
    D = {}; %return cell array, #rows = #digits, each row contains window identifiers per digit
    indexes = []; %column positions of ZCR_peaks when value exists
    j = 1;
    for i=1:length(ZCR_peaks)
        if(not(isnan(ZCR_peaks(i))))
            indexes(j) = i;
            j = j + 1;
        end
    end
    %--------------------------------------------------------------
    diff = []; %difference between ZCR peaks to detect voiced areas
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

function Digits = get_digits(DWE, DWZ, Segments)
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

function [iso, labels] = template_digits()
    iso = {9};
    labels = {9};
    for i=1:9
        filename = sprintf('%i.wav',i);
        file = fullfile('.\isolation',filename);
        x = audioread(file);
        [rows, columns] = size(x);
        x = reshape(x, 1, rows*columns);
        [~,name,~] = fileparts(file);
        iso{i,:} = x;                               %{x, name};
        labels{i,:} = name;
    end
end