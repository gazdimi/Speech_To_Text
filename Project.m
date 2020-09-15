clear, clc;

[x, fs] = audioread('Samples\6_numbers.wav');
time = (1:length(x))*1/fs;  %time in seconds

Hd = humanEar(); filteredSignal = filter(Hd, x);
figure("Name", 'Original signal vs Filtered signal'); plot(time, x);
hold on; plot(time, filteredSignal), xlabel('Time'), ylabel('Amplitude');
hold on;

E = short_time_energy(filteredSignal(:,1), 2048, 1024); %4096, 2048
DE = 10*log10(E); %energy in decibel
E_peaks = points_of_interest(DE);
%figure, plot(DE, '-r'); hold on; plot(E_peaks, 'linestyle', 'none', 'marker','*'); hold on;

figure("Name", 'Short time energy in decibel perceptible by humans');
plot(95 + 10*log10(E));                          %energy in decibel
hold on;

ZCR = zcr(filteredSignal(:,1), 2048, 1024); %4096, 2048
DZCR = 10*log10(ZCR);   %zero crossing rate in decibel
ZCR_peaks = points_of_interest(DZCR);
%figure, plot(DZCR, '-b'); hold on; plot(ZCR_peaks, 'linestyle', 'none', 'marker','*'); hold on;

indexes = []; %column positions of ZCR_peaks when value exists
j = 1;
for i=1:length(ZCR_peaks)
     if(not(isnan(ZCR_peaks(i))))
         indexes(j) = i;
         j = j + 1;
     end
end
diff = []; %difference between ZCR peaks to detect voiced areas
for i=1:length(indexes)-1
    diff(i)=indexes(i+1)-indexes(i);
end
%S={{1,2,3,7,8,10}; {4,5,6}; {7,8}}; | S{3,1} | S{3,1}{:} example with cell
%array
digits = {};
j=1;
digits{1,1}{1} = indexes(1);
k = 2;
for i=2:length(indexes)
    if(diff(i-1)>round(max(diff)/2))
        j = j + 1;
        k = 1;
    end 
    digits{j,1}{k} = indexes(i);
    k = k + 1;
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