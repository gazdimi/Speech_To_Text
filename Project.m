clear, clc;

[x, fs] = audioread('Samples\6_numbers.wav');
time = (1:length(x))*1/fs;  %time in seconds

Hd = humanEar(); filteredSignal = filter(Hd, x);
figure("Name", 'Original signal vs Filtered signal'); plot(time, x);
hold on; plot(time, filteredSignal), xlabel('Time'), ylabel('Amplitude');
hold on;

E = short_time_energy(filteredSignal(:,1), 4096, 2048);

figure("Name", 'Short time energy plot');
plot(95 + 10*log10(E));                          %energy in decibel
hold on;

ZCR = zcr(filteredSignal(:,1), 4096, 2048);

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
        ZCR = [ ZCR sum( abs(sgn(odd)-sgn(even)))/(2*(N-1))];
        m = m + 1;
    end
end

function s = sgn(x)
    s = 1*(x>=0) + (-1)*(x<0);
end