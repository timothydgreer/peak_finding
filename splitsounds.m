clc
close all
clear
tic
%read in the wav file
[d, sr]=wavread('INSERT YOUR FILE HERE!');
s = 2^8+1000; % set s to be arbitrarily high.
x = 1; %user input, ask for value between 1.5=double and 0 (not inclusive)
% Goes here if stereo audio file
if length(d(2,:))==2
        dl=d(:,1); % Left channel
        dr=d(:,2); % Right channel
        % Create spectrograms of the left and right channel
        [Sly Slf Slt SlP] = (spectrogram(dl,hanning(s),(2-x)*length(hanning(s))/2,s));
        [Sry Srf Srt SrP] = (spectrogram(dr,hanning(s),(2-x)*length(hanning(s))/2,s));
        % Plots just the left channel. Easy to do for right channel as well.
        surf(Slt/(Slt(length(Slt))*44100/length(dl)),Slf*21106/pi,10*log10(abs(SlP)),'EdgeColor','none');
        axis xy; axis tight; colormap(jet); view(0,90);
        xlabel('Time');
        ylabel('Frequency (Hz)');
        hold off
        % Make a temporary array of peaks for the left channel 
        peaksl = findpeaks(abs(real(Sly(:,5))));
        % Find the highest peak
        maxmuml = max(peaksl)
        % Filter out peaks that are less than 10% of the highest peak
        peaksl = peaksl(find(peaksl>=maxmuml*.1))
        % Make an array that will store a spectrogram into each entry
        temparrayl = zeros(length(peaksl),length(Sly(:,1)),length(Sly(1,:)));
        
        %go through all frequencies
        for j = 3:length(peaksl) 
            tempo = 0;
            %find the index of the peak
            ind = find(abs(real(Sly(:,5)))==peaksl(j));
            % then filter all sound but 2 rows up and down from the peak
            for i = ind-2:ind+2
                if(abs(real(Sly(i,5)))>=0.5*peaksl(j))
                    tempo = tempo+1;
                    tempend = i;
                end
            end
            %save to temp array
            temparrayl(j,tempend-tempo+1:tempend,:) = Sly(tempend-tempo+1:tempend,:);
        end
        %Go through the same process for the right channel
        figure
        surf(Srt/(Srt(length(Srt))*44100/length(dr)),Srf*21106/pi,10*log10(abs(SrP)),'EdgeColor','none');
        axis xy; axis tight; colormap(jet); view(0,90);
        xlabel('Time');
        ylabel('Frequency (Hz)');
        hold off
        tempr = zeros(length(Sry(:,1)),length(Sry(1,:)));
        peaksr = findpeaks(abs(real(Sry(:,5))));
        maxmumr = max(peaksr)
        peaksr = peaksr(find(peaksr>=maxmumr*.1))
        temparrayr = zeros(length(peaksr),length(Sry(:,1)),length(Sry(1,:)));
        for j = 3:length(peaksr) 
            tempo = 0;
            ind = find(abs(real(Sry(:,5)))==peaksr(j));
            for i = ind-2:ind+2
                if(abs(real(Sry(i,5)))>=0.5*peaksr(j))
                    tempo = tempo+1;
                    tempend = i;
                end
            end
            temparrayr(j,tempend-tempo+1:tempend,:) = Sry(tempend-tempo+1:tempend,:);
        end
        % Now want to pick one peak to output
        l1 = zeros(length(Sly(:,1)),length(Sly(1,:)));
        r1 = zeros(length(Sry(:,1)),length(Sry(1,:)));
        l1(:,:) = temparrayl(1,:,:);
        r1(:,:) = temparrayr(1,:,:);
        
        % ispecgram the peaks 
        Zl = ispecgram(l1,s,0);
        Zr = ispecgram(r1,s,0);
        %Make the sound file
        Z = [Zl,Zr];
        Q = audioplayer(Z, sr);
        toc
        play(Q)
else
        % Mono case. It is the same process as above, but assumes that
        % there is one channel instead of 2. For documentation, see the
        % left channel comments above
        [Sy Sf St SP] = (spectrogram(d,hanning(s),(2-x)*length(hanning(s))/2,s));
        surf(St/(St(length(St))*44100/length(d)),Sf*21106/pi,10*log10(abs(SP)),'EdgeColor','none');
        axis xy; axis tight; colormap(jet); view(0,90);
        xlabel('Time');
        ylabel('Frequency (Hz)');    
        temp = zeros(length(Sy(:,1)),length(Sy(1,:)));
        peaks = findpeaks(abs(real(Sy(:,1))));
        maxmum = max(peaks)
        peaks = peaks(find(peaks>=maxmum*.1))
        temparray = zeros(length(peaks),length(Sy(:,1)),length(Sy(1,:)));
        for j = 3:length(peaks) 
            tempo = 0;
            ind = find(abs(real(Sy(:,1)))==peaks(j));
            for i = ind-1:ind+1
                if(abs(real(Sy(i,1)))>=0.5*peaks(j))
                    tempo = tempo+1;
                    tempend = i;
                end
            end
            temparray(j,tempend-tempo+1:tempend,:) = Sy(tempend-tempo+1:tempend,:);
        end
        q1 = zeros(length(Sy(:,1)),length(Sy(1,:)));
        q1(:,:) = temparray(1,:,:);
        surf(St/(St(length(St))*44100/length(d)),Sf*21106/pi,10*log10(abs(q3)),'EdgeColor','none');
        axis xy; axis tight; colormap(jet); view(0,90);
        xlabel('Time');
        ylabel('Frequency (Hz)');        
        Z = ispecgram(q1,s,0);
        Q = audioplayer(Z, sr);
        toc
        play(Q)
end
