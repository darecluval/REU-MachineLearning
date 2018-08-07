%% 1 Detache Down - Energy Approach
DetacheDown = audioread('Dmajor_detache_down.wav');
 
%Variables
fs = 44100;
dt = 1/fs;
[s,~] = size(DetacheDown);
secs =  s / fs;
 
%Tolerance Level
tolerance = 0.60;
 
%Notes of interest
notes = [294 330 370 392 440 494 554 587];
 
%***********************Hanning entire Audiofile************************%
 
%Hanning and fft of the entire window
h = hanning(s);
DDh = (DetacheDown .* h);
DD = abs(fft(DDh));
 
%Amp vs. Freq Graph
plot(DD(1:s/2)); 
 
%Table that has interest values and precise frequency
Notes = zeros(numel(notes), 2);
 
%Finds the peaks and saves their frequency into Notes[]
for i = 1:numel(notes)
    first = round(notes(i)*secs-100);
    last = round(notes(i)*secs+100);
    [val, idx] = max(DD(first:last));
    y1 = find(DD==val);
    Notes(i, 1) = val;
    Notes(i, 2) = y1(1)/secs;
end
 
%Allocates memory for big static arrays
Co = zeros(s, numel(notes));
Corr = zeros(s+2000, 1);
Freq = zeros(s+100, 1);
N = zeros(3, 75);
 
%Sets note of interest
n = 1;
 
%*******************Cross Correlation of each note**********************%
 
%Fills cross correlation vector
while n <= numel(notes)
    
    %Will only look at frequencies of interest
    if Notes(n, 1) >= 1
        
        %Beginning Harmonics
        k = 1;
        
        while k <= 6 && (k*Notes(n,2) <= fs/2)
     
            %Creates a cosine wave matching the frequency
            T = 1/Notes(n, 2);
            intv = 0:dt:T+dt;
            t = 0:dt:10/round(Notes(n, 2));
            c = (cos(2*pi*k*round(Notes(n, 2))*t))';
        
            %Normalizes the cosine wave to mean = 0 and sum(energy) = 1
            C = sum(c(:).^2);
            c(:) = c(:) ./ sqrt(C);
        
            %Variables
            [lc, ~] = size(c);
            first = 1;
            last = lc;
            corr = zeros(lc, 1);
            DNorm = zeros(lc, 1);
            DetacheDown(end+1:end+lc) = 0;
        
            %Manual xcorr
            while first <= s
            
                %Normalizes the interval to sum(energy) = 1
                DNorm(:) = DetacheDown(first:last) - ...
                    mean(DetacheDown(first:last));
                DN = sum(DNorm.^2);
                DNorm(:) = DNorm(:) ./ sqrt(DN);
            
                %Determines the correlation coefficient
                corr(:) = DNorm(:) .* c(:);
                Corr(first,1) = sum(corr(:));
            
                %Increments to next window
                first = first + 1;
                last = last + 1;
            end
        
            %Saves the energy of the note (no harmonics)
            Co(:,n) = Co(:,n) + Corr(1:s).^2;
            
            %Progresses to next harmonic
            k = k + 1;
        end
        
        %Determines the max corr value of the note+harmonics
        [Max, idx] = max(Co(:,n));
        %*******************75 Percent Threshold************************%
        
        for i = 1:s
            if Co(i,n) < tolerance * Max
               Co(i,n) = 0;
            end
        end
        
    %Empties notes not of interest in the .wav
    else 
        Co(:,n) = 0;
    end
    n = n + 1;
end
 
i = 1;
 
while i <= s
    n = 1;
    while n <= numel(notes)
        while max(Co(i,:)) == 0 && i < s
            i = i + 1;
        end
        if Co(i,n) ~= max(Co(i,:))
            Co(i,n) = 0;
        else
            Co(i,n) = Notes(n,2);
        end
        n = n + 1;
    end
    Freq(i,1) = max(Co(i,:));
    i = i + 1;
end
 
%*****************************Clean Gaps********************************%
 
%Variables
i = 1;
j = 1;
count = 1;
 
%Cleans smallest gaps, observed to be 150 samples 
while i <= s
    
    %Skips big intervals of 0
    if Freq(i, 1) == 0
        i = i + 1; 
    else
        if i ~= s
            
            %If the next interval drops
            if Freq(i+1, 1) ~= Freq(i,1)
                z = 1;
                
                %Checks within 75 samples for a match
                while z < 18500 
                    if (Freq(i+z, 1) == Freq(i,1))
                        Freq(i:i+z,1) = Freq(i, 1);
                        z = 18501;
                    else 
                        z = z + 1;
                    end
                end
                i = i + 1;
                
            %Moves on
            else 
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end  
end
 
Freq(end+1,1) = 1;
[s,~] = size(Freq);
 
%Creates an matrix of intervals
while j <= s
    
   while Freq(j,1) == 0 && j<s
       j = j + 1;
   end
   
   if j < s-1
        %Saves the beginning of the interval 
        start = j;
        no = Freq(j,1);
        
        %Moves through the interval
        while Freq(j,1) == Freq(j+1, 1)
            j = j + 1;
        end
    
        stop = j;
    
        if j < s && (stop-start >= 7500)
        
            %Saves the note's interval
            N(:,count) = [no start stop];
            count = count + 1; 
        end
   end
    j = j + 1;
end
 
%Variables
[~, End] = size(N);
count = 1;
i = 1;
int = 128;
int2 = int/2;
int4 = int/4;
 
%**********Saves the parsed notes in dynamically named files************%
while i < End && N(1,i) ~= 0
    
    %Skips empty intervals
    while N(1, i) == 0 && i ~= End
        i = i + 1;
    end
    
    %Beginning of interval
    start = N(2, i);
    stop = N(3, i);
    
    %Creates wav and jpeg file names
    basename = round(N(1, i));
    filename = ['1_' num2str(count) '.wav'];
    imagename = ['1_' num2str(count) '.jpeg'];
    
    %Saves audio snippit
    audiowrite(filename, DetacheDown(start:stop), fs);
    save filename
    
    %Hanning Window
    w = hanning(int);
    
    %Variables
    first = 1;
    last = int;
    
    %Determines dimensions of the image
    length = stop-start+1;
    numWindows = ceil((length-int) / int4);
    newEnd = (numWindows * int4) + int;
    imageInt = zeros(newEnd, 1);
    imageInt(1:stop-start+1) = DetacheDown(start:stop);
    image = zeros(int2, numWindows);
  
    window = 1;
    %Creates spectrogram from the given note
    while window <= numWindows
        interval = imageInt(first:last);
        x = (interval .* w);
        x = abs(fft(x));
        z = x(1:int2);
        
        image(:, window) = z';
        
        window = window + 1;
        first = first + int4;
        last = last + int4;
    end
    
    %Scale image to 200x200 and SAVE
    m = max(image);
    image = image ./ m;
    Image = imresize(image, [200 200]);
    Image = flipud(Image);
    imwrite(Image, imagename);
    save imagename;
    
    count = count + 1;
    i = i + 1;
end
 
 
%End
done = 'detache down done'
 
%% 2 Detache Up
DetacheUp = audioread('Dmajor_detache_up.wav');
 
%Variables
fs = 44100;
dt = 1/fs;
[s,~] = size(DetacheUp);
secs =  s / fs;
 
%Tolerance Level
tolerance = 0.60;
 
%Notes of interest
notes = [294 330 370 392 440 494 554 587];
 
%***********************Hanning entire Audiofile************************%
 
%Hanning and fft of the entire window
h = hanning(s);
DDh = (DetacheUp .* h);
DD = abs(fft(DDh));
 
%Amp vs. Freq Graph
plot(DD(1:s/2)); 
 
%Table that has interest values and precise frequency
Notes = zeros(numel(notes), 2);
 
%Finds the peaks and saves their frequency into Notes[]
for i = 1:numel(notes)
    first = round(notes(i)*secs-100);
    last = round(notes(i)*secs+100);
    [val, idx] = max(DD(first:last));
    y1 = find(DD==val);
    Notes(i, 1) = val;
    Notes(i, 2) = y1(1)/secs;
end
 
%Allocates memory for big static arrays
Co = zeros(s, numel(notes));
Corr = zeros(s+2000, 1);
Freq = zeros(s+100, 1);
N = zeros(3, 75);
 
%Sets note of interest
n = 1;
 
%*******************Cross Correlation of each note**********************%
 
%Fills cross correlation vector
while n <= numel(notes)
    
    %Will only look at frequencies of interest
    if Notes(n, 1) >= 1
        
        %Beginning Harmonics
        k = 1;
        
        while k <= 6 && (k*Notes(n,2) <= fs/2)
     
            %Creates a cosine wave matching the frequency
            T = 1/Notes(n, 2);
            intv = 0:dt:T+dt;
            t = 0:dt:10/round(Notes(n, 2));
            c = (cos(2*pi*k*round(Notes(n, 2))*t))';
        
            %Normalizes the cosine wave to mean = 0 and sum(energy) = 1
            C = sum(c(:).^2);
            c(:) = c(:) ./ sqrt(C);
        
            %Variables
            [lc, ~] = size(c);
            first = 1;
            last = lc;
            corr = zeros(lc, 1);
            DNorm = zeros(lc, 1);
            DetacheUp(end+1:end+lc) = 0;
        
            %Manual xcorr
            while first <= s
            
                %Normalizes the interval to sum(energy) = 1
                DNorm(:) = DetacheUp(first:last) - ...
                    mean(DetacheUp(first:last));
                DN = sum(DNorm.^2);
                DNorm(:) = DNorm(:) ./ sqrt(DN);
            
                %Determines the correlation coefficient
                corr(:) = DNorm(:) .* c(:);
                Corr(first,1) = sum(corr(:));
            
                %Increments to next window
                first = first + 1;
                last = last + 1;
            end
        
            %Saves the energy of the note (no harmonics)
            Co(:,n) = Co(:,n) + Corr(1:s).^2;
            
            %Progresses to next harmonic
            k = k + 1;
        end
        
        %Determines the max corr value of the note+harmonics
        [Max, idx] = max(Co(:,n));
        %*******************75 Percent Threshold************************%
        
        for i = 1:s
            if Co(i,n) < tolerance * Max
               Co(i,n) = 0;
            end
        end
        
    %Empties notes not of interest in the .wav
    else 
        Co(:,n) = 0;
    end
    n = n + 1;
end
 
i = 1;
 
while i <= s
    n = 1;
    while n <= numel(notes)
        while max(Co(i,:)) == 0 && i < s
            i = i + 1;
        end
        if Co(i,n) ~= max(Co(i,:))
            Co(i,n) = 0;
        else
            Co(i,n) = Notes(n,2);
        end
        n = n + 1;
    end
    Freq(i,1) = max(Co(i,:));
    i = i + 1;
end
 
%*****************************Clean Gaps********************************%
 
%Variables
i = 1;
j = 1;
count = 1;
 
%Cleans smallest gaps, observed to be 150 samples 
while i <= s
    
    %Skips big intervals of 0
    if Freq(i, 1) == 0
        i = i + 1; 
    else
        if i ~= s
            
            %If the next interval drops
            if Freq(i+1, 1) ~= Freq(i,1)
                z = 1;
                
                %Checks within 75 samples for a match
                while z < 18500 
                    if (Freq(i+z, 1) == Freq(i,1))
                        Freq(i:i+z,1) = Freq(i, 1);
                        z = 18501;
                    else 
                        z = z + 1;
                    end
                end
                i = i + 1;
                
            %Moves on
            else 
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end  
end
 
Freq(end+1,1) = 1;
[s,~] = size(Freq);
 
%Creates an matrix of intervals
while j <= s
    
   while Freq(j,1) == 0 && j<s
       j = j + 1;
   end
   
   if j < s-1
        %Saves the beginning of the interval 
        start = j;
        no = Freq(j,1);
        
        %Moves through the interval
        while Freq(j,1) == Freq(j+1, 1)
            j = j + 1;
        end
    
        stop = j;
    
        if j < s && (stop-start >= 7500)
        
            %Saves the note's interval
            N(:,count) = [no start stop];
            count = count + 1; 
        end
   end
    j = j + 1;
end
 
%Variables
[~, End] = size(N);
count = 1;
i = 1;
int = 128;
int2 = int/2;
int4 = int/4;
 
%**********Saves the parsed notes in dynamically named files************%
while i < End && N(1,i) ~= 0
    
    %Skips empty intervals
    while N(1, i) == 0 && i ~= End
        i = i + 1;
    end
    
    %Beginning of interval
    start = N(2, i);
    stop = N(3, i);
    
    %Creates wav and jpeg file names
    basename = round(N(1, i));
    filename = ['2_' num2str(count) '.wav'];
    imagename = ['2_' num2str(count) '.jpeg'];
    
    %Saves audio snippit
    audiowrite(filename, DetacheUp(start:stop), fs);
    save filename
    
    %Hanning Window
    w = hanning(int);
    
    %Variables
    first = 1;
    last = int;
    
    %Determines dimensions of the image
    length = stop-start+1;
    numWindows = ceil((length-int) / int4);
    newEnd = (numWindows * int4) + int;
    imageInt = zeros(newEnd, 1);
    imageInt(1:stop-start+1) = DetacheUp(start:stop);
    image = zeros(int2, numWindows);
    
    window = 1;
    %Creates spectrogram from the given note
    while window <= numWindows
        interval = imageInt(first:last);
        x = (interval .* w);
        x = abs(fft(x));
        z = x(1:int2);
        
        image(:, window) = z';
        
        window = window + 1;
        first = first + int4;
        last = last + int4;
    end
    
    %Scale image to 200x200 and SAVE
    m = max(image);
    image = image ./ m;
    Image = imresize(image, [200 200]);
    Image = flipud(Image);
    imwrite(Image, imagename);
    save imagename;
    
    count = count + 1;
    i = i + 1;
end
 
 
%End
done = 'detache up done'
 
%% 3 Slur 2 - Energy Approach
Slur2 = audioread('Dmajor_slur2.wav');
 
%Variables
fs = 44100;
dt = 1/fs;
[s,~] = size(Slur2);
secs =  s / fs;
 
%Tolerance Level
tolerance = 0.60;
 
%Notes of interest
notes = [294 330 370 392 440 494 554 587];
 
%***********************Hanning entire Audiofile************************%
 
%Hanning and fft of the entire window
h = hanning(s);
DDh = (Slur2 .* h);
DD = abs(fft(DDh));
 
%Amp vs. Freq Graph
plot(DD(1:s/2)); 
 
%Table that has interest values and precise frequency
Notes = zeros(numel(notes), 2);
 
%Finds the peaks and saves their frequency into Notes[]
for i = 1:numel(notes)
    first = round(notes(i)*secs-100);
    last = round(notes(i)*secs+100);
    [val, idx] = max(DD(first:last));
    y1 = find(DD==val);
    Notes(i, 1) = val;
    Notes(i, 2) = y1(1)/secs;
end
 
%Allocates memory for big static arrays
Co = zeros(s, numel(notes));
Corr = zeros(s+2000, 1);
Freq = zeros(s+100, 1);
N = zeros(3, 75);
 
%Sets note of interest
n = 1;
 
%*******************Cross Correlation of each note**********************%
 
%Fills cross correlation vector
while n <= numel(notes)
    
    %Will only look at frequencies of interest
    if Notes(n, 1) >= 1
        
        %Beginning Harmonics
        k = 1;
        
        while k <= 6 && (k*Notes(n,2) <= fs/2)
     
            %Creates a cosine wave matching the frequency
            T = 1/Notes(n, 2);
            intv = 0:dt:T+dt;
            t = 0:dt:10/round(Notes(n, 2));
            c = (cos(2*pi*k*round(Notes(n, 2))*t))';
        
            %Normalizes the cosine wave to mean = 0 and sum(energy) = 1
            C = sum(c(:).^2);
            c(:) = c(:) ./ sqrt(C);
        
            %Variables
            [lc,~] = size(c);
            first = 1;
            last = lc;
            corr = zeros(lc, 1);
            DNorm = zeros(lc, 1);
            Slur2(end+1:end+lc) = 0;
        
            %Manual xcorr
            while first <= s
            
                %Normalizes the interval to sum(energy) = 1
                DNorm(:) = Slur2(first:last) - ...
                    mean(Slur2(first:last));
                DN = sum(DNorm.^2);
                DNorm(:) = DNorm(:) ./ sqrt(DN);
            
                %Determines the correlation coefficient
                corr(:) = DNorm(:) .* c(:);
                Corr(first,1) = sum(corr(:));
            
                %Increments to next window
                first = first + 1;
                last = last + 1;
            end
        
            %Saves the energy of the note (no harmonics)
            Co(:,n) = Co(:,n) + Corr(1:s).^2;
            
            %Progresses to next harmonic
            k = k + 1;
        end
        
        %Determines the max corr value of the note+harmonics
        [Max, idx] = max(Co(:,n));
        %*******************75 Percent Threshold************************%
        
        for i = 1:s
            if Co(i,n) < tolerance * Max
               Co(i,n) = 0;
            end
        end
        
    %Empties notes not of interest in the .wav
    else 
        Co(:,n) = 0;
    end
    n = n + 1;
end
 
i = 1;
 
while i <= s
    n = 1;
    while n <= numel(notes)
        while max(Co(i,:)) == 0 && i < s
            i = i + 1;
        end
        if Co(i,n) ~= max(Co(i,:))
            Co(i,n) = 0;
        else
            Co(i,n) = Notes(n,2);
        end
        n = n + 1;
    end
    Freq(i,1) = max(Co(i,:));
    i = i + 1;
end
 
%*****************************Clean Gaps********************************%
 
%Variables
i = 1;
j = 1;
count = 1;
 
%Cleans smallest gaps, observed to be 150 samples 
while i <= s
    
    %Skips big intervals of 0
    if Freq(i, 1) == 0
        i = i + 1; 
    else
        if i ~= s
            
            %If the next interval drops
            if Freq(i+1, 1) ~= Freq(i,1)
                z = 1;
                
                %Checks within 75 samples for a match
                while z < 18500 
                    if (Freq(i+z, 1) == Freq(i,1))
                        Freq(i:i+z,1) = Freq(i, 1);
                        z = 18501;
                    else 
                        z = z + 1;
                    end
                end
                i = i + 1;
                
            %Moves on
            else 
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end  
end
 
Freq(end+1,1) = 1;
[s,~] = size(Freq);
 
%Creates an matrix of intervals
while j <= s
    
    while Freq(j,1) == 0 && j<s
        j = j + 1;
    end
   
    if j < s-1
        %Saves the beginning of the interval 
        start = j;
        no = Freq(j,1);
        
        %Moves through the interval
        while Freq(j,1) == Freq(j+1, 1)
            j = j + 1;
        end
    
        stop = j;
    
        if j < s && (stop-start >= 7500)
        
            %Saves the note's interval
            N(:,count) = [no start stop];
            count = count + 1; 
        end
    end
    j = j + 1;
end
 
%Variables
[~, End] = size(N);
count = 1;
i = 1;
int = 128;
int2 = int/2;
int4 = int/4;
 
%**********Saves the parsed notes in dynamically named files************%
while i < End && N(1,i) ~= 0
    
    %Skips empty intervals
    while N(1, i) == 0 && i ~= End
        i = i + 1;
    end
    
    %Beginning of interval
    start = N(2, i);
    stop = N(3, i);
    
    %Creates wav and jpeg file names
    basename = round(N(1, i));
    filename = ['3_' num2str(count) '.wav'];
    imagename = ['3_' num2str(count) '.jpeg'];
    
    %Saves audio snippit
    audiowrite(filename, Slur2(start:stop), fs);
    save filename
    
    %Hanning Window
    w = hanning(int);
    
    %Variables
    first = 1;
    last = int;
    
    %Determines dimensions of the image
    length = stop-start+1;
    numWindows = ceil((length-int) / int4);
    newEnd = (numWindows * int4) + int;
    imageInt = zeros(newEnd, 1);
    imageInt(1:stop-start+1) = Slur2(start:stop);
    image = zeros(int2, numWindows);
    
    window = 1;
    %Creates spectrogram from the given note
    while window <= numWindows
        interval = imageInt(first:last);
        x = (interval .* w);
        x = abs(fft(x));
        z = x(1:int2);
        
        image(:, window) = z';
        
        window = window + 1;
        first = first + int4;
        last = last + int4;
    end
    
    %Scale image to 200x200 and SAVE
    m = max(image);
    image = image ./ m;
    Image = imresize(image, [200 200]);
    Image = flipud(Image);
    imwrite(Image, imagename);
    save imagename;
    
    count = count + 1;
    i = i + 1;
end
 
 
%End
done = 'slur 2 done'
%% 4 Slur 4 - Energy Approach
Slur4 = audioread('Dmajor_slur4.wav');
 
%Variables
fs = 44100;
dt = 1/fs;
[s,~] = size(Slur4);
secs =  s / fs;
 
%Tolerance Level
tolerance = 0.60;
 
%Notes of interest
notes = [294 330 370 392 440 494 554 587];
 
%***********************Hanning entire Audiofile************************%
 
%Hanning and fft of the entire window
h = hanning(s);
DDh = (Slur4 .* h);
DD = abs(fft(DDh));
 
%Amp vs. Freq Graph
plot(DD(1:s/2)); 
 
%Table that has interest values and precise frequency
Notes = zeros(numel(notes), 2);
 
%Finds the peaks and saves their frequency into Notes[]
for i = 1:numel(notes)
    first = round(notes(i)*secs-100);
    last = round(notes(i)*secs+100);
    [val, idx] = max(DD(first:last));
    y1 = find(DD==val);
    Notes(i, 1) = val;
    Notes(i, 2) = y1(1)/secs;
end
 
%Allocates memory for big static arrays
Co = zeros(s, numel(notes));
Corr = zeros(s+2000, 1);
Freq = zeros(s+100, 1);
N = zeros(3, 75);
 
%Sets note of interest
n = 1;
 
%*******************Cross Correlation of each note**********************%
 
%Fills cross correlation vector
while n <= numel(notes)
    
    %Will only look at frequencies of interest
    if Notes(n, 1) >= 1
        
        %Beginning Harmonics
        k = 1;
        
        while k <= 6 && (k*Notes(n,2) <= fs/2)
     
            %Creates a cosine wave matching the frequency
            T = 1/Notes(n, 2);
            intv = 0:dt:T+dt;
            t = 0:dt:10/round(Notes(n, 2));
            c = (cos(2*pi*k*round(Notes(n, 2))*t))';
        
            %Normalizes the cosine wave to mean = 0 and sum(energy) = 1
            C = sum(c(:).^2);
            c(:) = c(:) ./ sqrt(C);
        
            %Variables
            [lc, ~] = size(c);
            first = 1;
            last = lc;
            corr = zeros(lc, 1);
            DNorm = zeros(lc, 1);
            Slur4(end+1:end+lc) = 0;
        
            %Manual xcorr
            while first <= s
            
                %Normalizes the interval to sum(energy) = 1
                DNorm(:) = Slur4(first:last) - ...
                    mean(Slur4(first:last));
                DN = sum(DNorm.^2);
                DNorm(:) = DNorm(:) ./ sqrt(DN);
            
                %Determines the correlation coefficient
                corr(:) = DNorm(:) .* c(:);
                Corr(first,1) = sum(corr(:));
            
                %Increments to next window
                first = first + 1;
                last = last + 1;
            end
        
            %Saves the energy of the note (no harmonics)
            Co(:,n) = Co(:,n) + Corr(1:s).^2;
            
            %Progresses to next harmonic
            k = k + 1;
        end
        
        %Determines the max corr value of the note+harmonics
        [Max, idx] = max(Co(:,n));
        %*******************75 Percent Threshold************************%
        
        for i = 1:s
            if Co(i,n) < tolerance * Max
               Co(i,n) = 0;
            end
        end
        
    %Empties notes not of interest in the .wav
    else 
        Co(:,n) = 0;
    end
    n = n + 1;
end
 
i = 1;
 
while i <= s
    n = 1;
    while n <= numel(notes)
        while max(Co(i,:)) == 0 && i < s
            i = i + 1;
        end
        if Co(i,n) ~= max(Co(i,:))
            Co(i,n) = 0;
        else
            Co(i,n) = Notes(n,2);
        end
        n = n + 1;
    end
    Freq(i,1) = max(Co(i,:));
    i = i + 1;
end
 
%*****************************Clean Gaps********************************%
 
%Variables
i = 1;
j = 1;
count = 1;
 
%Cleans smallest gaps, observed to be 150 samples 
while i <= s
    
    %Skips big intervals of 0
    if Freq(i, 1) == 0
        i = i + 1; 
    else
        if i ~= s
            
            %If the next interval drops
            if Freq(i+1, 1) ~= Freq(i,1)
                z = 1;
                
                %Checks within 75 samples for a match
                while z < 18500 
                    if (Freq(i+z, 1) == Freq(i,1))
                        Freq(i:i+z,1) = Freq(i, 1);
                        z = 18501;
                    else 
                        z = z + 1;
                    end
                end
                i = i + 1;
                
            %Moves on
            else 
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end  
end
 
Freq(end+1,1) = 1;
[s,~] = size(Freq);
 
%Creates an matrix of intervals
while j <= s
    
   while Freq(j,1) == 0 && j<s
       j = j + 1;
   end
   
   if j < s-1
        %Saves the beginning of the interval 
        start = j;
        no = Freq(j,1);
        
        %Moves through the interval
        while Freq(j,1) == Freq(j+1, 1)
            j = j + 1;
        end
    
        stop = j;
    
        if j < s && (stop-start >= 7500)
        
            %Saves the note's interval
            N(:,count) = [no start stop];
            count = count + 1; 
        end
   end
    j = j + 1;
end
 
%Variables
[~, End] = size(N);
count = 1;
i = 1;
int = 128;
int2 = int/2;
int4 = int/4;
 
%**********Saves the parsed notes in dynamically named files************%
while i < End && N(1,i) ~= 0
    
    %Skips empty intervals
    while N(1, i) == 0 && i ~= End
        i = i + 1;
    end
    
    %Beginning of interval
    start = N(2, i);
    stop = N(3, i);
    
    %Creates wav and jpeg file names
    basename = round(N(1, i));
    filename = ['4_' num2str(count) '.wav'];
    imagename = ['4_' num2str(count) '.jpeg'];
    
    %Saves audio snippit
    audiowrite(filename, Slur4(start:stop), fs);
    save filename
    
    %Hanning Window
    w = hanning(int);
    
    %Variables
    first = 1;
    last = int;
    
    %Determines dimensions of the image
    length = stop-start+1;
    numWindows = ceil((length-int) / int4);
    newEnd = (numWindows * int4) + int;
    imageInt = zeros(newEnd, 1);
    imageInt(1:stop-start+1) = Slur4(start:stop);
    image = zeros(int2, numWindows);
    
    window = 1;
    %Creates spectrogram from the given note
    while window <= numWindows
        interval = imageInt(first:last);
        x = (interval .* w);
        x = abs(fft(x));
        z = x(1:int2);
        
        image(:, window) = z';
        
        window = window + 1;
        first = first + int4;
        last = last + int4;
    end
    
    %Scale image to 200x200 and SAVE
    m = max(image);
    image = image ./ m;
    Image = imresize(image, [200 200]);
    Image = flipud(Image);
    imwrite(Image, imagename);
    save imagename;
    
    count = count + 1;
    i = i + 1;
end
 
 
%End
done = 'slur 4 done'
 
 %% 5 Slur 8 - Energy Approach
Slur8 = audioread('Dmajor_slur8.wav');
 
%Variables
fs = 44100;
dt = 1/fs;
[s,~] = size(Slur8);
secs =  s / fs;
 
%Tolerance Level
tolerance = 0.65;
 
%Notes of interest
notes = [294 330 370 392 440 494 554 587];
 
%***********************Hanning entire Audiofile************************%
 
%Hanning and fft of the entire window
h = hanning(s);
DDh = (Slur8 .* h);
DD = abs(fft(DDh));
 
%Amp vs. Freq Graph
plot(DD(1:s/2)); 
 
%Table that has interest values and precise frequency
Notes = zeros(numel(notes), 2);
 
%Finds the peaks and saves their frequency into Notes[]
for i = 1:numel(notes)
    first = round(notes(i)*secs-100);
    last = round(notes(i)*secs+100);
    [val, idx] = max(DD(first:last));
    y1 = find(DD==val);
    Notes(i, 1) = val;
    Notes(i, 2) = y1(1)/secs;
end
 
%Allocates memory for big static arrays
Co = zeros(s, numel(notes));
Corr = zeros(s+2000, 1);
Freq = zeros(s+100, 1);
N = zeros(3, 75);
 
%Sets note of interest
n = 1;
 
%*******************Cross Correlation of each note**********************%
 
%Fills cross correlation vector
while n <= numel(notes)
    
    %Will only look at frequencies of interest
    if Notes(n, 1) >= 1
        
        %Beginning Harmonics
        k = 1;
        
        while k <= 6 && (k*Notes(n,2) <= fs/2)
     
            %Creates a cosine wave matching the frequency
            T = 1/Notes(n, 2);
            intv = 0:dt:T+dt;
            t = 0:dt:10/round(Notes(n, 2));
            c = (cos(2*pi*k*round(Notes(n, 2))*t))';
        
            %Normalizes the cosine wave to mean = 0 and sum(energy) = 1
            C = sum(c(:).^2);
            c(:) = c(:) ./ sqrt(C);
        
            %Variables
            [lc, ~] = size(c);
            first = 1;
            last = lc;
            corr = zeros(lc, 1);
            DNorm = zeros(lc, 1);
            Slur8(end+1:end+lc) = 0;
        
            %Manual xcorr
            while first <= s
            
                %Normalizes the interval to sum(energy) = 1
                DNorm(:) = Slur8(first:last) - ...
                    mean(Slur8(first:last));
                DN = sum(DNorm.^2);
                DNorm(:) = DNorm(:) ./ sqrt(DN);
            
                %Determines the correlation coefficient
                corr(:) = DNorm(:) .* c(:);
                Corr(first,1) = sum(corr(:));
            
                %Increments to next window
                first = first + 1;
                last = last + 1;
            end
        
            %Saves the energy of the note (no harmonics)
            Co(:,n) = Co(:,n) + Corr(1:s).^2;
            
            %Progresses to next harmonic
            k = k + 1;
        end
        
        %Determines the max corr value of the note+harmonics
        [Max, idx] = max(Co(:,n));
        %*******************75 Percent Threshold************************%
        
        for i = 1:s
            if n == 1
                if Co(i,n) < tolerance * Max
                    Co(i,n) = 0;
                end
            else
                if Co(i,n) < .60 * Max
                    Co(i,n) = 0;
                end
            end
        end
        
    %Empties notes not of interest in the .wav
    else 
        Co(:,n) = 0;
    end
    n = n + 1;
end
 
i = 1;
 
while i <= s
    n = 1;
    while n <= numel(notes)
        while max(Co(i,:)) == 0 && i < s
            i = i + 1;
        end
        if Co(i,n) ~= max(Co(i,:))
            Co(i,n) = 0;
        else
            Co(i,n) = Notes(n,2);
        end
        n = n + 1;
    end
    Freq(i,1) = max(Co(i,:));
    i = i + 1;
end
 
%*****************************Clean Gaps********************************%
 
%Variables
i = 1;
j = 1;
count = 1;
 
%Cleans smallest gaps, observed to be 150 samples 
while i <= s
    
    %Skips big intervals of 0
    if Freq(i, 1) == 0
        i = i + 1; 
    else
        if i ~= s
            
            %If the next interval drops
            if Freq(i+1, 1) ~= Freq(i,1)
                z = 1;
                
                %Checks within 75 samples for a match
                while z < 18500 
                    if (Freq(i+z, 1) == Freq(i,1))
                        Freq(i:i+z,1) = Freq(i, 1);
                        z = 18501;
                    else 
                        z = z + 1;
                    end
                end
                i = i + 1;
                
            %Moves on
            else 
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end  
end
 
Freq(end+1,1) = 1;
[s,~] = size(Freq);
 
%Creates an matrix of intervals
while j <= s
    
   while Freq(j,1) == 0 && j<s
       j = j + 1;
   end
   
   if j < s-1
        %Saves the beginning of the interval 
        start = j;
        no = Freq(j,1);
        
        %Moves through the interval
        while Freq(j,1) == Freq(j+1, 1)
            j = j + 1;
        end
    
        stop = j;
    
        if j < s && (stop-start >= 7500)
        
            %Saves the note's interval
            N(:,count) = [no start stop];
            count = count + 1; 
        end
   end
    j = j + 1;
end
 
%Variables
[~, End] = size(N);
count = 1;
i = 1;
int = 128;
int2 = int/2;
int4 = int/4;
 
%**********Saves the parsed notes in dynamically named files************%
while i < End && N(1,i) ~= 0
    
    %Skips empty intervals
    while N(1, i) == 0 && i ~= End
        i = i + 1;
    end
    
    %Beginning of interval
    start = N(2, i);
    stop = N(3, i);
    
    %Creates wav and jpeg file names
    basename = round(N(1, i));
    filename = ['5_' num2str(count) '.wav'];
    imagename = ['5_' num2str(count) '.jpeg'];
    
    %Saves audio snippit
    audiowrite(filename, Slur8(start:stop), fs);
    save filename
    
    %Hanning Window
    w = hanning(int);
    
    %Variables
    first = 1;
    last = int;
    
    %Determines dimensions of the image
    length = stop-start+1;
    numWindows = ceil((length-int) / int4);
    newEnd = (numWindows * int4) + int;
    imageInt = zeros(newEnd, 1);
    imageInt(1:stop-start+1) = Slur8(start:stop);
    image = zeros(int2, numWindows);
    
    window = 1;
    %Creates spectrogram from the given note
    while window <= numWindows
        interval = imageInt(first:last);
        x = (interval .* w);
        x = abs(fft(x));
        z = x(1:int2);
        
        image(:, window) = z';
        
        window = window + 1;
        first = first + int4;
        last = last + int4;
    end
    
    %Scale image to 200x200 and SAVE
    m = max(image);
    image = image ./ m;
    Image = imresize(image, [200 200]);
    Image = flipud(Image);
    imwrite(Image, imagename);
    save imagename;
    
    count = count + 1;
    i = i + 1;
end
 
 
%End
done = 'slur 8 done'
 
%% 6 Staccato - Energy Approach
Staccato = audioread('Dmajor_staccato.wav');
 
%Variables
fs = 44100;
dt = 1/fs;
[s,~] = size(Staccato);
secs =  s / fs;
 
%Tolerance Level
tolerance = 0.65;
 
%Notes of interest
notes = [294 330 370 392 440 494 554 587];
 
%***********************Hanning entire Audiofile************************%
 
%Hanning and fft of the entire window
h = hanning(s);
DDh = (Staccato .* h);
DD = abs(fft(DDh));
 
%Amp vs. Freq Graph
plot(DD(1:s/2)); 
 
%Table that has interest values and precise frequency
Notes = zeros(numel(notes), 2);
 
%Finds the peaks and saves their frequency into Notes[]
for i = 1:numel(notes)
    first = round(notes(i)*secs-100);
    last = round(notes(i)*secs+100);
    [val, idx] = max(DD(first:last));
    y1 = find(DD==val);
    Notes(i, 1) = val;
    Notes(i, 2) = y1(1)/secs;
end
 
%Allocates memory for big static arrays
Co = zeros(s, numel(notes));
Corr = zeros(s+2000, 1);
Freq = zeros(s+100, 1);
N = zeros(3, 75);
 
%Sets note of interest
n = 1;
 
%*******************Cross Correlation of each note**********************%
 
%Fills cross correlation vector
while n <= numel(notes)
    
    %Will only look at frequencies of interest
    if Notes(n, 1) >= 1
        
        %Beginning Harmonics
        k = 1;
        
        while k <= 6 && (k*Notes(n,2) <= fs/2)
     
            %Creates a cosine wave matching the frequency
            T = 1/Notes(n, 2);
            intv = 0:dt:T+dt;
            t = 0:dt:10/round(Notes(n, 2));
            c = (cos(2*pi*k*round(Notes(n, 2))*t))';
        
            %Normalizes the cosine wave to mean = 0 and sum(energy) = 1
            C = sum(c(:).^2);
            c(:) = c(:) ./ sqrt(C);
        
            %Variables
            [lc, ~] = size(c);
            first = 1;
            last = lc;
            corr = zeros(lc, 1);
            DNorm = zeros(lc, 1);
            Staccato(end+1:end+lc) = 0;
        
            %Manual xcorr
            while first <= s
            
                %Normalizes the interval to sum(energy) = 1
                DNorm(:) = Staccato(first:last) - ...
                    mean(Staccato(first:last));
                DN = sum(DNorm.^2);
                DNorm(:) = DNorm(:) ./ sqrt(DN);
            
                %Determines the correlation coefficient
                corr(:) = DNorm(:) .* c(:);
                Corr(first,1) = sum(corr(:));
            
                %Increments to next window
                first = first + 1;
                last = last + 1;
            end
        
            %Saves the energy of the note (no harmonics)
            Co(:,n) = Co(:,n) + Corr(1:s).^2;
            
            %Progresses to next harmonic
            k = k + 1;
        end
        
        %Determines the max corr value of the note+harmonics
        [Max, idx] = max(Co(:,n));
        %*******************75 Percent Threshold************************%
        
        for i = 1:s
            if n == 1
                if Co(i,n) < tolerance * Max
                    Co(i,n) = 0;
                end
            else
                if Co(i,n) < .60 * Max
                    Co(i,n) = 0;
                end
            end
        end
        
    %Empties notes not of interest in the .wav
    else 
        Co(:,n) = 0;
    end
    n = n + 1;
end
 
i = 1;
 
while i <= s
    n = 1;
    while n <= numel(notes)
        while max(Co(i,:)) == 0 && i < s
            i = i + 1;
        end
        if Co(i,n) ~= max(Co(i,:))
            Co(i,n) = 0;
        else
            Co(i,n) = Notes(n,2);
        end
        n = n + 1;
    end
    Freq(i,1) = max(Co(i,:));
    i = i + 1;
end
 
%*****************************Clean Gaps********************************%
 
%Variables
i = 1;
j = 1;
count = 1;
 
%Cleans smallest gaps, observed to be 150 samples 
while i <= s
    
    %Skips big intervals of 0
    if Freq(i, 1) == 0
        i = i + 1; 
    else
        if i ~= s
            
            %If the next interval drops
            if Freq(i+1, 1) ~= Freq(i,1)
                z = 1;
                
                %Checks within 75 samples for a match
                while z < 18500 
                    if (Freq(i+z, 1) == Freq(i,1))
                        Freq(i:i+z,1) = Freq(i, 1);
                        z = 18501;
                    else 
                        z = z + 1;
                    end
                end
                i = i + 1;
                
            %Moves on
            else 
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end  
end
 
Freq(end+1,1) = 1;
[s,~] = size(Freq);
 
%Creates an matrix of intervals
while j <= s
    
   while Freq(j,1) == 0 && j<s
       j = j + 1;
   end
   
   if j < s-1
        %Saves the beginning of the interval 
        start = j;
        no = Freq(j,1);
        
        %Moves through the interval
        while Freq(j,1) == Freq(j+1, 1)
            j = j + 1;
        end
    
        stop = j;
    
        if j < s && (stop-start >= 7500)
        
            %Saves the note's interval
            N(:,count) = [no start stop];
            count = count + 1; 
        end
   end
    j = j + 1;
end
 
%Variables
[~, End] = size(N);
count = 1;
i = 1;
int = 128;
int2 = int/2;
int4 = int/4;
 
%**********Saves the parsed notes in dynamically named files************%
while i < End && N(1,i) ~= 0
    
    %Skips empty intervals
    while N(1, i) == 0 && i ~= End
        i = i + 1;
    end
    
    %Beginning of interval
    start = N(2, i);
    stop = N(3, i);
    
    %Creates wav and jpeg file names
    basename = round(N(1, i));
    filename = ['6_' num2str(count) '.wav'];
    imagename = ['6_' num2str(count) '.jpeg'];
    
    %Saves audio snippit
    audiowrite(filename, Staccato(start:stop), fs);
    save filename
    
    %Hanning Window
    w = hanning(int);
    
    %Variables
    first = 1;
    last = int;
    
    %Determines dimensions of the image
    length = stop-start+1;
    numWindows = ceil((length-int) / int4);
    newEnd = (numWindows * int4) + int;
    imageInt = zeros(newEnd, 1);
    imageInt(1:stop-start+1) = Staccato(start:stop);
    image = zeros(int2, numWindows);
    
    window = 1;
    %Creates spectrogram from the given note
    while window <= numWindows
        interval = imageInt(first:last);
        x = (interval .* w);
        x = abs(fft(x));
        z = x(1:int2);
        
        image(:, window) = z';
        
        window = window + 1;
        first = first + int4;
        last = last + int4;
    end
    
    %Scale image to 200x200 and SAVE
    m = max(image);
    image = image ./ m;
    Image = imresize(image, [200 200]);
    Image = flipud(Image);
    imwrite(Image, imagename);
    save imagename;
    
    count = count + 1;
    i = i + 1;
end
 
 
%End
done = 'staccato done'
 
%% 7 Staccato Hooked - Energy Approach
StaccatoHooked = audioread('Dmajor_staccato_hooked.wav');
 
%Variables
fs = 44100;
dt = 1/fs;
[s,~] = size(StaccatoHooked);
secs =  s / fs;
 
%Tolerance Level
tolerance = 0.70;
 
%Notes of interest
notes = [294 330 370 392 440 494 554 587];
 
%***********************Hanning entire Audiofile************************%
 
%Hanning and fft of the entire window
h = hanning(s);
DDh = (StaccatoHooked .* h);
DD = abs(fft(DDh));
 
%Amp vs. Freq Graph
plot(DD(1:s/2)); 
 
%Table that has interest values and precise frequency
Notes = zeros(numel(notes), 2);
 
%Finds the peaks and saves their frequency into Notes[]
for i = 1:numel(notes)
    first = round(notes(i)*secs-100);
    last = round(notes(i)*secs+100);
    [val, idx] = max(DD(first:last));
    y1 = find(DD==val);
    Notes(i, 1) = val;
    Notes(i, 2) = y1(1)/secs;
end
 
%Allocates memory for big static arrays
Co = zeros(s, numel(notes));
Corr = zeros(s+2000, 1);
Freq = zeros(s+100, 1);
N = zeros(3, 75);
 
%Sets note of interest
n = 1;
 
%*******************Cross Correlation of each note**********************%
 
%Fills cross correlation vector
while n <= numel(notes)
    
    %Will only look at frequencies of interest
    if Notes(n, 1) >= 1
        
        %Beginning Harmonics
        k = 1;
        
        while k <= 6 && (k*Notes(n,2) <= fs/2)
     
            %Creates a cosine wave matching the frequency
            T = 1/Notes(n, 2);
            intv = 0:dt:T+dt;
            t = 0:dt:10/round(Notes(n, 2));
            c = (cos(2*pi*k*round(Notes(n, 2))*t))';
        
            %Normalizes the cosine wave to mean = 0 and sum(energy) = 1
            C = sum(c(:).^2);
            c(:) = c(:) ./ sqrt(C);
        
            %Variables
            [lc, ~] = size(c);
            first = 1;
            last = lc;
            corr = zeros(lc, 1);
            DNorm = zeros(lc, 1);
            StaccatoHooked(end+1:end+lc) = 0;
        
            %Manual xcorr
            while first <= s
            
                %Normalizes the interval to sum(energy) = 1
                DNorm(:) = StaccatoHooked(first:last) - ...
                    mean(StaccatoHooked(first:last));
                DN = sum(DNorm.^2);
                DNorm(:) = DNorm(:) ./ sqrt(DN);
            
                %Determines the correlation coefficient
                corr(:) = DNorm(:) .* c(:);
                Corr(first,1) = sum(corr(:));
            
                %Increments to next window
                first = first + 1;
                last = last + 1;
            end
        
            %Saves the energy of the note (no harmonics)
            Co(:,n) = Co(:,n) + Corr(1:s).^2;
            
            %Progresses to next harmonic
            k = k + 1;
        end
        
        %Determines the max corr value of the note+harmonics
        [Max, idx] = max(Co(:,n));
        %*******************75 Percent Threshold************************%
        
        for i = 1:s
            if Co(i,n) < tolerance * Max
               Co(i,n) = 0;
            end
        end
        
    %Empties notes not of interest in the .wav
    else 
        Co(:,n) = 0;
    end
    n = n + 1;
end
 
i = 1;
 
while i <= s
    n = 1;
    while n <= numel(notes)
        while max(Co(i,:)) == 0 && i < s
            i = i + 1;
        end
        if Co(i,n) ~= max(Co(i,:))
            Co(i,n) = 0;
        else
            Co(i,n) = Notes(n,2);
        end
        n = n + 1;
    end
    Freq(i,1) = max(Co(i,:));
    i = i + 1;
end
 
%*****************************Clean Gaps********************************%
 
%Variables
i = 1;
j = 1;
count = 1;
 
%Cleans smallest gaps, observed to be 150 samples 
while i <= s
    
    %Skips big intervals of 0
    if Freq(i, 1) == 0
        i = i + 1; 
    else
        if i ~= s
            
            %If the next interval drops
            if Freq(i+1, 1) ~= Freq(i,1)
                z = 1;
                
                %Checks within 75 samples for a match
                while z < 18500 
                    if (Freq(i+z, 1) == Freq(i,1))
                        Freq(i:i+z,1) = Freq(i, 1);
                        z = 18501;
                    else 
                        z = z + 1;
                    end
                end
                i = i + 1;
                
            %Moves on
            else 
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end  
end
 
Freq(end+1,1) = 1;
[s,~] = size(Freq);
 
%Creates an matrix of intervals
while j <= s
    
   while Freq(j,1) == 0 && j<s
       j = j + 1;
   end
   
   if j < s-1
        %Saves the beginning of the interval 
        start = j;
        no = Freq(j,1);
        
        %Moves through the interval
        while Freq(j,1) == Freq(j+1, 1)
            j = j + 1;
        end
    
        stop = j;
    
        if j < s && (stop-start >= 7500)
        
            %Saves the note's interval
            N(:,count) = [no start stop];
            count = count + 1; 
        end
   end
    j = j + 1;
end
 
%Variables
[~, End] = size(N);
count = 1;
i = 1;
int = 128;
int2 = int/2;
int4 = int/4;
 
%**********Saves the parsed notes in dynamically named files************%
while i < End && N(1,i) ~= 0
    
    %Skips empty intervals
    while N(1, i) == 0 && i ~= End
        i = i + 1;
    end
    
    %Beginning of interval
    start = N(2, i);
    stop = N(3, i);
    
    %Creates wav and jpeg file names
    basename = round(N(1, i));
    filename = ['7_' num2str(count) '.wav'];
    imagename = ['7_' num2str(count) '.jpeg'];
    
    %Saves audio snippit
    audiowrite(filename, StaccatoHooked(start:stop), fs);
    save filename
    
    %Hanning Window
    w = hanning(int);
    
    %Variables
    first = 1;
    last = int;
    
    %Determines dimensions of the image
    length = stop-start+1;
    numWindows = ceil((length-int) / int4);
    newEnd = (numWindows * int4) + int;
    imageInt = zeros(newEnd, 1);
    imageInt(1:stop-start+1) = StaccatoHooked(start:stop);
    image = zeros(int2, numWindows);
    
    window = 1;
    %Creates spectrogram from the given note
    while window <= numWindows
        interval = imageInt(first:last);
        x = (interval .* w);
        x = abs(fft(x));
        z = x(1:int2);
        
        image(:, window) = z';
        
        window = window + 1;
        first = first + int4;
        last = last + int4;
    end
    
    %Scale image to 200x200 and SAVE
    m = max(image);
    image = image ./ m;
    Image = imresize(image, [200 200]);
    Image = flipud(Image);
    imwrite(Image, imagename);
    save imagename;
    
    count = count + 1;
    i = i + 1;
end
 
 
%End
done = 'staccato hooked done'
 
%% 8 Portato Hooked - Energy Approach
PortatoHooked = audioread('Dmajor_portato_hooked.wav');
 
%Variables
fs = 44100;
dt = 1/fs;
[s,~] = size(PortatoHooked);
secs =  s / fs;
 
%Tolerance Level
tolerance = 0.75;
 
%Notes of interest
notes = [294 330 370 392 440 494 554 587];
 
%***********************Hanning entire Audiofile************************%
 
%Hanning and fft of the entire window
h = hanning(s);
DDh = (PortatoHooked .* h);
DD = abs(fft(DDh));
 
%Amp vs. Freq Graph
plot(DD(1:s/2)); 
 
%Table that has interest values and precise frequency
Notes = zeros(numel(notes), 2);
 
%Finds the peaks and saves their frequency into Notes[]
for i = 1:numel(notes)
    first = round(notes(i)*secs-100);
    last = round(notes(i)*secs+100);
    [val, idx] = max(DD(first:last));
    y1 = find(DD==val);
    Notes(i, 1) = val;
    Notes(i, 2) = y1(1)/secs;
end
 
%Allocates memory for big static arrays
Co = zeros(s, numel(notes));
Corr = zeros(s+2000, 1);
Freq = zeros(s+100, 1);
N = zeros(3, 75);
 
%Sets note of interest
n = 1;
 
%*******************Cross Correlation of each note**********************%
 
%Fills cross correlation vector
while n <= numel(notes)
    
    %Will only look at frequencies of interest
    if Notes(n, 1) >= 1
        
        %Beginning Harmonics
        k = 1;
        
        while k <= 6 && (k*Notes(n,2) <= fs/2)
     
            %Creates a cosine wave matching the frequency
            T = 1/Notes(n, 2);
            intv = 0:dt:T+dt;
            t = 0:dt:10/round(Notes(n, 2));
            c = (cos(2*pi*k*round(Notes(n, 2))*t))';
        
            %Normalizes the cosine wave to mean = 0 and sum(energy) = 1
            C = sum(c(:).^2);
            c(:) = c(:) ./ sqrt(C);
        
            %Variables
            [lc, ~] = size(c);
            first = 1;
            last = lc;
            corr = zeros(lc, 1);
            DNorm = zeros(lc, 1);
            PortatoHooked(end+1:end+lc) = 0;
        
            %Manual xcorr
            while first <= s
            
                %Normalizes the interval to sum(energy) = 1
                DNorm(:) = PortatoHooked(first:last) - ...
                    mean(PortatoHooked(first:last));
                DN = sum(DNorm.^2);
                DNorm(:) = DNorm(:) ./ sqrt(DN);
            
                %Determines the correlation coefficient
                corr(:) = DNorm(:) .* c(:);
                Corr(first,1) = sum(corr(:));
            
                %Increments to next window
                first = first + 1;
                last = last + 1;
            end
        
            %Saves the energy of the note (no harmonics)
            Co(:,n) = Co(:,n) + Corr(1:s).^2;
            
            %Progresses to next harmonic
            k = k + 1;
        end
        
        %Determines the max corr value of the note+harmonics
        [Max, idx] = max(Co(:,n));
        %*******************75 Percent Threshold************************%
        
        for i = 1:s
            if n == 1
                if Co(i,n) < tolerance * Max
                    Co(i,n) = 0;
                end
            else
                if Co(i,n) < .65 * Max
                    Co(i,n) = 0;
                end
            end
        end
        
    %Empties notes not of interest in the .wav
    else 
        Co(:,n) = 0;
    end
    n = n + 1;
end
 
i = 1;
 
while i <= s
    n = 1;
    while n <= numel(notes)
        while max(Co(i,:)) == 0 && i < s
            i = i + 1;
        end
        if Co(i,n) ~= max(Co(i,:))
            Co(i,n) = 0;
        else
            Co(i,n) = Notes(n,2);
        end
        n = n + 1;
    end
    Freq(i,1) = max(Co(i,:));
    i = i + 1;
end
 
%*****************************Clean Gaps********************************%
 
%Variables
i = 1;
j = 1;
count = 1;
 
%Cleans smallest gaps, observed to be 150 samples 
while i <= s
    
    %Skips big intervals of 0
    if Freq(i, 1) == 0
        i = i + 1; 
    else
        if i ~= s
            
            %If the next interval drops
            if Freq(i+1, 1) ~= Freq(i,1)
                z = 1;
                
                %Checks within 75 samples for a match
                while z < 18500 
                    if (Freq(i+z, 1) == Freq(i,1))
                        Freq(i:i+z,1) = Freq(i, 1);
                        z = 18501;
                    else 
                        z = z + 1;
                    end
                end
                i = i + 1;
                
            %Moves on
            else 
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end  
end
 
Freq(end+1,1) = 1;
[s,~] = size(Freq);
 
%Creates an matrix of intervals
while j <= s
    
   while Freq(j,1) == 0 && j<s
       j = j + 1;
   end
   
   if j < s-1
        %Saves the beginning of the interval 
        start = j;
        no = Freq(j,1);
        
        %Moves through the interval
        while Freq(j,1) == Freq(j+1, 1)
            j = j + 1;
        end
    
        stop = j;
    
        if j < s && (stop-start >= 7500)
        
            %Saves the note's interval
            N(:,count) = [no start stop];
            count = count + 1; 
        end
   end
    j = j + 1;
end
 
%Variables
[~, End] = size(N);
count = 1;
i = 1;
int = 128;
int2 = int/2;
int4 = int/4;
 
%**********Saves the parsed notes in dynamically named files************%
while i < End && N(1,i) ~= 0
    
    %Skips empty intervals
    while N(1, i) == 0 && i ~= End
        i = i + 1;
    end
    
    %Beginning of interval
    start = N(2, i);
    stop = N(3, i);
    
    %Creates wav and jpeg file names
    basename = round(N(1, i));
    filename = ['8_' num2str(count) '.wav'];
    imagename = ['8_' num2str(count) '.jpeg'];
    
    %Saves audio snippit
    audiowrite(filename, PortatoHooked(start:stop), fs);
    save filename
    
    %Hanning Window
    w = hanning(int);
    
    %Variables
    first = 1;
    last = int;
    
    %Determines dimensions of the image
    length = stop-start+1;
    numWindows = ceil((length-int) / int4);
    newEnd = (numWindows * int4) + int;
    imageInt = zeros(newEnd, 1);
    imageInt(1:stop-start+1) = PortatoHooked(start:stop);
    image = zeros(int2, numWindows);
    
    window = 1;
    %Creates spectrogram from the given note
    while window <= numWindows
        interval = imageInt(first:last);
        x = (interval .* w);
        x = abs(fft(x));
        z = x(1:int2);
        
        image(:, window) = z';
        
        window = window + 1;
        first = first + int4;
        last = last + int4;
    end
    
    %Scale image to 200x200 and SAVE
    m = max(image);
    image = image ./ m;
    Image = imresize(image, [200 200]);
    Image = flipud(Image);
    imwrite(Image, imagename);
    save imagename;
    
    count = count + 1;
    i = i + 1;
end
 
 
%End
done = 'portato hooked done'
 
%% 9 Portato Separate - Energy Approach
PortatoSeparate = audioread('Dmajor_portato_separate.wav');
 
%Variables
fs = 44100;
dt = 1/fs;
[s,~] = size(PortatoSeparate);
secs =  s / fs;
 
%Tolerance Level
tolerance = 0.60;
 
%Notes of interest
notes = [294 330 370 392 440 494 554 587];
 
%***********************Hanning entire Audiofile************************%
 
%Hanning and fft of the entire window
h = hanning(s);
DDh = (PortatoSeparate .* h);
DD = abs(fft(DDh));
 
%Amp vs. Freq Graph
plot(DD(1:s/2)); 
 
%Table that has interest values and precise frequency
Notes = zeros(numel(notes), 2);
 
%Finds the peaks and saves their frequency into Notes[]
for i = 1:numel(notes)
    first = round(notes(i)*secs-100);
    last = round(notes(i)*secs+100);
    [val, idx] = max(DD(first:last));
    y1 = find(DD==val);
    Notes(i, 1) = val;
    Notes(i, 2) = y1(1)/secs;
end
 
%Allocates memory for big static arrays
Co = zeros(s, numel(notes));
Corr = zeros(s+2000, 1);
Freq = zeros(s+100, 1);
N = zeros(3, 75);
 
%Sets note of interest
n = 1;
 
%*******************Cross Correlation of each note**********************%
 
%Fills cross correlation vector
while n <= numel(notes)
    
    %Will only look at frequencies of interest
    if Notes(n, 1) >= 1
        
        %Beginning Harmonics
        k = 1;
        
        while k <= 6 && (k*Notes(n,2) <= fs/2)
     
            %Creates a cosine wave matching the frequency
            T = 1/Notes(n, 2);
            intv = 0:dt:T+dt;
            t = 0:dt:10/round(Notes(n, 2));
            c = (cos(2*pi*k*round(Notes(n, 2))*t))';
        
            %Normalizes the cosine wave to mean = 0 and sum(energy) = 1
            C = sum(c(:).^2);
            c(:) = c(:) ./ sqrt(C);
        
            %Variables
            [lc, ~] = size(c);
            first = 1;
            last = lc;
            corr = zeros(lc, 1);
            DNorm = zeros(lc, 1);
            PortatoSeparate(end+1:end+lc) = 0;
        
            %Manual xcorr
            while first <= s
            
                %Normalizes the interval to sum(energy) = 1
                DNorm(:) = PortatoSeparate(first:last) - ...
                    mean(PortatoSeparate(first:last));
                DN = sum(DNorm.^2);
                DNorm(:) = DNorm(:) ./ sqrt(DN);
            
                %Determines the correlation coefficient
                corr(:) = DNorm(:) .* c(:);
                Corr(first,1) = sum(corr(:));
            
                %Increments to next window
                first = first + 1;
                last = last + 1;
            end
        
            %Saves the energy of the note (no harmonics)
            Co(:,n) = Co(:,n) + Corr(1:s).^2;
            
            %Progresses to next harmonic
            k = k + 1;
        end
        
        %Determines the max corr value of the note+harmonics
        [Max, idx] = max(Co(:,n));
        %*******************75 Percent Threshold************************%
        
        for i = 1:s
            if Co(i,n) < tolerance * Max
               Co(i,n) = 0;
            end
        end
        
    %Empties notes not of interest in the .wav
    else 
        Co(:,n) = 0;
    end
    n = n + 1;
end
 
i = 1;
 
while i <= s
    n = 1;
    while n <= numel(notes)
        while max(Co(i,:)) == 0 && i < s
            i = i + 1;
        end
        if Co(i,n) ~= max(Co(i,:))
            Co(i,n) = 0;
        else
            Co(i,n) = Notes(n,2);
        end
        n = n + 1;
    end
    Freq(i,1) = max(Co(i,:));
    i = i + 1;
end
 
%*****************************Clean Gaps********************************%
 
%Variables
i = 1;
j = 1;
count = 1;
 
%Cleans smallest gaps, observed to be 150 samples 
while i <= s
    
    %Skips big intervals of 0
    if Freq(i, 1) == 0
        i = i + 1; 
    else
        if i ~= s
            
            %If the next interval drops
            if Freq(i+1, 1) ~= Freq(i,1)
                z = 1;
                
                %Checks within 75 samples for a match
                while z < 18500 
                    if (Freq(i+z, 1) == Freq(i,1))
                        Freq(i:i+z,1) = Freq(i, 1);
                        z = 18501;
                    else 
                        z = z + 1;
                    end
                end
                i = i + 1;
                
            %Moves on
            else 
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end  
end
 
Freq(end+1,1) = 1;
[s,~] = size(Freq);
 
%Creates an matrix of intervals
while j <= s
    
   while Freq(j,1) == 0 && j<s
       j = j + 1;
   end
   
   if j < s-1
        %Saves the beginning of the interval 
        start = j;
        no = Freq(j,1);
        
        %Moves through the interval
        while Freq(j,1) == Freq(j+1, 1)
            j = j + 1;
        end
    
        stop = j;
    
        if j < s && (stop-start >= 7500)
        
            %Saves the note's interval
            N(:,count) = [no start stop];
            count = count + 1; 
        end
   end
    j = j + 1;
end
 
%Variables
[~, End] = size(N);
count = 1;
i = 1;
int = 128;
int2 = int/2;
int4 = int/4;
 
%**********Saves the parsed notes in dynamically named files************%
while i < End && N(1,i) ~= 0
    
    %Skips empty intervals
    while N(1, i) == 0 && i ~= End
        i = i + 1;
    end
    
    %Beginning of interval
    start = N(2, i);
    stop = N(3, i);
    
    %Creates wav and jpeg file names
    basename = round(N(1, i));
    filename = ['9_' num2str(count) '.wav'];
    imagename = ['9_' num2str(count) '.jpeg'];
    
    %Saves audio snippit
    audiowrite(filename, PortatoSeparate(start:stop), fs);
    save filename
    
    %Hanning Window
    w = hanning(int);
    
    %Variables
    first = 1;
    last = int;
    
    %Determines dimensions of the image
    length = stop-start+1;
    numWindows = ceil((length-int) / int4);
    newEnd = (numWindows * int4) + int;
    imageInt = zeros(newEnd, 1);
    imageInt(1:stop-start+1) = PortatoSeparate(start:stop);
    image = zeros(int2, numWindows);
    
    window = 1;
    %Creates spectrogram from the given note
    while window <= numWindows
        interval = imageInt(first:last);
        x = (interval .* w);
        x = abs(fft(x));
        z = x(1:int2);
        
        image(:, window) = z';
        
        window = window + 1;
        first = first + int4;
        last = last + int4;
    end
    
    %Scale image to 200x200 and SAVE
    m = max(image);
    image = image ./ m;
    Image = imresize(image, [200 200]);
    Image = flipud(Image);
    imwrite(Image, imagename);
    save imagename;
    
    count = count + 1;
    i = i + 1;
end
 
 
%End
done = 'portato separate done'
 
%% 10 Spiccato - Energy Approach
Spiccato = audioread('Dmajor_spiccato.wav');
 
%Variables
fs = 44100;
dt = 1/fs;
[s,~] = size(Spiccato);
secs =  s / fs;
 
%Tolerance Level
tolerance = 0.85 ;
 
%Notes of interest
notes = [294 330 370 392 440 494 554 587];
 
%***********************Hanning entire Audiofile************************%
 
%Hanning and fft of the entire window
h = hanning(s);
DDh = (Spiccato .* h);
DD = abs(fft(DDh));
 
%Amp vs. Freq Graph
plot(DD(1:s/2)); 
 
%Table that has interest values and precise frequency
Notes = zeros(numel(notes), 2);
 
%Finds the peaks and saves their frequency into Notes[]
for i = 1:numel(notes)
    first = round(notes(i)*secs-100);
    last = round(notes(i)*secs+100);
    [val, idx] = max(DD(first:last));
    y1 = find(DD==val);
    Notes(i, 1) = val;
    Notes(i, 2) = y1(1)/secs;
end
 
%Allocates memory for big static arrays
Co = zeros(s, numel(notes));
Corr = zeros(s+2000, 1);
Freq = zeros(s+100, 1);
N = zeros(3, 75);
 
%Sets note of interest
n = 1;
 
%*******************Cross Correlation of each note**********************%
 
%Fills cross correlation vector
while n <= numel(notes)
    
    %Will only look at frequencies of interest
    if Notes(n, 1) >= 1
        
        %Beginning Harmonics
        k = 1;
        
        while k <= 6 && (k*Notes(n,2) <= fs/2)
     
            %Creates a cosine wave matching the frequency
            T = 1/Notes(n, 2);
            intv = 0:dt:T+dt;
            t = 0:dt:10/round(Notes(n, 2));
            c = (cos(2*pi*k*round(Notes(n, 2))*t))';
        
            %Normalizes the cosine wave to mean = 0 and sum(energy) = 1
            C = sum(c(:).^2);
            c(:) = c(:) ./ sqrt(C);
        
            %Variables
            [lc, ~] = size(c);
            first = 1;
            last = lc;
            corr = zeros(lc, 1);
            DNorm = zeros(lc, 1);
            Spiccato(end+1:end+lc) = 0;
        
            %Manual xcorr
            while first <= s
            
                %Normalizes the interval to sum(energy) = 1
                DNorm(:) = Spiccato(first:last) - ...
                    mean(Spiccato(first:last));
                DN = sum(DNorm.^2);
                DNorm(:) = DNorm(:) ./ sqrt(DN);
            
                %Determines the correlation coefficient
                corr(:) = DNorm(:) .* c(:);
                Corr(first,1) = sum(corr(:));
            
                %Increments to next window
                first = first + 1;
                last = last + 1;
            end
        
            %Saves the energy of the note (no harmonics)
            Co(:,n) = Co(:,n) + Corr(1:s).^2;
            
            %Progresses to next harmonic
            k = k + 1;
        end
        
        %Determines the max corr value of the note+harmonics
        [Max, idx] = max(Co(:,n));
        %*******************75 Percent Threshold************************%
        
        for i = 1:s
            if n == 1
                if Co(i,n) < tolerance * Max
                    Co(i,n) = 0;
                end
            else
                if Co(i,n) < .60 * Max
                    Co(i,n) = 0;
                end
            end
        end
        
    %Empties notes not of interest in the .wav
    else 
        Co(:,n) = 0;
    end
    n = n + 1;
end
 
i = 1;
 
while i <= s
    n = 1;
    while n <= numel(notes)
        while max(Co(i,:)) == 0 && i < s
            i = i + 1;
        end
        if Co(i,n) ~= max(Co(i,:))
            Co(i,n) = 0;
        else
            Co(i,n) = Notes(n,2);
        end
        n = n + 1;
    end
    Freq(i,1) = max(Co(i,:));
    i = i + 1;
end
 
%*****************************Clean Gaps********************************%
 
%Variables
i = 1;
j = 1;
count = 1;
 
%Cleans smallest gaps, observed to be 150 samples 
while i <= s
    
    %Skips big intervals of 0
    if Freq(i, 1) == 0
        i = i + 1; 
    else
        if i ~= s
            
            %If the next interval drops
            if Freq(i+1, 1) ~= Freq(i,1)
                z = 1;
                
                %Checks within 75 samples for a match
                while z < 18500 
                    if (Freq(i+z, 1) == Freq(i,1))
                        Freq(i:i+z,1) = Freq(i, 1);
                        z = 18501;
                    else 
                        z = z + 1;
                    end
                end
                i = i + 1;
                
            %Moves on
            else 
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end  
end
 
Freq(end+1,1) = 1;
[s,~] = size(Freq);
 
%Creates an matrix of intervals
while j <= s
    
   while Freq(j,1) == 0 && j<s
       j = j + 1;
   end
   
   if j < s-1
        %Saves the beginning of the interval 
        start = j;
        no = Freq(j,1);
        
        %Moves through the interval
        while Freq(j,1) == Freq(j+1, 1)
            j = j + 1;
        end
    
        stop = j;
    
        if j < s && (stop-start >= 7500)
        
            %Saves the note's interval
            N(:,count) = [no start stop];
            count = count + 1; 
        end
   end
    j = j + 1;
end
 
%Variables
[~, End] = size(N);
count = 1;
i = 1;
int = 128;
int2 = int/2;
int4 = int/4;
 
%**********Saves the parsed notes in dynamically named files************%
while i < End && N(1,i) ~= 0
    
    %Skips empty intervals
    while N(1, i) == 0 && i ~= End
        i = i + 1;
    end
    
    %Beginning of interval
    start = N(2, i);
    stop = N(3, i);
    
    %Creates wav and jpeg file names
    basename = round(N(1, i));
    filename = ['10_' num2str(count) '.wav'];
    imagename = ['10_' num2str(count) '.jpeg'];
    
    %Saves audio snippit
    audiowrite(filename, Spiccato(start:stop), fs);
    save filename
    
    %Hanning Window
    w = hanning(int);
    
    %Variables
    first = 1;
    last = int;
    
    %Determines dimensions of the image
    length = stop-start+1;
    numWindows = ceil((length-int) / int4);
    newEnd = (numWindows * int4) + int;
    imageInt = zeros(newEnd, 1);
    imageInt(1:stop-start+1) = Spiccato(start:stop);
    image = zeros(int2, numWindows);
    
    window = 1;
    %Creates spectrogram from the given note
    while window <= numWindows
        interval = imageInt(first:last);
        x = (interval .* w);
        x = abs(fft(x));
        z = x(1:int2);
        
        image(:, window) = z';
        
        window = window + 1;
        first = first + int4;
        last = last + int4;
    end
    
    %Scale image to 200x200 and SAVE
    m = max(image);
    image = image ./ m;
    Image = imresize(image, [200 200]);
    Image = flipud(Image);
    imwrite(Image, imagename);
    save imagename;
    
    count = count + 1;
    i = i + 1;
end
 
 
%End
done = 'spiccato done'
 
%% 11 Spiccato Fast - Energy Approach
SpiccatoFast = audioread('Dmajor_spiccato_fast.wav');
 
%Variables
fs = 44100;
dt = 1/fs;
[s,~] = size(SpiccatoFast);
secs =  s / fs;
 
%Tolerance Level
tolerance = 0.925;
 
%Notes of interest
notes = [294 330 370 392 440 494 554 587];
 
%***********************Hanning entire Audiofile************************%
 
%Hanning and fft of the entire window
h = hanning(s);
DDh = (SpiccatoFast .* h);
DD = abs(fft(DDh));
 
%Amp vs. Freq Graph
plot(DD(1:s/2)); 
 
%Table that has interest values and precise frequency
Notes = zeros(numel(notes), 2);
 
%Finds the peaks and saves their frequency into Notes[]
for i = 1:numel(notes)
    first = round(notes(i)*secs-100);
    last = round(notes(i)*secs+100);
    [val, idx] = max(DD(first:last));
    y1 = find(DD==val);
    Notes(i, 1) = val;
    Notes(i, 2) = y1(1)/secs;
end
 
%Allocates memory for big static arrays
Co = zeros(s, numel(notes));
Corr = zeros(s+2000, 1);
Freq = zeros(s+100, 1);
N = zeros(3, 75);
 
%Sets note of interest
n = 1;
 
%*******************Cross Correlation of each note**********************%
 
%Fills cross correlation vector
while n <= numel(notes)
    
    %Will only look at frequencies of interest
    if Notes(n, 1) >= 1
        
        %Beginning Harmonics
        k = 1;
        
        while k <= 6 && (k*Notes(n,2) <= fs/2)
     
            %Creates a cosine wave matching the frequency
            T = 1/Notes(n, 2);
            intv = 0:dt:T+dt;
            t = 0:dt:10/round(Notes(n, 2));
            c = (cos(2*pi*k*round(Notes(n, 2))*t))';
        
            %Normalizes the cosine wave to mean = 0 and sum(energy) = 1
            C = sum(c(:).^2);
            c(:) = c(:) ./ sqrt(C);
        
            %Variables
            [lc, ~] = size(c);
            first = 1;
            last = lc;
            corr = zeros(lc, 1);
            DNorm = zeros(lc, 1);
            SpiccatoFast(end+1:end+lc) = 0;
        
            %Manual xcorr
            while first <= s
            
                %Normalizes the interval to sum(energy) = 1
                DNorm(:) = SpiccatoFast(first:last) - ...
                    mean(SpiccatoFast(first:last));
                DN = sum(DNorm.^2);
                DNorm(:) = DNorm(:) ./ sqrt(DN);
            
                %Determines the correlation coefficient
                corr(:) = DNorm(:) .* c(:);
                Corr(first,1) = sum(corr(:));
            
                %Increments to next window
                first = first + 1;
                last = last + 1;
            end
        
            %Saves the energy of the note (no harmonics)
            Co(:,n) = Co(:,n) + Corr(1:s).^2;
            
            %Progresses to next harmonic
            k = k + 1;
        end
        
        %Determines the max corr value of the note+harmonics
        [Max, idx] = max(Co(:,n));
        %*******************75 Percent Threshold************************%
        
        for i = 1:s
            if n == 1
                if Co(i,n) < tolerance * Max
                    Co(i,n) = 0;
                end
            else
                if Co(i,n) < .55 * Max
                    Co(i,n) = 0;
                end
            end
        end
        
    %Empties notes not of interest in the .wav
    else 
        Co(:,n) = 0;
    end
    n = n + 1;
end
 
i = 1;
 
while i <= s
    n = 1;
    while n <= numel(notes)
        while max(Co(i,:)) == 0 && i < s
            i = i + 1;
        end
        if Co(i,n) ~= max(Co(i,:))
            Co(i,n) = 0;
        else
            Co(i,n) = Notes(n,2);
        end
        n = n + 1;
    end
    Freq(i,1) = max(Co(i,:));
    i = i + 1;
end
 
%*****************************Clean Gaps********************************%
 
%Variables
i = 1;
j = 1;
count = 1;
 
%Cleans smallest gaps, observed to be 150 samples 
while i <= s
    
    %Skips big intervals of 0
    if Freq(i, 1) == 0
        i = i + 1; 
    else
        if i ~= s
            
            %If the next interval drops
            if Freq(i+1, 1) ~= Freq(i,1)
                z = 1;
                
                %Checks within 75 samples for a match
                while z < 18500 
                    if (Freq(i+z, 1) == Freq(i,1))
                        Freq(i:i+z,1) = Freq(i, 1);
                        z = 18501;
                    else 
                        z = z + 1;
                    end
                end
                i = i + 1;
                
            %Moves on
            else 
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end  
end
 
Freq(end+1,1) = 1;
[s,~] = size(Freq);
 
%Creates an matrix of intervals
while j <= s
    
   while Freq(j,1) == 0 && j<s
       j = j + 1;
   end
   
   if j < s-1
        %Saves the beginning of the interval 
        start = j;
        no = Freq(j,1);
        
        %Moves through the interval
        while Freq(j,1) == Freq(j+1, 1)
            j = j + 1;
        end
    
        stop = j;
    
        if j < s && (stop-start >= 7500)
        
            %Saves the note's interval
            N(:,count) = [no start stop];
            count = count + 1; 
        end
   end
    j = j + 1;
end
 
%Variables
[~, End] = size(N);
count = 1;
i = 1;
int = 128;
int2 = int/2;
int4 = int/4;
 
%**********Saves the parsed notes in dynamically named files************%
while i < End && N(1,i) ~= 0
    
    %Skips empty intervals
    while N(1, i) == 0 && i ~= End
        i = i + 1;
    end
    
    %Beginning of interval
    start = N(2, i);
    stop = N(3, i);
    
    %Creates wav and jpeg file names
    basename = round(N(1, i));
    filename = ['11_' num2str(count) '.wav'];
    imagename = ['11_' num2str(count) '.jpeg'];
    
    %Saves audio snippit
    audiowrite(filename, SpiccatoFast(start:stop), fs);
    save filename
    
    %Hanning Window
    w = hanning(int);
    
    %Variables
    first = 1;
    last = int;
    
    %Determines dimensions of the image
    length = stop-start+1;
    numWindows = ceil((length-int) / int4);
    newEnd = (numWindows * int4) + int;
    imageInt = zeros(newEnd, 1);
    imageInt(1:stop-start+1) = SpiccatoFast(start:stop);
    image = zeros(int2, numWindows);
    
    window = 1;
    %Creates spectrogram from the given note
    while window <= numWindows
        interval = imageInt(first:last);
        x = (interval .* w);
        x = abs(fft(x));
        z = x(1:int2);
        
        image(:, window) = z';
        
        window = window + 1;
        first = first + int4;
        last = last + int4;
    end
    
    %Scale image to 200x200 and SAVE
    m = max(image);
    image = image ./ m;
    Image = imresize(image, [200 200]);
    Image = flipud(Image);
    imwrite(Image, imagename);
    save imagename;
    
    count = count + 1;
    i = i + 1;
end
 
 
%End
done = 'spiccato fast done'
 
%% 12 Ricochet Fast - Energy Approach
RicochetFast = audioread('Dmajor_ricochet_fast.wav');
 
%Variables
fs = 44100;
dt = 1/fs;
[s,~] = size(RicochetFast);
secs =  s / fs;
 
%Tolerance Level
tolerance = 0.875;
 
%Notes of interest
notes = [294 330 370 392 440 494 554 587];
 
%***********************Hanning entire Audiofile************************%
 
%Hanning and fft of the entire window
h = hanning(s);
DDh = (RicochetFast .* h);
DD = abs(fft(DDh));
 
%Amp vs. Freq Graph
plot(DD(1:s/2)); 
 
%Table that has interest values and precise frequency
Notes = zeros(numel(notes), 2);
 
%Finds the peaks and saves their frequency into Notes[]
for i = 1:numel(notes)
    first = round(notes(i)*secs-100);
    last = round(notes(i)*secs+100);
    [val, idx] = max(DD(first:last));
    y1 = find(DD==val);
    Notes(i, 1) = val;
    Notes(i, 2) = y1(1)/secs;
end
 
%Allocates memory for big static arrays
Co = zeros(s, numel(notes));
Corr = zeros(s+2000, 1);
Freq = zeros(s+100, 1);
N = zeros(3, 75);
 
%Sets note of interest
n = 1;
 
%*******************Cross Correlation of each note**********************%
 
%Fills cross correlation vector
while n <= numel(notes)
    
    %Will only look at frequencies of interest
    if Notes(n, 1) >= 1
        
        %Beginning Harmonics
        k = 1;
        
        while k <= 6 && (k*Notes(n,2) <= fs/2)
     
            %Creates a cosine wave matching the frequency
            T = 1/Notes(n, 2);
            intv = 0:dt:T+dt;
            t = 0:dt:10/round(Notes(n, 2));
            c = (cos(2*pi*k*round(Notes(n, 2))*t))';
        
            %Normalizes the cosine wave to mean = 0 and sum(energy) = 1
            C = sum(c(:).^2);
            c(:) = c(:) ./ sqrt(C);
        
            %Variables
            [lc, ~] = size(c);
            first = 1;
            last = lc;
            corr = zeros(lc, 1);
            DNorm = zeros(lc, 1);
            RicochetFast(end+1:end+lc) = 0;
        
            %Manual xcorr
            while first <= s
            
                %Normalizes the interval to sum(energy) = 1
                DNorm(:) = RicochetFast(first:last) - ...
                    mean(RicochetFast(first:last));
                DN = sum(DNorm.^2);
                DNorm(:) = DNorm(:) ./ sqrt(DN);
            
                %Determines the correlation coefficient
                corr(:) = DNorm(:) .* c(:);
                Corr(first,1) = sum(corr(:));
            
                %Increments to next window
                first = first + 1;
                last = last + 1;
            end
        
            %Saves the energy of the note (no harmonics)
            Co(:,n) = Co(:,n) + Corr(1:s).^2;
            
            %Progresses to next harmonic
            k = k + 1;
        end
        
        %Determines the max corr value of the note+harmonics
        [Max, idx] = max(Co(:,n));
        %*******************75 Percent Threshold************************%
        
        for i = 1:s
            if n == 1
                if Co(i,n) < tolerance * Max
                    Co(i,n) = 0;
                end
            else
                if Co(i,n) < .55 * Max
                    Co(i,n) = 0;
                end
            end
        end
        
    %Empties notes not of interest in the .wav
    else 
        Co(:,n) = 0;
    end
    n = n + 1;
end
 
i = 1;
 
while i <= s
    n = 1;
    while n <= numel(notes)
        while max(Co(i,:)) == 0 && i < s
            i = i + 1;
        end
        if Co(i,n) ~= max(Co(i,:))
            Co(i,n) = 0;
        else
            Co(i,n) = Notes(n,2);
        end
        n = n + 1;
    end
    Freq(i,1) = max(Co(i,:));
    i = i + 1;
end
 
%*****************************Clean Gaps********************************%
 
%Variables
i = 1;
j = 1;
count = 1;
 
%Cleans smallest gaps, observed to be 150 samples 
while i <= s
    
    %Skips big intervals of 0
    if Freq(i, 1) == 0
        i = i + 1; 
    else
        if i ~= s
            
            %If the next interval drops
            if Freq(i+1, 1) ~= Freq(i,1)
                z = 1;
                
                %Checks within 75 samples for a match
                while z < 18500 
                    if (Freq(i+z, 1) == Freq(i,1))
                        Freq(i:i+z,1) = Freq(i, 1);
                        z = 18501;
                    else 
                        z = z + 1;
                    end
                end
                i = i + 1;
                
            %Moves on
            else 
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end  
end
 
Freq(end+1,1) = 1;
[s,~] = size(Freq);
 
%Creates an matrix of intervals
while j <= s
    
   while Freq(j,1) == 0 && j<s
       j = j + 1;
   end
   
   if j < s-1
        %Saves the beginning of the interval 
        start = j;
        no = Freq(j,1);
        
        %Moves through the interval
        while Freq(j,1) == Freq(j+1, 1)
            j = j + 1;
        end
    
        stop = j;
    
        if j < s && (stop-start >= 7500)
        
            %Saves the note's interval
            N(:,count) = [no start stop];
            count = count + 1; 
        end
   end
    j = j + 1;
end
 
%Variables
[~, End] = size(N);
count = 1;
i = 1;
int = 128;
int2 = int/2;
int4 = int/4;
 
%**********Saves the parsed notes in dynamically named files************%
while i < End && N(1,i) ~= 0
    
    %Skips empty intervals
    while N(1, i) == 0 && i ~= End
        i = i + 1;
    end
    
    %Beginning of interval
    start = N(2, i);
    stop = N(3, i);
    
    %Creates wav and jpeg file names
    basename = round(N(1, i));
    filename = ['12_' num2str(count) '.wav'];
    imagename = ['12_' num2str(count) '.jpeg'];
    
    %Saves audio snippit
    audiowrite(filename, RicochetFast(start:stop), fs);
    save filename
    
    %Hanning Window
    w = hanning(int);
    
    %Variables
    first = 1;
    last = int;
    
    %Determines dimensions of the image
    length = stop-start+1;
    numWindows = ceil((length-int) / int4);
    newEnd = (numWindows * int4) + int;
    imageInt = zeros(newEnd, 1);
    imageInt(1:stop-start+1) = RicochetFast(start:stop);
    image = zeros(int2, numWindows);
    
    window = 1;
    %Creates spectrogram from the given note
    while window <= numWindows
        interval = imageInt(first:last);
        x = (interval .* w);
        x = abs(fft(x));
        z = x(1:int2);
        
        image(:, window) = z';
        
        window = window + 1;
        first = first + int4;
        last = last + int4;
    end
    
    %Scale image to 200x200 and SAVE
    m = max(image);
    image = image ./ m;
    Image = imresize(image, [200 200]);
    Image = flipud(Image);
    imwrite(Image, imagename);
    save imagename;
    
    count = count + 1;
    i = i + 1;
end
 
 
%End
done = 'ricochet fast done'
 
%% 13 Ricochet - Energy Approach
Ricochet = audioread('Dmajor_ricochet.wav');
 
%Variables
fs = 44100;
dt = 1/fs;
[s,~] = size(Ricochet);
secs =  s / fs;
 
%Tolerance Level
tolerance = 0.825;
 
%Notes of interest
notes = [294 330 370 392 440 494 554 587];
 
%***********************Hanning entire Audiofile************************%
 
%Hanning and fft of the entire window
h = hanning(s);
DDh = (Ricochet .* h);
DD = abs(fft(DDh));
 
%Amp vs. Freq Graph
plot(DD(1:s/2)); 
 
%Table that has interest values and precise frequency
Notes = zeros(numel(notes), 2);
 
%Finds the peaks and saves their frequency into Notes[]
for i = 1:numel(notes)
    first = round(notes(i)*secs-100);
    last = round(notes(i)*secs+100);
    [val, idx] = max(DD(first:last));
    y1 = find(DD==val);
    Notes(i, 1) = val;
    Notes(i, 2) = y1(1)/secs;
end
 
%Allocates memory for big static arrays
Co = zeros(s, numel(notes));
Corr = zeros(s+2000, 1);
Freq = zeros(s+100, 1);
N = zeros(3, 75);
 
%Sets note of interest
n = 1;
 
%*******************Cross Correlation of each note**********************%
 
%Fills cross correlation vector
while n <= numel(notes)
    
    %Will only look at frequencies of interest
    if Notes(n, 1) >= 1
        
        %Beginning Harmonics
        k = 1;
        
        while k <= 6 && (k*Notes(n,2) <= fs/2)
     
            %Creates a cosine wave matching the frequency
            T = 1/Notes(n, 2);
            intv = 0:dt:T+dt;
            t = 0:dt:10/round(Notes(n, 2));
            c = (cos(2*pi*k*round(Notes(n, 2))*t))';
        
            %Normalizes the cosine wave to mean = 0 and sum(energy) = 1
            C = sum(c(:).^2);
            c(:) = c(:) ./ sqrt(C);
        
            %Variables
            [lc, ~] = size(c);
            first = 1;
            last = lc;
            corr = zeros(lc, 1);
            DNorm = zeros(lc, 1);
            Ricochet(end+1:end+lc) = 0;
        
            %Manual xcorr
            while first <= s
            
                %Normalizes the interval to sum(energy) = 1
                DNorm(:) = Ricochet(first:last) - ...
                    mean(Ricochet(first:last));
                DN = sum(DNorm.^2);
                DNorm(:) = DNorm(:) ./ sqrt(DN);
            
                %Determines the correlation coefficient
                corr(:) = DNorm(:) .* c(:);
                Corr(first,1) = sum(corr(:));
            
                %Increments to next window
                first = first + 1;
                last = last + 1;
            end
        
            %Saves the energy of the note (no harmonics)
            Co(:,n) = Co(:,n) + Corr(1:s).^2;
            
            %Progresses to next harmonic
            k = k + 1;
        end
        
        %Determines the max corr value of the note+harmonics
        [Max, idx] = max(Co(:,n));
        %*******************75 Percent Threshold************************%
        
        for i = 1:s
            if n == 1
                if Co(i,n) < tolerance * Max
                    Co(i,n) = 0;
                end
            else
                if Co(i,n) < .65 * Max
                    Co(i,n) = 0;
                end
            end
        end
        
    %Empties notes not of interest in the .wav
    else 
        Co(:,n) = 0;
    end
    n = n + 1;
end
 
i = 1;
n = 1;
 
while i <= s
    n = 1;
    while n <= numel(notes)
        while max(Co(i,:)) == 0 && i < s
            i = i + 1;
        end
        if Co(i,n) ~= max(Co(i,:))
            Co(i,n) = 0;
        else
            Co(i,n) = Notes(n,2);
        end
        n = n + 1;
    end
    Freq(i,1) = max(Co(i,:));
    i = i + 1;
end
 
%*****************************Clean Gaps********************************%
 
%Variables
i = 1;
j = 1;
count = 1;
 
%Cleans smallest gaps, observed to be 150 samples 
while i <= s
    
    %Skips big intervals of 0
    if Freq(i, 1) == 0
        i = i + 1; 
    else
        if i ~= s
            
            %If the next interval drops
            if Freq(i+1, 1) ~= Freq(i,1)
                z = 1;
                
                %Checks within 75 samples for a match
                while z < 18500 
                    if (Freq(i+z, 1) == Freq(i,1))
                        Freq(i:i+z,1) = Freq(i, 1);
                        z = 18501;
                    else 
                        z = z + 1;
                    end
                end
                i = i + 1;
                
            %Moves on
            else 
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end  
end
 
Freq(end+1,1) = 1;
[s,~] = size(Freq);
 
%Creates an matrix of intervals
while j <= s
    
   while Freq(j,1) == 0 && j<s
       j = j + 1;
   end
   
   if j < s-1
        %Saves the beginning of the interval 
        start = j;
        no = Freq(j,1);
        
        %Moves through the interval
        while Freq(j,1) == Freq(j+1, 1)
            j = j + 1;
        end
    
        stop = j;
    
        if j < s && (stop-start >= 7500)
        
            %Saves the note's interval
            N(:,count) = [no start stop];
            count = count + 1; 
        end
   end
    j = j + 1;
end
 
%Variables
[~, End] = size(N);
count = 1;
i = 1;
int = 128;
int2 = int/2;
int4 = int/4;
 
%**********Saves the parsed notes in dynamically named files************%
while i < End && N(1,i) ~= 0
    
    %Skips empty intervals
    while N(1, i) == 0 && i ~= End
        i = i + 1;
    end
    
    %Beginning of interval
    start = N(2, i);
    stop = N(3, i);
    
    %Creates wav and jpeg file names
    basename = round(N(1, i));
    filename = ['13_' num2str(count) '.wav'];
    imagename = ['13_' num2str(count) '.jpeg'];
    
    %Saves audio snippit
    audiowrite(filename, Ricochet(start:stop), fs);
    save filename
    
    %Hanning Window
    w = hanning(int);
    
    %Variables
    first = 1;
    last = int;
    
    %Determines dimensions of the image
    length = stop-start+1;
    numWindows = ceil((length-int) / int4);
    newEnd = (numWindows * int4) + int;
    imageInt = zeros(newEnd, 1);
    imageInt(1:stop-start+1) = Ricochet(start:stop);
    image = zeros(int2, numWindows);
    
    window = 1;
    %Creates spectrogram from the given note
    while window <= numWindows
        interval = imageInt(first:last);
        x = (interval .* w);
        x = abs(fft(x));
        z = x(1:int2);
        
        image(:, window) = z';
        
        window = window + 1;
        first = first + int4;
        last = last + int4;
    end
    
    %Scale image to 200x200 and SAVE
    m = max(image);
    image = image ./ m;
    Image = imresize(image, [200 200]);
    Image = flipud(Image);
    imwrite(Image, imagename);
    save imagename;
    
    count = count + 1;
    i = i + 1;
end
 
 
%End
done = 'ricochet done'
