function [a_wave, I_a, b_wave, I_b] = ab_detection(d2y, Ts)

L   = length(d2y);

% %Cancellation of b wave
% d2y_ab = d2y;
% 
% d2y(d2y<0) = 0;
for i = 1:L
    C(i) = max(0,d2y(i));
end

%Squaring
Z = (C).^2;

%Generating blocks of interest

%Calculating Moving Averages

%1st Moving Average
W1_time = 0.200; % [s]
W1      = round(W1_time/Ts); %Index

% round to nearest odd integer.
idx = mod(W1,2)<1;
W1 = floor(W1);
W1(idx) = W1(idx)+1;

% MA_peak = zeros(1,L);

% for i = 1 + (W1-1)/2 : L - (W1-1)/2
%     MA_peak(i) = 1/W1 * (sum(Z(i-(W1-1)/2:i+(W1-1)/2)));
% end

MA_peak = movmean(Z,W1);

%2nd Moving Average
W2_time = 0.667; % [s]
W2      = W2_time/Ts; %Index

% round to nearest odd integer.
idx = mod(W2,2)<1;
W2 = floor(W2);
W2(idx) = W2(idx)+1;

% MA_beat = zeros(1,L);

% for j = 1 + (W2-1)/2 : L - (W2-1)/2
%     MA_beat(j) = 1/W2 * (sum(Z(j-(W2-1)/2:j+(W2-1)/2)));
% end

MA_beat = movmean(Z,W2);

%Thresholding
avg_Z = nanmean(Z);
alpha = 0.02 * avg_Z;

THR_1 = MA_beat + alpha;

BlocksOfInterest = zeros(1,L);

for k = 1:length(MA_peak)
    if MA_peak(k) > THR_1(k)
        BlocksOfInterest(k) = 0.1;
    else
        BlocksOfInterest(k) = 0;
    end
end

Blocks = (BlocksOfInterest~=0);

THR_2 = W1;
% THR_2 = 20;

%Counting the number of blocks of interest
f_ones      = find(diff([0,Blocks,0]==1));
p           = f_ones(1:2:end-1);  % Start indices
cons_Blocks = f_ones(2:2:end) - p;  % Consecutive ones’ counts


B      = length(cons_Blocks);
% a_wave = zeros(1,B);
% I_a    = zeros(1,B);
index  = zeros(2,B);

for m = 1 : B
    
    if cons_Blocks(m) > THR_2 || cons_Blocks(m) == THR_2
        
        index_start = f_ones(m+m-1); %block start
        index_end   = f_ones(m+m) - 1; %block end
        
        [~, I_a(m)] = max(Z(index_start:index_end)); %the blocks where there is an a wave
        I_a(m)      = I_a(m) + index_start;
        
        a_wave(m) = d2y(I_a(m));
%         try
%             a_wave(m) = d2y_ab(I_a(m));
%         catch 
%             a_wave(m) = d2y_ab(I_a(m) - 1);
%         end
       
    end
    
end

%Finding the b wave

I_a(I_a==0) = [];

for m = 1:length(I_a) - 1
    if I_a(m) ~= 0 && I_a(m) ~= 1 
        [b_wave(m), I_b(m)] = min(d2y(I_a(m):I_a(m+1)));
        I_b(m) = I_b(m) + I_a(m);
    else
        I_a(m) = 1;
        [b_wave(m), I_b(m)] = min(d2y_ab);
    end
end

%Don't know the reason why the point are shifted one sample away, both PPGf
%and PPGw
I_a = I_a - 1;
I_b = I_b - 1;

end

