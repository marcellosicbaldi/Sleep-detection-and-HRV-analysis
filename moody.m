function SNR = moody(pulse)

% deltaT_s = 50; % [ms]
% deltaT   = deltaT_s * fs; %[samples]

% Determine Signal
maxPulse = max(pulse);
minPulse = min(pulse);

S = maxPulse - minPulse; 

% Determine Noise
avg     = mean(pulse);
avg_vec = avg * ones(size(pulse));

N     = sqrt(sum((pulse - avg_vec).^2)/length(pulse));

SNR = 10 * log10(S/N);

end

