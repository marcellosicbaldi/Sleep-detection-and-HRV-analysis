function [PPG_filt] = PPG_preproc(PPG,fs)

% By Vandercastele et al. 2018
% Band Pass Filtering
ordr = 5; 
ft_lp = 12; % [Hz]
Wn_lp = ft_lp/(fs/2); 
ft_hp = 0.5; % [Hz]
Wn_hp = ft_hp/(fs/2); 

Wn_bp = [Wn_hp,Wn_lp];

[b,a]    = butter(ordr,Wn_bp);
PPG_filt = filtfilt(b,a,PPG);


end

