% Use basic and diagnostic quality classifier
clear all
close all
clc

load('basic_classifier.mat');
load('diagnostic_classifier.mat'); 

% Put in parent_dir the path to the folder to be processed (for single
% subject)
parent_dir = 'C:\Users\Alice\OneDrive - Alma Mater Studiorum Università di Bologna\Perso nal_Health_Systems_Lab\Tutor\WearableSensors2021\Tesina\Rec\1';

PPGraw = importdata('Rec/7/BVP.csv');
ACCraw = importdata('Rec/7/ACC.csv');

% PPGraw = importdata('Rec/7/BVP.csv');
% ACCraw = importdata('Rec/7/ACC.csv');

% PPGraw = load(fullfile(parent_dir,'BVP.csv'));
% ACCraw = load(fullfile(parent_dir,'ACC.csv'));

template_struct = load('Template_v3.m','-mat');
template = template_struct.template;

start = datetime(PPGraw(1),'convertfrom','posixtime','timezone','Europe/Rome');

PPG    = PPGraw(3:end);
N_PPG  = length(PPG);
fs_PPG = PPGraw(2);
Ts_PPG = 1/fs_PPG;
t_PPGs = (0:1:N_PPG - 1) * Ts_PPG;
t_PPG  = start + seconds(t_PPGs);

ACC    = ACCraw(3:end,:);
N_ACC  = length(ACC);
fs_ACC = ACCraw(2,1);
Ts_ACC = 1/fs_ACC;
t_ACCs = (0:1:N_ACC - 1) * Ts_ACC;
t_ACC  = start + seconds(t_ACCs);

g = 9.8; %[m/s^2]
ACC_g(:,1) = ACC(:,1) / (64);
ACC_g(:,2) = ACC(:,2) / (64);
ACC_g(:,3) = ACC(:,3) / (64);

% Modulus from the three accelerometer components
ACC_mod = sqrt(ACC_g(:,1).^2 + ACC_g(:,2).^2 + ACC_g(:,3).^2);

PPG_filt = PPG_preproc(PPG,fs_PPG);
L        = length(PPG_filt);

[sys_peak, I_peak, sys_foot, I_foot] = ab_detection(PPG_filt,Ts_PPG);

I_foot_acc = round((I_foot .* Ts_PPG) .* fs_ACC);

for i = 1:length(I_foot) - 3
        
        pulse_analysis_nozscore = PPG_filt(I_foot(i):I_foot(i+1));
        
        pulse_analysis = zscore(PPG_filt(I_foot(i):I_foot(i+1)));
        P = length(pulse_analysis);
        
        pulse_analysis(pulse_analysis == 0) = 0.01;
        %In such this way, log(0) is replaced by log(0.01), that is not -Inf
        
        % Signal Quality
        
        % By Moody's algorithm
        SNR(i) = moody(pulse_analysis);
       
        % Accelerometer peak2peak amplitude
        pp_amp(i) = peak2peak(ACC_mod(I_foot_acc(i):I_foot_acc(i+1)));
   
        % Spectral Coherence, MSE, and PSNR
        if length(template) < length(pulse_analysis)
            N = length(pulse_analysis);
            template_2 = interp1(1:numel(template),template,linspace(1,numel(template),numel(pulse_analysis)));
            R_mat = corrcoef(template_2,pulse_analysis);
            R_corr(i) = R_mat(1,2);
            Cxy(i) = mean(mscohere(template_2,pulse_analysis,[],[],[],64));
            MSE(i) = mean((template_2 - pulse_analysis').^2);
            PSNR(i) = 10*log10(max(template_2)^2/MSE(i));
        elseif length(template) > length(pulse_analysis)
            N = length(template);
            pulse_analysis_2 = interp1(1:numel(pulse_analysis),pulse_analysis,linspace(1,numel(pulse_analysis),numel(template)));
            R_mat = corrcoef(template,pulse_analysis_2);
            R_corr(i) = R_mat(1,2);
            Cxy(i) = mean(mscohere(template,pulse_analysis_2,[],[],[],64));
            MSE(i) = mean((template - pulse_analysis_2).^2);
            PSNR(i) = 10*log10(max(template)/MSE(i));
        else
            N = length(pulse_analysis); % = length(template);
            R_mat = corrcoef(template,pulse_analysis);
            R_corr(i) = R_mat(1,2);
            Cxy(i) = mean(mscohere(template,pulse_analysis,[],[],[],64));
            MSE(i) = mean((template - pulse_analysis').^2);
            PSNR(i) = 10*log10(max(template)^2/MSE(i));
        end
        
        % By Jang18: Correlation coefficient between two adjancent pulses
        if i == 1
            Sig_sim(i) = 1;
        else
            pulse_analysis_prev = PPG_filt(I_foot(i-1):I_foot(i));
            if length(pulse_analysis_prev) < length(pulse_analysis)
                pulse_analysis_prev_2 = interp1(1:numel(pulse_analysis_prev),pulse_analysis_prev,linspace(1,numel(pulse_analysis_prev),numel(pulse_analysis)));
                Sig_sim_mat = corrcoef(pulse_analysis_prev_2,pulse_analysis);
                Sig_sim(i) = Sig_sim_mat(1,2);
            elseif length(pulse_analysis_prev) > length(pulse_analysis)
                pulse_analysis_2 = interp1(1:numel(pulse_analysis),pulse_analysis,linspace(1,numel(pulse_analysis),numel(pulse_analysis_prev)));
                Sig_sim_mat = corrcoef(pulse_analysis_prev,pulse_analysis_2);
                Sig_sim(i) = Sig_sim_mat(1,2);
            else
                Sig_sim_mat = corrcoef(pulse_analysis_prev,pulse_analysis);
                Sig_sim(i) = Sig_sim_mat(1,2);
            end
        end
        
        % Skewness
        mean_sd = mean(pulse_analysis) / std(pulse_analysis);
        vec_mean_sd = mean_sd * ones(size(pulse_analysis));
        
        S(i) = 1/P * sum((pulse_analysis - vec_mean_sd).^3)';
        
        % Kurtosis
        Kur(i) = 1/P * sum((pulse_analysis - vec_mean_sd).^4);
                
        % Relative Power
        [PSD,F] = pwelch(pulse_analysis,[],[],[],fs_PPG);
        
        f_1Hz = 1; %Hz
        f1    = 2.25; %Hz
        f2    = 8; %Hz
        
        index_f1Hz = find(F == f_1Hz);
        index_f1   = find(F == f1);
        index_f2   = find(F == f2);
        
        RelP(i) = sum(PSD(index_f1Hz:index_f1))/sum(PSD(1:index_f2));
        
        % Statistical parameters
        StdSig(i) = std(pulse_analysis_nozscore);
        
        MedianSig(i) = median(pulse_analysis);

        % By Sukor11: morphological features
       trough_depth(i) = abs(pulse_analysis(1) - pulse_analysis(end));
        
        pos_foot(i,1) = I_foot(i);
        pos_foot(i,2) = I_foot(i+1);
        
        % Updated 18.05
        % # of detected peaks
        [pks] = findpeaks(pulse_analysis);
        if isempty(pks)
            Npeaks(i) = 0; 
        else
            Npeaks(i) = length(pks); 
        end
        
        % Derivative zerocrossing rate
        der = diff(pulse_analysis);
        P_der = length(der);
        der_zerocross = length(find (der < 0));
        
        Z_der(i) = 1/P_der * sum(der_zerocross);
end
disp('end loop')
FeaturesNames_basic = cellstr({'Peak2peak ACC','Signal Similarity','Skewness','Trough depth','MedianPulse','StdPulse noZ','Npeaks','Z der rate','SigCoh','PSNR'});
Matrix_basic = [pp_amp', Sig_sim', S', trough_depth', MedianSig', StdSig', Npeaks', Z_der', Cxy', PSNR'];
Table_basic = array2table(Matrix_basic,'VariableNames',FeaturesNames_basic);

FeaturesNames_diagnostic = {'Peak2peak ACC','Signal Similarity','Kurtosis','Rel Power','Skewness','MedianPulse','StdPulse noZ','SNR Moody','Npeaks','Z der rate','MSE','PSNR'};
Matrix_diagnostic = [pp_amp', Sig_sim', Kur', RelP', S', MedianSig', StdSig', SNR', Npeaks', Z_der', MSE', PSNR'];
Table_diagnostic = array2table(Matrix_diagnostic,'VariableNames',FeaturesNames_diagnostic);

%% Preprocess features - Basic classifier
% Box-Cox transformation 
% To make the features normally-distributed
% Requirements: data > 0 
min_features = zeros(size(Matrix_basic,2),1); 
New_Matrix = zeros(size(Matrix_basic)); 

var_trans = zeros(size(Matrix_basic,2),1); 

Norm_Matrix = zeros(size(Matrix_basic)); 
lambda = zeros(size(Matrix_basic,2),1);

Matrix_min = zeros(size(Matrix_basic));

for i = 1:size(Matrix_basic,2)
    
    min_features(i) = min(Matrix_basic(:,i)); 
    
    if min_features(i) < 0 || min_features(i) == 0
        New_Matrix(:,i) = Matrix_basic(:,i) + abs(min_features(i)) + eps;
        var_trans(i) = 1; 
    else
        New_Matrix(:,i) = Matrix_basic(:,i); 
    end
    
    [Norm_Matrix(:,i), lambda(i)] = boxcox(New_Matrix(:,i)); 
    
    if min_features(i) < 0 || min_features(i) == 0 
        Matrix_min(:,i) = Norm_Matrix(:,i) - abs(min_features(i)) - eps;
    else
        Matrix_min(:,i) = Norm_Matrix(:,i);
    end
    
end

Matrix_basic_zscored = zscore(Matrix_min); 

%% Preprocess features - Diagnostic classifier
% Box-Cox transformation 
% To make the features normally-distributed
% Requirements: data > 0 
min_features = zeros(size(Matrix_diagnostic,2),1); 
New_Matrix = zeros(size(Matrix_diagnostic)); 

var_trans = zeros(size(Matrix_diagnostic,2),1); 

Norm_Matrix = zeros(size(Matrix_diagnostic)); 
lambda = zeros(size(Matrix_diagnostic,2),1);

Matrix_min = zeros(size(Matrix_diagnostic));

for i = 1:size(Matrix_diagnostic,2)
    
    min_features(i) = min(Matrix_diagnostic(:,i)); 
    
    if min_features(i) < 0 || min_features(i) == 0
        New_Matrix(:,i) = Matrix_diagnostic(:,i) + abs(min_features(i)) + eps;
        var_trans(i) = 1; 
    else
        New_Matrix(:,i) = Matrix_diagnostic(:,i); 
    end
    
    [Norm_Matrix(:,i), lambda(i)] = boxcox(New_Matrix(:,i)); 
    
    if min_features(i) < 0 || min_features(i) == 0 
        Matrix_min(:,i) = Norm_Matrix(:,i) - abs(min_features(i)) - eps;
    else
        Matrix_min(:,i) = Norm_Matrix(:,i);
    end
    
end
Matrix_diagnostic_zscored = zscore(Matrix_min); 
Table_diagnostic_zscored = array2table(Matrix_diagnostic_zscored,'VariableNames',FeaturesNames_diagnostic);

%% Basic classifier
% y_basic = 1 if the pulse can be used to estimate the HR, 0 otherwise
y_basic_double = basic_classifier(Matrix_basic_zscored'); 
y_basic = logical(round(y_basic_double)); 

% Find the foot indices of pulses that can be used to estimate the HR
ind_basic = find(y_basic == 1); 

% E.g.: Interbeat Interval 
IBI = (I_foot(ind_basic+1) - I_foot(ind_basic)) .* Ts_PPG;

%% Diagnostic classifier
% y_diagnostic = 1 if the pulse can be used for morphological analysis, 0 otherwise
y_diagnostic = diagnostic_classifier.predictFcn(Table_diagnostic_zscored); 

ind_diagnostic = find(y_diagnostic == 1); 

for i = 1:length(ind_diagnostic) - 1
    
    % E.g.: Crest Time
    CT(i) = (I_peak(ind_diagnostic(i)+1) - I_foot(ind_diagnostic(i))) * Ts_PPG;
    % + 1 in I_peak because the ab_detection algorithm detect the I_peak
    % first and then the I_foot
    
end


