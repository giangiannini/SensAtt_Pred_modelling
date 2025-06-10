
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Beta Bernoulli Toy Model %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear

% generate (n) random observations from probabilistic state (state)
state = 0.75;
n = 200; 

% initialize parameters
a = 1; 
b = 1; 
learning = 1; 
forgetting = 0.01; 

% set up other options
plot_sequence = 1; 
plot_indices = 1; 

% start the code
observations = rand(1,n);
observations(observations <= state) = 1;
observations(observations ~= 1) = 0; 
theta = linspace(0,1,100); 
binWidth = diff(theta);
binWidth = [binWidth(end),binWidth];
for i = 1:length(observations)
    %priors
    disp([a b])
    old_a = a;
    old_b = b;
    Y_prior{i} = betapdf(theta, old_a, old_b);

    % Likelihood
    Likelihood{i} = betapdf(theta, observations(i)+1, 1-observations(i)+1);

    % Update posterior
    a = (1-forgetting)*old_a + learning*observations(i); 
    b = (1-forgetting)*old_b + learning*(1 - observations(i)); 
    Y_posterior{i} = betapdf(theta, a, b);

    % Entropy
    temp = Y_prior{i}.*log(Y_prior{i}./binWidth);
    temp(isnan(temp)) = 0;
    entropy(i) = -sum(temp);

    % Update naive posterior
    a_naive = 1 + observations(i); 
    b_naive = 1 + (1 - observations(i)); 
    Y_posterior_naive{i} = betapdf(theta, a_naive, b_naive);

    % Predictive Surprise (Shannon surprise)
    post_pred = [old_b/(old_a+old_b), old_a/(old_a+old_b)];
    post_pred = post_pred(observations(i)+1);
    predictive_surprise(i) = -log(post_pred);

    % Bayesian surprise
    %bayesian_surprise(i) = KLDiv(Y_prior{i}, Y_posterior{i}); 
    bayesian_surprise(i) = KL_DIV([old_a old_b], [a b]); 

    % Confidence-corrected
    % KL(posterior, naive_posterior)
    % naive posterior: flat prior (without any observation) updated with
    % the new observation
    %confidence_surprise(i) = KLDiv(Y_prior{i}, Y_posterior_naive{i}); 
    confidence_surprise(i) = KL_DIV([old_a old_b], [a_naive b_naive]); 

    %PE
    [mean_beta var_beta] = betastat(old_a, old_b); 
%     mean_beta = old_a/(old_a + old_b); 
%     var_beta = (old_a*old_b)/((old_a+old_b+1)*(old_a+old_b)^2); 
% %     var_beta = (old_a*old_b)/(((old_a+old_b)^2)*(old_a+old_b+1));
%     SD_beta = sqrt((old_a*old_b)/(((old_a+old_b)^2)*(old_a+old_b+1)));
    PE(i) = observations(i) - mean_beta; 
    wPE(i) = (observations(i) - mean_beta)/var_beta;
%     sdwPE(i) = (observations(i) - mean_beta)/std(Y_prior{i}); 
    klPE(i) = KL_DIV([old_a old_b], [observations(i)+1 1-observations(i)+1]);

%     subplot(5,10,i)
%     %plot(var_beta(i))
%     hold on
%     plot(var_function(i))
%     %plot(Y_posterior{i})
%     hold off
end
% legend({'Posterior', 'Prior', 'Likelihood'})

% figure; 
% plot(var_beta)
% hold on
% plot(var_function)
% legend({'manual calculation' 'function', })

% figure; 
% plot(theta, Y_posterior{1}, '-b')
% hold on
% plot(theta, Y_posterior{10}, '-r')
% plot(theta, Y_posterior{30}, '-g')
% plot(theta, Y_posterior{50}, '-k')
% legend({'1 observation', '10 observations', '30 observations', '50 observations'})
% title('Beta distribution over observation randomly generated from state "0.25"')

% figure;
% observations1 = observations; 
% observations1(observations1 == 0) = nan; 
% plot(observations, 'xb');
% hold on
% plot(abs(PE))
% plot(abs(wPE))
% plot(abs(sdwPE))
% bayesian_surprise(bayesian_surprise==inf) = 1;
% confidence_surprise(confidence_surprise==inf) = 1;
% predictive_surprise(predictive_surprise==inf) = 1;
% plot(bayesian_surprise)
% plot(predictive_surprise)
% plot(confidence_surprise)
% title('symbolic train of 50 trials')
% legend({'Observations', 'PE', 'wPE', 'sdwPW', 'bs', 'ps', 'cs'})

figure;
subplot(6,1,1)
observations1 = observations; 
observations1(observations1 == 0) = nan; 
plot(observations, 'xb');
title('Observations')
subplot(6,1,2)
plot(abs(PE))
title('PE')
subplot(6,1,3)
plot(abs(wPE))
title('wPE')
% bayesian_surprise(bayesian_surprise==inf) = 1;
% confidence_surprise(confidence_surprise==inf) = 1;
% predictive_surprise(predictive_surprise==inf) = 1;
subplot(6,1,4)
plot(bayesian_surprise)
title('Bayesian Surprise')
subplot(6,1,5)
plot(predictive_surprise)
title('Predictive Surprise')
subplot(6,1,6)
plot(confidence_surprise)
title('Condifence Surprise')
% title('symbolic train of 50 trials')
% legend({'Observations', 'PE', 'wPE', 'sdwPW', 'bs', 'ps', 'cs'})

% figure;
% observations1 = observations; 
% observations1(observations1 == 0) = nan; 
% plot(observations, 'xb');
% hold on
% plot(abs(PE))
% plot(abs(wPE))
% plot(abs(sdwPE))
% bayesian_surprise(bayesian_surprise==inf) = 1;
% confidence_surprise(confidence_surprise==inf) = 1;
% predictive_surprise(predictive_surprise==inf) = 1;
% plot(zscore(bayesian_surprise))
% plot(zscore(predictive_surprise))
% plot(zscore(confidence_surprise))
% title('Zscored -- symbolic train of 50 trials')
% legend({'Observations', 'PE', 'wPE', 'sdwPW', 'bs', 'ps', 'cs'})
% 

%% RUN REAL SIMULATIONS (to modify)5state = 0.25; % %
trials = 50;
sequences = 1000; 
observations = rand(sequences,trials);
observations(observations <= state) = 1;
observations(observations ~= 1) = 0; 

Y_posterior = []; 
Y_prior = []; 
PE = []; 
wPE = []; 
predictive_surprise = []; 
bayesian_surprise = []; 
confidence_surprise = []; 

for i = 1:size(observations,1)
    a = 1; 
    b = 1; 
    theta = linspace(0,1,100); 
    binWidth = diff(theta);
    binWidth = [binWidth(end),binWidth];
    for y = 1:size(observations,2)
        old_a = a;
        old_b = b; 
        Y_prior{i,y} = betapdf(theta, old_a, old_b);
    
        Likelihood{i,y} = betapdf(theta, observations(i,y)+1, 1-observations(i,y)+1);
    
        % Update posterior
        a = (1-forgetting)*old_a + learning*observations(i); 
        b = (1-forgetting)*old_b + learning*(1 - observations(i)); 
        Y_posterior{i,y} = betapdf(theta, a, b); 
    
        % Entropy
        temp = Y_prior{i,y}.*log(Y_prior{i,y}./binWidth);
        temp(isnan(temp)) = 0;
        entropy(i,y) = -sum(temp);
    
        % Update naive posterior
        a_naive = 1 + observations(i,y); 
        b_naive = 1 + (1 - observations(i,y)); 
        Y_posterior_naive{i,y} = betapdf(theta, a_naive, b_naive);
    
        % Predictive Surprise (Shannon surprise)
        post_pred = [old_b/(old_a+old_b), old_a/(old_a+old_b)];
        post_pred = post_pred(observations(i,y)+1);
        predictive_surprise(i,y) = -log(post_pred);
    
        % Bayesian surprise
        %bayesian_surprise(i) = KLDiv(Y_prior{i}, Y_posterior{i}); 
        bayesian_surprise(i,y) = KL_DIV([old_a old_b], [a b]); 
    
        % Confidence-corrected
        % KL(posterior, naive_posterior)
        % naive posterior: flat prior (without any observation) updated with
        % the new observation
        %confidence_surprise(i) = KLDiv(Y_prior{i}, Y_posterior_naive{i}); 
        confidence_surprise(i,y) = KL_DIV([old_a old_b], [a_naive b_naive]); 
    
        %PE
        mean_beta = old_a/(old_a + old_b); 
        var_beta = (old_a*old_b)/((old_a+old_b)^2*(old_a+old_b+1));
        SD_beta = sqrt((old_a*old_b)/((old_a+old_b)^2*(old_a+old_b+1)));
        PE(i,y) = observations(i,y) - mean_beta; 
        wPE(i,y) = (observations(i,y) - mean_beta)/SD_beta; 
    end
end

matrix0 = observations; 
matrix0(matrix0 == 1) = nan;
matrix0(matrix0 == 0) = 1; 
matrix1 = observations; 
matrix1(matrix1 == 0) = nan; 



%This plot is to have a better understanding of how different measures
%might capture different aspects of learning the probabilistic underlying
%state. 
%As one notices, prediction error early on captures differences between
%having one stimulation or the other administered. 
figure;
%PE
subplot(5,1,1)
shadedErrorBar(1:size(observations,2), mean(abs(PE.*matrix1), 'omitnan'), var(abs(PE.*matrix1), 'omitnan'), 'or-')
hold on
shadedErrorBar(1:size(observations,2), mean(abs(PE.*matrix0), 'omitnan'), var(abs(PE.*matrix0), 'omitnan'), 'ob-')
title('PE')
%wPE
subplot(5,1,2)
shadedErrorBar(1:size(observations,2), mean(abs(wPE.*matrix1), 'omitnan'), var(abs(wPE.*matrix1), 'omitnan'), 'or-')
hold on
shadedErrorBar(1:size(observations,2), mean(abs(wPE.*matrix0), 'omitnan'), var(abs(wPE.*matrix0), 'omitnan'), 'ob-')
title('wPE')
%Predictive surprise
subplot(5,1,3)
shadedErrorBar(1:size(observations,2), mean(abs(predictive_surprise.*matrix1), 'omitnan'), var(abs(predictive_surprise.*matrix1), 'omitnan'), 'or-')
hold on
shadedErrorBar(1:size(observations,2), mean(abs(predictive_surprise.*matrix0), 'omitnan'), var(abs(predictive_surprise.*matrix0), 'omitnan'), 'ob-')
title('Predictive surprise')
%Bayesian surprise
subplot(5,1,4)
shadedErrorBar(1:size(observations,2), mean(abs(bayesian_surprise.*matrix1), 'omitnan'), var(abs(bayesian_surprise.*matrix1), 'omitnan'), 'or-')
hold on
shadedErrorBar(1:size(observations,2), mean(abs(bayesian_surprise.*matrix0), 'omitnan'), var(abs(bayesian_surprise.*matrix0), 'omitnan'), 'ob-')
title('Bayesian surprise')
%Confidence surprise
subplot(5,1,5)
shadedErrorBar(1:size(observations,2), mean(abs(confidence_surprise.*matrix1), 'omitnan'), var(abs(confidence_surprise.*matrix1), 'omitnan'), 'or-')
hold on
shadedErrorBar(1:size(observations,2), mean(abs(confidence_surprise.*matrix0), 'omitnan'), var(abs(confidence_surprise.*matrix0), 'omitnan'), 'ob-')
title('Confidence corr surprise')






