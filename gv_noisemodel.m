% Modelling script to estimate noise by using linear amplifier model. It takes into account both the
% consistency and accuracy of double-pass in CFN (contrast, faces, numbers) experiments. Makes
% predictions for finely sampled internal noise levels and then finds the
% best fitting simulation to fit each subjects data using the Pythagorean
% theorem.
% 
% Seems to work with linear noise levels way better than log
%
% GV 3/6/16

clear all
close all

tic

load('NFCdata.mat')

nsims = 100000; % number of simulations to do
intnoiselevelsdB = -40:0.1:40; % specify either this or intnoiselevels
subjectToPlot = 20; % pick a subject to plot as an example

if isempty(intnoiselevelsdB) == 0
    intnoiselevels = d_dBtoPercent(intnoiselevelsdB);
else
    intnoiselevels = 0:0.1:15; % use this if don't want noise to be logarithmic
end


fig = 0; % just helps keep count of the figures when plotting

for exp = 1:3 % refers to contrast, faces and numbers experiments; needs to be changed to 1:3 when script is finished
    disp(strcat('Modelling data for experiment ',num2str(exp)))
    switch exp
        case 1 % contrast experiment
            extnoiseSD = d_dBtoPercent(12); % stimulus noise standard deviation in linear units
            target_conds = [0 d_dBtoPercent(12)]; % experiment target conditions (no target and 12dB target in linear units)
            pedestal = 0;
        case 2 % faces exp
            extnoiseSD = d_dBtoPercent(24);
            target_conds = [0 d_dBtoPercent(24)];
            pedestal = d_dBtoPercent(30);
        case 3 % number summation exp
            extnoiseSD = 10; % standard deviation from numbers experiment
            target_conds = [0 10]; % makes the means of the numbers 50 and 60
            pedestal = 50;
    end
            
            
    correct = zeros(length(intnoiselevels),length(target_conds)); % matrix of zeros to put scores in; dimensions: experimental conds, internal noise levels
    consistent = zeros(length(intnoiselevels),length(target_conds));

    for a = 1:length(intnoiselevels)
        intnoiseSD = intnoiselevels(a);
    for cond = 1:length(target_conds) % 1:2 target conditions
        target = target_conds(cond);
        for n = 1:nsims

            int1 = extnoiseSD .* randn + pedestal + target;
            int2 = extnoiseSD .* randn + pedestal;

            % does the randn internal noise value need to be the same for both responses or different?
            resp11 = intnoiseSD .* randn + int1; % first interval of first pass
            resp12 = intnoiseSD .* randn + int2; % second interval of first pass
            resp21 = intnoiseSD .* randn + int1; % first interval of second pass
            resp22 = intnoiseSD .* randn + int2; % second interval of second pass        
            % calculating correct answer for pass 1:
            if resp11 > resp12
               correct(a,cond) = correct(a,cond) + 1; % if answer is correct, add 1 to the appropriate condition and int noise level
               pass1 = 1;
            else
               pass1 = 0;
            end
            % calculating correct answer for pass 2:
            if resp21 > resp22
               correct(a,cond) = correct(a,cond) + 1; % if answer is correct, add 1 to the appropriate condition and int noise level
               pass2 = 1;
            else
               pass2 = 0;
            end

            if pass1 == pass2
               consistent(a,cond) = consistent(a,cond) + 1;
            end
        end
    end
    end

    simcorr = correct ./ (2*nsims); % dividing by 2 times nsims because of two passes
    simcons = consistent ./ nsims;

    % plotting simulation spread:
    fig = fig + 1;
    figure(fig)
    hold on;
    for a = 1:length(intnoiselevels)
        plot([simcons(a,1) simcons(a,2)],[simcorr(a,1), simcorr(a,2)], 'k-', 'Color', [rand rand rand], 'LineWidth', 1);
    end
    
    if exp == 1
        plot([Cpcons0(subjectToPlot) Cpcons2(subjectToPlot)],[Cpcorr0(subjectToPlot), Cpcorr2(subjectToPlot)], 'ko-', 'LineWidth', 2);
    elseif exp == 2
        plot([Fpcons0(subjectToPlot) Fpcons2(subjectToPlot)],[Fpcorr0(subjectToPlot), Fpcorr2(subjectToPlot)], 'ko-', 'LineWidth', 2);
    else
        plot([Npcons0(subjectToPlot) Npcons2(subjectToPlot)],[Npcorr0(subjectToPlot), Npcorr2(subjectToPlot)], 'ko-', 'LineWidth', 2);
    end
    plot([0 1],[0.5 0.5],'k:');
    plot([0.5 0.5],[0 1],'k:');
    protpcpay = 0.5:0.01:1;
    protpcpax = 0.5 + 2.*(protpcpay-0.5).^2;      % from Klein & Levi, 2009
    plot(protpcpax,protpcpay,'k-');

    axis square;
    axis([0 1 0 1]);
    xlabel('Proportion consistent', 'FontSize', 18);
    ylabel('Proportion correct', 'FontSize', 18);

% finding the best fitting values for each subject:
for sub = 1:length(Cpcons0)
    for cond = 1:2
        
        switch exp
            case 1
                switch cond
                    case 1
                        currcons = Cpcons0(sub); % data points to compare to simulated values
                        currcorr = Cpcorr0(sub);
                    case 2
                        currcons = Cpcons2(sub);
                        currcorr = Cpcorr2(sub);
                end
            case 2
                switch cond
                    case 1
                        currcons = Fpcons0(sub); % data points to compare to simulated values
                        currcorr = Fpcorr0(sub);
                    case 2
                        currcons = Fpcons2(sub);
                        currcorr = Fpcorr2(sub);
                end
            case 3
                switch cond
                    case 1
                        currcons = Npcons0(sub); % data points to compare to simulated values
                        currcorr = Npcorr0(sub);
                    case 2
                        currcons = Npcons2(sub);
                        currcorr = Npcorr2(sub);
                end
        end

        for comp = 1:length(simcons) % for comparison in number of noise levels
            x = simcons(comp,cond); % simulated value of consistency to compare to
            y = simcorr(comp,cond);
            err(cond,comp) = sqrt((abs(x - currcons)).^2 + (abs(y - currcorr)).^2); % Pythagorean theorem to find the distance between points
        end   
    end % end of condition loop
    meanerr = mean(err,1); % taking the mean error over the two conditions so that we select the noise level which reflects performance in both conds
    [error i] = min(meanerr);
    if isempty(intnoiselevelsdB) == 0 % if we are doing the simulations in dB then index the dB variables
        subnoise(exp,sub) = intnoiselevelsdB(i);
    else % otherwise, index the linear noise level variable
        subnoise(exp,sub) = intnoiselevels(i);
    end
end % end of fitting subject data loop
end % end of experiments loop
intnoise = subnoise'

toc
