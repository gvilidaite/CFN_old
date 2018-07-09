function gv_fit_dipper_better

close all
clear all
% GV 10/2/16
% this version fits a model for each observer x number of times and then
% picks the one with the smallest error
% uses 2 free parameters


col = 'r';
no_of_fits = 50;
no_of_sims = 10000; % number of simulations
delete_outliers = 0;
outlierstd = 3;
dipperlevelsdB = [-18 6 18 30];
toplot = 1; % 1 = show a figure for each observer; 0 = nope

load('NFCdata.mat')

dippers = Cdippers;
AQ = AQmatrix(:,2);

no_of_observers = length(dippers);

if delete_outliers == 1     
    
    % finding the outliers
    del = 0;
    for x = 1:4
        themean(x) = mean(dippers(:,x));
        thestd(x) = std(dippers(:,x));
        therow = 0;
        for y = dippers(1:end,x)'
            therow = therow + 1;
            if (y < (themean(x)-(thestd(x).*outlierstd))) | (y > (themean(x) + (thestd(x).*outlierstd)))
                del = del + 1;
                deleterows3(del) = therow;
                
            end
        end
    end
    deleterows3 = unique(deleterows3)
end

no_of_observers_left = length(AQ)

dippers(no_of_observers_left+1,:) = mean(dippers,1); % adding the mean dipper on the end of the matrix so I can do the fit and later take it off the matrix

for f = 1:length(dippers)
    dipper = dippers(f,:);
    figure(f);
    %set(gcf,'Position',[100 50 950 900], 'color', [1 1 1]);
    %%%%%%%%%%%%%%% dipper %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on
    plot(dipperlevelsdB,dipper,'ko','LineWidth', 2, 'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',7);
    % calculating parameters until error (diperr) is less than 1:
    k = 0;
    for m = 1:no_of_fits
        k = k + 1;
        ap = fminsearch(@getdippererror, rand(1,4), [], dipperlevelsdB, dipper, 0); % putting in 4 random numbers as starting parameters

        lap = 10.^(ap/20); % antilog that doesn't let parameters be negative
        K.Z(m) = lap(1);
        K.n(m) = lap(2);
        K.q = lap(3);
        K.p = K.q + lap(4);

        K.q = 2;
        K.p = 2.4; % changing these parameters to canonical values is better for
        % comparisons between subjects and between measures (e.g. with noise masking)
        %K = [K.Z, K.n, K.q, K.p];
        dippermodel(m,:) = getdippermodel(ap, [-21:3:36], 0);
        diperrtemp(m) = getdippererror(ap, dipperlevelsdB, dipper, 0);
    end
    [diperr(f) ind] = min(diperrtemp);
    Z_param(f) = K.Z(ind);
    n_param(f) = K.n(ind);
    dippermodels(f,:) = dippermodel(ind,:);
    
if toplot == 1
    plot([-21:3:36], dippermodel(ind,:), 'k-', 'LineWidth', 2, 'color', col);
    xlabel('Pedestal contrast (dB)', 'FontSize', 12);
    ylabel('Threshold (dB)', 'FontSize', 12);
    title('Pedestal Masking', 'FontSize', 15);
end
    % simulating dipper data for noise masking and double-pass predictions
    [simthresh, simcorrect, simconsistent] = dodippersimulations(ap,dipperlevelsdB,no_of_sims,0); % use 10000 trials or so
    
    sim_dp_cor(f,:,:) = simcorrect;
    sim_dp_con(f,:,:) = simconsistent;
end


meandipper = dippers(end,:)
dippers(end,:) = [];
meandipperfit = dippermodels(end,:);
dippermodels(end,:) = [];
dipperSE = std(dippers,1,1)./sqrt(length(Z_param));
diperr

% save('AQ_dipper.mat', 'n_param', 'Z_param', 'AQ');
% save('dipper_splithalf', 'Z_param', 'n_param', 'diperr', 'dippers','dippermodels','meandipper','dipperSE','simthresh', 'sim_dp_cor', 'sim_dp_con','meandipperfit');


end
%--------------------------------------------------------------------------
function dippermodel = getdippermodel(ap, contrasts, type)
% this function produces a predicted dipper curve without touching the data
%
% ap = parameters (give four numbers)
% contrasts = vector of contrast levels in LINEAR VALUES                                               
% type = put in 0 if contrast levels are in dB OR 1 if contrast levels are linear
%
% noisemodel comes out in dB


ap = 10.^(ap/20); % antilog that doesn't let parameters be negative
K.Z = ap(1);
K.n = ap(2);
K.q = ap(3);
K.p = K.q + ap(4);

    K.q = 2;
    K.p = 2.4; % changing these parameters to canonical values is better for

dipperthresholds = getdipperthresh(K, contrasts, type);
dippermodel = dipperthresholds;  % e.g. if threshold is 10, log10(10) is 1, so noisemodel = 20

end
%--------------------------------------------------------------------------
function dipperthresholds = getdipperthresh(K, pedlevels, type)
% this function is used by getdippermodel and finds the thresholds for it
if (type == 0)
    pedlevels = 10.^(pedlevels./20);  % if contrasts are in dB convert them into linear units
    pedlevels(1) = 0;
end

null = (pedlevels.^K.p)./(K.Z + pedlevels.^K.q); % null interval (pedestal only condition)


for n = 1:length(pedlevels) % does this for each pedestal level
    target(n) = 0; % internal noise level or target contrast ?????????
    resp(n) = -1; % target + pedestal response set to a negative number ???
    while ((resp(n)-null(n)) < K.n)
    % loop set to terminate when the difference between the response to target
    % + pedestal becomes more than the noise parameter -- WHY???
    % this loop estimates the threshold with an error of 0.1
    % the threshold is the value of 'target'
        target(n) = target(n) + 0.1; % increasing the target contrast by small amounts each time
        resp(n) = ((pedlevels(n)+target(n)).^K.p)./(K.Z + (pedlevels(n)+target(n)).^K.q);
    end
    while ((resp(n)-null(n)) > K.n)
    % loop set to terminate when resp-null is less than the noise parameter
    % estimates threshold with and error of 0.001
        target(n) = target(n) - 0.001;
        resp(n) = ((pedlevels(n)+target(n)).^K.p)./(K.Z + (pedlevels(n)+target(n)).^K.q);
    end
end
dipperthresholds = 20.*log10(target);


end
%--------------------------------------------------------------------------
function rms = getdippererror(p, contrasts, thresholds, type)
% this function checks how good is the fit between the predicted curve and
% the data
% p = parameters (p1 = internal noise; p2 = external noise)
% contrasts = vector of contrast levels in LINEAR VALUES
% thresholds = the threshold values from the experiment for each of the
% contrast levels
% type = put in 0 if contrast levels are in dB OR 1 if contrast levels are linear
 


% uses the above function to get model predictions first:
dippermodel = getdippermodel(p, contrasts, type);
rms = sqrt(mean((dippermodel-thresholds).^2));

end
%--------------------------------------------------------------------------
function [simthresh, simcorrect, simconsistent] = dodippersimulations(ap,maskC,ntrials,type)
% this function is modified from the noise masking data simulation subfunction
% it takes the parameters of best fit from the dipper
% model and simulates a bunch of data. These data can be used both for
% plotting a simulated noise masking function and for working out double
% pass consistency and accuracy
%
% the simulated double pass measures can be used to compare with our actual
% double pass experiment measures
%
% simthresh = simulated thresholds for plotting a simulated noise masking
% function
% simcorrect and simconsistent = double pass measures
% k = parameters from model fit
% maskC = contrast levels
% ntrials = how many trials we want to simulate
% type = put in 0 if contrast levels are in dB OR 1 if contrast levels are linear

if (type == 0)
    maskC = 10.^(maskC./20);  % if contrasts are in dB convert them into linear units
elseif (type == 1)
    maskC = maskC;
end


ap = 10.^(ap/20); % antilog that doesn't let parameters be negative
K.Z = ap(1);
K.n = ap(2);
K.q = ap(3);
K.p = K.q + ap(4);
% 
    K.q = 2;
    K.p = 2.4; % changing these parameters to canonical values is better for

targetdB = -18:3:39; % target contrast levels for which we will simulate the data
targetC = 10.^(targetdB./20);% change to linear units

for m = 1:length(maskC)  % loop to simulate data for each pairing of mask and target levels
    for n = 1:length(targetC)
        ncorrect = 0; % this is where we add up all correct and consistent results
        nconsistent = 0;
        for o = 1:ntrials
        % simulating the stimuli that we 'showed' for each trial:
        int1 = randn.*maskC(m) + targetC(n); % interval-1-of-a-trial = (random-noise * mask-contrast) + (target-contrast * external-noise-parameter)
        int2 = randn.*maskC(m); % interval-1-of-a-trial = (random-noise * mask-contrast)
        
        sint1 = sign(int1);
        sint2 = sign(int2);
        int1 = abs(int1);
        int2 = abs(int2);
        
        int1 = sint1.*(int1.^K.p)./(K.Z + int1.^K.q);
        int2 = sint2.*(int2.^K.p)./(K.Z + int2.^K.q);

        
        
        % calculating the 'neural' response to the trial:
        resp1a = int1 + K.n.*randn; % 1a = interval 1 of pass 1; 1b = interval 2 of pass 1; etc..
        resp1b = int2 + K.n.*randn; % each interval of each pass will have slightly different responses due to random amounts of internal noise (p(1))
        resp2a = int1 + K.n.*randn; 
        resp2b = int2 + K.n.*randn;


            if resp1a>resp1b  % comparing the 'neural' responses to get a choice between intervals
                ncorrect = ncorrect + 1;
                pass1 = 1;
            else
                pass1 = 0;
            end
            if resp2a>resp2b % same for the 2nd pass
                ncorrect = ncorrect + 1;
                pass2 = 1;
            else
                pass2 = 0;
            end
            if pass1==pass2 % comparing the two passes for consistent/inconsistent
                nconsistent = nconsistent + 1;
            end


    end
        simcorrect(m,n) = ncorrect/(ntrials*2);
        simconsistent(m,n) = nconsistent/ntrials;
    end
end
% for m = 1:length(maskC)
%     a = find(simcorrect(m,:)<1);  % finds values of less than 1 and saves info about their location into a. This is just to make function interp1 work properly
%     simcorrect(m,a);
%     simthresh(m) = interp1(simcorrect(m,a), targetdB(a), 0.81625) % indice 'a' tells the function which items had values of less than 1
%                  % interp1 looks up the y value for a given x value
% 
% end
simthresh = [];

end