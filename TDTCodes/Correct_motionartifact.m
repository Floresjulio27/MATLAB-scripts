function [expCorrected] =  Correct_motionartifact(expCorrected,Params)    

    % Normalize using "isosbestic" channel
    % rescale 465 or 560 signal to remove motion artifact using 405

    % as there may be spikes, which will totally change the polynominal fit
    % first conduct a median filter to remove potential spikes
    filt_window = 5;%5;
    Medfit.F405 = medfilt1(expCorrected.F405, filt_window*Params.DataFs);
    Medfit.F465 = medfilt1(expCorrected.F465, filt_window*Params.DataFs);
    Medfit.F560 = medfilt1(expCorrected.F560, filt_window*Params.DataFs);

    Residual.F405 = expCorrected.F405-Medfit.F405;
    Residual.F465 = expCorrected.F465-Medfit.F465;
    Residual.F560 = expCorrected.F560-Medfit.F560;

    spike_thresh = 3;
    % removeSpike.F405 = Residual.F405 >= spike_thresh*std(Medfit.F405);
    % removeSpike.F465 = Residual.F465 >= spike_thresh*std(Medfit.F465);
    % removeSpike.F560 = Residual.F560 >= spike_thresh*std(Medfit.F560);

    removeSpike.F405 = Residual.F405 >= spike_thresh*(max(Medfit.F405)-min(Medfit.F405));
    removeSpike.F465 = Residual.F465 >= spike_thresh*(max(Medfit.F465)-min(Medfit.F465));
    removeSpike.F560 = Residual.F560 >= spike_thresh*(max(Medfit.F560)-min(Medfit.F560));

    if any(removeSpike.F405)
        disp('Spikes in 405 channel');
        expCorrected.F405(removeSpike.F405) = Medfit.F405(removeSpike.F405);
    end
    if any(removeSpike.F465)
        disp('Spikes in 465 channel');
        expCorrected.F465(removeSpike.F465) = Medfit.F465(removeSpike.F465);
    end
    if any(removeSpike.F560)
        disp('Spikes in 560 channel');
        expCorrected.F560(removeSpike.F560) = Medfit.F560(removeSpike.F560);
    end
    