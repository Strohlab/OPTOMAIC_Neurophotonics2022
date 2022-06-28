function [signal_corr]=DC_Baseline(signal,ws,ss)
% Correct Ca2+-2Photon-Data to bring it to its baseline
% Input: uncorrected signal-vector

 % create time-axis
 time=[1:numel(signal)]';
 % correct signal using the matlab-function msbackadj for signal with peaks
 % to do: test different properties for this function
 signal_corr=DC_msbackadj(time,signal,'WindowSize',ws,'StepSize',ss);
 % calculate an aditional shift using the standard deviation and
 % subtract this shift partially to consider peaks
 % to do: calculate weighting-factor based on the curve (no. of peaks...)
 % based perhaps on this: 
 % signal_corr=signal_corr-std(signal_corr*0.05);
 
 