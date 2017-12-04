
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% All rights reserved. Jonas Gros, Yasuyuki Zushi, and J. Samuel Arey.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% -----------------------------------------------------------------------------
% ALL PROGRAM PARAMETERS THAT SHOULD BE SET BY THE USER ARE SHOWN HERE.

% INSTRUMENT PARAMETERS
    
    % Precision for rounding of M/S values (e.g. 0.01)
    Precision=0.01;
    
    % Threshold value of signal intensity
    Intthreshold = 0;
    
    % The number of pixels for 2nd dimension (Target chromatogram)
    NbPix2ndD_target=160;
    
    % The number of pixels for 2nd dimension (Reference chromatogram)
    NbPix2ndD_Ref=160;
    
    % Value for drift mass spectrum. if it is set to "1", all the m/z 
    % values will be added "1" (e.g., m/z 300 => m/z 301)
    driftMS = 0;

% MODEL CHOICE PARAMETERS

    % What is the typical width of a peak in first and second dimension
    % (in units of pixels by default; unit must match the unit selected below): 
    typical_peak_width = [1,5];
     
    % Model choice ('normal' or 'DualSibson'):
    model_choice = 'normal';

% INPUT/OUTPUT PARAMETERS

    % Set plot_flag to a value of 0 to suppress plots.
    % Set to a value of 1 to see "normal" level of plotting (DEFAULT).
    plot_flag = 1;

    % Set the output file path
    output_path = 'users/output/';

    % Set the input file path
    input_path = 'users/input/';

    % Name of input-output files:
    % Reference chromatogram:
    Reference_chromatogram_file = 'Test_reference_chromatogram.cdf';

    % Target chromatogram:
    Target_chromatogram_file = 'Test_target_chromatogram.cdf';

    % Positions of alignment points in the Reference chromatogram:
    Reference_alignment_pts_file = 'Alignment_pts_Reference.csv';

    % Positions of alignment points in the Target chromatogram:
    Target_alignment_pts_file = 'Alignment_pts_Target.csv';

    % Set Matlab console output level. Choose: 'minimal', 'normal', or 'verbose'.
    prompt_output = 'normal';

% UNITS
    
    % units for "typical_peak_with" and for positions of the alignment 
    % points ('time' means (minutes (1st dimension) and seconds (2nd
    % dimension), 'pixel' means pixels as defined in the documentation:
    units = 'pixel';
    
% -----------------------------------------------------------------------------
% Do not modify the lines below.

cd('..');

addpath model_code

cd('model_code')

run main_code;

cd('../users');



