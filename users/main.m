
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% � All rights reserved. Jonas Gros, Yasuyuki Zushi, and J. Samuel Arey.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% -----------------------------------------------------------------------------
% ALL PROGRAM PARAMETERS THAT SHOULD BE SET BY THE USER ARE SHOWN HERE.

% INSTRUMENT PARAMETERS
    
    % Precision for rounding of M/S values (e.g. 0.01):
    Precision=0.01;
    
    % YASUYUKI PLEASE COMMENT:
    Intthresholdafter = 0;
    % YASUYUKI PLEASE COMMENT:
    Intthresholdbefore = 0; 
    % YASUYUKI PLEASE COMMENT:
    NbPix2ndD_target=160;
    % YASUYUKI PLEASE COMMENT:
    NbPix2ndD_Ref=160;
    % YASUYUKI PLEASE COMMENT:
    driftMS = 0;


% MODEL CHOICE PARAMETERS

    % What is the typical width of a peak in first and second dimension
    % (in units of pixels) 
    typical_peak_width = [1,5];

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
    Reference_chromatogram_file = 'LakeOntario_NormalOvenGCxGC1424rounded_intthr0.01_MSdrifted0.118.cdf';

    % Target chromatogram:
    Target_chromatogram_file = 'LakeOntario_HighTempOvenGCxGC1354rounded_intthr0.01_MSdrifted0.354.cdf';

    % Positions of alignment points in the Reference chromatogram:
    Reference_alignment_pts_file = 'Alignment_pts_Reference.csv';

    % Positions of alignment points in the Target chromatogram:
    Target_alignment_pts_file = 'Alignment_pts_Target.csv';

    % Set Matlab console output level. Choose: 'minimal', 'normal', or 'verbose'.
    prompt_output = 'normal';

% -----------------------------------------------------------------------------
% Do not modify the lines below.

cd('..');

addpath model_code

cd('model_code')

run main_code;

cd('../users');



