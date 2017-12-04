% main_code.m script called from "main.m" by user.
% 
% *** Do not modify this file. Normally the user should not need to adjust *** 
% *** anything in this script.                                             ***
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% © All rights reserved. Jonas Gros, Yasuyuki Zushi, and J. Samuel Arey.
% See academic license terms stated in file:
% LICENSE.txt
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

% addpath model_code

% Import chromatograms:


% A few error checks:
if ~exist(['../',input_path,Target_chromatogram_file],'file')
    error(['Target chromatogram file (',['../users/input/',Target_chromatogram_file],' ) does not exist.'])
end
if ~exist(['../',input_path,Reference_chromatogram_file],'file')
    error(['Reference chromatogram file (',['../users/input/',Reference_chromatogram_file],' ) does not exist.'])
end
if ~exist(['../',output_path],'dir')
    warning(['Indicated output directory, ''',...
        output_path,''' does not exist... Creating it.'])
    mkdir(['../',output_path])
end
% Check if allowed to write to the output_path:
[stat,CcC] = fileattrib(['../',output_path]); 
if ~CcC.UserWrite
    error(['Matlab is not allowed to write to the output_path folder...',...
        ' Try displacing the folder! (e.g. to ''Desktop'' on a Windows computer)'])
end

ChromatoTarget = matlab_cdf_open_function([strrep(pwd,'model_code',''),input_path,...
    Target_chromatogram_file],Intthreshold,driftMS);
ChromatoRef = matlab_cdf_open_function([strrep(pwd,'model_code',''),input_path,...
    Reference_chromatogram_file],Intthreshold,driftMS);
    
MPeriod_target = NbPix2ndD_target/ChromatoTarget.SamRate;
MPeriod_Ref = NbPix2ndD_Ref/ChromatoRef.SamRate;


% Import the positions of the alignment points:
% (Including a check whether the file do actually exist)
if ~exist([strrep(pwd,'model_code',''),input_path,...
    Reference_alignment_pts_file],'file')
    if exist([strrep(pwd,'model_code',''),output_path,...
    Reference_alignment_pts_file],'file')
        warning(['The file containing reference alignment points',...
            ' was found in the output_path and not input_path.',...
            ' Loading the file in output_path...'])
        Reference_peaks = importdata([strrep(pwd,'model_code',''),output_path,...
        Reference_alignment_pts_file]);
    else
        error(['The file containing reference alignment points (',...
            [strrep(pwd,'model_code',''),input_path,...
        Reference_alignment_pts_file],' ) does not exist.'])
    end
else
    Reference_peaks = importdata([strrep(pwd,'model_code',''),input_path,...
        Reference_alignment_pts_file]);
end

if ~exist([strrep(pwd,'model_code',''),input_path,...
    Target_alignment_pts_file],'file')
    if exist([strrep(pwd,'model_code',''),output_path,...
    Target_alignment_pts_file],'file')
        warning(['The file containing target alignment points',...
            ' was found in the output_path and not input_path.',...
            ' Loading the file in output_path...'])
        Target_peaks = importdata([strrep(pwd,'model_code',''),output_path,...
        Target_alignment_pts_file]);
    else
        error(['The file containing target alignment points (',...
            [strrep(pwd,'model_code',''),input_path,...
        Target_alignment_pts_file],' ) does not exist.'])
    end
else
    Target_peaks = importdata([strrep(pwd,'model_code',''),input_path,...
        Target_alignment_pts_file]);
end

if strcmpi(units,'time')
    % Store the positions in time units:
    Target_peaks_time = Target_peaks;
    Reference_peaks_time = Reference_peaks;
    % COnvert them to pixel units:
    Target_peaks = Time2Pix(Target_peaks, NbPix2ndD_target/ChromatoTarget.SamRate, ChromatoTarget.SamRate, ChromatoTarget.RTini);
    Reference_peaks = Time2Pix(Reference_peaks, NbPix2ndD_Ref/ChromatoRef.SamRate, ChromatoRef.SamRate, ChromatoRef.RTini);
    typical_peak_width = Time2Pix(typical_peak_width, NbPix2ndD_Ref/ChromatoRef.SamRate, ChromatoRef.SamRate, ChromatoRef.RTini);
end

% TIC is here used as the Reference/Target --> a first "normal" alignment
% is done for the TIC (using the code normal code for non-MS data), and
% some results are re-used to align all MS-ions
% (like Displacement, and local deformation correction
% (Deform_output))
Ref_TIC = reshape([ChromatoRef.MStotint(1:(length(ChromatoRef.MStotint)/NbPix2ndD_Ref)*NbPix2ndD_Ref);...
    zeros(NbPix2ndD_Ref-mod(length(ChromatoRef.MStotint),NbPix2ndD_Ref),1)],...
    NbPix2ndD_Ref,[]);

Target_TIC = reshape([ChromatoTarget.MStotint(1:(length(ChromatoTarget.MStotint)/NbPix2ndD_target)*NbPix2ndD_target);...
    zeros(NbPix2ndD_target-mod(length(ChromatoTarget.MStotint),NbPix2ndD_target),1)],...
    NbPix2ndD_target,[]);

%  % Round m/z values:
[ChromatoTarget.MSvaluebox,ChromatoTarget.MSintbox] = MSdataRound_v2(ChromatoTarget.MSvaluebox,...
    ChromatoTarget.MSintbox,Precision);

tic
[AlignedMSvaluebox,AlignedMSintbox,~,~,Aligned,Displacement,Deform_output] = Align2DChrom_MS_v5(Ref_TIC,Target_TIC,Reference_peaks,Target_peaks,ChromatoTarget.MSvaluebox,ChromatoTarget.MSintbox,NbPix2ndD_Ref,'Peak_widths',typical_peak_width,'model_choice',model_choice);
toc

Alignedscannum = ChromatoTarget.scannum;
Alignedflag = ChromatoTarget.flag;
Alignedmedmsmax = ChromatoTarget.medmsmax;
Alignedmedmsmin = ChromatoTarget.medmsmin;
Alignedscantime = ChromatoTarget.scantime;

%clear ChromatoTarget

Mat=AlignedMSvaluebox';
Mat2=AlignedMSintbox';

Vect=Mat(:);

AlignedMSvalueboxLine = Vect(Vect~=0);

Vecti=Mat2(:);

AlignedMSintboxLine = Vecti(Vect~=0);
% MSintbox as vector without zeros, except one for each empty line:

% Memory is a limiting factor:
clear Vect Vecti Mat Mat2

Alignedeachscannum= sum(AlignedMSvaluebox~=0,2);
Alignedionid = cumsum(Alignedeachscannum);

AlignedMStotint = sum(AlignedMSintbox,2);
toc

TICaligned = sum(AlignedMSintbox,2);
TIC_2Daligned=reshape([TICaligned;zeros(NbPix2ndD_target-mod(length(AlignedMStotint),NbPix2ndD_target),1)],...
     NbPix2ndD_target,[]);

if plot_flag 
    
    set(figure,'name','Total Ion Chromatograms (TICs)')
    subplot(3,1,1); 
    if ~strcmpi(units,'time')
        plotChromato(Ref_TIC);
        hold on
        plot(Reference_peaks(:,1),Reference_peaks(:,2),'ok','linewidth',2)
    else
        plotChromato(Ref_TIC,'MP',NbPix2ndD_target/ChromatoTarget.SamRate,'SR',ChromatoTarget.SamRate,'acquisition_delay',ChromatoTarget.RTini*60);
        hold on
        plot(Reference_peaks_time(:,1),Reference_peaks_time(:,2),'ok','linewidth',2)
    end
    c_a_xis = caxis*2; caxis(c_a_xis); title('\bfReference')
    ylabel('')
    xlabel('')
    
    subplot(3,1,2); 
    if ~strcmpi(units,'time')
        plotChromato(Target_TIC);
        hold on; 
        plot(Target_peaks(:,1),Target_peaks(:,2),'ok','linewidth',2)
    else
        plotChromato(Target_TIC,'MP',NbPix2ndD_target/ChromatoTarget.SamRate,'SR',ChromatoTarget.SamRate,'acquisition_delay',ChromatoTarget.RTini*60);
        hold on; 
        plot(Target_peaks_time(:,1),Target_peaks_time(:,2),'ok','linewidth',2)
    end
    caxis(c_a_xis); 
    title('\bfTarget'); 
    xlabel('')
    
    
    subplot(3,1,3); 
    if ~strcmpi(units,'time')
        plotChromato(TIC_2Daligned);
        hold on
        plot(Reference_peaks(:,1),Reference_peaks(:,2),'ok','linewidth',2)
    else
        plotChromato(TIC_2Daligned,'MP',NbPix2ndD_target/ChromatoTarget.SamRate,'SR',ChromatoTarget.SamRate,'acquisition_delay',ChromatoTarget.RTini*60);
        hold on
        plot(Reference_peaks_time(:,1),Reference_peaks_time(:,2),'ok','linewidth',2)
    end
    caxis(c_a_xis); title('\bfAligned')
    ylabel('')
    
    set(figure,'name','Difference chromatograms')
    subplot(3,1,2); 
    if ~strcmpi(units,'time')
        diffChromato(Ref_TIC,Target_TIC);
    else
        diffChromato(Ref_TIC,Target_TIC,'MP',NbPix2ndD_target/ChromatoTarget.SamRate,'SR',ChromatoTarget.SamRate,'acquisition_delay',ChromatoTarget.RTini*60);
    end 
    c_a_xis = caxis; 
    title('\bfReference - Target'); 
    xlabel('');
    
    subplot(3,1,3);  
    if ~strcmpi(units,'time')
        diffChromato(Ref_TIC,TIC_2Daligned);
    else
        diffChromato(Ref_TIC,TIC_2Daligned,'MP',NbPix2ndD_target/ChromatoTarget.SamRate,'SR',ChromatoTarget.SamRate,'acquisition_delay',ChromatoTarget.RTini*60);
    end 
    caxis(c_a_xis); 
    title('\bfReference - Aligned'); 
    ylabel('')
    
end

if strcmpi(prompt_output,'minimal')
    save_flag = 1;
else
    disp('Save aligned chromatogram? (1 = yes, 0 = no)')
    save_flag = input('');
end

if save_flag
    dot_indices = strfind(Target_chromatogram_file,'.');
    Output_file_name = [strrep(pwd,'model_code',''),output_path,...
        Target_chromatogram_file(1:dot_indices(end)-1),'_ALIGNED.cdf'];
    
    tic
    % Write cdf output file:
    cdf_write_function(Output_file_name ,Alignedscantime,AlignedMSintboxLine,Alignedflag,Alignedscannum,Alignedmedmsmax,Alignedmedmsmin,Alignedionid,Alignedeachscannum,AlignedMStotint,AlignedMSvalueboxLine)

    zu = toc;
    disp(' ')
    disp(['Done in ',num2str(zu),' seconds.'])
    disp(' ')
end


