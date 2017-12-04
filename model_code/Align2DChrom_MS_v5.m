function [AlignedMSvaluebox,AlignedMSintbox,Alignedeachscannum,Alignedmedionid,Aligned,Displacement,Deform_output] = Align2DChrom_MS_v5(Ref,Other,Peaks_Ref,Peaks_Other,MSvaluebox,MSintbox,NbPix2ndD,varargin)

% This file is a two-dimensional chromatograms alignment function, for MS
% data. It is an extension of the non-MS code of Gros et al. 2012.
% Inputs allowed
% are the same as for the normal (non MS) function. But includes some
% additional ones (described below).
% Refer to the normal (non-MS) code for more
% explanations about input & options that are the same.
% The main differences is that:
% Ref (the reference chromatogram) should be of size [m,n], and
% Other (the target chromatogram) should be of size [m,n] and will be used
% for the first alignment which determines the local
% Deformation correction & Displacement which will be later applied to 
% each ion. It could be the TIC, but might also be any ion chromatogram, as
% long as it contains enough common peaks with the corresponding Reference
% chromatogram.
% Peaks_Other should therefore refer to the position of alignment points in
% this chromatogram.
% MSvaluebox: it is the MS ions m/z value.
% MSintbox: it is the corresponding intensity value.
% NbPix2ndD: it is the number of pixels per 2nd column.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% See license terms stated in LICENSE.txt
% Authors : Jonas Gros, Yasuyuki Zushi, and J. Samuel Arey.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

% First, if ever the user used some options, give their values to the
% corresponding variables - It is just the same options as in the non-MS
% function):
% (Really interesting code starts around line 76).
for i=1:2:length(varargin)
    if strcmpi(varargin{i},'Interp_meth')
        Interp_meth=varargin{i+1};
    elseif strcmpi(varargin{i},'PowerFactor')
        PowerFactor=varargin{i+1};
    elseif strcmpi(varargin{i},'PlotVar')
        PlotVar=varargin{i+1};
    elseif strcmpi(varargin{i},'SaveVar')
        SaveVar=varargin{i+1};
    elseif strcmpi(varargin{i},'Peak_widths')
        Peak_widths=varargin{i+1};
        elseif strcmpi(varargin{i},'InterPixelInterpMeth')
        InterPixelInterpMeth=varargin{i+1};
    elseif strcmpi(varargin{i},'model_choice')
        model_choice=varargin{i+1};
    else
        error(['The option ''' varargin{i} ''' is unknown'])
    end
end

% Keep track of the amount of time needed for the program to run.
t_initial = now;
disp('alignment start time:')
disp(datestr(now));

% If Interp_meth was not given, set it to 'natural-neighbor':
if ~exist('Interp_meth','var')
     Interp_meth='natural-neighbor';
end

% If PowerFactor was not given, set its value to a default value of 2:
if ~exist('PowerFactor','var')
    PowerFactor=2;
end

% If Peak_widths was not given, it will be considered that no peak width
% correction has to be added (id est that distances will be considered in
% units of pixels), this is achieved by giving the same value to the two
% elements of Peak_widths:
if ~exist('Peak_widths','var')
     Peak_widths=[1,1];
end

% If InterPixelInterpMeth was not given, set it to cubic:
if ~exist('InterPixelInterpMeth','var')
%InterPixelInterpMeth='cubic';
InterPixelInterpMeth='linear';
end

% We don't want to plot each ion chromatogram:
     PlotVar=0;

%  Initialize the matrix for the aligned chromatogram:
Aligned=zeros(size(Other));

% First alignment (using the normal non-MS code):
tic % tic - toc, measures time elapsed
if strcmpi(model_choice,'normal')
[Aligned1,Displacement,Deform_output]=alignChromato(Ref,...
    squeeze(Other),Peaks_Ref,Peaks_Other,...
    'DisplacementInterpMeth',Interp_meth,'PowerFactor',PowerFactor,...
    'Peak_widths',Peak_widths,'InterPixelInterpMeth',InterPixelInterpMeth);
elseif strcmpi(model_choice,'DualSibson')
    [Aligned1,Displacement,Deform_output]=alignChromato_with_SibsonInterp_aslo_1stD(Ref,...
    squeeze(Other),Peaks_Ref,Peaks_Other,...
    'DisplacementInterpMeth',Interp_meth,'PowerFactor',PowerFactor,...
    'Peak_widths',Peak_widths,'InterPixelInterpMeth',InterPixelInterpMeth);
else
    error('unknown ''model_choice'' parameter value. Please choose either ''normal'' or ''DualSibson''') 
end
Aligned(:,:)=Aligned1;
disp('1st alignment in:')
toc

% We will compute the aligned chromatogram if taking each time the
% pixel at a different corner around the right position, than we will 
% perform some interpolation. 

% % Bilinear interpolation:
%
%     *c                       *d
%
%      |
%      u
%      |       
%      _       .p
%      |
%      t
%      |
%      *a  -r-  |    -  s  -   *d
% 
%  --> p = (r*b+s*a)*u + (r*d+s*c)*t
% 
% EXPLANATION: You want to interpolate the value at the position p in the
% Target chromatogram (that you need to put somewhere in the Aligned
% chromatogram). Therefore, you can interpolate a value based on the values
% of the four neighboring pixels (a,b,c,d). In the formula, p refers to
% interpolated intensity value, a,b,c,d refer to actual intensity values,
% and r,s,t,u refer to distances in pixels (where r+s = u+t = 1).

% So now we want to compute r,s,t,u, and also a,b,c,d for each pixel of the
% MS data.
% First, we intitialize some vectors:
MSpixelsInds=(1:size(MSvaluebox,1))'; % Pixel nbs, (1, 2, 3, etc.)
AlignedInds1=zeros(size(MSpixelsInds));  % a
AlignedInds2=zeros(size(MSpixelsInds));  % b
AlignedInds3=zeros(size(MSpixelsInds));  % c
AlignedInds4=zeros(size(MSpixelsInds));  % d
AlignedMSvaluebox1=zeros(size(MSvaluebox));  % a
AlignedMSintbox1=zeros(size(MSintbox));  % a
AlignedMSvaluebox2=zeros(size(MSvaluebox));  % b
AlignedMSintbox2=zeros(size(MSintbox));  % b
AlignedMSvaluebox3=zeros(size(MSvaluebox));  % c
AlignedMSintbox3=zeros(size(MSintbox));  % c
AlignedMSvaluebox4=zeros(size(MSvaluebox));  % d
AlignedMSintbox4=zeros(size(MSintbox));  % d
Interp_distr=zeros(size(MSpixelsInds));  % r
Interp_dists=zeros(size(MSpixelsInds));  % s
Interp_distt=zeros(size(MSpixelsInds));  % t
Interp_distu=zeros(size(MSpixelsInds));  % u
% Vectors for the "Deform_output":
    Defm1=zeros(size(MSpixelsInds));
    Defm2=zeros(size(MSpixelsInds));
    Defm3=zeros(size(MSpixelsInds));
    Defm4=zeros(size(MSpixelsInds));

% Convert displacement in 2-D chromatogram to displacement in 1-D vector:
AlignedInds1=MSpixelsInds; % Initialize with 0 displacement.
% AlignedInds will contain for each pixel in the aligned chromatogram its
% corresponding position in the target chromatogram.
AlignedInds2=MSpixelsInds; 
AlignedInds3=MSpixelsInds; 
AlignedInds4=MSpixelsInds; 

% Corner a (floor-floor)
% We will just need these 2 variables to refer to the good places in 2-D
% outputs from the non-MS alignment code):
FrstDFlag=1; % 1st D pixel count
ScndDFlag=1; % 2nd D pixel count
for ht=1:length(AlignedInds1)
%     We take position in the aligned chromatogram, and just add
%     displacements (add 2nd D displacement, and (1st D displacement *
%     number of pixels per modulation))
    AlignedInds1(ht)=AlignedInds1(ht)+floor(Displacement(ScndDFlag,FrstDFlag,1))+(NbPix2ndD*...
        floor(Displacement(ScndDFlag,FrstDFlag,2)));
%     Computes r,s,t,u:
    Interp_distr(ht)=Displacement(ScndDFlag,FrstDFlag,2)-...
        floor(Displacement(ScndDFlag,FrstDFlag,2));
    Interp_distt(ht)=Displacement(ScndDFlag,FrstDFlag,1)-...
        floor(Displacement(ScndDFlag,FrstDFlag,1));
     Interp_dists(ht)=-Displacement(ScndDFlag,FrstDFlag,2)+...
        ceil(Displacement(ScndDFlag,FrstDFlag,2));
    Interp_distu(ht)=-Displacement(ScndDFlag,FrstDFlag,1)+...
        ceil(Displacement(ScndDFlag,FrstDFlag,1));
% Extracts the local deformation correction stored in "Deform_output":
    Defm1(ht)=(Deform_output(ScndDFlag,FrstDFlag,1))*...
        (Deform_output(ScndDFlag,FrstDFlag,2))/4;
if ht~=NbPix2ndD*round(ht/NbPix2ndD) % If in a 2nd D column
ScndDFlag=ScndDFlag+1; % Go to next 2nd D pixel
else % i.e. if on the top of a 2nd D column, going to go to the next column at next iteration.
FrstDFlag=FrstDFlag+1; % Go to next column
ScndDFlag=1; % Go to first 2nd D pixel of the next column.
end
end

% Just correct r,s,t,u for cases when it is 0. (i.e. when it is on a border
% of the square "a,b,c,d" (see figure above)). In this case, just put it to
% 0.5 (it's because we will then just multiply values a,b,c,d by these
% numbers and make the sum to get the estimated value. That's why 0.5 seems
% the most indicated value in this case (0.5*0.5 = 0.25. 0.25, 4 times
% makes 1)).
Interp_distr(Interp_distr==0)=0.5;
Interp_dists(Interp_dists==0)=0.5;
Interp_distt(Interp_distt==0)=0.5;
Interp_distu(Interp_distu==0)=0.5;

% Pixels cannot come from outside the chromatogram. Label all of them with
% '0'
AlignedInds1(AlignedInds1>max(MSpixelsInds))=0;
AlignedInds1(AlignedInds1<0)=0;
sum(sum(AlignedInds1==0)) % Display the number of pixels "outside of chromatogram"
%  Create a vector to have the indices in the Aligned chromatogram, that
%  gives the position of the pixels that should come from inside the
%  chromatogram. These are the pixels that will be given some value during
%  alignments (others will remain at 0).
LpInds=1:length(AlignedInds1);
LpInds2=LpInds((AlignedInds1')~=0);
clear LpInds
for ht=LpInds2
%     Just keep the right MS ion m/z values:
    AlignedMSvaluebox1(ht,:)=MSvaluebox(AlignedInds1(ht),:);
%     Take the MS intensity value, but multiplied by s and u, and also
%     correct for local deformation (Deform_output, stored in Defm1).
        AlignedMSintbox1(ht,:)=MSintbox(AlignedInds1(ht),:)*Interp_dists(ht)*Interp_distu(ht).*Defm1(ht);
end
clear LpInds2 AlignedInds1

% Just do the same for corners b,c,d (same principle, refer to corner a for
% detailled comments).
% Corner b (ceil-floor)
FrstDFlag=1;
ScndDFlag=1;
for ht=1:length(AlignedInds2)
    AlignedInds2(ht)=AlignedInds2(ht)+floor(Displacement(ScndDFlag,FrstDFlag,1))+(NbPix2ndD*...
        ceil(Displacement(ScndDFlag,FrstDFlag,2)));
    Defm2(ht)=(Deform_output(ScndDFlag,FrstDFlag,1))*...
        (Deform_output(ScndDFlag,FrstDFlag,2))/4;
if ht~=NbPix2ndD*round(ht/NbPix2ndD)
ScndDFlag=ScndDFlag+1;
else
FrstDFlag=FrstDFlag+1;
ScndDFlag=1;
end
end
% Pixels cannot come from outside the chromatogram. Label all of them with
% '0'
AlignedInds2(AlignedInds2>max(MSpixelsInds))=0;
AlignedInds2(AlignedInds2<0)=0;
sum(sum(AlignedInds2==0))
LpInds=1:length(AlignedInds2);
LpInds2=LpInds((AlignedInds2')~=0);
clear LpInds
for ht=LpInds2
    AlignedMSvaluebox2(ht,:)=MSvaluebox(AlignedInds2(ht),:);
        AlignedMSintbox2(ht,:)=MSintbox(AlignedInds2(ht),:)*Interp_distr(ht)*Interp_distu(ht).*Defm2(ht);
end
clear LpInds2 AlignedInds2

% Corner c (floor-ceil)
FrstDFlag=1;
ScndDFlag=1;
for ht=1:length(AlignedInds3)
    AlignedInds3(ht)=AlignedInds3(ht)+ceil(Displacement(ScndDFlag,FrstDFlag,1))+(NbPix2ndD*...
        floor(Displacement(ScndDFlag,FrstDFlag,2)));
    Defm3(ht)=(Deform_output(ScndDFlag,FrstDFlag,1))*...
        (Deform_output(ScndDFlag,FrstDFlag,2))/4;
if ht~=NbPix2ndD*round(ht/NbPix2ndD)
ScndDFlag=ScndDFlag+1;
else
FrstDFlag=FrstDFlag+1;
ScndDFlag=1;
end
end
% Pixels cannot come from outside the chromatogram. Label all of them with
% '0'
AlignedInds3(AlignedInds3>max(MSpixelsInds))=0;
AlignedInds3(AlignedInds3<0)=0;
sum(sum(AlignedInds3==0))
LpInds=1:length(AlignedInds3);
LpInds2=LpInds((AlignedInds3')~=0);
clear LpInds
for ht=LpInds2
    AlignedMSvaluebox3(ht,:)=MSvaluebox(AlignedInds3(ht),:);
        AlignedMSintbox3(ht,:)=MSintbox(AlignedInds3(ht),:)*Interp_dists(ht)*Interp_distt(ht).*Defm3(ht);
end
clear LpInds2 AlignedInds3

% Corner d (ceil-ceil)
FrstDFlag=1;
ScndDFlag=1;
for ht=1:length(AlignedInds4)
    AlignedInds4(ht)=AlignedInds4(ht)+ceil(Displacement(ScndDFlag,FrstDFlag,1))+(NbPix2ndD*...
        ceil(Displacement(ScndDFlag,FrstDFlag,2)));
    Defm4(ht)=(Deform_output(ScndDFlag,FrstDFlag,1))*...
        (Deform_output(ScndDFlag,FrstDFlag,2))/4;
if ht~=NbPix2ndD*round(ht/NbPix2ndD)
ScndDFlag=ScndDFlag+1;
else
FrstDFlag=FrstDFlag+1;
ScndDFlag=1;
end
end
% Pixels cannot come from outside the chromatogram. Label all of them with
% '0'
AlignedInds4(AlignedInds4>max(MSpixelsInds))=0;
AlignedInds4(AlignedInds4<0)=0;
sum(sum(AlignedInds4==0))
LpInds=1:length(AlignedInds4);
LpInds2=LpInds((AlignedInds4')~=0);
clear LpInds
for ht=LpInds2
    AlignedMSvaluebox4(ht,:)=MSvaluebox(AlignedInds4(ht),:);
        AlignedMSintbox4(ht,:)=MSintbox(AlignedInds4(ht),:)*Interp_distr(ht)*Interp_distt(ht).*Defm4(ht);
end
clear LpInds2 AlignedInds4

disp('Computed 4 interpolated values')
toc

% First step, we will just put in matrices:
% [AlignedMSintbox1,AlignedMSintbox2,AlignedMSintbox3,AlignedMSintbox4].
% Here we initialize those vectors:
AlignedMSvalueboxI = zeros(size(MSvaluebox,1),size(MSvaluebox,2)*4);
AlignedMSintboxI = zeros(size(MSintbox,1),size(MSintbox,2)*4);

% Put the 4 corners interpolated values in the matrices.
AlignedMSvalueboxI=[AlignedMSvaluebox1,...
    AlignedMSvaluebox2,...
    AlignedMSvaluebox3,...
    AlignedMSvaluebox4];
AlignedMSintboxI=[AlignedMSintbox1,...
    AlignedMSintbox2,...
    AlignedMSintbox3,...
    AlignedMSintbox4];
% But, here we may have the same ion m/z several times for each pixel. 
% We want to have each m/z only once, and sum the corresponding intensity
% values.

% First step, remove useful zeros:
% Just temporarilly put the biggest value possible instead of the zeros:
AlignedMSvalueboxI(AlignedMSvalueboxI==0)=intmax;
% Sort the values (the big values corresponding to zeros go at the end)
% (sorted by ascending order):
[AlignedMSvalueboxI,IX]=sort(AlignedMSvalueboxI,2,'ascend');
% Re-put zeros instead of the big values:
AlignedMSvalueboxI(AlignedMSvalueboxI==intmax)=0;
% Just put the values in AlignedMSintboxI in the corresponding order:
for j = 1:size(AlignedMSintboxI,1)
    AlignedMSintboxI(j,:) = AlignedMSintboxI(j,IX(j,:)); 
end

% Just find the max of non-zero values (to remove useless zeros):
MaxNotZero=max(sum(AlignedMSintboxI~=0,2));
% 
% Remove useless zeros:
AlignedMSvalueboxII=AlignedMSvalueboxI(:,1:MaxNotZero);
AlignedMSintboxII=AlignedMSintboxI(:,1:MaxNotZero);
% size(AlignedMSvalueboxII)

% % % AlignedMSvaluebox=AlignedMSvalueboxII;
% % % AlignedMSintbox=AlignedMSintboxII;

% Clear no more useful variables from memory:
clear AlignedMSintboxI AlignedMSvalueboxI

% We will now just , for each pixel, only keep once each m/z value, and sum
% up the corresponding intensity values in case the m/z value was present
% more than once:
A = zeros(size(AlignedMSvalueboxII));
B = zeros(size(AlignedMSvalueboxII));
for kt=1:size(AlignedMSvalueboxII,1)
    B(kt,1) = AlignedMSvalueboxII(kt,1);
    A(kt,1) = AlignedMSintboxII(kt,1);
    Cnt = 1;
    for rr=1:size(AlignedMSvalueboxII,2)
        if AlignedMSvalueboxII(kt,rr) == B(kt,Cnt)
            A(kt,Cnt) = A(kt,Cnt) + AlignedMSintboxII(kt,rr);
        else
            Cnt = Cnt+1;
            B(kt,Cnt) = AlignedMSvalueboxII(kt,rr);
            A(kt,Cnt) = AlignedMSintboxII(kt,rr);
        end
        
    end

end

% We will again remove useless zeros:
MaxNotZero2 = max(sum(AlignedMSvalueboxII~=0,2));

MaxNotZero2 = max(sum(B~=0,2));
B = B(:,1:MaxNotZero2);
A = A(:,1:MaxNotZero2);
AlignedMSvaluebox = B;
AlignedMSintbox = A;

% disp('Data condensed')
% toc

Alignedeachscannum = size(AlignedMSintbox,2)*ones(size(AlignedMSintbox,1),1);

Alignedmedionid = cumsum(Alignedeachscannum);

disp('Alignment completed in (s):')
disp((now-t_initial)*24*3600)




