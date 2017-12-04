function [Chromato] = matlab_cdf_open_function(FileName, Intthreshold,driftMS)

% Function to open a cdf file containing a MS chromatogram. The output is
% a structure.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% See license terms stated in LICENSE.txt
% Authors : Jonas Gros, Yasuyuki Zushi, and J. Samuel Arey.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

cdf.dat1 = netcdf.open(FileName,'NOWRITE');

varidflag = netcdf.inqVarID(cdf.dat1,'flag_count');
varidscantime = netcdf.inqVarID(cdf.dat1,'scan_acquisition_time');
varidscannum = netcdf.inqVarID(cdf.dat1,'actual_scan_number');
varidmedmsmax = netcdf.inqVarID(cdf.dat1,'mass_range_max');
varidmedmsmin = netcdf.inqVarID(cdf.dat1,'mass_range_min');
varidionid = netcdf.inqVarID(cdf.dat1,'scan_index');
varideachscannum = netcdf.inqVarID(cdf.dat1,'point_count');
varidMStotint = netcdf.inqVarID(cdf.dat1,'total_intensity');
varidMSvalue = netcdf.inqVarID(cdf.dat1,'mass_values');
varidMSint = netcdf.inqVarID(cdf.dat1,'intensity_values');

flag = netcdf.getVar(cdf.dat1,varidflag);
scantime = netcdf.getVar(cdf.dat1,varidscantime);
scannum = netcdf.getVar(cdf.dat1,varidscannum);
medmsmax = netcdf.getVar(cdf.dat1,varidmedmsmax);
medmsmin = netcdf.getVar(cdf.dat1,varidmedmsmin);
ionid = netcdf.getVar(cdf.dat1,varidionid);
eachscannum = netcdf.getVar(cdf.dat1,varideachscannum);
MStotint = netcdf.getVar(cdf.dat1,varidMStotint);
MSvalue = netcdf.getVar(cdf.dat1,varidMSvalue);
MSint = netcdf.getVar(cdf.dat1,varidMSint);

MSvalue = MSvalue + driftMS;



Timepara = scantime(abs(eachscannum)<intmax );
RTini = min( Timepara )/60;
RTruntime = max( Timepara )/60-RTini;
SamRate = 1/(mean( Timepara(2:length(Timepara))-Timepara(1: (length(Timepara)-1) ) ));



pixelnum = length(scannum);
maxscannum = max(eachscannum);
minscannum = min(eachscannum);
MSdatabox = zeros(pixelnum, maxscannum+1 );
MSvaluebox = MSdatabox;
MSintbox = MSdatabox;

for inum=1:pixelnum
if (abs(eachscannum(inum))<intmax)
    initial = sum(eachscannum(1:inum))-eachscannum(inum)+1;
    acqrange = min(initial,(initial+eachscannum(inum)-1)):max(initial,(initial+eachscannum(inum)-1));
if(eachscannum(inum)==0)
    acqrange =[];
end
remainrep = zeros(1,maxscannum-eachscannum(inum)+1 );
MSvaluebox(inum,:) = [MSvalue( acqrange )', remainrep ];
MSintbox(inum,:) = [MSint( acqrange )',remainrep];
end
end

MSvaluebox(MSintbox < Intthreshold) =0;
MSintbox(MSintbox < Intthreshold)=0;


Chromato.flag = flag;
Chromato.scantime = scantime;
Chromato.scannum = scannum;
Chromato.medmsmax = medmsmax;
Chromato.medmsmin = medmsmin;
Chromato.eachscannum = eachscannum;
Chromato.MStotint = MStotint;
Chromato.MSvalue = MSvalue;
Chromato.MSint = MSint;
Chromato.Timepara = Timepara;
Chromato.RTini = RTini;
Chromato.RTruntime = RTruntime;
Chromato.SamRate = SamRate;
%Chromato.threshold = Intthreshold;
Chromato.MSthrethold = Intthreshold;
Chromato.pixelnum = pixelnum;
Chromato.maxscannum = maxscannum;
Chromato.minscannum = minscannum;
Chromato.MSdatabox = MSdatabox;
Chromato.MSvaluebox = MSvaluebox;
Chromato.MSintbox = MSintbox;
Chromato.ionid = ionid;
% Chromato.remainrep = remainrep;




