function  [] = cdf_write_function(FileName,Alignedscantime,AlignedMSintboxLine,Alignedflag,Alignedscannum,Alignedmedmsmax,Alignedmedmsmin,Alignedionid,Alignedeachscannum,AlignedMStotint,AlignedMSvalueboxLine)

% Function to write a MS chromatogram to a cdf file.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% See license terms stated in LICENSE.txt
% Authors : Jonas Gros, Yasuyuki Zushi, and J. Samuel Arey.
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%

ncnew = netcdf.create( FileName,'64BIT_OFFSET'  );



scan_number = netcdf.defDim(ncnew, 'scan_number', length(Alignedscantime) ) ;
point_number = netcdf.defDim(ncnew, 'point_number',length(AlignedMSintboxLine) );


vardimfirst = netcdf.defVar(ncnew,  'flag_count', 'NC_INT',scan_number );
vardimA = netcdf.defVar(ncnew, 'scan_acquisition_time', 'NC_DOUBLE',scan_number );
vardimB = netcdf.defVar(ncnew, 'actual_scan_number','NC_INT',scan_number );
vardimC = netcdf.defVar(ncnew, 'mass_range_max', 'NC_DOUBLE',scan_number );
vardimD = netcdf.defVar(ncnew, 'mass_range_min', 'NC_DOUBLE',scan_number );
vardimE = netcdf.defVar(ncnew, 'scan_index', 'NC_INT',scan_number );
vardimF = netcdf.defVar(ncnew, 'point_count', 'NC_INT',scan_number );
vardimG = netcdf.defVar(ncnew, 'total_intensity', 'NC_DOUBLE',scan_number );
vardimH = netcdf.defVar(ncnew, 'mass_values','NC_DOUBLE', point_number );
vardimI = netcdf.defVar(ncnew, 'intensity_values', 'NC_DOUBLE', point_number );

netcdf.endDef(ncnew)


netcdf.putVar(ncnew, vardimfirst, Alignedflag )
netcdf.putVar(ncnew, vardimA, Alignedscantime )
netcdf.putVar(ncnew, vardimB, Alignedscannum )
netcdf.putVar(ncnew, vardimC, Alignedmedmsmax )
netcdf.putVar(ncnew, vardimD, Alignedmedmsmin )
netcdf.putVar(ncnew, vardimE, Alignedionid )
netcdf.putVar(ncnew, vardimF, Alignedeachscannum )
netcdf.putVar(ncnew, vardimG, AlignedMStotint )
netcdf.putVar(ncnew, vardimH, AlignedMSvalueboxLine )
netcdf.putVar(ncnew, vardimI, AlignedMSintboxLine )


netcdf.close(ncnew)



