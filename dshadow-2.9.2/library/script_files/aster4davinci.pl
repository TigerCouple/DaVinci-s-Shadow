#!/usr/bin/perl 
#######################################################################
# 6/7/10 Added support for TIR SurfaceRadianceTIR and SkyIrradianceTIR
# 9/21/10 Added name/value pairs for projection parameters and UTM zones
#######################################################################

use strict;
use Getopt::Long;
use Math::Trig;

use POSIX qw(tmpnam);

use constant EARTH_FLATTENING => 1.0 - (6356.7523142/6378.1370); # (a-b/a) = 1 - b/a


# There are 15 bands in all: VNIR(1,2,3N,3B), SWIR(4,5,6,7,8,9), 
# TIR(10,11,12,13,14). Use the band ordering to extract only the desired bands
my %bandHash = ("1" => 1, "2" => 2, "3N" => 3, "3B"  => 4, "4"   => 5, "5"   => 6, "6" => 7, 
		"7" => 8, "8" => 9, "9"  => 10, "10" => 11, "11" => 12, "12" => 13, 
		"13" => 14, "14" => 15);

my $tmp;
my $bandPos;
my $lastSensor = "";
my $skipCnt = 0;
my $instrument;
my %keyWords;
my ($name, $value);
my $outputType;
my $band;
my ($year, $month, $day, $hour, $min, $sec);

my ($vnirPointingAngle, $swirPointingAngle, $tirPointingAngle, 
    $sceneCloudCoverage, $asterSceneId, $orbitNumber, $asterGDSDataGranuleId, 
    $solarAzimuth, $solarElevation, 
    $ulLat, $ulLon, $urLat, $urLon, $lrLat, $lrLon, $llLat, $llLon, 
    $cLat, $cLon, $ulLat_ac, $urLat_ac, $lrLat_ac, $llLat_ac, $cLat_ac, 
    $mapOrientAngle, $dataBrowseID, $band1Offset, $band2Offset, $band3NOffset, 
    $band3BOffset, $band4Offset, $band5Offset, $band6Offset, $band7Offset, $band8Offset, 
    $band9Offset, $band10Offset, $band11Offset, $band12Offset, $band13Offset, $band14Offset, 
    $band1CC, $band2CC, $band3NCC, $band3BCC, $band4CC, $band5CC, $band6CC, $band7CC, $band8CC,
    $band9CC, $band10CC, $band11CC, $band12CC, $band13CC, $band14CC, $utc, $band1_path, $band2_path,
    $band3N_path, $band3B_path, $band4_path, $band5_path, $band6_path, $band7_path, $band8_path, 
    $band9_path, $band10_path, $band11_path, $band12_path, $band13_path, $band14_path, 
    $bandSurfTemp_path,
    $day_night_flag, $param_name, $short_name, $data_plane_description, $band_scale_factors,
    $rti_start_line, $rti_stop_line, $rti_start_pixel, $rti_stop_pixel, $rti_aerosol_od_src,
    $rti_m_pres_src, $rti_molecular_od_src, $rti_junge_src, $rti_ssa_src, $rti_lut_vers,
    $rti_modtran_vers, $rti_m_o3_src, $rti_m_co2_src, $rti_m_h2o_src, $rti_m_temp_src,
    $rti_dem_src, $aerosol_od_src, $junge_src, $molecular_od_src, $aerosol_ssa_src,
    $modtran_o3_src, $modtran_o3_uncert, $modtran_h2o_src, $modtran_h2o_uncert,
    $modtran_temp_src, $modtran_temp_uncert, $modtran_pres_src, $modtran_pres_uncert,
    $dem_src, $dem_min_elevation, $dem_max_elevation, $dem_elevation_uncert, 
    $aster_map_projection, $map_projection_name, 
    $proj_param_band1, $proj_param_band2, $proj_param_band3N, $proj_param_band3B,
    $proj_param_band4, 
    $proj_param_band5, $proj_param_band6, $proj_param_band7, $proj_param_band8,
    $proj_param_band9, $proj_param_band10, $proj_param_band11, $proj_param_band12, 
    $proj_param_band13, $proj_param_band14,
    $utm_zone_band1, $utm_zone_band2, $utm_zone_band3N, $utm_zone_band3B, $utm_zone_band4, 
    $utm_zone_band5, $utm_zone_band6, $utm_zone_band7, $utm_zone_band8, 
    $utm_zone_band9, $utm_zone_band10, $utm_zone_band11, $utm_zone_band12, 
    $utm_zone_band13, $utm_zone_band14);

my $imgURL;
my $i;
my $destDir  = "";

my $gdalPath = "/usr/bin";

my $hdf      = "";
my $raw      = 0;
my $outFile;
my $cmdStr;
my $tStr;
my $tmpFile = POSIX::tmpnam();
my $gdalOutput = "gdalinfo.txt";
my ($imgCnt, $hdrCnt, $verbose, $help) = (0,0,0,0);

my $HDF_HDR;
my $GDAL_OUT;
my $ASTER_HDF;

my $usage   = "Usage: $0 [-v] [-h] [-gdal path] -d destDir  -s HDF\n\t" .
    "-v    verbose output\n\t" .
    "-h    this help\n\t" .
    "-d    destination path for output images\n\t" .
    "-gdal path to gdal utilities. Default: '/usr/bin'\n\t" .
    "-s    source path to ASTER HDF.\n";

my $success = GetOptions("-v"      => \$verbose, # Verbose output
			 "-gdal:s" => \$gdalPath,
			 "-d:s"    => \$destDir,
			 "-s:s"    => \$hdf,
                         "-h"      => \$help);

if(!$success || $help || $hdf eq "" || $destDir eq "") {
    print "$usage exiting...\n";
    exit(1);
}

if(!(-e $hdf)) {
    print "File not found\n";
    print "$usage exiting...\n";
    exit(1);
}

# Give the user a complete gdalinfo dump 
$cmdStr = "$gdalPath/gdalinfo -nogcp $hdf > $destDir/$gdalOutput";

system($cmdStr) == 0 
    or die "system $cmdStr failed: $?";

# Filter out everything but name = value pairs
open($HDF_HDR, "/bin/grep = $destDir/gdalinfo.txt |");

while(<$HDF_HDR>) {
    chomp($_);

    $_ =~ s/(\s*)\"//g; # Lose double quotes and trailing spaces
    $_ =~ s/^\s+//g;    # Lose leading spaces
    ($name, $value) = split("=", $_);		
			 
    if($name eq "SENSORNAME") {
	$lastSensor = $value;
    }
    elsif($name eq "ASTEROBSERVATIONMODE") {
	($instrument, $value) = split(", ", $value);		
	
	if($instrument eq "VNIR1") {
	    $name  = $instrument;
	}
	elsif($instrument eq "VNIR2") {
	    $name  = $instrument;
	}
	elsif($instrument eq "SWIR") {
	    $name  = $instrument;
	}
	elsif($instrument eq "TIR") {
	    $name  = $instrument;
	}
    }
    elsif($name eq "POINTINGANGLE") {
	if($lastSensor eq "VNIR") {
	    $name = "VNIR_POINTING_ANGLE";
	}
	elsif($lastSensor eq "SWIR") {
	    $name = "SWIR_POINTING_ANGLE";
	}
	elsif($lastSensor eq "TIR") {
	    $name = "TIR_POINTING_ANGLE";
	}
	else {
	    print STDERR "Error, illegal sensor name: $lastSensor\n";
	    exit(1);
	}
    }
    elsif($name eq "PLANE1_DESCRIPTION") {
	$name = "DATA_PLANE_DESCRIPTION";
    }

    $keyWords{$name} = $value;
}

close($HDF_HDR);

# Now scan just for the subdatasets
open($HDF_HDR, "/bin/grep -e SUBDATASET_'[0-9]*'_NAME $destDir/gdalinfo.txt |");

while(<$HDF_HDR>) {

    chomp($_);
    $_ =~ s/(\s*)\"//g; # Lose double quotes and trailing spaces"
    $_ =~ s/^\s+//g;    # Lose leading spaces

    ($name, $value) = split("=", $_);		
    $keyWords{$name} = $_;

    $bandPos = -1;

    $outputType = "";

    if(index($value, ":ImageData") > -1) {
	$band = substr($value, index($value, ":ImageData") + 10);

	$bandPos = $bandHash{$band};
    }

    elsif(index($value, ":SurfaceRadianceTIR:Band") > -1) {
	$band = substr($value, index($value, ":SurfaceRadianceTIR:Band") 
			+ length(":SurfaceRadianceTIR:Band"));

	# Radiance can be stored as a signed integer. Override so we can use PGM
	$outputType = " -ot UInt16 ";
	$bandPos = $bandHash{$band};
    }
#    elsif(index($value, ":SkyIrradianceTIR:Band") > -1) {
#	$band = substr($value, index($value, ":SkyIrradianceTIR:Band") 
#			+ length(":SkyIrradianceTIR:Band"));
#
#	# Radiance can be stored as a signed integer. Override so we can use PGM
#	$outputType = " -ot UInt16 ";
#	$bandPos = $bandHash{$band};
#    }
    elsif(index($value, ":SurfaceRadianceSWIR:Band") > -1) {
	$band = substr($value, index($value, ":SurfaceRadianceSWIR:Band") 
			+ length(":SurfaceRadianceSWIR:Band"));

	# Radiance can be stored as a signed integer. Override so we can use PGM
	$outputType = " -ot UInt16 ";
	$bandPos = $bandHash{$band};
    }
    elsif(index($value, ":SurfaceRadianceVNIR:Band") > -1) {
	$band = substr($value, index($value, ":SurfaceRadianceVNIR:Band") 
			+ length(":SurfaceRadianceVNIR:Band"));

	# Radiance can be stored as a signed integer. Override so we can use PGM
	$outputType = " -ot UInt16 ";
	$bandPos = $bandHash{$band};
    }
    elsif(index($value, ":SurfaceReflectanceVNIR:Band") > -1) {
	$band = substr($value, index($value, ":SurfaceReflectanceVNIR:Band") 
			+ length(":SurfaceReflectanceVNIR:Band"));

	# Refectance is stored as a signed integer. Override so we can use PGM
	$outputType = " -ot UInt16 ";
	$bandPos = $bandHash{$band};
    }
    elsif(index($value, ":SurfaceReflectanceSWIR:Band") > -1) {
	$band = substr($value, index($value, ":SurfaceReflectanceSWIR:Band") 
			+ length(":SurfaceReflectanceSWIR:Band"));

	# Refectance is stored as a signed integer. Override so we can use PGM
	$outputType = " -ot UInt16 ";
	$bandPos = $bandHash{$band};
    }
    elsif(index($value, ":SurfaceEmissivity:Band") > -1) {
	$band = substr($value, index($value, ":SurfaceEmissivity:Band") 
			+ length(":SurfaceEmissivity:Band"));

	$bandPos = $bandHash{$band};
    }
    elsif(index($value, ":SurfaceKineticTemperature:KineticTemperature") > -1) {
	$band    = "SurfTemp";
	$bandPos = 0;
    }

    if($bandPos != -1) { # Extract the data

	$outFile = "$destDir/band_" . $band . ".pgm"; # Destination of output

	$cmdStr = "$gdalPath/gdal_translate $outputType -of PNM  -sds \"$value\" $outFile > /dev/null";

	system($cmdStr) == 0 
	    or die "system $cmdStr failed: $?";
	    
	$keyWords{"band$band" . "_path"} = $outFile;
    }
    else {
                                #commented to not show annoying messages
#	print STDERR "No subdataset matches $value. Skipping...\n\n";
    }
}

close($HDF_HDR);

# Transfer parameters out of the hash
$sceneCloudCoverage = $keyWords{"SCENECLOUDCOVERAGE"};
$vnirPointingAngle  = $keyWords{"VNIR_POINTING_ANGLE"};	
$swirPointingAngle  = $keyWords{"SWIR_POINTING_ANGLE"};
$tirPointingAngle   = $keyWords{"TIR_POINTING_ANGLE"};
$asterSceneId       = $keyWords{"ASTERSCENEID"};
$orbitNumber        = $keyWords{"ORBITNUMBER"};

$asterGDSDataGranuleId = $keyWords{"IDOFASTERGDSDATAGRANULE"};

$tmp = $keyWords{"SOLARDIRECTION"};
($solarAzimuth, $solarElevation) = split(", ", $tmp);			    

$tmp = $keyWords{"UPPERLEFT"};
($ulLat, $ulLon) = split(", ", $tmp);

$tmp = $keyWords{"UPPERRIGHT"};
($urLat, $urLon) = split(", ", $tmp);

$tmp = $keyWords{"LOWERLEFT"};
($llLat, $llLon) = split(", ", $tmp);

$tmp = $keyWords{"LOWERRIGHT"};
($lrLat, $lrLon) = split(", ", $tmp);

$tmp = $keyWords{"SCENECENTER"};
($cLat, $cLon) = split(", ", $tmp);

# Flip from -180, +180 to 0 to 360
if($ulLon < 0) { $ulLon += 360.0; }
if($llLon < 0) { $llLon += 360.0; }
if($urLon < 0) { $urLon += 360.0; }
if($lrLon < 0) { $lrLon += 360.0; }

# Covert to Aerocentric latitude
$ulLat_ac = &ographic2ocentric($ulLat);
$urLat_ac = &ographic2ocentric($urLat);
$llLat_ac = &ographic2ocentric($llLat);
$lrLat_ac = &ographic2ocentric($lrLat);
$cLat_ac  = &ographic2ocentric($cLat);

$mapOrientAngle = $keyWords{"MAPORIENTATIONANGLE"};
$dataBrowseID   = $keyWords{"IDOFASTERGDSDATABROWSE"};

# Extract band offset
$band1Offset  = &extractField("OFFSET1", %keyWords, "");
$band2Offset  = &extractField("OFFSET2", %keyWords, "");
$band3NOffset = &extractField("OFFSET3N", %keyWords, "");
$band3BOffset = &extractField("OFFSET3B", %keyWords, "");
$band4Offset  = &extractField("OFFSET4", %keyWords, "");
$band5Offset  = &extractField("OFFSET5", %keyWords, "");
$band6Offset  = &extractField("OFFSET6", %keyWords, "");
$band7Offset  = &extractField("OFFSET7", %keyWords, "");
$band8Offset  = &extractField("OFFSET8", %keyWords, "");
$band9Offset  = &extractField("OFFSET9", %keyWords, "");
$band10Offset = &extractField("OFFSET10", %keyWords, "");
$band11Offset = &extractField("OFFSET11", %keyWords, "");
$band12Offset = &extractField("OFFSET12", %keyWords, "");
$band13Offset = &extractField("OFFSET13", %keyWords, "");
$band14Offset = &extractField("OFFSET14", %keyWords, "");

# Extract Correction Coefficients
$band1CC  = &extractField("INCL1",  %keyWords, "");
$band2CC  = &extractField("INCL2",  %keyWords, "");
$band3NCC = &extractField("INCL3N", %keyWords, "");
$band3BCC = &extractField("INCL3B", %keyWords, "");
$band4CC  = &extractField("INCL4",  %keyWords, "");
$band5CC  = &extractField("INCL5",  %keyWords, "");
$band6CC  = &extractField("INCL6",  %keyWords, "");
$band7CC  = &extractField("INCL7",  %keyWords, "");
$band8CC  = &extractField("INCL8",  %keyWords, "");
$band9CC  = &extractField("INCL9",  %keyWords, "");
$band10CC = &extractField("INCL10", %keyWords, "");
$band11CC = &extractField("INCL11", %keyWords, "");
$band12CC = &extractField("INCL12", %keyWords, "");
$band13CC = &extractField("INCL13", %keyWords, "");
$band14CC = &extractField("INCL14", %keyWords, "");

$rti_start_line       = &extractField("RTI_START_LINE", %keyWords, "");
$rti_stop_line        = &extractField("RTI_STOP_LINE", %keyWords, "");
$rti_start_pixel      = &extractField("RTI_START_PIXEL", %keyWords, "");
$rti_stop_pixel       = &extractField("RTI_STOP_PIXEL", %keyWords, "");
$rti_aerosol_od_src   = &extractField("RTI_AEROSOL_OD_SRC", %keyWords, "");
$rti_m_pres_src       = &extractField("RTI_M_PRES_SRC", %keyWords, "");
$rti_molecular_od_src = &extractField("RTI_MOLECULAR_OD_SRC", %keyWords, "");
$rti_junge_src        = &extractField("RTI_JUNGE_SRC", %keyWords, "");
$rti_ssa_src          = &extractField("RTI_SSA_SRC", %keyWords, "");
$rti_lut_vers         = &extractField("RTI_LUT_VERS", %keyWords, "");
$rti_modtran_vers     = &extractField("RTI_MODTRAN_VERS", %keyWords, "");
$rti_m_o3_src         = &extractField("RTI_M_O3_SRC", %keyWords, "");
$rti_m_co2_src        = &extractField("RTI_M_CO2_SRC", %keyWords, "");
$rti_m_h2o_src        = &extractField("RTI_M_H2O_SRC", %keyWords, "");
$rti_m_temp_src       = &extractField("RTI_M_TEMP_SRC", %keyWords, "");
$rti_dem_src          = &extractField("RTI_DEM_SRC", %keyWords, "");
$aerosol_od_src       = &extractField("AEROSOL_OD_SRC", %keyWords, "");
$junge_src            = &extractField("JUNGE_SRC", %keyWords, "");
$molecular_od_src     = &extractField("MOLECULAR_OD_SRC", %keyWords, "");
$aerosol_ssa_src      = &extractField("AEROSOL_SSA_SRC", %keyWords, "");
$modtran_o3_src       = &extractField("MODTRAN_O3_SRC", %keyWords, "");
$modtran_o3_uncert    = &extractField("MODTRAN_O3_UNCERT", %keyWords, "");
$modtran_h2o_src      = &extractField("MODTRAN_H2O_SRC", %keyWords, "");
$modtran_h2o_uncert   = &extractField("MODTRAN_H2O_UNCERT", %keyWords, "");
$modtran_temp_src     = &extractField("MODTRAN_TEMP_SRC", %keyWords, "");
$modtran_temp_uncert  = &extractField("MODTRAN_TEMP_UNCERT", %keyWords, "");
$modtran_pres_src     = &extractField("MODTRAN_PRES_SRC", %keyWords, "");
$modtran_pres_uncert  = &extractField("MODTRAN_PRES_UNCERT", %keyWords, "");
$dem_src              = &extractField("DEM_SRC", %keyWords, "");
$dem_min_elevation    = &extractField("DEM_MIN_ELEVATION", %keyWords, "");
$dem_max_elevation    = &extractField("DEM_MAX_ELEVATION", %keyWords, "");
$dem_elevation_uncert = &extractField("DEM_ELEVATION_UNCERT", %keyWords, "");

$aster_map_projection  = &extractField("ASTERMapProjection", %keyWords, "");
$map_projection_name   = &extractField("MAPPROJECTIONNAME", %keyWords, "");

$proj_param_band1   = &extractField("PROJECTIONPARAMETERS1", %keyWords, "");
$proj_param_band2   = &extractField("PROJECTIONPARAMETERS2", %keyWords, "");
$proj_param_band3N  = &extractField("PROJECTIONPARAMETERS3N", %keyWords, "");
$proj_param_band3B  = &extractField("PROJECTIONPARAMETERS3B", %keyWords, "");
$proj_param_band4   = &extractField("PROJECTIONPARAMETERS4", %keyWords, "");
$proj_param_band5   = &extractField("PROJECTIONPARAMETERS5", %keyWords, "");
$proj_param_band6   = &extractField("PROJECTIONPARAMETERS6", %keyWords, "");
$proj_param_band7   = &extractField("PROJECTIONPARAMETERS7", %keyWords, "");
$proj_param_band8   = &extractField("PROJECTIONPARAMETERS8", %keyWords, "");
$proj_param_band9   = &extractField("PROJECTIONPARAMETERS9", %keyWords, "");
$proj_param_band10  = &extractField("PROJECTIONPARAMETERS10", %keyWords, "");
$proj_param_band11  = &extractField("PROJECTIONPARAMETERS11", %keyWords, "");
$proj_param_band12  = &extractField("PROJECTIONPARAMETERS12", %keyWords, "");
$proj_param_band13  = &extractField("PROJECTIONPARAMETERS13", %keyWords, "");
$proj_param_band14  = &extractField("PROJECTIONPARAMETERS14", %keyWords, "");

$utm_zone_band1  = &extractField("UTMZONECODE1", %keyWords, "");
$utm_zone_band2  = &extractField("UTMZONECODE2", %keyWords, "");
$utm_zone_band3N = &extractField("UTMZONECODE3N", %keyWords, "");
$utm_zone_band3B = &extractField("UTMZONECODE3B", %keyWords, "");
$utm_zone_band4  = &extractField("UTMZONECODE4", %keyWords, "");
$utm_zone_band5  = &extractField("UTMZONECODE5", %keyWords, "");
$utm_zone_band6  = &extractField("UTMZONECODE6", %keyWords, "");
$utm_zone_band7  = &extractField("UTMZONECODE7", %keyWords, "");
$utm_zone_band8  = &extractField("UTMZONECODE8", %keyWords, "");
$utm_zone_band9  = &extractField("UTMZONECODE9", %keyWords, "");
$utm_zone_band10 = &extractField("UTMZONECODE10", %keyWords, "");
$utm_zone_band11 = &extractField("UTMZONECODE11", %keyWords, "");
$utm_zone_band12 = &extractField("UTMZONECODE12", %keyWords, "");
$utm_zone_band13 = &extractField("UTMZONECODE13", %keyWords, "");
$utm_zone_band14 = &extractField("UTMZONECODE14", %keyWords, "");

# Extract image paths
$band1_path   = &extractField("band1_path", %keyWords, "");
$band2_path   = &extractField("band2_path", %keyWords, "");
$band3N_path  = &extractField("band3N_path", %keyWords, "");
$band3B_path  = &extractField("band3B_path", %keyWords, "");
$band4_path   = &extractField("band4_path", %keyWords, "");
$band5_path   = &extractField("band5_path", %keyWords, "");
$band6_path   = &extractField("band6_path", %keyWords, "");
$band7_path   = &extractField("band7_path", %keyWords, "");
$band8_path   = &extractField("band8_path", %keyWords, "");
$band9_path   = &extractField("band9_path", %keyWords, "");
$band10_path  = &extractField("band10_path", %keyWords, "");
$band11_path  = &extractField("band11_path", %keyWords, "");
$band12_path  = &extractField("band12_path", %keyWords, "");
$band13_path  = &extractField("band13_path", %keyWords, "");
$band14_path  = &extractField("band14_path", %keyWords, "");
$bandSurfTemp_path = &extractField("bandSurfTemp_path", %keyWords, "");

$day_night_flag = &extractField("DAYNIGHTFLAG", %keyWords, "");
$param_name     = &extractField("PARAMETERNAME", %keyWords, "");
$short_name     = &extractField("SHORTNAME", %keyWords, "");

$band_scale_factors     = &extractField("BANDSCALEFACTORS", %keyWords, "");
$data_plane_description = &extractField("DATA_PLANE_DESCRIPTION", %keyWords, "");

my @tArray = ($vnirPointingAngle, $swirPointingAngle, $tirPointingAngle, 
	      $sceneCloudCoverage, $asterSceneId, 
	      $orbitNumber, $asterGDSDataGranuleId, $solarAzimuth, $solarElevation, 
	      $ulLat, $ulLon, $urLat, $urLon, $lrLat, $lrLon, $llLat, $llLon, 
	      $cLat, $cLon, $ulLat_ac, $urLat_ac, $lrLat_ac, $llLat_ac, $cLat_ac, 
	      $mapOrientAngle, $dataBrowseID, $band1Offset, $band2Offset, $band3NOffset, 
	      $band3BOffset, $band4Offset, 
	      $band5Offset, $band6Offset, $band7Offset, $band8Offset, $band9Offset, $band10Offset, 
	      $band11Offset, $band12Offset, $band13Offset, $band14Offset, $band1CC, $band2CC, 
	      $band3NCC, $band3BCC, $band4CC, $band5CC, $band6CC, $band7CC, $band8CC, $band9CC, 
	      $band10CC, $band11CC, $band12CC, $band13CC, $band14CC,
	      $band1_path, $band2_path, $band3N_path, $band3B_path, $band4_path, $band5_path,
	      $band6_path, $band7_path, $band8_path, $band9_path, $band10_path, $band11_path,
	      $band12_path, $band13_path, $band14_path,
	      $bandSurfTemp_path, $day_night_flag, $param_name, $short_name, $data_plane_description,
	      $rti_start_line, $rti_stop_line, $rti_start_pixel, $rti_stop_pixel, 
	      $rti_aerosol_od_src, $rti_m_pres_src,$rti_molecular_od_src, $rti_junge_src, $rti_ssa_src, 
	      $rti_lut_vers, $rti_modtran_vers, $rti_m_o3_src, $rti_m_co2_src, $rti_m_h2o_src, 
	      $rti_m_temp_src, $rti_dem_src, $aerosol_od_src, $junge_src, $molecular_od_src, 
	      $aerosol_ssa_src, $modtran_o3_src, $modtran_o3_uncert, $modtran_h2o_src, 
	      $modtran_h2o_uncert, $modtran_temp_src, $modtran_temp_uncert, $modtran_pres_src, 
	      $modtran_pres_uncert, $dem_src, $dem_min_elevation, $dem_max_elevation, $dem_elevation_uncert,
	      $aster_map_projection, $map_projection_name, $proj_param_band1, $proj_param_band2, 
	      $proj_param_band3N, $proj_param_band3B, $proj_param_band4, $proj_param_band5, 
	      $proj_param_band6, $proj_param_band7, $proj_param_band8, $proj_param_band9, $proj_param_band10,
	      $proj_param_band11, $proj_param_band12, $proj_param_band13, $proj_param_band14,
	      $utm_zone_band1, $utm_zone_band2, $utm_zone_band3N, $utm_zone_band3B, 
	      $utm_zone_band4, $utm_zone_band5, $utm_zone_band6, $utm_zone_band7, 
	      $utm_zone_band8, $utm_zone_band9, $utm_zone_band10, $utm_zone_band11, 
	      $utm_zone_band12, $utm_zone_band13, $utm_zone_band14);

my @asterHdrs = qw(vnir_pointing_angle swir_pointing_angle tir_pointing_angle 
		   scene_cloud_coverage aster_sceneid orbit_number aster_gds_granule_id 
		   solar_azimuth solar_elevation ul_lat ul_lon ur_lat ur_lon lr_lat 
		   lr_lon ll_lat ll_lon clat clon ul_lat_ac ur_lat_ac lr_lat_ac 
		   ll_lat_ac clat_ac map_orient_angle data_browse_id 
		   band1_offset band2_offset band3n_offset band3b_offset band4_offset 
		   band5_offset band6_offset band7_offset band8_offset band9_offset 
		   band10_offset band11_offset band12_offset band13_offset band14_offset
		   band1_corr_coeff band2_corr_coeff band3n_corr_coeff band3b_corr_coeff 
		   band4_corr_coeff band5_corr_coeff band6_corr_coeff band7_corr_coeff 
		   band8_corr_coeff band9_corr_coeff band10_corr_coeff band11_corr_coeff 
		   band12_corr_coeff band13_corr_coeff band14_corr_coeff 
		   band1_path band2_path band3N_path band3B_path band4_path band5_path
		   band6_path band7_path band8_path band9_path band10_path band11_path
		   band12_path band13_path band14_path bandSurfTemp_path day_night_flag 
		   param_name short_name data_plane_description 
		   rti_start_line rti_stop_line rti_start_pixel rti_stop_pixel 
		   rti_aerosol_od_src rti_m_pres_src rti_molecular_od_src rti_junge_src 
		   rti_ssa_src rti_lut_vers rti_modtran_vers rti_m_o3_src rti_m_co2_src 
		   rti_m_h2o_src rti_m_temp_src rti_dem_src aerosol_od_src junge_src 
		   molecular_od_src aerosol_ssa_src modtran_o3_src modtran_o3_uncert 
		   modtran_h2o_src modtran_h2o_uncert modtran_temp_src modtran_temp_uncert 
		   modtran_pres_src modtran_pres_uncert dem_src dem_min_elevation 
		   dem_max_elevation dem_elevation_uncert aster_map_projection 
		   map_projection_name
		   proj_param_band1 proj_param_band2 proj_param_band3N proj_param_band3B proj_param_band4 proj_param_band5 
		   proj_param_band6 proj_param_band7 proj_param_band8 proj_param_band9 proj_param_band10
		   proj_param_band11 proj_param_band12 proj_param_band13 proj_param_band14
		   utm_zone_band1 utm_zone_band2 utm_zone_band3N utm_zone_band3B 
		   utm_zone_band4 utm_zone_band5 utm_zone_band6 utm_zone_band7 
		   utm_zone_band8 utm_zone_band9 utm_zone_band10 utm_zone_band11 
		   utm_zone_band12 utm_zone_band13 utm_zone_band14);

for($i = 0; $i < $#tArray; $i++) {
    if($tArray[$i] ne "") { print "$asterHdrs[$i]\t$tArray[$i]\n"; }
}

# Dump the scale factors, if any
if($band_scale_factors ne "") { 

    @tArray = split(",", $band_scale_factors);

    for($i = 0; $i < 15; $i++) {
	if($i < 2) {
	    $band = $i + 1;
	}
	elsif($i == 2) {
	    $band = "3N";
	}
	elsif($i == 3) {
	    $band = "3B";
	}
	else {
	    $band = $i;
	}

	print "scale_factor_$band\t$tArray[$i]\n"; 
    }
}

print "gdalinfo_file\t$destDir/$gdalOutput\n";

exit(0);


sub extractField {
    my ($key, %keyWords, $default) = @_;

    if(exists $keyWords{"$key"}) {
	return $keyWords{"$key"};
    }
    else {
	return $default;
    }
}

sub ographic2ocentric {
    my ($lat) = @_;

    my $g2c = (1.0 - EARTH_FLATTENING) * (1.0 - EARTH_FLATTENING);

    return rad2deg(atan($g2c * tan(deg2rad($lat))));
}
