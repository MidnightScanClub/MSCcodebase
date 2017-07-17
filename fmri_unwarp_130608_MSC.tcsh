#!/bin/tcsh -f

set program = $0:t
if ( $?FSLDIR == 0 ) then
	echo "Error: FSLDIR environment variable needs to be set"
	exit 1
endif

if ( $?RELEASE == 0 ) then
	echo "Error: RELEASE environment variable needs to be set"
	exit 1
endif
if ( $#argv < 5 ) goto USAGE
##############################################################
# extract first character to insure first argument is a switch
##############################################################
if ( `echo $argv[1] | head -c1` != "-" ) then
	echo $program "error: First argument is not a method switch"; goto USAGE
endif
####################
# get method and run
####################
set method = `echo $argv[1] | sed 's/^.//'`; shift
if ( $method == "map" ) then
###############
# check options
###############
	echo $#argv
	if ( $#argv != 8) then
		echo $program "error: Incorrect number of options for own field map"; goto USAGE
	endif
###################
# process variables
###################
	set patid = $1;
	set epi = $2;
	if ($epi:e == "img") set epi = $epi:r;  if ($epi:e == "4dfp") set epi = $epi:r;
	if (! -e $epi.4dfp.img || ! -e $epi.4dfp.ifh) then
		echo $epi not found
		exit -1
	endif
	set mag = $3; set phase = $4;
	set dwell = $5; set te = $6; set ped = $7; set delta = $8;
#######################
# make output directory
#######################
	if (! -d unwarp) mkdir unwarp
####################
# convert epi to nii
####################
	nifti_4dfp -n $epi unwarp/$epi:t
################
# run epi_unwarp
################
	#epi_unwarp -m $mag -p $phase -e unwarp/$epi:t -dwell $dwell -te $te -delta $delta -r $patid -dir $ped -nomask 
	/data/nil-bluearc/GMT/Laumann/MSC/Scripts/epi_unwarp_MSC -m $mag -p $phase -e unwarp/$epi:t -dwell $dwell -te $te -delta $delta -r $patid -dir $ped -nomask
	if ($status) exit $status
#####################################################################
# create link so we have a shift image that matches the other methods
#####################################################################
	ln -s $cwd/unwarp/${patid}_epi_0_unwarp_shift_warp.nii.gz $cwd/unwarp/${epi:t}_uwrp_shift_warp.nii.gz
#############################
# convert result back to 4dfp
#############################
	gunzip -f     unwarp/${patid}_epi_unwarped.nii.gz
	nifti_4dfp -4 unwarp/${patid}_epi_unwarped unwarp/${epi:t}_uwrp
	exit $status

else if ( $method == "mean" ) then
	
	#Check options
	if ( $#argv != 5 ) then
		echo "Error: Incorrect number of options for mean field map. See usage."
		exit 1
	endif
	
	#Make variable names
	set epi = $1:r:r; set mean = $2:r:r; set mean_t4 = $3; set dwell = `echo "$4/1000" | bc -l`; set ped = $5
	
	#Make output folder
	if ( ! -d unwarp ) mkdir unwarp
	
	#Transform the mean field map to the epi space
	$RELEASE/t4_inv $mean_t4 unwarp/${mean:t}_to_${epi:t}_t4
	$RELEASE/t4img_4dfp unwarp/${mean:t}_to_${epi:t}_t4 $mean unwarp/${mean:t}_on_${epi:t} -O$epi

	#Convert mean and epi to nii
	foreach image ( unwarp/${mean:t}_on_${epi:t} $epi )
		$RELEASE/nifti_4dfp -n $image unwarp/$image:t
	end

	#Undistort EPI with new field map
	$FSLDIR/bin/fugue --loadfmap=unwarp/${mean:t}_on_${epi:t} --dwell=$dwell --in=unwarp/$epi:t -u unwarp/${epi:t}_uwrp --unwarpdir=$ped --saveshift=unwarp/${epi:t}_uwrp_shift_map
	
	#Convert the shiftwarp to an absolute warp
	$FSLDIR/bin/convertwarp -s unwarp/${epi:t}_uwrp_shift_map -r unwarp/$epi:t --shiftdir=$ped -o unwarp/${epi:t}_uwrp_shift_warp

	#Convert the result back to 4dfp land
	gunzip -f unwarp/${epi:t}_uwrp.nii.gz
	$RELEASE/nifti_4dfp -4 unwarp/${epi:t}_uwrp unwarp/${epi:t}_uwrp 
	exit $status

else if ( $method == "basis" ) then
###############
# check options
###############
	if ( $#argv != 9 && $#argv != 10 ) then
		echo $program "error: Incorrect number of options for basis functions"; goto USAGE
	endif
###################
# recover variables
###################
	set epi = $1:r:r; set t2w = $2:r:r; set mean = $3:r:r; set basis = $4:r:r
	set t2w_t4 = $5; set mean_t4 = $6; set dwell = $7; set ped = $8; set nbasis = $9
	
#######################
# make output directory
#######################
	if ( ! -d unwarp ) mkdir unwarp
#################
# create t2w mask
#################
	if ( $#argv == 10 ) then
		set t2w_mask = $10:r:r
	else
		nifti_4dfp -n $t2w unwarp/$t2w:t
		$FSLDIR/bin/bet unwarp/$t2w:t unwarp/${t2w:t}_brain -m -n -f 0.4 -R
		gunzip -f unwarp/${t2w:t}_brain_mask.nii.gz
		nifti_4dfp -4 unwarp/${t2w:t}_brain_mask unwarp/${t2w:t}_brain_mask 
		set t2w_mask = unwarp/${t2w:t}_brain_mask
	endif
	
	t4_inv $mean_t4 unwarp/${mean:t}_to_${epi:t}_t4
#################
# resample FMmean
#################
	t4img_4dfp unwarp/${mean:t}_to_${epi:t}_t4 $mean unwarp/${mean:t}_on_${epi:t} -O$epi
	if ($status) exit $status
#####################################################
# resample FMbases in epi space (only needed volumes)
#####################################################
	ln -s		$basis.4dfp.img			unwarp/${basis:t}.4dfp.img
	ln -s		$basis.4dfp.img.rec		unwarp/${basis:t}.4dfp.img.rec
	/bin/cp	-f	$basis.4dfp.ifh			unwarp/${basis:t}.4dfp.ifh
	echo "matrix size [4] := "$nbasis >>		unwarp/${basis:t}.4dfp.ifh
	t4img_4dfp unwarp/${mean:t}_to_${epi:t}_t4	unwarp/${basis:t} unwarp/${basis:t}_on_${epi:t} -O$epi
	if ($status) exit $status
####################
# convert epi to nii
####################
	nifti_4dfp -n $epi unwarp/$epi:t
	
###################################
# make parameter file for basis_opt
###################################
	echo "set epi = unwarp/$epi:t"				>! unwarp/${epi:t}_basis_opt.params
	echo "set basis = unwarp/${basis:t}_on_${epi:t}"	>> unwarp/${epi:t}_basis_opt.params
	echo "set mean = unwarp/${mean:t}_on_${epi:t}"		>> unwarp/${epi:t}_basis_opt.params
	echo "set dir = $ped"					>> unwarp/${epi:t}_basis_opt.params
	echo "set phase = unwarp/${epi:t}_basis_opt_phase"	>> unwarp/${epi:t}_basis_opt.params
	echo "set t2 = $t2w"					>> unwarp/${epi:t}_basis_opt.params
	echo "set t2_mskt = $t2w_mask"				>> unwarp/${epi:t}_basis_opt.params
	echo "set t4 = $t2w_t4"					>> unwarp/${epi:t}_basis_opt.params
	echo "set unwarp_Tyler = /data/nil-bluearc/benzinger2/Tyler/test/unwarp_eta.tcsh" >> unwarp/${epi:t}_basis_opt.params
##################################
# make a weight file for basis_opt
##################################
	echo "Echo Spacing: $dwell"				>! unwarp/${epi:t}_basis_opt.weights
	echo "Weights: $nbasis"					>> unwarp/${epi:t}_basis_opt.weights
###############################################
# run basis_opt with created params and weights
###############################################
	basis_opt unwarp/${epi:t}_basis_opt.params unwarp/${epi:t}_basis_opt.weights -n3 -e0.01
####################################
# extract the optimized echo spacing
####################################
	set opt_dwell = `cat unwarp/${epi:t}_basis_opt.weights | awk '/Echo Spacing/{printf("%.7f", $NF/1000);}'`
########################################
# undistort epi with optimized field map
########################################
	$FSLDIR/bin/fugue --loadfmap=unwarp/${epi:t}_basis_opt_phase --dwell=$opt_dwell --in=unwarp/$epi:t -u unwarp/${epi:t}_uwrp --unwarpdir=$ped --saveshift=unwarp/${epi:t}_uwrp_shift_map
###########################################
# convert the shiftwarp to an absolute warp
###########################################
	$FSLDIR/bin/convertwarp -s unwarp/${epi:t}_uwrp_shift_map -r unwarp/$epi:t --shiftdir=$ped -o unwarp/${epi:t}_uwrp_shift_warp
#################################
# convert the result back to 4dfp
#################################
	gunzip -f unwarp/${epi:t}_uwrp.nii.gz
	nifti_4dfp -4 unwarp/${epi:t}_uwrp unwarp/${epi:t}_uwrp 
	exit $status
else
	echo $program "error: Method $method not recognized"; goto USAGE
endif

USAGE:
echo "${program}: distortion correction wrapper script for fMRI preprocessing"
echo "Usage:"
echo "measured  field fap: ${program} -map	<patid> <epi> <mag> <phase> <dwell> <te> <ped> <delta>"
echo "mean      field map: ${program} -mean	<epi> <FMmean> <epi_to_atl_t4> <dwell> <ped>"
echo "basis_opt field map: ${program} -basis	<epi> <t2w> <FMmean> <FMbases> <epi_to_t2w_t4> <epi_to_atl_t4> <dwell> <ped> <nbasis> [t2w brain mask]"
echo "N.B.:	{} indicates required NIfTI images"
echo "N.B.:	<te> is the <epi> echo time, expressed in msec (typically, ~30 msec at 3T)"
echo "N.B.:	<dwell> refers to the <epi> = 1/(BandwidthPerPixelPhaseEncode*#PhaseEncodes), expressed in msec"
echo "N.B.:	typical <dwell> values fall in the range .3 to .4 msec in 3T epi"
echo "N.B.:	with option -basis, basis_opt optimizes the <dwell> value (aka, echo spacing) by default"
echo "N.B.:	<delta> is the phase map sequence difference between echoes, expressed in msec (typically, 2.46 in Siemens data)"
echo "N.B.:	<ped> is the phase encode direction (typically, y-)"

exit 1
