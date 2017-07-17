#!/bin/csh
#$Header: /data/petsun4/data1/solaris/csh_scripts/RCS/cross_bold_pp_130702.csh,v 1.6 2014/03/12 22:29:20 avi Exp $
#$Log: cross_bold_pp_130702.csh,v $
# Revision 1.6  2014/03/12  22:29:20  avi
# code working for measured field maps
#
# Revision 1.5  2014/02/22  07:40:34  avi
# $min_frames instead of hard coded 240
#
# Revision 1.4  2014/02/21  07:00:51  avi
# correct minor bug in computation of ${patid}_func_vols
#
# Revision 1.3  2014/02/21  04:39:14  avi
# handle cross-day data (but t2w must be in current atlas directory)
#
# Revision 1.2  2013/11/08  06:23:59  avi
# optional pre-blur (parameter $anat_aveb) on func_vols run_dvar_4dfp
#
# Revision 1.1  2013/11/08  05:38:02  avi
# Initial revision
#

set idstr = '$Id: cross_bold_pp_130702.csh,v 1.6 2014/03/12 22:29:20 avi Exp $'
echo $idstr
set program = $0; set program = $program:t

if (${#argv} < 1) then
	echo "usage:	"$program" params_file [instructions_file]"
	exit 1
endif
set prmfile = $1
echo "prmfile="$prmfile

if (! -e $prmfile) then
	echo $program": "$prmfile not found
	exit -1
endif
source $prmfile
set instructions = ""
if (${#argv} > 1) then
	set instructions = $2
	if (! -e $instructions) then
		echo $program": "$instructions not found
		exit -1
	endif
	cat $instructions
	source $instructions
endif

##########
# check OS
##########
set OS = `uname -s`
if ($OS != "Linux") then
	echo $program must be run on a linux machine
	exit -1
endif

if ($target:h != $target) then
	set tarstr = -T$target
else
	set tarstr = $target
endif

@ runs = ${#irun}
if ($runs != ${#fstd}) then
	echo "irun fstd mismatch - edit "$prmfile
	exit -1
endif

if (! ${?scrdir}) set scrdir = ""
@ usescr = `echo $scrdir | awk '{print length ($1)}'`
if ($usescr) then 
	if (! -e $scrdir) mkdir $scrdir
	if ($status) exit $status
endif
set sourcedir = $cwd
if (! ${?sorted}) @ sorted = 0

if (! ${?MB}) @ MB = 0			# skip slice timing correction and debanding
set MBstr = _faln_dbnd
if ($MB) set MBstr = ""

set squeezestr = ""
if (${?sx}) then
	set squeezestr = $squeezestr" -sx"$sx
endif
if (${?sy}) then
	set squeezestr = $squeezestr" -sy"$sy
endif

if (! ${?E4dfp}) @ E4dfp = 0
if (${E4dfp}) then
	echo "4dfp files have been pre-generated. Option E4dfp set with value $E4dfp. Skipping dcm_to_4dfp"
endif

if (! ${?use_anat_ave}) @ use_anat_ave = 0
if ($use_anat_ave) then
	set epi_anat = $patid"_anat_ave"
else
	set epi_anat = $patid"_func_vols_ave"
endif
if (! ${?min_frames}) @ min_frames = 240
if (! ${?day1_patid}) set day1_patid = "";

if (${?goto_UNWARP}) goto UNWARP

date

#goto REGISTER

###################
# process BOLD data
###################
if (${?epi_zflip}) then
	if ($epi_zflip) set zflip = "-z"
else
	set zflip = ""
endif
set interleave = ""
if (${?Siemens_interleave}) then
	if ($Siemens_interleave) set interleave = "-N"
endif
@ err = 0
@ k = 1
while ($k <= $runs)
	if (! $E4dfp) then
		if ($usescr) then		# test to see if user requested use of scratch disk
			if (-e bold$irun[$k]) /bin/rm bold$irun[$k]	# remove existing link
			if (! -d $scrdir/bold$irun[$k]) mkdir $scrdir/bold$irun[$k]
			ln -s $scrdir/bold$irun[$k] bold$irun[$k]
		else
			if (! -d bold$irun[$k]) mkdir bold$irun[$k]
		endif
	endif
	pushd bold$irun[$k]
	set y = $patid"_b"$irun[$k]${MBstr}
	if (-e $y.4dfp.img && -e $y.4dfp.ifh) goto POP
	if (! $E4dfp) then
		if ($sorted) then
			echo		dcm_to_4dfp -q -b study$fstd[$k] $inpath/study$fstd[$k]
			if ($go)	dcm_to_4dfp -q -b study$fstd[$k] $inpath/study$fstd[$k]
		else
			echo		dcm_to_4dfp -q -b study$fstd[$k] $inpath/$dcmroot.$fstd[$k]."*"
			if ($go)	dcm_to_4dfp -q -b study$fstd[$k] $inpath/$dcmroot.$fstd[$k].*
		endif
	endif
	endif
	echo		unpack_4dfp -V study$fstd[$k] $patid"_b"$irun[$k] -nx$nx -ny$ny $squeezestr $zflip
	if ($go)	unpack_4dfp -V study$fstd[$k] $patid"_b"$irun[$k] -nx$nx -ny$ny $squeezestr $zflip
	if ($status) then
		@ err++
		/bin/rm $patid"_b"$irun[$k]*
		goto POP
	endif
	echo		/bin/rm  study$fstd[$k]."*"
	if ($go)	/bin/rm  study$fstd[$k].*

	if ($MB) goto POP
	echo		frame_align_4dfp $patid"_b"$irun[$k] $skip -TR_vol $TR_vol -TR_slc $TR_slc -d $epidir $interleave
	if ($go)	frame_align_4dfp $patid"_b"$irun[$k] $skip -TR_vol $TR_vol -TR_slc $TR_slc -d $epidir $interleave

	echo		deband_4dfp -n$skip $patid"_b"$irun[$k]"_faln"
	if ($go)	deband_4dfp -n$skip $patid"_b"$irun[$k]"_faln"
	if ($status)	exit $status

	if ($economy > 2) then
		echo		/bin/rm $patid"_b"$irun[$k].4dfp."*"
		if ($go)	/bin/rm $patid"_b"$irun[$k].4dfp.*
	endif
	if ($economy > 3) then
		echo		/bin/rm $patid"_b"$irun[$k]"_faln".4dfp."*"
		if ($go)	/bin/rm $patid"_b"$irun[$k]"_faln".4dfp.*
	endif
POP:
	popd	# out of bold$irun[$k]
	@ k++
end
if ($err) then
	echo $program": one or more BOLD runs failed preliminary processing"
	exit -1
endif
if ($epi2atl == 2) goto ATL

if (-e  $patid"_xr3d".lst)	/bin/rm $patid"_xr3d".lst;	touch $patid"_xr3d".lst
if (-e  $patid"_anat".lst)	/bin/rm $patid"_anat".lst;	touch $patid"_anat".lst
@ k = 1
while ($k <= $runs)
	echo bold$irun[$k]/$patid"_b"$irun[$k]${MBstr} >>		$patid"_xr3d".lst
	echo bold$irun[$k]/$patid"_b"$irun[$k]${MBstr}_xr3d_norm >>	$patid"_anat".lst
	@ k++
end
/bin/cp $patid"_anat".lst $patid"_func_vols".lst

echo cat	$patid"_xr3d".lst
cat		$patid"_xr3d".lst
echo		cross_realign3d_4dfp -n$skip -qv$normode -l$patid"_xr3d".lst
if ($go)	cross_realign3d_4dfp -n$skip -qv$normode -l$patid"_xr3d".lst
if ($status)	exit $status

date
#################################
# compute mode 1000 normalization
#################################
@ k = 1
while ($k <= $runs)
	pushd bold$irun[$k]
	echo 		normalize_4dfp $patid"_b"$irun[$k]${MBstr}"_r3d_avg" -h
	if ($go)	normalize_4dfp $patid"_b"$irun[$k]${MBstr}"_r3d_avg" -h
	if ($economy > 4 && $epi2atl == 0) then
		echo		/bin/rm $patid"_b"$irun[$k]$MBstr.4dfp."*"
		if ($go)	/bin/rm $patid"_b"$irun[$k]$MBstr.4dfp.*
	endif
	popd	# out of bold$irun[$k]
	@ k++
end

date
###############################
# apply mode 1000 normalization
###############################
@ k = 1
while ($k <= $runs)
	pushd bold$irun[$k]
	set file = $patid"_b"$irun[$k]${MBstr}"_r3d_avg_norm".4dfp.img.rec
	set f = 1.0; if (-e $file) set f = `head $file | awk '/original/{print 1000/$NF}'`
	echo		scale_4dfp $patid"_b"$irun[$k]${MBstr}"_xr3d" $f -anorm
	if ($go)	scale_4dfp $patid"_b"$irun[$k]${MBstr}"_xr3d" $f -anorm
	echo		/bin/rm $patid"_b"$irun[$k]${MBstr}"_xr3d".4dfp."*"
	if ($go)	/bin/rm $patid"_b"$irun[$k]${MBstr}"_xr3d".4dfp.*
	popd	# out of bold$irun[$k]
	@ k++
end

date
###################
# movement analysis
###################
if (! -d movement) mkdir movement
@ k = 1
while ($k <= $runs)
	echo		mat2dat bold$irun[$k]/"*_xr3d".mat -RD -n$skip
	if ($go)	mat2dat bold$irun[$k]/*"_xr3d".mat -RD -n$skip
	echo		/bin/mv bold$irun[$k]/"*_xr3d.*dat"	movement
	if ($go)	/bin/mv bold$irun[$k]/*"_xr3d".*dat	movement
	@ k++
end

date
if (! -d atlas) mkdir atlas

######################################
# make EPI first frame (anatomy) image
######################################
echo cat	$patid"_anat".lst
cat		$patid"_anat".lst
echo		paste_4dfp -p1 $patid"_anat".lst	$patid"_anat_ave"
if ($go)	paste_4dfp -p1 $patid"_anat".lst	$patid"_anat_ave"
echo		ifh2hdr	-r2000				$patid"_anat_ave"
if ($go)	ifh2hdr	-r2000				$patid"_anat_ave"
echo		/bin/mv $patid"_anat*" atlas
if ($go)	/bin/mv $patid"_anat"* atlas


#######################################
# make func_vols_ave using actmapf_4dfp
#######################################
echo	conc_4dfp ${patid}_func_vols -l${patid}_func_vols.lst
	conc_4dfp ${patid}_func_vols -l${patid}_func_vols.lst
if ($status) exit $status
cat			${patid}_func_vols.conc
echo		/bin/mv	${patid}_func_vols."*" atlas
if ($go)	/bin/mv	${patid}_func_vols.*   atlas

pushd atlas		# into atlas
if (! ${?anat_aveb}) set anat_aveb = 0.
foreach k (1 2)		# iterate
	if (! -e ${patid}_func_vols.format) then
		set  format = `conc2format ${patid}_func_vols.conc $skip`
		echo $format
		echo $format >! ${patid}_func_vols.format
		echo hello
	else
		set format = `cat ${patid}_func_vols.format`
		set str = `format2lst -e $format | gawk '{k=0;l=length($1);for(i=1;i<=l;i++)if(substr($1,i,1)=="x")k++;}END{print k, l;}'`
		echo $str[1] out of $str[2] frames fail dvar criterion $anat_avet 
	#	@ j = $str[2] - $str[1]; if ($j < $min_frames) exit 1	# require at least $min_frames with dvar < $anat_avet to proceed
	endif
	echo	actmapf_4dfp "$format" ${patid}_func_vols.conc -aave
		actmapf_4dfp "$format" ${patid}_func_vols.conc -aave
	echo		ifh2hdr	-r2000	${patid}_func_vols_ave
	if ($go)	ifh2hdr	-r2000	${patid}_func_vols_ave
	if ($status) exit $status
	echo	zero_lt_4dfp 500 ${patid}_func_vols_ave
		zero_lt_4dfp 500 ${patid}_func_vols_ave
	if ($status) exit $status
	echo	run_dvar_4dfp ${patid}_func_vols.conc -m${patid}_func_vols_avez -n$skip -x$anat_avet -b$anat_aveb
		run_dvar_4dfp ${patid}_func_vols.conc -m${patid}_func_vols_avez -n$skip -x$anat_avet -b$anat_aveb
	if ($status) exit $status
end


REGISTER:
#pushd atlas
#################################################
# Run epi2t2w2mpr2atl_4dfp, register BOLD to atlas space
#################################################
set t2w = $patid"_t2w"
echo		epi2t2w2mpr2atl2_4dfp ${epi_anat} $t2w $patid"_mpr1" useold $tarstr
if ($go)	epi2t2w2mpr2atl2_4dfp ${epi_anat} $t2w $patid"_mpr1" useold $tarstr
if ($status) exit $status

EPI_to_ATL:
if (! $use_anat_ave && $day1_patid == "") then
	/bin/rm ${patid}_anat_ave_to_${target:t}_t4
	ln -s ${patid}_func_vols_ave_to_${target:t}_t4 ${patid}_anat_ave_to_${target:t}_t4
endif
if ($day1_patid != "") then
set echo
	ln -s ${patid}_anat_ave_to_${target:t}_t4 ${patid}_func_vols_ave_to_${target:t}_t4
	t4_mul ${patid}_anat_ave_to_${day1_patid}_anat_ave_t4 ${day1_patid}_anat_ave_to_${day1_patid}_t2w_t4
	/bin/rm						${patid}_func_vols_ave_to_${day1_patid}_t2w_t4
	ln -s ${patid}_anat_ave_to_${day1_patid}_t2w_t4 ${patid}_func_vols_ave_to_${day1_patid}_t2w_t4
endif

########################################################################
# make atlas transformed epi_anat and t2w in 111 222 and 333 atlas space
########################################################################
set t4file = ${patid}_anat_ave_to_${target:t}_t4
foreach O (111 222 333)
	echo		t4img_4dfp $t4file  ${epi_anat}	${epi_anat}_on_${target:t}_$O -O$O
	if ($go)	t4img_4dfp $t4file  ${epi_anat}	${epi_anat}_on_${target:t}_$O -O$O
	echo		ifh2hdr	 -r2000			${epi_anat}_on_${target:t}_$O
	if ($go)	ifh2hdr	 -r2000			${epi_anat}_on_${target:t}_$O
end
if ($status) exit $status

if ($day1_patid != "") goto SKIPT2W
set t4file = ${t2w}_to_${target:t}_t4
foreach O (111 222 333)
	echo		t4img_4dfp $t4file  ${t2w}	${t2w}_on_${target:t}_$O -O$O
	if ($go)	t4img_4dfp $t4file  ${t2w}	${t2w}_on_${target:t}_$O -O$O
	echo		ifh2hdr	 -r1000			${t2w}_on_${target:t}_$O
	if ($go)	ifh2hdr	 -r1000			${t2w}_on_${target:t}_$O
end
if ($status) exit $status
SKIPT2W:
/bin/rm *t4% >& /dev/null
popd		# out of atlas

#exit

UNWARP:
set t2w = $patid"_t2w"
################################
# convert dcm field maps to nii
################################
echo $gfm
if ($#gfm != 2) then
	echo $patid does not have gre_field_mapping
	exit 1
endif

dcm_to_4dfp -d 0018 0081 -b fieldmap_mag study$gfm[1]
echo fieldmap_mag.4dfp.img > temp.lst
echo fieldmap_mag_1.4dfp.img >> temp.lst
paste_4dfp -a temp.lst fieldmap_mags
nifti_4dfp -n fieldmap_mags ${patid}_amp
rm fieldmap_mag.4dfp.* fieldmap_mag_1.4dfp.* fieldmap_mags.4dfp.*

dcm_to_4dfp -b fieldmap_pha study$gfm[2] 
nifti_4dfp -n fieldmap_pha ${patid}_pha
rm fieldmap_pha.4dfp.*

#dcm2nii -d N -e N -f N -g N -i N	study$gfm[1]; if ($status) exit $status
#mv					study$gfm[1]/grefieldmapping*Slices.nii ${patid}_amp.nii
#dcm2nii -d N -e N -f N -g N -i N	study$gfm[2]; if ($status) exit $status
#mv					study$gfm[2]/grefieldmapping*Slices.nii ${patid}_pha.nii
#exit

################################
# compute fMRI distortion unwarp
################################
if ($status) exit $status
if ($day1_patid != "") then
	set patid1	= $day1_patid
else
	set patid1	= $patid
endif


##############################################################
# setting arguments for unwarp with measured field maps
##############################################################
set uwrp_args  = (-map $patid atlas/${epi_anat} ${patid}_amp.nii ${patid}_pha.nii $dwell $TE_vol $ped $delta) # This will change

##################################
# compute field mapping correction
##################################
set x = ${uwrp_cmnd:t}; set x = $x:r
set log		= ${patid}_$x.log
date				>! $log
echo  $uwrp_cmnd $uwrp_args
echo	$uwrp_cmnd $uwrp_args	>> $log
echo hello
	$uwrp_cmnd $uwrp_args	>> $log
if ($status) exit $status

if (! -e unwarp/${epi_anat}_uwrp_to_${t2w}_t4) then
	niftigz_4dfp -n atlas/$t2w atlas/$t2w
	bet atlas/$t2w atlas/${t2w}_brain -m -f 0.4 -R
	niftigz_4dfp -4 atlas/${t2w}_brain_mask atlas/${t2w}_brain_mask -N
	@ mode = 8192 + 2048 + 3
	/bin/cp atlas/${epi_anat}_to_${t2w}_t4 unwarp/${epi_anat}_uwrp_to_${t2w}_t4
	imgreg_4dfp atlas/$t2w atlas/${t2w}_brain_mask unwarp/${epi_anat}_uwrp none unwarp/${epi_anat}_uwrp_to_${t2w}_t4 $mode >! unwarp/${epi_anat}_uwrp_to_${t2w}.log
	if ($status) exit $status
endif

echo	t4_mul	unwarp/${epi_anat}_uwrp_to_${t2w}_t4 atlas/${t2w}_to_${target:t}_t4 unwarp/${epi_anat}_uwrp_to_${target:t}_t4
	t4_mul	unwarp/${epi_anat}_uwrp_to_${t2w}_t4 atlas/${t2w}_to_${target:t}_t4 unwarp/${epi_anat}_uwrp_to_${target:t}_t4
set t4file =    unwarp/${epi_anat}_uwrp_to_${target:t}_t4
foreach O (111 222 333)
	echo		t4img_4dfp $t4file  unwarp/${epi_anat}_uwrp	unwarp/${epi_anat}_uwrp_on_${target:t}_$O -O$O
	if ($go)	t4img_4dfp $t4file  unwarp/${epi_anat}_uwrp	unwarp/${epi_anat}_uwrp_on_${target:t}_$O -O$O
	echo		ifh2hdr	 -r2000					unwarp/${epi_anat}_uwrp_on_${target:t}_$O
	if ($go)	ifh2hdr	 -r2000					unwarp/${epi_anat}_uwrp_on_${target:t}_$O
end

if ($status) exit $status
echo	t4img_4dfp unwarp/${epi_anat}_uwrp_to_${t2w}_t4 unwarp/${epi_anat}_uwrp 	unwarp/${epi_anat}_uwrp_on_${t2w} -Oatlas/${t2w}
	t4img_4dfp unwarp/${epi_anat}_uwrp_to_${t2w}_t4 unwarp/${epi_anat}_uwrp 	unwarp/${epi_anat}_uwrp_on_${t2w} -Oatlas/${t2w}
if ($status) exit $status
ifh2hdr -r2000										unwarp/${epi_anat}_uwrp_on_${t2w}

ATL:
#################################
# one step resample unwarped fMRI
#################################
#if (! $epi2atl) exit 0
#set x = ${rsam_cmnd:t}; set x = $x:r
#set log		= ${patid}_$x.log
#date						>! $log
#echo	$rsam_cmnd $prmfile $instructions	>> $log
#	$rsam_cmnd $prmfile $instructions	>> $log
#if ($status) exit $status

####################################################################
# remake single resampled 333 atlas space fMRI volumetric timeseries
####################################################################
#set lst = ${patid}${MBstr}_xr3d_uwrp_atl.lst
#if (-e $lst) /bin/rm $lst
#touch $lst
#@ k = 1
#while ($k <= $#irun)
#	echo bold$irun[$k]/${patid}_b$irun[$k]${MBstr}_xr3d_uwrp_atl.4dfp.img >> $lst
#	@ k++
#end
#conc_4dfp ${lst:r}.conc -l$lst
#if ($status) exit $status
#set format = `cat atlas/${patid}_func_vols.format`
#if ($status) exit $status
#actmapf_4dfp $format	${patid}${MBstr}_xr3d_uwrp_atl.conc -aave
#if ($status) exit $status
#ifh2hdr -r2000 		${patid}${MBstr}_xr3d_uwrp_atl_ave
#mv			${patid}${MBstr}_xr3d_uwrp_atl_ave.4dfp.*	atlas
#var_4dfp -sf$format	${patid}${MBstr}_xr3d_uwrp_atl.conc
#ifh2hdr -r20		${patid}${MBstr}_xr3d_uwrp_atl_sd1
#mv			${patid}${MBstr}_xr3d_uwrp_atl_sd1*		atlas
#mv			${patid}${MBstr}_xr3d_uwrp_atl.conc*		atlas

#echo $program complete
exit 0
