#!/bin/csh
## Created by Avi Snyder and Tim Laumann 2013/2014
set basedir = /net/nil-bluearc/GMT/Laumann/MSC
set subject = $1
set subdir = ${basedir}/${subject}/Functionals
set outputdir = ${basedir}/${subject}/meanfield
set sesnums = `cat ${basedir}/${subject}/${subject}_func_sessions.txt`

goto MAKE_MOVIE



####################################
# Make movie to find bad fieldmaps
####################################
MAKE_MOVIE:
echo MAKE_MOVIE
mkdir ${outputdir}
set lst = ${outputdir}/sub.lst
if ( -e $lst ) /bin/rm $lst
touch $lst

foreach ses ( ${sesnums} )
	
	set sesfielddir = ${subdir}/${ses}/unwarp

	nifti_4dfp -4 ${sesfielddir}/${ses}_mag.nii ${sesfielddir}/${ses}_mag -N
	
	echo ${sesfielddir}/${ses}_mag >> $lst

end

paste_4dfp -ap1 $lst ${outputdir}/allsubs_mag_movie
#exit


####################################
# Register mag images to each other
####################################

REGISTER:
echo REGISTER
set modes = (0 0 0)
@ modes[1] = 3072 + 3 + 256
@ modes[2] = 2048 + 3 + 256
@ modes[3] = $modes[2] + 8192
@ n = $#sesnums
@ i = 1
while ($i <= $n)
	@ j = 1
	while ($j <= $n)
		if ($j != $i) then
			set t4file = $outputdir/$sesnums[$j]_mag_to_$sesnums[$i]_mag_t4
			@ k = 1
			while ($k <= $#modes)
			echo	imgreg_4dfp ${subdir}/$sesnums[$i]/unwarp/$sesnums[$i]_mag none ${subdir}/$sesnums[$j]/unwarp/$sesnums[$j]_mag none $t4file $modes[$k]
				imgreg_4dfp ${subdir}/$sesnums[$i]/unwarp/$sesnums[$i]_mag none ${subdir}/$sesnums[$j]/unwarp/$sesnums[$j]_mag none $t4file $modes[$k]
				if ($status) exit $status
				@ k++
			end
		endif
		@ j++
	end
	@ i++
end
#exit

#####################################
# Resolve t4 files from all sessions
#####################################

RESOLVE:
echo RESOLVE
@ n = $#sesnums
set images = ()
@ i = 1
while ($i <= $n)
	set images = ($images $sesnums[$i]_mag)
	@ i++
end
echo $images
cd $outputdir
t4_resolve $images -oresolved_mag -s >! ../mag_t4_resolve.log
echo "status="$status
#exit

##################################################
# Create list of mag images and their transforms
##################################################
RESOLVE_LIST:
echo RESOLVE_LIST
@ n = $#sesnums
cd $outputdir
set lst = mag_resolved.lst
if (-e $lst) /bin/rm $lst
touch $lst
@ i = 1
while ($i <= $n)
	echo $subdir/$sesnums[$i]/unwarp/$sesnums[$i]_mag"	t4="$sesnums[$i]_mag_to_resolved_mag_t4 >> $lst
	@ i++
end
cat $lst
#exit

#####################################################
# Create mag_resolved from mag_resolved.lst
#####################################################
RESOLVE_MAG_VOL:
echo RESOLVE_MAG_VOL
cd ${outputdir}
t4imgs_4dfp ${outputdir}/mag_resolved.lst mag_resolved -O${subdir}/${sesnums[1]}/unwarp/${sesnums[1]}_mag 
#exit

#####################################################
# Compute transform from mag_resolved to TRIO_Y_NDC
#####################################################
MAG_TO_TRIO:
echo MAG_TO_TRIO
set T1dir = ${basedir}/${subject}/T1
set T1_vol = ${subject}_mpr_debias_avgT
set T1_vol_mask = ${subject}_mpr1T_debias_bet
set T1_to_TRIO_t4 = ${subject}_mpr1T_debias_to_TRIO_Y_NDC_t4
cd ${outputdir}
imgreg_4dfp ${T1dir}/${T1_vol} ${T1dir}/${T1_vol_mask} mag_resolved none mag_resolved_to_T1_t4 2311
t4_mul mag_resolved_to_T1_t4 ${T1dir}/${T1_to_TRIO_t4} mag_resolved_to_TRIO_Y_NDC_t4
t4img_4dfp mag_resolved_to_TRIO_Y_NDC_t4 mag_resolved mag_resolved_on_TRIO_Y_NDC_111 -O111
#exit

#####################################
# Convert phase image to 4dfp
#####################################
PHASE:
echo PHASE
@ n = $#sesnums
@ i = 1
while ($i <= $n)
	pushd $subdir/$sesnums[$i]/unwarp
	/data/nil-bluearc/raichle/lin64-tools/niftigz_4dfp -4 $sesnums[$i]_phase_rad_unwrap_secs.nii.gz $sesnums[$i]_phase_rad_unwrap_secs -N
	if ($status) exit $status
	popd
	@ i++
end
#exit

##################################################
# Remove intensity scaling from resolved t4 files
##################################################
PHASE1:
echo PHASE1
@ n = $#sesnums
cd $outputdir
set lst = phase_rad_unwrap_secs_resolved.lst
if (-e $lst) /bin/rm $lst
touch $lst
@ i = 1
while ($i <= $n)
	sed '/scale:/d' $sesnums[$i]_mag_to_resolved_mag_t4 >! $sesnums[$i]_mag_to_resolved_mag_noscale_t4
	echo $subdir/$sesnums[$i]/unwarp/$sesnums[$i]_phase_rad_unwrap_secs"	t4="$sesnums[$i]_mag_to_resolved_mag_noscale_t4 >> $lst
	@ i++
end
cat $lst
#exit

###############################################################################################################################
# Creating transform from session fieldmap to TRIO_Y_NDC via mag_resolved, create list of phase images with their transforms
###############################################################################################################################
PHASE2:
echo PHASE2
@ n = $#sesnums
cd $outputdir
set lst = phase_rad_unwrap_secs_resolved_on_atlas.lst
if (-e $lst) /bin/rm $lst
touch $lst
@ i = 1
while ($i <= $n)
	t4_mul $sesnums[$i]_mag_to_resolved_mag_noscale_t4 mag_resolved_to_TRIO_Y_NDC_t4 $sesnums[$i]_mag_to_TRIO_Y_NDC_t4
	echo $subdir/$sesnums[$i]/unwarp/$sesnums[$i]_phase_rad_unwrap_secs"	t4="$sesnums[$i]_mag_to_TRIO_Y_NDC_t4 >> $lst
	@ i++
end
cat $lst
#exit

####################################################################################################
# Create mean field map in atlas space from phase_rad_unwrap_secs_resolved_on_atlas.lst
####################################################################################################
MEANFIELD_ATLAS:
echo MEANFIELD_ATLAS
cd ${outputdir}
t4imgs_4dfp phase_rad_unwrap_secs_resolved_on_atlas.lst phase_rad_unwrap_secs_resolved_on_TRIO_Y_NDC_111 -O111

exit




