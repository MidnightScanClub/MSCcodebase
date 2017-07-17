#!/bin/csh
## TOL, Version 1, 09/2014
set basedir = /data/nil-bluearc/GMT/Ortega
#set basedir = /data/nil-bluearc/GMT/Laumann/MSC/
set REFDIR = /data/petsun43/data1/atlas
set scriptdir = /data/nil-bluearc/GMT/Laumann/MSC/Scripts
set subjects = CIMT001T #(MSC04 MSC05 MSC06 MSC07 MSC08 MSC09 MSC10) #`cat $1`

foreach subject ( $subjects )

set subdir = $basedir/${subject}/
set sesnums = `cat ${basedir}/${subject}/${subject}_func_sessions.txt`

goto T1_PROC_nodebias


DCM_SORT_STRUCT:
################
# Sort dcm files
################
set sesnums_orig = `cat ${basedir}/${subject}/${subject}_struct_sessions_orig.txt`
set sesnums = `cat ${basedir}/${subject}/${subject}_struct_sessions.txt`
set origdir = ${basedir}/rawdata/${subject}/STRUCTURALS
#set origdir = /data/heisenberg/data1/MSC/rawdata/${subject}/STRUCTURALS
set structdir = ${subdir}/Structurals
mkdir ${structdir}
set k = 1

while ( $k <= $#sesnums) 
	
	set ses_orig = $sesnums_orig[$k]
	set ses = $sesnums[$k]
	mkdir ${structdir}/$ses
	pushd ${structdir}/$ses
	echo /data/nil-bluearc/raichle/lin64-tools/pseudo_dcm_sort.csh -i -s ${origdir}/${ses_orig}/scans
	/data/nil-bluearc/raichle/lin64-tools/pseudo_dcm_sort.csh -i -s ${origdir}/${ses_orig}/scans
	mv scans.studies.txt ${ses}.studies.txt
	popd
	@ k++
end

exit

DCM_SORT_FUNC:
################
# Sort dcm files
################
set sesnums_orig = `cat ${basedir}/${subject}/${subject}_func_sessions_orig.txt`
set sesnums = `cat ${basedir}/${subject}/${subject}_func_sessions.txt`
set origdir = /data/nil-bluearc/GMT/MSC/rawdata/${subject}/FUNCTIONALS
#set origdir = /data/heisenberg/data1/MSC/rawdata/${subject}/FUNCTIONALS
set funcdir = ${subdir}/Functionals
mkdir ${funcdir}
set k = 1

while ( $k <= $#sesnums) 
	
	set ses_orig = $sesnums_orig[$k]
	set ses = $sesnums[$k]
	mkdir ${funcdir}/$ses
	pushd ${funcdir}/$ses
	dcm_sort  ${origdir}/${ses_orig}
	#pseudo_dcm_sort.csh -i -s ${origdir}/${ses_orig}/SCANS
	#echo /data/nil-bluearc/raichle/lin64-tools/pseudo_dcm_sort.csh -i -s ${origdir}/${ses_orig}/scans
	#/data/nil-bluearc/raichle/lin64-tools/pseudo_dcm_sort.csh -i -s ${origdir}/${ses_orig}/scans
	mv scans.studies.txt ${ses}.studies.txt
	popd
	@ k++
end

exit

STRUCT_PARAMS:
################
# create struct params
################
set sesnums = `cat ${basedir}/${subject}/${subject}_struct_sessions.txt`
pushd ${basedir}/${subject}/Structurals

set T1_label = 
set T2_label = 
set MRA_label = 
set MRV_cor_label = 
set MRV_sag_label = 
foreach k ( $sesnums )
	pushd $k	
	set T1dcm = `cat $k.studies.txt | grep 'T1_.8x' | grep ' 224' | awk -F " " '{print $1}'`
	set T1dcm = `echo $T1dcm | tr -d '\n'`
	set T1num = $#T1dcm
	set t = 1
	while ( $t <= $T1num )
		if ( $T1dcm[$t] > 0 ) then 
			set T1_label = `echo ${T1_label} $k/study$T1dcm[$t]`			
		endif
		@ t++
	end

	set T2dcm = `cat $k.studies.txt | grep 'T2_.8x' | grep ' 224' | awk -F " " '{print $1}'`
	set T2dcm = `echo $T2dcm | tr -d '\n'`
	set T2num = $#T2dcm
	set t = 1
	while ( $t <= $T2num )
		if ( $T2dcm[$t] > 0 ) then 
			set T2_label = `echo ${T2_label} $k/study$T2dcm[$t]`			
		endif
		@ t++
	end

	set MRAdcm = `cat $k.studies.txt | grep 'MRA_whole' | grep ' 268' | awk -F " " '{print $1}'`
	set MRAdcm = `echo $MRAdcm | tr -d '\n'`
	set MRAnum = $#MRAdcm
	set t = 1
	while ( $t <= $MRAnum )
		if ( $MRAdcm[$t] > 0 ) then 
			set MRA_label = `echo ${MRA_label} $k/study$MRAdcm[$t]`			
		endif
		@ t++
	end

	set MRV_cor_dcm = `cat $k.studies.txt | grep 'MRV_cor' | grep ' 128' | awk -F " " '{print $1}'`
	set MRV_cor_dcm = `echo $MRV_cor_dcm | tr -d '\n'`
	set MRV_cor_num = $#MRV_cor_dcm
	set t = 1
	while ( $t <= $MRV_cor_num )
		if ( $MRV_cor_dcm[$t] > 0 ) then 
			set MRV_cor_label = `echo ${MRV_cor_label} $k/study$MRV_cor_dcm[$t]`			
		endif
		@ t++
	end

	set MRV_sag_dcm = `cat $k.studies.txt | grep 'MRV_sag' | grep ' 120' | awk -F " " '{print $1}'`
	set MRV_sag_dcm = `echo $MRV_sag_dcm | tr -d '\n'`
	set MRV_sag_num = $#MRV_sag_dcm
	set t = 1
	while ( $t <= $MRV_sag_num )
		if ( $MRV_sag_dcm[$t] > 0 ) then 
			set MRV_sag_label = `echo ${MRV_sag_label} $k/study$MRV_sag_dcm[$t]`			
		endif
		@ t++
	end
	popd
end

echo "set T1    = ( ${T1_label} )" > ${subject}.structparams
echo "set T2    = ( ${T2_label} )" >> ${subject}.structparams
echo "set MRA    = ( ${MRA_label} )" >> ${subject}.structparams
echo "set MRV_cor    = ( ${MRV_cor_label} )" >> ${subject}.structparams
echo "set MRV_sag    = ( ${MRV_sag_label} )" >> ${subject}.structparams
#exit


STRUCT_DCM:
################
# Convert dcm to 4dfp
################
set structtype = ( T1 T2 MRA MRV_cor MRV_sag )

pushd ${subdir}
source ${subdir}/Structurals/${subject}.structparams
foreach struct ( $structtype )
	mkdir $struct
end

set k = 1
pushd $structtype[1]
while ( $k <= $#T1 )
	set structscan = $T1[$k]	
	dcm_to_4dfp -b ${subject}_mpr$k ../Structurals/$structscan
	@ k++
end
popd
set k = 1
pushd $structtype[2]
while ( $k <= $#T2 )
	set structscan = $T2[$k]
	dcm_to_4dfp -b ${subject}_t2w$k ../Structurals/$structscan
	@ k++
end
popd

set k = 1
pushd $structtype[3]
while ( $k <= $#MRA )
	set structscan = $MRA[$k]
	dcm_to_4dfp -t T ../Structurals/$structscan
	rename analyze ${subject}_${structtype[3]}$k *
#	dcm2nii -d N -e N -f N -g N -i N ../Structurals/$structscan
#	mv ../Structurals/$structscan/MRAwholebrain.nii ./${subject}_${structtype[3]}$k.nii
#	rm ../Structurals/$structscan/*MRAwhole*
#	nifti_4dfp -4 ./${subject}_${structtype[3]}$k ./${subject}_${structtype[3]}$k -N
#	rm ./${subject}_${structtype[3]}$k.nii
	@ k++
end
popd

set k = 1
pushd $structtype[4]
while ( $k <= $#MRV_cor )
	set structscan = $MRV_cor[$k]
	dcm_to_4dfp -t T ../Structurals/$structscan
	rename analyze ${subject}_${structtype[4]}$k *
#	dcm2nii -d N -e N -f N -g N -i N ../Structurals/$structscan
#	mv ../Structurals/$structscan/MRVcoronalIpat.nii ./${subject}_${structtype[4]}$k.nii
#	rm ../Structurals/$structscan/*MRVcoronal*
#	nifti_4dfp -4 ./${subject}_${structtype[4]}$k ./${subject}_${structtype[4]}$k -N
#	rm ./${subject}_${structtype[4]}$k.nii
	@ k++
end
popd

set k = 1
pushd $structtype[5]
while ( $k <= $#MRV_sag)
	set structscan = $MRV_sag[$k]
	dcm_to_4dfp -t T ../Structurals/$structscan
	rename analyze ${subject}_${structtype[5]}$k *
#	dcm2nii -d N -e N -f N -g N -i N ../Structurals/$structscan
#	mv ../Structurals/$structscan/MRVsagsinus.nii ./${subject}_${structtype[5]}$k.nii
#	rm ../Structurals/$structscan/*MRVsag*
#	nifti_4dfp -4 ./${subject}_${structtype[5]}$k ./${subject}_${structtype[5]}$k -N
#	rm ./${subject}_${structtype[5]}$k.nii
	@ k++
end
popd
#exit

T1_PROC:
###############
# Register T1 to atlas, debias, and average
###############

set structdir = ${basedir}/${subject}/T1
source ${subdir}/Structurals/${subject}.structparams
pushd ${structdir}
set T1num = $#T1
#goto HERE1
# Transform from Sagittal to Transverse
set k = 1
while ( $k <= $T1num )
	S2T_4dfp ${subject}_mpr${k} ${subject}_mpr${k}T
	niftigz_4dfp -n ${subject}_mpr${k}T ${subject}_mpr${k}T
	@ k++
end

# Debias and convert back to 4dfp
set k = 1
while ( $k <= $T1num )
	echo ${scriptdir}/apply_debias.csh ${subject}_mpr${k}T
	${scriptdir}/apply_debias.csh ${subject}_mpr${k}T
	niftigz_4dfp -4 ${subject}_mpr${k}T_debias ${subject}_mpr${k}T_debias -N
	@ k++
end

# Mask first T1 for registration	
echo bet2 ${subject}_mpr1T_debias ${subject}_mpr1T_debias_bet
bet2 ${subject}_mpr1T_debias.nii.gz ${subject}_mpr1T_debias_bet
niftigz_4dfp -4 ${subject}_mpr1T_debias_bet ${subject}_mpr1T_debias_bet -N

#HERE1:
# Register first T1 to atlas
set modes	= (0 0 0 0 0)
@ modes[1]	= 1024 + 256 + 3
@ modes[2]	= 1024 + 256 + 3
@ modes[3]	= 3072 + 256 + 7
@ modes[4]	= 2048 + 256 + 7
@ modes[5]	= 2048 + 256 + 7
set t4file = ${subject}_mpr1T_debias_to_TRIO_Y_NDC_t4
set ref = /data/petsun43/data1/atlas/TRIO_Y_NDC
set mpr = ${subject}_mpr1T_debias
set mpr_mask = ${subject}_mpr1T_debias_bet
set log = ${subject}_mpr1T_to_TRIO_Y_NDC.log
@ k = 1
while ( $k <= $#modes )
	imgreg_4dfp $ref none $mpr $mpr_mask $t4file $modes[$k] >> $log
	@ k++
end

# Average T1s
set k = 1
set T1scans = ()
while ( $k <= $T1num )
	set T1scans = ( ${T1scans} ${subject}_mpr${k}T_debias )
	@ k++
end
avgmpr_4dfp ${T1scans} ${subject}_mpr_debias_avgT useold -T/data/petsun43/data1/atlas/TRIO_Y_NDC
t4imgs_4dfp -s ${subject}_mpr_debias_avgT.lst ${subject}_mpr_debias_avgT -O${subject}_mpr1T_debias

#exit
#goto T2_PROC

T1_PROC_nodebias:
###############
# Register T1 to atlas, debias, and average
###############
set structdir = ${basedir}/${subject}/T1
source ${subdir}/Structurals/${subject}.structparams
pushd ${structdir}
set T1num = $#T1

# Mask first T1 for registration	
echo bet2 ${subject}_mpr1T ${subject}_mpr1T_bet
bet2 ${subject}_mpr1T.nii.gz ${subject}_mpr1T_bet
niftigz_4dfp -4 ${subject}_mpr1T_bet ${subject}_mpr1T_bet -N

# Register first T1 to atlas
set modes	= (0 0 0 0 0)
@ modes[1]	= 1024 + 256 + 3
@ modes[2]	= 1024 + 256 + 3
@ modes[3]	= 3072 + 256 + 7
@ modes[4]	= 2048 + 256 + 7
@ modes[5]	= 2048 + 256 + 7
set t4file = ${subject}_mpr1T_to_TRIO_Y_NDC_t4
set ref = /data/petsun43/data1/atlas/TRIO_Y_NDC
set mpr = ${subject}_mpr1T
set mpr_mask = ${subject}_mpr1T_bet
set log = ${subject}_mpr1T_to_TRIO_Y_NDC.log
@ k = 1
while ( $k <= $#modes )
	imgreg_4dfp $ref none $mpr $mpr_mask $t4file $modes[$k] >> $log
	@ k++
end

# Average T1s
set k = 1
set T1scans = ()
while ( $k <= $T1num )
	set T1scans = ( ${T1scans} ${subject}_mpr${k}T )
	@ k++
end
avgmpr_4dfp ${T1scans} ${subject}_mpr_avgT useold -T/data/petsun43/data1/atlas/TRIO_Y_NDC
t4imgs_4dfp -s ${subject}_mpr_avgT.lst ${subject}_mpr_avgT -O${subject}_mpr1T

goto T2_PROC_nodebias
exit

T2_PROC:
###############
# Register T2 to T1 and average
###############

set structdir = ${basedir}/${subject}/T2
source ${subdir}/Structurals/${subject}.structparams
pushd ${structdir}
set T2num = $#T2
#goto HERE2
ln -s ../T1/${subject}_mpr1T_debias_to_TRIO_Y_NDC_t4 .
foreach e ( img img.rec ifh hdr )
	ln -s ../T1/${subject}_mpr1T_debias.4dfp.$e .
end
# Transform from sagittal to transverse
@ k = 1
while ( $k <= $T2num )
	S2T_4dfp ${subject}_t2w${k} ${subject}_t2w${k}T
	niftigz_4dfp -n ${subject}_t2w${k}T ${subject}_t2w${k}T
	@ k++
end

# Debias and convert back to 4dfp
set k = 1
while ( $k <= $T2num )
	echo ${scriptdir}/apply_debias.csh ${subject}_t2w${k}T
	${scriptdir}/apply_debias.csh ${subject}_t2w${k}T
	niftigz_4dfp -4 ${subject}_t2w${k}T_debias ${subject}_t2w${k}T_debias -N
	@ k++
end

#HERE2:
# Register T2 to T1
t2w2mpr_4dfp ${subject}_mpr1T_debias ${subject}_t2w1T_debias -T/data/petsun43/data1/atlas/TRIO_Y_NDC
set k = 1
set T2scans = ()
while ( $k <= $T2num )
	set T2scans = ( ${T2scans} ${subject}_t2w${k}T_debias )
	@ k++
end
avgmpr_4dfp ${T2scans} ${subject}_t2w_debias_avgT useold -T/data/petsun43/data1/atlas/TRIO_Y_NDC
t4imgs_4dfp -s ${subject}_t2w_debias_avgT.lst ${subject}_t2w_debias_avgT -O${subject}_t2w1T_debias
cp ${subject}_t2w1T_debias_to_TRIO_Y_NDC_t4 ${subject}_t2w_debias_avgT_to_TRIO_Y_NDC_t4
cp ${subject}_t2w1T_debias_to_${subject}_mpr1T_debias_t4 ${subject}_t2w_debias_avgT_to_${subject}_mpr1T_debias_t4
t4img_4dfp ${subject}_t2w_debias_avgT_to_${subject}_mpr1T_debias_t4 ${subject}_t2w_debias_avgT ${subject}_t2w_debias_avgT_on_${subject}_mpr1T_debias -O${subject}_mpr1T_debias.4dfp.ifh

# Create masked T2
niftigz_4dfp -n ${subject}_t2w_debias_avgT ${subject}_t2w_debias_avgT
bet2 ${subject}_t2w_debias_avgT ${subject}_t2w_debias_avgT_bet
niftigz_4dfp -4 ${subject}_t2w_debias_avgT_bet ${subject}_t2w_debias_avgT_bet
exit

T2_PROC_nodebias:
###############
# Register T2 to T1 and average
###############

set structdir = ${basedir}/${subject}/T2
source ${subdir}/Structurals/${subject}.structparams
pushd ${structdir}
echo ${structdir}
set T2num = $#T2

#goto HERE2
ln -s ../T1/${subject}_mpr1T_to_TRIO_Y_NDC_t4 .
foreach e ( img img.rec ifh hdr )
	ln -s ../T1/${subject}_mpr1T.4dfp.$e .
end

# Register T2 to T1
t2w2mpr_4dfp ${subject}_mpr1T ${subject}_t2w1T -T/data/petsun43/data1/atlas/TRIO_Y_NDC
set k = 1
set T2scans = ()
while ( $k <= $T2num )
	set T2scans = ( ${T2scans} ${subject}_t2w${k}T )
	@ k++
end
avgmpr_4dfp ${T2scans} ${subject}_t2w_avgT useold -T/data/petsun43/data1/atlas/TRIO_Y_NDC
t4imgs_4dfp -s ${subject}_t2w_avgT.lst ${subject}_t2w_avgT -O${subject}_t2w1T
cp ${subject}_t2w1T_to_TRIO_Y_NDC_t4 ${subject}_t2w_avgT_to_TRIO_Y_NDC_t4
cp ${subject}_t2w1T_to_${subject}_mpr1T_t4 ${subject}_t2w_avgT_to_${subject}_mpr1T_t4
t4img_4dfp ${subject}_t2w_avgT_to_${subject}_mpr1T_t4 ${subject}_t2w_avgT ${subject}_t2w_avgT_on_${subject}_mpr1T -O${subject}_mpr1T.4dfp.ifh

end
exit

MRA_PROCESS:
################
# Process MRA data
################
set structdir = ${basedir}/${subject}/MRA
pushd ${structdir}

source ${subdir}/Structurals/${subject}.structparams
set MRA_num = ( $#MRA )
echo $MRA_num
ln -s ../T1/${subject}_mpr1T_debias_to_TRIO_Y_NDC_t4 .
foreach e ( img img.rec ifh hdr )
	ln -s ../T1/${subject}_mpr1T_debias.4dfp.$e .
end

# Debias and convert back to 4dfp
set k = 1
while ( $k <= ${MRA_num} )
	niftigz_4dfp -n ${subject}_MRA${k} ${subject}_MRA${k}
	echo ${scriptdir}/apply_debias.csh ${subject}_MRA${k}
	${scriptdir}/apply_debias.csh ${subject}_MRA${k}
	niftigz_4dfp -4 ${subject}_MRA${k}_debias ${subject}_MRA${k}_debias -N
	@ k++
end

# Register MRA to T1
t2w2mpr_4dfp ${subject}_mpr1T_debias ${subject}_MRA1_debias -T/data/petsun43/data1/atlas/TRIO_Y_NDC
set structstring = ()
@ k = 1
while ( $k <= ${MRA_num} )
	set structstring = "$structstring ${subject}_MRA${k}_debias.4dfp.img"
	@ k++
end
avgmpr_4dfp ${structstring} ${subject}_MRA_debias_avg useold -T/data/petsun43/data1/atlas/TRIO_Y_NDC
#t4imgs_4dfp -s ${subject}_MRA_debias_avg.lst ${subject}_MRA_debias_avg -O${subject}_MRA1_debias
t4imgs_4dfp ${subject}_MRA_debias_avg.lst ${subject}_MRA_debias_avg -O${subject}_MRA1_debias # Without cubic spline interpolation
cp ${subject}_MRA1_debias_to_TRIO_Y_NDC_t4 ${subject}_MRA_debias_avg_to_TRIO_Y_NDC_t4
cp ${subject}_MRA1_debias_to_${subject}_mpr1T_debias_t4 ${subject}_MRA_debias_avg_to_${subject}_mpr1T_debias_t4
t4img_4dfp ${subject}_MRA_debias_avg_to_${subject}_mpr1T_debias_t4 ${subject}_MRA_debias_avg ${subject}_MRA_debias_avg_on_${subject}_mpr1T_debias -O${subject}_mpr1T_debias.4dfp.ifh

# Create masked T2
niftigz_4dfp -n ${subject}_MRA_debias_avg ${subject}_MRA_debias_avg
bet2 ${subject}_MRA_debias_avg ${subject}_MRA_debias_avg_bet
niftigz_4dfp -4 ${subject}_MRA_debias_avg_bet ${subject}_MRA_debias_avg_bet

#fslmaths ${outputdir}/${MRAvol}_debias_bet -thr 150 ${outpudir}/${MRAvol}_debias_bet_thresh150
popd
#exit 

MRV_PROCESS:
################
# Process MRV data
################

set MRV_dirs = ( MRV_cor MRV_sag )
source ${subdir}/Structurals/${subject}.structparams
set MRV_num = ( $#MRV_cor $#MRV_sag )
@ c = 1

foreach MRV_dir ( $MRV_dirs )
	set structdir = ${basedir}/${subject}/${MRV_dir}
	pushd ${structdir}
	
	ln -s ../T1/${subject}_mpr1T_to_TRIO_Y_NDC_t4 .
	foreach e ( img img.rec ifh hdr )
		ln -s ../T1/${subject}_mpr1T.4dfp.$e .
	end
#	SKIP:
#	set k = 1
#	set orient = `grep orientation ${subject}_${MRV_dir}${k}.4dfp.ifh | awk -F ":=" '{print $2}'`
#echo $orient
#	if ( $orient == 3 ) then
#		C2T_4dfp ${subject}_${MRV_dir}${k} ${subject}_${MRV_dir}${k}
#	else if ( $orient == 4 ) then
#		S2T_4dfp ${subject}_${MRV_dir}${k} ${subject}_${MRV_dir}${k}
#	endif
#	exit

	# Debias and convert back to 4dfp
	set k = 1
	while ( $k <= ${MRV_num[$c]} )
		niftigz_4dfp -n ${subject}_${MRV_dir}${k} ${subject}_${MRV_dir}${k}
		echo ${scriptdir}/apply_debias.csh ${subject}_${MRV_dir}${k}
		${scriptdir}/apply_debias.csh ${subject}_${MRV_dir}${k}	
		niftigz_4dfp -4 ${subject}_${MRV_dir}${k}_debias ${subject}_${MRV_dir}${k}_debias -N
		@ k++
	end

	# Register MRV to T1
	t2w2mpr_4dfp ${subject}_mpr1T ${subject}_${MRV_dir}1 -T/data/petsun43/data1/atlas/TRIO_Y_NDC
	set structstring = ()
	@ k = 1
	while ( $k <= ${MRV_num[$c]} )
		set structstring = "$structstring ${subject}_${MRV_dir}${k}"
		@ k++
	end
	avgmpr_4dfp ${structstring} ${subject}_${MRV_dir}_avg useold -T/data/petsun43/data1/atlas/TRIO_Y_NDC
	#t4imgs_4dfp -s ${subject}_${MRV_dir}_debias_avg.lst ${subject}_${MRV_dir}_debias_avg -O${subject}_${MRV_dir}1_debias
	t4imgs_4dfp ${subject}_${MRV_dir}_avg.lst ${subject}_${MRV_dir}_avg -O${subject}_${MRV_dir}1 # Without cubic spline interpolation
	cp ${subject}_${MRV_dir}1_to_TRIO_Y_NDC_t4 ${subject}_${MRV_dir}_avg_to_TRIO_Y_NDC_t4
	cp ${subject}_${MRV_dir}1_to_${subject}_mpr1T_t4 ${subject}_${MRV_dir}_avg_to_${subject}_mpr1T_t4
	t4img_4dfp ${subject}_${MRV_dir}_avg_to_${subject}_mpr1T_t4 ${subject}_${MRV_dir}_avg ${subject}_${MRV_dir}_avg_on_${subject}_mpr1T -O${subject}_mpr1T_debias.4dfp.ifh

	# Create masked T2
	niftigz_4dfp -n ${subject}_${MRV_dir}_avg ${subject}_${MRV_dir}_avg
	bet2 ${subject}_${MRV_dir}_avg ${subject}_${MRV_dir}_avg_bet
	niftigz_4dfp -4 ${subject}_${MRV_dir}_avg_bet ${subject}_${MRV_dir}_avg_bet

	#fslmaths ${outputdir}/${${MRV_dir}vol}_debias_bet -thr 150 ${outpudir}/${${MRV_dir}vol}_debias_bet_thresh150
	popd
	@ c++
end
#exit 

FUNC_PARAMS:
################
# create params
################
set sesnums = `cat ${basedir}/${subject}/${subject}_func_sessions.txt`

set runtypes = ( RSFC Memory_word Memory_scene Memory_face Motor Mixed )
set runlengths = ( 818 121 121 121 104 192)
set runnames = ( 1 mem_word mem_scene mem_face motor mixed )
set runnum = $#runtypes

foreach k ( $sesnums )
	
	pushd $subdir/Functionals/$k

	set r = 1
	set fstd_dcms = 
	set irun_label = 
	
	set bolddcm = `cat $k.studies.txt | grep "$runtypes[$r]" | grep " $runlengths[$r]" | awk -F " " '{print $1}'` 
	set bolddcm = `echo $bolddcm | tr -d '\n'`
	if ( $bolddcm > 0) then 
		set fstd_dcms = `echo $bolddcm`
		set irun_label = `echo "${runnames[$r]}"`
	endif
	@ r++

	while ( $r <= $runnum )
		set bolddcm = `cat $k.studies.txt | grep "$runtypes[$r]" | grep " $runlengths[$r]" | awk -F " " '{print $1}'`
		set bolddcm = `echo $bolddcm | tr -d '\n'`
		set boldnum = $#bolddcm
		if ( $boldnum > 1 ) then
			set t = 1
			while ( $t <= $boldnum )
				if ( $bolddcm[$t] > 0) then 
					set fstd_dcms = `echo $fstd_dcms $bolddcm[$t]`
					set irun_label = `echo $irun_label "$runnames[$r]${t}"`
				endif	
				@ t++
			end
		else
			if ( $bolddcm > 0) then 
				set fstd_dcms = `echo $fstd_dcms $bolddcm`
				set irun_label = `echo $irun_label "$runnames[$r]"`
			endif		
		endif
		@ r++
	end
	
	
	echo "set patid    = $k" > $k.params
	echo "set irun     = ($irun_label)" >> $k.params
	echo "set fstd     = ($fstd_dcms)" >> $k.params
	cat $k.studies.txt | awk 'BEGIN{n=0;};$3~/field/{s[n]=$1;n++;}END{printf("set gfm\t\t= (");for(i=0;i<n;i++)printf(" %d",s[i]);printf(")\n");}' >> $k.params
	echo "set boldruns = (1)" > $k.fcparams	
	popd
	
	
end
exit	

#FUNC_PARAMS_old:
################
# create params
################
#set sesnums = `cat ${basedir}/${subject}/${subject}_func_sessions.txt`
#
#foreach k ( $sesnums )
#	
#	pushd $subdir/Functionals/$k
#	
#	set fstd_dcms = 
#	set irun_label = 
#	set restdcm = `cat $k.studies.txt | grep 'RSFC' | grep ' 818' | awk -F " " '{print $1}'` 
#	set restdcm = `echo $restdcm | tr -d '\n'`
#	if ( $restdcm > 0) then 
#		set fstd_dcms = `echo $restdcm`
#		set irun_label = `echo "1"`
#	endif
#
#	set mem_word_dcm = `cat $k.studies.txt | grep 'Memory_word' | grep ' 121' | awk -F " " '{print $1}'`
#	set mem_word_dcm = `echo $mem_word_dcm | tr -d '\n'`
#	set mem_word_num = $#mem_word_dcm
#	if ( $mem_word_num > 1 ) then
#		set t = 1
#		while ( $t <= $mem_word_num )
#			if ( $mem_word_dcm[$t] > 0) then 
#				set fstd_dcms = `echo $fstd_dcms $mem_word_dcm`
#				set irun_label = `echo $irun_label mem_word$t`
#			endif	
#			@ t++
#		end
#	else
#		if ( $mem_word_dcm > 0) then 
#			set fstd_dcms = `echo $fstd_dcms $mem_word_dcm`
#			set irun_label = `echo $irun_label mem_word`
#		endif		
#	endif
#
#	set mem_scene_dcm = `cat $k.studies.txt | grep 'Memory_scene' | grep ' 121' | awk -F " " '{print $1}'` 
#	set mem_scene_dcm = `echo $mem_scene_dcm | tr -d '\n'`
#	set mem_scene_num = $#mem_scene_dcm	
#	if ( $mem_scene_num > 1 ) then
#		set t = 1
#		while ( $t <= $mem_scene_num )
#			if ( $mem_scene_dcm[$t] > 0) then 
#				set fstd_dcms = `echo $fstd_dcms $mem_scene_dcm`
#				set irun_label = `echo $irun_label mem_scene$t`
#			endif	
#			@ t++
#		end
#	else
#		if ( $mem_scene_dcm > 0) then 
#				set fstd_dcms = `echo $fstd_dcms $mem_scene_dcm`
#				set irun_label = `echo $irun_label mem_scene`
#		endif	
#	endif
#
#	set mem_face_dcm = `cat $k.studies.txt | grep 'Memory_face' | grep ' 121' | awk -F " " '{print $1}'`
#	set mem_face_dcm = `echo $mem_face_dcm | tr -d '\n'`
#	set mem_face_num = $#mem_face_dcm
#	if ( $mem_face_num > 1) then
#		set t = 1
#		while ( $t <= $mem_face_num )
#			if ( $mem_face_dcm[$t] > 0) then 
#				set fstd_dcms = `echo $fstd_dcms $mem_face_dcm`
#				set irun_label = `echo $irun_label mem_face$t`
#			endif
#			@ t++
#		end
#	else
#		if ( $mem_face_dcm > 0) then 
#				set fstd_dcms = `echo $fstd_dcms $mem_face_dcm`
#				set irun_label = `echo $irun_label mem_face`
#		endif
#	endif
#
#	set motordcm = `cat $k.studies.txt | grep 'Motor' | grep ' 104' | awk -F " " '{print $1}'`
#	set motordcm = `echo $motordcm | tr -d '\n'`
#	set motornum = $#motordcm
#	set t = 1
#	while ( $t <= $motornum )
#		if ( $motordcm[$t] > 0 ) then 
#			set fstd_dcms = `echo $fstd_dcms $motordcm[$t]`
#			set irun_label = `echo $irun_label motor$t`
#		endif
#		@ t++
#	end
#
#	set mixeddcm = `cat $k.studies.txt | grep 'Mixed' | grep ' 192' | awk -F " " '{print $1}'`
#	set mixeddcm = `echo $mixeddcm | tr -d '\n'`
#	set mixednum = $#mixeddcm
#	set t = 1
#	while ( $t <= $mixednum )
#		if ( $mixeddcm[$t] > 0 ) then 
#			set fstd_dcms = `echo $fstd_dcms $mixeddcm[$t]`
#			set irun_label = `echo $irun_label mixed$t`
#		endif
#		@ t++
#	end
#
#	echo "set patid    = $k" > $k.params
#	echo "set irun     = ($irun_label)" >> $k.params
#	echo "set fstd     = ($fstd_dcms)" >> $k.params
#	#echo "set gfm      = ($gfm_runs)" >> test_sub0$k.params
#	cat $k.studies.txt | awk 'BEGIN{n=0;};$3~/field/{s[n]=$1;n++;}END{printf("set gfm\t\t= (");for(i=0;i<n;i++)printf(" %d",s[i]);printf(")\n");}' >> $k.params
#	echo "set boldruns = (1)" > $k.fcparams	
#	popd
#	
#	
#end
#exit	

CREATE_ATLAS_LINK:
###############################################
# Generic preprocessing for dcm_to_4dfp etc...
###############################################


foreach k ( $sesnums )
	pushd $subdir/Functionals/$k	
	mkdir atlas
	pushd atlas
	ln -sf $basedir/${subject}/T1/${subject}_mpr1T_debias_to_TRIO_Y_NDC_t4 ./${k}_mpr1_to_TRIO_Y_NDC_t4
	ln -sf $basedir/${subject}/T2/${subject}_t2w_debias_avgT_to_${subject}_mpr1T_debias_t4 ./${k}_t2w_to_${k}_mpr1_t4
	ln -sf $basedir/${subject}/T2/${subject}_t2w_debias_avgT_to_TRIO_Y_NDC_t4 ./${k}_t2w_to_TRIO_Y_NDC_t4
#	ln -s $basedir/${subject}/T2/${subject}_t2w1_to_TRIO_Y_NDC_t4 ./${k}_t2w_to_TRIO_Y_NDC_t4
	foreach e ( img ifh img.rec hdr )
		ln -sf $basedir/${subject}/T1/${subject}_mpr_debias_avgT_111_t88.4dfp.$e ./${k}_mpr_n1_111_t88.4dfp.$e
		ln -sf $basedir/${subject}/T1/${subject}_mpr_debias_avgT_222_t88.4dfp.$e ./${k}_mpr_n1_222_t88.4dfp.$e
		ln -sf $basedir/${subject}/T1/${subject}_mpr_debias_avgT_333_t88.4dfp.$e ./${k}_mpr_n1_333_t88.4dfp.$e
		ln -sf $basedir/${subject}/T2/${subject}_t2w_debias_avgT_111_t88.4dfp.$e ./${k}_t2w_111_t88.4dfp.$e
		ln -sf $basedir/${subject}/T2/${subject}_t2w_debias_avgT_222_t88.4dfp.$e ./${k}_t2w_222_t88.4dfp.$e
		ln -sf $basedir/${subject}/T2/${subject}_t2w_debias_avgT_333_t88.4dfp.$e ./${k}_t2w_333_t88.4dfp.$e
		ln -sf $basedir/${subject}/T1/${subject}_mpr_debias_avgT.4dfp.$e ./${k}_mpr1.4dfp.$e
		ln -sf $basedir/${subject}/T2/${subject}_t2w_debias_avgT.4dfp.$e ./${k}_t2w.4dfp.$e
	end
	popd	
	popd
	#@ k++
end
#exit

GENERIC_PREPROCESS:
###############################################
# Generic preprocessing for dcm_to_4dfp etc...
###############################################

foreach k ( $sesnums )
	
	pushd $subdir/Functionals/$k
	
	#$RELEASE/generic_cross_bold_pp_090115 $k.params ${basedir}/MSC.params	

	/net/nil-bluearc/GMT/Laumann/MSC/Scripts/MSC_cross_bold_pp_130702.csh ${k}.params $basedir/MSC.params	

	popd
	
end
exit

RUN_DVAR_4dfp:
#######################################
# run_dvar_4dfp individually on each run
#######################################

foreach k ( $sesnums )

	set patid = $k
	pushd ${subdir}/Functionals/${patid}
	source ${patid}.params
	foreach	r ( $irun )
		pushd ./bold$r/
		echo ${patid}_b${r}_faln_dbnd_xr3d_norm > ${patid}_b${r}.lst
		conc_4dfp ${patid}_b${r}_faln_dbnd_xr3d_norm -l${patid}_b${r}.lst
		run_dvar_4dfp ${patid}_b${r}_faln_dbnd_xr3d_norm.conc -m../atlas/${patid}_func_vols_avez -n0 -b10 -x8
		popd
	end
	popd
end
#exit






APPLY_REGISTER_UWRP_FIRST_SESSION:
###########################################################
# Register func vol ave unwarp to first session epi and resample bold data
###########################################################
set T_epi      = ${subdir}/Functionals/$sesnums[1]/unwarp/$sesnums[1]_func_vols_ave_uwrp
set T_epi_mask = ${subdir}/Functionals/$sesnums[1]/unwarp/$sesnums[1]_func_vols_ave_uwrp_mskt
set U          = ${subdir}/Functionals/$sesnums[1]/unwarp/$sesnums[1]_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4
source ${basedir}/MSC.params

# generate mask for first session
pushd ${subdir}/Functionals/$sesnums[1]/unwarp
echo msktgen_4dfp $sesnums[1]_func_vols_ave_uwrp_mskt -T$REFDIR/TRIO_Y_NDC
msktgen_4dfp $sesnums[1]_func_vols_ave_uwrp -T$REFDIR/TRIO_Y_NDC
popd

# remake single resampled 333 atlas space fMRI volumetric timeseries for first session
set patid = $sesnums[1]
pushd ${subdir}/Functionals/${patid}
source ${patid}.params
$rsam_cmnd ${patid}.params ${basedir}/MSC.params

set MBstr = _faln_dbnd
set lst = ${patid}${MBstr}_xr3d_uwrp_atl.lst
if (-e $lst) /bin/rm $lst
touch $lst
@ k = 1
while ($k <= $#irun)
	echo bold$irun[$k]/${patid}_b$irun[$k]${MBstr}_xr3d_uwrp_atl.4dfp.img >> $lst
	@ k++
end
conc_4dfp ${lst:r}.conc -l$lst
if ($status) exit $status
set format = `cat atlas/${patid}_func_vols.format`
if ($status) exit $status
actmapf_4dfp $format	${patid}${MBstr}_xr3d_uwrp_atl.conc -aave
if ($status) exit $status
ifh2hdr -r2000 		${patid}${MBstr}_xr3d_uwrp_atl_ave
mv			${patid}${MBstr}_xr3d_uwrp_atl_ave.4dfp.*	unwarp
var_4dfp -sf$format	${patid}${MBstr}_xr3d_uwrp_atl.conc
ifh2hdr -r20		${patid}${MBstr}_xr3d_uwrp_atl_sd1
mv			${patid}${MBstr}_xr3d_uwrp_atl_sd1*		unwarp
mv			${patid}${MBstr}_xr3d_uwrp_atl.conc*		unwarp

# Register epi to first session epi, and resample BOLD to atlas
set modes = (0 0 0 0)
@ modes[1] = 2048 + 3 + 256 
@ modes[2] = 2048 + 3 + 256 + 4
@ modes[3] = 2048 + 3 + 256 + 4
@ modes[4] = $modes[3]

@ n = $#sesnums
@ i = 2
while ( $i <= $n )
	set patid = $sesnums[$i]
	cd ${subdir}/Functionals/${patid}
	source ${patid}.params
	pushd unwarp	# into unwarp
	set t4file = ${patid}_func_vols_ave_uwrp_to_$sesnums[1]_func_vols_ave_uwrp_t4
	if ($status) exit $status
	set log =    ${patid}_func_vols_ave_uwrp_to_$sesnums[1]_func_vols_ave_uwrp.log
	date >! $log
	@ k = 1
	while ($k <= $#modes)
	echo	imgreg_4dfp ${T_epi} ${T_epi_mask} ${patid}_func_vols_ave_uwrp none $t4file $modes[$k] >> $log
		imgreg_4dfp ${T_epi} ${T_epi_mask} ${patid}_func_vols_ave_uwrp none $t4file $modes[$k] >> $log
		if ($status) exit $status
		@ k++
	end
	t4_mul $t4file $U ${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4
	t4img_4dfp ${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4 ${patid}_func_vols_ave_uwrp ${patid}_func_vols_ave_uwrp_on_TRIO_Y_NDC_111 -O111
	t4img_4dfp ${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4 ${patid}_func_vols_ave_uwrp ${patid}_func_vols_ave_uwrp_on_TRIO_Y_NDC_333 -O333
	if ($status) exit $status
	popd		# out of unwarp
	$rsam_cmnd ${patid}.params ${basedir}/MSC.params

	# remake single resampled 333 atlas space fMRI volumetric timeseries	
	set MBstr = _faln_dbnd
	set lst = ${patid}${MBstr}_xr3d_uwrp_atl.lst
	if (-e $lst) /bin/rm $lst
	touch $lst
	@ k = 1
	while ($k <= $#irun)
		echo bold$irun[$k]/${patid}_b$irun[$k]${MBstr}_xr3d_uwrp_atl.4dfp.img >> $lst
		@ k++
	end
	conc_4dfp ${lst:r}.conc -l$lst
	if ($status) exit $status
	set format = `cat atlas/${patid}_func_vols.format`
	if ($status) exit $status
	actmapf_4dfp $format	${patid}${MBstr}_xr3d_uwrp_atl.conc -aave
	if ($status) exit $status
	ifh2hdr -r2000 		${patid}${MBstr}_xr3d_uwrp_atl_ave
	mv			${patid}${MBstr}_xr3d_uwrp_atl_ave.4dfp.*	unwarp
	var_4dfp -sf$format	${patid}${MBstr}_xr3d_uwrp_atl.conc
	ifh2hdr -r20		${patid}${MBstr}_xr3d_uwrp_atl_sd1
	mv			${patid}${MBstr}_xr3d_uwrp_atl_sd1*		unwarp
	mv			${patid}${MBstr}_xr3d_uwrp_atl.conc*		unwarp

	@ i++
end
#exit

MEAN_FIELD_MAP_MAKER:
###########################################################
# Make mean field map
###########################################################
${scriptdir}/meanfield_maker_MSC.csh ${subject}
exit

APPLY_REGISTER_UWRPMEAN_FIRST_SESSION:
###########################################################
# Apply mean distortion correction, register func vol ave unwarp to first session and resample BOLD data
###########################################################
source ${basedir}/MSC.params
set FMmean	= $basedir/${subject}/meanfield/phase_rad_unwrap_secs_resolved_on_TRIO_Y_NDC_111

# Apply mean distortion correction to first session, register to t2w, and resample BOLD to atlas
set patid = $sesnums[1]
pushd ${subdir}/Functionals/${patid}
source ${patid}.params

set T		= $basedir/${subject}/T2/${subject}_t2w_debias_avgT
set T_mask	= $basedir/${subject}/T2/${subject}_t2w_debias_avgT_bet
set U		= $basedir/${subject}/T2/${subject}_t2w_debias_avgT_to_TRIO_Y_NDC_t4
$uwrp_cmnd -mean atlas/${patid}_func_vols_ave $FMmean atlas/${patid}_func_vols_ave_to_TRIO_Y_NDC_t4 ${dwell} ${ped} 
if ($status) exit $status
pushd unwarp	# into unwarp
set t4file = ${patid}_func_vols_ave_uwrpmean_to_${patid}_t2w_t4
cp  ../atlas/${patid}_func_vols_ave_to_${patid}_t2w_t4 $t4file
if ($status) exit $status
set log =    ${patid}_func_vols_ave_uwrpmean_to_t2w_avgT.log
date >! $log
set modes = (0 0 0 0)
@ modes[1] = 3072 + 3 + 256
@ modes[2] = 2048 + 3 + 256
@ modes[3] = 2048 + 3 + 256
@ modes[4] = $modes[3]
@ j = 1
while ($j <= $#modes)
echo	imgreg_4dfp $T $T_mask ${patid}_func_vols_ave_uwrp none $t4file $modes[$j] >> $log
	imgreg_4dfp $T $T_mask ${patid}_func_vols_ave_uwrp none $t4file $modes[$j] >> $log
	if ($status) exit $status
	@ j++
end
t4_mul $t4file $U 	${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4
t4img_4dfp 		${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4 ${patid}_func_vols_ave_uwrp ${patid}_func_vols_ave_uwrp_on_TRIO_Y_NDC_333 -O333
t4img_4dfp 		${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4 ${patid}_func_vols_ave_uwrp ${patid}_func_vols_ave_uwrp_on_TRIO_Y_NDC_111 -O111
msktgen_4dfp ${patid}_func_vols_ave_uwrp -T$REFDIR/TRIO_Y_NDC # Create masked func_vols_ave_unwarp for subsequent registration
if ($status) exit $status
popd		# out of unwarp
$rsam_cmnd ${patid}.params ${basedir}/MSC.params
if ($status) exit $status
if ( -d unwarp_mean ) /bin/rm -r unwarp_mean
mv unwarp unwarp_mean

# remake single resampled 333 atlas space fMRI volumetric timeseries	
set MBstr = _faln_dbnd
set lst = ${patid}${MBstr}_xr3d_uwrp_atl.lst
if (-e $lst) /bin/rm $lst
touch $lst
@ k = 1
while ($k <= $#irun)
	echo bold$irun[$k]/${patid}_b$irun[$k]${MBstr}_xr3d_uwrp_atl.4dfp.img >> $lst
	@ k++
end
conc_4dfp ${lst:r}.conc -l$lst
if ($status) exit $status
set format = `cat atlas/${patid}_func_vols.format`
if ($status) exit $status
actmapf_4dfp $format	${patid}${MBstr}_xr3d_uwrp_atl.conc -aave
if ($status) exit $status
ifh2hdr -r2000 		${patid}${MBstr}_xr3d_uwrp_atl_ave
mv			${patid}${MBstr}_xr3d_uwrp_atl_ave.4dfp.*	unwarp_mean
var_4dfp -sf$format	${patid}${MBstr}_xr3d_uwrp_atl.conc
ifh2hdr -r20		${patid}${MBstr}_xr3d_uwrp_atl_sd1
mv			${patid}${MBstr}_xr3d_uwrp_atl_sd1*		unwarp_mean
mv			${patid}${MBstr}_xr3d_uwrp_atl.conc*		unwarp_mean

## Apply mean distortion correction to each session, register epi to first session epi, and resample BOLD to atlas
set T_epi      = ${subdir}/Functionals/$sesnums[1]/unwarp_mean/$sesnums[1]_func_vols_ave_uwrp
set T_epi_mask = ${subdir}/Functionals/$sesnums[1]/unwarp_mean/$sesnums[1]_func_vols_ave_uwrp_mskt
set U          = ${subdir}/Functionals/$sesnums[1]/unwarp_mean/$sesnums[1]_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4

set modes = (0 0 0 0)
@ modes[1] = 2048 + 3 + 256 
@ modes[2] = 2048 + 3 + 256 + 4
@ modes[3] = 2048 + 3 + 256 + 4
@ modes[4] = $modes[3]

@ n = $#sesnums
@ i = 2
while ( $i <= $n )
	set patid = $sesnums[$i]
	cd ${subdir}/Functionals/${patid}
	source ${patid}.params
	$uwrp_cmnd -mean atlas/${patid}_func_vols_ave $FMmean atlas/${patid}_func_vols_ave_to_TRIO_Y_NDC_t4 ${dwell} ${ped} 
	if ($status) exit $status
	pushd unwarp	# into unwarp
	set t4file = ${patid}_func_vols_ave_uwrp_to_$sesnums[1]_func_vols_ave_uwrp_t4
	if ($status) exit $status
	set log =    ${patid}_func_vols_ave_uwrpmean_to_$sesnums[1]_func_vols_ave_uwrp.log
	date >! $log
	@ k = 1
	while ($k <= $#modes)
	echo	imgreg_4dfp ${T_epi} ${T_epi_mask} ${patid}_func_vols_ave_uwrp none $t4file $modes[$k] >> $log
		imgreg_4dfp ${T_epi} ${T_epi_mask} ${patid}_func_vols_ave_uwrp none $t4file $modes[$k] >> $log
		if ($status) exit $status
		@ k++
	end
	t4_mul $t4file $U ${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4
	t4img_4dfp ${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4 ${patid}_func_vols_ave_uwrp ${patid}_func_vols_ave_uwrp_on_TRIO_Y_NDC_111 -O111
	t4img_4dfp ${patid}_func_vols_ave_uwrp_to_TRIO_Y_NDC_t4 ${patid}_func_vols_ave_uwrp ${patid}_func_vols_ave_uwrp_on_TRIO_Y_NDC_333 -O333
	if ($status) exit $status
	popd		# out of unwarp
	$rsam_cmnd ${patid}.params ${basedir}/MSC.params
	if ($status) exit $status
	if ( -d unwarp_mean ) /bin/rm -r unwarp_mean
	mv unwarp unwarp_mean

	# remake single resampled 333 atlas space fMRI volumetric timeseries	
	set MBstr = _faln_dbnd
	set lst = ${patid}${MBstr}_xr3d_uwrp_atl.lst
	if (-e $lst) /bin/rm $lst
	touch $lst
	@ k = 1
	while ($k <= $#irun)
		echo bold$irun[$k]/${patid}_b$irun[$k]${MBstr}_xr3d_uwrp_atl.4dfp.img >> $lst
		@ k++
	end
	conc_4dfp ${lst:r}.conc -l$lst
	if ($status) exit $status
	set format = `cat atlas/${patid}_func_vols.format`
	if ($status) exit $status
	actmapf_4dfp $format	${patid}${MBstr}_xr3d_uwrp_atl.conc -aave
	if ($status) exit $status
	ifh2hdr -r2000 		${patid}${MBstr}_xr3d_uwrp_atl_ave
	mv			${patid}${MBstr}_xr3d_uwrp_atl_ave.4dfp.*	unwarp_mean
	var_4dfp -sf$format	${patid}${MBstr}_xr3d_uwrp_atl.conc
	ifh2hdr -r20		${patid}${MBstr}_xr3d_uwrp_atl_sd1
	mv			${patid}${MBstr}_xr3d_uwrp_atl_sd1*		unwarp_mean
	mv			${patid}${MBstr}_xr3d_uwrp_atl.conc*		unwarp_mean

	@ i++
end
exit

SURFACE_CREATION:
###########################################################
# Surface processing and post freesurfer pipeline
###########################################################
set freesurfbin = /data/gizmo/data1/freesurfer5.3/bin
set freesurfdir = /net/nil-bluearc/GMT/Laumann/MSC/fs5.3_native_default
${freesurfbin}/recon-all -all -sd ${freesurfdir} -s ${subject} -i ${basedir}/${subject}/T1/${subject}_mpr_debias_avgT.nii.gz # -custom-tal-atlas TRIO_Y_NDC_as_mni_average_305
exit

CREATE_RIBBON:
###########################################################
# Create ribbon and create masks of voxels with high variability
###########################################################

# Resample structural volume for ribbon creation
pushd ${basedir}/${subject}/T1
echo niftigz_4dfp -n ${subject}_mpr_debias_avgT_333_t88 ${subject}_mpr_debias_avgT_333_t88
niftigz_4dfp -n ${subject}_mpr_debias_avgT_333_t88 ${subject}_mpr_debias_avgT_333_t88
#rm ${subject}_mpr_debias_avgT_333_t88.4dfp.* ${subject}_mpr_debias_avgT_111_t88.4dfp.*

# Create ribbon in 333 space
${scriptdir}/create_ribbon.csh ${subject}

# Create masks of variable voxels for each session
${scriptdir}/RibbonVolumetoSurfaceMapping_MSC_333.csh ${subject}
popd
end
exit

