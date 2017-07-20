#!/bin/csh
#set basedir = /data/heisenberg/data1/MSC
set basedir = /net/nil-bluearc/GMT/Laumann/MSC
set subject = $1
set fsdir = /net/nil-bluearc/GMT/Laumann/MSC/fs5.3_native/FREESURFER_fs_LR/${subject}/7112b_fs_LR/Ribbon
set neighsmooth = 5
set factor = 0.5
set sessions = $2 #`cat ${basedir}/${subject}/${subject}_func_sessions.txt` #sub014
set subdir = ${basedir}/${subject}

foreach session ( ${sessions} )
	
	set sesdir = ${subdir}/Functionals/${session}
	pushd ${sesdir}
	source ${session}.params	

	foreach r ( $irun )
		pushd ${sesdir}/bold${r}
		set preproc_runfunc = ${session}_b${r}_faln_dbnd_xr3d_uwrp_atl
	        set outputdir = ${sesdir}/bold${r}/goodvoxels_indiv
		echo ${outputdir}
		mkdir ${outputdir}
		set format = `cat ${session}_b${r}_faln_dbnd_xr3d_norm.format`
		actmapf_4dfp "${format}" -amean ${preproc_runfunc}
        	var_4dfp -f"${format}" -s ${preproc_runfunc} 
		mv ${preproc_runfunc}_mean.4dfp.* ${outputdir}
		mv ${preproc_runfunc}_sd1.4dfp.* ${outputdir}

		pushd ${outputdir}
		niftigz_4dfp -n ${preproc_runfunc}_mean ${session}_mean
		niftigz_4dfp -n ${preproc_runfunc}_sd1 ${session}_sd1
		rm ${preproc_runfunc}_mean.4dfp.*
		rm ${preproc_runfunc}_sd1.4dfp.*

	
		fslmaths ${outputdir}/${session}_sd1 -div ${outputdir}/${session}_mean ${outputdir}/${session}_cov
	
		fslmaths ${outputdir}/${session}_cov -mas ${fsdir}/ribbon_333.nii.gz ${outputdir}/${session}_cov_ribbon
	
		fslmaths ${outputdir}/${session}_cov_ribbon -div `fslstats ${outputdir}/${session}_cov_ribbon -M` ${outputdir}/${session}_cov_ribbon_norm
		fslmaths ${outputdir}/${session}_cov_ribbon_norm -bin -s $neighsmooth ${outputdir}/${session}_SmoothNorm
		fslmaths ${outputdir}/${session}_cov_ribbon_norm -s $neighsmooth -div ${outputdir}/${session}_SmoothNorm -dilD ${outputdir}/${session}_cov_ribbon_norm_s${neighsmooth}
		fslmaths ${outputdir}/${session}_cov -div `fslstats ${outputdir}/${session}_cov_ribbon -M` -div ${outputdir}/${session}_cov_ribbon_norm_s${neighsmooth} -uthr 1000 ${outputdir}/${session}_cov_norm_modulate
		fslmaths ${outputdir}/${session}_cov_norm_modulate -mas ${fsdir}/ribbon_333.nii.gz ${outputdir}/${session}_cov_norm_modulate_ribbon
	
		set STD = `fslstats ${outputdir}/${session}_cov_norm_modulate_ribbon -S`
		set MEAN = `fslstats ${outputdir}/${session}_cov_norm_modulate_ribbon -M`
	
		set Lower = `echo "${MEAN} - (${STD} * ${factor})" | bc -l`
		
		set Upper = `echo "${MEAN} + (${STD} * ${factor})" | bc -l`
	
		fslmaths ${outputdir}/${session}_mean -bin ${outputdir}/${session}_mask
		fslmaths ${outputdir}/${session}_cov_norm_modulate -thr $Upper -bin -sub ${outputdir}/${session}_mask -mul -1 ${outputdir}/${session}_goodvoxels
		popd
		popd
	end
	popd
end
