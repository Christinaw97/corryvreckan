executable=../bin/corry
prealignment_conf=mimosa_prealignment.conf
alignment_conf=mimosa_alignment.conf
prune="-o AlignmentTrackChi2.prune_tracks=true -o AlignmentTrackChi2.max_track_chi2ndof=10"
#${executable} -c ${prealignment_conf} -o Tracking4D.spatial_cut_abs=1000um,1000um -o Tracking4D.min_hits_on_track=3 -o histogram_file="out_MIMOSA.root"
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=1000um,1000um -o Tracking4D.min_hits_on_track=3 -o histogram_file="out_MIMOSA_1.root"
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=1000um,1000um -o Tracking4D.min_hits_on_track=3 -o histogram_file="out_MIMOSA_2.root"
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=1000um,1000um -o Tracking4D.min_hits_on_track=3 -o histogram_file="out_MIMOSA_3.root"
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=1000um,1000um -o Tracking4D.min_hits_on_track=3 -o histogram_file="out_MIMOSA_4.root"
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=1000um,1000um -o Tracking4D.min_hits_on_track=3 -o histogram_file="out_MIMOSA_5.root"
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=1000um,1000um -o Tracking4D.min_hits_on_track=3 -o histogram_file="out_MIMOSA_6.root"
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=500um,500um -o Tracking4D.min_hits_on_track=3 -o histogram_file="out_MIMOSA_7.root"
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=500um,500um -o Tracking4D.min_hits_on_track=3 -o histogram_file="out_MIMOSA_8.root"
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=200um,200um -o Tracking4D.min_hits_on_track=3 -o histogram_file="out_MIMOSA_9.root"
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=200um,200um -o Tracking4D.min_hits_on_track=3 -o histogram_file="out_MIMOSA_10.root"
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=200um,200um -o Tracking4D.min_hits_on_track=6 -o histogram_file="out_MIMOSA_11.root"
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=200um,200um -o Tracking4D.min_hits_on_track=6 -o histogram_file="out_MIMOSA_12.root"
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=200um,200um -o Tracking4D.min_hits_on_track=6 -o histogram_file="out_MIMOSA_13.root" ${prune}
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=200um,200um -o Tracking4D.min_hits_on_track=6 -o histogram_file="out_MIMOSA_14.root" ${prune}
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=150um,150um -o Tracking4D.min_hits_on_track=6 -o histogram_file="out_MIMOSA_15.root" ${prune}
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=150um,150um -o Tracking4D.min_hits_on_track=6 -o histogram_file="out_MIMOSA_16.root" ${prune}
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=100um,100um -o Tracking4D.min_hits_on_track=6 -o histogram_file="out_MIMOSA_17.root" ${prune}
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=100um,100um -o Tracking4D.min_hits_on_track=6 -o histogram_file="out_MIMOSA_18.root" ${prune}

prune="-o AlignmentTrackChi2.prune_tracks=true -o AlignmentTrackChi2.max_track_chi2ndof=5"
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=50um,50um -o Tracking4D.min_hits_on_track=6 -o histogram_file="out_MIMOSA_19.root" ${prune}
${executable} -c ${alignment_conf} -o Tracking4D.spatial_cut_abs=50um,50um -o Tracking4D.min_hits_on_track=6 -o histogram_file="out_MIMOSA_20.root" ${prune}
