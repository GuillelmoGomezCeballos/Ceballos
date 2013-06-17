awk '{printf("hadd histo_%s_all_noskim.root histo_%s_*_noskim.root\n",$5,$5)}' submitter_all.csh > merge_histo.csh 
