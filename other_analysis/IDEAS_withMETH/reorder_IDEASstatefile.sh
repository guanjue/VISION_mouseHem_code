time tail -n+2 run_IDEAS.state > run_IDEAS.noheader.state

time head -1 run_IDEAS.state > run_IDEAS.reordered.state

time paste ~/group/projects/vision/withMETH/check_order/bin.tab.mm10.200bp.noblacklist.bed.output.bed run_IDEAS.noheader.state | awk -F '\t' -v OFS=' ' '{print $4,$1,$2,$3,$5}' | awk -F ' ' -v OFS=' ' '{print $1,$2,$3,$4, $9,$10,$11,$12,$13,$14,$15,$16}' | awk -F 'R' -v OFS=' ' '{print $2}' | sort -k1,1n |  awk -F '\t' -v OFS=' ' '{print "R"$1}' >> run_IDEAS.reordered.state


