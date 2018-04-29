with(colData(BS.unstranded)[,1:2], message(sprintf("%s\t%s\n",Sample,Infection))) 

message(sprintf("%s\n", gsub("^(.*)_R1_merged_val_1.fq.gz_(bismark_bt2_pe).deduplicated.(CpG_report.txt.gz)", "\\1_\\2_\\3", list.files("extract_refined/Merged/"))))

colData(BS.unstranded)



finalFilenames <- colData(BS.unstranded)[,"Filename"]
message(sprintf("%s\n", gsub("^.*_([ATCG]{6})_.*","\\1", finalFilenames)))
