# SV-Det
sh a01multi.SV.Call.sh [-h] [-s] <FI: ref.fa>  <FI: tumor.bam> <FI:normal.bam>
perl a02multi.SV.filt.pl -i <FI: SV vcf> -o <STR:output directory> -t <STR:tool >  -sam <STR:sample >
perl a03multi.SV.overlap.pl  -i <input directory> -o <STR:out.dir>  -sam <STR:sample> -n <STR: overlap tools number >  -bp <STR: distance to consider one SV>
