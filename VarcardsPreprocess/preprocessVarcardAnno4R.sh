
#input Varcard annotation file, header NOT added
i=$1
#arguement 2, group name, can be Nomarl or PID
group=$2
#1. add group info, 2. remove useless fields, 
#3. remove Ensembl ID etc useless gene IDs, leave only gene symbl
#4. for gene symbols like "RAG1;RAG1", clean it, clean only the first one
#5. remove (%xxx.xxx) in score values
cat header $i|awk -v group=$group '{if(NR==1){print $0"\tGroup"}else{print $0"\t"group}}'|cut -f-5,8-|awk '{print $1"_"$2"\t"$0}'|cut -f1,5-|perl -pe 's/\|([\S]+)\t/\t/'|cut -f-6,8-|perl -nle 'my @F = split /\t/;$F[3] = (split /[\;\,]/,$F[3])[0];my $out = join "\t", @F;print "$out";'|perl -pe 's/\(\%[\d\.]+\)//g' >$i.CLIN_SIG;
#remove duplicate records
perl getUniq.varCards4R.pl $i.CLIN_SIG
