date > log.txt;
indir="../gbff";
outdir="../antismash_output_strict_known_clust";
counter=0;
date >> log.txt;
for file in `ls ../gbff/`;do 
    echo "Initializing $file" >> log.txt;
    if [ $counter -eq 7 ]; then
        sleep 7m;
        counter=0;
    fi
    counter=$((counter + 1));
    screen -dm run_antismash $indir/$file $outdir --cpu 5 --hmmdetection-strictness strict --genefinding-tool prodigal --cb-knownclusters
    #run_antismash $indir/$file $outdir -v --cpu 15 --hmmdetection-strictness strict --genefinding-tool prodigal --cb-knownclusters
    date >> log.txt;
done
date >> log.txt;

