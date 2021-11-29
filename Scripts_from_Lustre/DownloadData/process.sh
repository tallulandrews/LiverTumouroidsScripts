
for i in `seq 2 54`;
do 
	bsub -R"select[mem>10000] rusage[mem=10000]" -M10000 -o buildsce.%J.out -e buildsce.%J.err /software/R-3.4.2/bin/Rscript Process_data.R $i
done
