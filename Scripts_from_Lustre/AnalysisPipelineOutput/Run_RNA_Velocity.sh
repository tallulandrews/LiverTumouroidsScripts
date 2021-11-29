bsub -R"select[mem>80000] rusage[mem=80000] span[hosts=1]" -M80000 -n5 -o velotest.out -e velotest.err /software/R-3.4.2/bin/Rscript /nfs/users/nfs_t/ta6/Collaborations/LiverOrganoids/Laura_Pipeline/RNA_velocity.R SecondTrio

