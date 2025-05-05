# Summary
This workshop offers a balanced blend of theoretical foundation and hands-on experience in single-cell RNA sequencing (scRNA-seq) data analysis using the Seurat R package. It is designed for researchers new to single-cell data or looking to deepen their practical skills in cell clustering, visualization, and interpretation.

##### The container image file ("scRNA-seq.sif") can be found at: https://drive.google.com/drive/folders/18vXOcOPEPUGW85fM6ZbCxKQJQuHnanQD?usp=sharing

##### test_cluseter.sh

This is the job script we are going to use in the test cluster @spring2025-uofa.c3.ca

    sbatch test_cluster.sh

##### alliance_cluster.sh

This is the job script you can use in the Alliance cluster. Please change the account and server domain before submitting the job.

    sbatch alliance_cluseter.sh

##### local_linux.sh

This is the bash script that you run in your local computer

    cd /path/to/scRNA-seq.sif
    bash local_linux.sh
    
