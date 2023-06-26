# Run the pipeline on slurm using singularity
snakemake --unlock  --rerun-incomplete --cores 10 --local-cores 10
snakemake \
--snakefile Snakefile \
--configfile config.yaml \
--printshellcmds \
--use-singularity \
--singularity-args "--bind '$PWD/../','/scicore/home/zavolan/GROUP/StefanieCLIP/input'" \
--cluster-config cluster.json \
--cores 500 \
--local-cores 10 \
--jobs 20 \
--cluster "sbatch \
	--cpus-per-task={cluster.threads} \
	--mem={cluster.mem} \
	--qos={cluster.queue} \
	--time={cluster.time} \
	--output={cluster.out}" \
&>> bclip_integrator_complex.log
