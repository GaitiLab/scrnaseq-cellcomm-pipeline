# Notes

```sh
docker build --tag cell2cell . --platform linux/amd64

docker run --rm -ti --platform linux/amd64 cell2cell python3

```

```sh
docker build --tag cell2cell . --platform linux/amd64

docker run --rm -ti --platform linux/amd64 cpdb python3

```


```sh

nextflow run ${PWD} -profile conda,slurm -params-file "nf-params.yml" --outdir "output" -c "gaitilab.config"
```

****
```sh 
export APPTAINERENV_APPEND_PATH="/opt/conda/bin:/mnt/bin"

apptainer run --pwd ${PWD} --no-home --cleanenv -c --bind ${PWD}:/mnt containers/scrnaseqcellcomm.sif 000_get_metadata.R --input_file "example_data.rds"
```

^\S+\.csv$