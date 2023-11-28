

config.json: Snakegg.py samples_PED.csv .env
	python $<

.pipe.done:  Snakefile config.json
	snakemake --profile dmp -j 50

all: .pipe.done

