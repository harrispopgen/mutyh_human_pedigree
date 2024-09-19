

.preparegenome.done: scripts/fixpolfile.sh
	@echo "Preparing genome for mispolarization analysis"
	./scripts/fixpolfile.sh
	@touch .preparegenome.done

myhped_ms6.csv: scripts/script.py .preparegenome.done
	@echo "Running mispolarization analysis"
	./scripts/script.py

all: myhped_ms6.csv
