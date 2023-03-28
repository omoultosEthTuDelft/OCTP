CC = gcc
TARGET = frdf frdf_original frdf_omp

DUMPFILE = frdf.dump
LOGFILE = frdf.log
OUTPUTFILE = frdf.output
MINTIMESTEP = 0.0001

all: $(TARGET)

run: frdf
	./frdf $(DUMPFILE) $(LOGFILE) $(OUTPUTFILE) $(MINTIMESTEP)

frdf: force_rdf.c
	$(CC) -g -Wall force_rdf.c -lm -o frdf -Wno-unused-result -O3

frdf_omp: force_rdf_omp.c
	$(CC) -g -Wall force_rdf_omp.c -lm -o frdf_omp -Wno-unused-result -O3 -fopenmp

frdf_original: force_rdf_original.c
	$(CC) -Wall force_rdf_original.c -lm -o frdf_original -Wno-unused-result -O3

# Path: makefile
clean:
	rm -f $(TARGET)
