## PForLA - Path Finding f*or* Long read Assembly

This application is designed and developed for CS557 class in Bilkent University.

- To run the application, first you need to generate a DBG with dbgh5 tools provided with gatb-core library. For example;

```
dbgh5 -in Read1.fq,Read2.fq -kmer-size 19 -out graph -out-dir . -abundance-min 2
```

- Then, using generated graph, you can run long read path finding by

```
python3 graph.py <graph.h5> <long_reads.fasta> <output_file.txt> <kmer_size> <use_first_n_reads>
```

The program will output to stdout; its path coverage and also total length of long read data that it processed.

Found paths will be written in specified output file.