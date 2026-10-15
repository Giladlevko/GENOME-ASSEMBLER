De Bruijn Graph Genome Assembler
A high-performance genome assembler written in C++, designed to reconstruct genomic sequences from short paired-end reads using de Bruijn graphs. The assembler supports error-free and error-prone reads, read preprocessing, graph simplification, and repeat-aware contig stitching to produce highly contiguous genome assemblies.
🚀 Key Features
- Read Preprocessing: Cleans sequencing reads before graph construction to reduce sequencing errors and improve graph quality.
- De Bruijn Graph Construction: Builds directed de Bruijn graphs from short sequencing reads using k-mers.
- Graph Simplification:
  - Tip Removal: Prunes erroneous dead-end branches introduced by sequencing errors.
  - Bubble Detection: Identifies and resolves local graph bubbles caused by sequencing errors or sequence variation.
- Repeat-Aware Multi-Mapping: Allows graph nodes to map to multiple contigs using a std::vector of contig indices, preserving connections through repetitive regions.
- Contig Stitching and Scaffolding: Validates contig-link pairs and combines compatible paths while enforcing structural constraints to reduce the risk of misassemblies.
- Synthetic Dataset Generation: Includes a C++ test generator for simulating sequencing reads and evaluating assembly performance.
🧬 Assembly Pipeline
The assembler processes sequencing data through several stages:
1. Read Cleaning: Preprocesses sequencing reads to reduce errors before graph construction.
2. De Bruijn Graph Construction: Builds a directed graph from the cleaned reads using k-mers.
3. Graph Simplification: Removes erroneous tips and resolves bubbles to reduce graph complexity.
4. Repeat-Aware Traversal: Preserves multiple contig associations at repeat nodes to avoid unnecessarily fragmenting the assembly.
5. Contig Stitching: Validates and connects compatible contig links to produce longer assembled sequences.
📊 Benchmark Results
Tested on bacterial reference genomes, including Escherichia coli and N. deltocephalinicola, under varying read-length and sequencing-error conditions.
Metric	Reported Result
E. coli NGA50	Approximately 133.6 kb
Genome fraction	Approximately 99%
Read lengths	100 bp
Simulated substitution error rate	1%
Simulated sequencing coverage	30×


The assembler is designed to maintain assembly continuity while limiting structural errors across datasets with different read lengths and error profiles.
Benchmark values should be interpreted in the context of the reference genome, dataset, and evaluation methodology used.
🔍 Technical Highlight: Repeat-Aware Graph Traversal
Repetitive genomic sequences present a major challenge for de Bruijn graph assemblers. A graph node may correspond to a sequence that occurs in multiple genomic locations, making it difficult to determine which path a traversal should follow.
A strict one-to-one node-to-contig mapping can discard ambiguous repeat nodes or prevent paths from being connected, potentially fragmenting an otherwise contiguous assembly.
This assembler addresses the problem through a multi-mapping index strategy:
1. Multiple Contig Associations: Each graph node can store a std::vector of contig indices rather than a single index.
2. Preserved Repeat Connections: Nodes shared by multiple paths retain their associations instead of being automatically discarded.
3. Constraint-Based Stitching: Traversal and contig stitching use structural constraints to identify compatible connections through repetitive regions.
This approach is intended to preserve graph connectivity and improve contig continuity, including in datasets containing sequencing errors.
🧪 Synthetic Test Generator
The repository includes test_generator.cpp, a C++ utility for generating synthetic sequencing datasets for debugging, validation, and benchmarking.
Generator Capabilities
- Genome source: Loads reference sequences from N. deltocephalinicola or E. coli, or generates random genomes with a target GC content of 50%.
- Read simulation: Generates reads targeting a length of 100 bp, with lengths varying between 50 and 150 bp.
- Coverage: Targets 30× sequencing coverage.
- Error injection: Introduces random nucleotide substitutions at a target rate of 1%.
- Paired-end reads: Generates read pairs using insert distances between 100 and 300 bp.
- Input formatting: Writes paired-end reads using | delimiters and includes fragment-length information.
Generated Files
File	Description
genome_assembly/reference_gene.txt	Ground-truth reference sequence used for evaluation
genome_assembly/test_inputs.txt	Simulated sequencing reads supplied to the assembler


The reference sequence provides a ground truth against which the resulting assembly can be evaluated.
⚙️ Building and Running
Prerequisites
- A C++ compiler supporting C++11, such as GCC or Clang.
- A terminal or command-line environment.
1. Generate a Test Dataset
Compile and run the synthetic dataset generator:
g++ test_generator.cpp -o test_generator
./test_generator
This generates the reference genome and simulated sequencing reads in the genome_assembly/ directory.
2. Compile the Assembler
g++ genome_assembler.cpp -o genome_assembler
3. Run the Assembler
./genome_assembler genome_assembly/test_inputs.txt genome_assembly/contig_output.txt
The assembler reads the sequencing data and writes the resulting contigs to genome_assembly/contig_output.txt.
🛠️ Implementation Focus
The project explores several core genome-assembly challenges:
- Constructing and traversing de Bruijn graphs.
- Reducing sequencing errors through read preprocessing.
- Identifying and removing graph structures associated with sequencing errors.
- Preserving connectivity through repetitive regions.
- Connecting compatible contigs while respecting structural constraints.
- Generating synthetic datasets for repeatable testing and evaluation.
📁 Project Structure
.
├── main.cpp              # Genome assembler implementation
├── test_generator.cpp    # Synthetic sequencing dataset generator
└── genome_assembly/
    ├── reference_gene.txt
    ├── test_inputs.txt
    └── contig_output.txt
The directory structure above assumes the files are organized as shown.
🎯 Project Goals
The goal of this project is to explore the algorithmic and systems-level challenges of genome assembly, particularly read preprocessing, graph simplification, repeat handling, and contig reconstruction. By implementing these components in C++, the project combines graph algorithms, sequence processing, and performance-conscious data structures in a practical bioinformatics application.