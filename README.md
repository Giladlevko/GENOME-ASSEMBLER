
# De Bruijn Graph Genome Assembler

A C++ genome assembler that reconstructs genomic sequences from short paired-end reads using de Bruijn graphs. The project combines read preprocessing, graph simplification, repeat-aware traversal, and contig stitching to produce contiguous genome assemblies.

## 🚀 Key Features

### - Read Preprocessing: 

- Cleans sequencing reads before graph construction to reduce sequencing errors and improve graph quality.

### - De Bruijn Graph Construction: 

- Builds directed de Bruijn graphs from sequencing reads using k-mers.

### - Graph Simplification:

- Tip Removal: Prunes erroneous dead-end branches introduced by sequencing errors.
- Bubble Detection: Identifies and resolves local graph bubbles associated with sequencing errors or sequence variation.

### - Contig Stitching:

- Validates and connects compatible contig links while enforcing structural constraints to reduce the risk of misassemblies.

### - Synthetic Dataset Generation:

- Includes a C++ test generator for generating sequencing reads for debugging and benchmarking.

## 🧬 Assembly Pipeline

1. Read Cleaning: Preprocesses sequencing reads to reduce errors before graph construction.
2. Graph Construction: Builds a directed de Bruijn graph from the cleaned reads.
3. Graph Simplification: Removes erroneous tips and resolves bubbles to reduce graph complexity.
4. Repeat-Aware Traversal: Preserves multiple contig associations at repeat nodes to avoid unnecessary fragmentation.
5. Contig Stitching: Validates and connects compatible contig links to produce longer assembled sequences.

## 📊 Benchmark Results

Tested on bacterial reference genomes, including Escherichia coli and N. deltocephalinicola, under varying read-length and sequencing-error conditions.

|Metric | Reported Result |
| :--- | :--- |
| E. coli NGA50	| Approximately 178.2 kb |
| Largest alignment | Approximately 486.9 kb (>10% of the genome) |
| Genome fraction | Approximately 98.5% |
| Duplication rate | 1.0 |
| Missasemblies | 0 |
| Local missasemblies | 0 |
| Read lengths | 100 bp |
| Simulated substitution error rate	| 1% |
| Simulated sequencing coverage | 30× |
| Peak memory usage | <700 MB |


The assembler is designed to maintain assembly continuity while limiting structural errors across datasets with different read lengths and error profiles.
Benchmark values depend on the reference genome, dataset, and evaluation methodology.

## 🔍 Technical Highlight: Repeat-Aware Graph Traversal

Repeated sequences can make genome assembly difficult because the same graph node or link may correspond to multiple locations in the genome. Following an ambiguous path can lead to incorrect contig joins.

Rather than forcing connections through repeats, the assembler uses link-selection rules to identify reliable connections and leaves unresolved repeat boundaries unjoined.

### 1. Reciprocal Link Validation

- Requires mutual agreement between connected nodes: Node A must select Node B as its preferred continuation, and Node B must select Node A in return.
- Rejects one-sided connections that do not satisfy this requirement.

### 2. Directional Link Selection

- **Forward pass:** Selects links using a competitive ratio. The best candidate must have a weight greater than 1.3 times that of the second best candidate.
- **Backward pass:** Checks whether the selected edge is also the best edge in the opposite direction.

### 3. Filtering Invalid Edges

- Removes toxic overlaps and excessively wide negative gap connections during matrix construction, preventing these links from being considered during contig stitching.

These rules allow reliable connections to be chained together while leaving ambiguous repeat regions unresolved.

On the evaluated *E. coli* assembly, the assembler achieved **zero detected misassemblies**, according to QUAST.

## 🧪 Synthetic Test Generator

The repository includes test_generator.cpp, a C++ utility for generating synthetic sequencing datasets for debugging, validation, and benchmarking.

### Configuring the Generator

The generator's simulation parameters are defined in `test_generator.cpp`. To customize a dataset, modify the relevant values in the source file and recompile the generator before running it.

The default configuration targets 30× coverage, 50-150 bp reads, a 1% substitution error rate, and paired-end insert distances between 100 and 300 bp.

### Generated Files

| File | Description |
| --- | --- |
| genome_assembly/reference_gene.txt | Ground-truth reference sequence used for evaluation |
| genome_assembly/test_inputs.txt | Simulated sequencing reads supplied to the assembler |


The reference sequence provides ground truth for evaluating the resulting assembly.

## ⚙️ Building and Running

### Prerequisites

- A C++ compiler, such as GCC or Clang, supporting C++11.
- A terminal or command-line environment.

### 1. Generate a Test Dataset

Compile and run the synthetic dataset generator:

```bash
g++ test_generator.cpp -o test_generator
./test_generator
```

This generates the reference genome and simulated sequencing reads in the genome_assembly/ directory.

### 2. Compile the Assembler

```bash
g++ main.cpp -o main
```

### 3. Run the Assembler

```bash
./main genome_assembly/test_inputs.txt genome_assembly/contig_output.txt
```

The assembler reads the input data and writes the resulting contigs to genome_assembly/contig_output.txt.

## 🛠️ Implementation Focus

The project explores several core genome-assembly challenges:
- Constructing and traversing de Bruijn graphs.
- Reducing sequencing errors through read preprocessing.
- Identifying and removing graph structures associated with sequencing errors.
- Preserving connectivity through repetitive regions.
- Connecting compatible contigs while respecting structural constraints.
- Generating synthetic datasets for repeatable testing and evaluation.

## 📁 Project Structure
<pre>
.
├── main.cpp
├── test_generator.cpp
├── source/
    ├── common.hpp
    ├── de_bruijn_graph.hpp
    ├── genome_assembler.hpp
    └── kmer_bit_map.hpp
└── genome_assembly/
    ├── reference_gene.txt
    ├── test_inputs.txt
    └── contig_output.txt

</pre>

## 🎯 Project Goals

This project explores algorithmic and systems-level challenges in genome assembly, particularly read preprocessing, graph simplification, repeat handling, and contig reconstruction. Implemented in C++, it combines graph algorithms, sequence processing, and performance-conscious data structures in a practical bioinformatics application.