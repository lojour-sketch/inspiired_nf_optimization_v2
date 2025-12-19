# GENOME_INDEXING

## Description

This process indexes the reference genome using bowtie2, to be able to do the alignment a posteriori.

---

### Parameter configuration

- **runMode**: This is the run mode of STAR. This time we use the `--runMode genomeGenerate` option, which creates the indexed genome.  
- **genomeDir**: Directory that will contain the indexed genome.  
- **genomeFastaFiles**: Path to the reference genome fasta file.  
- **runThreadN**: Number of threads to use for indexing.  

---

## Tools

| Tool | Description                                   | Homepage |
|------|-----------------------------------------------|----------|
| STAR | Alignment tool for mapping DNA sequences     | [GitHub](https://github.com/alexdobin/STAR) |

---

## Input

| Name           | Type   | Description                       |
|----------------|--------|-----------------------------------|
| `genome_name`    | string | Name of the genome                |
| `refGemomeFile`  | file   | Path to the reference genome fasta file |

---

## Output

| Name              | Type      | Description                                      |
|------------------|-----------|--------------------------------------------------|
| `genome_name`       | string    | Name of the genome                               |
| `indexedGenomeDir`  | directory | Directory that will contain the indexed genome  |

---

## Authors

- @liberentaizp
