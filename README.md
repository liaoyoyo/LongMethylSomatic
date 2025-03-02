# LongMethylSomatic v1.0.0


## !!! Attention, still needs time to grow !!!

---

## **Program Objective**
This program is designed to detect the methylation status around somatic mutation sites and analyze whether the methylation ratio at adjacent CpG sites is consistent between mutant and non-mutant reads. Through this analysis, we can determine the association between somatic mutations and methylation patterns, providing insights into the relationship between genomic variations and epigenetics.

---

## **Input**
The program requires the following input files:

1. **Normal BAM File** (`-n`)
   - BAM file from normal tissue containing read alignment information and methylation tags.
2. **Tumor BAM File** (`-t`)
   - BAM file from tumor tissue, also containing read alignment information and methylation tags.
3. **Somatic VCF File** (`-v`)
   - A somatic mutation VCF file generated by tools like ClairS, containing mutation site details such as chromosome, position, reference allele, and mutant allele.
4. **Reference Genome File** (`-r`)
   - A reference genome file in FASTA format (e.g., hg38) used for alignment and positional reference.
5. **Output File Path** (`-o`)
   - The path for saving the analysis results.

---

## **Output**
The program generates a text file containing the following analysis results:

- **Somatic Mutation Sites**  
  - Example format: `chr1:1000 (A>G)`, indicating the chromosome, position, and allele change of the mutation.
- **Adjacent CpG Sites**  
  - Example format: `chr1:950`, representing a CpG site near the mutation site used for methylation analysis.
- **Methylation Ratio of Mutant Reads**  
  - Example format: `90% methylated (18/20)`, indicating the methylation percentage and count of methylated reads vs. total reads in mutant reads.
- **Methylation Ratio of Non-Mutant Reads**  
  - Example format: `10% methylated (2/20)`, showing the methylation ratio in non-mutant reads.
- **Conclusion**  
  - Example format: `Significant methylation discrepancy` or `No significant methylation difference`, determined based on the methylation ratio comparison.

---

## **Workflow**
The program follows these main steps:

1. **Parsing Command-Line Arguments**  
   - Reads input file paths from command-line arguments (Normal BAM, Tumor BAM, VCF, Reference Genome, Output File).
   - Displays usage instructions and terminates execution if required parameters are missing or incorrect.

2. **Loading the Somatic VCF File**  
   - Reads the VCF file and extracts information for each somatic mutation site (e.g., chromosome, position, reference, and mutant alleles).
   - Stores mutation sites in an internal data structure for processing.

3. **Opening BAM Files**  
   - Loads Normal BAM and Tumor BAM files, including headers and indexes, to enable fast lookup of reads in specific regions.

4. **Analyzing Mutation Sites**  
   - For each somatic mutation site:
     - **Extract Reads**: Retrieves reads from the Tumor BAM file within a predefined range (e.g., ±50 bp around the mutation site).
     - **Classify Reads**: Categorizes reads into mutant reads (matching the mutant allele in the VCF) and non-mutant reads.
     - **Calculate Methylation Ratios**: Computes the methylation status at adjacent CpG sites separately for mutant and non-mutant reads.
     - **Assess Consistency**: Compares methylation ratios between mutant and non-mutant reads to determine methylation pattern differences.

5. **Generating Output**  
   - Organizes the analysis results in a structured text format.
   - Writes the results to the specified output file.

---

## **Summary**
This program integrates information from somatic VCF and BAM files to analyze the methylation status of CpG sites near mutation sites. It evaluates whether methylation patterns differ between mutant and non-mutant reads, providing a systematic and efficient approach for large-scale genomic data analysis. The results contribute to research on the association between mutations and epigenetics.

