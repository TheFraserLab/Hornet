######
Hornet
######

A pipeline for unbiased allele-specific read mapping based on the
`WASP pipeline <https://github.com/bmvdgeijn/WASP>`_ by the
`Pritchard lab <http://pritchardlab.stanford.edu/home.html>`_.

Introduction
############

Hornet is a simple set of tools for unbiased allele-specific read mapping
that is originally described in
`van de Geijn B, McVicker G, Gilad Y, Pritchard JK. "WASP: allele-specific software for robust discovery of molecular quantitative trait loci" <https://www.nature.com/articles/nmeth.3582>

Overview
########

Reads are mapped normally using a mapper chosen by the user (must output
BAM or SAM format).  Then mapped reads that overlap single nucleotide
polymorphisms (SNPs) are identified. For each read that overlaps a SNP, its
genotype is swapped with that of the other allele and the read is re-mapped.
Re-mapped reads that fail to map to exactly the same location in the genome are
discarded.

Step 1:
-------

Map the fastq files using your favorite mapper/options and filter for quality using a cutoff of your choice

Example:
~~~~~~~~

.. code:: shell
 
	STAR \
            --genomeDir genome_dir/STAR \
            --outFileNamePrefix analysis_dir/STAR1/ \
            --outSAMattributes MD NH \
            --outSAMtype BAM Unsorted \
            --outTmpDir analysis_dir/STAR1/STARtmp \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --readFilesIn input.R1.fastq.gz input.R2.fastq.gz \
	    --outFilterMultimapNmax 1

Step 2:
-------

Use `find_intersecting_snps` to identify reads that may have mapping biases. 
It is recommended to remove duplicate reads prior to this step to save time 
provessing reads that will ultimately be discarded.

Usage::

	find_intersecting_snps.py [-p][-C][-P][-s] <input.bam> <SNP_directory>
	   -p indicates that reads are paired-end (default is single)
	   -C limits the list of SNPs to those found on a specific chromosome
	   -P indicates that the list of SNPs are phased; if so, only 2 reads will generated 
	      per input read that overlaps a SNP (as in a hybrid)
	   -s indicates that input data are coordinate sorted
	   <input.bam> is the bamfile from the initial mapping process
	   <SNP_directory> is the directory containing the SNPs segregating within the
	      sample in question (which need to be checked for mappability issues).  This directory
	      should contain coordinate sorted files of SNPs separated by chromosome and named:
	         chr<#>.snps.txt.gz
	      These files should contain 3 columns: position RefAllele AltAllele


Output::

	input.keep.bam - bamfile with reads that did not intersect SNPs or indels and therefore can
	   be kept without remapping
	input.to.remap.bam - bamfile with original reads that overlapped SNPs that need to be remapped
	input.remap.fq.gz - fastq file containing the reads with the new variants to remap. If the
	    paired-end option is used two files ending with .fq1.gz and .fq2.gz will be output.
	    
	Run statistics: run statistics will be reported once this step is complete, and include:
	  - Total input reads
	  - Reads with no SNPs
	  - Reads overlapping SNPs
	  - Total SNPs covered (refers to the number of instances where a read was found to overlap a SNP)
	  - Reference SNP matches (refers to the % of 'Total SNPs covered' that match EITHER allele in the SNP file; 
	  	should be as close to 100% as possible if SNPs have been called correctly)
	  - Non-reference SNP matches (should be a low %, caused by mutations or sequencing errors)
	  - Reads dropped [INDEL]
	  - Reads dropped [too many SNPs]
	  - Reads dropped [multivalent SNPs] (indicates a single read had SNPs matching both the ref and alt alleles)
	
	Note: Reads that overlap indels are currently excluded and will not be present in any of the 'remap' files
	or the input.keep.bam file. For this reason the total number of reads will not add up to the
	number of reads provided in the input.sort.bam file.


Example (after removing duplicates):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: shell

   python find_intersecting_snps.py \
        -p -P -s\
        rmdup.bam SNP_directory

Step 3
------
Map the input.remap.fq.gz file(s) using the same mapping arguments used in Step 1. 
The arguments should be exactly the same as those in Step 1 EXCEPT for arguments that
directly modify the reads that are used by the aligner. For example the read trimming
arguments to bowtie (-3 and -5 arguments) should be used in Step 1 ONLY because
they modify the reads that are output by bowtie. If you used 2-pass mapping in STAR 
where splice junctions were included, you can use the same SJ.out.tab file for remapping 
in this step.

Example:
~~~~~~~~

.. code:: shell

	STAR \
            --genomeDir genome_dir/STAR \
            --outFileNamePrefix analysis_dir/remap/ \
            --outSAMattributes MD NH \
            --outSAMtype BAM Unsorted \
            --outTmpDir analysis_dir/remap/STARtmp \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --readFilesIn rmdup.remap.fq1.gz rmdup.remap.fq2.gz \
	    --outFilterMultimapNmax 1
	    
	    mv analysis_dir/remap/remapAligned.out.bam rmdup.remap.bam

Step 4
------
Use filter_remapped_reads.py to retrieve reads that remapped correctly.
The remapped bam file MUST be sorted by read name for this step. The read names
contain the original mapping information, which is needed in this step to 
determine whether the reads mapped to the same location. 

Usage::

	filter_remapped_reads.py [-p] <to.remap.bam> <remapped_reads.bam> <output.bam>
	   -p option indicates that the reads are paired-end
	   <to.remap.bam> output from find_intersecting_snps.py which contains
	      the original aligned reads that were remapped
	   <remapped_reads.bam> output from the second mapping step (Step 3)
	   <output.bam> file where reads that are kept after remapping are stored
	   

Example:
~~~~~~~~

.. code:: shell

 samtools sort -n rmdup.remap.bam -o rmdup.remap.sort.bam

 python filter_remapped_reads.py \
        -p \
        rmdup.to.remap.bam rmdup.remap.sort.bam \
        rmdup.remap.kept.bam

At the end of the pipeline, rmdup.remap.keep.bam and rmdup.remap.keep.bam
can be sorted and merged for a complete set of mappability filtered aligned reads. 
The merged file should then be indexed:

.. code:: shell

  samtools sort rmdup.remap.kept.bam -o rmdup.remap.kept.sort.bam
  samtools sort rmdup.keep.bam -o rmdup.keep.sort.bam
  samtools merge rmdup.remap.kept.merged.bam rmdup.keep.sort.bam rmdup.remap.kept.sort.bam
  samtools index rmdup.remap.kept.merged.bam



Dependencies
############

Hornet is writte in python and will work with python 2.6+. It requires
`numpy <http://www.numpy.org>`_, `scipy <http://www.scipy.org>`_, and
`pysam <https://github.com/pysam-developers/pysam>`_.

It also depends on `argparse <https://code.google.com/p/argparse/>`_,
which is included by default in newer versions of python (>= 2.7).

Installation
############

.. code:: shell
   pip install https://github.com/TheFraserLab/Hornet/tarball/master

Dependencied will be installed automatically.
