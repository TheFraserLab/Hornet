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
`van de Geijn B, McVicker G, Gilad Y, Pritchard JK. "WASP: allele-specific software for robust discovery of molecular quantitative trait loci" <http://biorxiv.org/content/early/2014/11/07/011221>`_

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

  tophat --no-coverage-search -o ${LANE_NAME}_out Sequences/hg18_norand ${LANE_NAME}.fastq.gz
  samtools view -b -q 10 ${LANE_NAME}_out/accepted_hits.bam > ${LANE_NAME}_out/accepted_hits.quality.bam

Step 2:
-------

Use `find_intersecting_snps` to identify reads that may have mapping biases

Usage::

	find_intersecting_snps.py [-p] <input.bam> <SNP_file_directory>
	   -p indicates that reads are paired-end (default is single)
	   -m changes the maximum window to search for SNPs.  The default is
	      100,000 base pairs.  Reads or read pairs that span more than this distance
	      (usually due to splice junctions) will be thrown out.  Increasing this window
	      allows for longer junctions, but may increase run time and memory requirements.
	   <input.bam> is the bamfile from the initial mapping process
	   <SNP_file_directory> is the directory containing the SNPs segregating within the
	      sample in question (which need to be checked for mappability issues).  This directory
	      should contain sorted files of SNPs separated by chromosome and named:
	         chr<#>.snps.txt.gz
	      These files should contain 3 columns: position RefAllele AltAllele


Output::

	input.sort.bam - Sorted bamfile of the original input
	input.keep.bam - bamfile with reads that did not intersect SNPs or indels and therefore can
	   be kept without remapping
	input.to.remap.bam - bamfile with original reads that overlapped SNPs that need to be remapped
	input.to.remap.num.gz - the number of variants of the original read that must be remapped
	input.remap.fq.gz - fastq file containing the reads with the new variants to remap. If the
	    paired-end option is used two files ending with .fq1.gz and .fq2.gz will be output.
	
	Note: Reads that overlap indels are currently excluded and will not be present in any of the 'remap' files
	or the input.keep.bam file. For this reason the total number of reads will not add up to the
	number of reads provided in the input.sort.bam file.

To make the snps.txt.gz files, from vcf files separated by chromosome:

.. code:: bash

    mkdir snps
    for i in chr*.vcf.gz; do j=$(echo $i | sed 's/\..*//'); pigz -dc $i | grep -v "^#" | awk '{if ((length($4) == 1) && (length($5) == 1)) printf ("%s\t%s\t%s\n", $2, $4, $5)}' | pigz > ${j}.snps.txt.gz; done


From a single vcf:

.. code:: bash

    mkdir snps
    pigz -dc data.vcf.gz | grep -v "^#" | awk '{if ((length($4) == 1) && (length($5) == 1)) printf ("%s\t%s\t%s\n", $2, $4, $5) | "pigz > snps/"$1".snps.txt.gz"}'

pigz just parallelized gzip, if you don't have it, substitute `pigz` with `gzip`.

Example:
~~~~~~~~

.. code:: shell

   find_intersecting_snps ${LANE_NAME}_out/accepted_hits.quality.bam SNP_files/

Step 3
-----
Map the input.remap.fq.gz using the same mapping arguments used in Step 1. Note that
the arguments should be exactly the same as those in Step 1 EXCEPT for arguments that
directly modify the reads that are used by the aligner. For example the read trimming
arguments to bowtie (-3 and -5 arguments) should be used in Step 1 ONLY because
they modify the reads that are output by bowtie.

Example:
~~~~~~~~

.. code:: shell

  tophat --no-coverage-search -o ${LANE_NAME}_out_remap hg18_norand ${LANE_NAME}_out/accepted_hits.quality.remap.fq.gz
  samtools view -b -q 10 ${LANE_NAME}_out_remap/accepted_hits.bam > ${LANE_NAME}_out_remap/accepted_hits.quality.bam


Step 4
------
Use filter_remapped_reads.py to retrieve reads that remapped correctly

Usage::

	filter_remapped_reads.py [-p] <to.remap.bam> <remapped_reads.bam> <output.bam> <to.remap.num.gz>
	   -p option indicates that the reads are paired-end
	   <to.remap.bam> output from find_intersecting_snps.py which contains
	      the original aligned reads that were remapped
	   <remapped_reads.bam> output from the second mapping step (Step 3)
	   <output.bam> file where reads that are kept after remapping are stored
	   <to.remap.num.gz> is the file from find_intersecting_snps.py which contains
	      the number of remapped sequences

Example:
~~~~~~~~

.. code:: shell

  filter_remapped_reads ${LANE_NAME}_out/accepted_hits.quality.to.remap.bam ${LANE_NAME}_out_remap/accepted_hits.quality.bam ${LANE_NAME}.remap.keep.bam ${LANE_NAME}_out/accepted_hits.quality.to.remap.num.gz

At the end of the pipeline, ${LANE_NAME}.keep.bam and ${LANE_NAME}.remap.keep.bam
can be merged for a complete set of mappability filtered aligned reads. The merged
file should then be sorted and indexed:

.. code:: shell

  samtools merge ${LANE_NAME}.keep.merged.bam ${LANE_NAME}.keep.bam ${LANE_NAME}.remap.keep.bam
  samtools sort ${LANE_NAME}.keep.merged.bam ${LANE_NAME}.keep.merged.sorted
  samtools index ${LANE_NAME}.keep.merged.sorted.bam

Step 5 (Optional)
-----------------

Filter duplicate reads. Programs such as samtools rmdup introduce bias when
they filter duplicate reads because they retain the read with the highest score
(which usually matches the reference).

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

Testing
#######

To run the tests, execute `py.test` from within this directory.
