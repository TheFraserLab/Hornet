"""filter_remapped_reads

This is a rewrite of the official WASP code, which no longer requires
a .num.gz file. This code is compatable with the rewritten version of
the first step (find_intersecting_snps), which no longer generates the
.num.gz file. The bam file of remapped reads must be sorted by read name.

"""
from __future__ import print_function
import sys
import gzip
import time

import pysam


def run(to_remap_bam, remap_bam, keep_bam, is_paired_end):
    """Core function."""

    sys.stderr.write(
        "{} Starting run\n".format(
            time.strftime(("%b %d ") + time.strftime("%I:%M:%S"))
        )
    )

    to_remap_bam = pysam.Samfile(to_remap_bam, "rb")
    remap_bam = pysam.Samfile(remap_bam, "rb")
    keep_bam = pysam.Samfile(keep_bam, "wb", template=to_remap_bam)

    # List of correctly mapped reads, represented by read ID
    # (e.g. 42 in @42:chr1:14536:3)
    correct_maps = []

    # List of how many alternate versions of each read there are
    # (e.g. 3 in @42:chr1:14536:3)
    nums = []

    end_of_file = False
    prev_read_num = 0

    # Keep track of remapped read index
    counter = 1
    skip = 0

    # Get a list of reads that remapped correctly
    remap_read = next(remap_bam)

    while not end_of_file:

        chrm = remap_read.qname.strip().split(":")[1]

        if remap_read.is_reverse:
            if is_paired_end:
                pos = int(remap_read.qname.strip().split(":")[3])
            else:
                pos = int(remap_read.qname.strip().split(":")[2])

        else:
            pos = int(remap_read.qname.strip().split(":")[2])

        read_num = int(remap_read.qname.strip().split(":")[0])
        num = int(remap_read.qname.strip().split(":")[-1])

        # For each original read, keep track of how many alternate reads were mapped
        if read_num != prev_read_num:

            # If none of the alternate versions of the read got mapped, set num to 1
            if read_num != counter:
                skipped = read_num - counter
                nums += ([1]*skipped)
                counter = read_num

            if is_paired_end:
                nums.append(num*2)
            else:
                nums.append(num)

            prev_read_num = read_num
            counter += 1

        if (remap_read.tid != -1 and remap_read.pos == pos and
            remap_bam.getrname(remap_read.tid) == chrm):

            # Throw out the remapped read if it remapped with a deletion
            dels = 0
            for cig in remap_read.cigar:
                if not cig[0] in (0, 3, 4):
                    dels += 1

            if dels == 0:
                correct_maps.append(read_num)

        try:
            remap_read = next(remap_bam)
        except:
            end_of_file = True

    correct_maps.sort()

    # Original aligned reads
    orig_read = next(to_remap_bam)

    # Number of different reads generated from the original read
    orig_num = int(nums[0])

    # Line number of the remap_bam file (if single end data) or read pair
    # number if paired end data.
    line_num = 1

    # Index for walking through correct_maps.
    map_indx = 0

    # Number of correctly mapped reads for the current read (pair).
    correct = 0

    # Total number of correctly mapped read (pairs).
    total_correct = 0

    end_of_file = False
    while (not end_of_file and
           (map_indx < len(correct_maps)) and
           (line_num <= correct_maps[-1])):
        if line_num != correct_maps[map_indx]:

            # If we saw the correct number of remaps for the last read, keep it.
            if correct == orig_num:
                total_correct += 1
                keep_bam.write(orig_read)

                # If the data is paired end, write out the paired read.
                if is_paired_end:
                    try:
                        orig_read = next(to_remap_bam)
                    except:
                        raise ValueError(
                            "File ended unexpectedly (no pair found)."
                        )
                    keep_bam.write(orig_read)
            else:
                try:
                    second_read = next(to_remap_bam)
                except:
                    end_of_file=True
                    break
            try:
                orig_read = next(to_remap_bam)
                orig_num = nums[line_num]
            except StopIteration:
                end_of_file = True
            line_num += 1
            correct = 0
        else:
            correct += 1
            map_indx += 1

    if correct == orig_num:
        total_correct += 1
        keep_bam.write(orig_read)

        # If the data is paired end, write out the paired read.
        if is_paired_end:
            try:
                orig_read = next(to_remap_bam)
            except:
                sys.stderr.write("File ended unexpectedly (no pair found).")
                exit()
            keep_bam.write(orig_read)
    sys.stderr.write(
        "{} Finished\n\n".format(
            time.strftime(("%b %d ") + time.strftime("%I:%M:%S"))
        )
    )
    sys.stderr.write(
        ("RUN STATISTICS:\n\tTotal remapped read (pair)s: {}\n\t"
         "Read (pair)s remapped to the correct position: {} {:.2%}\n")
        .format(
            line_num, total_correct, total_correct/float(line_num)
        )
    )

def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", action='store_true', dest='is_paired_end',
                        default=False, help=('Indicates that reads are '
                                             'paired-end (default is single).'))
    h = ('to.remap.bam file from find_intersecting_snps.py.')
    parser.add_argument("to_remap_bam", help=h)
    parser.add_argument("remap_bam", help='Remapped bam file.')
    parser.add_argument("keep_bam", help=(
        'File to write correctly remapped reads to.'
    ))

    options = parser.parse_args()

    run(options.to_remap_bam, options.remap_bam, options.keep_bam,
        options.is_paired_end)

if __name__ == '__main__':
    main()
