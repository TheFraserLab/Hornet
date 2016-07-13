import gzip
from sys import argv

if __name__ == "__main__":
    fastq = gzip.open(argv[1], 'rt')
    base, _ = argv[1].rsplit('.remap', 1)
    num = gzip.open(base+'.to.remap.num.gz', 'wt')
    line = next(fastq)
    last = ''
    while line:
        curr, counts = line.strip().rsplit(':', 1)
        if curr != last:
            counts = int(counts)
            num.write('{}\n'.format(2*counts))
            num.flush()
        last = curr
        seq = next(fastq)
        header_again = next(fastq)
        quals = next(fastq)
        line = next(fastq)




