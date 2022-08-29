#!/usr/bin/env python3

"""
Encode or decode a genotype matrix extracted from a VCF file.

--------------
The VCF format
--------------

A VCF file consists of two main parts: the header and the actual variant records. Each variant record is stored in a
single line. Each variant record is separated into annotations, which can either have site level or sample level.

======================
Site-level annotations
======================

The first 8 fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) represent the properties observed at the variant site.
Hence, they are referred to as site-level annotations. When multiple samples are represented in a VCF file, some of the
site-level annotations represent a summary or average of the values obtained for that site from different samples.

The following site-level annotations are explicitly relevant in the scope of this program (explanations partly taken
from https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it/):

- CHROM: the contig (usually a chromosome) on which the variant occurs
- POS: the 1-based genomic coordinate on the contig of the variant. Note that for deletions the position given is the
  base preceding the event.
- REF: the reference allele observed
- ALT: the alternative allele(s) observed in a sample, or a set of samples. Note that REF and ALT are given on the
  forward strand. For insertions, the ALT allele includes the inserted sequence as well as the base preceding the
  insertion. For deletions, the ALT allele is the base before the deletion.

========================
Sample-level annotations
========================

Sample-level annotations are contained in the 9th column (FORMAT) and in the sample columns (10th and beyond).
Sample-level annotations are tag-value pairs. The tags are recorded in the FORMAT field. The values are then recorded in
corresponding order in each sample column.

The sample-level annotation with the tag GT is of importance. Its value indicates the genotype of a sample at a variant
site. For a diploid organism, its value indicates the two alleles carried by the sample, encoded by a 0 for the REF
allele, 1 for the first ALT allele, 2 for the second ALT allele, etc. When there is for example a single ALT allele
(which is by far the most common case), the value will be either:

- 0/0 - the sample is homozygous reference
- 0/1 - the sample is heterozygous, carrying one copy of each the REF and ALT alleles
- 1/1 - the sample is homozygous alternate

For non-diploids, the same pattern applies; in the haploid case there will be just a single number; for polyploids there
will be more, e.g. 4 numbers for a tetraploid organism.

Below is an example variant record::

    CHROM	POS	REF	ALT	FORMAT	S1
    1		1	T	G	GT	0/1

Here, a sample named S1 is observed with the allele combination T/G at position 1 on chromosome 1.

-------------------
The genotype matrix
-------------------

Suppose a set of m variant records. Each variant record contains GT annotations for n samples named S1, S2, ..., SN. The
genotype matrix G can be obtained by extracting only the GT values from the variant records.

Below are m=3 example variant records, each containing GT annotations for n=2 samples named S1 and S2::

    CHROM	POS	REF	ALT	FORMAT	S1	S2
    1		1	T	G	GT	0/1	0/0
    1		1	T	G	GT	0/1	0/0
    1		1	T	G	GT	0/1	0/0

The corresponding m-by-n genotype matrix G is as follows::

        0/1 0/0
    G = 0/1 0/0
        0/1 0/0

-------
Phasing
-------

For each diploid individual homologous chromosomes are labelled paternal and maternal. If the genotypes of the parents
are not known the labelling is arbitrary. Consider two SNPs: if it is possible to assign, for heterozygous calls, which
allele is on the paternal chromosome and which one is on the maternal chromosome, then it is possible to phase the
alleles at the two SNP sites. This suite of “ordered” SNPs then constitutes a haplotype.

Phasing refers to the separation of a consensus sequence into individual sequence strands to identify which variants are
grouped together or in phase. In the VCF format, for multiploid data the separator indicates whether the data are
phased (|) or unphased (/).

In the following example the two alleles C and G occur on the same chromosome in the sample named S1::

    CHROM	POS	REF	ALT	FORMAT	S1
    1		2	C	T	GT	0|1
    1		3	A	G	GT	1|0

-----------------------------
Coding of the genotype matrix
-----------------------------

The coding of a genotype matrix G comprises the following steps:
1. Row-wise splitting of the genotype matrix G into blocks that are subsequently processed separately. In this process,
   blocks may be constructed in such a way that they contain only a certain class of genomic variation (e.g. SNPs,
   indels or structural variants). An index must be maintained to reconstruct the original row order.
2. Splitting of the genotype matrix G into an allele matrix A and a phasing matrix P
3. Optional binarization of the allele matrix A (this process yields either bit planes B_q or a binary allele matrix C)
4. Optional row- and column-wise sorting of the allele matrix A, or the bit planes B_q, or the binary allele matrix C,
   and of the phasing matrix P
5. Entropy coding of the allele matrix A, or the bit planes B_q, or the binary allele matrix C, and of the phasing
   matrix P
"""

import sys
import logging as log

from .encoder import Encoder
from .decoder import Decoder
from .settings import PROGRAM_NAME, PROGRAM_DESC, LOG_LEVELS


def run(args, run_as_module: bool):
    #? Basic log config
    format_string = '[' + PROGRAM_NAME + '] [%(asctime)s] [%(levelname)-8s] --- [%(processName)-11s] [%(filename)15s] [%(funcName)20s] %(message)s'
    log.basicConfig(format=format_string, level=log.INFO)

    #? Log level
    try:
        log_level = LOG_LEVELS[args.log_level]
    except KeyError:
        log.warning("invalid log level '%s' (using fall-back log level 'info')", args.log_level)
        log_level = log.INFO
    logger = log.getLogger()
    logger.setLevel(log_level)

    #? Print banner
    log.info('********************************************************************************')
    log.info('    {}'.format(PROGRAM_NAME))
    log.info('    {}'.format(PROGRAM_DESC))
    log.info('********************************************************************************')

    #? Log the command line
    if not run_as_module:
        log.debug('command line: %s', ' '.join(sys.argv))

    #? Log the arguments
    log.debug('arguments: ')
    for arg in vars(args):
        log.debug('  %-16s: %s', arg, getattr(args, arg))

    if args.mode == "decode":
        if args.pos is None and args.samples is None:
            decoder = Decoder(args.input, args.output)
            decoder.decode_all()
            
        else:
            decoder = Decoder(args.input, args.output)
            decoder.random_access(
                args.pos,
                args.samples
            )

    elif args.mode == "encode":
        encoder = Encoder(
            args.input, 
            args.output, 
            binarization_name=args.binarization,
            axis=args.axis,
            sort_rows=args.sort_rows,
            sort_cols=args.sort_cols,
            dist=args.dist,
            solver=args.solver,
            preset_mode=args.preset_mode,
            codec_name=args.encoder,
            block_size=args.block_size,
            num_threads=args.num_threads
        )
        
        encoder.run()
    
    elif args.mode is None:
        raise RuntimeError("Please enter the mode.")
        
    else:
        raise NotImplementedError("Mode is not available:{}".format(args.mode))