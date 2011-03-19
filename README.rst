
Illumina reads adapter screening utilities
=================================================
This contains several utility programs that removes the adapters and low quality bases from
Illumina reads. 

:Author: Haibao Tang (`tanghaibao <http://github.com/tanghaibao>`_)
:Email: htang@jcvi.org
:License: `BSD <http://creativecommons.org/licenses/BSD/>`_

.. contents ::

Installation
-------------
The program depends on the excellent `SeqAn library <http://www.seqan.de/>`_.
Please download the library, and place ``seqan/`` in the same folder, and run::

    make


trimReads
----------
Functionality emulates `cutadapt <http://code.google.com/p/cutadapt/>`_.
The adapter sequences are identified through `Waterman-Eggert` algorithm
implemented in `SeqAn <http://www.seqan.de/>`_. The quality trimming are a
simple algorithm that takes the quality values, deduct a user specified cutoff,
and then finds the `max-sum segment
<http://en.wikipedia.org/wiki/Maximum_subarray_problem>`_. This method
guarantees that the average base quality is higher than the user cutoff. 

There are other options to cut adapters, including `cutadapt
<http://code.google.com/p/cutadapt/>`_ and `FASTX_TOOLKIT
<http://hannonlab.cshl.edu/fastx_toolkit/>`_. The main advantage of this program:

* Accepts an adapter FASTA file
* Fast, robust and flexibility
* Qual/adapter trimming in one step
* Can trim both 5`- and 3`- end
* Can trim adapters multiple times

Just run::

    trimReads

to see a list of program options::

    Illumina reads trimming utility
    Author: Haibao Tang <htang@jcvi.org>

    Usage: trimReads [options] fastqfile

      -h, --help              displays this help message
      -o, --outfile           Output file name, uses Sanger encoding for quality. (default replace suffix with .trimmed.fastq)
      -f, --adapterfile       FASTA formatted file containing the adapters for removal  (default adapters.fasta)
      -s, --score             Minimum score to call adapter match. Default scoring scheme for +1 match, -3 for mismatch/gapOpen/gapExtension. (default 15)
      -n, --times             Try to remove the adapters at most COUNT times. Useful when an adapter gets appended multiple times. (default 4)
      -q, --quality-cutoff    Trim low-quality regions below quality cutoff. The algorithm is similar to the one used by BWA by finding a max-sum segment within the quality string. Set it to 0 to skip quality trimming.  (default 20)
      -m, --minimum-length    Discard trimmed reads that are shorter than LENGTH. (default 64)
      -Q, --quality-encoding  Read quality encoding for input file. 64 for Illumina, 33 for Sanger.  (default 64)

Find a list of adapters to remove (more will slow down search), default is ``adapters.fasta``. When ready::

    trimReads test.fastq

to get a trimmed file `test.trimmed.fastq`. To turn off the quality trimming, just set ``-q`` to ``0``::

    trimReads -q 0 test.fastq


sortPairedReads
----------------
This program sorts all read pairs into three sets:

* Adapter set: the pairs with either /1 or /2 match adapters (in most cases
both will match). These are fragments up to 1X read length.
* Overlap set: the pairs with /1 and /2 having dovetail overlap. These are
fragments up to 2X read length.
* Clean set: survived the above two searches.

The reason for this sorting is to get rid of the short fragments (set 1 and set
2) commonly in the Illumina PE library. Some libraries are worse than others.
The goal is to input the mated library within nominal insert size ranges.

Just run::

    sortPairedReads

to see a list of program options::

    Sort pairs of Illumina reads
    Author: Haibao Tang <htang@jcvi.org>

    Usage: sortPairedReads [options] fastqfile1 fastqfile2

      -h, --help               displays this help message
      -f, --adapterfile        FASTA formatted file containing the adapters for removal  (default adapters.fasta)
      -s, --adapterMatchScore  Minimum score to call adapter match. Default scoring scheme for +1 match, -3 for mismatch/gapOpen/gapExtension. (default 15)
      -t, --endMatchScore      Minimum score to call dovetail match. Default scoring scheme for +1 match, -3 for mismatch/gapOpen/gapExtension. (default 20)
      -Q, --quality-encoding   Read quality encoding for input file. 64 for Illumina, 33 for Sanger.  (default 64)
      -v, --verbose            Print alignments for debugging  (default 0)
     
For any given two fastq files, the output contains 4 files: ``fastqfile1.adapters.fastq`` (set 1),
``fastqfile1.overlap.fastq`` (set 2), ``fastqfile1.clean.fastq`` and
``fastqfile2.clean.fastq`` (set 3). For genome assembler inputs, I recommend
discard set 1, treat set 2 as **unmated**, and treat set 3 as **mated**.

