
Adapter and quality trimming for Illumina reads
=================================================
This is a utility program that removes the adapters and low quality bases from
Illumina reads. Functionality emulates `cutadapt <http://code.google.com/p/cutadapt/>`_.

:Author: Haibao Tang (`tanghaibao <http://github.com/tanghaibao>`_)
:Email: htang@jcvi.org
:License: `BSD <http://creativecommons.org/licenses/BSD/>`_

.. contents ::

Algorithm
----------
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


Installation
-------------
The program depends on the excellent `SeqAn library <http://www.seqan.de/>`_.
Please download the library, and place ``seqan/`` in the same folder, and run::

    make


Usage
------
Just run::

    trimReads

to see a list of program options::

    Illumina reads trimming utility
    Author: Haibao Tang <htang@jcvi.org>

    Usage: trimReads [options] fastqfile

      -h, --help              displays this help message
      -o, --outfile           Output file name, uses Sanger encoding for quality. (default replace suffix with .trimmed.fastq)
      -f, --adapterfile       FASTA formatted file containing the adapters for removal [default: `adapters.fasta`] (default adapters.fasta)
      -s, --score             Minimum score to call adapter match. Default scoring scheme for +1 match, -3 for mismatch/gapOpen/gapExtension. (default 15)
      -n, --times             Try to remove the adapters at most COUNT times. Useful when an adapter gets appended multiple times. (default 4)
      -q, --quality-cutoff    Trim low-quality regions below quality cutoff. The algorithm is similar to the one used by BWA by finding a max-sum segment within the quality string. (default 20)
      -m, --minimum-length    Discard trimmed reads that are shorter than LENGTH. (default 64)
      -Q, --quality-encoding  Read quality encoding for input file. 64 for Illumina, 33 for Sanger. Output will always be Sanger encoding. (default 64)

Find a list of adapters to remove (more will slow down search), default is ``adapters.fasta``. When ready::

    trimReads test.fastq

to get a trimmed file `test.trimmed.fastq`.

