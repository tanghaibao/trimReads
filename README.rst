
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
<http://en.wikipedia.org/wiki/Maximum_subarray_problem>`_.


Installation
-------------
The program depends on the excellent `SeqAn library <http://www.seqan.de/>`_.
Please download the library, and place ``seqan/`` in the same folder, and run::

    make


Usage
------
Just run::

    trimReads

to see a list of program options. And use::

    trimReads test.fastq

to get a trimmed file `test.trimmed.fastq`.

