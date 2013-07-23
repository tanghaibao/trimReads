/*
 * Trim artificial base pairs within fastq reads.
 *
 * Author: Haibao Tang <htang@jcvi.org>
 * Date: 01/30/2011
 * License: BSD
 *
 * Trim given adapters using local alignments. Trimmed regions will be given
 * low phred quality (1) and then perform quality trim if asked (use -q 0 to
 * turn off).
 *
 * All quality values will deduct a CUTOFF value (specified by the user) and
 * the max sum segment within the quality string will then be used as the final
 * `trimmed` region.
 *
 * Input is a fastq file, output is a trimmed fastq file.
 */

#define SEQAN_PROFILE // enable time measurements
#include <iostream>
#include <fstream>

#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

using namespace seqan;
using namespace std;

struct Options
{
    int score;
    int times;
    int quality_cutoff;
    int minimum_length;
    int quality_encoding;
    bool discard_adapter_reads;
};

const int SangerOffset = 33;
const int IlluminaOffset = 64;
// score, times, quality_cutoff, minimum_length, quality_encoding,
// discard_adapter_reads
Options DEFAULTS = {15, 4, 20, 30, IlluminaOffset, false}, opts;

int qualityTrim(Dna5String &seq, CharString &qual,
                unsigned &trimStart, unsigned &trimEnd,
                int deduction)
{
    /* Maximum subarray problem, described in wiki
     * using Kadane's algorithm
     *
     * http://en.wikipedia.org/wiki/Maximum_subarray_problem
     *
     */
    unsigned seqlen = length(seq);
    assert (seqlen==length(qual));

    int qv;
    int maxSum = 0; // we are not interested in negative sum
    trimStart = trimEnd = 0;
    int currentMaxSum = 0;
    int currentStartIndex = 0;

    for (unsigned j = 0; j < seqlen; j++)
    {
        // Here we deduct user specified cutoff
        qv = (int)(ordValue(qual[j]) - deduction);
        // Any base with qv higher than cutoff get +1 score, and any base with
        // qv lower than cutoff get -1 score. The target is to get a max-sum
        // segment for the 'modified' score array.
        qv = (qv > 0) - (qv < 0);

        currentMaxSum += qv;

        if (currentMaxSum > maxSum)
        {
            maxSum = currentMaxSum;
            trimStart = currentStartIndex;
            trimEnd = j;
        }
        else if (currentMaxSum < 0)
        {
            currentMaxSum = 0;
            currentStartIndex = j + 1;
        }
    }

    return maxSum;
}

int toSangerQuality(CharString &qual, CharString &qualSanger, int offsetDiff)
{
    for (unsigned j = 0; j < length(qual); j++)
        append(qualSanger, qual[j] + offsetDiff);

    return 0;
}

int main (int argc, char const * argv[])
{
    SEQAN_PROTIMESTART(loadTime);
    ArgumentParser p("trimReads");
    addOption(p, ArgParseOption("o", "outfile",
                                   "Output file name. "
                                   "(default replace suffix with .trimmed.fastq).",
                                   ArgParseArgument::OUTPUTFILE, "OUT"));
    addOption(p, ArgParseOption("f", "adapterfile",
                                   "FASTA formatted file containing the adapters for removal.",
                                   ArgParseArgument::STRING, "STRING"));
    addOption(p, ArgParseOption("s", "score",
                                   "Minimum score to call adapter match. "
                                   "Default scoring scheme for +1 match, "
                                   "-3 for mismatch/gapOpen/gapExtension.",
                                   ArgParseArgument::INTEGER, "INT"));
    addOption(p, ArgParseOption("q", "quality-cutoff",
                                   "Trim low-quality regions below quality cutoff. "
                                   "The algorithm is similar to the one used by BWA "
                                   "by finding a max-sum segment within the quality string. "
                                   "Set it to 0 to skip quality trimming. ",
                                   ArgParseArgument::INTEGER, "INT"));
    addOption(p, ArgParseOption("m", "minimum-length",
                                   "Discard trimmed reads that are shorter than LENGTH.",
                                   ArgParseArgument::INTEGER, "INT"));
    addOption(p, ArgParseOption("Q", "quality-encoding",
                                   "Read quality encoding for input file. 64 for Illumina, "
                                   "33 for Sanger. ",
                                   ArgParseArgument::INTEGER, "INT"));
    addOption(p, ArgParseOption("d", "discard-adapter-reads",
                                   "Discard reads with adapter sequences rather than trim."));

    setDefaultValue(p, "adapterfile", "adapters.fasta");
    setDefaultValue(p, "score", DEFAULTS.score);
    setDefaultValue(p, "quality-cutoff", DEFAULTS.quality_cutoff);
    setDefaultValue(p, "minimum-length", DEFAULTS.minimum_length);
    setDefaultValue(p, "quality-encoding", DEFAULTS.quality_encoding);

    addArgument(p, ArgParseArgument(ArgParseArgument::INPUTFILE, "IN"));

    setShortDescription(p, "Illumina reads trimming utilities");
    setVersion(p, "1.1");
    setDate(p, "July 2013");

    addUsageLine(p, "[options] fastqfile");
    addDescription(p, "Author: Haibao Tang <htang@jcvi.org>");

    ArgumentParser::ParseResult res = parse(p, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;

    CharString infile, outfile, adapterfile;
    getArgumentValue(infile, p, 0);

    if (isSet(p, "outfile"))
    {
        getOptionValue(outfile, p, "outfile");
    }
    else
    {
        // replace the .suffix with .trimmed.fastq
        CharString ts = ".trimmed.fastq";
        for (unsigned i = 0; i < length(infile); i++)
        {
            if (infile[i] == '.') break;
            appendValue(outfile, infile[i]);
        }
        outfile += ts;
    }

    opts = DEFAULTS;
    getOptionValue(adapterfile, p, "adapterfile");
    getOptionValue(opts.score, p, "score");
    getOptionValue(opts.quality_cutoff, p, "quality-cutoff");
    getOptionValue(opts.minimum_length, p, "minimum-length");
    getOptionValue(opts.quality_encoding, p, "quality-encoding");
    if (isSet(p, "discard-adapter-reads")) opts.discard_adapter_reads = true;

    SequenceStream adapter_strm(toCString(adapterfile));
    SequenceStream seq_strm(toCString(infile));
    // Adapter library
    StringSet<CharString> adapterNames, adapters;
    CharString id, qual;
    Dna5String seq;
    readAll(adapterNames, adapters, adapter_strm);

    unsigned nadapters = length(adapters);
    unsigned seqCount = 0;
    unsigned tooShorts = 0;

    String<unsigned> startpos(nadapters); // keep track of adapter positions in concat string
    String<unsigned> adapterCounts(nadapters); // keep track of adapter counts

    // Concatenate all adapters
    unsigned pos = 0;
    startpos[0] = pos;
    CharString adapterdb = adapters[0];

    CharString Ns = "XXXXX"; // use X's to break alignments across adapters
    for (unsigned i = 1; i < nadapters; i++)
    {
        append(adapterdb, Ns);
        startpos[i] = length(adapterdb);
        append(adapterdb, adapters[i]);
    }

    for (unsigned i = 0; i < nadapters; i++)
        adapterCounts[i] = 0;

    Score<int> scoring(1, -3, -3, -3); // harsh penalty for mismatch and indel
    Align<CharString > ali;
    resize(rows(ali), 2);

    int deduction = opts.quality_encoding + opts.quality_cutoff;

    ofstream fout(toCString(outfile));
    unsigned int trimmed_reads = 0, discarded_adapter_reads = 0;

    if (!isGood(seq_strm))
    {
        std::cerr << "ERROR: Could not open " << infile << endl;
        return 1;
    }
    while (!atEnd(seq_strm))
    {
        seqCount ++;
        readRecord(id, seq, qual, seq_strm);
        //cerr << seqCount << "\t" << seq << endl;

        assignSource(row(ali, 0), adapterdb);
        assignSource(row(ali, 1), seq);

        int score = localAlignment(ali, scoring);
        bool containAdapters = (score >= opts.score);

        if (containAdapters)
        {
            unsigned clipStart = clippedBeginPosition(row(ali, 1));
            unsigned clipEnd = clippedEndPosition(row(ali, 1));

            //cout << "Score = " << score << endl;
            //cout << ali << endl;

            for (unsigned j = clipStart; j < clipEnd; j++)
                // mark the adapter region with qual of 1
                qual[j] = (char) (1 + opts.quality_encoding);

            // find out which adapters generated the alignment
            unsigned dbStart = clippedBeginPosition(row(ali, 0));
            unsigned idx;
            for (idx = 0; idx < nadapters; idx++)
            {
                if (startpos[idx] > dbStart) break;
            }
            idx--; // bisect startpos
            adapterCounts[idx]++;
        }

        unsigned trimStart = 0, trimEnd = length(seq) - 1;
        // results are [trimStart, trimEnd] (inclusive on ends)
        if (opts.quality_cutoff > 0)
        {
            qualityTrim(seq, qual, trimStart, trimEnd, deduction);

            if ((trimEnd - trimStart + 1) < (unsigned) opts.minimum_length)
            {
                tooShorts++;
                continue;
            }
        }

        if (opts.discard_adapter_reads && containAdapters)
        {
            discarded_adapter_reads++;
            continue;
        }

        fout << "@" << id << endl;
        fout << infix(seq, trimStart, trimEnd+1) << endl;
        fout << "+" << endl;
        fout << infix(qual, trimStart, trimEnd+1) << endl;
        trimmed_reads++;
    }

    fout.close();

    for (unsigned i = 0; i < nadapters; i++)
    {
        cerr << "[" << i <<"] " << adapterNames[i] << " found "
             << adapterCounts[i] << " times" << endl;
    }

    // Write a report of the trimming
    cerr << endl;
    cerr << "A total of " << tooShorts << " too short (trimmed length < "
         << opts.minimum_length << ") reads removed." << endl;
    if (opts.discard_adapter_reads)
    {
        cerr << "A total of " << discarded_adapter_reads
             << " adapter reads discarded." << endl;
    }
    cerr << "A total of " << trimmed_reads << " trimmed reads are written to `"
         << outfile << "`." << endl;
    cerr << "Processed " << seqCount << " sequences took " << SEQAN_PROTIMEDIFF(loadTime)
         << " seconds." << endl << endl;

    return 0;
}

