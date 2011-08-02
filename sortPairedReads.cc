/*
 * Sort Illumina read pairs into adapter set, overlap set and clean set.
 *
 * Author: Haibao Tang <htang@jcvi.org>
 * Date: 01/30/2011
 * License: BSD
 *
 * Input is a pair of fastq file, output is three sets of files:
 *
 * - Set 1: Adapter set (.adapters.fastq) contain read pairs where either /1 or
 *   /2 matches to given adapter.
 * - Set 2: Overlap set (.overlap.fastq) contain read pairs where /1 and /2 has
 *   dovetail overlap.
 * - Set 3: Surviving pairs (.clean.fastq).
 * 
 * This routine will detect the pairs of reads that contain adapters and throw
 * away both ends. Also detect the pairs of reads that have dovetail alignments
 * and push them into fragments.
 *
 * Paired end libraries often have a mixture of various fragment sizes and will
 * confuse some genome assemblers like Celera Assemblers and might make sense to
 * separate the small fragments.
 */

#define SEQAN_PROFILE // enable time measurements
#include <iostream>
#include <fstream>
#include <cassert>

#include <seqan/file.h>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/misc/misc_cmdparser.h>

using namespace seqan;
using namespace std;

struct Options
{
    int adapterMatchScore;
    int endMatchScore;
    int quality_encoding;
    bool verbose;
    bool nooverlap;
};

const int SangerOffset = 33;
const int IlluminaOffset = 64;
Options DEFAULTS = {15, 20, IlluminaOffset, false, false}, opts;

void replaceSuffix(CharString &infile, CharString &outfile, CharString &newSuffix)
{
    // replace the suffix of infile
    outfile = infile;
    unsigned suffix_idx;
    for (suffix_idx = length(infile)-1; suffix_idx > 0; suffix_idx--)
        if (infile[suffix_idx] == '.') break;
    if (suffix_idx==0) suffix_idx = length(outfile);
    suffix(outfile, suffix_idx) = newSuffix;
}

void writeFastQ(ofstream &fout, CharString &id, Dna5String &seq, CharString &qual)
{
    fout << "@" << id << endl << seq << endl;
    fout << "+" << endl << qual << endl;
}

// Dovetail alignment (or end-to-end alignment) between two read sequences
bool alignReads(CharString &id1, Dna5String &seq1, CharString &id2, Dna5String seq2, 
        Score<int>& scoring, int endMatchScore, bool verbose)
{
    AlignConfig<true,true,true,true> ac;
    Align<Dna5String > ali;
    resize(rows(ali), 2);

    reverseComplement(seq2);

    assignSource(row(ali, 0), seq1);
    assignSource(row(ali, 1), seq2);

    int score = globalAlignment(ali, scoring, ac, NeedlemanWunsch());

    if (score < endMatchScore) return false;

    if (verbose)
    {
        cerr << "[D] " << id1 << " vs " << id2 << ": " 
             << "score=" << score << endl << ali << endl;
    }

    return true;
}

// Adapters are called through Local alignment between the adapter database and
// the read sequence
bool alignAdapters(CharString &id, Dna5String &seq, CharString &qual,
        Align<Dna5String> ali, Score<int> &scoring, int adapterMatchScore,
        String<unsigned> &adapterCounts, String<unsigned> &startpos, 
        int quality_encoding, unsigned nadapters, bool verbose)
{
    assignSource(row(ali, 1), seq);
    LocalAlignmentFinder<> finder(ali);

    //int score = localAlignment(ali, scoring, WatermanEggert()); 
    int score = localAlignment(ali, finder, scoring, adapterMatchScore);
    if (score == 0) return false;

    score = getScore(finder);

    if (verbose)
    {
        cerr << "[A] " << id << " vs adapters: " 
             << "score=" << score << endl << ali << endl;
    }

    unsigned clipStart = clippedBeginPosition(row(ali, 1));
    unsigned clipEnd = clippedEndPosition(row(ali, 1));

    for (unsigned j = clipStart; j < clipEnd; j++)
        qual[j] = (char) (1 + quality_encoding);

    // find out which adapters generated the alignment
    unsigned dbStart = clippedBeginPosition(row(ali, 0));
    unsigned idx;
    for (idx = 0; idx < nadapters; idx++)
    {
        if (startpos[idx] > dbStart) break;
    }
    //cerr << dbStart << " " << startpos[idx] << " " << idx << endl; 
    idx--; // bisect startpos
    adapterCounts[idx]++;

    return true;
}

int main (int argc, char const * argv[])
{
    SEQAN_PROTIMESTART(loadTime);
    CommandLineParser p("sortPairedReads");
    addOption(p, CommandLineOption('O', "nooverlap",
                                   "Turn off overlapping reads detection, "
                                   "and do not create .overlap.fastq file.",
                                   OptionType::Bool, DEFAULTS.nooverlap));
    addOption(p, CommandLineOption('f', "adapterfile",
                                   "FASTA formatted file containing the adapters for removal ",
                                   OptionType::String, "adapters.fasta"));
    addOption(p, CommandLineOption('s', "adapterMatchScore",
                                   "Minimum score to call adapter match. "
                                   "Default scoring scheme for +1 match, "
                                   "-3 for mismatch/gapOpen/gapExtension.",
                                   OptionType::Int, DEFAULTS.adapterMatchScore));
    addOption(p, CommandLineOption('t', "endMatchScore",
                                   "Minimum score to call dovetail match. "
                                   "Default scoring scheme for +1 match, "
                                   "-3 for mismatch/gapOpen/gapExtension.",
                                   OptionType::Int, DEFAULTS.endMatchScore));
    addOption(p, CommandLineOption('Q', "quality-encoding",
                                   "Read quality encoding for input file. 64 for Illumina, "
                                   "33 for Sanger. ",
                                   OptionType::Int, DEFAULTS.quality_encoding));
    addOption(p, CommandLineOption('v', "verbose",
                                   "Print alignments for debugging ",
                                   OptionType::Bool, DEFAULTS.verbose));

    addTitleLine(p, "Sort pairs of Illumina reads");
    addTitleLine(p, "Author: Haibao Tang <htang@jcvi.org>");
    addUsageLine(p, "[options] fastqfile1 fastqfile2");

    if (!parse(p, argc, argv))
        return 1;

    String<CharString> args = getArgumentValues(p);
    if (length(args) != 2 || isSetShort(p, 'h'))
    {
        help(p, cerr);
        return 0;
    }

    CharString infile1 = args[0]; // /1
    CharString infile2 = args[1]; // /2
    CharString adapterfile; // the input adapter sequences

    CharString outfile1, outfile2; // sequences that are good
    CharString hasadaptersfile; // contain adapters, short fragments (< read size)
    CharString fragfile; // two ends overlap, short fragments (< 2X read size)

    CharString outSuffix = ".clean.fastq";
    CharString hasadaptersSuffix = ".adapters.fastq";
    CharString fragSuffix = ".overlap.fastq";

    replaceSuffix(infile1, outfile1, outSuffix);
    replaceSuffix(infile2, outfile2, outSuffix);
    replaceSuffix(infile1, hasadaptersfile, hasadaptersSuffix);
    replaceSuffix(infile1, fragfile, fragSuffix);

    opts = DEFAULTS;
    getOptionValueLong(p, "adapterfile", adapterfile);
    getOptionValueLong(p, "adapterMatchScore", opts.adapterMatchScore);
    getOptionValueLong(p, "endMatchScore", opts.endMatchScore);

    if (isSetShort(p, 'v')) opts.verbose = true;
    if (isSetShort(p, 'O')) opts.nooverlap = true;

    MultiSeqFile multiSeqFile1, multiSeqFile2, adapterFile;
    if (argc < 2 || !open(multiSeqFile1.concat, toCString(infile1), OPEN_RDONLY) 
                 || !open(multiSeqFile2.concat, toCString(infile2), OPEN_RDONLY)
                 || !open(adapterFile.concat, toCString(adapterfile), OPEN_RDONLY))
        return 1;

    AutoSeqFormat format;
    // Guess the format of the input, although currently only fastq is supported
    guessFormat(multiSeqFile1.concat, format);
    split(multiSeqFile1, format);
    split(multiSeqFile2, format);
    split(adapterFile, Fasta());

    unsigned seqCount = length(multiSeqFile1);
    unsigned seqCount2 = length(multiSeqFile2);
    assert (seqCount == seqCount2);

    unsigned nadapters = length(adapterFile);

    // Adapter library
    StringSet<CharString> adapterNames;
    StringSet<Dna5String> adapters;
    CharString id;
    Dna5String seq;

    for (unsigned i = 0; i < nadapters; i++)
    {
        assignSeqId(id, adapterFile[i], Fasta());   // read sequence id
        assignSeq(seq, adapterFile[i], Fasta());    // read sequence
        appendValue(adapterNames, id);
        appendValue(adapters, seq);
    }

    String<unsigned> startpos(nadapters); // keep track of adapter positions in concat string
    String<unsigned> adapterCounts(nadapters); // keep track of adapter counts

    // Concatenate all adapters
    unsigned pos = 0;
    startpos[0] = pos;
    Dna5String adapterdb = adapters[0];

    Dna5String Ns = "NNNNN"; // use N's to break alignments across adapters
    for (unsigned i = 1; i < nadapters; i++)
    {
        append(adapterdb, Ns);
        startpos[i] = length(adapterdb);
        append(adapterdb, adapters[i]);
    }

    for (unsigned i = 0; i < nadapters; i++)
        adapterCounts[i] = 0;

    Score<int> scoring(1, -3, -3, -3); // harsh penalty for mismatch and indel
    Align<Dna5String > ali;
    resize(rows(ali), 2);
    assignSource(row(ali, 0), adapterdb);

    CharString id1, qual1;
    CharString id2, qual2;
    Dna5String seq1, seq2;

    ofstream fout1(toCString(outfile1));
    ofstream fout2(toCString(outfile2));
    ofstream fouta(toCString(hasadaptersfile));
    ofstream foutf;
    if (! opts.nooverlap)
        foutf.open(toCString(fragfile));

    for (unsigned i = 0; i < seqCount; i++)
    {
        assignSeqId(id1, multiSeqFile1[i], format);   // read sequence id
        assignSeq(seq1, multiSeqFile1[i], format);    // read sequence
        assignQual(qual1, multiSeqFile1[i], format);  // read ascii quality values

        assignSeqId(id2, multiSeqFile2[i], format);   // read sequence id
        assignSeq(seq2, multiSeqFile2[i], format);    // read sequence
        assignQual(qual2, multiSeqFile2[i], format);  // read ascii quality values

        bool r1hasAdapter = alignAdapters(id1, seq1, qual1,
                ali, scoring, opts.adapterMatchScore, 
                adapterCounts, startpos,
                opts.quality_encoding, nadapters, opts.verbose);

        bool r2hasAdapter = alignAdapters(id2, seq2, qual2, 
                ali, scoring, opts.adapterMatchScore, 
                adapterCounts, startpos,
                opts.quality_encoding, nadapters, opts.verbose);

        // Here is the sorting logic, first check existence of adapters
        // then check the /1 and /2 read overlap
        //
        if (r1hasAdapter || r2hasAdapter)
        {
            writeFastQ(fouta, id1, seq1, qual1);
            writeFastQ(fouta, id2, seq2, qual2);
            continue;
        }

        if ((! opts.nooverlap) && alignReads(id1, seq1, id2, seq2,
                scoring, opts.endMatchScore, opts.verbose))
        {
            writeFastQ(foutf, id1, seq1, qual1);
            writeFastQ(foutf, id2, seq2, qual2);
        }
        else
        {
            writeFastQ(fout1, id1, seq1, qual1);
            writeFastQ(fout2, id2, seq2, qual2);
        }
    }

    // Report the adapter matches
    for (unsigned i = 0; i < nadapters; i++)
    {
        cerr << "[" << i <<"] " << adapterNames[i] << " found "
             << adapterCounts[i] << " times" << endl;
    }

    cerr << endl;
    cerr << "Processed " << seqCount << " sequences took " << SEQAN_PROTIMEDIFF(loadTime)
         << " seconds." << endl << endl;

    return 0;
}

