/**
 *
 */
package org.theseed.bins;

import org.theseed.counters.CountMap;
import org.theseed.genome.Contig;
import org.theseed.sequence.Sequence;

/**
 * This is a contig filter for binning.  It provides methods for creating a Bin object from a contig
 * sequence and for determining if a bin should be discarded, used for binning only, or used for the
 * SOUR protein search.
 *
 * @author Bruce Parrello
 *
 */
public class ContigFilter {

    // FIELDS
    /** parameter object containing filter parameters */
    private BinParms parms;

    /**
     * Initialize a contig filter.
     *
     * @param binParms	binning parameters
     */
    public ContigFilter(BinParms binParms) {
        this.parms = binParms;
    }

    /**
     * Create a bin from a contig sequence.
     *
     * @param contig		contig sequence to parse for length and coverage
     * @param counters		count map to record reasons for bin rejection
     *
     * @return a bin object for the contig
     */
    public Bin computeBin(Sequence contig, CountMap<String> counters) {
        // Compute the sequence length.
        int len = contig.length();
        double coverage = Contig.computeCoverage(contig);
        // We have the coverage and the length.  Build the bin.
        Bin retVal = new Bin(contig.getLabel(), len, coverage);
        // Compute what to do with the bin.
        Bin.Status status;
        if (len < this.parms.getBinLenFilter()) {
            status = Bin.Status.BAD;
            counters.count("contig-bad-TooShort");
        } else if (coverage < this.parms.getBinCovgFilter()) {
            status = Bin.Status.BAD;
            counters.count("contig-bad-LowCoverage");
        } else if (this.isAmbiguous(contig.getSequence(), parms.getXLimit())) {
            status = Bin.Status.BAD;
            counters.count("contig-bad-AmbiguityRun");
        } else if (len < this.parms.getLenFilter()) {
            status = Bin.Status.NORMAL;
            counters.count("contig-normal-TooShort");
        } else if (coverage < this.parms.getCovgFilter()) {
            status = Bin.Status.NORMAL;
            counters.count("contig-normal-LowCoverage");
        } else {
            status = Bin.Status.SEED_USABLE;
            counters.count("contig-seed-usable");
        }
        // Save the status and return the bin.
        retVal.setStatus(status);
        return retVal;
    }

    /**
     * Verify there are no long ambiguity-character runs in this sequence.
     *
     * @param sequence	DNA sequence to check
     * @param xLimit	maximum limit of an ambiguity run
     *
     * @return TRUE if a bad run was found, else FALSE
     */
    protected boolean isAmbiguous(String sequence, int xLimit) {
        final int n = sequence.length();
        // This will track the number of ambiguity characters in the current group.
        int runLength = 0;
        // When we hit a good character, we set run-length to 0.  Otherwise we increment it, and fail
        // when the run gets too long.
        for (int i = 0; i < n && runLength < xLimit; i++) {
            switch (sequence.charAt(i)) {
            case 'a':
            case 'A':
            case 'c':
            case 'C':
            case 'g':
            case 'G':
            case 't':
            case 'T':
            case 'u':
            case 'U':
                runLength = 0;
                break;
            default:
                runLength++;
            }
        }
        return (runLength >= xLimit);
    }

}
