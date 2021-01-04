/**
 *
 */
package org.theseed.sequence;

import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This class computes the tetramer profile for a DNA sequence.  A tetramer is represented by 8 bits (2 per nucleotide), and a
 * reverse compliment has the same value as the original string.  Because some sequences are palindromes, this nets us 136 different
 * possible tetramers.
 *
 * @author Bruce Parrello
 *
 */
public class TetramerProfile {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(TetramerProfile.class);

    /** normalized tetramer profile */
    private double[] profile;
    /** map of tetramer numbers to tetramer indices */
    protected static final int[] TETRA_MAP = createTetramerMap();
    /** number of distinct tetramer indices */
    protected static int TETRA_SIZE;

    /**
     * @return the index number (from 0 to 256) of a tetramer at the specified location, or -1 if the location is invalid
     *
     * @param dna		DNA sequence containing the tetramer
     * @param loc		location (1-based) of the tetramer in the sequence
     */
    public static int tetraNum(String dna, int loc) {
        int retVal = 0;
        int n = loc + 3;
        if (n  > dna.length())
            retVal = -1;
        else for (int i = loc - 1; i < n && retVal >= 0; i++) {
            retVal <<= 2;
            switch (dna.charAt(i)) {
            case 'A':
            case 'a':
                break;
            case 'C':
            case 'c':
                retVal += 1;
                break;
            case 'G':
            case 'g':
                retVal += 2;
                break;
            case 'T':
            case 't':
            case 'U':
            case 'u':
                retVal += 3;
                break;
            default :
                retVal = -1;
            }
        }
        return retVal;
    }

    /**
     * @return the reverse complement of a tetramer index number
     *
     * @param num	tetramer index to reverse complement
     */
    protected static int revTetraNum(int num) {
        return ~((num >> 6) | ((num & 0x3) << 6) | ((num & 0xC) << 2) | ((num & 0x30) >> 2)) & 0xFF;
    }

    /**
     * @return a map from tetramer numbers to the actual tetramer index.
     */
    private static int[] createTetramerMap() {
        int[] retVal = new int[256];
        int idx = 0;
        for (int i = 0; i < 256; i++) {
            int revI = revTetraNum(i);
            if (revI >= i) {
                retVal[i] = idx;
                retVal[revI] = idx;
                idx++;
            }
        }
        TETRA_SIZE = idx;
        return retVal;
    }

    /**
     * Construct the tetramer profile for a DNA sequence.
     *
     * @param dna	dna sequence whose profile is desired
     */
    public TetramerProfile(String dna) {
        this.init(dna, 1, dna.length());
    }

    /**
     * Create a tetramer profile for a subsequence.
     *
     * @param sequence		DNA sequence to use
     * @param pos			starting position (1-based)
     * @param len			length of desired profile
     */
    public TetramerProfile(String sequence, int pos, int len) {
        init(sequence, pos, len);
    }

    /**
     * Initialize this tetramer profile from a subsequence.
     *
     * @param sequence		DNA sequence to use
     * @param pos			starting position (1-based)
     * @param len			length of desired profile
     */
    protected void init(String sequence, int pos, int len) {
        this.profile = new double[TETRA_SIZE];
        Arrays.fill(this.profile, 0.0);
        // Compute the number of tetramers in this sequence.
        int n = pos + len - 1;
        if (n > sequence.length()) n = sequence.length();
        n -= 3;
        // Count each valid one.
        for (int i = pos; i <= n; i++) {
            int tetraNum = tetraNum(sequence, i);
            if (tetraNum >= 0) {
                int tetraIdx = TETRA_MAP[tetraNum];
                this.profile[tetraIdx]++;
            }
        }
        // Normalize the vector.
        if (len > 3) {
            double nDouble = (double) (len - 3);
            for (int i = 0; i < TETRA_SIZE; i++)
                this.profile[i] /= nDouble;
        }
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + Arrays.hashCode(this.profile);
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof TetramerProfile)) {
            return false;
        }
        TetramerProfile other = (TetramerProfile) obj;
        if (!Arrays.equals(this.profile, other.profile)) {
            return false;
        }
        return true;
    }

    @Override
    public String toString() {
        return Arrays.stream(this.profile).mapToObj(x -> Double.toString(x)).collect(Collectors.joining("\t"));
    }

    /**
     * @return a header line for tetramer profiles
     */
    public static String headers() {
        return IntStream.range(0, TETRA_SIZE).mapToObj(i -> String.format("T%03d", i)).collect(Collectors.joining("\t"));
    }


}
