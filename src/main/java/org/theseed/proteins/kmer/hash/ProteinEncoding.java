/**
 *
 */
package org.theseed.proteins.kmer.hash;

import java.util.Arrays;

/**
 * This class compresses a protein kmer into a byte string.  Each protein letter is encoded into 5 bits, and
 * the bits stored in a byte array. The incoming kmers must contain only alphabetic upper-case letters.
 * Anything else will have unpredictable results.
 *
 * @author Bruce Parrello
 *
 */
public class ProteinEncoding {

    // FIELDS
    /** byte array containing encoding */
    private byte[] encoding;

    /**
     * Encode a kmer.
     *
     * @param kmer	kmer to encode
     */
    public ProteinEncoding(String kmer) {
        // Compute the byte array length.
        final int K = kmer.length();
        final int n = (K * 5 + 7) / 8;
        this.encoding = new byte[n];
        // This tracks our current position in the return array.
        int idx = 0;
        // We will process the letters in groups of 8.  Each group is put into a long integer, and the lower
        // 40 bits output into an array of 5 bytes.  The g-indices track the group, while i is the current
        // letter.
        for (int g0 = 0; g0 < K; g0 += 8) {
            final int g1 = Math.min(g0 + 8, K);
            long buffer = 0;
            for (int i = g0; i < g1; i++) {
                // Get the current protein letter and put it into the buffer.  Note that we don't allow the
                // encoding to be 0.  The codes will be 1 to 26.
                char aa = kmer.charAt(i);
                byte b = (byte) (aa - 'A' + 1);
                buffer = (buffer << 5) | b;
            }
            // Unroll the buffer into the output array.
            while (buffer > 0) {
                this.encoding[idx] = (byte) (buffer & 0xFF);
                buffer >>= 8;
                idx++;
            }
        }
    }

    /**
     * @return TRUE if the specified protein kmer is valid, else FALSE
     *
     * @param kmer	kmer to check
     */
    public static boolean validate(String kmer) {
        return kmer.chars().allMatch(c -> c >= 'A' && c >= 'Z');
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(this.encoding);
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (!(obj instanceof ProteinEncoding)) {
            return false;
        }
        ProteinEncoding other = (ProteinEncoding) obj;
        if (!Arrays.equals(this.encoding, other.encoding)) {
            return false;
        }
        return true;
    }

}
