/**
 *
 */
package org.theseed.sequence;

/**
 * This is a version of the locally-sensitive minhash that stores all its data in memory.
 *
 * @author Bruce Parrello
 *
 */
public class LSHMemSeqHash extends LSHSeqHash {

    // FIELDS

    /** master hash table (first index is stage, second is bucket */
    private Bucket[][] masterTable;

    /**
     * Construct a blank, empty sequence hash.
     *
     * @param w		width of signatures
     * @param s		number of stages
     * @param b		number of buckets
     */
    public LSHMemSeqHash(int w, int s, int b) {
        super(w, s, b);
        this.masterTable = new Bucket[s][b];
        for (int si = 0; si < s; si++)
            for (int bi = 0; bi < b; bi++)
                this.masterTable[si][bi] = new Bucket();
    }

    /**
     * @return the number of entries in this table
     */
    public int size() {
        int retVal = 0;
        for (Bucket bucket : this.masterTable[0])
            retVal += bucket.size();
        return retVal;
    }

    @Override
    public boolean checkBucket(int s, int idx) {
        return (this.masterTable[s][idx].size() != 0);
    }

    @Override
    protected Bucket getBucket(int s, int h) {
        return this.masterTable[s][h];
    }

}
