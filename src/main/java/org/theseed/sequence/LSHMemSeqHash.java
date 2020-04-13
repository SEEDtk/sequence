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

    @Override
    public boolean checkBucket(int s, int idx) {
        return (this.masterTable[s][idx].size() != 0);
    }
    @Override
    protected void addToBucket(int s, int idx, Sketch sketch) {
        this.masterTable[s][idx].add(sketch);
    }

    @Override
    protected Bucket getBucket(int s, int h) {
        return this.masterTable[s][h];
    }

}
