/**
 *
 */
package org.theseed.reports;

/**
 * This simple object describes the BLAST results.
 */
public class BlastInfo {

    // FIELDS
    /** BLAST parameters */
    private String parms;
    /** query count */
    private int queriesIn;
    /** number of misses */
    private int missCount;
    /** number of hits */
    private int hitCount;

    /**
     * Construct the information object.
     *
     * @param parms		parameters for the BLAST
     * @param queries	number of incoming queries
     * @param misses	number of queries without a hit
     * @param hits		total number of hits
     */
    public BlastInfo(String parms, int queries, int misses, int hits) {
        this.parms = parms;
        this.queriesIn = queries;
        this.missCount = misses;
        this.hitCount = hits;
    }

    /**
     * @return the BLAST parameters
     */
    public String getParms() {
        return parms;
    }

    /**
     * @return the input query count
     */
    public int getQueriesIn() {
        return queriesIn;
    }

    /**
     * @return the query miss count
     */
    public int getMissCount() {
        return missCount;
    }

    /**
     * @return the hit count
     */
    public int getHitCount() {
        return hitCount;
    }


}
