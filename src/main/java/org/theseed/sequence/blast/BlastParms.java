/**
 *
 */
package org.theseed.sequence.blast;

import org.theseed.utils.Parms;

/**
 * This class manages the parameters for a Blast operation.  It is a normal Parms object,
 * but with fluent options for setting the common BLAST parameters.  In addition, a couple
 * of parameters required by post-processing are supported.
 *
 * @author Bruce Parrello
 */
public class BlastParms extends Parms implements Cloneable {

    // FIELDS
    /** minimum acceptable length for a result, as a percent of the query length */
    private double pctLenOfQuery;
    /** minimum acceptable length for a result, as a percent of the subject length */
    private double pctLenOfSubject;
    /** minimum acceptable percent identity */
    private double pctIdentity;
    /** minimum acceptable query-scaled bit score */
    private double minQueryBitScore;
    /** minimum acceptable query identity fraction */
    private double minQueryIdentity;
    /** TRUE if we should turn off filtering for blastn (-dust no) */
    private boolean dustOff;
    /** mismatch penalty for blastn (-penalty XX) */
    private int penalty;
    /** match reward for blastn (-reward XX) */
    private int reward;
    /** minimum absolute length for a match */
    private int minMatchLen;

    /**
     * Add an option parameter.  This calls through to the super-class, but returns an
     * instance of this class.
     *
     * @param name	parameter name
     */
    @Override
    public BlastParms set(String option) {
        super.set(option);
        return this;
    }

    /**
     * Add a string parameter.  This calls through to the super-class, but returns an
     * instance of this class.
     *
     * @param name		parameter name
     * @param value		parameter value
     */
    @Override
    public BlastParms set(String option, String value) {
        super.set(option, value);
        return this;
    }

    /**
     * Add a floating-point parameter.  This calls through to the super-class, but returns an
     * instance of this class.
     *
     * @param name		parameter name
     * @param value		parameter value
     */
    @Override
    public BlastParms set(String option, double value) {
        super.set(option, value);
        return this;
    }

    /**
     * Add an integer parameter.  This calls through to the super-class, but returns an
     * instance of this class.
     *
     * @param name		parameter name
     * @param value		parameter value
     */
    @Override
    public BlastParms set(String option, int value) {
        super.set(option, value);
        return this;
    }

    /**
     * Set the defaults for the post-processing parms.
     */
    @Override
    protected void setDefaults() {
        this.pctLenOfQuery = 0.0;
        this.pctLenOfSubject = 0.0;
        this.pctIdentity = 0.0;
        this.minQueryBitScore = 0.0;
        this.minQueryIdentity = 0.0;
        this.minMatchLen = 0;
        this.penalty = 0;
        this.dustOff = false;
        this.reward = 0;
    }

    /**
     * Specify the number of threads (the default is 1).
     *
     * @param threads	proposed thread count
     */
    public BlastParms num_threads(int threads) {
        return this.set("-num_threads", threads);
    }

    /**
     * Specify the genetic code for the database (the default is 1).
     *
     * @param gc	genetic code to use to translate database sequences
     */
    public BlastParms db_gencode(int gc) {
        return this.set("-db_gencode", gc);
    }

    /**
     * Specify the x-dropoff value (in bits).
     *
     * @param drop	gap for X dropoff
     */
    public BlastParms xdrop_gap(double drop) {
        return this.set("-xdrop_gap", drop);
    }

    /**
     * Specify the genetic code for the query sequences (the default is 1).
     *
     * @param gc	genetic code to use to translate database sequences
     */
    public BlastParms query_gencode(int gc) {
        return this.set("-query_gencode", gc);
    }

    /**
     * Specify the maximum e-value for results (the default is 10).
     *
     * @param evalue	maximum probability for an error
     */
    public BlastParms maxE(double evalue) {
        return this.set("-evalue", evalue);
    }

    /**
     * Specify the maximum number of results to keep per query sequence (the default is 500).
     *
     * @param max		maximum number of results per query
     */
    public BlastParms maxPerQuery(int max) {
        return this.set("-max_target_seqs", max);
    }

    /**
     * Specify the minimum percent identity.
     *
     * @param pct	minimum percent identity
     */
    public BlastParms minPercent(double pct) {
        this.pctIdentity = pct;
        return this;
    }

    /**
     * Denote we want filtering off for blastn.
     */
    public BlastParms dustOff() {
        this.dustOff = true;
        return this;
    }

    /**
     * Specify a blastn mismatch penalty.
     *
     * @param pen	penalty value (usually negative)
     */
    public BlastParms penalty(int pen) {
        this.penalty = pen;
        return this;
    }

    /**
     * Specify a blastn match reward.
     *
     * @param rwd	reward value (usually positive)
     */
    public BlastParms reward(int rwd) {
        this.reward = rwd;
        return this;
    }

    /**
     * Add special parameters for blastn.  This is always called on a clone of the client's original
     * parms object.
     */
    public void customizeForBlastn() {
        if (this.dustOff)
            this.set("-dust", "no");
        if (this.penalty != 0)
            this.set("-penalty", Integer.toString(this.penalty));
        if (this.reward != 0)
            this.set("-reward", Integer.toString(this.reward));
    }

    /**
     * @return a copy of this object.
     */
    @Override
    public BlastParms clone() {
        BlastParms retVal = new BlastParms();
        this.copyValues(retVal);
        retVal.pctLenOfQuery = this.pctLenOfQuery;
        retVal.pctLenOfSubject = this.pctLenOfSubject;
        retVal.pctIdentity = this.pctIdentity;
        retVal.minQueryBitScore = this.minQueryBitScore;
        retVal.minQueryIdentity = this.minQueryIdentity;
        return retVal;
    }

    /**
     * @return the pctLenOfQuery
     */
    public double getPctLenOfQuery() {
        return pctLenOfQuery;
    }

    /**
     * @param pctLenOfQuery the pctLenOfQuery to set
     */
    public BlastParms pctLenOfQuery(double pctLenOfQuery) {
        this.pctLenOfQuery = pctLenOfQuery;
        return this;
    }

    /**
     * @return the pctLenOfSubject
     */
    public double getPctLenOfSubject() {
        return pctLenOfSubject;
    }

    /**
     * @return the minimum match length
     */
    public int getMinLen() {
        return this.minMatchLen;
    }

    /**
     * @param minLen	the minimum match length to accept
     */
    public BlastParms minLen(int minLen) {
        this.minMatchLen = minLen;
        return this;
    }

    /**
     * @param pctLenOfSubject the pctLenOfSubject to set
     */
    public BlastParms pctLenOfSubject(double pctLenOfSubject) {
        this.pctLenOfSubject = pctLenOfSubject;
        return this;
    }

    /**
     * @param minQueryIdentity the minimum acceptable query identity to set
     */
    public BlastParms minQueryIdentity(double minQueryIdentity) {
        this.minQueryIdentity = minQueryIdentity;
        return this;
    }

    /**
     * @param minQueryBitScore the minimum acceptable query-scaled bit score to set
     */
    public BlastParms minQueryBitScore(double minQueryBitScore) {
        this.minQueryBitScore = minQueryBitScore;
        return this;
    }

    /**
     * @return TRUE if the specified blast hit satisfies the post-processing parameters
     */
    public boolean acceptable(BlastHit result) {
        boolean retVal = (result.getQueryPercentMatch() >= this.pctLenOfQuery
                && result.getSubjectPercentMatch() >= this.pctLenOfSubject
                && result.getPercentIdentity() >= this.pctIdentity
                && result.getQueryBitScore() >= this.minQueryBitScore
                && result.getQueryIdentity() >= this.minQueryIdentity
                && result.getAlignLen() >= this.minMatchLen);
        return retVal;
    }

    /**
     * @return the pctIdentity
     */
    public double getPctIdentity() {
        return pctIdentity;
    }

    @Override
    public String toString() {
        return super.toString() + (pctLenOfQuery > 0.0 ? " --pctLenOfQuery=" + pctLenOfQuery : "") +
                (pctLenOfSubject > 0.0 ? " --pctLenOfSubject=" + pctLenOfSubject : "") +
                (pctIdentity > 0.0 ? " --pctIdentity=" + pctIdentity : "");
    }

    /**
     * @return the minQueryBitScore
     */
    public double getMinQueryBitScore() {
        return minQueryBitScore;
    }

    /**
     * @return the minQueryIdentity
     */
    public double getMinQueryIdentity() {
        return minQueryIdentity;
    }

}
