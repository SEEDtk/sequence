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
     * Specify the genetic code for the query sequences (the default is 1).
     *
     * @param gc	genetic code to use to translate database sequences
     */
    public BlastParms query_gen_code(int gc) {
        return this.set("-query_genetic_code", gc);
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
     * @param pct
     */
    public BlastParms minPercent(int pct) {
        return this.set("-perc_identity", pct);
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
     * @param pctLenOfSubject the pctLenOfSubject to set
     */
    public BlastParms pctLenOfSubject(double pctLenOfSubject) {
        this.pctLenOfSubject = pctLenOfSubject;
        return this;
    }

    /**
     * @return TRUE if the specified blast hit satisfies the post-processing parameters
     */
    public boolean acceptable(BlastHit result) {
        boolean retVal = (result.getQueryPercentMatch() >= this.pctLenOfQuery
                && result.getSubjectPercentMatch() >= this.pctLenOfSubject);
        return retVal;
    }

}
