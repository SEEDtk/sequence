/**
 *
 */
package org.theseed.sequence.blast;

import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.theseed.locations.Location;

/**
 * This class represents a BLAST hit.
 *
 * @author Bruce Parrello
 *
 */
public class BlastHit {

    // FIELDS
    /** query sequence definition (comment) */
    private String queryDef;
    /** total query sequence length */
    private int queryLen;
    /** location of matching query sequence */
    private Location queryLoc;
    /** matching query sequence */
    private String querySeq;
    /** subject sequence definition (comment) */
    private String subjectDef;
    /** total subject sequence length */
    private int subjectLen;
    /** location of matching subject sequence */
    private Location subjectLoc;
    /** matching subject sequence */
    private String subjectSeq;
    /** match bit-score */
    private int bitScore;
    /** identity count */
    private int numIdentical;
    /** gap count */
    private int numGap;
    /** error value */
    private double evalue;
    /** number of positive (non-identical) matches */
    private int positives;
    /** alignment length */
    private int alignLen;
    /** query sequence type */
    private BlastDB.Type queryType;
    /** subject sequence type */
    private BlastDB.Type subjectType;
    /** required output format */
    public static final String OUT_FORMAT = "6 qseqid qlen qstart qend qseq stitle slen sstart send " +
                                            // 0      1    2      3    4    5      6    7      8
                                              "sseq bitscore nident gaps evalue positive length";
                                            // 9    10       11     12   13     14       15

    /**
     * Create a BLAST hit from an input line.
     *
     * @param line		input line to parse
     * @param qMap		map of query IDs to query definitions
     * @param qType		query type
     * @param sType		subject type
     */
    public BlastHit(String line, Map<String, String> qMap, BlastDB.Type queryType, BlastDB.Type subjectType) {
        String[] fields = StringUtils.split(line, '\t');
        String qid = fields[0];
        this.queryDef = qMap.getOrDefault(qid, "");
        this.queryLen = Integer.valueOf(fields[1]);
        this.queryLoc = Location.create(qid, Integer.valueOf(fields[2]), Integer.valueOf(fields[3]));
        this.querySeq = fields[4];
        String[] pieces = StringUtils.split(fields[5], " ", 2);
        this.subjectDef = pieces[1];
        this.subjectLen = Integer.valueOf(fields[6]);
        this.subjectLoc = Location.create(pieces[0], Integer.valueOf(fields[7]), Integer.valueOf(fields[8]));
        this.subjectSeq = fields[9];
        this.bitScore = Integer.valueOf(fields[10]);
        this.numIdentical = Integer.valueOf(fields[11]);
        this.numGap = Integer.valueOf(fields[12]);
        this.evalue = Double.valueOf(fields[13]);
        this.positives = Integer.valueOf(fields[14]);
        this.alignLen = Integer.valueOf(fields[15]);
        this.queryType = queryType;
        this.subjectType = subjectType;
    }

    /**
     * @return the query sequence description (comment)
     */
    public String getQueryDef() {
        return queryDef;
    }

    /**
     * @return the total length of the query sequence
     */
    public int getQueryLen() {
        return queryLen;
    }

    /**
     * @return the location of the matching portion of the query sequence
     */
    public Location getQueryLoc() {
        return queryLoc;
    }

    /**
     * @return the matching portion of the query sequence, in alignment format
     */
    public String getQuerySeq() {
        return querySeq;
    }

    /**
     * @return the subject sequence description (comment)
     */
    public String getSubjectDef() {
        return subjectDef;
    }

    /**
     * @return the total length of the subject sequence
     */
    public int getSubjectLen() {
        return subjectLen;
    }

    /**
     * @return the location of the matching portion of the subject sequence
     */
    public Location getSubjectLoc() {
        return subjectLoc;
    }

    /**
     * @return the matching portion of the subject sequence, in alignment format
     */
    public String getSubjectSeq() {
        return subjectSeq;
    }

    /**
     * @return the bit score of this match, which is the internal estimate of the quality of the match,
     * normalized with respect to the properties of the scoring system (but not the query length)
     */
    public int getBitScore() {
        return bitScore;
    }

    /**
     * @return the number of identical character positions
     */
    public int getNumIdentical() {
        return numIdentical;
    }

    /**
     * @return the number of gap positions
     */
    public int getNumGap() {
        return numGap;
    }

    /**
     * @return the probability that the match is in error
     */
    public double getEvalue() {
        return evalue;
    }

    /**
     * @return the number of non-identical matching positions
     */
    public int getPositives() {
        return positives;
    }

    /**
     * @return the match length as a percent of the query length
     */
    public double getQueryPercentMatch() {
        double retVal = 0.0;
        if (this.queryLen > 0) {
            double qLen = (double) this.queryLen;
            if (this.queryType == BlastDB.Type.DNA && this.subjectType == BlastDB.Type.PROTEIN)
                qLen /= 3.0;
            retVal = (this.positives + this.numIdentical) * 100 / qLen;
        }
        return retVal;
    }

    /**
     * @return the match length as a percent of the subject length
     */
    public double getSubjectPercentMatch() {
        double retVal = 0.0;
        if (this.subjectLen > 0) {
            double sLen = (double) this.subjectLen;
            if (this.subjectType == BlastDB.Type.DNA && this.queryType == BlastDB.Type.PROTEIN)
                sLen /= 3.0;
            retVal = (this.positives + this.numIdentical) * 100 / sLen;
        }
        return retVal;
    }

    /**
     * @return the alignment length
     */
    public int getAlignLen() {
        return alignLen;
    }

    /**
     * @return the percent identity
     */
    public double getPercentIdentity() {
        double retVal = 0.0;
        if (this.alignLen > 0)
            retVal = (double) (this.numIdentical * 100) / this.alignLen;
        return retVal;
    }

    /**
     * @return the percent similarity
     */
    public double getPercentSimilarity() {
        double retVal = 0.0;
        if (this.alignLen > 0)
            retVal = (double) (this.numIdentical + this.positives) * 100 / this.alignLen;
        return retVal;
    }

    /**
     * @return the query sequence ID
     */
    public String getQueryId() {
        return this.queryLoc.getContigId();
    }

    /**
     * @return the subject sequence ID
     */
    public String getSubjectId() {
        return this.subjectLoc.getContigId();
    }

}
