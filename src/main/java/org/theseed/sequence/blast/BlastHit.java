/**
 *
 */
package org.theseed.sequence.blast;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
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

    /**
     * Comparator for comparing blast hits by location in the query sequence
     * @author parrello
     *
     */
    public static class ByQueryLoc implements Comparator<BlastHit> {

        @Override
        public int compare(BlastHit o1, BlastHit o2) {
            int retVal = o1.getQueryLoc().compareTo(o2.getQueryLoc());
            if (retVal == 0) {
                retVal = o1.getSubjectId().compareTo(o2.getSubjectId());
                if (retVal == 0) {
                    retVal = o2.getNumSimilar() - o1.getNumSimilar();
                    if (retVal == 0) {
                        retVal = Double.compare(o1.getEvalue(), o2.getEvalue());
                        if (retVal == 0) {
                            retVal = o2.getNumIdentical() - o1.getNumIdentical();
                            if (retVal == 0) {
                                retVal = o1.getSubjectLoc().compareTo(o2.getSubjectLoc());
                            }
                        }
                    }
                }
            }
            return retVal;
        }

    }

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
    private double bitScore;
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
    private boolean queryIsProtein;
    /** subject sequence type */
    private boolean subjectIsProtein;
    /** required output format */
    public static final String OUT_FORMAT = "6 qseqid qlen qstart qend qseq stitle slen sstart send " +
                                            // 0      1    2      3    4    5      6    7      8
                                              "sseq bitscore nident gaps evalue positive length";
                                            // 9    10       11     12   13     14       15
    /** headers for output file */
    public static final String PRINT_HEADER =
            "q_id\tq_comment\tq_location\tq_length\t" +
            "s_id\ts_comment\ts_location\ts_length\t" +
            "bit_score\tn_ident\tn_gap\tn_positive\t" +
            "e_value\talign_len\tq_sequence\ts_sequence";

    /**
     * Create a BLAST hit from an input line.
     *
     * @param line		input line to parse
     * @param qMap		map of query IDs to query definitions
     * @param qType		query type
     * @param sType		subject type
     */
    public BlastHit(String line, Map<String, String> qMap, boolean queryIsProtein, boolean subjectIsProtein) {
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
        this.bitScore = Double.valueOf(fields[10]);
        this.numIdentical = Integer.valueOf(fields[11]);
        this.numGap = Integer.valueOf(fields[12]);
        this.evalue = Double.valueOf(fields[13]);
        this.alignLen = Integer.valueOf(fields[15]);
        this.positives = Integer.valueOf(fields[14]) - this.numIdentical;
        this.queryIsProtein = queryIsProtein;
        this.subjectIsProtein = subjectIsProtein;
    }

    /**
     * Sort a list of blast hits by query sequence ID, and map each such ID to a list of
     * the results that represent matches.  For each subject match, only keep the longest.
     *
     * @param results	list of BLAST results to sort
     *
     * @return a map from each query ID to a list of the longest match for each subject
     */
    public static Map<String, List<BlastHit>> sort(List<BlastHit> results) {
        Map<String, List<BlastHit>> retVal = new HashMap<String, List<BlastHit>>(results.size());
        for (BlastHit result : results) {
            String queryId = result.getQueryId();
            List<BlastHit> matches = retVal.get(queryId);
            if (matches == null) {
                matches = new ArrayList<BlastHit>(5);
                matches.add(result);
                retVal.put(queryId, matches);
            } else {
                // Here we need to insure this is the longest match for the result.
                int i = 0;
                String subjectId = result.getSubjectId();
                boolean keep = true;
                while (i < matches.size() && keep) {
                    BlastHit match = matches.get(i);
                    if (! subjectId.contentEquals(match.getSubjectId())) {
                        i++;
                    } else if (result.getNumIdentical() <= match.getNumIdentical()) {
                        keep = false;
                        i++;
                    } else  {
                        matches.remove(i);
                    }
                }
                if (keep)
                    matches.add(result);
            }
        }
        return retVal;
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
    public double getBitScore() {
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
            if (! queryIsProtein && subjectIsProtein)
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
            if (! subjectIsProtein && queryIsProtein)
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
     * @return the number of similar positions
     */
    public int getNumSimilar() {
        return (this.numIdentical + this.positives);
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

    /**
     * @return a tab-delimited print line for this result (with no LF)
     */
    public String getPrintLine() {
        return String.format("%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t" +
                "%4.2f\t%d\t%d\t%d\t%6.3e\t%d\t%s\t%s",
                this.getQueryId(), this.queryDef, this.queryLoc.toString(),
                this.queryLen, this.getSubjectId(), this.subjectDef,
                this.subjectLoc.toString(), this.subjectLen, this.bitScore,
                this.numIdentical, this.numGap, this.positives, this.evalue,
                this.alignLen, this.querySeq, this.subjectSeq);
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + alignLen;
        long temp;
        temp = Double.doubleToLongBits(bitScore);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(evalue);
        result = prime * result + (int) (temp ^ (temp >>> 32));
        result = prime * result + numGap;
        result = prime * result + numIdentical;
        result = prime * result + positives;
        result = prime * result + ((queryDef == null) ? 0 : queryDef.hashCode());
        result = prime * result + (queryIsProtein ? 1231 : 1237);
        result = prime * result + queryLen;
        result = prime * result + ((queryLoc == null) ? 0 : queryLoc.hashCode());
        result = prime * result + ((querySeq == null) ? 0 : querySeq.hashCode());
        result = prime * result + ((subjectDef == null) ? 0 : subjectDef.hashCode());
        result = prime * result + (subjectIsProtein ? 1231 : 1237);
        result = prime * result + subjectLen;
        result = prime * result + ((subjectLoc == null) ? 0 : subjectLoc.hashCode());
        result = prime * result + ((subjectSeq == null) ? 0 : subjectSeq.hashCode());
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        BlastHit other = (BlastHit) obj;
        if (alignLen != other.alignLen)
            return false;
        if (Double.doubleToLongBits(bitScore) != Double.doubleToLongBits(other.bitScore))
            return false;
        if (Double.doubleToLongBits(evalue) != Double.doubleToLongBits(other.evalue))
            return false;
        if (numGap != other.numGap)
            return false;
        if (numIdentical != other.numIdentical)
            return false;
        if (positives != other.positives)
            return false;
        if (queryDef == null) {
            if (other.queryDef != null)
                return false;
        } else if (!queryDef.equals(other.queryDef))
            return false;
        if (queryIsProtein != other.queryIsProtein)
            return false;
        if (queryLen != other.queryLen)
            return false;
        if (queryLoc == null) {
            if (other.queryLoc != null)
                return false;
        } else if (!queryLoc.equals(other.queryLoc))
            return false;
        if (querySeq == null) {
            if (other.querySeq != null)
                return false;
        } else if (!querySeq.equals(other.querySeq))
            return false;
        if (subjectDef == null) {
            if (other.subjectDef != null)
                return false;
        } else if (!subjectDef.equals(other.subjectDef))
            return false;
        if (subjectIsProtein != other.subjectIsProtein)
            return false;
        if (subjectLen != other.subjectLen)
            return false;
        if (subjectLoc == null) {
            if (other.subjectLoc != null)
                return false;
        } else if (!subjectLoc.equals(other.subjectLoc))
            return false;
        if (subjectSeq == null) {
            if (other.subjectSeq != null)
                return false;
        } else if (!subjectSeq.equals(other.subjectSeq))
            return false;
        return true;
    }

}
