/**
 *
 */
package org.theseed.sequence.blast;

import java.util.ArrayList;
import java.util.Arrays;
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
     * Class representing data on a sequence.  There is one of these for the query
     * and one for the subject.
     */
    public static class SeqData {
        /** sequence definition (comment) */
        private String seqDef;
        /** total sequence length */
        private int seqLen;
        /** location of matching sequence */
        private Location seqLoc;
        /** matching sequence in alignment form */
        private String seqAlignment;
        /** TRUE if this is a protein sequence */
        private boolean isProt;

        /**
         * Construct a sequence data object.
         *
         * @param def			definition (comment)
         * @param len			total length
         * @param loc			location of the matching part
         * @param alignment		alignment string
         * @param isProt		TRUE if this sequence is a protein
         */
        public SeqData(String def, int len, Location loc, String alignment, boolean isProt) {
            this.seqDef = def;
            this.seqLen = len;
            this.seqLoc = loc;
            this.seqAlignment = alignment;
            this.isProt = isProt;
        }

        /**
         * @return the sequence ID
         */
        public String getId() {
            return this.seqLoc.getContigId();
        }

        /**
         * @return the sequence definition (comment)
         */
        public String getDef() {
            return seqDef;
        }

        /**
         * @return the total sequence length
         */
        public int getLen() {
            return seqLen;
        }

        /**
         * @return the location of the matching part
         */
        public Location getLoc() {
            return seqLoc;
        }

        /**
         * @return the alignment string
         */
        public String getAlignment() {
            return seqAlignment;
        }

        /**
         * @return TRUE if this is a protein sequence
         */
        public boolean isProtein() {
            return isProt;
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((seqDef == null) ? 0 : seqDef.hashCode());
            result = prime * result + seqLen;
            result = prime * result + ((seqLoc == null) ? 0 : seqLoc.hashCode());
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
            SeqData other = (SeqData) obj;
            if (seqDef == null) {
                if (other.seqDef != null)
                    return false;
            } else if (!seqDef.equals(other.seqDef))
                return false;
            if (seqLen != other.seqLen)
                return false;
            if (seqLoc == null) {
                if (other.seqLoc != null)
                    return false;
            } else if (!seqLoc.equals(other.seqLoc))
                return false;
            return true;
        }

    }

    /**
     * Comparator for comparing blast hits by location in a target sequence
     */
    public static class ByLoc implements Comparator<BlastHit> {

        /** type to sort on */
        int type;
        /** other type */
        int other;

        /**
         * Construct a comparator to sort by the location in a target sequence,
         * and then by the ID of the other sequence.
         *
         * @param t		type whose location is to be sorted on (query or subject)
         */
        public ByLoc(int t) {
            this.type = t;
            this.other = 1 - t;
        }

        @Override
        public int compare(BlastHit o1, BlastHit o2) {
            int retVal = o1.seqs[type].getLoc().compareTo(o2.seqs[type].getLoc());
            if (retVal == 0) {
                retVal = o1.seqs[other].getId().compareTo(o2.seqs[other].getId());
                if (retVal == 0) {
                    retVal = o2.getNumSimilar() - o1.getNumSimilar();
                    if (retVal == 0) {
                        retVal = Double.compare(o1.getEvalue(), o2.getEvalue());
                        if (retVal == 0) {
                            retVal = o2.getNumIdentical() - o1.getNumIdentical();
                            if (retVal == 0) {
                                retVal = o1.seqs[other].getLoc().compareTo(o2.seqs[other].getLoc());
                            }
                        }
                    }
                }
            }
            return retVal;
        }

    }

    // FIELDS
    /** index of the query sequence data */
    public static final int QUERY = 0;
    /** index of the subject sequence data */
    public static final int SUBJECT = 1;
    /** sequence data for both sequences */
    private SeqData[] seqs;
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
        // Split the input line into fields.
        String[] fields = StringUtils.split(line, '\t');
        // Get the query ID and parse the subject title into ID and comment.
        String qid = fields[0];
        String[] pieces = StringUtils.split(fields[5], " ", 2);
        // Build the sequence data objects.
        this.seqs = new SeqData[] {
                new SeqData(qMap.getOrDefault(qid, ""), Integer.valueOf(fields[1]),
                        Location.create(qid, Integer.valueOf(fields[2]), Integer.valueOf(fields[3])),
                        fields[4], queryIsProtein),
                new SeqData(pieces[1], Integer.valueOf(fields[6]),
                        Location.create(pieces[0], Integer.valueOf(fields[7]), Integer.valueOf(fields[8])),
                        fields[9], subjectIsProtein)
        };
        this.bitScore = Double.valueOf(fields[10]);
        this.numIdentical = Integer.valueOf(fields[11]);
        this.numGap = Integer.valueOf(fields[12]);
        this.evalue = Double.valueOf(fields[13]);
        this.alignLen = Integer.valueOf(fields[15]);
        this.positives = Integer.valueOf(fields[14]) - this.numIdentical;
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
     * @return the sequence data for the sequence of the specified type
     *
     * @param type	type of data-- QUERY or SUBJECT
     */
    public SeqData getData(int type) {
        return seqs[type];
    }

    /**
     * @return the query sequence description (comment)
     */
    public String getQueryDef() {
        return seqs[QUERY].getDef();
    }

    /**
     * @return the total length of the query sequence
     */
    public int getQueryLen() {
        return seqs[QUERY].getLen();
    }

    /**
     * @return the location of the matching portion of the query sequence
     */
    public Location getQueryLoc() {
        return seqs[QUERY].getLoc();
    }

    /**
     * @return the matching portion of the query sequence, in alignment format
     */
    public String getQuerySeq() {
        return seqs[QUERY].getAlignment();
    }

    /**
     * @return the subject sequence description (comment)
     */
    public String getSubjectDef() {
        return seqs[SUBJECT].getDef();
    }

    /**
     * @return the total length of the subject sequence
     */
    public int getSubjectLen() {
        return seqs[SUBJECT].getLen();
    }

    /**
     * @return the location of the matching portion of the subject sequence
     */
    public Location getSubjectLoc() {
        return seqs[SUBJECT].getLoc();
    }

    /**
     * @return the matching portion of the subject sequence, in alignment format
     */
    public String getSubjectSeq() {
        return seqs[SUBJECT].getAlignment();
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
     * @return the match length as a percent of the specified sequence length
     *
     * @param type	denominator sequence type-- QUERY or SUBJECT
     */
    public double getPercentMatch(int type) {
        double retVal = 0.0;
        if (this.seqs[type].getLen() > 0) {
            double dLen = (double) this.seqs[type].getLen();
            if (! seqs[type].isProtein() && seqs[1-type].isProtein())
                dLen /= 3.0;
            retVal = (this.getPositives() + this.getNumIdentical()) * 100 / dLen;
        }
        return retVal;
    }

    /**
     * @return the match length as a percent of the query length
     */
    public double getQueryPercentMatch() {
        return this.getPercentMatch(QUERY);
    }

    /**
     * @return the match length as a percent of the subject length
     */
    public double getSubjectPercentMatch() {
        return this.getPercentMatch(SUBJECT);
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
        return this.seqs[QUERY].getLoc().getContigId();
    }

    /**
     * @return the subject sequence ID
     */
    public String getSubjectId() {
        return this.seqs[SUBJECT].getLoc().getContigId();
    }

    /**
     * @return a tab-delimited print line for this result (with no LF)
     */
    public String getPrintLine() {
        return String.format("%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t" +
                "%4.2f\t%d\t%d\t%d\t%6.3e\t%d\t%s\t%s",
                this.getQueryId(), this.getQueryDef(), this.getQueryLoc().toString(),
                this.getQueryLen(), this.getSubjectId(), this.getSubjectDef(),
                this.getSubjectLoc().toString(), this.getSubjectLen(), this.getBitScore(),
                this.getNumIdentical(), this.getNumGap(), this.getPositives(),
                this.getEvalue(), this.getAlignLen(), this.getQuerySeq(),
                this.getSubjectSeq());
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
        result = prime * result + Arrays.hashCode(seqs);
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
        if (!Arrays.equals(seqs, other.seqs))
            return false;
        return true;
    }


}
