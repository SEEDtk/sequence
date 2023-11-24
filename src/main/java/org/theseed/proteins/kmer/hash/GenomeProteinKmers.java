/**
 *
 */
package org.theseed.proteins.kmer.hash;

import java.util.List;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This object tracks the protein kmers in a genome.  For each kmer, we list the MD5s that contained it.
 * Each MD5 is then associated with a count that contains the best annotation.
 *
 * @author Bruce Parrello
 *
 */
public class GenomeProteinKmers {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(GenomeProteinKmers.class);
    /** protein kmer hash map for this genome, maps each MD5 to its best annotation (count and string) */
    private ProteinKmerHashMap<Proposal> protMap;
    /** minimum similarity threshold */
    private double minSim;
    /** default annotation */
    private static final String DEFAULT_ANNOTATION = "";

    /**
     * This class contains a similarity score and a proposed annotation.
     */
    public static class Proposal {

        /** similarity score */
        private double sim;
        /** proposed annotation */
        private String annotation;

        /**
         * Construct an empty proposal.
         */
        protected Proposal() {
            this.annotation = DEFAULT_ANNOTATION;
            this.sim = 0.0;
        }

        /**
         * Construct a proposal with a specified default annotation.
         *
         * @param anno		default annotation to use
         */
        public Proposal(String anno) {
            this.annotation = anno;
            this.sim = 0.0;
        }

        /**
         * Merge a new proposal into this one.
         *
         * @param score		score for the new proposal
         * @param anno		annotation for the new proposal
         */
        public void merge(double score, String anno) {
            if (this.sim < score) {
                this.sim = score;
                this.annotation = anno;
            }
        }

        /**
         * @return the similarity score
         */
        public double getSim() {
            return this.sim;
        }

        /**
         * @return the proposed annotation
         */
        public String getAnnotation() {
            return this.annotation;
        }

    }

    /**
     * Construct a blank, empty genome protein kmer structure.
     *
     * @param K		kmer size
     * @param min	minimum similarity score
     */
    public GenomeProteinKmers(int K, double min) {
        this.protMap = new ProteinKmerHashMap<Proposal>(K);
        this.minSim = min;
    }

    /**
     * Add a protein to the hash map.
     *
     * @param fid		feature ID
     * @param prot		protein sequence to add
     * @param anno		default annotation
     *
     * @return the protein MD5
     */
    public String addProtein(String fid, String prot, String anno) {
        String retVal = this.protMap.addProtein(prot, new Proposal(anno));
        return retVal;
    }

    /**
     * Process a proposed protein annotation.  The source protein sequence is matched against the genome, and
     * each protein to which it is sufficiently close is checked.  If the new sequence is more similar, its
     * annotation is kept.
     *
     * @param protein		source protein sequence
     * @param annotation	annotation for this sequence
     *
     * @return the number of close proteins found
     */
    public int processProposal(String protein, String annotation) {
        // Find the genome proteins close to this one.
        List<ProteinKmerHashMap<Proposal>.Result> closeProteins = this.protMap.findClose(protein, this.minSim);
        for (var closeProtein : closeProteins) {
            // Get the proposal for the close protein in the genome.
            Proposal protProposal = closeProtein.getValue();
            // Keep this annotation if the source protein is closer.
            protProposal.merge(closeProtein.getSimValue(), annotation);
        }
        return closeProteins.size();
    }

    /**
     * @return the number of kmers in the hash
     */
    public int getKmerCount() {
        return this.protMap.getKmerCount();
    }

    /**
     * @return the proposal for a protein based on its MD5, or NULL if none
     *
     * @param md5	MD5 for the protein in question
     */
    public Proposal getProposal(String md5) {
        return this.protMap.getByMd5(md5);
    }
}
