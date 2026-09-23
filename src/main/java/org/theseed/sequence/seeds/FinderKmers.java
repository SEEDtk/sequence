package org.theseed.sequence.seeds;

import java.util.Map;
import java.util.TreeMap;

import org.theseed.sequence.ProteinKmers;

/**
 * This class manages kmers for a Finder-based distance computation. For each genome, it stores a hash of role IDs to protein kmers.
 * Methods are provided to compute the distance between two genomes based on their finder kmers.
 * 
 * @author Bruce Parrello
 */
public class FinderKmers {

    // FIELDS
    /** logging facility */
    private static final org.slf4j.Logger log = org.slf4j.LoggerFactory.getLogger(FinderKmers.class);
    /** ID of the genome this kmer set belongs to */
    private final String genomeId;
    /** hash of role IDs to protein kmers for this genome */
    private final Map<String, ProteinKmers> roleKmers;

    // CONSTRUCTOR
    /**
     * Construct a new FinderKmers object for a given genome. This initializes the data structures
     * but does not populate the kmer map.
     * 
     * @param genome_id    ID of the genome this kmer set belongs to
     */
    public FinderKmers(String genome_id) {
        this.genomeId = genome_id;
        this.roleKmers = new TreeMap<>();
    }

    /**
     * @return the ID of the genome this kmer set belongs to
     */
    public String getGenomeId() {
        return this.genomeId;
    }

    /**
     * Add the specified protein sequence to the kmer map. On very rare occasions, we will get two
     * sequences for the same role, in which case either the original protein has split into two pieces,
     * or we have two very similar proteins for the same role. In either case, we merge the kmers.
     *
     * @param roleId        ID of the role the sequence belongs to
     * @param proteinSeq    protein sequence to add
     */
    public void addProteinSequence(String roleId, String proteinSeq) {
        ProteinKmers proteinKmers = new ProteinKmers(proteinSeq);
        ProteinKmers oldKmers = this.roleKmers.get(roleId);
        if (oldKmers != null)
            oldKmers.merge(proteinKmers);
        else
            this.roleKmers.put(roleId, proteinKmers);
    }

    /**
     * Compute the closeness between this genome and another genome based on their finder kmers.
     *
     * @param other    the other FinderKmers object to compare against
     * 
     * @return the similarity between this genome and the other genome, expressed as a value from 0 to 1
     */
    public double computeCloseness(FinderKmers other) {
        int commonKmers = 0;
        int totalKmers = 0;
        // First we find the common kmers between this genome and the other genome.
        // We also total the kmers for this genome.
        for (String roleId : this.roleKmers.keySet()) {
            ProteinKmers thisKmers = this.roleKmers.get(roleId);
            totalKmers += thisKmers.size();
            ProteinKmers otherKmers = other.roleKmers.get(roleId);
            if (otherKmers != null) {
                int sim = thisKmers.rawSimilarity(otherKmers);
                commonKmers += sim;
            }
        }
        // Now we add the kmers from the other genome to the total kmers.
        for (ProteinKmers otherKmers : other.roleKmers.values()) {
            totalKmers += otherKmers.size();
        }
        // Subtract the common kmers from the total to avoid double-counting.
        totalKmers -= commonKmers;
        // If the total kmers is 0, the closeness is 0. This can only happen if both of the genomes are horrendously
        // incomplete. Otherwise, the closeness is computed by the fraction of common kmers over total kmers.
        double retVal = 0.0;
        if (totalKmers > 0)
            retVal = commonKmers / (double) totalKmers;
        return retVal;
    }

    /**
     * Compute the distance between this genome and another genome based on their finder kmers.
     * This is just the inverse of the similarity.
     *
     * @param other    the other FinderKmers object to compare against
     * 
     * @return the distance between this genome and the other genome
     */
    public double computeDistance(FinderKmers other) {
        return 1.0 - this.computeCloseness(other);
    }


    /**
     * @return the number of kmers in this object
     */
    public int size() {
        int totalKmers = 0;
        for (ProteinKmers kmers : this.roleKmers.values()) {
            totalKmers += kmers.size();
        }
        return totalKmers;
    }

    /**
     * @return the number of roles in this object
     */
    public int roleCount() {
        return this.roleKmers.size();
    }

}
