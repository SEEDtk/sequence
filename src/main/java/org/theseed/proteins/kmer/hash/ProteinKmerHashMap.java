/**
 *
 */
package org.theseed.proteins.kmer.hash;

import java.io.UnsupportedEncodingException;
import java.security.NoSuchAlgorithmException;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.sequence.MD5Hex;
import org.theseed.sequence.ProteinKmers;

/**
 * This object is a hash map for proteins using a kmer lookup.  It's goal is to find a close protein
 * rather than an exact match.  The basic map is keyed by protein MD5, but there is also a secondary
 * map.  This map is keyed by kmer, and each kmer maps to a set of MD5s.  When a new protein comes in,
 * its kmers are all hashed, and the MD5 with the most hits is the closest.  The jaccard distance and
 * the associated value are both returned.
 */
public class ProteinKmerHashMap<T> {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ProteinKmerHashMap.class);
    /** protein MD5 computer */
    private MD5Hex md5Engine;
    /** MD5 map */
    private HashMap<String, Md5Entry> md5Map;
    /** kmer hash */
    private HashMap<String, Set<String>> kmerMap;
    /** kmer size */
    private int kmerSize;
    /** result returned when no protein is close */
    public final Result NOT_FOUND = new Result();

    /**
     * This object contains the kmer count and associated value for a protein.
     */
    public class Md5Entry {

        /** kmer count */
        private int kmerCount;
        /** associated value */
        private T value;

        /**
         * Construct a new MD5 entry.
         *
         * @param count		kmer count
         * @param val		associated value
         */
        protected Md5Entry(int count, T val) {
            this.kmerCount = count;
            this.value = val;
        }

        /**
         * @return the number of kmers in the protein
         */
        public int getKmerCount() {
            return this.kmerCount;
        }

        /**
         * @return the value associated with the protein
         */
        public T getValue() {
            return this.value;
        }

    }

    /**
     * This object represents a result from a close-protein query.
     */
    public class Result {

        /** MD5 of closest protein */
        private String md5;
        /** number of kmers in common */
        private int simCount;
        /** jaccard similarity to closest protein */
        private double simValue;
        /** associated value */
        private T value;

        /**
         * Construct a blank result.
         */
        protected Result() {
            this.md5 = "";
            this.simCount = 0;
            this.simValue = 0.0;
            this.value = null;
        }

        /**
         * Construct a result for a close protein.
         *
         * @param hitCount	number of kmers in common
         * @param origKmers	number of kmers in trial protein
         * @param md5		closest MD5
         */
        protected Result(int hitCount, int origKmers, String md5) {
            // Save the MD5.
            this.md5 = md5;
            // Get the MD5 entry from the main map and compute the scores.
            Md5Entry entry = ProteinKmerHashMap.this.md5Map.get(md5);
            if (entry == null)
                throw new IllegalArgumentException("\"" + md5 + "\" not found in protein kmer map.");
            this.simCount = hitCount;
            this.simValue = hitCount / (double) (origKmers + entry.getKmerCount() - hitCount);
            // Save the associated value.
            this.value = entry.getValue();
        }

        /**
         * @return the md5 for the closest protein
         */
        public String getMd5() {
            return this.md5;
        }

        /**
         * @return the number of kmers in common between the proteins
         */
        public int getSimCount() {
            return this.simCount;
        }

        /**
         * @return the jaccard similarity score between the proteins
         */
        public double getSimValue() {
            return this.simValue;
        }

        /**
         * @return the value associated with the closest protein
         */
        public T getValue() {
            return this.value;
        }

        /**
         * @return TRUE if this is a not-found result, else FALSE
         */
        public boolean isEmpty() {
            return this.simCount == 0;
        }

    }

    /**
     * Construct a blank, empty protein kmer hash.
     *
     * @param K		kmer size
     */
    public ProteinKmerHashMap(int K) {
        this.kmerSize = K;
        this.kmerMap = new HashMap<>();
        this.md5Map = new HashMap<>();
        try {
            this.md5Engine = new MD5Hex();
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException("Missing MD5 algorithm on this installation: " + e.getMessage());
        }
    }

    /**
     * @return the kmer size
     */
    public int getKmerSize() {
        return this.kmerSize;
    }

    /**
     * Add a protein to the kmer hash.
     *
     * @param prot		protein to add
     * @param value		associated value
     */
    public void addProtein(String prot, T value) {
        String protUC = StringUtils.upperCase(prot);
        // Get the MD5 for this protein.
        String md5;
        try {
            md5 = this.md5Engine.checksum(protUC);
        } catch (UnsupportedEncodingException e) {
            throw new IllegalArgumentException("Invalid characters in protein string.");
        }
        // Compute the protein kmers.
        ProteinKmers kmers = new ProteinKmers(protUC, this.kmerSize);
        // Put this protein into the hash.
        Md5Entry entry = new Md5Entry(kmers.size(), value);
        this.md5Map.put(md5, entry);
        // Now loop through the kmers, putting this MD5 into each kmer's protein set.
        for (String kmer : kmers) {
            Set<String> md5Set = this.kmerMap.computeIfAbsent(kmer, x -> new TreeSet<String>());
            md5Set.add(md5);
        }
    }

    /**
     * Find the closest result to a specified protein.
     *
     * @param prot		protein to search for
     *
     * @return a result containing stats about the protein found and its associated value
     */
    public Result findClosest(String prot) {
        String protUC = StringUtils.upperCase(prot);
        // TODO use the hash
        return null;
    }

}
