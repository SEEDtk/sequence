/**
 *
 */
package org.theseed.proteins.kmer.hash;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.security.NoSuchAlgorithmException;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.counters.CountMap;
import org.theseed.io.TabbedLineReader;
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
    private Map<String, Md5Entry> md5Map;
    /** kmer hash */
    private Map<String, Set<String>> kmerMap;
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
        this.setup(K);
        this.kmerMap = new ConcurrentHashMap<>();
        this.md5Map = new ConcurrentHashMap<>();
    }

    /**
     * Initialize the support fields.
     *
     * @param K		kmer size
     */
    protected void setup(int K) {
        this.kmerSize = K;
        try {
            this.md5Engine = new MD5Hex();
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException("Missing MD5 algorithm on this installation: " + e.getMessage());
        }
    }

    /**
     * Construct a blank, empty protein kmer hash using a size estimate.
     *
     * @param K		kmer size
     * @param est	estimated number of proteins.
     */
    public ProteinKmerHashMap(int K, int est) {
        this.setup(K);
        int hashSize = est * 4 / 3 + 1;
        int fullSize = hashSize * K;
        this.md5Map = new ConcurrentHashMap<>(hashSize);
        this.kmerMap = new ConcurrentHashMap<>(fullSize);
    }

    /**
     * Load a protein kmer hash from a tab-delimited file.
     *
     * @param inFile	input file name
     * @param K			kmer size
     * @param pCol		index (1-based) or name of the protein column
     * @param vCol		index (1-based) or name of the value column
     *
     * @return a protein kmer hash with string values loaded from the file
     *
     * @throws IOException
     */
    public static ProteinKmerHashMap<String> load(File inFile, int K, String pCol, String vCol)
            throws IOException {
        // Create the hash map.
        ProteinKmerHashMap<String> retVal = new ProteinKmerHashMap<String>(K);
        // Open the input file.
        log.info("Loading proteins from {} with kmer size {}.", inFile, K);
        try (TabbedLineReader inStream = new TabbedLineReader(inFile)) {
            // Get the input column indices.
            int pColIdx = inStream.findField(pCol);
            int vColIdx = inStream.findField(vCol);
            // Initialize the counters.
            int skipCount = 0;
            int protCount = 0;
            // Loop through the input lines.
            long lastMsg = System.currentTimeMillis();
            for (var line : inStream) {
                String protein = line.get(pColIdx);
                if (protein.length() < K)
                    skipCount++;
                else {
                    String value = line.get(vColIdx);
                    retVal.addProtein(protein, value);
                    protCount++;
                    if (log.isInfoEnabled() && lastMsg - System.currentTimeMillis() >= 10000) {
                        log.info("{} proteins loaded, {} skipped.", protCount, skipCount);
                        lastMsg = System.currentTimeMillis();
                    }
                }
            }
            log.info("{} proteins processed, {} skipped, {} proteins stored, {} kmers.",
                    protCount, skipCount, retVal.md5Map.size(), retVal.kmerMap.size());
        }
        return retVal;
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
            synchronized (md5Set) {
                md5Set.add(md5);
            }
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
        ProteinKmers kmers = new ProteinKmers(protUC, this.kmerSize);
        // This will track the hits for each MD5.
        CountMap<String> hitCounts = new CountMap<String>();
        // Loop through the kmers, counting them.
        for (String kmer : kmers) {
            Set<String> md5Set = this.kmerMap.get(kmer);
            if (md5Set != null) {
                for (String md5 : md5Set)
                    hitCounts.count(md5);
            }
        }
        Result retVal;
        // Insure we have at least one hit.
        if (hitCounts.size() == 0)
            retVal = NOT_FOUND;
        else {
            int kCount = kmers.size();
            var bestCount = hitCounts.getBestEntry();
            retVal = new Result(bestCount.getCount(), kCount, bestCount.getKey());
        }
        return retVal;
    }

    /**
     * @return the number of kmers in the map
     */
    public int getKmerCount() {
        return this.kmerMap.size();
    }

}
