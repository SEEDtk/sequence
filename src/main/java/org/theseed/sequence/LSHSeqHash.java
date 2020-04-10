/**
 *
 */
package org.theseed.sequence;

import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

import org.theseed.counters.QualityCountMap;
import org.theseed.sequence.Bucket.Result;

/**
 * This is a locality-sensitive hash that maps sequences to strings using a Mash/MinHash approach.  The sequence
 * is converted into a sketch using the Mash algorithm.  Each sketch is converted into an array of bucket numbers.
 * If the sequence is being inserted, it is placed into the matching bucket at each stage.  If the sequence is being
 * used for a search, we examine the matching bucket at each stage and do direct comparisons to find the closest
 * sequences.
 *
 * @author Bruce Parrello
 *
 */
public class LSHSeqHash  {

    /** large prime number for hashing */
    protected static final long LARGE_PRIME =  433494437;

    // FIELDS

    /** number of stages */
    private int stages;

    /** number of objects stored */
    private int size;

    /** width of each signature */
    private int width;

    /** number of buckets per stage */
    private int buckets;

    /** master hash table (first index is stage, second is bucket */
    private Bucket[][] masterTable;

    /**
     * Construct a blank, empty sequence hash.
     *
     * @param w		width of signatures
     * @param s		number of stages
     * @param b		number of buckets
     */
    public LSHSeqHash(int w, int s, int b) {
        this.stages = s;
        this.size = 0;
        this.masterTable = new Bucket[s][b];
        for (int si = 0; si < s; si++)
            for (int bi = 0; bi < b; bi++)
                this.masterTable[si][bi] = new Bucket();
        this.width = w;
        this.buckets = b;
    }

    /**
     * @return the distance between two signatures
     *
     * @param signature		first signature
     * @param otherSig		second signature
     */
    protected double sketchDistance(int[] signature, int[] otherSig) {
        return SequenceKmers.signatureDistance(signature, otherSig);
    }

    /**
     * @return the estimated distance between two sequences using sketches
     *
     * @param kmers		kmers for the first sequence
     * @param other		kmers for the second sequence
     */
    public double testSketches(SequenceKmers kmers, SequenceKmers other) {
        int[] signature = kmers.hashSet(this.width);
        int[] otherSig = other.hashSet(this.width);
        return sketchDistance(signature, otherSig);
    }

    /**
     * Add a sequence to the hash.
     *
     * @param seqKmers	kmer object for the sequence
     * @param target	string representing the target for the sequence
     */
    public void add(SequenceKmers seqKmers, String target) {
        // Compute the sketch for the new sequence.
        int[] signature = seqKmers.hashSet(this.width);
        Sketch entry = new Sketch(signature, target);
        // Add it to the hash.
        this.add(entry);
    }

    /**
     * Add a sketch to the hash.
     *
     * @param sketch	sketch containing a sequence signature and name
     */
    public void add(Sketch sketch) {
        // Calculate the buckets.
        int[] buckets = this.hashSignature(sketch.getSignature());
        // Place this target in the appropriate buckets.
        for (int i = 0; i < this.stages; i++)
            this.masterTable[i][buckets[i]].add(sketch);
        // Record the new entry.
        this.size++;
    }

    /**
     * Add a sketch to the hash.
     */

    /**
     * Hash a signature.  The signature is divided in s stages (or bands). Each stage is hashed to
     * one of the b buckets.
     *
     * @param signature		signature to hash
     *
     * @return a vector of s integers (between 0 and b-1)
     */
    public int[] hashSignature(final int[] signature) {
        // Create an accumulator for each stage
        int[] hash = new int[stages];
        // Number of signature elements per stage
        int rows = Math.max(1, signature.length / stages);
        for (int i = 0; i < signature.length; i++) {
            int stage = Math.min(i / rows, stages - 1);
            hash[stage] = (int)
                    ((hash[stage] + (long) signature[i] * LARGE_PRIME)
                    % buckets);
        }
        return hash;
    }


    /**
     * Search for the N closest sequences to an input sequence
     *
     * @param seqKmers	the search sequence of interest
     * @param n			the maximum number of results to return
     * @param maxDist	the maximum allowable distance
     *
     * @return the N closest sequence results in a sorted set
     */
    public SortedSet<Bucket.Result> getClosest(SequenceKmers seqKmers, int n,
            double maxDist) {
        // Compute the signature for the search sequence.
        int[] signature = seqKmers.hashSet(this.width);
        SortedSet<Bucket.Result> retVal = search(n, maxDist, signature);
        return retVal;
    }

    /**
     * Search for the N closest sequences to an input signature.
     *
     * @param n				maximum number of sequences to return
     * @param maxDist		maximum acceptable signature distance
     * @param signature		signature to find
     *
     * @return a sorted set of results found
     */
    private SortedSet<Bucket.Result> search(int n, double maxDist, int[] signature) {
        // This will be the return set.
        SortedSet<Bucket.Result> retVal = new TreeSet<Bucket.Result>();
        // Get the list of buckets to search.
        int[] buckets = this.hashSignature(signature);
        // Check each bucket for close sequences.
        for (int i = 0; i < this.stages; i++) {
            Bucket bucket = this.masterTable[i][buckets[i]];
            bucket.search(retVal, n, maxDist, signature);
        }
        return retVal;
    }

    /**
     * @return the number of entries in this table
     */
    public int size() {
        return this.size;
    }

    /**
     * Compute the quality of this hash.  In order to use this function, it is presumed that every
     * name (sequence ID) represents a cluster.  The quality of a cluster is the fraction of its
     * sequences that appear in at least one stage with another sequence in the same cluster.
     *
     * @return a QualityCountMap containing the number of good and bad sequences for each named cluster
     */
    public QualityCountMap<String> getQualityData() {
        QualityCountMap<String> retVal = new QualityCountMap<String>();
        // Each incoming sequence appears as a sketch, and a copy of that sketch will appear
        // once in every stage.  We go through the buckets in the first stage.  For each,
        // we compute its buckets in all the stages using hashSignature.  Then we check
        // each bucket and stop if we find a different sketch for the same cluster.
        for (Bucket bucket : this.masterTable[0]) {
            for (Sketch sketch : bucket) {
                int[] bucketList = this.hashSignature(sketch.getSignature());
                String name = sketch.getName();
                boolean matchFound = false;
                // Check each bucket for a match.
                for (int i = 0; i < this.stages && ! matchFound; i++) {
                    Bucket testBucket = this.masterTable[i][bucketList[i]];
                    Iterator<Sketch> iter = testBucket.iterator();
                    while (iter.hasNext() && ! matchFound) {
                        Sketch other = iter.next();
                        matchFound = (other != sketch && name.contentEquals(other.getName()));
                    }
                }
                if (matchFound)
                    retVal.setGood(name);
                else
                    retVal.setBad(name);
            }
        }
        return retVal;
    }

    /**
     * @return the quality ratio (good/total) for this hash
     */
    public double getQuality() {
        double retVal = 1.0;
        int good = 0;
        int bad = 0;
        QualityCountMap<String> counts = this.getQualityData();
        for (String cluster : counts.allKeys()) {
            good += counts.good(cluster);
            bad += counts.bad(cluster);
        }
        if (good > 0)
            retVal = good / (double) (good + bad);
        else if (bad > 0)
            retVal = 0.0;
        return retVal;
    }

    /**
     * @return a set of sketches close to an incoming sketch
     *
     * @param sketch	sketch whose neighbors are desired
     * @param maxDist	the maximum acceptable sketch distance
     */
    public SortedSet<Result> getClose(Sketch sketch, double maxDist) {
        return this.search(Integer.MAX_VALUE, maxDist, sketch.getSignature());
    }


}
