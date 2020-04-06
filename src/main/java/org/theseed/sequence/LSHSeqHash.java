/**
 *
 */
package org.theseed.sequence;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import info.debatty.java.lsh.LSH;
import info.debatty.java.lsh.MinHash;

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
public class LSHSeqHash extends LSH {

    /**
     * serialization class ID
     */
    private static final long serialVersionUID = -4410879652283658775L;

    /**
     * This subclass contains a sequence key and its connected object.
     */
    private class Entry {
        /** sketch for this entry's sequence */
        private int[] sketch;
        /** target object associated with the sequence */
        private String name;

        /**
         * Construct a hash entry.
         *
         * @param seq		kmer object to use for comparison
         * @param string	string associated with the sequence that generated the kmers
         */
        public Entry(int[] sk, String string) {
            this.sketch = sk;
            this.name = string;
        }

        /**
         * @return the distance between another sequence and this entry's sequence
         *
         * @param otherSketch	sketch of the other sequence
         */
        public double distance(int[] otherSketch) {
            return LSHSeqHash.this.sketchDistance(this.sketch, otherSketch);
        }

        /**
         * @return the sequence target
         */
        public String getName() {
            return name;
        }

    }

    /**
     * This subclass forms a bucket.
     */
    private static class Bucket implements Iterable<Entry> {
        private List<Entry> entries;

        /**
         * Create an empty bucket.
         */
        public Bucket() {
            this.entries = new LinkedList<Entry>();
        }

        /**
         * Add an entry to a bucket.
         *
         * @param entry		entry to add
         */
        public void add(Entry entry) {
            this.entries.add(entry);
        }

        @Override
        public Iterator<Entry> iterator() {
            return entries.iterator();
        }

    }

    /**
     * This subclass returns a result.  It contains the distance and the target (sequence ID),
     * and is sorted by distance followed by sequence ID.
     */
    public static class Result implements Comparable<Result> {

        private double distance;
        private String target;

        /**
         * Create a new result object.
         *
         * @param distance		distance to the search sequence
         * @param target		target (sequence ID) of the found sequence
         */
        private Result(double distance, String target) {
            this.distance = distance;
            this.target = target;
        }

        /**
         * Merge a new result into a result set.
         *
         * @param results	sorted result set to merge into
         * @param maxSize	maximum number of items allowed in result set
         * @param dist		distance from search sequence to found sequence
         * @param id		target (sequence ID) of found sequence
         *
         * @return TRUE if the new result is accepted, else FALSE
         */
        public static boolean merge(SortedSet<Result> results, int maxSize, double dist, String id) {
            boolean retVal = false;
            Result next = new Result(dist, id);
            if (results.size() < maxSize) {
                results.add(next);
                retVal = true;
            } else {
                Result last = results.last();
                if (last.compareTo(next) > 0) {
                    results.remove(last);
                    results.add(next);
                    retVal = true;
                }
            }
            return retVal;
        }

        @Override
        public int compareTo(Result o) {
            int retVal = Double.compare(this.distance, o.distance);
            if (retVal == 0)
                retVal = this.target.compareTo(o.target);
            return retVal;
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            long temp;
            temp = Double.doubleToLongBits(distance);
            result = prime * result + (int) (temp ^ (temp >>> 32));
            result = prime * result + ((target == null) ? 0 : target.hashCode());
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
            Result other = (Result) obj;
            if (Double.doubleToLongBits(distance) != Double.doubleToLongBits(other.distance))
                return false;
            if (target == null) {
                if (other.target != null)
                    return false;
            } else if (!target.equals(other.target))
                return false;
            return true;
        }

        /**
         * @return the distance from the search sequence to this result
         */
        public double getDistance() {
            return distance;
        }

        /**
         * @return the target (sequence ID) of this result
         */
        public String getTarget() {
            return target;
        }

    }

    // FIELDS

    /** number of stages */
    private int stages;

    /** number of objects stored */
    private int size;

    /** master hash table (first index is stage, second is bucket */
    private Bucket[][] masterTable;

    /** minhash object for creating sketches */
    private MinHash minHasher;

    /**
     * Construct a blank, empty sequence hash.
     *
     * @param w		width of signatures
     * @param s		number of stages
     * @param b		number of buckets
     * @param seed	random number seed
     */
    public LSHSeqHash(int w, int s, int b, long seed) {
        super(s, b);
        this.stages = s;
        this.size = 0;
        this.masterTable = new Bucket[s][b];
        for (int si = 0; si < s; si++)
            for (int bi = 0; bi < b; bi++)
                this.masterTable[si][bi] = new Bucket();
        this.minHasher = new MinHash(w, SequenceKmers.DICT_SIZE, seed);
    }

    /**
     * @return a random number for use as a seed.
     */
    public static long randomSeed() {
        return (System.currentTimeMillis() % LSH.LARGE_PRIME) * 13 + 1;
    }

    /**
     * @return the distance between two sketches
     *
     * @param sketch		first sketch
     * @param otherSketch	second sketch
     */
    protected double sketchDistance(int[] sketch, int[] otherSketch) {
        return 1.0 - minHasher.similarity(sketch, otherSketch);
    }

    /**
     * @return the estimated distance between two sequences using sketches
     *
     * @param kmers		kmers for the first sequence
     * @param other		kmers for the second sequence
     */
    public double testSketches(SequenceKmers kmers, SequenceKmers other) {
        int[] sketch = this.minHasher.signature(kmers.hashSet());
        int[] otherSketch = this.minHasher.signature(other.hashSet());
        return sketchDistance(sketch, otherSketch);
    }

    /**
     * Add a sequence to the hash.
     *
     * @param seqKmers	kmer object for the sequence
     * @param target	string representing the target for the sequence
     */
    public void add(SequenceKmers seqKmers, String target) {
        // Compute the sketch and the bucket array for the new sequence.
        int[] sketch = this.minHasher.signature(seqKmers.hashSet());
        int[] buckets = this.hashSignature(sketch);
        // Place this target in the appropriate buckets.
        Entry entry = new Entry(sketch, target);
        for (int i = 0; i < this.stages; i++)
            this.masterTable[i][buckets[i]].add(entry);
        // Record the new entry.
        this.size++;
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
    public SortedSet<Result> getClosest(SequenceKmers seqKmers, int n,
            double maxDist) {
        SortedSet<Result> retVal = new TreeSet<Result>();
        // Compute the sketch and the bucket array for the search sequence.
        int[] sketch = this.minHasher.signature(seqKmers.hashSet());
        int[] buckets = this.hashSignature(sketch);
        // Check each bucket for close sequences.
        for (int i = 0; i < this.stages; i++) {
            for (Entry entry : this.masterTable[i][buckets[i]]) {
                double dist = entry.distance(sketch);
                if (dist <= maxDist)
                    Result.merge(retVal, n, dist, entry.getName());
            }
        }
        return retVal;
    }

    /**
     * @return the number of entries in this table
     */
    public int size() {
        return this.size;
    }

}
