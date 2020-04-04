/**
 *
 */
package org.theseed.sequence;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

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
public class SeqHash extends LSH {

    /**
     * This subclass contains a sequence key and its connected object.
     */
    private class Entry {
        /** sequence for this entry */
        private SequenceKmers sequence;
        /** target object associated with the sequence */
        private String name;

        /**
         * Construct a hash entry.
         *
         * @param seq		kmer object to use for comparison
         * @param string	string associated with the sequence that generated the kmers
         */
        public Entry(SequenceKmers seq, String string) {
            this.sequence = seq;
            this.name = string;
        }

        /**
         * @return the distance between another sequence and this entry's sequence
         */
        public double distance(SequenceKmers other) {
            return sequence.distance(other);
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
    private class Bucket implements Iterable<Entry> {
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
     * @param s		number of stages
     * @param b		number of buckets
     * @param seed	random number seed
     */
    public SeqHash(int s, int b, long seed) {
        super(s, b);
        this.stages = s;
        this.size = 0;
        this.masterTable = new Bucket[s][b];
        this.minHasher = new MinHash(s, SequenceKmers.DICT_SIZE, seed);
    }


    //TODO use mh.signature(kmerSet.hashSet()) to get the sequence signature
    //TODO use this.hashSignature(sequence signature) to compute buckets per stage

}
