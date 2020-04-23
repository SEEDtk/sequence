/**
 *
 */
package org.theseed.sequence.hash;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * This class contains a list of sketches, forming a bucket.  The bucket
 * is searched to find sequences close to a search sequence.
 *
 * @author Bruce Parrello
 *
 */
public class Bucket implements Iterable<Sketch>, Serializable {

    // FIELDS
    /** serialization identifier */
    private static final long serialVersionUID = 1890888260042916972L;
    /** list of sketches in this bucket */
    private List<Sketch> entries;
    /** TRUE if this bucket has been modified since the last load */
    private transient boolean modified;
    /** time of last access */
    private transient long timeStamp;
    /** activity counter */
    private static transient long timeCounter = 0;



    /**
     * This nested class returns a result.  It contains the distance and the target (sequence ID),
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

    /**
     * Create an empty bucket.
     */
    public Bucket() {
        this.entries = new ArrayList<Sketch>(100);
        this.modified = false;
        this.timeStamp = 0;
    }

    /**
     * Add an entry to a bucket.
     *
     * @param entry		entry to add
     */
    public void add(Sketch entry) {
        this.entries.add(entry);
        this.modified = true;
        this.timeStamp = ++timeCounter;
    }

    @Override
    public Iterator<Sketch> iterator() {
        return entries.iterator();
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

    /**
     * Search for close sequences in this bucket and merge them into a result set.
     *
     * @param results	result set for output
     * @param n			maximum result set size
     * @param maxDist	maximum acceptable sketch distance
     * @param signature	signature to search for
     */
    public void search(SortedSet<Result> results, int n, double maxDist, int[] signature) {
        for (Sketch entry : this) {
            double dist = entry.distance(signature);
            if (dist <= maxDist)
                merge(results, n, dist, entry.getName());
        }
        this.timeStamp = ++timeCounter;
    }

    /**
     * Search for close sketches in this bucket.
     *
     * @param n			maximum result set size
     * @param maxDist	maximum acceptable sketch distance
     * @param signature	sketch (signature) to search for
     *
     * @return a sorted result set containing the closest sequences given the parameters
     */
    public SortedSet<Result> search(int n, double maxDist, int[] signature) {
        SortedSet<Result> retVal = new TreeSet<Result>();
        this.search(retVal, n, maxDist, signature);
        this.timeStamp = ++timeCounter;
        return retVal;
    }

    /**
     * Search for sketches with a given name.
     *
     * @param name		name (target) to search for
     *
     * @return a set of the sketches with the same name
     */
    public List<Sketch> search(String name) {
        List<Sketch> retVal = new ArrayList<Sketch>();
        for (Sketch sketch : this)
            if (name.contentEquals(sketch.getName()))
                retVal.add(sketch);
        this.timeStamp = ++timeCounter;
        return retVal;
    }

    /**
     * @return the number of sketches in this bucket
     */
    public int size() {
        return this.entries.size();
    }

    /**
     * Save this bucket to a file.
     *
     * @param outFile	file to which the bucket should be saved
     *
     * @throws IOException
     */
    public void save(File outFile) throws IOException {
        try (FileOutputStream outStream = new FileOutputStream(outFile)) {
            ObjectOutputStream objStream = new ObjectOutputStream(outStream);
            objStream.writeObject(this);
            this.modified = false;
        }
    }

    /**
     * Load a bucket from a file.
     *
     * @param inFile	file from which the bucket should be loaded
     *
     * @return the bucket loaded from the file
     *
     * @throws IOException
     * @throws ClassNotFoundException
     */
    public static Bucket load(File inFile) throws IOException, ClassNotFoundException {
        Bucket retVal = null;
        try (FileInputStream inStream = new FileInputStream(inFile)) {
            ObjectInputStream objStream = new ObjectInputStream(inStream);
            retVal = (Bucket) objStream.readObject();
            // The modification flag is transient, so we need to set it here.
            retVal.modified = false;
        }
        return retVal;
    }

    /**
     * @return the sketch at the specified index in the list
     *
     * @param i		index of the desired sketch
     */
    public Sketch get(int i) {
        return this.entries.get(i);
    }

    /**
     * @return the list of sketches after the specified position
     *
     * @param i		index of the postion after which to begin the return list
     */
    public List<Sketch> after(int i) {
        return this.entries.subList(i + 1, this.size());
    }

    /**
     * @return TRUE if this bucket is modified
     */
    public boolean isModified() {
        return this.modified;
    }

    /**
     * This is a comparison operation used to sort buckets for deletion
     * from the cache.  Buckets are compared by timestamp.  We want to sort the
     * oldest buckets first, so we go for the lowest time stamp.  If two buckets
     * are the same age, then they are usually unmodified.  In this case, we
     * prefer the shortest bucket.  Some buckets may compare equal, but that is
     * ok for us.
     */
    public int compareTo(Bucket o) {
        int retVal = Long.compare(this.timeStamp, o.timeStamp);
        if (retVal == 0) {
            retVal = this.size() - o.size();
        }
        return retVal;
    }

}
