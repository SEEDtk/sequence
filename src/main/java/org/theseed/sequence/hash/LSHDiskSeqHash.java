/**
 *
 */
package org.theseed.sequence.hash;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UncheckedIOException;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.io.FileUtils;

/**
 * This is a version of the locally-sensitive sequence hash that stores its data on disk.  The
 * main data item here is the bucket.  Each bucket is stored as a sequential disk file, and is
 * read in on demand.  A cache of buckets is kept in memory, with a least-recently-used
 * algorithm determining which to replace when memory is full.  The hash is stored on disk
 * as a directory, with each stage being a subdirectory and each bucket a file in that directory.
 *
 * Note that under this scheme, a single sketch may have multiple copies in memory due to the
 * caching.
 *
 * @author Bruce Parrello
 *
 */
public class LSHDiskSeqHash extends LSHSeqHash implements AutoCloseable {

    // FIELDS

    /** name for the control file */
    private static final String CONTROL_FILE_NAME = "lsh_hash.txt";

    /** name of the directory containing the bucket files */
    private File mainDir;

    /** hash of buckets in memory, by file name */
    private Map<File, Bucket> bucketCache;

    /** optimal cache size */
    private static int cacheLimit = 1000;

    /** recommended kmer size */
    private int kmerSize;

    /** comparator for bucket cache sort */
    private static EntryCompare COMPARATOR = new EntryCompare();

    protected LSHDiskSeqHash(int w, int s, int b, int K, File dir) {
        super(w, s, b);
        this.kmerSize = K;
        this.mainDir = dir;
        this.bucketCache = new HashMap<File, Bucket>(s * b);
    }

    /**
     * Create a new disk hash.
     *
     * @param w		width of each sequence sketch
     * @param s		number of stages
     * @param b		number of buckets per stage
     * @param K		kmer size for sequences
     * @param dir	directory in which to build the hash
     *
     * @throws IOException
     */
    public static LSHDiskSeqHash create(int w, int s, int b, int K, File dir) throws IOException {
        // Insure the directory exists.
        if (! dir.isDirectory()) {
            log.info("Creating hash directory {}.", dir);
            FileUtils.forceMkdir(dir);
        } else {
            log.info("Erasing hash directory {}.", dir);
            FileUtils.cleanDirectory(dir);
        }
        // Create the control file.
        File controlFile = new File(dir, CONTROL_FILE_NAME);
        try (PrintWriter writer = new PrintWriter(controlFile)) {
            writer.format("%d %d %d %d%n", w, s, b, K);
        }
        // Create the object.
        LSHDiskSeqHash retVal = new LSHDiskSeqHash(w, s, b, K, dir);
        // Create the stage directories.
        log.info("Creating {} stage directories in {}.", s, dir);
        for (int i = 0; i < s; i++) {
            FileUtils.forceMkdir(retVal.stageDir(i));
        }
        return retVal;
    }

    /**
     * Load an existing disk hash.
     *
     * @param dir	directory containing the hash
     *
     * @throws IOException
     */
    public static LSHDiskSeqHash load(File dir) throws IOException {
        // Read the control file.
        File controlFile = new File(dir, CONTROL_FILE_NAME);
        // The return value will be built in here.
        LSHDiskSeqHash retVal;
        try (Scanner scanner = new Scanner(controlFile)) {
            // Read the parameters.
            int w = scanner.nextInt();
            int s = scanner.nextInt();
            int b = scanner.nextInt();
            int K = scanner.nextInt();
            // Create the object.
            retVal = new LSHDiskSeqHash(w, s, b, K, dir);
        }
        return retVal;
    }

    /**
     * @return the stage directory for the specified stage
     *
     * @param s		relevant stage index
     */
    private File stageDir(int s) {
        return new File(this.mainDir, String.format("stage%04d", s));
    }

    @Override
    public boolean checkBucket(int s, int h) {
        boolean retVal;
        // Get the file name for this bucket.
        File bucketFile = this.getBucketFile(s, h);
        // Check for a memory copy.
        Bucket bucket = this.bucketCache.get(bucketFile);
        if (bucket != null)
            // Memory copy.  Check the number of sketches.
            retVal = (bucket.size() > 0);
        else
            // Disk copy.  Check for a bucket file.
            retVal = bucketFile.exists();
        return retVal;
    }

    @Override
    protected Bucket getBucket(int s, int h) {
        // Get the file name for this bucket.
        File bucketFile = this.getBucketFile(s, h);
        // Check for a memory copy.
        Bucket retVal = this.bucketCache.get(bucketFile);
        if (retVal == null) {
            // Here we need to get the bucket from disk.
            // First, make sure there's room in the cache.
            this.makeRoom(cacheLimit - 1);
            // Make sure the bucket exists.
            if (bucketFile.exists()) {
                try {
                    log.debug("Loading bucket from {} into cache.", bucketFile);
                    retVal = Bucket.load(bucketFile);
                } catch (Exception e) {
                    throw new RuntimeException("Error loading bucket file " + bucketFile + ": " + e.toString());
                }
            } else {
                retVal = new Bucket();
            }
            // Save the bucket in the cache.
            this.bucketCache.put(bucketFile, retVal);
        }
        return retVal;
    }

    /**
     * Insure the cache has no more than the specified number of buckets
     * in it.
     *
     * @param limit		maximum number of acceptable buckets.
     */
    private void makeRoom(int limit) {
        // Determine how many buckets we have to get rid of.
        int deletions = this.bucketCache.size() - limit;
        if (deletions > 0) {
            // We want get a set of the buckets that need to be deleted.
            SortedSet<Map.Entry<File, Bucket>> oldBuckets = new TreeSet<Map.Entry<File,Bucket>>(COMPARATOR);
            // Loop through the buckets in the hash.
            for (Map.Entry<File, Bucket> entry : this.bucketCache.entrySet())
                if (oldBuckets.size() < deletions)
                    oldBuckets.add(entry);
                else {
                    Map.Entry<File, Bucket> last = oldBuckets.last();
                    if (COMPARATOR.compare(entry, last) < 0) {
                        oldBuckets.remove(last);
                        oldBuckets.add(entry);
                }
            }
            // Now delete the old buckets.
            for (Map.Entry<File, Bucket> entry : oldBuckets) {
                Bucket bucket = entry.getValue();
                File file = entry.getKey();
                if (bucket.isModified()) {
                    this.saveBucket(file, bucket);
                }
                log.debug("Removing bucket for {} from cache.", file);
                this.bucketCache.remove(file);
            }
        }
    }

    /**
     * This is a comparator for comparing entries in the bucket cache.  We
     * compare them by bucket and then file name.
     */
    private static class EntryCompare implements Comparator<Map.Entry<File,Bucket>> {

        @Override
        public int compare(Map.Entry<File, Bucket> o1, Map.Entry<File, Bucket> o2) {
            int retVal = o1.getValue().compareTo(o2.getValue());
            if (retVal == 0)
                retVal = o1.getKey().compareTo(o2.getKey());
            return retVal;
        }

    }

    /**
     * @return the file name of the specified bucket
     *
     * @param s		relevant stage
     * @param h		bucket index in the stage
     */
    private File getBucketFile(int s, int h) {
        return new File(this.stageDir(s), String.format("bucket%06d.ser", h));
    }

    /**
     * @return the optimal number of buckets for the cache
     */
    public static int getCacheLimit() {
        return cacheLimit;
    }

    /**
     * Set a new cache limit size.  If the cache is shrunk, the actual
     * shrinkage will have no immediate effect.
     *
     * @param cacheLimit the new cache limit to set
     */
    public static void setCacheLimit(int cacheLimit) {
        LSHDiskSeqHash.cacheLimit = cacheLimit;
    }

    /**
     * Save this entire hash to disk.
     *
     * @throws IOException
     */
    public void save() throws IOException {
        // We write all the modified buckets to the appropriate bucket files.
        for (Map.Entry<File,Bucket> bucketEntry : this.bucketCache.entrySet()) {
            Bucket bucket = bucketEntry.getValue();
            if (bucket.isModified())
                this.saveBucket(bucketEntry.getKey(), bucket);
        }
    }

    /**
     * Save a bucket to a file.
     *
     * @param file		file name for the bucket
     * @param bucket	bucket to save
     */
    private void saveBucket(File file, Bucket bucket) {
        try {
            // If the bucket is empty, we have to delete the file.
            if (bucket.size() != 0) {
                log.debug("Updating bucket in {}.", file);
                bucket.save(file);
            } else if (file.exists()) {
                log.debug("Emptying bucket in {}.", file);
                FileUtils.forceDelete(file);
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    /**
     * @return the main directory
     */
    public File getMainDir() {
        return mainDir;
    }

    /**
     * @return the kmer size to use
     */
    public int getKmerSize() {
        return kmerSize;
    }

    @Override
    public void close() throws Exception {
        // Save this hash.
        this.save();
    }

}
