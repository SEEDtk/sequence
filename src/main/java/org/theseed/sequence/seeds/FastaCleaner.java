/**
 *
 */
package org.theseed.sequence.seeds;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.sequence.DnaKmers;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.ProteinKmers;
import org.theseed.sequence.Sequence;

/**
 * This object cleans a FASTA file.  The file is copied and any sequence with an ambiguity character is removed.
 * There are subclasses for protein and DNA files.
 *
 * @author Bruce Parrello
 *
 */
public abstract class FastaCleaner {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(FastaCleaner.class);

    /**
     * Clean a FASTA file.
     *
     * @param inFile	name of the file to clean
     *
     * @throws IOException
     */
    public void clean(File inFile) throws IOException {
        // Get the incoming file's parent directory.
        File parent = inFile.getParentFile();
        // Create a temporary file to copy into.
        File tempFile = File.createTempFile("temp", ".fasta", parent);
        // Copy all the clean sequences into the output file.
        boolean done = false;
        try (FastaInputStream inStream = new FastaInputStream(inFile);
                FastaOutputStream outStream = new FastaOutputStream(tempFile)) {
            int inCount = 0;
            int outCount = 0;
            for (Sequence seq : inStream) {
                inCount++;
                this.normalize(seq);
                if (this.isClean(seq.getSequence())) {
                    outStream.write(seq);
                    outCount++;
                }
                if (log.isInfoEnabled() && inCount % 5000 == 0)
                    log.info("{} sequences read, {} copied.", inCount, outCount);
            }
            log.info("{} sequences read from {}, {} copied.", inCount, inFile, outCount);
            done = true;
        } finally {
            if (done) {
                // The copy worked.  Delete the old file and rename the new one.
                log.info("Renaming temp file {} to {}.", tempFile, inFile);
                FileUtils.forceDelete(inFile);
                FileUtils.moveFile(tempFile, inFile);
            } else {
                // The copy failed.  Delete the temp file.
                log.info("Copy failed.  Deleting {}.", tempFile);
                FileUtils.forceDelete(tempFile);
            }
        }
    }

    /**
     * @return TRUE if the sequence has no ambiguity characters, else FALSE
     *
     * @param sequence	sequence to check
     */
    protected abstract boolean isClean(String sequence);

    /**
     * Normalize the specified sequence to the proper case.  The sequence
     * text is modified in place.
     *
     * @param seq	sequence to normalize
     */
    protected abstract void normalize(Sequence seq);

    /**
     * This subclass cleans DNA FASTA files.
     */
    public static class Dna extends FastaCleaner {

        @Override
        protected boolean isClean(String sequence) {
            return DnaKmers.isClean(sequence);
        }

        @Override
        protected void normalize(Sequence seq) {
            // DNA is lower-case
            seq.setSequence(seq.getSequence().toLowerCase());
        }

    }

    /**
     * This subclass cleans protein FASTA files.
     */
    public static class Protein extends FastaCleaner {

        @Override
        protected boolean isClean(String sequence) {
            return ProteinKmers.isClean(sequence);
        }

        @Override
        protected void normalize(Sequence seq) {
            // Proteins are upper-case
            seq.setSequence(seq.getSequence().toUpperCase());
        }

    }
}
