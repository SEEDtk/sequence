/**
 *
 */
package org.theseed.sequence;

import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastParms;

/**
 * This is the base class for sequence streams.
 *
 * @author Bruce Parrello
 *
 */
public interface SequenceStream extends Iterable<Sequence> {

    /**
     * @return TRUE if this is a protein stream, else FALSE (DNA stream)
     */
    public boolean isProtein();

    /**
     * BLAST this stream against the specified database with the specified parameters.
     *
     * @param blastDB	target BLAST database
     * @param parms		BLAST parameters to use
     *
     * @return a list of the BLAST hits
     */
    public List<BlastHit> blast(BlastDB blastDB, BlastParms parms) throws IOException, InterruptedException;

    /**
     * @return an iterator for looping through batches of data.
     *
     * @param batchSize		number of sequences to return per batch
     */
    public Iterator<SequenceStream> batchIterator(int batchSize);

}
