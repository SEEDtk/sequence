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
 * This is a base class that provides a list of protein sequences.
 *
 * @author Bruce Parrello
 *
 */
public abstract class ProteinStream implements SequenceStream {

    @Override
    public abstract Iterator<Sequence> iterator();

    @Override
    public boolean isProtein() {
        return true;
    }

    /**
     * BLAST this stream against the specified database with the specified parameters.
     *
     * @param blastDB	target BLAST database
     * @param parms		BLAST parameters to use
     *
     * @return a list of the BLAST hits
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public List<BlastHit> blast(BlastDB blastDB, BlastParms parms) throws IOException, InterruptedException {
        return blastDB.blast(this, parms);
    }

    /**
     * @return an iterator for looping through batches of data.
     *
     * @param source		source stream to loop through
     * @param batchSize		number of sequences to return
     */
    public Iterator<SequenceDataStream> batchIterator(int batchSize) {
        SequenceDataStream buffer = new ProteinDataStream(batchSize);
        return new BatchStreamIterator(this, buffer, batchSize);
    }


}
