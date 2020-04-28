/**
 *
 */
package org.theseed.sequence;

import java.util.Iterator;

/**
 * This is an iterator that batvh
 * @author parrello
 *
 */
public class BatchStreamIterator implements Iterator<SequenceStream> {

    // FIELDS
    /** iterator through the underlying sequence stream */
    private Iterator<Sequence> iter;
    /** buffer stream to contain each batch */
    private SequenceDataStream batch;
    /** batch size */
    int batchSize;

    /**
     * Create an iterator through batches of a source stream.
     *
     * @param source		source stream to iterate through
     * @param buffer		buffer stream to contain each batch
     * @param batchSize		size of each batch
     */
    public BatchStreamIterator(SequenceStream source, SequenceDataStream buffer, int batchSize) {
        this.iter = source.iterator();
        this.batch = buffer;
        this.batchSize = batchSize;
    }

    @Override
    public boolean hasNext() {
        return iter.hasNext();
    }

    @Override
    public SequenceStream next() {
        SequenceStream retVal = null;
        if (this.hasNext()) {
            // Here we have another batch.  Fill it in.
            this.batch.clear();
            while (this.iter.hasNext() && this.batch.size() < this.batchSize)
                this.batch.add(this.iter.next());
            retVal = this.batch;
        }
        return retVal;
    }

}
