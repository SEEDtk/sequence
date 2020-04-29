/**
 *
 */
package org.theseed.sequence;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.stream.Stream;

/**
 * This is a protein stream created from an in-memory list.
 *
 * @author Bruce Parrello
 *
 */
public class ProteinDataStream extends ProteinStream implements SequenceDataStream {

    // FIELDS
    private Collection<Sequence> sequences;

    /**
     * Create a protein stream from an in-memory list.
     *
     * @param sequences		collection of sequences to use
     */
    public ProteinDataStream(Collection<Sequence> sequences) {
        this.sequences = sequences;
    }

    /**
     * Create an empty protein data stream for a specified batch size.
     *
     * @param batchSize		batch size to use
     */
    public ProteinDataStream(int batchSize) {
        this.sequences = new ArrayList<Sequence>(batchSize);
    }

    @Override
    public Iterator<Sequence> iterator() {
        return sequences.iterator();
    }

    /**
     * Clear the data stream.
     */
    public void clear() {
        this.sequences.clear();
    }

    /**
     * Add a sequence to the stream.
     *
     * @param seq	sequence to add
     */
    public SequenceDataStream add(Sequence seq) {
        this.sequences.add(seq);
        return this;
    }

    /**
     * @return the number of sequences in the stream
     */
    public int size() {
        return this.sequences.size();
    }

    @Override
    public Stream<Sequence> stream() {
        return this.sequences.stream();
    }

}
