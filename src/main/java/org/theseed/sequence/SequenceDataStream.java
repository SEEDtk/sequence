/**
 *
 */
package org.theseed.sequence;

import java.util.stream.Stream;

/**
 * This is the interface for an in-memory sequence stream.  It allows clearing and rebuilding the stream.
 *
 * @author Bruce Parrello
 *
 */
public interface SequenceDataStream extends SequenceStream {

    /**
     * Clear the data stream.
     */
    public void clear();

    /**
     * Add a sequence to the stream.
     *
     * @param seq	sequence to add
     */
    public SequenceDataStream add(Sequence seq);

    /**
     * @return the number of sequences in the stream
     */
    public int size();

    /**
     * @return a java stream of these sequences
     */
    public Stream<Sequence> stream();

}
