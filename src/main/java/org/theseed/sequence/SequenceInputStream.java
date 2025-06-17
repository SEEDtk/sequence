/**
 *
 */
package org.theseed.sequence;


/**
 * This is a sequence stream attached to a file.  It includes closability.
 *
 * @author Bruce Parrello
 *
 */
public interface SequenceInputStream extends SequenceStream, AutoCloseable {

    public void close();

}
