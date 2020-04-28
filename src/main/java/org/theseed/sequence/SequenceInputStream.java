/**
 *
 */
package org.theseed.sequence;

import java.io.Closeable;

/**
 * This is a sequence stream attached to a file.  It includes closability.
 *
 * @author Bruce Parrello
 *
 */
public interface SequenceInputStream extends SequenceStream, Closeable, AutoCloseable {

    public void close();

}
