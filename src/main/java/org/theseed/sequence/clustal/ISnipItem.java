/**
 *
 */
package org.theseed.sequence.clustal;

import org.theseed.locations.Location;

/**
 * This interface represents an object that can function as a snip item.
 *
 * @author Bruce Parrello
 *
 */
public interface ISnipItem {

    /**
     * @return the location string for this snip
     *
     * @loc		location of the original sequence
     */
    public String getLocString(Location loc);

    /**
     * @return the length this snip consumes in the original sequence
     */
    public int getLen();

    /**
     * @return the characters of this snip
     */
    public String getChars();

    /**
     * @return the offset of this snip in the unaligned sequence
     */
    public int getOffset();


}
