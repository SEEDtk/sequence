/**
 *
 */
package org.theseed.sequence.clustal;

import org.theseed.locations.Location;
import org.theseed.sequence.ExtendedProteinRegion;

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
     * @return the location for this snip
     *
     * @loc		location of the original sequence
     */
    public Location getLoc(Location loc);

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

    /**
     * @return TRUE if this snip is significant (different from the base and the wild strains)
     */
    public boolean isSignificant();

    /**
     * @return TRUE if the differing portion of the snip is inside the contig
     *
     * @param base		base snip for comparison
     * @param region	region containing the snip
     */
    public boolean isReal(String base, ExtendedProteinRegion region);

}
