/**
 *
 */
package org.theseed.sequence.clustal;

import org.theseed.locations.Location;

/**
 * This represents an object that can function as a snip position.  A real snip position corresponds to an aligned sequence, while
 * a virtual snip position corresponds to a missing sequence.  They are implemented very differently.
 *
 * @author Bruce Parrello
 *
 */
public interface ISnipPosition {

    /**
     * Check the specified position in the alignment.
     *
     * @param iPos		position in the alignment to check; it should be the next physical position
     * @param chars		characters at the current position in the base sequences
     * @param inside	TRUE if we are accumulating a snip, FALSE if we are between snips
     *
     * @return TRUE if the current position differs from the base sequences, else FALSE
     */
    public boolean check(int iPos, CharSequence chars, boolean inside);

    /**
     * @return the snip item just found and reset for the next snip
     */
    public ISnipItem export();

    /**
     * @return the location of this position's region
     */
    public Location getLocation();

    /**
     * @return the ID of this position's feature
     */
    public String getFid();

    /**
     * @return the difference count
     */
    public int getDiffCount();

}
