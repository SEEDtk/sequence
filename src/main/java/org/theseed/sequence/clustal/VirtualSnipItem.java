/**
 *
 */
package org.theseed.sequence.clustal;

import org.theseed.locations.Location;

/**
 * A virtual snip item represents a snip in a missing sequence.  This is a sequence in a genome that has no region that aligns
 * with the base sequence.
 *
 * @author Bruce Parrello
 *
 */
public class VirtualSnipItem implements ISnipItem {

    /**
     * Construct a virtual snip item.
     */
    public VirtualSnipItem() {
    }

    @Override
    public String getLocString(Location loc) {
        return "";
    }

    @Override
    public int getLen() {
        return 0;
    }

    @Override
    public String getChars() {
        return "";
    }

    @Override
    public int getOffset() {
        return 0;
    }

    @Override
    public boolean isSignificant() {
        return false;
    }

}
