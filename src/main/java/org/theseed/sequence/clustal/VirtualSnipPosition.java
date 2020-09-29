/**
 *
 */
package org.theseed.sequence.clustal;

import org.theseed.locations.Location;

/**
 * This object represents a snip position in a missing sequence.  It supports all of the characteristics of a real position, but
 *
 * @author Bruce Parrello
 *
 */
public class VirtualSnipPosition implements ISnipPosition {

    /**
     * Initialize a virtual snip position.
     */
    public VirtualSnipPosition() {
    }

    @Override
    public boolean check(int iPos, CharSequence chars, boolean inside) {
        return false;
    }

    @Override
    public ISnipItem export() {
        ISnipItem retVal = new VirtualSnipItem();
        return retVal;
    }

    @Override
    public Location getLocation() {
        return null;
    }

    @Override
    public String getFid() {
        return "missing";
    }

    @Override
    public int getDiffCount() {
        return 0;
    }

}
