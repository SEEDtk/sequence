/**
 *
 */
package org.theseed.sequence.clustal;

import org.theseed.locations.Location;

/**
 * This class represents a single snip instance.  We think of a snip as being a set of subsequences that differ from the base
 * sequences.  This object represents a single one of those subsequences.
 *
 * @author Bruce Parrello
 *
 */
public class RealSnipItem implements ISnipItem {

    // FIELDS
    /** position in the original sequence, as an offset from the beginning */
    private int offset;
    /** length in the original sequence */
    private int len;
    /** characters of the snip */
    private String text;

    /**
     * Construct a snip item.
     *
     * @param text		characters of the snip
     * @param offset	position in the original sequence, as a 0-based offset
     * @param len		length in the original sequence (may be zero)
     */
    public RealSnipItem(String text, int offset, int len) {
        this.offset = offset;
        this.text = text;
        this.len = len;
    }

    @Override
    public String getLocString(Location loc) {
        return String.format("%s_%d%c%d", loc.getContigId(), loc.offsetPoint(this.offset), loc.getDir(), this.len);
    }

    @Override
    public int getLen() {
        return this.len;
    }

    @Override
    public String getChars() {
        return this.text;
    }

    @Override
    public int getOffset() {
        return this.offset;
    }

}
