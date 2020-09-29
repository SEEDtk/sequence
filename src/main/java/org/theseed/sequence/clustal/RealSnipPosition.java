/**
 *
 */
package org.theseed.sequence.clustal;

import org.apache.commons.lang3.StringUtils;
import org.theseed.locations.Location;
import org.theseed.sequence.ExtendedProteinRegion;

/**
 * This object represents a current position in an aligned sequence.  It is manipulated by the
 * SnipIterator to move through the alignment collecting snips.
 *
 * @author Bruce Parrello
 */
public class RealSnipPosition implements ISnipPosition {

    // FIELDS
    /** region represented by this snip */
    private ExtendedProteinRegion region;
    /** sequence being aligned */
    private String sequence;
    /** offset in sequence location of current position */
    private int offset;
    /** length in sequence location of current snip */
    private int len;
    /** snip currently being assembled */
    private StringBuilder snip;
    /** last character processed */
    private char lastChar;
    /** TRUE if the last character was different, else FALSE */
    private boolean lastDiff;
    /** offset of the first character of the snip */
    private int snipOffset;
    /** number of different characters */
    private int diffCount;

    /**
     * Construct a new snip position.
     *
     * @param region	region being aligned
     * @param sequence	aligned version of region sequence
     */
    public RealSnipPosition(ExtendedProteinRegion region, String sequence) {
        this.region = region;
        this.sequence = sequence;
        // Position at the beginning.
        this.offset = 0;
        this.len = 0;
        this.snip = new StringBuilder(20);
        this.lastChar = '.';
        this.lastDiff = false;
        this.diffCount = 0;
    }

    @Override
    public boolean check(int iPos, CharSequence chars, boolean inside) {
        int oldOffset = this.offset;
        char ch = this.sequence.charAt(iPos);
        // Push the offset past this character.
        if (ch != '-') this.offset++;
        // Return TRUE if this character is different.
        boolean retVal = ! StringUtils.contains(chars, ch);
        // If we are inside a snip, process the previous character.
        if (inside) {
            this.snip.append(this.lastChar);
            if (this.lastChar != '-') this.len++;
            if (lastDiff) this.diffCount++;
        } else {
            // Otherwise, save the old offset.  If a snip starts, this will be its offset.
            this.snipOffset = oldOffset;
        }
        // Buffer the current character.
        this.lastChar = ch;
        this.lastDiff = retVal;
        // Return TRUE if this character is different.
        return retVal;
    }

    @Override
    public ISnipItem export() {
        ISnipItem retVal = new RealSnipItem(this.snip.toString(), this.snipOffset, this.len);
        // Clear the snip buffer for the next snip.
        resetSnip();
        // Return the snip descriptor.
        return retVal;
    }

    /**
     * Reset the snip buffer for the next snip.
     */
    private void resetSnip() {
        this.snip.setLength(0);
        this.len = 0;
        this.diffCount = 0;
    }

    @Override
    public Location getLocation() {
        return this.region.getFullLocation();
    }

    @Override
    public String getFid() {
        return this.region.getLabel();
    }

    @Override
    public int getDiffCount() {
        return this.diffCount;
    }

}
