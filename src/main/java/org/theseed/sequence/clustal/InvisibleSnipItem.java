/**
 *
 */
package org.theseed.sequence.clustal;

/**
 * This object represents a snip item where there is no change from the base sequence.
 *
 * @author Bruce Parrello
 *
 */
public class InvisibleSnipItem extends RealSnipItem {


    /**
     * Construct a snip item.
     *
     * @param text		characters of the snip
     * @param offset	position in the original sequence, as a 0-based offset
     * @param len		length in the original sequence (may be zero)
     */
    public InvisibleSnipItem(String text, int offset, int len) {
        super(text, offset, len);
    }

    @Override
    public String getChars() {
        return "";
    }

}
