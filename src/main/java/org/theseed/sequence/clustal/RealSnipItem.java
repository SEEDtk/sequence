/**
 *
 */
package org.theseed.sequence.clustal;

import org.theseed.locations.Location;
import org.theseed.sequence.ExtendedProteinRegion;

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
    /** amino acid name map */
    private static final String[] AA_MAP = aaMap();

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

    /**
     * @return an array that can be used to map amino acid letters to names
     */
    private static String[] aaMap() {
        String[] retVal = new String[128];
        retVal['A'] = "Alanine"; retVal['R'] = "Arginine"; retVal['N'] = "Asparagine"; retVal['D'] = "Aspartic Acid";
        retVal['C'] = "Cysteine"; retVal['E'] = "Glutamic Acid"; retVal['Q'] = "Glutamine"; retVal['G'] = "Glycine";
        retVal['H'] = "Histidine"; retVal['I'] = "Isoleucine"; retVal['L'] = "Leucine"; retVal['K'] = "Lysine";
        retVal['M'] = "Methionine"; retVal['F'] = "Phenylalanine"; retVal['P'] = "Proline"; retVal['S'] = "Serine";
        retVal['T'] = "Threonine"; retVal['W'] = "Tryptophan"; retVal['Y'] = "Tyrosine"; retVal['V'] = "Valine";
        retVal['X'] = "Unknown"; retVal['-'] = "gap"; retVal[' '] = "upstream"; retVal['*'] = "stop";
        return retVal;
    }

    /**
     * Compute an array that describes the amino acid at each location of this snip.  There is one array item for each
     * snip character in the alignment.  Upstream locations are "upstream".  Instream gap locations are "gap".
     *
     * @param region	extended protein region for this row of the alignment
     *
     * @return an array containing the amino acid name for each location
     */
    public String[] getProteinMap(ExtendedProteinRegion region) {
        String snip = this.getChars();
        String[] retVal = new String[snip.length()];
        String protein = region.getProteinTranslation();
        // Determine the offset, relative to the start codon, of the snip's first character.
        int offsetI = this.offset - region.getUpstreamDistance();
        // Loop through the snip.
        for (int i = 0; i < snip.length(); i++) {
            // A negative offset puts us upstream.
            if (offsetI < 0)
                retVal[i] = AA_MAP[' '];
            else {
                int protI = offsetI / 3;
                // We have a stop if we are past the end.
                if (protI >= protein.length())
                    retVal[i] = AA_MAP['*'];
                else {
                    // Here we are in the protein.
                    String name = AA_MAP[protein.charAt(protI)];
                    if (name == null)
                        retVal[i] = AA_MAP['X'];
                    else
                        retVal[i] = name;
                }
            }
            if (snip.charAt(i) != '-') offsetI++;
        }
        return retVal;
    }

    @Override
    public String getLocString(Location loc) {
        return String.format("%s_%d%c%d", loc.getContigId(), loc.offsetPoint(this.offset), loc.getDir(), this.len);
    }

    @Override
    public Location getLoc(Location loc) {
        return loc.subLocation(this.offset, this.len);
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

    @Override
    public boolean isSignificant() {
        return true;
    }

    @Override
    public boolean isReal(String base, ExtendedProteinRegion region) {
        int offset2 = this.offset;
        boolean retVal = false;
        // Loop until we find a real difference.
        for (int i = 0; ! retVal && i < base.length(); i++) {
            char c = this.text.charAt(i);
            char cBase = base.charAt(i);
            if (c != cBase)
                retVal = ! region.isVirtual(offset2);
            if (c != '-') offset2++;
        }
        return retVal;
    }


}
