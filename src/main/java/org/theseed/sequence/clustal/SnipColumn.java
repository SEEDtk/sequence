/**
 *
 */
package org.theseed.sequence.clustal;

import org.theseed.locations.Location;

/**
 * This object represents a column of snips.  It is the primary output object of the snip iterator.
 *
 * @author Bruce Parrello
 *
 */
public class SnipColumn {

    // FIELDS
    /** width of this column, with respect to the alignment */
    private int width;
    /** array of feature IDs being aligned */
    private String[] fids;
    /** array of locations being aligned */
    private Location[] locs;
    /** array of snip items */
    private ISnipItem[] snips;

    /**
     * Construct a snip column.
     *
     * @param width		aligned width
     * @param fids		array of feature IDs
     * @param locs		array of aligned locations
     * @param snips		array of snip items
     */
    public SnipColumn(int width, String[] fids, Location[] locs, ISnipItem[] snips) {
        this.width = width;
        this.fids = fids;
        this.locs = locs;
        this.snips = snips;
    }

    /**
     * @return the snip for the specified row
     */
    public String getSnip(int row) {
        return this.snips[row].getChars();
    }

    /**
     * @return the full snip item for the specified row
     */
    public ISnipItem getItem(int row) {
        return this.snips[row];
    }

    /**
     * @return the location string for the specified row
     */
    public Location getLoc(int row) {
        return this.snips[row].getLoc(this.locs[row]);
    }

    /**
     * @return the location string for the specified row
     */
    public String getLocString(int row) {
        return this.snips[row].getLocString(this.locs[row]);
    }

    /**
     * @return the feature ID for the specified row
     */
    public String getFid(int row) {
        return this.fids[row];
    }

    /**
     * @return the alignment width of this column
     */
    public int getWidth() {
        return this.width;
    }

    /**
     * @return the offset of the snip item for the specified row
     *
     * @param row	row of interest
     */
    public int getOffset(int row) {
        return this.snips[row].getOffset();
    }

    /**
     * @return TRUE if the snip in the specified row is significant
     *
     * @param row	row of interest
     */
    public boolean isSignificant(int row) {
        return this.snips[row].isSignificant();
    }

    /**
     * @return the number of rows
     */
    public int getRows() {
        return this.snips.length;
    }

}
