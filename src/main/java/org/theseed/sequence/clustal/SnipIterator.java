/**
 *
 */
package org.theseed.sequence.clustal;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.locations.Location;
import org.theseed.sequence.ExtendedProteinRegion;
import org.theseed.sequence.RegionList;
import org.theseed.sequence.Sequence;

/**
 * This object iterates through the snip columns in an alignment.  The alignment must be between multiple extended regions.
 * One region acts as the base, while zero or more others act as alternate bases.  The aligned regions are processed in
 * parallel, one position at a time, tracking whether we are inside or outside a snip.  Each time we transition from inside
 * a snip to outside, the snips are exported to produce a snip column for output.
 *
 * @author Bruce Parrello
 *
 */
public class SnipIterator implements Iterator<SnipColumn> {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(SnipIterator.class);
    /** list of wild sequences (including the base) */
    private List<String> wildSequences;
    /** aligned-sequence snip positions (base is first) */
    private ISnipPosition[] alignedPositions;
    /** array of locations being aligned */
    private Location[] alignedLocations;
    /** array of feature IDs being aligned */
    private String[] alignedFids;
    /** current position in the alignment */
    private int iPos;
    /** next column to return, or NULL if we are at the end */
    private SnipColumn next;
    /** width an the alignment */
    private int width;

    /**
     * Create an iterator for a single alignment result.
     *
     * @param regions	list of regions being aligned
     * @param alignment	list of aligned sequences
     * @param wildSet	IDs of the wild genomes (including the base)
     * @param genomeIds	IDs of all the display genomes, in order; the base should be first
     */
    public SnipIterator(RegionList regions, List<Sequence> alignment, Set<String> wildSet, List<String> genomeIds) {
        // Create the arrays.
        this.wildSequences = new ArrayList<String>(wildSet.size());
        this.alignedPositions = new ISnipPosition[genomeIds.size()];
        this.alignedLocations = new Location[genomeIds.size()];
        this.alignedFids = new String[genomeIds.size()];
        // This next section is slow, but the number of aligned sequences will be small, so we can put up with it.
        for (int i = 0; i < genomeIds.size(); i++) {
            String genomeId = genomeIds.get(i);
            ExtendedProteinRegion region = (ExtendedProteinRegion) findGenome(regions, genomeId);
            if (region == null) {
                // This genome did not participate, so make its position virtual.
                this.alignedPositions[i] = new VirtualSnipPosition();
                this.alignedLocations[i] = null;
                this.alignedFids[i] = "missing";
            } else {
                // Here we have a participating genome.
                String alignString = findGenome(alignment, genomeId).getSequence();
                this.alignedPositions[i] = new RealSnipPosition(region, alignString);
                this.alignedLocations[i] = region.getFullLocation();
                this.alignedFids[i] = region.getLabel();
            }
        }
        // Now extract the wild sequences.
        for (Sequence seq : alignment) {
            String genomeId = Feature.genomeOf(seq.getLabel());
            if (wildSet.contains(genomeId))
                this.wildSequences.add(seq.getSequence());
        }
        // Denote we're at the beginning.
        this.iPos = 0;
        this.width = this.wildSequences.get(0).length();
        // Look for the first column.
        this.findNext();
    }

    /**
     * Find the next snip column.
     */
    private void findNext() {
        // Loop until we find a column with a visible snip.
        boolean done = false;
        while (! done) {
            boolean inside = false;
            // First, we skip to the start of the next snip.  At this point, we are positioned immediately after a known "outside"
            // column.
            while (this.iPos < this.width && ! inside) {
                // Check each snip position.  If we find one that differs, turn on the INSIDE flag.
                inside = this.checkColumn(false);
                // Move to the next column.
                this.iPos++;
            }
            if (! inside) {
                // Here we have run off the end and there is no next column.
                this.next = null;
                done = true;
            } else {
                int aligned = 0;
                // Now we have to find the end of the snip.
                while (this.iPos < this.width && inside) {
                    // Check each snip position, if none differ, turn off the INSIDE flag.
                    inside = checkColumn(true);
                    // Move to the next column.
                    this.iPos++;
                    aligned++;
                }
                // Now we export the snips to create the result column.  We need to start with the base snip.
                ISnipItem[] snips = new ISnipItem[this.alignedPositions.length];
                RealSnipPosition basePosition = (RealSnipPosition) this.alignedPositions[0];
                snips[0] = basePosition.export();
                // Loop through the aligned snips.
                for (int i = 1; i < this.alignedPositions.length; i++) {
                    ISnipPosition position = this.alignedPositions[i];
                    snips[i] = position.export();
                }
                this.next = new SnipColumn(aligned, this.alignedFids, this.alignedLocations, snips);
                done = true;
            }
        }
    }

    /**
     * Check the current column for changes and accumulate snips.
     *
     * @param inside	TRUE if we are inside a snip, else FALSE
     *
     * @return TRUE if any sequence has changed from the wild values, else FALSE
     */
    boolean checkColumn(boolean inside) {
        String base = this.wildSequences.stream().map(x -> x.substring(this.iPos, this.iPos + 1)).collect(Collectors.joining());
        boolean retVal = false;
        for (ISnipPosition pos : this.alignedPositions)
            if (pos.check(this.iPos, base, inside)) retVal = true;
        return retVal;
    }

    /**
     * @return the sequence in the specified list containing a feature from the specified genome
     *
     * @param seqList	list of sequences to search
     * @param genomeId	ID of desired genome
     */
    private static Sequence findGenome(List<? extends Sequence> seqList, String genomeId) {
        Sequence retVal = null;
        for (int i = 0; retVal == null && i < seqList.size(); i++) {
            Sequence seq = seqList.get(i);
            if (Feature.genomeOf(seq.getLabel()).contentEquals(genomeId))
                retVal = seq;
        }
        return retVal;
    }

    @Override
    public boolean hasNext() {
        return (this.next != null);
    }

    @Override
    public SnipColumn next() {
        SnipColumn retVal = this.next;
        this.findNext();
        return retVal;
    }

    /**
     * This class is a snip iterable, allowing snip iteration in for-each blocks.
     */
    public static class Run implements Iterable<SnipColumn> {

        /** underlying iterator */
        private SnipIterator iter;

        /**
         * Create an iteration run for a single alignment result.
         *
         * @param regions	list of regions being aligned
         * @param alignment	list of aligned sequences
         * @param wildSet	IDs of the wild genomes (including the base)
         * @param genomeIds	IDs of all the display genomes, in order; the base should be first
         */
        public Run(RegionList regions, List<Sequence> alignment, Set<String> wildSet, List<String> genomeIds) {
            this.iter = new SnipIterator(regions, alignment, wildSet, genomeIds);
        }

        @Override
        public Iterator<SnipColumn> iterator() {
            return this.iter;
        }

    }

}
