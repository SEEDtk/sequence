/**
 *
 */
package org.theseed.sequence;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Optional;

import org.theseed.genome.Genome;
import org.theseed.proteins.Function;
import org.theseed.proteins.FunctionMap;

/**
 * This object represents a list of extended protein regions.  It contains some useful constructors as well as a method to
 * find the closest region to a given sequence.
 *
 * @author Bruce Parrello
 *
 */
public class RegionList extends ArrayList<ExtendedProteinRegion> {

    // FIELDS
    /** serialization ID */
    private static final long serialVersionUID = -1683733432282843932L;

    /**
     * Create an empty region list with the default capacity.
     */
    public RegionList() {
        super();
    }

    /**
     * Create an empty region list with the specified initial capacity.
     *
     * @param initialCapacity	starting capacity of the array list
     */
    public RegionList(int initialCapacity) {
        super(initialCapacity);
    }

    /**
     * Create a region list from the specified list collection.
     *
     * @param c		collection of regions to use for initializing the list
     */
    public RegionList(Collection<? extends ExtendedProteinRegion> c) {
        super(c);
    }

    /**
     * Create a list of extended regions for a genome. Each extended region will contain a protein and its upstream
     * region on the strand, restrained by a maximum.
     *
     * @param genome			genome of interest
     * @param upstreamLimit		maximum upstream value
     */
    public RegionList(Genome genome, int upstreamLimit) {
        super(4000);
        Iterator<ExtendedProteinRegion> iter = new ExtendedProteinRegion.GenomeIterator(genome, upstreamLimit);
        while (iter.hasNext())
            this.add(iter.next());
    }

    /**
     * Create a map of function IDs to extended region lists.  We expect most such lists to be extremely small, except for the ones
     * relating to hypothetical or ambiguous proteins.
     *
     * @param functionMap		function ID definition map; this will be modified
     * @param genome			genome of interest
     * @param upstreamLimit		maximum upstream value
     */
    public static Map<String, RegionList> createMap(FunctionMap functionMap, Genome genome, int upstreamLimit) {
        Map<String, RegionList> retVal = new HashMap<String, RegionList>(4000);
        // Get an iterator through the genome.
        Iterator<ExtendedProteinRegion> iter = new ExtendedProteinRegion.GenomeIterator(genome, upstreamLimit);
        while (iter.hasNext()) {
            ExtendedProteinRegion region = iter.next();
            // Compute the function for this region's feature.
            String funString = region.getFeature().getPegFunction();
            Function fun = functionMap.findOrInsert(funString);
            // Get the function's region list and add this region.
            RegionList rList = retVal.computeIfAbsent(fun.getId(), x -> new RegionList(5));
            rList.add(region);
        }
        return retVal;
    }

    /**
     * @return the closest region in this list to the specified sequence, or NULL if there is none
     *
     * @param seq		sequence whose closest region is desired
     * @param maxDist	maximum permissible distance
     */
    public ExtendedProteinRegion getClosest(Sequence seq, double maxDist) {
        ExtendedProteinRegion retVal = null;
        double bestDist = maxDist;
        DnaKmers kmers = new DnaKmers(seq.getSequence());
        for (ExtendedProteinRegion region : this) {
            double newDist = region.getDistance(kmers);
            // We use <= so that we override NULL if a region is found at the max distance.
            if (newDist <= bestDist) {
                bestDist = newDist;
                retVal = region;
            }
        }
        return retVal;
    }

    /**
     * Save the sequences to a file.
     *
     * @param fileName	name of the output file
     *
     * @throws IOException
     */
    public void save(File fileName) throws IOException {
        try (FastaOutputStream outStream = new FastaOutputStream(fileName)) {
            outStream.write(this);
        }
    }

    /**
     * @return the region (if any) for the specified feature
     *
     * @param fid	ID of the feature of interest
     */
    public ExtendedProteinRegion get(String fid) {
        Optional<ExtendedProteinRegion> found = this.stream().filter(x -> x.getFeature().getId().contentEquals(fid)).findFirst();
        ExtendedProteinRegion retVal = null;
        if (found.isPresent()) retVal = found.get();
        return retVal;
    }


}
