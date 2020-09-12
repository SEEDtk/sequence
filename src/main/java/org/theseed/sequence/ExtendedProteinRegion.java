/**
 *
 */
package org.theseed.sequence;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;

/**
 * An extended protein region is a protein-encoding gene along with its upstream region.  It contains location information,
 * the dna sequence of the entire region and the source feature from which the extended region was built.
 *
 * @author Bruce Parrello
 *
 */
public class ExtendedProteinRegion extends Sequence {

    // FIELDS
    /** full location (including upstream) */
    private Location fullLocation;
    /** original feature */
    private Feature sourceFeature;
    /** cached copy of DNA kmers for distance check */
    private DnaKmers kmers;

    /**
     * Create an extended protein region from a feature.
     *
     * @param feat		source feature
     * @param upstream	length of upstream region
     */
    public ExtendedProteinRegion(Feature feat, int upstream) {
        this.sourceFeature = feat;
        // Determine the length of the contig containing the feature.
        Genome parent = feat.getParent();
        Location featLoc = feat.getLocation();
        Contig contig = parent.getContig(featLoc.getContigId());
        this.fullLocation = featLoc.expandUpstream(upstream, contig.length());
        // Get the DNA.
        this.setSequence(parent.getDna(this.fullLocation));
        // Store the label and comment.
        this.setLabel(feat.getId());
        this.setComment(feat.getFunction());
        // Denote there are no kmers cached.
        this.kmers = null;
    }

    /**
     * @return the full location
     */
    public Location getFullLocation() {
        return this.fullLocation;
    }

    /**
     * @return the source feature
     */
    public Feature getFeature() {
        return this.sourceFeature;
    }

    /**
     * @return the protein sequence
     */
    public String getProteinTranslation() {
        return this.sourceFeature.getProteinTranslation();
    }

    /**
     * @return the kmer distance to another sequence
     *
     * @param seq	other sequence to check
     */
    public double getDistance(Sequence seq) {
        if (this.kmers == null)
            this.kmers = new DnaKmers(this.getSequence());
        DnaKmers other = new DnaKmers(seq.getSequence());
        return this.kmers.distance(other);
    }

    /**
     * Create a list of extended regions for this genome. Each extended region will contain a protein and its upstream
     * region on the strand, restrained by a maximum.
     *
     * @param genome	genome of interest
     * @param limit		maximum upstream value
     */
    public static List<ExtendedProteinRegion> getGenomeExtendedProteins(Genome genome, int limit) {
        List<ExtendedProteinRegion> retVal = new ArrayList<ExtendedProteinRegion>(4000);
        // Get the features sorted by strand.
        List<Feature> features = genome.getPegs().stream().sorted(new Feature.StrandComparator()).collect(Collectors.toList());
        // Loop through the features.
        for (int i = 0; i < features.size(); i++) {
            // Get the current feature.
            Feature feat = features.get(i);
            // Get the previous feature.  This is strand-related.
            Location loc = feat.getLocation();
            Contig contig = genome.getContig(loc.getContigId());
            // Compute the upstream distance.  This requires knowing the strand.
            int upstream;
            if (loc.getDir() == '-') {
                int upstreamEdge = contig.length() + 1;
                int i0 = i + 1;
                if (i0 < features.size()) {
                    Location loc0 = features.get(i0).getLocation();
                    if (loc0.getContigId().contentEquals(contig.getId()))
                        upstreamEdge = loc0.getEnd();
                }
                upstream = upstreamEdge - loc.getEnd();
            } else {
                int upstreamEdge = 0;
                int i0 = i - 1;
                if (i0 > 0) {
                    Location loc0 = features.get(i0).getLocation();
                    if (loc0.getContigId().contentEquals(contig.getId()))
                        upstreamEdge = loc0.getEnd();
                }
                upstream = loc.getBegin() - upstreamEdge;
            }
            // We want the whole number of base pairs between the two edges.
            upstream--;
            // Constrain the upstream distance.
            if (upstream > limit) upstream = limit;
            if (upstream < 0) upstream = 0;
            // Create the region and save it.
            ExtendedProteinRegion region = new ExtendedProteinRegion(feat, upstream);
            retVal.add(region);
        }
        return retVal;
    }
}
