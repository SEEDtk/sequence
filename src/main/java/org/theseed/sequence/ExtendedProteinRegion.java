/**
 *
 */
package org.theseed.sequence;

import java.util.Iterator;
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
    /** upstream distance */
    private int upstreamDistance;
    /** length of the contig */
    private int contigLen;

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
        this.contigLen = contig.length();
        this.fullLocation = featLoc.expandUpstream(upstream, this.contigLen);
        // Save the upstream distance.
        this.upstreamDistance = upstream;
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
        DnaKmers other = new DnaKmers(seq.getSequence());
        return this.getDistance(other);
    }

    /**
     * @return the kmer distance to another sequence
     *
     * @param oKmers	kmers of the other sequence to check
     */
    public double getDistance(DnaKmers oKmers) {
        if (this.kmers == null)
            this.kmers = new DnaKmers(this.getSequence());
        return this.kmers.distance(oKmers);
    }

    /**
     * @return the length of the part of this region upstream from the protein
     */
    public int getUpstreamDistance() {
        return this.upstreamDistance;
    }

    /**
     * @return the upstream DNA in this region
     */
    public String getUpstreamDna() {
        return this.sequence.substring(0, this.upstreamDistance);
    }

    /**
     * @return TRUE if the specified offset is at the edge of the contig
     */
    public boolean isVirtual(int offset) {
        int offsetPoint = this.fullLocation.offsetPoint(offset);
        return (offsetPoint <= 1 || offsetPoint >= this.contigLen);
    }

    /**
     * This is an iterator for all the extended protein regions in a genome.
     */
    public static class GenomeIterator implements Iterator<ExtendedProteinRegion> {

        /** maximim upstream region size */
        private int limit;
        /** sorted list of pegs in the genome */
        private List<Feature> features;
        /** next feature to process */
        private int iPos;
        /** genome of interest */
        private Genome genome;

        /**
         * Construct an iterator for a specified genome with a specified limit on the upstream region size.
         *
         * @param genome	genome of interest
         * @param limit		maximum upstream region size
         */
        public GenomeIterator(Genome genome, int limit) {
            this.limit = limit;
            this.genome = genome;
            this.features = genome.getPegs().stream().sorted(new Feature.StrandComparator()).collect(Collectors.toList());
            this.iPos = 0;
        }

        @Override
        public boolean hasNext() {
            return (this.iPos < this.features.size());
        }

        @Override
        public ExtendedProteinRegion next() {
            // Get the current feature.
            Feature feat = features.get(this.iPos);
            // Get the previous feature.  This is strand-related.
            Location loc = feat.getLocation();
            Contig contig = genome.getContig(loc.getContigId());
            // Compute the upstream distance.  This requires knowing the strand.
            int upstream;
            if (loc.getDir() == '-') {
                int upstreamEdge = contig.length() + 1;
                int i0 = this.iPos + 1;
                if (i0 < features.size()) {
                    Location loc0 = features.get(i0).getLocation();
                    if (loc0.getContigId().contentEquals(contig.getId()))
                        upstreamEdge = loc0.getEnd();
                }
                upstream = upstreamEdge - loc.getBegin();
            } else {
                int upstreamEdge = 0;
                int i0 = this.iPos - 1;
                if (i0 >= 0) {
                    Location loc0 = features.get(i0).getLocation();
                    if (loc0.getContigId().contentEquals(contig.getId()))
                        upstreamEdge = loc0.getEnd();
                }
                upstream = loc.getBegin() - upstreamEdge;
            }
            // Position on the next feature.
            this.iPos++;
            // We want the whole number of base pairs between the two edges.
            upstream--;
            // Constrain the upstream distance.
            if (upstream > limit) upstream = limit;
            if (upstream < 0) upstream = 0;
            // Create the region.
            ExtendedProteinRegion retVal = new ExtendedProteinRegion(feat, upstream);
            return retVal;
        }

    }

}
