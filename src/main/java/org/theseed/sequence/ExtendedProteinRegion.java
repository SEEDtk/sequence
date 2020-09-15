/**
 *
 */
package org.theseed.sequence;

import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Contig;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.proteins.DnaTranslator;

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
    /** DNA translator */
    private DnaTranslator xlator;

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
        // Save the upstream distance.
        this.upstreamDistance = upstream;
        // Get the DNA.
        this.setSequence(parent.getDna(this.fullLocation));
        // Store the label and comment.
        this.setLabel(feat.getId());
        this.setComment(feat.getFunction());
        // Denote there are no kmers cached.
        this.kmers = null;
        // Denote there is no translator cached.
        this.xlator = null;
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
     * @return TRUE if the protein translation is changed with the specified snip inserted
     *
     * @param offset	offset in the full location at which to insert the snip
     * @param snip		text of the snip to insert
     * @param rLen		length of the region being replaced
     */
    public boolean isChanged(int offset, String snip, int rLen) {
        // Remove all the gaps from the snip.
        String compacted = StringUtils.remove(snip, '-');
        int len = compacted.length();
        // If the snip is all gaps, we record no change and stop.
        boolean retVal = false;
        if (len > 0) {
            // If the upstream region is changed, we record a change and stop.
            retVal = (offset < this.upstreamDistance);
            if (! retVal) {
                // Here there are no modifications to the upstream region, so we may have an invisible modification
                // to the protein.  Plug in the snip to the protein part.
                String pSeq = this.sequence.substring(this.upstreamDistance, offset) + compacted + this.sequence.substring(offset + rLen);
                // Get a translator.
                if (this.xlator == null) {
                    int gc = this.sourceFeature.getParent().getGeneticCode();
                    this.xlator = new DnaTranslator(gc);
                }
                // Now we want to isolate the snip part of the protein.  Compute the first altered codon.
                int start = offset - this.upstreamDistance;
                start -= start % 3;
                // Compute the position in the protein translation of the starting codon.
                int aaPos = start / 3;
                // Insure the translation length is a whole number of codons.
                len += 2; len -= len % 3;
                // Now we translate.  Note we have to convert "start" to a 1-based position.
                String newProt;
                if (start == 0) {
                    // Because this is a PEG, the first codon gets special treatment.
                    newProt = this.xlator.pegTranslate(pSeq, start + 1, len);
                } else {
                    // Here we're a middle codon, so we do a pure translate.
                    newProt = this.xlator.translate(pSeq, start + 1, len);
                }
                String oldProt = this.getProteinTranslation();
                // Compare the translations to check for a change.
                retVal = ! oldProt.regionMatches(aaPos, newProt, 0, newProt.length());
            }
        }
        return retVal;
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
                upstream = upstreamEdge - loc.getEnd();
            } else {
                int upstreamEdge = 0;
                int i0 = this.iPos - 1;
                if (i0 > 0) {
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
