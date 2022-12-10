/**
 *
 */
package org.theseed.sequence.seeds;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.locations.Location;
import org.theseed.p3api.Criterion;
import org.theseed.p3api.P3Connection;
import org.theseed.p3api.P3Connection.Table;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.ProteinInputStream;
import org.theseed.sequence.ProteinStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastParms;
import org.theseed.sequence.blast.DnaBlastDB;

import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This is a utility object that allows the use of BLAST to find seed proteins.  It provides methods to
 * create a protein database for a small set of proteins and methods to find the proteins in a set
 * of contig sequences.
 *
 * The object allows the maintenance of two databases:  a protein query FASTA file that can be used to
 * scan contig sequences for likely seed protein locations, and a gigantic DNA FASTA file that can be
 * used to find the closest genome for each seed.
 *
 * In the protein query file, the sequence label is the FIG ID of the protein, and the comment is the
 * protein role ID.
 *
 * In the DNA query file, the sequence label is the FIG ID of the protein, and the comment contains
 * the taxon ID and name for the species.
 *
 * @author Bruce Parrello
 *
 */
public class ProteinFinder {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ProteinFinder.class);
    /** name of directory containing the finder files */
    private File finderDir;
    /** role definition map for roles of interest */
    private RoleMap roleMap;
    /** parameter object for protein blasts */
    private BlastParms protParms;
    /** minimum e-value for protein blast */
    private int maxGapProt;
    /** minimum fraction hit length for protein blast */
    private double minFracProt;
    /** name to give to the role file */
    private static final String ROLE_FILE_NAME = "roles.for.finder";
    /** name to give to the protein FASTA */
    private static final String PROTEIN_FILE_NAME = "seedprot.fa";
    /** batch size for sequence queries */
    private static final int BATCH_SIZE = 400;
    /** minimum e-value for protein blast */
    private static double MIN_E_VAL_PROT = 1e-20;
    /** maximum gap size for protein blast */
    private static int MAX_GAP_PROT = 600;
    /** minimum fraction match length for protein blast (all pieces combined) */
    private static double MIN_FRAC_PROT = 0.5;
    /** genetic code for protein blast */
    private static int GC_PROT = 11;

    /**
     * Create a new protein-finder object.
     *
     * @param finderDir		directory containing the finder files
     * @param roleFile		role definition file to use
     *
     * @throws IOException
     */
    public ProteinFinder(File dir, File roleFile) throws IOException {
        // Insure the directory exists.
        if (! dir.isDirectory()) {
            log.info("Creating protein-finder directory {}.", dir);
            FileUtils.forceMkdir(dir);
        }
        File realRoleFile = new File(dir, ROLE_FILE_NAME);
        if (! realRoleFile.equals(roleFile)) {
            // Here we must copy the role file into the finder directory.
            log.info("Copying role file {} to {}.", roleFile, realRoleFile);
            FileUtils.copyFile(roleFile, realRoleFile);
        }
        // Initialize the finder.
        this.setup(dir);
    }

    /**
     * Create a new protein-finder object from a populated directory.
     *
     * @param finderDir		directory containing the finder files
     * @param roleFile		role definition file to use
     *
     * @throws IOException
     */
    public ProteinFinder(File dir) throws IOException {
        // Validate the directory.
        if (! dir.isDirectory())
            throw new FileNotFoundException("Finder directory " + dir + " is not found or invalid.");
        File roleFile = new File(dir, ROLE_FILE_NAME);
        if (! roleFile.exists())
            throw new FileNotFoundException("No role file found in finder directory " + dir + ".");
        // Initialize the finder.
        this.setup(dir);
    }

    /**
     * Set up this finder.  We basically load the role map and initialize the blast parms.
     *
     * @param dir		directory containing the finder files
     *
     * @throws IOException
     */
    private void setup(File dir) throws IOException {
        this.finderDir = dir;
        // Load the role map.
        File roleFile = new File(dir, ROLE_FILE_NAME);
        this.roleMap = RoleMap.load(roleFile);
        log.info("{} roles loaded from {}.", this.roleMap.size(), roleFile);
        // Set the tuning parameters from the defaults.
        this.maxGapProt = MAX_GAP_PROT;
        this.minFracProt = MIN_FRAC_PROT;
        // Set up the blast parameters.
        this.protParms = new BlastParms().maxE(MIN_E_VAL_PROT);
    }

    /**
     * Create the protein query file, if it does not already exist.
     *
     *  @param genomes		genome source containing the genomes to scan
     *
     *  @throws IOException
     */
    public void createProteinFile(GenomeSource genomes) throws IOException {
        // Our strategy is to loop through each genome's pegs, writing the ones
        // with interesting roles.
        File protFile = new File(this.finderDir, PROTEIN_FILE_NAME);
        if (! protFile.exists()) {
            log.info("Writing protein file to {}.", protFile);
            try (var protStream = new FastaOutputStream(protFile)) {
                int count = 0;
                // Loop through the genomes.
                for (Genome genome : genomes) {
                    log.info("Processing genome {}.", genome);
                    // Loop through the pegs with interesting roles.
                    Iterator<Feature> iter = genome.new Pegs();
                    while (iter.hasNext()) {
                        var peg = iter.next();
                        String prot = peg.getProteinTranslation();
                        if (! StringUtils.isBlank(prot)) {
                            // Now check the roles.
                            var roles = peg.getUsefulRoles(this.roleMap);
                            if (roles.size() > 0) {
                                // This is an interesting peg.  Write it out.
                                Sequence seq = new Sequence(peg.getId(), roles.get(0).getId(), prot);
                                protStream.write(seq);
                                count++;
                            }
                        }
                    }
                }
                log.info("{} protein sequences written to {}.", count, protFile);
            }
        }
    }

    /**
     * This is a dinky little class that describes a hit by the protein blast.
     * The hits are sorted by start location within query protein ID.  This insures
     * that if we merge hits, any overlapping merges will occur before any gapped
     * merges.  Since the merge never changes the start location, it doesn't affect
     * the ordering.
     */
    protected static class ProtHit implements Comparable<ProtHit> {

        /** location of combined hits */
        private Location loc;
        /** ID of the role hit */
        private String roleId;
        /** length of the query sequence that matched, converted to nucleotides */
        private int protLen;
        /** number of gap characters between hits */
        private int gapSize;
        /** length of relevant query sequence */

        /**
         * Create a protein hit object from a blast hit.
         *
         * @param hit	blast hit to convert
         */
        protected ProtHit(BlastHit hit) {
            // The role ID is in the query comment.
            this.roleId = hit.getQueryDef();
            // Remember the length of the matching query sequence.
            this.protLen = hit.getQueryLen() * 3;
            // Save the location in the contig hit.
            this.loc = hit.getSubjectLoc();
            // Gap size starts at 0.
            this.gapSize = 0;
        }

        /**
         * @return the length of the hit
         */
        public int getLength() {
            // The length is the matching part.  This is the total length of the location minus any gaps.
            return this.loc.getLength() - this.gapSize;
        }

        /**
         * @return the match fraction of the hit (length / protein-length)
         */
        public double getMatchFraction() {
            return this.getLength() / (double) this.protLen;
        }

        /**
         * @return the relevant role ID
         */
        public String getRole() {
            return this.roleId;
        }

        /**
         * @return the location of the protein hit
         */
        public Location getLocation() {
            return this.loc;
        }

        /**
         * This is used to determine if two hits are mergeable.  The first hit must have a begin-location less
         * than the second, which is always the case if we are processing hits in order.
         *
         * @param o		other hit to compare
         *
         * @return the gap between hits, which is negative if they overlap, and MAXVALUE if they are incompatible
         */
        protected int gapBetween(ProtHit o) {
            return this.loc.getDownstreamGap(o.loc);
        }

        /**
         * Merge another hit into this one.  This is always done after a call to "gapBetween".
         *
         * @param o		other hit to merge
         * @param gap	gap between hits
         */
        protected void merge(ProtHit o, int gap) {
            this.loc.merge(o.loc);
            if (gap > 0)
                this.gapSize += gap;
        }

        @Override
        public int compareTo(ProtHit o) {
            int retVal = this.roleId.compareTo(o.roleId);
            if (retVal == 0) {
                // The first part of the subject location is the sequence ID.
                retVal = this.loc.getContigId().compareTo(o.loc.getContigId());
                if (retVal == 0) {
                    // We don't care which strand sorts first, as long as it's consistent.
                    retVal = this.loc.getDir() - o.loc.getDir();
                    if (retVal == 0) {
                        // Here we use the downstream shift.
                        retVal = -this.loc.getDownstreamShift(o.loc);
                        // Finally, longer hits come first.  This insures we get overlaps whenever possible
                        // instead of gaps.
                        if (retVal == 0)
                            retVal = o.loc.getLength() - this.loc.getLength();
                    }
                }
            }
            return retVal;
        }

    }

    /**
     * Find the seed proteins in a FASTA file.  This method performs a protein-query BLAST against a DNA
     * database built from the file.  Adjacent hits may be combined if they are close enough together.
     *
     * @param	dnaFile		file of DNA sequences to search
     *
     * @return a map containing the set of protein locations found for each role ID
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public Map<String, List<Location>> findSeedProteins(File dnaFile) throws IOException, InterruptedException {
        var retVal = new HashMap<String, List<Location>>(this.roleMap.size() * 2);
        // Perform the BLAST to find the protein hits.
        DnaBlastDB db = DnaBlastDB.createOrLoad(dnaFile, GC_PROT);
        File protFile = new File(this.finderDir, PROTEIN_FILE_NAME);
        ProteinStream prots = new ProteinInputStream(protFile);
        var hits = db.blast(prots, this.protParms);
        // Only proceed if there is at least one hit.
        if (hits.size() <= 0)
            log.warn("No seed proteins found in {}.", dnaFile);
        else {
            log.info("{} hits returned from protein blast against {}.", hits.size(), dnaFile);
            // Sort the hits into a single sorted set as protein hits.
            var protHits = hits.stream().map(x -> new ProtHit(x)).collect(Collectors.toCollection(TreeSet<ProtHit>::new));
            Iterator<ProtHit> iter = protHits.iterator();
            // Now we loop through the protein hits, merging them whenever possible.  This variable contains
            // the hit we are merging into.
            ProtHit curr = iter.next();
            while (iter.hasNext()) {
                ProtHit next = iter.next();
                int gap = curr.gapBetween(next);
                if (gap <= this.maxGapProt) {
                    // Here the next hit is close enough to merge.
                    curr.merge(next, gap);
                    iter.remove();
                } else
                    curr = next;
            }
            log.info("{} hits remaining after merging.", protHits.size());
            // Now store the hits that are long enough in the return map.
            int stored = 0;
            for (ProtHit protHit : protHits) {
                if (protHit.getMatchFraction() >= this.minFracProt) {
                    // Here we have a good hit.  Add it to the return map.
                    List<Location> locList = retVal.computeIfAbsent(protHit.getRole(), x -> new ArrayList<Location>());
                    locList.add(protHit.getLocation());
                    stored++;
                }
            }
            log.info("{} hits for {} roles stored in seed-protein hit map.", stored, retVal.size());
        }
        return retVal;
    }

    /**
     * Create the DNA FASTA files for each of the seed proteins.  Files that already exist will be skipped.
     *
     * @param refMap	set of acceptable genomes to use
     *
     * @throws IOException
     */
    public void createDnaFiles(Map<String, Integer> refMap) throws IOException {
        // Loop through the roles. For each role, we build a DNA FASTA file by polling PATRIC.
        P3Connection p3 = new P3Connection();
        // Our first task is to get the species name for each reference genome species.
        Map<Integer, String> speciesMap = this.getSpeciesMap(p3, refMap);
        // Note that each role may have multiple aliases.  We need to query each one.  PATRIC stores roles
        // by description in the product field.
        int roleCount = 0;
        for (String role : roleMap.keySet()) {
            // Open a FASTA file for the role.
            File roleFastaFile = this.getDnaFastaFileName(role);
            // Only proceed if the file doesn't already exist.
            if (! roleFastaFile.exists()) {
                try (FastaOutputStream roleFastaStream = new FastaOutputStream(roleFastaFile)) {
                    log.info("Processing {}.", roleFastaFile);
                    // This will count the number of features output.
                    int outCount = 0;
                    for (Role roleDescriptor : roleMap.getAllById(role)) {
                        String roleName = P3Connection.clean(roleDescriptor.getName());
                        log.info("Searching for {}.  This will take several minutes.", roleDescriptor);
                        var featureRecords = p3.query(Table.FEATURE, "patric_id,product,na_sequence_md5,genome_id",
                                Criterion.EQ("product", roleName), Criterion.EQ("annotation", "PATRIC"));
                        log.info("{} features found.", featureRecords.size());
                        // We will accumulate a batch of features to process in here.
                        List<JsonObject> featureBatch = new ArrayList<JsonObject>(BATCH_SIZE);
                        for (var feature : featureRecords) {
                            // Verify the genome ID and the MD5.
                            String genomeId = P3Connection.getString(feature, "genome_id");
                            String md5 = P3Connection.getString(feature, "na_sequence_md5");
                            if (refMap.containsKey(genomeId) && ! StringUtils.isBlank(md5)) {
                                // Verify the role.
                                String function = P3Connection.getString(feature, "product");
                                var roles = Feature.usefulRoles(this.roleMap, function);
                                if (roles.stream().anyMatch(x -> x.getId().contentEquals(role))) {
                                    // Here the feature is worth keeping.  Insure there is room in the batch.
                                    if (featureBatch.size() >= BATCH_SIZE) {
                                        outCount += this.writeBatch(p3, roleFastaStream, featureBatch, refMap, speciesMap);
                                        featureBatch.clear();
                                        log.info("{} features written to {}.", outCount, roleFastaFile);
                                    }
                                    featureBatch.add(feature);
                                }
                            }
                        }
                        // Write out any residual features.
                        if (featureBatch.size() > 0)
                            outCount += this.writeBatch(p3, roleFastaStream, featureBatch, refMap, speciesMap);
                    }
                    log.info("{} features written to {}.", outCount, roleFastaFile);
                    roleCount++;
                }
            }
        }
        log.info("{} role FASTA files written.", roleCount);
    }

    /**
     * Write out the feature records in a batch.  We need to retrieve the DNA sequences and connect the species
     * name.
     *
     * @param p3				PATRIC connection for making queries
     * @param roleFastaStream	output FASTA stream
     * @param featureBatch		feature records to write
     * @param refMap			map of genome IDs to species IDs
     * @param speciesMap		map of species IDs to names
     *
     * @return the number of sequences actually written
     *
     * @throws IOException
     */
    private int writeBatch(P3Connection p3, FastaOutputStream roleFastaStream, List<JsonObject> featureBatch,
            Map<String, Integer> refMap, Map<Integer, String> speciesMap) throws IOException {
        // First, get the nucleotide sequences.
        Map<String, String> seqMap = getSequenceMap(p3, featureBatch);
        // Now loop through the batch, writing features.
        int retVal = 0;
        for (JsonObject feature : featureBatch) {
            String fid = P3Connection.getString(feature, "patric_id");
            String genomeId = P3Connection.getString(feature, "genome_id");
            int taxId = refMap.get(genomeId);
            String comment = String.format("%d\t%s", taxId, speciesMap.get(taxId));
            String sequence = seqMap.get(P3Connection.getString(feature, "na_sequence_md5"));
            if (! StringUtils.isBlank(sequence)) {
                roleFastaStream.write(new Sequence(fid, comment, sequence));
                retVal++;
            }
        }
        return retVal;
    }

    /**
     * Compute the map of sequence MD5s to sequence string for all the features in a batch.
     * @param p3			PATRIC connection for making queries
     * @param featureBatch	feature records to write
     *
     * @return a map from all the needed sequence MD5s to actual sequences
     */
    private Map<String, String> getSequenceMap(P3Connection p3, List<JsonObject> featureBatch) {
        Set<String> seq_md5s = featureBatch.stream().map(x -> P3Connection.getString(x, "na_sequence_md5")).collect(Collectors.toSet());
        var sequences = p3.getRecords(Table.SEQUENCE, seq_md5s, "md5,sequence");
        Map<String, String> retVal = sequences.entrySet().stream().collect(Collectors.toMap(x -> x.getKey(),
                x -> P3Connection.getString(x.getValue(), "sequence")));
        return retVal;
    }

    /**
     * @return the DNA FASTA file name for the specified role
     *
     * @param role	ID of the role in question
     */
    private File getDnaFastaFileName(String role) {
        return new File(this.finderDir, role + ".fna");
    }

    /**
     * Create a map of species taxon IDs to species names.  We'll need this to build the FASTA.
     *
     * @param p3		connection to PATRIC
     * @param refMap	reference genome map
     *
     * @return a map from species taxon IDs to species names
     */
    private Map<Integer, String> getSpeciesMap(P3Connection p3, Map<String, Integer> refMap) {
        // Get all the species IDs as strings.
        Set<String> species = refMap.values().stream().map(x -> x.toString()).collect(Collectors.toSet());
        // Create the return hash.
        var retVal = new HashMap<Integer, String>(species.size() * 4 / 3 + 1);
        var taxRecords = p3.getRecords(Table.TAXONOMY, species, "taxon_id,taxon_name");
        for (var taxRecordEntry : taxRecords.entrySet()) {
            int speciesId = Integer.valueOf(taxRecordEntry.getKey());
            retVal.put(speciesId, P3Connection.getString(taxRecordEntry.getValue(), "taxon_name"));
        }
        log.info("{} species names loaded.", retVal.size());
        return retVal;
    }

    /**
     * This is a small class that contains a closest-genome ID, a species taxon ID, and a species name.
     * It is used to indicate a hit against a reference genome by a contig.  We sort these by taxonomic
     * ID and then contig ID.
     */
    public static class DnaHit implements Comparable<DnaHit> {

        /** ID of the query contig */
        private String contigId;
        /** ID of the proposed closest genome */
        private String refId;
        /** taxonomic ID of the species */
        private int taxId;
        /** name of the species */
        private String name;
        /** total hit score */
        private double score;
        /** high hit score */
        private double hiScore;
        /** total match length */
        private int matchLen;

        /**
         * Construct a new DNA hit.  We save the score and the taxonomy information from the hit.
         *
         * @param hit	BLAST hit from which to construct this object
         */
        protected DnaHit(BlastHit hit) {
            // The contig ID is stored as the query comment.
            this.contigId = hit.getQueryDef();
            // The tax ID and name are parsed from the subject comment.
            String[] pieces = StringUtils.split(hit.getSubjectDef(), '\t');
            this.taxId = Integer.valueOf(pieces[0]);
            this.name = pieces[1];
            // The reference genome ID is in the subject ID.
            this.refId = Feature.genomeOf(hit.getSubjectId());
            // The score is the bit score.
            this.score = hit.getBitScore();
            this.hiScore = this.score;
            // Save the match length.
            this.matchLen = hit.getAlignLen();
        }

        /**
         * Merge another blast hit into this one.  If it is for the same contig and species, we keep the best hit
         * and add the score. It should compare equal:  that is, it should have the same contig ID and species taxon ID.
         * If it is a better hit, we keep its reference-genome selection.
         *
         * @param hit	DNA hit to merge into this one
         */
        protected void merge(DnaHit other) {
            this.score += other.hiScore;
            this.matchLen += other.matchLen;
            if (other.hiScore > this.hiScore) {
                // Here the new hit is better.
                this.hiScore = other.hiScore;
                this.refId = other.refId;
            }
        }

        /**
         * @return the ID of the query contig for this hit
         */
        public String getContigId() {
            return this.contigId;
        }
        /**
         * @return the ID of the proposed reference genome
         */
        public String getRefId() {
            return this.refId;
        }
        /**
         * @return the species ID
         */
        public int getTaxId() {
            return this.taxId;
        }
        /**
         * @return the species name
         */
        public String getName() {
            return this.name;
        }

        /**
         * @return the score
         */
        public double getScore() {
            return this.score;
        }

        @Override
        public int compareTo(DnaHit o) {
            int retVal = this.taxId - o.taxId;
            if (retVal == 0)
                retVal = this.contigId.compareTo(o.contigId);
            return retVal;
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + ((this.contigId == null) ? 0 : this.contigId.hashCode());
            result = prime * result + this.taxId;
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) {
                return true;
            }
            if (!(obj instanceof DnaHit)) {
                return false;
            }
            DnaHit other = (DnaHit) obj;
            if (this.contigId == null) {
                if (other.contigId != null) {
                    return false;
                }
            } else if (!this.contigId.equals(other.contigId)) {
                return false;
            }
            if (this.taxId != other.taxId) {
                return false;
            }
            return true;
        }

        /**
         * This is an alternate sorter that sorts the highest-scoring object to the front.
         */
        public static class ScoreSorter implements Comparator<DnaHit> {

            /**
             * @return a negative number if o1 has a higher score, a positive number if o2 has a higher score
             */
            @Override
            public int compare(DnaHit o1, DnaHit o2) {
                int retVal = Double.compare(o2.score, o1.score);
                if (retVal == 0) {
                    retVal = o2.matchLen - o1.matchLen;
                }
                return retVal;
            }

            /**
             * @return the best hit
             *
             * @param o1	first hit to compare
             * @param o2	second hit to compare
             */
            public DnaHit max(DnaHit o1, DnaHit o2) {
                int cmp = this.compare(o1, o2);
                return (cmp <= 0 ? o1 : o2);
            }

        }

    }

    /**
     * Use the DNA FASTA files to determine the closest genome to each hit in a protein hit role map.
     *
     * @param hitMap	map of role IDs to hit locations from "findSeedProteins"
     *
     * @return a map from contig IDs to closest-genome IDs
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public Map<String, DnaHit> findRefGenomes(Map<String, Collection<Location>> protHitMap) throws IOException, InterruptedException {
        // This is our working map.  For each contig, we remember the best hit against each species, and accumulate the scores
        // of other hits into it.  At the end, each contig will have a set of votesIt collects all the hits for each contig.
        var votingMap = new HashMap<String, Map<Integer, DnaHit>>(100);
        for (var protHitsEntry : protHitMap.entrySet()) {
            // Determine the role and the name of the BLAST database for it.
            String roleId = protHitsEntry.getKey();
            var locs = protHitsEntry.getValue();
            log.info("Processing {} hits for role {}.", locs.size(), roleId);
            File roleFileName = this.getDnaFastaFileName(roleId);
            DnaBlastDB db = DnaBlastDB.createOrLoad(roleFileName, GC_PROT);
            // TODO blast for this role and convert hits to DnaHits, keeping merged best hit for each species
        }
        // TODO collate votes to find the refID for each contig
        return null;
    }

}
