package org.theseed.sequence.seeds;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.theseed.p3api.KeyBuffer;
import org.theseed.p3api.P3CursorConnection;
import org.theseed.p3api.SolrFilter;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;

import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This class represents a batch of FinderKmers objects. It can be used to manage multiple genomes
 * and perform batch operations on their finder-kmers. More importantly, it loads finder-kmers from
 * the BV-BRC database.
 */
public class FinderKmerBatch {

    // FIELDS
    /** logging facility */
    private static final org.slf4j.Logger log = org.slf4j.LoggerFactory.getLogger(FinderKmerBatch.class);
    /** map of genome IDs to their corresponding FinderKmers objects */
    private final Map<String, FinderKmers> finderMap;
    /** number of validated sequences read from the database (for logging during load) */
    private int seqsIn;
    /** number of invalid sequences read from the database (for logging during load) */
    private int seqsInvalid;
    /** number of sequences kept (for logging during load) */
    private int seqsKept;

    // CONSTRUCTOR
    /**
     * Construct an empty FinderKmerBatch object.
     */
    public FinderKmerBatch() {
        this.finderMap = new HashMap<>();
    }

    /**
     * Load a batch of FinderKmers objects from the database.
     * 
     * @param roleMap   a role definition map containing the roles to use
     * @param p3        a cursor connection for the BV-BRC database
     * @param genomes   the set of genome IDs to load from the database
     * 
     * @throws IOException
     */
    public void loadBatch(RoleMap roleMap, P3CursorConnection p3, Collection<String> genomes) throws IOException {
        // Clear the counters.
        int rolesIn = 0;
        this.seqsIn = 0;
        this.seqsInvalid = 0;
        this.seqsKept = 0;
        // Loop through the roles and process each one.
        for (Role role : roleMap.objectValues()) {
            rolesIn++;
            String roleId = role.getId();
            String roleName = role.getName();
            log.info("Processing role #{} {}: {}.", rolesIn, roleId, roleName);
            // Build a filter for this role.
            SolrFilter findRole = SolrFilter.EQ("product", roleName);
            Collection<SolrFilter> criteria = Arrays.asList(findRole);
            p3.getRecords("feature", P3CursorConnection.MAX_LIMIT, 200, "genome_id", genomes,
                    "patric_id,genome_id,product,aa_sequence", criteria, x -> this.processRecord(x, role));
            log.info("Processed {} sequences. Role #{} {} in progress. {} sequences kept.", this.seqsIn, rolesIn, role, this.seqsKept);
        }
        log.info("{} sequences read, {} kept, {} invalid, {} roles processed.", this.seqsIn, this.seqsKept, this.seqsInvalid, rolesIn);
    }

    /**
     * Add a FinderKmers object to the batch.
     *
     * @param finder    the FinderKmers object to add
     */
    public void addFinderKmers(FinderKmers finder) {
        this.finderMap.put(finder.getGenomeId(), finder);
    }

    /**
     * Get the FinderKmers object for a specific genome ID.
     *
     * @param genomeId    the genome ID to look up
     *
     * @return the corresponding FinderKmers object, or null if not found
     */
    public FinderKmers getFinderKmers(String genomeId) {
        return this.finderMap.get(genomeId);
    }

    /**
     * @return the number of genomes in the batch
     */
    public int size() {
        return this.finderMap.size();
    }

    /**
     * @return the set of genome IDs for the batch
     */
    public Set<String> genomeIds() {
        return this.finderMap.keySet();
    }

    /**
     * This method processes a single feature record. We verify that the product matches the given role. If it does, we 
     * add the sequence to this finder-kmer batch.
     * 
     * @param record        JSON object containing the feature record
     * @param role          descriptor for the  role being processed
     */
    private void processRecord(JsonObject record, Role role) {
        // Denote we've read another sequence.
        this.seqsIn++;
        // Validate the product against the role name.
        String product = KeyBuffer.getString(record, "product");
        if (! role.matches(product))
            this.seqsInvalid++;
        else {
            // Add the sequence to this batch.
            String genomeId = KeyBuffer.getString(record, "genome_id");
            String protein = KeyBuffer.getString(record, "aa_sequence");
            FinderKmers finder = this.finderMap.computeIfAbsent(genomeId, x -> new FinderKmers(x));
            finder.addProteinSequence(role.getId(), protein);
            this.seqsKept++;
        }
    }

    /**
     * @return a list of all the finder-kmer objects in this batch
     */
    public List<FinderKmers> getFinders() {
        return this.finderMap.values().stream().collect(Collectors.toList());
    }

}