package org.theseed.sequence.seeds;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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
    /** batch counter for representative computation */
    private int batchCount;
    /** genome counter for representative computation */
    private int genomeCount;

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

    /**
     * Create a representative subset of the genomes in a genome ID stream and output them as a new FinderKmerBatch.
     * The representatives have the property that no two are within the specified closeness threshold of each other.
     * The incoming stream of genome IDs is processed in batches to improve performance when retrieving from the
     * database.
     * 
     * @param genomeIds    stream of genome IDs to consider for the representative subset
     * @param roleMap      role definition file for roles to use in selecting proteins
     * @param batchSize    size of the batches to process
     * @param closeness    the closeness threshold for selecting representative genomes
     *
     * @return a new FinderKmerBatch containing the representative genomes
     */
    public static FinderKmerBatch createRepresentativeSubset(Stream<String> genomeIds, RoleMap roleMap, int batchSize, double closeness) {
        // We will build the representative set in here.
        FinderKmerBatch retVal = new FinderKmerBatch();
        // We process the genomes a batch at a time. After a batch is encountered, we check for representatives and then
        // build the next one.
        Set<String> genomeBatch = new HashSet<>(batchSize * 3);
        // Connect to the database.
        P3CursorConnection p3 = new P3CursorConnection();
        // This counts the batches and genomes.
        retVal.batchCount = 0;
        retVal.genomeCount = 0;
        // Now process the stream.
        genomeIds.forEach(genomeId -> {
            retVal.genomeCount++;
            genomeBatch.add(genomeId);
            if (genomeBatch.size() >= batchSize) {
                retVal.batchCount++;
                log.info("Processing batch number {} with {} genomes.", retVal.batchCount, genomeBatch.size());
                retVal.processBatch(genomeBatch, roleMap, p3, closeness);
                genomeBatch.clear();
            }
        });
        // Process any remaining genomes in the last batch.
        if (! genomeBatch.isEmpty()) {
            log.info("Processing final batch with {} genomes.", genomeBatch.size());
            retVal.processBatch(genomeBatch, roleMap, p3, closeness);
            retVal.batchCount++;
        }
        log.info("Processed a total of {} batches with {} representatives found out of {} genomes.", retVal.batchCount, retVal.finderMap.size(), retVal.genomeCount);
        return retVal;
    }

    /**
     * This method is used to process a batch of genomes and update the representative set.
     *
     * @param genomeBatch  the set of genome IDs in the current batch
     * @param roleMap      role definition file for roles to use in selecting proteins
     * @param p3           connection to the BV-BRC database
     * @param closeness    the closeness threshold for selecting representative genomes
     */
    private void processBatch(Set<String> genomeBatch, RoleMap roleMap, P3CursorConnection p3, double closeness) {
        // Create a new FinderKmerBatch from the incoming genomes.
        FinderKmerBatch newBatch = new FinderKmerBatch();
        try {
            newBatch.loadBatch(roleMap, p3, genomeBatch);
        } catch (IOException e) {
            // We uncheck the IO exception so we can use this method in streams.
            throw new UncheckedIOException(e);
        }
        // Now loop through the batch. Any genome not close to an existing representative should be added to the representative set.
        for (FinderKmers genome : newBatch.finderMap.values()) {
            // We need to check this genome against the representatives already in this batch. If the
            // new genome is not close to any of them, we add it.
            Iterator<FinderKmers> newIter = this.finderMap.values().iterator();
            boolean isClose = false;
            while (! isClose && newIter.hasNext()) {
                FinderKmers repKmers = newIter.next();
                double repCloseness = repKmers.computeCloseness(genome);
                isClose = (repCloseness <= closeness);
            }
            if (! isClose) {
                this.finderMap.put(genome.getGenomeId(), genome);
            }
        }
    }

}