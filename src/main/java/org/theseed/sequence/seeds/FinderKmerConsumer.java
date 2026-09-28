package org.theseed.sequence.seeds;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Stream;

import org.theseed.p3api.P3CursorConnection;
import org.theseed.proteins.RoleMap;

/**
 * This object computes similarities (closeness) for finder-kmer batches and runs them through a caller-specified function.
 * The client can opt to do a dense analysis of every pair of genomes or a random sampling for efficiency.
 * 
 */
public class FinderKmerConsumer {

    // FIELDS
    /** logging facility */
    private static final org.slf4j.Logger log = org.slf4j.LoggerFactory.getLogger(FinderKmerConsumer.class);
    /** batch counter for genome stream processing */
    private int batchCounter;
    /** processing comment for log messages */
    private final String reason;

    /**
     * This enum defines the type of sampling to be performed. DENSE will process every pair of genomes, while RANDOM will
     * process pairs in a ring. Thus, DENSE is quadratic with respect to the batch size, while RANDOM is linear.
     */
    public enum SamplingType {
        /** process every possible pair of genomes */
        DENSE {
            @Override
            protected void processBatch(FinderKmerBatch batch, Consumer<Double> closeConsumer) {
                List<FinderKmers> finders = batch.getFinders();
                for (int i = 0; i < finders.size(); i++) {
                    for (int j = i + 1; j < finders.size(); j++) {
                        double closeness = finders.get(i).computeCloseness(finders.get(j));
                        closeConsumer.accept(closeness);
                    }
                }
            }
        },
        /** process pairs in a ring (linear with respect to batch size) */
        RANDOM {
            @Override
            protected void processBatch(FinderKmerBatch batch, Consumer<Double> closeConsumer) {
                // Shuffle the kmer objects to randomize the ring.
                List<FinderKmers> finders = batch.getFinders();
                Collections.shuffle(finders);
                for (int i = 1; i < finders.size(); i++) {
                    int j = i - 1;
                    double closeness = finders.get(i).computeCloseness(finders.get(j));
                    closeConsumer.accept(closeness);
                }
                double closeness = finders.get(0).computeCloseness(finders.get(finders.size() - 1));
                closeConsumer.accept(closeness);
            }
        };

        /**
         * Process the closeness values for the given list of FinderKmers.
         *
         * @param batch         the batch of FinderKmers to analyze
         * @param closeConsumer a consumer to receive the closeness values as they are computed
         */
        protected abstract void processBatch(FinderKmerBatch batch, Consumer<Double> closeConsumer);
    }

    /**
     * Construct a FinderKmerConsumer object.
     * 
     * @param myReason  a string describing the reason for creating this consumer
     */
    public FinderKmerConsumer(String myReason) {
        this.batchCounter = 0;
        this.reason = myReason;
    }

    /**
     * Process the closeness values of the genomes in a stream of genome IDs.
     * 
     * @param genomeStream  the stream of genome IDs to analyze
     * @param samplingType  the type of sampling to use
     * @param batchSize     the size of the batches to use when processing the genome stream
     * @param p3            the connection to use for database access
     * @param roleMap       role definitions to use
     * @param closeConsumer a consumer to receive the closeness values as they are computed
     * @param reason        a string describing the reason for processing this batch
     * 
     */
    public static void process(Stream<String> genomeStream, SamplingType samplingType, int batchSize, P3CursorConnection p3,
            RoleMap roleMap, Consumer<Double> closeConsumer, String reason) {
        Set<String> genomeSet = new HashSet<>(batchSize * 3);
        FinderKmerConsumer processor = new FinderKmerConsumer(reason);
        // Process the genome stream in batches of the specified size.
        genomeStream.forEach(genomeId -> {
            genomeSet.add(genomeId);
            if (genomeSet.size() >= batchSize) {
                processor.consumeBatch(roleMap, genomeSet, p3, samplingType, closeConsumer);
                genomeSet.clear();
            }
        });
        // Process the residual batch.
        if (! genomeSet.isEmpty()) {
            processor.consumeBatch(roleMap, genomeSet, p3, samplingType, closeConsumer);
        }
    }

    /**
     * Process data from a new batch of genome IDs.
     * 
     * @param roleMap       role definitions to use
     * @param genomeSet     the set of genome IDs in the new batch
     * @param p3            the connection to use for database access
     * @param samplingType  the type of sampling to use
     */
    protected void consumeBatch(RoleMap roleMap, Set<String> genomeSet, P3CursorConnection p3, SamplingType samplingType, 
        Consumer<Double> closeConsumer) {
        this.batchCounter++;
        log.info("Loading batch {} for {}", this.batchCounter, this.reason);
        FinderKmerBatch batch = new FinderKmerBatch();
        // We uncheck the IO exception to make it easier to use in streams without having to catch it explicitly.
        try {
            batch.loadBatch(roleMap, p3, genomeSet);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        log.info("Analyzing batch {} of {} genomes for {}.", this.batchCounter, genomeSet.size(), this.reason);
        samplingType.processBatch(batch, closeConsumer);
    }

}
