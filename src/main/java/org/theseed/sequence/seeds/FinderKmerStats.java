package org.theseed.sequence.seeds;

import java.util.stream.Stream;

import org.apache.commons.statistics.descriptive.DoubleStatistics;
import org.apache.commons.statistics.descriptive.Statistic;
import org.theseed.p3api.P3CursorConnection;
import org.theseed.proteins.RoleMap;

/**
 * This object computes similarities (closeness) for FinderKmer batches and outputs distribution statistics about them in 
 * the form of a DoubleStatistics object. The client can opt to do a dense analysis of every pair of genomes or a random sampling for
 * efficiency.
 * 
 * FinderKmerStats
 */
public class FinderKmerStats extends FinderKmerConsumer {

    // FIELDS
    /** summary statistics for the closeness */ 
    private final DoubleStatistics closeStats;

    /**
     * Construct a FinderKmerStats object. This is a consumer whose reason is to compute statistics.
     */
    public FinderKmerStats() {
        super("statistics");
        // Initialize the statistics object.
        this.closeStats = DoubleStatistics.of(
                    Statistic.MIN,
                    Statistic.MAX,
                    Statistic.MEAN,
                    Statistic.STANDARD_DEVIATION,
                    Statistic.SKEWNESS
                );
    }

    /**
     * Compute the statistics for a batch of FinderKmers using the specified sampling type. The
     * statistics can be retrieved from this object.
     *
     * @param batch        the batch of FinderKmers to analyze
     * @param samplingType  the type of sampling to use
     * 
     * @return a FinderKmerStats object for the closeness of the FinderKmers in the current batch
     */
    public static FinderKmerStats compute(FinderKmerBatch batch, SamplingType samplingType) {
        FinderKmerStats retVal = new FinderKmerStats();
        // Compute the statistics for the current batch using the specified sampling type.
        samplingType.processBatch(batch, retVal.closeStats::accept);
        return retVal;
    }

    /**
     * Compute the summary statistics for the closeness of the genomes in a stream of genome IDs.
     * 
     * @param genomeStream  the stream of genome IDs to analyze
     * @param samplingType  the type of sampling to use
     * @param batchSize     the size of the batches to use when processing the genome stream
     * @param p3            the connection to use for database access
     * @param roleMap       role definitions to use
     * 
     * @return a FinderKmerStats object containing the computed statistics for the genome stream
     * 
     */
    public static FinderKmerStats compute(Stream<String> genomeStream, SamplingType samplingType, int batchSize, P3CursorConnection p3,
            RoleMap roleMap) {
        FinderKmerStats retVal = new FinderKmerStats();
        process(genomeStream, samplingType, batchSize, p3, roleMap, x -> retVal.closeStats.accept(x), "statistics");
        return retVal;
    }

    /**
     * @return the minimum closeness
     */
    public double getMin() {
        return this.closeStats.getAsDouble(Statistic.MIN);
    }
    
    /**
     * @return the maximum closeness
     */
    public double getMax() {
        return this.closeStats.getAsDouble(Statistic.MAX);
    }

    /**
     * @return the mean closeness
     */
    public double getMean() {
        return this.closeStats.getAsDouble(Statistic.MEAN);
    }
    
    /**
     * @return the standard deviation of the closeness
     */
    public double getStdDev() {
        return this.closeStats.getAsDouble(Statistic.STANDARD_DEVIATION);
    }

    /**
     * @return the skewness of the closeness
     */
    public double getSkewness() {
        return this.closeStats.getAsDouble(Statistic.SKEWNESS);
    }

}
