package org.theseed.sequence.seeds;

import java.util.Collections;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

/**
 * This object computes similarities (closeness) for FinderKmer batches and outputs distribution statistics about them in 
 * the form of a SummaryStatistics object. The client can opt to do a dense analysis of every pair of genomes or a random sampling for
 * efficiency.
 * 
 * FinderKmerStats
 */
public class FinderKmerStats {

    // FIELDS
    /** summary statistics for the closeness */ 
    private SummaryStatistics closeStats;

    /**
     * This enum defines the type of sampling to be performed. DENSE will process every pair of genomes, while RANDOM will
     * process pairs in a ring. Thus, DENSE is quadratic with respect to the batch size, while RANDOM is linear.
     */
    public enum SamplingType {
        /** process every possible pair of genomes */
        DENSE {
            @Override
            protected void computeStats(SummaryStatistics stats, FinderKmerBatch batch) {
                List<FinderKmers> finders = batch.getFinders();
                for (int i = 0; i < finders.size(); i++) {
                    for (int j = i + 1; j < finders.size(); j++) {
                        double closeness = finders.get(i).computeCloseness(finders.get(j));
                        stats.addValue(closeness);
                    }
                }
            }
        },
        /** process pairs in a ring (linear with respect to batch size) */
        RANDOM {
            @Override
            protected void computeStats(SummaryStatistics stats, FinderKmerBatch batch) {
                // Shuffle the kmer objects to randomize the ring.
                List<FinderKmers> finders = batch.getFinders();
                Collections.shuffle(finders);
                for (int i = 1; i < finders.size(); i++) {
                    int j = i - 1;
                    double closeness = finders.get(i).computeCloseness(finders.get(j));
                    stats.addValue(closeness);
                }
                double closeness = finders.get(0).computeCloseness(finders.get(finders.size() - 1));
                stats.addValue(closeness);
            }
        };

        /**
         * Compute the summary statistics for the given list of FinderKmers.
         *
         * @param stats       the summary statistics object to populate with closeness values
         * @param batch       the batch of FinderKmers to analyze
         */
        protected abstract void computeStats(SummaryStatistics stats, FinderKmerBatch batch);
    }

    /**
     * Construct a FinderKmerStats object for the specified batch.
     * 
     * @param batch     FinderKmerBatch to analyze
     */
    public FinderKmerStats(FinderKmerBatch batch) {
        this.closeStats = new SummaryStatistics();
    }

    /**
     * Compute the summary statistics for the current batch of FinderKmers using the specified sampling type.
     *
     * @param batch        the batch of FinderKmers to analyze
     * @param samplingType  the type of sampling to use (DENSE or RANDOM)
     * 
     * @return the summary statistics for the closeness of the FinderKmers in the current batch
     */
    public SummaryStatistics compute(FinderKmerBatch batch, SamplingType samplingType) {
        this.closeStats = new SummaryStatistics();
        samplingType.computeStats(this.closeStats, batch);
        return this.closeStats;
    }

}
