package org.theseed.sequence.seeds;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.closeTo;
import static org.hamcrest.Matchers.equalTo;
import static org.hamcrest.Matchers.greaterThan;
import static org.hamcrest.Matchers.greaterThanOrEqualTo;
import static org.hamcrest.Matchers.is;
import static org.hamcrest.Matchers.lessThan;
import static org.hamcrest.Matchers.lessThanOrEqualTo;
import static org.hamcrest.Matchers.not;
import static org.hamcrest.Matchers.nullValue;
import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
import org.theseed.p3api.P3CursorConnection;
import org.theseed.proteins.RoleMap;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;

/**
 * Test class for FinderKmer.
 * 
 * FinderKmerTest
 */
public class FinderKmerTest {

    private static final Logger log = LoggerFactory.getLogger(FinderKmerTest.class);

    @Test
    public void testFinderKmer() throws IOException {
        // Create the three FinderKmer objects.
        FinderKmers finder1 = new FinderKmers("1169293.3");
        FinderKmers finder2 = new FinderKmers("357441.63");
        FinderKmers finder3 = new FinderKmers("224308.43");
        Map<String, FinderKmers> finderMap = Map.of(
                "1169293.3", finder1,
                "357441.63", finder2,
                "224308.43", finder3
        );
        // Now fill them in from the test FASTA file.
        try (FastaInputStream fastaStream = new FastaInputStream(new File("data", "sample.fasta"))) {
            for (Sequence seq : fastaStream) {
                String genome_id = Feature.genomeOf(seq.getLabel());
                FinderKmers finder = finderMap.get(genome_id);
                assertThat(genome_id, finder, not(nullValue(FinderKmers.class)));
                finder.addProteinSequence(seq.getComment(), seq.getSequence());
            }
        }
        // At this point, we have finder-kmer objects for three genomes. Check the distances.
        double dist11 = finder1.computeDistance(finder1);
        double dist12 = finder1.computeDistance(finder2);
        double dist13 = finder1.computeDistance(finder3);
        double dist23 = finder2.computeDistance(finder3);
        double close12 = finder1.computeCloseness(finder2);
        assertThat(dist11, equalTo(0.0));
        assertThat(dist12, greaterThan(0.0));
        assertThat(dist13, greaterThan(0.0));
        assertThat(dist23, greaterThan(dist12));
        assertThat(close12, lessThan(1.0));
        assertThat(dist12 + close12, closeTo(1.0, 1e-6));
        assertThat(dist12 + dist23, greaterThan(dist13));
        double dist21 = finder2.computeDistance(finder1);
        assertThat(dist21, equalTo(dist12));
    }

    @Test
    public void testFinderKmerBatch() throws IOException {
        // Get the list of genome IDs.
        Set<String> genomeIds = TabbedLineReader.readSet(new File("data", "finderTest.tbl"), "genome_id");
        // Load the role map.
        RoleMap roleMap = RoleMap.load(new File("data", "roles.for.finder"));
        // Connect to the database.
        P3CursorConnection p3 = new P3CursorConnection();
        // Create a FinderKmerBatch and load the genomes.
        FinderKmerBatch batch = new FinderKmerBatch();
        batch.loadBatch(roleMap, p3, genomeIds);
        assertThat(batch.size(), equalTo(genomeIds.size()));
        // Compute the distances between the first genome and all the others. Memorize the genome furthest away.
        Iterator<String> genomeIter = genomeIds.iterator();
        String furthestGenomeId = null;
        double maxDist = -1.0;
        if (! genomeIter.hasNext())
            assertThat("No genomes found", false, equalTo(true));
        String firstGenomeId = genomeIter.next();
        FinderKmers firstFinder = batch.getFinderKmers(firstGenomeId);
        assertThat(firstFinder.getGenomeId(), equalTo(firstGenomeId));
        while (genomeIter.hasNext()) {
            String otherGenomeId = genomeIter.next();
            FinderKmers otherFinder = batch.getFinderKmers(otherGenomeId);
            assertThat(otherFinder.getGenomeId(), equalTo(otherGenomeId));
            double dist = firstFinder.computeDistance(otherFinder);
            double close = firstFinder.computeCloseness(otherFinder);
            assertThat(otherGenomeId, dist + close, closeTo(1.0, 1e-6));
            if (dist > maxDist) {
                maxDist = dist;
                furthestGenomeId = otherGenomeId;
            }
        }
        // Check that we found a furthest genome.
        assertThat("Furthest not found.", furthestGenomeId, not(nullValue(String.class)));
        log.info("Furthest genome ID is {}, {} from {}.", furthestGenomeId, maxDist, firstGenomeId);
        assertThat(maxDist, greaterThan(0.0));
        assertThat(furthestGenomeId, not(equalTo(firstGenomeId)));
        // Now verify the triangle inequality between the first genome, the furthest genome, and all the others.
        FinderKmers furthestFinder = batch.getFinderKmers(furthestGenomeId);
        assertThat(furthestFinder.getGenomeId(), equalTo(furthestGenomeId));
        for (String otherGenomeId : genomeIds) {
            if (! otherGenomeId.equals(firstGenomeId) && ! otherGenomeId.equals(furthestGenomeId)) {
                FinderKmers otherFinder = batch.getFinderKmers(otherGenomeId);
                double dist1 = firstFinder.computeDistance(otherFinder);
                double dist1a = otherFinder.computeDistance(firstFinder);
                double dist2 = otherFinder.computeDistance(furthestFinder);
                double sum = dist1 + dist2;
                assertThat(otherGenomeId, dist1a, equalTo(dist1));
                assertThat(otherGenomeId, sum, greaterThanOrEqualTo(maxDist));
            }
        }

        // Compute the FinderKmerStats for the batch using both sampling types.
        FinderKmerStats denseStats = FinderKmerStats.compute(batch, FinderKmerStats.SamplingType.DENSE);
        FinderKmerStats randomStats = FinderKmerStats.compute(batch, FinderKmerStats.SamplingType.RANDOM);
        assertThat(denseStats.getMax(), greaterThanOrEqualTo(randomStats.getMax()));
        assertThat(denseStats.getMin(), lessThanOrEqualTo(randomStats.getMin()));
        assertThat(denseStats.getMin(), lessThanOrEqualTo(denseStats.getMean()));
        assertThat(denseStats.getMax(), greaterThanOrEqualTo(denseStats.getMean()));
        assertThat(randomStats.getMin(), lessThanOrEqualTo(randomStats.getMean()));
        assertThat(randomStats.getMax(), greaterThanOrEqualTo(randomStats.getMean()));
        // Perform a batched run.
        FinderKmerStats batchedStats = FinderKmerStats.compute(genomeIds.stream(), FinderKmerStats.SamplingType.DENSE, 10, p3,roleMap);
        assertThat(batchedStats.getMax(), lessThanOrEqualTo(denseStats.getMax()));
        assertThat(batchedStats.getMin(), greaterThanOrEqualTo(denseStats.getMin()));
        log.info("Random min, max, mean, sdev, skew: {}, {}, {}, {}, {}", randomStats.getMin(), randomStats.getMax(), randomStats.getMean(), randomStats.getStdDev(), randomStats.getSkewness());
        log.info("Dense min, max, mean, sdev, skew: {}, {}, {}, {}, {}", denseStats.getMin(), denseStats.getMax(), denseStats.getMean(), denseStats.getStdDev(), denseStats.getSkewness());
        log.info("Batched run min, max, mean, sdev, skew: {}, {}, {}, {}, {}", batchedStats.getMin(), batchedStats.getMax(), batchedStats.getMean(), batchedStats.getStdDev(), batchedStats.getSkewness());
        // Get the stats for the big set to give us a clue as to the closeness to use for the representative subset.
        genomeIds = TabbedLineReader.readSet(new File("data", "random.genomes.tbl"), "1");
        FinderKmerStats bigStats = FinderKmerStats.compute(genomeIds.stream(), FinderKmerStats.SamplingType.DENSE, 100, p3, roleMap);
        log.info("Big set min, max, mean, sdev, skew: {}, {}, {}, {}, {}", bigStats.getMin(), bigStats.getMax(), bigStats.getMean(), bigStats.getStdDev(), bigStats.getSkewness());
    }

    /**
     * Test the FinderKmerStats computation on streams.
     * 
     * @throws IOException
     */
    @Test
    public void testFinderKmerStatsStreams() throws IOException {
        P3CursorConnection p3 = new P3CursorConnection();
        Set<String> genomeIdSet = TabbedLineReader.readSet(new File("data", "random.genomes.tbl"), "1");
        Stream<String> genomeIds = genomeIdSet.stream();
        RoleMap roleMap = RoleMap.load(new File("data", "roles.for.finder"));
        FinderKmerStats stats = FinderKmerStats.compute(genomeIds, FinderKmerStats.SamplingType.DENSE, 200, p3, roleMap);
        assertThat(stats.getMin(), lessThanOrEqualTo(stats.getMax()));
        assertThat(stats.getMean(), greaterThanOrEqualTo(stats.getMin()));
        assertThat(stats.getMean(), lessThanOrEqualTo(stats.getMax()));
        assertThat(stats.getStdDev(), greaterThanOrEqualTo(0.0));
        log.info("Stream min, max, mean, sdev, skew: {}, {}, {}, {}, {}", stats.getMin(), stats.getMax(), stats.getMean(), stats.getStdDev(), stats.getSkewness());
    }

    /**
     * Test representative-subset creation in FinderKmerBatch.
     * 
     * @throws IOException 
     */
    @Test
    public void testRepresentativeSubset() throws IOException {
        P3CursorConnection p3 = new P3CursorConnection();
        // Create a stream of genome IDs for testing.
        Set<String> genomeIdSet = TabbedLineReader.readSet(new File("data", "random.genomes.tbl"), "1");
        Stream<String> genomeIds = genomeIdSet.stream();
        // Read in the role map.
        RoleMap roleMap = RoleMap.load(new File("data", "roles.for.finder"));
        // Create the representative subset.
        FinderKmerBatch repSubset = FinderKmerBatch.createRepresentativeSubset(genomeIds, p3, roleMap, 100, 0.6);
        assertThat(repSubset, not(nullValue(FinderKmerBatch.class)));
        assertThat(repSubset.genomeIds().size(), greaterThan(0));
        for (String genomeId : repSubset.genomeIds())
            assertThat(genomeId, genomeIdSet.contains(genomeId));
        FinderKmerStats  repStats = FinderKmerStats.compute(repSubset, FinderKmerStats.SamplingType.DENSE);
        assertThat(repStats.getMax(), lessThan(0.6));
        // Build a batch of the genomes not in the representative set. This will take a lot of time and memory.
        Set<String> nonRepGenomeIds = new HashSet<>(genomeIdSet);
        nonRepGenomeIds.removeAll(repSubset.genomeIds());
        FinderKmerBatch nonRepBatch = new FinderKmerBatch();
        nonRepBatch.loadBatch(roleMap, p3, nonRepGenomeIds);
        // Insure every genome ID in the non-representative set is close to at least one representative.
        for (String genomeId : nonRepGenomeIds) {
            FinderKmers nonRepKmers = nonRepBatch.getFinderKmers(genomeId);
            boolean found = repSubset.getFinders().stream().anyMatch(finder -> finder.computeCloseness(nonRepKmers) > 0.6);
            assertThat("Non-representative genome " + genomeId + " is not close to any representative", found, is(true));
        }
        log.info("All representative-genome tests passed.");
        log.info("Max closeness for repStats = {}. {} reps out of {} genomes.", repStats.getMax(), repSubset.size(), genomeIdSet.size());
    }

    /**
     * Test multiple representative-subset creation in FinderKmerBatch.
     */
    @Test
    public void testMultipleRepresentativeSubsets() throws IOException {
        P3CursorConnection p3 = new P3CursorConnection();
        // Create a stream of genome IDs for testing.
        Set<String> genomeIdSet = TabbedLineReader.readSet(new File("data", "random.genomes.tbl"), "1");
        Stream<String> genomeIds = genomeIdSet.stream();
        // Read in the role map.
        RoleMap roleMap = RoleMap.load(new File("data", "roles.for.finder"));
        // Create multiple representative subsets.
        FinderKmerStats stats = new FinderKmerStats();
        double closenessThresholds[] = {0.6, 0.8, 0.9};
        FinderKmerBatch[] repSubsets = FinderKmerBatch.createRepresentativeSubsets(genomeIds, p3, roleMap, 100, stats, closenessThresholds);
        assertThat(repSubsets, not(nullValue(FinderKmerBatch[].class)));
        assertThat(repSubsets.length, is(3));
        // Insure every representative subset contains only valid genome IDs.
        for (FinderKmerBatch repSubset : repSubsets) {
            assertThat(repSubset.genomeIds().size(), greaterThan(0));
            for (String genomeId : repSubset.genomeIds())
                assertThat(genomeId, genomeIdSet.contains(genomeId));
        }
        // Insure that each subset is smaller than the next one.
        for (int i = 0; i < repSubsets.length - 1; i++)
            assertThat(Integer.toString(i), repSubsets[i].size(), lessThanOrEqualTo(repSubsets[i + 1].size()));
        // Verify that each representative set is indeed representative, i.e., no two genomes within the set are too close.
        for (int i = 0; i < repSubsets.length; i++) {
            FinderKmerBatch repSubset = repSubsets[i];
            double closenessThreshold = closenessThresholds[i];
            for (String genomeId1 : repSubset.genomeIds()) {
                FinderKmers kmer1 = repSubset.getFinderKmers(genomeId1);
                for (String genomeId2 : repSubset.genomeIds()) {
                    if (! genomeId1.equals(genomeId2)) {
                        FinderKmers kmer2 = repSubset.getFinderKmers(genomeId2);
                        assertThat("Genomes " + genomeId1 + " and " + genomeId2 + " are too close in rep subset " + i, 
                                kmer1.computeCloseness(kmer2), lessThanOrEqualTo(closenessThreshold));
                    }
                }
            }
        }
        // Make sure we have statistics.
        assertThat(stats, not(nullValue(FinderKmerStats.class)));
        assertThat(stats.getMean(), greaterThan(0.0));
        assertThat(stats.getMin(), lessThanOrEqualTo(stats.getMax()));
        assertThat(stats.getStdDev(), greaterThanOrEqualTo(0.0));
        assertThat(stats.getSkewness(), not(Double.NaN));
        // Now verify the statistics.
        genomeIds = genomeIdSet.stream();
        FinderKmerStats testStats = FinderKmerStats.compute(genomeIds, FinderKmerStats.SamplingType.DENSE, 100, p3, roleMap);
        assertThat(testStats, not(nullValue(FinderKmerStats.class)));
        assertThat(testStats.getMean(), closeTo(stats.getMean(), 1e-5));
        assertThat(testStats.getMin(), closeTo(stats.getMin(), 1e-5));
        assertThat(testStats.getMax(), closeTo(stats.getMax(), 1e-5));
        assertThat(testStats.getStdDev(), closeTo(stats.getStdDev(), 1e-5));
        assertThat(testStats.getSkewness(), closeTo(stats.getSkewness(), 1e-5));
        log.info("All multiple representative-subset tests passed.");

    }

}
