package org.theseed.sequence.seeds;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.closeTo;
import static org.hamcrest.Matchers.equalTo;
import static org.hamcrest.Matchers.greaterThan;
import static org.hamcrest.Matchers.greaterThanOrEqualTo;
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
        FinderKmerStats stats = new FinderKmerStats();
        SummaryStatistics denseStats = stats.compute(batch, FinderKmerStats.SamplingType.DENSE);
        SummaryStatistics randomStats = stats.compute(batch, FinderKmerStats.SamplingType.RANDOM);
        assertThat(denseStats.getMax(), greaterThanOrEqualTo(randomStats.getMax()));
        assertThat(denseStats.getMin(), lessThanOrEqualTo(randomStats.getMin()));
        assertThat(denseStats.getMin(), lessThanOrEqualTo(denseStats.getMean()));
        assertThat(denseStats.getMax(), greaterThanOrEqualTo(denseStats.getMean()));
        assertThat(randomStats.getMin(), lessThanOrEqualTo(randomStats.getMean()));
        assertThat(randomStats.getMax(), greaterThanOrEqualTo(randomStats.getMean()));
        // Perform a batched run.
        SummaryStatistics batchedStats = FinderKmerStats.compute(genomeIds.stream(), FinderKmerStats.SamplingType.DENSE, 10, roleMap);
        assertThat(batchedStats.getMax(), lessThanOrEqualTo(denseStats.getMax()));
        assertThat(batchedStats.getMin(), greaterThanOrEqualTo(denseStats.getMin()));
        log.info("Random min, max, mean, sdev: {}, {}, {}, {}", randomStats.getMin(), randomStats.getMax(), randomStats.getMean(), randomStats.getStandardDeviation());
        log.info("Dense min, max, mean, sdev: {}, {}, {}, {}", denseStats.getMin(), denseStats.getMax(), denseStats.getMean(), denseStats.getStandardDeviation());
        log.info("Batched run min, max, mean, sdev: {}, {}, {}, {}", batchedStats.getMin(), batchedStats.getMax(), batchedStats.getMean(), batchedStats.getStandardDeviation());
        // Get the stats for the big set to give us a clue as to the closeness to use for the representative subset.
        genomeIds = TabbedLineReader.readSet(new File("data", "random.genomes.tbl"), "1");
        SummaryStatistics bigStats = FinderKmerStats.compute(genomeIds.stream(), FinderKmerStats.SamplingType.DENSE, 100, roleMap);
        log.info("Big set min, max, mean, sdev: {}, {}, {}, {}", bigStats.getMin(), bigStats.getMax(), bigStats.getMean(), bigStats.getStandardDeviation());
    }

    /**
     * Test representative-subset creation in FinderKmerBatch.
     * 
     * @throws IOException 
     */
    @Test
    public void testRepresentativeSubset() throws IOException {
        // Create a stream of genome IDs for testing.
        Set<String> genomeIdSet = TabbedLineReader.readSet(new File("data", "random.genomes.tbl"), "1");
        Stream<String> genomeIds = genomeIdSet.stream();
        // Read in the role map.
        RoleMap roleMap = RoleMap.load(new File("data", "roles.for.finder"));
        // Create the representative subset.
        FinderKmerBatch repSubset = FinderKmerBatch.createRepresentativeSubset(genomeIds, roleMap, 100, 0.6);
        assertThat(repSubset, not(nullValue(FinderKmerBatch.class)));
        assertThat(repSubset.genomeIds().size(), greaterThan(0));
        for (String genomeId : repSubset.genomeIds())
            assertThat(genomeId, genomeIdSet.contains(genomeId));
        FinderKmerStats stats = new FinderKmerStats();
        SummaryStatistics repStats = stats.compute(repSubset, FinderKmerStats.SamplingType.DENSE);
        assertThat(repStats.getMin(), greaterThan(0.6));
    }

}
