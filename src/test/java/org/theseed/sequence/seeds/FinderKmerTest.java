package org.theseed.sequence.seeds;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.closeTo;
import static org.hamcrest.Matchers.equalTo;
import static org.hamcrest.Matchers.greaterThan;
import static org.hamcrest.Matchers.lessThan;
import static org.hamcrest.Matchers.not;
import static org.hamcrest.Matchers.nullValue;
import org.junit.jupiter.api.Test;
import org.theseed.genome.Feature;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;

/**
 * Test class for FinderKmer.
 * 
 * FinderKmerTest
 */
public class FinderKmerTest {

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
                assertThat(genome_id, finder,not(nullValue(FinderKmers.class)));
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

}
