/**
 *
 */
package org.theseed.sequence;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.theseed.io.LineReader;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.stats.Shuffler;
import org.junit.jupiter.api.Test;
import static org.hamcrest.Matchers.*;
import static org.hamcrest.MatcherAssert.assertThat;

/**
 * Test sorting of blast hits.
 *
 * @author Bruce Parrello
 *
 */
public class TestHitSorts  {

    @Test
    public void testSorts() throws IOException {
        Shuffler<BlastHit> results0 = new Shuffler<BlastHit>(15);
        Map<String, String> qMap = new HashMap<String, String>();
        qMap.put("q1", "qtitle 1");
        qMap.put("q2", "qtitle 2");
        try (LineReader testStream = new LineReader(new File("data", "results.txt"))) {
            for (String line : testStream) {
                BlastHit result = new BlastHit(line, qMap, true, true);
                results0.add(result);
            }
        }
        // Scramble and sort by length.
        results0.shuffle(results0.size());
        results0.sort(new BlastHit.Longest());
        for (int i = 0; i < results0.size() - 1; i++)
            assertThat(Integer.toString(i), results0.get(i).getAlignLen(), greaterThanOrEqualTo(results0.get(i+1).getAlignLen()));
        // Scramble and sort by bit score.
        results0.shuffle(results0.size());
        results0.sort(new BlastHit.ByBitScore());
        for (int i = 0; i < results0.size() - 1; i++)
            assertThat(Integer.toString(i), results0.get(i).getBitScore(), greaterThanOrEqualTo(results0.get(i+1).getBitScore()));
    }
}
