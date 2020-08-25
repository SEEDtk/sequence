/**
 *
 */
package org.theseed.sequence;

import junit.framework.TestCase;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.SortedSet;

import org.theseed.counters.CountMap;
import org.theseed.sequence.hash.Bucket;
import org.theseed.sequence.hash.LSHDiskSeqHash;
import org.theseed.sequence.hash.Sketch;

public class DiskHashTest extends TestCase {


    public void testLSHDiskSeqHash() throws IOException {
        ProteinKmers.setKmerSize(8);
        LSHDiskSeqHash.setCacheLimit(100);
        String[] prots = new String[] { SeqTest.p1, SeqTest.p2, SeqTest.p3, SeqTest.p4, SeqTest.p5,
                "MLALDDLLKNPYETFRIADELRRELVGDTVTYVVNRNINFTDICINDCKFCSFRNRKKYLLSLDEIKQKVEEAVEFGCTELCIQGGLLPDADLDFYLSILQAVRDV",
                "MNAFTRMNGHYSEQTDLKAVLAKASPDVRTILEKALAGEELTQPEAIVLFETEGADYSAVLKTADAVRQQRCGDEASFIVTRNINFTNVCYMGCSFCNFSVAKDDA",
                "MTIIGDSGRQIGPLTDAEAVELLGADGADLEDLCARADALRRDLVGDTLTFVVNRNLDTERVGAGTDESRERVRALVAEAAGLGATEICMQGPLPAGAPRDGYLDL"
        };
        File hashFile = new File("data", "hash");
        LSHDiskSeqHash hash = LSHDiskSeqHash.create(100, 30, 100, 8, hashFile);
        SeqTest.createTestingHash(prots, hash);
        assertThat(hash.getKmerSize(), equalTo(8));
        int[] foundCount = new int[prots.length];
        for (int i = 0; i < prots.length; i++) {
            ProteinKmers seqKmers = new ProteinKmers(prots[i]);
            SortedSet<Bucket.Result> found = hash.getClosest(seqKmers, 5, 0.9);
            String type = String.format("p%d", i);
            for (Bucket.Result result : found) {
                assertThat(result.getDistance(), lessThanOrEqualTo(0.9));
                assertThat(result.getTarget(), equalTo(type));
            }
            foundCount[i] = found.size();
        }
        CountMap<String> skCounts = new CountMap<String>();
        for (Sketch sk : hash.sketches()) {
            skCounts.count(sk.getName());
        }
        hash.save();
        hash = LSHDiskSeqHash.load(hashFile);
        assertThat(hash.getKmerSize(), equalTo(8));
        for (int i = 0; i < prots.length; i++) {
            ProteinKmers seqKmers = new ProteinKmers(prots[i]);
            SortedSet<Bucket.Result> found = hash.getClosest(seqKmers, 5, 0.9);
            String type = String.format("p%d", i);
            for (Bucket.Result result : found) {
                assertThat(result.getDistance(), lessThanOrEqualTo(0.9));
                assertThat(result.getTarget(), equalTo(type));
            }
            assertThat(found.size(), equalTo(foundCount[i]));
        }
        CountMap<String> skCounts2 = new CountMap<String>();
        for (Sketch sk : hash.sketches()) {
            skCounts2.count(sk.getName());
        }
        for (String key : skCounts.keys()) {
            assertThat(key, skCounts2.getCount(key), equalTo(skCounts.getCount(key)));
        }
    }
}
