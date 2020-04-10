/**
 *
 */
package org.theseed.sequence;

import junit.framework.TestCase;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import org.theseed.io.TabbedLineReader;

/**
 * @author Bruce Parrello
 *
 */
public class BucketTest extends TestCase {

    private String p1 = "MDIQITHQVTEFDKEELLAGLRSYNAQFVDFSKNGQLGVYCRNESGEMVGGLIADRKGPWLCIDYLWVSESARNCGLGSKLMAMAEKEGLRKGCAHGLVD";
    private String p2 = "MSKDEISYQILYRYSLEKLYSTLTRRVDNVLSFALIFLGVGVTINVGSPFILGPGIVGIAILKRVLRFGTRSAQADRQSRAWLKLFNTQHRFPSDKTLFL";
    private String p3 = "MELQLMLNHFFERVRKDANFNAFLIDLEYNNIAYYIYFVATGNVKIITHAGHFISIKSNRKLIKVNSTPNTQLIKLTSDKHFSGEHSYEKYCTDLATAGV";
    private String p4 = "MTNITLSTQHYRIHRSDVEPVKEKTTEKDIFAKSITAVRNSFISLSTSLSDRFSLHQQTDIPTTHFHRGSASEGRAVLTSKTVKDFMLQKLNSLDIKGNA";
    private String p5 = "MSKDPAYARQTCEAILSAVYSNNKDHCCKLLISKGVSITPFLKEIGEAAQNAGLPGEIKNGVFTPGGAGANPFVVPLIASASIKYPHMFINHNQQVSFKA";

    /**
     * test sketch operations
     */
    public void testSketch() {
        String p12 = p1 + p2;
        ProteinKmers kmers1 = new ProteinKmers(p12);
        int[] signature = kmers1.hashSet(50);
        Sketch sketch1 = new Sketch(p12, "p12", 50);
        assertTrue(Arrays.equals(signature, sketch1.getSignature()));
        assertThat(sketch1.getName(), equalTo("p12"));
        String p21 = p2 + p1;
        Sketch sketch2 = new Sketch(p21, "p21", 50);
        double dist = sketch2.distance(sketch1);
        double dist2 = sketch2.distance(signature);
        assertThat(dist, equalTo(dist2));
        assertThat(dist, closeTo(0.148, 0.001));
        Sketch sketch1a = new Sketch(signature, "p12a");
        assertTrue(sketch1a.isSameSignature(sketch1));
    }

    /**
     * test bucket search
     */
    public void testSearch() {
        Bucket testBucket = new Bucket();
        Sketch sketch1 = new Sketch(p1, "g1", 360);
        Sketch sketch2 = new Sketch(p2, "g1", 360);
        Sketch sketch3 = new Sketch(p3, "p3", 360);
        Sketch sketch4 = new Sketch(p4, "p4", 360);
        Sketch sketch5 = new Sketch(p5, "p5", 360);
        testBucket.add(sketch1);
        assertThat(testBucket.size(), equalTo(1));
        testBucket.add(sketch2);
        testBucket.add(sketch3);
        testBucket.add(sketch4);
        testBucket.add(sketch5);
        List<Sketch> found = testBucket.search("g1");
        assertThat(found.size(), equalTo(2));
        for (Sketch fSketch : found)
            assertThat(fSketch.getName(), equalTo("g1"));
        found = testBucket.search("p3");
        assertThat(found.size(), equalTo(1));
        assertTrue(found.get(0) == sketch3);
        assertThat(testBucket.get(2).getName(), equalTo("p3"));
        List<Sketch> subList = testBucket.after(1);
        assertThat(subList.size(), equalTo(3));
        assertThat(subList.get(2).getName(), equalTo("p5"));
    }

    /**
     * test bucket save and load
     *
     * @throws IOException
     * @throws ClassNotFoundException
     */
    public void testBucketFiles() throws IOException, ClassNotFoundException {
        // The plan is to read a file of proteins and store them as sketches with their MD5 identifier.  We will
        // write out the sketches as a bucket file, then read them back in and verify that all the sketches are
        // still there.  This will be the case if every MD5 identifier is found, and the object found has an
        // identical sketch.
        // This will be the original bucket.
        Bucket original = new Bucket();
        assertThat(original.size(), equalTo(0));
        // Read the test file and fill the bucket.
        try (TabbedLineReader inStream = new TabbedLineReader(new File("src/test", "prot.fams.tbl"))) {
            // Column 2 is the MD5, 3 is the protein.
            for (TabbedLineReader.Line line : inStream) {
                String md5 = line.get(2);
                String prot = line.get(3);
                ProteinKmers kmers = new ProteinKmers(prot);
                Sketch sketch = new Sketch(kmers.hashSet(360), md5);
                original.add(sketch);
            }
        }
        // Now we write out the bucket.
        File saveFile = new File("src/test", "bucket.ser");
        original.save(saveFile);
        // Read in a copy.
        Bucket saved = Bucket.load(saveFile);
        assertThat(saved.size(), equalTo(original.size()));
        // Veryify the sketches are the same.
        for (Sketch oSketch : original) {
            List<Sketch> found = saved.search(oSketch.getName());
            assertThat(found.size(), equalTo(1));
            Sketch fSketch = found.get(0);
            assertThat(fSketch.getName(), equalTo(oSketch.getName()));
            assertTrue(fSketch.isSameSignature(oSketch));
        }
        // Verify that we can add a new sketch.

    }

}
