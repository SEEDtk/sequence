/**
 *
 */
package org.theseed.sequence;

import junit.framework.TestCase;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.List;
import org.theseed.io.TabbedLineReader;

/**
 * @author Bruce Parrello
 *
 */
public class BucketTest extends TestCase {

    /**
     * test bucket search
     */
    public void testSearch() {
        Bucket testBucket = new Bucket();
        String p1 = "MDIQITHQVTEFDKEELLAGLRSYNAQFVDFSKNGQLGVYCRNESGEMVGGLIADRKGPWLCIDYLWVSESARNCGLGSKLMAMAEKEGLRKGCAHGLVD";
        String p2 = "MSKDEISYQILYRYSLEKLYSTLTRRVDNVLSFALIFLGVGVTINVGSPFILGPGIVGIAILKRVLRFGTRSAQADRQSRAWLKLFNTQHRFPSDKTLFL";
        String p3 = "MELQLMLNHFFERVRKDANFNAFLIDLEYNNIAYYIYFVATGNVKIITHAGHFISIKSNRKLIKVNSTPNTQLIKLTSDKHFSGEHSYEKYCTDLATAGV";
        String p4 = "MTNITLSTQHYRIHRSDVEPVKEKTTEKDIFAKSITAVRNSFISLSTSLSDRFSLHQQTDIPTTHFHRGSASEGRAVLTSKTVKDFMLQKLNSLDIKGNA";
        String p5 = "MSKDPAYARQTCEAILSAVYSNNKDHCCKLLISKGVSITPFLKEIGEAAQNAGLPGEIKNGVFTPGGAGANPFVVPLIASASIKYPHMFINHNQQVSFKA";
        ProteinKmers kmer1 = new ProteinKmers(p1);
        ProteinKmers kmer2 = new ProteinKmers(p2);
        ProteinKmers kmer3 = new ProteinKmers(p3);
        ProteinKmers kmer4 = new ProteinKmers(p4);
        ProteinKmers kmer5 = new ProteinKmers(p5);
        Sketch sketch1 = new Sketch(kmer1.hashSet(360), "g1");
        Sketch sketch2 = new Sketch(kmer2.hashSet(360), "g1");
        Sketch sketch3 = new Sketch(kmer3.hashSet(360), "p3");
        Sketch sketch4 = new Sketch(kmer4.hashSet(360), "p4");
        Sketch sketch5 = new Sketch(kmer5.hashSet(360), "p5");
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
    }

}
