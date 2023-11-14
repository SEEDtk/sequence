package org.theseed.proteins.kmer.hash;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import org.junit.jupiter.api.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.TabbedLineReader;
import org.theseed.sequence.ProteinKmers;

class ProteinKmerHashTest {

    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ProteinKmerHashTest.class);


    @Test
    void testProteinKmerHash() throws IOException {
        File protFile = new File("data", "proteins.tbl");
        ProteinKmerHashMap<String> kmerMap = new ProteinKmerHashMap<String>(8);
        Map<String, ProteinKmers> md5Check = new HashMap<String, ProteinKmers>();
        try (TabbedLineReader protStream = new TabbedLineReader(protFile)) {
            for (var line : protStream) {
                String md5 = line.get(0);
                String prot = line.get(1);
                String annotation = line.get(2);
                md5Check.put(md5, new ProteinKmers(prot, 8));
                kmerMap.addProtein(prot, annotation);
            }
        }
        File testFile = new File("data", "testProteins.tbl");
        try (TabbedLineReader testStream = new TabbedLineReader(testFile)) {
            for (var line : testStream) {
                String md5 = line.get(0);
                String prot = line.get(1);
                String annotation = line.get(2);
                // Compute the expected result.
                ProteinKmers newKmers = new ProteinKmers(prot, 8);
                int expHits = 0;
                String expMd5 = "<none>";
                for (var md5Entry : md5Check.entrySet()) {
                    int hits = newKmers.similarity(md5Entry.getValue());
                    if (hits > expHits) {
                        expHits = hits;
                        expMd5 = md5Entry.getKey();
                    }
                }
                var result = kmerMap.findClosest(prot);
                log.info("Results for protein {} with annotation {}.", md5, annotation);
                log.info("Hit count {}, similarity {}.", result.getSimCount(), result.getSimValue());
                assertThat(md5, result.getSimCount(), equalTo(expHits));
                if (! result.isEmpty()) {
                    log.info("Value is {}.", result.getValue());
                    assertThat(result.getMd5(), equalTo(expMd5));
                }
            }
        }
    }

}
