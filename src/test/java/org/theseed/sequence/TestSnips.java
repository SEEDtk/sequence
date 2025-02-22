/**
 *
 */
package org.theseed.sequence;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.sequence.clustal.ClustalPipeline;
import org.theseed.sequence.clustal.RealSnipItem;
import org.theseed.sequence.clustal.SnipColumn;
import org.theseed.sequence.clustal.SnipIterator;

import org.junit.jupiter.api.Test;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 * @author Bruce Parrello
 *
 */
public class TestSnips {

    @Test
    public void testSnipItems() {
        RealSnipItem is = new RealSnipItem("aa-", 100, 2);
        assertThat(is.getLen(), equalTo(2));
        assertThat(is.getChars(), equalTo("aa-"));
        Location loc1 = Location.create("contig1", 1000, 2000);
        Location loc2 = Location.create("contig2", 2000, 1000);
        assertThat(is.getLocString(loc1), equalTo("contig1_1100+2"));
        assertThat(is.getLocString(loc2), equalTo("contig2_1900-2"));
        is = new RealSnipItem("----", 100, 0);
        assertThat(is.getLocString(loc1), equalTo("contig1_1100+0"));
        assertThat(is.getLocString(loc2), equalTo("contig2_1900-0"));
    }

    @Test
    public void testSnipAnalysis() throws IOException {
        Genome gto = new Genome(new File("data", "1313.5684.gto"));
        Feature feat = gto.getFeature("fig|1313.5684.peg.2088");
        ExtendedProteinRegion region = new ExtendedProteinRegion(feat, 100);
        RealSnipItem is = new RealSnipItem("aaa--atg-caacta--ttg-a", 97, 15);
        String[] aaMap = is.getProteinMap(region);
        assertThat(aaMap, arrayContaining("upstream", "upstream", "upstream", "Methionine", "Methionine", "Methionine", "Methionine", "Methionine",
                "Serine", "Serine", "Serine", "Serine", "Threonine", "Threonine", "Threonine", "Isoleucine", "Isoleucine",
                "Isoleucine", "Isoleucine", "Isoleucine", "Glutamic Acid", "Glutamic Acid"));
        is = new RealSnipItem("taaataa---", 1137, 7);
        aaMap = is.getProteinMap(region);
        assertThat(aaMap, arrayContaining("Glutamine", "Phenylalanine", "Phenylalanine", "Phenylalanine", "Lysine", "Lysine",
                "Lysine", "stop", "stop", "stop"));
    }

    @Test
    public void testSnipIterator() throws IOException, InterruptedException {
        RegionList phesRegions = new RegionList();
        Genome gto = new Genome(new File("data", "1313.5684.gto"));
        Feature feat = gto.getFeature("fig|1313.5684.peg.2088");
        ExtendedProteinRegion region = new ExtendedProteinRegion(feat, 100);
        phesRegions.add(region);
        gto = new Genome(new File("data", "1313.7001.gto"));
        feat = gto.getFeature("fig|1313.7001.peg.897");
        region = new ExtendedProteinRegion(feat, 95);
        phesRegions.add(region);
        gto = new Genome(new File("data", "1313.7090.gto"));
        feat = gto.getFeature("fig|1313.7090.peg.2070");
        region = new ExtendedProteinRegion(feat, 105);
        phesRegions.add(region);
        gto = new Genome(new File("data", "1313.5593.gto"));
        feat = gto.getFeature("fig|1313.5593.peg.1558");
        region = new ExtendedProteinRegion(feat, 98);
        phesRegions.add(region);
        gto = new Genome(new File("data", "360106.5.gto"));
        feat = gto.getFeature("fig|360106.5.peg.1206");
        region = new ExtendedProteinRegion(feat, 10);
        phesRegions.add(region);
        File tempFile = new File("data/temp", "temp.fa");
        phesRegions.save(tempFile);
        ClustalPipeline pipeline = new ClustalPipeline(tempFile);
        List<Sequence> alignment = pipeline.run();
        List<String> genomes = Arrays.asList("1313.5684", "1313.7001", "1313.7090", "1313.5593", "1313.6795");
        Set<String> wildSet = new TreeSet<String>();
        wildSet.add("1313.5684");
        wildSet.add("360106.5");
        SnipIterator.Run snipRun = new SnipIterator.Run(phesRegions, alignment, wildSet, genomes);
        for (SnipColumn snipCol : snipRun) {
            assertThat(snipCol.getRows(), equalTo(5));
            // Verify that the snips are properly located.
            String baseText = snipCol.getSnip(0);
            assertThat(snipCol.getWidth(), equalTo(baseText.length()));
            for (int row = 1; row < genomes.size(); row++) {
                String newText = snipCol.getSnip(row);
                if (! snipCol.isSignificant(row)) {
                    assertThat(newText, emptyString());
                } else {
                    assertThat(newText.length(), equalTo(baseText.length()));
                    assertThat(Feature.genomeOf(snipCol.getFid(row)), equalTo(genomes.get(row)));
                    int offset = snipCol.getOffset(row);
                    String compacted = StringUtils.remove(newText, '-');
                    if (! compacted.isEmpty())
                        assertThat(phesRegions.get(row).getSequence().substring(offset, offset + compacted.length()), equalTo(compacted));
                }
            }
        }
    }
}
