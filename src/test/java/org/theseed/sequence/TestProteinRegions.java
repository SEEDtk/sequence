/**
 *
 */
package org.theseed.sequence;

import org.junit.jupiter.api.Test;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.locations.Location;
import org.theseed.proteins.Function;
import org.theseed.proteins.FunctionMap;

/**
 * @author Bruce Parrello
 *
 */
public class TestProteinRegions {

    /**
     * test ExtendedProteinRegion
     *
     * @throws IOException
     */
    @Test
    public void testProteinRegions() throws IOException {
        Genome gto = new Genome(new File("data", "1313.7001.gto"));
        Feature feat = gto.getFeature("fig|1313.7001.peg.326");
        ExtendedProteinRegion region = new ExtendedProteinRegion(feat, 733);
        assertThat(region.getComment(), equalTo("hypothetical protein"));
        assertThat(region.getFeature(), equalTo(feat));
        Location loc = region.getFullLocation();
        assertThat(loc.getLeft(), equalTo(69444));
        assertThat(loc.getRight(), equalTo(70791));
        assertThat(loc.getDir(), equalTo('+'));
        assertThat(loc.getContigId(), equalTo("1313.7001.con.0011"));
        assertThat(gto.getDna(loc), equalTo(region.getSequence()));
        assertThat(region.getProteinTranslation(), equalTo(feat.getProteinTranslation()));
        feat = gto.getFeature("fig|1313.7001.peg.325");
        region = new ExtendedProteinRegion(feat, 500);
        loc = region.getFullLocation();
        assertThat(region.getFeature(), equalTo(feat));
        assertThat(loc.getLeft(), equalTo(69473));
        assertThat(loc.getLength(), equalTo(feat.getLocation().getLength() + 500));
        assertThat(loc.getDir(), equalTo('-'));
        assertThat(gto.getDna(loc), equalTo(region.getSequence()));
        assertThat(region.getProteinTranslation(), equalTo(feat.getProteinTranslation()));
        feat = gto.getFeature("fig|1313.7001.peg.897");
        region = new ExtendedProteinRegion(feat, 318);
        gto = new Genome(new File("data", "1313.7090.gto"));
        feat = gto.getFeature("fig|1313.7090.peg.2070");
        ExtendedProteinRegion region2 = new ExtendedProteinRegion(feat, 308);
        double dist = region.getDistance(region2);
        assertThat(dist, closeTo(0.1821, 0.0001));
    }

    /**
     * test ExtendedProteinRegion lists
     *
     * @throws IOException
     */
    public void testRegionLists() throws IOException {
        Genome gto = new Genome(new File("data", "1313.7001.gto"));
        RegionList regions = new RegionList(gto, 500);
        assertThat(regions.size(), equalTo(gto.getPegs().size()));
        for (ExtendedProteinRegion region : regions) {
            Feature feat = region.getFeature();
            assertThat(feat.getParent(), equalTo(gto));
            assertThat(feat.getProteinLength(), greaterThan(0));
            Location rLoc = region.getFullLocation();
            Location fLoc = feat.getLocation();
            int upstream = rLoc.getLength() - fLoc.getLength();
            assertThat(feat.getId(), upstream, lessThanOrEqualTo(500));
            assertThat(feat.getId(), upstream, greaterThanOrEqualTo(0));
            if (feat.getId().contentEquals("fig|1313.7001.peg.896"))
                assertThat(upstream, equalTo(12));
            else if (feat.getId().contentEquals("fig|1313.7001.peg.326"))
                assertThat(upstream, equalTo(500));
            else if (feat.getId().contentEquals("fig|1313.7001.peg.4"))
                assertThat(upstream, equalTo(0));
            else if (feat.getId().contentEquals("fig|1313.7001.peg.1215"))
                assertThat(upstream, equalTo(2));
            assertThat(regions.get(feat.getId()), sameInstance(region));
        }
        File tempFile = new File("data/temp", "regionList.fa");
        regions.save(tempFile);
        try (FastaInputStream inSeqs = new FastaInputStream(tempFile)) {
            int i = 0;
            for (Sequence seq : inSeqs) {
                ExtendedProteinRegion region = regions.get(i);
                assertThat(seq.getLabel(), equalTo(region.getLabel()));
                assertThat(seq.getComment(), equalTo(region.getComment()));
                assertThat(seq.getSequence(), equalTo(region.getSequence()));
                i++;
            }
        }
    }

    /**
     * test the function maps
     */
    public void testFunctionMap() throws IOException {
        Genome gto = new Genome(new File("data", "1313.7001.gto"));
        FunctionMap funMap = new FunctionMap();
        Map<String, RegionList> rMap = RegionList.createMap(funMap, gto, 500);
        for (Map.Entry<String, RegionList> rEntry : rMap.entrySet()) {
            String funId = rEntry.getKey();
            for (ExtendedProteinRegion region : rEntry.getValue()) {
                Function function = funMap.getByName(region.getFeature().getPegFunction());
                assertThat(function.getId(), equalTo(funId));
            }
        }
        RegionList regions = new RegionList(gto, 500);
        // Verify that each region in the genome is the closest one to its own function.
        for (ExtendedProteinRegion region : regions) {
            Function fun = funMap.getByName(region.getFeature().getPegFunction());
            RegionList funRegions = rMap.get(fun.getId());
            ExtendedProteinRegion closest = funRegions.getClosest(region, 0.1);
            assertThat(closest.getFeature(), equalTo(region.getFeature()));
        }
        // Verify that the region list sort works.
        List<RegionList> sorted = rMap.values().stream().sorted().collect(Collectors.toList());
        assertThat(sorted.size(), equalTo(rMap.size()));
        for (int i = 1; i < sorted.size(); i++) {
            Location loc0 = sorted.get(i - 1).get(0).getFullLocation();
            Location loc1 = sorted.get(i).get(0).getFullLocation();
            String label = loc0.toString() + " <-> " + loc1.toString();
            assertThat(label, loc0, lessThanOrEqualTo(loc1));
            assertThat(label, loc0.getContigId(), lessThanOrEqualTo(loc1.getContigId()));
        }
    }

}
