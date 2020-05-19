/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.RoleMap;

/**
 * A protein profile directory contains Psi-BLAST input files (all with suffix ".smp") plus a "_map.tbl" file.
 * The latter file contains all the role IDs and descriptions in the directory.  This class manages the directory
 * and provides a method to produce BLAST hits of all of the protein profiles against a BLAST database.
 *
 * @author Bruce Parrello
 *
 */
public class ProteinProfiles {

    // FIELDS
    /** role ID mapping to parse the functions */
    private RoleMap roleMap;
    /** cluster ID to function mapping */
    private Map<String, String> clusterMap;
    /** directory containing the profiles */
    private File profileDir;
    /** number of profiles hit in last run */
    private int hitCount;

    /**
     * Construct a new protein profile.
     *
     * @param profileDir	profile directory
     *
     * @throws IOException
     */
    public ProteinProfiles(File profileDir) throws IOException {
        this.profileDir = profileDir;
        this.roleMap = RoleMap.load(new File(profileDir, "_roles.tbl"));
        this.clusterMap = new HashMap<String, String>();
        try (TabbedLineReader clusterStream = new TabbedLineReader(new File(profileDir, "_map.tbl"), 2)) {
            for (TabbedLineReader.Line line : clusterStream)
                this.clusterMap.put(line.get(0), line.get(1));
        }
        // Verify that all the clusters exist.
        for (String cluster : this.clusterMap.keySet()) {
            File clusterFile = this.getFile(cluster);
            if (! clusterFile.canRead())
                throw new FileNotFoundException("Role " + cluster + " found in profile map but not in directory.");
        }
    }

    /**
     * @return the file containing the profile for the specified cluster
     *
     * @param cluster	ID of cluster whose profile is desired
     */
    private File getFile(String cluster) {
        return new File(this.profileDir, cluster + ".smp");
    }

    /**
     * @return a map listing the profile hits against each subject sequence in the specified blast database
     *
     * @param blastDB	blast database to process against the profiles
     * @param parms		BLAST parameters to use
     *
     */
    public Map<String, List<BlastHit>> profile(BlastDB blastDB, BlastParms parms) {
        Map<String, List<BlastHit>> retVal = new HashMap<String, List<BlastHit>>();
        // We run the BLASTs in parallel.  The "map" blasts the cluster, "flatMap"
        // converts a stream of lists to a stream of blast hits, and "collect"
        // converts the stream back to a single list.
        List<BlastHit> results = this.clusterMap.keySet().parallelStream()
                .map(x -> blastDB.psiBlast(this.getFile(x), parms, this.clusterMap))
                .flatMap(List::stream).collect(Collectors.toList());
        // This will track the profiles hit.
        Set<String> profilesHit = new HashSet<String>(results.size());
        // Put the results in the result map.
        for (BlastHit result : results) {
            String seqId = result.getSubjectId();
            profilesHit.add(result.getQueryId());
            List<BlastHit> list = retVal.computeIfAbsent(seqId, x -> new ArrayList<BlastHit>(10));
            list.add(result);
        }
        // Count the profiles hit.
        this.hitCount = profilesHit.size();
        return retVal;
    }

    /**
     * @return the number of profiles
     */
    public int size() {
        return this.roleMap.size();
    }

    /**
     * @return the number of profiles hit during the last run
     */
    public int getHitCount() {
        return hitCount;
    }

    /**
     * @return the role map
     */
    public RoleMap roleMap() {
        return this.roleMap;
    }


}
