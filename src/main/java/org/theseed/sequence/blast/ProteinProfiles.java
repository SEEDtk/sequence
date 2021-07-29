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

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.Role;
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
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ProteinProfiles.class);
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
        setup(profileDir, null);
    }

    /**
     * Construct a new protein profile.
     *
     * @param profileDir	profile directory
     * @param roleFilter	set of roles to include (all other roles will be excluded)
     *
     * @throws IOException
     */
    public ProteinProfiles(File profileDir, Set<String> roleFilter) throws IOException {
        setup(profileDir, roleFilter);
    }
    /**
    /**
     * Initialize a new protein profile.
     *
     * @param profileDir	profile directory
     * @param roleFilter	set of roles to include (all other roles will be excluded)
     *
     * @throws IOException
     * @throws FileNotFoundException
     */
    private void setup(File profileDir, Set<String> roleFilter) throws IOException, FileNotFoundException {
        this.profileDir = profileDir;
        this.roleMap = RoleMap.load(new File(profileDir, "_roles.tbl"));
        // If we are filtering, we must remove the roles we are not using.
        if (roleFilter != null) {
            // This is a two-pass process so that we are not modifying the key set while streaming through it.
            Set<String> badRoles = this.roleMap.keySet().stream().filter(x -> ! roleFilter.contains(x)).collect(Collectors.toSet());
            badRoles.stream().forEach(x -> this.roleMap.remove(x));
            log.info("{} roles kept after filtering.", this.roleMap.size());
        }
        this.clusterMap = new HashMap<String, String>();
        try (TabbedLineReader clusterStream = new TabbedLineReader(new File(profileDir, "_map.tbl"), 2)) {
            int kept = 0;
            int discarded = 0;
            for (TabbedLineReader.Line line : clusterStream) {
                String roleName = line.get(1);
                String roleId = null;
                String[] roleNames = Feature.rolesOfFunction(roleName);
                for (String roleName0 : roleNames) {
                    Role role = this.roleMap.getByName(roleName0);
                    if (role != null)
                        roleId = role.getId();
                }
                if (roleId != null) {
                    this.clusterMap.put(line.get(0), roleName);
                    kept++;
                } else
                    discarded++;
            }
            log.info("{} profiles kept, {} removed by filter.", kept, discarded);
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
