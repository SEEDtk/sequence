/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
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
    /** role ID to cluster list mapping */
    private Map<String, Set<String>> roleClusterMap;
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
        // Initialize the role-cluster map and the main cluster map.
        this.roleClusterMap = new HashMap<String, Set<String>>(this.roleMap.size());
        this.clusterMap = new HashMap<String, String>();
        // Read in the cluster map.
        try (TabbedLineReader clusterStream = new TabbedLineReader(new File(profileDir, "_map.tbl"), 2)) {
            int kept = 0;
            int discarded = 0;
            for (TabbedLineReader.Line line : clusterStream) {
                String roleName = line.get(1);
                String[] roleNames = Feature.rolesOfFunction(roleName);
                Set<String> roleIds = new TreeSet<String>();
                for (String roleName0 : roleNames) {
                    Role role = this.roleMap.getByName(roleName0);
                    if (role != null)
                        roleIds.add(role.getId());
                }
                if (! roleIds.isEmpty()) {
                    String clusterId = line.get(0);
                    this.clusterMap.put(clusterId, roleName);
                    kept++;
                    // Add this cluster to the relevant role maps.
                    for (String roleId : roleIds) {
                        Set<String> clusterIds = this.roleClusterMap.computeIfAbsent(roleId, x -> new TreeSet<String>());
                        clusterIds.add(clusterId);
                    }
                } else
                    discarded++;
            }
            log.info("{} profiles kept, {} removed by filter, {} roles with profiles.", kept, discarded, this.roleClusterMap.size());
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
     * @return the set of roles covered by profiles
     */
    public Set<String> getRoles() {
        return this.roleClusterMap.keySet();
    }

    /**
     * Blast all the profiles.
     *
     * @param blastDB	blast database to process against the profiles
     * @param parms		BLAST parameters to use
     *
     * @return a map listing the profile hits against each subject sequence in the specified blast database
     *
     */
    public Map<String, List<BlastHit>> profile(BlastDB blastDB, BlastParms parms) {
        Set<String> clusterSet = this.clusterMap.keySet();
        return blastProfiles(clusterSet, blastDB, parms);
    }

    /**
     * Blast the profiles for the selected roles.
     *
     * @param blastDB	blast database to process against the profiles
     * @param parms		BLAST parameters to use
     *
     * @return a map listing the profile hits against each subject sequence in the specified blast database
     */
    public Map<String, Set<BlastHit>> profile(Collection<String> roles, BlastDB blastDB, BlastParms parms) {
        Map<String, Set<BlastHit>> retVal = new HashMap<String, Set<BlastHit>>(roles.size());
        // Create a map of cluster IDs to role IDs.
        Map<String, Set<String>> clusterRoles = new HashMap<String, Set<String>>(roles.size() * 5);
        for (String role : roles) {
            Set<String> roleClusters = this.roleClusterMap.get(role);
            if (roleClusters != null) {
                for (String cluster : roleClusters) {
                    Set<String> rolesForCluster = clusterRoles.computeIfAbsent(cluster, x -> new TreeSet<String>());
                    rolesForCluster.add(role);
                }
            }
        }
        // Now we have a map of the roles for each cluster.  This will be used to convert the output to a hash keyed on
        // role ID.  Get the hits.
        Map<String, List<BlastHit>> hitMap = this.blastProfiles(clusterRoles.keySet(), blastDB, parms);
        // Convert this map to a role-based map.  For each cluster, we add all the cluster's hits to each of the cluster's
        // roles.
        for (List<BlastHit> hitList : hitMap.values()) {
            for (BlastHit hit : hitList) {
                String cluster = hit.getQueryLoc().getContigId();
                Set<String> rolesFound = clusterRoles.get(cluster);
                for (String role : rolesFound) {
                    Set<BlastHit> hitSet = retVal.computeIfAbsent(role, x -> new HashSet<BlastHit>());
                    hitSet.add(hit);
                }
            }
        }
        // Return the role-to-hit map.
        return retVal;
    }

    /**
     * Blast the specified clusters against a blast database.
     *
     * @param clusterSet	set of cluster IDs for the clusters to blast
     * @param blastDB		target blast database
     * @param parms			blast parameters to use
     *
     * @return a map listing the profile hits against each subject sequence in the specified blast database
     */
    private Map<String, List<BlastHit>> blastProfiles(Set<String> clusterSet, BlastDB blastDB, BlastParms parms) {
        Map<String, List<BlastHit>> retVal = new HashMap<String, List<BlastHit>>();
        // We run the BLASTs in parallel.  The "map" blasts the cluster, "flatMap"
        // converts a stream of lists to a stream of blast hits, and "collect"
        // converts the stream back to a single list.
        List<BlastHit> results = clusterSet.parallelStream()
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
