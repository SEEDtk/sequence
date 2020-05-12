/**
 *
 */
package org.theseed.sequence.blast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.theseed.io.TabbedLineReader;

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
    /** protein ID to description map */
    private Map<String, String> qMap;
    /** directory containing the profiles */
    private File profileDir;

    /**
     * Construct a new protein profile.
     *
     * @param profileDir	profile directory
     *
     * @throws IOException
     */
    public ProteinProfiles(File profileDir) throws IOException {
        this.profileDir = profileDir;
        this.qMap = new HashMap<String, String>();
        // Open the map file.
        try (TabbedLineReader mapStream = new TabbedLineReader(new File(profileDir, "_map.tbl"), 2)) {
            for (TabbedLineReader.Line line : mapStream) {
                // Put the role information in the map.
                String role = line.get(0);
                String description = line.get(1);
                this.qMap.put(role, description);
                // Verify that the role file exists.
                File roleFile = this.getFile(role);
                if (! roleFile.canRead())
                    throw new FileNotFoundException("Role " + role + " found in profile map but not in directory.");
            }
        }
    }

    /**
     * @return the file containing the profile for the specified role ID
     *
     * @param role	role whose profile is desired
     */
    private File getFile(String role) {
        return new File(this.profileDir, role + ".smp");
    }

    /**
     * @return a map listing the profile hits against each subject sequence in the specified DNA blast database
     *
     * @param dnaDB		DNA blast database to process against the profiles
     * @param parms		BLAST parameters to use
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public Map<String, List<BlastHit>> profile(DnaBlastDB profiler, BlastParms parms) throws IOException, InterruptedException {
        Map<String, List<BlastHit>> retVal = new HashMap<String, List<BlastHit>>();
        for (String role : this.qMap.keySet()) {
            // Blast this role.
            List<BlastHit> results = profiler.psiBlast(this.getFile(role), parms, this.qMap);
            // Put the results in the result map.
            for (BlastHit result : results) {
                String seqId = result.getSubjectId();
                List<BlastHit> list = retVal.computeIfAbsent(seqId, x -> new ArrayList<BlastHit>(10));
                list.add(result);
            }
        }
        return retVal;
    }

}
