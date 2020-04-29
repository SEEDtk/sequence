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

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.LineReader;
import org.theseed.sequence.DnaStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.ProteinStream;
import org.theseed.sequence.Sequence;
import org.theseed.sequence.SequenceStream;

/**
 * This object represents a blast database stored in a file.  It can be used to reference an
 * existing blast database, or create one from a FASTA file or a genome object.
 *
 * @author Bruce Parrello
 *
 */
public abstract class BlastDB {

    // FIELDS

    /** log stream */
    protected static Logger log = LoggerFactory.getLogger(BlastDB.class);
    /** path to BLAST+ software */
    protected static String BLAST_PATH = System.getenv("BLAST_PATH");
    /** file name of the blast database */
    private File dbFile;
    /** type of last BLAST run */
    private String blastType;
    /** parm string of last BLAST run */
    private String blastParms;

    /**
     * Create a blast database for a FASTA file.  This should only be used for blast databases
     * that already exist.
     *
     * @param fastaFile		the FASTA file to be used
     */
    public static BlastDB load(File fastaFile) throws IOException, InterruptedException {
        BlastDB retVal = null;
        File testFile = new File(fastaFile.getPath() + ".psq");
        if (testFile.exists()) {
            retVal = new ProteinBlastDB(fastaFile);
        } else {
            testFile = new File(fastaFile.getPath() + ".nsq");
            if (testFile.exists()) {
                retVal = new DnaBlastDB(fastaFile);
            } else {
                throw new FileNotFoundException("No support files found for " + fastaFile + ".");
            }
        }
        if (testFile.lastModified() < fastaFile.lastModified()) {
            log.info("BLAST database expired-- recreating.");
            retVal.createDb();
        }
        return retVal;
    }


    /**
     * Build the database support files using makeblastdb.
     *
     * @throws InterruptedException
     * @throws IOException
     */
    protected void createDb() throws IOException, InterruptedException {
        // Create the command.
        ProcessBuilder makeProcess = new ProcessBuilder(new File(BLAST_PATH, "makeblastdb").getPath(),
                "-in", this.dbFile.getAbsoluteFile().toString(), "-dbtype", this.getDbParm());
        // Merge STDERR into STDOUT, then start the process.
        makeProcess.redirectErrorStream(true);
        log.info("Running blast database make for {} and type {}.", this.dbFile, this.getDbParm());
        Process process = makeProcess.start();
        // Get a reader for the process output.
        try (LineReader reader = new LineReader(process.getInputStream())) {
            for (String line : reader)
                log.debug("*   {}", line);
            int exitCode = process.waitFor();
            if (exitCode != 0) {
                String commandString = StringUtils.join(makeProcess.command(), " ");
                log.error("Failing makeBlastDB command with exit code {} was {}.", exitCode, commandString);
                throw new RuntimeException("Exit code " + exitCode + " from BlastDB creation.");
            }
            log.info("Blast database created.");
        }
    }

    /**
     * @return the database creation type parameter
     */
    protected abstract String getDbParm();

    /**
     * @return the database type
     */
    public abstract boolean isProtein();

    /**
     * BLAST protein sequences against this database.
     *
     * @param proteins	list of protein sequences
     * @param parms		BLAST parameters
     *
     * @return a list of BLAST matches
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public abstract List<BlastHit> blast(ProteinStream proteins, BlastParms parms)
            throws IOException, InterruptedException;

    /**
     * BLAST DNA sequences against this database.
     *
     * @param contigs	list of DNA sequences
     * @param parms		BLAST parameters
     *
     * @return a list of BLAST matches
     *
     * @throws InterruptedException
     * @throws IOException
     */
    public abstract List<BlastHit> blast(DnaStream contigs, BlastParms parms)
            throws IOException, InterruptedException;

    /**
     * Run BLAST against the specified sequences with the specified parameters.
     *
     * @param blastProgram	name of the BLAST program to run
     * @param seqs			list of query sequences
     * @param myParms		current blast parameters
     * @param qProt			TRUE if the query sequence is protein
     *
     * @return a list of BlastHit objects representing the hits
     *
     * @throws IOException
     * @throws InterruptedException
     */
    protected List<BlastHit> runBlast(String blastProgram, SequenceStream seqs, BlastParms myParms)
            throws IOException, InterruptedException {
        List<BlastHit> retVal = new ArrayList<BlastHit>();
        // Save the command specs.
        this.blastType = blastProgram;
        this.blastParms = myParms.toString();
        // Set up the parameters and form the command.
        myParms.set("-outfmt", BlastHit.OUT_FORMAT);
        myParms.set("-db", this.dbFile.getAbsolutePath());
        List<String> command = new ArrayList<String>();
        command.add(new File(BLAST_PATH, blastProgram).getPath());
        command.addAll(myParms.get());
        ProcessBuilder blastCommand = new ProcessBuilder(command);
        // Start the BLAST.
        log.info("Running BLAST command {}.", blastProgram);
        Process blastProcess = blastCommand.start();
        // Write the sequences to the BLAST process and save a map of sequence IDs to comments.
        Map<String, String> qMap = new HashMap<String, String>();
        try (FastaOutputStream queryStream = new FastaOutputStream(blastProcess.getOutputStream())) {
            for (Sequence seq : seqs) {
                // Store this sequence's comment in the map.
                qMap.put(seq.getLabel(), seq.getComment());
                // Write the sequence to the BLAST process.
                queryStream.write(seq);
            }
        }
        log.info("{} query sequences submitted.", qMap.size());
        // Read back the results.
        try (LineReader blastReader = new LineReader(blastProcess.getInputStream())) {
            int count = 0;
            for (String line : blastReader) {
                BlastHit result = new BlastHit(line, qMap, seqs.isProtein(), this.isProtein());
                count++;
                if (myParms.acceptable(result))
                    retVal.add(result);
            }
            log.info("{} results from BLAST, {} returned.", count, retVal.size());
        }
        // Clean up the process.
        int exitCode = blastProcess.waitFor();
        if (exitCode != 0) {
            // We have an error. Output the error log.
            log.error("Output from BLAST error log follows.");
            try (LineReader logReader = new LineReader(blastProcess.getErrorStream())) {
                for (String line : logReader)
                    log.error("   {}", line);
            }
            throw new RuntimeException("Blast process failed with exit code " + exitCode);
        }
        return retVal;
    }

    /**
     * @return the file name of this database
     */
    public File getFile() {
        return this.dbFile;
    }

    /**
     * Specify the FASTA file for this database.
     *
     * @param file	file containing the main FASTA sequences
     */
    protected void setFile(File file) {
        this.dbFile = file;
    }

    /**
     * Specify that this FASTA database be deleted on exit.
     */
    public void deleteOnExit() {
        this.dbFile.deleteOnExit();
        cleanOnExit();
    }


    /**
     * Specify that this FASTA database's support files be deleted on exit.
     */
    public void cleanOnExit() {
        String[] suffixes = this.getSuffixes();
        String path = this.dbFile.getPath();
        for (String suffix : suffixes) {
            File serviceFile = new File(path + suffix);
            serviceFile.deleteOnExit();
        }
    }


    /**
     * @return a list of the suffixes
     */
    protected abstract String[] getSuffixes();


    /**
     * @return the type of the last BLAST run
     */
    public String getBlastType() {
        return blastType;
    }


    /**
     * @return the parameter string for the last BLAST run
     */
    public String getBlastParms() {
        return blastParms;
    }

}
