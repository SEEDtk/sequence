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
import org.theseed.io.MarkerFile;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.Sequence;

/**
 * This object represents a blast database stored in a file.  It can be used to reference an
 * existing blast database, or create one from a FASTA file or a genome object.
 *
 * @author Bruce Parrello
 *
 */
public class BlastDB {

    /**
     * sequence type of the blast database
     */
    public enum Type {
        DNA("nucl", "blastn", "tblastn"),
        PROTEIN("prot", "blastx", "blastp");

        /** parameter for dbtype option */
        private String dbParm;
        /** program to blast from DNA */
        private String dnaProg;
        /** program to blast from protein */
        private String protProg;

        private Type(String parm, String fromDna, String fromProt) {
            this.dbParm = parm;
            this.dnaProg = fromDna;
            this.protProg = fromProt;
        }

        /** return the dbType parameter to create a database of this type */
        public String getDbParm() {
            return this.dbParm;
        }

        /** return the program to use for a DNA subject sequence */
        public String getDnaBlaster() {
            return new File(BLAST_PATH, this.dnaProg).toString();
        }
        /** return the program to use for a protein subject sequence */
        public String getProteinBlaster() {
            return new File(BLAST_PATH, this.protProg).toString();
        }
    }

    // FIELDS

    /** log stream */
    protected static Logger log = LoggerFactory.getLogger(BlastDB.class);
    /** path to BLAST+ software */
    protected static String BLAST_PATH = System.getenv("BLAST_PATH");
    /** file name of the blast database */
    private File dbFile;
    /** type of the blast database */
    private Type type;
    /** genetic code for DNA database */
    private int geneticCode;

    /** genetic code for protein database */
    public static final int PROTEIN = 0;

    /**
     * Create a blast database from a FASTA file.  The database support files will be put in
     * the same directory.
     *
     * @param fastaFile		the FASTA file to be used
     * @param geneticCode	genetic code, or 0 for a protein database
     *
     * @throws IOException
     * @throws InterruptedException
     */
    public BlastDB(File fastaFile, int geneticCode) throws IOException, InterruptedException {
        this.dbFile = fastaFile;
        this.type = (geneticCode == 0 ? BlastDB.Type.PROTEIN : BlastDB.Type.DNA);
        this.geneticCode = geneticCode;
        this.createDb();
    }

    /**
     * Create a blast database for a FASTA file.  This should only be used for blast databases
     * that already exist.
     *
     * @param fastaFile		the FASTA file to be used
     */
    public BlastDB(File fastaFile) throws IOException, InterruptedException {
        this.dbFile = fastaFile;
        File testFile = new File(fastaFile.getPath() + ".psq");
        if (testFile.exists()) {
            log.info("Using protein BLAST database {}.", fastaFile);
            this.type = Type.PROTEIN;
        } else {
            testFile = new File(fastaFile.getPath() + ".nsq");
            if (testFile.exists()) {
                log.info("Using DNA BLAST database {}.", fastaFile);
                this.type = Type.DNA;
                this.geneticCode = MarkerFile.readInt(new File(fastaFile.getPath() + ".gc"));
            } else {
                throw new FileNotFoundException("No support files found for " + fastaFile + ".");
            }
        }
        if (testFile.lastModified() < fastaFile.lastModified()) {
            log.info("BLAST database expired-- recreating.");
            this.createDb();
        }
    }


    /**
     * Build the database support files using makeblastdb.
     *
     * @throws InterruptedException
     * @throws IOException
     */
    private void createDb() throws IOException, InterruptedException {
        // Create the command.
        ProcessBuilder makeProcess = new ProcessBuilder(new File(BLAST_PATH, "makeblastdb").toString(),
                "-in", this.dbFile.getAbsoluteFile().toString(), "-dbtype", this.type.getDbParm());
        // Merge STDERR into STDOUT, then start the process.
        makeProcess.redirectErrorStream(true);
        log.info("Running blast database make for {} and type {}.", this.dbFile, this.type);
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
            // Store the genetic code if this is a DNA database.
            if (this.type == BlastDB.Type.DNA)
                MarkerFile.write(new File(this.dbFile.getPath() + ".gc"), this.geneticCode);
            log.info("Blast database created.");
        }
    }

    /**
     * @return the database type
     */
    public Type getType() {
        return type;
    }

    /**
     * @return the genetic code (0 for a protein database)
     */
    public int getGeneticCode() {
        return this.geneticCode;
    }

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
     * @throws CloneNotSupportedException
     */
    public List<BlastHit> blastProteins(List<Sequence> proteins, BlastParms parms)
            throws CloneNotSupportedException, IOException, InterruptedException {
        String blastProgram = this.type.getProteinBlaster();
        BlastParms parms2 = parms.clone();
        if (this.geneticCode != 0) parms2.db_gencode(this.geneticCode);
        List<BlastHit> retVal = this.runBlast(blastProgram, proteins, parms2, Type.PROTEIN);
        return retVal;
    }

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
     * @throws CloneNotSupportedException
     */
    public List<BlastHit> blastDna(List<Sequence> contigs, BlastParms parms)
            throws CloneNotSupportedException, IOException, InterruptedException {
        String blastProgram = this.type.getDnaBlaster();
        List<BlastHit> retVal = this.runBlast(blastProgram, contigs, parms.clone(), Type.DNA);
        return retVal;
    }

    /**
     * Run BLAST against the specified sequences with the specified parameters.
     *
     * @param blastProgram	name of the BLAST program to run
     * @param seqs			list of query sequences
     * @param myParms			current blast parameters
     * @param qType			query sequence type
     *
     * @return a list of BlastHit objects representing the hits
     *
     * @throws CloneNotSupportedException
     * @throws IOException
     * @throws InterruptedException
     */
    private List<BlastHit> runBlast(String blastProgram, List<Sequence> seqs, BlastParms myParms, Type qType)
            throws CloneNotSupportedException, IOException, InterruptedException {
        List<BlastHit> retVal = new ArrayList<BlastHit>();
        // Set up the parameters and form the command.
        myParms.set("-outfmt", BlastHit.OUT_FORMAT);
        myParms.set("-db", this.dbFile.getAbsolutePath());
        List<String> command = new ArrayList<String>();
        command.add(blastProgram);
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
                BlastHit result = new BlastHit(line, qMap, qType, this.type);
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

}
