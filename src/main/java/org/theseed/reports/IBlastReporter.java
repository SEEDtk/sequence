/**
 *
 */
package org.theseed.reports;

import org.theseed.sequence.blast.BlastHit;

/**
 * @author Bruce Parrello
 *
 */
public interface IBlastReporter {

    /**
     * Record a hit by the blast.
     *
     * @param hit	blast hit to record
     */
    void recordHit(BlastHit hit);

    /**
     * Write the blast report.
     *
     * @param title			title to display
     * @param blastInfo		description of the BLAST
     */
    void writeReport(String string, BlastInfo blastInfo);

    /**
     * @return the number of sequences hit
     *
     * This is the number of queries for sort-type QUERY, subjects for sort-type SUBJECT.
     */
    int getSequencesHit();
}
