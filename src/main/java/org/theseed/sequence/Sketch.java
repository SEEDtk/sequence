/**
 *
 */
package org.theseed.sequence;

import java.io.Serializable;

/**
 * This object represents a kmer signature.  Each signature is associated with an identifying name; the name serves
 * as a key for accessing the named object-- it could be a protein family, a genome, a feature, or a role.
 *
 * @author Bruce Parrello
 *
 */
public class Sketch implements Serializable {

    // FIELDS

    private static final long serialVersionUID = -1880805874119392439L;

    /** signature for this entry's sequence */
    private int[] signature;
    /** target object associated with the sequence */
    private String name;
    /** numeric ID of this sketch */
    private int idNum;

    /**
     * Construct a hash entry.
     *
     * @param seq		kmer object to use for comparison
     * @param string	string associated with the sequence that generated the kmers
     * @param idNum		unique ID for this sketch
     */
    public Sketch(int[] sk, String string, int idNum) {
        this.signature = sk;
        this.name = string;
        this.idNum = idNum;
    }

    /**
     * @return the distance between another sequence and this entry's sequence
     *
     * @param otherSig	signature of the other sequence
     */
    public double distance(int[] otherSig) {
        return SequenceKmers.signatureDistance(this.signature, otherSig);
    }

    /**
     * @return the sequence target
     */
    public String getName() {
        return name;
    }

    /**
     * @return the ID number of this sketch
     */
    public int getIdNum() {
        return idNum;
    }

    /**
     * @return the signature for this sketch
     */
    public int[] getSignature() {
        return this.signature;
    }

}
