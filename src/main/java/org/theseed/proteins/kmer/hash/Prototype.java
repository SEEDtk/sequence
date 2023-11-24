/**
 *
 */
package org.theseed.proteins.kmer.hash;

import java.util.Map;

/**
 * This is a lightweight object that contains a protein sequence and an associated annotation.  It is used
 * to hold a single annotation file record in memory.
 *
 * @author Bruce Parrello
 *
 */
public class Prototype implements Map.Entry<String, String> {

    // FIELDS
    /** protein sequence */
    private String protein;
    /** annotation */
    private String annotation;

    /**
     * Construct a new annotation prototype.
     *
     * @param prot		protein sequence
     * @param anno		relevant annotation
     */
    public Prototype(String prot, String anno) {
        this.protein = prot;
        this.annotation = anno;
    }

    @Override
    public String getKey() {
        return this.protein;
    }

    @Override
    public String getValue() {
        return this.annotation;
    }

    @Override
    public String setValue(String value) {
        String retVal = this.annotation;
        this.annotation = value;
        return retVal;
    }

    /**
     * @return the protein sequence
     */
    public String getProtein() {
        return this.protein;
    }

    /**
     * @return the associated annotation
     */
    public String getAnnotation() {
        return this.annotation;
    }

}
