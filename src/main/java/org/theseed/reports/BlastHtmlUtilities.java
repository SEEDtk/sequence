/**
 *
 */
package org.theseed.reports;

import static j2html.TagCreator.*;

import org.theseed.locations.Location;
import org.theseed.sequence.blast.BlastDB;
import org.theseed.sequence.blast.BlastHit;
import org.theseed.sequence.blast.BlastHit.SeqData;

import j2html.tags.ContainerTag;

/**
 * This class contains useful utilities for BLAST reporting in HTML.
 *
 * @author Bruce Parrello
 *
 */
public class BlastHtmlUtilities {

    /**
     * @return a description of the color usage on a BLAST report
     *
     * @param colorType		color type scheme
     */
    public static ContainerTag showColorInfo(BlastDB.ColorType colorType) {
        return div().withClass("pod").with(p("Color is determined by " + colorType.description() + "."))
                .with(table(tr(
                        th("100%").withStyle("background-color: " + Color.BLUE.html() + "; color: white"),
                        th("90% to 99%").withStyle("background-color: " + Color.DARK_GREEN.html() + ";"),
                        th("70% to 89%").withStyle("background-color: " + Color.ORANGE.html() + ";"),
                        th("50% to 69%").withStyle("background-color: " + Color.RED.html() + ";"),
                        th("0% to 49%").withStyle("background-color: " + Color.DARK_GRAY.html() + "; color: white")
                        )));
    }

    /**
     * @return HTML reporting the basic information about the BLAST
     *
     * @param blastInfo		blast information object
     * @param sequencesHit	number of sequences hit
     * @param sortType		sort type
     * @param removed 		number of hits removed
     */
    public static ContainerTag showBlastInfo(BlastInfo blastInfo, int sequencesHit, BlastDB.SortType sortType, int removed) {
        ContainerTag retVal = ul(
                li(blastInfo.getParms()),
                li(String.format("%d queries produced %d hits.", blastInfo.getQueriesIn(),
                        blastInfo.getHitCount())),
                li(String.format("%d queries had no hits.", blastInfo.getMissCount())),
                li(String.format("%d %s had hits.", sequencesHit,
                        sortType.getPlural()))
                );
        if (removed > 0) {
            retVal.with(li(String.format("%d redundant hits are not shown.", removed)));
        }
        return retVal;
    }

    /**
     * @return the HTML hit arrow for a hit.
     *
     * @param target	target (hitting) sequence
     * @param anchor	anchor (hit) sequence
     * @param hit		blast hit to turn into an arrow
     * @param color		color to give to the arrow
     */
    public static HtmlHitSequence getHitArrow(SeqData target, SeqData anchor, BlastHit hit, Color color) {
        Location hitLoc = target.getLoc();
        char dir = (hitLoc.getDir() != anchor.getLoc().getDir() ? '-' : '+');
        String label = String.format("[e=%4.2e, ident=%4.1f%%, gap=%d, loc=(%d,%d)/%d] %s",
                hit.getEvalue(), hit.getPercentIdentity(), hit.getNumGap(), hitLoc.getBegin(),
                hitLoc.getEnd(), target.getLen(), target.getDef());
        HtmlHitSequence hitDescriptor = new HtmlHitSequence(target.getId(), label,
                anchor.getLoc(), dir, color);
        return hitDescriptor;
    }

}
