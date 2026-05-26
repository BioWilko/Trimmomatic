package org.usadellab.trimmomatic.trim;

import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * MAXAMBIG:<maxFraction>
 *
 * Drops the read if the fraction of ambiguous bases (N) exceeds maxFraction.
 * maxFraction is a value between 0.0 (no Ns allowed) and 1.0 (all Ns allowed).
 *
 * Example: MAXAMBIG:0.1  drops reads with more than 10% N bases.
 */
public class MaxAmbigTrimmer extends AbstractSingleRecordTrimmer {
    private final double maxFraction;

    public MaxAmbigTrimmer(String args) {
        maxFraction = Double.parseDouble(args);
        if (maxFraction < 0.0 || maxFraction > 1.0)
            throw new IllegalArgumentException(
                    "MAXAMBIG fraction must be between 0.0 and 1.0, got: " + maxFraction);
    }

    @Override
    public FastqRecord processRecord(FastqRecord in) {
        int len = in.getLength();
        if (len == 0)
            return null;

        String seq = in.getSequence();
        int nCount = 0;
        for (int i = 0; i < len; i++) {
            if (seq.charAt(i) == 'N')
                nCount++;
        }

        return ((double) nCount / len) <= maxFraction ? in : null;
    }
}
