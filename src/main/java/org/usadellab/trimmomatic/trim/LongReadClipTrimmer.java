package org.usadellab.trimmomatic.trim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * Long-read adapter clipping using Hamming distance (V0.42, no indels).
 *
 * Scans the 3' end of each read for every loaded adapter sequence. For each
 * adapter, overlap lengths from {@code minOverlap} up to
 * {@code min(readLen, adapterLen)} are tested, scanning longest-first so the
 * most-confident (longest) hit takes precedence. At most
 * {@code floor(overlap × maxErrorRate)} base mismatches are allowed; N bases
 * in either the read or the adapter are treated as wildcards and never
 * contribute to the mismatch count.
 *
 * Adapters are loaded from a FASTA file; both the forward orientation and its
 * reverse complement are added automatically (duplicates are skipped).
 *
 * Step syntax:
 * <pre>
 *   LONGREADCLIP:&lt;fasta&gt;:&lt;maxErrorRate&gt;[:&lt;minOverlap&gt;]
 * </pre>
 *
 * Examples:
 * <pre>
 *   LONGREADCLIP:adapters.fa:0.10:10
 *   LONGREADCLIP:/path/to/ONT_adapters.fasta:0.15
 * </pre>
 *
 * Recommended values:
 * <ul>
 *   <li>maxErrorRate 0.10 – moderate stringency; allows 1 mismatch per 10 bp</li>
 *   <li>maxErrorRate 0.15 – relaxed; suitable for noisier platforms (PacBio CLR)</li>
 *   <li>minOverlap  10   – default; rejects spurious 1–9 bp partial matches</li>
 * </ul>
 */
public class LongReadClipTrimmer extends AbstractSingleRecordTrimmer {

    private static final int DEFAULT_MIN_OVERLAP = 10;

    private final List<String> adapters;
    private final float        maxErrorRate;
    private final int          minOverlap;

    // ------------------------------------------------------------------
    // Construction

    public LongReadClipTrimmer(String args) throws IOException {
        // Windows paths contain a drive-letter colon (e.g. C:\...) which must not be
        // confused with the argument separator.  Parse numeric parameters from the
        // RIGHT of the colon-split array; the file path is everything to the left.
        String[] parts = args.split(":");
        if (parts.length < 2) {
            throw new IllegalArgumentException(
                    "LONGREADCLIP requires at least <fasta>:<maxErrorRate>");
        }

        int    lastIdx        = parts.length - 1;
        int    tempMinOverlap = DEFAULT_MIN_OVERLAP;
        float  tempErrorRate;
        int    fileEndIdx;

        // Try to peel minOverlap (integer) off the rightmost part, then maxErrorRate.
        try {
            tempMinOverlap = Integer.parseInt(parts[lastIdx]);
            if (lastIdx < 1) throw new IllegalArgumentException(
                    "LONGREADCLIP: missing maxErrorRate");
            tempErrorRate = Float.parseFloat(parts[lastIdx - 1]);
            fileEndIdx    = lastIdx - 1; // path is parts[0 .. fileEndIdx-1]
        } catch (NumberFormatException ignored) {
            // Last part is not an integer — try it as maxErrorRate directly.
            tempErrorRate = Float.parseFloat(parts[lastIdx]);
            fileEndIdx    = lastIdx; // path is parts[0 .. fileEndIdx-1]
        }

        // Reconstruct the file path (parts[0] through parts[fileEndIdx-1] joined by ':').
        StringBuilder pathBuilder = new StringBuilder();
        for (int i = 0; i < fileEndIdx; i++) {
            if (i > 0) {
                pathBuilder.append(':');
            }
            pathBuilder.append(parts[i]);
        }
        if (pathBuilder.length() == 0) {
            throw new IllegalArgumentException("LONGREADCLIP: missing FASTA file path");
        }

        File fastaFile = new File(pathBuilder.toString());
        maxErrorRate   = tempErrorRate;
        minOverlap     = tempMinOverlap;

        adapters = new ArrayList<>();
        loadAdapters(fastaFile);

        if (adapters.isEmpty()) {
            throw new IllegalArgumentException(
                    "LONGREADCLIP: no adapter sequences found in " + fastaFile);
        }
    }

    private void loadAdapters(File fastaFile) throws IOException {
        try (BufferedReader br = new BufferedReader(new FileReader(fastaFile))) {
            String line;
            StringBuilder seq = new StringBuilder();
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (line.startsWith(">")) {
                    if (seq.length() > 0) {
                        addAdapter(seq.toString());
                        seq.setLength(0);
                    }
                } else if (!line.isEmpty()) {
                    seq.append(line);
                }
            }
            if (seq.length() > 0) {
                addAdapter(seq.toString());
            }
        }
    }

    private void addAdapter(String seq) {
        String upper = seq.toUpperCase();
        adapters.add(upper);
        String rc = reverseComplement(upper);
        if (!rc.equals(upper)) {
            adapters.add(rc);
        }
    }

    // ------------------------------------------------------------------
    // Sequence utilities

    private static String reverseComplement(String seq) {
        int len = seq.length();
        StringBuilder sb = new StringBuilder(len);
        for (int i = len - 1; i >= 0; i--) {
            sb.append(complement(seq.charAt(i)));
        }
        return sb.toString();
    }

    private static char complement(char c) {
        return switch (c) {
            case 'A' -> 'T';
            case 'T' -> 'A';
            case 'G' -> 'C';
            case 'C' -> 'G';
            default  -> 'N';
        };
    }

    // ------------------------------------------------------------------
    // Trimming logic

    @Override
    public FastqRecord processRecord(FastqRecord in) {
        int readLen = in.getLength();
        if (readLen == 0) {
            return null;
        }

        String seq     = in.getSequence();
        int    trimPos = readLen; // sentinel: no trim found yet

        for (String adapter : adapters) {
            int adapterLen = adapter.length();
            int maxOverlap = Math.min(readLen, adapterLen);

            if (maxOverlap < minOverlap) {
                continue;
            }

            // Scan from longest overlap downward so the first hit is the most confident.
            for (int overlap = maxOverlap; overlap >= minOverlap; overlap--) {
                int readStart = readLen - overlap;
                if (readStart >= trimPos) {
                    // Any hit here would be to the right of an already-found trim point.
                    break;
                }

                int allowedMismatches = (int) (overlap * maxErrorRate);
                int mismatches        = 0;
                boolean exceeded      = false;

                for (int i = 0; i < overlap; i++) {
                    char rc = seq.charAt(readStart + i);
                    char ac = adapter.charAt(i);
                    // N in either strand is a wildcard — not counted as a mismatch.
                    if (rc != ac && ac != 'N' && rc != 'N') {
                        mismatches++;
                        if (mismatches > allowedMismatches) {
                            exceeded = true;
                            break;
                        }
                    }
                }

                if (!exceeded) {
                    trimPos = readStart; // update leftmost trim point for this adapter
                    break;              // take longest (most confident) match per adapter
                }
            }
        }

        if (trimPos == 0) {
            return null; // adapter found at position 0 — drop entire read
        }
        if (trimPos == readLen) {
            return in;   // no adapter detected
        }
        return new FastqRecord(in, 0, trimPos);
    }
}
