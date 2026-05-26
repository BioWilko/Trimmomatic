package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.File;
import java.io.FileWriter;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * Comprehensive tests for LongReadClipTrimmer.
 *
 * Covered here:
 *   - Basic 3' adapter clipping (perfect match)
 *   - Partial overlap at 3' end is detected and trimmed
 *   - Adapter not present → read returned unchanged
 *   - Partial overlap shorter than minOverlap → no trim
 *   - Mismatches within the allowed error rate → trimmed
 *   - Mismatches exceeding the error rate → no trim
 *   - Reverse complement of adapter detected
 *   - Palindrome adapter not double-counted
 *   - Entire read is adapter → dropped (null)
 *   - minOverlap=1 trims a single-base overlap
 *   - N bases in the read or adapter are wildcards (not counted as mismatches)
 *   - View records are trimmed correctly
 *   - Quality string is trimmed in sync with the sequence
 *   - Multiple adapters in FASTA: all checked, leftmost trim wins
 *   - Multi-line FASTA sequences assembled correctly
 *   - Empty FASTA file → IllegalArgumentException
 *   - Missing argument string → Exception
 */
public class LongReadClipTrimmerTest {

    @TempDir
    Path tempDir;

    // ------------------------------------------------------------------
    // Helpers

    private File writeFasta(String... entries) throws Exception {
        File f = tempDir.resolve("adapters.fasta").toFile();
        try (FileWriter fw = new FileWriter(f)) {
            for (String e : entries) {
                fw.write(e);
            }
        }
        return f;
    }

    private File singleAdapterFasta(String seq) throws Exception {
        return writeFasta(">adapter\n" + seq + "\n");
    }

    private FastqRecord makeRecord(String seq) {
        return new FastqRecord("r1", seq, "", "I".repeat(seq.length()), 33);
    }

    private FastqRecord makeRecord(String seq, String qual) {
        return new FastqRecord("r1", seq, "", qual, 33);
    }

    private LongReadClipTrimmer trimmer(File fasta, float errorRate, int minOverlap) throws Exception {
        return new LongReadClipTrimmer(fasta.getPath() + ":" + errorRate + ":" + minOverlap);
    }

    private LongReadClipTrimmer trimmer(File fasta, float errorRate) throws Exception {
        return new LongReadClipTrimmer(fasta.getPath() + ":" + errorRate);
    }

    // ------------------------------------------------------------------
    // Basic trimming

    @Test
    public void testPerfectFullAdapterMatch_3prime() throws Exception {
        File fasta = singleAdapterFasta("TTTTTTTT");
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 8);
        FastqRecord rec = makeRecord("ACGTACGTTTTTTTTT"); // 8 payload + 8 adapter
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTACGT", result.getSequence());
    }

    @Test
    public void testPartialAdapterOverlap_3prime() throws Exception {
        // Adapter 12-T's; read ends with 8 T's (payload has no T's so there is no
        // ambiguity about where the adapter starts).
        // Overlap=8 >= minOverlap=5, 0 mismatches → trimmed.
        File fasta = singleAdapterFasta("TTTTTTTTTTTT");
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 5);
        FastqRecord rec = makeRecord("GCGCGCGCTTTTTTTT"); // 8-G/C payload + 8 T's
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("GCGCGCGC", result.getSequence());
    }

    @Test
    public void testNoAdapterPresent_readUnchanged() throws Exception {
        File fasta = singleAdapterFasta("TTTTTTTT");
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 5);
        FastqRecord rec = makeRecord("ACGTACGTACGTACGT");
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTACGTACGTACGT", result.getSequence());
    }

    @Test
    public void testOverlapBelowMinOverlap_noTrim() throws Exception {
        // Read ends with 4 T's; minOverlap=5 → no trim.
        File fasta = singleAdapterFasta("TTTTTTTT");
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 5);
        // "ACGTACGTACGTTTTT" ends in 4 T's (overlap with TTTTTTTT = 4 < minOverlap=5).
        // Wait - last 5 chars = "CTTTT" which has C≠T mismatch. Last 4 = "TTTT" vs first 4 of adapter = "TTTT" → match of 4, below minOverlap=5.
        FastqRecord rec = makeRecord("ACGTACGTACGTTTTTT"); // ends in 5 T's, len=17; overlap=5 >= 5? yes
        // Actually: overlap=5, readStart=12, adapter[0..4]="TTTTT", read[12..16]="TTTTT" → match.
        // Let me use fewer T's: 3 T's → overlap=3 < 5 → no trim.
        FastqRecord rec2 = makeRecord("ACGTACGTACGTTTT"); // ends in 3 T's
        FastqRecord result = t.processRecord(rec2);
        assertNotNull(result);
        assertEquals("ACGTACGTACGTTTT", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Mismatch tolerance

    @Test
    public void testOneMismatch_withinErrorRate_trimmed() throws Exception {
        // Adapter TTTTTTTT. floor(8 * 0.15) = 1 allowed mismatch.
        File fasta = singleAdapterFasta("TTTTTTTT");
        LongReadClipTrimmer t = trimmer(fasta, 0.15f, 8);
        FastqRecord rec = makeRecord("ACGTACGTTTTTATTT"); // 1 mismatch at position 4
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTACGT", result.getSequence());
    }

    @Test
    public void testOneMismatch_exceedsErrorRate_noTrim() throws Exception {
        // floor(8 * 0.10) = 0 allowed mismatches.
        File fasta = singleAdapterFasta("TTTTTTTT");
        LongReadClipTrimmer t = trimmer(fasta, 0.10f, 8);
        FastqRecord rec = makeRecord("ACGTACGTTTTTATTT"); // 1 mismatch
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTACGTTTTTATTT", result.getSequence());
    }

    @Test
    public void testTwoMismatches_allowedAt25pct() throws Exception {
        // floor(8 * 0.25) = 2 allowed mismatches.
        // Adapter TTTTTTTT. Tail "TTATATTT" has exactly 2 mismatches (pos 2 and pos 4).
        File fasta = singleAdapterFasta("TTTTTTTT");
        LongReadClipTrimmer t = trimmer(fasta, 0.25f, 8);
        FastqRecord rec = makeRecord("ACGTACGTTTATATTT"); // 2 mismatches (A at pos 2,4 of tail)
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTACGT", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Reverse complement detection

    @Test
    public void testReverseComplementAdapterDetected() throws Exception {
        // Adapter AAAA; RC = TTTT. Read ends with 4 T's.
        File fasta = singleAdapterFasta("AAAA");
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 4);
        FastqRecord rec = makeRecord("ACGTACGTTTTT"); // ends in 4 T's
        // Actually: last 4 chars are "TTTT" = RC(AAAA). Trim at readLen-4=8.
        // But "ACGTACGTTTTT" is length 12; last 4 = "TTTT" → trim at 8.
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTACGT", result.getSequence());
    }

    @Test
    public void testPalindromeAdapter_addedOnce_stillWorks() throws Exception {
        // AATT → RC = AATT (palindrome). Loaded only once, works correctly.
        File fasta = singleAdapterFasta("AATT");
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 4);
        FastqRecord rec = makeRecord("GCGCAATT");
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("GCGC", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Entire read is adapter → dropped

    @Test
    public void testEntireReadIsAdapter_dropped() throws Exception {
        File fasta = singleAdapterFasta("TTTTTTTT");
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 8);
        FastqRecord rec = makeRecord("TTTTTTTT");
        assertNull(t.processRecord(rec));
    }

    // ------------------------------------------------------------------
    // minOverlap = 1

    @Test
    public void testMinOverlapOne_singleBaseMatch_trimmed() throws Exception {
        File fasta = singleAdapterFasta("T");
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 1);
        FastqRecord rec = makeRecord("ACGCACGCT"); // ends with T
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGCACGC", result.getSequence());
    }

    // ------------------------------------------------------------------
    // N bases as wildcards

    @Test
    public void testNInRead_isWildcard() throws Exception {
        // Adapter TTTTTTTT. Read ends with TNNNNNNN (T + 7 N's). All N's are wildcards.
        File fasta = singleAdapterFasta("TTTTTTTT");
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 8);
        FastqRecord rec = makeRecord("ACGTACGTTNNNNNNN");
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTACGT", result.getSequence());
    }

    @Test
    public void testNInAdapter_isWildcard() throws Exception {
        // Adapter NNNNNNNN (all N). Matches anything.
        File fasta = singleAdapterFasta("NNNNNNNN");
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 8);
        FastqRecord rec = makeRecord("ACGTACGTACGTACGT");
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTACGT", result.getSequence());
    }

    // ------------------------------------------------------------------
    // View records

    @Test
    public void testViewRecord_trimmedCorrectly() throws Exception {
        // Base: 8 padding + "GCGCTTTTTTTT" (4-base payload + 8-T adapter)
        FastqRecord base = makeRecord("XXXXXXXXGCGCTTTTTTTT");
        FastqRecord view = new FastqRecord(base, 8, 12); // "GCGCTTTTTTTT"
        File fasta = singleAdapterFasta("TTTTTTTT");
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 8);
        FastqRecord result = t.processRecord(view);
        assertNotNull(result);
        assertEquals("GCGC", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Quality string is trimmed in sync

    @Test
    public void testQualityTrimmedInSync() throws Exception {
        // Sequence "GCGCTTTT", quality "IIIIJJJJ". After adapter trim: "GCGC" / "IIII".
        File fasta = singleAdapterFasta("TTTT");
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 4);
        FastqRecord rec = makeRecord("GCGCTTTT", "IIIIJJJJ");
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("GCGC",    result.getSequence());
        assertEquals("IIII",    result.getQuality());
    }

    // ------------------------------------------------------------------
    // Multiple adapters

    @Test
    public void testMultipleAdapters_correctAdapterTrimmed() throws Exception {
        File fasta = writeFasta(">adp1\nAAAA\n", ">adp2\nTTTT\n");
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 4);
        FastqRecord rec = makeRecord("ACGTTTTT"); // ends with 4 T's
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGT", result.getSequence());
    }

    @Test
    public void testMultipleAdapters_leftmostTrimWins() throws Exception {
        // Adapter G = GGGGGGGG (8 G's), Adapter C = CCCC (4 C's).
        // Read: "XXGGGGGGGG" (10 chars). Adapter G matches 8 G's at readStart=2.
        // Adapter C (RC=GGGG): "GGGG" would match the last 4 of the G run at readStart=6.
        // Leftmost trim wins: readStart=2.
        File fasta = writeFasta(">adpG\nGGGGGGGG\n", ">adpC\nCCCC\n");
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 4);
        FastqRecord rec = makeRecord("XXGGGGGGGG");
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("XX", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Multi-line FASTA

    @Test
    public void testMultiLineFastaAssembled() throws Exception {
        File fasta = writeFasta(">adapter\nTTTT\nTTTT\n"); // sequence = TTTTTTTT
        LongReadClipTrimmer t = trimmer(fasta, 0.0f, 8);
        FastqRecord rec = makeRecord("ACGTACGTTTTTTTTT");
        FastqRecord result = t.processRecord(rec);
        assertNotNull(result);
        assertEquals("ACGTACGT", result.getSequence());
    }

    // ------------------------------------------------------------------
    // Error cases

    @Test
    public void testEmptyFasta_throwsIllegalArgument() throws Exception {
        File fasta = tempDir.resolve("empty.fasta").toFile();
        fasta.createNewFile();
        assertThrows(IllegalArgumentException.class,
                () -> new LongReadClipTrimmer(fasta.getPath() + ":0.1:10"));
    }

    @Test
    public void testMissingArg_throwsException() {
        assertThrows(Exception.class,
                () -> new LongReadClipTrimmer("noarguments"));
    }
}
