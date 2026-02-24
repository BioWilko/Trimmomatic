package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerPriorityTest {

	private String makeSequence(int len) {
		StringBuilder sb = new StringBuilder(len);
		String pattern = "ACGT";
		for (int i = 0; i < len; i++)
			sb.append(pattern.charAt(i % 4));
		return sb.toString();
	}

	private String makeQuality(char c, int len) {
		StringBuilder sb = new StringBuilder(len);
		for (int i = 0; i < len; i++)
			sb.append(c);
		return sb.toString();
	}

	@Test
	public void testPalindromeVsSimplePriority() {
		Logger logger = mock(Logger.class);
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		String fwdAdapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		String revAdapter = "TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"; // RC of fwdAdapter
		trimmer.addPrefixPair(fwdAdapter, revAdapter);
		
		// Add a simple adapter that appears EARLIER in the read than the palindrome match
		String simpleAdapter = "TGCATGCATGCATGCATGCA";
		trimmer.addClippingSeq(trimmer.new IlluminaMediumClippingSeq(simpleAdapter), true, false);
		
		// Construct read:
		// [Fragment 40bp] [SimpleAdapter] [More Fragment] [Palindrome Adapter]
		// Palindrome match would clip at end of fragment (say 100bp).
		// Simple match clips at 40bp.
		// Result should be 40bp.
		
		String part1 = makeSequence(40);
		String part2 = makeSequence(40);
		
		// Read 1: part1 + SimpleAdapter + part2 + FwdAdapter
		String seq1 = part1 + simpleAdapter + part2 + fwdAdapter;
		
		// Read 2: Just needs to trigger palindrome match. 
		// Let's make it a perfect read-through of (part1 + simple + part2).
		// Length = 40 + 20 + 40 = 100.
		// RC of that.
		// Then RevAdapter.
		String insert = part1 + simpleAdapter + part2;
		// Note: reverseComplement helper needed, but for this test we can just assume palindrome triggers if we set it up right.
		// Actually, simpler: Just rely on the fact that processRecords checks both.
		// We don't need a valid Read 2 for Simple Clipping to work on Read 1.
		// But to trigger Palindrome check, we need a pair.
		// Let's just pass a dummy Read 2 that WON'T trigger palindrome, to verify Simple works in PE mode.
		// Wait, the test is about PRIORITY. So we DO need palindrome to trigger.
		// Let's construct a valid palindrome pair.
		
		// Actually, constructing a valid palindrome pair that ALSO has a simple adapter embedded is complex.
		// Alternative: Read 1 has Simple Adapter at 40.
		// Palindrome logic detects "insert size" of 100.
		// So Palindrome says "Keep 100". Simple says "Keep 40".
		// Trimmer should keep min(100, 40) = 40.
		
		// We need a helper for RC.
		StringBuilder sb = new StringBuilder(insert.length());
		for (int i = insert.length() - 1; i >= 0; i--) {
			char c = insert.charAt(i);
			switch (c) { case 'A': sb.append('T'); break; case 'C': sb.append('G'); break; case 'G': sb.append('C'); break; case 'T': sb.append('A'); break; default: sb.append(c); break; }
		}
		String rcInsert = sb.toString();
		
		String seq2 = rcInsert + revAdapter;
		
		FastqRecord rec1 = new FastqRecord("r1", seq1, "", makeQuality('I', seq1.length()), 33);
		FastqRecord rec2 = new FastqRecord("r2", seq2, "", makeQuality('I', seq2.length()), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec1, rec2 });
		
		assertNotNull(result);
		// Should clip at 40 (Simple) instead of 100 (Palindrome)
		assertEquals(40, result[0].getSequence().length());
	}
}