package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerPEvsSETest {

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

	private String reverseComplement(String seq) {
		StringBuilder sb = new StringBuilder(seq.length());
		for (int i = seq.length() - 1; i >= 0; i--) {
			char c = seq.charAt(i);
			sb.append(switch (c) {
				case 'A' -> 'T';
				case 'C' -> 'G';
				case 'G' -> 'C';
				case 'T' -> 'A';
				default -> c;
			});
		}
		return sb.toString();
	}

	@Test
	public void testPalindromeAdapterBehaviorInSE() {
		Logger logger = mock(Logger.class);
		// Setup trimmer with a prefix pair (Palindrome mode adapters)
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		String fwdAdapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		String revAdapter = "TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"; // RC of fwdAdapter
		
		// Add ONLY as a pair. This simulates loading from a file where these are defined as /1 and /2
		// and thus removed from the simple clipping sets.
		trimmer.addPrefixPair(fwdAdapter, revAdapter);
		
		// Construct a read that looks like a perfect palindrome read-through
		String fragment = makeSequence(50);
		String rcFragment = reverseComplement(fragment);
		
		String seq1 = fragment + fwdAdapter + makeSequence(150 - 50 - fwdAdapter.length());
		String seq2 = rcFragment + revAdapter + makeSequence(150 - 50 - revAdapter.length());
		
		String qual = makeQuality('I', 150);
		
		FastqRecord rec1 = new FastqRecord("r1", seq1, "", qual, 33);
		FastqRecord rec2 = new FastqRecord("r2", seq2, "", qual, 33);
		
		// 1. Process as Paired End (both records)
		// Should clip because palindrome logic is triggered.
		FastqRecord[] resultPE = trimmer.processRecords(new FastqRecord[] { rec1, rec2 });
		assertEquals(50, resultPE[0].getSequence().length());
		// Reverse read dropped by default
		assertNotNull(resultPE);
		assertEquals(null, resultPE[1]);
		
		// 2. Process as Single End (one record at a time)
		// Should NOT clip because the adapter is only known as part of a pair, 
		// and simple clipping sets are empty.
		
		FastqRecord[] resultSE1 = trimmer.processRecords(new FastqRecord[] { rec1 });
		assertEquals(150, resultSE1[0].getSequence().length()); // No clip
	}
}