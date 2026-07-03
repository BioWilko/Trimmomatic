package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaPalindromeTest {

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
	public void testPalindromeLowQualityMismatch() {
		Logger logger = mock(Logger.class);
		// minPalindromeLikelihood = 30
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		String fwdAdapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		String revAdapter = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
		trimmer.addPrefixPair(fwdAdapter, revAdapter);

		// 50bp fragment
		String fragment = makeSequence(50);
		String rcFragment = reverseComplement(fragment);
		
		// Construct reads with mismatches in the overlap region
		// Overlap is 50bp.
		// Perfect match score: 50 * 0.602 = 30.1
		// Threshold is 30.
		// Even a small penalty will drop it below 30.
		
		// Read 1: Fragment + FwdAdapter
		String seq1 = fragment + fwdAdapter + makeSequence(100 - fwdAdapter.length());
		// Read 2: RCFragment + RevAdapter (but with a mismatch in the fragment part)
		
		// Introduce mismatch in the fragment part of read 2
		// Change last base of rcFragment
		String rcFragmentMut = rcFragment.substring(0, 49) + (rcFragment.charAt(49) == 'A' ? 'T' : 'A');
		String seq2 = rcFragmentMut + revAdapter + makeSequence(100 - revAdapter.length());
		
		// High quality mismatch (Q40) -> Penalty 4.
		// Score approx 30.1 - 4 = 26.1 < 30. Should NOT clip.
		
		FastqRecord rec1 = new FastqRecord("r1", seq1, "", makeQuality('I', 150), 33);
		FastqRecord rec2 = new FastqRecord("r2", seq2, "", makeQuality('I', 150), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec1, rec2 });
		
		// Should not clip
		assertEquals(150, result[0].getSequence().length());
		assertEquals(150, result[1].getSequence().length());
	}
}