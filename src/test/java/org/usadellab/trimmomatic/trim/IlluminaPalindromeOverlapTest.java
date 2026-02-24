package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaPalindromeOverlapTest {

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
			switch (c) {
			case 'A': sb.append('T'); break;
			case 'C': sb.append('G'); break;
			case 'G': sb.append('C'); break;
			case 'T': sb.append('A'); break;
			default: sb.append(c); break;
			}
		}
		return sb.toString();
	}

	@Test
	public void testPalindromeOverlapBoundary() {
		Logger logger = mock(Logger.class);
		// Adapters are 33bp. Two adapters contribute ~40 to the score (66 * 0.602).
		// We want to test fragment overlap boundary.
		// Set threshold to 65.
		// Case 1: 40bp fragment. Total 106bp. Score 63.8 < 65. No Clip.
		// Case 2: 50bp fragment. Total 116bp. Score 69.8 > 65. Clip.
		
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 65, 10, 8, false);
		
		String fwdAdapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
		String revAdapter = "TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"; // RC of fwdAdapter
		trimmer.addPrefixPair(fwdAdapter, revAdapter);

		// Case 1: 40bp overlap (should NOT clip)
		String frag40 = makeSequence(40);
		String rcFrag40 = reverseComplement(frag40);
		
		String seq1_40 = frag40 + fwdAdapter + makeSequence(150 - 40 - fwdAdapter.length());
		String seq2_40 = rcFrag40 + revAdapter + makeSequence(150 - 40 - revAdapter.length());
		
		FastqRecord rec1_40 = new FastqRecord("r1", seq1_40, "", makeQuality('I', 150), 33);
		FastqRecord rec2_40 = new FastqRecord("r2", seq2_40, "", makeQuality('I', 150), 33);
		
		FastqRecord[] res40 = trimmer.processRecords(new FastqRecord[] { rec1_40, rec2_40 });
		assertEquals(150, res40[0].getSequence().length());

		// Case 2: 50bp overlap (should clip)
		String frag50 = makeSequence(50);
		String rcFrag50 = reverseComplement(frag50);
		
		String seq1_50 = frag50 + fwdAdapter + makeSequence(150 - 50 - fwdAdapter.length());
		String seq2_50 = rcFrag50 + revAdapter + makeSequence(150 - 50 - revAdapter.length());
		
		FastqRecord rec1_50 = new FastqRecord("r1", seq1_50, "", makeQuality('I', 150), 33);
		FastqRecord rec2_50 = new FastqRecord("r2", seq2_50, "", makeQuality('I', 150), 33);
		
		FastqRecord[] res50 = trimmer.processRecords(new FastqRecord[] { rec1_50, rec2_50 });
		assertEquals(50, res50[0].getSequence().length());
	}
}