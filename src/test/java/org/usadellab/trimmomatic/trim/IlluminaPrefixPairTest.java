package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaPrefixPairTest {

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
	public void testAsymmetricAdapterLengths() {
		Logger logger = mock(Logger.class);
		// Lower palindrome threshold to 10.
		// Fragment 6bp + Adapter 10bp + Prefix 10bp = 26bp match. Score ~15.6.
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 10, 10, 8, false);
		
		// Adapter 1: 20bp "AAAAAAAAAAAAAAAAAAAA"
		String ad1 = "AAAAAAAAAAAAAAAAAAAA";
		// Adapter 2: 10bp "TTTTTTTTTT" (RC of the suffix of ad1)
		String ad2 = "TTTTTTTTTT";
		
		// The logic in IlluminaPrefixPair truncates the start of the longer adapter.
		// So ad1 becomes "AAAAAAAAAA" (last 10 bases). P2 is "TTTTTTTTTT".
		// P1 (AAAA...) matches RC(P2) (AAAA...).
		
		trimmer.addPrefixPair(ad1, ad2);
		
		// Fragment: 6bp (CCCCCC). RC: GGGGGG.
		// We need 6bp to ensure the 16-mer seed (Frag + Adapter) fits within the available 10bp Adapter prefix.
		String frag = "CCCCCC";
		String rcFrag = "GGGGGG";
		
		// Read 1: CCCCC + AAAAAAAAAAAAAAAAAAAA...
		// Read 2: GGGGG + TTTTTTTTTT...
		
		String seq1 = frag + ad1 + makeSequence(50);
		String seq2 = rcFrag + ad2 + makeSequence(50);
		
		FastqRecord rec1 = new FastqRecord("r1", seq1, "", makeQuality('I', seq1.length()), 33);
		FastqRecord rec2 = new FastqRecord("r2", seq2, "", makeQuality('I', seq2.length()), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec1, rec2 });
		
		assertNotNull(result);
		// Should clip to fragment length (6)
		assertEquals(6, result[0].getSequence().length());
		assertEquals(frag, result[0].getSequence());
	}
}