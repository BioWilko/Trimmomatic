package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerInvalidCharTest {

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
	public void testInvalidCharInAdapter() {
		Logger logger = mock(Logger.class);
		// seedMaxMiss = 0 (Strict)
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 0, 30, 10, 8, false);
		
		// Adapter with 'Z'. 'Z' is not A, C, G, T.
		// Use 21bp length with Z in the middle.
		// This ensures that all 16-mer seeds overlap the Z and fail to match (since Z packs to 0).
		// Even if scoring were triggered, the Z splits the match into two 10bp blocks (score 6.0 each),
		// which when merged with penalty (-4) gives 8.0 < 10. No Clip.
		String adapter = "AAAAAAAAAAZAAAAAAAAAA"; // 21bp
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);
		
		// Read matches the adapter, but has 'A' where adapter has 'Z'.
		String seq = "AAAAAAAAAAAAAAAAAAAAA" + makeSequence(129);
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('I', 150), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		// Should NOT clip because 'Z' vs 'A' is a mismatch (both pack to 0? No, A packs to 1/2/4/8).
		// packCh('Z') returns 0. packCh('A') returns non-zero.
		// So bitwise check fails.
		assertEquals(150, result[0].getSequence().length());
	}
}