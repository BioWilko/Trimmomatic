package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerMediumPartialTest {

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
	public void testMediumPartialOverlap() {
		Logger logger = mock(Logger.class);
		// minSequenceLikelihood = 10.
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		// Medium adapter: 20bp
		String adapter = "AGATCGGAAGAGCACACGTC";
		trimmer.addClippingSeq(trimmer.new IlluminaMediumClippingSeq(adapter), true, false);
		
		// Read ends with the first 18bp of the adapter.
		// 18bp match score: 18 * 0.602 = 10.8 > 10. Should clip.
		// Read: 132bp DNA + 18bp Adapter prefix.
		
		String dna = makeSequence(132);
		String adapterPart = adapter.substring(0, 18);
		String seq = dna + adapterPart;
		
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('I', 150), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		
		assertNotNull(result);
		// Should clip the 18bp adapter part
		assertEquals(132, result[0].getSequence().length());
		assertEquals(dna, result[0].getSequence());
	}
}