package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerRobustnessExtendedTest {

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
	public void testNullReverseRecord() {
		Logger logger = mock(Logger.class);
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		String adapter = "AGATCGGAAGAGCACACGTC";
		trimmer.addClippingSeq(trimmer.new IlluminaMediumClippingSeq(adapter), true, false); // Forward only
		
		String seq = makeSequence(100) + adapter + makeSequence(30);
		FastqRecord rec = new FastqRecord("r1", seq, "", makeQuality('I', 150), 33);
		
		// Input: [rec, null]
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec, null });
		
		assertNotNull(result);
		assertEquals(2, result.length);
		assertNotNull(result[0]);
		assertEquals(100, result[0].getSequence().length());
		assertNull(result[1]);
	}

	@Test
	public void testBothNullRecords() {
		Logger logger = mock(Logger.class);
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { null, null });
		
		assertNotNull(result);
		assertEquals(2, result.length);
		assertNull(result[0]);
		assertNull(result[1]);
	}
}