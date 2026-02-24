package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerDirectionalTest {

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
	public void testDirectionalClipping() {
		Logger logger = mock(Logger.class);
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		String fwdAdapter = "AGATCGGAAGAGCACACGTC"; // 20bp non-repetitive
		String revAdapter = "AGATCGGAAGAGCGTCGTGT"; // 20bp non-repetitive
		
		// Add fwdAdapter as Forward-only (true, false)
		trimmer.addClippingSeq(trimmer.new IlluminaMediumClippingSeq(fwdAdapter), true, false);
		// Add revAdapter as Reverse-only (false, true)
		trimmer.addClippingSeq(trimmer.new IlluminaMediumClippingSeq(revAdapter), false, true);
		
		// Construct reads
		// Read 1 (Forward): Contains revAdapter (should NOT clip)
		String seq1 = makeSequence(100) + revAdapter + makeSequence(30);
		// Read 2 (Reverse): Contains fwdAdapter (should NOT clip)
		String seq2 = makeSequence(100) + fwdAdapter + makeSequence(30);
		
		FastqRecord rec1 = new FastqRecord("r1", seq1, "", makeQuality('I', 150), 33);
		FastqRecord rec2 = new FastqRecord("r2", seq2, "", makeQuality('I', 150), 33);
		
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec1, rec2 });
		
		assertEquals(150, result[0].getSequence().length());
		assertEquals(150, result[1].getSequence().length());
		
		// Now construct reads that match the direction
		// Read 3 (Forward): Contains fwdAdapter (SHOULD clip)
		String seq3 = makeSequence(100) + fwdAdapter + makeSequence(30);
		// Read 4 (Reverse): Contains revAdapter (SHOULD clip)
		String seq4 = makeSequence(100) + revAdapter + makeSequence(30);
		
		FastqRecord rec3 = new FastqRecord("r3", seq3, "", makeQuality('I', 150), 33);
		FastqRecord rec4 = new FastqRecord("r4", seq4, "", makeQuality('I', 150), 33);
		
		FastqRecord[] result2 = trimmer.processRecords(new FastqRecord[] { rec3, rec4 });
		
		assertEquals(100, result2[0].getSequence().length());
		assertEquals(100, result2[1].getSequence().length());
	}
	
	@Test
	public void testCommonClipping() {
		Logger logger = mock(Logger.class);
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		String commonAdapter = "AGATCGGAAGAGCACACGTC";
		// Add as Common (true, true)
		trimmer.addClippingSeq(trimmer.new IlluminaMediumClippingSeq(commonAdapter), true, true);
		
		String seq = makeSequence(100) + commonAdapter + makeSequence(30);
		FastqRecord rec = new FastqRecord("r", seq, "", makeQuality('I', 150), 33);
		
		// Should clip whether passed as forward or reverse
		assertEquals(100, trimmer.processRecords(new FastqRecord[] { rec })[0].getSequence().length());
		assertEquals(100, trimmer.processRecords(new FastqRecord[] { null, rec })[1].getSequence().length());
	}
}