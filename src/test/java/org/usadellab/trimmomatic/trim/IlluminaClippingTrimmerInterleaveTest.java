package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerInterleaveTest {

	private String makeSequence(int len) {
		StringBuilder sb = new StringBuilder(len);
		for (int i = 0; i < len; i++)
			sb.append('C');
		return sb.toString();
	}

	private String makeQuality(char c, int len) {
		StringBuilder sb = new StringBuilder(len);
		for (int i = 0; i < len; i++)
			sb.append(c);
		return sb.toString();
	}

	@Test
	public void testInterleaveOffsets() {
		Logger logger = mock(Logger.class);
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		// 32bp adapter (Long clipping seq)
		String adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTC"; 
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);
		
		// Test offsets 0, 1, 2, 3 to cover all interleave phases
		for (int offset = 0; offset < 4; offset++) {
			// Construct read: [Offset padding] + [Adapter] + [Rest]
			// We want to clip the adapter.
			// If offset is 0, adapter starts at 0.
			// If offset is 1, adapter starts at 1.
			
			String padding = makeSequence(offset);
			// Ensure padding doesn't match adapter start (A)
			if (padding.endsWith("A")) padding = padding.substring(0, padding.length()-1) + "C";
			
			String seq = padding + adapter + makeSequence(50);
			// Pad to 150bp if needed, but length doesn't strictly matter as long as it's long enough
			
			FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('I', seq.length()), 33);
			
			FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
			
			assertNotNull(result);
			
			if (offset == 0) {
				assertNull(result[0], "Record should be dropped (null) when offset is 0");
			} else {
				// Should clip at the start of the adapter (which is at 'offset')
				// So remaining length should be 'offset'
				assertEquals(offset, result[0].getSequence().length(), "Failed at offset " + offset);
				assertEquals(padding, result[0].getSequence());
			}
		}
	}
}