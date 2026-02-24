package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerEdgeCasesTest {

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
	public void testAdapterLongerThanRead() {
		Logger logger = mock(Logger.class);
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		// Adapter 200bp
		String adapter = makeSequence(200);
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter), true, false);
		
		// Read 150bp. Contains the start of the adapter.
		String seq = adapter.substring(0, 150);
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('I', 150), 33);
		
		// Should clip everything (return null) because the read is entirely adapter
		FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
		assertNotNull(result);
		if (result.length > 0) {
			assertNull(result[0]);
		}
	}

	@Test
	public void testAdapterWithNs() {
		Logger logger = mock(Logger.class);
		IlluminaClippingTrimmer trimmer = new IlluminaClippingTrimmer(logger, 2, 30, 10, 8, false);
		
		// Adapter 25bp with Ns.
		// Ns in adapter match Ns in read in the seed check (0^0=0), but score 0 in quality check.
		String adapter25 = "ACGTNACGTNACGTNACGTNACGTN"; // 25bp
		trimmer.addClippingSeq(trimmer.new IlluminaLongClippingSeq(adapter25), true, false);
		
		// Read matches adapter exactly (including Ns)
		String seq25 = "ACGTNACGTNACGTNACGTNACGTN" + makeSequence(125);
		FastqRecord rec25 = new FastqRecord("head", seq25, "", makeQuality('I', 150), 33);
		
		FastqRecord[] result25 = trimmer.processRecords(new FastqRecord[] { rec25 });
		
		// Score calculation:
		// 20 matches * 0.6 = 12.
		// 5 Ns * 0 = 0.
		// Total 12 > 10 (threshold).
		// Should clip. Since adapter is at start, read is dropped.
		
		assertNotNull(result25);
		if (result25.length > 0) {
			assertNull(result25[0]);
		}
	}
}