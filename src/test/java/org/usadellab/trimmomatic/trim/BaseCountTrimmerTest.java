package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class BaseCountTrimmerTest {

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
	public void testBaseCountKeep() {
		// 150bp, 2 Ns. Range 0-5.
		BaseCountTrimmer trimmer = new BaseCountTrimmer("N:0:5");
		String seq = makeSequence(148) + "NN";
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('!', 150), 33);
		assertNotNull(trimmer.processRecord(rec));
	}

	@Test
	public void testBaseCountDropMax() {
		// 150bp, 6 Ns. Range 0-5.
		BaseCountTrimmer trimmer = new BaseCountTrimmer("N:0:5");
		String seq = makeSequence(144) + "NNNNNN";
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('!', 150), 33);
		assertNull(trimmer.processRecord(rec));
	}

	@Test
	public void testBaseCountDropMin() {
		// 150bp, 0 Ns. Range 1-5.
		BaseCountTrimmer trimmer = new BaseCountTrimmer("N:1:5");
		String seq = makeSequence(150); // No Ns
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('!', 150), 33);
		assertNull(trimmer.processRecord(rec));
	}

	@Test
	public void testBaseCountGC_Boundary() {
		// 150bp. ACGT repeating.
		// 150 / 4 = 37 remainder 2 (AC).
		// A: 38, C: 38, G: 37, T: 37.
		// GC count = 38 + 37 = 75.

		// Test exact match
		BaseCountTrimmer trimmer = new BaseCountTrimmer("GC:75:75");
		String seq = makeSequence(150);
		FastqRecord rec = new FastqRecord("head", seq, "", makeQuality('!', 150), 33);
		assertNotNull(trimmer.processRecord(rec));

		// Test just outside
		trimmer = new BaseCountTrimmer("GC:76:100");
		assertNull(trimmer.processRecord(rec));
	}
}