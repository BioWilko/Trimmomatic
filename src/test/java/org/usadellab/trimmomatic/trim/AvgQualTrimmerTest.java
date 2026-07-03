package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class AvgQualTrimmerTest {

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
	public void testAvgQualKeep() {
		// 150bp, Q30 ('?'), Threshold 20
		AvgQualTrimmer trimmer = new AvgQualTrimmer(20);
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", makeQuality('?', 150), 33);
		assertNotNull(trimmer.processRecord(rec));
	}

	@Test
	public void testAvgQualDrop() {
		// 150bp, Q10 ('+'), Threshold 20
		AvgQualTrimmer trimmer = new AvgQualTrimmer(20);
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", makeQuality('+', 150), 33);
		assertNull(trimmer.processRecord(rec));
	}

	@Test
	public void testAvgQualBoundaryExact() {
		// 150bp, Q20 ('5'), Threshold 20. Total = 20*150. Required < 20*150 is false. Keep.
		AvgQualTrimmer trimmer = new AvgQualTrimmer(20);
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", makeQuality('5', 150), 33);
		assertNotNull(trimmer.processRecord(rec));
	}

	@Test
	public void testAvgQualBoundaryDrop() {
		// 150bp. 149 bases Q20 ('5'), 1 base Q19 ('4'). Avg < 20. Drop.
		AvgQualTrimmer trimmer = new AvgQualTrimmer(20);
		String qual = makeQuality('5', 149) + "4";
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", qual, 33);
		assertNull(trimmer.processRecord(rec));
	}

	@Test
	public void testEmptyRead() {
		// 0bp. Total 0. 0 < 0 is false. Keep.
		AvgQualTrimmer trimmer = new AvgQualTrimmer(20);
		FastqRecord rec = new FastqRecord("head", "", "", "", 33);
		assertNotNull(trimmer.processRecord(rec));
	}
}