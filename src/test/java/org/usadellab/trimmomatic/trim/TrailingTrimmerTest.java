package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class TrailingTrimmerTest {

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
	public void testTrailingNoTrim() {
		// Threshold 20. All Q30 ('?').
		TrailingTrimmer trimmer = new TrailingTrimmer(20);
		String seq = makeSequence(150);
		String qual = makeQuality('?', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		FastqRecord result = trimmer.processRecord(rec);
		assertNotNull(result);
		assertEquals(150, result.getSequence().length());
		assertEquals(seq, result.getSequence());
	}

	@Test
	public void testTrailingTrimEnd() {
		// Threshold 20. First 100 Q30 ('?'), last 50 Q10 ('+').
		TrailingTrimmer trimmer = new TrailingTrimmer(20);
		String seq = makeSequence(150);
		String qual = makeQuality('?', 100) + makeQuality('+', 50);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		FastqRecord result = trimmer.processRecord(rec);
		assertNotNull(result);
		assertEquals(100, result.getSequence().length());
		assertEquals(seq.substring(0, 100), result.getSequence());
	}

	@Test
	public void testTrailingDropAllLow() {
		// Threshold 20. All Q10 ('+'). Should drop.
		TrailingTrimmer trimmer = new TrailingTrimmer(20);
		String seq = makeSequence(150);
		String qual = makeQuality('+', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		assertNull(trimmer.processRecord(rec));
	}

	@Test
	public void testTrailingDropSingleBaseHigh() {
		// Threshold 20. First base Q30 ('?'), rest Q10 ('+').
		// Due to implementation (i > 0), TrailingTrimmer drops reads if only the first base survives.
		TrailingTrimmer trimmer = new TrailingTrimmer(20);
		String seq = makeSequence(150);
		String qual = "?" + makeQuality('+', 149);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		assertNull(trimmer.processRecord(rec));
	}

	@Test
	public void testTrailingKeepTwoBasesHigh() {
		// Threshold 20. First 2 bases Q30 ('?'), rest Q10 ('+').
		TrailingTrimmer trimmer = new TrailingTrimmer(20);
		String seq = makeSequence(150);
		String qual = "??" + makeQuality('+', 148);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		FastqRecord result = trimmer.processRecord(rec);
		assertNotNull(result);
		assertEquals(2, result.getSequence().length());
	}
}