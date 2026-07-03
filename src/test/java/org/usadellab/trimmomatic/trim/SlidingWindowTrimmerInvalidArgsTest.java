package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class SlidingWindowTrimmerInvalidArgsTest {

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
	public void testWindowSizeZero() {
		// Window size 0.
		// Logic analysis: totalRequiredQuality becomes 0.
		// The sliding window loop condition (total < totalRequired) is never met (0 < 0 is false).
		// So it keeps the full length initially.
		// Then it enters the backtracking loop: while(lastBase < required).
		// This effectively makes it behave like a TrailingTrimmer!
		
		SlidingWindowTrimmer trimmer = new SlidingWindowTrimmer(0, 20);
		
		// 150bp read: 100bp High Qual, 50bp Low Qual
		String seq = makeSequence(150);
		String qual = makeQuality('I', 100) + makeQuality('!', 50);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		
		FastqRecord result = trimmer.processRecord(rec);
		assertNotNull(result);
		// Should trim the 50bp low quality tail
		assertEquals(100, result.getSequence().length());
	}
}