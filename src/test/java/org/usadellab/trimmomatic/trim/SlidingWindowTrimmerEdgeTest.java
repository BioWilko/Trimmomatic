package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class SlidingWindowTrimmerEdgeTest {

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
	public void testWindowLargerThanRead() {
		// Window 151, Read 150.
		SlidingWindowTrimmer trimmer = new SlidingWindowTrimmer("151:20");
		
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150); // High quality
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		
		// Should drop because window cannot fit.
		assertNull(trimmer.processRecord(rec));
	}

	@Test
	public void testWindowEqualsRead() {
		// Window 150, Read 150.
		SlidingWindowTrimmer trimmer = new SlidingWindowTrimmer("150:20");
		
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150); // Avg 40 > 20.
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		
		// Should keep.
		assertNotNull(trimmer.processRecord(rec));
	}
	
	@Test
	public void testWindowEqualsReadFail() {
		// Window 150, Read 150.
		SlidingWindowTrimmer trimmer = new SlidingWindowTrimmer("150:20");
		FastqRecord rec = new FastqRecord("head", makeSequence(150), "", makeQuality('!', 150), 33);
		// Should drop.
		assertNull(trimmer.processRecord(rec));
	}
}