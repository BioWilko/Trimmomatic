package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class MaximumInformationTrimmerTest {

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
	public void testMaxInfoPerfect() {
		// Target 150, Strictness 0.5. Perfect read.
		MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("150:0.5");
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150); // High quality Q40
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		FastqRecord result = trimmer.processRecord(rec);
		assertNotNull(result);
		assertEquals(150, result.getSequence().length());
	}

	@Test
	public void testMaxInfoTrimTail() {
		// Target 150, Strictness 0.8.
		// 100bp High Qual, 50bp Low Qual.
		MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("150:0.8");
		String seq = makeSequence(150);
		String qual = makeQuality('I', 100) + makeQuality('!', 50); // 100 Q40, 50 Q0
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		FastqRecord result = trimmer.processRecord(rec);
		assertNotNull(result);
		// Should trim off the low quality tail.
		// Exact length depends on the scoring formula, but should be around 100.
		assertTrue(result.getSequence().length() < 140);
		assertTrue(result.getSequence().length() >= 90);
	}

	@Test
	public void testMaxInfoDropGarbage() {
		// Target 150, Strictness 0.9.
		// All Low Qual.
		MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("150:0.9");
		String seq = makeSequence(150);
		String qual = makeQuality('!', 150); // Q0
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		// The trimmer finds the maximum information point. Even if quality is low, 
		// it returns the best cut point found.
		FastqRecord result = trimmer.processRecord(rec);
		assertNotNull(result);
		assertTrue(result.getSequence().length() < 150);
	}

	@Test
	public void testMaxInfoStrictness0() {
		// Strictness 0 -> Purely length based (mostly).
		// Should keep as much as possible even if quality is low, as long as it adds "unique" info.
		// But max info formula is complex.
		// Coverage weighting = length^(1-0) = length.
		// Quality weighting = prob^0 = 1.
		// So score accumulates length score.
		
		MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("150:0");
		String seq = makeSequence(150);
		String qual = makeQuality('!', 150); // Low quality
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		FastqRecord result = trimmer.processRecord(rec);
		assertNotNull(result);
		// With strictness 0, quality doesn't penalize. It should keep the whole read if it's long enough.
		assertEquals(150, result.getSequence().length());
	}

	@Test
	public void testMaxInfoStrictness1() {
		// Strictness 1 -> Quality is paramount.
		// Coverage weighting = length^(1-1) = 1.
		// Quality weighting = prob^1 = prob.
		// If quality is low, prob is low (or negative in log space?).
		// The implementation uses log space.
		// qualProbTmp[i] = Math.log(...) * strictness.
		// If strictness is 1, low quality bases have large negative scores.
		
		MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("150:1.0");
		String seq = makeSequence(150);
		// 100bp High, 50bp Low
		String qual = makeQuality('I', 100) + makeQuality('!', 50); 
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		FastqRecord result = trimmer.processRecord(rec);
		assertNotNull(result);
		// Should trim the low quality part aggressively.
		assertTrue(result.getSequence().length() <= 100);
	}
}