package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class GeneralOvertrimTest {

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
	public void testSlidingWindowBoundary() {
		// Window 4, Threshold 20 (Phred 33 '5').
		// We want a window with average exactly 20.
		// 20 * 4 = 80 total quality.
		// Quals: 20, 20, 20, 20. ('5', '5', '5', '5')
		
		SlidingWindowTrimmer trimmer = new SlidingWindowTrimmer("4:20");
		
		String seq = makeSequence(100);
		// 50 bases High (Q30), 4 bases Q20, rest High.
		String qual = makeQuality('?', 50) + makeQuality('5', 4) + makeQuality('?', 46);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord res = trimmer.processRecord(rec);
		
		assertNotNull(res);
		// Should NOT trim. The window of Q20s has average 20. 20 < 20 is False.
		assertEquals(100, res.getSequence().length());
	}

	@Test
	public void testSlidingWindowLowBaseHighAvg() {
		// Window 4, Threshold 20.
		// Window: 10, 30, 30, 30. Sum 100. Avg 25.
		// Should NOT trim despite the Q10 base at the start of the window.
		
		SlidingWindowTrimmer trimmer = new SlidingWindowTrimmer("4:20");
		
		String seq = makeSequence(100);
		// 50 High, 1 Low (Q10='+'), 3 High (Q30='?'), rest High.
		String qual = makeQuality('?', 50) + "+" + makeQuality('?', 3) + makeQuality('?', 46);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord res = trimmer.processRecord(rec);
		
		assertNotNull(res);
		assertEquals(100, res.getSequence().length());
	}

	@Test
	public void testLeadingBoundary() {
		// Leading 20.
		// First base Q20 ('5'). Should keep.
		LeadingTrimmer trimmer = new LeadingTrimmer(20);
		
		String seq = makeSequence(100);
		String qual = "5" + makeQuality('?', 99);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord res = trimmer.processRecord(rec);
		
		assertEquals(100, res.getSequence().length());
	}

	@Test
	public void testTrailingBoundary() {
		// Trailing 20.
		// Last base Q20 ('5'). Should keep.
		TrailingTrimmer trimmer = new TrailingTrimmer(20);
		
		String seq = makeSequence(100);
		String qual = makeQuality('?', 99) + "5";
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord res = trimmer.processRecord(rec);
		
		assertEquals(100, res.getSequence().length());
	}

	@Test
	public void testAvgQualBoundary() {
		// AvgQual 20.
		// Read all Q20. Avg 20. Should keep.
		AvgQualTrimmer trimmer = new AvgQualTrimmer(20);
		
		String seq = makeSequence(100);
		String qual = makeQuality('5', 100);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord res = trimmer.processRecord(rec);
		
		assertNotNull(res);
	}

	@Test
	public void testMaxInfoPerfectRead() {
		// MaxInfo: Target 50, Strictness 0.5.
		// Read 100bp, Perfect Quality.
		// Should keep all 100bp because information content increases with length
		// and there is no quality penalty to drag the score down.
		
		MaximumInformationTrimmer trimmer = new MaximumInformationTrimmer("50:0.5");
		
		String seq = makeSequence(100);
		String qual = makeQuality('I', 100);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord res = trimmer.processRecord(rec);
		
		assertNotNull(res);
		assertEquals(100, res.getSequence().length());
	}
	
	@Test
	public void testHeadCropExact() {
		// HeadCrop 10.
		// Ensure it doesn't crop 11.
		HeadCropTrimmer trimmer = new HeadCropTrimmer("10");
		FastqRecord rec = new FastqRecord("head", makeSequence(100), "", makeQuality('I', 100), 33);
		FastqRecord[] res = trimmer.processRecords(new FastqRecord[]{rec});
		
		assertEquals(90, res[0].getSequence().length());
	}
	
	@Test
	public void testCropExact() {
		// Crop 50.
		// Ensure it doesn't crop to 49.
		CropTrimmer trimmer = new CropTrimmer(50);
		FastqRecord rec = new FastqRecord("head", makeSequence(100), "", makeQuality('I', 100), 33);
		FastqRecord res = trimmer.processRecord(rec);
		
		assertEquals(50, res.getSequence().length());
	}
}