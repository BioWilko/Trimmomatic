package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.trim.HeadCropTrimmer;
import org.usadellab.trimmomatic.trim.LeadingTrimmer;

public class FastqRecordHeadPosTest {

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
	public void testHeadPosAccumulation() {
		// 150bp read
		String seq = makeSequence(150);
		// First 10 bases low quality, rest high
		String qual = makeQuality('!', 10) + makeQuality('I', 140);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		assertEquals(0, rec.getHeadPos());
		
		// 1. HeadCrop 5. headPos -> 5.
		HeadCropTrimmer headCrop = new HeadCropTrimmer("5");
		FastqRecord rec1 = headCrop.processRecords(new FastqRecord[] { rec })[0];
		assertEquals(5, rec1.getHeadPos());
		
		// 2. Leading 20. The first 5 bases of rec1 are still low quality (indices 5-9 of original).
		// LeadingTrimmer should trim those 5. headPos -> 5 + 5 = 10.
		LeadingTrimmer leading = new LeadingTrimmer(20);
		FastqRecord rec2 = leading.processRecord(rec1);
		assertEquals(10, rec2.getHeadPos());
	}
}