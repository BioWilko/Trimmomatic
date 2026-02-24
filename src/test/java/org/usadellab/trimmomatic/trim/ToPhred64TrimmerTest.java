package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class ToPhred64TrimmerTest {

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
	public void testConvertPhred33To64() {
		// Input Phred33. 'I' = 73. 73 + 31 = 104 ('h', Q40 in Phred64).
		ToPhred64Trimmer trimmer = new ToPhred64Trimmer("");
		String seq = makeSequence(150);
		String qual33 = makeQuality('I', 150);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual33, 33);
		FastqRecord result = trimmer.processRecord(rec);

		assertNotNull(result);
		assertEquals(64, result.getPhredOffset());
		assertEquals(makeQuality('h', 150), result.getQuality());
	}

	@Test
	public void testAlreadyPhred64() {
		// Input Phred64. Should remain unchanged.
		ToPhred64Trimmer trimmer = new ToPhred64Trimmer("");
		String seq = makeSequence(150);
		String qual64 = makeQuality('h', 150);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual64, 64);
		FastqRecord result = trimmer.processRecord(rec);

		assertEquals(rec, result);
	}
}