package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class ToPhred33TrimmerTest {

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
	public void testConvertPhred64To33() {
		// Input Phred64. 'h' = 104. 104 - 31 = 73 ('I', Q40 in Phred33).
		ToPhred33Trimmer trimmer = new ToPhred33Trimmer("");
		String seq = makeSequence(150);
		String qual64 = makeQuality('h', 150);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual64, 64);
		FastqRecord result = trimmer.processRecord(rec);

		assertNotNull(result);
		assertEquals(33, result.getPhredOffset());
		assertEquals(makeQuality('I', 150), result.getQuality());
	}

	@Test
	public void testAlreadyPhred33() {
		// Input Phred33. Should remain unchanged.
		ToPhred33Trimmer trimmer = new ToPhred33Trimmer("");
		String seq = makeSequence(150);
		String qual33 = makeQuality('I', 150);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual33, 33);
		FastqRecord result = trimmer.processRecord(rec);

		assertEquals(rec, result);
	}
}