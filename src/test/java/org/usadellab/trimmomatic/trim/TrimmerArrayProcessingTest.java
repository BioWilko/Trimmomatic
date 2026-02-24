package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;

public class TrimmerArrayProcessingTest {

	// Anonymous subclass to test the base class logic
	private Trimmer identityTrimmer = new AbstractSingleRecordTrimmer() {
		@Override
		public FastqRecord processRecord(FastqRecord in) {
			return in; // Pass through
		}
	};

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

	private FastqRecord makeRecord() {
		// 150bp read
		return new FastqRecord("head", makeSequence(150), "", makeQuality('I', 150), 33);
	}

	@Test
	public void testProcessNullArray() {
		assertNull(identityTrimmer.processRecords(null));
	}

	@Test
	public void testProcessEmptyArray() {
		FastqRecord[] input = new FastqRecord[0];
		FastqRecord[] output = identityTrimmer.processRecords(input);
		assertNotNull(output);
		assertEquals(0, output.length);
	}

	@Test
	public void testProcessArrayWithNulls() {
		// Input: [Record, null]
		FastqRecord[] input = new FastqRecord[] { makeRecord(), null };
		
		FastqRecord[] output = identityTrimmer.processRecords(input);
		
		assertNotNull(output);
		assertEquals(2, output.length);
		assertNotNull(output[0]);
		assertNull(output[1]);
	}
	
	@Test
	public void testProcessArrayNormal() {
		FastqRecord[] input = new FastqRecord[] { makeRecord(), makeRecord() };
		FastqRecord[] output = identityTrimmer.processRecords(input);
		assertEquals(2, output.length);
	}
}