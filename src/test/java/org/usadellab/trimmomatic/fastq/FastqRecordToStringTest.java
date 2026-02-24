package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertNotNull;

import org.junit.jupiter.api.Test;

public class FastqRecordToStringTest {

	@Test
	public void testToString() {
		FastqRecord rec = new FastqRecord("head", "ACGT", "comm", "IIII", 33);
		String str = rec.toString();
		assertNotNull(str);
		// FastqRecord doesn't override toString(), so it uses Object.toString()
		// Just ensuring it doesn't crash.
	}
}