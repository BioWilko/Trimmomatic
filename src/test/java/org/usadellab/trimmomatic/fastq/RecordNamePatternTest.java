package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;

import org.junit.jupiter.api.Test;

public class RecordNamePatternTest {

	@Test
	public void testCasava13() {
		// Format: @HWUSI-EAS100R:6:73:941:1973#0/1
		// Should extract: HWUSI-EAS100R:6:73:941:1973#0
		String header = "HWUSI-EAS100R:6:73:941:1973#0/1";
		assertEquals("HWUSI-EAS100R:6:73:941:1973#0", RecordNamePattern.CASAVA_13.canonicalizeOne(header));
	}

	@Test
	public void testCasava18() {
		// Format: @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
		// Should extract: EAS139:136:FC706VJ:2:2104:15343:197393
		String header = "EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG";
		assertEquals("EAS139:136:FC706VJ:2:2104:15343:197393", RecordNamePattern.CASAVA_18.canonicalizeOne(header));
	}

	@Test
	public void testNoMatch() {
		String header = "SimpleName";
		assertNull(RecordNamePattern.CASAVA_13.canonicalizeOne(header));
	}
}