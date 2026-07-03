package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

public class IlluminaBitwiseExtendedTest {

	@Test
	public void testPackN() {
		// 'N' is not A, C, G, or T.
		// packCh returns 0.
		// unpack converts 0 to "[0]".
		long[] packed = IlluminaClippingTrimmer.packSeqExternal("N");
		String unpacked = IlluminaClippingTrimmer.unpack(packed[0]);
		
		assertTrue(unpacked.startsWith("[0]"));
	}
}