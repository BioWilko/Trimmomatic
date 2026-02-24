package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

public class IlluminaBitwiseTest {

	private String makeSequence(int len) {
		StringBuilder sb = new StringBuilder(len);
		String pattern = "ACGT";
		for (int i = 0; i < len; i++)
			sb.append(pattern.charAt(i % 4));
		return sb.toString();
	}

	@Test
	public void testPackUnpack150bp() {
		String seq = makeSequence(150);
		long[] packed = IlluminaClippingTrimmer.packSeqExternal(seq);
		
		assertEquals(150, packed.length);
		
		// Check first 16 bases (index 0)
		String first16 = IlluminaClippingTrimmer.unpack(packed[0]);
		assertEquals(seq.substring(0, 16), first16);
		
		// Check middle 16 bases (index 50 -> bases 50..65)
		String mid16 = IlluminaClippingTrimmer.unpack(packed[50]);
		assertEquals(seq.substring(50, 66), mid16);
		
		// Check last full 16-mer (index 134 -> bases 134..149)
		String last16 = IlluminaClippingTrimmer.unpack(packed[134]);
		assertEquals(seq.substring(134, 150), last16);
	}
	
	@Test
	public void testPackPadding() {
		// Sequence "A" (length 1). packSeqExternal returns array of length 1.
		// out[0] contains 'A' followed by 15 0s (which unpack to [0])
		long[] packed = IlluminaClippingTrimmer.packSeqExternal("A");
		String unpacked = IlluminaClippingTrimmer.unpack(packed[0]);
		
		assertTrue(unpacked.startsWith("A"));
		assertTrue(unpacked.endsWith("[0]"));
	}
	
	@Test
	public void testCalcSingleMask() {
		// Length 16 -> all 1s (full mask)
		assertEquals(0xFFFFFFFFFFFFFFFFL, IlluminaClippingTrimmer.calcSingleMask(16));
		
		// Length 1 -> top 4 bits set (0xF000...)
		// 16-1 = 15. 15*4 = 60 shift.
		long expected = 0xFFFFFFFFFFFFFFFFL << 60;
		assertEquals(expected, IlluminaClippingTrimmer.calcSingleMask(1));
	}

	@Test
	public void testPackSeqInternalReverse() {
		// Sequence: "ACGT..."
		// Reverse Complement: "ACGT" -> "ACGT" (A->T, C->G, G->C, T->A)
		// Wait: RC of "ACGT" is "ACGT".
		// A -> T, C -> G, G -> C, T -> A.
		// Rev: T G C A.
		// "ACGT" -> T G C A.
		
		// Let's use "AAAA" -> "TTTT"
		String seq = "AAAAAAAAAAAAAAAA"; // 16 As
		long[] packed = IlluminaClippingTrimmer.packSeqInternal(seq, true);
		
		// Should be 1 long (16 bases)
		// 16 - 15 = 1 entry in output array for the sliding window start at 0?
		// packSeqInternal returns array of size seq.length() - 15.
		assertEquals(1, packed.length);
		
		String unpacked = IlluminaClippingTrimmer.unpack(packed[0]);
		assertEquals("TTTTTTTTTTTTTTTT", unpacked);
	}
}