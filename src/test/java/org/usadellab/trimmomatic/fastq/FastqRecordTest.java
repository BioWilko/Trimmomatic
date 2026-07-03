package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import org.junit.jupiter.api.Test;

public class FastqRecordTest {

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
	public void testPhred33QualityConversion() {
		// 150bp read
		String seq = makeSequence(150);
		// Mix of low ('!'=0) and high ('I'=40) quality
		String qual = makeQuality('!', 75) + makeQuality('I', 75);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		int[] quals = rec.getQualityAsInteger(true);
		
		assertEquals(150, quals.length);
		for(int i=0; i<75; i++) assertEquals(0, quals[i]);
		for(int i=75; i<150; i++) assertEquals(40, quals[i]);
	}

	@Test
	public void testPhred64QualityConversion() {
		// 150bp read
		String seq = makeSequence(150);
		// Mix of low ('@'=0) and high ('h'=40) quality for Phred64
		String qual = makeQuality('@', 75) + makeQuality('h', 75);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 64);
		int[] quals = rec.getQualityAsInteger(true);
		
		assertEquals(150, quals.length);
		for(int i=0; i<75; i++) assertEquals(0, quals[i]);
		for(int i=75; i<150; i++) assertEquals(40, quals[i]);
	}

	@Test
	public void testSubRecordCreation() {
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);
		FastqRecord original = new FastqRecord("head", seq, "", qual, 33);
		
		// Slice middle 50bp (offset 50, length 50)
		FastqRecord sub = new FastqRecord(original, 50, 50);
		
		assertEquals("head", sub.getName());
		assertEquals(50, sub.getSequence().length());
		assertEquals(seq.substring(50, 100), sub.getSequence());
		assertEquals(qual.substring(50, 100), sub.getQuality());
		assertEquals(33, sub.getPhredOffset());
	}
	
	@Test
	public void testSubRecordOutOfBounds() {
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);
		FastqRecord original = new FastqRecord("head", seq, "", qual, 33);
		
		// Try to slice starting completely beyond end: Offset 200 > 150
		// This should fail with StringIndexOutOfBoundsException
		assertThrows(StringIndexOutOfBoundsException.class, () -> {
			new FastqRecord(original, 200, 50); 
		});
	}

	@Test
	public void testEquals() {
		String seq = makeSequence(150);
		String qual1 = makeQuality('I', 150);
		String qual2 = makeQuality('J', 150);

		FastqRecord rec1 = new FastqRecord("head", seq, "", qual1, 33);
		FastqRecord rec2 = new FastqRecord("head", seq, "", qual1, 33);
		FastqRecord rec3 = new FastqRecord("head", seq, "", qual2, 33);

		assertEquals(rec1, rec2);
		assertNotEquals(rec1, rec3);
	}

	@Test
	public void testConstructorLengthMismatch() {
		String seq = makeSequence(150);
		String qual = makeQuality('I', 149); // One base short

		assertThrows(RuntimeException.class, () -> {
			new FastqRecord("head", seq, "", qual, 33);
		});
	}

	@Test
	public void testQualityWithNs() {
		// 150bp read with Ns
		String seq = makeSequence(140) + "NNNNNNNNNN";
		String qual = makeQuality('I', 150);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		int[] quals = rec.getQualityAsInteger(true);
		
		assertEquals(150, quals.length);
		for(int i=0; i<140; i++) assertEquals(40, quals[i]);
		for(int i=140; i<150; i++) assertEquals(0, quals[i]);
	}

	@Test
	public void testNestedSlicing() {
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);
		FastqRecord original = new FastqRecord("head", seq, "", qual, 33);
		
		FastqRecord slice1 = new FastqRecord(original, 10, 140);
		assertEquals(10, slice1.getHeadPos());
		
		FastqRecord slice2 = new FastqRecord(slice1, 10, 130);
		assertEquals(20, slice2.getHeadPos());
		assertEquals(130, slice2.getSequence().length());
	}

	@Test
	public void testTrimNegativeHead() {
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		
		assertThrows(RuntimeException.class, () -> {
			new FastqRecord(rec, -1, 100);
		});
	}

	@Test
	public void testTrimClamping() {
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		
		// Request 200 bases from offset 0, should clamp to 150
		FastqRecord slice = new FastqRecord(rec, 0, 200);
		assertEquals(150, slice.getSequence().length());
	}

	@Test
	public void testBarcodeLabelPreservation() {
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);
		FastqRecord original = new FastqRecord("head", seq, "", qual, 33);
		original.setBarcodeLabel("BC1");
		
		FastqRecord slice = new FastqRecord(original, 10, 100);
		assertEquals("BC1", slice.getBarcodeLabel());
	}

	@Test
	public void testStaticMakeWithQuality() {
		// Create record with integer quality 30 (Phred+33 = 63 = '?')
		FastqRecord rec = FastqRecord.make("seq1", "ACGT", 30);
		assertEquals("seq1", rec.getName());
		assertEquals("ACGT", rec.getSequence());
		assertEquals("????", rec.getQuality());
		assertEquals(33, rec.getPhredOffset());
	}

	@Test
	public void testStaticMakeDefaultQuality() {
		// Default quality is 40 (Phred+33 = 73 = 'I')
		FastqRecord rec = FastqRecord.make("seq1", "ACGT");
		assertEquals("ACGT", rec.getSequence());
		assertEquals("IIII", rec.getQuality());
	}

	@Test
	public void testQualityAsIntegerNoZeroNs() {
		// 150bp read with Ns
		String seq = makeSequence(140) + "NNNNNNNNNN";
		// Quality 'I' (40) everywhere
		String qual = makeQuality('I', 150);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		// zeroNs = false -> should return actual quality even for Ns
		int[] quals = rec.getQualityAsInteger(false);
		
		assertEquals(150, quals.length);
		for(int i=0; i<150; i++) assertEquals(40, quals[i]);
	}

	@Test
	public void testQualityAsIntegerZeroNs() {
		// 150bp read with Ns
		String seq = makeSequence(140) + "NNNNNNNNNN";
		// Quality 'I' (40) everywhere
		String qual = makeQuality('I', 150);
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		// zeroNs = true -> Ns should be 0
		int[] quals = rec.getQualityAsInteger(true);
		
		assertEquals(150, quals.length);
		for(int i=0; i<140; i++) assertEquals(40, quals[i]);
		for(int i=140; i<150; i++) assertEquals(0, quals[i]);
	}

	@Test
	public void testGetRecordLength() {
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);
		String name = "header";
		String comment = "comment";
		
		FastqRecord rec = new FastqRecord(name, seq, comment, qual, 33);
		
		// Length = name.len + seq.len + comment.len + qual.len + 10 (constants)
		// 6 + 150 + 7 + 150 + 10 = 323
		assertEquals(323, rec.getRecordLength());
	}

	@Test
	public void testSetPhredOffset() {
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150); // 'I' = 73.
		
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);
		// At 33, 'I' is 40.
		assertEquals(40, rec.getQualityAsInteger(true)[0]);
		
		rec.setPhredOffset(64);
		// At 64, 'I' (73) is 9.
		assertEquals(9, rec.getQualityAsInteger(true)[0]);
	}

	@Test
	public void testConstructorFromBase() {
		FastqRecord base = new FastqRecord("name", "ACGT", "comment", "IIII", 33);
		base.setBarcodeLabel("BC1");
		
		String newSeq = "AC";
		String newQual = "II";
		
		FastqRecord newRec = new FastqRecord(base, newSeq, newQual, 64);
		
		assertEquals("name", newRec.getName());
		assertEquals("comment", newRec.getComment());
		assertEquals("AC", newRec.getSequence());
		assertEquals("II", newRec.getQuality());
		assertEquals(64, newRec.getPhredOffset());
		assertEquals("BC1", newRec.getBarcodeLabel());
	}

	@Test
	public void testHashCode() {
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);
		
		FastqRecord rec1 = new FastqRecord("head", seq, "", qual, 33);
		FastqRecord rec2 = new FastqRecord("head", seq, "", qual, 33);
		
		assertEquals(rec1.hashCode(), rec2.hashCode());
	}

	@Test
	public void testConstructorNulls() {
		// Should throw NPE if sequence or quality is null because length() is called
		assertThrows(NullPointerException.class, () -> {
			new FastqRecord("head", null, "", "qual", 33);
		});
	}

	@Test
	public void testPhredOffsetIndependence() {
		FastqRecord base = new FastqRecord("head", "ACGT", "", "IIII", 33);
		FastqRecord sub = new FastqRecord(base, 0, 4);
		
		// Change offset of base
		base.setPhredOffset(64);
		
		// Sub should still be 33
		assertEquals(33, sub.getPhredOffset());
		// Base should be 64
		assertEquals(64, base.getPhredOffset());
	}
}