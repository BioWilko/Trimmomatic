package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

import org.junit.jupiter.api.Test;

public class FastqSerializerTest {

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
	public void testWriteRecord() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		FastqSerializer serializer = new FastqSerializer();
		serializer.open(baos);

		// 150bp read
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);
		FastqRecord rec = new FastqRecord("seq1", seq, "comment", qual, 33);

		serializer.writeRecord(rec);
		serializer.close();

		String output = baos.toString();
		String expected = "@seq1\n" + seq + "\n+comment\n" + qual + "\n";
		
		assertEquals(expected, output);
	}

	@Test
	public void testWriteMultipleRecords() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		FastqSerializer serializer = new FastqSerializer();
		serializer.open(baos);

		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);
		
		FastqRecord rec1 = new FastqRecord("seq1", seq, "comment1", qual, 33);
		FastqRecord rec2 = new FastqRecord("seq2", seq, "comment2", qual, 33);

		serializer.writeRecord(rec1);
		serializer.writeRecord(rec2);
		serializer.close();

		String output = baos.toString();
		String expected = "@seq1\n" + seq + "\n+comment1\n" + qual + "\n" +
				          "@seq2\n" + seq + "\n+comment2\n" + qual + "\n";
		
		assertEquals(expected, output);
	}
}