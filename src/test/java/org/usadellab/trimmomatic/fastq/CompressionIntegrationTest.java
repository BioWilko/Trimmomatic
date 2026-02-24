package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

public class CompressionIntegrationTest {

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
	public void testGzipRoundTrip(@TempDir File tempDir) throws IOException {
		File gzFile = new File(tempDir, "test.fastq.gz");
		
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		// Write
		FastqSerializer serializer = new FastqSerializer();
		serializer.open(gzFile);
		try {
			serializer.writeRecord(rec);
		} finally {
			serializer.close();
		}

		// Verify file exists and has content
		assertTrue(gzFile.exists());
		assertTrue(gzFile.length() > 0);

		// Read
		FastqParser parser = new FastqParser(33);
		parser.open(gzFile);
		
		try {
			assertTrue(parser.hasNext());
			FastqRecord readRec = parser.next();
			
			assertEquals(rec.getName(), readRec.getName());
			assertEquals(rec.getSequence(), readRec.getSequence());
			assertEquals(rec.getQuality(), readRec.getQuality());
		} finally {
			parser.close();
		}
	}

	@Test
	public void testBzip2RoundTrip(@TempDir File tempDir) throws IOException {
		File bz2File = new File(tempDir, "test.fastq.bz2");
		
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150);
		FastqRecord rec = new FastqRecord("head", seq, "", qual, 33);

		// Write
		FastqSerializer serializer = new FastqSerializer();
		serializer.open(bz2File);
		try {
			serializer.writeRecord(rec);
		} finally {
			serializer.close();
		}

		// Read
		FastqParser parser = new FastqParser(33);
		parser.open(bz2File);
		
		try {
			assertTrue(parser.hasNext());
			FastqRecord readRec = parser.next();
			
			assertEquals(rec.getName(), readRec.getName());
			assertEquals(rec.getSequence(), readRec.getSequence());
		} finally {
			parser.close();
		}
	}
}