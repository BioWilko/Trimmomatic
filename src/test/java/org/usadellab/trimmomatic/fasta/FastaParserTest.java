package org.usadellab.trimmomatic.fasta;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

public class FastaParserTest {

	private String makeSequence(int len) {
		StringBuilder sb = new StringBuilder(len);
		String pattern = "ACGT";
		for (int i = 0; i < len; i++)
			sb.append(pattern.charAt(i % 4));
		return sb.toString();
	}

	@Test
	public void testParseSingleLine(@TempDir File tempDir) throws IOException {
		File fastaFile = new File(tempDir, "test_single.fa");
		String seq150 = makeSequence(150);
		
		try (FileWriter writer = new FileWriter(fastaFile)) {
			writer.write(">seq1\n" + seq150 + "\n");
		}

		FastaParser parser = new FastaParser();
		try {
			parser.parse(fastaFile);
	
			assertTrue(parser.hasNext());
			FastaRecord rec1 = parser.next();
			assertEquals("seq1", rec1.getName());
			assertEquals(seq150, rec1.getSequence());
			
			assertFalse(parser.hasNext());
		} finally {
			parser.close();
		}
	}

	@Test
	public void testParseMultiLine(@TempDir File tempDir) throws IOException {
		File fastaFile = new File(tempDir, "test_multi.fa");
		String seq150 = makeSequence(150);
		// Split into 50bp lines
		String part1 = seq150.substring(0, 50);
		String part2 = seq150.substring(50, 100);
		String part3 = seq150.substring(100, 150);
		
		try (FileWriter writer = new FileWriter(fastaFile)) {
			writer.write(">seq1\n" + part1 + "\n" + part2 + "\n" + part3 + "\n");
		}

		FastaParser parser = new FastaParser();
		try {
			parser.parse(fastaFile);
	
			assertTrue(parser.hasNext());
			FastaRecord rec1 = parser.next();
			assertEquals("seq1", rec1.getName());
			assertEquals(seq150, rec1.getSequence());
			
			assertFalse(parser.hasNext());
		} finally {
			parser.close();
		}
	}
	
	@Test
	public void testParseMultipleRecords(@TempDir File tempDir) throws IOException {
		File fastaFile = new File(tempDir, "test_multiple.fa");
		String seqA = makeSequence(150);
		String seqB = makeSequence(150); // Same content but distinct record
		
		try (FileWriter writer = new FileWriter(fastaFile)) {
			writer.write(">recA\n" + seqA + "\n>recB\n" + seqB + "\n");
		}

		FastaParser parser = new FastaParser();
		try {
			parser.parse(fastaFile);
	
			assertTrue(parser.hasNext());
			assertEquals("recA", parser.next().getName());
			assertTrue(parser.hasNext());
			assertEquals("recB", parser.next().getName());
			assertFalse(parser.hasNext());
		} finally {
			parser.close();
		}
	}

	@Test
	public void testNonExistentFile(@TempDir File tempDir) {
		File missingFile = new File(tempDir, "missing.fa");
		FastaParser parser = new FastaParser();
		// Should throw IOException if file doesn't exist
		assertThrows(IOException.class, () -> parser.parse(missingFile));
	}
}