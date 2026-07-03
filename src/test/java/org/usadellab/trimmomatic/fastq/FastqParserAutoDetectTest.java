package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

public class FastqParserAutoDetectTest {

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
	public void testDetectPhred33(@TempDir File tempDir) throws IOException {
		// Phred33 detection looks for chars in range [33, 58] ('!' to ':')
		// and NO chars in range [80, 104] ('P' to 'h')
		
		File fastqFile = new File(tempDir, "test_phred33.fastq");
		String seq = makeSequence(150);
		// Use '!' (33) and '5' (53) which are in the Phred33 detection bin
		String qual = makeQuality('!', 75) + makeQuality('5', 75);

		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\n" + seq + "\n+\n" + qual + "\n");
		}

		FastqParser parser = new FastqParser(0); // 0 triggers auto-detect logic in open()
		parser.open(fastqFile);
		
		try {
			assertEquals(33, parser.determinePhredOffset());
		} finally {
			parser.close();
		}
	}

	@Test
	public void testDetectPhred64(@TempDir File tempDir) throws IOException {
		// Phred64 detection looks for chars in range [80, 104] ('P' to 'h')
		// and NO chars in range [33, 58] ('!' to ':')
		
		File fastqFile = new File(tempDir, "test_phred64.fastq");
		String seq = makeSequence(150);
		// Use 'P' (80) and 'h' (104) which are in the Phred64 detection bin
		String qual = makeQuality('P', 75) + makeQuality('h', 75);

		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\n" + seq + "\n+\n" + qual + "\n");
		}

		FastqParser parser = new FastqParser(0);
		parser.open(fastqFile);
		try {
			assertEquals(64, parser.determinePhredOffset());
		} finally {
			parser.close();
		}
	}

	@Test
	public void testDetectAmbiguous(@TempDir File tempDir) throws IOException {
		// If qualities fall outside both specific bins (e.g. only 'I' (73)), 
		// or mix them in a way that violates exclusivity, it returns 0.
		
		File fastqFile = new File(tempDir, "test_ambiguous.fastq");
		String seq = makeSequence(150);
		String qual = makeQuality('I', 150); // 'I' is 73, between 58 and 80.

		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\n" + seq + "\n+\n" + qual + "\n");
		}

		FastqParser parser = new FastqParser(0);
		parser.open(fastqFile);
		try {
			assertEquals(0, parser.determinePhredOffset());
		} finally {
			parser.close();
		}
	}

	@Test
	public void testDetectMixedBins(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "test_mixed.fastq");
		String seq = makeSequence(150);
		// '!' (33) from Phred33 bin, 'h' (104) from Phred64 bin
		String qual = makeQuality('!', 75) + makeQuality('h', 75);

		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\n" + seq + "\n+\n" + qual + "\n");
		}

		FastqParser parser = new FastqParser(0);
		parser.open(fastqFile);
		
		try {
			assertEquals(0, parser.determinePhredOffset());
		} finally {
			parser.close();
		}
	}
}