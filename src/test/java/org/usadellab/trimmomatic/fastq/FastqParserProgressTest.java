package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

public class FastqParserProgressTest {

	private String makeSequence(int len) {
		StringBuilder sb = new StringBuilder(len);
		String pattern = "ACGT";
		for (int i = 0; i < len; i++)
			sb.append(pattern.charAt(i % 4));
		return sb.toString();
	}

	@Test
	public void testProgressUpdates(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "progress.fastq");
		String seq = makeSequence(150);
		String qual = makeSequence(150); // Dummy quality
		
		// Write 10000 records to ensure file size > buffer size (32KB)
		try (FileWriter writer = new FileWriter(fastqFile)) {
			for (int i = 0; i < 10000; i++) {
				writer.write("@seq" + i + "\n" + seq + "\n+\n" + qual + "\n");
			}
		}

		FastqParser parser = new FastqParser(33);
		parser.open(fastqFile);

		try {
			int initialProgress = parser.getProgress();
			assertTrue(initialProgress >= 0 && initialProgress <= 100);
	
			// Read half
			for (int i = 0; i < 5000; i++) {
				parser.next();
			}
			
			int midProgress = parser.getProgress();
			assertTrue(midProgress > initialProgress, "Progress should increase");
			assertTrue(midProgress <= 100);
		} finally {
			parser.close();
		}
	}
}