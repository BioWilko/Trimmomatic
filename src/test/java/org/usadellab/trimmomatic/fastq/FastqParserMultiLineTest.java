package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

public class FastqParserMultiLineTest {

	@Test
	public void testMultiLineSequenceFail(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "multiline.fastq");
		try (FileWriter writer = new FileWriter(fastqFile)) {
			// Multi-line sequence is valid FASTQ but Trimmomatic parser expects 4 lines per record
			writer.write("@seq1\nACGT\nACGT\n+\nIIII\nIIII\n");
		}

		FastqParser parser = new FastqParser(33);
		
		// Should fail because it reads "ACGT" (line 3) as the comment line, 
		// and it doesn't start with '+'
		try {
			assertThrows(RuntimeException.class, () -> {
				parser.open(fastqFile);
			});
		} finally {
			parser.close();
		}
	}
}