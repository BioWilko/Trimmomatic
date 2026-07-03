package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNull;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

public class FastqParserIteratorTest {

	@Test
	public void testNextAtEOF(@TempDir File tempDir) throws IOException {
		File fastqFile = new File(tempDir, "single.fastq");
		try (FileWriter writer = new FileWriter(fastqFile)) {
			writer.write("@seq1\nACGT\n+\nIIII\n");
		}

		FastqParser parser = new FastqParser(33);
		parser.open(fastqFile);

		try {
			// 1. hasNext should be true
			// 2. next should return record
			if (parser.hasNext()) {
				parser.next();
			}
			
			// 3. hasNext should be false
			assertFalse(parser.hasNext());
			
			// 4. next should return null (based on implementation: return record; nextRecord = parseOne();)
			// When hasNext is false, nextRecord is null.
			// Calling next() returns null and tries to parse again (which returns null).
			assertNull(parser.next());
		} finally {
			parser.close();
		}
	}
}