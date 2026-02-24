package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.mock;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.UUID;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.util.Logger;

public class IlluminaClippingTrimmerLoadTest {

	@Test
	public void testLoadFromAbsoluteFile() throws IOException {
		// Create a dummy adapter file
		// Cannot use absolute path on Windows due to ':' separator in args. Use CWD.
		File adapterFile = new File("adapters_" + UUID.randomUUID() + ".fa");
		try (FileWriter writer = new FileWriter(adapterFile)) {
			writer.write(">PrefixPE/1\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n");
			writer.write(">PrefixPE/2\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n");
		}

		Logger logger = mock(Logger.class);
		
		// Args: file:seed:palindrome:simple:min:keep
		String args = adapterFile.getName() + ":2:30:10:8:true";
		
		try {
			IlluminaClippingTrimmer trimmer = IlluminaClippingTrimmer.makeIlluminaClippingTrimmer(logger, args);
			assertNotNull(trimmer);
		} finally {
			if (!adapterFile.delete()) {
				System.gc();
				if (!adapterFile.delete())
					adapterFile.deleteOnExit();
			}
		}
	}

	@Test
	public void testLoadMissingFile() {
		Logger logger = mock(Logger.class);
		// Use a non-existent path
		String args = "/path/to/nonexistent/file.fa:2:30:10";
		
		// Should throw IOException wrapped in RuntimeException or directly IOException depending on caller,
		// but makeIlluminaClippingTrimmer declares throws IOException.
		// However, the logic inside might catch and log, or throw. 
		// Looking at code: "if (fastaFile == null || !fastaFile.exists()) throw new IOException..."
		
		assertThrows(IOException.class, () -> {
			IlluminaClippingTrimmer.makeIlluminaClippingTrimmer(logger, args);
		});
	}

	@Test
	public void testLoadEmptyFile() throws IOException {
		File emptyFile = new File("empty_" + UUID.randomUUID() + ".fa");
		emptyFile.createNewFile();
		
		Logger logger = mock(Logger.class);
		String args = emptyFile.getName() + ":2:30:10:8:true";
		
		try {
			IlluminaClippingTrimmer trimmer = IlluminaClippingTrimmer.makeIlluminaClippingTrimmer(logger, args);
			assertNotNull(trimmer);
			
			// Should process records without error (and do nothing)
			FastqRecord rec = new FastqRecord("head", "ACGT", "", "IIII", 33);
			FastqRecord[] result = trimmer.processRecords(new FastqRecord[] { rec });
			
			assertEquals(1, result.length);
			assertEquals("ACGT", result[0].getSequence());
		} finally {
			if (!emptyFile.delete()) {
				System.gc();
				if (!emptyFile.delete())
					emptyFile.deleteOnExit();
			}
		}
	}

	@Test
	public void testLoadWindowsLineEndings() throws IOException {
		File adapterFile = new File("windows_" + UUID.randomUUID() + ".fa");
		try (FileWriter writer = new FileWriter(adapterFile)) {
			// \r\n line endings
			writer.write(">seq1\r\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCA\r\n");
		}

		Logger logger = mock(Logger.class);
		String args = adapterFile.getName() + ":2:30:10:8:true";
		
		try {
			IlluminaClippingTrimmer trimmer = IlluminaClippingTrimmer.makeIlluminaClippingTrimmer(logger, args);
			assertNotNull(trimmer);
		} finally {
			if (!adapterFile.delete()) {
				System.gc();
				if (!adapterFile.delete())
					adapterFile.deleteOnExit();
			}
		}
	}
}