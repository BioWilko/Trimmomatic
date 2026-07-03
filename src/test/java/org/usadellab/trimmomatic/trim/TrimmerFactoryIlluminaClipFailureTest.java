package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.mock;

import java.io.File;
import java.io.IOException;
import java.util.UUID;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.util.Logger;

public class TrimmerFactoryIlluminaClipFailureTest {

	@Test
	public void testMissingArguments() throws IOException {
		// Use CWD to avoid Windows absolute path colon issues
		File adapterFile = new File("adapters_missing_" + UUID.randomUUID() + ".fa");
		adapterFile.createNewFile();
		
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		
		// Missing simple clip threshold
		// ILLUMINACLIP:file:seed:palindrome
		String args = "ILLUMINACLIP:" + adapterFile.getName() + ":2:30";
		
		try {
			assertThrows(ArrayIndexOutOfBoundsException.class, () -> factory.makeTrimmer(args));
		} finally {
			if (!adapterFile.delete()) {
				System.gc();
				if (!adapterFile.delete())
					adapterFile.deleteOnExit();
			}
		}
	}

	@Test
	public void testInvalidIntegerArgs() throws IOException {
		File adapterFile = new File("adapters_invalid_" + UUID.randomUUID() + ".fa");
		adapterFile.createNewFile();
		
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		
		// Invalid seed mismatches (float instead of int)
		String args = "ILLUMINACLIP:" + adapterFile.getName() + ":2.5:30:10";
		
		try {
			assertThrows(NumberFormatException.class, () -> factory.makeTrimmer(args));
		} finally {
			if (!adapterFile.delete()) {
				System.gc();
				if (!adapterFile.delete())
					adapterFile.deleteOnExit();
			}
		}
	}

	@Test
	public void testMissingAdapterFile() {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		String args = "ILLUMINACLIP:nonexistent_file_" + UUID.randomUUID() + ".fa:2:30:10";
		assertThrows(IOException.class, () -> factory.makeTrimmer(args));
	}
}